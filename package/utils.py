#!/usr/bin/env python
# -*- coding: utf-8 -*-
#==========================================================================
# 
#--------------------------------------------------------------------------
# Copyright (c) 2022 Light Conversion, UAB
# All rights reserved.
# www.lightcon.com
#==========================================================================
from cProfile import label
import matplotlib.pyplot as plt
from scipy.optimize import minimize
import numpy as np
import sys
import json


def get_wavelength_calibration(wlc, wlc_sample, descriptions):
    wlc_sample_norm = [None, None]

    wl = [wlc[0][:,0], wlc[0][:,0]]

    plt.figure('Calibration report')
    plt.clf()
    axes = [plt.subplot(211), plt.subplot(212)]

    def get_shifted_sample(par):
        wl_new = wl[0] + np.poly1d(par)(np.arange(len(wl[0])))
        return (wl_new, np.interp(wl_new, wl[0], wlc_sample_norm[1]))

    def fun_min(par):
        wl_new, wlc_new = get_shifted_sample(par)
        wl_range = [np.max([wl[0][0], wl_new[0]]), np.min([wl[0][-1], wl_new[-1]])]
        return np.sum(np.abs([wlc_sample_norm[0][i] - wlc_new[i] for i,wl in enumerate(wl[0])]))

    for i in [0,1]:
        wlc_sample_norm[i] = wlc_sample[i][:,1]/wlc[i][:,1]

    par0 = [0.0, 0.0]
    out1 = minimize(fun_min, par0, method = 'Nelder-Mead', options = {'maxiter' : 10000, 'fatol': 1.0e-10} )

    par1 = [0.0, 0.0] + list(out1.x)
    par1[0] = 0.0
    out = minimize(fun_min, par1, method = 'Nelder-Mead', options = {'maxiter' : 10000, 'fatol': 1.0e-15} )
    
    print (out)

    wl_shifted, wlc_shifted = get_shifted_sample(out.x)
    wl_range = [np.max([wl[0][0], wl_shifted[0]]), np.min([wl[0][-1], wl_shifted[-1]])]

    axes[0].plot(wl[0], wlc_sample_norm[0], label = 'Normalized WLSc with ref sample (det {:})'.format(descriptions[0]))
    axes[0].plot(wl[0], wlc_shifted, color='C1', linestyle='dashed', 
    label = 'Calibrated normalized WLSc with ref sample  (det {:})'.format(descriptions[1]) )

    axes[1].plot(wl[0], wlc_sample_norm[0] - wlc_sample_norm[1], label='before calibration')
    axes[1].plot(wl[0], wlc_sample_norm[0] - wlc_shifted, label='after calibration')
    
    for ax in axes:
        ax.set_xlim(wl_range)
        ax.set_xlabel('Wavelength, nm')
        ax.grid(True)
        ax.legend()

    axes[0].set_ylabel('Voltage, V')
    axes[1].set_ylabel('Difference (det {:} - det {:}) voltage, V'.format(descriptions[0], descriptions[1]))
    
    return  {"polynomial" : list(out.x) }

def print_carpet_view_header(fnames, spectra_per_acquisition, wavelength_axes):
    for i,fname in enumerate(fnames):
        with open(fname, 'w') as f: 
            f.write('Pump-probe, Not referenced, {:} shots per spectrum,\n'.format(spectra_per_acquisition) +
                    'Background: PumpPower, SmplBckPumped, RefBckPumped, SmplBckUnpumped, RefBckUnpumped,\n' + 
                    'Measurement: PumpPower, SmplMeasPumped, RefMeasPumped, SmplMeasUnpumped, RefMeasUnpumped\n' + 
                    'Wavelength: ' + '\t'.join(['{:.8e}'.format(wl) for wl in wavelength_axes[i]]) + '\n')

def print_carpet_view_line(fname, measurement_type, scan, pump, delays, intensities, data_lines):
    with open(fname, 'a') as f:
        f.write(measurement_type + ', scan {:}, Delay {:.8e} {:.8e}, Intensity {:.8e} {:.8e}'.format(scan, delays[0], delays[1], intensities[0], intensities[1]) + '\n')
        f.write('Pump={:.8e}\n'.format(pump))
        for line in data_lines:
            f.write('\t'.join(['{:.8e}'.format(item) for item in line]) + '\n')


def print_carpet_view_measurements(fnames, info, is_background, i_run):
    for i, fname in enumerate(fnames):        
        measurement_type = 'Background' if is_background else 'Measurement'
        sig_key = 'out_bckg' if is_background else 'out_signal'

        print_carpet_view_line(
            fname,
            measurement_type,
            i_run+1,
            info['pump_high'],
            [info['delay'], 0.0], 
            [0.0, 0.0], 
            [
                info[sig_key][i]['pumped'],
                info[sig_key][i]['pumped'],
                info[sig_key][i]['not_pumped'],
                info[sig_key][i]['not_pumped']
                ])        
def average_in_range(signal, length = 256):
    '''Calculates mean value of whole spectral range (returns array of length defined by 'length')'''
    return [np.average(signal[i::length]) for i in range(length)]

def average_in_place(signal, length = 256):
    '''Calculates mean values in segments of given length of a 1d array'''
    out = [list([np.average(signal[i*length:(i*length + length)])])*length for i in range(int(len(signal)/length))]    
    return [el for row in out for el in row]

def get_pumped_notpumped(spectrum, pd_integral, pumped_uncertainty, not_pumped_uncertainty, wavelength_axis = None, wavelength_axis2 = None, datapoints_per_spectrum = 256, scale = False):
    '''Separates pumped and not pumped spectra from continuous array according to given pd_integral value and uncertainties'''
    gap = np.max(pd_integral) - np.min(pd_integral)
    pumped_floor = np.max(pd_integral) - gap * pumped_uncertainty
    not_pumped_ceil = np.min(pd_integral) + gap * not_pumped_uncertainty
    
    pumped_signal_idx = [1 if pd_integral[i*datapoints_per_spectrum] >= pumped_floor else 0 for i in range(int(len(spectrum)/datapoints_per_spectrum))]
    not_pumped_signal_idx = [1 if pd_integral[i*datapoints_per_spectrum] <= not_pumped_ceil else 0 for i in range(int(len(spectrum)/datapoints_per_spectrum))]
    
    pumped_signal = np.zeros(datapoints_per_spectrum)
    not_pumped_signal = np.zeros(datapoints_per_spectrum)
    
    for i in range(len(pumped_signal_idx)):
        if pumped_signal_idx[i] != 0:
            pumped_signal = pumped_signal + np.array(spectrum[i*datapoints_per_spectrum:((i+1)*datapoints_per_spectrum)])
        if not_pumped_signal_idx[i] != 0:
            not_pumped_signal = not_pumped_signal + np.array(spectrum[i*datapoints_per_spectrum:((i+1)*datapoints_per_spectrum)])
            
    pumped_signal = pumped_signal / np.sum(pumped_signal_idx)
    not_pumped_signal = not_pumped_signal / np.sum(not_pumped_signal_idx)
    
    
    if wavelength_axis is not None and wavelength_axis2 is not None:
        return {'pumped': np.interp(wavelength_axis2, wavelength_axis, pumped_signal), 
                'not_pumped': np.interp(wavelength_axis2, wavelength_axis, not_pumped_signal), 
                'pumped_count': np.sum(pumped_signal_idx), 
                'not_pumped_count': np.sum(not_pumped_signal_idx)}
    else:
        return {'pumped': pumped_signal, 'not_pumped': not_pumped_signal, 'pumped_count': np.sum(pumped_signal_idx), 'not_pumped_count': np.sum(not_pumped_signal_idx)}