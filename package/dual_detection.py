#!/usr/bin/env python
# -*- coding: utf-8 -*-
#==========================================================================
# Harpia REST API Interface example
#--------------------------------------------------------------------------
# Copyright (c) 2021 Light Conversion, UAB
# All rights reserved.
# www.lightcon.com
#==========================================================================
# This example requires `lightcon` package installed into the current, active
# Python environment. Use `pip install lightcon` to install the package.
# 
# In this example, uprocessed raw data is read from Harpia and the neccessary separation 
# into pumped/not-pumped data sets is performed pump-probe spectrum is calculated. 
# The background signal is calculated and separated analogously. This example
# follow the same logic and mathematics which is done by Harpia Service App and
# calculated spectrum is retrieved through harpia.pump_probe_spectrum(). The same
# approach given here can be applied to calculating multi-pulse experiments when
# AuxiliarySignal is available.
     
from lightcon.harpia import Harpia
import lightcon.style

import matplotlib.pyplot as plt
import numpy as np
import sys
import json

#lightcon.style.apply_style()

def average_in_place(signal, length = 256):
    '''Calculates mean values in segments of given length of a 1d array'''
    out = [list([np.average(signal[i*length:(i*length + length)])])*length for i in range(int(len(signal)/length))]    
    return [el for row in out for el in row]

def get_pumped_notpumped(spectrum, pd_integral, pumped_uncertainty, not_pumped_uncertainty, datapoints_per_spectrum = 256, scale = False):
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
    
    
    if scale and scale_wl_axis:
        return {'pumped': np.interp(wavelength_axis2, wavelength_axis, pumped_signal), 
                'not_pumped': np.interp(wavelength_axis2, wavelength_axis, not_pumped_signal), 
                'pumped_count': np.sum(pumped_signal_idx), 
                'not_pumped_count': np.sum(not_pumped_signal_idx)}
    else:
        return {'pumped': pumped_signal, 'not_pumped': not_pumped_signal, 'pumped_count': np.sum(pumped_signal_idx), 'not_pumped_count': np.sum(not_pumped_signal_idx)}
    
if __name__ == "__main__":
    ip_address = '127.0.0.1'
    number_of_spectra = 200
    pumped_uncertainty = 0.15
    not_pumped_uncertainty = 0.15
    out_signal = [None, None]
    out_bckg = [None, None]
    lines_wlc = [None, None]
    lines_pp = [None, None]
    is_plotted = False
    scale_wl_axis = True
    with open("calibration.json", "r") as f:
        scale_wl_poly = json.loads(f.read())['polynomial']
    
    def plot_pump_probe_spectra(outs, outs_bckg, descriptions = None):
        global is_plotted    
        
        if not descriptions:
            descriptions = ['detector {:}'.format(i+1) for i,_ in enumerate(outs)]
            
        fig = plt.figure('PP spectrum' + (' aligned' if scale_wl_axis else ''), figsize = (10, 8))
        
        if not is_plotted:
            plt.clf()
            
        actual_delay = harpia.delay_line_actual_delay()
        
        ax1 = plt.subplot(211)
        ax2 = plt.subplot(212)
        
        if not is_plotted:
            for i in [0,1]:
                lines_wlc[i], = ax1.plot(wavelength_axis, outs[i]['not_pumped'], label='WLC {:}'.format(descriptions[i]))
            for i in [0,1]:
                lines_pp[i], = ax2.plot(wavelength_axis, 
                        -1000.0 * np.log10((outs[i]['pumped'] - outs_bckg[i]['pumped'])/(outs[i]['not_pumped'] - outs_bckg[i]['not_pumped'])))
                ax2.set_ylim([-5, 5])
        
            ax1.set_title('WLC spectrum')
            ax1.set_xlabel('Wavelength (nm)')
            ax1.set_ylabel('Voltage (V)')
            ax1.set_xlim(wl_range)
            ax1.grid()
            
        
            ax2.set_title('Pump-probe spectrum')
            ax2.set_xlabel('Wavelength (nm)')
            ax2.set_ylabel('$\Delta A$ (mOD)')
            ax2.set_xlim(wl_range)
            ax2.grid()
            plt.tight_layout()
            
        for i, out in enumerate(outs):
            lines_wlc[i].set_xdata(wavelength_axis)
            lines_wlc[i].set_ydata(out['not_pumped'])
    
        
        for i, out in enumerate(outs):
            out_bckg = outs_bckg[i]
            lines_pp[i].set_xdata(wavelength_axis)
            lines_pp[i].set_ydata(-1000.0 * np.log10((out['pumped'] - out_bckg['pumped'])/(out['not_pumped'] - out_bckg['not_pumped']))) 
            lines_pp[i].set_label('{:}, $T={:.2f}$ ps'.format(descriptions[i], actual_delay))
            
        ax1.legend()
        ax2.legend()
        is_plotted = True
        fig.canvas.draw()
        plt.pause(0.1)
#        print ('Calculated pump-probe spectrum at {:.3f} ps delay by averaging {:} pumped and {:} not pumped spectra'.format(harpia.delay_line_actual_delay(), out['pumped_count'], out['not_pumped_count']))

    
    # initialize connection to Harpia 
    harpia = Harpia(ip_address)
    
    # check if connection successful
    if not harpia.connected:
        sys.exit("Could not connect to Harpia")
        
    # ensure chopper is started
    harpia.chopper_start()
    
    # number of datapoints per spectrum
    datapoints_per_spectrum = harpia.datapoints_per_spectrum()
    
    # prepare for background signal measurement
    harpia.close_all_shutters()
    harpia.open_pump_shutter()
    
    # we will be collecting 5x more data for background signal
    harpia.set_spectra_per_acquisition(5*number_of_spectra)
    data_bckg = harpia.raw_signal()
    
    # measure pump-probe signals
    harpia.open_all_shutters()

    # restore number of spectra per acquisition
    harpia.set_spectra_per_acquisition(number_of_spectra)    
    
    # total number of datapoints per channel
    number_of_datapoints = number_of_spectra * harpia.datapoints_per_spectrum()
    
    wavelength_axis = np.array(harpia.wavelength_axis())
    
    if scale_wl_axis:
        wavelength_axis2 = wavelength_axis + np.poly1d(scale_wl_poly)(np.arange(len(wavelength_axis)))
        wl_range = [np.max([wavelength_axis2[0], wavelength_axis[0]]), np.min([wavelength_axis2[-1], wavelength_axis[-1]])]        
    else:
        wl_range = [wavelength_axis[0], wavelength_axis[-1]]

    while True:
        try:
            data = harpia.raw_signal()
            
            # separate measured pump-probe signal into averaged pumped/not-pumped spectra 
            # for detector 1
            out_signal[0] = get_pumped_notpumped(data['SpectrometerSignal'], average_in_place(data['PumpPhotodetectorSignal']), pumped_uncertainty, not_pumped_uncertainty)
            out_signal[1] = get_pumped_notpumped(data['AuxiliarySignal'], average_in_place(data['PumpPhotodetectorSignal']), pumped_uncertainty, not_pumped_uncertainty, scale = True)
            
            # separate measured background signal into averaged pumped/not-pumped spectra
            out_bckg[0] = get_pumped_notpumped(data_bckg['SpectrometerSignal'], average_in_place(data_bckg['PumpPhotodetectorSignal']), pumped_uncertainty, not_pumped_uncertainty)
            out_bckg[1] = get_pumped_notpumped(data_bckg['AuxiliarySignal'], average_in_place(data_bckg['PumpPhotodetectorSignal']), pumped_uncertainty, not_pumped_uncertainty, scale = True)
            
            # plot pump-probe spectrum
            plot_pump_probe_spectra(out_signal, out_bckg, ['V', 'H'])
        except KeyboardInterrupt:
            plot_pump_probe_spectra(out_signal, out_bckg, ['V', 'H'])
            # close all shutters
            harpia.close_all_shutters()
            sys.exit("")
        except:
            pass
    
    
    # close all shutters
    harpia.close_all_shutters()