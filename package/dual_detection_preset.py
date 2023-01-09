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
import pickle
import time
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
                'not_pumped_count': np.sum(not_pumped_signal_idx),
                'pump_high' : np.max(pd_integral)}
    else:
        return {'pumped': pumped_signal, 'not_pumped': not_pumped_signal, 'pumped_count': np.sum(pumped_signal_idx), 'not_pumped_count': np.sum(not_pumped_signal_idx), 'pump_high' : np.max(pd_integral)}
    
    
if __name__ == "__main__":
    ip_address = '10.0.8.5'
    sample_name = '3_degraded'
    number_of_spectra = 2000
    pumped_uncertainty = 0.15
    not_pumped_uncertainty = 0.15
    number_of_runs = 1
    cuts = [505.626, 764.700, 707.884]
    out_signal = [None, None]
    out_bckg = [None, None]
    lines_wlc = [None, None]
    lines_pp = [None, None]
    lines_cuts = [None, None]
    lines_cuts[0] = [None] * len(cuts)
    lines_cuts[1] = [None] * len(cuts)
    is_plotted = False
    scale_wl_axis = True
    # scale_wl_poly = [ 4.30664195e-05,  1.68839201e-03, -4.11486921e-01]
    scale_wl_poly = [ 3.54640396e-05,  8.26729511e-03, -3.10629795e+00 ]    # 03-07 calibration

    delays = np.concatenate((np.arange(-5, -1., 0.2), np.arange(-1, 1, 0.03), np.arange(1, 2, 0.1), np.arange(2, 5, 0.2), np.arange(5, 10, 0.5), np.arange(10, 20, 0.75), np.arange(20, 50, 1.5), np.arange(50, 100, 5)))

    def plot_pump_probe_spectra(outs, outs_bckg, descriptions = None, figname = 'PP'):
        global is_plotted    
        global data
        
        if not descriptions:
            descriptions = ['detector {:}'.format(i+1) for i,_ in enumerate(outs)]
            
        fig = plt.figure(figname, figsize = (10, 8))
        
        if not is_plotted:
            plt.clf()
            
        actual_delay = harpia.delay_line_actual_delay()
        
        ax1 = plt.subplot(311)
        ax2 = plt.subplot(312)
        ax3 = plt.subplot(313)
        
        if not is_plotted:
            for i in [0,1]:
                lines_wlc[i], = ax1.plot(wavelength_axis, outs[i]['not_pumped'], label='{:}'.format(descriptions[i]))
            for i in [0,1]:
                lines_pp[i], = ax2.plot(wavelength_axis, 
                        -1000.0 * np.log10((outs[i]['pumped'] - outs_bckg[i]['pumped'])/(outs[i]['not_pumped'] - outs_bckg[i]['not_pumped'])))
                ax2.set_ylim([-5, 5])
            for iic, ic in enumerate(cuts_idx):
                for i in [0,1]:
                    sp_key = ['spectra', 'spectra_aux'][i]
                    bckg_key = ['background_signals', 'background_signals_aux'][i]
                    linestyle = ['solid', 'dashed'][i]
                    x = [data['delays'][j] for j,sp in enumerate(data['spectra']) if sp]
                    y = [-1000.0 * np.log10((data[sp_key][j]['pumped'][ic] - data[bckg_key]['pumped'][ic])/(data[sp_key][j]['not_pumped'][ic] - data[bckg_key]['not_pumped'][ic])) for j,sp in enumerate(data['spectra']) if sp]
                    lines_cuts[i][iic], = ax3.plot(x, y, color = 'C{:}'.format(iic), linestyle = linestyle, label = '{:.0f} nm, {:}'.format(cuts[iic], descriptions[i]))
                    
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
            
            ax3.set_title('Pump-probe cuts')
            ax3.set_xlabel('Delay (ps)')
            ax3.set_ylabel('$\Delta A$ (mOD)')
            ax3.grid()
            ax3.legend()
            
            plt.tight_layout()
            
        for i, out in enumerate(outs):
            lines_wlc[i].set_xdata(wavelength_axis)
            lines_wlc[i].set_ydata(out['not_pumped'])
    
        
        for i, out in enumerate(outs):
            out_bckg = outs_bckg[i]
            lines_pp[i].set_xdata(wavelength_axis)
            lines_pp[i].set_ydata(-1000.0 * np.log10((out['pumped'] - out_bckg['pumped'])/(out['not_pumped'] - out_bckg['not_pumped']))) 
            lines_pp[i].set_label('{:}, $T={:.2f}$ ps'.format(descriptions[i], actual_delay))
            
        for iic, ic in enumerate(cuts_idx):
            for i in [0,1]:
                sp_key = ['spectra', 'spectra_aux'][i]
                bckg_key = ['background_signals', 'background_signals_aux'][i]
                x = [data['delays'][j] for j, sp in enumerate(data['spectra']) if sp]
                y = [-1000.0 * np.log10((data[sp_key][j]['pumped'][ic] - data[bckg_key]['pumped'][ic])/(data[sp_key][j]['not_pumped'][ic] - data[bckg_key]['not_pumped'][ic])) for j,sp in enumerate(data['spectra']) if sp]
                
                lines_cuts[i][iic].set_xdata(x)
                lines_cuts[i][iic].set_ydata(y)
                ax3.set_xlim([np.min([ax3.get_xlim()[0], np.min(x)]), np.max([ax3.get_xlim()[1], np.max(x)])])
                ax3.set_ylim([np.min([ax3.get_ylim()[0], np.min(y)]), np.max([ax3.get_ylim()[1], np.max(y)])])
            
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
    
    wavelength_axis = np.array(harpia.wavelength_axis())


    # total number of datapoints per channel
    number_of_datapoints = number_of_spectra * harpia.datapoints_per_spectrum()
    
    wavelength_axis = np.array(harpia.wavelength_axis())
    
    if scale_wl_axis:
        wavelength_axis2 = wavelength_axis + np.poly1d(scale_wl_poly)(np.arange(len(wavelength_axis)))
        wl_range = [np.max([wavelength_axis2[0], wavelength_axis[0]]), np.min([wavelength_axis2[-1], wavelength_axis[-1]])]        
    else:
        wl_range = [wavelength_axis[0], wavelength_axis[-1]]
        
    # separate measured background signal into averaged pumped/not-pumped spectra
    out_bckg[0] = get_pumped_notpumped(data_bckg['SpectrometerSignal'], average_in_place(data_bckg['PumpPhotodetectorSignal']), pumped_uncertainty, not_pumped_uncertainty)
    out_bckg[1] = get_pumped_notpumped(data_bckg['AuxiliarySignal'], average_in_place(data_bckg['PumpPhotodetectorSignal']), pumped_uncertainty, not_pumped_uncertainty, scale = True)    
        
    cuts_idx = [np.argmin(np.abs(wavelength_axis - cut_wl)) for cut_wl in cuts]
    
    data = {}
    data['wavelenghts'] = np.array(wavelength_axis)
    data['delays'] = delays
    data['number_of_runs'] = 1
    data['spectra_per_acquisition'] = number_of_spectra
    data['background_signals'] = out_bckg[0]
    data['background_signals_aux'] = out_bckg[1]
    data['spectra'] = [None] * len(delays)
    data['spectra_aux'] = [None] * len(delays)
    
    figname = 'PP, sample {:} pump {:.1f}%'.format(sample_name, harpia.pump_vndf_transmittance())
    timestamp = time.strftime('%Y%m%d_%H_%M_%S')
    dump_filename = 'pp-int_'+timestamp+'.p'

    for i, delay in enumerate(data['delays']):
        try:            
            harpia.set_delay_line_target_delay(delay)
            
            # restore number of spectra per acquisition
            harpia.set_spectra_per_acquisition(number_of_spectra * (1.0 if np.abs(delay) > 0.5 else 1.0))    
                
            # total number of datapoints per channel
            number_of_datapoints = number_of_spectra * harpia.datapoints_per_spectrum()
            
            raw = harpia.raw_signal()            
            # separate measured pump-probe signal into averaged pumped/not-pumped spectra 
            # for detector 1
            out_signal[0] = get_pumped_notpumped(raw['SpectrometerSignal'], average_in_place(raw['PumpPhotodetectorSignal']), pumped_uncertainty, not_pumped_uncertainty)
            out_signal[1] = get_pumped_notpumped(raw['AuxiliarySignal'], average_in_place(raw['PumpPhotodetectorSignal']), pumped_uncertainty, not_pumped_uncertainty, scale = True)
            
            data['spectra'][i] = out_signal[0]
            data['spectra_aux'][i] = out_signal[1]            
            
            # plot pump-probe spectrum
            plot_pump_probe_spectra(out_signal, out_bckg, ['V', 'H'], figname = figname)
        except KeyboardInterrupt:
            plot_pump_probe_spectra(out_signal, out_bckg, ['V', 'H'], figname = figname)
            # close all shutters
            harpia.close_all_shutters()
            sys.exit("")
    
    plt.savefig(figname+ ' ' + timestamp + '.pdf')
    
    pickle.dump(data, open(dump_filename, 'wb' ))
    
    # close all shutters
    harpia.close_all_shutters()
    
    #%%
# =============================================================================
# MEASUREMENT SLICES TO CarpetView
# =============================================================================


def print_measurement_lines(measurement_type, scan, pump, delays, intensities, data_lines):
    f.write(measurement_type + ', scan {:}, Delay {:.8e} {:.8e}, Intensity {:.8e} {:.8e}'.format(scan, delays[0], delays[1], intensities[0], intensities[1]) + '\n')
    f.write('Pump={:.8e}\n'.format(pump))
    for line in data_lines:
        f.write('\t'.join(['{:.8e}'.format(item) for item in line]) + '\n')

for fname, bckg_key, sig_key in [(figname+ ' ' + timestamp + '_det1.dat', 'background_signals', 'spectra'),
                                 (figname+ ' ' + timestamp + '_det2.dat', 'background_signals_aux', 'spectra_aux')]:
    
    f = open(fname, 'w')    
    
    f.write('Pump-probe, Not referenced, {:} shots per spectrum,\n'.format(data['spectra_per_acquisition']) +
            'Background: PumpPower, SmplBckPumped, RefBckPumped, SmplBckUnpumped, RefBckUnpumped,\n' + 
            'Measurement: PumpPower, SmplMeasPumped, RefMeasPumped, SmplMeasUnpumped, RefMeasUnpumped\n' + 
            'Wavelength: ' + '\t'.join(['{:.8e}'.format(wl) for wl in data['wavelenghts']]) + '\n')
    
    for i_run in np.arange(number_of_runs):
        print_measurement_lines('Background', i_run+1, data[bckg_key]['pump_high'], [data['delays'][0], 0.0], [0.0, 0.0], 
                                [data[bckg_key]['pumped'],
                                 data[bckg_key]['pumped'],
                                 data[bckg_key]['not_pumped'],
                                 data[bckg_key]['not_pumped']])
        
        for j,delay in enumerate(data['delays']):
            print_measurement_lines('Measurement', i_run+1, data[sig_key][j]['pump_high'],[data['delays'][j], 0.0], [0.0, 0.0], 
                                    [data[sig_key][j]['pumped'],
                                     data[sig_key][j]['pumped'],
                                     data[sig_key][j]['not_pumped'],
                                     data[sig_key][j]['not_pumped'],])
            
    print('Exported CarpetView-compatible data to {:}'.format(fname))
    f.close()
    
    fname_matrix = fname[:-4]+'_matrix.txt'
    
    mat = [-1000.0 * np.log10((data[sig_key][i]['pumped'] - data[bckg_key]['pumped'])/(data[sig_key][i]['not_pumped'] - data[bckg_key]['not_pumped'])) for i,_ in enumerate(data['delays'])]
    mat = np.hstack([np.array(data['delays'], ndmin=2).transpose(), mat])
    mat = np.vstack([np.concatenate([[0], np.array(wavelength_axis)]), mat])
    
    np.savetxt(fname_matrix, mat, header = 'XAxisTitle Wavelength (nm)\nYAxisTitle Delay (ps)', comments='')
    print('Exported matrix data to {:}'.format(fname_matrix))