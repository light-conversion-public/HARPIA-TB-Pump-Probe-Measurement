#!/usr/bin/env python
# -*- coding: utf-8 -*-
#==========================================================================
# Harpia REST API Interface example
#--------------------------------------------------------------------------
# Copyright (c) 2021 Light Conversion, UAB
# All rights reserved.
# www.lightcon.com
#==========================================================================
     
import lclauncher

import os
import sys    
import time
import numpy as np
import lightcon.style
import json
from utils import *

import matplotlib
matplotlib.use('Qt5Agg')

from PyQt5.QtWidgets import *
from PyQt5.QtGui import QIcon
from PyQt5.QtCore import pyqtSlot, QObject, QThread, pyqtSignal
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure

lightcon.style.apply_style()

sys.path.append(os.path.dirname(os.path.realpath(sys.argv[0])))
os.chdir(os.path.dirname(os.path.realpath(sys.argv[0])))

# if connections to devices are used, they are initiated here:
connections = lclauncher.establish_connections()

# initialize and connect to HARPIA
harpia = connections.get_connection('harpia')

# check if connection successful
if not harpia:
    sys.exit("Could not connect to Harpia")
    
# =============================================================================
# Parameters
# =============================================================================
settings = None

def load_settings():
    global settings
    with open("./package/settings.json", "r") as f:
        settings = json.loads(f.read())

def save_settings():
    global settings
    with open("./package/settings.json", "w") as f:
        json.dump(settings, f)

load_settings()

pumped_uncertainty = settings.get('pumped_uncertainty') or 0.15
not_pumped_uncertainty = settings.get('not_pumped_uncertainty') or 0.15

descriptions = settings.get('descriptions') or ["H", "V"]

wlc = [None, None]
wlc_sample = [None, None]

# =============================================================================
#  Parameters - END
# =============================================================================



class Worker(QObject):    
    finished = pyqtSignal()
    progress = pyqtSignal(dict)
    out_signal = [None, None]
    out_bckg = [None, None]
    
    def stop(self):
        self.is_running = False

    def measure_preset(self):
        self.is_running = True     
        
        self.actions_before_measurement()

        preset = harpia._get('/Basic/CurrentPreset')

        if preset.get('IsMeasurementTA'):
            number_of_runs = preset.get('NumberOfRuns') or 1
            delays = preset.get('DelayTimes') or []

            print_carpet_view_header(fnames, self.spectra_per_acquisition, [self.wavelength_axis, self.wavelength_axis2])

            for irun in np.arange(number_of_runs):
                harpia.set_delay_line_target_delay(preset['DelayTimes'][0])

                self.measure_background()


                for idelay, delay in enumerate(preset['DelayTimes']):
                    if not self.is_running:
                        break
                    harpia.set_delay_line_target_delay(delay)
                    info = self.measure_once()
                    if info is not None:
                        info['for_preset'] = True

                        if idelay==0:
                            print_carpet_view_measurements(fnames, info, True, irun)

                        print_carpet_view_measurements(fnames, info, False, irun)

                        self.progress.emit(info)

                if not self.is_running:
                        break

        print('stop')
        
        self.actions_after_measurement()
        
        self.finished.emit()

    def measure_continuously(self):
        self.is_running = True     
        
        self.actions_before_measurement()
                
        while self.is_running:      
            info = self.measure_once()
            if info is not None:
                self.progress.emit(info)
        
        print('stop')
        
        self.actions_after_measurement()
        
        self.finished.emit()
        
    def measure_once(self):
        try:
            data = harpia.raw_signal()
            
            out_signal=[get_pumped_notpumped(data['SpectrometerSignal'], average_in_place(data['PumpPhotodetectorSignal']), pumped_uncertainty, not_pumped_uncertainty),
                        get_pumped_notpumped(data['AuxiliarySignal'], average_in_place(data['PumpPhotodetectorSignal']), pumped_uncertainty, not_pumped_uncertainty, wavelength_axis=self.wavelength_axis, wavelength_axis2=self.wavelength_axis2)]
            out_bckg = [get_pumped_notpumped(self.data_bckg['SpectrometerSignal'], average_in_place(self.data_bckg['PumpPhotodetectorSignal']), pumped_uncertainty, not_pumped_uncertainty),
                        get_pumped_notpumped(self.data_bckg['AuxiliarySignal'], average_in_place(self.data_bckg['PumpPhotodetectorSignal']), pumped_uncertainty, not_pumped_uncertainty, wavelength_axis=self.wavelength_axis, wavelength_axis2=self.wavelength_axis2)]
                    
            return {
                'delay': harpia.delay_line_actual_delay(),
                'out_signal': out_signal,
                'out_bckg': out_bckg,
                'number_of_spectra' : self.spectra_per_acquisition,
                'pump_high' : np.max(data['PumpPhotodetectorSignal']),
                'spectra' : [-1000.0 * np.log10((np.abs(out_signal[i]['pumped'] - out_bckg[i]['pumped'])/(out_signal[i]['not_pumped'] - out_bckg[i]['not_pumped']))) for i in [0,1]],
                #'spectra' : [np.random.randn(256), np.random.randn(256)],
                'wavelength_axis' : self.wavelength_axis,
                'wavelength_axis2' : self.wavelength_axis2,
                'wavelength_range' : self.wavelength_range
                }
            
        except Exception as e:
            print (e)
            self.actions_after_measurement()
        except:
            pass
        
    def measure_background(self):
        # prepare for background signal measurement
        harpia.close_all_shutters()
        harpia.open_pump_shutter()
        
        # we will be collecting 5x more data for background signal
        harpia.set_spectra_per_acquisition(5*self.spectra_per_acquisition)
        self.data_bckg = harpia.raw_signal()
        
        # measure pump-probe signals
        harpia.open_all_shutters()
    
        # restore number of spectra per acquisition
        harpia.set_spectra_per_acquisition(self.spectra_per_acquisition)    


    def actions_before_measurement(self):
        # ensure chopper is started
        harpia.chopper_start()
        
        # number of datapoints per spectrum
        datapoints_per_spectrum = harpia.datapoints_per_spectrum()
        
        self.spectra_per_acquisition = harpia.spectra_per_acquisition()

        self.measure_background()
        
        # total number of datapoints per channel
        self.number_of_datapoints = self.spectra_per_acquisition * harpia.datapoints_per_spectrum()
        
        self.wavelength_axis = np.array(harpia.wavelength_axis())        
        
        self.scale_wl_poly = settings['calibration']['polynomial']
            
        self.wavelength_axis2 = self.wavelength_axis + np.poly1d(self.scale_wl_poly)(np.arange(len(self.wavelength_axis)))
        self.wavelength_range = [np.max([self.wavelength_axis2[0], self.wavelength_axis[0]]), np.min([self.wavelength_axis2[-1], self.wavelength_axis[-1]])]            
        
    def actions_after_measurement(self):
        harpia.close_all_shutters()
        
class MplCanvas(FigureCanvasQTAgg):
    lines_wlc = [None, None]
    lines_pp = [None, None]
    
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        lightcon.style.apply_style()
        self.fig = Figure(figsize=(width, height), dpi=dpi)                
        
        self.ax1 = self.fig.add_subplot(211)
        self.ax2 = self.fig.add_subplot(212)
                
        plt.ion()
        
        for i in [0,1]:
            self.lines_wlc[i], = self.ax1.plot([],[])
            self.lines_pp[i], = self.ax2.plot([],[])
            
        self.ax1.set_title('WLSc spectrum')
        self.ax1.set_xlabel('Wavelength (nm)')
        self.ax1.set_ylabel('Voltage (V)')
        self.ax2.legend()
        self.ax1.grid(True)        
    
        self.ax2.set_title('Pump-probe spectrum')
        self.ax2.set_xlabel('Wavelength (nm)')
        self.ax2.set_ylabel('$\Delta A$ (mOD)')
        self.ax2.grid(True)
        self.ax2.legend()
        self.fig.tight_layout()
                
        super(MplCanvas, self).__init__(self.fig)        

class MainWindow(QMainWindow):
    sc = None
    canDraw = True
    worker = None
    calibration_done = [False,False]
    working_directory = settings.get('working_directory') or r"C:\Users\Public\Desktop"
    
    def __init__(self, title):
        super().__init__()
        self.title = title
        self.left = 100
        self.top = 100
        self.width = 1200
        self.height = 900
        self.initUI()
        
    def initUI(self):
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)
        
        self.measure_continuously_button = QPushButton('START/STOP\nCONTINUOUS ACQUISITION', self)
        self.measure_continuously_button.setToolTip('Start measuring')
        # self.measure_continuously_button.setFixedHeight(75)
        self.measure_continuously_button.clicked.connect(self.measure_continuously_button_on_click)
                
        self.measure_preset_button = QPushButton('START/STOP\nPRESET MEASUREMENT', self)
        self.measure_preset_button.setToolTip('Start measuring preset')
        self.measure_preset_button.clicked.connect(self.measure_preset_button_on_click)
        
        self.working_directory_line = QLineEdit()
        self.working_directory_line.setText(self.working_directory)
        
        self.wlc_no_sample_button = QPushButton('Measure WlSc without ref', self)
        self.wlc_no_sample_button.setToolTip('Measure WLSc without reference sample')
        self.wlc_no_sample_button.clicked.connect(self.wlc_no_sample_button_on_click)
        
        self.wlc_sample_button = QPushButton('Measure WlSc with ref', self)
        self.wlc_sample_button.setToolTip('Measure WLSc with reference sample')
        self.wlc_sample_button.clicked.connect(self.wlc_sample_button_on_click)
        
        self.wl_calibrate_button = QPushButton('CALIBRATE', self)
        self.wl_calibrate_button.setToolTip('Calibrate wavelength axis')
        # self.wl_calibrate_button.setEnabled(False)
        self.wl_calibrate_button.clicked.connect(self.wl_calibrate_button_on_click)
        
        self.sc = [MplCanvas(self, width=6, height=3, dpi=100)]

        outerLayout = QHBoxLayout()
        topLayout = QVBoxLayout()            
        bottomLayout = QVBoxLayout()
        
        presetLayout = QVBoxLayout()
        presetLayout.addWidget(self.measure_preset_button)
        presetLayout.addWidget(QLabel("Working directory:"))
        presetLayout.addWidget(self.working_directory_line)
        
        calibrationLayout = QVBoxLayout()        
        calibrationLayout.addWidget(self.wlc_no_sample_button)
        calibrationLayout.addWidget(self.wlc_sample_button)
        calibrationLayout.addWidget(self.wl_calibrate_button)
        calibrationLayout.addStretch()
        
        calibrationGroupBox = QGroupBox("Wavelength axis calibration")
        calibrationGroupBox.setFixedHeight(130)

        presetGroupBox = QGroupBox("Preset measurement")
        

        topLayout.addWidget(self.measure_continuously_button, 1)
        topLayout.addWidget(presetGroupBox)
        topLayout.addWidget(calibrationGroupBox)
        topLayout.addStretch()
        # topLayout.addWidget(self.reset_button, 1)
        # topLayout.addWidget(self.clear_button, 1)

        calibrationGroupBox.setFlat(True)
        calibrationGroupBox.setLayout(calibrationLayout)
        
        presetGroupBox.setFlat(True)
        presetGroupBox.setLayout(presetLayout)
        presetGroupBox.setFixedHeight(130)

        bottomLayout.addWidget(self.sc[0])        
        
        outerLayout.addLayout(topLayout,1)
        outerLayout.addLayout(bottomLayout,3)
        
        widget = QWidget()
        widget.setLayout(outerLayout)

        self.setCentralWidget(widget)  
        self.show()
        
    def measure_continuously_task(self):        
        if self.worker:
            if self.worker.is_running:
                self.worker.is_running = False
                return
        
        self.thread = QThread()
        self.worker = Worker()
        self.worker.moveToThread(self.thread)
        self.thread.started.connect(self.worker.measure_continuously)
        self.worker.finished.connect(self.thread.quit)
        self.worker.finished.connect(self.worker.deleteLater)
        self.thread.finished.connect(self.thread.deleteLater)
        self.worker.progress.connect(self.updatePlots)
        self.thread.start()
        # Final resets       
        # self.thread.finished.connect(self.addToPlots)

    def measure_preset_task(self):        
        if self.worker:
            if self.worker.is_running:
                self.worker.is_running = False
                return
                
        global fnames
        global descriptions
        timestamp = time.strftime('%Y%m%d_%H%M')
        fnames = [os.path.join(self.working_directory_line.text(), timestamp+"_det" + descriptions[i] + ".dat") for i in [0,1]]

        settings['working_directory'] = self.working_directory_line.text()
        save_settings()

        self.thread = QThread()
        self.worker = Worker()
        self.worker.moveToThread(self.thread)
        self.thread.started.connect(self.worker.measure_preset)
        self.worker.finished.connect(self.thread.quit)
        self.worker.finished.connect(self.worker.deleteLater)
        self.thread.finished.connect(self.thread.deleteLater)
        self.worker.progress.connect(self.updatePlots)
        self.thread.start()
        # Final resets       
        # self.thread.finished.connect(self.addToPlots)
    def updatePlots(self, info): 
        if (not self.canDraw): 
            print("Cannot draw")
            return()
        self.canDraw = False
        
        spectra = info['spectra']
        out_signal = info['out_signal']
        out_bckg = info['out_bckg']
        wavelength_axis = info['wavelength_axis']
        wavelength_axis2 = info['wavelength_axis2']
        wavelength_range = info['wavelength_range']        
        
        for i, _ in enumerate(out_signal):
            self.sc[0].lines_wlc[i].set_xdata(wavelength_axis)
            self.sc[0].lines_wlc[i].set_ydata(out_signal[i]['not_pumped'])                


            self.sc[0].lines_pp[i].set_xdata(wavelength_axis)
            self.sc[0].lines_pp[i].set_ydata(spectra[i])            

            self.sc[0].lines_pp[i].set_label('{:}, $T={:.2f}$ ps'.format(descriptions[i], info.get('delay') or 0.0))
            
            
        self.sc[0].ax1.set_xlim(wavelength_range)
        self.sc[0].ax1.set_ylim([0, np.nanmax([np.nanmax(out_signal[0]['not_pumped']), np.nanmax(out_signal[1]['not_pumped'])])])
        self.sc[0].ax1.legend()
            
        self.sc[0].ax2.set_xlim(wavelength_range)
        self.sc[0].ax2.set_ylim([np.nanmin([np.nanmin(spectra[0]), np.nanmin(spectra[1])]),np.max([np.nanmax(spectra[0]), np.nanmax(spectra[1])])])
        self.sc[0].ax2.legend()
                
        self.sc[0].draw()
        app.processEvents()
        self.canDraw = True

    def measure_for_calibration(self,is_without_sample):
        harpia.close_all_shutters()
        harpia.open_probe_shutter()
        
        data = harpia.raw_signal()

        harpia.close_all_shutters()

        wl = np.array(harpia.wavelength_axis())
        
        fname_prefix = 'no_sample_' if is_without_sample else 'sample_'

        wlc = [np.vstack((wl, average_in_range(data['SpectrometerSignal']))), np.vstack((wl,average_in_range(data['AuxiliarySignal'])))]
        np.savetxt('./package/' + fname_prefix + '0_WLSc.txt', np.transpose(wlc[0]), delimiter='\t', header="Wavelength (nm)\tDetector signal (V)")
        np.savetxt('./package/' + fname_prefix + '1_WLSc auxiliary.txt', np.transpose(wlc[1]), delimiter='\t', header="Wavelength (nm)\tDetector signal (V)")
        
    @pyqtSlot()
    def wlc_no_sample_button_on_click(self):
        self.measure_for_calibration(True)
        self.calibration_done[0] = True
        self.wl_calibrate_button.setEnabled(np.all(self.calibration_done))
        pass
    
    @pyqtSlot()
    def wlc_sample_button_on_click(self):
        self.measure_for_calibration(False)
        self.calibration_done[1] = True
        self.wl_calibrate_button.setEnabled(np.all(self.calibration_done))
        pass
    
    @pyqtSlot()
    def wl_calibrate_button_on_click(self):
        wlc = [np.loadtxt('./package/no_sample_0_WLSc.txt', skiprows=1, delimiter='\t'),np.loadtxt('./package/no_sample_1_WLSc auxiliary.txt', skiprows=1, delimiter='\t')]
        wlc_sample = [np.loadtxt('./package/sample_0_WLSc.txt', skiprows=1, delimiter='\t'),np.loadtxt('./package/sample_1_WLSc auxiliary.txt', skiprows=1, delimiter='\t')]

        calibration = get_wavelength_calibration(wlc, wlc_sample, descriptions)
        settings['calibration'] = calibration
        save_settings()
        pass
        
    @pyqtSlot()
    def measure_continuously_button_on_click(self):
        # camera.enable_beam_profiler()
        self.measure_continuously_task()

    @pyqtSlot()
    def measure_preset_button_on_click(self):        
        self.measure_preset_task()

    @pyqtSlot()
    def reset_button_on_click(self):
        pass
        # camera.enable_beam_profiler()
        # self.resetLongTask()
        # self.disableButtons()
        
app = QApplication([])
w = MainWindow('HARPIA Dual detection')
app.exec_()
        