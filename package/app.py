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
settings = {}

def load_settings():
    global settings
    with open("./package/settings.json", "a+") as f:
        pass

    with open("./package/settings.json", "r") as f:
        settings = json.loads(f.read() or "{}")

def save_settings():
    global settings
    with open("./package/settings.json", "w") as f:
        json.dump(settings, f)

load_settings()

pumped_uncertainty = settings.get('pumped_uncertainty') or 0.15
not_pumped_uncertainty = settings.get('not_pumped_uncertainty') or 0.15

descriptions = settings.get('descriptions') or ["TA"]

wlc = [None]
wlc_sample = [None]

# =============================================================================
#  Parameters - END
# =============================================================================



class Worker(QObject):    
    finished = pyqtSignal()
    progress = pyqtSignal(dict)
    out_signal = [None]
    out_bckg = [None]
    
    def stop(self):
        self.is_running = False

    def measure_preset(self):
        self.is_running = True     
        
        self.actions_before_measurement()

        preset = harpia._get('/Basic/CurrentPreset')

        if preset.get('IsMeasurementSWTA') or preset.get('IsMeasurementTA'):
            number_of_runs = preset.get('NumberOfRuns') or 1
            delays = preset.get('DelayTimes') or []

            print_carpet_view_header(fnames, self.spectra_per_acquisition, [self.wavelength_axis])

            for irun in np.arange(number_of_runs):
                harpia.harpiatb_set_delay_line_target_delay(preset['DelayTimes'][0])

                self.measure_background()


                for idelay, delay in enumerate(preset['DelayTimes']):
                    if not self.is_running:
                        break
                    harpia.harpiatb_set_delay_line_target_delay(delay)
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
            
            out_signal=[get_pumped_notpumped(data['SpectrometerSignal'], average_in_place(data['PumpPhotodetectorSignal'], length = self.datapoints_per_spectrum), pumped_uncertainty, not_pumped_uncertainty, datapoints_per_spectrum = self.datapoints_per_spectrum)]
            out_bckg = [get_pumped_notpumped(self.data_bckg['SpectrometerSignal'], average_in_place(self.data_bckg['PumpPhotodetectorSignal'], length = self.datapoints_per_spectrum), pumped_uncertainty, not_pumped_uncertainty, datapoints_per_spectrum = self.datapoints_per_spectrum)]
                    
            return {
                'delay': harpia.harpiatb_delay_line_actual_delay(),
                'out_signal': out_signal,
                'out_bckg': out_bckg,
                'number_of_spectra' : self.spectra_per_acquisition,
                'pump_high' : np.max(data['PumpPhotodetectorSignal']),
                'spectra' : [-1000.0 * np.log10((np.abs(out_signal[i]['pumped'] - out_bckg[i]['pumped'])/(out_signal[i]['not_pumped'] - out_bckg[i]['not_pumped']))) for i in [0]],
                #'spectra' : [np.random.randn(256), np.random.randn(256)],
                'wavelength_axis' : self.wavelength_axis,
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
        harpia.open_third_beam_shutter()
    
        # restore number of spectra per acquisition
        harpia.set_spectra_per_acquisition(self.spectra_per_acquisition)    

    def actions_before_measurement(self):
        # ensure chopper is started
        harpia.chopper_start()
        
        # number of datapoints per spectrum
        datapoints_per_spectrum = harpia.datapoints_per_spectrum()
        
        self.spectra_per_acquisition = harpia.spectra_per_acquisition()

        self.measure_background()
        
        self.datapoints_per_spectrum = harpia.datapoints_per_spectrum()
        # total number of datapoints per channel
        self.number_of_datapoints = self.spectra_per_acquisition * self.datapoints_per_spectrum
        
        self.wavelength_axis = np.array(harpia.wavelength_axis())        
        
        self.wavelength_range = [self.wavelength_axis.min() - 1, self.wavelength_axis.max() + 1]
        
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
        
        for i in [0]:
            self.lines_wlc[i], = self.ax1.plot([],[], '.-')
            self.lines_pp[i], = self.ax2.plot([],[], '.-')
            
        self.ax1.set_title('WLSc spectrum')
        self.ax1.set_xlabel('Delay (ps)')
        self.ax1.set_ylabel('Voltage (V)')
        #self.ax2.legend()
        self.ax1.grid(True)        
    
        self.ax2.set_title('Pump-probe spectrum')
        self.ax2.set_xlabel('Delay (ps)')
        self.ax2.set_ylabel('$\Delta A$ (mOD)')
        self.ax2.grid(True)
        #self.ax2.legend()
        self.fig.tight_layout()
                
        super(MplCanvas, self).__init__(self.fig)        

class MainWindow(QMainWindow):
    sc = None
    canDraw = True
    worker = None
    calibration_done = [False,False]
    working_directory = settings.get('working_directory') or r"C:\Users\Public\Desktop"
    plot_data = {}
    
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
        self.measure_continuously_button.clicked.connect(self.measure_continuously_button_on_click)
                
        self.measure_preset_button = QPushButton('START/STOP\nPRESET MEASUREMENT', self)
        self.measure_preset_button.setToolTip('Start measuring preset')
        self.measure_preset_button.clicked.connect(self.measure_preset_button_on_click)
        
        self.working_directory_line = QLineEdit()
        self.working_directory_line.setText(self.working_directory)
        
        self.sc = [MplCanvas(self, width=6, height=3, dpi=100)]

        outerLayout = QHBoxLayout()
        topLayout = QVBoxLayout()            
        bottomLayout = QVBoxLayout()
        
        presetLayout = QVBoxLayout()
        presetLayout.addWidget(self.measure_preset_button)
        presetLayout.addWidget(QLabel("Working directory:"))
        presetLayout.addWidget(self.working_directory_line)
        
        presetGroupBox = QGroupBox("Preset measurement")

        topLayout.addWidget(self.measure_continuously_button, 1)
        topLayout.addWidget(presetGroupBox)
        topLayout.addStretch()
        
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

        self.y_lims = [[0,-np.inf], [np.inf,-np.inf]]
        self.plot_data = {}

        self.thread = QThread()
        self.worker = Worker()
        self.worker.moveToThread(self.thread)
        self.thread.started.connect(self.worker.measure_continuously)
        self.worker.finished.connect(self.thread.quit)
        self.worker.finished.connect(self.worker.deleteLater)
        self.thread.finished.connect(self.thread.deleteLater)
        self.worker.progress.connect(self.updatePlots)
        self.thread.start()

    def measure_preset_task(self):        
        if self.worker:
            if self.worker.is_running:
                self.worker.is_running = False
                return
                
        global fnames
        global descriptions
        timestamp = time.strftime('%Y%m%d_%H%M')
        fnames = [os.path.join(self.working_directory_line.text(), timestamp + ".dat") for i in [0]]

        settings['working_directory'] = self.working_directory_line.text()
        save_settings()

        self.plot_data = {}
        self.y_lims = [[0,-np.inf], [np.inf,-np.inf]]

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
        
        delay = info["delay"]
        spectra = info['spectra']
        out_signal = info['out_signal']
        out_bckg = info['out_bckg']
        wavelength_axis = info['wavelength_axis']
        wavelength_range = info['wavelength_range']        

        self.plot_data[delay] = (out_signal[0]['not_pumped'], spectra[0])

        x = [delay for delay in self.plot_data]
        y1 = [self.plot_data[item][0] for item in self.plot_data]
        y2 = [self.plot_data[item][1] for item in self.plot_data]

        y1_sorted = [y for _, y in sorted(zip(x, y1), key=lambda pair: pair[0])]
        y2_sorted = [y for _, y in sorted(zip(x, y2), key=lambda pair: pair[0])]
        x_sorted = sorted(x)

        delay_range = [x_sorted[0] - 1, x_sorted[-1] + 1]
        
        self.sc[0].lines_wlc[0].set_xdata(x_sorted)
        self.sc[0].lines_wlc[0].set_ydata(y1_sorted)

        self.sc[0].lines_pp[0].set_xdata(x_sorted)
        self.sc[0].lines_pp[0].set_ydata(y2_sorted)            

        self.sc[0].lines_pp[0].set_label('{:}, $T={:.2f}$ ps'.format(descriptions[0], info.get('delay') or 0.0))
            
        self.y_lims[0][0] = 0
        if np.nanmax(y1_sorted) > self.y_lims[0][1]:
            self.y_lims[0][1] = np.nanmax(y1_sorted)
            
        if np.nanmin(y2_sorted) < self.y_lims[1][0]:
            self.y_lims[1][0] = np.nanmin(y2_sorted)
        if np.nanmax(y2_sorted) > self.y_lims[1][1]:
            self.y_lims[1][1] = np.nanmax(y2_sorted)

        try:
            self.sc[0].ax1.set_xlim(delay_range)    
            self.sc[0].ax1.set_ylim(self.y_lims[0])
            #self.sc[0].ax1.legend()
                
            self.sc[0].ax2.set_xlim(delay_range)
            self.sc[0].ax2.set_ylim(self.y_lims[1])
            #self.sc[0].ax2.legend()
        except ValueError:
            pass

                
        self.sc[0].draw()
        app.processEvents()
        self.canDraw = True

    @pyqtSlot()
    def measure_continuously_button_on_click(self):
        self.measure_continuously_task()

    @pyqtSlot()
    def measure_preset_button_on_click(self):        
        self.measure_preset_task()
        
app = QApplication([])
w = MainWindow('HARPIA-TB Pump-Probe Measurement')
app.exec_()
        