# -- coding: utf-8 --

# Copyright 2018 Olivier Scholder <o.scholder@gmail.com>

"""
This is a standalone script which allows the user to perfome cross-section on ToF-SIMS images on different channels
"""

import sys
from PyQt5.QtWidgets import QApplication, QMainWindow, QTableWidgetItem, QFileDialog, QWidget, QAction, QProgressBar, QStatusBar
from PyQt5.QtGui import QPixmap
from PyQt5.QtCore import Qt, QCoreApplication
from slicer import Ui_slicer
import os
import re
import pySPM
import numpy as np
from pySPM.utils import CDF
from scipy import optimize as opt

class SlicerApp(QWidget):
    def __init__(self, filename=None):
        super(QWidget, self).__init__()
        self.ui = Ui_slicer()
        self.ui.setupUi(self)
        
        self.canvas = self.ui.mpl.canvas
        self.fig = self.ui.mpl.canvas.fig
        self.initPlotLayout()
        self.level = None
        if filename is None:
            if len(sys.argv) < 2:
                self.path, _ = QFileDialog.getOpenFileName(self,"ITA Image", "","(*.ITA)")
            else:
                # If an argument is sent to the script, the first argument will be used as a Path. Very usefull for debugging the script without having to selectr the folder each time with window dialog
                self.path = sys.argv[1]
        else:
            self.path = filename
        if not os.path.exists(self.path):
            raise Exception("File \"{}\" is not found".format(self.path))
        if os.path.exists(self.path+".level.npy"):
            self.level = np.load(self.path+".level.npy")
        self.curs = [0,0,0]
        self.volume = None
        self.ITA = pySPM.ITA(self.path) 
        for i,x in enumerate(self.ITA.get_masses()):
            self.ui.peakList.setRowCount(i+1)
            self.ui.peakList.setItem(i, 0, QTableWidgetItem(x['assign']))
            self.ui.peakList.setItem(i, 1, QTableWidgetItem("{:.2f}u".format((x['cmass']))))
            self.ui.peakList.setItem(i, 2, QTableWidgetItem("{:.2f}u".format(x['umass'] - x['lmass'])))
        self.ui.peakList.show()
        self.ui.cmap.currentIndexChanged.connect(self.plot)
        self.ui.prof1daxis.currentIndexChanged.connect(self.plot)
        self.ui.peakList.cellClicked.connect(self.load_channel)
        self.canvas.mpl_connect('button_press_event', self.on_pick)
        self.flatAction = QAction("Flatten substrate from this channel")
        self.flatAction.triggered.connect(self.flatten)
        self.ui.correction.stateChanged.connect(self.plot)
        self.ui.peakList.addAction(self.flatAction)
        self.ui.status.setText("IDLE")
    
    def flatline(self, y):
        for x in range(self.volume.shape[1]):
            popt, pcov = opt.curve_fit(CDF, np.arange(self.volume.shape[2]), self.volume[y,x,:], (10, self.volume.shape[2]/2, 1))
            self.level[y, x] = popt[1]
         
    def flatten(self):
        from scipy import optimize as opt
        self.ui.status.setText("Start the flattening...")
        self.level = np.zeros(self.volume.shape[:2])
        self.ui.pb.setMaximum(self.volume.shape[0])
        for y in range(self.volume.shape[0]):
            self.ui.pb.setValue(y)
            self.flatline(y)
            QCoreApplication.processEvents()
            self.ax.clear()
            self.ax.imshow(self.level)
            self.canvas.draw()
        self.ui.pb.setValue(0)
        self.ui.status.setText("Flattening finished")
        np.save(self.path+".level", self.level)
        
    def load_channel(self, row, col):
        self.ui.status.setText("Loading channel...")
        id = row
        vol = []
        for i in range(self.ITA.Nscan):
            x = self.ITA.getImage(id, i)
            vol.append(x)
        self.volume = np.stack([x for x in vol], axis=2)
        if not self.level is None:
            self.corrected = np.zeros(self.volume.shape)
            z = np.arange(self.ITA.Nscan)
            self.ui.pb.setMaximum(self.level.shape[0])
            for y in np.arange(self.level.shape[0]):
                self.ui.pb.setValue(y)
                for x in np.arange(self.level.shape[1]):
                    dz = int(-self.level[y,x] + np.max(self.level))
                    self.corrected[y,x,dz:] = self.volume[y,x,:self.volume.shape[2]-dz]
            self.ui.pb.setValue(0)
        self.plot()
        self.ui.status.setText("IDLE")
        
    def plot(self):
        if self.ui.correction.isChecked():
            A = self.corrected
        else:
            A = self.volume
        cmap = self.ui.cmap.currentText()
        self.axXY.clear()
        self.axXZ.clear()
        self.axYZ.clear()
        self.ax.clear()
        if self.volume is None:
            return
        self.axXY.imshow(A[:,:,self.curs[2]],cmap=cmap)
        self.axXZ.imshow(A[self.curs[1],:,:].T,cmap=cmap)
        self.axYZ.imshow(A[:,self.curs[0],:].T,cmap=cmap)
        self.axXY.axhline(self.curs[1])
        self.axXY.axvline(self.curs[0])
        self.axXZ.axhline(self.curs[2])
        self.axXZ.axvline(self.curs[0])
        self.axYZ.axhline(self.curs[2])
        self.axYZ.axvline(self.curs[1])
        self.axXY.set_title("XY")
        self.axXZ.set_title("XZ")
        self.axYZ.set_title("YZ")
        if self.ui.prof1daxis.currentText() in 'XYZ':
            i = 'XYZ'.index(self.ui.prof1daxis.currentText())
            self.ax.set_xlabel("XYZ"[i])
            if i==0:
                self.ax.plot(A[:,self.curs[0],self.curs[2]])
            elif i==1:
                self.ax.plot(A[self.curs[1],:,self.curs[2]])
            elif i==2:
                self.ax.plot(A[self.curs[1],self.curs[0],:])
        self.canvas.draw()
        
    def on_pick(self, event):
        if not event.inaxes in [self.axYZ,self.axXZ,self.axXY]:
            return
        x = event.xdata
        y = event.ydata
        axis = [self.axYZ,self.axXZ,self.axXY].index(event.inaxes)
        xdata = int(x)
        ydata = int(y)
        
        if event.inaxes == self.axXY:
            self.curs[0] = xdata
            self.curs[1] = ydata
        elif event.inaxes == self.axYZ:
            self.curs[1] = xdata
            self.curs[2] = ydata
        elif event.inaxes == self.axXZ:
            self.curs[0] = xdata
            self.curs[2] = ydata
        else:
            print("Click event not handled")
        self.curs = [np.clip(0,self.curs[i],self.volume.shape[i]-1) for i in range(3)]
        self.plot()
        
    def initPlotLayout(self):
        """
        Setup the plotting layout.
        """
        self.axXY = self.fig.add_subplot(2,2,1)
        self.axYZ = self.fig.add_subplot(2,2,2)
        self.axXZ = self.fig.add_subplot(2,2,3)
        self.ax = self.fig.add_subplot(2,2,4)
        self.fig.tight_layout()

app = QApplication(sys.argv)
window = SlicerApp()
window.show()
sys.exit(app.exec_())