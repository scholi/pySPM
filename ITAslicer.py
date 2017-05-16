import sys
from PyQt5.QtWidgets import QApplication, QMainWindow, QTableWidgetItem, QFileDialog, QWidget
from PyQt5.QtGui import QPixmap
from PyQt5.QtCore import Qt
from slicer import Ui_slicer
import os
import re
import pySPM
import numpy as np

class SlicerApp(QWidget):
    def __init__(self, filename=None):
        super(QWidget, self).__init__()
        self.ui = Ui_slicer()
        self.ui.setupUi(self)
        
        self.canvas = self.ui.mpl.canvas
        self.fig = self.ui.mpl.canvas.fig
        self.initPlotLayout()
        
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
        self.curs = [0,0,0]
        self.volume = None
        self.ITA = pySPM.ITA(self.path) 
        for i,x in enumerate(self.ITA.get_masses()):
            self.ui.peakList.setRowCount(i+1)
            self.ui.peakList.setItem(i, 0, QTableWidgetItem(x['assign']))
            self.ui.peakList.setItem(i, 1, QTableWidgetItem("{:.2f}u".format((x['cmass']))))
            self.ui.peakList.setItem(i, 2, QTableWidgetItem("{:.2f}u".format(x['umass'] - x['lmass'])))
        self.ui.peakList.show()
        self.ui.peakList.cellClicked.connect(self.load_channel)
        self.canvas.mpl_connect('button_press_event', self.on_pick)
        
    def load_channel(self, row, col):
        id = row
        vol = []
        for i in range(self.ITA.Nscan):
            x = self.ITA.getImage(id,i)
            vol.append(x)
        self.volume = np.stack([x for x in vol],axis=2)
        self.plot()
        
    def plot(self):
        self.axXY.clear()
        self.axXZ.clear()
        self.axYZ.clear()
        if self.volume == None:
            return
        self.axXY.imshow(self.volume[:,:,self.curs[2]])
        self.axXZ.imshow(self.volume[self.curs[1],:,:].T)
        self.axYZ.imshow(self.volume[:,self.curs[0],:].T)
        self.axXY.axhline(self.curs[1])
        self.axXY.axvline(self.curs[0])
        self.axXZ.axhline(self.curs[2])
        self.axXZ.axvline(self.curs[0])
        self.axYZ.axhline(self.curs[2])
        self.axYZ.axvline(self.curs[1])
        self.axXY.set_title("XY")
        self.axXZ.set_title("XZ")
        self.axYZ.set_title("YZ")
        self.canvas.draw()
        
    def on_pick(self, event):
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
window = SlicerApp(r"C:\Users\ols\ownCloud\ToFSIMS\Micropat_auto__ @micropat2 (-)_2.ita")
window.show()
sys.exit(app.exec_())