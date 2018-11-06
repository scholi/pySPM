#!python
import sys
import os
from PyQt5.QtWidgets import QApplication, QWidget, QVBoxLayout, QMainWindow, \
    QAction, qApp, QFileDialog, QComboBox
from pySPM import ITM
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar

def DualPlot(ax, col1='C0', col2='C1'):
    axb = ax.twinx()
    axb.spines['left'].set_color(col1)
    axb.spines['right'].set_color(col2)
    ax.yaxis.label.set_color(col1)
    axb.yaxis.label.set_color(col2)
    ax.tick_params(axis='y', colors=col1)
    axb.tick_params(axis='y', colors=col2)
    return axb
    
class Plotter(QMainWindow):
    def plot(self, ev=None, scans=3):
        if self.cwd is not None:
            self.figure.clf()
            ax = self.figure.add_subplot(1,1,1)
            
            filename = os.path.join(self.cwd, self.fileDrop.currentText())
            I = ITM(filename)
            I.show_stability(ax)
            del I
            self.canvas.draw()
            
    def open(self, dirs=None):
        if dirs is None:
            dirs = QFileDialog.getExistingDirectory(self)
        if dirs:
            self.cwd = dirs
            self.fileDrop.clear()
            for x in os.listdir(dirs):
                if x.endswith(".itm") or x.endswith(".ITM"):
                    self.fileDrop.addItem(x)

    def __init__(self):
        QWidget.__init__(self)
        self.cwd = None
        self.figure = Figure()
        self.canvas = FigureCanvas(self.figure)
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.markers = None
        self.draggable = None
        self.img = None
        self.msize = 6

        openAction = QAction('&Open', self)        
        openAction.setShortcut('Ctrl+O')
        openAction.setStatusTip('Open folder')
        openAction.triggered.connect(self.open)
        
        menubar = self.menuBar()
        fileMenu = menubar.addMenu('&File')
        fileMenu.addAction(openAction)
        
        self.fileDrop = QComboBox()
        
        layout = QVBoxLayout()
        layout.addWidget(self.fileDrop)
        layout.addWidget(self.toolbar)
        layout.addWidget(self.canvas)
        
        window = QWidget()
        window.setLayout(layout);
        self.setCentralWidget(window)
        
        self.fileDrop.currentIndexChanged.connect(self.plot)
        
        self.show()

def main():
    app = QApplication(sys.argv)
    a = Plotter()
    a.open()
    sys.exit(app.exec_())

if __name__ == '__main__':
    main()