import sys
import os
from PyQt5.QtWidgets import QApplication, QWidget, QVBoxLayout, QMainWindow, \
    QAction, qApp, QFileDialog, QComboBox
from minITM import ITM
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
            print(I.Nscan)
            t,MeasData,s = I.get_meas_data()
            
            ax.plot(t, MeasData*1e6)
            ax.set_xlabel("Time [s]");
            ax.set_ylabel("Emission Current [$\\mu$A]")
            if scans:
                assert type(scans) is int
                lim = ax.get_ylim()
                axs = ax.twiny()
                axs.set_xticks(.5*s[:-1:scans]+.5*s[1::scans])
                axs.set_xticklabels([str(i+1) for i in range(0,I.Nscan,scans)])
                colors = [i%2 for i in range(0,I.Nscan,scans)]
                for i,tick in enumerate(axs.xaxis.get_ticklabels()):
                    tick.set_color(["black","green"][colors[i]])
                axs.set_xlim(ax.get_xlim())
                for i in range(1,I.Nscan-1,2):
                    ax.fill_between([s[i],s[i+1]],*lim,color='green',alpha=.1)
                axs.set_xlabel("Scan number")
                axs.set_xlim(ax.get_xlim())
                axs.set_ylim(*lim)
            axb = DualPlot(ax)
            t, MeasData,_ = I.get_meas_data("Instrument.LMIG.Suppressor")
            axb.plot(t, MeasData, 'C1')
            axb.set_ylabel("Suppressor Voltage [V]")
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

app = QApplication(sys.argv)
a = Plotter()
a.open()  #r"Z:\Olivier\180328_Xue")
sys.exit(app.exec_())