# Copyright 2018 Olivier Scholder <o.scholder@gmail.com>

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from PyQt5.QtCore import QSize
from PyQt5.QtWidgets import QSizePolicy, QVBoxLayout, QWidget


class MplCanvas(FigureCanvas):
    def __init__(self):
        self.fig = Figure()
        FigureCanvas.__init__(self, self.fig)
        FigureCanvas.setSizePolicy(self, QSizePolicy.Expanding, QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)

    def sizeHint(self):
        w, h = self.get_width_height()
        return QSize(w, h)


class MplWidget(QWidget):
    def __init__(self, parent=None):
        QWidget.__init__(self, parent)
        self.canvas = MplCanvas()
        self.mpl_toolbar = NavigationToolbar(self.canvas, self)
        layout = QVBoxLayout()
        self.setLayout(layout)
        layout.addWidget(self.mpl_toolbar)
        layout.addWidget(self.canvas)
