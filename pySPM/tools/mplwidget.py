from PyQt5.QtWidgets import QWidget, QVBoxLayout, QSizePolicy
from PyQt5.QtCore import QSize
import matplotlib as mpl
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure

class MplCanvas(FigureCanvas):
	def __init__(self):
		self.fig = Figure()
		FigureCanvas.__init__(self, self.fig)
		FigureCanvas.setSizePolicy(self,
			QSizePolicy.Expanding,
			QSizePolicy.Expanding)
		FigureCanvas.updateGeometry(self)

	def sizeHint(self):
		w, h = self.get_width_height()
		return QSize(w,h)

class MplWidget(QWidget):
		def __init__(self, parent = None):
			QWidget.__init__(self, parent)
			self.canvas = MplCanvas()
			self.mpl_toolbar = NavigationToolbar(self.canvas, self)
			layout = QVBoxLayout()
			self.setLayout(layout)
			layout.addWidget(self.mpl_toolbar)
			layout.addWidget(self.canvas)
