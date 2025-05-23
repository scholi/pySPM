# Form implementation generated from reading ui file 'spectraviewer.ui'
#
# Created by: PyQt5 UI code generator 5.9
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtWidgets

from .mplwidget import MplWidget


class Ui_SpectraViewer:
    def setupUi(self, SpectraViewer):
        SpectraViewer.setObjectName("SpectraViewer")
        SpectraViewer.resize(1148, 690)
        self.centralwidget = QtWidgets.QWidget(SpectraViewer)
        self.centralwidget.setObjectName("centralwidget")
        self.gridLayout = QtWidgets.QGridLayout(self.centralwidget)
        self.gridLayout.setObjectName("gridLayout")
        self.pushButton_2 = QtWidgets.QPushButton(self.centralwidget)
        sizePolicy = QtWidgets.QSizePolicy(
            QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding
        )
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.pushButton_2.sizePolicy().hasHeightForWidth())
        self.pushButton_2.setSizePolicy(sizePolicy)
        self.pushButton_2.setMaximumSize(QtCore.QSize(18, 16777215))
        self.pushButton_2.setObjectName("pushButton_2")
        self.gridLayout.addWidget(self.pushButton_2, 1, 2, 4, 1)
        self.lab_sf = QtWidgets.QLabel(self.centralwidget)
        self.lab_sf.setObjectName("lab_sf")
        self.gridLayout.addWidget(self.lab_sf, 2, 0, 1, 2)
        self.tableMassCal = QtWidgets.QTableWidget(self.centralwidget)
        sizePolicy = QtWidgets.QSizePolicy(
            QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding
        )
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.tableMassCal.sizePolicy().hasHeightForWidth())
        self.tableMassCal.setSizePolicy(sizePolicy)
        self.tableMassCal.setSelectionMode(QtWidgets.QAbstractItemView.SingleSelection)
        self.tableMassCal.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
        self.tableMassCal.setColumnCount(4)
        self.tableMassCal.setObjectName("tableMassCal")
        self.tableMassCal.setRowCount(0)
        item = QtWidgets.QTableWidgetItem()
        self.tableMassCal.setHorizontalHeaderItem(0, item)
        item = QtWidgets.QTableWidgetItem()
        self.tableMassCal.setHorizontalHeaderItem(1, item)
        item = QtWidgets.QTableWidgetItem()
        self.tableMassCal.setHorizontalHeaderItem(2, item)
        item = QtWidgets.QTableWidgetItem()
        self.tableMassCal.setHorizontalHeaderItem(3, item)
        self.tableMassCal.horizontalHeader().setDefaultSectionSize(64)
        self.gridLayout.addWidget(self.tableMassCal, 1, 0, 1, 2)
        self.lab_k0 = QtWidgets.QLabel(self.centralwidget)
        self.lab_k0.setObjectName("lab_k0")
        self.gridLayout.addWidget(self.lab_k0, 3, 0, 1, 2)
        self.pushButton = QtWidgets.QPushButton(self.centralwidget)
        self.pushButton.setObjectName("pushButton")
        self.gridLayout.addWidget(self.pushButton, 4, 0, 1, 2)
        self.mpl = MplWidget(self.centralwidget)
        sizePolicy = QtWidgets.QSizePolicy(
            QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding
        )
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.mpl.sizePolicy().hasHeightForWidth())
        self.mpl.setSizePolicy(sizePolicy)
        self.mpl.setStyleSheet("")
        self.mpl.setObjectName("mpl")
        self.gridLayout.addWidget(self.mpl, 1, 3, 4, 1)
        self.lab_m0 = QtWidgets.QLabel(self.centralwidget)
        self.lab_m0.setObjectName("lab_m0")
        self.gridLayout.addWidget(self.lab_m0, 0, 3, 1, 1)
        self.show_mass = QtWidgets.QCheckBox(self.centralwidget)
        self.show_mass.setEnabled(False)
        self.show_mass.setCheckable(True)
        self.show_mass.setChecked(True)
        self.show_mass.setObjectName("show_mass")
        self.gridLayout.addWidget(self.show_mass, 0, 0, 1, 2)
        SpectraViewer.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(SpectraViewer)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 1148, 21))
        self.menubar.setObjectName("menubar")
        SpectraViewer.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(SpectraViewer)
        self.statusbar.setObjectName("statusbar")
        SpectraViewer.setStatusBar(self.statusbar)

        self.retranslateUi(SpectraViewer)
        QtCore.QMetaObject.connectSlotsByName(SpectraViewer)

    def retranslateUi(self, SpectraViewer):
        _translate = QtCore.QCoreApplication.translate
        SpectraViewer.setWindowTitle(_translate("SpectraViewer", "MainWindow"))
        self.pushButton_2.setText(_translate("SpectraViewer", "«"))
        self.lab_sf.setText(_translate("SpectraViewer", "sf"))
        item = self.tableMassCal.horizontalHeaderItem(0)
        item.setText(_translate("SpectraViewer", "Element"))
        item = self.tableMassCal.horizontalHeaderItem(1)
        item.setText(_translate("SpectraViewer", "Mass"))
        item = self.tableMassCal.horizontalHeaderItem(2)
        item.setText(_translate("SpectraViewer", "Time"))
        item = self.tableMassCal.horizontalHeaderItem(3)
        item.setText(_translate("SpectraViewer", "Delta [u]"))
        self.lab_k0.setText(_translate("SpectraViewer", "k0"))
        self.pushButton.setText(_translate("SpectraViewer", "Delete Element"))
        self.lab_m0.setText(_translate("SpectraViewer", "Center Mass:"))
        self.show_mass.setText(_translate("SpectraViewer", "Show masses"))
