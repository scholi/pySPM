# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'C:\Users\ols\Dropbox\Projects\ITA_slicer\slicer.ui'
#
# Created by: PyQt5 UI code generator 5.6
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_slicer(object):
    def setupUi(self, slicer):
        slicer.setObjectName("slicer")
        slicer.resize(802, 518)
        self.gridLayout = QtWidgets.QGridLayout(slicer)
        self.gridLayout.setObjectName("gridLayout")
        self.splitter = QtWidgets.QSplitter(slicer)
        self.splitter.setOrientation(QtCore.Qt.Horizontal)
        self.splitter.setObjectName("splitter")
        self.peakList = QtWidgets.QTableWidget(self.splitter)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.peakList.sizePolicy().hasHeightForWidth())
        self.peakList.setSizePolicy(sizePolicy)
        self.peakList.setObjectName("peakList")
        self.peakList.setColumnCount(3)
        self.peakList.setRowCount(0)
        item = QtWidgets.QTableWidgetItem()
        self.peakList.setHorizontalHeaderItem(0, item)
        item = QtWidgets.QTableWidgetItem()
        self.peakList.setHorizontalHeaderItem(1, item)
        item = QtWidgets.QTableWidgetItem()
        self.peakList.setHorizontalHeaderItem(2, item)
        self.mpl = MplWidget(self.splitter)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.mpl.sizePolicy().hasHeightForWidth())
        self.mpl.setSizePolicy(sizePolicy)
        self.mpl.setMinimumSize(QtCore.QSize(500, 500))
        self.mpl.setObjectName("mpl")
        self.gridLayout.addWidget(self.splitter, 0, 0, 1, 1)
        self.mpl.raise_()
        self.peakList.raise_()

        self.retranslateUi(slicer)
        QtCore.QMetaObject.connectSlotsByName(slicer)

    def retranslateUi(self, slicer):
        _translate = QtCore.QCoreApplication.translate
        slicer.setWindowTitle(_translate("slicer", "Slicer"))
        item = self.peakList.horizontalHeaderItem(0)
        item.setText(_translate("slicer", "Name"))
        item = self.peakList.horizontalHeaderItem(1)
        item.setText(_translate("slicer", "center mass"))
        item = self.peakList.horizontalHeaderItem(2)
        item.setText(_translate("slicer", "Δ mass"))

from mplwidget import MplWidget
