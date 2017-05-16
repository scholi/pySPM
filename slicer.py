# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'C:\Users\ols\Dropbox\Python\pySPM\slicer.ui'
#
# Created by: PyQt5 UI code generator 5.6
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_slicer(object):
    def setupUi(self, slicer):
        slicer.setObjectName("slicer")
        slicer.resize(802, 546)
        self.gridLayout = QtWidgets.QGridLayout(slicer)
        self.gridLayout.setObjectName("gridLayout")
        self.peakList = QtWidgets.QTableWidget(slicer)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.peakList.sizePolicy().hasHeightForWidth())
        self.peakList.setSizePolicy(sizePolicy)
        self.peakList.setContextMenuPolicy(QtCore.Qt.ActionsContextMenu)
        self.peakList.setObjectName("peakList")
        self.peakList.setColumnCount(3)
        self.peakList.setRowCount(0)
        item = QtWidgets.QTableWidgetItem()
        self.peakList.setHorizontalHeaderItem(0, item)
        item = QtWidgets.QTableWidgetItem()
        self.peakList.setHorizontalHeaderItem(1, item)
        item = QtWidgets.QTableWidgetItem()
        self.peakList.setHorizontalHeaderItem(2, item)
        self.gridLayout.addWidget(self.peakList, 0, 0, 1, 1)
        self.mpl = MplWidget(slicer)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.mpl.sizePolicy().hasHeightForWidth())
        self.mpl.setSizePolicy(sizePolicy)
        self.mpl.setMinimumSize(QtCore.QSize(500, 500))
        self.mpl.setObjectName("mpl")
        self.gridLayout.addWidget(self.mpl, 0, 1, 1, 1)
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.status = QtWidgets.QLabel(slicer)
        self.status.setObjectName("status")
        self.horizontalLayout.addWidget(self.status)
        self.correction = QtWidgets.QCheckBox(slicer)
        self.correction.setLayoutDirection(QtCore.Qt.RightToLeft)
        self.correction.setObjectName("correction")
        self.horizontalLayout.addWidget(self.correction)
        self.gridLayout.addLayout(self.horizontalLayout, 1, 0, 1, 1)
        self.pb = QtWidgets.QProgressBar(slicer)
        self.pb.setEnabled(True)
        self.pb.setProperty("value", 0)
        self.pb.setInvertedAppearance(False)
        self.pb.setObjectName("pb")
        self.gridLayout.addWidget(self.pb, 1, 1, 1, 1)

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
        self.status.setText(_translate("slicer", "Loading..."))
        self.correction.setText(_translate("slicer", "Apply Correction"))

from mplwidget import MplWidget
