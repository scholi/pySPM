# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'timer_display.ui'
#
# Created by: PyQt5 UI code generator 5.9
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_ToF_Timer(object):
    def setupUi(self, ToF_Timer):
        ToF_Timer.setObjectName("ToF_Timer")
        ToF_Timer.resize(466, 108)
        self.centralWidget = QtWidgets.QWidget(ToF_Timer)
        self.centralWidget.setObjectName("centralWidget")
        self.verticalLayout = QtWidgets.QVBoxLayout(self.centralWidget)
        self.verticalLayout.setContentsMargins(11, 11, 11, 11)
        self.verticalLayout.setSpacing(6)
        self.verticalLayout.setObjectName("verticalLayout")
        self.progressBar = QtWidgets.QProgressBar(self.centralWidget)
        self.progressBar.setProperty("value", 0)
        self.progressBar.setObjectName("progressBar")
        self.verticalLayout.addWidget(self.progressBar)
        self.label = QtWidgets.QLabel(self.centralWidget)
        self.label.setObjectName("label")
        self.verticalLayout.addWidget(self.label)
        self.label_2 = QtWidgets.QLabel(self.centralWidget)
        self.label_2.setObjectName("label_2")
        self.verticalLayout.addWidget(self.label_2)
        self.label_3 = QtWidgets.QLabel(self.centralWidget)
        self.label_3.setObjectName("label_3")
        self.verticalLayout.addWidget(self.label_3)
        ToF_Timer.setCentralWidget(self.centralWidget)

        self.retranslateUi(ToF_Timer)
        QtCore.QMetaObject.connectSlotsByName(ToF_Timer)

    def retranslateUi(self, ToF_Timer):
        _translate = QtCore.QCoreApplication.translate
        ToF_Timer.setWindowTitle(_translate("ToF_Timer", "ToF_Timer"))
        self.label.setText(_translate("ToF_Timer", "Remaining time: 00:00:00"))
        self.label_2.setText(_translate("ToF_Timer", "Scans: 0/0"))
        self.label_3.setText(_translate("ToF_Timer", "Analysis time: 0 s"))

