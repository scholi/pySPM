# -- coding: utf-8 --

# Copyright 2018 Olivier Scholder <o.scholder@gmail.com>

"""
This module allows to display a small GUI in order to display a table of key/values.

It is used by the class ITM.show_values(gui=True).
"""

from PyQt5.QtGui import QStandardItemModel, QStandardItem
from PyQt5.QtWidgets import QWidget, QApplication, QTreeView, QVBoxLayout


class GUI_values(QWidget):
    def __init__(self, data):
        QWidget.__init__(self)
        self.treeView = QTreeView()
        self.model = QStandardItemModel()
        self.addItems(self.model, data)
        self.treeView.setModel(self.model)
        layout = QVBoxLayout()
        layout.addWidget(self.treeView)
        self.setLayout(layout)

    def addItems(self, parent, elements):
        for k in sorted(elements.keys()):
            item = QStandardItem(k)
            parent.appendRow(item)
            if type(elements[k]) == dict:
                self.addItems(item, elements[k])
            else:
                child = QStandardItem(str(elements[k]))
                item.appendRow(child)


def show_values(data):
    app = QApplication([])
    G = GUI_values(data)
    G.show()
    app.exec_()
