#!/usr/bin/env python


import sys, os, random
from PyQt4 import QtGui, QtCore



class ApplicationWindow(QtGui.QMainWindow):
    def __init__(self):
        QtGui.QMainWindow.__init__(self)

        self.main_widget = QtGui.QWidget(self)


        self.main_widget.setFocus()
        self.setCentralWidget(self.main_widget)



#qApp = QtGui.QApplication(sys.argv)

aw = ApplicationWindow()
aw.show()
#sys.exit(qApp.exec_())
#qApp.exec_()
