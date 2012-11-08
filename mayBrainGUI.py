# -*- coding: utf-8 -*-
"""
Created on Wed Oct 31 21:00:09 2012

@author: hypocrates
"""

# First, and before importing any Enthought packages, set the ETS_TOOLKIT
# environment variable to qt4, to tell Traits that we will use Qt.


#import os
#os.environ['ETS_TOOLKIT'] = 'qt4'
import networkutils_binary_camgrid_0_2 as mb

# To be able to use PySide or PyQt4 and not run in conflicts with traits,
# we need to import QtGui and QtCore from pyface.qt
from pyface.qt import QtGui, QtCore
#import sip
#sip.setapi('QString', 2)


import sys
#from PyQt4 import QtCore, QtGui
import mayBrainUI as ui


class mayBrainGUI(QtGui.QMainWindow):
    def __init__(self, parent=None):
        
        QtGui.QMainWindow.__init__(self, parent)
        # set up UI
        self.ui = ui.Ui_MainWindow()
        self.ui.setupUi(self)

        # set up function
        self.brains = {}
        self.plot = mb.plotObj()
        
        # link up buttons and functions
        self.connectSlots()
        
        
        
    def connectSlots(self):
        ''' connect buttons and functions '''

        # quit app        
        QtCore.QObject.connect(self.ui.quitButton, QtCore.SIGNAL('clicked()'), app.quit)
        
        # tab 1: general settings
        QtCore.QObject.connect(self.ui.spatialFnameButton, QtCore.SIGNAL('clicked()'), self.getSpatialFilename)
        QtCore.QObject.connect(self.ui.adjFnameButton, QtCore.SIGNAL('clicked()'), self.getAdjFilename)
        QtCore.QObject.connect(self.ui.skullFnameButton, QtCore.SIGNAL('clicked()'), self.getSkullFilename)
        QtCore.QObject.connect(self.ui.propsFnameButton, QtCore.SIGNAL('clicked()'), self.getPropsFilename)                
        
        
    ## =============================================

    ## Functions to connect slots and Maybrain function
        
    # Basic info, files and loading
    def getSpatialFilename(self):
        ''' open a dialog to get the spatial filename'''
        f = QtGui.QFileDialog.getOpenFileName(self, 'Choose spatial file', '.')
        self.ui.spatialFilename.setText(f)
        
    def getAdjFilename(self):
        ''' open a dialog to get the adjacency filename '''
        f = QtGui.QFileDialog.getOpenFileName(self, 'Choose adjacency file', '.')
        self.ui.adjFilename.setText(f)
        
    def getSkullFilename(self):
        ''' get a skull filename '''
        f = QtGui.QFileDialog.getOpenFileName(self, 'Choose skull file', '.')
        self.ui.skullFilename.setText(f)
        
    def getPropsFilename(self):
        f = QtGui.QFileDialog.getOpenFileName(self, 'Choose properties file', '.')
        self.ui.propsFilename.setText(f)
        
    def loadBrain(self):
        ''' load a spatial info file '''
        f = str(self.ui.fileNameBox.text())
        
        br = mb.brainObj()
        try:
            thold = self.ui.thresholdValue.test()
        except:
            print('threshold value not recovered, is it a ')
        br.readAdjFile(f, )
        
        
        

        
if __name__ == "__main__":

    # using instance allows Mayavi to run alongside without any problems.
    app = QtGui.QApplication.instance()
    ex = mayBrainGUI()
    ex.show()
    sys.exit(app.exec_())

    

    # create and show app    
#    app = QtGui.QApplication(sys.argv)
#    myapp = mayBrainGUI()
##    myapp.show()
#    sys.exit(app.exec_())        