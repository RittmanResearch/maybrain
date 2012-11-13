# -*- coding: utf-8 -*-
"""
Created on Wed Oct 31 21:00:09 2012

@author: hypocrates
"""

# First, and before importing any Enthought packages, set the ETS_TOOLKIT
# environment variable to qt4, to tell Traits that we will use Qt.


#import os
#os.environ['ETS_TOOLKIT'] = 'qt4'
import mayBrainTools as mb

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
        self.ui.brainPlot.setEnabled(0)
        self.ui.skullPlot.setEnabled(0)

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
        QtCore.QObject.connect(self.ui.brainLoad, QtCore.SIGNAL('clicked()'), self.loadBrain)
        QtCore.QObject.connect(self.ui.skullLoad, QtCore.SIGNAL('clicked()'), self.loadSkull)
        QtCore.QObject.connect(self.ui.brainPlot, QtCore.SIGNAL('clicked()'), self.plotBrain)
        QtCore.QObject.connect(self.ui.skullPlot, QtCore.SIGNAL('clicked()'), self.plotSkull)
        
        
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
        ''' load a brain using given filenames '''
        # get adjacency filename
        f = str(self.ui.adjFilename.text())
        
        # get threshold
        try:
            thold = float(self.ui.thresholdValue.text())
        except:
            print('threshold value not recovered, is it a float?', self.ui.thresholdValue.text())
        
        # get spatial info file
        g = str(self.ui.spatialFilename.text())

        # create brain and load adjacency file
#        try:
            # create brain object if it doesn't exist
        if not('mainBrain' in self.brains):
            br = mb.brainObj()
            self.brains['mainBrain'] = br
        else:
            br = self.brains['mainBrain']
        # read in files
        br.readAdjFile(f, thold)
        br.readSpatialInfo(g)
                
        # enable plot button
        self.ui.brainPlot.setEnabled(True)
        
#        except:
#            print('could not create brain object and load in files')
            
            
        
    def loadSkull(self):
        ''' load a skull file '''
        # get filename
        f = str(self.ui.skullFilename.text())

        try:        
            # create brain object if it doesn't exist
            if not('mainBrain' in self.brains):
                br = mb.brainObj()
                self.brains['mainBrain'] = br
            else:
                br = self.brains['mainBrain']
            # read in file
            br.importSkull(f)
            
            self.ui.skullPlot.setEnabled(True)
            
        except:
            print('problem loading skull file')           
        

    ## =============================================
    
    ## plotting functions
    
    def plotBrain(self):
        ''' plot the brain object '''
        
        # check if brain object exists
        if not('mainBrain' in self.brains):
            return

        # plot the brain
        try:            
            self.plot.plotBrain(self.brains['mainBrain'], label = 'mainBrain')
        except:
            print('problem plotting brain, have files been loaded?')
                
    
    def plotSkull(self):
        ''' plot the skull '''
        
        # check if brain object exists
        if not('mainBrain' in self.brains):
            return
            
        # plot
        try:
            self.plot.plotSkull(self.brains['mainBrain'], label = 'mainSkull')
        except:
            print('could not plot skull, has file been loaded?')


    ## =============================================

    ## subbrain functions
    def plotSubBrain(self):
        ''' plot a subBrain with certain properties '''
        a=1
        # create the subBrain
        
        
        # pass info to different parts of the program
        
        
        # do the plot

                        
        

        
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