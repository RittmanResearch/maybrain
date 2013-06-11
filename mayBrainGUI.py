# -*- coding: utf-8 -*-
"""
Created on Wed Oct 31 21:00:09 2012

@author: hypocrates


Some code for QTreeWidget stuff
 1 # Let's do something interesting: load the database contents
 2 # into our task list widget
 3 for task in todo.Task.query().all():
 4     tags=','.join([t.name for t in task.tags])
 5     item=QtGui.QTreeWidgetItem([task.text,str(task.date),tags])
 6     item.task=task
 7     if task.done:
 8         item.setCheckState(0,QtCore.Qt.Checked)
 9     else:
10         item.setCheckState(0,QtCore.Qt.Unchecked)
11     self.ui.list.addTopLevelItem(item)

change log:
    2/1/2013
    - several bug fixes, including 'no selected plot' error, opacity adjustment throwing an error
    - plotBrain button becomes a replot after loading and plotting once

    9/6/13
    - added ability to make subplots

to do:

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

from numpy import log10

import sys
#from PyQt4 import QtCore, QtGui
import mayBrainUI as ui
from os import path


class mayBrainGUI(QtGui.QMainWindow):
    def __init__(self, parent=None):
        
        QtGui.QMainWindow.__init__(self, parent)
        # set up UI
        self.ui = ui.Ui_MainWindow()
        self.ui.setupUi(self)
        self.ui.adjPlot.setEnabled(0)
        self.ui.skullPlot.setEnabled(0)

        # set up function
        self.brains = {}
        self.plot = mb.plotObj()
        
        # link up buttons and functions
        self.connectSlots()
        
        # plot selected in ui.plotTree
        self.selectedPlot = None
        
        # local variables
        self.lastFolder = None # last folder viewed by user
        self.brainName = ['brain', 0]
        
        
        
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
        QtCore.QObject.connect(self.ui.skullPlot, QtCore.SIGNAL('clicked()'), self.plotSkull)
        QtCore.QObject.connect(self.ui.plotTree, QtCore.SIGNAL('itemClicked(QTreeWidgetItem*,int)'), self.setPlotValues)
        self.ui.opacitySlider.valueChanged.connect(self.setOpacity)
        self.ui.opacityBox.valueChanged.connect(self.setOpacityBox)
        self.ui.visibleCheckBox.clicked.connect(self.setVisibility)
        self.ui.redSlider.valueChanged.connect(self.setColourRed)
        self.ui.redValueBox.valueChanged.connect(self.setColourRedDial)
        self.ui.greenSlider.valueChanged.connect(self.setColourGreen)
        self.ui.greenValueBox.valueChanged.connect(self.setColourGreenDial)
        self.ui.blueSlider.valueChanged.connect(self.setColourBlue)
        self.ui.blueValueBox.valueChanged.connect(self.setColourBlueDial)
        self.ui.subPlotButtonNodes.clicked.connect(self.plotSubBrainNodes)
        self.ui.subPlotButtonEdges.clicked.connect(self.plotSubBrainEdges)
        self.ui.propsLoad.clicked.connect(self.addProperties)
        
        
#        QtCore.QObject.connect(self.ui.opacitySlider, QtCore.SIGNAL('mouseReleaseEvent()'), self.setOpacity)
        
        
    ## =============================================

    ## Functions to connect slots and Maybrain function
        
    # Basic info, files and loading
    def getSpatialFilename(self):
        ''' open a dialog to get the spatial filename'''
        f = QtGui.QFileDialog.getOpenFileName(self, 'Choose spatial file', directory = self.lastFolder)
        self.ui.spatialFilename.setText(f)
        self.lastFolder = path.dirname(f)
        
    def getAdjFilename(self):
        ''' open a dialog to get the adjacency filename '''
        f = QtGui.QFileDialog.getOpenFileName(self, 'Choose adjacency file', directory = self.lastFolder)
        self.ui.adjFilename.setText(f)
        self.lastFolder = path.dirname(f)
        
    def getSkullFilename(self):
        ''' get a skull filename '''
        f = QtGui.QFileDialog.getOpenFileName(self, 'Choose skull file', directory = self.lastFolder)
        self.ui.skullFilename.setText(f)
        self.lastFolder = path.dirname(f)
        
    def getPropsFilename(self):
        f = QtGui.QFileDialog.getOpenFileName(self, 'Choose properties file', directory = self.lastFolder)
        self.ui.propsFilename.setText(f)
        self.lastFolder = path.dirname(f)
        
    def loadBrain(self):
        ''' load a brain using given filenames '''
        self.ui.adjPlot.setEnabled(False)        
        
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
        self.brainName[1] = self.brainName[1] + 1
        brainName = self.brainName[0] + str(self.brainName[1])
        if not(brainName in self.brains):
            br = mb.brainObj()
            self.brains[brainName] = br
            # add to box for subplotting
            self.ui.subPlotBrain.addItem(brainName)
        else:
            br = self.brains[brainName]
        # read in files
        br.readAdjFile(f, threshold = thold) # need to allow different methods here
        br.readSpatialInfo(g)
                
        # enable plot button
        try:
            QtCore.QObject.disconnect(self.ui.adjPlot, QtCore.SIGNAL('clicked()'), self.rePlotBrain)
        except:
            pass
        QtCore.QObject.connect(self.ui.adjPlot, QtCore.SIGNAL('clicked()'), self.plotBrain)
        self.ui.adjPlot.setEnabled(True)
        
     
        
#        except:
#            print('could not create brain object and load in files')
            
            
        
    def loadSkull(self):
        ''' load a skull file '''
        # get filename
        f = str(self.ui.skullFilename.text())
        # get brain name
        brainName = self.brainName[0] + str(self.brainName[1])

        try:        
            # create brain object if it doesn't exist
            if not(brainName in self.brains):
                br = mb.brainObj()
                self.brains[brainName] = br
            else:
                br = self.brains[brainName]
            # read in file
            br.importSkull(f)
            
            self.ui.skullPlot.setEnabled(True)
            
        except:
            print('problem loading skull file')           
        

    ## =============================================
    
    ## plotting functions
    
    def plotBrain(self, labelIn = None):
        ''' plot the brain object '''
        
        # sort label for brain object
        if not(labelIn):
            label = self.brainName[0] + str(self.brainName[1])
        else:
            label = labelIn
        # check if brain object exists
        if not(label in self.brains):
            return

        # plot the brain
        try:            
            self.plot.plotBrain(self.brains[label], label = label)

            # add to tree view        
            QtGui.QTreeWidgetItem(self.ui.plotTree, ['brainNode', label])
            QtGui.QTreeWidgetItem(self.ui.plotTree, ['brainEdge', label])

        except:
            print('problem plotting brain, have files been loaded?')
            
        # change plot button to replot
        QtCore.QObject.disconnect(self.ui.adjPlot, QtCore.SIGNAL('clicked()'), self.plotBrain)
        QtCore.QObject.connect(self.ui.adjPlot, QtCore.SIGNAL('clicked()'), self.rePlotBrain)
    
    def rePlotBrain(self):
        ''' plot brain with altered threhold '''        
        
        if not('mainBrain' in self.brains):
            return

        # remove old plots
        self.plot.brainEdgePlots['mainBrain'].remove()
        self.plot.brainNodePlots['mainBrain'].remove()
            
#        try:
        # get new threshold
        threshold = float(self.ui.thresholdValue.text())
        # replot
        br = self.brains['mainBrain']
        br.adjMatThresholding(tVal = threshold)
        self.plot.plotBrain(br, label = 'mainBrain')       
#        except:
#            print('problem plotting brain, is threshold correct?')
            
            
    def plotSkull(self):
        ''' plot the skull '''
        
        # check if brain object exists
        if not('mainBrain' in self.brains):
            return
            
        # plot
        try:
            self.plot.plotSkull(self.brains['mainBrain'], label = 'skull')
            
            # add to treeview
            QtGui.QTreeWidgetItem(self.ui.plotTree, ['skull', 'skull'])     
            
        except:
            print('could not plot skull, has file been loaded?')
            
        


    ## =============================================

    ## subbrain functions
    def plotSubBrainNodes(self):
        self.plotSubBrain('nodes')
        
    def plotSubBrainEdges(self):
        label = self.plotSubBrain('edges')
        
        # toggle visibility of nodes
        # THIS DOESN'T WORK FOR SOME REASON...
        self.selectedPlot = label
        self.selectedPlotType = 'brainNode'
        self.setVisibility()
        
    
    def plotSubBrain(self, nodesOrEdges = 'nodes'):
        ''' plot a subBrain with certain properties '''
        
        # get values from ui
        prop = str(self.ui.subPlotProp.text())
        val = str(self.ui.subPlotValue.text())
        val2 = str(self.ui.subPlotValue_2.text())
        sb = str(self.ui.subPlotBrain.currentText())
        
        # need to add check for self.ui.subPlotValue_2 box, to define a range of values
        if val2!='':
            newName = sb+prop+val+val2
            val = {'min':float(val), 'max':float(val2)}            
        else:            
            # create the subBrain
            newName = sb+prop+val
        
        if nodesOrEdges == 'nodes':
            self.brains[newName] = self.brains[sb].makeSubBrain(prop, val)
        else:
            self.brains[newName] = self.brains[sb].makeSubBrainEdges(prop, val)
        
        # do the plot
        self.plotBrain(labelIn = newName)
        
        # add to list of plots for subBrains
        self.ui.subPlotBrain.addItem(newName)
        
        return newName
        
        
        
        
    ## =============================================
    
    ## setting and altering plot properties e.g. visibility, opacity
    
    def addProperties(self):
        ''' add properties to a brain from file '''
        
        # get properties filename from GUI
        fname = str(self.ui.propsFilename.text())        
        
        # call function from mayBrainTools to add properties to brain
        brain = (self.ui.subPlotBrain.currentText())
        self.brains[brain].nodePropertiesFromFile(fname)
        
    
    def setPlotValues(self, item):
        ''' update all values for sliders and stuff of a specific plot '''
        
        # set the selected plot
        self.selectedPlot = str(item.text(1))
        self.selectedPlotType = str(item.text(0))
        
        # set values related to plot on sliders
        props = ['opacity', 'visibility', 'colour']
        
        for p in props:
            # get the property
            v = self.plot.getPlotProperty(self.selectedPlotType, p, self.selectedPlot)
            
            # pass values to gui
            if p == 'visibility':
                if v == 1:
                    v = 2
                self.ui.visibleCheckBox.setCheckState(v)
            elif p == 'opacity':
                v = log10(9.*v+1.)*100.
                print(v)
                self.ui.opacitySlider.setValue(v)
            elif p == 'colour':
                self.ui.redSlider.setValue(v[0]*100)
                self.ui.redValueBox.setValue(v[0])
                self.ui.greenSlider.setValue(v[1]*100)
                self.ui.greenValueBox.setValue(v[1])
                self.ui.blueSlider.setValue(v[2]*100)
                self.ui.blueValueBox.setValue(v[2])
                
    def setOpacity(self):
        ''' set the opacity for the currently selected plot from the slider '''

        v = float(self.ui.opacitySlider.value())
        v = (10**(v/100.)-1.)/9.
        # update value in box
        self.ui.opacityBox.setValue(v)
        # update plot
        self.plot.changePlotProperty(self.selectedPlotType, 'opacity', self.selectedPlot, v)
        

    def setOpacityBox(self):
        ''' set the opacity from a box '''
        
        v = self.ui.opacityBox.value()
        # set value in slider
        vs = log10(9.*v+1.)*100.
        self.ui.opacitySlider.setValue(vs)
        # change value in plot
        self.plot.changePlotProperty(self.selectedPlotType, 'opacity', self.selectedPlot, v)

        
    def setVisibility(self):
        ''' toggle visibility from checkbox '''        
        
        v = self.ui.visibleCheckBox.checkState()
        if v==0:
            v=False
        elif v==2:
            v = True
        self.plot.changePlotProperty(self.selectedPlotType, 'visibility', self.selectedPlot, value = v)
        
    def setColourRed(self):
        ''' change colours from sliders '''     
        # get old values
        v = self.plot.getPlotProperty(self.selectedPlotType, 'colour', self.selectedPlot)        
        # get new red value
        r = float(self.ui.redSlider.value()) * 0.01
        v1 = (r, v[1], v[2])        
        # change value in plot
        self.plot.changePlotProperty(self.selectedPlotType, 'colour', self.selectedPlot, value = v1)        
        # set dial value
        self.ui.redValueBox.setValue(r)
        
    def setColourRedDial(self):
        ''' change red colour from dial '''
        # get old values
        v = self.plot.getPlotProperty(self.selectedPlotType, 'colour', self.selectedPlot)
        # get new red value
        r = self.ui.redValueBox.value()
        v1 = (r, v[1], v[2])        
        # change value in plot
        self.plot.changePlotProperty(self.selectedPlotType, 'colour', self.selectedPlot, value = v1)        
        # set slider value
        self.ui.redSlider.setValue(int(r*100.))
        

    def setColourGreen(self):
        ''' change colours from sliders '''        
        # get old values
        v = self.plot.getPlotProperty(self.selectedPlotType, 'colour', self.selectedPlot)        
        # get new green value
        g = float(self.ui.greenSlider.value()) * 0.01
        v1 = (v[0], g, v[2])        
        # change value in plot
        self.plot.changePlotProperty(self.selectedPlotType, 'colour', self.selectedPlot, value = v1)        
        # set dial value
        self.ui.greenValueBox.setValue(g)
        
    def setColourGreenDial(self):
        ''' change green colour from dial '''
        # get old values
        v = self.plot.getPlotProperty(self.selectedPlotType, 'colour', self.selectedPlot)
        # get new green value
        g = self.ui.greenValueBox.value()
        v1 = (v[0], g, v[2])        
        # change value in plot
        self.plot.changePlotProperty(self.selectedPlotType, 'colour', self.selectedPlot, value = v1)        
        # set slider value
        self.ui.greenSlider.setValue(int(g*100.))
        
    def setColourBlue(self):
        ''' change colours from sliders '''        
        # get old values
        v = self.plot.getPlotProperty(self.selectedPlotType, 'colour', self.selectedPlot)        
        # get new blue value
        b = float(self.ui.blueSlider.value()) * 0.01
        v1 = (v[0], v[1], b)        
        # change value in plot
        self.plot.changePlotProperty(self.selectedPlotType, 'colour', self.selectedPlot, value = v1)        
        # set dial value
        self.ui.blueValueBox.setValue(b)
        
    def setColourBlueDial(self):               
        ''' change blue colour from dial '''
        # get old values
        v = self.plot.getPlotProperty(self.selectedPlotType, 'colour', self.selectedPlot)
        # get new blue value
        b = self.ui.blueValueBox.value()
        v1 = (v[0], v[1], b) 
        # change value in plot
        self.plot.changePlotProperty(self.selectedPlotType, 'colour', self.selectedPlot, value = v1)        
        # set slider value
        self.ui.blueSlider.setValue(int(b*100.))
                

        
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