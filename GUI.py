# -*- coding: utf-8 -*-
"""

GUI for Maybrain

"""

# First, and before importing any Enthought packages, set the ETS_TOOLKIT
# environment variable to qt4, to tell Traits that we will use Qt.


#import os
#os.environ['ETS_TOOLKIT'] = 'qt4'

# To be able to use PySide or PyQt4 and not run in conflicts with traits,
# we need to import QtGui and QtCore from pyface.qt
from pyface.qt import QtGui, QtCore
#import sip
#sip.setapi('QString', 2)

from numpy import log10
import maybrain as mb

import sys
import mayBrainUI as ui
from os import path


class mayBrainGUI(QtGui.QMainWindow):
    
    #### Initialisation and seutp
    
    def __init__(self, parent=None):
#        app=QtGui.QApplication([])  # is this required to get it working properly?
#        app.exec_()
        
        QtGui.QMainWindow.__init__(self, parent)
        # set up UI
        self.ui = ui.Ui_Maybrain()
        self.ui.setupUi(self)
        # disable some buttons
        self.ui.adjPlot.setEnabled(0)
        self.ui.skullPlot.setEnabled(0)

        # set up function
        self.brains = {} # dictionary to hold brain objects
        self.plot = mb.mbplot.plotObj() # start plot object
        self.highlights = {} # dictionary of highlights, each entry contains [brainName, highlightName]
        
        # link up buttons and functions
        self.connectSlots()
        
        # plot selected in ui.plotTree
        self.selectedPlot = None
        
        # local variables
        self.lastFolder = '' # last folder viewed by user
        self.brainName = ['brain', 0] # used to auto-create labels for brains if no user input given (currently only supports 1 brain)
        self.currentBrainName = None # currently used brain
        self.plName = ['plot', 0]
        self.bgName = ['bg', 0] # default background name prefix and index
        self.isoName = ['iso', 0] # default isosurface name and index
        self.brInds = 0 # counter for number of brains
        
        
    def connectSlots(self):
        ''' connect buttons and functions '''

        # quit app        
        QtCore.QObject.connect(self.ui.quitButton, QtCore.SIGNAL('clicked()'), app.quit)
        
        # tab 1: loading files and global plotting
        QtCore.QObject.connect(self.ui.spatialFnameButton, QtCore.SIGNAL('clicked()'), self.getSpatialFilename)
        QtCore.QObject.connect(self.ui.adjFnameButton, QtCore.SIGNAL('clicked()'), self.getAdjFilename)
        QtCore.QObject.connect(self.ui.skullFnameButton, QtCore.SIGNAL('clicked()'), self.getSkullFilename)
        QtCore.QObject.connect(self.ui.propsFnameButton, QtCore.SIGNAL('clicked()'), self.getPropsFilename)                
        QtCore.QObject.connect(self.ui.brainLoad, QtCore.SIGNAL('clicked()'), self.loadBrain)
        QtCore.QObject.connect(self.ui.skullLoad, QtCore.SIGNAL('clicked()'), self.loadSkull)
        QtCore.QObject.connect(self.ui.isoLoad, QtCore.SIGNAL('clicked()'), self.loadIso)   
        QtCore.QObject.connect(self.ui.skullPlot, QtCore.SIGNAL('clicked()'), self.plotSkull)
        QtCore.QObject.connect(self.ui.isoPlot, QtCore.SIGNAL('clicked()'), self.plotIso)
        QtCore.QObject.connect(self.ui.clearFigButton, QtCore.SIGNAL('clicked()'), self.clearPlot)
        
        # tab 2: highlights
        self.ui.hlApplyButton.clicked.connect(self.makeHighlight)
        self.ui.propsLoad.clicked.connect(self.addProperties)
        
        # tab 3: plotting        
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
        self.ui.brainSelect.currentIndexChanged.connect(self.selectActiveBrain)
        
        # tab 4: metrics
        self.ui.degreeButton.clicked.connect(self.getDegree)

        
#        QtCore.QObject.connect(self.ui.opacitySlider, QtCore.SIGNAL('mouseReleaseEvent()'), self.setOpacity)

    ## =============================================
    
    #### General functions

    def selectActiveBrain(self):
        ''' changes the active brain when selected in the box at the top of the GUI (brainSelect)'''

        print(str(self.ui.brainSelect.currentText()))
        self.activeBrain = self.brains[str(self.ui.brainSelect.currentText())]
        
        
        
    ## =============================================

    #### Functions to connect slots and Maybrain functions
        
    # Basic info, files and loading
    def getSpatialFilename(self):
        ''' open a dialog to get the spatial filename'''
        f = QtGui.QFileDialog.getOpenFileName(self, 'Choose spatial file', directory = self.lastFolder)
        self.ui.spatialFilename.setText(f)
        self.lastFolder = path.dirname(str(f))
        
    def getAdjFilename(self):
        ''' open a dialog to get the adjacency filename '''
        f = QtGui.QFileDialog.getOpenFileName(self, 'Choose adjacency file', directory = self.lastFolder)
        self.ui.adjFilename.setText(f)
        self.lastFolder = path.dirname(str(f))
        
    def getSkullFilename(self):
        ''' get a skull filename '''
        f = QtGui.QFileDialog.getOpenFileName(self, 'Choose skull file', directory = self.lastFolder)
        self.ui.skullFilename.setText(f)
        self.lastFolder = path.dirname(str(f))

    def getIsoFilename(self):
        ''' get an isosurface filename '''
        f = QtGui.QFileDialog.getOpenFileName(self, 'Choose isosurface file', directory = self.lastFolder)
        self.ui.isoFilename.setText(f)
        self.lastFolder = path.dirname(str(f))
        
    def getPropsFilename(self):
        f = QtGui.QFileDialog.getOpenFileName(self, 'Choose properties file', directory = self.lastFolder)
        self.ui.propsFilename.setText(f)
        self.lastFolder = path.dirname(str(f))
        
    def findBrainName(self):
        ''' construct a new brain name if necessary '''

        # Look for name entered by user
        try:
            brainName = str(self.ui.brainName.text())
        except:
            print("invalid name")
            brainName = ''
            
        # get rid of unwanted characters
        brainName = mb.extraFns.stripString(brainName)
            
        if brainName == '':
        
            brainName = self.brainName[0] + str(self.brainName[1])
            self.brainName[1] = self.brainName[1] + 1
            
        # check if it's been used previously
        brainUsed = False
        if brainName in self.brains:
            brainUsed = True
                        
        # returns the name and whether the name has already been used to label a brain
        return brainName, brainUsed


    def findPlotName(self):
        ''' construct a new brain name if necessary '''

        # Look for name entered by user
        try:
            plName = str(self.ui.plotName.text())
        except:
            print("invalid name")
            plName = ''
            
        # get rid of unwanted characters
        plName = mb.extraFns.stripString(plName)
        print(plName)
            
        if plName == '':
        
            plName = self.plName[0] + str(self.plName[1])
            self.plName[1] = self.plName[1] + 1
            
#        # check if it's been used previously
#        plotUsed = 0
#        if (plName in self.plot.brainNodePlots)|(plName in self.plot.brainEdgePlots):
#            plotUsed = 1
            
        self.currentPlotName = plName
            
        # returns the name and whether the name has already been used to label a brain
        return plName #, plotUsed
        
        
    def findBGname(self):
        ''' get the name for a background '''
        
        try:
            bgName = str(self.ui.bgName.text())
        except:
            print("invalid background name")
            bgName = ''
            
        # remove unwanted characters
        bgName = mb.extraFns.stripString(bgName)
        
        if bgName == '':
            bgName = self.bgName[0] + str(self.bgName[1])
            self.bgName[1] = self.bgName[1] + 1
            
        return bgName
        
        
    def findIsoName(self):
        ''' get the name for an isosurface '''
        
        try:
            isoName = str(self.ui.isoName.text())
        except:
            print("invalid isosurface name")
            isoName = ''
            
        # remove unwanted characters
        isoName = mb.extraFns.stripString(isoName)
        
        if isoName == '':
            isoName = self.isoName[0] + str(self.isoName[1])
            self.isoName[1] = self.isoName[1] + 1
            
        return isoName        
        

        
    def loadBrain(self):
        ''' load a brain using given filenames '''
        self.ui.adjPlot.setEnabled(False)        
        
        # get adjacency filename
        adj = str(self.ui.adjFilename.text())
        # get threshold
        thType, thVal = self.getThresholdType()
        # get spatial info file
        coords = str(self.ui.spatialFilename.text())
        
        # make name for brain
        brName, brainUsedBool = self.findBrainName()
        self.currentBrainName = brName
        
        # create and add to list of brains
        if brainUsedBool:
            # case where brain name exists already
            br = self.brains[brName]
        else:
            # make a new brain
            br = mb.brainObj()
                        
        # add properties
        br.importAdjFile(adj)
        br.importSpatialInfo(coords)
        br.applyThreshold(thresholdType = thType, value = thVal)
        
        self.brains[brName] = brma
        
        # add to brains selected for highlighting
        if not(brainUsedBool):
            self.ui.brainSelect.addItem(brName)

        # make sure the active brain is the desired one - this is a bit of a hack, but seems to work
        brInd = len(self.brains.keys()) - self.brains.keys().index(brName)-1
        self.ui.brainSelect.setCurrentIndex(brInd)
        self.selectActiveBrain()
                   
        # enable plot button
#        try:
#            QtCore.QObject.disconnect(self.ui.adjPlot, QtCore.SIGNAL('clicked()'), self.rePlotBrain)
#        except:
#            pass
        QtCore.QObject.connect(self.ui.adjPlot, QtCore.SIGNAL('clicked()'), self.plotBrain)
        self.ui.adjPlot.setEnabled(True)
     
        
    def loadSkull(self):
        ''' load a skull file '''
        # get filename
        f = str(self.ui.skullFilename.text())
        # get brain name
        brainName, brUsedBool = self.findBrainName()
    
        # create brain object if it doesn't exist
        if not(brUsedBool):
            br = mb.brainObj()
            self.brains[brainName] = br
        else:
            br = self.brains[brainName]
        # read in file
        br.importSkull(f)
        
        # enable plot button
        self.ui.skullPlot.setEnabled(True)
        
    def loadIso(self):
        ''' load a skull file '''
        # get filename
        f = str(self.ui.isoFilename.text())
        # get brain name
        brainName, brUsedBool = self.findBrainName()
    
        # create brain object if it doesn't exist
        if not(brUsedBool):
            br = mb.brainObj()
            self.brains[brainName] = br
        else:
            br = self.brains[brainName]
        # read in file
        br.importISO(f)
        
        # enable plot button
        self.ui.isoPlot.setEnabled(True)        
        
        
    def getThresholdType(self):
        ''' read the value in the threshold dropdown box to determine thresholding method '''
        
        thType = str(self.ui.tholdDropdown.currentText())
        value = float(self.ui.thresholdValue.text())
        
        print(thType, value)
        
#        edgePC = None 
#        totalEdges = None
#        tVal = None
        
        if thType=='Value':
            thType = 'tVal'
            
        elif thType=='%':
            thType = 'edgePC'
            
        elif thType=='# edges':
            thType = 'totalEdges'
        else:
            print('type not recognised')
            thType = 'tVal'
            
        return thType, value

    ## =============================================
    
    #### plotting functions
    
    def plotBrain(self):
        ''' plot the entire brain object '''
        
        brName = self.getActiveBrainName()
        
#        # sort label for brain object
#        brName, brUsed = self.findBrainName()
#        if not(brUsed):
#            brName = self.currentBrainName
            
        # get threshold
        thType, thVal = self.getThresholdType()
        # apply threshold
        self.brains[brName].applyThreshold(thType, thVal)
                
        # get plot name
        plname = self.findPlotName()
                
        # plot the brain
        self.plot.plotBrain(self.brains[brName], opacity = 0.2, edgeOpacity = None, label=plname, plotHighlights = False)

        # add to tree view
        self.readAvailablePlots() 
            
            
    def plotSkull(self):
        ''' plot the skull '''
        
        # get current brain name
        brName = self.getActiveBrainName()    
        
        # get label
        bgLabel = self.findBGname()
        
        # plot
        self.plot.plotSkull(self.brains[brName], label = bgLabel)
            
        # add to treeview
        self.readAvailablePlots()

         
    def plotIso(self):
        ''' plot an isosurface '''
        
        # get current brain name
        brName = self.getActiveBrainName()
        
        # get label
        isolabel = self.findIsoName
        
        # plot
        self.plot.plotIsosurface(self.brains[brName], label = isolabel)
        
    
        
    def readAvailablePlots(self):
        ''' read in the available plots from the plotting object '''
        
        # clear current values
        self.ui.plotTree.clear()
        
        # add new values
        lists = [self.plot.brainEdgePlots, self.plot.brainNodePlots, self.plot.skullPlots, self.plot.isosurfacePlots]
        labels = ['edges', 'nodes', 'skull', 'isosurf']
        
        for ls in range(len(lists)):
            names = lists[ls].keys()
            names.sort()
            lab = labels[ls]
            for p in names:
                QtGui.QTreeWidgetItem(self.ui.plotTree, [lab, p])
        
#        for p in self.plot.brainEdgePlots:
#            QtGui.QTreeWidgetItem(self.ui.plotTree, ['brainEdge', p])
#        for p in self.plot.brainNodePlots:            
#            QtGui.QTreeWidgetItem(self.ui.plotTree, ['brainNode', p]) 
#        for s in self.plot.skullPlots:
#            QtGui.QTreeWidgetItem(self.ui.plotTree, ['skull', s])
#        for s in self.plot.isosurfacePlots:
#            QtGui.QTreeWidgetItem(self.ui.plotTree, ['isosurf', s])


    def clearPlot(self):
        ''' clear all current plots'''
        
        # delete plots from figure
        for p in self.plot.brainNodePlots:
            try:
                self.plot.brainNodePlots[p].remove()
            except ValueError:
                continue
        for p in self.plot.brainEdgePlots:
            try:
                self.plot.brainEdgePlots[p].remove()
            except ValueError:
                continue
        for p in self.plot.skullPlots:
            try:
                self.plot.skullPlots.remove()
            except ValueError:
                continue                
        for p in self.plot.isosurfacePlots:
            try:
                self.plot.isosurfacePlots.remove()
            except ValueError:
                continue                
        
        
        # remove plots from list in plot object
        self.plot.brainEdgePlots = {}
        self.plot.brainNodePlots = {}
        self.plot.skullPlots = {}
        self.plot.isosurfacePlots = {}

        # clear the mayavi figure
        mb.mbplot.mlab.clf()
        
        # redisplay list
        self.readAvailablePlots()


    ## =============================================
    
    #### Highlight functions
    
    def makeHighlight(self):
        ''' make a highlight using settings from ui '''
        
        # get current brain object
        br = self.selectActiveBrain()
        br = self.activeBrain
        
        # get settings from ui
        propName = str(self.ui.hlProp.currentText())
        relation = self.getRelation()
        try:
            val1 = float(self.ui.hlValue1.text())
        except:
            val1 = str(self.ui.hlValue1.text())
            val1 = mb.extraFns.stripString(val1)
        try:
            val2 = float(self.ui.hlValue2.text())
        except:
            val2 = str(self.ui.hlValue2.text())
            val2 = mb.extraFns.stripString(val2)
        # *** should change the name to hlName ***
        label = str(self.ui.subPlotName.text())
        # NEEDS MODIFICATION
        if label=='':
            label = None
        mode = str(self.ui.hlNodesOrEdgesBox.currentText())
        red = float(self.ui.hlRedSpin.value())        
        green = float(self.ui.hlGreenSpin.value())
        blue = float(self.ui.hlBlueSpin.value())
        # *** maybe change to a box??
        opacity = float(self.ui.hlOpacityValue.text())
        
        # get value correct
        if relation in ['in()', 'in[)', 'in(]', 'in[]']:
            val = [val1, val2]
        else:
            val = val1
        
        # create the highlight object
        br.highlightFromConds(propName, relation, val, label=label, mode=mode, colour = (red, green, blue), opacity=opacity)
        
        # plot
        # bodge for cases when labels change:
#        if label in br.highlights:
        print('highlights: ', br.highlights)
        self.plot.plotBrainHighlights(br, highlights=[label])
        try:
            print(br.highlights['green'].edgeIndices)
        except:
            print('it aint green')
#        else:
#            self.plot.plotBrainHighlights(br, highlights = [label+'_edges'])
#            self.plot.plotBrainHighlights(br, highlights = [label+'_nodes'])
        
        # add to list of plots
        self.readAvailablePlots()
        
        
    def getRelation(self):
        ''' retrieve information from the relation-box and translate into maybraintools language '''
        # get current value
        val = str(self.ui.hlRelationBox.currentText())
        
        # translate        
        if val=='=':
            outval = 'eq'
        elif val=='<':
            outval = 'lt'
        elif val=='>':
            outval = 'gt'
        elif val=='<=':
            outval = 'leq'
        elif val=='>=':
            outval = 'geq'
        elif val=='contains text':
            outval = 'contains'
        else:
            # covers cases in + 2 brackets
            outval = val
        
        return outval
    
#    def plotHighlight(self):
#        ''' plot the selected highlight '''
#        # Get highlight name
##        brainName, hlName = 
#        
#        # Do the plot
#        self.plot.plotHighlights(self.br[brainName], highlights=[hlName])
#        
#        QtGui.QTreeWidgetItem(self.ui.plotTree, ['brainNode', label])
#        QtGui.QTreeWidgetItem(self.ui.plotTree, ['brainEdge', label])

        
    ## =============================================
    
    ## setting and altering plot properties e.g. visibility, opacity
    
    def addProperties(self):
        ''' add properties to a brain from file '''
        
        # get properties filename from GUI
        fname = str(self.ui.propsFilename.text())        
        
        # find the active brain
        brain, usedBool = self.findBrainName()
        if not(usedBool):
            brain = self.currentBrainName
        
        # call function from mayBrainTools to add properties to brain
        prop = self.brains[brain].importProperties(fname)
        
        # add to plot properties box
        self.ui.hlProp.addItem(str(prop))
        
        for e in self.brains[brain].G.edges(data = True):
            print(e)
             
    
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
        
        
    def getActiveBrainName(self):
        
        # get current brain object
        brName = str(self.ui.brainSelect.currentText())
        
#        brName, nameUsedBool = self.findBrainName()
        if not(brName in self.brains):
            brName = self.currentBrainName
        else:
            self.currentBrainName = brName
        
        return brName
                
    ## =============================================
    
    #### Analysis Functions
    
    def getDegree(self):
        ''' Find the degree of the active brain and save to file '''

        # get the active brain and find nodes and degrees        
        self.activeBrain = self.brains[str(self.ui.brainSelect.currentText())]
        nodes = self.activeBrain.G.nodes()
        degree = self.activeBrain.G.degree().values()
        
        # get the output file ready
        fname = self.ui.degreeFilename.text()
        f = open(fname, 'w')
        f.write('node\tdegree')
        print(degree)
        
        for d in range(len(degree)):
            print(nodes[d], degree[d])
            f.write('\n'+str(nodes[d]) + '\t' + str(degree[d]))


        f.close()

        
if __name__ == "__main__":

    # using instance allows Mayavi to run alongside without any problems.
#    app = QtGui.QApplication.instance()
#    ex = mayBrainGUI()
#    ex.show()
#    sys.exit(app.exec_())

    

    # create and show app    
    app = QtGui.QApplication(sys.argv)
    myapp = mayBrainGUI()
    myapp.show()
    sys.exit(app.exec_())        