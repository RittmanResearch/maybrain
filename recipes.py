# -*- coding: utf-8 -*-
"""
Created on Tue Jan  7 22:38:39 2014

@author: hypocrates

A series of recipes for maybrain, intended to act as a quick start

"""

from mayBrainTools import brainObj
from mbplot import plotObj

def loadFiles(adjname, coordname):
    ''' create a brain object and load adjacency matrix and coordinates. Returns a brain object. '''

    brain = brainObj()
    brain.importAdjFile(adjname)
    brain.importSpatialInfo(coordname)
    
    return brain
    
def loadAndThreshold(adjname, coordname, threshold):
    ''' create a brain object using the given files, apply a threshold value
        and create the networkx objects. Returns a brain object. '''
        
    brain = brainObj()
    brain.importAdjFile(adjname)
    brain.importSpatialInfo(coordname)
    
    brain.applyThreshold(tVal = threshold)
    
    return brain
    
def loadAndPlot(adjname, coordname, threshold, opacity = 1.0):
    ''' create a brain object from the given values, apply a threshold value
        then plot. Returns brain, plot. '''
        
    # make brain object
    brain = brainObj()
    
    # load date from files
    brain.importAdjFile(adjname)
    brain.importSpatialInfo(coordname)
    
    # apply threshold
    brain.applyThreshold(tVal = threshold)
    
    # plot and show
    plt = plotObj()
    plt.plotBrain(brain, opacity=opacity)
    
    return brain, plt