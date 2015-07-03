# -*- coding: utf-8 -*-
"""

Some unit tests for Maybrain recipes - loading of data

"""

# add the local path to python executable path
import sys
try:
    sys.path.append("/home/galileo/Dropbox/smart laptop/maybrain/maybrain")
except:
    pass

import unittest
import maybrain as mb
from maybrain import recipes


class TestSequenceFunctions(unittest.TestCase):

    def setUp(self):
        ''' load inital parameters and/or data for testing'''
        self.fnameAdj = 'data/3d_grid_adj.txt'
        self.fnameCo = 'data/3d_grid_coords.txt'
        
        self.brain = mb.brainObj()

    #### loading things 
    def test_loadAdj(self):
        ''' test basic loading '''       
        
        self.brain.importAdjFile(self.fnameAdj)
        
    
    #### Different types of thresholding
    def test_loadAndThresholdVal(self):
        ''' load files and threshold '''
        
        # load brain
        self.brain.importAdjFile(self.fnameAdj)
        
        # threshold by absolute value
        self.brain.applyThreshold(tVal = 0.5)
        
        # threshold by perecentage of edges
        self.brain.applyThreshold(edgePC = 50)
                
        # threshold by number of edges
        self.brain.applyThreshold(totalEdges = 2)
        
    def test_loadAndPlot(self):
        ''' load, threhsold and plot '''
        
        br, plt = recipes.loadAndPlot(self.fnameAdj, self.fnameCo, 0.5, opacity = 0.2)
        
        plt.show()

        

if __name__ == '__main__':
    unittest.main()