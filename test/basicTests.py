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
        self.fnameProp = 'data/3d_grid_properties.txt'
        
        self.brain = mb.brainObj()

    #### loading things 
    def test_loadAdj(self):
        ''' test basic loading '''       
        
        self.brain.import_adj_file(self.fnameAdj)
        
    
    #### Different types of thresholding
    def test_loadAndThresholdVal(self):
        ''' load files and threshold '''
        
        # load brain
        self.brain.import_adj_file(self.fnameAdj)
        
        # threshold by absolute value
        self.brain.apply_threshold(tVal = 0.5)
        
        # threshold by perecentage of edges
        self.brain.apply_threshold(edgePC = 50)
                
        # threshold by number of edges
        self.brain.apply_threshold(totalEdges = 2)
        
    def test_loadAndPlot(self):
        ''' load, threhsold and plot '''
        
        br, plt = recipes.loadAndPlot(self.fnameAdj, self.fnameCo, 0.5, opacity = 0.2)
        
        plt.show()
        
    def test_addProperties(self):
        ''' add properties and highlights from a file and plot '''
        
        br = recipes.loadAndThreshold(self.fnameAdj, self.fnameCo, 0.5)
        
        br.import_properties(self.fnameProp)
        
        # highlight nodes with x value greater than 5
        br.highlight_from_conds('x', 'gt', 0.5, label ='x1', mode ='node', colour = (0.5, 0.5, 0.), opacity = 0.5)
        
        # highlight edges labelled green 
        br.highlight_from_conds('colour', 'eq', 'green', label ='green', mode ='edge', colour = (0., 1., 0.), opacity = 0.5)
        
#        br.highlight_from_conds(prop, rel, val, label = None, mode = 'edge', colour = (1.,0.,0.), opacity = 1.0)
#        
#        br.highlight_from_conds(prop, rel, val, label = None, mode = 'edge', colour = (1.,0.,0.), opacity = 1.0)
#        
#        br.highlight_from_conds(prop, rel, val, label = None, mode = 'edge', colour = (1.,0.,0.), opacity = 1.0)
        
        
        
        

        

if __name__ == '__main__':
    unittest.main()