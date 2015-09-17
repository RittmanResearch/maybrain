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
import maybrain.recipes as mb


class TestSequenceFunctions(unittest.TestCase):

    def setUp(self):
        ''' load inital parameters and/or data for testing'''
        self.fnameAdj = 'data/3d_grid_adj.txt'
        self.fnameCo = 'data/3d_grid_coords.txt'

    def test_load(self):
        ''' test basic loading '''
        
        mb.loadFiles(self.fnameAdj, self.fnameCo)
        
    def test_loadAndThreshold(self):
        ''' load files and threshold '''
        
        mb.loadAndThreshold(self.fnameAdj, self.fnameCo, 0.5)
        
    def test_loadAndPlot(self):
        ''' load, threhsold and plot '''
        
        br, plt = mb.loadAndPlot(self.fnameAdj, self.fnameCo, 0.5, opacity = 0.2)
        
        plt.show()

        

if __name__ == '__main__':
    unittest.main()