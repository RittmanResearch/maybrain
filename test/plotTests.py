# -*- coding: utf-8 -*-
"""
Created on Wed Oct  1 07:13:02 2014

@author: galileo

Test module for plotting functions of Maybrain
Basic plotting is included in recipeTests, this module tests more advanced
options:

 - 

"""

import sys
try:
    sys.path.append("/home/galileo/Documents/programming/maybrain/maybrain")
except:
    pass

import unittest
import maybrain as mb
import networkx as nx

class TestSequenceFunctions(unittest.TestCase):

    def setUp(self):
        ''' load inital parameters and/or data for testing'''
        self.fnameAdj = 'data/3d_grid_adj.txt'
        self.fnameCo = 'data/3d_grid_coords.txt'
        
        # initialise brain and plot objects
        self.a = mb.brainObj()
        self.b = mb.plotObj()
        
        # some data
        self.edges = [(0,2),(0,4),(0,6)]
        
        self.locations = {0:(1,1,1), 2:(2,2,2), 4:(3,3,3), 6:(4,4,4)}


    def test_plotSelectedCoords(self):

        # make an empty graph
#        a = 1
        self.a.G = nx.Graph()
        print('run')
        
        # populate graph with edges and nodes
        self.a.G.add_edges_from(self.edges)
        for n in self.locations:
            self.a.G.node[n]['xyz'] = self.locations[n]
                
        
        # plot everything
        self.b.plotBrainCoords(self.a)
        
        # this line will only plot two points and give a printed error 
        # for the others (1 and 3)
        nodes = [0,1,2,3]
        print(nodes)
        self.b.plotBrainCoords(self.a, nodes=nodes)
        
        # plot selected nodes
        # this one will work fine
        nodes = [0,2,6]
        print(nodes)
        self.b.plotBrainCoords(self.a, nodes=nodes)
        
        
    def test_plotOnlyEdges(self):
        
        # clear current plot and brain
        self.b.clear()
        
        # load more data (in a fresh brain)
        self.a = mb.loadAndThreshold(self.fnameAdj, self.fnameCo, 0.5)
        
        # plot the edges
        self.b.plotBrainEdges(self.a)
        

    def test_plotOnlyCoords(self):
        
        # clear current plot and brain
        self.b.clear()
        
        # load more data (in a fresh brain)
        self.a = mb.loadAndThreshold(self.fnameAdj, self.fnameCo, 0.5)
        
        # plot the edges
        self.b.plotBrainCoords(self.a)

        
if __name__ == '__main__':
    unittest.main()        
