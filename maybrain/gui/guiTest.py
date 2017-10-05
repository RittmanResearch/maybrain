# -*- coding: utf-8 -*-
"""
Created on Tue Sep 23 07:19:08 2014

@author: galileo
"""

# add the local path to python executable path
#try:
#    import sys
#    sys.path.append("/home/galileo/Dropbox/Share/maybrain/august 2014 dev/dev2-2")
#except:
#    pass

import unittest
#import maybrain as mb
from . import mayBrainGUI as gui
import sys
from pyface.qt import QtGui, QtCore

class TestSequenceFunctions(unittest.TestCase):

    def setUp(self):
        ''' load inital parameters and/or data for testing'''

        # parameters
        self.fnameAdj = 'data/3d_grid_adj.txt'
        self.fnameCo = 'data/3d_grid_coords.txt'


        # start gui        
        self.app = gui.runGui()

    def test_loadData(self):
        ''' test inputting of basic data '''
        
        # enter filenames
        self.app.close()
        
        
        # load files
        
        
        # plot 
        
        # quit application
        
        
        
    def test_highlights(self):
        ''' test the ability to highlight some part of a brain '''
        
        # load files (as in test_loadData)
        a=1
        
        # set higlighting properties
        
        # plot highlight
        
        
        # change some property






if __name__ == '__main__':
    unittest.main()