# -*- coding: utf-8 -*-
"""
Created on Tue Jan 29 14:43:01 2013

@author: martyn

Initialisation file for maybrain



Module description
==========================
Maybrain is a module for visualizing brain connectivity data.

Some practical details
--------------------------
Modules are imported in this file to avoid having to call the 
submodule each time.

To avoid having to import mayavi you can use, e.g.

  ``from mayavin import brainObj``



"""

#from mayBrainTools import *
#from mbplot import *
#import writeFns
#import mayBrainExtraFns
#from recipes import *
#
#
#write = writeFns.writeFns()

# imports to avoid having to import submodules each time
#from recipes import *
from . import brainObj
#from plot import plotObj
from .extraFns import *
#from . import mbplot

#import plot
#import plotObjs
