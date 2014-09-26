# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19 22:48:02 2012

@author: hypocrates

problem in makeSubBrain 

"""

import mayBrainTools as nb
from mayavi import mlab
from numpy import random
import demoVideoAnimate as animate

# initialise brain
a = nb.brainObj()
plot = nb.plotObj()

# read in files for adjacency, spatial info and skull
a.readAdjFile("test data/wave_cor_mat_level_4d_400.txt", 0.5)
a.readSpatialInfo("test data/parcel_400_xyz.txt")
a.importSkull('test data/avg152T1_LR_nifti.nii')

# plot brain and skull
plot.plotBrain(a, opacity = 0.2, label='brain')
contourVals = range(0, 255, 75)
plot.plotSkull(a, contourVals = contourVals, opacity = 1.0, label='skull')

# assign colours at random
l = len(a.G.nodes())
cols = ['red', 'green', 'blue']

for ind in range(l):
    n = int(len(cols)*random.rand())
    col = cols[n]
    a.G.node[ind]['colour'] = col
#    a.G.node[ind]['X'] = a.G.node[ind]['xyz'][0]

# get sub brains of each colour
redBrain = a.makeSubBrain('colour', 'red')
blueBrain = a.makeSubBrain('colour', 'blue')
greenBrain = a.makeSubBrain('colour', 'green')

# plot brains of each colour
plot.plotBrain(redBrain, col=(1, 0, 0), opacity = 1., label = 'red')
plot.plotBrain(greenBrain, col=(0, 1, 0), opacity = 1., label = 'green')
plot.plotBrain(blueBrain, col=(0, 0, 1), opacity = 1., label = 'blue')

# animate the image
animate.anim(plot, 1, folder = 'demoVideo/')
mlab.show()