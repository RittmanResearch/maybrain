# -*- coding: utf-8 -*-
"""
Created on Sat Oct 20 14:54:34 2012

Animates the image and saves a snapshot of each view to file. To turn into video
with mencoder, use mencoder "mf://demoVideo/im %04d.png" -mf fps=10 -o demoVideo.avi -ovc lavc -lavcopts vcodec=msmpeg4v2:vbitrate=500

To use this script:
    1. change the anim function to do what you require the animation to
    2. import this script to your main script (you may want to save it under a different filename)
    3. call the anim(...) function at the end of the main script
    4. call mlab.show() after calling anim(...)


"""

import mayBrainTools as nb
from mayavi import mlab
from numpy import arange

## ================================================================================

# Functions to animate a plot
# The first function here can be edited to carry out the required animation

# function for animation of plot
@mlab.animate
def anim(obj, saveBool, folder = ''):
    ''' Animate the object as specified. If saveBool, each frame is saved to a 
        folder specified as a global variable in this script '''
        
    def makefname(folder, count, length):
        ''' make a filename, all files are called 'im' and saved in folder specified at the top of this script '''
        
        n = str(count)
        zeros = (length - len(n)) * '0'
        n = zeros + n
        print n
        fname = folder + 'im ' + n + '.png'
                
        return fname
        
        
    def postChange(fig, saveFig, count):
        ''' standard stuff that happens after the image is changed. Need to pass the current figure
            and a boolean (saveFig) to determine whether to save or not, plus count to determine figure
            name '''
        fig.scene.render()
        if saveFig:
            fname = makefname(folder, count, 4)
            mlab.savefig(fname)
            count = count + 1
        
        return count        
    
    # initialise things
    filecount = 0    # count the number of images saved
    f = mlab.gcf()   # get the current mayavi figure
    r = arange(1., 0., -0.06) # a range of values for opacity

    # set initial views
    
    # turn off red brain
    obj.change_plot_property('brainEdge', 'visibility', 'red')
    obj.change_plot_property('brainNode', 'visibility', 'red')

    # turn off green brain
    obj.change_plot_property('brainEdge', 'visibility', 'green')
    obj.change_plot_property('brainNode', 'visibility', 'green')

    # turn off blue brain
    obj.change_plot_property('brainEdge', 'visibility', 'blue')
    obj.change_plot_property('brainNode', 'visibility', 'blue')



    # skull fading
    for s in r:
        if s<0.:
            s = 0.
        obj.change_plot_property('skull', 'opacity', 'skull', value = s)
        filecount = postChange(f, saveBool, filecount)
        yield

    obj.change_plot_property('skull', 'opacity', 'skull', value = 0.01)

    # brain nodes and edges fading
    for s in r:
        if s<0.:
            s = 0.
        obj.change_plot_property('brainNode', 'opacity', 'brain', value = s)
        obj.change_plot_property('brainEdge', 'opacity', 'brain', value = s)
        filecount = postChange(f, saveBool, filecount)
        yield        

    # fix properties
    obj.change_plot_property('brainNode', 'opacity', 'brain', value = 1.)
    obj.change_plot_property('brainEdge', 'opacity', 'brain', value = 0.1)
    obj.change_plot_property('skull', 'visibility', 'skull')
    
    # highlight each colour
    for ind in ['red', 'blue', 'green']:
        # switch on
        obj.change_plot_property('brainEdge', 'visibility', ind)
        obj.change_plot_property('brainNode', 'visibility', ind)
        obj.change_plot_property('brainEdge', 'opacity', ind, value = 0.1)
        obj.change_plot_property('brainNode', 'opacity', ind, value = 1.)
        for x in range(10):
            filecount = postChange(f, saveBool, filecount)
            yield
        # switch off
        obj.change_plot_property('brainEdge', 'visibility', ind)
        obj.change_plot_property('brainNode', 'visibility', ind)
        postChange(f, saveBool, filecount)
        yield
        
    # fix properties
    obj.change_plot_property('brainNode', 'opacity', 'brain', value = 1.)
    obj.change_plot_property('brainEdge', 'opacity', 'brain', value = 0.1)
    obj.change_plot_property('skull', 'visibility', 'skull')
    obj.change_plot_property('skull', 'opacity', 'skull', value = 0.1)
    
    
    # spin it round
    for x in range(60):
        f.scene.camera.azimuth(1)
        filecount = postChange(f, saveBool, filecount)
        yield        
#    
#    while 1:
#        f.scene.camera.azimuth(1)
#        filecount = postChange(f, 0, filecount)
#        yield
        



