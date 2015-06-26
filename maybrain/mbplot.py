# -*- coding: utf-8 -*-
"""
Plot functions for Maybrain

The plotObj class creates an object to which many different plots are added. 

Examples of plots that can be added are brain networks with highlights, isosurfaces
and other surfaces from nibabel (see maybraintools.py). Various properties of the
plot can be changed using this class as well.

The plots are all made using mayavi (http://code.google.com/p/maybrain/), and any 
mayavi functions can be accessed by calling the mlab function with the appropriate 
plot. Plots are all added to a dictionary (one dictionary per data type, e.g. one 
for node plots, self.brainNodePlots). In addition, Mayavi has a powerful gui editor. 
We plot with this so you can take advantage of its many features to customize your 
plot once it has been made.

You need to run mlab.show to see the finished plots. 

Released under some kind of GNU/open source license. Basically, do what you want
but give us and the guys at Mayavi at least some credit and make sure your code 
is available for others to use.

More details at http://code.google.com/p/maybrain/.

"""


from mayavi import mlab
from numpy import array,repeat,max,power


#!! plotObj taken from dev version. Legacy code removed    
class plotObj():
    ''' classes that plot various aspects of a brain object '''
    
    
    def __init__(self, bg=None):
        
        # initialise mayavi figure
        self.startMayavi(bg)
        
    def startMayavi(self,bg):
        ''' initialise the Mayavi figure for plotting '''        
        
        # start the engine
        from mayavi.api import Engine
        self.engine = Engine()
        self.engine.start()
        
        if not bg:
            bg=(1., 1., 1.)
        
        # create a Mayavi figure
        self.mfig = mlab.figure(bgcolor = bg,
                                fgcolor = (0, 0, 0),
                                engine = self.engine,
                                size=(1500, 1500))
        
        # holders for plot objects
        self.brainNodePlots = {}
        self.brainEdgePlots = {}
        self.skullPlots = {}
        self.isosurfacePlots = {}

        # autolabel for plots
        self.labelNo = 0         
        
        
    def coordsToList(self, brain, nodeList=None):
        ''' get coordinates from lists in the brain object, possibly using only
        the indices given by nodeList and edgeList '''
        
        # select some rows if necessary
        if not nodeList:
            nodeList = brain.G.nodes()        
        
        # get coordinates
        # make into array for easy output                    
        coords = array([brain.G.node[x]['xyz'] for x in nodeList])
            
#        # get nodes from networkx object
#        else:
#            for x in brain.G.nodes():
#                # put into list
#                coords.append(brain.G.node[x]['xyz'])
                    
#        coords = array(coords)                
        
        # return x, y and z coordinates
        return coords[:, 0], coords[:, 1], coords[:, 2]
        
        
    def edgesToList(self, brain):
        ''' Turn the edges of a brain into coordinates '''
                
        # intialise output
        x1 = []
        x2 = []
        
        y1 = []
        y2 = []
        
        z1 = []
        z2 = []
        
        s = []
        
        for e in brain.G.edges(data = True):
            # get coord and vector from each edge
            p1 = brain.G.node[e[0]]['xyz']
            p2 = brain.G.node[e[1]]['xyz']
            
            x1.append(p1[0])
            x2.append(p2[0]-p1[0]) 
            
            y1.append(p1[1])
            y2.append(p2[1]-p1[1])
            
            z1.append(p1[2])
            z2.append(p2[2]-p1[2])
            
            # set scalar value as edge weight
            s.append(e[2]['weight'])
        
        return x1, y1, z1, x2, y2, z2, s

        

    def plotBrain(self, brain, opacity = 1.0, edgeOpacity = None, label = 'plot'):
        ''' plot all the coords, edges and highlights in a brain '''     
               
        # sort out the edge opacity
        if not(edgeOpacity):
            edgeOpacity = opacity

        # plot coords
        coords = self.coordsToList(brain)
        self.plotCoords(coords, opacity = opacity, label=label)
        
        # plot  edges
        ex1, ey1, ez1, ux, uy, yz, s = self.edgesToList(brain)    
        self.plotEdges(ex1, ey1, ez1, ux, uy, yz, s, col = (0.,0.,0.), opacity = opacity, label=label)  
        
        # plot the highlights
        self.plotBrainHighlights(brain)
        
    def plotBrainCoords(self, brain, nodes=None, opacity = 1.0, label = 'coordplot', sizeList=None, col=(0.,0.,0.),
                        sf=None, sfRange=None):
        ''' plot all coordinates in a brain '''
        
        coords = self.coordsToList(brain, nodeList=nodes)
        self.plotCoords(coords, opacity = opacity, label=label, col=col, sizeList=sizeList, sf=sf, sfRange=sfRange)
        
        
    def plotBrainEdges(self, brain, opacity = 1.0, label = 'edgeplot', col=None, cm ='GnBu', lw=2., scalars=None):
        ''' plot all edges within a brain 
        lw = line width
        #### Not working currently - plots in a different window #### --Really, seems OK to me!
        ''' 
        
        ex1, ey1, ez1, ux, uy, yz, s = self.edgesToList(brain)
        if scalars:
            s=scalars
        self.plotEdges(ex1, ey1, ez1, ux, uy, yz, s, lw=lw, col = col, cm=cm, opacity = opacity, label=label)  


    def plotBrainHighlights(self, brain, highlights = [], labelPre = ''):    
        ''' plot all or some of the highlights in a brain
            labelPre allow for a standard prefix (this is used by the GUI) '''

        if highlights == []:
            highlights = brain.highlights
        
        # plot highlights (subsets of edges or points)
        for h in highlights:
            label = labelPre + h
            try:
                ho = brain.highlights[h]
            except:
                print('highlight not found: ' + h)
                continue
                
            ex1, ey1, ez1, ux, uy, yz, s = ho.getEdgeCoordsToPlot(brain)
            
            # case where edge opacity is not separately defined
            if not(ho.edgeOpacity):
                ho.edgeOpacity = ho.opacity
            
            # plot the edges
            if not(len(ex1)==0):
                self.plotEdges(ex1, ey1, ez1, ux, uy, yz, s, col = ho.colour, opacity = ho.edgeOpacity, label=label)
            
            # plot the nodes
            hp = ho.nodeIndices
            if not(len(hp))==0:
                x, y, z = ho.getCoords(brain)
                self.plotCoords((x,y,z), col = ho.colour, opacity = ho.opacity, label=label)        
        
        
    def plotCoords(self, coords, col = (1.,1.,1.), opacity = 1., label='plot', sizeList=None, sf=1., sfRange=None):
        ''' plot the coordinates of a brain object
            "absoluteScaling" is an option to use and absolute rather than a relative scaling, particularly useful for multiple plots
           
        '''
        if sizeList==None:
            # note that scalar value is currently set to 1.
            ptdata = mlab.pipeline.scalar_scatter(coords[0], coords[1], coords[2],
                                                  figure = self.mfig)
            sf = 1.
            
        else:
            try:
                float(sizeList)
                sizeList = repeat(sizeList, len(coords[0]))
            except:
                pass

            if not sf:
                sf = 5./power(max(sizeList), 1/3)
                print "sf calculated as: "+str(sf)
                
            ptdata = mlab.pipeline.scalar_scatter(coords[0], coords[1], coords[2],
                                                      sizeList, figure = self.mfig)
             
        p = mlab.pipeline.glyph(ptdata, color = col, opacity = opacity, scale_factor=sf,
                                scale_mode="scalar")
        
        if sfRange:
            print "Adjusting glyph range"
            p.glyph.glyph.range = array(sfRange)

        
        self.brainNodePlots[label] = p
        print(label, p)
        
        
    def plotEdges(self, ex1, ey1, ez1, ux, uy, uz, s, lw=2., col = None, opacity = 1., cm = 'GnBu', label='plot'):
        ''' plot some edges
        
            ec has the order [x0, x1, y0, y1, z0, z1]
            s is the scalars - used to determine the node colour
            cm = 'GnBu'  colormap for plotting scalars
            col is a triple of numbers in the interval (0,1), or None
            lw is the line width
            
            
        '''
        
        plotmode = '2ddash' # could be modified later
        
        # add data to mayavi
#        edata = mlab.pipeline.vector_scatter(ex1, ey1, ez1, ux, uy, yz, scalars = s, figure = self.mfig)
        # plot
        v = mlab.quiver3d(ex1, ey1, ez1, ux, uy, uz, scalars = s, line_width=lw, opacity=opacity, mode = plotmode, color = col, scale_factor = 1., scale_mode = 'vector', colormap = cm)
        if not col:
            v.glyph.color_mode = 'color_by_scalar'
            
        self.brainEdgePlots[label] = v

            
#    #!! old version of plotbrain removed here
#
#    def plotBrainNodes(self, brain, nodes = None, col = (0.,0.,0.,), opacity = 1., label=None):
#        ''' plot the nodes using Mayavi TO BE DEPRECATED'''
#        
#        # sort out keywords
#        if not nodes:
#            nodeList = brain.G.nodes()
#        else:
#            nodeList = nodes
#            
#        if not(label):
#            label = self.getAutoLabel()       
#            
#        # turn nodes into lists for plotting
#        xn, yn, zn = self.getNodesList(brain, nodeList=nodeList)
#        
#        # plot nodes
#        s = mlab.points3d(xn, yn, zn, scale_factor = self.nodesf, color = col, opacity = opacity)
#        self.brainNodePlots[label] = s
#
#
    #!! old plotBrainEdges and plotSubset removed here       
            
    def getCoords(self, brain, edge):
        ''' get coordinates from nodes and return a coordinate and a vector '''
        
        c1 = brain.G.node[edge[0]]["xyz"]
        c2 = brain.G.node[edge[1]]["xyz"]
        
        diff = [c2[0]-c1[0], c2[1]-c1[1], c2[2]-c1[2]]    
        
        return c1, diff           
    
    
    def plotSkull(self, brain, label = None, contourVals = [], opacity = 0.1, cmap='Spectral'):
        ''' plot the skull using Mayavi '''
        
        if not(label):
            label = self.getAutoLabel()        
        
        if contourVals == []:            
            s = mlab.contour3d(brain.background, opacity = opacity, colormap=cmap)
        else:
            s = mlab.contour3d(brain.background, opacity = opacity, contours = contourVals, colormap=cmap)
            
        # get the object for editing
        self.skullPlots[label] = s
        
        
    def plotIsosurface(self, brain, label = None, contourVals = [], opacity = 0.1, cmap='autumn'):
        ''' plot an isosurface using Mayavi, almost the same as skull plotting '''
        
        if not(label):
            label = self.getAutoLabel()        
        
        if contourVals == []:            
            s = mlab.contour3d(brain.iso, opacity = opacity, colormap=cmap)
        else:
            s = mlab.contour3d(brain.iso, opacity = opacity, contours = contourVals, colormap=cmap)
            
        # get the object for editing
        self.isosurfacePlots[label] = s
        
    def plotParcels(self, brain, label = None, contourVals = [], opacity = 0.5, cmap='autumn'):
        ''' plot an isosurface using Mayavi, almost the same as skull plotting '''
        
        if not(label):
            label = self.getAutoLabel()        
        
        if contourVals == []:            
            s = mlab.contour3d(brain.parcels, opacity = opacity, colormap=cmap)
        else:
            s = mlab.contour3d(brain.parcels, opacity = opacity, contours = contourVals, colormap=cmap)
            
        # get the object for editing
        self.isosurfacePlots[label] = s
            
    def changePlotProperty(self, plotType, prop, plotLabel, value = 0.):
        ''' change a specified property prop of a plot of type plotType, index is used if multiple plots of
            the same type have been made. Value is used by some properties.
            
            Allowed plotTypes: skull, brainNode, brainEdge
            Allowed props: opacity, visibility, colour
            
            This is basically a shortcut to mayavi visualisation functions 
            
        '''
        

        try:            
            # get plot
            if plotType == 'skull':
                plot = self.skullPlots[plotLabel]
            elif plotType == 'nodes':
                plot = self.brainNodePlots[plotLabel]
            elif plotType == 'edges':
                plot = self.brainEdgePlots[plotLabel]
            else:
                print('plotType not recognised: ' + plotType)
                return
        except:
            # quietly go back if the selected plot doesn't exist
            return
        
        # change plot opacity
        if prop == 'opacity':
            try:
                plot.actor.property.opacity = value
            except:
                print('opacity value not recognised, should be a float', value)
        
        # toggle plot visibility
        elif prop == 'visibility':
            if type(value)!=bool:
                if plot.actor.actor.visibility:
                    plot.actor.actor.visibility = False
                else:
                    plot.actor.actor.visibility = True  
                
            else:
                plot.actor.actor.visibility = value
        # change plot colour
        elif prop == 'colour':
            try:
                plot.actor.property.color = value            
            except:
                print('colour not recognised, should be a triple of values between 0 and 1', value)
                
        else:
            print('property not recognised')
            
        
    def getPlotProperty(self, plotType, prop, plotLabel):
        ''' return the value of a given property for a given plot ''' 
        
        # get plot
        if plotType == 'skull':
            plot = self.skullPlots[plotLabel]
        elif plotType == 'nodes':
            plot = self.brainNodePlots[plotLabel]
        elif plotType == 'edges':
            plot = self.brainEdgePlots[plotLabel]
        else:
            print('plotType not recognised: ' + plotType )
            return
            
        if prop == 'opacity':
            value = plot.actor.property.opacity
        elif prop == 'visibility':
            value =  plot.actor.actor.visibility
        elif prop == 'colour':
            value = plot.actor.property.color
        else:
            print('property not recognised')
            return
            
        return value
        

    def getAutoLabel(self):
        ''' generate an automatic label for a plot object if none given '''
        
        # get index of label
        num = str(self.labelNo)
        num = '0' * (4-len(num)) + num
        
        # make label and print
        label = 'plot ' + num
        print('automatically generated label: '+ label)
        
        # iterate label index
        self.labelNo = self.labelNo + 1
        
        return label