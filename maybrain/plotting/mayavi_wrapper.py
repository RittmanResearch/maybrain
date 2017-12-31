# -*- coding: utf-8 -*-
"""
Plot functions for Maybrain using Mayavi

The MayaviWrapper class creates an object to which many different plots are added.

Examples of plots that can be added are brain networks with highlights, isosurfaces
and other surfaces from nibabel. Various properties of the
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


"""

from mayavi import mlab
from numpy import array, repeat, max, power
from maybrain import constants as ct


class MayaviWrapper:
    """ classes that plot various aspects of a brain object """

    def __init__(self, bg=None):
        """ initialise the Mayavi figure for plotting """

        # start the engine
        # from mayavi.api import Engine
        # self.engine = Engine()
        # self.engine.start()

        if not bg:
            bg = (1., 1., 1.)

        # create a Mayavi figure
        self.mfig = mlab.figure(bgcolor=bg,
                                fgcolor=(0, 0, 0),
                                # engine = self.engine,
                                size=(1500, 1500))

        # holders for plot objects
        self.brainNodePlots = {}
        self.brainEdgePlots = {}
        self.skullPlots = {}
        self.isosurfacePlots = {}

        # record order of plots
        self.plotKeys = []

        # autolabel for plots
        self.labelNo = 0

    def show(self):
        """ show the final plot """
        mlab.show()

    def coords_to_list(self, brain, node_list=None):
        """ get coordinates from lists in the brain object, possibly using only
        the indices given by node_list"""

        # select some rows if necessary
        if not node_list:
            node_list = brain.G.nodes()

            # get coordinates
            # make into array for easy output
        #        coords = array([brain.G.nodes[x][ct.XYZ] for x in nodeList])

        #        # get nodes from networkx object
        #        else:
        coords = []
        for x in brain.G.nodes():
            # put into list
            try:
                coords.append(brain.G.nodes[x][ct.XYZ])
            except KeyError:
                print(('node ' + str(x) + ' not found in function coords_to_list'))

        coords = array(coords)

        # return x, y and z coordinates
        return coords[:, 0], coords[:, 1], coords[:, 2]

    def edges_to_list(self, brain):
        """ Turn the edges of a brain into coordinates """

        # initialise output
        x1 = []
        x2 = []

        y1 = []
        y2 = []

        z1 = []
        z2 = []

        s = []

        for e in brain.G.edges(data=True):
            # get coord and vector from each edge
            p1 = brain.G.nodes[e[0]][ct.XYZ]
            p2 = brain.G.nodes[e[1]][ct.XYZ]

            x1.append(p1[0])
            x2.append(p2[0] - p1[0])

            y1.append(p1[1])
            y2.append(p2[1] - p1[1])

            z1.append(p1[2])
            z2.append(p2[2] - p1[2])

            # set scalar value as edge weight
            s.append(e[2][ct.WEIGHT])

        return x1, y1, z1, x2, y2, z2, s

    def plot_brain(self, brain, opacity=1.0, edge_opacity=None, label='plot', plot_highlights=True):
        """ plot all the coords, edges and highlights in a brain """

        # sort out the edge opacity
        if not edge_opacity:
            edge_opacity = opacity

        # plot coords
        coords = self.coords_to_list(brain)
        self.plot_coords(coords, opacity=opacity, label=label)

        # plot  edges
        ex1, ey1, ez1, ux, uy, yz, s = self.edges_to_list(brain)
        self.plot_edges(ex1, ey1, ez1, ux, uy, yz, s, col=(0., 0., 0.), opacity=opacity, label=label)

        # plot the highlights
        if plot_highlights:
            self.plot_brain_highlights(brain)

    def plot_brain_coords(self, brain, nodes=None, opacity=1.0, label='coordplot', size_list=None, col=(0., 0., 0.),
                          sf=None, sf_range=None):
        """ plot all coordinates in a brain """

        if nodes:
            coords = self.coords_to_list(brain, node_list=nodes)
        else:
            coords = self.coords_to_list(brain)
        self.plot_coords(coords, opacity=opacity, label=label, col=col, size_list=size_list, sf=sf, sf_range=sf_range)

    def plot_brain_edges(self, brain, opacity=1.0, label='edgeplot', col=None, cm='GnBu', lw=2., scalars=None):
        """ plot all edges within a brain
        lw = line width
        #### Not working currently - plots in a different window #### --Really, seems OK to me!
        """

        ex1, ey1, ez1, ux, uy, yz, s = self.edges_to_list(brain)
        if scalars:
            s = scalars
        self.plot_edges(ex1, ey1, ez1, ux, uy, yz, s, lw=lw, col=col, cm=cm, opacity=opacity, label=label)

    def plot_brain_highlights(self, brain, highlights=[], labelpre=''):
        """ plot all or some of the highlights in a brain
            labelpre allow for a standard prefix (this is used by the GUI) """

        if highlights == []:
            highlights = brain.highlights

        # plot highlights (subsets of edges or points)
        for h in highlights:
            label = labelpre + h
            try:
                ho = brain.highlights[h]
            except:
                print(('highlight not found: ' + h))
                continue

            # get edge data                
            ex1, ey1, ez1, ux, uy, yz, s = ho.getEdgeCoordsToPlot(brain)
            # get node data
            hp = ho.nodeIndices

            #            # case where there are node and edges
            #            if (len(ex1)<0 | len(ex1)<0):
            #                print('highlight ' + h + ' has nodes and edges, label will change')
            #                labelEdge = label + '_edge'
            #                labelNode = label + '_node'
            #            else:
            #                labelEdge = label
            #                labelNode = label

            # case where edge opacity is not separately defined
            if not ho.edgeOpacity:
                ho.edgeOpacity = ho.opacity

            # plot the edges
            if not (len(ex1) == 0):
                self.plot_edges(ex1, ey1, ez1, ux, uy, yz, s, col=ho.colour, opacity=ho.edgeOpacity, label=label)

            # plot the nodes
            if not (len(hp) == 0):
                x, y, z = ho.get_coords(brain)
                self.plot_coords((x, y, z), col=ho.colour, opacity=ho.opacity, label=label)

    def plot_coords(self, coords, col=(1., 1., 1.), opacity=1., label='plot', size_list=None, sf=1., sf_range=None):
        """ plot the coordinates of a brain object
            "absoluteScaling" is an option to use and absolute rather than a relative scaling, particularly useful
            for multiple plots
        """

        # remove old version, if any
        if label in self.brainNodePlots:
            self.brainNodePlots[label].remove()

        if size_list is None:
            # note that scalar value is currently set to 1.
            ptdata = mlab.pipeline.scalar_scatter(coords[0], coords[1], coords[2],
                                                  figure=self.mfig)
            sf = 1.

        else:
            try:
                float(size_list)
                size_list = repeat(size_list, len(coords[0]))
            except:
                pass

            if not sf:
                sf = 5. / power(max(size_list), 1 / 3)
                print("sf calculated as: " + str(sf))

            ptdata = mlab.pipeline.scalar_scatter(coords[0], coords[1], coords[2],
                                                  size_list, figure=self.mfig)

        self.brainNodePlots[label] = mlab.pipeline.glyph(ptdata, color=col, opacity=opacity, scale_factor=sf,
                                                         scale_mode="scalar")

        if sf_range:
            print("Adjusting glyph range")
            self.brainNodePlots[label].glyph.glyph.range = array(sf_range)

            # record label for order

            #        self.plotKeys.append(label)

            #        self.brainNodePlots[label] = p
            #        print(label, p)

    def plot_edges(self, ex1, ey1, ez1, ux, uy, uz, s, lw=2., col=None, opacity=1., cm='GnBu', label='plot'):
        """ plot some edges

            ec has the order [x0, x1, y0, y1, z0, z1]
            s is the scalars - used to determine the node colour
            cm = 'GnBu'  colormap for plotting scalars
            col is a triple of numbers in the interval (0,1), or None
            lw is the line width


        """

        if label in self.brainEdgePlots:
            self.brainEdgePlots[label].remove()

        plotmode = '2ddash'  # could be modified later

        # add data to mayavi
        #        edata = mlab.pipeline.vector_scatter(ex1, ey1, ez1, ux, uy, yz, scalars = s, figure = self.mfig)
        # plot
        self.brainEdgePlots[label] = mlab.quiver3d(ex1, ey1, ez1, ux, uy, uz, scalars=s, line_width=lw, opacity=opacity,
                                                   mode=plotmode, color=col, scale_factor=1., scale_mode='vector',
                                                   colormap=cm)
        if not col:
            self.brainEdgePlots[label].glyph.color_mode = 'color_by_scalar'

    def get_coords(self, brain, edge):
        """ get coordinates from nodes and return a coordinate and a vector """

        c1 = brain.G.nodes[edge[0]][ct.XYZ]
        c2 = brain.G.nodes[edge[1]][ct.XYZ]

        diff = [c2[0] - c1[0], c2[1] - c1[1], c2[2] - c1[2]]

        return c1, diff

    def plot_skull(self, brain, label=None, contour_vals=[], opacity=0.1, cmap='Spectral'):
        """ plot the skull using Mayavi """

        if not label:
            label = self.get_auto_label()

            # remove old version
        if label in self.skullPlots:
            self.skullPlots[label].remove()

        if contour_vals == []:
            self.skullPlots[label] = mlab.contour3d(brain.background, opacity=opacity, colormap=cmap)
        else:
            self.skullPlots[label] = mlab.contour3d(brain.background, opacity=opacity, contours=contour_vals,
                                                    colormap=cmap)

            # get the object for editing
            #        self.skullPlots[label] = s

    def plot_isosurface(self, brain, label=None, contour_vals=[], opacity=0.1, cmap='autumn'):
        """ plot an isosurface using Mayavi, almost the same as skull plotting """

        if not label:
            label = self.get_auto_label()

        if label in self.isosurfacePlots:
            self.isosurfacePlots[label].remove()

        if contour_vals == []:
            self.isosurfacePlots[label] = mlab.contour3d(brain.iso, opacity=opacity, colormap=cmap)
        else:
            self.isosurfacePlots[label] = mlab.contour3d(brain.iso, opacity=opacity, contours=contour_vals,
                                                         colormap=cmap)

            # get the object for editing
            #        self.isosurfacePlots[label] = s

    def plot_parcels(self, brain, label=None, contour_vals=[], opacity=0.5, cmap='autumn'):
        """ plot an isosurface using Mayavi, almost the same as skull plotting """

        if not label:
            label = self.get_auto_label()

        if contour_vals == []:
            self.isosurfacePlots[label] = mlab.contour3d(brain.parcels, opacity=opacity, colormap=cmap)
        else:
            self.isosurfacePlots[label] = mlab.contour3d(brain.parcels, opacity=opacity, contours=contour_vals,
                                                         colormap=cmap)

            # get the object for editing
            #        self.isosurfacePlots[label] = s

    def change_plot_property(self, plot_type, prop, plot_label, value=0.):
        """ change a specified property prop of a plot of type plot_type, index is used if multiple plots of
            the same type have been made. Value is used by some properties.

            Allowed plotTypes: skull, brainNode, brainEdge
            Allowed props: opacity, visibility, colour

            This is basically a shortcut to mayavi visualisation functions

        """

        try:
            # get plot
            if plot_type == 'skull':
                plot = self.skullPlots[plot_label]
            elif plot_type == 'nodes':
                plot = self.brainNodePlots[plot_label]
            elif plot_type == 'edges':
                plot = self.brainEdgePlots[plot_label]
            else:
                print(('plot_type not recognised: ' + plot_type))
                return
        except:
            # quietly go back if the selected plot doesn't exist
            return

        # change plot opacity
        if prop == 'opacity':
            try:
                plot.actor.property.opacity = value
            except:
                print(('opacity value not recognised, should be a float', value))

        # toggle plot visibility
        elif prop == 'visibility':
            if type(value) != bool:
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
                print(('colour not recognised, should be a triple of values between 0 and 1', value))

        else:
            print('property not recognised')

    def get_plot_property(self, plot_type, prop, plot_label):
        """ return the value of a given property for a given plot """

        # get plot
        if plot_type == 'skull':
            plot = self.skullPlots[plot_label]
        elif plot_type == 'nodes':
            plot = self.brainNodePlots[plot_label]
        elif plot_type == 'edges':
            plot = self.brainEdgePlots[plot_label]
        else:
            print(('plot_type not recognised: ' + plot_type))
            return

        if prop == 'opacity':
            value = plot.actor.property.opacity
        elif prop == 'visibility':
            value = plot.actor.actor.visibility
        elif prop == 'colour':
            value = plot.actor.property.color
        else:
            print('property not recognised')
            return

        return value

    def get_auto_label(self):
        """ generate an automatic label for a plot object if none given """

        # get index of label
        num = str(self.labelNo)
        num = '0' * (4 - len(num)) + num

        # make label and print
        label = 'plot ' + num
        print(('automatically generated label: ' + label))

        # iterate label index
        self.labelNo = self.labelNo + 1

        return label

    def clear(self):
        """ clear current plot """

        mlab.clf()

        # need to clear other things here
        self.brainNodePlots = {}
        self.brainEdgePlots = {}
        self.skullPlots = {}
        self.isosurfacePlots = {}

        # autolabel for plots
        self.labelNo = 0
