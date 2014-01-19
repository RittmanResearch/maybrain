# -*- coding: utf-8 -*-
"""
Created on Fri Mar 16 18:10:07 2012

@author: -

to do:
    - improve memory retention by instant read of file, rather than reloading data
    - write mask for mayavi data??
    - split plotting into separate file (and rename things)
    - set coord scalar value in plotting
    - edge properties from file
    
    
Changes for v0.2:
    - highlight object rather than sub-brain
    - don't change/create networkx or Mayavi objects until needed
    - allow for colouring of edges based on strength

"""

import string,os,csv,community
from shutil import move
import networkx as nx
import numpy as np
from networkx.drawing import *
from networkx.algorithms import centrality
from networkx.algorithms import components
import random
from numpy import shape, fill_diagonal, array, where, zeros, sqrt, sort, min, max, ones, isnan
from mayavi import mlab
from string import split
import nibabel as nb
#from mayavi.core.ui.api import MlabSceneModel, SceneEditor

class brainObj:
    """
    A class that defines a brain network created from an adjacency matrix and spatial inforamtion 
    with certain properties. The added extras are:
        - a list of hub nodes
        - spatial xyz information for each nodes (if input files available)
        - actual length information for each edge
        - layout for plotting
        - counter for iterations of any subsequent process
        
    """
    
    def __init__(self):
        ''' 
        Initialise the brain model.
    
        '''        
        
        # create an empty graph
        self.G = nx.Graph()
        self._directed = False # is this a directed graph or not?
        self.iter = None # not sure where this is used. It is often defined, but never actually used!
        
        # initialise global variables
        self.coords = [] # coordinates of points - should be the same length as dims of adj matrix
        self.adjMat = None # adjacency matrix, containing weighting of edges. Should be square.
        self.adjInds = [] # list of points used from coord/adj data
        self.edgeInds = [] # list of active edges (ordered list of pairs) - only updated if networkx not used
        self.threshold = 0 # value of threshold for including edges
        self.networkxUsed = False # changed to true when self.makeNetwork is run
        
        # need to define the following, what do they do???
        self.hubs = []
        self.lengthEdgesRemoved = None
        self.bigconnG = None
        
        self.nbskull = None
        self.skull = None
        self.skullHeader = None
        
        self.nbiso = None
        self.iso = None
        self.isoHeader = None
        
        self.labelNo = 0 # index for autolabeling of highlights
        self.highlights = {} # highlights items consist of a list contating a name, a highlight object and a colour

    ## ================================================

    ## File inputs        

           

    def importAdjFile(self, fname, delimiter = None):
        ''' get the adjacency data from a file and return as an array '''
        
        # open file
        f = open(fname,"rb")
        reader = f.readlines()        

        # get line that data starts in 
        startLine = 0
        for line in reader:
            if 'begins line' in str(line):
                lstr = str(line)
                whereLabel = lstr.find('begins line')
                startLine = int(lstr[whereLabel + 12])-1
                break
                
        # get data and convert to lists of floats
        linesStr = reader[startLine:]
        lines = []
        for l in linesStr:
            while l[-1] in ('\n', '\t'):
                l = l[:-1]
            lines.append(map(float, [v if v != "NA" else np.nan for v in split(l, sep=delimiter)]))                


        # close file                
        f.close()

        # set adjacency matrix
        self.adjMat = array(lines)
      
      
    def importSpatialInfo(self, fname):
        ''' add 3D coordinate information for each node from a given file '''
        
        try:
            f = open(fname,"rb")
        except IOError, error:
            (errorno, errordetails) = error
            print "Couldn't find 3D position information"
            print "Problem with opening file: "+errordetails                
            return
            
        self.coords = []
        self.anatLabels = []
        
        # get data from file
        lines = f.readlines()
        nodeCount=0
        for line in lines:
            l = split(line)
#            self.G.node[nodeCount]['anatlabel'] = l[0]
#            self.G.node[nodeCount]['xyz'] = (float(l[1]),float(l[2]),float(l[3]))
            self.coords.append([float(l[1]),float(l[2]),float(l[3])])
            self.anatLabels.append(l[0])
            nodeCount+=1
                 

    def importSkull(self, fname):
        ''' Import a file for skull info using nbbabel
            gives a 3D array with data range 0 to 255 for test data
            could be 4d??
            defines an nibabel object, plus ndarrays with data and header info in
        
        '''        
        
        self.nbskull = nb.load(fname)
        self.skull = self.nbskull.get_data()
        self.skullHeader = self.nbskull.get_header()        
                
    def importISO(self, fname):
        ''' Import a file for isosurface info using nbbabel
            gives a 3D array with data range 0 to 255 for test data
            could be 4d??
            defines an nibabel object, plus ndarrays with data and header info in
        
        '''
        self.nbiso = nb.load(fname)
        self.iso = self.nbiso.get_data()
        self.isoHeader = self.nbiso.get_header()
        
    def parcels(self, nodeList):
        """
        Plots 3D parcels specified in the nodeList. This function assumes the parcellation template has
        been loaded to the brain using brain.importISO.
        """        
        zeroArr = zeros(self.iso.shape)
                
        for n in nodeList:
            nArr = np.ma.masked_where(self.iso !=n, self.iso)
            nArr.fill_value=0.0
            zeroArr = zeroArr + nArr
            zeroArr.mask=None
            
        self.parcelList = np.ma.masked_values(zeroArr, 0.0)

    def exportParcelsNii(self, outname='brain'):
        """
        This function saves the parcelList as a nifti file. It requires the
        brain.parcels function has been run first.
        """
        N = nb.Nifti1Image(self.parcelList, self.nbiso.get_affine(), header=self.isoHeader)
        nb.save(N, outname+'.nii')
        
    def importNodeProperties(self, filename):
        ''' add node properties from a file. first lines should contain the property 
            name and the following lines tabulated node indices and property value e.g.:
                
            colour
            1 red
            2 white                
        ''' 
        
        f = open(filename)
        data = f.readlines()
        
        # check that there are enough lines to read
        if len(data)<2:
            print('no data in properties file')
            return
        
        prop = data[0][:-1]       
        
        nodes = []
        propVals = []
        for l in data[1:]:
            try:
                vals = split(l)
                nodes.append(int(vals[0]))
                propVals.append(vals[1])
            except:
                pass

        print(nodes)
        print(propVals)
        print(prop)
        
        self.addNodeProperties(prop, nodes, propVals)
            
        
    def addNodeProperties(self, propertyName, nodeList, propList):
        ''' add properties to nodes, reading from a list of nodes and a list of 
            corresponding properties '''
        
        for ind in range(len(nodeList)):
            n = nodeList[ind]
            p = propList[ind]
            try:
                self.G.node[n][propertyName] = p
            except:
                print('property assignment failed: ' + propertyName + ' ' + str(n) + ' ' + str(p))

    def addEdgeProperty(self, propertyName, edgeList, propList):
        ''' add a property to a selection of edges '''
        
        for ind in range(len(edgeList)):
            e = edgeList[ind]
            p = propList[ind]
            try:
                self.G.edge[e[0]][e[1]][propertyName] = p
            except:
                print('edge property assignment failed: ' + propertyName + ' ' + str(e) + ' ' + str(p))
                    

    
    ## ===========================================================================
    
    ## newtworkx functions    
    
    def createNetworkX(self):
        ''' create a new networkx graph from input data (coords and edges) '''
        
        if self._directed:
            self.G = nx.DiGraph()
        else:
            self.G = nx.Graph()
        
        # make nodes
        try:
            self.G.add_nodes_from(range(len(self.coords)))
        except AttributeError:
            print('brain has no coordinates')
            return
        
        # add coordinates to nodes
        print("adding coords")
        ind = 0
        for c in self.coords:
            self.G.node[ind]['xyz'] = (float(c[0]),float(c[1]),float(c[2]))            
            
            ind = ind + 1

        # add anatomy labels if present
        try:
            print("adding coords")
            ind = 0
            for c in self.coords:
                self.G.node[ind]['anatlabel'] = self.anatLabels[ind]
            ind = ind + 1
        except AttributeError:
            print('brain has no anatlabels')
        
        # make edges
        for e in self.edgeInds:
            # add edge
            self.G.add_edge(e[0],e[1], weight = self.adjMat[e[0],e[1]])
            
        # change networkx bool
        self.networkxUsed = True
        self.edgeInds = None
    

    ## ===========================================================================    
    
    
    ## Functions to alter the brain

    def applyThreshold(self, edgePC = None, totalEdges = None, tVal = -1.1, rethreshold=False):

        #### Determine threshold value
        
        ## case 1: edgePC - percent of edges are shown

        # get the number of edges to link
        if edgePC:
            n = shape(self.adjMat)[0]
            # note this works for undirected graphs because it is applied to the whole adjacency matrix
            if self._directed:
                edgeNum = int(round(edgePC/100. * n * (n-1) * 0.5))                
            else:
                edgeNum = int(round(edgePC/100. * n * (n-1) ))

            self.edgePC=edgePC # !! why  is this necessary?
            print(edgeNum)
            
        ## case 2: totalEdges - give a fixed number of edges
        elif totalEdges:
            # allow a fixed number of edges
            edgeNum = totalEdges
            if not(self._directed):
                edgeNum =  edgeNum * 2

        ## case 3: absolute threshold - show all edges with weighting above threshold
        else:
            edgeNum = -1
        
        
        ##### get threshold value
        if edgeNum>=0:
            # cases 1 and 2 from above - get threshold value for given number of edges
#            if rethreshold:
#                weights = [self.G[v[0]][v[1]]['weight'] for v in self.G.edges()]
#            else:
            weights = self.adjMat.flatten()
            weights.sort()
            while (str(weights[-1]) == 'nan'):
                weights = weights[:-1]
#            weights = [v for v in self.adjMat.flatten() if not str(v)=="nan"]
#            weights.sort()

            # case where all edges are included
            if edgeNum > len(weights):
                self.threshold = weights[0]

            # case where only some edges are included
            else:
                self.threshold = weights[-edgeNum]
#            try:
#                threshold = weights[-edgeNum]
#            except IndexError:
#                threhsold = weights[0]
        else:
            self.threshold  = tVal
        
            
        ##### carry out thresholding on adjacency matrix
        boolMat = self.adjMat>=self.threshold
        try:
            fill_diagonal(boolMat, 0)
        except:
            for x in range(len(boolMat[0,:])):
                boolMat[x,x] = 0
                
        es = where(boolMat) # lists of where edges should be

        # exclude duplicates if not directed 
        if not(self._directed):
            newes1 = []
            newes2 = []
            for ii in range(len(es[1])):
                if es[0][ii]<=es[1][ii]:
                    newes1.append(es[0][ii])
                    newes2.append(es[1][ii])
                    
            es = (newes1, newes2)                   
        
        # add edges to networkx or self.edgeInds
        if self.networkxUsed:
            # could improve the next few lines by spotting when you're only adding edges            
            
            # remove previous edges
            self.G.remove_edges_from(self.G.edges())
            # add new edges
            for ii in range(len(es[0])):
                self.G.add_edge(es[0][ii],es[1][ii], weight = self.adjMat[es[0][ii],es[1][ii]])
            
        else:
            self.edgeInds = []            
            for ii in range(len(es[0])):
                eadd = [es[0][ii], es[1][ii]]
                self.edgeInds.append(eadd)
            
            self.edgeInds.sort()

    def getEdgeWeights(self):
        ''' returns a list of the weighting of each edge '''

        l = []        
        for e in self.edgeInds:
            e.append(self.adjMat(e[0], e[1]))
            
        return l


#    def adjMatThresholding(self, edgePC = None, totalEdges = None, tVal = -1.1, rethreshold=False, doPrint=True):
#        ''' LEGACY FUNCTION - DO NOT USE
#        
#        apply thresholding to the adjacency matrix. This can be done in one of
#            three ways (in order of decreasing precendence):
#                edgePC - this percentage of nodes will be linked
#                totalEdges - this number of nodes will be linked
#                tVal - give an absolute threshold value. Pairs of nodes with a corresponding adj matrix
#                       value greater than this will be linked. Defaults to -1.1 to produce a fully connected graph.
#        '''
#
#        # check if adjacency matrix is there
#        if self.adjMat == None:
#            print("No adjacency matrix. Please load one.")
#            return
#
#        if not rethreshold:
#            # remove existing edges
#            self.G.remove_edges_from(self.G.edges())
#
#        nodecount = len(self.G.nodes())
#            
#        if rethreshold:
#            edgesToRemove = [v for v in self.G.edges() if self.G[v[0]][v[1]]['weight'] < threshold]
#            self.G.remove_edges_from(edgesToRemove)
#        else:
#            
#                
#                if not(self.G.has_edge(node1, node2)):
#                    self.G.add_edge(node1, node2, weight = self.adjMat[node1, node2])        
                
    def highlightFromConds(self, prop, rel, val, label = None, mode = 'edge', colour = (1.,0.,0.)):
        ''' Creates a highlight by asking if the propertly prop is related to val by rel 
        
            type can be 'edge' or 'nodes', to filter for edges or nodes (coordinates) 
            with the given property

            rel can be one of the following strings:
                geq - greater than or equal to
                leq - less than or equal to
                gt - strictly greater than
                lt - stricctly less than
                eq - equal to (i.e. exactly)
                in(), in[), in(], in[] - within an interval, in this case val is a list of two numbers
                contains - val is a string
        '''

        # check filter mode
        if not(mode in ['edge', 'node']):
            print('filter mode not recognised')
            return
        # check if label given
        if not label:
            label = self.getAutoLabel()
            
        # make a highlight object
        h = highlightObj()
        h.colour = colour
        
        # extract lists from edges        
        if mode == 'edge':
            h.edges = []
            ind = -1
            for e in self.G.edges(data = True):
                ind = ind +1
                try:
                    d = self.G.edge[e[0]][e[1]][prop]
                except:
                    continue           
                
                # match properties
                boolval = self.propCompare(d, rel, val)
                
                # save data in highlight
                if boolval:
                    h.edges.append((e[0],e[1])) 
                    
        
        # extract lists from nodes
        elif mode == 'node':
            h.points = []
            for c in range(len(self.G.nodes())):
                # get property

                # special treatment for 'x', 'y' and 'z'
                if prop=='x':
                    d = self.G.node[c]['xyz'][0]
                elif prop=='y':
                    d = self.G.node[c]['xyz'][1]
                elif prop=='z':
                    d = self.G.node[c]['xyz'][2]
                else:
                    # any other property
                    try:
                        d = self.G.node[c][prop]
                    except:
                        continue
                
                # test property against criteria
                boolval = self.propCompare(d, rel, val)
                
                # add to highlight if good
                if boolval:
                    h.points.append(c)
                    
        # add highlight to dictionary
        self.highlights[label] = h
                
        

    def propCompare(self, d, rel, val):
        ''' compare d relative to val, used by highlightFromConds
        
        geq - greater than or equal to
                leq - less than or equal to
                gt - strictly greater than
                lt - stricctly less than
                eq - equal to (i.e. exactly)
                in(), in[), in(], in[] - with
                contains '''
                
        
        if rel == 'eq':
            b = d == val
        elif rel == 'gt':
            b = d > val
        elif rel == 'lt':
            b = d < val
        elif rel == 'leq':
            b = d <= val
        elif rel == 'geq':
            b = d >=val
        elif rel == 'in()':
            b = (d>val[0]) & (d<val[1])
        elif rel == 'in[)':
            b = (d>=val[0]) & (d<val[1])
        elif rel == 'in(]':
            b = (d>val[0]) & (d<=val[1])
        elif rel == 'in[]':
            b = (d>=val[0]) & (d<=val[1])
        elif rel == 'contains':
            b = d in val
        
        return b
            
                            
                
    def makeHighlight(self, edgeInds, coordInds, col, label = None):
        ''' create a highlight object from some edge indices '''
        
        h = highlightObj()
        
        h.edgeIndices = edgeInds
        h.points = coordsInds
        h.colour = col
        
        if not(label):
            label = self.getAutoLabel()
        
        self.highlights[label] = h
   
   
    ## =============================================================
    
    ## Analysis functions

    def reconstructAdjMat(self):
        ''' redefine the adjacency matrix from the edges and weights '''
        n = len(self.G.nodes())
        adjMat = zeros([n,n])
        
        for e in self.G.edges():
            try:
                w = self.G.edge[e[0]][e[1]]['weight']
                adjMat[e[0], e[1]] = w
                adjMat[e[1], e[0]] = w
            except:
                #print("no weight found for edge " + str(e[0]) + " " + str(e[1]) + ", skipped" )
                adjMat[e[0], e[1]] = np.nan

        self.adjMat = adjMat
        
        return adjMat
        
    def updateAdjMat(self, edge):
        ''' update the adjacency matrix for a single edge '''
        
        try:
            w = self.G.edge[edge[0]][edge[1]]['weight']
            self.adjMat[edge[0], edge[1]] = w
            self.adjMat[edge[1], edge[0]] = w
        except:
            print("no weight found for edge " + str(edge[0]) + " " + str(edge[1]) + ", skipped" )
            self.adjMat[edge[0], edge[1]] = np.nan
            self.adjMat[edge[1], edge[0]] = np.nan
        
    
    def findSpatiallyNearest(self, nodeList):
        # find the spatially closest node as no topologically close nodes exist
        print "Finding spatially closest node"
        if isinstance(nodeList, list):
            duffNode = random.choice(nodeList)
        else:
            duffNode = nodeList
            
        nodes = [v for v in self.G.nodes() if v!=duffNode]
        nodes = [v for v in nodes if not v in nodeList]

        shortestnode = (None, None)
        for node in nodes:
            try:
                distance = np.linalg.norm(np.array(self.G.node[duffNode]['xyz'] - np.array(self.G.node[node]['xyz'])))
            except:
                print "Finding the spatially nearest node requires x,y,z values"
                
            if shortestnode[0]:
                if distance < shortestnode[1]:
                    if self.G.degree(node) > 0:
                        shortestnode = (node, distance)
            
            else:
                if self.G.degree(node) > 0:
                    shortestnode = (node, distance)
                    
        return shortestnode[0]
        
    
    def findSpatiallyNearestNew(self, nodeList, threshold=1.):
        ''' find the spatially nearest nodes to each node within a treshold 
        
        Comment - did you have something in mind for this?
        '''
        
        randNode = random.choice(nodeList)
            
        nodes = [v for v in self.G.nodes() if v!=randNode]
        nodes = [v for v in nodes if not v in nodeList]
        
        # get list of node positions
        xyzList = []
        count = 0
        for node in nodes:
            xyzList.append([count] + list(self.G.node[node]['xyz']))
            count = count  + 1

        # cut down in x,y and z coords
        xyz0 = self.G.node[randNode]['xyz']
        xyzmax = [0, xyz0[0] + threshold, xyz0[1] + threshold, xyz0[2] + threshold]
        xyzmin = [0, xyz0[0] - threshold, xyz0[1] - threshold, xyz0[2] - threshold]
        
        # check so that you don't get an empty answer
        count = 0
        countmax = 10
        newxyzList = []
        while (newxyzList==[]) & (count <countmax):           
            # see if it's close to orig point            
            for l in xyzList:
                cond = 1
                # check x, y and z coords
                for ind in [1,2,3]:
                    cond = (l[ind]>xyzmin[ind]) & (l[ind]<xyzmax[ind])
                    
                    if cond==0:
                        break
                    
                # append to new list if close
                if cond:
                    newxyzList.append(l)
                    cond = False
                    
            # increase threshold for next run, if solution is empty
            threshold = threshold * 2

        if newxyzList == []:
            print('unable to find a spatially nearest node')
            return -1, 0
        
        # find shortest distance
        
        # find distances
        dists = []
        print 'newxyzlist'
        print newxyzList
        for l in newxyzList:
            d = sqrt((l[1]-xyz0[0])**2 + (l[2]-xyz0[1])**2 + (l[3]-xyz0[2]**2))
            dists = dists + [(d, l[0])]
        
        print('presort')
        print(dists)
        # sort distances
        dtype = [('d', float), ('ind', int)]
        dists = array(dists, dtype = dtype)
        dists = sort(dists, order = ['d', 'ind'])
        print('postsort')
        print(dists)
        
        # get shortest node
        nodeIndex = dists[0][1]        
        closestNode = self.G.node[nodeIndex]

        return nodeIndex, closestNode        
        
        
    def findLinkedNodes(self):
        ''' give each node a list containing the linked nodes '''
        
        for l in self.G.edges():
            
            # add to list of connecting nodes for each participating node
            try:
                self.G.node[l[0]]['linkedNodes'] = self.G.node[l[0]]['linkedNodes'] + [l[1]]
            except:
                self.G.node[l[0]]['linkedNodes'] = [l[1]]
                        
            try:
                self.G.node[l[1]]['linkedNodes'] = self.G.node[l[1]]['linkedNodes'] + [l[0]]
            except:
                self.G.node[l[1]]['linkedNodes'] = [l[0]]
    
    def hubHelper(self, node):
        hubscore = self.betweenessCentrality[node] + self.closenessCentrality[node] + self.degrees[node]
        return(hubscore)

    def hubIdentifier(self, weighted=False, assign=False):
        """ 
        define hubs by generating a hub score, based on the sum of normalised scores for:
            betweenness centrality
            closeness centrality
            degree
        
        hubs are defined as nodes 2 standard deviations above the mean hub score
        
        defines self.hubs
        
        if assign is true, then each node's dictionary is assigned a hub score
        
        Changelog 7/12/12:
            - added possibility of weighted measures
        """
        
        self.hubs = []
        
    #    get centrality measures
        if weighted:
            self.betweenessCentrality = np.array((centrality.betweenness_centrality(self.G, weight='weight').values()))
            
            ###### this next line doesn't work! Because of negative weights #######
            self.closenessCentrality = np.array((centrality.closeness_centrality(self.G, distance=True).values()))
            self.degrees = np.array((nx.degree(self.G, weight='weight').values()))
            
            
        else:
            self.betweenessCentrality = np.array((centrality.betweenness_centrality(self.G).values()))
            self.closenessCentrality = np.array((centrality.closeness_centrality(self.G).values()))
            self.degrees = np.array((nx.degree(self.G).values()))
            
        self.betweenessCentrality /= np.sum(self.betweenessCentrality)
        self.closenessCentrality /=  np.sum(self.closenessCentrality)        
        self.degrees /= np.sum(self.degrees)
        
        
        # deprecated code follows:
#    #    come normalised measures for each node to generate a hub score
#        hubScores = []
#        for node in self.G.nodes():
#            if weighted:
#                self.G.node[node]['hubScore'] = betweenessCentrality[node]/sum_betweenness + closenessCentrality[node]/sum_closeness + degrees[node]/sum_degrees
#            else:
#                self.G.node[node]['hubScore'] = betweenessCentrality[node]/sum_betweenness + closenessCentrality[node]/sum_closeness + degrees[node]/sum_degrees
#                
#            
#            hubScores.append(self.G.node[node]['hubScore'])

        hubScores = map(self.hubHelper, range(len(self.G.nodes())))
        
        if assign:
            for n,node in enumerate(self.G.nodes()):
                self.G.node['hubscore'] = hubScores[n]
            
    #   find standard deviation of hub score
        upperLimit = np.mean(np.array(hubScores)) + 2*np.std(np.array(hubScores))
    
    #   identify nodes as hubs if 2 standard deviations above hub score
        
        self.hubs = [n for n,v in enumerate(hubScores) if v > upperLimit]
                
    def psuedohubIdentifier(self):
        """ 
        define hubs by generating a hub score, based on the sum of normalised scores for:
            betweenness centrality
            closeness centrality
            degree
            
        hubs are the two 5% connected nodes
        """
        self.hubs = []
        # get degree
        degrees = nx.degree(self.G)
        sum_degrees = np.sum(degrees.values())
        
    #    get centrality measures
        betweenessCentrality = centrality.betweenness_centrality(self.G)
        sum_betweenness = np.sum(betweenessCentrality.values())
        
        closenessCentrality = centrality.closeness_centrality(self.G)
        sum_closeness = np.sum(closenessCentrality.values())
        
    #   calculate the length of 5% of nodes
        numHubs = len(self.G.nodes()) * 0.05
        if numHubs < 1:
            numHubs = 1
            
    #    combine normalised measures for each node to generate a hub score
        hubScores = {}
        for node in self.G.nodes():
            self.G.node[node]['hubScore'] = betweenessCentrality[node]/sum_betweenness + closenessCentrality[node]/sum_closeness + degrees[node]/sum_degrees
            
    #   Check if hub scores is more than those previously collected, and if so replace it in the list of hubs            
            if len(hubScores.keys()) < numHubs:
                hubScores[node] = self.G.node[node]['hubScore']
                
            else:
                minhubScore = np.min(hubScores.values())
                if self.G.node[node]['hubScore'] > minhubScore:
                    minKeyList = [v for v in hubScores.keys() if hubScores[v] == minhubScore]
                    for key in minKeyList:
                        del(hubScores[key])
                    
                    hubScores[node] = self.G.node[node]['hubScore']
                    
    
    #   identify nodes as hubs if 2 standard deviations above hub score
        for node in hubScores.keys():
            self.hubs.append(node)
    
    def clusters(self, drawpos=False):
        """
        Defines clusters using community detection algorithm and adjust layout for pretty presentations
        
        """
        try:
            clusters = community.best_partition(self.G)
#            clusters = community.generate_dendogram(self.G)[0]
                
            # add community and degenerating attributes for each node
            for node in self.G.nodes():
                self.G.node[node]['cluster']=clusters[node] # adds a community attribute to each node
            self.clusternames = set(clusters.values())
        
        except:
            print "Can not assign clusters"
            print "Setting all nodes to cluster 0"
            for node in self.G.nodes():
                self.G.node[node]['cluster'] = 0
            self.clusternames = [0]
            clusters = None
    
#        # set layout position for plotting
#        xy = (0,400)
#        try:
#            angle = 360/len(self.clusternames)
#        except:
#            angle = 180
#        
#        points = {0:xy}
#        for n in range(1,len(self.clusternames)):
#            x = points[n-1][0]
#            y = points[n-1][1]
#            points[n] = (x*np.cos(angle)-y*np.sin(angle),x*np.sin(angle)+y*np.cos(angle))
#        
#        self.pos = {}
#        
#        for clust in self.clusternames:
#            clusternodes = [v for v in self.G.nodes() if self.G.node[v]['cluster']==clust]
#            clusteredges = [v for v in self.G.edges(clusternodes) if v[0] in clusternodes and v[1] in clusternodes]
#            
#            subgraph = nx.Graph()
#            subgraph.add_nodes_from(clusternodes)
#            subgraph.add_edges_from(clusteredges)
#            
#            if drawpos:
#                centre = points[clust]
#                
#                clusterpos = nx_agraph.graphviz_layout(subgraph,prog='neato')
#               
#                for node in clusternodes:
#                    self.pos[node] = (clusterpos[node][0]+centre[0],clusterpos[node][1]+centre[1])
    
        # calculate modularity
        if clusters:
            self.modularity = community.modularity(clusters,self.G)
        else:
            self.modularity = 0
            
    def degenerate(self, weightloss=0.1, edgesRemovedLimit=1, weightLossLimit=None, toxicNodes=None, riskEdges=None, spread=False, updateAdjmat=True):
        ''' remove random edges from connections of the toxicNodes set, or from the riskEdges set. This occurs either until edgesRemovedLimit
        number of edges have been removed (use this for a thresholded weighted graph), or until the weight loss
        limit has been reached (for a weighted graph). For a binary graph, weight loss should be set
        to 1.
        
        The spread option recruits connected nodes of degenerating edges to the toxic nodes list.
        
        By default this function will enact a random attack model, with a weight loss of 0.1 each iteration.
        '''
        
        if toxicNodes:
            nodeList = [v for v in toxicNodes]
        else:
            nodeList = []
            
        # set limit
        if weightLossLimit:
            limit = weightLossLimit
        
        else:
            limit = edgesRemovedLimit
        
        if not riskEdges:
            reDefineEdges=True
            # if no toxic nodes defined, select the whole graph
            if not nodeList:
                nodeList = self.G.nodes()
            
            # generate list of at risk edges
            riskEdges = nx.edges(self.G, nodeList)
        else:
            reDefineEdges=False
            
        # iterate number of steps
        self.lengthEdgesRemoved = []
        while limit>0:
            if not riskEdges:
                # find spatially closest nodes if no edges exist
                # is it necessary to do this for all nodes?? - waste of computing power,
                # choose node first, then calculated spatially nearest of a single node
                newNode = self.findSpatiallyNearest(nodeList)
                if newNode:
                    print "Found spatially nearest node"
                    nodeList.append(newNode)
                    riskEdges = nx.edges(self.G, nodeList)
                else:
                    print "No further edges to degenerate"
                    break
            # choose at risk edge to degenerate from           
            dyingEdge = random.choice(riskEdges)
                        
            # remove specified weight from edge
            w = self.G[dyingEdge[0]][dyingEdge[1]]['weight']
            
            if np.absolute(w) < weightloss:
                loss = w
                self.G[dyingEdge[0]][dyingEdge[1]]['weight'] = 0.
            
            elif w>0:
                loss = weightloss
                self.G[dyingEdge[0]][dyingEdge[1]]['weight'] -= weightloss
                
                
            else:
                loss = weightloss
                self.G[dyingEdge[0]][dyingEdge[1]]['weight'] += weightloss
                
            try:
                self.lengthEdgesRemoved.append(np.linalg.norm( np.array((self.G.node[dyingEdge[0]]['xyz'])) - np.array((self.G.node[dyingEdge[1]]['xyz']))  ))
            except:
                pass                
            
            # update the adjacency matrix (essential if robustness is to be calculated)            
            if updateAdjmat:
                self.updateAdjMat(dyingEdge)
                            
            # add nodes to toxic list if the spread option is selected
            if spread:
                for node in dyingEdge:
                    if not node in nodeList:
                        nodeList.append(node)
            # remove edge if below the graph threshold
            if self.G[dyingEdge[0]][dyingEdge[1]]['weight'] < self.threshold and self.threshold != -1:      # checks that the graph isn't fully connected and weighted, ie threshold = -1
                self.G.remove_edge(dyingEdge[0], dyingEdge[1])
                print ' '.join(["Edge removed:",str(dyingEdge[0]),str(dyingEdge[1])])
                if not weightLossLimit:
                    limit-=1
            
            if weightLossLimit:
                limit -= loss
                
            # redefine at risk edges
            if reDefineEdges:
                riskEdges = nx.edges(self.G, nodeList)
        
#        # Update adjacency matrix to reflect changes
#        self.reconstructAdjMat()
        
        print "Number of toxic nodes: "+str(len(nodeList))
        
        return nodeList
        
    def degenerateNew(self, weightloss=0.1, edgesRemovedLimit=1, weightLossLimit=None, riskNodes=None, riskEdges=None, spread=False):
        ''' remove random edges from connections of the riskNodes set, or from the riskEdges set. This occurs either until edgesRemovedLimit
        number of edges have been removed (use this for a thresholded weighted graph), or until the weight loss
        limit has been reached (for a weighted graph). For a binary graph, weight loss should be set
        to 1.
        
        The spread option recruits connected nodes of degenerating edges to the toxic nodes list.
        
        By default this function will enact a random attack model.
        
        What is the role of riskNodes??
        '''  
        
        # generate list of at-risk nodes
        nodeList = []
        if riskNodes:
            nodeList = [v for v in riskNodes]
        else:
            # if no toxic nodes defined, select a random node
            if nodeList==[]:
                nodeList = [random.choice(self.G.nodes())]
#                nodeList = nodeList + self.G.neighbors(nodeList[0])
#                nodeList = self.G.nodes()

        # generate list of at risk edges                    
        if not(riskEdges):
            riskEdges = nx.edges(self.G, nodeList)
            print('new risk edges', len(riskEdges))

        # set limit in terms of total weight lost or number of edges removed
        if weightLossLimit:
            limit = weightLossLimit        
        else:
            limit = edgesRemovedLimit

        # recording
        deadEdgesRec = [[]]
            
        # iterate number of steps
        while limit>0:
            print(len(nodeList), 'nodes left')
            if not riskEdges:
                # find spatially closest nodes if no edges exist
                # is it necessary to do this for all nodes?? - waste of computing power,
                # choose node first, then calculated spatially nearest of a single node
                newNode = self.findSpatiallyNearest(nodeList) # ugh? why for all nodes??
                if newNode:
                    print "Found spatially nearest node"
                    nodeList.append(newNode)
                    riskEdges = nx.edges(self.G, nodeList)
                else:
                    print "No further edges to degenerate"
                    break
            # choose at risk edge to degenerate from    
            # not sure a random choice is suitable here, should be weighted by the edge weight??
            dyingEdge = random.choice(riskEdges)
            print("edge selected ", dyingEdge)
                        
            # remove specified weight from edge
            w = self.G[dyingEdge[0]][dyingEdge[1]]['weight']            
            # get amount to remove
            if weightloss == 0:
                loss = w
            elif w>0:
                loss = weightloss                
            else:
                # seems a bit weird to have a negative case!!
                loss = -weightloss                
            # remove amount from edge
            new_w = max(w - loss, 0.)
            tBool = new_w<self.threshold
            self.G[dyingEdge[0]][dyingEdge[1]]['weight'] = new_w
            print("old and new weights", w, new_w)
            
            # add nodes to toxic list if the spread option is selected
            if spread:
                for node in dyingEdge:
                    if not (node in nodeList):
                        nodeList.append(node)
                        
            # remove edge if below the graph threshold
            if tBool & (self.threshold != -1):      # checks that the graph isn't fully connected and weighted, ie threshold = -1)                
                self.G.remove_edges_from([dyingEdge])
                riskEdges.pop(riskEdges.index(dyingEdge))
                print ("Edge removed: " + str(dyingEdge[0]) + ' ' + str(dyingEdge[1]) )
                
                # recording
                deadEdgesRec.append([dyingEdge])
        
            # iterate
            if weightLossLimit:
                limit -= loss
            else:
                limit = limit -1
                
            # redefine at risk edges
            # not very efficient, this line
            print(len(riskEdges))
#            riskEdges = nx.edges(self.G, nodeList)

        # Update adjacency matrix to reflect changes
        self.reconstructAdjMat()
        
        print "Number of toxic nodes: "+str(len(nodeList))
        
        return nodeList, deadEdgesRec
        


    def contiguousspread(self, edgeloss, largestconnectedcomp=False, startNodes = None):
        ''' degenerate nodes in a continuous fashion. Doesn't currently include spreadratio '''

        # make sure nodes have the linkedNodes attribute
        try:
            self.G.node[0]['linkedNodes']
        except:
            self.findLinkedNodes()
            
        # make sure all nodes have degenerating attribute
        try:
            self.G.node[0]['degenerating']
        except:
            for n in range(len(self.G.nodes())):
                self.G.node[n]['degenerating']=False 
        
        # start with a random node or set of nodes
        if not(startNodes):
            # start with one random node if none chosen
            toxicNodes = [random.randint(len(self.G.nodes))]
        else:
            # otherwise use user provided nodes
            toxicNodes = startNodes
        # make all toxic nodes degenerating
        for t in toxicNodes:
            self.G.node[t]['degenerating'] = True
                
        # put at-risk nodes into a list
        riskNodes = []
        for t in toxicNodes:
            l = self.G.node[t]['linkedNodes']
            newl = []
            # check the new indices aren't already toxic
            for a in l:
                if a in toxicNodes:
                    continue
                if self.G.node[a]['degenerating']:
                    continue
#                if not(a in toxicNodes)&(not(self.G.node[a]['degenerating'])):
                newl.append(a)

            riskNodes = riskNodes + newl


        
        # iterate number of steps
        toxicNodeRecord = [toxicNodes[:]]
        for count in range(edgeloss):
            # find at risk nodes
            ind = random.randint(0, len(riskNodes)-1)
            deadNode = riskNodes.pop(ind) # get the index of the node to be removed and remove from list
            # remove all instances from list
            while deadNode in riskNodes:
                riskNodes.remove(deadNode)
            
            # add to toxic list    
            toxicNodes.append(deadNode)
            # make it degenerate
            self.G.node[deadNode]['degenerating'] = True
            print('deadNode', deadNode)
            
            
            # add the new at-risk nodes
            l = self.G.node[deadNode]['linkedNodes']
            newl = []
            # check the new indices aren't already toxic
            for a in l:
                if a in toxicNodes:
                    continue
                if self.G.node[a]['degenerating']:
                    continue
                newl.append(a)
                
            riskNodes = riskNodes + newl
            
            toxicNodeRecord.append(toxicNodes[:])
            
            # check that there are any more nodes at risk
            if len(riskNodes)==0:
                break
            
#            print(toxicNodes)
            
        # Update adjacency matrix to reflect changes
        self.reconstructAdjMat()
            
        return toxicNodes, toxicNodeRecord              
            
            
    def neuronsusceptibility(self, edgeloss=1, largestconnectedcomp=False):
        """
        Models loss of edges according to a neuronal suceptibility model with the most highly connected nodes losing
        edges. Inputs are the number of edges to be lost each iteration and the number of iterations.
        """
        self.lengthEdgesRemoved = []
        edgesleft = edgeloss
        if not self.iter:
            self.iter = 0
        
        while edgesleft > 0:
            try:
                # redefine hubs
                self.hubIdentifier()
                
                if self.G.edges(self.hubs) == []:
                    print "No hub edges left, generating pseudohubs"
                    self.psuedohubIdentifier()
                                    
                edgetoremove = random.choice(self.G.edges(self.hubs))
                
                try: # records length of edge removal if spatial information is available
                    self.lengthEdgesRemoved.append(np.linalg.norm(np.array(self.G.node[edgetoremove[0]]['xyz']) - np.array(self.G.node[edgetoremove[1]]['xyz'])))
                    
                except:
                    pass
                    
                self.G.remove_edge(edgetoremove[0],edgetoremove[1])
                edgesleft -= 1
            
            except:
                if self.G.edges(self.hubs) == []:
                    print "No hub edges left, redefining hubs"
                    self.hubIdentifier()

                if self.G.edges(self.hubs) == []:
                    print "No hub edges left, generating pseudohubs"
                    self.psuedohubIdentifier()
                
                if self.G.edges(self.hubs) == []:
                    print "Still no hub edges left, exiting loop"
                    break
                    
                else:
                    continue
        
        self.iter += 1
        
        if largestconnectedcomp:
            self.bigconnG = components.connected.connected_component_subgraphs(self.G)[0]  # identify largest connected component
            
        # Update adjacency matrix to reflect changes
        self.reconstructAdjMat()                      
                      
        
    def percentConnected(self):
        '''
        This will only give the correct results 
        '''
        if self.directed:
            totalConnections = len(self.G.nodes()*(len(self.G.nodes())-1))
        else:
            totalConnections = len(self.G.nodes()*(len(self.G.nodes())-1)) / 2
        self.percentConnections = float(len(self.G.edges()))/float(totalConnections)
        
        return self.percentConnections
    
    def binarise(self):
        '''
            removes weighting from edges 
        '''
        for edge in self.G.edges():
            self.G.edge[edge[0]][edge[1]]['weight'] = 1
        
    def largestConnComp(self):
        self.bigconnG = components.connected.connected_component_subgraphs(self.G)[0]  # identify largest connected component
        
    def strnum(num, length=5):
        ''' convert a number into a string of a given length'''
        
        sn = str(num)
        lenzeros = length - len(sn)
        sn = lenzeros*'0' + sn
        
        return sn
        
    def checkrobustness(self, conVal, step):
        ''' Robustness is a measure that starts with a fully connected graph, \
        then reduces the threshold incrementally until the graph breaks up in \
        to more than one connected component. The robustness level is the \
        threshold at which this occurs. '''

        self.adjMatThresholding(edgePC = conVal)
        conVal -= step
        
        sgLenStart = len(components.connected.connected_component_subgraphs(self.G))
        #print "Starting sgLen: "+str(sgLenStart)
        sgLen = sgLenStart

        while(sgLen == sgLenStart and conVal > 0.):
            self.adjMatThresholding(edgePC = conVal)
            sgLen = len(components.connected.connected_component_subgraphs(self.G))  # identify largest connected component
            conVal -= step
            # print "New connectivity:" +str(conVal)+ " Last sgLen:" + str(sgLen)
        return conVal+ (2*step)


    def checkrobustnessNew(self, decs = 1, t0low = -1., edgePCBool=False):  #, minThr = None, maxThr = None
        ''' Robustness is a measure that starts with a fully connected graph,
        then reduces the threshold incrementally until the graph breaks up in
        to more than one connected component. The robustness level is the
        threshold at which this occurs. 
        
        t0low is the value of the lowest starting threshold
        step is the resolution of the threshold to find
        
        uses a simple 3 point interpolation algorithm 
        
        decs gives the number of decimal places for the output
        
        '''

        oldThr = self.threshold


#        if not(minThr):
#            minThr = min(self.adjMat)
#        if not(maxThr):
#            maxThr = max(self.adjMat)
#            
#            
#        
#        print('min max of adjmat')
#        print(min(self.adjMat), max(self.adjMat))
#                
        # new method
        # get starting threshold values
        if not edgePCBool:
            ths = [t0low, np.mean(np.array((self.threshold, t0low))), self.threshold]
        else:
            ths = [self.edgePC, np.mean(np.array(((self.edgePC), 1.))), 1.]
            
        ths.sort()
        
        if edgePCBool:
            ths.reverse()
            
        # get the connectedness for each threshold value
        sgLen=[]
        for n in range(3):
            if not edgePCBool:
                self.adjMatThresholding(tVal = ths[n])
            else:
                self.adjMatThresholding(edgePC = ths[n])
            sgLen.append(len(components.connected.connected_component_subgraphs(self.G)))
        
        # check to see if tlow0 is too high        
        if sgLen[2] == sgLen[0]:
            print('wrong boundaries for robustness checking')
            print(sgLen, ths)
            return
        
        contBool = 1
        count = 0
        maxCounts = 1000
        while contBool:            
            # compare connectedness of the three values, take the interval in which
            # a difference is found

            if sgLen[0]!=sgLen[1]:
                # case where there's a difference between 0th and 1st components
                
                # set thresholds
                ths = [ths[0], np.mean(np.array((ths[0], ths[1]))), ths[1]]
                
                # get corresponding connectedness
                if not edgePCBool:
                    self.adjMatThresholding(tVal = ths[1])
                else:
                    self.adjMatThresholding(edgePC = ths[1])
                sgLenNew = len(components.connected.connected_component_subgraphs(self.G))
                sgLen = [sgLen[0], sgLenNew, sgLen[1]]
                                
            elif sgLen[1]!=sgLen[2]:
                # case where there's a difference between 1st and 2nd components
                
                # get thresholds
                ths = [ths[1], np.mean(np.array((ths[1], ths[2]))), ths[2]]
                
                # get corresponding connectedness
                if not edgePCBool:
                    self.adjMatThresholding(tVal = ths[1])
                else:
                    self.adjMatThresholding(edgePC = ths[1])
                sgLenNew = len(components.connected.connected_component_subgraphs(self.G))
                
                sgLen = [sgLen[1], sgLenNew, sgLen[2]]                

            else:
                # case where all components are the same (condition satisfied)
#                print(sgLen, ths)
                contBool = 0
            
            # check if final condition satisfied
            if abs(ths[2]-ths[1])< (0.1**decs):
#                print(sgLen, ths)
                contBool = 0

#            print(sgLen, ths)
                
            # stop if too slow
            count = count+1
            if count == maxCounts:
                print('maximum iteration steps exceeded in robustness check')
                print( 'resolution reached: ', str( abs(ths[1]-ths[1]) ) )
                contBool = 0
                
        # check that something happened in the above loop
        if count == 1:
            print("starting threshold error in robustness check, is thrMin value correct?"+ str(t0low0))


        # reset threshold to old value
        if edgePCBool:
            outThs = self.threshold
        self.threshold = oldThr
        self.adjMatThresholding(tVal = self.threshold)

        # return the maximum threshold value where they were found to be the same
        fm = "{:."+str(decs)+"f}"
        if not edgePCBool:
            return fm.format(ths[2])
        else:
            return([fm.format(v) for v in [ths[2],outThs]])
            
    
    def robustness(self, outfilebase="brain", conVal=1.0, decPoints=3, append=True):
        """
        Function to calculate robustness.
        """
        # record starting threhold
        startthresh = self.threshold
        
        if append:
            writeMode="a"
        else:
            if os.path.exists(outfilebase+'_Robustness.txt'):
                move(outfilebase+'_Robustness.txt', outfilebase+'_Robustness.txt.old')
                print ' '.join(["Moving", outfilebase+'_Robustness.txt', "to", outfilebase+'_Robustness.txt.old' ]) 
            writeMode="w"
        
        # iterate through decimal points of connectivity 
        for decP in range(1,decPoints+1):
            step = float(1)/(10**decP)
            #print "Step is: " + str(step)
            conVal = self.checkrobustness(conVal, step)
        
        conVal = conVal-step
        
        if not os.path.exists(outfilebase+'_Robustness.txt'):
            
            log = open(outfilebase+'_Robustness.txt', writeMode)
            log.writelines('\t'.join(["RobustnessConnectivity","RobustnessThresh"])+'\n')
        else:
            log = open(outfilebase+'_Robustness.txt', "a")
        
        log.writelines('\t'.join([str(conVal), str(self.threshold)])+ '\n')
        log.close()
        
        # return graph to original threshold
        self.adjMatThresholding(tVal=startthresh)
        

    def getAutoLabel(self):
        ''' generate an automatic label for a highlight object if none given '''
        
        # get index of label
        num = str(self.labelNo)
        num = '0' * (4-len(num)) + num
        
        # make label and print
        label = 'plot ' + num
        print('automatically generated label: '+ label)
        
        # iterate label index
        self.labelNo = self.labelNo + 1
        
        return label        


class highlightObj():
    ''' object to hold information to highlight a subsection of a brainObj '''
    
    def __init__(self, points = [], edges = []):
        ''' Points refer to the indices of node coordinates in the brain object to which it
        is related. Edges is a set of pairs of coordinates of the same brain object '''
        
#        self._mode = 'pe' # p, e or pe for points only, edges only or points and edges
        self.points = points # indices of points used from a brain object
        self.edges = edges # list of ordered pairs of edges
        self.colour = (0.5, 0.5, 0.5)
        self.opacity = 1.0
        self.edgeOpacity = None
        
        
    def getEdgeCoordsToPlot(self, brain):
        ''' turn list of edges into lists of coordinates - note that the second set of coordinates are the 
        vector from the first to the second points of the edge '''
        
        # initialise outputs
        x1 = []
        x2 = []
        y1 = []
        y2 = []
        z1 = []
        z2 = []
        s = []
        
        print('\n' + 'get edges coords to plot')
        for e in self.edges:
            # get coordinates of edge indices
            p1 = brain.coords[e[0]]
            p2 = brain.coords[e[1]]
            print(p1, p2)
            
            x1.append(p1[0])
            x2.append(p2[0]-p1[0]) 
            
            y1.append(p1[1])
            y2.append(p2[1]-p1[1])
            
            z1.append(p1[2])
            z2.append(p2[2]-p1[2])
            
            s.append(brain.G.edge[e[0]][e[1]]['weight'])
            print(brain.G.edge[e[0]][e[1]]['weight'])
            
        return x1, y1, z1, x2, y2, z2, s
        
    
    def getCoords(self, brain):
        ''' turn indices of points into coordinates '''

        x = []
        y = []
        z = []        
        
        for p in self.points:
            pc = brain.coords[p]
            x.append(pc[0])
            y.append(pc[1])
            z.append(pc[2])
            
        return x, y, z
        
    
class plotObj():
    ''' classes that plot various aspects of a brain object '''
    
    
    def __init__(self):
        
        # initialise mayavi figure
        self.startMayavi()  
        
        self.nodesf = 0.5 # scale factor for nodes
        
        
    def startMayavi(self):
        ''' initialise the Mayavi figure for plotting '''        
        
        # start the engine
        from mayavi.api import Engine
        self.engine = Engine()
        self.engine.start()
        
        # create a Mayavi figure
        self.mfig = mlab.figure(bgcolor = (0, 0, 0), fgcolor = (1., 1., 1.), engine = self.engine, size=(1500, 1500))
        
        # holders for plot objects
        self.brainNodePlots = {}
        self.brainEdgePlots = {}
        self.skullPlots = {}
        self.isosurfacePlots = {}

        # autolabel for plots
        self.labelNo = 0         
        
    def coordsToList(self, brain, nodeList='all'):
        ''' get coordinates from lists in the brain object, possibly using only
        the indices given by nodeList and edgeList '''
        
        coords = array(brain.coords)
        if nodeList!='all':
            coords = coords[nodeList,:]

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
            p1 = brain.coords[e[0]]
            p2 = brain.coords[e[1]]
            
            x1.append(p1[0])
            x2.append(p2[0]-p1[0]) 
            
            y1.append(p1[1])
            y2.append(p2[1]-p1[1])
            
            z1.append(p1[2])
            z2.append(p2[2]-p1[2])
            
            # set scalar value as edge weight
            s.append(e[2]['weight'])
        
        return x1, y1, z1, x2, y2, z2, s
        
        
#                    # get coordinates of edge indices
#            p1 = brain.coords[e[0]]
#            p2 = brain.coords[e[1]]
#            print(p1, p2)
#            
#            x1.append(p1[0])
#            x2.append(p2[0]-p1[0]) 
#            
#            y1.append(p1[1])
#            y2.append(p2[1]-p1[1])
#            
#            z1.append(p1[2])
#            z2.append(p2[2]-p1[2])
#            
#            s.append(brain.G.edge[e[0]][e[1]]['weight'])
#            print(brain.G.edge[e[0]][e[1]]['weight'])
            

    def plotBrain(self, brain, opacity = 1.0, edgeOpacity = None, label = 'plot'):
        ''' plot all the coords, edges and highlights in a brain '''     
               
        print('brainbase')
        self.plotBrainBase(brain, opacity=opacity, edgeOpacity=edgeOpacity, label=label)
        
        print('brainhighlights')
        self.plotBrainHighlights(brain, opacity=edgeOpacity)
        

                
    def plotBrainBase(self, brain, opacity = 1.0, edgeOpacity = None, label='plot'):
        ''' plot all coordinates and edges in a brain '''

        if not(edgeOpacity):
            edgeOpacity = opacity

        # plot coords
        coords = self.coordsToList(brain)
        self.plotCoords(coords, opacity = opacity, label=label)
        
        
        # plot all edges
        ex1, ey1, ez1, ux, uy, yz, s = self.edgesToList(brain)

        self.plotEdges(ex1, ey1, ez1, ux, uy, yz, s, col = (1., 1., 1.), opacity = opacity, label=label)  


    def plotBrainHighlights(self, brain, highlights = [], opacity = 1.0, labelPre = ''):    
        ''' plot all or some of the highlights in a brain '''

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
            # plot the edges
            if not(len(ex1)==0):
                self.plotEdges(ex1, ey1, ez1, ux, uy, yz, s, col = ho.colour, opacity = ho.edgeOpacity, label=label)
            
            # plot the nodes
            hp = ho.points
            if not(len(hp))==0:
                x, y, z = ho.getCoords(brain)
                self.plotCoords((x,y,z), col = ho.colour, opacity = ho.opacity, label=label)        
        
        

    def plotCoords(self, coords, col = (1,1,1), opacity = 1., label='plot'):
        ''' plot the coordinates of a brain object '''
        
        # note that scalar value is currently set to x
        ptdata = mlab.pipeline.scalar_scatter(coords[0], coords[1], coords[2], figure = self.mfig)
        p = mlab.pipeline.glyph(ptdata, color = col, opacity = opacity, scale_factor = 0.1)
        
        self.brainNodePlots[label] = p
        
        
    def plotEdges(self, ex1, ey1, ez1, ux, uy, uz, s, col = None, opacity = 1., label='plot'):
        ''' plot some edges
        
            ec has the order [x0, x1, y0, y1, z0, z1]
            s is the scalars
            col is a triple of numbers in the interval (0,1), or None
            
        '''
        
        plotmode = '2ddash' # could be modified later
        cm = 'GnBu' # colormap for plotting
        
        # add data to mayavi
#        edata = mlab.pipeline.vector_scatter(ex1, ey1, ez1, ux, uy, yz, scalars = s, figure = self.mfig)
        
        # plot
        v = mlab.quiver3d(ex1, ey1, ez1, ux, uy, uz, scalars = s, line_width=1., opacity=opacity, mode = plotmode, color = col, scale_factor = 1., scale_mode = 'vector', colormap = cm)
        if not(col):
            v.glyph.color_mode = 'color_by_scalar'
            
        self.brainEdgePlots[label] = v

            
    def plotBrainOld(self, brain, label = None, nodes = None, edges = None, col = (1, 1, 1), opacity = 1., edgeCol = None):
        ''' CAN SOON BE DEPRECATED
        
            plot the nodes and edges using Mayavi
        
            brain is a brain object
            label is a label used to store the dictinoary
        
            THIS FUNCTION NEEDS SPLITTING UP SOME

            points and edges are stored as different data sets. Hopefully filters
            are applied to plot each of them. 
            
            col is the colour of the plot. Can be a triplet of values  in the interval [0,1],
            or 'vector' to colour by vector scalar.
        
        '''
        
        # sort out keywords
        if not nodes:
            # use all the nodes
            nodeList = range(len(brain.coords))
        else:
            # use a subset of nodes
            nodeList = nodes
        if not edges:
            # use all the edges available
            edgeList = range(len(brain.edgeInds))
        else:
            # use a subset of edges
            edgeList = edges
            
        # get output label for the plot
        if not(label):
            label = self.getAutoLabel()
                        
        # turn nodes into lists for plotting
#        xn, yn, zn, xe, ye, ze, xv, yv, zv, h = self.nodeToList(brain, nodeList=nodeList, edgeList=edgeList)
        xn, yn, zn, xe, ye, ze, xv, yv, zv, h = self.coordsToList(brain, nodeList=nodeList, edgeList=edgeList)
        s = ones(len(xn))
        

        # add point data to mayavi
        ptdata = mlab.pipeline.scalar_scatter(xn, yn, xn, s)

        # add vector data to mayavi
        edata = mlab.pipeline.vector_scatter(xe, xv, ye, yv, ze, zv, scalars = h, figure = self.mfig)
#
        # plot nodes
        mlab.pipeline.glyph(ptdata, color = col, opacity = opacity)
#
        # plot edges
        v = mlab.pipeline.vectors(edata, colormap='GnBu', line_width=1., opacity=opacity, mode = '2ddash', color = edgeCol)
        if not(edgeCol):
            v.glyph.color_mode = 'color_by_scalar'

        
        # plot nodes
#        s = mlab.points3d(xn, yn, zn, scale_factor = self.nodesf, color = col, opacity = opacity)
#        self.brainNodePlots[label] = s
        
        # plot edges
#        t = mlab.quiver3d(xe, ye, ze, xv, yv, zv, line_width = 1., mode = '2ddash', scale_mode = 'vector', scale_factor = 1., color = col, opacity = opacity)
#        self.brainEdgePlots[label] = t
        


    def plotBrainNodes(self, brain, nodes = None, col = (1, 1, 1), opacity = 1., label=None):
        ''' plot the nodes using Mayavi TO BE DEPRECATED'''
        
        # sort out keywords
        if not nodes:
            nodeList = brain.G.nodes()
        else:
            nodeList = nodes
            
        if not(label):
            label = self.getAutoLabel()            
            
        # turn nodes into lists for plotting
        xn, yn, zn = self.getNodesList(brain, nodeList=nodeList)
        
        # plot nodes
        s = mlab.points3d(xn, yn, zn, scale_factor = self.nodesf, color = col, opacity = opacity)
        self.brainNodePlots[label] = s


    def plotBrainEdges(self, brain, edges = None, col = (1, 1, 1), opacity = 1., label = 'plot'):
        ''' plot the nodes and edges using Mayavi TO BE DEPRECATED'''
        
        # sort out keywords
        if not edges:
            edgeList = brain.G.edges()
        else:
            edgeList = edges
            
        if not(label):
            label = self.getAutoLabel()            
            
        # turn nodes into lists for plotting
        xe, ye, ze, xv, yv, zv = self.getEdgesList(brain, edgeList=edgeList)
                
        # plot edges
        t = mlab.quiver3d(xe, ye, ze, xv, yv, zv, line_width = 1., mode = '2ddash', scale_mode = 'vector', scale_factor = 1., color = col, opacity = opacity)
        self.brainEdgePlots[label] = t
                       
                       
    def plotSubset(self, brain, hl, col):
        ''' plot a subset of nodes and edges. Nodes are plotted with colour 'col', a tuple of 3 numbers between 0 and 1, e.g. (0, 0.4, 0.6) 
            TO BE DEPRECATED '''
        

        # get indices of edges
        edgeInds = hl.getEdgeInds(brain)
            
        # plot the highlight
        self.plotBrain(brain, nodes = hl.points, edges = edgeInds, edgeCol = col)
        
            
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
            s = mlab.contour3d(brain.skull, opacity = opacity, colormap=cmap)
        else:
            s = mlab.contour3d(brain.skull, opacity = opacity, contours = contourVals, colormap=cmap)
            
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
            
            THIS IS GOING TO NEED SOME CHANGES
            
        '''
        

        try:            
            # get plot
            if plotType == 'skull':
                plot = self.skullPlots[plotLabel]
            elif plotType == 'brainNode':
                plot = self.brainNodePlots[plotLabel]
            elif plotType == 'brainEdge':
                plot = self.brainEdgePlots[plotLabel]
            else:
                print 'plotType not recognised'
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
        elif plotType == 'brainNode':
            plot = self.brainNodePlots[plotLabel]
        elif plotType == 'brainEdge':
            plot = self.brainEdgePlots[plotLabel]
        else:
            print('plotType not recognised')
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
    
        
