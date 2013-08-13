# -*- coding: utf-8 -*-
"""
Created on Fri Mar 16 18:10:07 2012

@author: -

Change log:
5/2/2013
    - critical bug fix: readSpatial info was assigning coordinates to the wrong nodes
    
7/2/2013
    - Changed readAdjFile to deal with non-whitespace delimiters through the 'delimiter' option
    - Changed readAdjFile to replace string "NA" with 0 from adjacency matrices
    - Added 'weighted' property to readAdjFile to allow a choice of binary or weighted graphs
    - Changed binarise to change the weight of each edge to 1 rather than remove all edge properties
    - Removed contiguousspread and contiguousspreadold functions (superceded by degenerate)
    - Removed the 'symmetric' argument and replaced with directed, which creates a directed graph if required.
        Note that if the graph is undirected, it can not have reciprocal connections eg (2,1) and (1,2)
    - Added clusterslices function to produce 2D images with nodes plotted ***this needs testing***
    - Removed a lost, lonely and wandering dictionary definition from readSpatialInfo
    - Added a function to import isosurfaces, enabling the plotting function to be useful
    - Simplified the spatial positioning function

"""

import os,community
from shutil import move
import networkx as nx
import numpy as np
from networkx.drawing import *
from networkx.algorithms import centrality
from networkx.algorithms import components
import random
from numpy import shape, fill_diagonal, array, where, zeros, sqrt, sort, min, max
#from mayavi import mlab
from string import split
#import nibabel as nb
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
        initialise the brain model. Arguments are:
    
        '''        
        
        # create an empty graph
        self.G = nx.Graph()
        self.iter = None # not sure where this is used. It is often defined, but never actually used!
        
        # initialise global variables
        self.adjMat = None # adjacency matrix, containing weighting of edges.
        self.threshold = 0
        
        self.hubs = []
        self.lengthEdgesRemoved = None
        self.bigconnG = None
        
        self.nbskull = None
        self.skull = None
        self.skullHeader = None
        
        self.nbiso = None
        self.iso = None
        self.isoHeader = None

    ## ================================================

    ## File inputs        

    def readAdjFile(self, fname, threshold = None, edgePC = None, totalEdges = None, directed = False, delimiter=None, weighted=True, NAval="nan", thresholdtype="global", MST=True, excludedNodes=None, delNanNodes=True):
        ''' load adjacency matrix with filename and threshold for edge definition 
        Comment - I don't seem to be able to construct an unweighted directed graph, ie with edge number N(N-1) for an NxN adjacency matrix 
        '''
        
        print('loading data from ' + fname)

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
            lines.append(map(float, [v if not NAval in v else np.nan for v in split(l, sep=delimiter)]))
        nodecount = len(lines)                

        # close file                
        f.close()

        # create adjacency matrix as an array
        # would it be more memory efficient to save the adjacency matrix as a memory map?
        self.adjMat = np.memmap("adjmat.mat",dtype="float32", shape=(len(lines), len(lines)), mode="w+")
        self.adjMat[:] = array(lines)       # array of data
        
        # set excluded node connectivity values to nan
        if excludedNodes:
            for en in excludedNodes:
                self.adjMat[:,en] = np.nan
                self.adjMat[en,:] = np.nan
        try:
            fill_diagonal(self.adjMat, np.nan)
        except: # needed for older version of numpy
            for x in range(len(self.adjMat[0,:])):
                self.adjMat[x,x] = np.nan
                    
        # check if it's diagonal
        sh = shape(self.adjMat)
        if sh[0]!=sh[1]:
            print("Note Bene: adjacency matrix not square.")
            print(sh)

        # create directed graph if necessary
        self.directed = directed
        if directed:
            self.G = nx.DiGraph()
            
        # create nodes on graph
        self.G.add_nodes_from(range(nodecount))  # creates one node for every line in the adjacency matrix
        
        # removes nodes with no connectivity (nan values) in the adjacency matrix
        if delNanNodes:
            for node in self.G.nodes():
                if np.all(np.isnan(self.adjMat[node,:])):
                    self.G.remove_node(node)

        # create edges by thresholding adjacency matrix
        self.adjMatThresholding(edgePC=None, totalEdges=None, MST=False)
        if(threshold or edgePC or totalEdges):
            if thresholdtype=="global":
                self.adjMatThresholding(edgePC, totalEdges, threshold, MST=MST)
    
            elif thresholdtype=="local":
                self.localThresholding(totalEdges, edgePC)
        
        if not weighted:
            self.binarise()
            
                
    
    def NNG(self, k):
        G = nx.Graph()
        nodes = range(len(self.adjMat[0]))
        
        G.add_nodes_from(nodes)
        
        for i in nodes:
            l = np.ma.masked_array(self.adjMat[i,:], mask=np.isnan(self.adjMat[i]))
            l.mask[i] = True
            
            for j in range(k):
                node = np.argmax(l)
                
                if not np.isnan(self.adjMat[i,node]):
                    G.add_edge(i,node)
                    
                l.mask[node] = True
        
        return(G)
        
    def minimum_spanning_edges(self, weight='weight', data=True):
        """Generate edges in a minimum spanning forest of an undirected 
        weighted graph.
    
        A minimum spanning tree is a subgraph of the graph (a tree)
        with the minimum sum of edge weights.  A spanning forest is a
        union of the spanning trees for each connected component of the graph.
    
        Parameters
        ----------
        G : NetworkX Graph
        
        weight : string
           Edge data key to use for weight (default 'weight').
    
        data : bool, optional
           If True yield the edge data along with the edge.
           
        Returns
        -------
        edges : iterator
           A generator that produces edges in the minimum spanning tree.
           The edges are three-tuples (u,v,w) where w is the weight.
        
        Examples
        --------
        >>> G=nx.cycle_graph(4)
        >>> G.add_edge(0,3,weight=2) # assign weight 2 to edge 0-3
        >>> mst=nx.minimum_spanning_edges(G,data=False) # a generator of MST edges
        >>> edgelist=list(mst) # make a list of the edges
        >>> print(sorted(edgelist))
        [(0, 1), (1, 2), (2, 3)]
    
        Notes
        -----
        Uses Kruskal's algorithm.
    
        If the graph edges do not have a weight attribute a default weight of 1
        will be used.
    
        Modified code from David Eppstein, April 2006
        http://www.ics.uci.edu/~eppstein/PADS/
        """
        # Modified code from David Eppstein, April 2006
        # http://www.ics.uci.edu/~eppstein/PADS/
        # Kruskal's algorithm: sort edges by weight, and add them one at a time.
        # We use Kruskal's algorithm, first because it is very simple to
        # implement once UnionFind exists, and second, because the only slow
        # part (the sort) is sped up by being built in to Python.
        from networkx.utils import UnionFind
        if self.G.is_directed():
            raise nx.NetworkXError(
                "Mimimum spanning tree not defined for directed graphs.")
    
        subtrees = UnionFind()
        edges = sorted(self.G.edges(data=True),key=lambda t: t[2][weight], reverse=True)
    #    print edges[0]    
    #    edges = [ v for v in edges if not isnan(v[2]) ]
        
        for u,v,d in edges:
            if subtrees[u] != subtrees[v]:
                if data:
                    yield (u,v,d)
                else:
                    yield (u,v)
                subtrees.union(u,v)
                
    def minimum_spanning_tree(self, weight='weight'):
        """Return a minimum spanning tree or forest of an undirected 
        weighted graph.
    
        A minimum spanning tree is a subgraph of the graph (a tree) with
        the minimum sum of edge weights.
    
        If the graph is not connected a spanning forest is constructed.  A
        spanning forest is a union of the spanning trees for each
        connected component of the graph.
    
        Parameters
        ----------
        G : NetworkX Graph
        
        weight : string
           Edge data key to use for weight (default 'weight').
    
        Returns
        -------
        G : NetworkX Graph
           A minimum spanning tree or forest. 
        
        Examples
        --------
        >>> G=nx.cycle_graph(4)
        >>> G.add_edge(0,3,weight=2) # assign weight 2 to edge 0-3
        >>> T=nx.minimum_spanning_tree(G)
        >>> print(sorted(T.edges(data=True)))
        [(0, 1, {}), (1, 2, {}), (2, 3, {})]
    
        Notes
        -----
        Uses Kruskal's algorithm.
    
        If the graph edges do not have a weight attribute a default weight of 1
        will be used.
        """
        T=nx.Graph(self.minimum_spanning_edges(weight="weight", data=True))
        # Add isolated nodes
        if len(T)!=len(self.G):
            T.add_nodes_from([n for n,d in self.G.degree().items() if d==0])
        # Add node and graph attributes as shallow copy
        for n in T:
            T.node[n] = self.G.node[n].copy()
        T.graph = self.G.graph.copy()
        return T
    
    def localThresholding(self, totalEdges=None, edgePC=None):
        nodecount = len(self.G.nodes())
        
        # get the number of edges to link
        if not edgePC == None:  # needs to be written this way in case edgePC is 0
            # find threshold as a percentage of total possible edges
            # note this works for undirected graphs because it is applied to the whole adjacency matrix
            edgeNum = int(edgePC * nodecount * (nodecount-1) / 2) 
            self.edgePC=edgePC
            
        elif totalEdges:
            # allow a fixed number of edges
            edgeNum = totalEdges
        else:
            edgeNum = -1

        k=1 # number of degrees for NNG
    
        # create minimum spanning tree
        T = self.minimum_spanning_tree(self.G)
        lenEdges = len(T.edges())
        if lenEdges > edgeNum:
            print "The minimum spanning tree already has: "+ lenEdges + " edges, select more edges."
        
        while lenEdges<edgeNum:
            print "NNG degree: "+str(k)
            # create nearest neighbour graph
            nng = self.NNG(k)
            
            # failsafe in case there are no more edges to add
            if len(nng.edges())==0:
                print "There are no edges in the nearest neighbour graph - check you have set the delimiter appropriately"
                break
            
            # remove edges from the NNG that exist already in the new graph/MST
            nng.remove_edges_from(T.edges())
            
            # add weights to NNG    
            for e in nng.edges():
                nng.edge[e[0]][e[1]] = self.G.edge[e[0]][e[1]]
            nng.edges(data=True)
            
            # get a list of edges from the NNG in order of weight
            edgeList = sorted(nng.edges(data=True), key=lambda t: t[2]['weight'], reverse=True)
    
            # add edges to graph in order of connectivity strength
            for edge in edgeList:
                T.add_edges_from([edge])
                lenEdges = len(T.edges())
                if lenEdges >= edgeNum:
                    break
            
            k+=1
        
        self.G = T

      
    def readSpatialInfo(self, fname):
        ''' add 3D coordinate information for each node from a given file '''
        
        try:
            f = open(fname,"rb")
        except IOError, error:
            (errorno, errordetails) = error
            print "Couldn't find 3D position information"
            print "Problem with opening file: "+errordetails                
            return
        
        # get data from file
        lines = f.readlines()
        nodeCount=0
        for line in lines:
            l = split(line)
            try:
                self.G.node[nodeCount]['anatlabel'] = l[0]
                self.G.node[nodeCount]['xyz'] = (float(l[1]),float(l[2]),float(l[3]))
            except:
                print "Node "+str(nodeCount)+" does not exist in graph"
            nodeCount+=1
        del(lines)
                 

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
        
    def inputNodeProperties(self, propertyName, nodeList, propList):
        ''' add properties to nodes, reading from a list of nodes and a list of corresponding properties '''
        
        for ind in range(len(nodeList)):
            n = nodeList[ind]
            p = propList[ind]
            self.G.node[n][propertyName] = p

    
    ## ===========================================================================
    
    ## Functions to alter the brain

    def adjMatThresholding(self, edgePC = None, totalEdges = None, tVal = -1.1, rethreshold=False, doPrint=True, MST=True):
        ''' apply thresholding to the adjacency matrix. This can be done in one of
            three ways (in order of decreasing precendence):
                edgePC - this percentage of nodes will be linked
                totalEdges - this number of nodes will be linked
                tVal - give an absolute threshold value. Pairs of nodes with a corresponding adj matrix
                       value greater than this will be linked. Defaults to -1.1 to produce a fully connected graph.
        '''

        # check if adjacency matrix is there
        if self.adjMat == None:
            print("No adjacency matrix. Please load one.")
            return

        if MST:
            # create minimum spanning tree
            T = self.minimum_spanning_tree(self.G)
 
        if not rethreshold:
            # remove existing edges
            self.G.remove_edges_from(self.G.edges())

        nodecount = len(self.G.nodes())
        print nodecount
            
        # get the number of edges to link
        if not edgePC == None:  # needs to be written this way in case edgePC is 0
            # find threshold as a percentage of total possible edges
            # note this works for undirected graphs because it is applied to the whole adjacency matrix
            edgeNum = int(edgePC * nodecount * (nodecount-1))
            self.edgePC=edgePC
            
        elif totalEdges:
            # allow a fixed number of edges
            edgeNum = totalEdges
        else:
            edgeNum = -1
            
        # get threshold
        if edgeNum>0:
            # get threshold value
            if rethreshold:
                weights = [self.G[v[0]][v[1]]['weight'] for v in self.G.edges()]
            else:
                weights = [v for v in self.adjMat.flatten() if not np.isnan(v)]
            weights.sort()

            try:
                threshold = weights[-edgeNum]
            except IndexError:
                print "Check you are not trying to apply a lower threshold than the current "+str(self.threshold)
                
                if not 'self.threshold' in locals():
                    threshold = -1                
                else:
                    return
            if doPrint:
                print("Threshold set at: "+str(threshold))
        elif edgeNum==0.:
            threshold=1.
        else:
            threshold  = tVal      
            
        self.threshold = threshold
            
        if rethreshold:
            edgesToRemove = [v for v in self.G.edges() if self.G[v[0]][v[1]]['weight'] < threshold]
            self.G.remove_edges_from(edgesToRemove)
        else:
            # carry out thresholding on adjacency matrix
            boolMat = self.adjMat>threshold
            
            edgeCos = where(boolMat) # lists of where edges should be
            
            # display coordinates (for testing only)
    #        for x in range(len(edgeCos[0])):
    #            print(edgeCos[0][x], edgeCos[1][x])
    
            if MST:
                # recreate minimum spanning tree
                self.G = T
            
            edgeCount = len(self.G.edges())
            
            if edgeNum == -1:
                edgeNum = len(self.G.nodes()) * (len(self.G.nodes())-1)
            
            # correct for undirected graph
            edgeNum = edgeNum/2
            
            print ' '.join([str(edgeNum), str(edgeCount) ])
            print(len(edgeCos[0]))
            for ind in range(len(edgeCos[0])):
                node1 = edgeCos[0][ind]
                node2 = edgeCos[1][ind]
                
                # check if edge doesn't already exist and weight isn't nan
                if not(self.G.has_edge(node1, node2)) and not(np.isnan(self.adjMat[node1, node2])):
                    self.G.add_edge(node1, node2, weight = self.adjMat[node1, node2])        
                
                    if MST:
                        edgeCount = len(self.G.edges())
                        if edgeCount >= edgeNum:
                            break
                    
    
    def makeSubBrain(self, propName, value):
        ''' separate out nodes and edges with a certain property '''

        # initialise new brain
        subBrain = brainObj()
        
        # option for making a direct copy
        if value == "any":
            subBrain.G = self.G
        # use some subvalues
        else:
        
            # get nodes from the current list
            acceptedNodes = []
    
            # sort cases for different values
            if type(value) == dict:
                try:
                    v1 = float(value['min'])
                    v2 = float(value['max'])
                except:
                    print 'min and max value not found in makeSubBrain'
    
                       
                for n in self.G.nodes(data = True):
                    try:
                        v = n[1][propName]
                        if (v>=v1)&(v<=v2):
        #                if n[1][propName] == value:
                            subBrain.G.add_nodes_from([n])
                            acceptedNodes.append(n[0])
                    except:
                        continue
                    
            else:
                # make into a list if need be
                if not(type(value))==list:
                    value = [value]
    
                # check nodes to see if property is true
                for n in self.G.nodes(data = True):
                    try:
    #                    if self.G.node[n][propName] in value:
                        if n[1][propName] in value:
                            subBrain.G.add_nodes_from([n])
                            acceptedNodes.append(n[0])
                    except:
                        continue
                
            # add edges from the current brain if both nodes are in the current brain
            for e in self.G.edges():
                if (e[0] in acceptedNodes) & (e[1] in acceptedNodes):
                    subBrain.G.add_edges_from([e])
        
        # construct adjacency matrix in new brain
        subBrain.reconstructAdjMat()
        
        # transfer skull and isosurface
        if self.skull != None:
            subBrain.nbskull = self.nbskull
            subBrain.skull = self.skull
            subBrain.skullHeader = self.skullHeader
        if self.iso != None:
            subBrain.nbiso = self.nbiso
            subBrain.iso = self.iso
            subBrain.isoHeader = self.isoHeader
        
        return subBrain
        
        
    def makeSubBrainEdges(self, propName, value):
        ''' separate out edges with a certain property

            value can be a a single value, or a defined range. The latter is 
            given as a dictionary e.g. {'min': 0, 'max':1}
        
        '''

        # initialise new brain
        subBrain = brainObj()
        subBrain.G.add_nodes_from(self.G.nodes(data=True))
        
        # get nodes from the current list
        acceptedEdges = []

        # sort cases for different values
        if type(value) == dict:
            try:
                v1 = float(value['min'])
                v2 = float(value['max'])
            except:
                print 'min and max value not found in makeSubBrainEdges'

            # check to see if property holds
            for n in self.G.edges(data = True):
                try:
                    v = n[2][propName]
                    if (v>=v1)&(v<=v2):
    #                if n[1][propName] == value:
                        subBrain.G.add_edges_from([n])
                        acceptedEdges.append(n[0])
                except:
                    continue
                
        else:
            # make into a list if need be
            if not(type(value))==list:
                value = [value]

            # check edges to see if property is true
            for n in self.G.edges(data = True):
                try:
#                    if self.G.node[n][propName] in value:
                    if n[2][propName] in value:
                        subBrain.G.add_nodes_from([n])
                        acceptedEdges.append(n[0])
                except:
                    continue

        # construct adjacency matrix in new brain
        subBrain.reconstructAdjMat()
        
        # transfer skull and isosurface
        if self.skull != None:
            subBrain.nbskull = self.nbskull
            subBrain.skull = self.skull
            subBrain.skullHeader = self.skullHeader
        if self.iso != None:
            subBrain.nbiso = self.nbiso
            subBrain.iso = self.iso
            subBrain.isoHeader = self.isoHeader            
                            
        return subBrain    
        
    def makeSubBrainIndList(self, indices):
        ''' Create a subbrain using the given list of nodes. Also taks edges those nodes are in '''

        subbrain = brainObj()        
        for ind in indices:
            n = self.G.nodes(data=True)[ind]
            subbrain.G.add_nodes_from([n])

        # add edges from the current brain if both nodes are in the current brain
        for e in self.G.edges():
            if (e[0] in indices) & (e[1] in indices):
                subbrain.G.add_edges_from([e])                      
            
        # construct adjacency matrix in new brain
        subbrain.reconstructAdjMat()
        
        # transfer skull and isosurface
        if self.skull != None:
            subBrain.nbskull = self.nbskull
            subBrain.skull = self.skull
            subBrain.skullHeader = self.skullHeader
        if self.iso != None:
            subBrain.nbiso = self.nbiso
            subBrain.iso = self.iso
            subBrain.isoHeader = self.isoHeader            
            
        # should also transfer the adjacency matrix here
        return subbrain
        
    def copy(self):
        ''' make an exact copy of this brain '''
        
        newbrain = self.makeSubBrain(None, 'any')
        newbrain.hubs = self.hubs
        newbrain.lengthEdgesRemoved = self.lengthEdgesRemoved
        newbrain.bigconnG = self.bigconnG

        return newbrain

    def randomiseGraph(self, largestconnectedcomp = False):
        self.G = nx.gnm_random_graph(len(self.G.nodes()), len(self.G.edges()))
        if largestconnectedcomp:
            self.bigconnG = components.connected.connected_component_subgraphs(self.G)[0]  # identify largest connected component
            
    def randomremove(self,edgeloss):
        if not self.iter:
            self.iter=0
        try:
            edges_to_remove = random.sample(self.G.edges(), edgeloss)
            self.G.remove_edges_from(edges_to_remove)
            
        except ValueError:
            print "No further edges left"
            
        self.iter += 1
        


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
#    #    combine normalised measures for each node to generate a hub score
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
        
        ## Update adjacency matrix to reflect changes
        #self.reconstructAdjMat()
        
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

        self.adjMatThresholding(edgePC = conVal, doPrint=False)
        conVal -= step
        
        sgLenStart = len(components.connected.connected_component_subgraphs(self.G))
        #print "Starting sgLen: "+str(sgLenStart)
        sgLen = sgLenStart

        while(sgLen == sgLenStart and conVal > 0.):
            self.adjMatThresholding(edgePC = conVal, doPrint=False)
            sgLen = len(components.connected.connected_component_subgraphs(self.G))  # identify largest connected component
            conVal -= step
            # print "New connectivity:" +str(conVal)+ " Last sgLen:" + str(sgLen)
        return conVal+ (2*step)


    def checkrobustnessNew(self, decs = 4, t0low = 1., edgePCBool=False):  #, minThr = None, maxThr = None
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
        if not edgePCBool and self.threshold:
            ths = [t0low, np.mean(np.array((self.threshold, t0low))), self.threshold]
        elif edgePCBool:
            ths = [self.edgePC, np.mean(np.array(((self.edgePC), 0.))), 0.]
        else:
            ths= [t0low, 0.5, 0.]
            
        ths.sort()
        
        if edgePCBool:
            ths.reverse()
            
        # get the connectedness for each threshold value
        sgLen=[]
        for n in range(3):
            if not edgePCBool:
                self.adjMatThresholding(tVal = ths[n], doPrint=False)
            else:
                self.adjMatThresholding(edgePC = ths[n], doPrint=False)
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
                    self.adjMatThresholding(tVal = ths[1], doPrint=False)
                else:
                    self.adjMatThresholding(edgePC = ths[1], doPrint=False)
                sgLenNew = len(components.connected.connected_component_subgraphs(self.G))
                sgLen = [sgLen[0], sgLenNew, sgLen[1]]
                                
            elif sgLen[1]!=sgLen[2]:
                # case where there's a difference between 1st and 2nd components
                
                # get thresholds
                ths = [ths[1], np.mean(np.array((ths[1], ths[2]))), ths[2]]
                
                # get corresponding connectedness
                if not edgePCBool:
                    self.adjMatThresholding(tVal = ths[1], doPrint=False)
                else:
                    self.adjMatThresholding(edgePC = ths[1], doPrint=False)
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
            print("starting threshold error in robustness check, is t0low0 value correct?"+ str(t0low0))


        # reset threshold to old value
        if edgePCBool:
            outThs = self.threshold
        self.threshold = oldThr
        self.adjMatThresholding(tVal = self.threshold, doPrint = False)

        # return the maximum threshold value where they were found to be the same
        fm = "{0:."+str(decs)+"f}"
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
        self.mfig = mlab.figure(bgcolor = (1,1,1), fgcolor = (0,0,0), engine = self.engine, size=(1500, 1500))
        
        # holders for plot objects
        self.brainNodePlots = {}
        self.brainEdgePlots = {}
        self.skullPlots = {}
        self.isosurfacePlots = {}

        # autolabel for plots
        self.labelNo = 0         
        
           
    def nodeToList(self, brain, nodeList=None, edgeList=None):
        ''' convert nodes and edges to list of coordinates '''
        if not nodeList:
            nodeList=brain.G.nodes()
        if not edgeList:
            edgeList=brain.G.edges()
        # get node coordinates into lists        
        nodesX = []
        nodesY = []
        nodesZ = []
        
        for n in nodeList:
            nodesX.append(brain.G.node[n]["xyz"][0])
            nodesY.append(brain.G.node[n]["xyz"][1])
            nodesZ.append(brain.G.node[n]["xyz"][2])
        
        # get edge coordinates and vectors into lists
        edgesCx = []
        edgesCy = []
        edgesCz = []
        edgesVx = []
        edgesVy = []
        edgesVz = []
        
        for e in edgeList:
            c, v = self.getCoords(brain, e)
            edgesCx.append(c[0])
            edgesCy.append(c[1])
            edgesCz.append(c[2])
            
            edgesVx.append(v[0])
            edgesVy.append(v[1])
            edgesVz.append(v[2])
            
            
        
        return nodesX, nodesY, nodesZ, edgesCx, edgesCy, edgesCz, edgesVx, edgesVy, edgesVz
        

    def getNodesList(self, brain, nodeList=None, edgeList=None):
        ''' convert nodes and edges to list of coordinates '''
        if not nodeList:
            nodeList=brain.G.nodes()
        # get node coordinates into lists        
        nodesX = []
        nodesY = []
        nodesZ = []
        
        for n in nodeList:
            nodesX.append(brain.G.node[n]["xyz"][0])
            nodesY.append(brain.G.node[n]["xyz"][1])
            nodesZ.append(brain.G.node[n]["xyz"][2])
                
        return nodesX, nodesY, nodesZ


    def getEdgesList(self, brain, edgeList=None):
        ''' and edges to list of coordinates '''
        if not edgeList:
            edgeList=brain.G.edges()
        
        # get edge coordinates and vectors into lists
        edgesCx = []
        edgesCy = []
        edgesCz = []
        edgesVx = []
        edgesVy = []
        edgesVz = []
        
        for e in edgeList:
            c, v = self.getCoords(brain, e)
            edgesCx.append(c[0])
            edgesCy.append(c[1])
            edgesCz.append(c[2])
            
            edgesVx.append(v[0])
            edgesVy.append(v[1])
            edgesVz.append(v[2])
            
            
        
        return edgesCx, edgesCy, edgesCz, edgesVx, edgesVy, edgesVz


    
    def plotBrain(self, brain, label = None, nodes = None, edges = None, col = (0, 0, 0), opacity = 1.):
        ''' plot the nodes and edges using Mayavi '''
        
        # sort out keywords
        if not nodes:
            nodeList = brain.G.nodes()
        else:
            nodeList = nodes
        if not edges:
            edgeList = brain.G.edges()
        else:
            edgeList = edges
            
        if not(label):
            label = self.getAutoLabel()
            
        # turn nodes into lists for plotting
        xn, yn, zn, xe, ye, ze, xv, yv, zv = self.nodeToList(brain, nodeList=nodeList, edgeList=edgeList)
        
        # plot nodes
        s = mlab.points3d(xn, yn, zn, scale_factor = self.nodesf, color = col, opacity = opacity)
        self.brainNodePlots[label] = s
        
        # plot edges
        t = mlab.quiver3d(xe, ye, ze, xv, yv, zv, line_width = 1., mode = '2ddash', scale_mode = 'vector', scale_factor = 1., color = col, opacity = opacity)
        self.brainEdgePlots[label] = t

    def plotBrainNodes(self, brain, nodes = None, col = (0.5, 0.5, 0.5), opacity = 1., label=None, sizeList=None):
        ''' plot the nodes and edges using Mayavi '''
        
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
        if sizeList:
            s = mlab.points3d(xn, yn, zn, sizeList, scale_factor = self.nodesf, color = col, opacity = opacity)

        else:
            s = mlab.points3d(xn, yn, zn, scale_factor = self.nodesf, color = col, opacity = opacity)
        self.brainNodePlots[label] = s


    def plotBrainEdges(self, brain, label = None, edges = None, col = (0, 0, 0), opacity = 1.):
        ''' plot the nodes and edges using Mayavi '''
        
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
                       
                       
    def plotSubset(self, brain, nodeIndices, edgeIndices, col):
        ''' plot a subset of nodes and edges. Nodes are plotted with colour 'col', a tuple of 3 numbers between 0 and 1, e.g. (0, 0.4, 0.6) '''
        
        nodeSubset = []
        edgeSubset = []
        
        for ind in nodeIndices:
            nodeSubset.append(brain.G.nodes()[ind])
        for ind in edgeIndices:
            edgeSubset.append(brain.G.edges()[ind])
            
        self.plotBrain(nodes = nodeSubset, edges = edgeSubset, col = col)
        
            
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
    
        