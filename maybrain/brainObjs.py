# -*- coding: utf-8 -*-
"""
Created on Fri Mar 16 18:10:07 2012

@authors: Timothy Rittman, Martyn Rittman

Documentation available at http://code.google.com/p/maybrain/

This is a merge version. Notes prepended with #!!

To do:
Move nibabel functions into another script
Add in link to brain connectivity toolbox
Date of conference: 3rd June
change skull to background
pysurfer??



"""
import networkx as nx
import numpy as np
#from networkx.drawing import *
from networkx.algorithms import centrality
from networkx.algorithms import components
import random
from string import split
import extraFns

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
    
    def __init__(self, directed=False):
        ''' 
        Initialise the brain model.
    
        '''        
        
        # create an empty graph
        self.directed = directed # is this a directed graph or not?
        
        # initialise global variables
        self.adjMat = None # adjacency matrix, containing weighting of edges. Should be square.
        self.threshold = 0 # value of threshold for including edges
        
        # need to define the following, what do they do???
#        self.edgePC = 5 # an arbitrary value for percentage of edges to plot. #!! is this necessary?
        self.hubs = []
        self.lengthEdgesRemoved = None
        self.bigconnG = None
        #!! following two added in merge
        self.dyingEdges = {}
        self.nodesRemoved = None
        
        # background info imported by nibabel
        self.nbbackground = None # the nibabel object
        self.background = None # coordinates of background
        self.backgroundHeader = None # header of nibabel data
        
        # isosurface information imported by nibabel
        self.nbiso = None # all the isosurface info, nibabel object
        self.iso = None # the isosurface
        self.isoHeader = None # header information
        
        self.labelNo = 0 # index for autolabeling of highlights
        self.highlights = {} # highlights items consist of a list contating a name, a highlight object and a colour

        # create a new networkX graph object
        if self.directed:
            # directional case
            self.G = nx.DiGraph()
        else:
            # non-directional case
            self.G = nx.Graph()

    ## ================================================

    ##### File inputs and outputs

    ### edges and nodes

    #!! readAdjFile removed           
    def importAdjFile(self, fname, delimiter = None, exclnodes=[], naVals=["NA"]):
        ''' get the adjacency data from a file and return as an array '''
        self.exclnodes=exclnodes
        
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
            lines.append(map(float, [v if v not in naVals else np.nan for v in split(l, sep=delimiter)]))

        # close file                
        f.close()

        # set adjacency matrix
        self.adjMat = np.array(lines)
        
        # add nodes
        self.G.add_nodes_from([v for v in range(len(lines)) if not v in exclnodes])

        # update adjacency matrix to null values of excluded nodes
        if exclnodes:
            for en in exclnodes:
                self.adjMat[:,en]=np.nan
                self.adjMat[en,:]=np.nan

    #!! After much deliberation, this function was kept from the master, but renamed
    def importSpatialInfo(self, fname, delimiter=None, convertMNI=False):
        ''' add 3D coordinate information for each node from a given file
            note that the graph object is recreated here by default
        '''
        # open file
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
            if delimiter:
                l = split(line,sep=delimiter)
            else:
                l = split(line)

            if convertMNI:
                l[1] = 45 - (float(l[1])/2)
                l[2] = 63 + (float(l[2])/2)
                l[3] = 36 + (float(l[3])/2)

            if nodeCount in self.G.nodes():
                self.G.node[nodeCount]["xyz"]=(float(l[1]),float(l[2]),float(l[3]))
                self.G.node[nodeCount]["anatlabel"]=l[0]

            nodeCount+=1
                        
    #!! in merge import Properties taken from dev2
    def importProperties(self, filename):
        ''' add properties from a file. first lines should contain the property 
            name and the following lines tabulated node indices and property value e.g.:
                
            colour
            1 red
            2 white    
            2   4   green
            
            if 2 indices are given, the property is applied to edges instead. Can mix nodes and edges in the same file.
            
            
        ''' 
        
        # load file from data
        f = open(filename)
        data = f.readlines()
        
        # check that there are enough lines to read
        if len(data)<2:
            print('no data in properties file')
            return
        
        # get the proprety name
        prop = data[0][:-1]       

        # initialise output        
        nodes = []
        propValsNodes = []
        edges = []
        propValsEdges = []
        
        # extract data from file
        for l in data[1:]:
            # determine if it should apply to edges or nodes (by length)
            try:
                print(l)
                value = split(l)
                if len(value)==2:
                    mode = 'nodes'
                elif len(value)==3:
                    mode =  'edges'
            except:
                print('couldn\'t parse data line. Please check property file', l)
            
            # get data
            if mode=='nodes':
                nodes.append(int(value[0]))
                propValsNodes.append(value[1])
                
            elif mode=='edges':
                edges.append([int(value[0]),int(value[1])])
                propValsEdges.append(value[2])
        
        # add data to brain
        if len(nodes)>0:
            self.addNodeProperties(prop, nodes, propValsNodes)
        if len(edges)>0:
            self.addEdgeProperty(prop, edges, propValsEdges)

        return mode, prop # required for GUI
                    
        
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

        
    ### supplementary structures

    #!! rename background to template
    def importBackground(self, fname):
        import nibabel as nb
        ''' Import a file for background info using nbbabel
            gives a 3D array with data range 0 to 255 for test data
            could be 4d??
            defines an nibabel object, plus ndarrays with data and header info in
        
        '''        
        
        self.nbbackground = nb.load(fname)
        self.background = self.nbbackground.get_data()
        self.backgroundHeader = self.nbbackground.get_header()        
                
    def importISO(self, fname):
        ''' Import a file for isosurface info using nibabel
            gives a 3D array with data range 0 to 255 for test data
            could be 4d??
            defines an nibabel object, plus ndarrays with data and header info in
        
        '''
        import nibabel as nb
        self.nbiso = nb.load(fname)
        self.iso = self.nbiso.get_data()
        self.isoHeader = self.nbiso.get_header()
        
    #!! in merge, parcels taken from master branch
    def parcels(self, nodeList):
        """
        Plots 3D parcels specified in the nodeList. This function assumes the 
        parcellation template has been loaded to the brain using brain.importISO.
        Note, values passed to this function should corresond with those in the 
        iso image, not necessarily the node values.
        """
        
        zeroArr = np.zeros(self.iso.shape)
        
        for n in nodeList:
            n = float(n)
            nArr = np.ma.masked_where(self.iso !=n, self.iso)
            nArr.fill_value=0.0
            zeroArr = zeroArr + nArr
            zeroArr.mask=None
            
        self.parcelList = np.ma.masked_values(zeroArr, 0.0)

    def exportParcelsNii(self, outname='brain', valueDict=None):
        """
        This function saves the parcelList as a nifti file. It requires the
        brain.parcels function has been run first.
        """
        import nibabel as nb
        if valueDict: # creates a numpy array based on the dictionary provided
            outMat = np.zeros(self.nbiso.get_data().shape, dtype="float64")
            for n in valueDict.keys():
                outMat[np.where(self.nbiso.get_data()==n+1)] = valueDict[n]
        else:
            outMat = self.parcelList
        
        N = nb.Nifti1Image(outMat, self.nbiso.get_affine(), header=self.isoHeader)

        nb.save(N, outname+'.nii')
                     
        

    ## ===========================================================================    
    
    
    ##### Functions to alter the brain

    ### adjacency matrix and edge thresholding

    #!! need to create minimum spanning tree option??
    #!! doPrint option removed in merge
    def applyThreshold(self, edgePC = None, totalEdges = None, tVal = None, rethreshold=False):
        ''' Treshold the adjacency matrix to determine which nodes are linked by edges. There are 
            three options:
            
            edgePC - a percentage of edges are taken
            totalEdges - give a number for the total edges found
            tVal - give a value for the threshold at which edges are chosen 
            
            rethreshold can be used if the threshold has already been applied --- NOT CURRENTLY A LIVE OPTION'''

        #### Determine threshold value
        weights = self.adjMat.flatten()
        
        ## cases 1 and 2: edgePC - percent of edges are shown
        ##                totalEdges - a number of edges is given
        
        # get the number of edges to link
        if tVal == None:  
            # find threshold as a percentage of total possible edges
            # note this works for undirected graphs because it is applied to the whole adjacency matrix
        
            # get a list of weights from the adjacency matrix
            if not self.directed:
                # only flatten the upper right part of the matrix
                weights = extraFns.undirectedFlatten(self.adjMat)
            else:
                weights = list(self.adjMat.flatten())
            print(weights)
            weights.sort()
            
            # remove NaNs from end
            while (str(weights[-1]) == 'nan'):
                weights = weights[:-1]
            nEdges = len(weights)
            
            # percentage case
            if not (edgePC == None):
                # get number of edges
                edgeNum = int(edgePC/100. * nEdges)
            else:
            # case where number of edges given
            
            # not sure if the following if is really needed. 
#                if not self.directed:
#                    edgeNum = totalEdges * 2
#                else:
                edgeNum = totalEdges
            
            # get threshold value
            if edgeNum > len(weights):
                # case where all edges are included
                self.threshold = weights[0]
            else:
                # case where some edges are included
                self.threshold = weights[-edgeNum]
            
            print edgeNum, self.threshold             
#            self.edgePC=edgePC # never used

        
        # case 3 - absolute threshold
        elif tVal:
            self.threshold = tVal
        
        # case 4 - weighted graph, no threshold specified
        else:
            self.threshold  = np.min(weights[~np.isnan(weights)])
        
        print self.threshold
            
        ##### carry out thresholding on adjacency matrix
        boolMat = self.adjMat>=self.threshold
        try:
            np.fill_diagonal(boolMat, 0)
        except:
            for x in range(len(boolMat[0,:])):
                boolMat[x,x] = 0
                
        es = np.where(boolMat) # lists of where edges should be

        # exclude duplicates if not directed 
        if not(self.directed):
            newes1 = []
            newes2 = []
            for ii in range(len(es[1])):
                if es[0][ii]<es[1][ii]:
                    newes1.append(es[0][ii])
                    newes2.append(es[1][ii])
                    
            es = (newes1, newes2)                   
        # add edges to networkx

        # could improve the next few lines by spotting when you're only adding edges            
        
        # remove previous edges
        self.G.remove_edges_from(self.G.edges())
        # add new edges
        for ii in range(len(es[0])):
            self.G.add_edge(es[0][ii],es[1][ii], weight = self.adjMat[es[0][ii],es[1][ii]])
        

    def reconstructAdjMat(self):
        ''' redefine the adjacency matrix from the edges and weights '''
        n = len(self.G.nodes())
        adjMat = np.zeros([n,n])
        
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
            

    #!! added from master is this in the right place?? 
    def localThresholding(self, totalEdges=None, edgePC=None, removeUnconnected=True):
        '''
        Threshold the association matrix by building from the minimum spanning
        tree and adding successive N-nearest neighbour degree graphs.        
        '''
        self.applyThreshold()
        if removeUnconnected:
            self.removeUnconnectedNodes()
        # get the number of edges to link
        if not edgePC == None:  # needs to be written this way in case edgePC is 0
            # find threshold as a percentage of total possible edges
            # note this works for undirected graphs because it is applied to the whole adjacency matrix
            weights = self.adjMat.flatten()
            weights.sort()
            while (str(weights[-1]) == 'nan'):
                weights = weights[:-1]
            nEdges = len(weights)
            
            # note this works for undirected graphs because it is applied to the whole adjacency matrix
            edgeNum = int(edgePC/100. * nEdges)
                
            # correct for diagonal if graph is undirected
            if not self.directed:
                edgeNum -= len([self.adjMat[x,x] for x in range(len(self.adjMat[0,:])) if not np.isnan(self.adjMat[x,x])])
            print edgeNum

            self.edgePC=edgePC
            
        elif totalEdges:
            # allow a fixed number of edges
            edgeNum = totalEdges
        else:
            edgeNum = -1

        k=1 # number of degrees for NNG
    
        # create minimum spanning tree
        T = self.minimum_spanning_tree(self)
        print nx.is_connected(T)
        lenEdges = len(T.edges())
        if lenEdges > edgeNum:
            print "The minimum spanning tree already has: "+ str(lenEdges) + " edges, select more edges."
        
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
                nng.edge[e[0]][e[1]]['weight'] = self.adjMat[e[0],e[1]]
            
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
        
        #!! are you sure you want to redefine G here??
        self.G = T
        
        
    def binarise(self):
        '''
            removes weighting from edges 
        '''
        for edge in self.G.edges():
            self.G.edge[edge[0]][edge[1]]['weight'] = 1        

    def removeUnconnectedNodes(self):
        '''
         Remove nodes with no connections
        '''
        nodeList = [v for v in self.G.nodes() if self.G.degree(v)==0]
        self.G.remove_nodes_from(nodeList)

    ### making highlights

                
    def highlightFromConds(self, prop, rel, val, label = None, mode = 'edge', colour = (1.,0.,0.), opacity = 1.0):
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
        h.opacity = opacity
        
        # extract lists from edges        
        if mode == 'edge':
            h.edgeIndices = []
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
                    h.edgeIndices.append((e[0],e[1])) 
                    
        
        # extract lists from nodes
        elif mode == 'node':
            h.nodeIndices = []
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
                    h.nodeIndices.append(c)

        
        # add highlight to dictionary
        self.highlights[label] = h
                
                
    def makeHighlight(self, edgeInds, coordInds, col, label = None):
        ''' create a highlight object from some edge indices '''
        
        h = highlightObj()
        
        h.edgeIndices = edgeInds
        h.nodeIndices = coordInds
        h.colour = col
        
        if not(label):
            label = self.getAutoLabel()
        
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
        else:
            print('relation not recognised: ' + rel )
        
        return b        


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


    ### edge removal/degeneration
            
    def degenerate(self, weightloss=0.1, edgesRemovedLimit=1, threshLimit=None,
                   pcLimit=None, weightLossLimit=None, toxicNodes=None,
                   riskEdges=None, spread=False, updateAdjmat=True,
                   distances=False, spatialSearch=False):
        ''' remove random edges from connections of the toxicNodes set, or from the riskEdges set. This occurs either until edgesRemovedLimit
        number of edges have been removed (use this for a thresholded weighted graph), or until the weight loss
        limit has been reached (for a weighted graph). For a binary graph, weight loss should be set
        to 1.
        
        The spread option recruits connected nodes of degenerating edges to the toxic nodes list.
        
        By default this function will enact a random attack model, with a weight loss of 0.1 each iteration.
        
        Weights are taken as absolute values, so the weight in any affected edge tends to 0.    
        
        Spread can either be False, or a number specifying the weight above which to add
        nodes within the list of at-risk nodes.
        
        '''
        
        if toxicNodes:
            nodeList = [v for v in toxicNodes]
        else:
            nodeList = []
            
        # set limit
        if weightLossLimit and pcLimit:
            print "You have asked for both a weight and percentage connectivity limit, using the percentage connectivity limit"
        
        if threshLimit:
            pcLimit = self.thresholdToPercentage(threshLimit)
        
        if pcLimit:
            lenNodes = len(self.G.nodes())
            lenEdges = len(self.G.edges())
            
            maxEdges = float(lenNodes * (lenNodes-1))
            if not self.G.is_directed():
                maxEdges = maxEdges / 2
                
            newEdgeNum = int(round(pcLimit * maxEdges))
            if newEdgeNum > lenEdges:
                print "The percentage threshold set is lower than the current graph, please choose a larger value"
            
            limit = lenEdges - newEdgeNum
            weightLossLimit = False
            
        elif weightLossLimit:
            limit = weightLossLimit
        
        else:
            limit = edgesRemovedLimit
        
        if not riskEdges:
            reDefineEdges=True
            # if no toxic nodes defined, select the whole graph
            if not nodeList:
                nodeList = self.G.nodes()
            
            # generate list of at risk edges
            riskEdges = [v for v in nx.edges(self.G, nodeList) if self.G.edge[v[0]][v[1]]['weight'] != 0.]
        else:
            reDefineEdges=False
            
        if spread:
            nodeList = []
            
        # iterate number of steps
        self.lengthEdgesRemoved = []
        
        # check if there are enough weights left
        riskEdgeWtSum = np.sum([self.G.edge[v[0]][v[1]]['weight'] for v in riskEdges])
        if limit > riskEdgeWtSum:
            print "Not enough weight left to remove"
            return nodeList
            
        
        while limit>0.:
            if not riskEdges and spatialSearch:
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
                loss = np.absolute(w)
                self.G[dyingEdge[0]][dyingEdge[1]]['weight'] = 0.
                riskEdges.remove(dyingEdge)
            
            elif w>0:
                loss = weightloss
                self.G[dyingEdge[0]][dyingEdge[1]]['weight'] -= weightloss
                
            else:
                loss = weightloss
                self.G[dyingEdge[0]][dyingEdge[1]]['weight'] += weightloss
            
            # record the edge length of edges lost
            if distances:
                self.dyingEdges[dyingEdge] = self.G[dyingEdge[0]][dyingEdge[1]]
                self.dyingEdges[dyingEdge]['distance'] =  np.linalg.norm( np.array((self.G.node[dyingEdge[0]]['xyz'])) - np.array((self.G.node[dyingEdge[1]]['xyz']))  )
            
            # update the adjacency matrix (essential if robustness is to be calculated)            
            if updateAdjmat:
                self.updateAdjMat(dyingEdge)
                            
            # add nodes to toxic list if the spread option is selected
            if spread:
                for node in dyingEdge:
                    if not node in nodeList and self.G.edge[dyingEdge[0]][dyingEdge[1]] > spread:
                        nodeList.append(node)
                        
            # remove edge if below the graph threshold
            if self.G[dyingEdge[0]][dyingEdge[1]]['weight'] < self.threshold and self.threshold != -1:      # checks that the graph isn't fully connected and weighted, ie threshold = -1
                self.G.remove_edge(dyingEdge[0], dyingEdge[1])
                riskEdges.remove(dyingEdge)
                print ' '.join(["Edge removed:",str(dyingEdge[0]),str(dyingEdge[1])])
                if not weightLossLimit:
                    limit-=1
            
            if weightLossLimit:
                limit -= loss
                
            # redefine at risk edges
            if reDefineEdges or spread:
                riskEdges = nx.edges(self.G, nodeList)

        
        
        ## Update adjacency matrix to reflect changes
        #self.reconstructAdjMat()
        
        print "Number of toxic nodes: "+str(len(nodeList))
        
        return nodeList

    #!! degenerateNew function removed

    def contiguousspread(self, edgeloss, startNodes = None):
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
            toxicNodes = [random.randint(0, len(self.G.nodes()))]
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
            
    #!! neuronsusceptibility removed            


    #!! added from master branch. Which other function uses this one? Can it be removed (randomly, or otherwise?)            
    def randomremove(self,edgeloss):
        if not self.iter:
            self.iter=0
        try:
            edges_to_remove = random.sample(self.G.edges(), edgeloss)
            self.G.remove_edges_from(edges_to_remove)
            
        except ValueError:
            print "No further edges left"
            
        self.iter += 1  

    ### other modifying functions

    #!! added direct from master branch        
    def randomiseGraph(self, largestconnectedcomp = False):
        self.G = nx.gnm_random_graph(len(self.G.nodes()), len(self.G.edges()))
        if largestconnectedcomp:
            self.bigconnG = components.connected.connected_component_subgraphs(self.G)[0]  # identify largest connected component


        
   
    # =============================================================
    
    ##### Analysis functions

    ### Minimum spanning tree functions

    #!! This and following functions added during merge. Need explanation.
    def NNG(self, k):
        ''' '''
        #!! docstring missing
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


    ### basic proximities
   
    def findSpatiallyNearest(self, nodeList, contra=False, midline=44.5, connected=True):
        # find the spatially closest node as no topologically close nodes exist
        if isinstance(nodeList, list):
            duffNode = random.choice(nodeList)
        else:
            duffNode = nodeList
            
        nodes = [v for v in self.G.nodes() if v!=duffNode]
        nodes = [v for v in nodes if not v in nodeList]

        shortestnode = (None, None)
        
        # get the contralaterally closest node if desired
        pos = [v for v in self.G.node[duffNode]['xyz']]
        if contra:
            if pos[0] < midline:
                pos[0] = midline + (midline - pos[0])
            else:
                pos[0] = midline + (pos[0] - midline)
        pos = tuple(pos)
            
        for node in nodes:
            try:
                distance = np.linalg.norm(np.array(pos - np.array(self.G.node[node]['xyz'])))
            except:
                print "Finding the spatially nearest node requires x,y,z values"
                
            if shortestnode[0]:
                if distance < shortestnode[1]:
                    if connected:
                        if self.G.degree(node) > 0:
                            shortestnode = (node, distance)
                    else:
                        shortestnode = (node, distance)
            
            else:
                if connected:
                    if self.G.degree(node) > 0:
                        shortestnode = (node, distance)
                else:
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
            d = np.sqrt((l[1]-xyz0[0])**2 + (l[2]-xyz0[1])**2 + (l[3]-xyz0[2]**2))
            dists = dists + [(d, l[0])]
        
        print('presort')
        print(dists)
        # sort distances
        dtype = [('d', float), ('ind', int)]
        dists = np.array(dists, dtype = dtype)
        dists = np.sort(dists, order = ['d', 'ind'])
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

                
    #!! new function from master
    def weightToDistance(self):
        #!! docstring added
        ''' convert weights to a positive distance '''
        edgeList = [self.G.edge[v[0]][v[1]]['weight'] for v in self.G.edges() ]
        
        # get the maximum edge value, plus any negative correction required
        # and a small correction to keep the values above zero
        # the correction is the inverse of the number of nodes - designed to keep
        # calculations of efficiency sensible
        eMax = np.max(edgeList) + 1/float(len(self.G.nodes()))
        
        for edge in self.G.edges():
                self.G.edge[edge[0]][edge[1]]["distance"] = eMax - self.G.edge[edge[0]][edge[1]]["weight"] # convert weights to a positive distance
                
        
    ### hubs

    def hubIdentifier(self, weighted=False, assign=False, recalc=True):
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
        #!! got an error because reCalc wasn't a string, should it be??
        if recalc or not 'hubscore' in self.G.node[0].keys():
            if weighted:
                self.betweenessCentrality = np.array((centrality.betweenness_centrality(self.G, weight='distance').values()))
                
                self.weightToDistance()
                self.closenessCentrality = np.array((centrality.closeness_centrality(self.G, distance="distance").values()))
                self.degrees = np.array((nx.degree(self.G, weight='weight').values()))
                          
            else:
                self.betweenessCentrality = np.array((centrality.betweenness_centrality(self.G).values()))
                self.closenessCentrality = np.array((centrality.closeness_centrality(self.G).values()))
                self.degrees = np.array((nx.degree(self.G).values()))
            
            # normalise values to mean 0, sd 1
            self.betweenessCentrality -= np.mean(self.betweenessCentrality)
            self.closenessCentrality -=  np.mean(self.closenessCentrality)        
            self.degrees -= np.mean(self.degrees)
            
            self.betweenessCentrality /= np.std(self.betweenessCentrality)
            self.closenessCentrality /=  np.std(self.closenessCentrality)        
            self.degrees /= np.std(self.degrees)

        hubScores = map(self.hubHelper, range(len(self.G.nodes())))
        
        for n,node in enumerate(self.G.nodes()):
                self.G.node[node]['hubscore'] = hubScores[n]
        else:
            hubScores = [ self.G.node[v]['hubscore'] for v in self.G.nodes() ]
            
        # find 2 standard deviations above mean hub score
        #!! sdT not defined
        upperLimit = np.mean(np.array(hubScores)) + sdT*np.std(np.array(hubScores))
    
        # identify nodes as hubs if 2 standard deviations above hub score
        self.hubs = [n for n in self.G.nodes() if self.G.node[n]['hubscore'] > upperLimit ]

                
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
       
       
    def hubHelper(self, node):
        #!! docstring missing
        hubscore = self.betweenessCentrality[node] + self.closenessCentrality[node] + self.degrees[node]
        return(hubscore)
    
    def copyHemisphere(self, hSphere="R", midline=44.5):
        """
        This copies all the nodes and attributes from one hemisphere to the other, deleting any pre-existing
        data on the contralateral side. Particularly useful when you only have data from a single 
        hemisphere available.        
        """
        if hSphere=="L":
            for node in self.G.nodes():
                if self.G.node[node]['xyz'][0] < midline:
                    self.G.add_node(str(node)+"R")
                    self.G.node[str(node)+"R"] = {v:w for v,w in self.G.node[node].iteritems()}
                    pos = self.G.node[node]['xyz']
                    pos = (midline + (midline - pos[0]), pos[1], pos[2])
                    self.G.node[str(node)+"R"]['xyz'] = pos
                else:
                    self.G.remove_node(node)
                    
        elif hSphere=="R":
            for node in self.G.nodes():
                if self.G.node[node]['xyz'][0] > midline:
                    self.G.add_node(str(node)+"L")
                    self.G.node[str(node)+"L"] = {v:w for v,w in self.G.node[node].iteritems()}
                    pos = self.G.node[node]['xyz']
                    self.G.node[str(node)+"L"]['xyz'] = (midline - (pos[0] - midline), pos[1], pos[2])
                else:
                    self.G.remove_node(node)


    ### other functions        

    def modularity(self, hierarchy=False, diagVal=0., nodesToExclude=None):
        '''
        Modularity function borrowed (after asking nicely!) from
        https://sites.google.com/site/bctnet/measures/list and converted from 
        matlab to python code.
        
        The main modification is to allow NA values in the association matrix.
        The code is being integrated in to maybrain: http://code.google.com/p/maybrain/
        
        The function only returns a hierarchical dictionary of matrices and
        modularities if hierarchy is True. Otherwise, labels are added to
        individual nodes and the modularity is assigned as 'Q', eg brain.Q
        '''

        
        W = self.adjMat.copy()
        n0 = len(W)                                # number of nodes
        
        W = np.ma.array(W, mask=False)    # convert to masked array
        W.mask = W.data
        W.mask = False
        W[np.isnan(self.adjMat)] = 0.
        
        h=0                                     # hierarchy index
        Ci = { h:np.ma.array(np.zeros(n0),mask=False, dtype=int) } # create dictionary of hierarchy assignments and blank arrays
        if nodesToExclude:
            Ci[h].mask = Ci[h].data
            Ci[h].mask = False
            for i in [int(v) for v in nodesToExclude]:
                Ci[h].mask[i] = True
                W.mask[i,:] = True
                W.mask[:,i] = True
        
        
        # change diagonals to d only in non-masked rows/columns and assign
        # initial values
        count = 0
        for i in range(n0):
            if np.ma.is_masked(Ci[h][i]):
                pass
            else:
                Ci[h][i] = int(count)
                count+=1
                W[i,i] = diagVal
                W.mask[i,i] = False
        
        Q = { h:-1 }
        
        # get rid of nan's
        W = W[np.invert(W.mask)]
        W.shape = np.repeat(np.sqrt(len(W)),2)
        n = len(W)

        s = np.sum(W)                           # weight of edges
    
        while 1:
            K = np.sum(W, axis=1)                   # node degree
            Km = K.copy()                            # module degree
            Knm = W.copy()                          # node-to-module degree
            
            M = np.array([v for v in range(n)])     # initial module assignments
            
            Nm = np.ones(n)                         # number of nodes in modules
    
            flag=True                               # flag for within network hierarchy search
            
            while flag:
                flag=False            
                nList = [v for v in range(n)]
                random.shuffle(nList)
                while nList:
                    i = nList.pop()
                    dQ = (Knm[i,:] - Knm[i,M[i]] + W[i,i]) - K[i] * (Km-Km[M[i]]+K[i]) /s  # algorithm condition
    #            dQ=(Knm(i,:)-Knm(i,M(i))+W(i,i)) - K(i).*(Km-Km(M(i))+K(i))/s;
    
                    dQ[M[i]]=0
                    
                    max_dQ = np.max(dQ)                # find maximal increase in modularity
                    
                    if max_dQ>0:                        # if maximal increase is positive
                        j = np.argmax(dQ)
                                            
                        Knm[:,j] = Knm[:,j]+W[:,i]      # change node-to-module degrees
                        Knm[:,M[i]] = Knm[:,M[i]]-W[:,i]
                        
                        Km[j] = Km[j]+K[i]
                        Km[M[i]] = Km[M[i]] - K[i]       # change module degrees
                        
                        Nm[j] += 1                      # change number of nodes in modules
                        Nm[M[i]] -= 1
                        
                        M[i]=j                          # reassign module
                        
                        flag=True
    
    
            x, M1 = np.unique(M, return_inverse=True)

            h+=1
            Ci[h] = np.ma.array(np.zeros(n0), dtype=int)
            
            for i in range(n):
                Ci[h][Ci[h-1]==i] = int(M[i])
            Ci[h].mask=Ci[0].mask.copy()
            
            n = len(x)                                 # new number of modules
            
            W1 = np.zeros((n,n))                            # new weighted matrix
    
            for i in range(n):
                for j in range(i,n):                          # pool weights of nodes in same module w=sum(sum(W(M1==i,M1==j)));
                    A = np.zeros(W.shape)
                    indRow = np.array([z for z,v in enumerate(M1) if v==i])
                    indCol = np.array([z for z,v in enumerate(M1) if v==j])
                    
                    for x in indRow:
                        for y in indCol:
                            A[x,y] = W[x,y]
    
                    w = np.sum(A)
    #                print w
    
                    W1[i,j] = w
                    W1[j,i] = w
                    
            W = W1.copy()
            del(W1)      
            
            Q[h] = np.sum(np.diagonal(W))/s - np.sum(np.sum(W/s, axis=0)**2)     # compute modularity
            if Q[h] <= Q[h-1]:                     # if modularity does not increase
                break
            
        for node in self.G.nodes():
            self.G.node[node]['module'] = Ci[h-1][node]
        
        self.Q = Q[h-1]
        
        # return hierarchy only if desired
        if hierarchy:
            return(Ci, Q)            
            
    #!! clusters function removed
            
    def thresholdToPercentage(self, threshold):
        '''
        Functional to convert a threshold to a percentage connectivity.
        
        As this is returns a ratio between nodes and edges, it doesn't matter
        whether the graph is directed (ie an asymmetric association matrix)
        '''
        lenNodes = float(len(self.G.nodes()))
        maxEdges = float(lenNodes) * (lenNodes-1)
        
        lenEdges = len(self.adjMat.flatten()[self.adjMat.flatten()>threshold])

        pc = lenEdges / maxEdges
        return(pc)

        
    def percentConnected(self):
        '''
        This returns the percentage of the total number of possible connections
        where an edge actually exists.
        '''
        if self.directed:
            totalConnections = len(self.G.nodes()*(len(self.G.nodes())-1))
        else:
            totalConnections = len(self.G.nodes()*(len(self.G.nodes())-1)) / 2
        self.percentConnections = float(len(self.G.edges()))/float(totalConnections)
        
        return self.percentConnections
        
        
    def largestConnComp(self):
        self.bigconnG = components.connected.connected_component_subgraphs(self.G)[0]  # identify largest connected component

                
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
        
    ### brain connectivity toolbox
    def makebctmat(self):
        """
        Create a matrix for use with brain connectivity toolbox measures.
        See https://pypi.python.org/pypi/bctpy
        
        Note that missing nodes are not included, so the matrix order in
        the resulting matrix may not match the node number in the maybrain
        networkx object
        """
        self.bctmat = np.zeros((len(self.G.nodes()),len(self.G.nodes())))
        nodeIndices = dict(zip(self.G.nodes(), range(len(self.G.nodes()))))
        for nn,x in enumerate(self.G.nodes()):
            for y in self.G.edge[x].keys():
                try:
                    self.bctmat[nn,nodeIndices[y]] = self.G.edge[x][y]['weight']
                except:
                    pass
        
    def assignbctResult(self, bctRes):
        ''' translate a maybrain connectome into a bct compatible format ''' 
        out = dict(zip(self.G.nodes(), bctRes))
        return(out) 

    ##### miscellaneous functions


    
    def strnum(num, length=5):
        ''' convert a number into a string of a given length'''
        
        sn = str(num)
        lenzeros = length - len(sn)
        sn = lenzeros*'0' + sn
        
        return sn        


class highlightObj():
    ''' object to hold information to highlight a subsection of a brainObj '''
    
    def __init__(self, nodes = [], edges = []):
        ''' Points refer to the indices of node coordinates in the brain object to which it
        is related. Edges is a set of pairs of coordinates of the same brain object '''
        
#        self._mode = 'pe' # p, e or pe for points only, edges only or points and edges
        self.nodeIndices = nodes # indices of points used from a brain object
        self.edgeIndices = edges # list of ordered pairs of edges
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
                
        for e in self.edgeIndices:
            # get coordinates of edge indices
            p1 = brain.G.node[e[0]]['xyz']
            p2 = brain.G.node[e[1]]['xyz']    
            
            x1.append(p1[0])
            x2.append(p2[0]-p1[0]) 
            
            y1.append(p1[1])
            y2.append(p2[1]-p1[1])
            
            z1.append(p1[2])
            z2.append(p2[2]-p1[2])
            
            s.append(brain.G.edge[e[0]][e[1]]['weight'])
            
        return x1, y1, z1, x2, y2, z2, s
        
    
    def getCoords(self, brain):
        ''' turn indices of points into coordinates '''

        x = []
        y = []
        z = []        
        
        for p in self.nodeIndices:
            pc = brain.G.node[p]['xyz']
            # pc = brain.G.node[p[1]]['xyz']
            x.append(pc[0])
            y.append(pc[1])
            z.append(pc[2])
            
        return x, y, z
        
    
    def hlPrint(self, brain):
        ''' print the contents of the highlight object '''
        
        # introduction
        print('# ================ #')
        print('\n'+'highlight information:')
        # colour
        print('colour', self.colour)
        # opacity
        print('opacity', self.opacity) 
        # edge opacity
        print('edge opacity', self.edgeOpacity)
        
        # nodes
        print('\n nodes')
        for n in self.nodeIndices:
            print(str(n), brain.G.node[n]['xyz'])
            
        # edges
        print('\n edges')
        for e in self.edgeIndices:
            print(str(e), brain.G.node[e[0]]['xyz'], brain.G.node[e[1]]['xyz'])
            
        print('\n # ================ # \n')
        
       