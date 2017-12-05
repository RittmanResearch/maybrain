# -*- coding: utf-8 -*-
"""
Created on Fri Mar 16 18:10:07 2012

@authors: Timothy Rittman, Martyn Rittman

Documentation available at https://github.com/RittmanResearch/maybrain

"""
import networkx as nx
import numpy as np
from networkx.algorithms import components
import random
from maybrain import extraFns
from maybrain import highlightObj


class brainObj:
    """
    A class that defines a brain network created from an adjacency matrix and spatial information 
    with certain properties. The added extras are:
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
        self.subject = None # identification of the subject to which this brain object belongs
        self.scan = None # information about the scan which generated this brain object
     
        # TODO: need to define the following, what do they do???
        self.dyingEdges = {}
        
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
        
        self.riskEdges = None
        
        # Properties in edges/nodes
        self.WEIGHT       = 'weight'
        self.LINKED_NODES = 'linkedNodes'
        self.XYZ          = 'xyz'
        self.ANAT_LABEL   = 'anatlabel'
        self.DISTANCE     = 'distance'

        # create a new networkX graph object
        if self.directed:
            self.G = nx.DiGraph()
        else:
            self.G = nx.Graph()

    ## ================================================
    ##### File inputs and outputs

    ### edges and nodes        
    def importAdjFile(self, fname, delimiter = None, exclnodes=[], naVals=["NA"]):
        '''
        Imports an adjacency matrix from a file.
        fname : file name
        delimiter : the delimiter of the values inside the matrix, like ","
        exclnodes : Nodes you don't want to load, in an array format (nodes no. starts from zero)
        naVals : How the "Not a Number" values are represented in the file
        '''
        self.exclnodes=exclnodes
        
        lines = []
        try:
            with open(fname, "r") as f:
                for l in f:
                    while l[-1] in ('\n', '\t'):
                        l = l[:-1]
                    lines.append(list(map(float, [v if v not in naVals else np.nan for v in l.split(sep=delimiter)])))
        except IOError as error:            
            error.strerror = 'Problem with opening file "' + fname + '": ' + error.strerror
            raise error

        # set adjacency matrix
        self.adjMat = np.array(lines)
        
        # add nodes
        self.G.add_nodes_from([v for v in range(len(lines)) if not v in exclnodes])

        # update adjacency matrix to null values of excluded nodes
        if exclnodes:
            for en in exclnodes:
                self.adjMat[:,en]=np.nan
                self.adjMat[en,:]=np.nan


    def importSpatialInfo(self, fname, delimiter=None, convertMNI=False):
        ''' 
        Add 3D coordinate information for each node from a given file. It needs to be called after importAdjFile()
        fname : file name
        delimiter : the delimiter of the values inside the matrix, like ","
        convertMNI : Whether you want to convert coordinates from voxel-wise to MNI space
        '''
        # open file
        try:
            f = open(fname,"r")
        except IOError as error:    
            error.strerror = 'Problem with opening 3D position file "' + fname + '": ' + error.strerror
            raise error   
               
        # get data from file
        lines = f.readlines()
        nodeCount=0
        for line in lines:
            l = line.split(sep=delimiter)

            if convertMNI:
                l[1] = 45 - (float(l[1])/2)
                l[2] = 63 + (float(l[2])/2)
                l[3] = 36 + (float(l[3])/2)

            if nodeCount in self.G.nodes():
                self.G.node[nodeCount][self.XYZ]=(float(l[1]),float(l[2]),float(l[3]))
                self.G.node[nodeCount][self.ANAT_LABEL]=l[0]

            nodeCount+=1

        f.close()
#TODO: When applying threshold, apply properties again.
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
                value = l.split()
                if len(value)==2:
                    mode = 'nodes'
                elif len(value)==3:
                    mode =  'edges'
            except:
                print(('couldn\'t parse data line. Please check property file', l))
                
            # make last entry of value a float or string (property value)
            try:
                value[-1] = float(value[-1])
            except ValueError:
                value[-1] = str(value[-1])
                value[-1] = extraFns.stripString(value[-1])
                print((value[-1]))
            
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

        return prop # required for GUI
                    
        
    def addNodeProperties(self, propertyName, nodeList, propList):
        ''' add properties to nodes, reading from a list of nodes and a list of 
            corresponding properties '''
        
        for ind in range(len(nodeList)):
            n = nodeList[ind]
            p = propList[ind]
            try:
                self.G.node[n][propertyName] = p
            except:
                print(('property assignment failed: ' + propertyName + ' ' + str(n) + ' ' + str(p)))


    def addEdgeProperty(self, propertyName, edgeList, propList):
        ''' add a property to a selection of edges '''
        
        for ind in range(len(edgeList)):
            e = edgeList[ind]
            p = propList[ind]
            try:
                self.G.edge[e[0]][e[1]][propertyName] = p
            except KeyError:
                print(('edge property assignment failed: ' + propertyName + ' ' + str(e) + ' ' + str(p)))

        
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
        ''' 
        Imports an isosurface info using nibabel
        gives a 3D array with data range 0 to 255 for test data
        defines an nibabel object, plus ndarrays with data and header info in
        
        fname: File Name with the isosurface
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
            # parcel files start from 1, zero is for background
            nArr = np.ma.masked_where(self.iso != (n+1), self.iso)
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
            for n in list(valueDict.keys()):
                outMat[np.where(self.nbiso.get_data()==n+1)] = valueDict[n]
        else:
            outMat = self.parcelList
        
        N = nb.Nifti1Image(outMat, self.nbiso.get_affine(), header=self.isoHeader)

        nb.save(N, outname+'.nii')
                     
        

    ## ===========================================================================    
    
    
    ##### Functions to alter the brain

    def applyThreshold(self, thresholdType = None, value = 0., useAbsolute=False):
        ''' 
        Treshold the adjacency matrix to determine which nodes are linked by edges.
        
        thresholdType : The type of threshold applied. Four options are available:
            "edgePC" -> retain the most connected edges as a percentage of the total possible number of edges. "value" must be between 0 and 100
            "totalEdges" -> retain the most strongly connected edges
            "tVal" -> retain edges with a weight greater or equal than value
            None -> all possible edges are created
        value : Value according to thresholdType
        useAbsolute : Thresholding by absolute value. For example, if this is set to False, a \
            weight of 1 is stronger than -1. If this is set to True, these values are equally \
            strong. This affects thresholding with "edgePC", "totalEdges" and "tVal". In case of \
            "tVal", it will threshold for weights >= abs(tVal) and <= -abs(tVal)
        '''
        
        # Controlling input
        if thresholdType not in ["edgePC", "totalEdges", "tVal", None]:
            raise TypeError("Not a valid thresholdType for applyThreshold()")
        if thresholdType == "edgePC" and (value < 0 or value > 100):
            raise TypeError("Invalid value for edgePC in applyThreshold()")
        

        # Creating the array with weights and edges
        upperValues = np.triu_indices(np.shape(self.adjMat)[0], k= 1)
        weights = []
        
        # Creating weights in which each element is (node1, node2, weight)
        for x in np.nditer(upperValues):
            if not np.isnan(self.adjMat[x[0]][x[1]]):
                weights.append((int(x[0]), 
                                int(x[1]),
                                self.adjMat[x[0]][x[1]]))
    
        # If directed, also add the lower down part of the adjacency matrix
        if self.directed:
            belowValues = np.tril_indices(np.shape(self.adjMat)[0], k= -1)
            for x in np.nditer(belowValues):
                if not np.isnan(self.adjMat[x[0]][x[1]]):
                    weights.append((int(x[0]),
                                    int(x[1]),
                                    self.adjMat[x[0]][x[1]]))

        # Filtering weights when with thresholdType of "edgePC" or "totalEdges"
        if thresholdType in ["edgePC", "totalEdges"]:
            # Sorting ascending
            if useAbsolute:
                weights.sort(key=lambda x: abs(x[2]))
            else:
                weights.sort(key=lambda x: x[2])
            
            # Getting the number of edges to include
            if thresholdType == 'edgePC':
                edgeNum = int( (value/100.) * len(weights))   
            else: #totalEdges
                edgeNum = int(value)
            
            # Removing weak edges
            if edgeNum <= 0:
                weights = []
            elif edgeNum > len(weights):
                pass #include all weights
            else:
                weights = weights[-edgeNum:]
            
        # remove previous edges
        self.G.remove_edges_from(self.G.edges())   
        
        # Adding the edges
        for e in weights:
            if thresholdType == 'tVal' and useAbsolute:
                if e[2] >= abs(value) or e[2] <= -abs(value):
                    self.G.add_edge(e[0], e[1], weight = e[2])
            elif thresholdType == 'tVal' and not useAbsolute:
                if e[2] >= value:
                    self.G.add_edge(e[0], e[1], weight = e[2])
            else: # None, edgePC, totalEdges
                self.G.add_edge(e[0], e[1], weight = e[2])
                        

    def reconstructAdjMat(self):
        '''
        It redefines the adjacency matrix from the edges' weights of G
        It assumes that size of adjMat is maintained
        '''
        self.adjMat[:] = np.nan        
        
        for e in self.G.edges():
            self.updateAdjMat(e)
                
    def updateAdjMat(self, edge):
        ''' 
        It updates the adjacency matrix by bringing the weight of an edge in G \
        to the adjacency matrix
        
        edge : The edge in G to bring to adjMat'''
        
        try:
            w = self.G.edge[edge[0]][edge[1]][self.WEIGHT]
            self.adjMat[edge[0], edge[1]] = w
            
            if not self.directed:
                self.adjMat[edge[1], edge[0]] = w
        except KeyError as error:
            import sys
            _, _, tb = sys.exc_info()
            raise KeyError(error, "Edge does not exist in G or doesn't have WEIGHT property").with_traceback(tb)
        except IndexError as error:
            import sys
            _, _, tb = sys.exc_info()
            raise IndexError("adjMat too small to have such an edge").with_traceback(tb).with_traceback(tb)

    def localThresholding(self, thresholdType=None, value=0.):
        '''
        Threshold the adjacency matrix by building from the minimum spanning tree (MST) and adding successive N-nearest neighbour degree graphs.
        Thus, if you want to have a local thresholding of N edges when the MST has more than N edges, thresholding will retain the MST
        It only works for undirected graphs
        
        thresholdType : The type of threshold applied. Three options are available:
            "edgePC" -> retain the most connected edges as a percentage of the total possible number of edges. "value" must be between 0 and 100
            "totalEdges" -> retain the most strongly connected edges
            None -> retain the minimum spanning tree
        value : Value according to thresholdType
        '''
        
        # Controlling input
        if thresholdType not in ["edgePC", "totalEdges", None]:
            raise TypeError("Not a valid thresholdType for localThresholding()")
        if self.directed:
            raise TypeError("localThresholding() not available for directed graphs")
        if thresholdType == "edgePC" and (value < 0 or value > 100):
            raise TypeError("Invalid value for edgePC for localThresholding()")
        
        
        # Putting all the edges in the G object for local thresholding
        self.applyThreshold()
        
        if not nx.is_connected(self.G):
            raise TypeError("Adjacency Matrix is not connected. Impossible to execute localThresholding()")
            
        # create minimum spanning tree
        T = nx.minimum_spanning_tree(self.G)

        if not(thresholdType):
            self.G = T
            return #Nothing else to do, just return
        elif thresholdType == 'edgePC':
            # find threshold as a percentage of total possible edges
            upperValues = np.triu_indices(np.shape(self.adjMat)[0], k= 1)
            # only flatten the upper right part of the matrix
            weights = np.array(self.adjMat[upperValues])
            
            # remove NaNs
            weights = weights[~np.isnan(weights)]
            # calculate percentage
            edgeNum = int(value/100. * len(weights))
        else: # 'totalEdges' option
            edgeNum = value

        
        lenEdges = len(T.edges())
        if lenEdges > edgeNum:
            print("Warning: The minimum spanning tree already has: "+ str(lenEdges) + " edges, select more edges.",
                  "Local Threshold will be applied by just retaining the Minimum Spanning Tree")
            self.localThresholding()
            return
 
        k=1 # number of degrees for NNG
        while lenEdges < edgeNum:
            # create nearest neighbour graph
            nng = self._NNG(k)
            
            # remove edges from the NNG that exist already in the new graph/MST
            nng.remove_edges_from(T.edges())
            
            # Ending condition. No more edges to add so break the cycle
            if len(nng.edges())==0:
                break
            
            # add weights to NNG
            for e in nng.edges():
                nng.edge[e[0]][e[1]][self.WEIGHT] = self.adjMat[e[0],e[1]]
            
            # get a list of edges from the NNG in order of weight
            edgeList = sorted(nng.edges(data=True), key=lambda t: t[2][self.WEIGHT], reverse=True)
            
            # add edges to graph in order of connectivity strength
            for edge in edgeList:
                T.add_edges_from([edge])
                lenEdges = len(T.edges())
                if lenEdges >= edgeNum:
                    break
            
            k+=1
        
        self.G = T
        
    def binarise(self):
        '''
        Removes weighting from edges by assigning a weight of 1 to the existing edges
        '''
        for edge in self.G.edges(data=True):
            edge[2][self.WEIGHT] = 1        

    def removeUnconnectedNodes(self):
        '''
        Removes nodes with no connections
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
        if not(mode in ['edge', 'node', 'node or edge']):
            print('filter mode not recognised')
            return
        # check if label given
        if not label:
            label = self._getAutoLabel()
            
        # make a highlight object
        h = highlightObj()
        h.colour = colour
        h.opacity = opacity
        h.edgeIndices = []
        h.nodeIndices = []
        
        print((prop, rel, val, label, mode))
        
        # extract lists from edges  
        if mode in ['edge', 'node or edge']:
            ind = -1
            for e in self.G.edges(data = True):
                ind = ind +1
                print((self.G.edge[e[0]][e[1]]))
                try:
                    d = self.G.edge[e[0]][e[1]][prop]
                except KeyError:
                    continue           
                
                # match properties
                boolval = self.propCompare(d, rel, val)
                print((d, rel, val, boolval))
                
                # save data in highlight
                if boolval:
                    h.edgeIndices.append((e[0],e[1])) 
                    
        
        # extract lists from nodes
        if mode in ['node', 'node or edge']:
            for c in range(len(self.G.nodes())):
                # get property

                # special treatment for 'x', 'y' and 'z'
                if prop=='x':
                    d = self.G.node[c][self.XYZ][0]
                elif prop=='y':
                    d = self.G.node[c][self.XYZ][1]
                elif prop=='z':
                    d = self.G.node[c][self.XYZ][2]
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
            print(('relation not recognised: ' + rel ))
        
        return b        


    def _getAutoLabel(self):
        ''' generate an automatic label for a highlight object if none given '''
        
        # get index of label
        num = str(self.labelNo)
        num = '0' * (4-len(num)) + num
        
        # make label and print
        label = 'plot ' + num
        print(('automatically generated label: '+ label))
        
        # iterate label index
        self.labelNo = self.labelNo + 1
        
        return label     

            
    def degenerate(self, weightloss=0.1, edgesRemovedLimit=1, threshLimit=None,
                   pcLimit=None, weightLossLimit=None, nodeList=[],
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
        
        # get rid of the existing list of edges if the node list is specified
        if nodeList:
            self.riskEdges = None
            
        # set limit
        if weightLossLimit and pcLimit:
            print("You have asked for both a weight and percentage connectivity limit, using the percentage connectivity limit")
        
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
                print("The percentage threshold set is lower than the current graph, please choose a larger value")
            
            limit = lenEdges - newEdgeNum
            weightLossLimit = False
            
        elif weightLossLimit:
            limit = weightLossLimit
        
        else:
            limit = edgesRemovedLimit
        
        if not self.riskEdges:
            reDefineEdges=True
            # if no toxic nodes defined, select the whole graph
            if not nodeList:
                nodeList = self.G.nodes()
            
            # generate list of at risk edges
            self.riskEdges = [v for v in nx.edges(self.G, nodeList) if self.G.edge[v[0]][v[1]][self.WEIGHT] != 0.]
        else:
            reDefineEdges=False
            
        if spread:
            nodeList = []
        
        # check if there are enough weights left
        riskEdgeWtSum = np.sum([self.G.edge[v[0]][v[1]][self.WEIGHT] for v in self.riskEdges])
        if limit > riskEdgeWtSum:
            print("Not enough weight left to remove")
            return nodeList
            
        while limit>0.:
            if not self.riskEdges and spatialSearch:
                # find spatially closest nodes if no edges exist
                # is it necessary to do this for all nodes?? - waste of computing power,
                # choose node first, then calculated spatially nearest of a single node
                newNode = self.findSpatiallyNearest(nodeList)
                if newNode:
                    print("Found spatially nearest node")
                    nodeList.append(newNode)
                    self.riskEdges = nx.edges(self.G, nodeList)
                else:
                    print("No further edges to degenerate")
                    break
            # choose at risk edge to degenerate from           
            dyingEdge = random.choice(self.riskEdges)            
                        
            # remove specified weight from edge
            w = self.G[dyingEdge[0]][dyingEdge[1]][self.WEIGHT]
            
            if np.absolute(w) < weightloss:
                loss = np.absolute(w)
                self.G.remove_edge(dyingEdge[0], dyingEdge[1])
                self.riskEdges.remove(dyingEdge)
                if not weightLossLimit:
                    limit-=1
            
            elif w>0:
                loss = weightloss
                self.G[dyingEdge[0]][dyingEdge[1]][self.WEIGHT] -= weightloss
                
            else:
                loss = weightloss
                self.G[dyingEdge[0]][dyingEdge[1]][self.WEIGHT] += weightloss
            
            # record the edge length of edges lost
            if distances:
                self.dyingEdges[dyingEdge] = self.G[dyingEdge[0]][dyingEdge[1]]
                self.dyingEdges[dyingEdge][self.DISTANCE] =  np.linalg.norm( np.array((self.G.node[dyingEdge[0]][self.XYZ])) - np.array((self.G.node[dyingEdge[1]][self.XYZ]))  )
            
            # update the adjacency matrix (essential if robustness is to be calculated)            
            if updateAdjmat:
                self.updateAdjMat(dyingEdge)
                            
            # add nodes to toxic list if the spread option is selected
            if spread:
                for node in dyingEdge:
                    if not node in nodeList and self.G.edge[dyingEdge[0]][dyingEdge[1]] > spread:
                        nodeList.append(node)
                        
            if weightLossLimit:
                limit -= loss
                
            # redefine at risk edges
            if reDefineEdges or spread:
                self.riskEdges = nx.edges(self.G, nodeList)        
        
        ## Update adjacency matrix to reflect changes
        #self.reconstructAdjMat()
        
        print("Number of toxic nodes: "+str(len(nodeList)))
        
        return nodeList

    def contiguousSpread(self, edgeloss, startNodes = None):
        ''' degenerate nodes in a continuous fashion. Doesn't currently include spreadratio '''

        # make sure nodes have the linkedNodes attribute
        try:
            self.G.node[0][self.LINKED_NODES]
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
            l = self.G.node[t][self.LINKED_NODES]
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
            print(('deadNode', deadNode))
            
            
            # add the new at-risk nodes
            l = self.G.node[deadNode][self.LINKED_NODES]
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
    
   
    # =============================================================
    
    ##### Analysis functions

    def _NNG(self, k):
        ''' Private method to help local thresholding by creating a k-nearest neighbour graph'''
        G = nx.Graph()
        nodes = list(range(len(self.adjMat[0])))
        
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
        pos = [v for v in self.G.node[duffNode][self.XYZ]]
        if contra:
            if pos[0] < midline:
                pos[0] = midline + (midline - pos[0])
            else:
                pos[0] = midline + (pos[0] - midline)
        pos = tuple(pos)
            
        for node in nodes:
            try:
                distance = np.linalg.norm(np.array(pos - np.array(self.G.node[node][self.XYZ])))
            except:
                print("Finding the spatially nearest node requires x,y,z values")
                
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
            xyzList.append([count] + list(self.G.node[node][self.XYZ]))
            count = count  + 1

        # cut down in x,y and z coords
        xyz0 = self.G.node[randNode][self.XYZ]
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
        print('newxyzlist')
        print(newxyzList)
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
        ''' 
        It gives to each node a list containing the linked nodes.
        If Graph is undirected, and there is an edge (1,2), node 1 is linked \
         to node 2, and vice-versa. If Graph is directed, thus node 1 is linked \
         to node 2, but not the other way around
        This property can be accessed through self.LINKED_NODES
        Be sure to call this method again if you threshold your brainObj again
        '''
    
        # Resetting all nodes from some past information (few edges might not \
        #  be able to reset this field in all nodes)
        for n in self.G.nodes(data=True):
            n[1][self.LINKED_NODES] = []
        
        for l in self.G.edges():            
            # add to list of connecting nodes for each participating node
            self.G.node[l[0]][self.LINKED_NODES].append(l[1])
            
            if not self.directed:
                self.G.node[l[1]][self.LINKED_NODES].append(l[0])

    def weightToDistance(self):
        '''
        It inverts all the edges' weights so they become equivalent to a distance measure. \ 
        With a weight, the higher the value the stronger the connection. With a distance, \ 
        the higher the value the "weaker" the connection. 
        In this case there is no measurement unit for the distance, as it is just \
        a conversion from the weights.
        The distances can be accessed in each node's property with self.DISTANCE
        '''
        edgeList = [v[2][self.WEIGHT] for v in self.G.edges(data=True) ]
        
        # get the maximum edge value, plus a small correction to keep the values above zero
        # the correction is the inverse of the number of nodes - designed to keep
        # calculations of efficiency sensible
        eMax = np.max(edgeList) + 1/float(self.G.number_of_nodes())
        
        for edge in self.G.edges():
            self.G.edge[edge[0]][edge[1]][self.DISTANCE] = eMax - self.G.edge[edge[0]][edge[1]][self.WEIGHT] # convert weights to a positive distance
                
    def copyHemisphere(self, hSphere="R", midline=44.5):
        """
        This copies all the nodes and attributes from one hemisphere to the other, deleting any pre-existing
        data on the contralateral side. Particularly useful when you only have data from a single 
        hemisphere available.        
        """
        if hSphere=="L":
            for node in self.G.nodes():
                if self.G.node[node][self.XYZ][0] < midline:
                    self.G.add_node(str(node)+"R")
                    self.G.node[str(node)+"R"] = {v:w for v,w in self.G.node[node].items()}
                    pos = self.G.node[node][self.XYZ]
                    pos = (midline + (midline - pos[0]), pos[1], pos[2])
                    self.G.node[str(node)+"R"][self.XYZ] = pos
                else:
                    self.G.remove_node(node)
                    
        elif hSphere=="R":
            for node in self.G.nodes():
                if self.G.node[node][self.XYZ][0] > midline:
                    self.G.add_node(str(node)+"L")
                    self.G.node[str(node)+"L"] = {v:w for v,w in self.G.node[node].items()}
                    pos = self.G.node[node][self.XYZ]
                    self.G.node[str(node)+"L"][self.XYZ] = (midline - (pos[0] - midline), pos[1], pos[2])
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
            
            
    def thresholdToPercentage(self, threshold):
        '''
        It returns a ratio between the edges on adjMat above a certain threshold value \
        and the total possible edges of adjMat.
        In an unidrected graph, the total possible edges are the upper right \
        part elements of adjMat different from np.nan. In a directed graph, the total \
        possible edges are all the elements of adjMat except the diagonal and \
        np.nan
        
        threshold : The threshold value
        '''
        upperValues = np.triu_indices(np.shape(self.adjMat)[0], k= 1)
        belowValues = np.tril_indices(np.shape(self.adjMat)[0], k= -1)
        if not self.directed:
            # only flatten the upper right part of the matrix
            weights = np.array(self.adjMat[upperValues])
        else:
            weights = np.concatenate((self.adjMat[upperValues], self.adjMat[belowValues]))
        # remove NaNs
        weights = weights[~np.isnan(weights)]
        
        maxEdges = len(weights)
        lenEdges = len(weights[weights>threshold])
        
        return lenEdges / maxEdges
        
    def percentConnected(self):
        '''
        This returns the ratio of the current number of edges in our `G` object \
        and the total number of possible connections.
        If N is the number of nodes in our `G` object, the total number of \
        possible connections is (N * (N - 1))/2 for an undirected graph, and \
        N * (N-1) for a directed graph.
        
        '''
        nodes = self.G.number_of_nodes()
        
        if self.directed:
            totalConnections = nodes * (nodes - 1)
        else:
            totalConnections = nodes * (nodes - 1) / 2
        return float(self.G.number_of_edges()) / float(totalConnections)
        
         
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
        nodeIndices = dict(list(zip(self.G.nodes(), list(range(len(self.G.nodes()))))))
        for nn,x in enumerate(self.G.nodes()):
            for y in list(self.G.edge[x].keys()):
                try:
                    self.bctmat[nn,nodeIndices[y]] = self.G.edge[x][y][self.WEIGHT]
                except:
                    pass
        
    def assignbctResult(self, bctRes):
        ''' translate a maybrain connectome into a bct compatible format ''' 
        out = dict(list(zip(self.G.nodes(), bctRes)))
        return(out)

