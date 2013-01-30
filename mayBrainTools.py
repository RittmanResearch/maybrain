# -*- coding: utf-8 -*-
"""
Created on Fri Mar 16 18:10:07 2012

@author: -

Change log:
    - adjacency thresholding moved into separate function
    - symmetric option for adjacency matrix
    - fixed bug: all nodes linked to themselves
    - contiguous spread functions
    - function to save adjacency matrix to file
        
24/01/2013 Modifications from Tim's files - 13-01
	- added degeneracy fn to brainObj
	- binarise and largest conn comp added
	- hub identifier can take weighted measures
        
"""

import string,os,csv,sys,itertools # ,community
import networkx as nx
import numpy as np
from networkx.drawing import *
from networkx.algorithms import centrality
from networkx.algorithms import components
from numpy import shape, random, fill_diagonal, array, where, tril, zeros
from mayavi import mlab
from string import split
import nibabel as nb
from os import path, remove

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
        		- drawpos: draw 2D position ??? What does this mean ??? removed from arguments for now
    
        '''        
        
        # create an empty graph
        self.G = nx.Graph()
        self.iter = None # not sure where this is used. It is often defined, but never actually used!
        
        # define global variables
        self.adjMat = None # adjacency matrix, containing weighting of edges.
        
        self.hubs = []
        self.lengthEdgesRemoved = None
        
    

    ## ================================================

    ## File inputs        

    def readAdjFile(self, fname, threshold = None, edgePC = None, totalEdges = None, symmetric = False):
        ''' load adjacency matrix with filename and threshold for edge definition '''
        
        print('loading data from' + fname)

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
            lines.append(map(float, split(l)))
        nodecount = len(lines)                

        # close file                
        f.close()

        # create adjacency matrix as an array
        self.adjMat = array(lines)       # array of data
        # zero below-diagonal values if symmetric
        if symmetric:
            self.adjMat = self.adjMat - tril(self.adjMat)
        sh = shape(self.adjMat)
        # check if it's diagonal
        if sh[0]!=sh[1]:
            print("Note Bene: adjacency matrix not square.")
            print(sh)

        # create nodes on graph
        self.G.add_nodes_from(range(nodecount))  # creates one node for every line in the adjacency matrix
        
        # create edges by thresholding adjacency matrix
        self.adjMatThresholding(edgePC, totalEdges, threshold)
      
      
    def readSpatialInfo(self, fname):
        ''' add 3D coordinate information for each node from a given file '''
        
        coords = {}

        # Try something slightly crazy if filename not defined
        if not fname:
            if os.path.exists("anat_labels.txt"):
                ROIs = {}
                anat = csv.reader(open('anat_labels.txt','rb'),skipinitialspace = True)
                lblsfile = csv.reader(open('ROI_labels.txt','rb'),delimiter = '\t')
                lbls = lblsfile.next()
                
                count = 1
                for line in anat:
                    ROIs[count] = [v for v in line]
                    count += 1
                
                x = string.split(open('ROI_xyz.txt').readlines()[0])
                y = string.split(open('ROI_xyz.txt').readlines()[1])
                z = string.split(open('ROI_xyz.txt').readlines()[2])
                
                for n in range(len(self.G.nodes())):
                    self.G.node[n]['anatlabel'] = ROIs[int(lbls[n])]
                    self.G.node[n]['xyz'] = (float(x[n]),float(y[n]),float(z[n]))

        else:
            # open file
            try:
                f = open(fname,"rb")
            except IOError, error:
                (errorno, errordetails) = error
                print "Couldn't find 3D position information"
                print "Problem with opening file: "+errordetails                
            
            # get data from file
            lines = f.readlines()
            for l in lines:
                line = l.split()
                if len(line) == 4:
                    # read in line
                    coords[line[0]] = {"x":float(line[1]), "y":float(line[2]), "z":float(line[3])}
                    ROIs = coords.keys()
                    ROIs.sort()
                    
                    # get x, y and z coords
                    x = [coords[ROI]["x"] for ROI in ROIs]
                    y = [coords[ROI]["y"] for ROI in ROIs]
                    z = [coords[ROI]["z"] for ROI in ROIs]
                    
            # add coords to nodes
            for n in range(len(self.G.nodes())):
                self.G.node[n]['anatlabel'] = ROIs[n]
                self.G.node[n]['xyz'] = (float(x[n]),float(y[n]),float(z[n]))                
                

    def importSkull(self, fname):
        ''' Import a file for skull info using nbbabel
            gives a 3D array with data range 0 to 255 for test data
            could be 4d??
            defines an nibabel object, plus ndarrays with data and header info in
        
        '''        
        
        self.nbskull = nb.load(fname)
        self.skull = self.nbskull.get_data()
        self.skullHeader = self.nbskull.get_header()        
                
                
                
    def inputNodeProperties(self, propertyName, nodeList, propList):
        ''' add properties to nodes, reading from a list of nodes and a list of corresponding properties '''
        
        for ind in range(len(nodeList)):
            n = nodeList[ind]
            p = propList[ind]
            self.G.node[n][propertyName] = p
            
            
    def outputAdjMatrix(self, filename, header = None):
        ''' output the adjacency matrix to file. Header is a string.'''
        
        if path.exists(filename):
            print("old adjacency file removed")
            remove(filename)
        
        # write header
        if header:
            f = open(filename, 'w+')
            f.write(header + '\n')
            f.close()                

        f = open(filename, 'a+')

        l = f.readlines()
        nl = len(l)
        
        # append line to say where data starts
        headerApp = 'data begins line '  + str(nl+2) + '\n'
        f.write(headerApp)
        
        # write data to file from adjacency array
        for line in self.adjMat:
            for num in line[:-1]:
                f.write(str(num) + '\t')
            f.write(str(line[-1]) + '\n')
            
        f.close()
            
        print("data written to " + filename)
        
        
    def outputEdges(self, filename, header = None, properties = []):
        ''' output the edges to file '''        
        
        if path.exists(filename):
            print("old edge file removed")
            remove(filename)   
            
        # open file and write header
        f = open(filename, 'w')
        if header:
            f.write(header + '\n')
            
        # write column headers
        line = 'x' + '\t' + 'y'
        for p in properties:
            line = line + '\t' + p
        line = line + '\n'
        f.write(line)
        
        for n in self.G.edges(data = True):
            # add coordinates
            line = str(n[0]) + '\t' + str(n[1]) 
            # add other properties
            for p in properties:
                try:
                    line = line + '\t' + str(n[2][p])
                except:
                    print(p, "omitted in edge output")
            line = line + '\n'
            # write out
            f.write(line)
        f.close()

        print("edges written to " + filename)            
        
    
    def outputEdgesMatrix(self, filename):
        ''' output the edge data as a boolean matrix '''
        
        if path.exists(filename):
            print("old edge matrix file removed")
            remove(filename)           
        
        n = self.G.number_of_nodes()
        mat = zeros([n,n], dtype = int)
        
        for ed in self.G.edges():
            x = ed[0]
            y = ed[1]
            if y>x:
                mat[x,y] = 1
            else:
                mat[y,x] = 1
            
        f = open(filename, 'w')
        for row in mat:
            for ch in row[:-1]:
                 f.write(str(ch) + '\t')
            f.write(str(row[-1]) + '\n')
        f.close()
                
        print("edges written to " + filename)            
    
    
    ## ===========================================================================
    
    ## Functions to alter the brain

    def adjMatThresholding(self, edgePC = None, totalEdges = None, tVal = None):
        ''' apply thresholding to the adjacency matrix. This can be done in one of
            three ways (in order of decreasing precendence):
                edgePC - this percentage of nodes will be linked
                totalEdges - this number of nodes will be linked
                tVal - give an absolute threshold value. Pairs of nodes with a corresponding adj matrix
                       value greater than this will be linked.
        '''
        
        # check if adjacency matrix is there
        if self.adjMat == None:
            print("No adjacency matrix. Please load one.")
            return
        
        # remove existing edges
        self.G.remove_edges_from(self.G.edges())
        nodecount = len(self.G)
            
        # check some method of thresholding was defined
        if (edgePC==None)&(totalEdges==None)&(tVal==None):
            print("No method of thresholding given. Please define edgePC, totalEdges or tVal.")
            return
        
        # get the number of edges to link
        if edgePC:
            # find threshold as a percentage of total possible edges
            edgeNum = int(edgePC * (nodecount * (nodecount-1)))
        elif totalEdges:
            # allow a fixed number of edges
            edgeNum = totalEdges
        else:
            edgeNum = -1
        
        # get threshold
        if edgeNum>=0:
            # get treshold value
            weights = self.adjMat.flatten()
            weights.sort()
            threshold = weights[-edgeNum]
            print("Threshold set at: "+str(threshold))
        else:
            threshold  = tVal        
            
        
        # carry out thresholding on adjacency matrix
        boolMat = self.adjMat>threshold
        fill_diagonal(boolMat, 0)
        edgeCos = where(boolMat) # lists of where edges should be
        
        for ind in range(len(edgeCos[0])):
            node1 = edgeCos[0][ind]
            node2 = edgeCos[1][ind]
            
            if not(self.G.has_edge(node1, node2)):
                self.G.add_edge(node1, node2, weight = self.adjMat[node1, node2])            
            
    
    def makeSubBrain(self, propName, value):
        ''' separate out nodes and edges with a certain property '''

        # initialise new brain
        subBrain = brainObj()
        
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
                
        return subBrain
        
        
    def makeSubBrainEdges(self, propName, value):
        ''' separate out edges with a certain property '''

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
                
                            
        return subBrain        


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
        


    def contiguousspread(self, edgeloss, largestconnectedcomp=False, startNodes = None):
        ''' remove edges in a continuous fashion. Doesn't currently include spreadratio '''

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
                if not(a in toxicNodes)&(not(self.G.node[a]['degenerating'])):
                    newl.append(a)
            riskNodes = riskNodes + newl
            
#        print(riskNodes)
        
        # iterate number of steps
        for count in range(edgeloss):
            # find at risk nodes
            ind = random.randint(0, len(riskNodes))
            deadNode = riskNodes.pop(ind) # get the index of the node to be removed and remove from list
            # remove all instances from list
            while deadNode in riskNodes:
                riskNodes.remove(deadNode)
            
            # add to toxic list    
            toxicNodes.append(deadNode)
            # make it degenerate
            self.G.node[deadNode]['degenerating'] = 'True'
            
            
            # add the new at-risk nodes
            l = self.G.node[deadNode]['linkedNodes']
            newl = []
            # check the new indices aren't already toxic
            for a in l:
                if not(a in toxicNodes)&(not(self.G.node[a]['degenerating'])):
                    newl.append(a)
            riskNodes = riskNodes + newl
            
            # check that there are any more nodes at risk
            if len(riskNodes)==0:
                break
            
#            print(toxicNodes)
            
        return toxicNodes
                
                
                
        
        


    def contiguousspreadOld(self, edgeloss, spreadratio, largestconnectedcomp=False):
        """
        removes edges in a contiguous fashion
        
        edges = number of edges to remove each time
        spreadratio = Ratio of mean degree to the recruitment for new nodes.
                      New nodes are added [to what??] when spreadratio*mean degree number of edges are lost 
                      Effectively a measure of the speed of diffusitivity.
        """
        self.lengthEdgesRemoved = []
        edgesleft = edgeloss
        if not self.iter:
            self.iter = 0
        # get group of degenerating nodes
        # if the ratio to the mean degree requires more nodes, add some more now
        # find mean node degree and ratio of mean degree
        meandegree = np.mean(self.G.degree().values())
        ratio_meandegree  = meandegree*spreadratio
        while edgesleft > 0:
            try:
                degennodes = [v for v in self.G.nodes() if self.G.node[v]['degenerating']]
                
            except:
                # if no degenerating nodes, pick a random node to start
                for node in self.G.nodes():
                    self.G.node[node]['degenerating'] = False
                
                connectednodes = [v for v in self.G.nodes() if nx.degree(self.G,v)!=0]
                
                self.startingnode = random.choice(connectednodes)
    
                self.G.node[self.startingnode]['degenerating'] = True
                degennodes = [self.startingnode]
                self.edgeloss = 0
                
                # number of nodes to add is a ratio of the mean degree
                numNodestoAdd = ratio_meandegree
                
                while numNodestoAdd > 0:
                    try:
                        neighbours = []
                        for node in degennodes:
                            neighbours.append( [ v for v in self.G.neighbors(node) ] )
                        
                        neighbours = list(itertools.chain.from_iterable(neighbours))  # converts list of lists in to list
                        neighbours = set(neighbours)
                        neighbours = [v for v in neighbours if not v in degennodes]
                        
                        newnode = random.choice(neighbours)
                        
                        self.G.node[newnode]['degenerating'] = True
                        degennodes.append(newnode)
                        
                        numNodestoAdd -= 1
                    
                    except:
                        nextnode = self.findSpatiallyNearest(degennodes)
                        degennodes.append(nextnode)
                        numNodestoAdd -= 1
    
                    
            # identify edges of all degenerating nodes
            degenedges = self.G.edges(degennodes)
            
            # check to see if there are any edges left in the graph if no degenerating edges exist
            if not degenedges:
                if not self.G.edges():
                    sys.exit("No edges left in graph")
                else:
                    pass
            
            try:
                degenedges_toremove = random.sample(degenedges, edgeloss)
            except:
                degenedges_toremove = None
                    
            # recruit neighbouring nodes to degenerate only if number of edges lost is greater than a ratio of the mean degree
            # or if there are no edges to degenerate                
            if (not degenedges_toremove) or (self.edgeloss >= ratio_meandegree):
                neighbours = []
                numnewrecruits = 1
                
                if self.edgeloss >= ratio_meandegree:
                    numnewrecruits = int(self.edgeloss / ratio_meandegree)
                    self.edgeloss = 0 # resets the counter
                    
                    for node in degennodes:
                        for v in self.G.neighbors_iter(node):
                            neighbours.append(v)
                    
                    neighbours = set(neighbours)            
                    neighbours = [v for v in neighbours if not v in degennodes]
                
                for i in range(numnewrecruits):
                    # first step, try and find nodes topologically connected
                    try:
                        newrecruit = random.sample(neighbours,numnewrecruits)
                        self.G.node[newrecruit]['degenerating'] = True
                    
                    # if no topologically connected nodes, find nearest spatial neighbour
                    except:
                        newrecruit = self.findSpatiallyNearest(degennodes)
                        if not newrecruit:
                            print "The graph has no non-degenerating nodes left with a non-zero degree"
                            sys.exit("No edges left in graph")
                            
                        try:
                            self.G.node[newrecruit]['degenerating'] = True
                        except:
                            continue
                            
                        # if we're looking for new neighbours purely to remove an edge, enter this loop
                        if not degenedges_toremove:
                            degenedges_toremove = random.sample(self.G.edges(self.G.neighbors(newrecruit)), 1)
                                        
            
            # remove 
            if edgesleft - len(degenedges_toremove) < 0:
                degenedges_toremove = [v for v in degenedges_toremove[0:edgesleft]]
            self.G.remove_edges_from(degenedges_toremove)
            self.lastnode = degenedges_toremove[0][1]
            self.edgeloss += edgeloss
            
    
            try: # records length of edge removal if spatial information is available
                for edge in degenedges_toremove:
                    print edge
                    print np.linalg.norm(np.array(self.G.node[edge[0]]['xyz']) - np.array(self.G.node[edge[1]]['xyz']))
                    self.lengthEdgesRemoved.append(np.linalg.norm(np.array(self.G.node[edge[0]]['xyz']) - np.array(self.G.node[edge[1]]['xyz'])))
    
            except:
                pass
    
                  
            self.iter += 1
            
            if largestconnectedcomp:
                self.bigconnG = components.connected.connected_component_subgraphs(self.G)[0]  # identify largest connected component
            edgesleft -= len(degenedges_toremove)        


    ## =============================================================
    
    ## Analysis functions

    
    def findSpatiallyNearestOld(self, degennodes):
        # find the spatially closest node as no topologically close nodes exist
        print "Finding spatially closest node"
        duffNode = random.choice(degennodes)
        nodes = [v for v in self.G.nodes() if v!=duffNode]
        nodes = [v for v in nodes if not v in degennodes]

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
        
    
    def findSpatiallyNearest(self, nodeList, threshold):
        ''' find the spatially nearest nodes to each node within a treshold '''
        
        a = 1
        
        
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
                                                

    def hubIdentifier(self, weighted=False):
        """ 
        define hubs by generating a hub score, based on the sum of normalised scores for:
            betweenness centrality
            closeness centrality
            degree
        
        hubs are defined as nodes 2 standard deviations above the mean hub score
        
        Changelog 7/12/12:
            - added possibility of weighted measures
        """
        
        self.hubs = []
        
    #    get centrality measures
        if weighted:
            betweenessCentrality = centrality.betweenness_centrality(self.G, weight='weight')
            closenessCentrality = centrality.closeness_centrality(self.G, distance=True)
            degrees = nx.degree(self.G, weight='weight')
            
            
        else:
            betweenessCentrality = centrality.betweenness_centrality(self.G)
            closenessCentrality = centrality.closeness_centrality(self.G)
            degrees = nx.degree(self.G)
            
        sum_degrees = np.sum(degrees.values())
        sum_betweenness = np.sum(betweenessCentrality.values())        
        sum_closeness = np.sum(closenessCentrality.values())
        
    #    combine normalised measures for each node to generate a hub score
        hubScores = []
        for node in self.G.nodes():
            if weighted:
                self.G.node[node]['hubScore'] = betweenessCentrality[node]/sum_betweenness + closenessCentrality[node]/sum_closeness + degrees[node]/sum_degrees
            else:
                self.G.node[node]['hubScore'] = betweenessCentrality[node]/sum_betweenness + closenessCentrality[node]/sum_closeness + degrees[node]/sum_degrees
                
            
            hubScores.append(self.G.node[node]['hubScore'])
            
    #   find standard deviation of hub score
        upperLimit = np.mean(np.array(hubScores)) + 2*np.std(np.array(hubScores))
    
    #   identify nodes as hubs if 2 standard deviations above hub score
        for node in self.G.nodes():
            if self.G.node[node]['hubScore'] > upperLimit:
                self.hubs.append(node)
                
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
    
        # set layout position for plotting
        xy = (0,400)
        try:
            angle = 360/len(self.clusternames)
        except:
            angle = 180
        
        points = {0:xy}
        for n in range(1,len(self.clusternames)):
            x = points[n-1][0]
            y = points[n-1][1]
            points[n] = (x*np.cos(angle)-y*np.sin(angle),x*np.sin(angle)+y*np.cos(angle))
        
        self.pos = {}
        
        for clust in self.clusternames:
            clusternodes = [v for v in self.G.nodes() if self.G.node[v]['cluster']==clust]
            clusteredges = [v for v in self.G.edges(clusternodes) if v[0] in clusternodes and v[1] in clusternodes]
            
            subgraph = nx.Graph()
            subgraph.add_nodes_from(clusternodes)
            subgraph.add_edges_from(clusteredges)
            
            if drawpos:
                centre = points[clust]
                
                clusterpos = nx_agraph.graphviz_layout(subgraph,prog='neato')
               
                for node in clusternodes:
                    self.pos[node] = (clusterpos[node][0]+centre[0],clusterpos[node][1]+centre[1])
    
        # calculate modularity
        if clusters:
            self.modularity = community.modularity(clusters,self.G)
        else:
            self.modularity = 0
            
    def degenerate(self, weightloss=0.1, edgesRemovedLimit=1, weightLossLimit=None, toxicNodes=None, riskEdges=None, spread=False):
        ''' remove random edges from connections of the toxicNodes set, or from the riskEdges set. This occurs either until edgesRemovedLimit
        number of edges have been removed (use this for a thresholded weighted graph), or until the weight loss
        limit has been reached (for a weighted graph). For a binary graph, weight loss should be set
        to 1.
        
        The spread option recruits connected nodes of degenerating edges to the toxic nodes list.
        
        By default this function will enact a random attack model, with a weight loss of 0.1 each iteration.
        '''        
        try:
            self.degenEdges
        except AttributeError:
            self.degenEdges = []
            
            
        # set limit
        if weightLossLimit:
            limit = weightLossLimit
        
        else:
            limit = edgesRemovedLimit
        
        if not riskEdges:
            # if no toxic nodes defined, select the whole graph
            if not toxicNodes:
                toxicNodes = self.G.nodes()
            
            # make all toxic nodes degenerating
            for t in toxicNodes:
                self.G.node[t]['degenerating'] = True
                          
             # generate list of at risk edges
            riskEdges = nx.edges(self.G, toxicNodes)
        
        # iterate number of steps
        while limit>0:
            if not riskEdges:
                # find spatially closest nodes if no edges exist
                # is it necessary to do this for all nodes?? - waste of computing power,
                # choose node first, then calculated spatially nearest of a single node
                newNode = self.findSpatiallyNearest(toxicNodes)
                if newNode:
                    print "Found spatially nearest node"
                    toxicNodes.append(newNode)
                    riskEdges = nx.edges(self.G, toxicNodes)
                else:
                    print "No further edges to degenerate"
                    break
            
            # choose at risk edge to degenerate from           
            dyingEdge = random.choice(riskEdges)            
            if not dyingEdge in self.degenEdges:
                self.degenEdges.append(dyingEdge)
            
            # remove specified weight from edge
            w = self.G[dyingEdge[0]][dyingEdge[1]]['weight']
            
            if np.absolute(w) < weightloss:
                loss = w
                self.G[dyingEdge[0]][dyingEdge[1]]['weight'] = 0
            
            elif w>0:
                loss = weightloss
                self.G[dyingEdge[0]][dyingEdge[1]]['weight'] -= weightloss
                
            else:
                loss = weightloss
                self.G[dyingEdge[0]][dyingEdge[1]]['weight'] += weightloss
            
            # add nodes to toxic list if the spread option is selected
            if spread:
                for node in dyingEdge:
                    if not node in toxicNodes:
                        toxicNodes.append(node)
            
            # remove edge if below the graph threshold
            if self.G[dyingEdge[0]][dyingEdge[1]]['weight'] < self.threshold and self.threshold != -1:      # checks that the graph isn't fully connected and weighted, ie threshold = -1
                self.G.remove_edge(dyingEdge[0], dyingEdge[1])
                print ' '.join(["Edge removed:",str(dyingEdge[0]),str(dyingEdge[1])])
                if not weightLossLimit:
                    limit-=1
            
            if weightLossLimit:
                limit -= loss
                
            # redefine at risk edges
            riskEdges = nx.edges(self.G, toxicNodes)
        
        print "Number of toxic nodes: "+str(len(toxicNodes))
                                     
        return toxicNodes
            
            
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
                      
        
    def percentConnected(self):
        totalConnections = len(self.G.nodes()*(len(self.G.nodes())-1))
        self.percentConnections = float(len(self.G.edges()))/float(totalConnections)
    
    def binarise(self):
        binEdges = [v for v in self.G.edges()]
        self.G.remove_edges_from(self.G.edges())
        self.G.add_edges_from(binEdges)
        
    def largestConnComp(self):
        self.bigconnG = components.connected.connected_component_subgraphs(self.G)[0]  # identify largest connected component
        

class plotObj():
    ''' classes that plot various aspects of a brain object '''
    
    def __init__(self):
        
        # initialise mayavi figure
        self.startMayavi()  
        
    
    
        
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


    
    def plotBrain(self, brain, label = None, nodes = None, edges = None, col = (1, 1, 1), opacity = 1.):
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
        s = mlab.points3d(xn, yn, zn, scale_factor = 0.5, color = col, opacity = opacity)
        self.brainNodePlots[label] = s
        
        # plot edges
        t = mlab.quiver3d(xe, ye, ze, xv, yv, zv, line_width = 1., mode = '2ddash', scale_mode = 'vector', scale_factor = 1., color = col, opacity = opacity)
        self.brainEdgePlots[label] = t

    def plotBrainNodes(self, brain, nodes = None, col = (1, 1, 1), opacity = 1.):
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
        s = mlab.points3d(xn, yn, zn, scale_factor = 0.2, color = col, opacity = opacity)
        self.brainNodePlots[label] = s


    def plotBrainEdges(self, brain, label = None, edges = None, col = (1, 1, 1), opacity = 1.):
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
    
    
    def plotSkull(self, brain, label = None, contourVals = [], opacity = 0.1):
        ''' plot the skull using Mayavi '''
        
        if not(label):
            label = self.getAutoLabel()        
        
        if contourVals == []:            
            s = mlab.contour3d(brain.skull, opacity = opacity)
        else:
            s = mlab.contour3d(brain.skull, opacity = opacity, contours = contourVals)
            
        # get the object for editing
        self.skullPlots[label] = s
        
    
    def changePlotProperty(self, plotType, prop, plotLabel, value = 0.):
        ''' change a specified property prop of a plot of type plotType, index is used if multiple plots of
            the same type have been made. Value is used by some properties.
            
            Allowed plotTypes: skull, brainNode, brainEdge
            Allowed props: opacity, visibility, colour
            
            
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
    
    

