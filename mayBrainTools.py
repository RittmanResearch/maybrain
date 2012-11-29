# -*- coding: utf-8 -*-
"""
Created on Fri Mar 16 18:10:07 2012

@author: -

Change log from networkutils_binary_cambrid_0_1.py:
    - split brain and plotting classes
    - can define a subBrain from properties of edges or nodes
    - animate function a bit more sophisticated (see animate.py)
    - fixed spelling error of 'visibility' in changePlotProperty
    - separated plotting of nodes and edges (backwards compatible, kept plotBrain function)
    - added colour to list of changeable properties of a plot
    - lists of plots in plotObj converted to dictionaries, with autolabel generation if none supplied
        
To do:
    - write an error function so that errors can be better handled when using a gui??
        
"""

import string,os,community,csv,sys,itertools,operator
import networkx as nx
import numpy as np
from networkx.drawing import *
from networkx.generators import random_graphs
from networkx.algorithms import cluster
from networkx.algorithms import centrality
from networkx.algorithms import components
from numpy import shape, random
from mayavi import mlab
from string import split
import nibabel as nb
import mayBrainExtraFns as fns

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
    
    def __init__(self, largestconnectedcomp=False, drawpos=False):
        ''' 
        initialise the brain model. Arguments are:
            - adjmat: adjacency matrix filename
            - threshold: threshold of adjacency for plotting edges
            - spatialinfoFile: filename for spatial info
            - largestconnectedcomp: compute the largest conneted component
            - drawpos: draw 2D position ??? What does this mean ???
    
        '''        
        
        # create an empty graph
        self.G = nx.Graph()
        self.iter = None
        self.largestconnectedcomp = largestconnectedcomp
        self.drawpos = drawpos
        
#        # read input files
#        self.readAdjFile(adjmat, threshold)
#        
#        # read spatial information
#        if spatialinfoFile!='':
#            self.readSpatialInfo(spatialinfoFile)
        
#        # identify largest connected component
#        if self.largestconnectedcomp:
#            self.bigconnG = components.connected.connected_component_subgraphs(self.G)[0]
#            
#        # draw 2D position
#        if self.drawpos:
#            self.pos = nx_agraph.graphviz_layout(self.G,prog='neato')
        
      
        

    def readAdjFile(self, fname, threshold):
        ''' load adjacency matrix with filename and threshold for edge definition '''
        
        print 'loading data from', fname

        # open file
        f = open(fname,"rb")
        reader = f.readlines()        

        # get line that data starts in 
        startLine = 0
        for line in reader:
            if 'begins line' in str(line):
                lstr = str(line)
                whereLabel = lstr.find('begins line')
                startLine = int(lstr[whereLabel + 12])
                break
                
        # get data and convert to lists of floats
        linesStr = reader[startLine:]
        lines = []
        for l in linesStr:
            lines.append(map(float, split(l)))
        nodecount = len(lines)                

        # close file                
        f.close()
        
        ## check somewhere here that data is square?? Does it need to be??

        # create nodes
        self.G.add_nodes_from(range(nodecount))  # creates one node for every line in the adjacency matrix
        
        # add edges (if above threshold)
        mat = np.array(lines)       # array of data
        boolMat = mat>threshold
        edgeCos = np.where(boolMat) # lists of where edges should be
        
        for ind in range(len(edgeCos[0])):
            node1 = edgeCos[0][ind]
            node2 = edgeCos[1][ind]
            
            if not(self.G.has_edge(node1, node2)):
                self.G.add_edge(node1, node2, weight = mat[node1, node2])
                
    def inputNodeProperties(self, propertyName, nodeList, propList):
        ''' add properties to nodes, reading from a list of nodes and a list of corresponding properties '''
        
        for ind in range(len(nodeList)):
            n = nodeList[ind]
            p = propList[ind]
            self.G.node[n][propertyName] = p
            
    
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


    def readSpatialInfo(self, fname):
        ''' add 3D coordinate information for each node from a given file '''
        
        coords = {}

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
                reader = csv.reader(f,delimiter=" ")
            except IOError, error:
                (errorno, errordetails) = error
                print "Couldn't find 3D position information"
                print "Problem with opening file: "+errordetails                
            
            # get data from file
            for line in reader:
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
                

            
#                except:
#                    print "Problem adding 3D position information - check the information in the files is correct"
#                    print "Ignoring the problem and carrying on"        
            


        
        
    def importSkull(self, fname):
        ''' Import a file for skull info using nbbabel
            gives a 3D array with data range 0 to 255 for test data
            could be 4d??
            defines an nibabel object, plus ndarrays with data and header info in
        
        '''        
        
        self.nbskull = nb.load(fname)
        self.skull = self.nbskull.get_data()
        self.skullHeader = self.nbskull.get_header()
        print(shape(self.skull))
        
        
            


    
    def findSpatiallyNearest(self, degennodes):
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
                                
    
    def hubIdentifier(self):
        """ 
        define hubs by generating a hub score, based on the sum of normalised scores for:
            betweenness centrality
            closeness centrality
            degree
        
        hubs are defined as nodes 2 standard deviations above the mean hub score
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
        
    #    combine normalised measures for each node to generate a hub score
        hubScores = []
        for node in self.G.nodes():
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
    
    def clusters(self):
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
            
            if self.drawpos:
                centre = points[clust]
                
                clusterpos = nx_agraph.graphviz_layout(subgraph,prog='neato')
               
                for node in clusternodes:
                    self.pos[node] = (clusterpos[node][0]+centre[0],clusterpos[node][1]+centre[1])
    
        # calculate modularity
        if clusters:
            self.modularity = community.modularity(clusters,self.G)
        else:
            self.modularity = 0
    
    def neuronsusceptibility(self, edgeloss=1):
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
        
        if self.largestconnectedcomp:
            self.bigconnG = components.connected.connected_component_subgraphs(self.G)[0]  # identify largest connected component

    def contiguousspread(self, edgeloss, spreadratio):
        """
        removes edges in a contiguous fashion
        
        edges = number of edges to remove each time
        spreadratio = Ratio of mean degree to the recruitment for new nodes.
                      New nodes are added when spreadratio*mean degree number of edges are lost 
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
            
            if self.largestconnectedcomp:
                self.bigconnG = components.connected.connected_component_subgraphs(self.G)[0]  # identify largest connected component
            edgesleft -= len(degenedges_toremove)
                
    def randomiseGraph(self):
        self.G = nx.gnm_random_graph(len(self.G.nodes()), len(self.G.edges()))
        if self.largestconnectedcomp:
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
        
        
    def percentConnected(self):
        totalConnections = len(self.G.nodes()*(len(self.G.nodes())-1))
        self.percentConnections = float(len(self.G.edges()))/float(totalConnections)
    

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
        s = mlab.points3d(xn, yn, zn, scale_factor = 0.2, color = col, opacity = opacity)
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
    
    

