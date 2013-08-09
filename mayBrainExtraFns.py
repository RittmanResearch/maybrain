# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19 22:25:22 2012

Functions that used to be at the end of networkutils_bin_camb

"""
import os
import numpy as np
import networkx as nx
#import csv
from networkx.algorithms import cluster
from networkx.algorithms import centrality
import random
from networkx.algorithms import components
#from matplotlib import pyplot as plt
from numpy import linalg as lg
from numpy import fill_diagonal

def efficiencyfunc(node, G):
    # nicked from nx.single_source_shortest_path_length to sum shortest lengths
    seen={}                  # level (number of hops) when seen in BFS
    level=0                  # the current level
    nextlevel={node:1}  # dict of nodes to check at next level
    while nextlevel:
        thislevel=nextlevel  # advance to next level
        nextlevel={}         # and start a new list (fringe)
        for v in thislevel:
            if v not in seen: 
                seen[v]=level # set the level of vertex v
                nextlevel.update(G[v]) # add neighbors of v
        level=level+1
    if sum(seen.values())>0:
        invpls = [1/float(v) for v in seen.values() if v!=0 ]
        invpl = np.sum(invpls)
    else:
        invpl = 0.
        
    return(invpl)
    
    
def globalefficiency(G):
    """
    A set of definitions to calculate global efficiency and local efficiency. Definitions are taken from Latora and Marchiori
    2001, Physical Review Letters. 
    """

    N = len(G.nodes()) # count nodes
    ssl = 0            # sum of inverse of the shortest path lengths
    
    ssl = [ efficiencyfunc(v, G) for v in G.nodes() ]
    ssl = np.sum(np.array((ssl)))
    
    if N>1:
        Geff = float(ssl) / (float(N)*(float(N-1)))
        return Geff
    else:
        "Number of nodes <1, can't calculate global efficiency"
        return None
    

def localefficiency(G, nodes=None):
    """
    Returns a dictionary of local efficiency values for each node in the graph.
    """        
    
    # if node list is not specified, use all nodes in the graph
    if not nodes:
        nodes = G.nodes()
    
    outDict={}
    for node in nodes:
        Ginodes = nx.neighbors(G,node)
        Giedges = G.edges(Ginodes)
        Gi = nx.Graph()
        Gi.add_nodes_from(Ginodes)
        Gi.add_edges_from(Giedges)
        
        if Gi.edges():
            ssl = efficiencyfunc(node,Gi)
            
            # correct by number of edges in local graph
            locEff = ssl / (len(Gi.nodes()) * (len(Gi.nodes())-1))
        
            outDict[node] = locEff
        else:
            outDict[node] = 0.
    
    return(outDict)
            
def nodalefficiency(G, nodes=None):
    """
    Returns a dictionary of nodal efficiency values for each node in the graph.
    """        
    
    # if node list is not specified, use all nodes in the graph
    if not nodes:
        nodes = G.nodes()
    
    outDict={}
    for node in nodes:
        nodEff = efficiencyfunc(node,G)
        outDict[node] = nodEff / (len(nodes) * (len(nodes)-1))
    
    return(outDict)
    

def smallworldparameters(brain, outfilebase = "brain", append=True, smallWorldness=False, SWiters=50):
    """A calculation for average cluster coefficient and average shortest path length (as defined in Humphries
    2008 http://www.plosone.org/article/info:doi/10.1371/journal.pone.0002051).
    """
    outfile = outfilebase+'_smallworldparameters'
    
    if not append and os.path.exists(outfile):
        print "Moving existing file to "+outfile+'.old'
        os.rename(outfile,outfile+'.old')
    
    if append and os.path.exists(outfile):
        f = open(outfile,"ab")
    
    else:
        f= open(outfile,"wb")
        if smallWorldness:
            line = 'ClusterCoeff,AvShortPathLength,SmallWorldness\n'
        else:
            line = 'ClusterCoeff,AvShortPathLength\n'
        f.writelines(line)
    
    cc_pl =  [None,None]
    
    # check biggest connected component is defined
    if not brain.bigconnG:
        brain.largestConnComp()

    # calculate parameters
    try:
        cc_pl = [cluster.average_clustering(brain.bigconnG),nx.average_shortest_path_length(brain.bigconnG)]
    
    except:
        if brain.bigconnG.nodes() == []:
            print "No nodes left in network, can not calculate path length, clustering coeff or smallworldness"
        
        else:
            print "No edges, can not calculate shortest path length"
            print "Writing clustering coefficient only"
            cc_pl = [cluster.average_clustering(brain.bigconnG),None]
    
    if smallWorldness:
        # peform multiple iterations for robustness
        randcc = []
        randpl = []
        for i in range(SWiters):
            # generate random brain
            randG = nx.Graph()
            randG.add_nodes_from([v for v in brain.G.nodes()])
            for n in range(len(brain.G.edges())):
                node1 = random.choice(randG.nodes())
                node2 = random.choice([v for v in randG.nodes() if v != node1])
                randG.add_edge(node1, node2)
            
            # find largest connected component
            randG = components.connected.connected_component_subgraphs(randG)[0]
            
            # take random brain metrics
            randcc.append(cluster.average_clustering(randG))
            randpl.append(nx.average_shortest_path_length(randG))
            del(randG)
        
        randcc_pl = (np.mean(randcc), np.mean(randpl))
        print randcc_pl
        print cc_pl
        
        # calculate small worldness
        sw = (cc_pl[0]/randcc_pl[0]) / (cc_pl[1]/randcc_pl[1])
        cc_pl.append(sw)
        
    f.writelines(','.join([str(v) for v in cc_pl])+'\n')
    f.close()
        
## Everything below this line will soon be deprecated. It's currently kept for backwards
## compatibility.     
    

def degreewrite(brain,outfilebase="brain", append=True):
    # node degrees
    outfile = outfilebase+'_degrees_nodes'

    try:
        if not append and os.path.exists(outfile):
            print "Moving existing file to "+outfile+'.old'
            os.rename(outfile,outfile+'.old')
    except:
        pass
    
    if append and os.path.exists(outfile):
        f = open(outfile,"ab")
        
    else:
        f= open(outfile,"wb")
        f.writelines(','.join([str(v) for v in brain.G.nodes()])+'\n')
    
    degDict = brain.G.degree()
    f.writelines(','.join([str(degDict[v]) for v in sorted(degDict.iterkeys()) ])+'\n')

    f.close()
    
    # hub degrees
    outfile = outfilebase+'_degrees_hubs'
    hubidfile = outfilebase+'_degrees_hubs_ids'

    try:
        if not append and os.path.exists(outfile):
            print "Moving existing file to "+outfile+'.old'
            try:
                os.rename(outfile,outfile+'.old')
                os.rename(hubidfile,hubidfile+'.old')
                print "Moving existing degrees file to "+outfile+'.old'
                print "Moving existing hub ID file to "+hubidfile+'.old'
                
            except IOError, error:
                (errorno, errordetails) = error
                print "Error moving files "+errordetails
                print "Carrying on anyway"
    except:
        pass
    
    if append and os.path.exists(outfile):
        f = open(outfile,"ab")
        g = open(hubidfile,"ab")
        
    else:
        f= open(outfile,"wb")
        g = open(hubidfile,"wb")
    
    deghubs = [hub for hub in brain.hubs if hub in brain.G] # hubs within largest connected graph component

    # write hub identifies to file   
    g.writelines(','.join([str(v) for v in brain.hubs])+'\n')
    
    try:
        degs = brain.G.degree(deghubs)
        f.writelines(','.join([ str(degs[v]) for v in sorted(degs.iterkeys()) ])+'\n')
    except:
        print "no hubs in largest connected component"
    f.close()
    g.close()
    
    # write degree histogram
    outfileHist = outfilebase+'_degreehistogram'
    if not append and os.path.exists(outfileHist):
        print "Moving existing file to "+outfileHist+'.old'
        os.rename(outfileHist,outfileHist+'.old')
    
    if append and os.path.exists(outfileHist):
        f = open(outfileHist,"ab")
    
    else:
        f= open(outfileHist,"wb")
        
    degreeHist = nx.degree_histogram(brain.G)
    
    f.writelines(','.join([str(v) for v in degreeHist])+'\n')
    f.close()
        

            

def betweennesscentralitywrite(brain,outfilebase = "brain", append=True):
    """ Calculates node and hub betwenness centralities. For hub centralities there are two files, one with the values in and another
    with the hub identities in corresponding rows.
    """
    
    ## betweenness centrality
    # node centrality
    outfile = outfilebase+'_betweenness_centralities_nodes'
    
    if not append and os.path.exists(outfile):
        print "Moving existing file to "+outfile+'.old'
        os.rename(outfile,outfile+'.old')
    
    if append and os.path.exists(outfile):
        f= open(outfile,"ab")
        
    else:
        f = open(outfile,"wb")
        f.writelines(','.join([str(v) for v in brain.G.nodes()])+'\n')
        
        
    centralities = centrality.closeness_centrality(brain.G)  # calculate centralities for largest connected component
    f.writelines(','.join([str(centralities[v]) for v in sorted(centralities.iterkeys())]) + '\n')                    # write out centrality values
    f.close()
    
    # hub centrality
    outfile = outfilebase+'_betweenness_centralities_hubs'
    hubidfile = outfilebase+'_betweenness_centralities_hubs_ids'
    
    if not append and os.path.exists(outfile):
        print "Moving existing file to "+outfile+'.old'
        try:
            os.rename(outfile,outfile+'.old')
            os.rename(hubidfile,hubidfile+'.old')
            print "Moving existing degrees file to "+outfile+'.old'
            print "Moving existing hub ID file to "+hubidfile+'.old'
            
        except IOError, error:
            (errorno, errordetails) = error
            print "Error moving files "+errordetails
            print "Carrying on anyway"
    
    if append and os.path.exists(outfile):
        f = open(outfile,"ab")
        g = open(hubidfile,"ab")
        
    else:
        f= open(outfile,"wb")
        g = open(hubidfile,"wb")
        
    centhubs = [hub for hub in brain.hubs if hub in brain.G] # hubs within largest connected graph component

    # write hub identifies to file      
    g.writelines(','.join([ str(v) for v in brain.hubs]) + '\n')
    
    hubcentralitieistowrite = dict((n,None) for n in brain.hubs) # empty dictionary to populate with centralities data

    for hub in centhubs:
        hubcentralitieistowrite[hub] = centralities[hub]
        
    f.writelines(','.join([ str(v) for v in sorted(hubcentralitieistowrite.iterkeys())]) +'\n')
    f.close()
    g.close()
    
def closenesscentralitywrite(brain,outfilebase = "brain", append=True):
    """ Calculates node and hub betwenness centralities. For hub centralities there are two files, one with the values in and another
    with the hub identities in corresponding rows.
    """
    ## closeness centrality
    # node centrality
    outfile = outfilebase+'_closeness_centralities_nodes'
    
    if not append and os.path.exists(outfile):
        print "Moving existing file to "+outfile+'.old'
        os.rename(outfile,outfile+'.old')
    
    if append and os.path.exists(outfile):
        f= open(outfile,"ab")
        
    else:
        f = open(outfile,"wb")
        f.writelines(','.join([str(v) for v in brain.G.nodes()])+'\n')
        
        
    centralities = centrality.closeness_centrality(brain.G)  # calculate centralities for largest connected component
    f.writelines(','.join([str(centralities[v]) for v in sorted(centralities.iterkeys())]) + '\n')                    # write out centrality values
    f.close()
    
    # hub centrality
    outfile = outfilebase+'_closeness_centralities_hubs'
    hubidfile = outfilebase+'_closeness_centralities_hubs_ids'
    
    if not append and os.path.exists(outfile):
        print "Moving existing file to "+outfile+'.old'
        try:
            os.rename(outfile,outfile+'.old')
            os.rename(hubidfile,hubidfile+'.old')
            print "Moving existing degrees file to "+outfile+'.old'
            print "Moving existing hub ID file to "+hubidfile+'.old'
            
        except IOError, error:
            (errorno, errordetails) = error
            print "Error moving files "+errordetails
            print "Carrying on anyway"
    
    if append and os.path.exists(outfile):
        f = open(outfile,"ab")
        g = open(hubidfile,"ab")
        
    else:
        f= open(outfile,"wb")
        g = open(hubidfile,"wb")
        
    centhubs = [hub for hub in brain.hubs if hub in brain.G] # hubs within largest connected graph component

    # write hub identifies to file       
    g.writelines(','.join([ str(v) for v in brain.hubs ])+'\n')
    
    hubcentralitieistowrite = dict((n,None) for n in brain.hubs) # empty dictionary to populate with centralities data

    for hub in centhubs:
        hubcentralitieistowrite[hub] = centralities[hub]
        
    f.writelines(','.join([str(hubcentralitieistowrite[v]) for v in sorted(hubcentralitieistowrite.iterkeys())])+'\n')
    f.close()
    g.close()

def GlobalEfficiencywrite(brain,outfilebase = "brain", append=True):
    """
    Writes to a file the output of global efficiencies
    """
    
    outfile = outfilebase+'_GlobalEfficiency'
    if not append and os.path.exists(outfile):
        print "Moving existing file to "+outfile+'.old'
        os.rename(outfile,outfile+'.old')
    
    if append and os.path.exists(outfile):
        f = open(outfile,"ab")
    
    else:
        f= open(outfile,"wb")

    f.writelines("{:0.5f}".format(globalefficiency(brain.G))+'\n')
    f.close()
    
def LocalEfficiencywrite(brain,outfilebase = "brain", append=True, nodes=None):
    """
    Writes to a file the output of local efficiencies
    """
    if not nodes:
        nodes = brain.G.nodes()
            
    outfile = outfilebase+'_LocalEfficiency'
    if not append and os.path.exists(outfile):
        print "Moving existing file to "+outfile+'.old'
        os.rename(outfile,outfile+'.old')
    
    if append and os.path.exists(outfile):
        f = open(outfile,"ab")
    
    else:
        f= open(outfile,"wb")
        
        f.writelines(','.join([str(v) for v in nodes])+'\n')
    
    lEffs = localefficiency(brain.G, nodes=nodes)
    f.writelines(','.join(["{:0.5f}".format(lEffs[v]) for v in nodes])+'\n')
    f.close()

def NodalEfficiencywrite(brain,outfilebase = "brain", append=True, nodes=None):
    """
    Writes to a file the output of local efficiencies
    """
    if not nodes:
        nodes = brain.G.nodes()

    outfile = outfilebase+'_NodalEfficiency'
    if not append and os.path.exists(outfile):
        print "Moving existing file to "+outfile+'.old'
        os.rename(outfile,outfile+'.old')
    
    if append and os.path.exists(outfile):
        f = open(outfile,"ab")
    
    else:
        f= open(outfile,"wb")
        f.writelines(','.join([str(v) for v in nodes])+'\n')

    nEffs = nodalefficiency(brain.G, nodes=nodes)
    f.writelines(','.join(["{:0.5f}".format(nEffs[v]) for v in nodes])+'\n')
    f.close()
  
  

def modularstructure(brain,outfilebase = "brain", redefine_clusters=True, append=True):
    """
    Writes to a file the output of cluster membership and modularity.
    """
    # redefine graph clusters
    if redefine_clusters == True:
        brain.clusters()
    
    # write out modularity
    outfile = outfilebase+'_modularity'
    
    if not append and os.path.exists(outfile):
        print "Moving existing file to "+outfile+'.old'
        os.rename(outfile,outfile+'.old')
    
    if append and os.path.exists(outfile):
        f = open(outfile,"ab")
    
    else:
        f= open(outfile,"wb")

    f.writelines("{:0.5f}".format(brain.modularity)+'\n')
    f.close()

def writeEdgeLengths(brain,outfilebase = "brain", append=True):
    outfile = outfilebase+'_edgeLengths'
    try:
        if not append and os.path.exists(outfile):
            print "Moving existing file to "+outfile+'.old'
            os.rename(outfile,outfile+'.old')
        
        if append and os.path.exists(outfile):
            f = open(outfile,"ab")
        
        else:
            f= open(outfile,"wb")
        
        try:
            f.writelines(','.join([str(v) for v in brain.lengthEdgesRemoved])+'\n')
        except TypeError:
            f.writelines('NA\n')
        f.close()
        
    except AttributeError:
        pass
        
    # mean edge lengths
    edgeLengths = []
    for edge in brain.G.edges():
        edgeLengths.append(np.linalg.norm(np.array(brain.G.node[edge[0]]['xyz']) - np.array(brain.G.node[edge[1]]['xyz'])))
        
    meanEdgeLengths = np.mean(np.absolute(edgeLengths))
    medianEdgeLengths = np.median(np.absolute(edgeLengths))
    
    hubEdgeLengths = []
    for edge in brain.G.edges(brain.hubs):
        hubEdgeLengths.append(np.linalg.norm(np.array(brain.G.node[edge[0]]['xyz']) - np.array(brain.G.node[edge[1]]['xyz'])))

    meanHubLengths = np.mean(np.absolute(hubEdgeLengths))
    medianHubEdgeLengths = np.median(np.absolute(hubEdgeLengths))
    
    outfile = outfilebase+'_meanEdgeLengths'
    
    if not append and os.path.exists(outfile):
        print "Moving existing file to "+outfile+'.old'
        os.rename(outfile,outfile+'.old')
    
    if append and os.path.exists(outfile):
        f = open(outfile,"ab")
            
    else:
        f= open(outfile,"wb")
        f.writelines("MeanEdgeLengths,MeanHubEdgeLengths,MedianEdgeLengths,MedianHubEdgeLengths\n")
        
    f.writelines(','.join([str(v) for v in [meanEdgeLengths, meanHubLengths, medianEdgeLengths, medianHubEdgeLengths]])+'\n')
    f.close()
    
def writeEdgeNumber(brain, outfilebase = "brain", append=True):
    outfile = outfilebase+'_EdgeNumber'
    
    
    if not append and os.path.exists(outfile):
        print "Moving existing file to "+outfile+'.old'
        os.rename(outfile,outfile+'.old')
    
    if append and os.path.exists(outfile):
        f = open(outfile,"ab")
        
    else:
        f= open(outfile,"wb")
        f.writelines("EdgeNumber\n")
        
    edgeNum = len(brain.G.edges())
    f.writelines(str(edgeNum)+'\n')
    f.close()
    
def histograms(brain, outfilebase="brain"):
    """ 
    Produces histograms of association matrix weights and edge lengths.
    Requires spatial information to have been loaded on to the graph.
    """
    # define lengths for each edge
    for edge in brain.G.edges():
        brain.G.edge[edge[0]][edge[1]]['length'] = abs(lg.norm(np.array(brain.G.node[edge[0]]['xyz']) - np.array(brain.G.node[edge[1]]['xyz'])))

    # get a list of weights
    weightList = brain.adjMat.copy()
    fill_diagonal(weightList, 0.)
    weightList = np.array([v for v in weightList.flatten() if not str(v)=="nan"])
    
    # plot the weights
    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    ax1.set_title("Weights")
    ax1.hist(weightList)
    
    # add a second axis
    ax2 = fig.add_axes([0.15, 0.1, 0.7, 0.3])    

    # get a list of lengths
    lengths = [brain.G.edge[v[0]][v[1]]['length'] for v in brain.G.edges()]
    lengths=np.array((lengths))
    
    # Title for lengths
    if not brain.threshold:
        title = "Lengths at no threshold"
        print "You have not set a threshold on the graph, so all edge lengths will be included giving a normal distribution"
    else:
        try:
            title = "Lengths at " + "{0:.1f}".format(brain.edgePC*100) + "% connectivity"
        except AttributeError:
            title = "Lengths at " + "{0.2f}".format(brain.threshold) + " threshold"
    ax2.set_title(title)
    
    # plot the lengths
    ax2.hist(lengths)
    
    # save the figure
    fig.savefig(outfilebase+"_histograms.png")
    
