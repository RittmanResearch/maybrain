# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19 22:25:22 2012

Functions that used to be at the end of networkutils_bin_camb

"""
from os import path,rename
import numpy as np
import networkx as nx
from matplotlib import pyplot as plt
from numpy import linalg as lg
from numpy import fill_diagonal

def efficiencyfunc(node, G, weight=None):
    pls = nx.shortest_path_length(G, source=node, weight=weight)
    del(pls[node])
                                        
    invpls = [1/float(v) for v in pls.values()]
    invpl = np.sum(invpls)
    
    return(invpl)
    
def globalefficiency(G, weight=None):
    """
    A set of definitions to calculate global efficiency and local efficiency.
    Definitions are taken from Latora and Marchiori 2001, Physical Review
    Letters. 
    """

    N = float(len(G.nodes())) # count nodes
    ssl = 0            # sum of inverse of the shortest path lengths
    
    ssl = [ efficiencyfunc(v, G, weight=weight) for v in G.nodes() ]
    ssl = np.sum(np.array((ssl)))
    
    if N>1:
        Geff = float(ssl) / (N*N-N)
        return Geff
    else:
        "Number of nodes <1, can't calculate global efficiency"
        return None

def localefficiency(G, nodes=None, weight=None):
    """
    Returns a dictionary of local efficiency values for each node in the graph.
    """    
    # if node list is not specified, use all nodes in the graph
    if not nodes:
        nodes = G.nodes()
    
    outDict={}

    for node in nodes:
        Ginodes = nx.neighbors(G,node)
        Ginodes.append(node)
        Giedges = G.edges(Ginodes)

        Giedges = [ e for e in Giedges if e[1] in Ginodes ]
        Giedges = [ e for e in Giedges if e[0] in Ginodes ]

        Gi = nx.Graph()
        Gi.add_edges_from(Giedges)
        for e in Gi.edges():
            Gi.edge[e[0]][e[1]] = G.edge[e[0]][e[1]]

        if Giedges:
            ssl = [ efficiencyfunc(v,
                                   Gi,
                                   weight=weight) for v in Gi.nodes() if not v==node]
            
            # correct by number of edges in local graph
            locEff = np.sum(ssl) / (len(Ginodes)-1)  # -1 because self connections are not expected
            outDict[node] = locEff
        else:
            outDict[node] =  0.
    
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
    
  
def edgeLengths(G, nodeWise=False):
    el =  dict(zip(G.edges(),
                            [ np.absolute(np.linalg.norm(np.array(G.node[edge[0]]['xyz']) - np.array(G.node[edge[1]]['xyz']))) for edge in G.edges()]))
    
    if nodeWise:
        eln = {}
        for n in G.nodes():
            eln[n] = np.sum([el[e] for e in el.keys() if n in e])
        return(eln)
    else:
        return(el)

def hubscore(G, bc=None, cc=None, degs=None, weighted=False, Gaussian=False):
    """ 
    define hubs by generating a hub score, based on the sum of normalised scores for:
        betweenness centrality
        closeness centrality
        degree
    
    hubs are defined as nodes 2 standard deviations above the mean hub score
    
    defines hubs
    
    if assign is true, then each node's dictionary is assigned a hub score
    
    sdT defines the number of standard deviations above the mean to define a
    node as a hub (default is 2)
    
    Changelog 7/12/12:
        - added possibility of weighted measures
    """
    
#    get centrality measures
    if weighted:
        if not all([bc, cc, degs]):
            if not 'distance' in G.node[G.nodes()[0]].keys():
                for edge in G.edges():
                    G.edge[edge[0]][edge[1]]["distance"] = 1.00001 - G.edge[edge[0]][edge[1]]["weight"] # convert weights to a positive distance

            bc = nx.centrality.betweenness_centrality(G, weight='distance')
            cc = nx.centrality.closeness_centrality(G, distance="distance")
            degs = nx.degree(G)
        
    else:
        if not all([bc, cc, degs]):
            bc = nx.centrality.betweenness_centrality(G)
            cc = nx.centrality.closeness_centrality(G)
            degs = nx.degree(G)
            
    bc = np.array([bc[v] for v in G.nodes()])
    cc = np.array([cc[v] for v in G.nodes()])
    degs = np.array([degs[v] for v in G.nodes()])
    
    if Gaussian:
        # normalise values to mean 0, sd 1
        bc -= np.mean(bc)
        cc -=  np.mean(cc)        
        degs -= np.mean(degs)
        
        bc /= np.std(bc)
        cc /=  np.std(cc)        
        degs /= np.std(degs)
    
    else:
        # rescale between 0 and 1
        bc -= np.min(bc)
        cc -=  np.min(cc)        
        degs -= np.min(degs)
        
        bc /= np.max(bc)
        cc /=  np.max(cc)        
        degs /= np.max(degs)

    hubscores = {node:bc[n]+cc[n]+degs[n] for n,node in enumerate(G.nodes())}
    return(hubscores)

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
    
def withinModuleDegree(G, ci, weight=None):
    """
    To calculate mean within module degree
    """
    moduleList = [v for v in set(ci.values())] # get module list
    withinDegDict = {}  # output dictionary of modules and mean within module degree
    modDict = {m:[v for v in ci.keys() if ci[v]==m] for m in moduleList} # sort nodes in to modules
    
    for n in G.nodes():
        m = ci[n] # iterate through modules
        
        eList = G.edges([n])
        eList = [e for e in eList if all([ci[e[0]]==m, ci[e[1]]==m])] # find edges exclusively within the module
        
        if weight:
            wts = np.sum([float(G.edge[e[0]][e[1]]['weight']) for e in eList])  # get weights/degree
        else:
            wts = float(len(eList))
        
        if len(modDict[m]) > 1:
            withinDegDict[n] = wts/(len(modDict[m])-1)  # mean of weight/degree, ie average degree within module
        else:
            withinDegDict[n] = 0.
            
    return(withinDegDict)
    
def writeResults(results, measure,
                 outfilebase="brain",
                 append=True,
                 propDict=None):
    """ 
    Function to write out results    
    """
    outFile = outfilebase+measure+'.txt'
    if not path.exists(outFile) or not append:
        if path.exists(outFile):
            rename(outFile, outFile+'.old')
        f = open(outFile, "w")
        writeHeadFlag=True

    else:
        f = open(outFile, "a")
        writeHeadFlag=False
    
    # check to see what form the results take
    if isinstance(results, (dict)):
        headers = [v for v in results.keys()]
        headers.sort()

        out = ' '.join([ str(results[v]) if v in results.keys() else 'NA' for v in headers] )
        
    elif isinstance(results, (list)):
        if writeHeadFlag:
            headers = ' '.join([str(v) for v in range(len(results))])
        out = ' '.join([str(v) for v in results])
    
    else:
        if writeHeadFlag:
            headers = [measure]
        out = str(results)

    # write headers
    if propDict:
        propHeaders = propDict.keys()
        propHeaders.sort()

    if writeHeadFlag:
        headers = ' '.join([str(v) for v in headers])
        if propDict:
            headers = ' '.join([' '.join(propHeaders), headers])
            
        f.writelines(headers+'\n')

    # add on optional extras
    if propDict:
        outProps = ' '.join([str(propDict[v]) for v in propHeaders])
        out = ' '.join([outProps, out])
        
    f.writelines(out+'\n')
    f.close()
    