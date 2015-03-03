# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19 22:25:22 2012

Functions that used to be at the end of networkutils_bin_camb

"""
from os import path,rename
import numpy as np
import networkx as nx
from networkx.algorithms import cluster, centrality, components
from networkx import configuration_model
from networkx import connected_component_subgraphs as ccs
import random
from numpy import linalg as lg

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

def localefficiency(G, nodes=None, weight=None, fc=False):
    """
    Returns a dictionary of local efficiency values for each node in the graph.
    fc is a bit of a bodge, specify this if the graph is fully connected to
    make it computationally possible to run.
    """    
    # if node list is not specified, use all nodes in the graph
    if not nodes:
        nodes = G.nodes()

    if weight and fc:
        outDict = {v:efficiencyfunc(v, G, weight=weight)/(len(G.nodes())-1) for v in G.nodes()}
                
    else:
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
    
#!! robustness function from master, checkrobustness removed
def robustness(G, iterLen=500, N=50):
    ''' a function to calculate robustness based on "Error and attack
    tolerance of complex networks" Albert et al Nature 2000 406:378-382
    
    The function calculates the rate of change in the size of the largest
    connected component (S) as nodes are randomly removed. The process runs
    iteratively and takes a mean. The gradient of S is smoothed to provide
    a more accurate measure by a sliding window.
    
    N = size of the sliding window for smoothing the gradient
    iterLen = number of iterations
    
    Note, this function is relatively slow compared to other metrics due to
    the multiple iterations.
    
    '''
    fList = np.zeros(iterLen)
    for i in range(iterLen):
        a = sorted(nx.connected.connected_component_subgraphs(G), key=len, reverse=True)[0]
        nList = [v for v in G.nodes()]
        random.shuffle(nList)
        nList = nList[:-1]
        mat = np.zeros((len(nList)))
        count = 0
        S = len(a.nodes())
        
        while nList and S>1:
            n = nList.pop()
            if n in a.nodes():
                a.remove_node(n)
                if not nx.is_connected(a): # only recalculate if the further fragmentation
                    a =sorted(nx.connected.connected_component_subgraphs(a), key=len, reverse=True)[0]
                    S = len(a.nodes())
                else:
                    S-=1
            mat[count] = S            
            count+=1
        
        g = np.gradient(mat)
        runMean = np.convolve(g, np.ones((N,))/N)[(N-1):]
        diffs = np.diff(runMean)
        nr = np.argmin(diffs)
        
        fList[i] = nr
    return(np.mean(fList) / len(G.nodes()))
  
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
        
def histograms(brain, outfilebase="brain"):
    from matplotlib import pyplot as plt
    from numpy import fill_diagonal
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
#    moduleList = [v for v in set(ci.values())] # get module list
    withinDegDict = {}  # output dictionary of nodes and mean within module degree
#    modDict = {m:[v for v in ci.keys() if ci[v]==m] for m in moduleList} # sort nodes in to modules
    
    for n in G.nodes():
        m = ci[n] # select module
        
        eList = G.edges([n])
        eList = [e for e in eList if all([ci[e[0]]==m, ci[e[1]]==m])] # find edges exclusively within the module
        
        if weight:
            wts = np.sum([float(G.edge[e[0]][e[1]]['weight']) for e in eList])  # get weights/degree
            wts = wts/float(len(eList))
        else:
            wts = float(len(eList))
        
        withinDegDict[n] = wts
        
#        if len(modDict[m]) > 1:
#            withinDegDict[n] = wts/(len(modDict[m])-1)  # mean of weight/degree, ie average degree within module
#        else:
#            withinDegDict[n] = 0.
            
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
    
#### deprecated function ###
def writeResultsOld(results, measure,
                 outfilebase="brain",
                 append=True,
                 edgePC=None,
                 threshold=None):
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
    if writeHeadFlag:
        headers = ' '.join([str(v) for v in headers])
        if edgePC:
            headers = ' '.join(['edgePC', headers])
        if threshold:
            headers = ' '.join(['threshold', headers])
        f.writelines(headers+'\n')

    # add on optional extras
    if edgePC:
        out = ' '.join([str(edgePC), out])
    if threshold:
        out = ' '.join([str(threshold), out])
    f.writelines(out+'\n')
    f.close()
    
def normalise(G, func, n=500, retNorm=True, inVal=None, printLCC=False):
    """
    This function normalises the function G by generating a series of n random
    graphs and averaging the results. If retNorm is specified, the normalised 
    value is returned, else a list of n values for a random graph are returned.
    """
    vals = []
    for i in range(n):
        rand = configuration_model(G.degree().values())
        if nx.components.number_connected_components(rand)>1:
            graphList = ccs(rand)
            rand = graphList.next()
            if printLCC:
                print "Selecting largest connected components of "+str(len(rand.nodes()))+" nodes"
        try:
            vals.append(func(rand)) # collect values using the function
        except KeyError: # exception raised if the spatial information is missing
            nx.set_node_attributes(rand, 'xyz', {rn:G.node[v]['xyz'] for rn,v in enumerate(G.nodes())})
            vals.append(func(rand)) # collect values using the function

    if retNorm: # return the normalised values
        if not inVal:
            inVal = func(G)
        return(inVal/np.mean(vals))
        
    else: # return a list of the values from the random graph
        return(vals)
        
def normaliseNodeWise(G, func, n=500, retNorm=True, inVal=None):
    """
    This function normalises the function G by generating a series of n random
    graphs and averaging the results. If retNorm is specified, the normalised 
    value is returned, else a list of n values for a random graph are returned.
    """
    nodesDict = {v:[] for v in G.nodes()}
    for i in range(n):
        rand = configuration_model(G.degree().values())
        rand = nx.Graph(rand) # convert to simple graph from multigraph
        nodeList = [v for v in rand.nodes() if rand.degree(v)==0]
        rand.remove_nodes_from(nodeList)
        
        try:
            res = func(rand) # collect values using the function
        except KeyError: # exception raised if the spatial information is missing
            print "Adding spatial info"
            nx.set_node_attributes(rand, 'xyz', {rn:G.node[v]['xyz'] for rn,v in enumerate(G.nodes())}) # copy across spatial information
            res = func(rand) # collect values using the function
            
        for x,node in enumerate(nodesDict):
            nodesDict[node].append(res[x])

    if retNorm: # return the normalised values
        if not inVal:
            inVal = func(G)
        for node in nodesDict:
            nodesDict[node] = inVal[node]/np.mean(nodesDict[node])
        return(nodesDict)
        
    else: # return a list of the values from the random graph
        return(nodesDict)

def extractCoordinates(template, outFile="ROI_xyz.txt"):
    import nibabel as nb
    import csv
    
    # load image
    f = nb.load(template)
    fData = f.get_data()
    
    # extract voxel coordinates
    I,J,K=fData.shape
    coords = np.zeros((I, J, K, 3))
    coords[...,0] = np.arange(I)[:,np.newaxis,np.newaxis]
    coords[...,1] = np.arange(J)[np.newaxis,:,np.newaxis]
    coords[...,2] = np.arange(K)[np.newaxis,np.newaxis,:]
    
    # convert voxel coordinates to mm values using affine values
    M = f.affine[:3,:3]
    abc = f.affine[:3,3]
    
    for x in range(len(coords[:,0,0,0])):
        for y in range(len(coords[0,:,0,0])):
            for z in range(len(coords[0,0,:,0])):
                coords[x,y,z,:] = M.dot(coords[x,y,z,:]) + abc
                
    # get unique values in the mask
    valList = np.unique(fData)
    valList = valList[valList!=0]
    
    out = open(outFile, "w")
    writer = csv.DictWriter(out,
                            fieldnames = ["Node", "x", "y", "z"],
                            delimiter=" "
                            )
    for v in valList:
        tempArr = np.zeros((I, J, K, 3), dtype=bool)
        tfArray = fData==v
        tempArr[...,0] = tfArray
        tempArr[...,1] = tfArray
        tempArr[...,2] = tfArray
        
        tempArr = coords[tempArr]
        tempArr = np.mean(tempArr.reshape([tempArr.shape[0]/3,3]),axis=0)
        outList = [str(int(v))]
        outList.extend(["{:.2f}".format(x) for x in tempArr])
        outDict = dict(zip(writer.fieldnames,
                           outList))
        writer.writerow(outDict)
        
    out.close()
