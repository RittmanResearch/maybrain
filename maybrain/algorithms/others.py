# TODO: The functions defined here need to be analysed, moved to other modules, or removed (for example if they are already implemented elsewhere
# That's why they are not included by default in the __init__ from the package


def efficiencyfunc(node, G, weight=None):
    pls = nx.shortest_path_length(G, source=node, weight=weight)
    del(pls[node])

    invpls = [1/float(v) for v in list(pls.values())]
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
        print("You have not set a threshold on the graph, so all edge lengths will be included giving a normal distribution")
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