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