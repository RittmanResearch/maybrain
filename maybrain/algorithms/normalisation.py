import numpy as np
import networkx as nx
from networkx import configuration_model
from random import shuffle


def normalise(G, func, n=500, retNorm=True, inVal=None, weight='weight'):
    """
    This function normalises the function G by generating a series of n random
    graphs and averaging the results. If retNorm is specified, the normalised
    value is returned, else a list of n values for a random graph are returned.
    """
    vals = []
    for i in range(n):
        rand = configuration_model(list(G.degree().values()))
        rand = nx.Graph(rand) # convert to simple graph from multigraph

        # check whether the graph is weighted
        # if so, the algorithm randomly reassigns weights from the input graph to the random graph node by node
        # this means that for each node the total weighted degree will be similar
        # nb this depends on the random graph having the same or fewer vertices for each node
        if isinstance(G.degree(weight=weight)[0], float):
            randDegs = rand.degree()
            for node in rand.nodes():
                edgeVals = [v[2][weight] for v in G.edges(G.nodes()[node],data=True)]
                shuffle(edgeVals)
                for en,edge in enumerate(rand.edges(node)):
                    rand.edge[edge[0]][edge[1]]['weight'] = edgeVals[en]

        # remove unconnected nodes
        nodeList = [rand.degree()[v] for v in rand.nodes() if rand.degree(v)==0]
        rand.remove_nodes_from(nodeList)

        # add spatial information if necessary
        if spatial:
            print("Adding spatial info")
            nx.set_node_attributes(rand, 'xyz', {rn:G.node[v]['xyz'] for rn,v in enumerate(G.nodes())}) # copy across spatial information

        try:
            vals.append(func(rand),weight='weight') # collect values using the function
        except:
            vals.append(func(rand)) # collect values using the function

    if retNorm: # return the normalised values
        if not inVal:
            inVal = func(G)
        return(inVal/np.mean(vals))

    else: # return a list of the values from the random graph
        return(vals)

def normaliseNodeWise(G, func, inVal={}, n=500, retNorm=True, distance=False, weight='weight', spatial=False):
    """
    This function normalises the function G by generating a series of n random
    graphs and averaging the results. If retNorm is specified, the normalised
    value is returned, else a list of n values for a random graph are returned.
    """
    nodesDict = {v:[] for v in G.nodes()}
    for i in range(n):
        rand = configuration_model(list(G.degree().values()))
        rand = nx.Graph(rand) # convert to simple graph from multigraph

        # check whether the graph is weighted
        # if so, the algorithm randomly reassigns weights from the input graph to the random graph node by node
        # this means that for each node the total weighted degree will be similar
        # nb this depends on the random graph having the same or fewer vertices for each node
        if isinstance(G.degree(weight=weight)[0], float):
            randDegs = rand.degree()
            for node in rand.nodes():
                edgeVals = [v[2][weight] for v in G.edges(G.nodes()[node],data=True)]
                shuffle(edgeVals)
                for en,edge in enumerate(rand.edges(node)):
                    rand.edge[edge[0]][edge[1]]['weight'] = edgeVals[en]

        # remove unconnected nodes
        nodeList = [rand.degree()[v] for v in rand.nodes() if rand.degree(v)==0]
        rand.remove_nodes_from(nodeList)

        #        if distance:
        #            edgeList = [rand.edge[v[0]][v[1]]['weight'] for v in rand.edges() ]
        #
        #            # get the maximum edge value, plus any negative correction required
        #            # and a small correction to keep the values above zero
        #            # the correction is the inverse of the number of nodes - designed to keep
        #            # calculations of efficiency sensible
        #            eMax = np.max(edgeList) + 1/float(len(rand.nodes()))
        #
        #            for edge in rand.edges():
        #                rand.edge[edge[0]][edge[1]]["distance"] = eMax - rand.edge[edge[0]][edge[1]]["weight"] # convert weights to a positive distance

        # add spatial information if necessary
        if spatial:
            print("Adding spatial info")
            nx.set_node_attributes(rand, 'xyz', {rn:G.node[v]['xyz'] for rn,v in enumerate(G.nodes())}) # copy across spatial information

        if distance:
            res = func(rand, distance=True) # collect values using the function
        else:
            try:
                # collect values using the function
                # nb it doesn't matter what the specified weight value is (ie weight or distance), the values are stored as 'weight' in the random graph
                res = func(rand, weight='weight')
            except:
                res = func(rand) # collect values using the function

        for x,node in enumerate(nodesDict):
            nodesDict[node].append(res[x])

    if retNorm: # return the normalised values
        for node in nodesDict:
            nodesDict[node] = inVal[node]/np.mean(nodesDict[node])
        return(nodesDict)

    else: # return a list of the values from the random graph
        return(nodesDict)
