from random import shuffle
import numpy as np
import networkx as nx
from networkx import configuration_model


from maybrain import constants as ct

#TODO: meter lista de atributos em vez de "spatial
#TODO indicar que considera que os edges no a têm que ter o attributo 'weight' ja atribuido para cada edge
#   |-> se calhar meter um teste a ver se o atributo existe, e lançar excepçao se nao?-> meter isto nos testes!
#TODO: Check inVal is not empty
#TODO: criar normalise geral para que isto tudo fique reduzido, e no inicio testar o tipo

def normalise(brain, func, n=500, retNorm=True, inVal=None, **kwargs):
    """
    This function normalises the function G by generating a series of n random
    graphs and averaging the results. If retNorm is specified, the normalised
    value is returned, else a list of n values for a random graph are returned.
    """
    if brain.directed:
        raise TypeError("normalise() not available for directed graphs")
            
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

def normalise_node_wise(brain, func, init_vals=None, n_iter=500, ret_normalised=True, 
                        node_attrs=None, edge_attrs=None, **kwargs):
    """
    It normalises measures taken from a brain by generating a series of n random graphs and
    averaging them.

    brain: an instance of the `Brain` class
    func: the function that calculates the measure in each node of brain
    init_vals: the initial measures calculated from brain.G that will be averaged. If this is None,
               this will be equal to func(brain.G, **kwargs)
    n_iter: number of iteratios that will be used to generate the random graphs
    ret_normalised: if True, the normalised measures are returned
    node_attrs: If the node attributes of brain.G are necessary to a correct calculation of func()
                in the random graph, just pass the attributes as a list of strings in this parameter
    edge_attrs: If the edge attributes of brain.G are necessary to a correct calculation of func()
                in the random graph, just pass the attributes as a list of strings in this parameter.
                `ct.WEIGHT` is already passed to the edges of the random graphs, so need to pass it
                in this parameter
    **kwargs: keyword arguments if you need to pass them to func()
    """
    if brain.directed:
        raise TypeError("normalise_node_wise() not available for directed graphs")
    if init_vals is None:
        init_vals = func(brain.G, **kwargs)
    if node_attrs is None:
        node_attrs = []
    if edge_attrs is None:
        edge_attrs = []
        
    nodes_dict = {v:[] for v in brain.G.nodes()}
    for _ in range(n_iter):
        rand = configuration_model(list(dict(brain.G.degree()).values()))
        rand = nx.Graph(rand) # convert to simple graph from multigraph

        # Reassigning weights from the input graph to the random graph node by node
        # this means that for each node the total weighted degree will be similar
        # nb this depends on the random graph having the same or fewer vertices for each node
        shift = 0 # Some nodes might not exist in brain.G. rand.nodes() will go [0, G.number_of_nodes()]
        for node in sorted(rand.nodes()):
            # rand.nodes() have all the values from 0 to brain.G.number_of_nodes() 
            #  while brain.G.nodes() might not have all the nodes
            #  Thus, a shift is necessary to map correctly
            while not brain.G.has_node(node + shift):
                shift += 1
            
            for attr in node_attrs:
                rand.nodes[node][attr] = brain.G.nodes[node+shift][attr]
            
            edge_vals = [v[2] for v in brain.G.edges(node+shift, data=True)]
            
            shuffle(edge_vals)
            for en, edge in enumerate(rand.edges(node)):
                rand.edges[edge[0], edge[1]][ct.WEIGHT] = edge_vals[en][ct.WEIGHT]
                for attr in edge_attrs:
                    rand.edges[edge[0], edge[1]][attr] = edge_vals[en][attr] 
        
        # Applying func() to the random graph
        res = func(rand, **kwargs)

        for x, node in enumerate(nodes_dict):
            nodes_dict[node].append(res[x])
    
    if ret_normalised: # update nodes_dict to return the normalised values
        for node in nodes_dict:
            nodes_dict[node] = init_vals[node]/np.mean(nodes_dict[node])
    
    # Returns a list of the values from the random graph
    return nodes_dict
