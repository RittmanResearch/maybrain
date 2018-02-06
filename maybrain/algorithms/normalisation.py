from random import shuffle
import numpy as np
import networkx as nx
from networkx import configuration_model

from maybrain import constants as ct
import numbers


def normalise_single(brain, func, init_val=None, n_iter=500, ret_normalised=True, 
              node_attrs=None, edge_attrs=None, **kwargs):
    """
    See `normalise()` method's documentation for explanation. This method just expects a single
    initial measure (init_val) to be averaged, instead of a dictionary
    """
    if brain.directed:
        raise TypeError("normalise_single() not available for directed graphs")
    if not isinstance(init_val, numbers.Real):
        raise TypeError("normalise_single() expects init_val to be a real number")
        
    return normalise(brain, func, init_vals=init_val, n_iter=n_iter, ret_normalised=ret_normalised, 
                     node_attrs=node_attrs, edge_attrs=edge_attrs, **kwargs)


def normalise_node_wise(brain, func, init_vals=None, n_iter=500, ret_normalised=True, 
              node_attrs=None, edge_attrs=None, **kwargs):
    """
    See `normalise()` method's documentation for explanation. This method just expects init_vals
    to be a dictionary with the measures, for each node, that will be averaged
    """
    if brain.directed:
        raise TypeError("normalise() not available for directed graphs")
    if not isinstance(init_vals, dict):
        raise TypeError("normalise_node_wise() expects init_vals to be a dictionary")
        
    return normalise(brain, func, init_vals=init_vals, n_iter=n_iter, ret_normalised=ret_normalised, 
                     node_attrs=node_attrs, edge_attrs=edge_attrs, **kwargs)

def normalise(brain, func, init_vals=None, n_iter=500, ret_normalised=True, 
              node_attrs=None, edge_attrs=None, **kwargs):
    """
    It normalises measures taken from a brain by generating a series of n random graphs and
    averaging them.

    brain: an instance of the `Brain` class
    func: the function that calculates the measure in each node of brain
    init_vals: the initial measures calculated from brain.G that will be averaged. 
                If this is None, this will be equal to func(brain.G, **kwargs)
                If this is a dictionary, this will be a normalisation node-wise
                Otherwise, it will be treated as a single measure to be averaged
    n_iter: number of iteratios that will be used to generate the random graphs
    ret_normalised: if True, the normalised measures are returned, otherwise a list of n values 
                    for a random graph is returned
    node_attrs: If the node attributes of brain.G are necessary to a correct calculation of func()
                in the random graph, just pass the attributes as a list of strings in this parameter
    edge_attrs: If the edge attributes of brain.G are necessary to a correct calculation of func()
                in the random graph, just pass the attributes as a list of strings in this parameter.
                `ct.WEIGHT` is already passed to the edges of the random graphs, so no need to pass 
                it in this parameter
    **kwargs: keyword arguments if you need to pass them to func()
    """
    if brain.directed:
        raise TypeError("normalise() not available for directed graphs")
    if init_vals is None:
        init_vals = func(brain.G, **kwargs)
    if node_attrs is None:
        node_attrs = []
    if edge_attrs is None:
        edge_attrs = []
    # Checking there is the attribute ct.WEIGHT in the edges
    try:
        if list(brain.G.edges(data=True))[0][2][ct.WEIGHT]:
            pass
    except KeyError as error:
        import sys
        _, _, tb = sys.exc_info()
        raise KeyError(error, "Edge doesn't have ct.WEIGHT property").with_traceback(tb)
    
    if isinstance(init_vals, dict):
        nodes_dict = {v:[] for v in brain.G.nodes()}
    else:
        vals = []
        
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

        if isinstance(init_vals, dict):
            for x, node in enumerate(nodes_dict):
                nodes_dict[node].append(res[x])
        else:
            vals.append(res)
    
    
    if isinstance(init_vals, dict):
        if ret_normalised: # update nodes_dict to return the normalised values
            for node in nodes_dict:
                nodes_dict[node] = init_vals[node]/np.mean(nodes_dict[node])
        return nodes_dict
    else:
        if ret_normalised:
            return init_vals / np.mean(vals)
        else:
            return vals
