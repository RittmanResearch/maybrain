from random import shuffle
import numpy as np
import networkx as nx

from maybrain import constants as ct
import numbers


class RandomGenerationError(Exception):
    """
    Exception raised for errors in the generation of random graphs
    """

    def __init__(self, *args, **kwargs):
        Exception.__init__(self, *args, **kwargs)


def generate_random_graph_from_degree(brain, throw_exception=False, node_attrs=None, edge_attrs=None):
    """
    It returns a graph with the same degree sequence as the brain specified as argument.

    This algorithm is an adaptation from Networkx's `degree_seq.configuration_model` in order to have a working version
    for nx.Graph (no selfloops and parallel edges). As it is an adaptation using random shuffle under the hood,
    sometimes it is not possible to have a new random graph with exactly the same degree sequence. The behaviour of
    what to do in this case is specified in `throw_exception` parameter.

    Parameters
    ----------
    brain: maybrain.brain.Brain
        An instance of the `Brain` class
    throw_exception: bool
        If set to True, when the algorithm doesn't find a random graph with the same degree sequence, it throws an
        exception. Usually, depending on the input brain, running this function again will probably solve the problem
        as it is highly dependent on random behaviour.
        If set to False, the returned graph might not have exaclty the same degree distribution as the original input
        brain, because some edges are ignored
    node_attrs: list of str
        If the node attributes of brain.G are necessary for a correct calculation of func()
        in the random graph, just pass the attributes as a list of strings in this parameter.
        This way, these node attributes from the original brain will be presented in the random graph's nodes.
    edge_attrs: list of str
        If the edge attributes of brain.G are necessary for a correct calculation of func()
        in the random graph, just pass the attributes as a list of strings in this parameter.
        This way, these edge attributes from the original brain will be shuffled in the random graph's edges. This
        shuffling is made in a way that each node in the random graph will have edges with similar attributes as in the
        original brain. Although the degree sequence is maintained, as we are creating a random graph, it is impossible
        to have the nodes with the same edges' attributes around it. This is made in a way that all attributes are
        maintained in the random graph, regardless of their exact place. In the case the random graph doesn't maintain
        the same degree distribution (check `throw_exception` parameter), obviously not all the edge's attributes will
        be in the random graph
    Returns
    -------
    new_g: nx.Graph
        A graph with the same, or almost the same, degree sequence of original `brain`

    Raises
    ------
    RandomGenerationError : Exception
        This depends on the value of `throw_exception` parameter
    """
    if node_attrs is None:
        node_attrs = []
    if edge_attrs is None:
        edge_attrs = []

    nodes = []
    new_g = nx.Graph()
    for deg in brain.G.degree():
        # Creating the list with the nodes according to degree
        for _ in range(deg[1]):
            nodes.append(deg[0])
        # Adding node just so random graph will have the same nodes as original one
        if deg[1] == 0:
            new_g.add_node(deg[0])

    shuffle(nodes)

    added = set()
    while len(nodes) >= 2:  # while there are 2 elements or more in nodes
        elem = nodes.pop()
        # Trying to find a pair to this node to form an edge
        for i, n in enumerate(nodes):
            edge = tuple(sorted((elem, n)))
            # if edge is not a selfloop and not already added, we can add to the brain
            if elem != n and edge not in added:
                added.add(edge)
                new_g.add_edge(edge[0], edge[1])
                del nodes[i]

                break  # break the for to come back to the while
            # no pair found for this node
            if throw_exception and i == len(nodes) - 1:  # no pair found for this node
                raise RandomGenerationError(str(elem) + " without pair")

    # Reassigning attributes from the input graph to the random graph node by node
    for node in sorted(new_g.nodes()):
        for attr in node_attrs:
            new_g.nodes[node][attr] = brain.G.nodes[node][attr]

    # Just return the graph straightaway if there are no edge attributes to handle
    if not edge_attrs:
        return new_g

    # Putting the attributes from the original edges into the random graph's edges
    edges = [v for v in brain.G.edges(data=True)]
    # Finding the edges that need the attributes
    # Copying the set is to make it easier to delete elements
    for edg in list(added):
        # Iterating until finding a similar one. Reversed for easier deletion
        for pos, edg_orig in reversed(list(enumerate(edges))):
            if edg_orig[0] == edg[0] or edg_orig[0] == edg[1] or edg_orig[1] == edg[0] or edg_orig[1] == edg[1]:
                # Putting the attributes in the edge
                for attr in edge_attrs:
                    new_g.edges[edg[0], edg[1]][attr] = edg_orig[2][attr]
                del edges[pos]
                added.remove(edg)
                break

    # Add the rest of the edges left randomly
    shuffle(edges)
    for edg in list(added):
        # Iterating until finding a similar one
        edg_orig = edges.pop()
        # Putting the attributes in the edge
        for attr in edge_attrs:
            new_g.edges[edg[0], edg[1]][attr] = edg_orig[2][attr]

    return new_g


def normalise_single(brain, func, init_val=None, n_iter=500, ret_normalised=True, exact_random=False,
                     node_attrs=None, edge_attrs=None, random_location=None, **kwargs):
    """
    See `normalise()` method's documentation for explanation. This method just expects a single
    initial measure (init_val) to be averaged, instead of a dictionary
    """
    if brain.directed:
        raise TypeError("normalise_single() not available for directed graphs")
    if not isinstance(init_val, numbers.Number):
        raise TypeError("normalise_single() expects init_val to be a number")

    return normalise(brain, func, init_vals=init_val, n_iter=n_iter,
                     ret_normalised=ret_normalised, exact_random=exact_random,
                     node_attrs=node_attrs, edge_attrs=edge_attrs, random_location=random_location, **kwargs)


def normalise_node_wise(brain, func, init_vals=None, n_iter=500, ret_normalised=True, exact_random=False,
                        node_attrs=None, edge_attrs=None, random_location=None, **kwargs):
    """
    See `normalise()` method's documentation for explanation. This method just expects init_vals
    to be a dictionary with the measures, for each node, that will be averaged
    """
    if brain.directed:
        raise TypeError("normalise() not available for directed graphs")
    if not isinstance(init_vals, dict):
        raise TypeError("normalise_node_wise() expects init_vals to be a dictionary")

    return normalise(brain, func, init_vals=init_vals, n_iter=n_iter,
                     ret_normalised=ret_normalised, exact_random=exact_random,
                     node_attrs=node_attrs, edge_attrs=edge_attrs, random_location=random_location, **kwargs)


def normalise(brain, func, init_vals=None, n_iter=500, ret_normalised=True, exact_random=False,
              node_attrs=None, edge_attrs=None, random_location=None, **kwargs):
    """
    It normalises measures taken from a brain by generating a series of n random graphs and averaging them.

    Previously generated random graphs can be specified using parameter `random_location`

    Parameters
    ----------
    brain: maybrain.brain.Brain
        An instance of the `Brain` class
    func
        the function that calculates the measure in each node of brain
    init_vals: dictionary or number
        the initial measures calculated from brain.G that will be averaged.
        If this is None, this will be equal to func(brain.G, **kwargs)
        If this is a dictionary, this will be a normalisation node-wise
        Otherwise, it will be treated as a single measure to be averaged
    n_iter: int
        number of iteratios that will be used to generate the random graphs
    ret_normalised: bool
        if True, the normalised measures are returned, otherwise a list of n values
        for a random graph is returned
    exact_random: bool
        This is passed to `throw_exception` argument in `algorithms.generate_random_graph_from_degree()`
    node_attrs: list of str
        If the node attributes of brain.G are necessary for a correct calculation of func()
        in the random graph, just pass the attributes as a list of strings in this parameter.
        This way, these node attributes from the original brain will be presented in the random graph's nodes.
        This is only used if random_location is None
    edge_attrs: list of str
        If the edge attributes of brain.G are necessary for a correct calculation of func()
        in the random graph, just pass the attributes as a list of strings in this parameter.
        `ct.WEIGHT` is already passed to the edges of the random graphs, so no need to pass it in this parameter.
        This way, these edge attributes from the original brain will be shuffled in the random graph's edges.
        This is only used if random_location is None
    random_location: str
        If the random graphs were previously generated, put here the location of them.
        Consider that for each iteration `i = 0...n_iter`, "i" will be added at the end of this path to get
        each random graph. Otherwise, `algorithms.generate_random_graph_from_degree()` will be used to
        create the random graphs
    kwargs
        Keyword arguments if you need to pass them to func()

    Returns
    -------
    vals
        either a dictionary or set, depending on init_vals and ret_normalised

    Raises
    ------
    TypeError: Exception
        If the graph is directed
    KeyError: Exception
        If the edges don't have ct.WEIGHT property
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
        raise KeyError(error, "Edge doesn't have constants.WEIGHT property").with_traceback(tb)

    if isinstance(init_vals, dict):
        nodes_dict = {v: [] for v in brain.G.nodes()}
    elif not isinstance(init_vals, numbers.Number):
        raise TypeError("normalise() expects init_vals to be a number or a dict")
    else:
        vals = []

    for i in range(n_iter):
        if random_location is not None:
            rand = nx.read_gpickle(random_location + str(i))
        else:
            edge_attrs.append(ct.WEIGHT)
            while True:
                try:
                    rand = generate_random_graph_from_degree(brain, throw_exception=exact_random,
                                                             node_attrs=node_attrs, edge_attrs=edge_attrs)
                    break  # if it reaches here, means randomiser didn't throw any exception, so break While
                except RandomGenerationError:
                    pass

        # Applying func() to the random graph
        res = func(rand, **kwargs)

        if isinstance(init_vals, dict):
            for node in nodes_dict:
                nodes_dict[node].append(res[node])
        else:
            vals.append(res)

    if isinstance(init_vals, dict):
        if ret_normalised:  # update nodes_dict to return the normalised values
            for node in nodes_dict:
                nodes_dict[node] = init_vals[node] / np.mean(nodes_dict[node])
        return nodes_dict
    else:
        if ret_normalised:
            return init_vals / np.mean(vals)
        else:
            return vals
