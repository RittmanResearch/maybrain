# -*- coding: utf-8 -*-

import numpy as np


def threshold_to_percentage(brain, threshold):
    """
    It returns a ratio between the edges on adjMat above a certain threshold value \
    and the total possible edges of adjMat.
    In an unidrected graph, the total possible edges are the upper right \
    part elements of adjMat different from np.nan. In a directed graph, the total \
    possible edges are all the elements of adjMat except the diagonal and \
    np.nan

    threshold : The threshold value
    """
    upper_values = np.triu_indices(np.shape(brain.adjMat)[0], k=1)
    below_values = np.tril_indices(np.shape(brain.adjMat)[0], k=-1)
    if not brain.directed:
        # only flatten the upper right part of the matrix
        weights = np.array(brain.adjMat[upper_values])
    else:
        weights = np.concatenate((brain.adjMat[upper_values], brain.adjMat[below_values]))
    # remove NaNs
    weights = weights[~np.isnan(weights)]

    max_edges = len(weights)
    len_edges = len(weights[weights > threshold])

    return len_edges / max_edges


def percent_connected(brain):
    """
    This returns the ratio of the current number of edges in our `G` object \
    and the total number of possible connections.
    If N is the number of nodes in our `G` object, the total number of \
    possible connections is (N * (N - 1))/2 for an undirected graph, and \
    N * (N-1) for a directed graph.

    """
    nodes = brain.G.number_of_nodes()

    if brain.directed:
        total_connections = nodes * (nodes - 1)
    else:
        total_connections = nodes * (nodes - 1) / 2
    return float(brain.G.number_of_edges()) / float(total_connections)
