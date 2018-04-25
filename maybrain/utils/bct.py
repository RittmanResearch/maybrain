"""
Module with utility functions to integrate maybrain with bctpy
"""
import numpy as np
import networkx as nx


def makebctmat(brain, nonedge=np.nan):
    """
    Create a matrix from brain.G for use with Brain Connectivity Toolbox measures
    (https://github.com/aestrivex/bctpy/)

    Note that missing nodes are not included, so the matrix order in the resulting matrix
    may not match the node number in the maybrain networkx object

    brain: an instance of the `Brain` class
    nonedge: the value to be put in the bct matrix when there is no edge presented

    Parameters
    ----------
    brain: maybrain.brain.Brain
        An instance of the `Brain` class
    nonedge
        The value to be put in the bct matrix when there is no edge presented
    Returns
    -------
    array: np.array
        A connectivity array ready to be used by bctpy
    """
    return np.copy(nx.to_numpy_matrix(brain.G, nonedge=nonedge))


def assignbctresult(brain, bct_res):
    """
    It assigns the results of a Brain Connectivity Toolbox matrix to the respective nodes in
    maybrain

    Parameters
    ----------
    brain: maybrain.brain.Brain
        An instance of the `Brain` class
    bct_res: np.array
        A matrix from bctpy with values to be assigned for each node

    Returns
    -------
    result: dict
        A dictionary in which the keys are the nodes of `brain` and the values are the ones presented in `bct_res`
    """
    return dict(list(zip(brain.G.nodes(), bct_res)))
