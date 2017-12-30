import numpy as np
from maybrain import constants as ct


def makebctmat(brain):
    """
    Create a matrix for use with brain connectivity toolbox measures.
    See https://pypi.python.org/pypi/bctpy

    Note that missing nodes are not included, so the matrix order in
    the resulting matrix may not match the node number in the maybrain
    networkx object
    """
    bctmat = np.zeros((len(brain.G.nodes()), len(brain.G.nodes())))
    node_indices = dict(list(zip(brain.G.nodes(), list(range(len(brain.G.nodes()))))))
    for nn, x in enumerate(brain.G.nodes()):
        for y in list(brain.G.edge[x].keys()):
            try:
                bctmat[nn, node_indices[y]] = brain.G.edge[x][y][ct.WEIGHT]
            except:
                pass
    return bctmat


def assignbctresult(brain, bct_res):
    """ translate a maybrain connectome into a bct compatible format """
    out = dict(list(zip(brain.G.nodes(), bct_res)))
    return out