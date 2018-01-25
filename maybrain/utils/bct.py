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
    """
    return np.copy(nx.to_numpy_matrix(brain.G, nonedge=nonedge))


def assignbctresult(brain, bct_res):
    """ 
    It assigns the results of a Brain Connectivity Toolbox matrix to the respective nodes in 
    maybrain
    """
    return dict(list(zip(brain.G.nodes(), bct_res)))
