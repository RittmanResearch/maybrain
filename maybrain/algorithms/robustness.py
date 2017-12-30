import random
import networkx as nx
import numpy as np


def robustness(G, iterLen=500, N=50):
    """
    A function to calculate robustness based on "Error and attack
    tolerance of complex networks" Albert et al Nature 2000 406:378-382

    The function calculates the rate of change in the size of the largest
    connected component (S) as nodes are randomly removed. The process runs
    iteratively and takes a mean. The gradient of S is smoothed to provide
    a more accurate measure by a sliding window.

    N = size of the sliding window for smoothing the gradient
    iterLen = number of iterations

    Note, this function is relatively slow compared to other metrics due to
    the multiple iterations.

    """
    fList = np.zeros(iterLen)
    for i in range(iterLen):
        a = sorted(nx.connected.connected_component_subgraphs(G), key=len, reverse=True)[0]
        nList = [v for v in G.nodes()]
        random.shuffle(nList)
        nList = nList[:-1]
        mat = np.zeros((len(nList)))
        count = 0
        S = len(a.nodes())

        while nList and S > 1:
            n = nList.pop()
            if n in a.nodes():
                a.remove_node(n)
                if not nx.is_connected(a):  # only recalculate if the further fragmentation
                    a = sorted(nx.connected.connected_component_subgraphs(a), key=len, reverse=True)[0]
                    S = len(a.nodes())
                else:
                    S -= 1
            mat[count] = S
            count += 1

        g = np.gradient(mat)
        runMean = np.convolve(g, np.ones((N,)) / N)[(N - 1):]
        diffs = np.diff(runMean)
        nr = np.argmin(diffs)

        fList[i] = nr
    return np.mean(fList) / len(G.nodes())
