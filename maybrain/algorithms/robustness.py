"""
Module to calculate robustness of a brain's graph
"""
import random
import networkx as nx
import numpy as np


def robustness(brain, iter_len=500, window_size=50):
    """
    A function to calculate robustness based on "Error and attack
    tolerance of complex networks" Albert et al. Nature 2000 406:378-382

    The function calculates the rate of change in the size of the largest
    connected component (S) as nodes are randomly removed. The process runs
    iteratively and takes a mean. The gradient of S is smoothed to provide
    a more accurate measure by a sliding window.

    Note, this function is relatively slow compared to other metrics due to
    the multiple iterations.

    Parameters
    ----------
    brain: maybrain.brain.Brain
        An instance of the `Brain` class
    iter_len: int
        number of iterations
    window_size: int
        size of the sliding window for smoothing the gradient

    """
    f_list = np.zeros(iter_len)
    for i in range(iter_len):
        a_max = max(nx.connected.connected_component_subgraphs(brain.G), key=len)
        n_list = [v for v in brain.G.nodes()]
        random.shuffle(n_list)
        n_list = n_list[:-1]
        mat = np.zeros(len(n_list))
        count = 0
        a_len = a_max.number_of_nodes()

        while n_list and a_len > 1:
            n = n_list.pop()
            if n in a_max.nodes():
                a_max.remove_node(n)
                if not nx.is_connected(a_max):  # only recalculate if the further fragmentation
                    a_max = max(nx.connected.connected_component_subgraphs(a_max), key=len)
                    a_len = a_max.number_of_nodes()
                else:
                    a_len -= 1
            mat[count] = a_len
            count += 1

        grad = np.gradient(mat)
        run_mean = np.convolve(grad, np.ones((window_size,)) / window_size)[(window_size - 1):]
        diffs = np.diff(run_mean)
        nr_arg = np.argmin(diffs)

        f_list[i] = nr_arg
    return np.mean(f_list) / brain.G.number_of_nodes()
