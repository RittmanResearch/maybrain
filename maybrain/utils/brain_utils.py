# -*- coding: utf-8 -*-
"""
Utility module with recipes to calculate some common measures
"""
import csv
from os import path, rename
import nibabel as nb
import numpy as np


def threshold_to_percentage(brain, threshold):
    """
    It returns a ratio between the edges on adjMat above a certain threshold value \
    and the total possible edges of adjMat.
    In an unidrected graph, the total possible edges are the upper right \
    part elements of adjMat different from np.nan. In a directed graph, the total \
    possible edges are all the elements of adjMat except the diagonal and \
    np.nan

    Parameters
    ----------
    brain: maybrain.brain.Brain
        An instance of the `Brain` class
    threshold: number
        The threshold value

    Returns
    -------
    ratio: float
        The final result
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

    Parameters
    ----------
    brain: maybrain.brain.Brain
        An instance of the `Brain` class
    """
    nodes = brain.G.number_of_nodes()

    if brain.directed:
        total_connections = nodes * (nodes - 1)
    else:
        total_connections = nodes * (nodes - 1) / 2
    return float(brain.G.number_of_edges()) / float(total_connections)


def write_results(results, measure,
                  outfilebase="brain",
                  append=True,
                  propdict=None):
    """
    Function to write out results
    """
    outfile = outfilebase + measure + '.txt'
    if not path.exists(outfile) or not append:
        if path.exists(outfile):
            rename(outfile, outfile + '.old')
        file = open(outfile, "w")
        writeheadflag = True

    else:
        file = open(outfile, "a")
        writeheadflag = False

    # check to see what form the results take
    if isinstance(results, dict):
        headers = list(results.keys())
        headers.sort()

        out = ' '.join([str(results[v]) if v in list(results.keys()) else 'NA' for v in headers])

    elif isinstance(results, (list)):
        if writeheadflag:
            headers = ' '.join([str(v) for v in range(len(results))])
        out = ' '.join([str(v) for v in results])

    else:
        if writeheadflag:
            headers = [measure]
        out = str(results)

    # write headers
    if propdict:
        propheaders = list(propdict.keys())
        propheaders.sort()

    if writeheadflag:
        headers = ' '.join([str(v) for v in headers])
        if propdict:
            headers = ' '.join([' '.join(propheaders), headers])

        file.writelines(headers + '\n')

    # add on optional extras
    if propdict:
        outprops = ' '.join([str(propdict[v]) for v in propheaders])
        out = ' '.join([outprops, out])

    file.writelines(out + '\n')
    file.close()


def extract_coordinates(template, outfile="ROI_xyz.txt"):
    """
    Legacy code
    """
    # load image
    file = nb.load(template)
    f_data = file.get_data()

    # extract voxel coordinates
    i, j, k = f_data.shape
    coords = np.zeros((i, j, k, 3))
    coords[..., 0] = np.arange(i)[:, np.newaxis, np.newaxis]
    coords[..., 1] = np.arange(j)[np.newaxis, :, np.newaxis]
    coords[..., 2] = np.arange(k)[np.newaxis, np.newaxis, :]

    # convert voxel coordinates to mm values using affine values
    mmm = file.affine[:3, :3]
    abc = file.affine[:3, 3]

    for x in range(len(coords[:, 0, 0, 0])):
        for y in range(len(coords[0, :, 0, 0])):
            for z in range(len(coords[0, 0, :, 0])):
                coords[x, y, z, :] = mmm.dot(coords[x, y, z, :]) + abc

    # get unique values in the mask
    val_list = np.unique(f_data)
    val_list = val_list[val_list != 0]

    out = open(outfile, "w")
    writer = csv.DictWriter(out,
                            fieldnames=["Node", "x", "y", "z"],
                            delimiter=" ")
    for val in val_list:
        temp_arr = np.zeros((i, j, k, 3), dtype=bool)
        tf_array = f_data == val
        temp_arr[..., 0] = tf_array
        temp_arr[..., 1] = tf_array
        temp_arr[..., 2] = tf_array

        temp_arr = coords[temp_arr]
        temp_arr = np.mean(temp_arr.reshape([temp_arr.shape[0] / 3, 3]), axis=0)
        outlist = [str(int(val))]
        outlist.extend(["{:.2f}".format(x) for x in temp_arr])
        outdict = dict(list(zip(writer.fieldnames,
                                outlist)))
        writer.writerow(outdict)

    out.close()
