# -*- coding: utf-8 -*-

import numpy as np
from os import path, rename


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


def writeResults(results, measure,
                 outfilebase="brain",
                 append=True,
                 propDict=None):
    """
    Function to write out results
    """
    outFile = outfilebase + measure + '.txt'
    if not path.exists(outFile) or not append:
        if path.exists(outFile):
            rename(outFile, outFile + '.old')
        f = open(outFile, "w")
        writeHeadFlag = True

    else:
        f = open(outFile, "a")
        writeHeadFlag = False

    # check to see what form the results take
    if isinstance(results, (dict)):
        headers = [v for v in list(results.keys())]
        headers.sort()

        out = ' '.join([str(results[v]) if v in list(results.keys()) else 'NA' for v in headers])

    elif isinstance(results, (list)):
        if writeHeadFlag:
            headers = ' '.join([str(v) for v in range(len(results))])
        out = ' '.join([str(v) for v in results])

    else:
        if writeHeadFlag:
            headers = [measure]
        out = str(results)

    # write headers
    if propDict:
        propHeaders = list(propDict.keys())
        propHeaders.sort()

    if writeHeadFlag:
        headers = ' '.join([str(v) for v in headers])
        if propDict:
            headers = ' '.join([' '.join(propHeaders), headers])

        f.writelines(headers + '\n')

    # add on optional extras
    if propDict:
        outProps = ' '.join([str(propDict[v]) for v in propHeaders])
        out = ' '.join([outProps, out])

    f.writelines(out + '\n')
    f.close()


def extractCoordinates(template, outFile="ROI_xyz.txt"):
    import nibabel as nb
    import csv

    # load image
    f = nb.load(template)
    fData = f.get_data()

    # extract voxel coordinates
    I, J, K = fData.shape
    coords = np.zeros((I, J, K, 3))
    coords[..., 0] = np.arange(I)[:, np.newaxis, np.newaxis]
    coords[..., 1] = np.arange(J)[np.newaxis, :, np.newaxis]
    coords[..., 2] = np.arange(K)[np.newaxis, np.newaxis, :]

    # convert voxel coordinates to mm values using affine values
    M = f.affine[:3, :3]
    abc = f.affine[:3, 3]

    for x in range(len(coords[:, 0, 0, 0])):
        for y in range(len(coords[0, :, 0, 0])):
            for z in range(len(coords[0, 0, :, 0])):
                coords[x, y, z, :] = M.dot(coords[x, y, z, :]) + abc

    # get unique values in the mask
    valList = np.unique(fData)
    valList = valList[valList != 0]

    out = open(outFile, "w")
    writer = csv.DictWriter(out,
                            fieldnames=["Node", "x", "y", "z"],
                            delimiter=" "
                            )
    for v in valList:
        tempArr = np.zeros((I, J, K, 3), dtype=bool)
        tfArray = fData == v
        tempArr[..., 0] = tfArray
        tempArr[..., 1] = tfArray
        tempArr[..., 2] = tfArray

        tempArr = coords[tempArr]
        tempArr = np.mean(tempArr.reshape([tempArr.shape[0] / 3, 3]), axis=0)
        outList = [str(int(v))]
        outList.extend(["{:.2f}".format(x) for x in tempArr])
        outDict = dict(list(zip(writer.fieldnames,
                                outList)))
        writer.writerow(outDict)

    out.close()
