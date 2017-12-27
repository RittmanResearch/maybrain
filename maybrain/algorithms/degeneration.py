# -*- coding: utf-8 -*-

"""
Functions for degeneration of a brain object
"""

import networkx as nx
import numpy as np
import random
from maybrain import constants as ct

def degenerate(brain, weight_loss=0.1, edges_removed_limit=1,
               thresh_limit=None, pc_limit=None, weight_loss_limit=None,
               node_list=[], risk_edges=None, spread=False,
               update_adj_mat=True, distances=False, spatial_search=False):
    '''
    Remove random edges from connections of the toxicNodes set, or from the
    riskEdges set. This occurs either until edgesRemovedLimit number of edges
    have been removed (use this for a thresholded weighted graph), or until the
    weight loss limit has been reached (for a weighted graph). For a binary
    graph, weight loss should be set to 1.

    The spread option recruits connected nodes of degenerating edges to the
    toxic nodes list.

    By default this function will enact a random attack model, with a weight
    loss of 0.1 each iteration.

    Weights are taken as absolute values, so the weight in any affected edge
    tends to 0.

    Spread can either be False, or a number specifying the weight above which
    to add nodes within the list of at-risk nodes.

    '''

    # get rid of the existing list of edges if the node list is specified
    if node_list:
        brain.riskEdges = None

    # set limit
    if weight_loss_limit and pc_limit:
        print("You have asked for both a weight and percentage"
              "connectivitylimit, using the percentage connectivity limit")

    if thresh_limit:
        pc_limit = brain.thresholdToPercentage(thresh_limit)

    if pc_limit:
        len_nodes = len(brain.G.nodes())
        len_edges = len(brain.G.edges())

        max_edges = float(len_nodes * (len_nodes-1))
        if not brain.G.is_directed():
            max_edges = max_edges / 2

        new_edge_num = int(round(pc_limit * max_edges))
        if new_edge_num > len_edges:
            print("The percentage threshold set is lower than the current"
                  "graph, please choose a larger value")

        limit = len_edges - new_edge_num
        weight_loss_limit = False

    elif weight_loss_limit:
        limit = weight_loss_limit

    else:
        limit = edges_removed_limit

    if not brain.riskEdges:
        reDefineEdges = True
        # if no toxic nodes defined, select the whole graph
        if not node_list:
            nodeList = brain.G.nodes()

        # generate list of at risk edges
        brain.riskEdges = [v for v in nx.edges(brain.G, nodeList)
                           if brain.G.edge[v[0]][v[1]][ct.WEIGHT] != 0.]
    else:
        reDefineEdges = False

    if spread:
        nodeList = []

    # check if there are enough weights left
    riskEdgeWtSum = np.sum([brain.G.edge[v[0]][v[1]][ct.WEIGHT]
                            for v in brain.riskEdges])
    if limit > riskEdgeWtSum:
        print("Not enough weight left to remove")
        return nodeList

    while limit > 0.:
        if not brain.riskEdges and spatial_search:
        # find spatially closest nodes if no edges exist
        # is it necessary to do this for all nodes?? - waste of computing power
        # choose node first, then calculated spatially nearest of a single node
            newNode = brain.findSpatiallyNearest(nodeList)
            if newNode:
                print("Found spatially nearest node")
                nodeList.append(newNode)
                brain.riskEdges = nx.edges(brain.G, nodeList)
            else:
                print("No further edges to degenerate")
                break
        # choose at risk edge to degenerate from
        dyingEdge = random.choice(brain.riskEdges)

        # remove specified weight from edge
        w = brain.G[dyingEdge[0]][dyingEdge[1]][ct.WEIGHT]

        if np.absolute(w) < weight_loss:
            loss = np.absolute(w)
            brain.G.remove_edge(dyingEdge[0], dyingEdge[1])
            brain.riskEdges.remove(dyingEdge)
            if not weight_loss_limit:
                limit -= 1

        elif w > 0:
            loss = weight_loss
            brain.G[dyingEdge[0]][dyingEdge[1]][ct.WEIGHT] -= weight_loss

        else:
            loss = weight_loss
            brain.G[dyingEdge[0]][dyingEdge[1]][ct.WEIGHT] += weight_loss

        # record the edge length of edges lost
        if distances:
            brain.dyingEdges[dyingEdge] = brain.G[dyingEdge[0]][dyingEdge[1]]
            brain.dyingEdges[dyingEdge][ct.DISTANCE] = \
                np.linalg.norm(np.array((brain.G.node[dyingEdge[0]][ct.XYZ]))
                               - np.array((brain.G.node[dyingEdge[1]][ct.XYZ])))

        # update the adjacency matrix (essential if robustness is to be calculated)
        if update_adj_mat:
            brain.updateAdjMat(dyingEdge)

        # add nodes to toxic list if the spread option is selected
        if spread:
            for node in dyingEdge:
                if (not node in nodeList and
                        brain.G.edge[dyingEdge[0]][dyingEdge[1]] > spread):
                    nodeList.append(node)

        if weight_loss_limit:
            limit -= loss

        # redefine at risk edges
        if reDefineEdges or spread:
            brain.riskEdges = nx.edges(brain.G, nodeList)

    print("Number of toxic nodes: "+str(len(nodeList)))

    return nodeList

def contiguousSpread(brain, edgeloss, startNodes=None):
    '''
    Degenerate nodes in a continuous fashion.
    Doesn't currently include spreadratio
    '''

    # make sure nodes have the linkedNodes attribute
    try:
        brain.G.node[0][ct.LINKED_NODES]
    except:
        brain.findLinkedNodes()

    # make sure all nodes have degenerating attribute
    try:
        brain.G.node[0]['degenerating']
    except:
        for n in range(len(brain.G.nodes())):
            brain.G.node[n]['degenerating'] = False

    # start with a random node or set of nodes
    if not startNodes:
        # start with one random node if none chosen
        toxicNodes = [random.randint(0, len(brain.G.nodes()))]
    else:
        # otherwise use user provided nodes
        toxicNodes = startNodes
    # make all toxic nodes degenerating
    for t in toxicNodes:
        brain.G.node[t]['degenerating'] = True

    # put at-risk nodes into a list
    riskNodes = []
    for t in toxicNodes:
        l = brain.G.node[t][ct.LINKED_NODES]
        newl = []
        # check the new indices aren't already toxic
        for a in l:
            if a in toxicNodes:
                continue
            if brain.G.node[a]['degenerating']:
                continue
#                if not(a in toxicNodes)&(not(brain.G.node[a]['degenerating'])):
            newl.append(a)

        riskNodes = riskNodes + newl

    # iterate number of steps
    toxicNodeRecord = [toxicNodes[:]]
    for count in range(edgeloss):
        # find at risk nodes
        ind = random.randint(0, len(riskNodes)-1)
        # get the index of the node to be removed and remove from list
        deadNode = riskNodes.pop(ind)
        # remove all instances from list
        while deadNode in riskNodes:
            riskNodes.remove(deadNode)

        # add to toxic list
        toxicNodes.append(deadNode)
        # make it degenerate
        brain.G.node[deadNode]['degenerating'] = True
        print(('deadNode', deadNode))


        # add the new at-risk nodes
        l = brain.G.node[deadNode][ct.LINKED_NODES]
        newl = []
        # check the new indices aren't already toxic
        for a in l:
            if a in toxicNodes:
                continue
            if brain.G.node[a]['degenerating']:
                continue
            newl.append(a)

        riskNodes = riskNodes + newl

        toxicNodeRecord.append(toxicNodes[:])

        # check that there are any more nodes at risk
        if len(riskNodes) == 0:
            break

#            print(toxicNodes)

    # Update adjacency matrix to reflect changes
    brain.reconstructAdjMat()

    return toxicNodes, toxicNodeRecord
