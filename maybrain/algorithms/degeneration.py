# -*- coding: utf-8 -*-

"""
Functions for degeneration of a brain object
"""
import random

import networkx as nx
import numpy as np

from maybrain import constants as ct


def degenerate(brain, weight_loss=0.1, edges_removed_limit=1,
               thresh_limit=None, pc_limit=None, weight_loss_limit=None,
               node_list=[], spread=False, update_adj_mat=True,
               distances=False, spatial_search=False):
    """
    Remove random edges from connections of the toxicNodes set, or from the
    risk_edges set. This occurs either until edgesRemovedLimit number of edges
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

    """

    # get rid of the existing list of edges if the node list is specified
    if node_list:
        brain.risk_edges = None

    # set limit
    if weight_loss_limit and pc_limit:
        print("You have asked for both a weight and percentage"
              "connectivitylimit, using the percentage connectivity limit")

    if thresh_limit:
        pc_limit = brain.threshold_to_percentage(thresh_limit)

    if pc_limit:
        len_nodes = len(brain.G.nodes())
        len_edges = len(brain.G.edges())

        max_edges = float(len_nodes * (len_nodes - 1))
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

    if not brain.risk_edges:
        redefine_edges = True
        # if no toxic nodes defined, select the whole graph
        if not node_list:
            node_list = brain.G.nodes()

        # generate list of at risk edges
        brain.risk_edges = [v for v in nx.edges(brain.G, node_list)
                            if brain.G.edge[v[0]][v[1]][ct.WEIGHT] != 0.]
    else:
        redefine_edges = False

    if spread:
        node_list = []

    # check if there are enough weights left
    risk_edge_wt_sum = np.sum([brain.G.edge[v[0]][v[1]][ct.WEIGHT]
                               for v in brain.risk_edges])
    if limit > risk_edge_wt_sum:
        print("Not enough weight left to remove")
        return node_list

    while limit > 0.:
        if not brain.risk_edges and spatial_search:
            # find spatially closest nodes if no edges exist
            # is it necessary to do this for all nodes?? - waste of computing power
            # choose node first, then calculated spatially nearest of a single node
            new_node = brain.find_spatially_nearest(node_list)
            if new_node:
                print("Found spatially nearest node")
                node_list.append(new_node)
                brain.risk_edges = nx.edges(brain.G, node_list)
            else:
                print("No further edges to degenerate")
                break
        # choose at risk edge to degenerate from
        dying_edge = random.choice(brain.risk_edges)

        # remove specified weight from edge
        wei = brain.G[dying_edge[0]][dying_edge[1]][ct.WEIGHT]

        if np.absolute(wei) < weight_loss:
            loss = np.absolute(wei)
            brain.G.remove_edge(dying_edge[0], dying_edge[1])
            brain.risk_edges.remove(dying_edge)
            if not weight_loss_limit:
                limit -= 1

        elif wei > 0:
            loss = weight_loss
            brain.G[dying_edge[0]][dying_edge[1]][ct.WEIGHT] -= weight_loss

        else:
            loss = weight_loss
            brain.G[dying_edge[0]][dying_edge[1]][ct.WEIGHT] += weight_loss

        # record the edge length of edges lost
        if distances:
            brain.dying_edges[dying_edge] = brain.G[dying_edge[0]][dying_edge[1]]
            brain.dying_edges[dying_edge][ct.DISTANCE] = \
                np.linalg.norm(np.array((brain.G.node[dying_edge[0]][ct.XYZ]))
                               - np.array((brain.G.node[dying_edge[1]][ct.XYZ])))

        # update the adjacency matrix (essential if robustness is to be calculated)
        if update_adj_mat:
            brain.update_adj_mat(dying_edge)

        # add nodes to toxic list if the spread option is selected
        if spread:
            for node in dying_edge:
                if not (node in node_list or not (brain.G.edge[dying_edge[0]][dying_edge[1]] > spread)):
                    node_list.append(node)

        if weight_loss_limit:
            limit -= loss

        # redefine at risk edges
        if redefine_edges or spread:
            brain.risk_edges = nx.edges(brain.G, node_list)

    print("Number of toxic nodes: " + str(len(node_list)))

    return node_list


def contiguous_spread(brain, edgeloss, startnodes=None):
    """
    Degenerate nodes in a continuous fashion.
    Doesn't currently include spreadratio
    """

    # make sure nodes have the linkedNodes attribute
    try:
        brain.G.node[0][ct.LINKED_NODES]
    except:
        brain.find_linked_nodes()

    # make sure all nodes have degenerating attribute
    try:
        brain.G.node[0]['degenerating']
    except:
        for n in range(len(brain.G.nodes())):
            brain.G.node[n]['degenerating'] = False

    # start with a random node or set of nodes
    if not startnodes:
        # start with one random node if none chosen
        toxic_nodes = [random.randint(0, len(brain.G.nodes()))]
    else:
        # otherwise use user provided nodes
        toxic_nodes = startnodes
    # make all toxic nodes degenerating
    for t in toxic_nodes:
        brain.G.node[t]['degenerating'] = True

    # put at-risk nodes into a list
    risk_nodes = []
    for t in toxic_nodes:
        l = brain.G.node[t][ct.LINKED_NODES]
        newl = []
        # check the new indices aren't already toxic
        for a in l:
            if a in toxic_nodes:
                continue
            if brain.G.node[a]['degenerating']:
                continue
            #                if not(a in toxic_nodes)&(not(brain.G.node[a]['degenerating'])):
            newl.append(a)

        risk_nodes = risk_nodes + newl

    # iterate number of steps
    toxic_node_record = [toxic_nodes[:]]
    for _ in range(edgeloss):
        # find at risk nodes
        ind = random.randint(0, len(risk_nodes) - 1)
        # get the index of the node to be removed and remove from list
        dead_node = risk_nodes.pop(ind)
        # remove all instances from list
        while dead_node in risk_nodes:
            risk_nodes.remove(dead_node)

        # add to toxic list
        toxic_nodes.append(dead_node)
        # make it degenerate
        brain.G.node[dead_node]['degenerating'] = True
        print(('dead_node', dead_node))

        # add the new at-risk nodes
        l = brain.G.node[dead_node][ct.LINKED_NODES]
        newl = []
        # check the new indices aren't already toxic
        for a in l:
            if a in toxic_nodes:
                continue
            if brain.G.node[a]['degenerating']:
                continue
            newl.append(a)

        risk_nodes = risk_nodes + newl

        toxic_node_record.append(toxic_nodes[:])

        # check that there are any more nodes at risk
        if len(risk_nodes) == 0:
            break

        #            print(toxic_nodes)

    # Update adjacency matrix to reflect changes
    brain.reconstruct_adj_mat()

    return toxic_nodes, toxic_node_record
