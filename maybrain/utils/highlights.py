from maybrain import constants as ct
import numpy as np


__number_of_highlights = 1  # Used for automatic labeling
highlights = {}  # Our highlights


def _get_auto_label():
    """ generate an automatic label for a highlight object """
    global __number_of_highlights, highlights

    # If number doesn't exist in highlights, use it
    if __number_of_highlights not in highlights:
        return __number_of_highlights
    # Otherwise, try until a key is available
    else:
        while __number_of_highlights in highlights:
            __number_of_highlights += 1
        return __number_of_highlights


class Highlight:
    """
    It holds information to highlight from a subsection of a Brain
    """
    def __init__(self, nodes=None, edges=None):
        if edges is None:
            edges = []
        if nodes is None:
            nodes = []

        self.nodes = nodes
        self.edges = edges


def make_highlight(edge_inds=None, nodes_inds=None, label=None):
    """
    It create a highlight object from some edge/nodes indices, and put it in the highlights global

    edge_inds: a list of tuples of edges
    nodes_inds: a list with nodes
    label: the name given to the highlight. If non is provided, it will be automatically generated
    """
    global highlights

    if not edge_inds and not nodes_inds:
        raise TypeError("Your are trying to create an empty highlight")
    if not all(isinstance(item, tuple) and len(item) == 2 for item in edge_inds):
        raise TypeError("The list of edges is not a list of tuples")

    h = Highlight(edges=edge_inds, nodes=nodes_inds)

    if not label:
        label = _get_auto_label()

    highlights[label] = h


def highlight_from_conds(brain, prop, rel, val, mode, label=None):
    """
    Creates a highlight by asking if the property `prop` is related to `val` by `rel`.

    brain: An instance of Brain to highlight
    prop: The property to look for in `brain`. If this value is equal to x/X/y/Y/z/Z, it will look out for the
          respective value from the property ct.XYZ
    rel: The relation between the property and value that is being looked for in the brain. It can be:
        'geq' - greater than or equal to
        'leq' - less than or equal to
        'gt' - strictly greater than
        'lt' - strictly less than
        'eq' - equal to (i.e. exactly)
        'in()', 'in[)', 'in(]', 'in[]' - within an interval, in this case val is a list of two numbers
                                         "[" and "]" means inclusive, "(" and ")" means exclusive
        'in' - in `val`
    val: The value, or range of values, to be compared to
    mode: Whether looking for nodes and/or edges. It can be 'node', 'edge', 'node|edge' or 'edge|node'
    label: Identification of the highlight for further access. If None, one will be automatically created
    """
    global highlights

    # check filter mode
    if not (mode in ['edge', 'node', 'node|edge', 'edge|node']):
        raise TypeError('Filter mode is not recognised')
    # check if label given
    if not label:
        label = _get_auto_label()

    # make an instance of Highlight
    h = Highlight()

    # extract lists from edges
    if mode in ['edge', 'node|edge', 'edge|node']:
        for e in brain.G.edges(data=True):
            try:
                val_to_compare = e[2][prop]
            except KeyError:
                continue  # property doesn't exist, continuing

            # add to highlight if good
            if _prop_compare(val_to_compare, rel, val):
                h.edges.append((e[0], e[1]))

    # extract lists from nodes
    if mode in ['node', 'node|edge', 'edge|node']:
        for node in brain.G.nodes(data=True):
            try:
                # special treatment for 'x', 'y' and 'z'
                if prop in ['x', 'X']:
                    val_to_compare = node[1][ct.XYZ][0]
                elif prop in ['y', 'Y']:
                    val_to_compare = node[1][ct.XYZ][1]
                elif prop in ['z', 'Z']:
                    val_to_compare = node[1][ct.XYZ][2]
                else:
                    # any other property
                    val_to_compare = node[1][prop]
            except KeyError:
                continue  # property doesn't exist, continuing

            # add to highlight if good
            if _prop_compare(val_to_compare, rel, val):
                h.nodes.append(node[0])

    # add highlight to dictionary
    highlights[label] = h


def _prop_compare(val_to_compare, rel, val):
    """
    Compare `val_to_compare` relative to `val`
    This is meant to be used just by highlight_from_conds()
    """
    d = val_to_compare  # making code more readily down here
    if rel == 'eq':
        b = d == val
    elif rel == 'gt':
        b = d > val
    elif rel == 'lt':
        b = d < val
    elif rel == 'leq':
        b = d <= val
    elif rel == 'geq':
        b = d >= val
    elif rel == 'in()':
        b = (d > val[0]) and (d < val[1])
    elif rel == 'in[)':
        b = (d >= val[0]) and (d < val[1])
    elif rel == 'in(]':
        b = (d > val[0]) and (d <= val[1])
    elif rel == 'in[]':
        b = (d >= val[0]) and (d <= val[1])
    elif rel == 'in':
        b = d in val
    else:
        raise TypeError('Relation not recognised: ' + str(rel))

    return b
