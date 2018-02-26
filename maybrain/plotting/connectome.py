import matplotlib.pyplot as plt
import numpy as np
import colorsys
import networkx as nx
from nilearn import plotting
from maybrain import constants as ct


def plot_connectome(brain, only_nodes=False, node_property=None, node_attributes=None, **kwargs):
    """ 
    Wrapper over `nilearn.plotting.plot_connectome` to plot the connectome of a brain.

    Brain's nodes should have `constants.XYZ` attribute (spatial information)

    Parameters
    ----------
    brain: maybrain.brain.Brain
        An instance of the `Brain` class
    only_nodes: bool
        If True, only nodes will be plotted (no edges)
    kwargs
        Keyword arguments if you need to pass them to nilearn's plot_connectome()
if node_property is, node_color is overwritten
    Returns
    -------
    display
        the display from nilearn (return from plot_connectome)

    Raises
    ------
    KeyError: Exception
        If the edges don't have ct.XYZ property
    """
    try:
        if list(brain.G.nodes(data=True))[0][1][ct.XYZ]:
            pass
    except KeyError as error:
        import sys
        _, _, tb = sys.exc_info()
        raise KeyError(error, "Node doesn't have constants.XYZ property").with_traceback(tb)

    if 'edge_cmap' not in kwargs:
        kwargs['edge_cmap'] = plt.get_cmap('YlGnBu')
    if 'node_color' not in kwargs:
        kwargs['node_color'] = 'red'

    connection_matrix = np.copy(nx.to_numpy_matrix(brain.G, nonedge=0))
    if only_nodes:
        connection_matrix[:] = 0

    # If node_property is defined, let's create the custom colours
    if node_property:
        palette = []

        # Creating a palette of colours based on the amount of
        for i in range(len(node_attributes) + 1):  # +1 to account for "other" regions
            rgb = colorsys.hsv_to_rgb(i / (len(node_attributes) + 1), 1.0, 1.0)
            palette.append('#%02x%02x%02x' % tuple([(round(255*x)) for x in rgb]))

        # Defining the colour for each node according to the attribute it has
        colours = []
        for n in brain.G.nodes(data=True):
            we = np.where(np.array(node_attributes) == n[1][node_property])
            if len(we[0]) > 0:
                colours.append(palette[we[0][0]])
            else:
                colours.append(palette[-1])  # not found, so another colour

        kwargs['node_color'] = colours

    return plotting.plot_connectome(connection_matrix,
                                    list(dict(brain.G.nodes(data=ct.XYZ)).values()),
                                    **kwargs)
