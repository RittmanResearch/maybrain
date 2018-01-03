import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from maybrain import brain as mbt
from maybrain import resources as rr
from maybrain import constants as ct

import numpy as np
import copy


def _get_ordered_array_and_labels(matrix, dummy_adj_file, hemi_prop, lobes_prop, anat_prop,
                                  hemi_name, lobes_name, anat_name):
    dummy_brain = mbt.Brain()
    dummy_brain.import_adj_file(dummy_adj_file)
    dummy_brain.apply_threshold()
    dummy_brain.import_properties(anat_prop)
    dummy_brain.import_properties(lobes_prop)
    dummy_brain.import_properties(hemi_prop)

    # Sorting the nodes, first by hemisphere, then by lobe, then by anatomical label
    nodes = copy.deepcopy(list(dummy_brain.G.nodes(data=True)))
    nodes.sort(key=lambda x: (x[1][hemi_name], x[1][lobes_name], x[1][anat_name]))

    # Getting the labels of the nodes and the numbers, now in ascending order
    labels = [x[1][anat_name] for x in nodes]
    permutation = [x[0] for x in nodes]

    # Copying adjMat, and ordering the rows and columns
    order = np.argsort(permutation)

    arr = np.copy(matrix)
    arr = arr[:, order]
    arr = arr[order, :]

    return arr, labels


def _plot_array(arr, title, labels, output_file):
    fig, ax = plt.subplots(figsize=(len(arr) / 50, len(arr) / 50))
    # Plotting
    im = ax.imshow(arr, interpolation='none')  # plotting the ordered adjacency matrix
    # setting the labels by anatlabels name instead of numbers
    ax.set_xticklabels(labels)
    ax.set_yticklabels(labels)
    ax.xaxis.tick_top()
    # Show the label every in every ticker
    ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(1))

    # label size should be small otherwise impossible to see
    plt.setp(ax.get_xticklabels(), visible=True)
    plt.setp(ax.get_yticklabels(), visible=True)
    ax.tick_params(axis='both', which='both', labelsize=1, length=0.2, width=0.1)

    ax.tick_params(axis='both', pad=0.9)
    ax.set_xticks(np.arange(len(arr) + 1))  # +1 creates an extra tick but otherwise plt would clip the image
    ax.set_yticks(np.arange(len(arr) + 1))
    # vertical labels
    for label in im.axes.xaxis.get_ticklabels():
        label.set_rotation(90)

    plt.colorbar(im)
    ax.set_title(title, y=1.08)

    # Removing outer lines because they hide part of the first line/column
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)

    # If outputfile is defined, fig is correctly closed, otherwise returned so others can add more information to it
    if output_file is not None:
        fig.savefig(output_file)
        plt.close(fig)
    else:
        return fig, ax


def plot_avg_matrix(brains, output_file=None, dummy_adj_file=rr.DUMMY_ADJ_FILE_500,
                    hemi_prop=rr.PROPERTIES_HEMISPHERES_500, hemi_name=ct.HEMISPHERE,
                    lobes_prop=rr.PROPERTIES_LOBES_500, lobes_name=ct.LOBE,
                    anat_prop=rr.PROPERTIES_ANATLABEL_500, anat_name=ct.ANAT_LABEL):
    """
    It uses matplotlib to plot the connection strength average matrix of all `adjMat`s in `brains`.
    The nodes are ordered first by hemisphere, then by lobes, and then by anatomical labels.

    For each edge, the final averaged connection strength ignores NaNs in the calculations.

    brains: a dictionary where the values are instances of `Brain`s. The keys don't matter.
    output_file: if you want to create a file. It then calls fig.savefig(output_file) from matplotlib
    dummy_adj_file: an adjacency matrix needed to create a dummy instance of Brain
    hemi_prop: hemisphere properties needed to create a dummy instance of Brain
    hemi_name: name to lookup for the hemisphere property
    lobes_prop: lobes properties needed to create a dummy instance of Brain
    lobes_name: name to lookup for the lobe property
    anat_prop: anatomical labels properties needed to create a dummy instance of Brain
    anat_name: name to lookup for the anatomical labels property
    """
    # With empty dictionary, nothing to do
    if not brains.values():
        raise TypeError("brains is empty, nothing can be done")
    # Check correct size of `brains`
    size_brain = None
    for brn in brains.values():
        if not size_brain:
            size_brain = brn.adjMat.shape[0]
        elif brn.adjMat.shape != (size_brain, size_brain):
            raise TypeError("The brains are not all with the same size")

    # The final matrix with the averaged connection strength
    avg_matrix = np.empty([size_brain, size_brain])
    # For each edge, it says how many brains contain that edge. Used to calculate the mean (we have to divide by the
    #  total number)
    elements = np.zeros([size_brain, size_brain], dtype=int)

    avg_matrix[:] = np.nan

    for brn in brains.values():
        for (x, y), value in np.ndenumerate(brn.adjMat):
            if np.isnan(value):
                continue

            if np.isnan(avg_matrix[x][y]):
                avg_matrix[x][y] = value
            else:
                avg_matrix[x][y] += value
            elements[x][y] += 1

    # Calculating the actual average
    for x in range(size_brain):
        for y in range(size_brain):
            if np.isnan(avg_matrix[x][y]):
                continue
            if elements[x][y] == 0:
                continue
            avg_matrix[x][y] /= elements[x][y]

    arr, labels = _get_ordered_array_and_labels(avg_matrix, dummy_adj_file=dummy_adj_file,
                                                hemi_prop=hemi_prop, lobes_prop=lobes_prop, anat_prop=anat_prop,
                                                hemi_name=hemi_name, lobes_name=lobes_name, anat_name=anat_name)

    return _plot_array(arr, title="Average Connection Strength Matrix", labels=labels, output_file=output_file)


def plot_strength_matrix(brain, title="", output_file=None, dummy_adj_file=rr.DUMMY_ADJ_FILE_500,
                         hemi_prop=rr.PROPERTIES_HEMISPHERES_500, hemi_name=ct.HEMISPHERE,
                         lobes_prop=rr.PROPERTIES_LOBES_500, lobes_name=ct.LOBE,
                         anat_prop=rr.PROPERTIES_ANATLABEL_500, anat_name=ct.ANAT_LABEL):
    """
    It uses matplotlib to plot the connection strength matrix of the `adjMat` of `brain`.
    The nodes are ordered first by hemisphere, then by lobes, and then by anatomical labels.

    brain: an instance of the `Brain` class
    title: title to appear on top of the matrix
    output_file: if you want to create a file. It then calls fig.savefig(output_file) from matplotlib
    dummy_adj_file: an adjacency matrix needed to create a dummy instance of Brain
    hemi_prop: hemisphere properties needed to create a dummy instance of Brain
    hemi_name: name to lookup for the hemisphere property
    lobes_prop: lobes properties needed to create a dummy instance of Brain
    lobes_name: name to lookup for the lobe property
    anat_prop: anatomical labels properties needed to create a dummy instance of Brain
    anat_name: name to lookup for the anatomical labels property
    """
    arr, labels = _get_ordered_array_and_labels(brain.adjMat, dummy_adj_file=dummy_adj_file,
                                                hemi_prop=hemi_prop, lobes_prop=lobes_prop, anat_prop=anat_prop,
                                                hemi_name=hemi_name, lobes_name=lobes_name, anat_name=anat_name)

    return _plot_array(arr, title=title, labels=labels, output_file=output_file)
