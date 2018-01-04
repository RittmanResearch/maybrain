# -*- coding: utf-8 -*-
"""
Module which contains the definition of Brain class.
"""
import networkx as nx
import numpy as np
import random
from maybrain import constants as ct


class Brain:
    """
    A class that defines a brain network created from an adjacency matrix and
    spatial information with certain properties
    """

    def __init__(self, directed=False):
        """
        Initialise the brain object.
        """

        # is this a directed graph or not?
        self.directed = directed
        # adjacency matrix, containing weighting of edges. Should be square.
        self.adjMat = None
        # identification of the subject to which this brain object belongs
        self.subject = None
        # information about the scan which generated this brain object
        self.scan = None

        # Used for degeneration module
        self.dying_edges = {}

        # background info imported by nibabel
        self.nbbackground = None  # the nibabel object
        self.background = None  # coordinates of background
        self.backgroundHeader = None  # header of nibabel data

        # isosurface information imported by nibabel
        self.nbiso = None  # all the isosurface info, nibabel object
        self.iso = None  # the isosurface
        self.isoHeader = None  # header information

        self.risk_edges = None

        # For the properties features
        self.nodeProperties = []
        self.edgeProperties = []
        self.update_properties_after_threshold = False

        # create a new networkX graph object
        if self.directed:
            self.G = nx.DiGraph()
        else:
            self.G = nx.Graph()

    def import_adj_file(self, fname, delimiter=None, nodes_to_exclude=[], na_vals=["NA"]):
        """
        Imports an adjacency matrix from a file.
        fname : file name
        delimiter : the delimiter of the values inside the matrix, like ","
        nodes_to_exclude : Nodes you don't want to load, in an array format (nodes no. starts from zero)
        naVals : How the "Not a Number" values are represented in the file
        """
        lines = []
        try:
            with open(fname, "r") as f:
                for l in f:
                    l = l.strip()
                    lines.append(list(map(float, [v if v not in na_vals else np.nan for v in l.split(sep=delimiter)])))
        except IOError as error:
            error.strerror = 'Problem with opening file "' + fname + '": ' + error.strerror
            raise error

        # set adjacency matrix
        self.adjMat = np.array(lines)

        # add nodes
        self.G.add_nodes_from([v for v in range(len(lines)) if v not in nodes_to_exclude])

        # update adjacency matrix to null values of excluded nodes
        if nodes_to_exclude:
            for en in nodes_to_exclude:
                self.adjMat[:, en] = np.nan
                self.adjMat[en, :] = np.nan

    def import_spatial_info(self, fname, delimiter=None, convert_mni=False):
        """
        Add 3D coordinate information for each node from a given file. It needs to be called after importAdjFile()
        fname : file name
        delimiter : the delimiter of the values inside the matrix, like ","
        convert_mni : Whether you want to convert coordinates from voxel-wise (2 mm) to MNI space
        """
        # open file
        try:
            f = open(fname, "r")
        except IOError as error:
            error.strerror = 'Problem with opening 3D position file "' + fname + '": ' + error.strerror
            raise error

        # get data from file
        lines = f.readlines()
        node_count = 0
        for line in lines:
            l = line.strip().split(sep=delimiter)

            if convert_mni:
                l[1] = 45 - (float(l[1]) / 2)
                l[2] = 63 + (float(l[2]) / 2)
                l[3] = 36 + (float(l[3]) / 2)

            if node_count in self.G.nodes():
                self.G.nodes[node_count][ct.XYZ] = (float(l[1]), float(l[2]), float(l[3]))
                self.G.nodes[node_count][ct.ANAT_LABEL] = l[0]

            node_count += 1

        f.close()

    def import_properties(self, filename):
        """
        Add properties from a file. First line should contain the property name and the following \
        lines spaced node indices and property value e.g.:

            colour
            1 red
            2 white
            2 4 green

        Not that if 2 indices are given, the property is applied to edges instead.
        Properties will be treated as strings.
        You can mix nodes and edges in the same file.
        """
        edges_p = []
        nodes_p = []
        with open(filename, 'r') as file:
            prop = file.readline().strip()
            for line in file:
                info = line.strip().split(' ')

                if len(info) == 2:
                    nodes_p.append([prop, int(info[0]), info[1]])
                elif len(info) == 3:
                    edges_p.append([prop, int(info[0]), int(info[1]), info[2]])
                else:
                    raise ValueError(
                        'Problem in parsing %s, it has an invalid structure at line %d' % (filename, file.tell()))

        self._add_properties(edges_p)
        self._add_properties(nodes_p)

        self.nodeProperties.extend(nodes_p)
        self.edgeProperties.extend(edges_p)

    def _add_properties(self, properties):
        """
        It receives a list with properties and add them to either the nodes or edges according to the structure.
        For nodes properties, the format is:
            [ [property_name, node_id, property_value], [property_name, node_id, property_value], ...]

        For edges properties, the format is:
            [ [property_name, edge1, edge2, property_value], [property_name, edge1, edge2, property_value], ...]
        """
        for prop in properties:
            try:
                if len(prop) == 3:  # nodes
                    self.G.nodes[prop[1]][prop[0]] = prop[2]
                elif len(prop) == 4:  # edges
                    self.G.edges[prop[1], prop[2]][prop[0]] = prop[3]
            except:
                print('Warning! Unable to process property %s' % prop)

    def import_background(self, fname):
        import nibabel as nb
        ''' Import a file for background info using nbbabel
            gives a 3D array with data range 0 to 255 for test data
            could be 4d??
            defines an nibabel object, plus ndarrays with data and header info in

        '''

        self.nbbackground = nb.load(fname)
        self.background = self.nbbackground.get_data()
        self.backgroundHeader = self.nbbackground.get_header()

    def import_iso(self, fname):
        """
        Imports an isosurface info using nibabel
        gives a 3D array with data range 0 to 255 for test data
        defines an nibabel object, plus ndarrays with data and header info in

        fname: File Name with the isosurface
        """
        import nibabel as nb
        self.nbiso = nb.load(fname)
        self.iso = self.nbiso.get_data()
        self.isoHeader = self.nbiso.get_header()

    def parcels(self, node_list):
        """
        Plots 3D parcels specified in the node_list. This function assumes the
        parcellation template has been loaded to the brain using brain.importISO.
        Note, values passed to this function should correspond with those in the
        iso image, not necessarily the node values.
        """

        zero_arr = np.zeros(self.iso.shape)

        for n in node_list:
            n = float(n)
            # parcel files start from 1, zero is for background
            n_arr = np.ma.masked_where(self.iso != (n + 1), self.iso)
            n_arr.fill_value = 0.0
            zero_arr = zero_arr + n_arr
            zero_arr.mask = None

        self.parcelList = np.ma.masked_values(zero_arr, 0.0)

    def export_parcels_nii(self, outname='brain', value_dict=None):
        """
        This function saves the parcelList as a nifti file. It requires the
        brain.parcels function has been run first.
        """
        import nibabel as nb
        if value_dict:  # creates a numpy array based on the dictionary provided
            out_mat = np.zeros(self.nbiso.get_data().shape, dtype="float64")
            for n in list(value_dict.keys()):
                out_mat[np.where(self.nbiso.get_data() == n + 1)] = value_dict[n]
        else:
            out_mat = self.parcelList

        n = nb.Nifti1Image(out_mat, self.nbiso.get_affine(), header=self.isoHeader)

        nb.save(n, outname + '.nii')

    def apply_threshold(self, threshold_type=None, value=0., use_absolute=False):
        """
        Threshold the adjacency matrix to determine which nodes are linked by edges.

        threshold_type : The type of threshold applied. Four options are available:
            "edgePC" -> retain the most connected edges as a percentage of the total possible number of edges.
                         "value" must be between 0 and 100
            "totalEdges" -> retain the most strongly connected edges
            "tVal" -> retain edges with a weight greater or equal than value
            None -> all possible edges are created
        value : Value according to threshold_type
        use_absolute : Thresholding by absolute value. For example, if this is set to False, a \
            weight of 1 is stronger than -1. If this is set to True, these values are equally \
            strong. This affects thresholding with "edgePC", "totalEdges" and "tVal". In case of \
            "tVal", it will threshold for weights >= abs(tVal) and <= -abs(tVal)
        """

        # Controlling input
        if threshold_type not in ["edgePC", "totalEdges", "tVal", None]:
            raise TypeError("Not a valid threshold_type for apply_threshold()")
        if threshold_type == "edgePC" and (value < 0 or value > 100):
            raise TypeError("Invalid value for edgePC in apply_threshold()")

        # Creating the array with weights and edges
        upper_values = np.triu_indices(np.shape(self.adjMat)[0], k=1)
        weights = []

        # Creating weights in which each element is (node1, node2, weight)
        for x in np.nditer(upper_values):
            if not np.isnan(self.adjMat[x[0]][x[1]]):
                weights.append((int(x[0]),
                                int(x[1]),
                                self.adjMat[x[0]][x[1]]))

        # If directed, also add the lower down part of the adjacency matrix
        if self.directed:
            below_values = np.tril_indices(np.shape(self.adjMat)[0], k=-1)
            for x in np.nditer(below_values):
                if not np.isnan(self.adjMat[x[0]][x[1]]):
                    weights.append((int(x[0]),
                                    int(x[1]),
                                    self.adjMat[x[0]][x[1]]))

        # Filtering weights when with threshold_type of "edgePC" or "totalEdges"
        if threshold_type in ["edgePC", "totalEdges"]:
            # Sorting ascending
            if use_absolute:
                weights.sort(key=lambda x: abs(x[2]))
            else:
                weights.sort(key=lambda x: x[2])

            # Getting the number of edges to include
            if threshold_type == 'edgePC':
                edgenum = int((value / 100.) * len(weights))
            else:  # totalEdges
                edgenum = int(value)

            # Removing weak edges
            if edgenum <= 0:
                weights = []
            elif edgenum > len(weights):
                pass  # include all weights
            else:
                weights = weights[-edgenum:]

        # remove previous edges
        self.G.remove_edges_from(list(self.G.edges()))

        # Adding the edges
        for e in weights:
            if threshold_type == 'tVal' and use_absolute:
                if e[2] >= abs(value) or e[2] <= -abs(value):
                    self.G.add_edge(e[0], e[1], weight=e[2])
            elif threshold_type == 'tVal' and not use_absolute:
                if e[2] >= value:
                    self.G.add_edge(e[0], e[1], weight=e[2])
            else:  # None, edgePC, totalEdges
                self.G.add_edge(e[0], e[1], weight=e[2])

        # Apply existing properties
        if self.update_properties_after_threshold:
            self._add_properties(self.nodeProperties)
            self._add_properties(self.edgeProperties)

    def reconstruct_adj_mat(self):
        """
        It redefines the adjacency matrix from the edges' weights of G
        It assumes that size of adjMat is maintained
        """
        self.adjMat[:] = np.nan

        for e in self.G.edges():
            self.update_adj_mat(e)

    def update_adj_mat(self, edge):
        """
        It updates the adjacency matrix by bringing the weight of an edge in G \
        to the adjacency matrix

        edge : The edge in G to bring to adjMat
        """

        try:
            w = self.G.edges[edge[0], edge[1]][ct.WEIGHT]
            self.adjMat[edge[0], edge[1]] = w

            if not self.directed:
                self.adjMat[edge[1], edge[0]] = w
        except KeyError as error:
            import sys
            _, _, tb = sys.exc_info()
            raise KeyError(error, "Edge does not exist in G or doesn't have WEIGHT property").with_traceback(tb)
        except IndexError:
            import sys
            _, _, tb = sys.exc_info()
            raise IndexError("adjMat too small to have such an edge").with_traceback(tb).with_traceback(tb)

    def local_thresholding(self, threshold_type=None, value=0.):
        """
        Threshold the adjacency matrix by building from the minimum spanning tree (MST) and adding successive N-nearest
        neighbour degree graphs.
        Thus, if you want to have a local thresholding of N edges when the MST has more than N edges, thresholding will
        retain the MST
        It only works for undirected graphs

        threshold_type : The type of threshold applied. Three options are available:
            "edgePC" -> retain the most connected edges as a percentage of the total possible number of edges. "value"
                        must be between 0 and 100
            "totalEdges" -> retain the most strongly connected edges
            None -> retain the minimum spanning tree
        value : Value according to threshold_type
        """

        # Controlling input
        if threshold_type not in ["edgePC", "totalEdges", None]:
            raise TypeError("Not a valid threshold_type for local_thresholding()")
        if self.directed:
            raise TypeError("local_thresholding() not available for directed graphs")
        if threshold_type == "edgePC" and (value < 0 or value > 100):
            raise TypeError("Invalid value for edgePC for local_thresholding()")

        # Putting all the edges in the G object for local thresholding
        self.apply_threshold()

        if not nx.is_connected(self.G):
            raise TypeError("Adjacency Matrix is not connected. Impossible to execute local_thresholding()")

        # create minimum spanning tree
        t = nx.minimum_spanning_tree(self.G)

        if not threshold_type:
            self.G = t
            return  # Nothing else to do, just return
        elif threshold_type == 'edgePC':
            # find threshold as a percentage of total possible edges
            upper_values = np.triu_indices(np.shape(self.adjMat)[0], k=1)
            # only flatten the upper right part of the matrix
            weights = np.array(self.adjMat[upper_values])

            # remove NaNs
            weights = weights[~np.isnan(weights)]
            # calculate percentage
            edgenum = int(value / 100. * len(weights))
        else:  # 'totalEdges' option
            edgenum = value

        len_edges = len(t.edges())
        if len_edges > edgenum:
            print("Warning: The minimum spanning tree already has: " + str(len_edges) + " edges, select more edges.",
                  "Local Threshold will be applied by just retaining the Minimum Spanning Tree")
            self.local_thresholding()
            return

        k = 1  # number of degrees for NNG
        while len_edges < edgenum:
            # create nearest neighbour graph
            nng = self._nng(k)

            # remove edges from the NNG that exist already in the new graph/MST
            nng.remove_edges_from(t.edges())

            # Ending condition. No more edges to add so break the cycle
            if len(nng.edges()) == 0:
                break

            # add weights to NNG
            for e in nng.edges():
                nng.edges[e[0], e[1]][ct.WEIGHT] = self.adjMat[e[0], e[1]]

            # get a list of edges from the NNG in order of weight
            edge_list = sorted(nng.edges(data=True), key=lambda t: t[2][ct.WEIGHT], reverse=True)

            # add edges to graph in order of connectivity strength
            for edge in edge_list:
                t.add_edges_from([edge])
                len_edges = len(t.edges())
                if len_edges >= edgenum:
                    break

            k += 1

        self.G = t

        # Apply existing properties
        if self.update_properties_after_threshold:
            self._add_properties(self.nodeProperties)
            self._add_properties(self.edgeProperties)

    def binarise(self):
        """
        Removes weighting from edges by assigning a weight of 1 to the existing edges
        """
        for edge in self.G.edges(data=True):
            edge[2][ct.WEIGHT] = 1

    def make_edges_absolute(self):
        """
        Makes all the edges in self.G absolute
        """
        for edge in self.G.edges(data=True):
            edge[2][ct.WEIGHT] = abs(edge[2][ct.WEIGHT])

    def remove_unconnected_nodes(self):
        """
        Removes nodes with no connections
        """
        node_list = [v for v in self.G.nodes() if self.G.degree(v) == 0]
        self.G.remove_nodes_from(node_list)

    def _nng(self, k):
        """ Private method to help local thresholding by creating a k-nearest neighbour graph"""
        g = nx.Graph()
        nodes = list(range(len(self.adjMat[0])))

        g.add_nodes_from(nodes)

        for i in nodes:
            l = np.ma.masked_array(self.adjMat[i, :], mask=np.isnan(self.adjMat[i]))
            l.mask[i] = True

            for j in range(k):
                node = np.argmax(l)

                if not np.isnan(self.adjMat[i, node]):
                    g.add_edge(i, node)

                l.mask[node] = True

        return g

    def find_spatially_nearest(self, node_list, contra=False, midline=44.5, connected=True, threshold=None):
        # find the spatially closest node as no topologically close nodes exist
        if isinstance(node_list, list):
            duff_node = random.choice(node_list)
        else:
            duff_node = node_list

        nodes = [v for v in self.G.nodes() if v != duff_node]
        nodes = [v for v in nodes if v not in node_list]

        shortestnode = (None, None)
        tlist = []  # list of closest nodes according to threshold
        if not threshold:
            tval = 0.

        # get the contralaterally closest node if desired
        pos = [v for v in self.G.nodes[duff_node][ct.XYZ]]
        if contra:
            if pos[0] < midline:
                pos[0] = midline + (midline - pos[0])
            else:
                pos[0] = midline + (pos[0] - midline)
        pos = tuple(pos)

        for node in nodes:
            try:
                distance = np.linalg.norm(np.array(pos - np.array(self.G.nodes[node]['xyz'])))
                if distance < tval:
                    tlist.append(node)
            except:
                print("Finding the spatially nearest node requires x,y,z values")

            if shortestnode[0]:
                if distance < shortestnode[1]:
                    if connected:
                        if self.G.degree(node) > 0:
                            shortestnode = (node, distance)
                    else:
                        shortestnode = (node, distance)

            else:
                if connected:
                    if self.G.degree(node) > 0:
                        shortestnode = (node, distance)
                else:
                    shortestnode = (node, distance)

        if threshold:
            return tlist
        else:
            return shortestnode[0]

    def find_linked_nodes(self):
        """
        It gives to each node a list containing the linked nodes.
        If Graph is undirected, and there is an edge (1,2), node 1 is linked \
         to node 2, and vice-versa. If Graph is directed, thus node 1 is linked \
         to node 2, but not the other way around
        This property can be accessed through ct.LINKED_NODES
        Be sure to call this method again if you threshold your brainObj again
        """

        # Resetting all nodes from some past information (few edges might not \
        #  be able to reset this field in all nodes)
        for n in self.G.nodes(data=True):
            n[1][ct.LINKED_NODES] = []

        for l in self.G.edges():
            # add to list of connecting nodes for each participating node
            self.G.nodes[l[0]][ct.LINKED_NODES].append(l[1])

            if not self.directed:
                self.G.nodes[l[1]][ct.LINKED_NODES].append(l[0])

    def weight_to_distance(self):
        """
        It inverts all the edges' weights so they become equivalent to a distance measure. \
        With a weight, the higher the value the stronger the connection. With a distance, \
        the higher the value the "weaker" the connection.
        In this case there is no measurement unit for the distance, as it is just \
        a conversion from the weights.
        The distances can be accessed in each node's property with ct.DISTANCE
        """
        edge_list = [v[2][ct.WEIGHT] for v in self.G.edges(data=True)]

        # get the maximum edge value, plus a small correction to keep the values above zero
        # the correction is the inverse of the number of nodes - designed to keep
        # calculations of efficiency sensible
        emax = np.max(edge_list) + 1 / float(self.G.number_of_nodes())

        for edge in self.G.edges():
            self.G.edges[edge[0], edge[1]][ct.DISTANCE] = emax - self.G.edges[edge[0], edge[1]][
                ct.WEIGHT]  # convert weights to a positive distance

    def copy_hemisphere(self, hsphere="R", midline=44.5):
        """
        This copies all the nodes and attributes from one hemisphere to the other, deleting any pre-existing
        data on the contralateral side. Particularly useful when you only have data from a single
        hemisphere available.
        """
        if hsphere == "L":
            for node in self.G.nodes():
                if self.G.nodes[node][ct.XYZ][0] < midline:
                    self.G.add_node(str(node) + "R")
                    self.G.nodes[str(node) + "R"] = {v: w for v, w in self.G.nodes[node].items()}
                    pos = self.G.nodes[node][ct.XYZ]
                    pos = (midline + (midline - pos[0]), pos[1], pos[2])
                    self.G.nodes[str(node) + "R"][ct.XYZ] = pos
                else:
                    self.G.remove_node(node)

        elif hsphere == "R":
            for node in self.G.nodes():
                if self.G.nodes[node][ct.XYZ][0] > midline:
                    self.G.add_node(str(node) + "L")
                    self.G.nodes[str(node) + "L"] = {v: w for v, w in self.G.nodes[node].items()}
                    pos = self.G.nodes[node][ct.XYZ]
                    self.G.nodes[str(node) + "L"][ct.XYZ] = (midline - (pos[0] - midline), pos[1], pos[2])
                else:
                    self.G.remove_node(node)
