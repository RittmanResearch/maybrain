import unittest

from maybrain import brain as mbt
from maybrain import constants as ct
from maybrain import utils
import networkx as nx
import numpy as np


class TestBrainObj(unittest.TestCase):
    """
    Test Brain class from maybrain
    """

    def setUp(self):
        self.a = mbt.Brain()
        self.SMALL_FILE = "test/data/3d_grid_adj.txt"
        self.SMALL_NEG_FILE = "test/data/3d_grid_adj_neg.txt"
        self.MODIF_FILE = "test/data/3d_grid_adj2.txt"
        self.COORD_FILE = "test/data/3d_grid_coords.txt"
        self.PROPS_FILE = "test/data/3d_grid_properties.txt"

    def test_constructor(self):
        self.assertIsInstance(self.a, mbt.Brain)
        self.assertIsInstance(self.a.G, nx.Graph)
        self.assertNotIsInstance(self.a.G, nx.DiGraph)
        self.assertFalse(self.a.directed)
        self.assertEqual(self.a.G.number_of_nodes(), 0)
        self.assertEqual(self.a.G.number_of_edges(), 0)
        self.assertEqual(self.a.adjMat, None)

    def test_importAdjFile(self):
        self.assertRaises(IOError, self.a.import_adj_file, "sdfasdf")
        self.a.import_adj_file(self.SMALL_FILE)
        self.assertEqual(self.a.adjMat.shape, (4, 4))
        self.assertEqual(self.a.adjMat[0][0], 0.802077230054)
        self.assertEqual(self.a.G.number_of_nodes(), 4)
        self.assertEqual(len(self.a.G.edges()), 0)

        b = mbt.Brain()
        b.import_adj_file(self.MODIF_FILE, delimiter=",", nodes_to_exclude=[2, 4])
        # Confirm general info
        self.assertEqual(b.adjMat.shape, (15, 15))
        self.assertEqual(b.G.number_of_nodes(), 13)
        self.assertEqual(b.G.number_of_edges(), 0)
        self.assertEqual(b.adjMat[0][0], 0)
        self.assertEqual(b.adjMat[1][0], 1.541495524150943430e-01)
        self.assertTrue(np.isnan(b.adjMat[6][0]))
        # Confirm excluded nodes
        self.assertTrue(all(np.isnan(x) for x in b.adjMat[:, 2]))
        self.assertTrue(all(np.isnan(x) for x in b.adjMat[:, 4]))
        self.assertTrue(all(np.isnan(x) for x in b.adjMat[2, :]))
        self.assertTrue(all(np.isnan(x) for x in b.adjMat[4, :]))

    def test_importSpatialInfo(self):
        self.assertRaises(FileNotFoundError, self.a.import_spatial_info, "sdfasdf")
        self.a.import_adj_file(self.SMALL_FILE)
        self.a.import_spatial_info(self.COORD_FILE)

        attrs = mbt.nx.get_node_attributes(self.a.G, ct.XYZ)
        self.assertEqual(self.a.G.number_of_nodes(), 4)
        self.assertEqual(self.a.G.number_of_edges(), 0)
        self.assertEqual(attrs[0][0], 0)
        self.assertEqual(attrs[0][1], 0)
        self.assertEqual(attrs[0][2], 0)
        self.assertEqual(attrs[3][0], 2.)
        self.assertEqual(attrs[3][1], 2.)
        self.assertEqual(attrs[3][2], 0)

        attrs2 = mbt.nx.get_node_attributes(self.a.G, ct.ANAT_LABEL)
        self.assertEqual(attrs2[0], '0')
        self.assertEqual(attrs2[3], '3')

    def test_applyThreshold(self):
        self.a.import_adj_file(self.MODIF_FILE, delimiter=",", nodes_to_exclude=[2, 4])
        self.a.import_spatial_info(self.COORD_FILE)  # making sure this doesn't influence the rest
        self.a.apply_threshold()
        all_edges = self.a.G.number_of_edges()
        # Although there are 3 NAs in the file, just the upper half of the matrix is considered
        degrees = self.a.G.degree()
        self.assertEqual(degrees[0], 11)
        self.assertRaises(KeyError, lambda: degrees[2])
        self.assertRaises(KeyError, lambda: degrees[4])
        self.assertEqual(degrees[5], 12)

        accepted_edges = ((0, 1), (0, 5), (1, 5), (3, 5))
        not_accepted = ((0, 0), (1, 1), (2, 1), (0, 2), (0, 3), (0, 4), (1, 2), (1, 3), (1, 4))
        self.assertTrue(all(e in self.a.G.edges() for e in accepted_edges))
        self.assertTrue(all(e not in self.a.G.edges() for e in not_accepted))
        self.assertEqual(self.a.G.number_of_edges(), 76)  # (15*15 - 15) /2 - 2 - 27

        # edgePC
        b = mbt.Brain()
        b.import_adj_file(self.MODIF_FILE, delimiter=",", nodes_to_exclude=[2, 4])
        for val in [True, False]:
            b.apply_threshold(threshold_type="edgePC", value=10.5, use_absolute=val)
            self.assertEqual(b.G.number_of_edges(), int(0.105 * all_edges))
            self.assertTrue((1, 12) in b.G.edges())
            b.apply_threshold(threshold_type="edgePC", value=0, use_absolute=val)
            self.assertEqual(b.G.number_of_edges(), 0)

        # totalEdges
        for val in [True, False]:
            b.apply_threshold(threshold_type="totalEdges", value=3, use_absolute=val)
            self.assertEqual(b.G.number_of_edges(), 3)
            self.assertTrue((1, 12) in b.G.edges())
            b.apply_threshold(threshold_type="totalEdges", value=1000, use_absolute=val)
            self.assertEqual(b.G.number_of_edges(), 76)
            b.apply_threshold(threshold_type="totalEdges", value=0, use_absolute=val)
            self.assertEqual(b.G.number_of_edges(), 0)

        # tVal
        b.apply_threshold(threshold_type="tVal", value=3)
        self.assertTrue(all(e[2][ct.WEIGHT] >= 3 for e in b.G.edges(data=True)))
        self.assertEqual(b.G.number_of_edges(), 0)
        b.apply_threshold(threshold_type="tVal", value=6.955292039622642530e-01)
        self.assertTrue(all(e[2][ct.WEIGHT] >= 6.955292039622642530e-01 for e in b.G.edges(data=True)))
        self.assertEqual(b.G.number_of_edges(), 1)
        b.apply_threshold(threshold_type="tVal", value=0.5)
        self.assertTrue(all(e[2][ct.WEIGHT] >= 0.5 for e in b.G.edges(data=True)))

        # directed
        c = mbt.Brain(directed=True)
        c.import_adj_file(self.MODIF_FILE, delimiter=",")
        for val in [True, False]:
            c.apply_threshold(use_absolute=val)
            all_edges = c.G.number_of_edges()
            self.assertEqual(c.G.number_of_edges(), 207)  # 15*15 - 15 -3NAs
            c.apply_threshold(threshold_type="edgePC", value=10.5, use_absolute=val)
            self.assertEqual(c.G.number_of_edges(), int(0.105 * all_edges))
            c.apply_threshold(threshold_type="totalEdges", value=76, use_absolute=val)
            self.assertEqual(c.G.number_of_edges(), 76)
            c.apply_threshold(threshold_type="totalEdges", value=10000, use_absolute=val)
            self.assertEqual(c.G.number_of_edges(), 207)
        c.apply_threshold(threshold_type="tVal", value=0.5)
        self.assertTrue(all(e[2][ct.WEIGHT] >= 0.5 for e in c.G.edges(data=True)))

        # Absolute thresholding
        c.import_adj_file(self.SMALL_NEG_FILE)
        c.apply_threshold(threshold_type="totalEdges", value=2, use_absolute=True)
        self.assertTrue(c.G.edges[1, 2][ct.WEIGHT] == -0.843798947781)
        self.assertTrue(c.G.edges[2, 0][ct.WEIGHT] == 0.858463902674)
        c.apply_threshold(threshold_type="tVal", value=0.5, use_absolute=True)
        self.assertTrue(all(e[2][ct.WEIGHT] >= 0.5 or e[2][ct.WEIGHT] <= -0.5
                            for e in c.G.edges(data=True)))
        b.import_adj_file(self.SMALL_NEG_FILE)
        b.apply_threshold(threshold_type="totalEdges", value=1, use_absolute=True)
        self.assertTrue(b.G.edges[1, 2][ct.WEIGHT] == -0.843798947781)
        b.apply_threshold(threshold_type="edgePC", value=20, use_absolute=True)
        self.assertTrue(b.G.edges[1, 2][ct.WEIGHT] == -0.843798947781)
        self.assertEqual(b.G.number_of_edges(), 1)

    def test_binarise(self):
        self.a.import_adj_file(self.MODIF_FILE, delimiter=",")
        self.a.apply_threshold()
        self.a.binarise()
        self.assertTrue(all(e[2][ct.WEIGHT] == 1 for e in self.a.G.edges(data=True)))

    def test_make_absolute(self):
        c = mbt.Brain(directed=True)
        c.import_adj_file(self.SMALL_NEG_FILE)
        c.apply_threshold()
        c.make_edges_absolute()
        self.assertTrue(all(e[2][ct.WEIGHT] >= 0 for e in c.G.edges(data=True)))

    def test_localThreshold(self):
        self.a.import_adj_file(self.MODIF_FILE, delimiter=",")
        self.a.apply_threshold()
        temp = nx.minimum_spanning_tree(self.a.G)

        # Normal
        self.a.local_thresholding()
        self.assertEqual(temp.number_of_edges(), self.a.G.number_of_edges())
        self.assertTrue(nx.is_connected(self.a.G))

        # totalEdges
        # normal totalEdges
        self.a.local_thresholding(threshold_type="totalEdges", value=20)
        self.assertEqual(self.a.G.number_of_edges(), 20)
        self.assertTrue(nx.is_connected(self.a.G))
        # short totalEdges
        self.a.local_thresholding(threshold_type="totalEdges", value=1)
        self.assertEqual(self.a.G.number_of_edges(), temp.number_of_edges())
        self.assertTrue(nx.is_connected(self.a.G))
        # bigger totalEdges
        self.a.local_thresholding(threshold_type="totalEdges", value=500000)
        self.assertTrue(nx.is_connected(self.a.G))

        # edgePC
        self.a.apply_threshold()
        all_edges = self.a.G.number_of_edges()

        self.a.local_thresholding(threshold_type="edgePC", value=100)
        self.assertEqual(self.a.G.number_of_edges(), all_edges)
        self.assertTrue(nx.is_connected(self.a.G))

        self.a.local_thresholding(threshold_type="edgePC", value=20)
        self.assertEqual(self.a.G.number_of_edges(), int(0.2 * all_edges))
        self.assertTrue(nx.is_connected(self.a.G))

    def test_AdjMat(self):
        self.a.import_adj_file(self.MODIF_FILE, delimiter=",")
        self.a.apply_threshold(threshold_type="totalEdges", value=1)
        self.a.reconstruct_adj_mat()
        # Size must be 2 because adjMat is undirected
        self.assertEqual(len([e for e in self.a.adjMat.flatten() if not np.isnan(e)]), 2)

        self.a.G.add_edge(2, 2, weight=23)
        self.a.G.add_edge(3, 3, weight=12)

        self.a.update_adj_mat((2, 2))
        self.assertEqual(self.a.adjMat[2][2], 23)
        self.assertTrue(np.isnan(self.a.adjMat[3][3]))

        self.a.reconstruct_adj_mat()
        self.assertEqual(self.a.adjMat[2][2], 23)
        self.assertEqual(self.a.adjMat[3][3], 12)

        # Error handling
        self.assertRaises(KeyError, self.a.update_adj_mat, (1, 1))
        self.a.G.add_edge(50, 50)  # dummy
        self.assertRaises(KeyError, self.a.update_adj_mat, (50, 50))
        self.a.G.edges[50, 50][ct.WEIGHT] = 12
        self.assertRaises(IndexError, self.a.update_adj_mat, (50, 50))

    def test_removeUnconnectedNodes(self):
        self.a.import_adj_file(self.MODIF_FILE, delimiter=",")
        self.a.apply_threshold(threshold_type="totalEdges", value=1)
        self.assertEqual(self.a.G.number_of_nodes(), 15)  # info from the file
        self.a.remove_unconnected_nodes()
        # There is only 2 nodes left connecting the only existing edge
        self.assertEqual(self.a.G.number_of_nodes(), 2)

    def test_percentages(self):
        self.a.import_adj_file(self.SMALL_FILE)
        self.a.apply_threshold(threshold_type="totalEdges", value=1)
        self.assertEqual(utils.threshold_to_percentage(self.a, 0.6), 3 / 6)  # from the file

        self.assertEqual(utils.percent_connected(self.a), 1 / 6)
        self.a.apply_threshold()
        self.assertEqual(utils.percent_connected(self.a), 1)

        # directed
        c = mbt.Brain(directed=True)
        c.import_adj_file(self.SMALL_FILE)
        self.assertEqual(utils.threshold_to_percentage(c, 0.6), 5 / 12)

        c.apply_threshold(threshold_type="totalEdges", value=1)
        self.assertEqual(utils.percent_connected(c), 1 / 12)
        c.apply_threshold()
        self.assertEqual(utils.percent_connected(c), 1)
        c.apply_threshold(threshold_type="totalEdges", value=0)
        self.assertEqual(utils.percent_connected(c), 0)

    def test_linkedNodes(self):
        self.a.import_adj_file(self.SMALL_FILE)
        self.a.apply_threshold()
        for n in self.a.G.nodes(data=True):
            self.assertRaises(KeyError, lambda: n[1][ct.LINKED_NODES])

        self.a.find_linked_nodes()
        for n in self.a.G.nodes(data=True):
            n[1][ct.LINKED_NODES]  # it should not raise any exception

        self.assertTrue(sorted(self.a.G.node[0][ct.LINKED_NODES]) == [1, 2, 3])
        self.assertTrue(sorted(self.a.G.node[1][ct.LINKED_NODES]) == [0, 2, 3])

        # directed
        c = mbt.Brain(directed=True)
        c.import_adj_file(self.SMALL_FILE)
        c.apply_threshold()
        c.find_linked_nodes()
        for n in c.G.nodes(data=True):
            n[1][ct.LINKED_NODES]  # it should not raise any exception

        self.assertTrue(sorted(c.G.node[0][ct.LINKED_NODES]) == [1, 2, 3])
        self.assertTrue(sorted(c.G.node[1][ct.LINKED_NODES]) == [0, 2, 3])

    def test_weightToDistance(self):
        self.a.import_adj_file(self.MODIF_FILE, delimiter=",")
        self.a.apply_threshold()
        self.a.weight_to_distance()

        edge_list = [v[2][ct.WEIGHT] for v in self.a.G.edges(data=True)]
        emax = np.max(edge_list) + 1 / float(self.a.G.number_of_nodes())

        self.assertTrue(all(emax == e[2][ct.DISTANCE] + e[2][ct.WEIGHT] for e in self.a.G.edges(data=True)))
        self.assertTrue(all(e[2][ct.DISTANCE] > 0 for e in self.a.G.edges(data=True)))

    def test_properties(self):
        self.a.import_adj_file(self.SMALL_FILE)
        self.a.apply_threshold()
        self.a.import_properties(self.PROPS_FILE)

        for e in range(2):
            # 2nd iteration
            if e == 1:
                self.a.apply_threshold(threshold_type="totalEdges", value=0)
                self.a.update_properties_after_threshold = True
                self.a.apply_threshold()

            self.assertRaises(KeyError, lambda: self.a.G.node[2]['colour'])
            self.assertRaises(KeyError, lambda: self.a.G.node[6])  # testing node 6 is not created
            self.assertTrue(self.a.G.node[1]['colour'], 'red')
            self.assertTrue(self.a.G.node[0]['colour'], 'blue')
            self.assertTrue(self.a.G.node[3]['colour'], 'red')
            # edges
            self.assertTrue(self.a.G.edges[0, 2]['colour'], 'red')
            self.assertTrue(self.a.G.edges[2, 0]['colour'], 'red')
            self.assertTrue(self.a.G.edges[1, 3]['colour'], 'green')

    def test_highlights(self):
        self.a.import_adj_file(self.SMALL_FILE)
        self.a.apply_threshold()
        self.a.import_properties(self.PROPS_FILE)

        utils.highlight_from_conds(self.a, 'colour', 'eq', 'blue', mode='node')
        self.assertEqual(utils.highlights[1].nodes, [0])
        self.assertEqual(utils.highlights[1].edges, [])

        utils.highlight_from_conds(self.a, ct.WEIGHT, 'geq', 0.84379894778099995, mode='edge', label=2)
        self.assertEqual(utils.highlights[2].nodes, [])
        self.assertEqual(utils.highlights[2].edges, [(1, 2)])

        utils.highlight_from_conds(self.a, ct.WEIGHT, 'gt', 0.84379894778099995, mode='edge')
        self.assertEqual(utils.highlights[3].nodes, [])
        self.assertEqual(utils.highlights[3].edges, [])

        utils.highlight_from_conds(self.a, ct.WEIGHT, 'leq', 0.203602458588, mode='edge', label=4)
        self.assertEqual(utils.highlights[4].nodes, [])
        self.assertEqual(sorted(utils.highlights[4].edges), [(0, 2), (0, 3)])

        utils.highlight_from_conds(self.a, ct.WEIGHT, 'lt', 0.203602458588, mode='edge', label=5)
        self.assertEqual(utils.highlights[5].nodes, [])
        self.assertEqual(utils.highlights[5].edges, [(0, 3)])

        utils.highlight_from_conds(self.a, ct.WEIGHT, 'in[]', [0.16390494700200001, 0.60080034391699999], mode='node|edge', label=5)
        self.assertEqual(utils.highlights[5].nodes, [])
        self.assertEqual(sorted(utils.highlights[5].edges), [(0, 1), (0, 2), (0, 3), (1, 3)])
        utils.highlight_from_conds(self.a, ct.WEIGHT, 'in(]', [0.16390494700200001, 0.60080034391699999], mode='edge|node', label=5)
        self.assertEqual(utils.highlights[5].nodes, [])
        self.assertEqual(sorted(utils.highlights[5].edges), [(0, 1), (0, 2), (1, 3)])
        utils.highlight_from_conds(self.a, ct.WEIGHT, 'in[)', [0.16390494700200001, 0.60080034391699999], mode='edge', label=5)
        self.assertEqual(utils.highlights[5].nodes, [])
        self.assertEqual(sorted(utils.highlights[5].edges), [(0, 2), (0, 3), (1, 3)])
        utils.highlight_from_conds(self.a, ct.WEIGHT, 'in()', [0.16390494700200001, 0.60080034391699999], mode='edge', label=5)
        self.assertEqual(utils.highlights[5].nodes, [])
        self.assertEqual(sorted(utils.highlights[5].edges), [(0, 2), (1, 3)])

        utils.highlight_from_conds(self.a, 'colour', 'in', ['red', 'blue'], mode='node|edge')
        self.assertEqual(sorted(utils.highlights[6].nodes), [0, 1, 3])
        self.assertEqual(utils.highlights[6].edges, [(0, 2)])

        utils.make_highlight(edge_inds=[(1, 2), (4, 5)], nodes_inds=[2, 4], label='custom')
        self.assertEqual(sorted(utils.highlights['custom'].nodes), [2,4])
        self.assertEqual(sorted(utils.highlights['custom'].edges), [(1, 2), (4, 5)])

        self.a.import_spatial_info(self.COORD_FILE)
        utils.highlight_from_conds(self.a, 'y', 'eq', 2, mode='node', label='custom')
        self.assertEqual(sorted(utils.highlights['custom'].nodes), [2, 3])
        self.assertEqual(utils.highlights['custom'].edges, [])

if __name__ == '__main__':
    unittest.main()
