import unittest

import networkx as nx

from maybrain import brain as mbt
from maybrain import resources as rt
import maybrain.plotting as mpt
import maybrain.algorithms as mba
import maybrain.constants as ct
import matplotlib.pyplot as plt


class TestPlottingAndAlgs(unittest.TestCase):
    """
    Test plotting and algorithms features, for now just smoke tests to see that plotting doesn't raise exceptions while executing
    """

    def setUp(self):
        self.a = mbt.Brain()
        self.SMALL_FILE = "test/data/3d_grid_adj.txt"
        self.SMALL_NEG_FILE = "test/data/3d_grid_adj_neg.txt"
        self.COORD_FILE = "test/data/3d_grid_coords.txt"
        self.PROPS_FILE = "test/data/3d_grid_properties_full.txt"

    def test_histograms(self):
        self.a.import_adj_file(self.SMALL_FILE)
        self.a.apply_threshold()
        fig, _ = mpt.plot_weight_distribution(self.a)
        plt.close(fig)

    def test_matrices(self):
        dummy_brain = mbt.Brain()
        dummy_brain.import_adj_file(rt.DUMMY_ADJ_FILE_500)
        dictionary = {0: dummy_brain, 1:dummy_brain}
        fig, _ = mpt.plot_avg_matrix(dictionary)
        plt.close(fig)

        fig, _ = mpt.plot_strength_matrix(dummy_brain)
        plt.close(fig)

    def test_robustness(self):
        self.a.import_adj_file(self.SMALL_FILE)
        self.a.apply_threshold()
        mba.robustness(self.a)

    def test_normalisation(self):
        self.a.import_adj_file(self.SMALL_NEG_FILE)
        self.a.apply_threshold()
        # Getting random graph
        while True:
            try:
                rand = mba.generate_random_graph_from_degree(self.a, throw_exception=True, edge_attrs=[ct.WEIGHT])
                break  # if it reaches here, means randomiser didn't throw any exception, so break while
            except mba.RandomGenerationError as error:
                pass
        orig_deg = nx.degree(self.a.G)
        for n in nx.degree(rand):
            self.assertEqual(n[1], orig_deg[n[0]])

        normalised = mba.normalise_node_wise(self.a,
                                             nx.degree,
                                             init_vals=dict(orig_deg),
                                             n_iter=5)
        self.assertTrue(all(i == 1 for i in normalised.values()))
        self.assertEqual(sum(dict(nx.degree(rand, weight=ct.WEIGHT)).values()),
                         sum(dict(nx.degree(self.a.G, weight=ct.WEIGHT)).values()))

    def test_connectome(self):
        self.a.import_adj_file(self.SMALL_FILE)
        self.a.import_spatial_info(self.COORD_FILE)
        self.a.import_properties(self.PROPS_FILE)
        display = mpt.plot_connectome(self.a)
        display.close()
        display = mpt.plot_connectome(self.a, only_nodes=True)
        display.close()
        display = mpt.plot_connectome(self.a, node_property="colour", node_attributes=["red", "green"])
        display.close()


if __name__ == '__main__':
    unittest.main()
