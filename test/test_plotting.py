import unittest

from maybrain import brain as mbt
from maybrain import resources as rt
import maybrain.plotting as mpt
import matplotlib.pyplot as plt


class TestPlotting(unittest.TestCase):
    """
    Test plotting features, for now just smoke tests to see that plotting doesn't raise exceptions while executing
    """

    def setUp(self):
        self.a = mbt.Brain()
        self.SMALL_FILE = "test/data/3d_grid_adj.txt"

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


if __name__ == '__main__':
    unittest.main()
