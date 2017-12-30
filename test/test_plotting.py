import unittest

from maybrain import brain as mbt
import maybrain.plotting as mpt
import matplotlib.pyplot as plt


class TestPlotting(unittest.TestCase):
    """
    Test plotting features, for now just execute and check everything goes ok
    """

    def setUp(self):
        self.a = mbt.Brain()
        self.SMALL_FILE = "data/3d_grid_adj.txt"

    def test_histograms(self):
        self.a.import_adj_file(self.SMALL_FILE)
        self.a.apply_threshold()
        fig, _ = mpt.plot_weight_distribution(self.a)
        plt.close(fig)

if __name__ == '__main__':
    unittest.main()
