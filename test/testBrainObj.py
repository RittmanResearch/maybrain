import unittest

from maybrain import brainObj as mbt
import networkx as nx

class TestBrainObj(unittest.TestCase):
    """
    Test brainObj class from maybrain
    """
    
    def setUp(self):
        self.a = mbt.brainObj() 

    def test_constructor(self):
        self.assertIsInstance(self.a, mbt.brainObj)
        self.assertIsInstance(self.a.G, nx.Graph)
        self.assertNotIsInstance(self.a.G, nx.DiGraph)

    def test_importAdjFile(self):
        self.assertEqual(self.a.importAdjFile("sdfasdf"), -1)
        self.a.importAdjFile("test/data/3d_grid_adj.txt")
        self.assertEqual(self.a.adjMat.shape, (4,4))
        self.assertEqual(self.a.adjMat[0][0], 0.802077230054)
        b = mbt.brainObj()
        b.importAdjFile("test/data/3d_grid_adj2.txt", delimiter=",")
        self.assertEqual(b.adjMat.shape, (15,15))
        self.assertEqual(b.adjMat[0][0], 0)

if __name__ == '__main__':
    unittest.main()
