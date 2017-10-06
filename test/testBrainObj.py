import unittest

from maybrain import brainObj as mbt
import networkx as nx
import numpy as np

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
        self.assertFalse(self.a.directed)

    def test_importAdjFile(self):
        self.assertEqual(self.a.importAdjFile("sdfasdf"), -1)
        self.a.importAdjFile("test/data/3d_grid_adj.txt")
        self.assertEqual(self.a.adjMat.shape, (4,4))
        self.assertEqual(self.a.adjMat[0][0], 0.802077230054)
        
        b = mbt.brainObj()
        b.importAdjFile("test/data/3d_grid_adj2.txt", delimiter=",", exclnodes = [2,4])
        # Confirm general info
        self.assertEqual(b.adjMat.shape, (15,15))
        self.assertEqual(b.adjMat[0][0], 0)
        self.assertEqual(b.adjMat[1][0], 1.541495524150943430e-01)
        self.assertTrue(np.isnan(b.adjMat[6][0]))
        #Confirm excluded nodes
        self.assertTrue(all(np.isnan(x) for x in b.adjMat[:,2]))
        self.assertTrue(all(np.isnan(x) for x in b.adjMat[:,4]))
        self.assertTrue(all(np.isnan(x) for x in b.adjMat[2,:]))
        self.assertTrue(all(np.isnan(x) for x in b.adjMat[4,:]))
    
    def test_importSpatialInfo(self):
        self.assertEqual(self.a.importSpatialInfo("sdfasdf"), -1)
        self.a.importAdjFile("test/data/3d_grid_adj.txt")
        self.a.importSpatialInfo("test/data/3d_grid_coords.txt")
       
        attrs = mbt.nx.get_node_attributes(self.a.G, "xyz")
        self.assertEqual(mbt.nx.number_of_nodes(self.a.G), 4)
        self.assertEqual(attrs[0][0], 0)
        self.assertEqual(attrs[0][1], 0)
        self.assertEqual(attrs[0][2], 0)
        self.assertEqual(attrs[3][0], 2.)
        self.assertEqual(attrs[3][1], 2.)
        self.assertEqual(attrs[3][2], 0)
        
        attrs2 = mbt.nx.get_node_attributes(self.a.G, "anatlabel")
        self.assertEqual(attrs2[0], '1')
        self.assertEqual(attrs2[3], '4')
        

if __name__ == '__main__':
    unittest.main()
