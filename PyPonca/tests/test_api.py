import unittest
import numpy as np

import pyponca
from pyponca import _pyponca

"""
    Test the API, not that its results are correct !
"""
class TestAPI(unittest.TestCase):
    def setUp(self):
        self.randomtensor = lambda shape, dtype: np.random.uniform(0, 1, shape).astype(dtype)

        self.N = 16
        self.pc2d = self.randomtensor((self.N, 2), np.float64)
        self.pc3d = self.randomtensor((self.N, 3), np.float64)
        self.pc2f = self.randomtensor((self.N, 2), np.float32)
        self.pc3f = self.randomtensor((self.N, 3), np.float32)

        self.pcs = [self.pc2d, self.pc3d, self.pc2f, self.pc3f]

    def test_pointcloud(self):
        """
            Test that we can build a PointCloud and that mangling works properly
        """
        self.assertIsInstance(pyponca.PointCloud(self.pc2d).object, _pyponca.PointCloud2dPN)
        self.assertIsInstance(pyponca.PointCloud(self.pc3d).object, _pyponca.PointCloud3dPN)
        self.assertIsInstance(pyponca.PointCloud(self.pc2f).object, _pyponca.PointCloud2fPN)
        self.assertIsInstance(pyponca.PointCloud(self.pc3f).object, _pyponca.PointCloud3fPN)

        # Test with normal
        self.assertIsInstance(pyponca.PointCloud(self.pc2d, self.pc2d).object, _pyponca.PointCloud2dPN)
        self.assertIsInstance(pyponca.PointCloud(self.pc3d, self.pc3d).object, _pyponca.PointCloud3dPN)
        self.assertIsInstance(pyponca.PointCloud(self.pc2f, self.pc2f).object, _pyponca.PointCloud2fPN)
        self.assertIsInstance(pyponca.PointCloud(self.pc3f, self.pc3f).object, _pyponca.PointCloud3fPN)

        # Test that we can't construct mismatched objects
        for i in range(len(self.pcs)):
            for j in range(len(self.pcs)): 
                if i != j:
                    with self.assertRaises((RuntimeError, TypeError)):
                        cls = pyponca.PointCloud(self.pcs[i], self.pcs[j])
                        print(cls._cls)

    def test_kdtree(self):
        """
            Test that we can build a KDtree and that mangling works properly
        """
        # Test we can build a kdtree
        self.assertIsInstance(pyponca.KDTree(self.pc2d).object, _pyponca.KdTree2dPN)
        self.assertIsInstance(pyponca.KDTree(self.pc3d).object, _pyponca.KdTree3dPN)
        self.assertIsInstance(pyponca.KDTree(self.pc2f).object, _pyponca.KdTree2fPN)
        self.assertIsInstance(pyponca.KDTree(self.pc3f).object, _pyponca.KdTree3fPN)

    def test_kdtree_query(self):
        """
            Test queries can be performed on kdtrees
        """
        for pc in self.pcs:
            kdtree = pyponca.KDTree(pc)
            kdtree.rangeNeighbors(pc[0], 1)
            kdtree.rangeNeighbors(pc, np.ones((pc.shape[0])))
    
    def test_compute(self):
        """
            Test that we can attach pointclouds and filter to a compute object
        """
        for co in _pyponca.ComputeObjectList:
            for pc in self.pcs:
                for flt in pyponca.Filters:
                    try:
                        ts = self.randomtensor((self.N), pc.dtype)

                        object = pyponca.__dict__[co]()
                        object.setNeighborFilter(pc, ts, flt)
                        object.addComputation(_pyponca.Computation.DRY)
                        result = object.compute(pc)
                    except NotImplementedError:
                        # We let notimplementederror that may come from classes that only support 3d
                        pass
    
    def test_compute_func(self):
        """
            Test for potential memory errors and settings in a classical loop
        """
        for co in _pyponca.ComputeObjectList:
            for pc in [self.pc3d, self.pc3f]:
                ts = self.randomtensor((self.N), pc.dtype)

                cls = pyponca.__dict__[co]

                # Helper function to try with different data
                def test(data):
                    obj = cls()
                    obj.setNeighborFilter(pc, ts)
                    obj.addComputation(_pyponca.Computation.DRY)
                    obj.compute(data)

                test(pc)
                test(pyponca.PointCloud(pc))
                test(pyponca.KDTree(pc))

if __name__ == "__main__":
    unittest.main()
