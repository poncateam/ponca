import unittest
import contextlib
import numpy as np

import pyponca
from pyponca import _pyponca

try:
    import torch
except ImportError:
    torch = None

"""
    Test the API, not that its results are correct !
"""
class TestAPI(unittest.TestCase):
    # Whether the tests run on the host, where every operation is available
    onHost = True

    def todevice(self, x):
        """
            Sends an array to the device the tests run on
        """
        return x

    def tohost(self, x):
        """
            Converts an array (stored on the device, or returned by pyponca) to a numpy array
        """
        return np.from_dlpack(x)

    def supports(self, co):
        """
            Tells if a compute object can run on the device the tests run on
        """
        return self.onHost or pyponca.__dict__[co].isCudaAvailable

    def assertRaisesUnless(self, supported):
        """
            Expects a RuntimeError when an operation is not supported on the device the tests run on
        """
        return contextlib.nullcontext() if supported else self.assertRaises(RuntimeError)

    def kdtree(self, pc, dense = True):
        """
            Builds a kdtree on the host from pc, then sends it to the device pc is stored on
        """
        return pyponca.KDTree(self.tohost(pc), dense).toDevice(pc)

    def setUp(self):
        self.randomtensor = lambda shape, dtype: self.todevice(np.random.uniform(0, 1, shape).astype(dtype))

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

    def test_pointcloud_todevice(self):
        """
            Test that point clouds can be sent to a device, and back to the host
        """
        for pc in self.pcs:
            hostpc = self.tohost(pc)
            devicetype, deviceid = (int(v) for v in pc.__dlpack_device__())

            # Without, then with normals
            for cloud in [pyponca.PointCloud(hostpc), pyponca.PointCloud(hostpc, hostpc)]:
                self.assertEqual((cloud.device_type, cloud.device_id), (1, 0)) # kDLCPU

                # With a device type or an array
                for devicecloud in [cloud.toDevice(pc), cloud.toDevice(devicetype, deviceid)]:
                    self.assertIsInstance(devicecloud, pyponca.PointCloud)
                    self.assertIsInstance(devicecloud.object, cloud._cls)
                    self.assertEqual((devicecloud.device_type, devicecloud.device_id), (devicetype, deviceid))

                    # From the device: to the same device, then back to the host
                    self.assertEqual(devicecloud.toDevice(pc).device_type, devicetype)
                    self.assertEqual(devicecloud.toDevice(1).device_type, 1)

                # The original point cloud is left on the host
                self.assertEqual(cloud.device_type, 1)

    def test_kdtree(self):
        """
            Test that we can build a KDtree and that mangling works properly
        """
        classes = [_pyponca.KdTree2dPN, _pyponca.KdTree3dPN, _pyponca.KdTree2fPN, _pyponca.KdTree3fPN]
        for pc, cls in zip(self.pcs, classes):
            # Trees can only be built on the host
            with self.assertRaisesUnless(self.onHost):
                pyponca.KDTree(pc)

            # Test we can build a kdtree
            hosttree = pyponca.KDTree(self.tohost(pc))
            self.assertIsInstance(hosttree.object, cls)
            self.assertEqual((hosttree.device_type, hosttree.device_id), (1, 0)) # kDLCPU

            # Then send it to the device, with a device type or an array
            devicetype, deviceid = (int(v) for v in pc.__dlpack_device__())
            for tree in [hosttree.toDevice(pc), hosttree.toDevice(devicetype, deviceid)]:
                self.assertIsInstance(tree, pyponca.KDTree)
                self.assertIsInstance(tree.object, cls)
                self.assertEqual((tree.device_type, tree.device_id), (devicetype, deviceid))

                # From the device: to the same device, then back to the host
                self.assertEqual(tree.toDevice(pc).device_type, devicetype)
                self.assertEqual(tree.toDevice(1).device_type, 1)

            # The original tree is left on the host
            self.assertEqual(hosttree.device_type, 1)

    def test_kdtree_roundtrip(self):
        """
            Test that a tree sent to the device and back to the host returns the same neighbors
        """
        for pc in self.pcs:
            hostpc   = self.tohost(pc)
            hosttree = pyponca.KDTree(hostpc)
            backtree = hosttree.toDevice(pc).toDevice(1)

            radii = np.full((self.N,), 0.5, dtype=hostpc.dtype)
            np.testing.assert_array_equal(np.from_dlpack(backtree.rangeNeighbors(hostpc, radii)),
                                          np.from_dlpack(hosttree.rangeNeighbors(hostpc, radii)))

    def test_kdtree_query(self):
        """
            Test queries can be performed on kdtrees. They are only performed on the host
        """
        for pc in self.pcs:
            radii = self.randomtensor((self.N,), self.tohost(pc).dtype)

            # Tree stored on the device the tests run on, then on the host
            for kdtree in [self.kdtree(pc), pyponca.KDTree(self.tohost(pc))]:
                with self.assertRaisesUnless(self.onHost):
                    kdtree.rangeNeighbors(pc[0], 1)
                with self.assertRaisesUnless(self.onHost):
                    kdtree.rangeNeighbors(pc, radii)

    def test_cuda_available(self):
        """
            Test that compute objects tell if they can run on CUDA devices
        """
        for co in pyponca.ComputeObjectList:
            isCudaAvailable = pyponca.__dict__[co].isCudaAvailable
            self.assertIsInstance(isCudaAvailable, bool)
            if not _pyponca.internal._WithCUDA:
                self.assertFalse(isCudaAvailable, co)

    def test_compute(self):
        """
            Test that we can attach pointclouds and filter to a compute object, with every kind of data
        """
        for co in pyponca.ComputeObjectList:
            for pc in self.pcs:
                for flt in pyponca.Filters:
                    ts = self.randomtensor((self.N,), self.tohost(pc).dtype)

                    for data in [pc, pyponca.PointCloud(pc, pc), self.kdtree(pc), self.kdtree(pc, False)]:
                        object = pyponca.__dict__[co]()
                        object.setNeighborFilter(pc, ts, flt, pc)
                        object.addComputation(pyponca.Computation.NONE)
                        try:
                            result = object.compute(data)
                        except NotImplementedError:
                            # We let notimplementederror that may come from classes that only support 3d
                            continue
                        except RuntimeError:
                            # Classes that can not run on the device
                            self.assertFalse(self.supports(co), co)
                            continue

                        self.assertTrue(self.supports(co), co)
                        self.assertEqual(self.tohost(result.fitStatus).shape, (self.N,))
                        self.assertEqual(result.fitStatus.__dlpack_device__(), pc.__dlpack_device__())

    def test_compute_func(self):
        """
            Test for potential memory errors and settings in a classical loop
        """
        for co in pyponca.ComputeObjectList:
            for pc in [self.pc3d, self.pc3f]:
                ts = self.randomtensor((self.N,), self.tohost(pc).dtype)

                cls = pyponca.__dict__[co]

                # Helper function to try with different data
                def test(data):
                    obj = cls()
                    obj.setNeighborFilter(pc, ts, normals = pc)
                    obj.addComputation(pyponca.Computation.NONE)
                    with self.assertRaisesUnless(self.supports(co)):
                        obj.compute(data)

                test(pc)
                test(pyponca.PointCloud(pc))
                test(self.kdtree(pc))

    def test_device_mismatch(self):
        """
            Test that host and device arrays can not be mixed
        """
        for pc in self.pcs:
            hostpc = self.tohost(pc)
            object = pyponca.MeanPosition()
            object.setNeighborFilter(hostpc, np.ones((self.N,), dtype=hostpc.dtype))
            object.addComputation(pyponca.Computation.MEAN_POSITION)
            with self.assertRaisesUnless(self.onHost):
                object.compute(pc)

    def test_no_normals(self):
        """
            Test that missing normals are set to zero
        """
        pc = self.pc3d
        object = pyponca.MeanNormal()
        # The radius covers the whole cloud, so that every fit has neighbors
        object.setNeighborFilter(pc, self.todevice(np.full((1,), 10, dtype=np.float64)))
        object.addComputation(pyponca.Computation.MEAN_NORMAL)
        result = self.tohost(object.compute(pc).results[0])
        np.testing.assert_array_equal(result, np.zeros((self.N, 3)))


@unittest.skipUnless(torch is not None and torch.cuda.is_available() and pyponca.WithCuda,
                     "Requires torch with a CUDA device and pyponca built with CUDA")
class TestAPICUDA(TestAPI):
    """
        Same tests, with every array sent to the GPU
    """
    onHost = False

    def todevice(self, x):
        return torch.as_tensor(x, device="cuda")

    def tohost(self, x):
        # Convert to numpy as we use the np testing suite. 
        return torch.from_dlpack(x).cpu().numpy()


if __name__ == "__main__":
    unittest.main()
