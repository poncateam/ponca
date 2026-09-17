import os
import json
import unittest
import warnings
import numpy as np
from pathlib import Path

import pyponca

def stridearray(x):
    """
        Returns a new array such that seld.stridearray(x)[3*(N-1)::-3, ...] == x 

        For now pyponca does not support non contiguous arrays in the second dimension !
    """
    N, D = x.shape

    y = np.empty((3 * N - 2, D), dtype=x.dtype)
    y[3 * (N - 1)::-3] = x

    # Check that the striding works as intented !
    np.testing.assert_allclose(x, y[3 * (N-1)::-3, :])
    return y


class TestAccuracy(unittest.TestCase):
    def setUp(self):
        with open(os.environ["TEST_CASES"], "rt") as file:
            self.expected = json.load(file)

        # Load main pointclouds
        N = self.expected["pointcloud"]["pos"]["N"]
        D = self.expected["pointcloud"]["pos"]["D"]
        self.pos     = np.array(self.expected["pointcloud"]["pos"]["data"]    ).reshape((N, D))
        self.normals = np.array(self.expected["pointcloud"]["normals"]["data"]).reshape((N, D))
        self.points  = pyponca.PointCloud(self.pos, self.normals)

        M = self.expected["analysisLocation"]["N"]
        self.analysisLocation = np.array(self.expected["analysisLocation"]["data"]).reshape((M, D))
        self.analysisScale    = np.array(self.expected["analysisScale"]   ["data"]).reshape((M))

    def runtestcase(self, pointcloud, method, function, input, signed, result):
        """
            Run a testcase using the provided inputs and expected results
        """
        if method in pyponca.ComputeObjectList:
            if function in pyponca.Computation.__members__: 
                object = pyponca.__dict__[method]()
                id = pyponca.Computation[function]

                object.setNeighborFilter(self.analysisLocation, self.analysisScale)
                object.addComputation(id, input)
                pyresult = np.from_dlpack(object.compute(pointcloud)[0].data)

                # The abs here is for signed / unsigned distances. This is sufficient, but we
                # might want to check from the class whether this is necessary
                if signed:
                    np.testing.assert_allclose(pyresult, result, err_msg=method)    
                else:
                    np.testing.assert_allclose(np.abs(pyresult), np.abs(result), err_msg=method)    
            else:
                warnings.warn(f"Function not found: {function}")
        else:
            warnings.warn(f"Method not found: {method}")

    def test_pointcloud(self):
        """
            Test that the binding returns the correct results
        """
        for run in self.expected["runs"]:
            method   = run["method"]
            function = run["function"]
            result   = run["result"]["data"]
            input    = np.array(run["input"]) if "input" in run else None
            sign     = bool(run["signed"]) if "signed" in run else False

            self.runtestcase(self.points, method, function, input, sign, result)

    def test_pointcloud_strided(self):
        """
            Test that the binding returns the correct results when given
            a strided array. See stridearray for more information about the stride
            in question
        """
        for run in self.expected["runs"]:
            method   = run["method"]
            function = run["function"]
            result   = run["result"]["data"]
            input    = np.array(run["input"]) if "input" in run else None
            sign     = bool(run["signed"]) if "signed" in run else False

            N = self.pos.shape[0]
            newpos     = stridearray(self.pos)[3*(N-1)::-3, :]
            newnormals = stridearray(self.normals)[3*(N-1)::-3, :]
            newpoints  = pyponca.PointCloud(newpos, newnormals)

            # Simpler stride here: simply a reverse
            if input is not None:
                input = input[::-1].ascontiguousarray()[::-1]

            self.runtestcase(newpoints, method, function, input, sign, result)

    def test_kdtree(self):
        """
            Test that the binding returns the correct results with a kdtree
        """
        for run in self.expected["runs"]:
            method   = run["method"]
            function = run["function"]
            result   = run["result"]["data"]
            input    = np.array(run["input"]) if "input" in run else None
            sign     = bool(run["signed"]) if "signed" in run else False

            # Dense then sparse
            self.runtestcase(pyponca.KDTree(self.points, False), method, function, input, sign, result)
            self.runtestcase(pyponca.KDTree(self.points, True) , method, function, input, sign, result)

    def test_kdtree_strided(self):
        """
            Test that the binding returns the correct results with a kdtree
        """
        for run in self.expected["runs"]:
            method   = run["method"]
            function = run["function"]
            result   = run["result"]["data"]
            input    = np.array(run["input"]) if "input" in run else None
            sign     = bool(run["signed"]) if "signed" in run else False

            N = self.pos.shape[0]
            newpos     = stridearray(self.pos)[3*(N-1)::-3, :]
            newnormals = stridearray(self.normals)[3*(N-1)::-3, :]
            newpoints  = pyponca.PointCloud(newpos, newnormals)

            # Simpler stride here: simply a reverse
            if input is not None:
                input = input[::-1].ascontiguousarray()[::-1]
    
            # Dense then sparse
            self.runtestcase(pyponca.KDTree(newpoints, False), method, function, input, sign, result)
            self.runtestcase(pyponca.KDTree(newpoints, True) , method, function, input, sign, result)


if __name__ == "__main__":
    unittest.main()
