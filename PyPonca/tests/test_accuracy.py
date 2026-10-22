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

    y = np.empty((3 * N - 2, 2 * D - 1), dtype=x.dtype)
    y[3 * (N - 1)::-3, 2 * (D - 1)::-2] = x

    # Check that the striding works as intented !
    np.testing.assert_allclose(x, y[3 * (N-1)::-3, 2 * (D - 1)::-2])
    return y

def readarray(testcase, name):
    """
        Reads an array from the test cases
    """
    data = None
    if name in testcase:
        data = np.asarray(testcase[name]["data"])
        data = data.reshape((testcase[name]["N"], testcase[name]["D"]))
        data = data.squeeze()
        data = data.astype(np.float64)
    
    return data



class TestAccuracy(unittest.TestCase):
    def setUp(self):
        with open(os.environ["TEST_CASES"], "rt") as file:
            content = file.read()
        content = content.replace("-nan", "NaN")
        content = content.replace("nan", "NaN")
        self.expected = json.loads(content)

        # Load main pointclouds
        self.pos     = readarray(self.expected["pointcloud"], "pos")
        self.normals = readarray(self.expected["pointcloud"], "normals")
        self.points  = pyponca.PointCloud(self.pos, self.normals)

        self.analysisLocation = readarray(self.expected, "analysisLocation")
        self.analysisScale    = readarray(self.expected, "analysisScale")
        self.analysisNormals  = readarray(self.expected, "analysisNormals")

    def runtestcase(self, pointcloud, method, function, input, status, signed, result):
        """
            Run a testcase using the provided inputs and expected results
        """
        if method in pyponca.ComputeObjectList:
            # Disable tests for CNC. It uses some random sampling within it's
            # triangle generation algorithm and there are no API to set seed for this. 
            # We can hope that the results are correct (due to binding code being
            # a for loop over) if it passes for all other cases. 
            if "CNC" in method:
                return 
            
            if function in pyponca.Computation.__members__: 
                object = pyponca.__dict__[method]()
                id = pyponca.Computation[function]

                object.setNeighborFilter(self.analysisLocation, self.analysisScale, normals = self.analysisNormals)
                object.addComputation(id, input)

                pyresult = object.compute(pointcloud)
                pystatus = np.from_dlpack(pyresult.fitStatus)
                pyresult = np.from_dlpack(pyresult.results[0])

                if not signed:
                    pyresult = np.abs(pyresult)
                    result = np.abs(result)

                # Test that status match (stable and undefined)
                np.testing.assert_allclose(pystatus, status)
                mask = (status == pyponca.FitResult.STABLE)
                np.testing.assert_allclose(pyresult[mask], result[mask], atol=1e-8, err_msg=f"Failed: {method} / {id} / {input.shape if input is not None else ""}")
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
            input    = readarray(run, "input")
            result   = readarray(run, "result")
            status   = readarray(run, "status")
            sign     = bool(run["signed"]) if "signed" in run else False

            self.runtestcase(self.points, method, function, input, status, sign, result)

    def test_pointcloud_strided(self):
        """
            Test that the binding returns the correct results when given
            a strided array. See stridearray for more information about the stride
            in question
        """
        for run in self.expected["runs"]:
            method   = run["method"]
            function = run["function"]
            input    = readarray(run, "input")
            result   = readarray(run, "result")
            status   = readarray(run, "status")
            sign     = bool(run["signed"]) if "signed" in run else False

            N = self.pos.shape[0]
            D = self.pos.shape[1]
            newpos     = stridearray(self.pos)[3*(N-1)::-3, 2 * (D - 1)::-2]
            newnormals = stridearray(self.normals)[3*(N-1)::-3, 2 * (D - 1)::-2]
            newpoints  = pyponca.PointCloud(newpos, newnormals)

            self.runtestcase(newpoints, method, function, input, status, sign, result)

    def test_kdtree(self):
        """
            Test that the binding returns the correct results with a kdtree
        """
        for run in self.expected["runs"]:
            method   = run["method"]
            function = run["function"]
            input    = readarray(run, "input")
            result   = readarray(run, "result")
            status   = readarray(run, "status")
            sign     = bool(run["signed"]) if "signed" in run else False

            # Dense then sparse
            self.runtestcase(pyponca.KDTree(self.points, False), method, function, input, status, sign, result)
            self.runtestcase(pyponca.KDTree(self.points, True) , method, function, input, status, sign, result)

    def test_kdtree_strided(self):
        """
            Test that the binding returns the correct results with a kdtree
        """
        for run in self.expected["runs"]:
            method   = run["method"]
            function = run["function"]
            input    = readarray(run, "input")
            result   = readarray(run, "result")
            status   = readarray(run, "status")
            sign     = bool(run["signed"]) if "signed" in run else False

            N = self.pos.shape[0]
            D = self.pos.shape[1]
            newpos     = stridearray(self.pos)[3*(N-1)::-3, 2 * (D - 1)::-2]
            newnormals = stridearray(self.normals)[3*(N-1)::-3, 2 * (D - 1)::-2]
            newpoints  = pyponca.PointCloud(newpos, newnormals)
    
            # Dense then sparse
            self.runtestcase(pyponca.KDTree(newpoints, False), method, function, input, status, sign, result)
            self.runtestcase(pyponca.KDTree(newpoints, True) , method, function, input, status, sign, result)


if __name__ == "__main__":
    unittest.main()
