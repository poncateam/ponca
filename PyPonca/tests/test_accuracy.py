import os
import json
import unittest
import warnings
import numpy as np
from pathlib import Path

import pyponca

try:
    import torch
except ImportError:
    torch = None

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
    # Whether the tests run on the host, where every method is available
    onHost = True

    def todevice(self, x):
        """
            Sends an array to the device the tests run on
        """
        return x

    def tohost(self, x):
        """
            Converts an array returned by pyponca to a numpy array
        """
        return np.from_dlpack(x)

    def strided(self, x):
        """
            Sends an array to the device the tests run on, as a non contiguous view.
            See stridearray for more information about the stride in question
        """
        N, D = x.shape
        return self.todevice(stridearray(x))[3 * (N - 1)::-3, 2 * (D - 1)::-2]

    def setUp(self):
        with open(os.environ["TEST_CASES"], "rt") as file:
            content = file.read()
        content = content.replace("-nan", "NaN")
        content = content.replace("nan", "NaN")
        self.expected = json.loads(content)

        # Load main pointclouds
        self.pos     = readarray(self.expected["pointcloud"], "pos")
        self.normals = readarray(self.expected["pointcloud"], "normals")
        self.points  = pyponca.PointCloud(self.todevice(self.pos), self.todevice(self.normals))

        self.analysisLocation = readarray(self.expected, "analysisLocation")
        self.analysisScale    = readarray(self.expected, "analysisScale")
        self.analysisNormals  = readarray(self.expected, "analysisNormals")

    def stridedpoints(self):
        """
            Returns the point cloud built from strided arrays. See stridearray
            for more information about the stride in question
        """
        return pyponca.PointCloud(self.strided(self.pos), self.strided(self.normals))

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

            # Methods that can not run on the device
            if not (self.onHost or pyponca.__dict__[method].isCudaAvailable):
                return

            if function in pyponca.Computation.__members__: 
                object = pyponca.__dict__[method]()
                id = pyponca.Computation[function]

                object.setNeighborFilter(self.todevice(self.analysisLocation), self.todevice(self.analysisScale),
                                         normals = self.todevice(self.analysisNormals))
                object.addComputation(id, self.todevice(input) if input is not None else None)

                pyresult = object.compute(pointcloud)
                pystatus = self.tohost(pyresult.fitStatus)
                pyresult = self.tohost(pyresult.results[0])

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

            self.runtestcase(self.stridedpoints(), method, function, input, status, sign, result)

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

            newpoints = self.stridedpoints()

            # Dense then sparse
            self.runtestcase(pyponca.KDTree(newpoints, False), method, function, input, status, sign, result)
            self.runtestcase(pyponca.KDTree(newpoints, True) , method, function, input, status, sign, result)


@unittest.skipUnless(torch is not None and torch.cuda.is_available() and pyponca.WithCuda,
                     "Requires torch with a CUDA device and pyponca built with CUDA")
class TestAccuracyCUDA(TestAccuracy):
    """
        Same tests, with every array sent to the GPU
    """
    onHost = False

    def todevice(self, x):
        return torch.as_tensor(x, device="cuda") if x is not None else None

    def tohost(self, x):
        # Convert to numpy as we use the np testing suite. 
        return torch.from_dlpack(x).cpu().numpy()

    def strided(self, x):
        # Torch does not support negative strides: stridearray(x[::-1, ::-1])[::3, ::2] == x
        return self.todevice(stridearray(x[::-1, ::-1]))[::3, ::2]


if __name__ == "__main__":
    unittest.main()
