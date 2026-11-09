# 
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#
from . import _pyponca
from enum import StrEnum

def _parseDevice(device, device_id):
    """
        Returns the (device type, device id) pair given to toDevice functions: device is either a DLPack
        device type (e.g. 2 for CUDA), or an object exposing __dlpack_device__ (e.g. a numpy array or a
        torch tensor), whose device is used
    """
    if hasattr(device, "__dlpack_device__"):
        device, device_id = device.__dlpack_device__()
    return int(device), int(device_id)

# Dispatcher for point cloud
class PointCloud:
    """
        Main ponintcloud object that serves as a dispatch
        class magnles the array and retrieve the corresponding
        Ponca point cloud structure
    """
    def __init__(self, *args):
        self.N = 0
        self.dim = 0
        self._mangledName = "unknown"
        self._cls = object

        if len(args) == 1 or len(args) == 2:
            if args[0].ndim != 2:
                raise RuntimeError("Only 2D arrays are supported !")
            
            if len(args) == 2:
                if args[0].shape != args[1].shape:
                    raise RuntimeError("Position and normal must have the same shape !")
            
            self.object = self._pyponcapointcloud(*args)
        else:
            raise RuntimeError("At most 2 arrays of shape (N, D) are supported")
    
    def _pyponcapointcloud(self, *args):
        mangledArray = _pyponca.internal._mangleArray(args[0])

        clsname = "PointCloud" + mangledArray + _pyponca.internal._PointName
        if clsname not in _pyponca.__dict__:
            raise NotImplementedError(f"PointCloud does not have specialization for {args[0].dtype} arrays")

        self._cls  = _pyponca.__dict__[clsname]
        self._mangledName = mangledArray + _pyponca.internal._PointName
        return self._cls(*args)

    def toDevice(self, device, device_id = 0):
        """
            Returns a copy of the point cloud (positions and normals, with their strides) stored on a device.
            The point cloud may be stored on any device. The new point cloud owns its arrays.

            see _parseDevice
        """
        # Not built with __init__, which expects arrays
        cloud = PointCloud.__new__(PointCloud)
        cloud.__dict__.update(self.__dict__)
        cloud.object = self.object.toDevice(*_parseDevice(device, device_id))
        return cloud

    @property
    def device_type(self):
        """ DLPack device type of the arrays of the point cloud """
        return self.object.device_type

    @property
    def device_id(self):
        """ Index of the device the arrays of the point cloud are stored on """
        return self.object.device_id

class KDTree:
    """
        Main kdtree object that serves as a dispatch
        class magnles the array and retrieve the corresponding
        Ponca kdtree structure
    """

    def __init__(self, pc, isDense = True):
        self.N = 0
        self.dim = 0
        self._mangledName = "unknown"
        self._cls = object

        if not isinstance(pc, PointCloud):
            mangledArray = _pyponca.internal._mangleArray(pc)
            clsname = "KdTree" + mangledArray + _pyponca.internal._PointName
            self._mangledName = mangledArray + _pyponca.internal._PointName
        else:
            clsname = "KdTree" + pc._mangledName
            self._mangledName = pc._mangledName
            pc = pc.object

        if clsname not in _pyponca.__dict__:
            raise NotImplementedError(f"PointCloud does not have specialization for {self._mangledName} arrays")

        self._cls = _pyponca.__dict__[clsname]
        self.object = self._cls(pc, isDense)

    def toDevice(self, device, device_id = 0):
        """
            Returns a copy of the tree, and of the point cloud it was built from, stored on a device.
            The tree may be stored on any device (trees are always built on the CPU). When it is already
            stored on that device, the returned tree shares its buffers.

            see _parseDevice
        """
        # We do not want to
        tree = KDTree.__new__(KDTree)
        tree.__dict__.update(self.__dict__)
        tree.object = self.object.toDevice(*_parseDevice(device, device_id))
        return tree

    def __getattr__(self, name):
        return getattr(self.object, name)
        
class Filters(StrEnum):
    """
        Filter mangling mapping
    """
    SMOOTH_WEIGHT = "SW"
    CONSTANT_WEIGHT = "CW"
    NO_WEIGHT = "NW"

# Following methods are dispatch functions

def __initComputeObject(self, name: str):
    """
        Ctor of to-be-dispatch compute object
    """
    self.name = name
    self.nf = Filters.SMOOTH_WEIGHT
    self._computationList = []

    self._object = None
    self._pos = None
    self._radius = None

# Creates the object
def __createComputeObject(self, data):
    """
        Helper function to dispatch to the correct underlying class
    """
    clsname = self.name + data._mangledName + str(self.nf)
    if clsname not in _pyponca.__dict__:
        # With filter failed, try without:
        clsname = self.name + data._mangledName
        if clsname not in _pyponca.__dict__:
            raise NotImplementedError(f"Compute object {self.name} do not have a specialization for {data._mangledName + str(self.nf)} point clouds")

    return _pyponca.__dict__[clsname]()
    
def __dispatchCompute(self, data):
    """
        Dispatch the attach function to the correct class
    """
    if isinstance(data, KDTree):
        pass
    elif not isinstance(data, PointCloud):
        # Convert array to pointcloud !
        data = PointCloud(data)

    self.object = __createComputeObject(self, data)

    # May fail with a pretty good error message generated by nanobind 
    # if dimension or type does not match 
    for computation in self._computationList:
        self.object.addComputation(*computation)

    self.object.setNeighborFilter(self._pos, self._radius, self._normal)
    return self.object.compute(data.object)

def __dispatchSetNeighborFilter(self, pos, radius, nf = Filters.SMOOTH_WEIGHT, normals = None):
    """
        Dispatch the setNeighborFilter call. Also accept the neighbor filter
        as input
    """
    # Stores information for further dispatch
    # We leave the __dispatchCompute function the responsibility
    # of calling setNeighborFilter, which saves codes as we do not
    # know yet if the underlying object will support these parameters
    self.nf = nf

    self._pos = pos
    self._radius = radius
    self._normal = normals


def __dispatchAddComputation(self, id, input = None):
    """
        Adds a computation to the object
    """
    self._computationList.append((id, input))

def __isCudaAvailable(name: str):
    """
        Tells if a compute object can run on CUDA devices. This does not depend
        on the specialization, so the first one found is used
    """
    for clsname, cls in _pyponca.__dict__.items():
        # Specializations are named after the object
        if clsname.startswith(name):
            return cls.isCudaAvailable
    return False

# Forward types to global module 
ComputeObjectList = _pyponca.ComputeObjectList
Computation = _pyponca.Computation
FitResult = _pyponca.FitResult
WithCuda = _pyponca.internal._WithCUDA

# Dynamically create object based on what was exposed by the binding. 
for co in ComputeObjectList:
    newco = type(co, (), {
        "__init__": lambda self, co=co: __initComputeObject(self, co),
        "compute": __dispatchCompute,
        "setNeighborFilter": __dispatchSetNeighborFilter,
        "addComputation": __dispatchAddComputation,
        "isCudaAvailable": __isCudaAvailable(co),
        # Dynamically find the function within the subobject
        # "__getattr__": __dispatchgetattr
    })

    # Register the new class within the pyponca namespace
    globals()[co] = newco