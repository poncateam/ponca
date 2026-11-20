# 
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#
from . import _pyponca

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

WithCuda = _pyponca.internal._WithCUDA
