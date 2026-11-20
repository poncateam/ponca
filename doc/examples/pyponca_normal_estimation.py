# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

# PyPonca documentation example (doc/src/ponca_pyponca.mdoc): normal estimation with a plane fit.
# The code between the snippet markers is shown in the documentation.

# [normal_estimation]
import numpy as np
import pyponca

points = np.random.uniform(-1, 1, (10000, 3))  # float64 point cloud of shape (N, 3)
radii  = np.array([0.1])                       # A single radius, broadcast to every center

fit = pyponca.CovariancePlaneFit()
fit.setNeighborFilter(points, radii)           # One fit centered on each point
fit.addComputation(pyponca.Computation.PRIMITIVE_GRADIENT)
fit.addComputation(pyponca.Computation.MEAN_POSITION)

output = fit.compute(pyponca.KDTree(points))

normals     = np.from_dlpack(output.results[0])  # Shape (N, 3)
barycenters = np.from_dlpack(output.results[1])  # Shape (N, 3)

# Fits with too few neighbors are UNDEFINED: only keep the stable ones
stable  = np.from_dlpack(output.fitStatus) == pyponca.FitResult.STABLE
normals = normals[stable]
# [normal_estimation]
