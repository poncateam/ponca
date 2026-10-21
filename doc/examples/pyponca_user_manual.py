# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

# Code fragments of the PyPonca user manual (doc/src/ponca_pyponca.mdoc), runnable as a whole.
# The code between each pair of snippet markers is shown in the documentation.

import numpy as np
import pyponca

pos     = np.random.uniform(-1, 1, (1000, 3))
normals = pos / np.linalg.norm(pos, axis=1, keepdims=True)

# [pointcloud]
cloud = pyponca.PointCloud(pos)           # Positions only
cloud = pyponca.PointCloud(pos, normals)  # Positions and normals
# [pointcloud]

# [kdtree]
tree = pyponca.KDTree(cloud)         # From a PointCloud (or a raw (N, D) array)
tree = pyponca.KDTree(cloud, False)  # Second argument selects the storage, see below
# [kdtree]

fit = pyponca.APSS()
fit.setNeighborFilter(pos, np.array([0.2]))

# [results]
fit.addComputation(pyponca.Computation.POTENTIAL)
fit.addComputation(pyponca.Computation.PRIMITIVE_GRADIENT)
output = fit.compute(tree)

potential, gradient = [np.from_dlpack(r) for r in output.results]
status = np.from_dlpack(output.fitStatus)
# [results]

# [fit_status]
stable   = status == pyponca.FitResult.STABLE  # Boolean mask of shape (M,)
gradient = gradient[stable]
# [fit_status]
