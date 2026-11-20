# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

# PyPonca documentation example (doc/src/ponca_pyponca.mdoc): APSS projection and curvature.
# The code between the snippet markers is shown in the documentation.

# [apss]
import numpy as np
import pyponca

pos     = np.random.normal(size=(5000, 3))
pos    /= np.linalg.norm(pos, axis=1, keepdims=True)  # Points on the unit sphere
normals = pos.copy()

cloud = pyponca.PointCloud(pos, normals)
tree  = pyponca.KDTree(cloud)

queries = pos[::10] * 1.05                            # Some locations near the surface
radii   = np.full(queries.shape[0], 0.2)

fit = pyponca.APSS()
fit.setNeighborFilter(queries, radii, pyponca.Filters.SMOOTH_WEIGHT)
fit.addComputation(pyponca.Computation.PROJECTION, queries)  # Project each query on its own fit
fit.addComputation(pyponca.Computation.GLS)

projection, gls = fit.compute(tree).results

projection = np.from_dlpack(projection)               # Shape (M, 3)
gls        = np.from_dlpack(gls)                      # Shape (M, 2 + 3): tau, kappa, eta
kappa      = gls[:, 1]                                # Curvature
# [apss]
