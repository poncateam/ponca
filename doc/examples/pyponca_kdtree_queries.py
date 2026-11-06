# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

# PyPonca documentation example (doc/src/ponca_pyponca.mdoc): neighbor queries with a KDTree.
# The code between the snippet markers is shown in the documentation.

import numpy as np
import pyponca

points = np.random.uniform(-1, 1, (10000, 3))

# [kdtree_queries]
tree = pyponca.KDTree(points)

# Single query: returns a 1D array of indices
ids = np.from_dlpack(tree.rangeNeighbors(points[0], 0.1))

# Multiple queries: returns a (M, K) array, padded with -1
ids = np.from_dlpack(tree.rangeNeighbors(points[:100], np.array([0.1])))
neighbors_of_3 = ids[3][ids[3] >= 0]
# [kdtree_queries]
