# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

# PyPonca documentation example (doc/src/ponca_pyponca.mdoc): using PyTorch tensors instead of NumPy arrays.
# The code between the snippet markers is shown in the documentation.

# [torch]
import torch
import pyponca

points = torch.rand((1000, 3), dtype=torch.float32)

fit = pyponca.SphereFit()
fit.setNeighborFilter(points, torch.tensor([0.2]))
fit.addComputation(pyponca.Computation.ALGEBRAIC_SPHERE)
sphere = torch.from_dlpack(fit.compute(points).results[0])  # Shape (1000, 2 + 3)
# [torch]
