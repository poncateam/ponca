/*
This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <Ponca/Ponca>

namespace nb = nanobind;
#include "Common/RegisterCommon.h"
#include "Fitting/RegisterFitting.h"
#include "SpatialPartitioning/RegisterSpatialPartitioning.h"

NB_MODULE(_pyponca, m)
{
    auto internal = m.def_submodule("internal");

#ifdef __CUDACC__
    internal.attr("_WithCUDA") = true;
#else
    internal.attr("_WithCUDA") = false;
#endif

    RegisterCommon(m, internal);
    RegisterFitting(m, internal);
    RegisterSpatialPartitioning(m, internal);
}

