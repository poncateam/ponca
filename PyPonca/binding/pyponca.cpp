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

    RegisterCommon(m, internal);
    RegisterFitting(m, internal);
    RegisterSpatialParitioning(m, internal);
}

