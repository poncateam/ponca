/*
This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/
#include <nanobind/nanobind.h>
#include "Mangling.h"
#include "pypoint.h"

namespace nb = nanobind;
using namespace nb::literals;

/**
 * \brief Register Mangling functions
 *
 * These function are made available in python code for coherence
 * between the two codes.
 *
 * \param m The module. Expected to be the internal module.
 */
void RegisterManglingUtils(nanobind::module_& m)
{
    namespace nb = nanobind;

    // No constraints on the array, should work with anything
    m.def("_mangleArray", [](const nb::ndarray<>& arr) { return MangleArray(arr); });
    m.attr("_PointName") = nb::cast(MangledPointName);
}

/**
 * \brief Register a PyPoncaPointCloud depending on the point type
 *
 * \param m The module to register the class within
 */
template <typename _Scalar, unsigned int _Dim>
void RegisterPyPointCloud(nb::module_& m)
{
    using PointCloud = PyPointCloud<_Scalar, _Dim>;
    using Point      = typename PointCloud::Point;

    const std::string pointName = ManglePoint<Point>() + "PN";
    const std::string className = "PointCloud" + pointName;

    // Only the constructor and the mangledName are exposed. Nothing else should
    // be necessary in client code.
    // MangledName is necessary for dispatching to the correct PointCloud type
    auto cls = nb::class_<PointCloud>(m, className.c_str());

    cls.def(nb::init<const PyVectorArray<Point>&>(), nb::arg("pos").noconvert());
    cls.def(nb::init<const PyVectorArray<Point>&, const PyVectorArray<Point>&>(), nb::arg("pos").noconvert(),
            nb::arg("normal").noconvert());
}

/**
 * \brief Register all point clouds
 *
 * \param m The module to register the class within
 */
void RegisterPointClouds(nb::module_& m)
{
    RegisterPyPointCloud<double, 3>(m);
    RegisterPyPointCloud<float, 3>(m);
    RegisterPyPointCloud<double, 2>(m);
    RegisterPyPointCloud<float, 2>(m);
}

/**
 * \brief Register class, functions and enums for the common module
 *
 * \param m The main module
 * \param internal An internal module reserved for the library
 */
void RegisterCommon(nb::module_& m, nb::module_& internal)
{
    RegisterManglingUtils(internal);
    RegisterPointClouds(m);
}

