/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/
#pragma once

#include <nanobind/nanobind.h>
#include "FittingUtils.h"

/**
 * \brief List of computation allowed
 */
enum class Computation
{
    DRY,
    POTENTIAL,
    PROJECTION
};

/**
 * \brief Register the computation enumeration defined above
 * 
 * \param m Nanobind module
 */
void RegisterComputations(nanobind::module_& m)
{
    nb::enum_<Computation>(m, "Computation")
        .value("DRY", Computation::DRY)
        .value("POTENTIAL", Computation::POTENTIAL)
        .value("PROJECTION", Computation::PROJECTION)
        .export_values();
}

/**
 * \brief Returns the expected shape given the input and the type of computation
 *
 * This function allows to pre allocate the result array.
 *
 * The first element of the returned can be any value and will be overidden by the
 * PyComputeObject instance and filled with the number of filter centers.
 *
 * This function is only called on host code.
 *
 * \param id The id of the computation
 * \param array The input array. The shape of this array may indicates no inputs or an array of inputs for each center.
 */
template <typename _PyCO>
inline std::vector<size_t> ComputeOutputDimension(Computation id, const nb::ndarray<>& array)
{
    using Point = typename _PyCO::Point;
    using CO    = typename _PyCO::ComputeObject;

    switch (id)
    {
    case Computation::DRY:
        return std::vector<size_t>{0, 0};
    case Computation::POTENTIAL:
        if constexpr (Ponca::ProvidesImplicitPrimitive<CO>)
        {
            using Result = std::remove_cv_t<decltype(CO{}.potential())>;
            return ComputeVectorDim<Point, Result>(array);
        }
        break;
    case Computation::PROJECTION:
        if constexpr (Ponca::ProvidesProjectionOperator<CO>)
        {
            using Result = std::remove_cv_t<decltype(CO{}.project({}))>;
            return ComputeVectorDim<Point, Result>(array);
        }
        break;
    default:
        break;
    }

    throw std::runtime_error("Invalid combination of object / computation (" + std::to_string(static_cast<int>(id)) +
                             ") requested");
    return {};
}

/**
 * \brief Perform the computation
 *
 * This function must be compatible with all device the code may want to run.
 * The out array is already preallocated to the correct size, as prescribed by
 * ComputeOutputDimension.
 *
 * \param id The id of the computation
 * \param object The fitted computed object
 * \param i The index of the center being processed
 * \param in The input array for this computation.
 * \param out The output array
 */
template <typename CO, typename Scalar = typename CO::Scalar>
inline void ExtractComputation(size_t id, CO& object, size_t i, const DeviceArrayView<Scalar>& in,
                                               DeviceArrayView<Scalar>& out)
{
    Computation computation = static_cast<Computation>(id);
    switch (computation)
    {
    case Computation::DRY:
        return;
    case Computation::POTENTIAL:
        if constexpr (Ponca::ProvidesImplicitPrimitive<CO>)
        {
            if (in.data != nullptr)
                RunVectorMethod<CO>(object, MAKE_VECTOR_FUNCTION(potential), in, out, i);
            else
                Write(out, i * out.stride[0], object.potential());
        }
        break;
    case Computation::PROJECTION:
        if constexpr (Ponca::ProvidesProjectionOperator<CO>)
            RunVectorMethod<CO>(object, MAKE_VECTOR_FUNCTION(project), in, out, i);
        break;
    };
}
