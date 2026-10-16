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
    NONE,
    MEAN_POSITION, 
    MEAN_NORMAL,
    POTENTIAL,
    PRIMITIVE_GRADIENT, 
    PROJECTION,
    PRINCIPAL_CURVATURE,
    ALGEBRAIC_SPHERE,
    GLS, 
    GLS_NORMALIZED,  
    QUADRATIC_HEIGHTFIELD_MAP
};

/**
 * \brief Register the computation enumeration defined above
 *
 * \param m Nanobind module
 */
void RegisterComputations(nanobind::module_& m)
{
    nb::enum_<Computation>(m, "Computation")
        .value("NONE", Computation::NONE)
        .value("MEAN_POSITION", Computation::MEAN_POSITION)
        .value("MEAN_NORMAL", Computation::MEAN_NORMAL)
        .value("POTENTIAL", Computation::POTENTIAL)
        .value("PRIMITIVE_GRADIENT", Computation::PRIMITIVE_GRADIENT)
        .value("PROJECTION", Computation::PROJECTION)
        .value("PRINCIPAL_CURVATURE", Computation::PRINCIPAL_CURVATURE)
        .value("ALGEBRAIC_SPHERE", Computation::ALGEBRAIC_SPHERE)
        .value("GLS", Computation::GLS)
        .value("GLS_NORMALIZED", Computation::GLS_NORMALIZED)
        .value("QUADRATIC_HEIGHTFIELD_MAP", Computation::QUADRATIC_HEIGHTFIELD_MAP)
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
    using Point      = typename _PyCO::Point;
    using Scalar     = typename _PyCO::Scalar;
    using VectorType = typename _PyCO::VectorType;
    using CO         = typename _PyCO::ComputeObject;

    switch (id)
    {
    case Computation::NONE:
        return std::vector<size_t>{0, 0};
    case Computation::MEAN_POSITION:
        if constexpr (Ponca::ProvidesMeanPosition<CO>)
        {
            if (array.data() != nullptr)
                throw std::invalid_argument("MEAN_POSITION does not take any data as input.");

            return ComputeVectorDim<Point, VectorType>(array);
        }
        break;
    case Computation::MEAN_NORMAL:
        if constexpr (Ponca::ProvidesMeanNormal<CO>)
        {
            if (array.data() != nullptr)
                throw std::invalid_argument("MEAN_NORMAL does not take any data as input.");

            return ComputeVectorDim<Point, VectorType>(array);
        }
        break;
    case Computation::POTENTIAL:
        if constexpr (Ponca::ProvidesImplicitPrimitive<CO>)
        {
            return ComputeVectorDim<Point, Scalar>(array);
        }
        break;
    case Computation::PRIMITIVE_GRADIENT:
        if constexpr (Ponca::ProvidesImplicitPrimitive<CO>)
        {
            return ComputeVectorDim<Point, VectorType>(array);
        }
        break;
    case Computation::PROJECTION:
        if constexpr (Ponca::ProvidesProjectionOperator<CO>)
        {
            if (array.data() == nullptr)
                throw std::invalid_argument("PROJECTION requires an input array");

            return ComputeVectorDim<Point, VectorType>(array);
        }
        break;
    case Computation::PRINCIPAL_CURVATURE:
        if constexpr (Ponca::ProvidesPrincipalCurvatures<CO>)
        {
            if (array.data() != nullptr)
                throw std::invalid_argument("PRINCIPAL_CURVATURE does not take any data as input.");

            // kmin, kmax, gauss, kmindir, kmaxdir
            return ComputeVectorDimGroup<Point, Scalar, Scalar, Scalar, VectorType, VectorType>(array);
        }
        break;
    case Computation::ALGEBRAIC_SPHERE:
        if constexpr (Ponca::ProvidesAlgebraicSphere<CO>)
        {
            if (array.data() != nullptr)
                throw std::invalid_argument("ALGEBRAIC_SPHERE does not take any data as input.");
                
            // Radius, Pratt Norm, center
            return ComputeVectorDimGroup<Point, Scalar, Scalar, VectorType>(array);     
        }
        break;
    case Computation::GLS:
    case Computation::GLS_NORMALIZED: // Same as GLS as far as dimensions are concerned
        if constexpr (Ponca::ProvidesGLSParam<CO>)
        {
            if (array.data() != nullptr)
                throw std::invalid_argument("GLS(_NORMALIZED) does not take any data as input.");
                
            // tau, kappa, eta
            return ComputeVectorDimGroup<Point, Scalar, Scalar, VectorType>(array);            
        }
        break;
    case Computation::QUADRATIC_HEIGHTFIELD_MAP:
        if constexpr(Ponca::ProvidesHeightField<CO>)
        {
            // TODO:
            throw std::runtime_error("QUADRATIC_HEIGHTFIELD_MAP Compute method not implemented yet");

            if (array.data() == nullptr)
                throw std::invalid_argument("QUADRATIC_HEIGHTFIELD_MAP requires input data.");
            
            if (array.shape(array.ndim() - 1) != 2)
                throw std::invalid_argument("");
                
            return ComputeVectorDim<Point, Scalar>(array);
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
PONCA_MULTIARCH inline void ExtractComputation(size_t id, CO& object, size_t i, const DeviceArrayView<Scalar>& in,
                                               DeviceArrayView<Scalar>& out)
{
    using Point = typename CO::DataPoint;

    // if constexprs are just here to allow code to compile. The actual guard is useless and already perfomed
    // in the ComputeOutputDimension.
    Computation computation = static_cast<Computation>(id);
    Scalar* base = out.data + i * out.stride[0];

    switch (computation)
    {
    case Computation::NONE:
        break;
    case Computation::MEAN_POSITION:
        if constexpr (Ponca::ProvidesMeanPosition<CO>)
        {
            Write(base, object.meanPosition().barycenter());
        }
        break;
    case Computation::MEAN_NORMAL:
        if constexpr (Ponca::ProvidesMeanNormal<CO>)
        {
            Write(base, object.meanNormal().meanNormalVector());
        }
        break;
    case Computation::POTENTIAL:
        if constexpr (Ponca::ProvidesImplicitPrimitive<CO>)
        {
            if (in.data != nullptr)
                RunVectorMethod<CO>(object, MAKE_VECTOR_FUNCTION(implicitPrimitive().potential), in, out, i);
            else
                Write(base, object.implicitPrimitive().potential());
        }
        break;
    case Computation::PRIMITIVE_GRADIENT:
        if constexpr (Ponca::ProvidesImplicitPrimitive<CO>)
        {
            if (in.data != nullptr)
                RunVectorMethod<CO>(object, MAKE_VECTOR_FUNCTION(implicitPrimitive().primitiveGradient), in, out, i);
            else
                Write(base, object.implicitPrimitive().primitiveGradient());
        }
        break;
    case Computation::PROJECTION:
        if constexpr (Ponca::ProvidesProjectionOperator<CO>)
        {
            RunVectorMethod<CO>(object, MAKE_VECTOR_FUNCTION(projectionOperator().project), in, out, i);
        }
        break;
    case Computation::PRINCIPAL_CURVATURE:
        if constexpr (Ponca::ProvidesPrincipalCurvatures<CO>)
        {
            Write(base + 0 * out.stride[1], object.curvatureTensor().kmin());
            Write(base + 1 * out.stride[1], object.curvatureTensor().kmax());
            Write(base + 2 * out.stride[1], object.curvatureTensor().GaussianCurvature());
            Write(base + 3 * out.stride[1], object.curvatureTensor().kminDirection());
            Write(base + (3 + Point::Dim) * out.stride[1], object.curvatureTensor().kmaxDirection());
        }
        break;
    case Computation::ALGEBRAIC_SPHERE:
        if constexpr (Ponca::ProvidesAlgebraicSphere<CO>)
        {
            Write(base + 0 * out.stride[1], object.algebraicSphere().radius());
            Write(base + 1 * out.stride[1], object.algebraicSphere().prattNorm());
            Write(base + 2 * out.stride[1], object.algebraicSphere().center());
        }
        break;
    case Computation::GLS:
        if constexpr (Ponca::ProvidesGLSParam<CO>)
        {
            Write(base + 0 * out.stride[1], object.glsParam().tau());
            Write(base + 1 * out.stride[1], object.glsParam().kappa());
            Write(base + 2 * out.stride[1], object.glsParam().eta());
        }
        break;
    case Computation::GLS_NORMALIZED: // Same as GLS as far as dimensions are concerned
        if constexpr (Ponca::ProvidesGLSParam<CO>)
        {
            Write(base + 0 * out.stride[1], object.glsParam().tau_normalized());
            Write(base + 1 * out.stride[1], object.glsParam().kappa_normalized());
            Write(base + 2 * out.stride[1], object.glsParam().eta_normalized());           
        }
        break;
    case Computation::QUADRATIC_HEIGHTFIELD_MAP:
        if constexpr(Ponca::ProvidesHeightField<CO>)
        {
            // Can not use run vector method here 
            // TODO
        }
        break;
    default:
        break;
    }
}
