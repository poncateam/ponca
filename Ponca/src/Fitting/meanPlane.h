#pragma once

#include "./defines.h"

#include PONCA_MULTIARCH_INCLUDE_STD(cmath)
#include PONCA_MULTIARCH_INCLUDE_STD(limits)

namespace Ponca
{

template < class DataPoint, class _WFunctor, typename T >
class MeanPlane : public T
{
    PONCA_FITTING_DECLARE_DEFAULT_TYPES
    PONCA_FITTING_DECLARE_MATRIX_TYPE

protected:
    enum
    {
        check = Base::PROVIDES_MEAN_POSITION && Base::PROVIDES_MEAN_NORMAL && Base::PROVIDES_PLANE,  /*!< \brief Requires PrimitiveBase */
        PROVIDES_MEAN_PLANE_BASIS                                       /*!< \brief Provides Algebraic Sphere */
    };

protected:
    MatrixType B; // local frame

public:
    PONCA_EXPLICIT_CAST_OPERATORS(MeanPlane,meanPlane)
    /*! \brief Set the scalar field values to 0 and reset the isNormalized() status

        \warning Set m_ul to Zero(), which leads to nans in OrientedSphere::normal()
    */
    PONCA_MULTIARCH inline void init(const VectorType& _basisCenter = VectorType::Zero())
    {
        Base::init(_basisCenter);
        B = MatrixType::Identity();
    }

    /*!
     * \brief Express a point in ambient space relatively to the tangent plane.
     *
     * Output vector is: [h, u, v]^T, where u, v are 2d coordinates on the plane,
     * and h the height of the sample.
     * \tparam ignoreTranslation must be set to true when passing vectors instead of points
     * \param _q Point coordinates expressed in ambient space
     * \return Point coordinates expressed in local tangent frame
     */
    template <bool ignoreTranslation = false>
    PONCA_MULTIARCH inline VectorType worldToTangentPlane(const VectorType &_q) const;

        /*!
     * \brief Transform a point from the tangent plane [h, u, v]^T to ambient space
     *
     * \tparam ignoreTranslation must be set to true when passing vectors instead of points
     * \param _q Point coordinates expressed in local tangent frame
     * \return Point coordinates expressed in ambient space
     */
    template <bool ignoreTranslation = false>
    PONCA_MULTIARCH inline VectorType tangentPlaneToWorld(const VectorType &_q) const;

}; //class MeanPlane

#include "meanPlane.hpp"

} // namespace Ponca