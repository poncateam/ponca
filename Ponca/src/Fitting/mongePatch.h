/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "./defines.h"

#include <Eigen/Dense>

#include "curvature.h"

namespace Ponca
{

/*!
 * \brief Extension to compute the best fit quadric on 3d points expressed as \f$f(u,v)=h\f$
 *
 * \note This procedure requires at least two passes, the first one for plane fitting,
 * the second one for quadric fitting.
 * \warning This class is valid only in 3D.
 *
 * \note This class mixes the primitive (MongePatch) and its fitting procedure.
 *       Could makes sense to split the two
 */
template < class DataPoint, class _NFilter, typename T>
class MongePatch : public CurvatureEstimatorBase<DataPoint, _NFilter, T>
{
PONCA_FITTING_DECLARE_DEFAULT_TYPES

protected:
    enum { Check = Base::PROVIDES_PLANE && Base::PROVIDES_TANGENT_PLANE_BASIS };

public:
    using SampleMatrix = Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>;
    using Vector6      = Eigen::Matrix<Scalar,6,1>;

protected:
    SampleMatrix m_A; /*!< \brief Quadric input samples */
    Vector6      m_x {Vector6::Zero()};      /*!< \brief Quadric parameters */
    Vector6      m_b {Vector6::Zero()};      /*!< \brief Observations */

    bool m_planeIsReady {false};
public:
    PONCA_EXPLICIT_CAST_OPERATORS(MongePatch,mongePatch)
    PONCA_FITTING_DECLARE_INIT_ADD_FINALIZE

    //! \brief Returns an estimate of the mean curvature
    PONCA_MULTIARCH [[nodiscard]] inline Scalar kMean() const;

    //! \brief Returns an estimate of the Gaussian curvature
    PONCA_MULTIARCH [[nodiscard]] inline Scalar GaussianCurvature() const;

    PONCA_MULTIARCH [[nodiscard]] inline Scalar evalUV(Scalar u, Scalar v) const {
      return h_uu()*u*u + h_vv()*v*v + h_uv()*u*v + h_u()*u + h_v()*v + h_c();
    }

    //! \brief Value of the scalar field at the evaluation point
    //! \see method `#isSigned` of the plane fit to check if the sign is reliable
    PONCA_MULTIARCH [[nodiscard]] inline Scalar potential(const VectorType& _q, const bool /*convertToLocalBasis*/ = true) const {
      VectorType x = Base::worldToTangentPlane(_q);
      return evalUV(*(x.data()+1),*(x.data()+2)) - *(x.data());
    }

    //! \brief Orthogonal projecting on the patch, such that h = f(u,v)
    PONCA_MULTIARCH [[nodiscard]] inline VectorType project (const VectorType& _q) const
    {
        VectorType x = Base::worldToTangentPlane(_q);
        *(x.data()) = evalUV(*(x.data()+1),*(x.data()+2));
        return Base::tangentPlaneToWorld(x);
    }

    PONCA_MULTIARCH [[nodiscard]] inline const Scalar & h_uu () const { return *(m_x.data()); }
    PONCA_MULTIARCH [[nodiscard]] inline const Scalar & h_vv () const { return *(m_x.data()+1); }
    PONCA_MULTIARCH [[nodiscard]] inline const Scalar & h_uv () const { return *(m_x.data()+2); }
    PONCA_MULTIARCH [[nodiscard]] inline const Scalar & h_u  () const { return *(m_x.data()+3); }
    PONCA_MULTIARCH [[nodiscard]] inline const Scalar & h_v  () const { return *(m_x.data()+4); }
    PONCA_MULTIARCH [[nodiscard]] inline const Scalar & h_c  () const { return *(m_x.data()+5); }

};

#include "mongePatch.hpp"

} //namespace Ponca
