/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "./defines.h"
#include "./curvature.h"
#include "./weingarten.h"
#include "./heightField.h"

#include PONCA_MULTIARCH_INCLUDE_STD(cmath)

#include <Eigen/Dense>

namespace Ponca
{
/*!
    \brief Monge Patch primitive, defined as \f$ \mathbf{x}(u,v)= (u,v,h(u,v)) \f$,
     with \f$h(u,v)\f$ defined by a Base class.

    \see MongePatchQuadraticFitImpl and MongePatchRestrictedQuadraticFitImpl

    \warning Internally, the patch is represented as (h(u,v), u, v) because
    CovariancePlaneFitImpl::worldToTangentPlane() wraps up coordinates in this order (height, u and v).

    This primitive provides:
    \verbatim PROVIDES_MONGE_PATCH, PROVIDES_FIRST_FUNDAMENTAL_FORM_COMPONENTS, PROVIDES_SECOND_FUNDAMENTAL_FORM_COMPONENTS \endverbatim

    This primitive requires:
    \verbatim PROVIDES_PLANE, PROVIDES_TANGENT_PLANE_BASIS, PROVIDES_HEIGHTFIELD \endverbatim
    */
    template < class DataPoint, class _NFilter, typename T >
    class MongePatchPrimitive : public T
    {
        PONCA_FITTING_DECLARE_DEFAULT_TYPES
        static_assert ( DataPoint::Dim == 3, "MongePatch is only valid in 3D");

    protected:
        enum {
            Check = Base::PROVIDES_PLANE &&
                    Base::PROVIDES_TANGENT_PLANE_BASIS &&
                    Base::PROVIDES_HEIGHTFIELD,      /*!< \brief Requires a heightfield function */
                    PROVIDES_MONGE_PATCH,            /*!< \brief Provides MongePatch API */
                    PROVIDES_FIRST_FUNDAMENTAL_FORM_COMPONENTS, /*!< \brief Provides first fundamental form */
                    PROVIDES_SECOND_FUNDAMENTAL_FORM_COMPONENTS  /*!< \brief Provides second fundamental form */
        };

    public:

        /*! \brief Default constructor */
        PONCA_MULTIARCH inline MongePatchPrimitive() : Base() { Base::init(); }

        PONCA_EXPLICIT_CAST_OPERATORS(MongePatchPrimitive,mongePatchPrimitive)

        //! \brief Value of the scalar field at the evaluation point
        //! \see method `#isSigned` of the fit to check if the sign is reliable
        PONCA_MULTIARCH [[nodiscard]] inline Scalar potential ( ) const
        {
            return Base::height(Scalar(0),Scalar(0));
        }

        //! \brief Value of the scalar field at a given point
        //! \see method `#isSigned` of the plane fit to check if the sign is reliable
        PONCA_MULTIARCH [[nodiscard]] inline Scalar potential(const VectorType& _q) const {
            const VectorType x = Base::worldToTangentPlane(_q);
            return Base::height(Base::getUFromLocalCoordinates(x),Base::getVFromLocalCoordinates(x)) - Base::getHFromLocalCoordinates(x);
        }

        //! \brief Orthogonal projection on the patch
        ///
        /// Given a point p and its local representation q, project the point such that \f$h(q.x(),q.y())-q.z()=0\f$
        PONCA_MULTIARCH [[nodiscard]] inline VectorType project (const VectorType& _q) const
        {
            VectorType x = Base::worldToTangentPlane(_q);
            Base::getHFromLocalCoordinates(x) = Base::height(Base::getUFromLocalCoordinates(x),Base::getVFromLocalCoordinates(x));
            return Base::tangentPlaneToWorld(x);
        }

        //! \brief Scalar field gradient direction at the basis center
        PONCA_MULTIARCH [[nodiscard]] inline VectorType primitiveGradient () const
        {
            return Base::tangentPlaneToWorld(primitiveGradientLocal());
        }

        //! \brief Scalar field gradient direction at \f$ \mathbf{q}\f$
        PONCA_MULTIARCH [[nodiscard]] inline VectorType primitiveGradient (const VectorType& _q) const
        {
            return Base::tangentPlaneToWorld(primitiveGradientLocal(Base::worldToTangentPlane(_q)), false);
        }

        //! \brief Scalar field gradient direction, both input and output vectors are expressed in the local basis
        PONCA_MULTIARCH [[nodiscard]] inline VectorType primitiveGradientLocal (const VectorType& _localQ = VectorType::Zero()) const
        {
            VectorType  uVect = Base::heightTangentULocal(_localQ);
            VectorType  vVect = Base::heightTangentVLocal(_localQ);
            VectorType  cross = uVect.cross(vVect);
            VectorType  crossN = uVect.normalized().cross(vVect.normalized());
            return Base::heightTangentULocal(_localQ).cross(Base::heightTangentVLocal(_localQ));
        }

        PONCA_MULTIARCH inline void firstFundamentalFormComponents (Scalar &E, Scalar &F, Scalar &G) const
        {
            PONCA_MULTIARCH_STD_MATH(sqrt);
            const VectorType fu {Base::dh_du(), Scalar(1), Scalar(0)};
            const VectorType fv {Base::dh_dv(), Scalar(0), Scalar(1)};
            E = fu.dot(fu);
            F = fu.dot(fv);
            G = fv.dot(fv);
        }

        PONCA_MULTIARCH inline void secondFundamentalFormComponents (Scalar &L, Scalar &M, Scalar &N) const
        {
            PONCA_MULTIARCH_STD_MATH(sqrt);
            const VectorType fuu {Base::d2h_duu(), Scalar(0), Scalar(0)};
            const VectorType fuv {Base::d2h_duv(), Scalar(0), Scalar(0)};
            const VectorType fvv {Base::d2h_dvv(), Scalar(0), Scalar(0)};
            const VectorType fn = primitiveGradientLocal();

            L = fuu.dot(fn);
            M = fuv.dot(fn);;
            N = fvv.dot(fn);;
        }
    }; //class MongePatchPrimitive

/*!
 * \brief Extension to compute the best fit quadric on 3d points expressed as \f$f(u,v)=h\f$
 *
 * \see MongePatchPrimitive
 * \note This procedure requires at least two passes, the first one for plane fitting,
 * the second one for quadric fitting.
 * \warning This class is valid only in 3D.
 */
    template < class DataPoint, class _NFilter, typename T>
    class MongePatchQuadraticFitImpl : public T
    {
    PONCA_FITTING_DECLARE_DEFAULT_TYPES

    protected:
        enum {
            Check = Base::PROVIDES_QUADRIC_HEIGHTFIELD &&
                    Base::PROVIDES_MONGE_PATCH
        };

    public:
        // we need to use dynamic matric to use ThinU, ThinV for solving
        using SampleMatrix                     = Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>;
        using QuadraticHeightFieldCoefficients = typename Base::HeightFieldCoefficients;

//protected:
    public: //temporary DO NOT MERGE THIS
        SampleMatrix                     m_A; /*!< \brief Quadric input samples */
        QuadraticHeightFieldCoefficients m_b {QuadraticHeightFieldCoefficients::Zero()}; /*!< \brief Observations */

        bool m_planeIsReady {false};
    public:
        PONCA_EXPLICIT_CAST_OPERATORS(MongePatchQuadraticFitImpl,mongePatchQuadraticFit)
        PONCA_FITTING_DECLARE_INIT_ADD_FINALIZE
    }; // MongePatchQuadraticFitImpl
/*!
 * \brief Extension to compute the best fit restricted quadric on 3d points expressed as \f$f(u,v)=h\f$
 *
 * \see MongePatchPrimitive
 * \note This procedure requires at least two passes, the first one for plane fitting,
 * the second one for quadric fitting.
 * \warning This class is valid only in 3D.
 */
    template < class DataPoint, class _NFilter, typename T>
    class MongePatchRestrictedQuadraticFitImpl : public T
    {
    PONCA_FITTING_DECLARE_DEFAULT_TYPES

    protected:
        enum {
            Check = Base::PROVIDES_RESTRICTED_QUADRIC_HEIGHTFIELD &&
                    Base::PROVIDES_MONGE_PATCH
        };

    public:
        // we need to use dynamic matric to use ThinU, ThinV for solving
        using SampleMatrix                     = Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>;
        using QuadraticHeightFieldCoefficients = typename Base::HeightFieldCoefficients;

protected:
        SampleMatrix                     m_A; /*!< \brief Quadric input samples */
        QuadraticHeightFieldCoefficients m_b {QuadraticHeightFieldCoefficients::Zero()}; /*!< \brief Observations */

        bool m_planeIsReady {false};
    public:
        PONCA_EXPLICIT_CAST_OPERATORS(MongePatchRestrictedQuadraticFitImpl,mongePatchQuadraticFit)
        PONCA_FITTING_DECLARE_INIT_ADD_FINALIZE
    }; // MongePatchRestrictedQuadraticFitImpl

template < class DataPoint, class _NFilter, typename T>
using MongePatchQuadraticFit =
        WeingartenCurvatureEstimator<DataPoint, _NFilter,
            CurvatureEstimator<DataPoint, _NFilter,
                FundamentalFormWeingartenEstimator<DataPoint, _NFilter,
                    MongePatchQuadraticFitImpl<DataPoint, _NFilter,
                        MongePatchPrimitive<DataPoint, _NFilter,
                            QuadraticHeightField<DataPoint, _NFilter,
                                HeightField<DataPoint, _NFilter,
                                    CovariancePlaneFit<DataPoint, _NFilter,T>>>>>>>>;

template < class DataPoint, class _NFilter, typename T>
using MongePatchRestrictedQuadraticFit =
        WeingartenCurvatureEstimator<DataPoint, _NFilter,
            CurvatureEstimator<DataPoint, _NFilter,
                FundamentalFormWeingartenEstimator<DataPoint, _NFilter,
                    MongePatchRestrictedQuadraticFitImpl<DataPoint, _NFilter,
                        MongePatchPrimitive<DataPoint, _NFilter,
                            RestrictedQuadraticHeightField<DataPoint, _NFilter,
                                HeightField<DataPoint, _NFilter,
                                    CovariancePlaneFit<DataPoint, _NFilter,T>>>>>>>>;

#include "mongePatch.hpp"

} //namespace Ponca
