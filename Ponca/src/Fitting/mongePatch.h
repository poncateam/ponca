/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "./defines.h"
#include "./curvature.h"

#include PONCA_MULTIARCH_INCLUDE_STD(cmath)

#include <Eigen/Dense>

namespace Ponca
{
/*!
    \brief Monge Patch primitive, defined as \f$ \mathbf{x}(u,v)= (u,v,h(u,v)) \f$,
     with \f$h(u,v)\f$ defined by a child class.

    \see MongePatchFit

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

        //! \brief Value of the scalar field at the evaluation point
        //! \see method `#isSigned` of the plane fit to check if the sign is reliable
        PONCA_MULTIARCH [[nodiscard]] inline Scalar potential(const VectorType& _q) const {
            VectorType x = Base::worldToTangentPlane(_q);
            return Base::height(*(x.data()+1),*(x.data()+2)) - *(x.data());
        }

        //! \brief Orthogonal projecting on the patch, such that \f$h = f(u,v)\f$
        PONCA_MULTIARCH [[nodiscard]] inline VectorType project (const VectorType& _q) const
        {
            VectorType x = Base::worldToTangentPlane(_q);
            *(x.data()) = height(*(x.data()+1),*(x.data()+2));
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
            return Base::tangentPlaneToWorld(primitiveGradientLocal(Base::worldToTangentPlane(_q)));
        }

        //! \brief Scalar field gradient direction, both input and output vectors are expressed in the local basis
        PONCA_MULTIARCH [[nodiscard]] inline VectorType primitiveGradientLocal (const VectorType& _localQ = VectorType::Zero()) const
        {
            return Base::heightTangentULocal(_localQ).cross(Base::heightTangentVLocal(_localQ));
        }

        PONCA_MULTIARCH inline void firstFundamentalFormComponents (Scalar &E, Scalar &F, Scalar &G) const
        {
            PONCA_MULTIARCH_STD_MATH(sqrt);
            VectorType fu {Scalar(0), Base::dh_u(), Scalar(0)};
            VectorType fv {Scalar(0), Scalar(0),    Base::dh_v()};
            E = fu.dot(fu);
            F = fu.dot(fv);
            G = fv.dot(fv);
        }

        PONCA_MULTIARCH inline void secondFundamentalFormComponents (Scalar &L, Scalar &M, Scalar &N) const
        {
            PONCA_MULTIARCH_STD_MATH(sqrt);
            VectorType fuu {Scalar(0), Base::dh_uu(), Scalar(0)};
            VectorType fuv {Scalar(0), Base::dh_uv(), Base::dh_uv()};
            VectorType fvv {Scalar(0), Scalar(0),     Base::dh_vv()};
            VectorType fn = primitiveGradientLocal();

            L = fuu.dot(fn);
            M = fuv.dot(fn);;
            N = fvv.dot(fn);;
        }
    }; //class MongePatchPrimitive

    /*!
    \brief Quadratic height field defined as \f$h(u,v)=h_{uu}u^2 + h_{vv}v^2 + h_u u + h_v v + h_c \f$

    \see MongePatchPrimitive

    This primitive provides:
    \verbatim PROVIDES_HEIGHTFIELD, PROVIDES_QUADRIC_HEIGHTFIELD \endverbatim
    */
    template < class DataPoint, class _NFilter, typename T >
    class QuadraticHeightField : public T
    {
    PONCA_FITTING_DECLARE_DEFAULT_TYPES
        using QuadraticHeightFieldCoefficients = Eigen::Matrix<Scalar,6,1>;
        static_assert ( DataPoint::Dim == 3, "QuadraticHeightField is only valid in 3D");

    protected:
        enum {
            PROVIDES_HEIGHTFIELD, /*!< \brief Provides generic heightfield API */
            PROVIDES_QUADRIC_HEIGHTFIELD /*!< \brief Provides quadric heightfield API */
        };
        QuadraticHeightFieldCoefficients  m_coeffs {QuadraticHeightFieldCoefficients::Zero()}; /*!< \brief Quadric parameters */

    public:

        /*! \brief Default constructor */
        PONCA_MULTIARCH inline QuadraticHeightField() : Base() { init(); }

        PONCA_EXPLICIT_CAST_OPERATORS(QuadraticHeightField,quadraticHeightField)

        /// \brief Set the scalar field values
        PONCA_MULTIARCH inline void setQuadric(const QuadraticHeightFieldCoefficients& coeffs) { m_coeffs = coeffs; }

        /// \brief Set the scalar field values to 0
        PONCA_MULTIARCH inline void init()
        {
            Base::init();
            m_coeffs.setZero();
        }

        /// \brief Tell if the plane as been correctly set.
        /// Used to set CONFLICT_ERROR_FOUND during fitting
        /// \return false when called straight after #init. Should be true after fitting
        PONCA_MULTIARCH [[nodiscard]] inline bool isValid() const{
            return ! m_coeffs.isApprox(QuadraticHeightFieldCoefficients::Zero());
        }

        PONCA_MULTIARCH [[nodiscard]] inline bool operator==(const QuadraticHeightField<DataPoint, _NFilter, T>& other) const{
            return m_coeffs.isApprox(other.m_params);
        }

        /*! \brief Comparison operator, convenience function */
        PONCA_MULTIARCH [[nodiscard]] inline bool operator!=(const QuadraticHeightField<DataPoint, _NFilter, T>& other) const{
            return ! ((*this) == other);
        }

        //! \brief Height value at local uv
        PONCA_MULTIARCH [[nodiscard]] inline Scalar height(Scalar u, Scalar v) const {
            return h_uu()*u*u + h_vv()*v*v + h_uv()*u*v + h_u()*u + h_v()*v + h_c();
        }

        PONCA_MULTIARCH [[nodiscard]] inline const Scalar & h_uu () const { return *(m_coeffs.data()); }
        PONCA_MULTIARCH [[nodiscard]] inline const Scalar & h_vv () const { return *(m_coeffs.data()+1); }
        PONCA_MULTIARCH [[nodiscard]] inline const Scalar & h_uv () const { return *(m_coeffs.data()+2); }
        PONCA_MULTIARCH [[nodiscard]] inline const Scalar & h_u  () const { return *(m_coeffs.data()+3); }
        PONCA_MULTIARCH [[nodiscard]] inline const Scalar & h_v  () const { return *(m_coeffs.data()+4); }
        PONCA_MULTIARCH [[nodiscard]] inline const Scalar & h_c  () const { return *(m_coeffs.data()+5); }

        /// \brief Partial derivative \f$ \frac{\delta h }{\delta u} (u,v) = 2h_{uu} u + h_{uv} v + h_u \f$
        PONCA_MULTIARCH [[nodiscard]] inline Scalar dh_u (Scalar u = Scalar(0), Scalar v = Scalar(0)) const
        { return Scalar(2)*h_uu() * u + h_uv()*v + h_u(); }

        /// \brief Partial derivative \f$ \frac{\delta h }{\delta v} (u,v) = 2h_{vv} v + h_{uv} u + h_v \f$
        PONCA_MULTIARCH [[nodiscard]] inline Scalar dh_v (Scalar u = Scalar(0), Scalar v = Scalar(0)) const
        { return Scalar(2)*h_vv() * v + h_uv()*u + h_v(); }

        /// \brief Second order partial derivative \f$ \frac{\delta h }{\delta u^2} (u,v) = 2h_{uu} \f$
        PONCA_MULTIARCH [[nodiscard]] inline Scalar dh_uu (Scalar u = Scalar(0), Scalar v = Scalar(0)) const
        { return Scalar(2)*h_uu(); }

        /// \brief Second order partial derivative \f$ \frac{\delta h }{\delta v} (u,v) = 2h_{vv} \f$
        PONCA_MULTIARCH [[nodiscard]] inline Scalar dh_vv (Scalar u = Scalar(0), Scalar v = Scalar(0)) const
        { return Scalar(2)*h_vv(); }

        /// \brief Second order partial derivative \f$ \frac{\delta }{\delta v} (\frac{\delta h }{\delta u}) (u,v) = h_{uv} \f$
        PONCA_MULTIARCH [[nodiscard]] inline Scalar dh_uv (Scalar u = Scalar(0), Scalar v = Scalar(0)) const
        { return h_uv(); }

        //! \brief Local tangent vector in the direction of \f$u\f$
        PONCA_MULTIARCH [[nodiscard]] inline VectorType heightTangentULocal (const VectorType& _localQ) const {
            PONCA_MULTIARCH_STD_MATH(sqrt);
            const Scalar tu = dh_u(*(_localQ.data()+1), *(_localQ.data()+2));
            return VectorType(tu, Scalar(0), sqrt(Scalar(1) - tu*tu) ); // we build a normalized vector, ie its norm = 1
        }

        //! \brief Local tangent vector in the direction of \f$v\f$
        PONCA_MULTIARCH [[nodiscard]] inline VectorType heightTangentVLocal (const VectorType& _localQ) const {
            PONCA_MULTIARCH_STD_MATH(sqrt);
            const Scalar tv = dh_v(*(_localQ.data()+1), *(_localQ.data()+2));
            return VectorType (Scalar(0), tv, sqrt(Scalar(1) - tv*tv)); // we build a normalized vector, ie its norm = 1
        }
    }; //class QuadraticHeightField

/*!
 * \brief Extension to compute the best fit quadric on 3d points expressed as \f$f(u,v)=h\f$
 *
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
    using SampleMatrix                     = Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>;
    using QuadraticHeightFieldCoefficients = typename Base::QuadraticHeightFieldCoefficients;

protected:
    SampleMatrix                     m_A; /*!< \brief Quadric input samples */
    QuadraticHeightFieldCoefficients m_b {QuadraticHeightFieldCoefficients::Zero()}; /*!< \brief Observations */

    bool m_planeIsReady {false};
public:
    PONCA_EXPLICIT_CAST_OPERATORS(MongePatchQuadraticFitImpl,mongePatchQuadraticFit)
    PONCA_FITTING_DECLARE_INIT_ADD_FINALIZE
}; // MongePatchQuadraticFitImpl

template < class DataPoint, class _NFilter, typename T>
using MongePatchQuadraticFit =
        FundamentalFormCurvatureEstimator<DataPoint, _NFilter,
            MongePatchQuadraticFitImpl<DataPoint, _NFilter,
                MongePatchPrimitive<DataPoint, _NFilter,
                    QuadraticHeightField<DataPoint, _NFilter,
                        CovariancePlaneFit<DataPoint, _NFilter,T>>>>>;

#include "mongePatch.hpp"

} //namespace Ponca
