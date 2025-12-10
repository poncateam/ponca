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
     with \f$h(u,v)\f$ defined by a Base class.

    \see MongePatchFit

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
    \brief Internal base classe for height fields.

     Provides protected convenience methods for heightfields

    This primitive provides:
    \verbatim PROVIDES_HEIGHTFIELD \endverbatim
    */
    template < class DataPoint, class _NFilter, typename T >
    class HeightField : public T
    {
        PONCA_FITTING_DECLARE_DEFAULT_TYPES
        static_assert ( DataPoint::Dim == 3, "HeightField is only valid in 3D");
    protected:
        enum {
            PROVIDES_HEIGHTFIELD /*!< \brief Provides generic heightfield API */
        };

        /// \brief get access to height from local coordinate vector
        PONCA_MULTIARCH [[nodiscard]] inline const Scalar& getHFromLocalCoordinates (const VectorType& _lq) const
        { return *(_lq.data()); }

        /// \brief get access to height from local coordinate vector
        PONCA_MULTIARCH [[nodiscard]] inline Scalar& getHFromLocalCoordinates (VectorType& _lq) const
        { return *(_lq.data()); }

        /// \brief get access to u from local coordinate vector
        PONCA_MULTIARCH [[nodiscard]] inline const Scalar& getUFromLocalCoordinates (const VectorType& _lq) const
        { return *(_lq.data()+1); }

        /// \brief get access to u from local coordinate vector
        PONCA_MULTIARCH [[nodiscard]] inline Scalar& getUFromLocalCoordinates (VectorType& _lq) const
        { return *(_lq.data()+1); }

        /// \brief get access to v from local coordinate vector
        PONCA_MULTIARCH [[nodiscard]] inline const Scalar& getVFromLocalCoordinates (const VectorType& _lq) const
        { return *(_lq.data()+2); }

        /// \brief get access to v from local coordinate vector
        PONCA_MULTIARCH [[nodiscard]] inline Scalar& getVFromLocalCoordinates (VectorType& _lq)
        { return *(_lq.data()+2); }
    };

    /*!
    \brief Quadratic height field defined as \f$h(u,v)=h_{uu}u^2 + h_{vv}v^2 + h_{uv}uv + h_u u + h_v v + h_c \f$

    \see MongePatchPrimitive

    This primitive provides:
    \verbatim PROVIDES_QUADRIC_HEIGHTFIELD \endverbatim
    */
    template < class DataPoint, class _NFilter, typename T >
    class QuadraticHeightField : public T
    {
    PONCA_FITTING_DECLARE_DEFAULT_TYPES
        using HeightFieldCoefficients = Eigen::Matrix<Scalar,6,1>;
        static_assert ( DataPoint::Dim == 3, "QuadraticHeightField is only valid in 3D");

    protected:
        enum {
            Check = Base::PROVIDES_HEIGHTFIELD,
            PROVIDES_QUADRIC_HEIGHTFIELD /*!< \brief Provides quadric heightfield API */
        };
        /// \brief Quadric parameters, stored as \f$[h_uu, h_vv, h_uv, h_u, h_v, h_c]\f$
        HeightFieldCoefficients  m_coeffs {HeightFieldCoefficients::Zero()};

    public:

        /*! \brief Default constructor */
        PONCA_MULTIARCH inline QuadraticHeightField() : Base() { init(); }

        PONCA_EXPLICIT_CAST_OPERATORS(QuadraticHeightField,quadraticHeightField)

        /// \brief Set the scalar field values
        PONCA_MULTIARCH inline void setQuadric(const HeightFieldCoefficients& coeffs) { m_coeffs = coeffs; }

        PONCA_MULTIARCH [[nodiscard]] inline const HeightFieldCoefficients& coeffs() const { return m_coeffs; }

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
            return ! m_coeffs.isApprox(HeightFieldCoefficients::Zero());
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
        PONCA_MULTIARCH [[nodiscard]] inline Scalar dh_du (Scalar u = Scalar(0), Scalar v = Scalar(0)) const
        { return Scalar(2)*h_uu() * u + h_uv()*v + h_u(); }

        /// \brief Partial derivative \f$ \frac{\delta h }{\delta v} (u,v) = 2h_{vv} v + h_{uv} u + h_v \f$
        PONCA_MULTIARCH [[nodiscard]] inline Scalar dh_dv (Scalar u = Scalar(0), Scalar v = Scalar(0)) const
        { return Scalar(2)*h_vv() * v + h_uv()*u + h_v(); }

        /// \brief Second order partial derivative \f$ \frac{\delta h }{\delta u^2} (u,v) = 2h_{uu} \f$
        PONCA_MULTIARCH [[nodiscard]] inline Scalar d2h_duu (Scalar u = Scalar(0), Scalar v = Scalar(0)) const
        { return Scalar(2)*h_uu(); }

        /// \brief Second order partial derivative \f$ \frac{\delta h }{\delta v} (u,v) = 2h_{vv} \f$
        PONCA_MULTIARCH [[nodiscard]] inline Scalar d2h_dvv (Scalar u = Scalar(0), Scalar v = Scalar(0)) const
        { return Scalar(2)*h_vv(); }

        /// \brief Second order partial derivative \f$ \frac{\delta }{\delta v} (\frac{\delta h }{\delta u}) (u,v) = h_{uv} \f$
        PONCA_MULTIARCH [[nodiscard]] inline Scalar d2h_duv (Scalar u = Scalar(0), Scalar v = Scalar(0)) const
        { return h_uv(); }

        //! \brief Local tangent vector in the direction of \f$u\f$
        PONCA_MULTIARCH [[nodiscard]] inline VectorType heightTangentULocal (const VectorType& _localQ) const {
            const Scalar tu = dh_du(Base::getUFromLocalCoordinates(_localQ), Base::getVFromLocalCoordinates(_localQ));
            return VectorType (tu, Scalar(1), Scalar(0)).normalized();
        }

        //! \brief Local tangent vector in the direction of \f$v\f$
        PONCA_MULTIARCH [[nodiscard]] inline VectorType heightTangentVLocal (const VectorType& _localQ) const {
            const Scalar tv = dh_dv(Base::getUFromLocalCoordinates(_localQ), Base::getVFromLocalCoordinates(_localQ));
            return VectorType (tv, Scalar(0), Scalar(1)).normalized();
        }
    }; //class QuadraticHeightField
    /*!
     \brief Quadratic height field defined as \f$h(u,v)=h_{uu}u^2 + h_{vv}v^2 + h_{uv}uv \f$

     \see MongePatchPrimitive, QuadraticHeightField, HeightField

     In contrast to QuadraticHeightField, this version does not hold the linear and constant terms are they are expected
     to be obtained from the support plane

     This primitive provides:
     \verbatim PROVIDES_QUADRIC_HEIGHTFIELD \endverbatim
    */
    template < class DataPoint, class _NFilter, typename T >
    class RestrictedQuadraticHeightField : public T
    {
    PONCA_FITTING_DECLARE_DEFAULT_TYPES
        using HeightFieldCoefficients = Eigen::Matrix<Scalar,3,1>;
        static_assert ( DataPoint::Dim == 3, "QuadraticHeightField is only valid in 3D");

    protected:
        enum {
            Check = Base::PROVIDES_HEIGHTFIELD,
            PROVIDES_RESTRICTED_QUADRIC_HEIGHTFIELD /*!< \brief Provides quadric heightfield API */
        };
        /// \brief Quadric parameters, stored as \f$[h_uu, h_vv, h_uv]\f$
        HeightFieldCoefficients  m_coeffs {HeightFieldCoefficients::Zero()};

    public:

        /*! \brief Default constructor */
        PONCA_MULTIARCH inline RestrictedQuadraticHeightField() : Base() { init(); }

        PONCA_EXPLICIT_CAST_OPERATORS(RestrictedQuadraticHeightField,quadraticHeightField)

        /// \brief Set the scalar field values
        PONCA_MULTIARCH inline void setQuadric(const HeightFieldCoefficients& coeffs) { m_coeffs = coeffs; }

        PONCA_MULTIARCH [[nodiscard]] inline const HeightFieldCoefficients& coeffs() const { return m_coeffs; }

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
            return ! m_coeffs.isApprox(HeightFieldCoefficients::Zero());
        }

        PONCA_MULTIARCH [[nodiscard]] inline bool operator==(const RestrictedQuadraticHeightField<DataPoint, _NFilter, T>& other) const{
            return m_coeffs.isApprox(other.m_params);
        }

        /*! \brief Comparison operator, convenience function */
        PONCA_MULTIARCH [[nodiscard]] inline bool operator!=(const RestrictedQuadraticHeightField<DataPoint, _NFilter, T>& other) const{
            return ! ((*this) == other);
        }

        //! \brief Height value at local uv
        PONCA_MULTIARCH [[nodiscard]] inline Scalar height(Scalar u, Scalar v) const {
            return h_uu()*u*u + h_vv()*v*v + h_uv()*u*v;
        }

        PONCA_MULTIARCH [[nodiscard]] inline const Scalar & h_uu () const { return *(m_coeffs.data()); }
        PONCA_MULTIARCH [[nodiscard]] inline const Scalar & h_vv () const { return *(m_coeffs.data()+1); }
        PONCA_MULTIARCH [[nodiscard]] inline const Scalar & h_uv () const { return *(m_coeffs.data()+2); }

        /// \brief Partial derivative \f$ \frac{\delta h }{\delta u} (u,v) = 2h_{uu} u + h_{uv} v\f$
        PONCA_MULTIARCH [[nodiscard]] inline Scalar dh_du (Scalar u = Scalar(0), Scalar v = Scalar(0)) const
        { return Scalar(2)*h_uu() * u + h_uv()*v; }

        /// \brief Partial derivative \f$ \frac{\delta h }{\delta v} (u,v) = 2h_{vv} v + h_{uv} u\f$
        PONCA_MULTIARCH [[nodiscard]] inline Scalar dh_dv (Scalar u = Scalar(0), Scalar v = Scalar(0)) const
        { return Scalar(2)*h_vv() * v + h_uv()*u ; }

        /// \brief Second order partial derivative \f$ \frac{\delta h }{\delta u^2} (u,v) = 2h_{uu} \f$
        PONCA_MULTIARCH [[nodiscard]] inline Scalar d2h_duu (Scalar u = Scalar(0), Scalar v = Scalar(0)) const
        { return Scalar(2)*h_uu(); }

        /// \brief Second order partial derivative \f$ \frac{\delta h }{\delta v} (u,v) = 2h_{vv} \f$
        PONCA_MULTIARCH [[nodiscard]] inline Scalar d2h_dvv (Scalar u = Scalar(0), Scalar v = Scalar(0)) const
        { return Scalar(2)*h_vv(); }

        /// \brief Second order partial derivative \f$ \frac{\delta }{\delta v} (\frac{\delta h }{\delta u}) (u,v) = h_{uv} \f$
        PONCA_MULTIARCH [[nodiscard]] inline Scalar d2h_duv (Scalar u = Scalar(0), Scalar v = Scalar(0)) const
        { return h_uv(); }

        //! \brief Local tangent vector in the direction of \f$u\f$
        PONCA_MULTIARCH [[nodiscard]] inline VectorType heightTangentULocal (const VectorType& _localQ) const {
            const Scalar tu = dh_du(Base::getUFromLocalCoordinates(_localQ), Base::getVFromLocalCoordinates(_localQ));
            return VectorType (tu, Scalar(1), Scalar(0)).normalized();
        }

        //! \brief Local tangent vector in the direction of \f$v\f$
        PONCA_MULTIARCH [[nodiscard]] inline VectorType heightTangentVLocal (const VectorType& _localQ) const {
            const Scalar tv = dh_dv(Base::getUFromLocalCoordinates(_localQ), Base::getVFromLocalCoordinates(_localQ));
            return VectorType (tv, Scalar(0), Scalar(1)).normalized();
        }
    }; //class RestrictedQuadraticHeightField

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
            CurvatureEstimatorBase<DataPoint, _NFilter,
                FundamentalFormWeingartenEstimator<DataPoint, _NFilter,
                    MongePatchQuadraticFitImpl<DataPoint, _NFilter,
                        MongePatchPrimitive<DataPoint, _NFilter,
                            QuadraticHeightField<DataPoint, _NFilter,
                                HeightField<DataPoint, _NFilter,
                                    CovariancePlaneFit<DataPoint, _NFilter,T>>>>>>>>;

template < class DataPoint, class _NFilter, typename T>
using MongePatchRestrictedQuadraticFit =
        WeingartenCurvatureEstimator<DataPoint, _NFilter,
            CurvatureEstimatorBase<DataPoint, _NFilter,
                FundamentalFormWeingartenEstimator<DataPoint, _NFilter,
                    MongePatchRestrictedQuadraticFitImpl<DataPoint, _NFilter,
                        MongePatchPrimitive<DataPoint, _NFilter,
                            RestrictedQuadraticHeightField<DataPoint, _NFilter,
                                HeightField<DataPoint, _NFilter,
                                    CovariancePlaneFit<DataPoint, _NFilter,T>>>>>>>>;

#include "mongePatch.hpp"

} //namespace Ponca
