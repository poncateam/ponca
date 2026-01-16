/*
This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "./defines.h"

#include <Eigen/Dense>

namespace Ponca
{
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
     \brief Quadratic height field defined as \f$h(u,v)=h_{uu}u^2 + h_{vv}v^2 + h_{uv}uv + h_c \f$

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
        using HeightFieldCoefficients = Eigen::Matrix<Scalar,4,1>;
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
            return h_uu()*u*u + h_vv()*v*v + h_uv()*u*v + h_c();
        }

        PONCA_MULTIARCH [[nodiscard]] inline const Scalar & h_uu () const { return *(m_coeffs.data()); }
        PONCA_MULTIARCH [[nodiscard]] inline const Scalar & h_vv () const { return *(m_coeffs.data()+1); }
        PONCA_MULTIARCH [[nodiscard]] inline const Scalar & h_uv () const { return *(m_coeffs.data()+2); }
        PONCA_MULTIARCH [[nodiscard]] inline const Scalar & h_c ()  const { return *(m_coeffs.data()+3); }

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

} //namespace Ponca

