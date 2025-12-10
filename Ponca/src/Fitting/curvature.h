/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "./defines.h"

#include <Eigen/Eigenvalues>
#include <Eigen/Core>

namespace Ponca
{
    template < class DataPoint, class _NFilter, typename T>
/**
 *
 * \brief Base class for any 3d curvature estimator: holds \f$k_{\min}\f$, \f$k_{\max}\f$ and associated vectors,
 * such that \f$ k_{\min} <= k_{\max} \f$
 */
    class CurvatureEstimatorBase : public T
    {
    PONCA_FITTING_DECLARE_DEFAULT_TYPES
    PONCA_FITTING_DECLARE_MATRIX_TYPE

    protected:
        enum
        {
            PROVIDES_PRINCIPAL_CURVATURES
        };

    private:
        /// \brief Minimal principal curvature
        Scalar m_kmin {0},
        /// \brief Maximal principal curvature
        m_kmax {0};
        /// \brief Direction associated to the minimal principal curvature
        VectorType m_vmin {VectorType::Zero()},
        /// \brief Direction associated to the maximal principal curvature
        m_vmax {VectorType::Zero()};

        /// \brief Internal state indicating if the curvature are set to default (false), or have been computed (true)
        bool m_isValid {false};

        static_assert ( DataPoint::Dim == 3, "CurvatureEstimatorBase is only valid in 3D");

    public:
        PONCA_EXPLICIT_CAST_OPERATORS(CurvatureEstimatorBase, curvatureEstimatorBase)
        PONCA_FITTING_DECLARE_INIT

        /// \brief Returns true if contains valid curvature values (and not default ones)
        PONCA_MULTIARCH [[nodiscard]] inline bool isValid() const { return m_isValid; }

        /**************************************************************************/
        /* Use results                                                            */
        /**************************************************************************/
        //! \brief Returns an estimate of the minimal principal curvature value
        PONCA_MULTIARCH [[nodiscard]] inline Scalar kmin() const { return m_kmin; }

        //! \brief Returns an estimate of the maximal principal curvature value
        PONCA_MULTIARCH [[nodiscard]] inline Scalar kmax() const { return m_kmax; }

        //! \brief Returns an estimate of the minimal principal curvature direction
        PONCA_MULTIARCH [[nodiscard]] inline VectorType kminDirection() const { return m_vmin; }

        //! \brief Returns an estimate of the maximal principal curvature direction
        PONCA_MULTIARCH [[nodiscard]] inline VectorType kmaxDirection() const { return m_vmax; }

        //! \brief Returns an estimate of the mean curvature
        PONCA_MULTIARCH [[nodiscard]] inline Scalar kMean() const { return (m_kmin + m_kmax)/Scalar(2);}

        //! \brief Returns an estimate of the Gaussian curvature
        PONCA_MULTIARCH [[nodiscard]] inline Scalar GaussianCurvature() const { return m_kmin * m_kmax;}

    protected:
        /// \brief Set curvature values. To be called in finalize() by child classes
        ///
        /// If the given parameters are such that \f$ k_{\min} > k_{\max} \f$, then this
        /// method swap the two curvature values and directions and store them such that
        /// \f$ k_{\min} <= k_{\max} \f$.
        ///
        PONCA_MULTIARCH inline void setCurvatureValues(Scalar kmin, Scalar kmax, const VectorType& vmin, const VectorType& vmax);
    };

    template < class DataPoint, class WFunctor, int /*DiffType*/, typename T>
    using CurvatureEstimatorBaseDer = CurvatureEstimatorBase<DataPoint, WFunctor, T>;



/*!
    \brief Compute principal curvatures from a base class providing fundamental forms

    \see MongePatchPrimitive


    This primitive provides:
    \verbatim PROVIDES_WEINGARTEN_MAP \endverbatim

    This primitive requires:
    \verbatim PROVIDES_TANGENT_PLANE_BASIS, PROVIDES_FIRST_FUNDAMENTAL_FORM_COMPONENTS, PROVIDES_SECOND_FUNDAMENTAL_FORM_COMPONENTS \endverbatim
    */
    template < class DataPoint, class _NFilter, typename T >
    class FundamentalFormWeingartenEstimator : public T
    {
    PONCA_FITTING_DECLARE_DEFAULT_TYPES
        using Matrix2 = Eigen::Matrix<Scalar, 2, 2>;
        static_assert ( DataPoint::Dim == 3, "FundamentalFormCurvatureEstimator is only valid in 3D");

    protected:
        enum {
            Check = Base::PROVIDES_FIRST_FUNDAMENTAL_FORM_COMPONENTS &&
                    Base::PROVIDES_SECOND_FUNDAMENTAL_FORM_COMPONENTS,
            PROVIDES_WEINGARTEN_MAP
        };
    public:
        PONCA_EXPLICIT_CAST_OPERATORS(FundamentalFormWeingartenEstimator, fundamentalFormWeingartenEstimator)

        /// \brief Construct and return the first fundamental form from the base class
        ///
        /// \return first fundamental form
        /// \see firstFundamentalForm(Matrix2& first) for an in-place version
        PONCA_MULTIARCH [[nodiscard]] inline Matrix2 firstFundamentalForm() const;

        /// \copybrief firstFundamentalForm()
        /// \tparam Matrix2Derived Input matrix type that must have same interface than Matrix2
        template<typename Matrix2Derived>
        PONCA_MULTIARCH inline void firstFundamentalForm(Matrix2Derived& first) const;

        /// \brief Construct and return the second fundamental form from the base class
        ///
        /// \return second fundamental form
        /// \see secondFundamentalForm(Matrix2& second) for an in-place version
        PONCA_MULTIARCH [[nodiscard]] inline Matrix2 secondFundamentalForm() const;

        /// \copybrief secondtFundamentalForm()
        /// \tparam Matrix2Derived Input matrix type that must have same interface than Matrix2
        template<typename Matrix2Derived>
        PONCA_MULTIARCH inline void secondFundamentalForm(Matrix2Derived &second) const;

        //! \brief Returns the Weingarten Map, defined as
        /// \FIXME fix documentation
        ///
        /// \return Weingarten Map
        /// \see weingartenMap(Matrix2& w) for an in-place version
        PONCA_MULTIARCH [[nodiscard]] inline Matrix2 weingartenMap() const;

        //! \copybrief weingartenMap()
        /// \tparam Matrix2Derived Input matrix type that must have same interface than Matrix2
        template<typename Matrix2Derived>
        PONCA_MULTIARCH inline void weingartenMap(Matrix2Derived &w) const;

        /// \brief Returns an estimate of the mean curvature directly from the fundamental forms
        PONCA_MULTIARCH inline Scalar kMean() const;

        //! \brief Returns an estimate of the Gaussian curvature directly from the fundamental forms
        PONCA_MULTIARCH inline Scalar GaussianCurvature() const;
    };


/*!
    \brief Compute principal curvatures from a base class providing fundamental forms

    This class extracts curvature information from the spatial derivatives of the normal field \f$ N \f$.
    It first assemble a 2x2 matrix representation of the shape operator


    This primitive provides:
    \verbatim PROVIDES_WEINGARTEN_MAP \endverbatim

    This primitive requires:
    \verbatim PROVIDES_TANGENT_PLANE_BASIS, PROVIDES_FIRST_FUNDAMENTAL_FORM_COMPONENTS, PROVIDES_SECOND_FUNDAMENTAL_FORM_COMPONENTS \endverbatim
    */
    template < class DataPoint, class _NFilter, int DiffType, typename T >
    class NormalDerivativeWeingartenEstimator : public T
    {
        PONCA_FITTING_DECLARE_DEFAULT_TYPES
        PONCA_FITTING_DECLARE_MATRIX_TYPE
        PONCA_FITTING_DECLARE_DEFAULT_DER_TYPES
        using Matrix2 = Eigen::Matrix<Scalar, 2, 2>;
        static_assert ( DataPoint::Dim == 3, "FundamentalFormCurvatureEstimator is only valid in 3D");

    protected:
        enum {
            Check = Base::PROVIDES_NORMAL_DERIVATIVE,
            PROVIDES_WEINGARTEN_MAP,
            PROVIDES_TANGENT_PLANE_BASIS
        };

    private:
        MatrixType m_tangentBasis {MatrixType::Zero()};

    public:
        PONCA_EXPLICIT_CAST_OPERATORS_DER(NormalDerivativeWeingartenEstimator, normalDerivativeWeingartenEstimator)
        PONCA_FITTING_DECLARE_FINALIZE

        //! \brief Returns the Weingarten Map, defined as
        /// \FIXME fix documentation
        ///
        /// \return Weingarten Map
        /// \see weingartenMap(Matrix2& w) for an in-place version
        PONCA_MULTIARCH [[nodiscard]] inline Matrix2 weingartenMap() const;

        //! \copybrief weingartenMap()
        /// \tparam Matrix2Derived Input matrix type that must have same interface than Matrix2
        template<typename Matrix2Derived>
        PONCA_MULTIARCH inline void weingartenMap(Matrix2Derived &w) const;


        /// \copydoc CovariancePlaneFitImpl::worldToTangentPlane
        PONCA_MULTIARCH inline VectorType worldToTangentPlane(const VectorType &_q,
                                                              bool _isPositionVector = true) const;

        /// \copydoc CovariancePlaneFitImpl::tangentPlaneToWorld
        PONCA_MULTIARCH inline VectorType tangentPlaneToWorld(const VectorType &_q,
                                                              bool _isPositionVector = true) const;
    };

/*!
    \brief Compute principal curvatures from a base class providing fundamental forms



    \warning To map from the tangent plane to ambiant space, we assume that the patch is represented as
    (h(u,v), u, v) because CovariancePlaneFitImpl::worldToTangentPlane() wraps up coordinates in this order (height, u and v).

    This primitive provides:
    \verbatim PROVIDES_PRINCIPAL_CURVATURES \endverbatim

    This primitive requires:
    \verbatim PROVIDES_TANGENT_PLANE_BASIS, PROVIDES_WEINGARTEN_MAP \endverbatim
    */
    template < class DataPoint, class _NFilter, typename T >
    class WeingartenCurvatureEstimatorBase : public T
    {
    PONCA_FITTING_DECLARE_DEFAULT_TYPES
        using Matrix2 = Eigen::Matrix<Scalar, 2, 2>;
        static_assert ( DataPoint::Dim == 3, "WeingartenCurvatureEstimator is only valid in 3D");

    protected:
        enum {
            Check = Base::PROVIDES_TANGENT_PLANE_BASIS &&  // required for tangentPlaneToWorld
                    Base::PROVIDES_PRINCIPAL_CURVATURES && // required curvature storage
                    Base::PROVIDES_WEINGARTEN_MAP,
            PROVIDES_PRINCIPAL_CURVATURES  /*!< \brief Provides curvature API */
        };

    public:
        PONCA_FITTING_DECLARE_FINALIZE
    };

    template < class DataPoint, class _NFilter, typename T >
    struct WeingartenCurvatureEstimator : public WeingartenCurvatureEstimatorBase<DataPoint, _NFilter, T>{
        PONCA_FITTING_DECLARE_DEFAULT_TYPES
        PONCA_EXPLICIT_CAST_OPERATORS(WeingartenCurvatureEstimator, weingartenCurvatureEstimator)
    };

    template < class DataPoint, class _NFilter, int DiffType, typename T >
    struct WeingartenCurvatureEstimatorDer : public WeingartenCurvatureEstimatorBase<DataPoint, _NFilter, T>{
        PONCA_FITTING_DECLARE_DEFAULT_TYPES
        PONCA_FITTING_DECLARE_DEFAULT_DER_TYPES
        PONCA_EXPLICIT_CAST_OPERATORS_DER(WeingartenCurvatureEstimatorDer, weingartenCurvatureEstimator)

    };

} //namespace Ponca

#include "curvature.hpp"
