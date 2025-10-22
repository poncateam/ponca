/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#pragma once

#include "./defines.h"

#include PONCA_MULTIARCH_INCLUDE_STD(cmath)
#include PONCA_MULTIARCH_INCLUDE_STD(limits)

#include <Eigen/Core>

namespace Ponca
{

/*!
    \brief Algebraic Sphere primitive

    Method published in \cite Guennebaud:2007:APSS

    An algebraic hyper-sphere is defined as the \f$0\f$-isosurface of the scalar field

    \f$ s_\mathbf{u}(\mathbf{x}) = \left[ 1 \; \mathbf{x}^T \; \mathbf{x}^T\mathbf{x}\right]^T \cdot \mathbf{u} \f$

    with \f$ \mathbf{u} \left[ u_c \; \mathbf{u_l} \; u_q\right]^T \f$ is the
    vector of the constant, linear and quadratic parameters.

    \note If internally the scalar fields are stored in a local frame defined
    by the evaluation position, the public methods involving a query (such as
    project, potential, gradient) have to be defined in global
    coordinates (e.g. you don't need to convert your query in the current locale
    frame).

    This primitive provides:
    \verbatim PROVIDES_ALGEBRAIC_SPHERE \endverbatim

    \todo Deal with planar case (_uq == 0) and what about _ul == 0 ?
*/


template < class DataPoint, class _NFilter, typename T >
class AlgebraicSphere : public T
{
    PONCA_FITTING_DECLARE_DEFAULT_TYPES

protected:
    enum
    {
        check = Base::PROVIDES_PRIMITIVE_BASE,  /*!< \brief Requires PrimitiveBase */
        PROVIDES_ALGEBRAIC_SPHERE               /*!< \brief Provides Algebraic Sphere */
    };

protected:
    //! \brief Is the implicit scalar field normalized using Pratt
    bool m_isNormalized;

// results
public:
    Scalar m_uc {0},       /*!< \brief Constant parameter of the Algebraic hyper-sphere */
           m_uq {0};       /*!< \brief Quadratic parameter of the Algebraic hyper-sphere */
    VectorType m_ul {VectorType::Zero()};   /*!< \brief Linear parameter of the Algebraic hyper-sphere */

public:
    PONCA_EXPLICIT_CAST_OPERATORS(AlgebraicSphere,algebraicSphere)

    /*! \brief Set the scalar field values to 0 and reset the isNormalized() status

        \warning Set m_ul to Zero(), which leads to nans in OrientedSphere::normal()
    */
    PONCA_MULTIARCH inline void init()
    {
        Base::init();

        m_uc = Scalar(0);
        m_ul = VectorType::Zero();
        m_uq = Scalar(0);

        m_isNormalized = false;
    }

    /// \brief Tell if the sphere as been correctly set.
    /// Used to set CONFLICT_ERROR_FOUND during fitting
    /// \return false when called straight after #init. Should be true after fitting
    PONCA_MULTIARCH inline bool isValid() const{
        return !( m_ul.isApprox(VectorType::Zero()) && m_uc == Scalar(0) && m_uq == Scalar(0) );
    }

    /// \brief Comparison operator \warning Assume that other shares the same basis \see changeBasis()
    PONCA_MULTIARCH inline bool operator==(const AlgebraicSphere<DataPoint, NeighborFilter, T>& other) const{
        PONCA_MULTIARCH_STD_MATH(pow);
        const Scalar epsilon        = Eigen::NumTraits<Scalar>::dummy_precision();
        const Scalar squaredEpsilon = epsilon*epsilon;
        return  pow(m_uc - other.m_uc, Scalar(2)) < squaredEpsilon &&
                pow(m_uq - other.m_uq, Scalar(2)) < squaredEpsilon &&
                 m_ul.isApprox(other.m_ul);
    }
    /// \brief Approximate operator \warning Assume that other shares the same basis \see changeBasis()
    template<typename Other>
    PONCA_MULTIARCH inline bool isApprox(const Other& other, const Scalar& epsilon = Eigen::NumTraits<Scalar>::dummy_precision()) const{
        PONCA_MULTIARCH_STD_MATH(pow);
        const Scalar squaredEpsilon = epsilon*epsilon;
        return  pow(m_uc - other.m_uc, Scalar(2)) < squaredEpsilon &&
                pow(m_uq - other.m_uq, Scalar(2)) < squaredEpsilon &&
                 m_ul.isApprox(other.m_ul);
    }


    /*! \brief Comparison operator, convenience function */
    PONCA_MULTIARCH inline bool operator!=(const AlgebraicSphere<DataPoint, NeighborFilter, T>& other) const{
        return ! ((*this) == other);
    }

    /*! \brief Express the scalar field relatively to a new basis 
    
    The sphere in dimension \f$d\f$ is parametrized in the original basis by \f$\mathbf{u} = \left[u_c, \mathbf{u}_l^T, u_q\right]^T \in \mathbb{R}^{d+2}\f$.
    The sphere equation is given in this basis by:

    \begin{equation}
    u_c + \mathbf{u}_l^T.\mathbf{x} + u_q \mathbf{x}^T.\mathbf{x}= 0 
    \end{equation}

    The same equation written in a different basis is:
    \begin{equation}
    u_c + \mathbf{u}_l^T.(\mathbf{x}-\mathbf{\Delta}) + u_q (\mathbf{x}-\mathbf{\Delta})^T.(\mathbf{x}-\mathbf{\Delta})= 0 
    \end{equation}
    where \f$\mathbf{\Delta}\f$ is the difference vector between the two basis: \f$\mathbf{\Delta} = \mathbf{b}_\mathrm{original}-\mathbf{b}_\mathrm{new}\f$.
    It corresponds to a translation of \f$\Delta\f$ of the sphere.
    This develops in:
    \begin{equation}
    \left[u_c - \mathbf{u}_l^T.\mathbf{\Delta} + u_q \mathbf{\Delta}^T.\mathbf{\Delta}\right] + \left[\mathbf{u}_l - 2 u_q \mathbf{\Delta}\right]^T.\mathbf{x} + u_q \mathbf{x}^T.\mathbf{x} = 0 
    \end{equation}

    By identification, the parametrization \f$\left[u_c', \mathbf{u}_l'^T, u_q'\right]^T\f$ of the sphere in the new basis is given by

    \f$u_c' = u_c - \mathbf{u}_l^T.\mathbf{\Delta} + u_q \mathbf{\Delta}^T.\mathbf{\Delta}\f$

    \f$\mathbf{u}_l' = \mathbf{u}_l - 2 u_q \mathbf{\Delta}\f$

    \f$u_q'=u_q\f$

    */
    PONCA_MULTIARCH inline void changeBasis(const VectorType& newbasis)
    {
        VectorType diff = Base::getNeighborFilter().evalPos() - newbasis;
        Base::setNeighborFilter({newbasis, Base::getNeighborFilter().evalScale()});
        Base::init();
        m_uc = m_uc - m_ul.dot(diff) + m_uq * diff.dot(diff);
        m_ul = m_ul - Scalar(2.)*m_uq*diff;
        //m_uq is not changed
        m_isNormalized = false;
        applyPrattNorm();
    }

    /*! \brief compute the Pratt norm of the implicit scalar field. */
    PONCA_MULTIARCH inline Scalar prattNorm() const
    {
        PONCA_MULTIARCH_STD_MATH(sqrt);
        return sqrt(prattNorm2());
    }

    /*! \brief compute the squared Pratt norm of the implicit scalar field. */
    PONCA_MULTIARCH inline Scalar prattNorm2() const
    {
        return m_ul.squaredNorm() - Scalar(4.) * m_uc * m_uq;
    }

    /*!
        \brief Normalize the scalar field by the Pratt norm
        \return false when the normalization fails (sphere is already normalized)
    */
    PONCA_MULTIARCH inline bool applyPrattNorm()
    {
        if (! m_isNormalized)
        {
            Scalar pn = prattNorm();
            m_uc /= pn;
            m_ul *= Scalar(1.)/pn;
            m_uq /= pn;

            m_isNormalized = true;
        }
        return true;
    }

    /*!
        \brief return the estimated radius of the sphere
        \warning return inf if the fitted surface is planar
    */
    PONCA_MULTIARCH inline Scalar radius() const
    {
        PONCA_MULTIARCH_STD_MATH(numeric_limits);
        if(isPlane())
            return numeric_limits<Scalar>::infinity(); // non-sense value

        PONCA_MULTIARCH_STD_MATH(sqrt);
        Scalar b = Scalar(1.)/m_uq;
        return Scalar(sqrt( ((Scalar(-0.5)*b)*m_ul).squaredNorm() - m_uc*b ));
    }

    /*!
        \brief return the estimated center of the sphere
        \warning return Vector inf if the fitted surface is planar
    */
    PONCA_MULTIARCH inline VectorType center() const
    {
        PONCA_MULTIARCH_STD_MATH(numeric_limits);
        if(isPlane())
            return VectorType::Constant(numeric_limits<Scalar>::infinity()); // non-sense value

        Scalar b = Scalar(1.)/m_uq;
        return Base::getNeighborFilter().convertToGlobalBasis((Scalar(-0.5)*b)*m_ul);
    }

    //! \brief State indicating when the sphere has been normalized
    PONCA_MULTIARCH inline bool isNormalized() const { return m_isNormalized; }

    //! \brief Value of the scalar field at the location \f$ \mathbf{q}\f$.
    //! \see method `#isSigned` of the fit to check if the sign is reliable
    PONCA_MULTIARCH inline Scalar potential (const VectorType& _q) const
    {
        // Turn to centered basis
        const VectorType lq = Base::getNeighborFilter().convertToLocalBasis(_q);
        return potentialLocal(lq);
    }

    //! \brief Value of the scalar field at the evaluation point
    //! \see method `#isSigned` of the fit to check if the sign is reliable
    PONCA_MULTIARCH inline Scalar potential() const { return m_uc; }

    /*!
       \brief Project a point on the algebraic hypersphere

        This projection is realized in closed-form: the algebraic hypersphere is converted
        to a geometrical representation (hyperplane or hypersphere), and _q is orthogonally
        projected on the primtive.
        \note This function is in most cases more accurate and faster than #projectDescent
     */
    template <typename NF = NeighborFilter, std::enable_if_t<NF::isLocal, int> = 0> // Enable project only if NF::isLocal
    PONCA_MULTIARCH inline VectorType project (const VectorType& _q) const;

    /*! \brief Approximation of the scalar field gradient at \f$ \mathbf{q}\f$
        \warning The gradient is not normalized by default */
    PONCA_MULTIARCH inline VectorType primitiveGradient (const VectorType& _q) const
    {
        // Turn to centered basis
        const VectorType lq = Base::getNeighborFilter().convertToLocalBasis(_q);
        return primitiveGradientLocal(lq);
    }

    /*! \brief Approximation of the scalar field gradient at the evaluation point
        \warning The gradient is not normalized by default */
    PONCA_MULTIARCH inline const VectorType& primitiveGradient () const { return m_ul; }

    /*!
        \brief Used to know if the fitting result to a plane
        \return true if finalize() have been called and the fitting result to a plane
    */
    PONCA_MULTIARCH inline bool isPlane() const
    {
        PONCA_MULTIARCH_STD_MATH(abs);
        bool bPlanar   = Eigen::internal::isApprox(m_uq,Scalar(0)); 
        bool bReady    = Base::isReady();
        return bReady && bPlanar;
    }

protected:
    /// \copydoc AlgebraicSphere::potential
    PONCA_MULTIARCH inline Scalar potentialLocal (const VectorType& _lq) const;

    /// \copydoc AlgebraicSphere::primitiveGradient
    PONCA_MULTIARCH inline VectorType primitiveGradientLocal (const VectorType& _lq) const;
}; //class AlgebraicSphere

#include "algebraicSphere.hpp"
}
