
/*
 Copyright (C) 2021 aniket agarwalla <aniketagarwalla37@gmail.com>
 
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "./defines.h"
#include <Eigen/Geometry>
#include <Eigen/Core>

namespace Ponca
{

/*!
    \brief A parametrized line is defined by an origin point \f$\mathbf{o}\f$ and a unit direction vector
    \f$\overrightarrow{\mathbf{d}}\f$ such that the line corresponds to the set
    \f$l(t)=\mathbf{o}+t\overrightarrow{\mathbf{d}}, t\in \mathbb{R}\f$.

    This class inherits Eigen::ParametrizedLine.

    This primitive requires the definition of n-dimensionnal vectors
    (VectorType) in Concept::PointConcept.

    This primitive provides:
    \verbatim PROVIDES_LINE \endverbatim
*/

template < class DataPoint, class _WFunctor, typename T >
class Line : public T,
             public Eigen::ParametrizedLine<typename DataPoint::Scalar, DataPoint::Dim >
{
PONCA_FITTING_DECLARE_DEFAULT_TYPES

public:
    /// \brief Specialization of Eigen::ParametrizedLine inherited by Ponca::Line
    using EigenBase = Eigen::ParametrizedLine<typename DataPoint::Scalar, DataPoint::Dim >;

protected:

    enum
    {
        check = Base::PROVIDES_PRIMITIVE_BASE,  /*!< \brief Requires PrimitiveBase */
        PROVIDES_LINE                           /*!< \brief Provides  Line */
    };

public:
    PONCA_EXPLICIT_CAST_OPERATORS(Line,line)

    /*!
     * \brief Set the scalar field values to 0 and reset the distance() and origin() status
    */
    PONCA_MULTIARCH inline void init()
    {
        Base::init();
        EigenBase::origin().setZero();
        EigenBase::direction().setZero();
    }

    /// \brief Tell if the line as been correctly set.
    /// Used to set CONFLICT_ERROR_FOUND during fitting
    /// \return false when called straight after #init. Should be true after fitting
    PONCA_MULTIARCH inline bool isValid() const{
        static const typename EigenBase::VectorType zeros = EigenBase::VectorType::Zero();
        return ! ( EigenBase::origin().isApprox(zeros) && EigenBase::direction().isApprox(zeros) ) ;
    }

    /*! \brief Comparison operator */
    PONCA_MULTIARCH inline bool operator==(const Line<DataPoint, WFunctor, T>& other) const{
        return EigenBase::isApprox(other);
    }

    /*! \brief Comparison operator, convenience function */
    PONCA_MULTIARCH inline bool operator!=(const Line<DataPoint, WFunctor, T>& other) const{
        return ! ((*this) == other);
    }

    /*! \brief Init the line from a direction and a position
       \param direction Orientation of the line, does not need to be normalized
       \param origin Position of the line
    */
    PONCA_MULTIARCH inline void setLine (const VectorType& origin,
                                         const VectorType& direction)
    {
        EigenBase* cc = static_cast<EigenBase*>(this);
        *cc = EigenBase(origin, direction);
    }

    //! \brief Value of the scalar field at the evaluation point
    //! \see method `#isSigned` of the fit to check if the sign is reliable
    PONCA_MULTIARCH inline Scalar potential ( ) const
    {
        // The potential is the distance from a point to the line
        return EigenBase::squaredDistance(VectorType::Zero());
    }

    /*!  \brief Value of the scalar field at the location \f$ \mathbf{q} \f$,
     * defined as the squared distance between \f$ \mathbf{q} \f$ and the line
     *  \see method `#isSigned` of the fit to check if the sign is reliable
     */
    PONCA_MULTIARCH inline Scalar potential (const VectorType& _q) const
    {
        // The potential is the distance from a point to the line
        return EigenBase::squaredDistance(Base::m_w.convertToLocalBasis(_q));
    }

    //! \brief Project a point on the line
    PONCA_MULTIARCH inline VectorType project (const VectorType& _q) const
    {
        // Project on the normal vector and add the offset value
        return Base::m_w.convertToGlobalBasis(EigenBase::projection(Base::m_w.convertToLocalBasis(_q)));
    }
}; //class Line


}
