
/*
 Copyright (C) 2021 aniket agarwalla <aniketagarwalla37@gmail.com>
 
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "./defines.h"
#include "./primitive.h" // PrimitiveBase
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

    \ingroup fitting
*/

template < class DataPoint, class _WFunctor, typename T = void  >
class Line : public PrimitiveBase<DataPoint, _WFunctor>,
             public Eigen::ParametrizedLine<typename DataPoint::Scalar, DataPoint::Dim >
{
private:

   typedef PrimitiveBase<DataPoint, _WFunctor> Base;

public:
    /// \brief Specialization of Eigen::ParametrizedLine inherited by Ponca::Line
    using EigenBase = Eigen::ParametrizedLine<typename DataPoint::Scalar, DataPoint::Dim >;

protected:

    enum
    {
        PROVIDES_LINE /*!< \brief Provides  Line */
    };

public:

    /*! \brief Scalar type inherited from DataPoint */
    typedef typename DataPoint::Scalar      Scalar;
    /*! \brief Vector type inherited from DataPoint */
    typedef typename DataPoint::VectorType  VectorType;
    /*! \brief Matrix type inherited from DataPoint */
    typedef typename DataPoint::MatrixType  MatrixType;
    /*! \brief Weight Function */
    typedef _WFunctor                       WFunctor;

private:

    /*! \brief Evaluation position (needed for centered basis) */
    VectorType m_p;

public:

    /*! \brief Default constructor */
    PONCA_MULTIARCH inline Line()
        : Base()
    {
        m_p = VectorType::Zero();
        resetPrimitive();
    }


    /*! \brief Explicit conversion to Line, to access methods potentially hidden by inheritage */
    PONCA_MULTIARCH inline
    Line<DataPoint, WFunctor, T>& line()
    { return * static_cast<Line<DataPoint, WFunctor, T>*>(this); }

    /*!
     * \brief Set the scalar field values to 0 and reset the distance() and origin() status
    */
    PONCA_MULTIARCH inline void resetPrimitive()
    {
        Base::resetPrimitive();
        EigenBase* cc = static_cast<EigenBase*>(this);
        *cc = EigenBase();
    }
    /*! \brief Comparison operator */
    PONCA_MULTIARCH inline bool operator==(const Line<DataPoint, WFunctor, T>& other) const{
        return EigenBase::isApprox(other);
    }

    /*! \brief Comparison operator, convenience function */
    PONCA_MULTIARCH inline bool operator!=(const Line<DataPoint, WFunctor, T>& other) const{
        return ! ((*this) == other);
    }
      /*! \brief Reading access to the basis center (evaluation position) */
    PONCA_MULTIARCH inline const VectorType& basisCenter () const { return m_p; }
    /*! \brief Writing access to the (evaluation position) */
    PONCA_MULTIARCH inline       VectorType& basisCenter ()       { return m_p; }


    /*! \brief Init the line from a direction and a position
       \param _dir Orientation of the line, does not need to be normalized
       \param _pos Position of the line
    */
    PONCA_MULTIARCH inline void setLine (const VectorType& origin,
                                         const VectorType& direction)
    {
        EigenBase* cc = static_cast<EigenBase*>(this);
        *cc = EigenBase(origin, direction);
    }

    /*! \brief Value of the scalar field at the evaluation point */
    PONCA_MULTIARCH inline Scalar potential ( ) const
    {
        // The potential is the distance from a point to the line
        return EigenBase::squaredDistance( m_p);
    }

    /*!  \brief Value of the scalar field at the location \f$ \mathbf{q} \f$,
     * defined as the squared distance between \f$ \mathbf{q} \f$ and the line
     */
    PONCA_MULTIARCH inline Scalar potential (const VectorType& _q) const
    {
        // The potential is the distance from a point to the line
        return EigenBase::squaredDistance(_q - m_p);
    }

    //! \brief Project a point on the line
    PONCA_MULTIARCH inline VectorType project (const VectorType& _q) const
    {
        // Project on the normal vector and add the offset value
        return EigenBase::projection(_q - m_p) + m_p;
    }
}; //class Line


}
