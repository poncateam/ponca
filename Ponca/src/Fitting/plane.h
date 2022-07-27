/*
 Copyright (C) 2014 Nicolas Mellado <nmellado0@gmail.com>

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "./defines.h"
#include <Eigen/Geometry>

namespace Ponca
{

/*!
    \brief Implicit hyperplane defined by an homogeneous vector \f$\mathbf{p}\f$.

    In n-dimensionnal space, the plane is defined as
    the \f$0\f$-isosurface of the scalar field

    \f$ s_\mathbf{u}(\mathbf{x}) =
    \left[ \mathbf{x}^T \; 1 \;\right]^T \cdot \mathbf{p} \f$.

    This class inherits Eigen::Hyperplane.

    This primitive requires the definition of n-dimensionnal vectors
    (VectorType) in Concept::PointConcept.

    This primitive provides:
    \verbatim PROVIDES_PLANE \endverbatim

    \ingroup fitting

*/
template < class DataPoint, class _WFunctor, typename T >
class Plane : public T,
              public Eigen::Hyperplane<typename DataPoint::Scalar, DataPoint::Dim >
{
private:
    using Base = T;

public:
    /// \brief Specialization of Eigen::Hyperplane inherited by Ponca::Plane
    using EigenBase = Eigen::Hyperplane<typename DataPoint::Scalar, DataPoint::Dim >;

protected:

    enum
    {
        check = Base::PROVIDES_PRIMITIVE_BASE,  /*!< \brief Requires PrimitiveBase */
        PROVIDES_PLANE                          /*!< \brief Provides a Plane primitive */
    };

public:
    using Scalar     = typename Base::Scalar;     /*!< \brief Inherited scalar type*/
    using VectorType = typename Base::VectorType; /*!< \brief Inherited vector type*/
    using WFunctor   = typename Base::WFunctor;   /*!< \brief Weight Function*/


public:

    /*! \brief Default constructor */
    PONCA_MULTIARCH inline Plane()
        : Base(), EigenBase()
    {
        init(VectorType::Zero());
    }

    PONCA_EXPLICIT_CAST_OPERATORS(Plane,compactPlane) //< \fixme To be removed, kept for compatibility only
    PONCA_EXPLICIT_CAST_OPERATORS(Plane,plane)

    /// \brief Set the scalar field values to 0
    PONCA_MULTIARCH inline void init(const VectorType& _basisCenter = VectorType::Zero())
    {
        Base::init(_basisCenter);
        EigenBase::coeffs().setZero();
    }
    }

    PONCA_MULTIARCH inline bool operator==(const Plane<DataPoint, WFunctor, T>& other) const{
        return EigenBase::isApprox(other);
    }

    /*! \brief Comparison operator, convenience function */
    PONCA_MULTIARCH inline bool operator!=(const Plane<DataPoint, WFunctor, T>& other) const{
        return ! ((*this) == other);
    }

    /* \brief Init the plane from a direction and a position
       \param _dir Orientation of the plane, does not need to be normalized
       \param _pos Position of the plane
    */
    PONCA_MULTIARCH inline void setPlane (const VectorType& _dir,
                                    const VectorType& _pos)
    {
        EigenBase* cc = static_cast<EigenBase*>(this);
        *cc = EigenBase(_dir.normalized(), _pos);
    }

    /*! \brief Value of the scalar field at the evaluation point */
    PONCA_MULTIARCH inline Scalar potential ( ) const
    {
        return EigenBase::signedDistance(VectorType::Zero());
    }

    //! \brief Value of the scalar field at the location \f$ \mathbf{q} \f$
    PONCA_MULTIARCH inline Scalar potential (const VectorType& _q) const
    {
        // The potential is the distance from the point to the plane
        return EigenBase::signedDistance(Base::m_w.convertToLocalBasis(_q) );
    }

    //! \brief Project a point on the plane
    PONCA_MULTIARCH inline VectorType project (const VectorType& _q) const
    {
        // Project on the normal vector and add the offset value
        return EigenBase::projection(Base::m_w.convertToLocalBasis(_q)) + Base::m_w.basisCenter();
    }

    //! \brief Scalar field gradient direction at the evaluation point
    PONCA_MULTIARCH inline VectorType primitiveGradient () const
    {
        // Uniform gradient defined only by the orientation of the plane
        return EigenBase::normal();
    }

    //! \brief Scalar field gradient direction at \f$ \mathbf{q}\f$
    PONCA_MULTIARCH inline VectorType primitiveGradient (const VectorType&) const
    {
        // Uniform gradient defined only by the orientation of the plane
        return EigenBase::normal();
    }


}; //class Plane

}
