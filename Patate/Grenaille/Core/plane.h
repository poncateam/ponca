/*
 Copyright (C) 2014 Nicolas Mellado <nmellado0@gmail.com>
 
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#ifndef _GRENAILLE_PLANE_
#define _GRENAILLE_PLANE_

namespace Grenaille
{

/*!
    \brief Implicit hyperplane defined by an homogeneous vector \f$\mathbf{p}\f$. 
    
    In n-dimensionnal space, the plane is defined as
    the \f$0\f$-isosurface of the scalar field
    
    \f$ s_\mathbf{u}(\mathbf{x}) = 
    \left[ \mathbf{x}^T \; 1 \;\right]^T \cdot \mathbf{p} \f$.
    
    
    This primitive requires the definition of n-dimensionnal vectors 
    (VectorType) and homogeneous n-dimensionnal vectors (HVectorType) in 
    Concept::PointConcept.
    
    This primitive provides: 
    \verbatim PROVIDES_PLANE \endverbatim
    
    
    \TODO Use Generic Primitive class and push fit state into it
    
    \note The first n-components of the plane must define a normalized vector
    
*/
template < class DataPoint, class _WFunctor, typename T = void  >
class Plane
{
protected:

    enum
    {
        PROVIDES_PLANE /*!< \brief Provides a Plane primitive */
    };

public:

    /*! \brief Scalar type inherited from DataPoint */
    typedef typename DataPoint::Scalar      Scalar;     
    /*! \brief Vector type inherited from DataPoint */
    typedef typename DataPoint::VectorType  VectorType;
    /*! \brief Homogeneous vector type inherited from DataPoint */
    typedef typename DataPoint::HVectorType HVectorType;
    /*! \brief Weight Function */
    typedef _WFunctor                       WFunctor;  
    
protected:
    
    //! \brief Represent the current state of the fit (finalize function 
    //! update the state)
    FIT_RESULT m_eCurrentState;

    //! \brief Give the number of neighbors
    int m_nbNeighbors;

// results
public:

    HVectorType m_p; /*!< \brief Homogeneous plane representation */
    
public:

    /*! \brief Default constructor */
    MULTIARCH inline Plane()
    {
        resetPrimitive();
    }    
    
    /*! \brief Set the scalar field values to 0 and reset the isNormalized() 
         status */
    MULTIARCH inline void resetPrimitive()
    {
        m_p = HVectorType::Zero();
      
        m_eCurrentState = UNDEFINED;
        m_nbNeighbors = 0;
    }
    
    /*! \brief Is the plane fitted an ready to use (finalize has been called)
    \warning The fit can be unstable (having neighbors between 3 and 6) */
    MULTIARCH inline bool isReady() const 
    { 
        return (m_eCurrentState == STABLE) || (m_eCurrentState == UNSTABLE); 
    }

    /*! \brief Is the plane fitted an ready to use (finalize has been called 
    and the result is stable, eq. having more than 6 neighbors) */
    MULTIARCH inline bool isStable() const { return m_eCurrentState == STABLE; }

    /*! \return the current test of the fit */
    MULTIARCH inline FIT_RESULT getCurrentState() const 
    { 
        return m_eCurrentState; 
    }
        
    //! \brief Value of the scalar field at the location \f$ \mathbf{q} \f$
    MULTIARCH inline Scalar potential (const VectorType& _q) const
    {
        // Project on the normal vector and add the offset value
        return    m_p.template head<DataPoint::Dim>().dot(_q) + 
               *( m_p.template tail<1>().data() );
    }
    
    //! \brief Project a point on the plane
    MULTIARCH inline VectorType project (const VectorType& _q) const
    {
        // The potential is the distance from the point to the plane
        return _q - potential(_q) * m_p.template head<DataPoint::Dim>();
    }
    
    //! \brief Scalar field gradient direction at \f$ \mathbf{q}\f$
    MULTIARCH inline VectorType primitiveGradient (const VectorType&) const
    {
        // Uniform gradient defined only by the orientation of the plane
        return m_p.template head<DataPoint::Dim>();
    }


}; //class Plane

}
#endif  // _GRENAILLE_PLANE_
