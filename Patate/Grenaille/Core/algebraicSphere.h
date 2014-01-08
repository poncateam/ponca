/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/. 
*/


#ifndef _GRENAILLE_ALGEBRAIC_SPHERE_
#define _GRENAILLE_ALGEBRAIC_SPHERE_

namespace Grenaille
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
    
    This fitting procedure provides: 
    \verbatim PROVIDES_ALGEBRAIC_SPHERE \endverbatim

    \todo Deal with planar case (_uq == 0) and what about _ul == 0 ?
   */
  template < class DataPoint, class _WFunctor, typename T = void  >
  class AlgebraicSphere {
  protected:
    enum
      {
        PROVIDES_ALGEBRAIC_SPHERE /*!< \brief Provides Algebraic Sphere */
      };

  public:
    /*! \brief Scalar type inherited from DataPoint*/
    typedef typename DataPoint::Scalar     Scalar;     
    /*! \brief Vector type inherited from DataPoint*/
    typedef typename DataPoint::VectorType VectorType;
    /*! \brief Weight Function*/
    typedef _WFunctor                      WFunctor;  
    
  private:    
    //! \brief Evaluation position (needed for centered basis)
    VectorType _p; 
    
  protected:
    //! \brief Is the implicit scalar field normalized using Pratt
    bool _isNormalized;
    
    //! \brief Represent the current state of the fit (finalize function update the state)
    FIT_RESULT _eCurrentState;

	//! \brief Give the number of neighbors
	int _nbNeighbors;

    // results
  public:
    Scalar _uc,       /*!< \brief Constant parameter of the Algebraic hyper-sphere */
           _uq;       /*!< \brief Quadratic parameter of the Algebraic hyper-sphere */
    VectorType _ul;   /*!< \brief Linear parameter of the Algebraic hyper-sphere */
    
  public:
    /*! \brief Default constructor */
    MULTIARCH inline AlgebraicSphere(){
      _p = VectorType::Zero();
      resetPrimitive();
    }    
    
    /*! \brief Set the scalar field values to 0 and reset the isNormalized() status */
    MULTIARCH inline void resetPrimitive() {
      _uc = Scalar(0.0);
      _ul = VectorType::Zero();
      _uq = Scalar(0.0);
      
      _isNormalized = false;
      _eCurrentState = UNDEFINED;
	  _nbNeighbors = 0;
    }
    
    /*! \brief Is the sphere fitted an ready to use (finalize has been called)
		\warning The fit can be unstable (having neighbors between 3 and 6) */
    MULTIARCH inline bool isReady() const { return (_eCurrentState == STABLE) || (_eCurrentState == UNSTABLE); }

	/*! \brief Is the sphere fitted an ready to use (finalize has been called and the result is stable, eq. having more than 6 neighbors) */
    MULTIARCH inline bool isStable() const { return _eCurrentState == STABLE; }

	/*! \return the current test of the fit */
	MULTIARCH inline FIT_RESULT getCurrentState() const { return _eCurrentState; }

    /*! \brief Reading access to the basis center (evaluation position) */
    MULTIARCH inline const VectorType& basisCenter () const { return _p; }
    /*! \brief Writing access to the (evaluation position) */
    MULTIARCH inline       VectorType& basisCenter ()       { return _p; }

    /*! \brief compute the Pratt norm of the implicit scalar field. */
    MULTIARCH inline Scalar prattNorm() const {
      MULTIARCH_STD_MATH(sqrt);
      return sqrt(prattNorm2());
    }
    
    /*! \brief compute the squared Pratt norm of the implicit scalar field. */
    MULTIARCH inline Scalar prattNorm2() const {
      return _ul.squaredNorm() - Scalar(4.) * _uc*_uq;
    }

    /*!
       \brief Normalize the scalar field by the Pratt norm
       \return false when the normalization fails (sphere is already normalized)
     */
    MULTIARCH inline bool applyPrattNorm() {
      if (! _isNormalized){
        Scalar pn = prattNorm();
        _uc /= pn;
        _ul *= Scalar(1.)/pn;
        _uq /= pn;

        _isNormalized = true;
      }
      return true;
    }
    
	/*! 
		\brief return the estimated radius of the sphere
		\return inf if the fitted surface is planar
		\warning return +inf if the fitted surface is planar
	*/
    MULTIARCH inline Scalar radius()
	{
	  if(isPlane())
	  {
	    return std::numeric_limits<Scalar>::infinity();
	  }

      MULTIARCH_STD_MATH(sqrt);
      Scalar b = 1./_uq;
      return sqrt( ((-0.5*b)*_ul).squaredNorm() - _uc*b );
    }
    
	/*! 
		\brief return the estimated center of the sphere
	*/
    MULTIARCH inline VectorType center()
	{
	  if(isPlane())
	  {
	    return VectorType(std::numeric_limits<Scalar>::infinity());
	  }

      Scalar b = 1./_uq;
      return (-0.5*b)*_ul + basisCenter();
    }
    
    //! \brief State indicating when the sphere has been normalized 
    MULTIARCH inline bool isNormalized() const { return _isNormalized; }
    
    //! \brief Value of the scalar field at the location \f$ \mathbf{q} \f$
    MULTIARCH inline Scalar potential (const VectorType& q) const;
    
    //! \brief Project a point on the sphere
    MULTIARCH inline VectorType project (const VectorType& q) const;
    
    //! \brief Approximation of the scalar field gradient at \f$ \mathbf{q} (not normalized) \f$
    MULTIARCH inline VectorType primitiveGradient (const VectorType& q) const;

	/*! 
		\brief Used to know if the fitting result to a plane
		\return true if finalize() have been called and the fitting result to a plane
	*/
	MULTIARCH inline bool isPlane() const
	{
		Scalar epsilon = Eigen::NumTraits<Scalar>::dummy_precision();
		bool bPlanar = Eigen::internal::isMuchSmallerThan(std::abs(_uq), 1., epsilon);
		bool bReady = isReady();

		if(bReady && bPlanar)
		{
			return true;
		}

		return false;
	}
	
  }; //class AlgebraicSphere


#include "algebraicSphere.hpp"

}
#endif  // _GRENAILLE_ALGEBRAIC_SPHERE_
