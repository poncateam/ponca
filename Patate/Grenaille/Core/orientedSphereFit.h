/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/. 
*/


#ifndef _GRENAILLE_ORIENTED_SPHERE_FIT_
#define _GRENAILLE_ORIENTED_SPHERE_FIT_

#include "internal.h"

namespace Grenaille
{
  /*!
    \brief Algebraic Sphere fitting procedure on oriented point sets
    
    An algebraic hyper-sphere is defined as the \f$0\f$-isosurface of the scalar field
    
    \f$ s_\mathbf{u}(\mathbf{x}) = \left[ 1 \; \mathbf{x}^T \; \mathbf{x}^T\mathbf{x}\right]^T \cdot \mathbf{u} \f$    
    
    with \f$ \mathbf{u} \left[ u_c \; \mathbf{u_l} \; u_q\right]^T \f$ is the 
    vector of the constant, linear and quadratic parameters.
    
    
    
    This fitting procedure provides: 
    \verbatim PROVIDES_ALGEBRAIC_SPHERE \endverbatim

    \todo Deal with planar case
   */
  template < class DataPoint, class _WFunctor, typename T = void >
  class OrientedSphereFit {
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
    
  protected:
    //! \brief Evaluation position (needed for centered basis)
    VectorType _p;
    
    // computation data
    VectorType _sumN, /*!< \brief Sum of the normal vectors */
               _sumP; /*!< \brief Sum of the relative positions */
    Scalar _sumDotPN, /*!< \brief Sum of the dot product betwen relative positions and normals */
           _sumDotPP, /*!< \brief Sum of the squared relative positions */
           _sumW;     /*!< \brief Sum of queries weight */
    
    WFunctor _w;      /*!< \brief Weight function (must inherits BaseWeightFunc) */

    //! Is the implicit scalar field normalized using Pratt
    bool _isNormalized;

    // results
  public:
    Scalar _uc,       /*!< \brief Constant parameter of the Algebraic hyper-sphere */
           _uq;       /*!< \brief Quadratic parameter of the Algebraic hyper-sphere */
    VectorType _ul;   /*!< \brief Linear parameter of the Algebraic hyper-sphere */
    
  public:
    /*! \brief Default constructor */
    MULTIARCH inline OrientedSphereFit(){ }

    /*! \brief Reading access to the evaluation position */
    MULTIARCH inline const VectorType& evalPos () const { return _p; }
    /*! \brief Writing access to the evaluation position */
    MULTIARCH inline       VectorType& evalPos ()       { return _p; }
    
    
    /**************************************************************************/
    /* Initialization                                                         */
    /**************************************************************************/
    /*! \copydoc FittingProcedureInterface::setWeightFunc() */
    MULTIARCH inline void setWeightFunc (const WFunctor& w) { _w  = w; }
    
    /*! \copydoc FittingProcedureInterface::init() */
    MULTIARCH inline void init (const VectorType& evalPos);
    

    /**************************************************************************/
    /* Processing                                                             */
    /**************************************************************************/
    /*! \copydoc FittingProcedureInterface::addNeighbor() */
    MULTIARCH inline void addNeighbor(const DataPoint &nei);
    
    /*! \copydoc FittingProcedureInterface::finalize() */
    MULTIARCH inline void finalize   ();


    /**************************************************************************/
    /* Use results                                                            */
    /**************************************************************************/
    /*! \brief compute the Pratt norm of the implicit scalar field. */
    MULTIARCH inline Scalar prattNorm() const {
      MULTIARCH_STD_MATH(sqrt);
      return sqrt(prattNorm2());
    }
    
    /*! \brief compute the squared Pratt norm of the implicit scalar field. */
    MULTIARCH inline Scalar prattNorm2() const {
      return _ul.squaredNorm() - Scalar(4.) * _uc*_uq;
    }

    //! Normalize the scalar field by the Pratt norm
    /*!
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
    
    //! State that indicates if the sphere has been normalized 
    MULTIARCH inline bool isNormalized() const { return _isNormalized; }
    
    //! Project a point on the sphere
    MULTIARCH inline VectorType project (VectorType q);

    
    //    MULTIARCH VectorType gradient(VectorType q, bool normalize = true);

  }; //class OrientedSphereFit


  namespace internal{

    enum {
      FitScaleDer = 0x01, /*!< \brief Flag indicating a scale differentiation. */
      FitSpaceDer = 0x02  /*!< \brief Flag indicating a space differentiation. */
    };

    /*! 
      \brief Internal generic class performing the Fit derivation 
      \inherit FittingExtensionInterface
      
      The differentiation can be done automatically in scale and/or space, by
      combining the enum values FitScaleDer and FitSpaceDer in the template 
      parameter Type.
      
      The differenciated values are stored in static arrays. The size of the
      arrays is computed with respect to the derivation type (scale and/or space)
      and the number of the dimension of the ambiant space.      
      By convention, the scale derivatives are stored at index 0 when Type 
      contains at least FitScaleDer. The size of these arrays can be known using
      derDimension(), and the differentiation type by isScaleDer() and 
      isSpaceDer().
      
    */
    template < class DataPoint, class _WFunctor, typename T, int Type>
    class OrientedSphereDer : public T{
    private:
      typedef T Base; /*!< \brief Generic base type */

    protected:
      enum
	{
	  Check = Base::PROVIDES_ALGEBRAIC_SPHERE, /*!< \brief Needs Algebraic Sphere */
	  PROVIDES_ALGEBRAIC_SPHERE_DERIVATIVE     /*!< \brief Provides Algebraic Sphere derivative*/
	};

    public:
      typedef typename Base::Scalar     Scalar;     /*!< \brief Inherited scalar type*/
      typedef typename Base::VectorType VectorType; /*!< \brief Inherited vector type*/
      typedef typename Base::WFunctor   WFunctor;   /*!< \brief Weight Function*/

#define GLS_DER_NB_DERIVATIVES(TYPE,DIM) ((TYPE & FitScaleDer) ? 1 : 0 ) + ((TYPE & FitSpaceDer) ? DIM : 0)
      /*! \brief Static array of scalars with a size adapted to the differentiation type */
      typedef FixedSizeArray <VectorType, GLS_DER_NB_DERIVATIVES(Type,DataPoint::Dim)> VectorArray;
      /*! \brief Static array of scalars with a size adapted to the differentiation type */
      typedef FixedSizeArray <Scalar,     GLS_DER_NB_DERIVATIVES(Type,DataPoint::Dim)> ScalarArray;
      
    private:
      // computation data
      VectorArray _dSumN,     /*!< \brief Sum of the normal vectors with differenciated weights */
                  _dSumP;     /*!< \brief Sum of the relative positions with differenciated weights*/
      ScalarArray _dSumDotPN, /*!< \brief Sum of the dot product betwen relative positions and normals with differenciated weights */
                  _dSumDotPP, /*!< \brief Sum of the squared relative positions with differenciated weights */
                  _dSumW;     /*!< \brief Sum of queries weight with differenciated weights */
      
    public:
      // results
      ScalarArray _dUc, /*!< \brief Derivative of the hyper-sphere constant term  */
                  _dUq; /*!< \brief Derivative of the hyper-sphere quadratic term */
      VectorArray _dUl; /*!< \brief Derivative of the hyper-sphere linear term    */
    
      /************************************************************************/
      /* Initialization                                                       */
      /************************************************************************/
      /*! \see FittingProcedureInterface::init() */
      MULTIARCH void init       (const VectorType &evalPos);
    
      /************************************************************************/
      /* Processing                                                           */
      /************************************************************************/
      /*! \see FittingProcedureInterface::addNeighbor() */
      MULTIARCH void addNeighbor(const DataPoint  &nei);
      /*! \see FittingProcedureInterface::finalize() */
      MULTIARCH void finalize   ();


    /**************************************************************************/
    /* Use results                                                            */
    /**************************************************************************/
      /*! \brief compute the square of the Pratt norm derivative for dimension d */
      MULTIARCH inline Scalar dprattNorm2(unsigned int d) const {
        return   Scalar(2.) * _dUl[d].dot(Base::_ul) 
               - Scalar(4.) * _dUc[d]*Base::_uq
	             - Scalar(4.) * Base::_uc*_dUq[d];}

      /*! \brief compute the Pratt norm derivative for the dimension d */
      MULTIARCH inline Scalar dprattNorm(unsigned int d) const {
	      MULTIARCH_STD_MATH(sqrt);
        return sqrt(dprattNorm2(d));
      }

      /*! \brief State specified at compilation time to differenciate the fit in scale */
      MULTIARCH inline bool isScaleDer() const {return Type & FitScaleDer;}
      /*! \brief State specified at compilation time to differenciate the fit in space */
      MULTIARCH inline bool isSpaceDer() const {return Type & FitSpaceDer;}
      /*! \brief Number of dimensions used for the differentiation */
      MULTIARCH inline unsigned int derDimension() const { return VectorArray::size();}


      //! Normalize the scalar field by the Pratt norm
      /*!
	      \warning Requieres that isNormalized() return false
	      \return false when the original sphere has already been normalized.
       */
      MULTIARCH inline bool applyPrattNorm();

    }; //class OrientedSphereFitDer

  }// namespace internal  

  /*!
    \brief Differentiation in scale of the OrientedSphereFit
    \inherit FittingExtensionInterface
    
    Requierement: 
    \verbatim PROVIDES_ALGEBRAIC_SPHERE \endverbatim
    Provide: 
    \verbatim PROVIDES_ALGEBRAIC_SPHERE_SCALE_DERIVATIVE \endverbatim
  */
  template < class DataPoint, class _WFunctor, typename T>
  class OrientedSphereScaleDer:public internal::OrientedSphereDer<DataPoint, _WFunctor, T, internal::FitScaleDer>
  {
  protected:
    /*! \brief Inherited class */
    typedef internal::OrientedSphereDer<DataPoint, _WFunctor, T, internal::FitScaleDer> Base;
    enum { PROVIDES_ALGEBRAIC_SPHERE_SCALE_DERIVATIVE };
  };


  /*!
    \brief Spatial differentiation of the OrientedSphereFit
    \inherit FittingExtensionInterface
    
    Requierement: 
    \verbatim PROVIDES_ALGEBRAIC_SPHERE \endverbatim
    Provide: 
    \verbatim PROVIDES_ALGEBRAIC_SPHERE_SPACE_DERIVATIVE \endverbatim
  */
  template < class DataPoint, class _WFunctor, typename T>
  class OrientedSphereSpaceDer:public internal::OrientedSphereDer<DataPoint, _WFunctor, T, internal::FitSpaceDer>
  {
  protected:
    /*! \brief Inherited class */
    typedef internal::OrientedSphereDer<DataPoint, _WFunctor, T, internal::FitSpaceDer> Base;
    enum {  PROVIDES_ALGEBRAIC_SPHERE_SPACE_DERIVATIVE };
  };


  /*!
    \brief Differentiation both in scale and space of the OrientedSphereFit
    \inherit FittingExtensionInterface
    
    Requierement: 
    \verbatim PROVIDES_ALGEBRAIC_SPHERE \endverbatim
    Provide: 
    \verbatim PROVIDES_ALGEBRAIC_SPHERE_SCALE_DERIVATIVE 
    PROVIDES_ALGEBRAIC_SPHERE_SPACE_DERIVATIVE
    \endverbatim
  */
  template < class DataPoint, class _WFunctor, typename T>
  class OrientedSphereScaleSpaceDer:public internal::OrientedSphereDer<DataPoint, _WFunctor, T, internal::FitSpaceDer | internal::FitScaleDer>
  {
  protected:
    /*! \brief Inherited class */
    typedef internal::OrientedSphereDer<DataPoint, _WFunctor, T, internal::FitSpaceDer | internal::FitScaleDer> Base;
    enum
      {
        PROVIDES_ALGEBRAIC_SPHERE_SCALE_DERIVATIVE,
        PROVIDES_ALGEBRAIC_SPHERE_SPACE_DERIVATIVE
      };
  };


#include "orientedSphereFit.hpp"

} //namespace Grenaille


#endif
