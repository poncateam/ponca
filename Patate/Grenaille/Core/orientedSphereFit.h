/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/. 
*/


#ifndef _GRENAILLE_ORIENTED_SPHERE_FIT_
#define _GRENAILLE_ORIENTED_SPHERE_FIT_

#include "internal.h"
#include "algebraicSphere.h"

namespace Grenaille
{
  /*!
    \brief Algebraic Sphere fitting procedure on oriented point sets
    
    Method published in \cite Guennebaud:2007:APSS.
    
    \inherit Concept::FittingProcedureConcept
    
    \see AlgebraicSphere
   */
  template < class DataPoint, class _WFunctor, typename T = void >
  class OrientedSphereFit : public AlgebraicSphere<DataPoint, _WFunctor>{
    private:
      typedef AlgebraicSphere<DataPoint, _WFunctor> Base; 
  public:
    /*! \brief Scalar type inherited from DataPoint*/
    typedef typename Base::Scalar     Scalar;     
    /*! \brief Vector type inherited from DataPoint*/
    typedef typename Base::VectorType VectorType;
    /*! \brief Weight Function*/
    typedef _WFunctor                 WFunctor;   
    
  protected:
    
    // computation data
    VectorType _sumN, /*!< \brief Sum of the normal vectors */
               _sumP; /*!< \brief Sum of the relative positions */
    Scalar _sumDotPN, /*!< \brief Sum of the dot product betwen relative positions and normals */
           _sumDotPP, /*!< \brief Sum of the squared relative positions */
           _sumW;     /*!< \brief Sum of queries weight */
    
    WFunctor _w;      /*!< \brief Weight function (must inherits BaseWeightFunc) */

    
  public:
    /*! \brief Default constructor */
    MULTIARCH inline OrientedSphereFit()
    : AlgebraicSphere<DataPoint,WFunctor>(){ }
    
    
    /**************************************************************************/
    /* Initialization                                                         */
    /**************************************************************************/
    /*! \copydoc Concept::FittingProcedureConcept::setWeightFunc() */
    MULTIARCH inline void setWeightFunc (const WFunctor& w) { _w  = w; }
    
    /*! \copydoc Concept::FittingProcedureConcept::init() */
    MULTIARCH inline void init (const VectorType& evalPos);
    

    /**************************************************************************/
    /* Processing                                                             */
    /**************************************************************************/
    /*! \copydoc Concept::FittingProcedureConcept::addNeighbor() */
    MULTIARCH inline void addNeighbor(const DataPoint &nei);
    
    /*! \copydoc Concept::FittingProcedureConcept::finalize() */
    MULTIARCH inline void finalize   ();
  }; //class OrientedSphereFit


  namespace internal{

    enum {
      FitScaleDer = 0x01, /*!< \brief Flag indicating a scale differentiation. */
      FitSpaceDer = 0x02  /*!< \brief Flag indicating a space differentiation. */
    };

    /*! 
      \brief Internal generic class performing the Fit derivation 
      \inherit Concept::FittingExtensionConcept
      
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
      typedef typename Base::MatrixType MatrixType; /*!< \brief Inherited Matrix type*/
      typedef typename Base::WFunctor   WFunctor;   /*!< \brief Weight Function*/

#define GLS_DER_NB_DERIVATIVES(TYPE,DIM) ((TYPE & FitScaleDer) ? 1 : 0 ) + ((TYPE & FitSpaceDer) ? DIM : 0)
#define GLS_DER_STORAGE_ORDER(TYPE)      ((TYPE & FitSpaceDer) ? Eigen::RowMajor : Eigen::ColMajor )

      /*! \brief Static array of scalars with a size adapted to the differentiation type */
      typedef Eigen::Matrix < Scalar, 
                              DataPoint::Dim, 
                              GLS_DER_NB_DERIVATIVES(Type,DataPoint::Dim),
                              GLS_DER_STORAGE_ORDER(Type) > VectorArray;
      /*! \brief Static array of scalars with a size adapted to the differentiation type */
      typedef Eigen::Matrix < Scalar, 
                              1, 
                              GLS_DER_NB_DERIVATIVES(Type,DataPoint::Dim)/*,
                              GLS_DER_STORAGE_ORDER(Type)*/ > ScalarArray;
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
      /*! \see Concept::FittingProcedureConcept::init() */
      MULTIARCH void init       (const VectorType &evalPos);
    
      /************************************************************************/
      /* Processing                                                           */
      /************************************************************************/
      /*! \see Concept::FittingProcedureConcept::addNeighbor() */
      MULTIARCH void addNeighbor(const DataPoint  &nei);
      /*! \see Concept::FittingProcedureConcept::finalize() */
      MULTIARCH void finalize   ();


    /**************************************************************************/
    /* Use results                                                            */
    /**************************************************************************/
      /*! \brief compute  the square of the Pratt norm derivative */
	  MULTIARCH inline ScalarArray dprattNorm2() const{
        return   Scalar(2.) * Base::_ul.transpose() * _dUl 
               - Scalar(4.) * Base::_uq * _dUc
	             - Scalar(4.) * Base::_uc * _dUq;      
      }
    
      /*! \brief compute the square of the Pratt norm derivative for dimension d */
      MULTIARCH inline Scalar dprattNorm2(unsigned int d) const {
        return   Scalar(2.) * _dUl.col(d).dot(Base::_ul) 
               - Scalar(4.) * _dUc.col(d)[0]*Base::_uq
	             - Scalar(4.) * _dUq.col(d)[0]*Base::_uc;}

      /*! \brief compute the Pratt norm derivative for the dimension d */
      MULTIARCH inline Scalar dprattNorm(unsigned int d) const {
	      MULTIARCH_STD_MATH(sqrt);
        return sqrt(dprattNorm2(d));
      }

      /*! \brief compute the Pratt norm derivative */
      MULTIARCH inline Scalar dprattNorm() const {
	      MULTIARCH_STD_MATH(sqrt);
        return dprattNorm2().array().sqrt();
      }

      /*! \brief State specified at compilation time to differenciate the fit in scale */
      MULTIARCH inline bool isScaleDer() const {return Type & FitScaleDer;}
      /*! \brief State specified at compilation time to differenciate the fit in space */
      MULTIARCH inline bool isSpaceDer() const {return Type & FitSpaceDer;}
      /*! \brief Number of dimensions used for the differentiation */
      MULTIARCH inline unsigned int derDimension() const { return GLS_DER_NB_DERIVATIVES(Type,DataPoint::Dim);}


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
    \inherit Concept::FittingExtensionConcept
    
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
    \inherit Concept::FittingExtensionConcept
    
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
    \inherit Concept::FittingExtensionConcept
    
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
