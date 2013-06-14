/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/. 
*/


#ifndef _GRENAILLE_GLS_
#define _GRENAILLE_GLS_


namespace Grenaille
{

  /*!
    \brief Growing Least Squares reparemetrization of the OrientedSphereFit
    \inherit FittingExtensionInterface
    
    Method published in \cite Mellado:2012:GLS
        
    This class assumes that the WeightFunc defines the accessor
    \code
    w.evalScale();
    \endcode
    in order to access to the evaluation scale, needed to compute the 
    scale invariant GLS reparametrization (all *_normalized methods).
    
    Computed values:
     - tau(), eta() and kappa(): the GLS descriptor 
     \f$ \left[ \tau \; \eta \; \kappa \right]\f$ 
     - tau_normalized(), eta_normalized() and kappa_normalized(): 
     the scale invariant GLS descriptor
     \f$ \left[ \frac{\tau}{t} \; \eta \; t\kappa \right]\f$ 
          
    Requierement: 
    \verbatim PROVIDES_ALGEBRAIC_SPHERE \endverbatim
    Provide: 
    \verbatim PROVIDES_GLS_PARAMETRIZATION \endverbatim
   */
  template < class DataPoint, class _WFunctor, typename T>
  class GLSParam : public T{
  private:
    typedef T Base;

  protected:
    enum
      {
        Check = Base::PROVIDES_ALGEBRAIC_SPHERE,
        PROVIDES_GLS_PARAMETRIZATION
      };

  public:
    typedef typename Base::Scalar     Scalar;     /*!< \brief Inherited scalar type*/
    typedef typename Base::VectorType VectorType; /*!< \brief Inherited vector type*/
    typedef typename Base::WFunctor   WFunctor;   /*!< \brief Weight Function*/


  protected:
    Scalar _t; /*!< \brief Evaluation scale. Needed to computed the normalized values*/
    Scalar _fitness; /*!< \brief Save the fitness value to avoid side effect with Pratt normalization*/

  public:
    /*! \brief Default constructor */
    MULTIARCH inline GLSParam() : _t(0) {}

    
    /**************************************************************************/
    /* Initialization                                                         */
    /**************************************************************************/
    /*! \copydoc FittingProcedureInterface::setWeightFunc() */
    MULTIARCH inline void setWeightFunc (const WFunctor& w){
      Base::setWeightFunc(w);
      _t = w.evalScale();
    }
    
    
    /**************************************************************************/
    /* Processing                                                             */
    /**************************************************************************/
    /*! \copydoc FittingProcedureInterface::finalize() */
    MULTIARCH inline void finalize   (){
      Base::finalize();
      _fitness = Scalar(1.) - Base::prattNorm2();
    }


    /**************************************************************************/
    /* Use results                                                            */
    /**************************************************************************/
    /*! \brief Compute and return \f$ \tau \f$ */
    MULTIARCH inline Scalar     tau()   const 
    {return Base::isNormalized() ? Base::_uc : Base::_uc / Base::prattNorm();}

    /*! \brief Compute and return \f$ \eta \f$ */
    MULTIARCH inline VectorType eta()   const 
    {return Base::_ul / Base::_ul.norm();}

    /*! \brief Compute and return \f$ \kappa \f$ */
    MULTIARCH inline Scalar     kappa() const 
    {return Scalar(2.) * (Base::isNormalized() ? Base::_uq : Base::_uq / Base::prattNorm());}
    
    /*! \brief Compute and return \f$ \frac{\tau}{t} \f$ */
    MULTIARCH inline Scalar     tau_normalized()   const {return tau()/_t;}
    
    /*! \brief Compute and return \f$ \eta \f$ */    
    MULTIARCH inline VectorType eta_normalized()   const {return eta();}
    
    /*! \brief Compute and return \f$ t \kappa \f$ */
    MULTIARCH inline Scalar     kappa_normalized() const {return kappa()*_t;}    
    
    /*! \brief Return the fitness, e.g. the pratt norm of the initial scalar field */
    MULTIARCH inline Scalar     fitness()          const {return _fitness;}

  }; //class GLSParam

  
  /*!
    \brief Differentiation of GLSParam
    \inherit FittingExtensionInterface

    Method published in \cite Mellado:2012:GLS
   */
  template < class DataPoint, class _WFunctor, typename T>
  class GLSDer : public T{
  private:
    typedef T Base;

  protected:
    enum
      {
        Check = Base::PROVIDES_GLS_PARAMETRIZATION,
        PROVIDES_GLS_DERIVATIVE
      };

  public:
    typedef typename Base::Scalar     Scalar;
    typedef typename Base::VectorType VectorType;
    typedef typename Base::WFunctor   WFunctor;

    typedef typename Base::VectorArray VectorArray;
    typedef typename Base::ScalarArray ScalarArray;   
    
    MULTIARCH inline ScalarArray dtau()   const;
    MULTIARCH inline VectorArray deta()   const;
    MULTIARCH inline ScalarArray dkappa() const;
    
    MULTIARCH inline ScalarArray dtau_normalized()   const;
    MULTIARCH inline VectorArray deta_normalized()   const;
    MULTIARCH inline ScalarArray dkappa_normalized() const;
  }; //class GLSScaleDer


  /*!
    \brief Extension to compute the Geometric Variation of GLSParam
    \inherit FittingExtensionInterface
    
    Method published in \cite Mellado:2012:GLS
    \todo Add more details
   */
  template < class DataPoint, class _WFunctor, typename T>
  class GLSGeomVar : public T{
  private:
    typedef T Base;

  protected:
    enum
      {
        Check = Base::PROVIDES_ALGEBRAIC_SPHERE_SCALE_DERIVATIVE & Base::PROVIDES_GLS_DERIVATIVE,
	      PROVIDES_GLS_GEOM_VAR
      };

  public:
    typedef typename Base::Scalar Scalar;

    /*!
      \brief Compute and return the Geometric Variation
     */
    MULTIARCH inline Scalar geomVar(Scalar wtau   = Scalar(1), 
  			                            Scalar weta   = Scalar(1),
	                        			    Scalar wkappa = Scalar(1)) const;
  };

  #include "gls.hpp"

} //namespace Grenaille

#endif
