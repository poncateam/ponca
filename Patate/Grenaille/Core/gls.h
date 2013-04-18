#ifndef _GRENAILLE_GLS_
#define _GRENAILLE_GLS_


namespace Grenaille
{

  /*!
    \brief Growing Least Squares reparemetrization of the OrientedSphereFit
    
    Requierement: 
    \verbatim PROVIDES_ALGEBRAIC_SPHERE \endverbatim
    Provide: 
    \verbatim PROVIDES_GLS_PARAMETRIZATION \endverbatim

    
    This class assumes that the WeightFunc define the accessor
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
    typedef typename Base::Scalar     Scalar;
    typedef typename Base::VectorType VectorType;
    typedef typename Base::WFunctor   WFunctor;


  protected:
    Scalar _t;

  public:
    MULTIARCH inline GLSParam() : _t(0) {}

    MULTIARCH inline void setWeightFunc (const WFunctor& w){
      Base::setWeightFunc(w);
      _t = w.evalScale();
    }

    MULTIARCH inline Scalar     tau()   const 
    {return Base::isNormalized() ? Base::_uc : Base::_uc / Base::prattNorm();}

    MULTIARCH inline VectorType eta()   const 
    {return Base::_ul * ( Scalar(1.) / Base::_ul.norm());}

    MULTIARCH inline Scalar     kappa() const 
    {return Scalar(2.) * Base::_uq / Base::prattNorm();}
    
    MULTIARCH inline Scalar     tau_normalized()   const {return tau()/_t;}
    MULTIARCH inline VectorType eta_normalized()   const {return eta();}
    MULTIARCH inline Scalar     kappa_normalized() const {return kappa()*_t;}
    
    MULTIARCH inline Scalar     fitness() const {
        return Scalar(1.) - Base::_ul.squaredNorm() - Scalar(4.) * Base::_uc * Base::_uq;}

  }; //class GLSParam

  
  /*!
    \brief Differentiation of GLSParam

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

    MULTIARCH void finalize (){
      Base::finalize();
      Base::applyPrattNorm();
    }    
    
    MULTIARCH inline ScalarArray dtau()   const;
    MULTIARCH inline VectorArray deta()   const;
    MULTIARCH inline ScalarArray dkappa() const;
    
    MULTIARCH inline ScalarArray dtau_normalized()   const;
    MULTIARCH inline VectorArray deta_normalized()   const;
    MULTIARCH inline ScalarArray dkappa_normalized() const;
  }; //class GLSScaleDer


  /*!
    \brief Extension to compute the Geometric Variation of the GLSParam
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

    MULTIARCH inline Scalar geomVar(Scalar wtau   = Scalar(1), 
				    Scalar weta   = Scalar(1),
				    Scalar wkappa = Scalar(1)) const;
  };

  #include "gls.hpp"

} //namespace Grenaille

#endif
