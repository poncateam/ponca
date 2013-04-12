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


    MULTIARCH inline Scalar     tau()   const 
    {return Base::_uc / Base::prattNorm();}

    MULTIARCH inline VectorType eta()   const 
    {return Base::_ul * ( Scalar(1.) / Base::_ul.norm());}

    MULTIARCH inline Scalar     kappa() const 
    {return Scalar(2.) * Base::_uq / Base::prattNorm();}
    
    MULTIARCH inline Scalar     tau_normalized()   const {return tau()/Base::_tmax;}
    MULTIARCH inline VectorType eta_normalized()   const {return eta();}
    MULTIARCH inline Scalar     kappa_normalized() const {return kappa()*Base::_tmax;}
    
    MULTIARCH inline Scalar     fitness() const {
        return Scalar(1.) - Base::_ul.squaredNorm() - Scalar(4.) * Base::_uc * Base::_uq;}

  }; //class GLSParam

} //namespace Grenaille

#endif
