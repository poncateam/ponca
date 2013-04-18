#ifndef _GRENAILLE_WEIGHT_FUNC_
#define _GRENAILLE_WEIGHT_FUNC_

namespace Grenaille{
  template <class DataPoint, typename Derived >
  class BaseWeightFunc {

  public:        
    typedef typename DataPoint::Scalar Scalar;
    typedef typename DataPoint::VectorType VectorType;
      
    MULTIARCH inline Scalar w(const VectorType& relativeQuery, 
			      const DataPoint&  attributes) const
    { return _der().w(relativeQuery, attributes); }   
       
    MULTIARCH inline VectorType spacedw(const VectorType& relativeQuery, 
				    const DataPoint&  attributes) const
    { return _der().spacedw(relativeQuery, attributes); }   
       
    MULTIARCH inline Scalar scaledw(const VectorType& relativeQuery, 
				    const DataPoint&  attributes) const
    { return _der().scaledw(relativeQuery, attributes); }

    MULTIARCH inline Scalar evalScale() const
    { return _der().evalScale(); }

  protected:
    MULTIARCH inline const Derived& _der() { return &static_cast<Derived*>(this); }    
  };// class BaseWeightFunc


  /*!
    \warning Assumes that the evaluation scale t is strictly positive
   */
  template <class DataPoint, class WeightKernel>
  class DistWeightFunc: public BaseWeightFunc<DataPoint, DistWeightFunc<DataPoint, WeightKernel> >{
  public:
    typedef typename DataPoint::Scalar Scalar;
    typedef typename DataPoint::VectorType VectorType;
    
    MULTIARCH inline DistWeightFunc(const Scalar& t = Scalar(1.)): _t(t) {}

    /*!
      Compute a weight using the norm of the query \f$ \mathbf{q} \f$ 
      (expressed in centered basis)
     */
    MULTIARCH inline Scalar w(const VectorType& q, 
			      const DataPoint&  /*attributes*/) const;
    
    
    /*!
      First order derivative in space (for each dimension \f$\mathsf{x})\f$:
      
      \f$ \frac{\delta \frac{\mathbf{q}_\mathsf{x}}{t}}{\delta \mathsf{x}} 
      \nabla w(\frac{\mathbf{q}_\mathsf{x}}{t}) 
      = \frac{ \nabla{w(\frac{\mathbf{q}_\mathsf{x}}{t})}}{t}  \f$
      
      where \f$ \mathbf{q}_\mathsf{x} \f$ represent the query coordinate in 
      the spatial dimension \f$ \mathsf{x}\f$ expressed in centered basis.
    */
    MULTIARCH inline VectorType spacedw(const VectorType& q, 
			     const DataPoint&  /*attributes*/) const;
       
    
    /*!
      First order derivative in scale t:
      
      \f$ \frac{\delta \frac{\mathbf{q}}{t}}{\delta t} 
      \nabla w(\frac{\mathbf{q}}{t}) 
      = - \frac{\mathbf{q}}{t^2} \nabla{w(\frac{\mathbf{q}}{t})} \f$
      
      where \f$ \mathbf{q} \f$ represent the query coordinate expressed in 
      centered basis.
    */
    MULTIARCH inline Scalar scaledw(const VectorType& q, 
			     const DataPoint&  /*attributes*/) const;

    MULTIARCH inline Scalar evalScale() const { return _t; }

  protected:
    Scalar       _t;
    WeightKernel _wk;
    

  };// class DistWeightFunc

#include "weightFunc.hpp"

}// namespace Grenaille


#endif // _GRENAILLE_WEIGHT_FUNC_
