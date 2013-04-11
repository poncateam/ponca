#ifndef _GRENAILLE_WEIGHT_KERNEL_
#define _GRENAILLE_WEIGHT_KERNEL_

/*!
  \file weightKernel Define 1D weight kernel functors
 */


namespace Grenaille{

  template <typename _Scalar>
  class ConstantWeightKernel{
  public:
    typedef _Scalar Scalar;
    
    // Init 
    MULTIARCH inline ConstantWeightKernel(const Scalar& value = Scalar(1.)) : _y(value){}
    MULTIARCH inline void setValue(const Scalar& value){ _y = value; } 

    // Functor
    MULTIARCH inline Scalar f  (const Scalar&) const { return _y; }
    MULTIARCH inline Scalar df (const Scalar&) const { return Scalar(0.); }
    MULTIARCH inline Scalar ddf(const Scalar&) const { return Scalar(0.); }

  private: 
    Scalar _y;
  };// class ConstantWeightKernel

  
  /*!
    \todo Add way to specify a degree parameter (in this class or another one)
   */
  template <typename _Scalar>
  class SmoothWeightKernel{
  public:
    typedef _Scalar Scalar;

    // Functor
    MULTIARCH inline Scalar f  (const Scalar& x) const { Scalar v = x*x - Scalar(1.); return v*v; }
    MULTIARCH inline Scalar df (const Scalar& x) const { return Scalar(4.)*x*(x*x-Scalar(1.)); }
    MULTIARCH inline Scalar ddf(const Scalar& x) const { return Scalar(12.)*x*x - Scalar(4.); }
  };//class SmoothWeightKernel

}// namespace Grenaille


#endif //_GRENAILLE_WEIGHT_KERNEL_
