#ifndef _GRENAILLE_ORIENTED_SPHERE_FIT_
#define _GRENAILLE_ORIENTED_SPHERE_FIT_


#ifndef __CUDA_ARCH__
#include <math.h>
#endif



namespace Grenaille
{
  template < class DataPoint, class _WFunctor, typename T = void >
  class OrientedSphereFit {
  protected:
    enum
      {
        PROVIDES_ALGEBRAIC_SPHERE
      };

  public:
    typedef typename DataPoint::Scalar     Scalar;
    typedef typename DataPoint::VectorType VectorType;
    typedef _WFunctor             WFunctor;
    
  protected:
    // Evaluation position (needed for centered basis)
    VectorType _p;
    
    // computation data
    VectorType _sumN, _sumP;
    Scalar _sumDotPN, _sumDotPP, _sumW;
    
    WFunctor _w;

    // results
  public:
    Scalar _uc, _uq;
    VectorType _ul;
    
  public:
    MULTIARCH inline OrientedSphereFit(){ }

    // getters
    MULTIARCH inline const VectorType& evalPos () const { return _p; }
    MULTIARCH inline       VectorType& evalPos ()       { return _p; }
    
    // init
    MULTIARCH inline void setWeightFunc (const WFunctor& w) { _w  = w; }
    MULTIARCH inline void init (const VectorType& evalPos);

    // processing
    //\todo Change impl
    MULTIARCH inline void addNeighbor(const DataPoint &nei);
    MULTIARCH inline void finalize   ();
    
    //! compute the squared Pratt norm of the implicit scalar field
    MULTIARCH inline Scalar prattNorm2() const {
      return _ul.squaredNorm() - Scalar(4.) * _uc*_uq;
    }
    
    //! compute the Pratt norm of the implicit scalar field
    MULTIARCH inline Scalar prattNorm() const {
      MULTIARCH_STD_MATH(sqrt);
      return sqrt(prattNorm2());
    }
    
    //! Project a point on the sphere
    MULTIARCH inline VectorType project (VectorType q);
    
    //    MULTIARCH VectorType gradient(VectorType q, bool normalize = true);

  }; //class OrientedSphereFit

#include "orientedSphereFit.hpp"

} //namespace Grenaille


#endif
