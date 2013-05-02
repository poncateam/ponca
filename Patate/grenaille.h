#ifndef _PATATE_GRENAILLE_


// Include common stuff
#include "common/defines.h"

// Include Grenaille Core components
#include "Grenaille/Core/basket.h"

#include "Grenaille/Core/weightKernel.h"
#include "Grenaille/Core/weightFunc.h"

#include "Grenaille/Core/rawSphereFit.h"
#include "Grenaille/Core/orientedSphereFit.h"
#include "Grenaille/Core/unorientedSphereFit.h"
#include "Grenaille/Core/gls.h"


// Include Grenaille Algorithms



/*!

  \brief Documentation

  Toto

  You must define a Point class:
  \code
  // Define our working data structure
class MyPoint{
public:
  enum {Dim = 3};
  typedef double Scalar;
  typedef Eigen::Matrix<Scalar, Dim, 1> VectorType;

  MULTIARCH inline MyPoint(const VectorType &pos    = VectorType::Zero(), 
		 const VectorType& normal = VectorType::Zero())
    : _pos(pos), _normal(normal) {}
    
  MULTIARCH inline const VectorType& pos()    const { return _pos; }  
  MULTIARCH inline const VectorType& normal() const { return _normal; }

  MULTIARCH inline VectorType& pos()    { return _pos; }  
  MULTIARCH inline VectorType& normal() { return _normal; }

  static inline MyPoint Random() { 
    return MyPoint (VectorType::Random(), VectorType::Random());
  };

private:
  VectorType _pos, _normal;
};
  \endcode

  and define
  \code
typedef DistWeightFunc<MyPoint,SmoothWeightKernel<Scalar> > WeightFunc; 
typedef Basket<MyPoint,WeightFunc,OrientedSphereFit, GLSParam, OrientedSphereScaleDer, GLSDer, GLSGeomVar> Fit;
  \endcode
 */
namespace Grenaille
{ } // End namespace Grenaille




/*!
  \example Grenaille/basic_cpu.cpp
  This is an example of how to use Grenaille to compute 
  the GLS Geometric variation on random data.

 */



#endif //_PATATE_GRENAILLE_
