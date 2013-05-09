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

  \brief Grenaille: the simplest way to fit and analyze point-set data

  This Patate is based on top of Eigen, and use advanced C++ template 
  programming to define generic fitting procedures and analysis.

  In this documentation we present an \ref grenaille_overview_sec of 
  the Grenaille programming concept, as well as a simple 
  \ref grenaille_howto_sec. For a complete example, please go to the
  \ref cpp/grenaille_basic_cpu.cpp "basic_cpu" example.

  \section grenaille_overview_sec Overview

  \section grenaille_howto_sec How to

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
{ 
} // End namespace Grenaille

#endif //_PATATE_GRENAILLE_
