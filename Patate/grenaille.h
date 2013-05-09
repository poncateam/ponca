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

  \brief Grenaille: the simplest way to fit and analyze point-sets efficiently.

  This Patate is based on top of Eigen, and use advanced C++ template 
  programming to define generic fitting procedures and analysis.

  In this documentation we present an \ref grenaille_overview_sec of 
  the Grenaille programming concept, as well as a simple 
  \ref grenaille_howto_sec. For a complete example, please go to the
  \ref cpp/grenaille_basic_cpu.cpp "basic_cpu" example.

  \section grenaille_overview_sec Overview

  In Grenaille, you have access to Fitting Kernels and Extensions you can 
  combine using the \ref Basket template class to generate a specialized 
  fitting or analysis kernel:
  \image html grenaille/concept.svg

  \subsection grenaille_overview_basket_sec Basket
  \ref Basket "Basket <class Point, class Weight, class Fit, Ext-set ...>" 
  is an helper template class that you will specialize to define 
  the structure of the targeted kernel. To do so, please consider each of its 
  template parameters:
    - <b>Point</b>: defines the type of data the kernel will be applied on (e.g. number
  of dimensions, attributes),
    - <b>Weight</b>: defines the weighting function that will be applied to collect
  neighbors for the fit,
    - <b>Fit</b>: the fitting kernel
    - <b>Ext-set</b>: set of extensions that will be applied over the fit to
  obtain the desired behavior.

  In practice, you will have to define a Point class (you could use the example below
  as basis), a weighting function (most of the time you could use those provided in 
  Grenaille), choose your fitting kernel and the associated extensions. The interface
  of your class Point is defined by the requirements of the weighting, fitting and 
  extension classes.

  An example of Basket, to fit spheres on oriented point-sets and compute the GLS 
  descriptor:
  \code
  typedef Basket<MyPoint,WeightFunc,OrientedSphereFit, GLSParam> Fit;
  \endcode

  \section grenaille_howto_sec How to

  \subsection grenaille_howto_point_sec Point

  You must implemente the following interface to write your Point class:
    - Define a public enum value called Dim, to define the number of dimensions:
  \code
enum {Dim = 3};
  \endcode
    - Define a public Scalar type:
  \code
typedef float Scalar; //use floating precision
  \endcode
    - Define a public Vector type. It is strongly recommended to use an Eigen type:
  \code
typedef Eigen::Matrix<Scalar, Dim, 1> VectorType;
  \endcode
    - Provide an implicit constructor,
    - Provide at least a position vector, with r/w accessors:
  \code
  MULTIARCH const VectorType& pos()const;
  MULTIARCH       VectorType& pos();
  \endcode
  
  Example of a minimal class:
  \code
class MyPoint{
public:
  enum {Dim = 3};
  typedef float Scalar;
  typedef Eigen::Matrix<Scalar, Dim, 1> VectorType;

  MULTIARCH inline MyPoint(const VectorType &pos    = VectorType::Zero())
    : _pos(pos) {}
    
  MULTIARCH inline const VectorType& pos()    const { return _pos; }  
  MULTIARCH inline VectorType& pos()    { return _pos; }  

private:
  VectorType _pos;
};
  \endcode

  The other requirements defined by the \ref Basket specialization are listed 
  in the documentation of each used class (Weigthing function, fit or 
  extension). The following sections present at least these requirements for 
  the common cases.

  \subsection grenaille_howto_fit_sec Fitting kernel
  \warning Documentation not written yet

  \subsection grenaille_howto_weight_sec Weighting function

  There is two components combined to weight samples:
    - a 1D kernel, defined in the interval \f$ \left[0 \; : \; 1 \right]\f$. 
  Actually their is two kernels included in the library: 
  \ref ConstantWeightKernel and \ref SmoothWeightKernel. Please go to the 
  related documentation pages for more details.
    - a specialized function that will use this kernel with respect to a given
  point structure. It must derive from \ref BaseWeightFunc. Actually only 
  \ref DistWeightFunc is provided, to compute a weight with respect to the
  distance between samples for a given scale \f$ t \f$ specified at runtime.

  An example of Weighting function with the previous MyPoint class:
  \code
typedef DistWeightFunc<MyPoint,SmoothWeightKernel<Scalar> > WeightFunc; 
  \endcode

  \subsection grenaille_howto_ext_sec Extensions
  \warning Documentation not written yet



 */
namespace Grenaille
{ 
} // End namespace Grenaille

#endif //_PATATE_GRENAILLE_
