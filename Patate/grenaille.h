/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/. 
*/


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

  A fitting kernel define at least four methods:
  \code
// init
MULTIARCH void setWeightFunc (const WeightFunc& w);
MULTIARCH void init (const VectorType& evalPos);

// processing
MULTIARCH void addNeighbor(const Point &nei);
MULTIARCH void finalize   ();
  \endcode

  Please go to the documentation of each fitting primitive to know their 
  requirements. For example, \ref OrientedSphereFit requires to define
  a normal vector and the associated accessors in the Point class:
  \code
class MyPoint{
public:
  enum {Dim = 3};
  typedef float Scalar;
  typedef Eigen::Matrix<Scalar, Dim, 1> VectorType;

  MULTIARCH inline MyPoint(const VectorType &pos    = VectorType::Zero(), 
		 const VectorType& normal = VectorType::Zero())
    : _pos(pos), _normal(normal) {}
    
  MULTIARCH inline const VectorType& pos()    const { return _pos; }  
  MULTIARCH inline const VectorType& normal() const { return _normal; }

  MULTIARCH inline VectorType& pos()    { return _pos; }  
  MULTIARCH inline VectorType& normal() { return _normal; }
  };

private:
  VectorType _pos, _normal;
};
  \endcode


  The fitting kernel should also provide a macro, indicating the type of the fit. For example,
  \ref OrientedSphereFit defines
  \code
PROVIDES_ALGEBRAIC_SPHERE
  \endcode


  \subsection grenaille_howto_compat_sec Basket compatibility check

  The provided macros are checked by all the extensions to ensure the compatibility 
  during the compilation and throw human-readable error message. For example:
  \code
Patate/Grenaille/Core/gls.h:91:7: error: ‘PROVIDES_GLS_PARAMETRIZATION’ is not a member of ...
  \endcode
  indicates that the current \ref requires an extension defining the GLS 
  parametrization (in practice: \ref GLSParam).

  This error is thrown when trying to compile the type
\code
typedef Basket<MyPoint,WeightFunc,OrientedSphereFit, OrientedSphereScaleDer, GLSDer> Fit;
\endcode

  Please go the the documentation of the extensions you plan to use
  to know their requirements and the macros they provide.

  \subsection grenaille_howto_ext_sec Extensions
  \warning This documentation is not yet ready

  \subsection grenaille_howto_weight_sec Weighting function

  There is two components combined to weight samples:
    - a 1D kernel, defined in the interval \f$ \left[0 \; : \; 1 \right]\f$,
  that inherits BaseWeightKernel. 
  Actually their is two kernels included in the library: 
  \ref ConstantWeightKernel and \ref SmoothWeightKernel. Please go to the 
  related documentation pages for more details.
    - a specialized function that will use this kernel with respect to a given
  point structure. It must derive from \ref BaseWeightFunc. Actually only 
  \ref DistWeightFunc is provided, to compute a weight with respect to the
  distance between samples for a given scale \f$ t \f$ specified at runtime.

  An example of Weighting function with the MyPoint class:
  \code
typedef DistWeightFunc<MyPoint,SmoothWeightKernel<Scalar> > WeightFunc; 
  \endcode



 */
namespace Grenaille
{ 
} // End namespace Grenaille

#endif //_PATATE_GRENAILLE_
