/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/
#if defined(_PONCA_INSTANTIATE_SCALESPACEDER) || defined(_PONCA_INSTANTIATE_ALL)
#    define _DT FitScaleSpaceDer
#    include "basketsdiff.h"
#    undef _DT
#endif
#if defined(_PONCA_INSTANTIATE_SPACEDER) || defined(_PONCA_INSTANTIATE_ALL)
#    define _DT FitSpaceDer
#    include "basketsdiff.h"
#    undef _DT
#endif

// Define the list of basket in this file:
// - _S: Scalar type
// - _D: Dimension
// - _P: Point type
// - _NF: Neighbor filter
//
// You may filter by dimension, but not by "type" as the
// preprocessor does not allows this. You may however define
// a custom macros and then check for this definition. This
// is what is done with '_Normal'. Please, do not forget to
// add them to the following list:
//
// Type information macros:
// - _Normal: Defined when point provides normals
//
// Wrap the type within the _PONCA_BASKET_X call:
//  - Name of the basket
//  - Description
//  - Actual basket class
//
// The name and description might be used to generate mappings
// or bindings to these class
//
// The _PONCA_BASKET_X macro will be defined elsewhere to generate
// the desired code
#if _D == 3
_PONCA_BASKET_X("MongePatchFit", "", Basket<_P, _NF, MongePatchQuadraticFit>)
_PONCA_BASKET_X("MongePatchRestrictedFit", "Provides curvature estimation without spatial differentiation",
                Basket<_P, _NF, MongePatchRestrictedQuadraticFit>)
#endif

#if _D == 2 || _D == 3
_PONCA_BASKET_X("CovariancePlaneFit", "", Basket<_P, _NF, CovariancePlaneFit>)
_PONCA_BASKET_X("SphereFit", "", Basket<_P, _NF, SphereFit, GLSParam>)
#    if defined(_Normal)
_PONCA_BASKET_X("MeanPlaneFit", "", Basket<_P, _NF, MeanPlaneFit>)
_PONCA_BASKET_X("OrientedSphereFit", "AKA Algebraic Point Set Surfaces (APSS)",
                Basket<_P, _NF, OrientedSphereFit, GLSParam>)
_PONCA_BASKET_X("UnorientedSphereFit", "AKA", Basket<_P, _NF, UnorientedSphereFit, GLSParam>)
#    endif
#endif
