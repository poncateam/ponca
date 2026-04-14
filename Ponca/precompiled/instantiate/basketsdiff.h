/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

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
// Wrap the type within the _PONCA_BASKET_DIFF_X call:
//  - Name of the basket
//  - Description
//  - Actual basket class
//
// The name and description might be used to generate mappings
// or bindings to these class
//
// The _PONCA_BASKET_DIFF_X macro will be defined elsewhere to generate
// the desired code

#if _D == 3 && defined(_Normal)
_PONCA_BASKET_DIFF_X(
    "Diff", "",
    BasketDiff<Basket<_P, _NF, OrientedSphereFit, GLSParam>, _DT, OrientedSphereDer, GLSDer, CurvatureEstimatorDer,
               NormalDerivativeWeingartenEstimator, WeingartenCurvatureEstimatorDer>)
#endif
