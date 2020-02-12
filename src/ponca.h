/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

// First inclue Eigen Core
#include <Eigen/Core>

// Include common stuff
#include "Ponca/defines.h"

// Include Ponca Core components
#include "Ponca/enums.h"
#include "Ponca/basket.h"

#include "Ponca/weightKernel.h"
#include "Ponca/weightFunc.h"

#include "Ponca/plane.h"
#include "Ponca/meanPlaneFit.h"
#include "Ponca/covariancePlaneFit.h"
#include "Ponca/mongePatch.h"

#include "Ponca/sphereFit.h"
#include "Ponca/orientedSphereFit.h"
#include "Ponca/mlsSphereFitDer.h"
#include "Ponca/curvature.h"
#include "Ponca/gls.h"

#include "Ponca/curvatureEstimation.h"

// not supported on cuda
#ifndef __CUDACC__
# include "Ponca/unorientedSphereFit.h"
#endif


// Include Ponca Algorithms
