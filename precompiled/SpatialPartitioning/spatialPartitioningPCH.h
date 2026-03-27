/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include <Ponca/src/Common/pointTypes.h>
#include <Ponca/SpatialPartitioning>

template <typename Scalar>
using Point = Ponca::PointPositionNormal<Scalar, 3>;

extern template class Ponca::KdTreeBase<Ponca::KdTreeDefaultTraits<Point<float>>>;
extern template class Ponca::KdTreeBase<Ponca::KdTreeDefaultTraits<Point<double>>>;
extern template class Ponca::KdTreeBase<Ponca::KdTreeDefaultTraits<Point<long double>>>;

extern template class Ponca::KnnGraphBase<Ponca::KnnGraphDefaultTraits<Point<float>>>;
extern template class Ponca::KnnGraphBase<Ponca::KnnGraphDefaultTraits<Point<double>>>;
extern template class Ponca::KnnGraphBase<Ponca::KnnGraphDefaultTraits<Point<long double>>>;
