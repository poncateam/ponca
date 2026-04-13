/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "spatialPartitioningPCH.h"

template class Ponca::KdTreeBase<Ponca::KdTreeDefaultTraits<Point<float>>>;
template class Ponca::KdTreeBase<Ponca::KdTreeDefaultTraits<Point<double>>>;
template class Ponca::KdTreeBase<Ponca::KdTreeDefaultTraits<Point<long double>>>;

template class Ponca::KnnGraphBase<Ponca::KnnGraphDefaultTraits<Point<float>>>;
template class Ponca::KnnGraphBase<Ponca::KnnGraphDefaultTraits<Point<double>>>;
template class Ponca::KnnGraphBase<Ponca::KnnGraphDefaultTraits<Point<long double>>>;
