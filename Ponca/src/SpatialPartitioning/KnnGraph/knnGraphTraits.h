/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#pragma once

#include <Eigen/Geometry>

namespace Ponca {

/*!
 * \brief The default traits type used by the kd-tree.
 */
template <typename _DataPoint>
struct KnnGraphDefaultTraits
{
    /*!
     * \brief The type used to store point data.
     *
     * Must provide `Scalar` and `VectorType` typedefs.
     *
     * `VectorType` must provide a `squaredNorm()` function returning a `Scalar`, as well as a
     * `maxCoeff(int*)` function returning the dimension index of its largest scalar in its output
     * parameter (e.g. 0 for *x*, 1 for *y*, etc.).
     */
    using DataPoint = _DataPoint;

private:
    using Scalar     = typename DataPoint::Scalar;
    using VectorType = typename DataPoint::VectorType;

public:
    /*!
     * \brief The type used to calculate node bounding boxes.
     *
     * Must provide `min()`, `max()`, and `center()` functions, all returning a `VectorType`.
     */
    using AabbType = Eigen::AlignedBox<Scalar, DataPoint::Dim>;

    // Containers
    using IndexType      = int;
    using PointContainer = std::vector<DataPoint>;
    using IndexContainer = std::vector<IndexType>;
};
} // namespace Ponca
