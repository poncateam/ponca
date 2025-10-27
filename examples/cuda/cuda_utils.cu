#pragma once

/*!
 * \brief Extract a vector from a flattened array of vectors
 *
 * \tparam DataPoint The DataPoint type containing the VectorType and the number of dimensions.
 * \param idx The id of the vector that needs to be extracted (the index in the flattened array corresponds to : idx * number of dimensions).
 * \param flattenedVectorArray The flattened vector array, that is of size : total number of vector * Dimension.
 * \return The vector that was extracted from the flattened vector array.
 */
template <typename DataPoint>
PONCA_MULTIARCH typename DataPoint::VectorType extractVectorFromFlattenedArray(
    const int idx,
    const typename DataPoint::Scalar * const flattenedVectorArray
) {
    using VectorType = typename DataPoint::VectorType;
    const int singleDimIndex = idx * DataPoint::Dim;
    VectorType v;

    for (int d = 0; d < DataPoint::Dim; ++d) {
        v.row(d) << flattenedVectorArray[singleDimIndex + d];
    }

    return v;
}

/*! \brief Kernel that binds each point buffer to the interlaced buffer of positions and normals in the Device.
 *
 * \tparam DataPoint The DataPoint type.
 * \tparam Scalar A scalar type (float, double...)
 * \param points The point cloud (needs to be bind to the interlacedArray on the device).
 * \param interlacedArray The interlaced buffer of positions and normals.
 * \param nbPoints The total number of points in the point cloud.
 */
template<typename DataPoint, typename Scalar>
__global__ void bindPointsKernel(
    DataPoint* points,
    Scalar* interlacedArray,
    const int nbPoints
) {
    // Get global index
    const unsigned int i = threadIdx.x + blockIdx.x * blockDim.x;

    // Skip when not in the point cloud
    if (i >= nbPoints) return;

    points[i].bind(interlacedArray);
}
