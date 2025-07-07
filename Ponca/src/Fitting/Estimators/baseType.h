#pragma once

#include <Eigen/Core>

namespace Ponca {
    /// Map a block to a Ponca point
    template <typename _Scalar>
    struct BlockPointAdapter {
    public:
        enum {Dim = 3};
        using Scalar     = _Scalar;
        using VectorType = Eigen::Matrix<Scalar, Dim, 1>;
        using MatrixType = Eigen::Matrix<Scalar, Dim, Dim>;

        /// \brief Map a vector as ponca Point
        PONCA_MULTIARCH inline BlockPointAdapter(VectorType v, VectorType n)
                : m_pos (v), m_nor (n) {}

        PONCA_MULTIARCH inline VectorType pos()    const { return m_pos; }
        PONCA_MULTIARCH inline VectorType normal() const { return m_nor; }

        PONCA_MULTIARCH inline BlockPointAdapter& operator=(const BlockPointAdapter& other) = default;

    private:
        VectorType m_pos;
        VectorType m_nor;
    };

    using Scalar = float;
    using PPAdapter = BlockPointAdapter<Scalar>;
    using VectorType = Eigen::Matrix<Scalar, 3, 1>;
    using SampleMatrixType = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
    using SampleVectorType = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
}