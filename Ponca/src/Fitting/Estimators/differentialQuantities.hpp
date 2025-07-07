#pragma once

#include <utility>
#include <vector>
#include <Eigen/Dense>

namespace Ponca {
    template <typename _Scalar>
    struct Quantity {
        using Scalar = _Scalar;
        using VectorType = Eigen::Matrix<Scalar, 3, 1>;

    public :
        Scalar k1 = 0;
        Scalar k2 = 0;
        Scalar mean = 0;
        Scalar gauss = 0;
        Scalar neighbor_count = 0;
        VectorType d1 = VectorType(1, 0, 0);
        VectorType d2 = VectorType(1, 0, 0);
        VectorType normal = VectorType(1, 0, 0);
        VectorType projection = VectorType::Zero();

        int non_stable = 0;

        Quantity() = default;

        Quantity(const Quantity& q) {
            this->k1 = q.k1;
            this->k2 = q.k2;
            this->mean = q.mean;
            this->gauss = q.gauss;
            this->d1 = q.d1;
            this->d2 = q.d2;

            this->neighbor_count = q.neighbor_count;

            this->projection = q.projection;
            this->normal = q.normal;

            this->non_stable = q.non_stable;
        }
    };
}