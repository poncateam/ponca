/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#pragma once

#include "defines.h"
#include PONCA_MULTIARCH_INCLUDE_STD(cmath)
#include "cncFormulaEigen.h"

#define DEFINE_CNC_FUNC(CNC_FUNC)                                             \
	template<bool differentOrder = false>                                     \
    inline Scalar CNC_FUNC () {                                               \
        return CNCEigen::CNC_FUNC(                                            \
			 points[0],  points[2-differentOrder],  points[1+differentOrder], \
			normals[0], normals[2-differentOrder], normals[1+differentOrder]  \
		);                                                                    \
    }


namespace Ponca
{

namespace internal {

// triangle storing indices of points
template < class DataPoint >
struct Triangle {
    typedef typename DataPoint::Scalar Scalar;
    typedef typename DataPoint::VectorType VectorType;
    typedef typename DataPoint::MatrixType MatrixType;

    std::array < VectorType, 3 > points;
    std::array < VectorType, 3 > normals;
    // Maybe need to store the normal too

    Triangle(const std::array < VectorType, 3 >& points, const std::array < VectorType, 3 >& normals) {
        this->points = points;
        this->normals = normals;
    }

    bool operator==(const Triangle& other) const {
        return (points[0] == other.points[0])
			&& (points[1] == other.points[1])
			&& (points[2] == other.points[2]);
    }

    bool operator!=(const Triangle& other) const {
        return !((*this) == other);
    }

	DEFINE_CNC_FUNC(mu0InterpolatedU)
	DEFINE_CNC_FUNC(mu1InterpolatedU)
	DEFINE_CNC_FUNC(mu2InterpolatedU)
	DEFINE_CNC_FUNC(muXYInterpolatedU)
};

} // namespace internal

/*!
    \brief CNC generation of triangles from a set of points
*/

enum class TriangleGenerationMethod {
    UniformGeneration
};

template < class P, class WeightFunc, TriangleGenerationMethod _method = TriangleGenerationMethod::UniformGeneration>
class CNC : BasketBase<P, WeightFunc> {
public:
    using DataPoint = P;
    using MatrixType = typename DataPoint::MatrixType;
    using Scalar     = typename DataPoint::Scalar;
    using VectorType = typename DataPoint::VectorType;
    typedef Eigen::VectorXd  DenseVector;
    typedef Eigen::MatrixXd  DenseMatrix;

protected:
    static const auto& randomInt = Eigen::internal::random<int>;
	// Basis
	VectorType _evalPointNormal   {VectorType::Zero()};

    //! \brief protected variables
    std::array < Scalar, 6 > _cos;
    std::array < Scalar, 6 > _sin;

    int _nb_vt {0}; // Number of valid triangles
    std::vector <internal::Triangle < DataPoint > > _triangles;
    Scalar _A {0}; // Area
    Scalar _H {0}; // Mean Curvatures
    Scalar _G {0}; // Gaussian Curvatures
    Scalar _T11 {0}; // T11
    Scalar _T12 {0}; // T12
    Scalar _T13 {0}; // T13
    Scalar _T22 {0}; // T22
    Scalar _T23 {0}; // T23
    Scalar _T33 {0}; // T33

    Scalar k1 {0};
    Scalar k2 {0};

    VectorType v1;
    VectorType v2;

    // Hexagram
    std::array< Scalar    ,    6 > _distance2;
    std::array< VectorType,    6 > _targets;

// results
public:
    /*!< \brief Parameters of the triangles */
    int _maxtriangles {100};
    Scalar _avgnormals {Scalar(0.5)};

    PONCA_FITTING_DECLARE_FINALIZE

    /*! \brief Set the scalar field values to 0 and reset the isNormalized() status

    */

    PONCA_MULTIARCH inline void init() {
        k1 = Scalar(0);
        k2 = Scalar(0);

        v1 = VectorType::Zero();
        v2 = VectorType::Zero();

        // Instantiate the parameters
        _maxtriangles = 100;
        for ( int j = 0; j < 6; j++ )
        {
            const Scalar a = j * M_PI / 3.0;
            _cos[ j ] = std::cos( a );
            _sin[ j ] = std::sin( a );
        }
    }

    template <typename PointContainer>
    PONCA_MULTIARCH inline FIT_RESULT compute( const PointContainer& points );

    template <typename PointContainer>
    PONCA_MULTIARCH inline std::enable_if_t<_method == TriangleGenerationMethod::UniformGeneration, bool> generateTriangles( const PointContainer& points );

    PONCA_MULTIARCH inline int getNumTriangles() const {
        return _nb_vt;
    }

    PONCA_MULTIARCH inline void getTriangles( std::vector<std::array<Scalar, 3>>& triangles ) {

        for (int i = 0; i < _triangles.size(); i++) {
            std::array <Scalar, 3> point0 = {_triangles[i].points[0][0], _triangles[i].points[0][1], _triangles[i].points[0][2]};
            std::array <Scalar, 3> point1 = {_triangles[i].points[1][0], _triangles[i].points[1][1], _triangles[i].points[1][2]};
            std::array <Scalar, 3> point2 = {_triangles[i].points[2][0], _triangles[i].points[2][1], _triangles[i].points[2][2]};
            triangles.push_back(point0);
            triangles.push_back(point1);
            triangles.push_back(point2);        
        }

    }

	void setEvalPointNormal(const VectorType& evalPointNormal) {
		_evalPointNormal = evalPointNormal;
	}

    bool operator==(const CNC& other) const
    {
        // We use the matrix to compare the fitting results
        return (_T11 == other._T11) && (_T12 == other._T12) && (_T13 == other._T13) && (_T22 == other._T22) && (_T23 == other._T23) && (_T33 == other._T33);
    }
    bool operator!=(const CNC& other) const
    {
        // We use the matrix to compare the fitting results
        return !(this == &other);
    }

    PONCA_MULTIARCH inline Scalar kmin() { return k1; }

    PONCA_MULTIARCH inline Scalar kmax() { return k2; }

    PONCA_MULTIARCH inline VectorType kminDirection() { return v1; }

    PONCA_MULTIARCH inline VectorType kmaxDirection() { return v2; }

    PONCA_MULTIARCH inline Scalar kMean() { return _H; }

    PONCA_MULTIARCH inline Scalar kgauss() { return _G; }
}; //class CNC

} // namespace Ponca

#include "cnc.hpp"