/**
Copyright (c) 2022
 Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr)
 Laboratory of Mathematics (CNRS, UMR 5807), University of Savoie, France,

All rights reserved.

*/

#pragma once

#include "defines.h"
#include PONCA_MULTIARCH_INCLUDE_STD(cmath)
#include "cncFormulaEigen.h"
#include <Ponca/src/Fitting/weightFunc.h>


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

    Triangle(DataPoint pointA, DataPoint pointB, DataPoint pointC) {
        points = {
            pointA.pos(),
            pointB.pos(),
            pointC.pos()
        };
        normals = {
            pointA.normal(),
            pointB.normal(),
            pointC.normal()
        };
    }
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

#define DEFINE_CNC_FUNC(CNC_FUNC, RETURN_TYPE)                                \
    template<bool differentOrder = false>                                     \
    inline RETURN_TYPE CNC_FUNC () {                                          \
        return CNCEigen<DataPoint>::CNC_FUNC(                                 \
            points[0] ,  points[2-differentOrder],  points[1+differentOrder], \
            normals[0], normals[2-differentOrder], normals[1+differentOrder]  \
        );                                                                    \
}

	DEFINE_CNC_FUNC(mu0InterpolatedU , Scalar)
	DEFINE_CNC_FUNC(mu1InterpolatedU , Scalar)
	DEFINE_CNC_FUNC(mu2InterpolatedU , Scalar)
	DEFINE_CNC_FUNC(muXYInterpolatedU, MatrixType)
};

} // namespace internal

/*!
    \brief CNC generation of triangles from a set of points
*/

enum TriangleGenerationMethod {
    UniformGeneration, HexagramGeneration, IndependentGeneration, AvgHexagramGeneration
};

template < class P, TriangleGenerationMethod _method = UniformGeneration>
class CNC : BasketBase<P, NoWeightFunc<P>> {
public:
    using DataPoint = P;
    using MatrixType = typename DataPoint::MatrixType;
    using Scalar     = typename DataPoint::Scalar;
    using VectorType = typename DataPoint::VectorType;
    typedef Eigen::VectorXd  DenseVector;
    typedef Eigen::MatrixXd  DenseMatrix;
protected:
	// Basis
    VectorType _evalPointNormal = VectorType::Zero();
    VectorType _evalPointPos = VectorType::Zero();

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

// results
public:
    PONCA_FITTING_DECLARE_FINALIZE

    /*! \brief Set the scalar field values to 0 and reset the isNormalized() status
    */
    PONCA_MULTIARCH inline void init() {
        k1 = Scalar(0);
        k2 = Scalar(0);
    }

    /*! \brief Compute function for STL-like containers */
    /*! Add neighbors stored in a PointContainer and call finalize at the end.*/
    template <typename PointContainer>
    PONCA_MULTIARCH inline FIT_RESULT compute( const PointContainer& points);

    /*! \brief Compute function to iterate over a subset of samples in a PointContainer  */
    /*! Add neighbors stored in a PointContainer and sampled using indices stored in ids.*/
    /*! \tparam IndexContainer An STL-like container storing the indices of the neighbors */
    /*! \tparam PointContainer An STL-like container storing the points */
    template <typename IndexContainer, typename PointContainer>
    PONCA_MULTIARCH inline FIT_RESULT computeWithIds( const IndexContainer& ids, const PointContainer& points );

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

	void setEvalPoint(const DataPoint& evalPoint) {
        _evalPointNormal = evalPoint.normal();
        _evalPointPos    = evalPoint.pos();
    }
    void setEvalPoint(const VectorType& evalPointNormal, const VectorType& evalPointPos) {
        _evalPointNormal = evalPointNormal;
        _evalPointPos    = evalPointPos;
    }

    bool operator==(const CNC& other) const {
        // We use the matrix to compare the fitting results
        return (_T11 == other._T11) && (_T12 == other._T12) && (_T13 == other._T13) && (_T22 == other._T22) && (_T23 == other._T23) && (_T33 == other._T33);
    }
    bool operator!=(const CNC& other) const {
        // We use the matrix to compare the fitting results
        return !(this == &other);
    }

    bool isApprox(const CNC& other, const Scalar& epsilon = Eigen::NumTraits<Scalar>::dummy_precision()) const {
        // Simply compare the kMean and kGauss results
        return std::abs(kMean()  - other.kMean())  < epsilon
            && std::abs(kGauss() - other.kGauss()) < epsilon;
    }

    PONCA_MULTIARCH inline Scalar kmin() const { return k1; }

    PONCA_MULTIARCH inline Scalar kmax() const { return k2; }

    PONCA_MULTIARCH inline VectorType kminDirection() const { return v1; }

    PONCA_MULTIARCH inline VectorType kmaxDirection() const { return v2; }

    PONCA_MULTIARCH inline Scalar kMean() const { return _H; }

    PONCA_MULTIARCH inline Scalar kGauss() const { return _G; }

    PONCA_MULTIARCH inline VectorType project(const VectorType& _q) const {return _q;}
    PONCA_MULTIARCH inline VectorType primitiveGradient() const {return VectorType::Zero();}
}; //class CNC

} // namespace Ponca

#include "cnc.hpp"