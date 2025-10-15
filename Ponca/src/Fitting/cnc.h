/**
Copyright (c) 2022
 Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr)
 Laboratory of Mathematics (CNRS, UMR 5807), University of Savoie, France,

All rights reserved.

*/

#pragma once

#include "defines.h"
#include "cncFormulaEigen.h"


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
class CNC : ComputeObject<CNC<P, _method>> {
protected:
    enum
    {
        PROVIDES_PRINCIPAL_CURVATURES
    };
public:
    using DataPoint = P;
    using MatrixType = typename DataPoint::MatrixType;
    using Scalar     = typename DataPoint::Scalar;
    using VectorType = typename DataPoint::VectorType;
    typedef Eigen::VectorXd  DenseVector;
    typedef Eigen::MatrixXd  DenseMatrix;
protected:
	// Basis
    VectorType m_evalPointNormal = VectorType::Zero();
    VectorType m_evalPointPos = VectorType::Zero();

    // Triangles used for the computation
    int m_nb_vt {0}; // Number of valid triangles
    std::vector <internal::Triangle < DataPoint > > m_triangles;

    // Results of the fit
    Scalar m_A {0}; // Area
    Scalar m_H {0}; // Mean Curvatures
    Scalar m_G {0}; // Gaussian Curvatures
    Scalar m_T11 {0}; // T11
    Scalar m_T12 {0}; // T12
    Scalar m_T13 {0}; // T13
    Scalar m_T22 {0}; // T22
    Scalar m_T23 {0}; // T23
    Scalar m_T33 {0}; // T33

    Scalar m_k1 {0};
    Scalar m_k2 {0};

    VectorType m_v1;
    VectorType m_v2;

public:
    PONCA_FITTING_DECLARE_FINALIZE

    /*! \brief Set the scalar field values to 0 and reset the isNormalized() status
    */
    PONCA_MULTIARCH inline void init() {
        m_k1 = Scalar(0);
        m_k2 = Scalar(0);
    }

    /*! \brief Compute function for STL-like containers */
    /*! Add neighbors stored in a PointContainer and call finalize at the end.*/
    template <typename PointContainer>
    PONCA_MULTIARCH inline FIT_RESULT compute( const PointContainer& points);

    /*! \brief Compute function to iterate over a subset of samples in a PointContainer  */
    /*! Add neighbors stored in a PointContainer and sampled using indices stored in ids.*/
    /*! \tparam IndexRange An STL-like container storing the indices of the neighbors */
    /*! \tparam PointContainer An STL-like container storing the points */
    template <typename IndexRange, typename PointContainer>
    PONCA_MULTIARCH inline FIT_RESULT computeWithIds( const IndexRange& ids, const PointContainer& points );

    /*! \brief Returns the number of fitted triangles  */
    PONCA_MULTIARCH inline int getNumTriangles() const {
        return m_nb_vt;
    }

	void setEvalPoint(const DataPoint& evalPoint) {
        m_evalPointNormal = evalPoint.normal();
        m_evalPointPos    = evalPoint.pos();
    }
    void setEvalPoint(const VectorType& evalPointNormal, const VectorType& evalPointPos) {
        m_evalPointNormal = evalPointNormal;
        m_evalPointPos    = evalPointPos;
    }

    bool operator==(const CNC& other) const {
        // We use the matrix to compare the fitting results
        return (m_T11 == other.m_T11) && (m_T12 == other.m_T12) && (m_T13 == other.m_T13) && (m_T22 == other.m_T22) && (m_T23 == other.m_T23) && (m_T33 == other.m_T33);
    }
    bool operator!=(const CNC& other) const {
        // We use the matrix to compare the fitting results
        return !(this == &other);
    }

    template<typename Fit>
    bool isApprox(const Fit& other, const Scalar& epsilon = Eigen::NumTraits<Scalar>::dummy_precision()) const {
        // Simply compare the kMean and kGauss results
        return std::abs(kMean()  - other.kMean())  < epsilon
            && std::abs(GaussianCurvature() - other.GaussianCurvature()) < epsilon;
    }

    PONCA_MULTIARCH inline Scalar kmin() const { return m_k1; }

    PONCA_MULTIARCH inline Scalar kmax() const { return m_k2; }

    PONCA_MULTIARCH inline VectorType kminDirection() const { return m_v1; }

    PONCA_MULTIARCH inline VectorType kmaxDirection() const { return m_v2; }

    PONCA_MULTIARCH inline Scalar kMean() const { return m_H; }

    PONCA_MULTIARCH inline Scalar GaussianCurvature() const { return m_G; }
}; //class CNC

} // namespace Ponca

#include "cnc.hpp"