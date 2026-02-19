/**
Copyright (c) 2022
 Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr)
 Laboratory of Mathematics (CNRS, UMR 5807), University of Savoie, France,

All rights reserved.

*/

#pragma once

#include "defines.h"
#include "cncFormulaEigen.h"
#include "weightFunc.h"
#include "weightKernel.h"


namespace Ponca
{

namespace internal {
/*!
 * Stores the three points and normals of the triangles and provides access to Corrected Normal Current formula
 *
 * \tparam DataPoint Type of input points.
 *
 * \see CNCEigen
 */
template < class DataPoint >
struct Triangle {
public:
    typedef typename DataPoint::Scalar Scalar;
    typedef typename DataPoint::VectorType VectorType;
    typedef typename DataPoint::MatrixType MatrixType;
protected:
    std::array < VectorType, 3 > m_points;
    std::array < VectorType, 3 > m_normals;
public:
    Triangle(DataPoint pointA, DataPoint pointB, DataPoint pointC) {
        m_points = {
            pointA.pos(),
            pointB.pos(),
            pointC.pos()
        };
        m_normals = {
            pointA.normal(),
            pointB.normal(),
            pointC.normal()
        };
    }

    Triangle(const std::array < VectorType, 3 >& points, const std::array < VectorType, 3 >& normals) {
        m_points = points;
        m_normals = normals;
    }

    /*! \brief Get the position of the point at the given index.
     *
     * @param index Index of one of the three vertices of the triangle (between 0 and 2)
     * @return The position of vertex
     */
    PONCA_MULTIARCH [[nodiscard]] VectorType& getPos(const int index) {
        return m_points[index];
    }

    PONCA_MULTIARCH [[nodiscard]] bool operator==(const Triangle& other) const {
        return (m_points[0] == other.m_points[0])
			&& (m_points[1] == other.m_points[1])
			&& (m_points[2] == other.m_points[2]);
    }

    PONCA_MULTIARCH [[nodiscard]] bool operator!=(const Triangle& other) const {
        return !((*this) == other);
    }

#define DEFINE_CNC_FUNC(CNC_FUNC, RETURN_TYPE)                                      \
    template<bool differentOrder = false>                                           \
    inline RETURN_TYPE CNC_FUNC () {                                                \
        return CNCEigen<DataPoint>::CNC_FUNC(                                       \
            m_points[0] ,  m_points[2-differentOrder],  m_points[1+differentOrder], \
            m_normals[0], m_normals[2-differentOrder], m_normals[1+differentOrder]  \
        );                                                                          \
    }

	DEFINE_CNC_FUNC(mu0InterpolatedU , Scalar)
	DEFINE_CNC_FUNC(mu1InterpolatedU , Scalar)
	DEFINE_CNC_FUNC(mu2InterpolatedU , Scalar)
	DEFINE_CNC_FUNC(muXYInterpolatedU, MatrixType)
};
} // namespace internal


/*!
 * \breif Generation method of the triangles for the Corrected Normal Current formula
 */
enum TriangleGenerationMethod {
    UniformGeneration, HexagramGeneration, IndependentGeneration, AvgHexagramGeneration
};

/*!
 * \brief Corrected Normal Current Fit type.
 *
 * This fitting method generates triangles from a set a points cloud and use a statistical formula to compute :
 * - The principal curvatures values and directions
 * - The mean and gaussian curvatures
 *
 * \see PROVIDES_PRINCIPAL_CURVATURES
*/
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
    using NeighborFilter = NeighborFilterStoreNormal<DataPoint, NoWeightFunc<DataPoint>>;
protected:
	// Basis
    NeighborFilter m_nFilter;

    // Triangles used for the computation
    int m_nb_vt {0}; // Number of valid triangles
    std::vector <internal::Triangle < DataPoint > > m_triangles;

    // Results of the fit
    Scalar m_A   {0}; // Area
    Scalar m_H   {0}; // Mean Curvatures
    Scalar m_G   {0}; // Gaussian Curvatures
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

    /*! \brief Represent the current state of the fit
     * (finalize function update the state)
     */
    FIT_RESULT m_eCurrentState {UNDEFINED};
public:
    PONCA_FITTING_DECLARE_FINALIZE

    //! \brief Set the scalar field values to 0 and reset the isNormalized() status
    PONCA_MULTIARCH inline void init() {
        m_eCurrentState = UNDEFINED;
        m_A  = Scalar(0);
        m_H  = Scalar(0);
        m_G  = Scalar(0);

        m_k1 = Scalar(0);
        m_k2 = Scalar(0);

        m_triangles.clear();
    }

    /*!
     * \brief Compute function for STL-like containers.
     * \tparam PointContainer An STL-like container storing the points
     */
    template <typename PointContainer>
    PONCA_MULTIARCH inline FIT_RESULT compute( const PointContainer& points);

    /*!
     * \brief Compute function that iterates over a subset of sampled points from an STL-Like container.
     * \tparam IndexRange An STL-like container storing the indices of the neighbors
     * \tparam PointContainer An STL-like container storing the points
     */
    template <typename IndexRange, typename PointContainer>
    PONCA_MULTIARCH inline FIT_RESULT computeWithIds( const IndexRange& ids, const PointContainer& points );

    /*!
     * \brief Get the number of triangles that were generated with the compute method.
     * \return Returns the number of triangles
     * \see CNC::compute
     */
    PONCA_MULTIARCH [[nodiscard]] inline size_t getNumTriangles() const {
        return static_cast<size_t>(m_nb_vt);
    }

    PONCA_FITTING_APIDOC_SETWFUNC
    PONCA_MULTIARCH inline void setNeighborFilter (const NeighborFilter& _nFilter) {
        m_nFilter  = _nFilter;
    }

    /*!
     * \brief Returns the triangles
     *
     * \return A pointer to the array of triangle that was generated during the CNC Fit
     */
    PONCA_MULTIARCH [[nodiscard]] std::vector< internal::Triangle<DataPoint> >& getTriangles() {
        return m_triangles;
    }

    //! \brief Comparison operator
    PONCA_MULTIARCH [[nodiscard]] bool operator==(const CNC& other) const {
        // We use the matrix to compare the fitting results
        return (m_eCurrentState == other.m_eCurrentState)
            && (kMean() == other.kMean())
            && (kmin() == other.kmin())
            && (kmax() == other.kmax())
            && (kminDirection() == other.kminDirection())
            && (kmaxDirection() == other.kmaxDirection())
            && (GaussianCurvature() == other.GaussianCurvature())
            && (m_T11 == other.m_T11) && (m_T12 == other.m_T12) && (m_T13 == other.m_T13)
            && (m_T22 == other.m_T22) && (m_T23 == other.m_T23)
            && (m_T33 == other.m_T33);
    }

    //! \brief Comparison operator, convenience function
    PONCA_MULTIARCH [[nodiscard]] bool operator!=(const CNC& other) const {
        // We use the matrix to compare the fitting results
        return !(this == &other);
    }

    //! \brief Approximate operator
    PONCA_MULTIARCH [[nodiscard]] bool isApprox(const CNC& other, const Scalar& epsilon = Eigen::NumTraits<Scalar>::dummy_precision()) const {
        PONCA_MULTIARCH_STD_MATH(abs);

        return (m_eCurrentState == other.m_eCurrentState)
            && (std::abs(kMean()  - other.kMean())  < epsilon)
            && (std::abs(GaussianCurvature() - other.GaussianCurvature()) < epsilon)
            && (std::abs(kmin() - other.kmin()) < epsilon)
            && (std::abs(kmax() - other.kmax()) < epsilon);
    }

    //! \brief Is the fitted primitive ready to use (finalize has been called and the result is stable)
    PONCA_MULTIARCH [[nodiscard]] inline bool isStable() const { return m_eCurrentState == STABLE; }

    //! \brief Returns an estimate of the minimal principal curvature value
    PONCA_MULTIARCH [[nodiscard]] inline Scalar kmin() const { return m_k1; }

    //! \brief Returns an estimate of the maximal principal curvature value
    PONCA_MULTIARCH [[nodiscard]] inline Scalar kmax() const { return m_k2; }

    //! \brief Returns an estimate of the minimal principal curvature direction
    PONCA_MULTIARCH [[nodiscard]] inline VectorType kminDirection() const { return m_v1; }

    //! \brief Returns an estimate of the maximal principal curvature direction
    PONCA_MULTIARCH [[nodiscard]] inline VectorType kmaxDirection() const { return m_v2; }

    //! \brief Returns an estimate of the mean curvature
    PONCA_MULTIARCH [[nodiscard]] inline Scalar kMean() const { return m_H; }

    //! \brief Returns an estimate of the Gaussian curvature
    PONCA_MULTIARCH [[nodiscard]] inline Scalar GaussianCurvature() const { return m_G; }
}; //class CNC

} // namespace Ponca

#include "cnc.hpp"
