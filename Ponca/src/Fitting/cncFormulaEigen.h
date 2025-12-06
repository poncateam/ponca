/**
 Copyright (c) 2020
 Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr)
 Laboratory of Mathematics (CNRS, UMR 5127), University of Savoie, France,
 David Coeurjolly (\c david.coeurjolly@liris.cnrs.fr)
 LIRIS (CNRS, UMR 5205), CNRS, France

All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
* Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.

* Neither the name of the <organization> nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT
HOLDER> BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

/**
 * @file cncFormulaEigen.h
 * @brief Adapted version of CorrectedNormalCurrentFormulaEigen.h for Ponca.
 *
 * @details
 * ## Original file :
 * https://github.com/JacquesOlivierLachaud/PointCloudCurvCNC/blob/main/CorrectedNormalCurrentFormulaEigen.h
 *
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr)
 * Laboratory of Mathematics (CNRS, UMR 5127), University of Savoie, France
 * @author David Coeurjolly (\c david.coeurjolly@liris.cnrs.fr)
 * LIRIS (CNRS, UMR 5205), CNRS, France
 *
 * @date 2022/06/05
 *
 * ## Modifications made :
 *
 * - Renamed file name to cncFormulaEigen.h
 * - Moved CNCEigen struct to Ponca::internal namespace
 * - Adapted CNCEigen to use the Ponca DataPoint Types
 *
 * @author Florian Auberval (\c florian.auberval@irit.fr)
 *
 * @date 2025/11/19
 */

#pragma once

#include <Eigen/Dense>

namespace Ponca::internal {

/**
 * @brief This class contains some stand-alone CorrectedNormalCurrent formulas for triangles,
 * using eigen as linear algebra backend.
 */
template<typename DataPoint>
struct CNCEigen {
	using MatrixType = typename DataPoint::MatrixType;
	using Scalar     = typename DataPoint::Scalar;
	using VectorType = typename DataPoint::VectorType;
	/// Small constant used to approximate zero.
	static constexpr Scalar epsilon = Eigen::NumTraits<Scalar>::epsilon();

	/// Represents a triangle on a sphere of radius one.
	struct SphericalTriangle {
		///Spherical point data type
		static bool isDegenerate(const VectorType& a, const VectorType& b, const VectorType& c) {
			Scalar d[3] = {(a - b).norm(), (a - c).norm(), (b - c).norm()};
			// Checks that the spherical triangle is small or thin.
			if ((d[0] < epsilon) || (d[1] < epsilon) || (d[2] < epsilon))
				return true;
			// Checks that the spherical triangle is flat.
			size_t m = 0;
			if (d[1] > d[m]) m = 1;
			if (d[2] > d[m]) m = 2;
			return (fabs(d[m] - d[(m + 1) % 3] - d[(m + 2) % 3]) < epsilon);
		}

		/// @return the polar triangle associated with this triangle.
		static void polarTriangle(
			const VectorType& a, const VectorType& b, const VectorType& c,
			VectorType& Ap     , VectorType& Bp     , VectorType& Cp
		) {
			Ap = b.cross(c);
			Bp = c.cross(a);
			Cp = a.cross(b);
			// Reorient points.
			if (Ap.dot(a) < 0.0) Ap = -Ap;
			if (Bp.dot(b) < 0.0) Bp = -Bp;
			if (Cp.dot(c) < 0.0) Cp = -Cp;
		}

		/// Returns the interior angles of the spherical triangle ABC.
		/// @param[out] alpha the interior angle at vertex A.
		/// @param[out] beta  the interior angle at vertex B.
		/// @param[out] gamma the interior angle at vertex C.
		static void interiorAngles(
			const VectorType& a, const VectorType& b, const VectorType& c,
		    Scalar& alpha      , Scalar& beta       , Scalar& gamma
		) {
			VectorType Ta, Tb, Tc;
			polarTriangle(a, b, c, Ta, Tb, Tc);
			Ta /= Ta.norm();
			Tb /= Tb.norm();
			Tc /= Tc.norm();
			if (Ta == VectorType::Zero() || Tb == VectorType::Zero() || Tc == VectorType::Zero())
				alpha = beta = gamma = Scalar(0.0);
			else
			{
				Scalar ca = std::max(Scalar(-1.0), std::min(Scalar(1.0), Tb.dot(Tc)));
				Scalar cb = std::max(Scalar(-1.0), std::min(Scalar(1.0), Tc.dot(Ta)));
				Scalar cc = std::max(Scalar(-1.0), std::min(Scalar(1.0), Ta.dot(Tb)));
				alpha = acos(ca);
				beta = acos(cb);
				gamma = acos(cc);
			}
		}

		/// @return the (unsigned) area of the spherical triangle (below 2pi).
		static Scalar area(const VectorType& a, const VectorType& b, const VectorType& c) {
			Scalar alpha, beta, gamma;
			if (isDegenerate(a, b, c)) return 0.0;
			interiorAngles(a, b, c, alpha, beta, gamma);
			return ((fabs(alpha) < epsilon)
				       || (fabs(beta) < epsilon)
				       || (fabs(gamma) < epsilon))
				       ? Scalar(0.0)
				       : 2.0 * M_PI - alpha - beta - gamma;
		}

		/// @return the (signed) area of the spherical triangle (below 2pi).
		static Scalar algebraicArea(const VectorType& a, const VectorType& b, const VectorType& c) {
			Scalar S = area(a, b, c);
			VectorType M = a + b + c;
			VectorType X = (b - a).cross(c - a);
			if (M.template lpNorm<1>() <= epsilon || X.template lpNorm<1>() <= epsilon) return 0.0;
			return M.dot(X) < 0.0 ? -S : S;
		}
	};


	// ---------------------- Main functions ----------------
public:
	/// @name Functions for computing measures
	/// @{

	/// Computes mu0 measure (area) of triangle abc given an interpolated
	/// corrected normal vector \a ua, \a \ub, \a uc.
	/// @param a any point
	/// @param b any point
	/// @param c any point
	/// @param ua the corrected normal vector at point a
	/// @param ub the corrected normal vector at point b
	/// @param uc the corrected normal vector at point c
	/// @param unit_u when 'true' considers that interpolated
	/// corrected normals should be made unitary, otherwise
	/// interpolated corrected normals may have smaller norms.
	/// @return the mu0-measure of triangle abc, i.e. its area.
	static Scalar mu0InterpolatedU(
		const VectorType& a,
	    const VectorType& b,
	    const VectorType& c,
		const VectorType& ua,
		const VectorType& ub,
		const VectorType& uc,
		bool unit_u = false
	) {
		// MU0=1/2*det( uM, B-A, C-A )
		//    =  1/2 < ( (u_A + u_B + u_C)/3.0 ) | (AB x AC ) >
		VectorType uM = (ua + ub + uc) / 3.0;
		if (unit_u) {
			auto uM_norm = uM.norm();
			uM = uM_norm == 0.0 ? uM : uM / uM_norm;
		}
		return 0.5 * ((b - a).cross(c - a)).dot(uM);
	}


	/// Computes mu1 measure (mean curvature) of triangle abc given an interpolated
	/// corrected normal vector \a ua, \a \ub, \a uc.
	/// @param a any point
	/// @param b any point
	/// @param c any point
	/// @param ua the corrected normal vector at point a
	/// @param ub the corrected normal vector at point b
	/// @param uc the corrected normal vector at point c
	/// @param unit_u when 'true' considers that interpolated
	/// corrected normals should be made unitary, otherwise
	/// interpolated corrected normals may have smaller norms.
	/// @return the mu1-measure of triangle abc, i.e. its mean curvature.
	static Scalar mu1InterpolatedU(
		const VectorType& a,
		const VectorType& b,
		const VectorType& c,
		const VectorType& ua,
		const VectorType& ub,
		const VectorType& uc,
		bool unit_u = false
    ) {
		// MU1=1/2( | uM u_C-u_B A | + | uM u_A-u_C B | + | uM u_B-u_A C |
		VectorType uM = (ua + ub + uc) / 3.0;
		if (unit_u) uM /= uM.norm();
		return 0.25 * (uM.cross(uc - ub).dot(a)
			+ uM.cross(ua - uc).dot(b)
			+ uM.cross(ub - ua).dot(c));
	}


	/// Computes mu2 measure (Gaussian curvature) of triangle abc given an interpolated
	/// corrected normal vector \a ua, \a \ub, \a uc.
	/// @param a any point
	/// @param b any point
	/// @param c any point
	/// @param ua the corrected normal vector at point a
	/// @param ub the corrected normal vector at point b
	/// @param uc the corrected normal vector at point c
	/// @param unit_u when 'true' considers that interpolated
	/// corrected normals should be made unitary, otherwise
	/// interpolated corrected normals may have smaller norms.
	/// @return the mu2-measure of triangle abc, i.e. its Gaussian curvature.
	static Scalar mu2InterpolatedU(
		const VectorType& a,
		const VectorType& b,
		const VectorType& c,
		const VectorType& ua,
		const VectorType& ub,
		const VectorType& uc,
		bool unit_u = false
	) {
		// Using non unitary interpolated normals give
		// MU2=1/2*det( uA, uB, uC )
		// When normals are unitary, it is the area of a spherical triangle.
		if (unit_u)
			return SphericalTriangle::algebraicArea(ua, ub, uc);
		else
			return 0.5 * (ua.cross(ub).dot(uc));
	}


	/// Computes muXY measure (anisotropic curvature) of triangle abc given an interpolated
	/// corrected normal vector \a ua, \a \ub, \a uc.
	/// @param a any point
	/// @param b any point
	/// @param c any point
	/// @param ua the corrected normal vector at point a
	/// @param ub the corrected normal vector at point b
	/// @param uc the corrected normal vector at point c
	/// @return the muXY-measure of triangle abc, i.e. its anisotropic curvature.
	static MatrixType muXYInterpolatedU(
		const VectorType& a,
		const VectorType& b,
		const VectorType& c,
		const VectorType& ua,
		const VectorType& ub,
		const VectorType& uc,
		bool unit_u = false
	) {
		MatrixType T = MatrixType::Zero();
		VectorType uM = (ua + ub + uc) / 3.0;
		if (unit_u) uM /= uM.norm();
		const VectorType uac = uc - ua;
		const VectorType uab = ub - ua;
		const VectorType ab = b - a;
		const VectorType ac = c - a;
		for (size_t i = 0; i < 3; ++i) {
			VectorType X = VectorType::Zero();
			X(i) = 1.0;
			for (size_t j = 0; j < 3; ++j) {
				// Since RealVector Y = RealVector::base( j, 1.0 );
				// < Y | uac > = uac[ j ]
				const Scalar tij =
					0.5 * uM.dot(uac[j] * X.cross(ab)
						- uab[j] * X.cross(ac));
				T(i, j) = tij;
			}
		}
		return T;
	}

	/// @}

	// ---------------------- Helper functions ----------------
public:
	/// @name Helper functions
	/// @{

	/// Computing principal curvatures k1 and k2 from tensor
	/// @param tensor The muXY integrated tensor
	/// @param area Area of the face
	/// @param N the normal vector
	/// @return a pair of principal directions.
	static std::pair<VectorType, VectorType> curvDirFromTensor(
		const MatrixType& tensor,
		const Scalar area,
		const VectorType& N
	) {
		auto Mt = tensor.transpose();
		auto M = tensor;
		M += Mt;
		M *= 0.5;
		const Scalar coef_N = 1000.0 * area;
		// Adding 1000 area n x n to anisotropic measure
		for (int j = 0; j < 3; j++)
			for (int k = 0; k < 3; k++)
				M(j, k) += coef_N * N[j] * N[k];

		Eigen::SelfAdjointEigenSolver<MatrixType> eigensolver(M);
		if (eigensolver.info() != Eigen::Success) abort();

		//SelfAdjointEigenSolver returns sorted eigenvalues, no
		//need to reorder the eigenvectors.
		assert(eigensolver.eigenvalues()(0) <= eigensolver.eigenvalues()(1));
		assert(eigensolver.eigenvalues()(1) <= eigensolver.eigenvalues()(2));
		VectorType v1 = eigensolver.eigenvectors().col(1);
		VectorType v2 = eigensolver.eigenvectors().col(0);
		return std::pair<VectorType, VectorType>(v1, v2);
	}

	/// Computing principal curvatures k1 and k2 from tensor
	/// @param tensor The muXY integrated tensor
	/// @param area Area of the face
	/// @param N the normal vector
	/// @return a pair of principal directions.
	static std::tuple<Scalar, Scalar, VectorType, VectorType> curvaturesFromTensor(
		const MatrixType& tensor,
		const Scalar area,
		const VectorType& N
	) {
		auto Mt = tensor.transpose();
		auto M = tensor;
		M += Mt;
		M *= 0.5;
		const Scalar coef_N = 1000.0 * area;
		// Adding 1000 area n x n to anisotropic measure
		for (int j = 0; j < 3; j++)
			for (int k = 0; k < 3; k++)
				M(j, k) += coef_N * N[j] * N[k];

		Eigen::SelfAdjointEigenSolver<MatrixType> eigensolver(M);
		if (eigensolver.info() == Eigen::Success) {
			// SelfAdjointEigenSolver returns sorted eigenvalues, no
			// need to reorder the eigenvectors.
			assert(eigensolver.eigenvalues()(0) <= eigensolver.eigenvalues()(1));
			assert(eigensolver.eigenvalues()(1) <= eigensolver.eigenvalues()(2));
			VectorType v1 = eigensolver.eigenvectors().col(0);
			VectorType v2 = eigensolver.eigenvectors().col(1);
			return std::make_tuple(-eigensolver.eigenvalues()(0),
			                       -eigensolver.eigenvalues()(1),
			                       v1, v2);
		} else {
			std::cerr << "Incorrect diagonalization for tensor " << M << std::endl;
			VectorType v1, v2;
			return std::make_tuple(Scalar(0.0), Scalar(0.0), v1, v2);
		}
	}

	/// @}
}; // struct CNCEigen

}
