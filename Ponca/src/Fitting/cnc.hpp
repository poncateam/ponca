/*
This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

namespace Ponca {

namespace internal
{
    /*!
        \internal
        \brief Class to generate a random integer from a presetted boundary
        \note Calling the () operator on this objet after the initialization of its boundary generates a random integer
    */
    class GetRandomInt {
    public:
        const int _nMin;
        const int _nMax;

        explicit GetRandomInt( const int nMax, const int nMin = 0 ) : _nMax(nMax), _nMin(nMin) { }

        /// \internal
        /// \brief Returns a random integer in bounds of : [ _nMin, _nMax ]
        int operator()() const {
            // random operator
            const int r = Eigen::internal::random<int>(_nMin, _nMax);
            if (_nMin > r || r > _nMax)
                throw std::runtime_error(
                    "Random index values must be in range :"
                    + std::to_string(_nMin) + " <= i <= " + std::to_string(_nMax)
                    + " But got result : " + std::to_string(r));
            return r;
        }
    };

    /*!
        \internal
        \brief Getting a random element from an STL-like container.
        Stores the container to then pick an element from it
        \note Calling the () operator on this objet after the initialization of its boundary picks a random element from the container
        \inherit GetRandomIndex
    */
    template<typename Container>
    class GetRandomElementFromContainer : GetRandomInt {
    private:
        Container& _elements;
    public:
        GetRandomElementFromContainer(Container& elements, const int nMax, const int nMin = 0) :
            GetRandomInt(nMax, nMin), _elements(elements) { }

        /// \internal
        /// \brief Returns a random elements from the container in the index range of : [ _nMin, _nMax ]
        /// \note Overloads the () operator to return an element picked from the container with the random value, instead of a random integer
        auto operator()() const {
            // random operator
            const int r = Eigen::internal::random<int>(_nMin, _nMax);
            if (_nMin > r || r > _nMax)
                throw std::runtime_error(
                    "Random index values must be in range :"
                    + std::to_string(_nMin) + " <= i <= " + std::to_string(_nMax)
                    + " But got result : " + std::to_string(r));
            return _elements[r];
        }
    };
}


template < class P, class W, TriangleGenerationMethod M>
template <typename PointContainer>
FIT_RESULT CNC<P, W, M>::compute( const PointContainer& points ) {
    // Random index from the size of the point container
    internal::GetRandomInt rdmIndex( points.size()-1 );
    generateTriangles( points, rdmIndex );

	return finalize();
}
template < class P, class W, TriangleGenerationMethod M>
template <typename IndexContainer, typename PointContainer>
FIT_RESULT CNC<P, W, M>::computeWithIds( const IndexContainer& ids, const PointContainer& points ) {
    // Getting a random index from an index container
    internal::GetRandomElementFromContainer rdmIndex( ids, ids.size()-1 );
    generateTriangles( points, rdmIndex );

    return finalize();
}

/// Generates the triangle used by the CNC Fit depending on the method (UniformGeneration)
template <class P, class W, TriangleGenerationMethod M>
template <typename PointContainer, typename RandomIndexGetter>
std::enable_if_t<M == TriangleGenerationMethod::UniformGeneration, bool>
CNC<P, W, M>::generateTriangles(
	const PointContainer& points,
    const RandomIndexGetter& rdmId
) {
    _nb_vt = 0; // Number of valid generated triangles

    for (int i = 0; i < _maxtriangles; ++i) {
        // Randomly select triangles
        int i1 = rdmId();
        int i2 = rdmId();
        int i3 = rdmId();
        if (i1 == i2 || i1 == i3 || i2 == i3) continue;

        std::array <VectorType, 3> positions  = {
			points[i1].pos(),
			points[i2].pos(),
			points[i3].pos()
		};
        std::array <VectorType, 3> normals = {
			points[i1].normal(),
			points[i2].normal(),
			points[i3].normal()
		};

        _triangles.push_back(internal::Triangle<P>(positions, normals));
        _nb_vt++;
    }
    return _nb_vt > 0;
}

template < class P, class W, TriangleGenerationMethod M>
FIT_RESULT CNC<P, W, M>::finalize( ) {
    _A = Scalar(0);
    _H = Scalar(0);
    _G = Scalar(0);

    MatrixType localT = MatrixType::Zero();

    for (int t = 0; t < _nb_vt; ++t) {

        // Simple estimation.
        Scalar tA = _triangles[t].mu0InterpolatedU();
        if (tA < -internal::CNCEigen<P>::epsilon) {
            _A     -= tA;
            _H     += _triangles[t].template mu1InterpolatedU<true>();
            _G     += _triangles[t].template mu2InterpolatedU<true>();
            localT += _triangles[t].template muXYInterpolatedU<true>();
        } else if (tA > internal::CNCEigen<P>::epsilon) {
            _A     += tA;
            _H     += _triangles[t].mu1InterpolatedU();
            _G     += _triangles[t].mu2InterpolatedU();
            localT += _triangles[t].muXYInterpolatedU();
        }

    } // end for t

    _T11 = localT(0,0);
    _T12 = 0.5 * (localT(0,1) + localT(1,0));
    _T13 = 0.5 * (localT(0,2) + localT(2,0));
    _T22 = localT(1,1);
    _T23 = 0.5 * (localT(1,2) + localT(2,1));
    _T33 = localT(2,2);

    MatrixType T;

    if (_A != Scalar(0)){
        T  << _T11, _T12, _T13,
              _T12, _T22, _T23,
              _T13, _T23, _T33;
        T  /= _A;
        _H /= _A;
        _G /= _A;
    } else {
        _H = Scalar(0);
        _G = Scalar(0);
    }

    std::tie (k2, k1, v2, v1) = internal::CNCEigen<P>::curvaturesFromTensor(T, 1.0, _evalPointNormal);

    return STABLE;
}
}