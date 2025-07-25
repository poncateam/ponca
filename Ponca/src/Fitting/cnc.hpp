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
        \brief Getting a random index
    */
    class GetRandomIndex {
    public:
        int _indexMax {0};

        GetRandomIndex(const int indexMax) : _indexMax(indexMax) { }

        //! \brief Returns a random index in bounds of : [ 0, indexMax ]
        int operator()() const {
            // random operator
            int r = Eigen::internal::random<int>(0, _indexMax);
            if (0 > r || r > _indexMax)
                throw std::runtime_error(
                "Random index values must be in range :"
                    "0 < i < " + std::to_string(_indexMax) +
                    " But got result : " + std::to_string(r));
            return r;
        }
    };

    /*!
        \brief Getting a random element from an STL-like container
        \inherit GetRandomIndex
    */
    template<typename Container>
    class GetRandomElementFromContainer : GetRandomIndex {
    private:
        Container& _ids;
    public:
        GetRandomElementFromContainer(const int indexMax, Container& ids) :
            GetRandomIndex(indexMax), _ids(ids) { }

        //! \brief Returns a random index from the index container
        int operator()() const {
            // random operator
            int r = Eigen::internal::random<int>(0, _indexMax);
            if (0 > r || r > _indexMax)
                throw std::runtime_error(
                "Random index values must be in range :"
                    "0 < i < " + std::to_string(_indexMax) +
                    " But got result : " + std::to_string(r));
            return _ids[r];
        }
    };
}


template < class P, class W, TriangleGenerationMethod M>
template <typename PointContainer>
FIT_RESULT CNC<P, W, M>::compute( const PointContainer& points ) {
    // Random index from the size of the point container
    internal::GetRandomIndex rdmId( points.size()-1 );
    generateTriangles( points, rdmId );

	return finalize();
}
template < class P, class W, TriangleGenerationMethod M>
template <typename IndexContainer, typename PointContainer>
FIT_RESULT CNC<P, W, M>::computeWithIds( const IndexContainer& ids, const PointContainer& points ) {
    // Getting a random index from an index container
    internal::GetRandomElementFromContainer rdmId( ids.size()-1, ids );
    generateTriangles( points, rdmId );

    return finalize();
}

/// Generates the triangle used by the CNC Fit depending on the method (UniformGeneration)
template <class P, class W, TriangleGenerationMethod M>
template <typename PointContainer, typename IndexRandomGetter>
std::enable_if_t<M == TriangleGenerationMethod::UniformGeneration, bool>
CNC<P, W, M>::generateTriangles(
	const PointContainer& points,
    const IndexRandomGetter& rdmId
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