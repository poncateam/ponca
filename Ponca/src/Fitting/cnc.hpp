/*
This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once
#include <iostream>
#include <ostream>

namespace Ponca
{
    namespace internal {
        /*!
            \internal
            \brief Class to generate or generate a random integer from a presetted boundary
            \note Calling the () operator on this objet after the initialization of its boundary generates a random integer
        */
        class BoundedIntRange {
        public:
            const int _nMin;
            const int _nMax;
            explicit BoundedIntRange( const int nMax, const int nMin = 0 ) : _nMin(nMin), _nMax(nMax) { }

            void verifyBounds(const int n) const {
                if (_nMin > n || n > _nMax)
                    throw std::runtime_error(
                        "Index values must be in range :"
                        + std::to_string(_nMin) + " <= i <= " + std::to_string(_nMax)
                        + " But got result : " + std::to_string(n));
            }

            /// \brief Can be overwridden to something else in children class
            [[nodiscard]] int get(const int n) const {
                verifyBounds(n);
                return n;
            }
            /// \internal
            /// \brief Returns a random integer in bounds of : [ _nMin, _nMax ]
            [[nodiscard]] int random() const {
                // random operator
                const int r = Eigen::internal::random<int>(_nMin, _nMax);
                verifyBounds(r);
                return r;
            }
            class Iterator {
                int _current;
            public:
                explicit Iterator(int start) : _current(start) {}

                int operator*() const { return _current; }
                Iterator& operator++() { ++_current; return *this; }

                bool operator!=(const Iterator& other) const {
                    return _current != other._current;
                }
            };

            [[nodiscard]] Iterator begin() const { return Iterator(_nMin); }
            [[nodiscard]] Iterator end() const { return Iterator(_nMax); }
        };

        /*!
            \internal
            \brief Getting a random element from an STL-like container.
            Stores the container to then pick an element from it
            \note Calling the () operator on this objet after the initialization of its boundary picks a random element from the container
            \inherit GetRandomIndex
        */
        template<typename Container>
        class ElementSampler : public BoundedIntRange {
        private:
            Container& _elements;
        public:
            ElementSampler(Container& elements, const int nMax, const int nMin = 0) :
                BoundedIntRange(nMax, nMin), _elements(elements) { }

            [[nodiscard]] int get(const int i) const {
                verifyBounds(i);
                return _elements[i];
            }
            /// \internal
            /// \brief Returns a random elements from the container in the index range of : [ _nMin, _nMax ]
            /// \note Overloads the () operator to return an element picked from the container with the random value, instead of a random integer
            [[nodiscard]] int random() const {
                // random operator
                const int r = Eigen::internal::random<int>(_nMin, _nMax);
                verifyBounds(r);
                return _elements[r]; // Returns the element of the STL-like array
            }
            class Iterator {
                const ElementSampler& _parent;
                int _current;
            public:
                Iterator(const ElementSampler& parent, int start)
                    : _parent(parent), _current(start) {}

                auto operator*() const { return _parent._elements[_current]; }
                Iterator& operator++() { ++_current; return *this; }

                bool operator!=(const Iterator& other) const {
                    return _current != other._current;
                }
            };

            Iterator begin() const { return Iterator(*this, _nMin); }
            Iterator end() const { return Iterator(*this, _nMax); }
        };


    /// Generates the triangle used by the CNC Fit depending on the method
    template <TriangleGenerationMethod Method, typename P>
    struct TriangleGenerator {
        template <typename PointContainer, typename IndexGetter>
        static bool generate(
            const PointContainer& points,
            const IndexGetter& getIndex,
            const P& evalPoint,
            std::vector<internal::Triangle<P>>& triangles)
        {
            static_assert(true, "Triangle generation method not implemented!");
            return false;
        }
    };

    /// Generates the triangle used by the CNC Fit depending on the method (UniformGeneration)
    template <typename P>
    struct TriangleGenerator<TriangleGenerationMethod::UniformGeneration, P> {
        template <typename PointContainer, typename IndexGetter>
        static int generate(
            const PointContainer& points,
            const IndexGetter& getIndex,
            const P& /*evalPoint*/,
            std::vector<internal::Triangle<P>>& triangles)
        {
            int maxtriangles {100};
            int nb_vt = 0; // Number of valid generated triangles

            for (int i = 0; i < maxtriangles; ++i) {
                // Randomly select triangles
                int i1 = getIndex.random();
                int i2 = getIndex.random();
                int i3 = getIndex.random();
                if (i1 == i2 || i1 == i3 || i2 == i3) continue;

                std::array <typename P::VectorType, 3> positions  = {
                    points[i1].pos(),
                    points[i2].pos(),
                    points[i3].pos()
                };
                std::array <typename P::VectorType, 3> normals = {
                    points[i1].normal(),
                    points[i2].normal(),
                    points[i3].normal()
                };

                triangles.push_back(internal::Triangle<P>(positions, normals));
                nb_vt++;
            }
            return nb_vt;
        }
    };

    /// Generates the triangle used by the CNC Fit depending on the method (HexagramGeneration)
    template <typename P>
    struct TriangleGenerator<TriangleGenerationMethod::HexagramGeneration, P> {
        using VectorType = typename P::VectorType;
        using Scalar = typename P::Scalar;
        // Hexagram
        static std::array< Scalar    ,    6 > _distance2;
        static std::array< VectorType,    6 > _targets;

        template <typename PointContainer, typename IndexGetter>
        static int generate(
            const PointContainer& points,
            const IndexGetter& getIndex,
            const P& evalPoint,
            std::vector<internal::Triangle<P>>& triangles)
        {
            Scalar avgnormals  = Scalar(0.5);
            // BIN
            VectorType c = evalPoint.pos();
            VectorType n = evalPoint.normal();
            VectorType a;
            a.setZero();

            int iSource = -1;
            Scalar avgd = Scalar(0);

            for ( int index : getIndex ) {
                avgd += ( points[ index ].pos() - c ).norm();
                a    += points[ index ].normal();
                // if avgd == 0 then it is the evalPoint
                if ( iSource == -1 && points[ index ].pos() == c  ) {
                    iSource = index;
                }
            }

            a /= a.norm();
            n = ( Scalar(1) - avgnormals ) * n + avgnormals * a;
            n /= n.norm();
            avgd /= getIndex._nMax;

            const int m = ( std::abs( n[0] ) > std::abs ( n[1] ))
                    ? ( ( std::abs( n[0] ) ) > std::abs( n[2] ) ? 0 : 2 )
                    : ( ( std::abs( n[1] ) ) > std::abs( n[2] ) ? 1 : 2 );
            const VectorType e =
                ( m == 0 ) ? VectorType( Scalar(0), Scalar(1), Scalar(0) ) :
                ( m == 1 ) ? VectorType( Scalar(0), Scalar(0), Scalar(1) ) :
                VectorType( Scalar(1), Scalar(0), Scalar(0) );

            VectorType u = n.cross( e );
            VectorType v = n.cross( u );
            u /= u.norm();
            v /= v.norm();

            std::array<int, 6> indices = {iSource, iSource, iSource, iSource, iSource, iSource};

            for ( int i = 0 ; i < 6 ; i++ ){
                _distance2 [ i ] = avgd * avgd;
                _targets   [ i ] = avgd * ( u * std::cos( i * M_PI / 3.0 ) + v * std::cos( i * M_PI / 3.0 ) );
            }

            for ( int index : getIndex ) {
                VectorType p = points[ index ].pos();
                if ( p == c ) {
                    std::cout << "p == c" << std::endl;
                    continue;
                }

                const VectorType d = p - c;
                for ( int j = 0 ; j < 6 ; j++ ){
                    const Scalar d2 = ( d - _targets[ j ]).squaredNorm();
                    if ( d2 < _distance2[ j ] ){
                        indices[ j ] = index;
                        _distance2[ j ] = d2;
                    }
                }
            }
            std::array <VectorType, 3> t1_points  = {points[indices[0]].pos(), points[indices[2]].pos(), points[indices[4]].pos()};
            std::array <VectorType, 3> t1_normals = {points[indices[0]].normal(), points[indices[2]].normal(), points[indices[4]].normal()};

            std::array <VectorType, 3> t2_points  = {points[indices[1]].pos(), points[indices[3]].pos(), points[indices[5]].pos()};
            std::array <VectorType, 3> t2_normals = {points[indices[1]].normal(), points[indices[3]].normal(), points[indices[5]].normal()};

            triangles.push_back(internal::Triangle<P>(t1_points, t1_normals));
            triangles.push_back(internal::Triangle<P>(t2_points, t2_normals));

            return 2;
        }
    };

    }
    template < class P, class W, TriangleGenerationMethod M>
    template <typename PointContainer>
    FIT_RESULT CNC<P, W, M>::compute( const PointContainer& points ) {
        if (M != TriangleGenerationMethod::UniformGeneration) {
            throw std::runtime_error("Used TriangleGenerationMethod UniformGeneration but forgot to precise the eval point");
        }
        auto p = points[0]; // Dummy point
        // Random index from the size of the point container
        internal::BoundedIntRange indexGetter( points.size()-1 );
        internal::TriangleGenerator<M, P>::generate( points, indexGetter, p, _triangles);

        return finalize();
    }

    template < class P, class W, TriangleGenerationMethod M>
    template <typename PointContainer>
    FIT_RESULT CNC<P, W, M>::compute( const PointContainer& points, const P& evalPoint ) {
        // Random index from the size of the point container
        internal::BoundedIntRange indexGetter( points.size()-1 );
        internal::TriangleGenerator<M, P>::generate( points, indexGetter, evalPoint, _triangles);

        return finalize();
    }

    template < class P, class W, TriangleGenerationMethod M>
    template <typename IndexContainer, typename PointContainer>
        FIT_RESULT CNC<P, W, M>::computeWithIds( const IndexContainer& ids, const PointContainer& points ) {
        if (M != TriangleGenerationMethod::UniformGeneration) {
            throw std::runtime_error("Used TriangleGenerationMethod UniformGeneration but forgot to precise the eval point");
        }
        auto p = points[0]; // Dummy point
        // Getting a random index from an index container
        internal::ElementSampler indexGetter( ids, ids.size()-1);
        internal::TriangleGenerator<M, P>::generate( points, indexGetter, p, _triangles);

        return finalize();
    }

    template < class P, class W, TriangleGenerationMethod M>
    template <typename IndexContainer, typename PointContainer>
    FIT_RESULT CNC<P, W, M>::computeWithIds( const IndexContainer& ids, const PointContainer& points, const P& evalPoint ) {
        // Getting a random index from an index container
        internal::ElementSampler indexGetter( ids, ids.size()-1 );
        internal::TriangleGenerator<M, P>::generate( points, indexGetter, evalPoint, _triangles);

        return finalize();
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
