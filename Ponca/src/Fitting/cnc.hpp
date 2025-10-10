/**
Copyright (c) 2022
 Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr)
 Laboratory of Mathematics (CNRS, UMR 5807), University of Savoie, France,

All rights reserved.

*/

#pragma once

#include <iostream>
#include <random>
#include <vector>
#include "boundedRange.h"

namespace Ponca::internal {
    /*!
        \internal
        \brief Generates the triangles used by the CNC Fit depending on the method.
            As an output, pushes every generated triangle into the "triangles" vector and returns the number of triangle that was pushed into the List.
        \note Needs to be implemented for each triangle generation method by specializing the template over the TriangleGenerationMethod
    */
    template <TriangleGenerationMethod Method, typename P>
    struct TriangleGenerator {
        using VectorType = typename P::VectorType;
        template <typename PointContainer, typename IndicesGetter>
        static int generate(
            const PointContainer& /*points*/,
            const IndicesGetter& /*indicesGetter*/,
            const VectorType& /*_evalPointPos*/, const VectorType& /*_evalPointNormal*/,
            std::vector<Triangle<P>>& /*triangles*/
        ) {
            static_assert(true, "Triangle generation method not implemented!");
            return 0;
        }
    };

    /// Generates the triangles used by the CNC Fit using UniformGeneration
    template <typename P>
    struct TriangleGenerator<UniformGeneration, P> {
    private:
        static constexpr int maxTriangles {100};
    public:
        using VectorType = typename P::VectorType;

        template <typename PointContainer, typename IndicesGetter>
        static int generate(
            const PointContainer& points,
            const IndicesGetter& indicesGetter,
            const VectorType& /*_evalPointPos*/, const VectorType& /*_evalPointNormal*/,
            std::vector<Triangle<P>>& triangles
        ) {
            int nb_vt = 0; // Number of valid generated triangles

            for (int i = 0; i < maxTriangles; ++i) {
                // Randomly select triangles
                int i1 = indicesGetter[Eigen::internal::random<int>(0, indicesGetter.size()-1)];
                int i2 = indicesGetter[Eigen::internal::random<int>(0, indicesGetter.size()-1)];
                int i3 = indicesGetter[Eigen::internal::random<int>(0, indicesGetter.size()-1)];
                if (i1 == i2 || i1 == i3 || i2 == i3) continue;

                triangles.push_back(internal::Triangle<P>(points[i1], points[i2], points[i3]));
                nb_vt++;
            }
            return nb_vt;
        }
    };

    /// Generates the triangles used by the CNC Fit using HexagramGeneration
    template <typename P>
    struct TriangleGenerator<HexagramGeneration, P> {
        using VectorType = typename P::VectorType;
        using Scalar = typename P::Scalar;

        template <typename PointContainer, typename IndicesGetter>
        static int generate(
            const PointContainer& points,
            const IndicesGetter& indicesGetter,
            const VectorType& _evalPointPos, const VectorType& _evalPointNormal,
            std::vector<Triangle<P>>& triangles
        ) {
            // Hexagram
            std::array< Scalar    ,    6 > _distance2;
            std::array< VectorType,    6 > _targets;

            Scalar avg_normal  = Scalar(0.5);
            // BIN
            VectorType c = _evalPointPos;
            VectorType n = _evalPointNormal;
            VectorType a;
            a.setZero();

            int iSource = -1;
            Scalar avg_d = Scalar(0);

            for ( int index : indicesGetter ) {
                avg_d += ( points[ index ].pos() - c ).norm();
                a     += points[ index ].normal();
                // if avg_d == 0 then it is the evalPoint
                if ( iSource == -1 && points[ index ].pos() == c  ) {
                    iSource = index;
                }
            }

            a     /= a.norm();
            n      = ( Scalar(1) - avg_normal ) * n + avg_normal * a;
            n     /= n.norm();
            avg_d /= indicesGetter.size();

            const int m = ( std::abs( n[0] ) > std::abs ( n[1] ))
                    ? ( ( std::abs( n[0] ) ) > std::abs( n[2] ) ? 0 : 2 )
                    : ( ( std::abs( n[1] ) ) > std::abs( n[2] ) ? 1 : 2 ) ;

            const VectorType e =
                ( m == 0 ) ? VectorType( Scalar(0), Scalar(1), Scalar(0) ) :
                ( m == 1 ) ? VectorType( Scalar(0), Scalar(0), Scalar(1) ) :
                             VectorType( Scalar(1), Scalar(0), Scalar(0) ) ;

            VectorType u = n.cross( e );
            VectorType v = n.cross( u );
            u /= u.norm();
            v /= v.norm();

            std::array<int, 6> indices = {iSource, iSource, iSource, iSource, iSource, iSource};

            for ( int i = 0 ; i < 6 ; i++ ) {
                _distance2 [ i ] = avg_d * avg_d;
                _targets   [ i ] = avg_d * ( u * std::cos( i * M_PI / 3.0 ) + v * std::sin( i * M_PI / 3.0 ) );
            }

            for ( int index : indicesGetter ) {
                VectorType p = points[ index ].pos();
                if ( p == c ) continue;

                const VectorType d = p - c;
                for ( int j = 0 ; j < 6 ; j++ ){
                    const Scalar d2 = ( d - _targets[ j ]).squaredNorm();
                    if ( d2 < _distance2[ j ] ){
                        indices[ j ] = index;
                        _distance2[ j ] = d2;
                    }
                }
            }
            triangles.push_back(internal::Triangle<P>(points[indices[0]], points[indices[2]], points[indices[4]]));
            triangles.push_back(internal::Triangle<P>(points[indices[1]], points[indices[3]], points[indices[5]]));

            return 2;
        }
    };

    /// Generates the triangles used by the CNC Fit using AvgHexagramGeneration
    template <typename P>
    struct TriangleGenerator<AvgHexagramGeneration, P> {
        using VectorType = typename P::VectorType;
        using Scalar = typename P::Scalar;

        template <typename PointContainer, typename IndicesGetter>
        static int generate(
            const PointContainer& points,
            const IndicesGetter& indicesGetter,
            const VectorType& _evalPointPos, const VectorType& _evalPointNormal,
            std::vector<Triangle<P>>& triangles
        ) {
            VectorType c = _evalPointPos;
            VectorType n = _evalPointNormal;
            Scalar avg_d = Scalar(0);
            VectorType a = VectorType::Zero();

            std::array< VectorType,6 > array_avg_normals{VectorType::Zero()};
            std::array< VectorType,6 > array_avg_points{VectorType::Zero()};
            std::array< int, 6 >    array_nb {0};

            std::array< VectorType,    6 > _targets;
            Scalar avg_normal  = Scalar(0.5);

            for ( int index : indicesGetter ) {
                avg_d += ( points[ index ].pos() - c ).norm();
                a     += points[ index ].normal();
            }

            a     /= a.norm();
            n      = ( Scalar(1) - avg_normal ) * n + avg_normal * a;
            n     /= n.norm();
            avg_d /= indicesGetter.size();

            const int m = ( std::abs( n[0] ) > std::abs ( n[1] ))
                    ? ( ( std::abs( n[0] ) ) > std::abs( n[2] ) ? 0 : 2 )
                    : ( ( std::abs( n[1] ) ) > std::abs( n[2] ) ? 1 : 2 );

            const VectorType e = ( m == 0 ) ? VectorType( 0, 1, 0 ) :
                                 ( m == 1 ) ? VectorType( 0, 0, 1 ) :
                                              VectorType( 1, 0, 0 ) ;
            VectorType u = n.cross( e );
            VectorType v = n.cross( u );
            u /= u.norm();
            v /= v.norm();

            for (int i = 0 ; i < 6 ; i++ ) {
                _targets[ i ]          = avg_d * ( u * std::cos( i * M_PI / 3.0 ) + v * std::sin( i * M_PI / 3.0 ) );
                array_avg_normals[ i ] = VectorType::Zero();
                array_avg_points[ i ]  = VectorType::Zero();
            }

            for (int index : indicesGetter) {
                if ( points[ index ].pos() == c )
                    continue;

                VectorType p = points[ index ].pos() - c;
                int best_k = 0;
                Scalar best_d2 = ( p - _targets[ 0 ] ).squaredNorm();
                for (int k = 1 ; k < 6 ; k++) {
                    const Scalar d2 = ( p - _targets[ k ] ).squaredNorm();
                    if ( d2 < best_d2 ) {
                        best_k = k;
                        best_d2 = d2;
                    }
                }
                array_avg_normals[ best_k ] += points[ index ].normal();
                array_avg_points [ best_k ] += points[ index ].pos();
                array_nb[ best_k ] += 1;
            }

            for (int i = 0 ; i < 6 ; i++) {
                if ( array_nb[ i ] == 0 ) {
                    array_avg_normals[ i ] = n;
                    array_avg_points [ i ] = c;
                } else {
                    array_avg_normals[ i ] /= array_avg_normals[ i ].norm();
                    array_avg_points [ i ] /= array_nb[ i ];
                }
            }

            std::array <VectorType, 3> t1_points  = { array_avg_points[0] , array_avg_points[2] , array_avg_points[4] };
            std::array <VectorType, 3> t1_normals = { array_avg_normals[0], array_avg_normals[2], array_avg_normals[4] };

            std::array <VectorType, 3> t2_points  = { array_avg_points[1] , array_avg_points[3] , array_avg_points[5] };
            std::array <VectorType, 3> t2_normals = { array_avg_normals[1], array_avg_normals[3], array_avg_normals[5] };

            triangles.push_back(internal::Triangle<P>(t1_points, t1_normals));
            triangles.push_back(internal::Triangle<P>(t2_points, t2_normals));

            return 2;
        }
    };

    /// Generates the triangles used by the CNC Fit using IndependentGeneration
    template <typename P>
    struct TriangleGenerator<IndependentGeneration, P> {
    private:
        static constexpr int maxTriangles {100};
    public:
        using VectorType = typename P::VectorType;
        using Scalar = typename P::Scalar;

        template <typename PointContainer, typename IndicesGetter>
        static int generate(
            const PointContainer& points,
            const IndicesGetter& indicesGetter,
            const VectorType& /*_evalPointPos*/, const VectorType& /*_evalPointNormal*/,
            std::vector<Triangle<P>>& triangles
        ) {
            int nb_vt = 0; // Number of valid generated triangles

            // Makes a new array to shuffle
            std::vector<int> indices(indicesGetter.size());
            for (int i = 0; i < indicesGetter.size() ; ++i)
                indices[i] = indicesGetter[i];

            // Shuffles the neighbors
            std::random_device rd;
            std::mt19937 rg(rd());
            std::shuffle(indices.begin(), indices.end(), rg);

            // Compute the triangles
            triangles.clear();
            for (const int max_triangles = std::min(maxTriangles, static_cast<int>(indicesGetter.size()) / 3); nb_vt < max_triangles-2; nb_vt++) {
                int i1 = indices[nb_vt];
                int i2 = indices[nb_vt+1];
                int i3 = indices[nb_vt+2];
                triangles.push_back(internal::Triangle<P>(points[i1], points[i2], points[i3]));
            }
            return nb_vt;
        }
    };
} // namespace Ponca::internal

namespace Ponca {
    template < class P, TriangleGenerationMethod M>
    template <typename PointContainer>
    FIT_RESULT CNC<P, M>::compute( const PointContainer& points ) {
        init();
        internal::BoundedIntRange indicesSample( points.size() ); // Provides an index iterator and randomizer based on the number of points
        _nb_vt = internal::TriangleGenerator<M, P>::generate( points, indicesSample, _evalPointPos, _evalPointNormal, _triangles);

        return finalize();
    }

    template < class P, TriangleGenerationMethod M>
    template <typename IndexContainer, typename PointContainer>
    FIT_RESULT CNC<P, M>::computeWithIds( const IndexContainer& ids, const PointContainer& points ) {
        init();
        _nb_vt = internal::TriangleGenerator<M, P>::generate( points, ids, _evalPointPos, _evalPointNormal, _triangles);
        return finalize();
    }

    template < class P, TriangleGenerationMethod M>
    FIT_RESULT CNC<P, M>::finalize( ) {
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
        if (_A != Scalar(0)) {
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
} // namespace Ponca
