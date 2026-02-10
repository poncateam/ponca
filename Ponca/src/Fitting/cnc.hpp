/**
Copyright (c) 2022
 Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr)
 Laboratory of Mathematics (CNRS, UMR 5807), University of Savoie, France,

All rights reserved.

*/

#pragma once

#include <random>

namespace Ponca::internal {
    /*! \internal
     *
     * \brief Generates the triangles used by the CNC Fit depending on the method.
     *
     * As an output, pushes every generated triangle into the `triangles` vector and returns the number of triangle that was pushed to the List.
     *
     * \note Needs to be implemented for each triangle generation method by specializing the template over the `TriangleGenerationMethod` enum.
     *
     * \see TriangleGenerationMethod
     */
    template <TriangleGenerationMethod Method, typename P>
    struct TriangleGenerator {
        using VectorType = typename P::VectorType;
        template <typename IndexRange, typename PointContainer, typename NeighborFilter>
        static FIT_RESULT generate(
            const IndexRange& /*ids*/,
            const PointContainer& /*points*/,
            const NeighborFilter& /*w*/,
            std::vector<Triangle<P>>& /*triangles*/
        ) {
            throw std::invalid_argument("Triangle generation method not implemented!");
        }
    };

    /*! \copydoc TriangleGenerator
     *
     * Generates the triangles using UniformGeneration
     */
    template <typename P>
    struct TriangleGenerator<UniformGeneration, P> {
    private:
        static constexpr int maxTriangles {100};
    public:
        using VectorType = typename P::VectorType;
        using Scalar     = typename P::Scalar;

        template <typename IndexRange, typename PointContainer, typename NeighborFilter>
        static FIT_RESULT generate(
            const IndexRange& ids,
            const PointContainer& points,
            const NeighborFilter& w,
            std::vector<Triangle<P>>& triangles
        ) {
            // Makes a new array
            std::vector<int> indices;
            for ( int index : ids ) {
                // Skip the points that are outside the kernel radius
                if (w(points[ index ]).first == Scalar(0.))
                    continue;
                indices.push_back(index);
            }
            if (indices.empty())
                return UNDEFINED;

            const int lastIndex = int(indices.size()) - 1;
            for (int i = 0; i < maxTriangles; ++i) {
                // Randomly select triangles
                int i1 = indices[Eigen::internal::random<int>(0, lastIndex)];
                int i2 = indices[Eigen::internal::random<int>(0, lastIndex)];
                int i3 = indices[Eigen::internal::random<int>(0, lastIndex)];
                if (i1 == i2 || i1 == i3 || i2 == i3) continue;

                triangles.push_back(internal::Triangle<P>(points[i1], points[i2], points[i3]));
            }
            return STABLE;
        }
    };

    /*! \copydoc TriangleGenerator
     *
     * Generates the triangles using IndependentGeneration
     */
    template <typename P>
    struct TriangleGenerator<IndependentGeneration, P> {
    private:
        static constexpr int maxTriangles {100};
    public:
        using VectorType = typename P::VectorType;
        using Scalar     = typename P::Scalar;

        template <typename IndexRange, typename PointContainer, typename NeighborFilter>
        static FIT_RESULT generate(
            const IndexRange& ids,
            const PointContainer& points,
            const NeighborFilter& w,
            std::vector<Triangle<P>>& triangles
        ) {
            // Makes a new array to shuffle
            std::vector<int> indices;
            for ( int index : ids ) {
                // Skip the points that are outside the kernel radius
                if (w(points[ index ]).first == Scalar(0.))
                    continue;
                indices.push_back(index);
            }
            if (indices.empty())
                return UNDEFINED;

            // Shuffles the neighbors
            std::random_device rd;
            std::mt19937 rg(rd());
            std::shuffle(indices.begin(), indices.end(), rg);

            // Compute the triangles
            triangles.clear();
            const int max_triangles = std::min(maxTriangles, int(ids.size() / 3));
            for (int nb_vt = 0; nb_vt < max_triangles-2; nb_vt++) {
                int i1 = indices[nb_vt];
                int i2 = indices[nb_vt+1];
                int i3 = indices[nb_vt+2];
                triangles.push_back(internal::Triangle<P>(points[i1], points[i2], points[i3]));
            }
            return STABLE;
        }
    };

    template<typename P>
    struct HexagramBase {
        using Scalar = typename P::Scalar;
        static constexpr Scalar avg_normal_coef {Scalar(0.5)};
    };


    /*! \copydoc TriangleGenerator
     *
     * Generates the triangles using HexagramGeneration
     */
    template <typename P>
    struct TriangleGenerator<HexagramGeneration, P> : protected HexagramBase<P> {
        using VectorType = typename P::VectorType;
        using Scalar     = typename P::Scalar;

        template <typename IndexRange, typename PointContainer, typename NeighborFilter>
        static FIT_RESULT generate(
            const IndexRange& ids,
            const PointContainer& points,
            const NeighborFilter& w,
            std::vector<Triangle<P>>& triangles
        ) {
            PONCA_MULTIARCH_STD_MATH(abs);
            // Compute normal and maximum distance.
            VectorType c = w.evalPos();
            VectorType n = w.evalNormal();
            VectorType a {VectorType::Zero()};
            Scalar avg_d = Scalar(0);

            bool isUndefined = true;
            for ( int index : ids ) {
                // Skip the points that are outside the kernel radius
                if (w(points[ index ]).first == Scalar(0.))
                    continue;
                auto p = points[ index ];
                avg_d += ( p.pos() - c ).norm();
                a     += p.normal();
                isUndefined = false;
            }
            if (isUndefined)
                return UNDEFINED;

            a     /= a.norm();
            n      = ( Scalar(1) - HexagramBase<P>::avg_normal_coef ) * n + HexagramBase<P>::avg_normal_coef * a;
            n     /= n.norm();
            avg_d /= ids.size();

            // Define basis for sector analysis.
            const int m = abs( n[0] ) > abs( n[1] )
                      ? ( abs( n[0] ) > abs( n[2] ) ? 0 : 2 )
                      : ( abs( n[1] ) > abs( n[2] ) ? 1 : 2 ) ;

            const VectorType e =
                ( m == 0 ) ? VectorType( 0, 1, 0 ) :
                ( m == 1 ) ? VectorType( 0, 0, 1 ) :
                             VectorType( 1, 0, 0 ) ;

            VectorType u = n.cross( e );
            VectorType v = n.cross( u );
            u /= u.norm();
            v /= v.norm();

            std::array<VectorType , 6> positions ;
            std::array<VectorType , 6> normals   ;
            std::array< Scalar    , 6 > distance2;
            std::array< VectorType, 6 > targets;

            for ( int i = 0 ; i < 6 ; i++ ) {
                distance2[ i ] = avg_d * avg_d;
                targets  [ i ] = avg_d * ( u * cos(i * M_PI / 3.0) + v * sin(i * M_PI / 3.0) );
                positions[ i ] = w.evalPos();
                normals  [ i ] = w.evalNormal();
            }

            // Compute closest points.
            for ( int index : ids ) {
                // Skip the points that are outside the kernel radius
                if (w(points[ index ]).first == Scalar(0.))
                    continue;

                VectorType p = points[ index ].pos();
                if ( p == c ) continue; // Skip the eval point
                const VectorType d = p - c;

                for ( int j = 0 ; j < 6 ; j++ ) {
                    const Scalar d2 = ( d - targets[ j ]).squaredNorm();
                    if ( d2 < distance2[ j ] ) {
                        distance2[ j ] = d2;
                        positions[ j ] = p;
                        normals  [ j ] = points[ index ].normal();
                    }
                }
            }
            triangles.push_back(internal::Triangle<P>({positions[0] , positions[2] , positions[4]}, {normals[0] , normals[2], normals[4]}));
            triangles.push_back(internal::Triangle<P>({positions[1] , positions[3] , positions[5]}, {normals[1] , normals[3], normals[5]}));

            return STABLE;
        }
    };

    /*! \copydoc TriangleGenerator
     *
     * Generates the triangles using AvgHexagramGeneration
     */
    template <typename P>
    struct TriangleGenerator<AvgHexagramGeneration, P> : protected HexagramBase<P> {
        using VectorType = typename P::VectorType;
        using Scalar     = typename P::Scalar;

        template <typename IndexRange, typename PointContainer, typename NeighborFilter>
        static FIT_RESULT generate(
            const IndexRange& ids,
            const PointContainer& points,
            const NeighborFilter& w,
            std::vector<Triangle<P>>& triangles
        ) {
            // Compute normal and maximum distance.
            VectorType c = w.evalPos();
            VectorType n = w.evalNormal();
            VectorType a = w.evalNormal();
            Scalar avg_d = Scalar(0);

            std::array< VectorType, 6 > targets;
            Scalar avg_normal  = Scalar(0.5);

            bool isUndefined = true;
            for ( int index : ids ) {
                if (w(points[ index ]).first == Scalar(0.))
                    continue; // Skip the points that are outside the kernel radius
                a     += points[ index ].normal();
                avg_d += ( points[ index ].pos() - c ).norm();
                isUndefined = false;
            }
            if (isUndefined)
                return UNDEFINED;

            a     /= a.norm();
            n      = ( Scalar(1) - HexagramBase<P>::avg_normal_coef ) * n + HexagramBase<P>::avg_normal_coef * a;
            n     /= n.norm();
            avg_d /= ids.size();

            // Define basis for sector analysis.
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

            // Initialize the average values
            std::array< VectorType, 6 > array_avg_normals;
            std::array< VectorType, 6 > array_avg_pos;
            std::array< int, 6 >    array_nb {};
            for (int i = 0 ; i < 6 ; i++ ) {
                targets[ i ]          = avg_d * ( u * cos(i * M_PI / 3.0) + v * sin(i * M_PI / 3.0) );
                array_avg_normals[ i ] = VectorType::Zero();
                array_avg_pos    [ i ] = VectorType::Zero();
            }

            // Compute closest points.
            for (int index : ids) {
                VectorType p = points[ index ].pos() - c;
                int best_k = 0;
                Scalar best_d2 = ( p - targets[ 0 ] ).squaredNorm();
                for (int k = 1 ; k < 6 ; k++) {
                    const Scalar d2 = ( p - targets[ k ] ).squaredNorm();
                    if ( d2 < best_d2 ) {
                        best_k = k;
                        best_d2 = d2;
                    }
                }
                array_avg_normals[ best_k ] += points[ index ].normal();
                array_avg_pos    [ best_k ] += points[ index ].pos();
                array_nb[ best_k ] += 1;
            }

            for (int i = 0 ; i < 6 ; i++) {
                if ( array_nb[ i ] == 0 ) {
                    array_avg_normals[ i ] = w.evalNormal();
                    array_avg_pos    [ i ] = w.evalPos();
                } else {
                    array_avg_normals[ i ] /= array_avg_normals[ i ].norm();
                    array_avg_pos    [ i ] /= array_nb[ i ];
                }
            }

            triangles.push_back(internal::Triangle<P>(
                { array_avg_pos[0] , array_avg_pos[2] , array_avg_pos[4] },
                { array_avg_normals[0], array_avg_normals[2], array_avg_normals[4] }
            ));
            triangles.push_back(internal::Triangle<P>(
                { array_avg_pos[1] , array_avg_pos[3] , array_avg_pos[5] },
                { array_avg_normals[1], array_avg_normals[3], array_avg_normals[5] }
            ));
            return STABLE;
        }
    };
} // namespace Ponca::internal

namespace Ponca {
    template < class P, TriangleGenerationMethod M>
    template <typename PointContainer>
    FIT_RESULT CNC<P, M>::compute( const PointContainer& points ) {
        init();
        std::vector<unsigned int> indicesSample(points.size());
        std::iota(indicesSample.begin(), indicesSample.end(), 0);

        m_eCurrentState = internal::TriangleGenerator<M, P>::generate( indicesSample, points, m_nFilter, m_triangles);
        if (m_eCurrentState != STABLE) return m_eCurrentState;
        m_nb_vt = m_triangles.size();
        return finalize();
    }

    template < class P, TriangleGenerationMethod M>
    template <typename IndexRange, typename PointContainer>
    FIT_RESULT CNC<P, M>::computeWithIds( const IndexRange& ids, const PointContainer& points ) {
        init();
        m_eCurrentState = internal::TriangleGenerator<M, P>::generate( ids, points, m_nFilter, m_triangles);
        if (m_eCurrentState != STABLE) return m_eCurrentState;
        m_nb_vt = m_triangles.size();
        return finalize();
    }

    template < class P, TriangleGenerationMethod M>
    FIT_RESULT CNC<P, M>::finalize( ) {
        m_A = Scalar(0);
        m_H = Scalar(0);
        m_G = Scalar(0);

        MatrixType localT = MatrixType::Zero();

        for (int t = 0; t < m_nb_vt; ++t) {
            // Simple estimation.
            Scalar tA = m_triangles[t].mu0InterpolatedU();
            if (tA < -internal::CNCEigen<P>::epsilon) {
                m_A     -= tA;
                m_H     += m_triangles[t].template mu1InterpolatedU<true>();
                m_G     += m_triangles[t].template mu2InterpolatedU<true>();
                localT  += m_triangles[t].template muXYInterpolatedU<true>();
            } else if (tA > internal::CNCEigen<P>::epsilon) {
                m_A     += tA;
                m_H     += m_triangles[t].mu1InterpolatedU();
                m_G     += m_triangles[t].mu2InterpolatedU();
                localT  += m_triangles[t].muXYInterpolatedU();
            }
        } // end for t

        m_T11 = localT(0,0);
        m_T12 = 0.5 * (localT(0,1) + localT(1,0));
        m_T13 = 0.5 * (localT(0,2) + localT(2,0));
        m_T22 = localT(1,1);
        m_T23 = 0.5 * (localT(1,2) + localT(2,1));
        m_T33 = localT(2,2);

        MatrixType T;
        if (m_A != Scalar(0)) {
            T  << m_T11, m_T12, m_T13,
                  m_T12, m_T22, m_T23,
                  m_T13, m_T23, m_T33;
            T  /= m_A;
            m_H /= m_A;
            m_G /= m_A;
        } else {
            m_H = Scalar(0);
            m_G = Scalar(0);
        }

        std::tie (m_k2, m_k1, m_v2, m_v1) = internal::CNCEigen<P>::curvaturesFromTensor(T, 1.0, m_nFilter.evalNormal());

        return STABLE;
    }
} // namespace Ponca
