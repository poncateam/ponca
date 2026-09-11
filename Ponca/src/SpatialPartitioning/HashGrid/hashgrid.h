/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/
#pragma once

#include <type_traits>
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numeric>
#include <memory>
#include <vector>
#include <cmath>

#include "../defines.h"
#include "../utils.h"
// #include "../../Common/concepts.h"

namespace Ponca
{
    template <typename Point, typename Int, typename IntContainer, typename CellContainer>
    struct HashGridBuffer
    {
        constexpr static std::uint8_t Dim = Point::Dim;
        using Coords                      = std::array<Int, Dim>;
        using VectorType                  = typename Point::VectorType;
        using Scalar                      = typename Point::Scalar;

        VectorType lowerleft{};
        VectorType cellSize{};
        Coords cellCount{};
        size_t totalCellCount{};

        IntContainer indices{};
        CellContainer cells{};
    };

    template <typename Point, std::uint32_t MaxK, typename Int = std::uint32_t>
    struct NearestNeighborQuery
    {
        struct ReturnType
        {
            std::array<Int, MaxK> indices;
            std::array<typename Point::Scalar, MaxK> dists;
            std::uint32_t count;
        };

        ReturnType MakeReturn() const { return ReturnType{}; }

        std::uint32_t k;
        typename Point::VectorType loc;
    };

    template <typename Int>
    struct HashGridCell
    {
        Int start;
        Int end;

        Int count() const { return end - start; }
        bool empty() const { return end == start; }
    };

    template <typename Point, typename Int = std::uint32_t, typename IntContainer = std::vector<Int>,
              typename CellContainer = std::vector<HashGridCell<Int>>, bool IsDense = true>
    struct HashGrid
    {
        static constexpr std::uint32_t Dim = Point::Dim;
        using BufferType                   = HashGridBuffer<Point, Int, IntContainer, CellContainer>;
        using Scalar                       = typename Point::Scalar;
        using Coords                       = typename BufferType::CellCount;
        using VectorType                   = typename BufferType::VectorType;

        HashGrid(BufferType&& buffer) : m_buffer(std::move(buffer)) {}

        template <typename PointContainer>
        HashGrid(const Coords& cellCount, const PointContainer& points)
        {
            VectorType lower, upper;
            GetBounds(points, lower, upper);

            m_buffer.totalCellCount = 1;
            for (size_t d = 0; d < Dim; ++d)
            {
                m_buffer.lowerleft[d] = lower[d];
                m_buffer.cellSize[d]  = (upper[d] - lower[d]) / cellCount[d];
                m_buffer.cellCount[d] = cellCount[d];
                m_buffer.totalCellCount *= m_buffer.cellCount;
            }

            build(points);
        }

        template <typename PointContainer>
        HashGrid(const VectorType& cellSize, const PointContainer& points)
        {
            VectorType lower, upper;
            GetBounds(points, lower, upper);

            m_buffer.totalCellCount = 1;
            for (size_t d = 0; d < Dim; ++d)
            {
                m_buffer.lowerleft[d] = lower[d];
                m_buffer.cellSize[d]  = cellSize[d];
                m_buffer.cellCount[d] = std::floor((upper[d] - lower[d]) / cellSize[d]);
                m_buffer.totalCellCount *= m_buffer.cellCount;
            }

            build(points);
        }

        template <typename PointContainer>
        HashGrid(const VectorType& lower, const VectorType& upper, const VectorType& cellSize,
                 const PointContainer& points)
        {
            m_buffer.totalCellCount = 1;
            for (size_t d = 0; d < Dim; ++d)
            {
                m_buffer.lowerleft[d] = lower[d];
                m_buffer.cellSize[d]  = cellSize[d];
                m_buffer.cellCount[d] = std::floor((upper[d] - lower[d]) / cellSize[d]);
                m_buffer.totalCellCount *= m_buffer.cellCount;
            }

            build(points);
        }

        template <typename PointContainer>
        HashGrid(const VectorType& lower, const VectorType& upper, const Coords& cellCount,
                 const PointContainer& points)
        {
            m_buffer.totalCellCount = 1;
            for (size_t d = 0; d < Dim; ++d)
            {
                m_buffer.lowerleft[d] = lower[d];
                m_buffer.cellSize[d]  = (upper[d] - lower[d]) / cellCount[d];
                m_buffer.cellCount[d] = cellCount[d];
                m_buffer.totalCellCount *= m_buffer.cellCount;
            }

            build(points);
        }

        const BufferType& GetBuffer() const { return m_buffer; }

        template <std::uint32_t MaxK>
        decltype(auto) Query(const NearestNeighborQuery<Point, MaxK, Int>& query)
        {
            const std::uint32_t k = std::min(MaxK, query.k);
            const auto maxSize    = std::max_element(m_buffer.cellSize.begin(), m_buffer.cellSize.end());
            auto ret              = query.MakeReturn();

            Int hash       = Hash(query.loc);
            Int level      = 0;
            Scalar maxDist = std::numeric_limits<typename Point::Scalar>::max();
            for (std::uint32_t i = 0; i < k; ++i)
                ret.dists[i] = std::numeric_limits<typename Point::Scalar>::max();

            while (maxDist > 1.5 * level * maxSize)
            {
                SearchRing(hash, level, [&](const auto& _, const HashGridCell<Int>& cell) {
                    for (Int i = cell.start; i < cell.end; ++cell)
                    {
                    }
                });
            }
        }

    private:
        template <typename Visitor>
        void SearchRing(Int centerHash, Int radius, Visitor&& visit)
        {
            const Coords center = ReverseHash(centerHash);

            for (Int axis = 0; axis < Dim; ++axis)
            {
                for (char side : {0, 1})
                {
                    Coords current{};
                    current[axis] = (side == 0) ? 0 : radius;

                    bool done = false;
                    while (!done)
                    {
                        // Bound checking without overflows
                        {
                            Coords frameCoords{};
                            bool isOk = true;
                            for (Int d = 0; d < Dim; ++d)
                            {
                                // Bound checking
                                const bool leftOk  = (center[d] >= radius);
                                const bool rightOk = (center[d] + radius < m_buffer.cellCount[d]);
                                isOk               = isOk && leftOk && rightOk;
                                if (leftOk && rightOk)
                                    frameCoords[d] = center[d] - radius + current[d];
                            }

                            if (isOk)
                                visit(frameCoords, m_buffer.cells[Hash(frameCoords)]);
                        }

                        for (Int i = 0; i < Dim; ++i)
                        {
                            if (i == axis)
                            {
                                done = done | (i == Dim - 1);
                                continue;
                            }

                            current[i]++;
                            if (current[i] < radius)
                            {
                                break;
                            }
                            else
                            {
                                current[i] = 0;
                                if (i == Dim - 1 || (i == Dim - 2 && axis == Dim - 1))
                                {
                                    done = true;
                                }
                            }
                        }
                    }
                }
            }
        }

        template <typename PointContainer>
        void build(const PointContainer& points)
        {
            m_buffer.indices.resize(points.size());
            if constexpr (IsDense)
                m_buffer.cellCount.resize(points.size());

            std::iota(m_buffer.indices.begin(), m_buffer.indices.end());
            std::sort(m_buffer.indices.begin(), m_buffer.indices.end(),
                      [&](Int i, Int j) { return Hash(points[i]) < Hash(points[j]); });

            Int currentHash = Hash(points[0]);

            m_buffer.cells[currentHash].start = 0;
            m_buffer.cells[currentHash].end   = 1;
            for (size_t i = 1; i < m_buffer.indices.size(); ++i)
            {
                Int hash = Hash(points[i]);
                if (hash == currentHash)
                {
                    m_buffer.cells[currentHash].end++;
                }
                else
                {
                    m_buffer.cells[hash].start = i;
                    m_buffer.cells[hash].end   = i + 1;

                    currentHash = hash;
                }
            }
        }

        PONCA_MULTIARCH Int Hash(const Point& p) const
        {
            Int idx    = 0;
            Int stride = 1;
            for (size_t d = 0; d < Dim; ++d)
            {
                const Int coord = (p[d] - m_buffer.lowerleft[d]) / m_buffer.cellSize[d];
                idx += coord * stride;
                stride *= m_buffer.cellCount[d];
            }

            return idx;
        }

        PONCA_MULTIARCH Int Hash(const Coords& p) const
        {
            Int idx    = 0;
            Int stride = 1;
            for (size_t d = 0; d < Dim; ++d)
            {
                idx += p[d] * stride;
                stride *= m_buffer.cellCount[d];
            }

            return idx;
        }

        PONCA_MULTIARCH Coords ReverseHash(Int hash) const
        {
            Coords coords{};
            for (size_t d = 0; d < Dim; ++d)
            {
                coords[d] = hash % m_buffer.cellCount[d];
                hash /= m_buffer.cellCount[d];
            }
            return coords;
        }

        template <typename PointContainer>
        void GetBounds(const PointContainer& container, VectorType& lower, VectorType& upper)
        {
            for (size_t d = 0; d < Point::Dim; ++d)
            {
                lower[d] = std::numeric_limits<Scalar>::max();
                upper[d] = std::numeric_limits<Scalar>::min();
            }

            for (size_t i = 0; i < container.size(); ++i)
            {
                const auto& point = container[i];
                for (size_t d = 0; d < Point::Dim; ++d)
                {
                    lower[d] = (lower[d] > point[d]) ? point[d] : lower[d];
                    upper[d] = (upper[d] > point[d]) ? point[d] : upper[d];
                }
            }
        }

    private:
        BufferType m_buffer;
    };

} // namespace Ponca
