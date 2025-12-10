#pragma once

namespace Ponca::internal
{
    /*!
        \internal
        \brief Parent class to manage integers inside a boundary. Can be iterated over or it
    */
    class BoundedIntRange {
    protected:
        const int m_nMin; // included in the bounded range
        const int m_nMax; // excluded from the bounded range
    public:
        /// \internal
        /// \brief The operations on this can be restricted to the bounds : [ _nMin, _nMax [
        /// \note nMin is optional and is equal to zero by default, but an nMax must be provided
        explicit BoundedIntRange( const int nMax, const int nMin = 0 ) : m_nMin(nMin), m_nMax(nMax) { }

        constexpr void verifyBounds(const int i) const {
            if (m_nMin > i || i >= m_nMax)
                throw std::runtime_error(
                    "Index values must be in range :"
                    + std::to_string(m_nMin) + " <= i < " + std::to_string(m_nMax)
                    + " But got result : " + std::to_string(i));
        }

        /// \internal
        /// \brief Simply verify that n is in bounds and returns it.
        constexpr int operator[](const int i) const {
            verifyBounds(i);
            return i;
        }

        /// \internal
        /// \brief Get the size of the range
        [[nodiscard]] int size() const {
            return m_nMax - m_nMin;
        }

        /// \internal
        /// \brief Makes the class iterable
        class Iterator {
            int m_current;
        public:
            explicit Iterator(const int start) : m_current(start) {}

            int operator*() const { return m_current; }
            Iterator& operator++() { ++m_current; return *this; }

            bool operator!=(const Iterator& other) const {
                return m_current != other.m_current;
            }
        };

        [[nodiscard]] Iterator begin() const { return Iterator(m_nMin); }
        [[nodiscard]] Iterator end() const { return Iterator(m_nMax); }
    };
}
