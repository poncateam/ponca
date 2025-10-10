#pragma once

namespace Ponca::internal
{
    /*!
        \internal
        \brief Parent class to manage integers inside a boundary. Can be iterated over or it
    */
    class BoundedIntRange {
    public:
        const int _nMin; // included
        const int _nMax; // excluded
        /// \internal
        /// \brief The operations on this can be restricted to the bounds : [ _nMin, _nMax [
        /// \note nMin is optional and is equal to zero by default, but an nMax must be provided
        explicit BoundedIntRange( const int nMax, const int nMin = 0 ) : _nMin(nMin), _nMax(nMax) { }

        constexpr void verifyBounds(const int n) const {
            if (_nMin > n || n >= _nMax)
                throw std::runtime_error(
                    "Index values must be in range :"
                    + std::to_string(_nMin) + " <= i < " + std::to_string(_nMax)
                    + " But got result : " + std::to_string(n));
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
            return _nMax - _nMin;
        }

        /// \internal
        /// \brief Makes the class iterable
        class Iterator {
            int _current;
        public:
            explicit Iterator(const int start) : _current(start) {}

            int operator*() const { return _current; }
            Iterator& operator++() { ++_current; return *this; }

            bool operator!=(const Iterator& other) const {
                return _current != other._current;
            }
        };

        [[nodiscard]] Iterator begin() const { return Iterator(_nMin); }
        [[nodiscard]] Iterator end() const { return Iterator(_nMax); }
    };
}