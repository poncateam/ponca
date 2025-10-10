#pragma once

namespace Ponca::internal
{
    /*!
        \internal
        \brief Parent class to manage integers inside a boundary. Can be iterated over or used to generates a random integer inside the boundary
        \note Calling the on this objet after the initialization of its boundary generates a random integer in between [_nMin, _nMax[
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
        /// \brief Simply verify that n is in bounds and returns it. Can be overwritten to something else in children class
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

    /*!
        \internal
        \brief Class that provides usefull operators to iterate over an STL-like container of ints (to manage indices map for instance).
        It will store a reference to the container, and can iterate over it or pick an element from it (random or not)
        \inherit BoundedIntRange
    */
    template<typename Container>
    class IndexMap : public BoundedIntRange {
    private:
        Container& _elements;
    public:
        /// \internal
        /// \brief The operations on the container can be restricted to the bounds : [ _nMin, _nMax [
        /// \note nMin is optional and is equal to zero by default, but an nMax must be provided
        IndexMap(Container& elements, const int nMax, const int nMin = 0) :
            BoundedIntRange(nMax, nMin), _elements(elements) { }

        /// \internal
        /// \brief Returns an elements from the integer container after verifying that the index i is in bounds
        constexpr int operator[](const int i) const {
            verifyBounds(i);
            return _elements[i];
        }

        /// \internal
        /// \brief Makes the class iterable
        class Iterator {
            const IndexMap& _parent;
            int _current;
        public:
            Iterator(const IndexMap& parent, const int start)
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
}