#pragma once

#include <algorithm>


/*!
 * \copybreif hasDuplicate(Iter first, Iter last)
 *
 * This function is restricted to random access iterators
 * (iterators that can access the container non-sequentially, by jumping around).
 *
 * \tparam RandomIt Random Access Iterator type.
 * \param first The beginning iterator (obtained using `std::begin(container)` on an STL-like container).
 * \param last The end iterator (obtained using `std::end(container)` on an STL-like container).
 * \return True if the container has a duplicate element inside it.
 */
template<class RandomIt>
bool randomAccessIteratorHasDuplicate(RandomIt first, RandomIt last)
{
    static_assert(std::is_base_of_v<
        std::random_access_iterator_tag,
        typename std::iterator_traits<RandomIt>::iterator_category
    >);

    return std::any_of(first, last, [&](const auto& cur)->bool
    {
        // next is the iterator pointing after the current value cur
        RandomIt next = first + (&cur - &(*first)) + 1;
        return std::find(next, last, cur) != last;
    });
}

/*!
 * \copybreif hasDuplicate(Iter first, Iter last)
 *
 * Works on any forward iterators
 * (iterators that can scan sequentially the container multiple time and edit its values).
 *
 * \tparam ForwardIt Forward Iterator type.
 * \param first The beginning iterator (obtained using `std::begin(container)` on an STL-like container).
 * \param last The end iterator (obtained using `std::end(container)` on an STL-like container).
 * \return True if the container has a duplicate element inside it.
 */
template<class ForwardIt>
bool forwardIteratorHasDuplicate(ForwardIt first, ForwardIt last)
{
    static_assert(std::is_base_of_v<
        std::forward_iterator_tag,
        typename std::iterator_traits<ForwardIt>::iterator_category
    >);

    for (ForwardIt it = first; it != last; ++it)
    {
        ForwardIt next = it;
        ++next;

        if (std::find(next, last, *it) != last)
            return true;
    }
    return false;
}

/*!
 * \brief Detects if the container (represented by the first / last iterators) contains duplicated elements.
 *
 * \tparam IteratorT An Iterator type (any Forward iterator).
 * \param first The beginning iterator (obtained using `std::begin(container)` on an STL-like container).
 * \param last The end iterator (obtained using `std::end(container)` on an STL-like container).
 * \return True if the container has a duplicate element inside it.
 */
template<typename IteratorT>
bool hasDuplicate(IteratorT first, IteratorT last)
{
    if constexpr (std::is_base_of_v<
        std::random_access_iterator_tag,
        typename std::iterator_traits<IteratorT>::iterator_category()
    >) return randomAccessIteratorHasDuplicate(first, last);

    return forwardIteratorHasDuplicate(first, last);
}

/*!
 * \brief Checks if the container contains duplicated items.
 *
 * Complexity = O(n^2)
 *
 * \tparam ContainerT An STL-like container.
 * \param container The container to check.
 * \return True if the container has a duplicate element inside it.
 */
template<class ContainerT>
bool hasDuplicate(ContainerT container)
{
    return hasDuplicate(std::begin(container), std::end(container));
}
