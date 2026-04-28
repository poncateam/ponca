/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

/*!
 * \file tests/src/common_containers.cpp
 * \brief Validate LimitedPriorityQueue, HashSet and BitSet
 */

#include "../common/testing.h"
#include "../common/testUtils.h"
#include <Ponca/src/Common/Containers/bitset.h>
#include <Ponca/src/Common/Containers/hashset.h>

#include <random>
#include <set>
#include <vector>

#include "Ponca/src/Common/Containers/limitedPriorityQueue.h"

using namespace Ponca;
using namespace std;

/*! \brief Builds a vector of indice that are shuffled
 * \return The number of total elements
 */
int makeShuffledIndexVector(vector<int>& indices, const int min, const int max)
{
    const int nbTotal = QUICK_TESTS ? 1 : Eigen::internal::random<int>(min, max-1)+1;
    indices.resize(nbTotal);

    std::iota(indices.begin(), indices.end(), 1);
    ranges::shuffle(indices, std::mt19937{std::random_device{}()});

    return nbTotal;
}

/*!
 * \brief Test the validity of Ponca Common Set of int classes by inserting random values
 * and comparing the results with standard std::set<int>
 *
 * \tparam IndexSet The Set of indices to evaluate
 * \tparam RandomFunctor The Random functor type
 * \param _maxIndex The maximum possible index (must be > 100)
 * \param _pickRandom The function to call when generating a random value
 */
template <typename IndexSet, typename RandomFunctor>
void testSetStandardCapabilities(const int _maxIndex, RandomFunctor _pickRandom)
{
    assert(_maxIndex > 100);
    IndexSet indexSet;

    // Verify that the set starts empty
    for (int i = 0; i < (QUICK_TESTS ? 1 : _maxIndex); ++i)
        VERIFY(!indexSet.contains(i));

    std::set<int> indexSetSTD; // Stores id
    IndexSet indexSetPonca;

    // Pick a number of elements to insert in the set
    const int nbInsertion = QUICK_TESTS ? 1 : Eigen::internal::random<int>(100, _maxIndex - 1);

    for (int i = 0; i < nbInsertion; ++i)
    {
        // Select a random index to add to the set
        const int idx = _pickRandom();
        VERIFY((indexSetSTD.contains(idx) == indexSetPonca.contains(idx))); // Check before insert
        VERIFY((indexSetSTD.insert(idx).second == indexSetPonca.insert(idx)));
        VERIFY((indexSetSTD.contains(idx) == indexSetPonca.contains(idx))); // Check after insert
    }

    for (int i = 0; i < nbInsertion; ++i)
    {
        VERIFY((indexSetSTD.contains(i) == indexSetPonca.contains(i)));
    }
}

/*!
 * \brief Test the validity of Ponca Common Set of int classes by randomly inserting positive values
 * and comparing the results with standard std::set<int>
 *
 * \tparam IndexSet The Set of indices to evaluate
 * \param _maxIndex The maximum possible index (must be > 100)
 */
template <typename IndexSet>
void testSetStandardCapabilities(const int _maxIndex)
{
    testSetStandardCapabilities<IndexSet>(_maxIndex, [&_maxIndex]()
    {
        return Eigen::internal::random<int>(0, _maxIndex - 1);
    });
}

/*
 * \brief Test the insert capabilities of the limited Set
 *
 * \param _maxIndex The maximum possible index
 * \param _setCapacity The maximum capacity of the limited set
 */
template <typename IndexSet>
void testLimitedSet(const int _maxIndex, const int _setCapacity)
{
    IndexSet indexSet;

    vector<int> indices; // To keep track that we insert each index only once
    const int nbTotalInsertion = // A number of insertion that exceed available total capacity
        makeShuffledIndexVector(indices, _setCapacity, _maxIndex);


    // Insert until we reach max capacity
    for (int i = 0; i < _setCapacity; ++i)
        VERIFY((indexSet.insert(indices[i])));

    // Test insert above capacity
    for (int i = _setCapacity; i < nbTotalInsertion; ++i)
        VERIFY(! (indexSet.insert(indices[i])));
}

template <int MAX_INSERT_SIZE>
void testLimitedQueue(const int _maxIndex, const int _setCapacity)
{
    LimitedPriorityQueue<int, MAX_INSERT_SIZE, std::greater<>> queue(_setCapacity);

    vector<int> indices; // To keep track that we insert each index only once
    const int nbTotalInsertion = QUICK_TESTS ? 1 : Eigen::internal::random<int>(_setCapacity, _maxIndex-1)+1;

    // Insert until we reach max capacity
    for (int i = 0; i < _setCapacity; ++i)
        VERIFY(queue.push(i));

    int i = 0;
    // Test insert above capacity
    for (int j = _setCapacity; j < nbTotalInsertion; ++j, ++i)
    {
        VERIFY(queue.bottom() == i);
        VERIFY(queue.push(j));
        // Verify that it removed the lowest element (bottom)
        VERIFY(queue.bottom() == i+1);
        // Verify the new top
        VERIFY(queue.top() == j);
    }
}

int main(const int argc, char** argv)
{
    if (!init_testing(argc, argv))
    {
        return EXIT_FAILURE;
    }
    // The upper limit of index range (index can't go higher than this number for BITSET)
    constexpr int MAX_INDEX = 10000;
    // The maximum size of a limited set (values won't be inserted after that)
    constexpr int MAX_INSERT_SIZE = 100;

    cout << "Check consistency of Set of Indices class" << endl;

    for (int i = 0; i < g_repeat; ++i)
    {
        CALL_SUBTEST((testSetStandardCapabilities<BitSet<MAX_INDEX>>(MAX_INDEX)));
        CALL_SUBTEST((testSetStandardCapabilities<HashSet<MAX_INDEX>>(MAX_INDEX, []()
        {
            // Also test storing negative, but not -1 as it's not allowed by the HashSet
            int x = -1;
            while (x == -1) {
                x = Eigen::internal::random<int>(-(MAX_INDEX-1), MAX_INDEX - 1);
            }
            return x;
        })));
        CALL_SUBTEST((testLimitedSet<HashSet<MAX_INSERT_SIZE>>(MAX_INDEX, MAX_INSERT_SIZE)));
        CALL_SUBTEST((testLimitedQueue<MAX_INSERT_SIZE>(MAX_INDEX, MAX_INSERT_SIZE)));
    }
}
