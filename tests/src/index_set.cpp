/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

/*!
 * \file tests/src/set.cpp
 * \brief Validate custom HashSet and BitSet
 */

#include "../common/testing.h"
#include "../common/testUtils.h"
#include <Ponca/src/Common/Containers/bitset.h>
#include <Ponca/src/Common/Containers/hashset.h>

#include <random>
#include <set>
#include <vector>

using namespace Ponca;
using namespace std;

/*!
 * \brief Test the validity of Ponca Common Set of int classes by randomly inserting values
 * and comparing the results with standard std::set<int>
 *
 * \tparam IndexSet The Set of indices to evaluate
 * \tparam RandomFunctor The Random functor type
 * \param _maxIndex The maximum possible index (must be > 100)
 * \param _pickRandom The function to call when generating a random value
 */
template <typename IndexSet, typename RandomFunctor>
void testSetStandardFeatures(const int _maxIndex, RandomFunctor _pickRandom)
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
void testSetStandardFeatures(const int _maxIndex)
{
    testSetStandardFeatures<IndexSet>(_maxIndex, [&_maxIndex]()
    {
        return Eigen::internal::random<int>(0, _maxIndex - 1);
    });
}

/*
 * \brief Test the insert capabilities of the limited Set
 *
 * \tparam IndexSet A Set of indices limited in size to evaluate
 * \param _maxIndex The maximum possible index
 * \param _setCapacity The maximum capacity of the limited set
 */
template <typename IndexSet>
void testLimitedSet(const int _maxIndex, const int _setCapacity)
{
    IndexSet indexSet;

    // Select a number of insertion that exceed available total capacity
    const int nbTotalInsertion = QUICK_TESTS ? 1 : Eigen::internal::random<int>(_setCapacity + 1, _maxIndex);

    // To keep track that we insert each index only once
    vector<int> indices(nbTotalInsertion);
    std::iota(indices.begin(), indices.end(), 1);
    ranges::shuffle(indices, std::mt19937{std::random_device{}()});

    // Insert until we reach max capacity
    for (int i = 0; i < _setCapacity; ++i)
        VERIFY((indexSet.insert(indices[i]) == true));

    // Select a number of insertion that exceed available total capacity
    for (int i = _setCapacity; i < nbTotalInsertion; ++i)
        VERIFY((indexSet.insert(indices[i]) == false));
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
    constexpr int MAX_INSERT_SIZE = MAX_INDEX;

    cout << "Check consistency of Set of Indices class" << endl;

    for (int i = 0; i < g_repeat; ++i)
    {
        CALL_SUBTEST((testSetStandardFeatures<BitSet<MAX_INDEX>>(MAX_INDEX)));
        CALL_SUBTEST((testSetStandardFeatures<HashSet<MAX_INSERT_SIZE>>(MAX_INDEX, []()
        {
            // Also test storing negative, but not -1 as it's not allowed by the HashSet
            int x = -1;
            while (x == -1) {
                x = Eigen::internal::random<int>(-(MAX_INDEX-1), MAX_INDEX - 1);
            }
            return x;
        })));
        CALL_SUBTEST((testLimitedSet<HashSet<MAX_INSERT_SIZE>>(MAX_INDEX, MAX_INSERT_SIZE)));
    }
}
