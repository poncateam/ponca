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
 * \brief Test the validity of Ponca Common Set of int classes by randomly inserting elements
 * and comparing the results with standard std::set<int>
 *
 * \tparam IndexSet The Set of indices to evaluate
 * \param _maxIndex The maximum possible index (must be > 100)
 */
template <typename IndexSet>
void testSetStandardFeatures(const int _maxIndex)
{
    // assert(_maxIndex >= 0);
    IndexSet indexSet;

    // Verify that the set starts empty
    for (int i = 0; i < (QUICK_TESTS ? 1 : _maxIndex); ++i)
    {
        VERIFY(!indexSet.contains(i));
    }

    std::set<int> indexSetSTD; // Stores id
    IndexSet indexSetPonca;

    // Pick a number of elements to insert in the set
    const int nbInsertion = QUICK_TESTS ? 1 : Eigen::internal::random<int>(100, _maxIndex-1);

    for (int i = 0; i < nbInsertion; ++i)
    {
        // Select a random index to add to the set
        const int idx = Eigen::internal::random<int>(0, _maxIndex-1);
        // Check before insert
        const bool valFound1 = indexSetSTD.contains(idx);
        const bool valFound2 = indexSetPonca.contains(idx);

        VERIFY((valFound1 == valFound2));
        // VERIFY((indexSetSTD.contains(idx) == indexSetPonca.contains(idx)));

        // VERIFY((flagSTD.insert(idx).second == flagPonca.insert(idx)));
        auto val1 = indexSetSTD.insert(idx).second;
        auto val2 = indexSetPonca.insert(idx);
        VERIFY((val1 == val2)); // Check insert result

        // Check after insert
        VERIFY((indexSetSTD.contains(idx) == indexSetPonca.contains(idx)));
    }

    for (int i = 0; i < nbInsertion; ++i)
    {
        VERIFY((indexSetSTD.contains(i) == indexSetPonca.contains(i)));
    }
}

/*
 * \brief Test the Ponca limited indice set capabilities
 *
 * \tparam IndexSet A Set of indices limited in size to evaluate
 * \param _maxIndex The maximum possible index
 * \param _setCapacity The maximum capacity of the limited set
 */
template <typename IndexSet>
void testLimitedSet(const int _maxIndex, const int _setCapacity)
{
    IndexSet indexSet;
    // Insert until we reach max capacity
    for (int i = 0; i < _setCapacity; ++i)
    {
        const int idx = Eigen::internal::random<int>(0, _maxIndex-1);
        VERIFY((indexSet.insert(idx) == true));
    }

    // Select a number of insertion that exceed available total capacity
    const int nbInsertion = QUICK_TESTS ? 1 :  Eigen::internal::random<int>(_setCapacity+1, _setCapacity+10000);
    for (int i = _setCapacity; i < nbInsertion; ++i)
    {
        const int idx = Eigen::internal::random<int>(0, _maxIndex-1);
        VERIFY((indexSet.insert(idx) == false));
    }
}

int main(const int argc, char** argv)
{
    if (!init_testing(argc, argv))
    {
        return EXIT_FAILURE;
    }

    constexpr int MAX_INDEX       = 10000;     //< The upper limit of index range (index can't go higher than this number for BITSET)
    constexpr int MAX_INSERT_SIZE = MAX_INDEX; //< The maximum size of a limited set (values won't be inserted after that)

    cout << "Check consistency of Set of Indices class" << endl;

    for (int i = 0; i < g_repeat; ++i)
    {
        CALL_SUBTEST(( testSetStandardFeatures<BitSet<MAX_INDEX>>(MAX_INDEX)));
        CALL_SUBTEST(( testSetStandardFeatures<HashSet<MAX_INSERT_SIZE>>(MAX_INDEX) ));
        CALL_SUBTEST(( testLimitedSet<HashSet<MAX_INSERT_SIZE>>(MAX_INDEX, MAX_INSERT_SIZE) ));
    }
}
