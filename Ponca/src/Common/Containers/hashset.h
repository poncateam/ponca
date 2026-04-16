/**
This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
 \author Auberval Florian
*/

#pragma once

#include "../defines.h"

namespace Ponca {
    /*! \brief An HashMap implementation, that stores the internal values in an array of signed values.
     *
     * The internal values stored must not be equal to -1
     * For the search, it's best case complexity is O(1), and the worst case complexity is O(n).
     * \tparam N The maximum size of the HashSet
     * \tparam T The value type stored in the HashSet : Must be a signed value type that can be used as an array index (int-like)
     */
    template <int N, typename T>
    class HashSet {
    public:
        PONCA_MULTIARCH void clear();
        PONCA_MULTIARCH bool insert(int key);

    private:
        static constexpr T EMPTY = T(-1);
        T table[N] = {};

        PONCA_MULTIARCH [[nodiscard]] static int hash(const int x)
        {
            return (x * 2654435761u) % N;
        }
    };

    template <int N, typename T>
    PONCA_MULTIARCH void HashSet<N, T>::clear()
    {
        for (int i = 0; i < N; i ++)
            table[i] = EMPTY;
    }

    template <int N, typename T>
    PONCA_MULTIARCH bool HashSet<N, T>::insert(const int key)
    {
        const int h = hash(key);

        for (int i = 0; i < N; ++i) {
            int idx = (h + i) % N;

            int& slot = table[idx];

            if (slot == EMPTY) {
                slot = key;
                return true;
            }
            if (slot == key) {
                return false;
            }
        }
        return false; // full
    }
}