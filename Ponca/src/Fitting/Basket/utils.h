/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once
#include <iostream>

namespace Ponca
{
    /**
     * \brief Helper class to check for status of required extensions
     *
     * It allows to cache whether the last call to finalize was
     * FIT_RESULT::STABLE or not.
     *
     * We encourage replacing `PONCA_FITTING_DECLARE_DEFAULT_TYPES`
     * with `PONCA_FITTING_DECLARE_MULTIPASS_TYPES`

     * And use the class as follows in addLocalNeighbor / finalize
     * \code
     *  if (!Status::ready()) return Base::finalize()
     *  else // Insert the code here
     * \endcode
     */
    template <class DataPoint, class _NFilter, typename T>
    class MultipassStatus : public T
    {
        PONCA_FITTING_DECLARE_DEFAULT_TYPES
    public:
        MultipassStatus() : m_ready(false) {}

        PONCA_MULTIARCH inline void init()
        {
            Base::init();
            m_ready = false;
        }

        PONCA_MULTIARCH [[nodiscard]] inline FIT_RESULT finalize()
        {
            FIT_RESULT res = Base::finalize();
            m_ready        = (res == STABLE);
            res            = (res == STABLE) ? NEED_OTHER_PASS : res;

            return Base::m_eCurrentState = res;
        }

        PONCA_MULTIARCH inline bool ready() const { return m_ready; }

    protected:
        bool m_ready;
    };
} // namespace Ponca
