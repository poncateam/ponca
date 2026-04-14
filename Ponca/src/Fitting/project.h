/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/
#pragma once

namespace Ponca
{
    /**
     * \brief Project a point using the primitive projection operator
     */
    struct SimpleProject
    {
    public:
        /**
         * \brief Project a point using the primitive projection operator
         * 
         * \tparam Fit The fitting object
         * 
         * \param f The fitting object
         * \param pos The position to project
         */
        template<typename Fit>
        typename Fit::VectorType project(const Fit& f, const typename Fit::VectorType& pos) const
        {
            return f.project(pos);
        }

        /**
         * \copydoc project
         */
        template <typename Fit>
        typename Fit::VectorType operator()(const Fit& f, const typename Fit::VectorType& pos) const
        {
            return project(f, pos);
        }

    private:
    };
}; // namespace Ponca
