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
        PONCA_MULTIARCH typename Fit::VectorType project(const Fit& _f, const typename Fit::VectorType& _pos) const
        {
            return _f.project(_pos);
        }

        /**
         * \copydoc project
         */
        template <typename Fit>
        PONCA_MULTIARCH typename Fit::VectorType operator()(const Fit& _f, const typename Fit::VectorType& _pos) const
        {
            return project(_f, _pos);
        }
    };
    
    /**
     * \brief Project a point using the gradinet of the sdf
     * 
     * \tparam Fit The fitting object
     * 
     * \param f The fitting object
     * \param pos The position to project
     */
    struct GradientDescentProject
    {
    public:
        GradientDescentProject(unsigned int _nbIter = 16) : nbIter(_nbIter) {}

        /**
         * \brief Project a point using the gradinet of the sdf
         * 
         * \tparam Fit The fitting object
         * 
         * \param f The fitting object
         * \param pos The position to project
         */
        template<typename Fit>
        PONCA_MULTIARCH typename Fit::VectorType project(const Fit& _f, const typename Fit::VectorType& _pos) const
        {
            PONCA_MULTIARCH_STD_MATH(min)
            using VectorType = typename Fit::VectorType;
            using Scalar = typename Fit::Scalar;

            VectorType grad;
            VectorType dir  = _f.primitiveGradient(_pos);
            Scalar ilg      = Scalar(1.) / dir.norm();
            dir             = dir * ilg;
            Scalar ad       = _f.potential(_pos);
            Scalar delta    = -ad * min(ilg, Scalar(1.));
            VectorType proj = _pos + dir * delta;

            for (unsigned int i = 0; i < nbIter; ++i)
            {
                grad  = _f.primitiveGradient(proj);
                ilg   = Scalar(1.) / grad.norm();
                delta = -_f.potential(proj) * min(ilg, Scalar(1.));
                proj += dir * delta;
            }
            return proj;
        }

        /**
         * \copydoc project
         */
        template <typename Fit>
        PONCA_MULTIARCH typename Fit::VectorType operator()(const Fit& f, const typename Fit::VectorType& pos) const
        {
            return project(f, pos);
        }
    public:
        unsigned int nbIter;
    };
}; // namespace Ponca
