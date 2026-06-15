/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/
#pragma once

#include "factorylist.h"
#include "typelist.h"

namespace Ponca
{
    namespace internal
    {
        /**
         * \brief List of compute object that can be iterated over
         *
         * Unlike a Factory, this object has non static methods and
         * stores a copy of object.
         *
         * \tparam _FactoryEntries A std::tuple of FactoryEntry
         */
        template <typename _FactoryEntries>
        struct ComputeObjectList
        {
            using FactoryEntries = _FactoryEntries;

            /**
             * \brief Returns a new object that is a filtered version of the current one
             *
             * \tparam Preds List of predicate to filter by
             */
            template <template <typename> class... Preds>
            static decltype(auto) Filter()
            {
                using namespace Factory;
                return ComputeObjectList<typename filterList<conjunction<Preds...>, FactoryEntries>::type>{};
            };

            /**
             * \brief Returns a new object that is a filtered version of the current one
             *
             * \tparam MethodIds List of ids to filter by methods
             */
            template <Method... MethodIds>
            static decltype(auto) Filter()
            {
                return Filter<MethodProvider<(unsigned int)MethodIds>::template pred...>();
            }

            /**
             * \brief Return an object by it's id
             *
             * In order for this method to work, it is expected that only one method
             * has the specified id !
             *
             * This method returns the ComputeObject, not a ComputeObjectList !
             *
             * \tparam id The id of the method
             */
            template <Method id>
            static decltype(auto) GetMethod()
            {
                using MethodList = decltype(Filter<MethodProvider<(unsigned int)id>::template pred>());

                static_assert(std::tuple_size_v<typename MethodList::FactoryEntries> == 1);
                return std::get<0>(typename MethodList::FactoryEntries{}).object;
            }

            /**
             * \brief Applies a function to each entry of the list
             *
             * This is expected to be a template function. Example
             * usage is:
             * ComputeObjectList.foreach([&](auto& x)
             * {
             *      // Code here should be valid for all possible
             *      // types !
             *
             *      // If you know the type, you may run specialized code
             *      using CurrentType = std::decay_t<decltype(x)>;
             *      if constexpr (std::is_same_v<CurrentType, TheExpectedType>)
             *      {
             *          // The code goes here.
             *      }
             * });
             *
             * \tparam Func The function to apply to each element
             */
            template <typename Func>
            void foreach (Func&& fun)
            {
                std::apply([&](auto&&... xs) { (fun(xs), ...); }, m_entries);
            }

        private:
            FactoryEntries m_entries;
        };

        /**
         * \brief Base class for factory
         *
         * This class behaves similarly to ComputeObjectList but do not
         * store any information. It is here to have plain static methods.
         *
         * \tparam _FactoryEntries A std::tuple<FactoryEntry>
         */
        template <typename _FactoryEntries>
        struct FactoryBase
        {
            using FactoryEntries        = _FactoryEntries;
            using FullComputeObjectList = ComputeObjectList<FactoryEntries>;

            /**
             * \copydoc ComputeObjectList::Filter
             */
            template <template <typename> class... Preds>
            static decltype(auto) Filter()
            {
                return FullComputeObjectList::template Filter<Preds...>();
            }

            /**
             * \copydoc ComputeObjectList::Filter
             */
            template <Method... MethodIds>
            static decltype(auto) Filter()
            {
                return FullComputeObjectList::template Filter<MethodIds...>();
            }

            /**
             * \copydoc ComputeObjectList::GetMethod
             */
            template <Method id>
            static decltype(auto) GetMethod()
            {
                return FullComputeObjectList::template GetMethod<id>();
            }

            /**
             * \copydoc ComputeObjectList::Foreach
             */
            template <typename Func>
            static void foreach (Func&& fun)
            {
                FullComputeObjectList{}.foreach (fun);
            }
        };
    } // namespace internal

    /**
     * \brief Default factory templated on Point, Filter and Derivative
     *
     * \tparam P Point
     * \tparam NF NeighborFilter
     * \tparam DerType Derivative Type
     */
    template <typename P, typename NF, int DerType>
    struct Factory : public internal::FactoryBase<FactoryGenericList<P, NF, DerType>>
    {
    };

    /**
     * \brief Default factory templated on Point, Filter and Derivative
     *
     * This version is specialized for 3D points
     *
     * \tparam P Point
     * \tparam NF NeighborFilter
     * \tparam DerType Derivative Type
     */
    template <typename P, typename NF, int DerType>
        requires Is3D<P>
    struct Factory<P, NF, DerType>
        : public internal::FactoryBase<
              internal::Factory::tuple_cat_t<FactoryGenericList<P, NF, DerType>, Factory3DList<P, NF, DerType>>>
    {
    };
} // namespace Ponca
