/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/
#pragma once

#include <tuple>
#include <cstddef>
#include <algorithm>

namespace Ponca
{
    namespace internal
    {
        // Another namespace layer to avoid future potential conflicts in the name choices
        namespace Factory
        {
            /**
             * \brief StringLiteral as template parameters
             *
             * We use the C++20 NTTP extension to allow to template by string literals.
             * This class will be used to provide a name to FactoryEntries
             */
            template <size_t N>
            struct StringLiteral
            {
                char value[N];
                constexpr StringLiteral(const char (&str)[N]) { std::copy_n(str, N, value); }
            };
        } // namespace Factory
    } // namespace internal

    // We leave this in the Ponca namespace as it can be used by user

    /**
     * \brief Factory entry
     *
     * This class stores as compile time constant:
     * * A name of the entry (with string)
     * * A type, presumably a ComputeObject but this is not mandatory
     * * An identifier
     *
     * Any instance of this class stores an instance of the associated
     * type.
     */
    template <internal::Factory::StringLiteral _Name, class Type, unsigned int _MethodId = 0>
    struct FactoryEntry
    {
        using type = Type;

    private:
        static constexpr internal::Factory::StringLiteral _name = _Name;

    public:
        static constexpr size_t MethodId  = _MethodId;
        static constexpr const char* name = _name.value;

        FactoryEntry(size_t i) : idx(i) {}

        // This index (which is the index within the list) is not part of the type
        // because this would make the code horrendous. The factory specification
        // would need to wait for the list to be complete and hence store "to be templated classes".
        // Filling this index can be done in 3 lines and will work flawlessly without
        // any runtime cost (most of it will be at compile time) and the cost of an int in memory...
        const size_t idx;
        Type object;
    };

    namespace internal
    {
        // Another namespace layer to avoid future potential conflicts in the name choices
        namespace Factory
        {
            // Convert of std::tuple<FactoryEntry> to a std::tuple<FactoryEntry::type>
            template <typename EntryList>
            struct FactoryEntryTypeExtractor;

            // Specialized version for std::tuple
            template <typename... Ts>
            struct FactoryEntryTypeExtractor<std::tuple<Ts...>>
            {
                using type = std::tuple<typename Ts::type...>;
            };

            // Base filtering template interface
            template <class Pred, typename... Ts>
            struct filter;

            // Base case nothing to filter
            template <class Pred>
            struct filter<Pred>
            {
                using type = std::tuple<>;
            };

            // Shortcut to concatenate std::tuple
            template <typename... Tuples>
            using tuple_cat_t = decltype(std::tuple_cat(std::declval<Tuples>()...));

            // General case, head conditionnaly adds current type to the list;
            template <class Pred, typename T, typename... Ts>
            struct filter<Pred, T, Ts...>
            {
                static constexpr bool keepT = Pred::template value<T>;
                using head                  = std::conditional_t<keepT, std::tuple<T>, std::tuple<>>;
                using tail                  = typename filter<Pred, Ts...>::type;
                using type                  = tuple_cat_t<head, tail>;
            };

            // Type to allow type indirection. ie being able to extract all
            // types from a tuple while letting the class being templated
            // with a tuple
            template <class Pred, typename... Ts>
            struct filterList;

            template <typename Pred, typename... Ts>
            struct filterList<Pred, std::tuple<Ts...>>
            {
                using type = typename filter<Pred, Ts...>::type;
            };

            // Combination of multiple predicates
            template <template <typename> class... Preds>
            struct conjunction
            {
                template <typename T>
                static constexpr bool value = (Preds<T>::value && ...);
            };

        }; // namespace Factory
    } // namespace internal
} // namespace Ponca

