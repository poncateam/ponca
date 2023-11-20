#pragma once

#include <type_traits>

namespace Ponca 
{
namespace internal
{

//
// see https://stackoverflow.com/a/62292282/5317819
// work in c++17
// compile-time check if a class T has static data member M:
//
//     constexpr bool has = has_member_impl<T>( [](auto&& obj)->decltype( obj.M ){} );
//
template<typename T, typename F>
constexpr auto has_member_impl(F&& f) -> decltype(f(std::declval<T>()), true) {
    return true;
}
template<typename>
constexpr bool has_member_impl(...) {
    return false; 
}

} // namespace internal
} // namespace Ponca

#define IS_PROVIDED(CAPABILITY) \
    Ponca::internal::has_member_impl<Base>( [](auto&& obj)->decltype(obj.PROVIDES_ ## CAPABILITY){} )
#define PROVIDES(CAPABILITY) \
    static constexpr bool PROVIDES_ ## CAPABILITY{true}
#define REQUIRES(CAPABILITY) \
    static_assert(IS_PROVIDED(CAPABILITY), "Ponca capability '" #CAPABILITY "' not provided by any base class")
