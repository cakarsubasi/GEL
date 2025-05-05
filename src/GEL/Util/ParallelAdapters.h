//
// Created by Cem Akarsubasi on 5/5/25.
//
// TODO: examples
// TODO: map2, map3
// TODO: fold

#ifndef GEL_UTIL_PARALLELADAPTERS_H
#define GEL_UTIL_PARALLELADAPTERS_H

#include <GEL/Util/ThreadPool.h>

namespace GEL::Util
{

/// Concepts for constraining the inputs of parallel adapters
namespace Concepts
{
    template <typename C>
    concept Collection = requires(const C& c, size_t i)
    {
        { c.size() } -> std::convertible_to<size_t>;
        { c.at(i) } -> std::same_as<typename C::value_type const&>;
        { c[i] } -> std::same_as<typename C::value_type const&>;
    };

    template <typename F, typename In>
    using OutputOf = std::invoke_result_t<F, In>;

    template <typename F, typename... In>
    using Function = std::is_invocable<F, In...>;

    template <typename F, typename T, typename U>
    concept UnaryFunction =
        requires(T a, F f)
        {
            { std::is_invocable<F, T>() }; //hello
            { f(a) } -> std::same_as<U>; // hello
        };

    static constexpr auto identity = []<typename T0>(T0 x) -> T0
    {
        return x;
    };
    static_assert(UnaryFunction<decltype(identity), int, int>);
    static_assert(UnaryFunction<decltype(identity), float, float>);
} // namespace Concepts

using namespace Concepts;

///
/// @tparam F a unary function type of the form T -> void
/// @tparam InputIt An input iterator type that yields elements of T
/// @tparam BoundsChecking Whether to perform bound checking in the iterators
/// @param pool The threadpool to use
/// @param f the function to map over
/// @param it the input iterator
///
template <typename F,
          typename InputIt,
          bool BoundsChecking = true>
auto parallel_foreach(ThreadPool& pool, F f, const InputIt& it) -> void
{
    const auto pool_size = pool.size();
    const auto work_size = it.size();
    const auto reduced_size = work_size / pool_size;

    for (auto i = 0; i < pool_size; ++i)
    {
        pool.addTask([i, reduced_size, work_size, &it, &f]
        {
            const auto max_size = std::min((i + 1) * reduced_size, work_size);
            for (auto j = i * reduced_size; j < max_size; ++j)
            {
                if constexpr (BoundsChecking)
                {
                    f(it.at(j));
                }
                else
                {
                    f(it[j]);
                }
            }
        });
    }
    pool.waitAll();
}

///
/// @tparam F a unary function type of the form T -> U
/// @tparam InputIt An input iterator type that yields elements of T
/// @tparam OutputIt An output iterator type that takes elements of type U; uses std::vector<U> if unspecified
/// @tparam BoundsChecking Whether to perform bound checking in the iterators
/// @param pool The threadpool to use
/// @param f the function to map over
/// @param it the input iterator
/// @param out the output iterator; std::vector<U> by default
/// @return the output iterator by value
///
template <typename F,
          typename InputIt,
          typename OutputIt = std::vector<OutputOf<F, typename InputIt::value_type>>,
          bool BoundsChecking = false>
    requires Collection<InputIt> &&
    Collection<OutputIt> &&
    UnaryFunction<F, typename InputIt::value_type, typename OutputIt::value_type>
auto parallel_map(
    ThreadPool& pool,
    const F& f,
    const InputIt& it,
    OutputIt&& out = std::vector<typename OutputIt::value_type>()
) -> std::remove_reference_t<decltype(out)>
{
    const auto pool_size = pool.size();
    const auto work_size = it.size();
    const auto reduced_size = work_size / pool_size;
    if (out.size() != work_size)
    {
        out.reserve(work_size);
        out.resize(work_size);
    }

    for (auto i = 0; i < pool_size; ++i)
    {
        pool.addTask([&out, i, reduced_size, work_size, &it, &f]
        {
            auto max_size = std::min((i + 1) * reduced_size, work_size);
            for (auto j = i * reduced_size; j < max_size; ++j)
            {
                if constexpr (BoundsChecking)
                {
                    out.at(j) = f(it.at(j));
                }
                else
                {
                    out[j] = f(it[j]);
                }
            }
        });
    }
    pool.waitAll();
    return out;
}

} // namespace GEL::Util


#endif //GEL_UTIL_PARALLELADAPTERS_H
