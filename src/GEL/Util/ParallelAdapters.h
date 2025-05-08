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
    concept ContiguousSizedCollection = requires(const C& c, size_t i)
    {
        { c.size() } -> std::convertible_to<size_t>;
        { c.begin() } -> std::random_access_iterator;
    };

    template <typename F, typename In>
    using OutputOf = std::invoke_result_t<F, In>;

    template <typename F, typename... In>
    using Function = std::is_invocable<F, In...>;

    template <typename F, typename T, typename U>
    concept UnaryFunction =
        requires(T a, F f)
        {
            { std::is_invocable<F, decltype(a)>() };
            { f(a) } -> std::convertible_to<U>;
        };

    static constexpr auto identity = []<typename T0>(T0 x) -> T0 {
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
requires
ContiguousSizedCollection<InputIt> &&
UnaryFunction<F, typename InputIt::value_type, void>
auto parallel_foreach(ThreadPool& pool, F&& f, const InputIt& it) -> void
{
    const auto pool_size = pool.size();
    const auto work_size = it.size();
    const auto reduced_size = work_size / pool_size + (work_size % pool_size != 0);
    if (work_size == 0) {
        return;
    }
    for (auto i = 0; i < pool_size; ++i) {
        pool.addTask([i, reduced_size, work_size, &it, &f] {
            const auto max_size = std::min((i + 1) * reduced_size, work_size);
            for (auto j = i * reduced_size; j < max_size; ++j) {
                if constexpr (BoundsChecking) {
                    f(it.at(j));
                } else {
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
    requires
    ContiguousSizedCollection<InputIt> &&
    ContiguousSizedCollection<OutputIt> &&
    UnaryFunction<F, typename InputIt::value_type, typename std::remove_reference_t<OutputIt>::value_type>
auto parallel_map(
    ThreadPool& pool,
    F&& f,
    const InputIt& it,
    OutputIt&& out = std::vector<typename OutputIt::value_type>()
) -> decltype(out)
{
    const auto pool_size = pool.size();
    const auto work_size = it.size();
    const auto reduced_size = work_size / pool_size + (work_size % pool_size != 0);
    if (work_size == 0) {
        return std::forward<OutputIt>(out);
    }
    if (out.size() != work_size) {
        out.reserve(work_size);
        out.resize(work_size);
    }

    for (auto i = 0; i < pool_size; ++i) {
        pool.addTask([&out, i, reduced_size, work_size, &it, &f] {
            auto max_size = std::min((i + 1) * reduced_size, work_size);
            for (auto j = i * reduced_size; j < max_size; ++j) {
                if constexpr (BoundsChecking) {
                    out.at(j) = f(it.at(j));
                } else {
                    out[j] = f(it[j]);
                }
            }
        });
    }
    pool.waitAll();
    return std::forward<OutputIt>(out);
}

/// @brief perform a filter
/// @tparam F a unary function predicate of the form T -> bool
/// @tparam InputIt An input iterator type that yields elements of T
/// @tparam OutputIt An output iterator type that takes elements of T
/// @tparam BoundsChecking Whether to perform bound checking in the iterators
/// @param pool The threadpool to use
/// @param p the predicate
/// @param it the input iterator
/// @param out the output iterator; std::vector<U> by default
/// @return the output iterator by value
///
template <typename F,
          typename InputIt,
          typename OutputIt = std::vector<typename InputIt::value_type>,
          bool BoundsChecking = false>
    requires
    ContiguousSizedCollection<InputIt> &&
    ContiguousSizedCollection<OutputIt> &&
    UnaryFunction<F, typename InputIt::value_type, bool>
auto parallel_filter(
    ThreadPool& pool,
    F&& p,
    const InputIt& it,
    OutputIt&& out = std::vector<typename OutputIt::value_type>()
) -> decltype(out)
{
    // strategy: we allocate out to the same size as it, every thread has its own counter, and when they are all done,
    // we perform a serial memcpy. memcpy is memory-bound, meaning there is not a lot to gain from parallelization.
    const auto pool_size = pool.size();
    const auto work_size = it.size();
    const auto reduced_size = work_size / pool_size + (work_size % pool_size != 0);
    if (work_size == 0) {
        return std::forward<OutputIt>(out);
    }
    if (out.size() != work_size) {
        out.reserve(work_size);
        // safety post-condition: we need to manually resize out to the correct size at the end
        out.resize(work_size);
    }
    // TODO: move this to a scratch space so we amortize its creation
    std::vector<size_t> counters(pool_size, 0);
    for (auto thread_id = 0; thread_id < pool_size; ++thread_id) {
        auto& thread_counter = counters[thread_id];
        thread_counter = 0;
        pool.addTask([&out, thread_id, reduced_size, work_size, &it, &p, &thread_counter] {
            auto max_size = std::min((thread_id + 1) * reduced_size, work_size);
            auto begin = thread_id * reduced_size;
            for (auto j = begin; j < max_size; ++j) {
                if constexpr (BoundsChecking) {
                    const auto& value = it.at(j);
                    if (p(value)) {
                        out.at(begin + thread_counter) = value;
                        thread_counter++;
                    }
                } else {
                    const auto& value = it[j];
                    if (p(value)) {
                        out[begin + thread_counter] = value;
                        thread_counter++;
                    }
                }
            }
        });
    }
    pool.waitAll();

    decltype(counters)::value_type total_size = 0;
    for (const auto& counter : counters) {
        total_size += counter;
    }

    auto end_of_chunk = out.begin() + counters.at(0);
    for (auto chunk = 1; chunk < pool_size; ++chunk) {
        auto chunk_size = counters[chunk];
        auto chunk_begin = out.begin() + chunk * reduced_size;
        auto chunk_end = chunk_begin + chunk_size;
        std::copy(chunk_begin, chunk_end, end_of_chunk);
        end_of_chunk += chunk_size;
    }
    out.resize(total_size);

    return std::forward<OutputIt>(out);
}

/// @brief perform a map and a filter operation simultaneously
/// @tparam F a unary function type of the form T -> std::optional<U>
/// @tparam InputIt An input iterator type that yields elements of T
/// @tparam OutputIt An output iterator type that takes elements of type U; uses std::vector<U> if unspecified
/// @tparam BoundsChecking Whether to perform bound checking in the iterators
/// @param pool The threadpool to use
/// @param f the function to map-filter over
/// @param it the input iterator
/// @param out the output iterator; std::vector<U> by default
/// @return the output iterator by value
///
template <typename F,
          typename InputIt,
          typename OutputIt = std::vector<typename OutputOf<F, typename InputIt::value_type>::value_type>,
          bool BoundsChecking = false>
    requires
    ContiguousSizedCollection<InputIt> &&
    ContiguousSizedCollection<OutputIt> &&
    UnaryFunction<F, typename InputIt::value_type, std::optional<typename std::remove_reference_t<OutputIt>::value_type>>
auto parallel_map_filter(
    ThreadPool& pool,
    F&& f,
    const InputIt& it,
    OutputIt&& out = std::vector<typename OutputIt::value_type>()
) -> decltype(out)
{
    // strategy: we allocate out to the same size as it, every thread has its own counter, and when they are all done,
    // we perform a serial memcpy. memcpy is memory-bound, meaning there is not a lot to gain from parallelization.
    const auto pool_size = pool.size();
    const auto work_size = it.size();
    const auto reduced_size = work_size / pool_size + (work_size % pool_size != 0);
    if (work_size == 0) {
        return std::forward<OutputIt>(out);
    }
    if (out.size() != work_size) {
        out.reserve(work_size);
        // safety post-condition: we need to manually resize out to the correct size at the end
        out.resize(work_size);
    }
    // TODO: move this to a scratch space so we amortize its creation
    std::vector<size_t> counters(pool_size, 0);
    for (auto thread_id = 0; thread_id < pool_size; ++thread_id) {
        auto& thread_counter = counters[thread_id];
        thread_counter = 0;
        pool.addTask([&out, thread_id, reduced_size, work_size, &it, &f, &thread_counter] {
            auto max_size = std::min((thread_id + 1) * reduced_size, work_size);
            auto begin = thread_id * reduced_size;
            for (auto j = begin; j < max_size; ++j) {
                if constexpr (BoundsChecking) {
                    const auto& input_value = it.at(j);
                    const auto&& result = f(input_value);
                    if (result.has_value()) {
                        auto&& inner = result.value();
                        out.at(begin + thread_counter) = inner;
                        thread_counter++;
                    }
                } else {
                    const auto& input_value = it[j];
                    const auto&& result = f(input_value);
                    if (result.has_value()) {
                        auto&& inner = result.value();
                        out[begin + thread_counter] = inner;
                        thread_counter++;
                    }
                }
            }
        });
    }
    pool.waitAll();

    decltype(counters)::value_type total_size = 0;
    for (const auto& counter : counters) {
        total_size += counter;
    }

    auto end_of_chunk = out.begin() + counters.at(0);
    for (auto chunk = 1; chunk < pool_size; ++chunk) {
        auto chunk_size = counters[chunk];
        auto chunk_begin = out.begin() + chunk * reduced_size;
        auto chunk_end = chunk_begin + chunk_size;
        std::copy(chunk_begin, chunk_end, end_of_chunk);
        end_of_chunk += chunk_size;
    }
    out.resize(total_size);

    return std::forward<OutputIt>(out);
}
} // namespace GEL::Util


#endif //GEL_UTIL_PARALLELADAPTERS_H
