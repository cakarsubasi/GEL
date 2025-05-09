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
    using OutputOfN = std::invoke_result_t<F, In...>;

    template <typename F, typename... In>
    using Function = std::is_invocable<F, In...>;

    template <typename F, typename T, typename U>
    concept UnaryFunction =
        requires(T a, F f)
        {
            { std::is_invocable<F, decltype(a)>() };
            { f(a) } -> std::convertible_to<U>;
        };

    template <typename F, typename T, typename U, typename V>
    concept BinaryFunction =
        requires(T a, U b, F f)
        {
            { std::is_invocable<F, decltype(a), decltype(b)>() };
            { f(a, b) } -> std::convertible_to<V>;
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

/// @brief Map over an iterator.
/// @tparam F a unary function type of the form T -> U
/// @tparam InputIt An input iterator type that yields elements of T
/// @tparam OutputIt An output iterator type that takes elements of type U; uses std::vector<U> if unspecified
/// @tparam BoundsChecking Whether to perform bound checking in the iterators
/// @param pool The threadpool to use
/// @param f the function to map over
/// @param it the input iterator
/// @param out the output iterator; std::vector<U> by default
/// @return the forwarded output iterator
///
template <typename F,
          typename InputIt,
          typename OutputIt = std::vector<OutputOfN<F, typename InputIt::value_type>>,
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

/// @brief Map over an iterator with the indexes
/// @tparam F a binary function type of the form (size_t, T) -> U
/// @tparam InputIt An input iterator type that yields elements of T
/// @tparam OutputIt An output iterator type that takes elements of type U; uses std::vector<U> if unspecified
/// @tparam BoundsChecking Whether to perform bound checking in the iterators
/// @param pool The threadpool to use
/// @param f the function to map over
/// @param it the input iterator
/// @param out the output iterator; std::vector<U> by default
/// @return the forwarded output iterator.
///
template <typename F,
          typename InputIt,
          typename OutputIt = std::vector<OutputOfN<F, typename InputIt::value_type>>,
          bool BoundsChecking = false>
    requires
    ContiguousSizedCollection<InputIt> &&
    ContiguousSizedCollection<OutputIt> &&
    BinaryFunction<F, size_t, typename InputIt::value_type, typename std::remove_reference_t<OutputIt>::value_type>
auto parallel_enumerate_map(
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
                    out.at(j) = f(j, it.at(j));
                } else {
                    out[j] = f(j, it[j]);
                }
            }
        });
    }
    pool.waitAll();
    return std::forward<OutputIt>(out);
}

/// @brief Maps over two input iterators at the same time. If the iterators are of different length,
/// iterates up to the end of the shorter iterator.
/// @tparam F a unary function type of the form (T, U) -> V
/// @tparam InputIt1 An input iterator type that yields elements of T
/// @tparam InputIt2 An input iterator type that yields elements of U
/// @tparam OutputIt An output iterator type that takes elements of type V; uses std::vector<U> if unspecified
/// @tparam BoundsChecking Whether to perform bound checking in the iterators
/// @param pool The threadpool to use
/// @param f the function to map over
/// @param it1 the first input iterator
/// @param it2 the second input iterator
/// @param out the output iterator; std::vector<U> by default
/// @return the output iterator by value
///
template <typename F,
          typename InputIt1,
          typename InputIt2,
          typename OutputIt = std::vector<OutputOf<F, typename InputIt1::value_type>>,
          bool BoundsChecking = false>
    requires
    ContiguousSizedCollection<InputIt1> &&
    ContiguousSizedCollection<InputIt2> &&
    ContiguousSizedCollection<OutputIt> &&
    BinaryFunction<F, typename InputIt1::value_type, typename InputIt2::value_type, typename std::remove_reference_t<
                       OutputIt>::value_type>
auto parallel_map2(
    ThreadPool& pool,
    F&& f,
    const InputIt1& it1,
    const InputIt2& it2,
    OutputIt&& out = std::vector<typename OutputIt::value_type>()
) -> decltype(out)
{
    const auto pool_size = pool.size();
    const auto work_size = [&] {
        return it1.size() > it2.size() ? it2.size() : it1.size();
    }();
    const auto reduced_size = work_size / pool_size + (work_size % pool_size != 0);
    if (work_size == 0) {
        return std::forward<OutputIt>(out);
    }
    if (out.size() != work_size) {
        out.reserve(work_size);
        out.resize(work_size);
    }

    for (auto i = 0; i < pool_size; ++i) {
        pool.addTask([&out, i, reduced_size, work_size, &it1, &it2, &f] {
            auto max_size = std::min((i + 1) * reduced_size, work_size);
            for (auto j = i * reduced_size; j < max_size; ++j) {
                if constexpr (BoundsChecking) {
                    out.at(j) = f(it1.at(j), it2.at(j));
                } else {
                    out[j] = f(it1[j], it2[j]);
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
    UnaryFunction<F, typename InputIt::value_type, std::optional<typename std::remove_reference_t<
                      OutputIt>::value_type>>
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
                    if (const auto&& result = f(input_value); result.has_value()) {
                        auto&& inner = result.value();
                        out.at(begin + thread_counter) = inner;
                        thread_counter++;
                    }
                } else {
                    const auto& input_value = it[j];
                    if (const auto&& result = f(input_value); result.has_value()) {
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

/// @brief perform a map and a filter operation simultaneously over an enumerated iterator
/// @tparam F a binary function type of the form (size_t, T) -> std::optional<U>
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
          typename OutputIt = std::vector<typename OutputOfN<F, size_t, typename InputIt::value_type>::value_type>,
          bool BoundsChecking = false>
    requires
    ContiguousSizedCollection<InputIt> &&
    ContiguousSizedCollection<OutputIt>
// &&
    //BinaryFunction<F, size_t, typename InputIt::value_type, std::optional<typename std::remove_reference_t<
     //                 OutputIt>::value_type>>
auto parallel_enumerate_map_filter(
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
                    if (const auto&& result = f(j, input_value); result.has_value()) {
                        auto inner = std::move(result.value());
                        out.at(begin + thread_counter) = inner;
                        thread_counter++;
                    }
                } else {
                    const auto& input_value = it[j];
                    if (const auto&& result = f(j, input_value); result.has_value()) {
                        auto inner = std::move(result.value());
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

/// @brief perform a map and a filter operation simultaneously over two input iterators
/// @tparam F a binary function type of the form (T, U) -> std::optional<V>
/// @tparam InputIt1 An input iterator type that yields elements of T
/// @tparam InputIt2 An input iterator type that yields elements of U
/// @tparam OutputIt An output iterator type that takes elements of type V; uses std::vector<V> if unspecified
/// @tparam BoundsChecking Whether to perform bound checking in the iterators
/// @param pool The threadpool to use
/// @param f the function to map-filter over
/// @param it1 the first input iterator
/// @param it2 the second input iterator
/// @param out the output iterator; std::vector<U> by default
/// @return the output iterator by value
///
template <typename F,
          typename InputIt1,
          typename InputIt2,
          typename OutputIt = std::vector<typename OutputOfN<
              F, typename InputIt1::value_type, typename InputIt2::value_type>::value_type>,
          bool BoundsChecking = false>
    requires
    ContiguousSizedCollection<InputIt1> &&
    ContiguousSizedCollection<InputIt2> &&
    ContiguousSizedCollection<OutputIt> &&
    BinaryFunction<F, typename InputIt1::value_type, typename InputIt2::value_type, std::optional<typename
                       std::remove_reference_t<OutputIt>::value_type>>
auto parallel_map2_filter(
    ThreadPool& pool,
    F&& f,
    const InputIt1& it1,
    const InputIt2& it2,
    OutputIt&& out = std::vector<typename OutputIt::value_type>()
) -> decltype(out)
{
    // strategy: we allocate out to the same size as it, every thread has its own counter, and when they are all done,
    // we perform a serial memcpy. memcpy is memory-bound, meaning there is not a lot to gain from parallelization.
    const auto pool_size = pool.size();
    const auto work_size = [&] {
        return it1.size() > it2.size() ? it2.size() : it1.size();
    }();
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
        pool.addTask([&out, thread_id, reduced_size, work_size, &it1, &f, &thread_counter, &it2] {
            auto max_size = std::min((thread_id + 1) * reduced_size, work_size);
            auto begin = thread_id * reduced_size;
            for (auto j = begin; j < max_size; ++j) {
                if constexpr (BoundsChecking) {
                    const auto& input_value1 = it1.at(j);
                    const auto& input_value2 = it2.at(j);
                    const auto&& result = f(input_value1, input_value2);
                    if (result.has_value()) {
                        auto&& inner = result.value();
                        out.at(begin + thread_counter) = inner;
                        thread_counter++;
                    }
                } else {
                    const auto& input_value1 = it1[j];
                    const auto& input_value2 = it2[j];
                    const auto&& result = f(input_value1, input_value2);
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

/// @brief perform a map operation over two input iterators and a filter operation through two output iterators
/// @tparam F a binary function type of the form (T, U) -> std::optional<std::pair<V, W>>
/// @tparam InputIt1 An input iterator type that yields elements of T
/// @tparam OutputIt1 An output iterator type that takes elements of type U; uses std::vector<U> if unspecified
/// @tparam BoundsChecking Whether to perform bound checking in the iterators
/// @param pool The threadpool to use
/// @param f the function to map-filter over
/// @param it1 the first input iterator
/// @param it2 the second input iterator
/// @param out1 the first output iterator; std::vector<V> by default
/// @param out2 the second output iterator; std::vector<W> by default
/// @return the pair of output iterators by value
///
template <typename F,
          typename InputIt1,
          typename InputIt2,
          // TODO: the default types for these seem to be causing some problems
          typename OutputIt1 =
          std::vector<
              typename OutputOfN<F, typename InputIt1::value_type, typename InputIt2::value_type>
              ::value_type::first_type>,
          typename OutputIt2 =
          std::vector<
              typename OutputOfN<F, typename InputIt1::value_type, typename InputIt2::value_type>
              ::value_type::second_type>,
          bool BoundsChecking = false>
    requires
    ContiguousSizedCollection<InputIt1> &&
    ContiguousSizedCollection<InputIt2> &&
    ContiguousSizedCollection<OutputIt1> &&
    ContiguousSizedCollection<OutputIt2>
// TODO: fix this
// &&
//     BinaryFunction<F, typename InputIt1::value_type, typename InputIt2::value_type ,std::optional<
//      std::pair<typename std::remove_reference_t<OutputIt1>::value_type, typename std::remove_reference_t<OutputIt2>::value_type>>>
auto parallel_map2_filter2(
    ThreadPool& pool,
    F&& f,
    const InputIt1& it1,
    const InputIt2& it2,
    OutputIt1&& out1 = std::vector<typename OutputIt1::value_type>(),
    OutputIt2&& out2 = std::vector<typename OutputIt2::value_type>()
) -> std::pair<decltype(out1), decltype(out2)>
{
    // strategy: we allocate out to the same size as it, every thread has its own counter, and when they are all done,
    // we perform a serial memcpy. memcpy is memory-bound, meaning there is not a lot to gain from parallelization.
    const auto pool_size = pool.size();
    const auto work_size = [&] {
        return it1.size() > it2.size() ? it2.size() : it1.size();
    }();
    const auto reduced_size = work_size / pool_size + (work_size % pool_size != 0);
    if (work_size == 0) {
        return std::make_pair(std::forward<OutputIt1>(out1), std::forward<OutputIt2>(out2));
    }
    if (out1.size() != work_size) {
        out1.reserve(work_size);
        // safety post-condition: we need to manually resize out to the correct size at the end
        out1.resize(work_size);
    }
    if (out2.size() != work_size) {
        out2.reserve(work_size);
        // safety post-condition: we need to manually resize out to the correct size at the end
        out2.resize(work_size);
    }
    // TODO: move this to a scratch space so we amortize its creation
    std::vector<size_t> counters(pool_size, 0);
    for (auto thread_id = 0; thread_id < pool_size; ++thread_id) {
        auto& thread_counter = counters[thread_id];
        thread_counter = 0;
        pool.addTask([&out1, thread_id, reduced_size, work_size, &it1, &f, &thread_counter, &it2, out2] {
            auto max_size = std::min((thread_id + 1) * reduced_size, work_size);
            auto begin = thread_id * reduced_size;
            for (auto j = begin; j < max_size; ++j) {
                if constexpr (BoundsChecking) {
                    const auto& input_value1 = it1.at(j);
                    const auto& input_value2 = it2.at(j);
                    const auto&& result = f(input_value1, input_value2);
                    if (result.has_value()) {
                        auto [inner1, inner2] = result.value();
                        out1.at(begin + thread_counter) = std::move(inner1);
                        out2.at(begin + thread_counter) = std::move(inner2);
                        thread_counter++;
                    }
                } else {
                    const auto& input_value1 = it1[j];
                    const auto& input_value2 = it2[j];
                    const auto&& result = f(input_value1, input_value2);
                    if (result.has_value()) {
                        auto [inner1, inner2] = result.value();
                        out1[begin + thread_counter] = std::move(inner1);
                        out2[begin + thread_counter] = std::move(inner2);
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

    auto fix_chunks = [pool_size](auto& target, const auto& counter, auto max_chunk_size, auto total_size) {
        auto end_of_chunk1 = target.begin() + counter.at(0);
        for (auto chunk = 1; chunk < pool_size; ++chunk) {
            auto this_chunk_size = counter[chunk];
            auto chunk_begin = target.begin() + chunk * max_chunk_size;
            auto chunk_end = chunk_begin + this_chunk_size;
            std::copy(chunk_begin, chunk_end, end_of_chunk1);
            end_of_chunk1 += this_chunk_size;
        }
        target.resize(total_size);
    };
    fix_chunks(out1, counters, reduced_size, total_size);
    fix_chunks(out2, counters, reduced_size, total_size);

    return std::make_pair(std::forward<OutputIt1>(out1), std::forward<OutputIt2>(out2));
}

/// @brief perform a map operation over two input iterators and a filter operation through two output iterators
/// @tparam F a ternary function type of the form (size_t, T, U) -> std::optional<std::pair<V, W>>
/// @tparam InputIt1 An input iterator type that yields elements of T
/// @tparam OutputIt1 An output iterator type that takes elements of type U; uses std::vector<U> if unspecified
/// @tparam BoundsChecking Whether to perform bound checking in the iterators
/// @param pool The threadpool to use
/// @param f the function to map-filter over
/// @param it1 the first input iterator
/// @param it2 the second input iterator
/// @param out1 the first output iterator; std::vector<V> by default
/// @param out2 the second output iterator; std::vector<W> by default
/// @return the pair of output iterators by value
///
template <typename F,
          typename InputIt1,
          typename InputIt2,
          // TODO: the default types for these seem to be causing some problems
          typename OutputIt1 =
          std::vector<
              typename OutputOfN<F, typename InputIt1::value_type, typename InputIt2::value_type>
              ::value_type::first_type>,
          typename OutputIt2 =
          std::vector<
              typename OutputOfN<F, typename InputIt1::value_type, typename InputIt2::value_type>
              ::value_type::second_type>,
          bool BoundsChecking = false>
    requires
    ContiguousSizedCollection<InputIt1> &&
    ContiguousSizedCollection<InputIt2> &&
    ContiguousSizedCollection<OutputIt1> &&
    ContiguousSizedCollection<OutputIt2>
// TODO: fix this requirement
// &&
//     BinaryFunction<F, typename InputIt1::value_type, typename InputIt2::value_type ,std::optional<
//      std::pair<typename std::remove_reference_t<OutputIt1>::value_type, typename std::remove_reference_t<OutputIt2>::value_type>>>
auto parallel_enumerate_map2_filter2(
    ThreadPool& pool,
    F&& f,
    const InputIt1& it1,
    const InputIt2& it2,
    OutputIt1&& out1 = std::vector<typename OutputIt1::value_type>(),
    OutputIt2&& out2 = std::vector<typename OutputIt2::value_type>()
)
//-> std::pair<decltype(out1), decltype(out2)>
{
    // strategy: we allocate out to the same size as it, every thread has its own counter, and when they are all done,
    // we perform a serial memcpy. memcpy is memory-bound, meaning there is not a lot to gain from parallelization.
    const auto pool_size = pool.size();
    const auto work_size = [&] {
        return it1.size() > it2.size() ? it2.size() : it1.size();
    }();
    const auto reduced_size = work_size / pool_size + (work_size % pool_size != 0);
    if (work_size == 0) {
    //    return std::make_pair(std::forward<OutputIt1&&>(out1), std::forward<OutputIt2&&>(out2));
        return;
    }
    if (out1.size() != work_size) {
        out1.reserve(work_size);
        // safety post-condition: we need to manually resize out to the correct size at the end
        out1.resize(work_size);
    }
    if (out2.size() != work_size) {
        out2.reserve(work_size);
        // safety post-condition: we need to manually resize out to the correct size at the end
        out2.resize(work_size);
    }
    // TODO: move this to a scratch space so we amortize its creation
    std::vector<size_t> counters(pool_size, 0);
    for (auto thread_id = 0; thread_id < pool_size; ++thread_id) {
        auto& thread_counter = counters[thread_id];
        thread_counter = 0;
        pool.addTask([&out1, thread_id, reduced_size, work_size, &it1, &f, &thread_counter, &it2, &out2] {
            auto max_size = std::min((thread_id + 1) * reduced_size, work_size);
            auto begin = thread_id * reduced_size;
            for (auto j = begin; j < max_size; ++j) {
                if constexpr (BoundsChecking) {
                    const auto& input_value1 = it1.at(j);
                    const auto& input_value2 = it2.at(j);
                    const auto result = f(j, input_value1, input_value2);
                    if (result.has_value()) {
                        auto [inner1, inner2] = result.value();
                        out1.at(begin + thread_counter) = std::move(inner1);
                        out2.at(begin + thread_counter) = std::move(inner2);
                        thread_counter++;
                    }
                } else {
                    const auto& input_value1 = it1[j];
                    const auto& input_value2 = it2[j];
                    auto result = f(j, input_value1, input_value2);
                    if (result.has_value()) {
                        auto [inner1, inner2] = std::move(result.value());
                        out1[begin + thread_counter] = inner1;
                        out2[begin + thread_counter] = inner2;
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

    auto fix_chunks = [pool_size](auto& target, const auto& counter, auto max_chunk_size, auto total_size) {
        auto end_of_chunk1 = target.begin() + counter.at(0);
        for (auto chunk = 1; chunk < pool_size; ++chunk) {
            auto this_chunk_size = counter[chunk];
            auto chunk_begin = target.begin() + chunk * max_chunk_size;
            auto chunk_end = chunk_begin + this_chunk_size;
            std::copy(chunk_begin, chunk_end, end_of_chunk1);
            end_of_chunk1 += this_chunk_size;
        }
        target.resize(total_size);
    };
    fix_chunks(out1, counters, reduced_size, total_size);
    fix_chunks(out2, counters, reduced_size, total_size);

    // return std::make_pair(std::forward<OutputIt1>(out1), std::forward<OutputIt2>(out2));
}


} // namespace GEL::Util

#endif //GEL_UTIL_PARALLELADAPTERS_H
