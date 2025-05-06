//
// Created by arkeo on 4/25/25.
//

#include <vector>
#include <functional>
#include <thread>
#include <iostream>
#include <numeric>
#include <algorithm>

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest.h>
#include <nanobench.h>

#include <GEL/Util/ParallelAdapters.h>

using GEL::Util::ThreadPool;
using GEL::Util::parallel_map;

bool test_empty()
{
    ThreadPool pool(2);
    const auto empty = parallel_map(pool, [](auto x) { return x; }, std::vector<int>{});
    return empty.empty();
}

static constexpr auto problem_size = 10000;

inline float non_trivial_function(const int x)
{
    float total = 0.0;
    for (int i = 0; i < 500; ++i) {
        total += 1.0f / static_cast<float>(i) * std::cos(static_cast<float>(x)/static_cast<float>(i+1));
    }
    return total;
}

bool test_oracle()
{
    constexpr auto size = 10;
    ThreadPool pool(3);
    std::vector<int> v(size);
    std::vector<float> vout(size);
    std::vector<float> vout2(size);
    std::iota(v.begin(), v.end(), 0);
    std::ranges::transform(v.cbegin(), v.cend(), vout.begin(), non_trivial_function);
    const auto right = parallel_map(pool, non_trivial_function, v, std::move(vout2));
    for (auto i = 0; i < 10; ++i) {
        std::cout << i << " : " << vout[i] << "\n";
        std::cout << i << " : " << right[i] << "\n";
        //CHECK_EQ(vout[i], right[i]);
    }
    CHECK_EQ(vout, right);
    return true;
}

TEST_CASE("parallel_map")
{
    CHECK(test_empty());
    test_oracle();
}



TEST_CASE("parallel_map_performance_singlethread")
{
    constexpr auto size = problem_size;
    ThreadPool pool(1);
    std::vector<int> v(size);
    std::vector<float> vout(size);
    std::iota(v.begin(), v.end(), 0);
    ankerl::nanobench::Bench().run("parallel_map_singlethreaded", [&]()
    {
        ankerl::nanobench::doNotOptimizeAway(
            parallel_map(pool, non_trivial_function, v, std::move(vout))
        );
    });
}

TEST_CASE("parallel_map_performance_2_threads")
{
    constexpr auto size = problem_size;
    ThreadPool pool(2);
    std::vector<int> v(size);
    std::vector<float> vout(size);
    std::iota(v.begin(), v.end(), 0);
    ankerl::nanobench::Bench().run("parallel_map_2_threads", [&]()
    {
        ankerl::nanobench::doNotOptimizeAway(
            parallel_map(pool, non_trivial_function, v, std::move(vout))
        );
    });
}

TEST_CASE("parallel_map_performance_4_threads")
{
    constexpr auto size = problem_size;
    ThreadPool pool(4);
    std::vector<int> v(size);
    std::vector<float> vout(size);
    std::iota(v.begin(), v.end(), 0);
    ankerl::nanobench::Bench().run("parallel_map_4_threads", [&]()
    {
        ankerl::nanobench::doNotOptimizeAway(
            parallel_map(pool, non_trivial_function, v, std::move(vout))
        );
    });
}

TEST_CASE("parallel_map_performance_8_threads")
{
    constexpr auto size = problem_size;
    ThreadPool pool(8);
    std::vector<int> v(size);
    std::vector<float> vout(size);
    std::iota(v.begin(), v.end(), 0);
    ankerl::nanobench::Bench().run("parallel_map_8_threads", [&]()
    {
        ankerl::nanobench::doNotOptimizeAway(
            parallel_map(pool, non_trivial_function, v, std::move(vout))
        );
    });
}

TEST_CASE("transform_performance")
{
    constexpr auto size = problem_size;
    std::vector<int> v(size);
    std::vector<float> vout(size);
    std::iota(v.begin(), v.end(), 0);
    ankerl::nanobench::Bench().run("std::transform", [&]()
    {
        ankerl::nanobench::doNotOptimizeAway(
            std::transform(v.cbegin(), v.cend(), vout.begin(), non_trivial_function)
        );
    });
}

TEST_CASE("ranges::transform_performance")
{
    constexpr auto size = problem_size;
    std::vector<int> v(size);
    std::vector<float> vout(size);
    std::iota(v.begin(), v.end(), 0);
    ankerl::nanobench::Bench().run("std::ranges::transform", [&]()
    {
        ankerl::nanobench::doNotOptimizeAway(
            std::ranges::transform(v.cbegin(), v.cend(), vout.begin(), non_trivial_function)
        );
    });
}

TEST_CASE("for_loop_performance")
{
    constexpr auto size = problem_size;
    std::vector<int> v(size);
    std::vector<float> vout(size);
    std::iota(v.begin(), v.end(), 0);
    auto for_loop = [&]()
    {
        for (auto i = 0; i < vout.size(); ++i)
        {
            vout[i] = non_trivial_function(i);
        }
        return 0;
    };

    ankerl::nanobench::Bench().run("for_loop", [&]()
    {
        ankerl::nanobench::doNotOptimizeAway(
            for_loop()
        );
    });
}

// int not_main()
// {
//
//     constexpr auto size = 10000;
//     std::vector<int> v(size);
//     for (int i = 0; i < size; ++i)
//     {
//         v[i] = i;
//     }
//
//     std::function<void(int,int)> lambda = [&](auto i, auto j)
//     {
//         v.at(i) = j;
//     };
//
//     auto lambda2 = [](int x) -> int
//     {
//         return x * 2 + 5;
//     };
//
//     //ParallelIterator<decltype(v)> pi(10);
//     ThreadPool pool(2);
//     std::vector<int> v2;
//     v2 = parallel_map(pool, lambda2, v, std::move(v2));
//
//     const auto v3 = parallel_map(pool, lambda2, v);
//     for (int i = 0; i < 50; ++i)
//     {
//         std::cout << i << " : " << v3.at(i) << "\n";
//     }
//     // std::vector<int> r2 = ParallelFor<int>(pool, v, [&](auto x) -> auto { return x; });
//     //pi.ForEach(v, lambda);
//     return 0;
// }