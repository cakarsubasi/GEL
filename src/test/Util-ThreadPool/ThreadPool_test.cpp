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

TEST_CASE("parallel_map")
{
    CHECK(test_empty());
}

static constexpr auto problem_size = 10000000;

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
            parallel_map(pool, [=](const int x){ return std::cos(static_cast<float>(x)/static_cast<float>(size)); }, v, std::move(vout))
        );
    });
}

auto f(const int x) { return std::cos(static_cast<float>(x)/static_cast<float>(problem_size)); };

TEST_CASE("parallel_map_performance_singlethread_static")
{
    constexpr auto size = problem_size;

    ThreadPool pool(1);
    std::vector<int> v(size);
    std::vector<float> vout(size);
    std::iota(v.begin(), v.end(), 0);
    ankerl::nanobench::Bench().run("parallel_map_singlethreaded_static", [&]()
    {
        ankerl::nanobench::doNotOptimizeAway(
            parallel_map(pool, f, v, std::move(vout))
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
            std::transform(v.cbegin(), v.cend(), vout.begin(), [=](const int x) {return std::cos(static_cast<float>(x)/static_cast<float>(size));})
        );
    });
}

TEST_CASE("transform_performance")
{
    constexpr auto size = problem_size;
    std::vector<int> v(size);
    std::vector<float> vout(size);
    std::iota(v.begin(), v.end(), 0);
    ankerl::nanobench::Bench().run("std::ranges::transform", [&]()
    {
        ankerl::nanobench::doNotOptimizeAway(
            std::ranges::transform(v.cbegin(), v.cend(), vout.begin(), [=](const int x) {return std::cos(static_cast<float>(x)/static_cast<float>(size));})
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
            vout[i] = std::cos(static_cast<float>(v[i])/static_cast<float>(size));
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