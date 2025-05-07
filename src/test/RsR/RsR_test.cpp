//
// Created by Cem Akarsubasi on 4/15/25.
//
#include <iostream>
#include <GEL/HMesh/RsR.h>
#include <GEL/HMesh/RsR2.h>
#include <GEL/HMesh/HMesh.h>
#include <GEL/Geometry/load.h>

using GEL::HMesh::RSR::point_cloud_to_mesh;

static constexpr auto file_name =
//        "../../../../data/PointClouds/Capital_A.obj";
        //    "../../../../data/bunny.obj";
//    "../../../../data/torus.obj";
    "../../../../data/PointClouds/bun_complete.obj";
static constexpr auto output_name =
        "bunny.obj";

constexpr auto is_euclidean = false;
constexpr auto k = 30;
constexpr auto genus = -1;
constexpr auto r = 20;
constexpr auto theta = 60;
constexpr auto n = 50;

auto manifold_is_identical(const HMesh::Manifold& left, const HMesh::Manifold& right) -> bool
{
    // This is a horrendous way of actually checking if two manifolds are identical,
    // but assuming we did not mess something up during construction, they should be
    // using identical IDs, which is good enough for quick regression analysis

    const auto left_edges = left.halfedges();
    const auto right_edges = right.halfedges();
    for (auto left_begin = left_edges.begin(), right_begin = right_edges.begin();
         left_begin != left_edges.end() && right_begin != right_edges.end() ;
         ++left_begin, ++right_begin)
    {
        if (*left_begin != *right_begin)
        {
            return false;
        }
    }
    return true;
}

auto test_new() -> HMesh::Manifold
{
    std::cout << "======================\n"
    << "Begin new function\n";
    TriMesh input;
    load(file_name, input);
    //assert(success);

    std::cout << "obj vertices: " << input.geometry.no_vertices() << "\n";
    std::cout << "obj normals: " << input.normals.no_vertices() << "\n";

    std::vector<Point> points;
    std::vector<Vector> normals;

    for (auto& point: input.geometry.vertices_copy())
    {
        points.emplace_back(point);
    }
    for (auto& point: input.normals.vertices_copy())
    {
        normals.emplace_back(point);
    }
    normals = {};
    GEL::HMesh::RSR::RsROpts opts;
    opts.isEuclidean = is_euclidean;
    opts.k = k;
    opts.genus = genus;
    opts.r = r;
    opts.theta = theta;
    opts.n = n;
    HMesh::Manifold output = point_cloud_to_mesh(points, normals, opts);
    // k: 70 is too large
    // r: needs isEuclidean false
    std::cout << output.positions.size() << "\n";
    const bool result = HMesh::obj_save(output_name, output);
    assert(result);
    return output;
}


auto test_old() -> HMesh::Manifold
{
    std::cout << "======================\n"
    << "Begin original function\n";
    TriMesh input;
    load(file_name, input);
    //assert(success);

    std::cout << "obj vertices: " << input.geometry.no_vertices() << "\n";
    std::cout << "obj normals: " << input.normals.no_vertices() << "\n";

    HMesh::Manifold output;
    std::vector<Point> points;
    std::vector<Vector> normals;

    for (auto& point: input.geometry.vertices_copy())
    {
        points.emplace_back(point);
    }
    for (auto& point: input.normals.vertices_copy())
    {
        normals.emplace_back(point);
    }
    normals = {};
    reconstruct_single(output, points, normals, is_euclidean,  genus, k, r, theta, n);
    // k: 70 is too large
    // r: needs isEuclidean false
    std::cout << output.positions.size() << "\n";
    //const bool result = HMesh::obj_save(output_name, output);

    return output;
}

int main()
{
    auto left = test_new();
    auto right = test_old();
    assert(manifold_is_identical(left, right));
    return 0;
}