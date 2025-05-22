//
// Created by Cem Akarsubasi on 4/15/25.
//

#include <GEL/HMesh/RsR.h>
#include <GEL/HMesh/RsR2.h>
#include <GEL/HMesh/HMesh.h>
#include <GEL/Util/RawObj.h>

using GEL::HMesh::RSR::point_cloud_to_mesh;
using namespace GEL::HMesh::RSR;

static constexpr auto file_name =
//        "../../../../data/PointClouds/Capital_A.obj";
//        "../../../../data/bunny.obj";
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

auto test_options()
{
    GEL::HMesh::RSR::RsROpts opts;
    opts.is_euclidean = is_euclidean;
    opts.k = k;
    opts.genus = genus;
    opts.r = r;
    opts.theta = theta;
    opts.n = n;
    return opts;
}

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
    auto input = GEL::Util::read_raw_obj(file_name);

    std::cout << "obj vertices: " << input.vertices.size() << "\n";
    std::cout << "obj normals: " << input.normals.size() << "\n";

    const auto opts = test_options();
    HMesh::Manifold output = point_cloud_to_mesh(input.vertices, input.normals, opts);
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
    auto input = GEL::Util::read_raw_obj(file_name);

    std::cout << "obj vertices: " << input.vertices.size() << "\n";
    std::cout << "obj normals: " << input.normals.size() << "\n";

    HMesh::Manifold output;

    reconstruct_single(output, input.vertices, input.normals, is_euclidean,  genus, k, r, theta, n);
    // k: 70 is too large
    // r: needs isEuclidean false
    std::cout << output.positions.size() << "\n";
    //const bool result = HMesh::obj_save(output_name, output);

    return output;
}

template <typename T>
auto indexed_select(const std::vector<T>& vec, const std::vector<size_t>& indices) -> std::vector<T>
{
    std::vector<T> result;
    result.reserve(indices.size());
    for (auto idx : indices)
    {
        result.push_back(vec[idx]);
    }
    return result;
}

auto test_collapse() -> void
{
    auto input = GEL::Util::read_raw_obj(file_name);
    auto collapse = GEL::HMesh::RSR::point_cloud_collapse(input.vertices, input.normals);
    auto new_vertices = indexed_select(input.vertices, collapse.important_points);
    auto new_normals = indexed_select(input.normals, collapse.important_points);
    auto obj = GEL::Util::RawObj(new_vertices, new_normals, std::vector<Vec2d>());
    auto mesh = point_cloud_to_mesh(new_vertices, new_normals, test_options());

    GEL::Util::write_raw_obj("bun_collapsed.obj", obj);
    HMesh::obj_save("bun_collapsed_reconstructed.obj", mesh);
}

int main()
{
    //auto left = test_new();
    try {
        test_collapse();
    } catch (const std::exception& e) {}

    //auto right = test_old();
    //assert(manifold_is_identical(left, right));
    return 0;
}