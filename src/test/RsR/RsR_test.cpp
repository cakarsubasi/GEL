//
// Created by Cem Akarsubasi on 4/15/25.
//

#include <GEL/HMesh/RsR.h>
#include <GEL/HMesh/RsR2.h>
#include <GEL/HMesh/HMesh.h>
#include <GEL/Util/RawObj.h>

#include <filesystem>
#include <ranges>
#include <array>
#include <string_view>

#include <nanobench.h>

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest.h>

using namespace HMesh::RSR;

static constexpr auto FILE_CAPITAL_A = "../../../../data/PointClouds/Capital_A.obj";
static constexpr auto FILE_BUNNY_SIMPLE = "../../../../data/bunny.obj";
static constexpr auto FILE_BUNNY_COMPLEX = "../../../../data/PointClouds/bun_complete.obj";
static constexpr auto FILE_THINGY = "../../../../data/thingy.obj";
static constexpr auto FILE_AS = "../../../../data/as.obj";

template <typename... Args>
constexpr auto make_array(Args... args)
{
    constexpr auto size = sizeof...(args);
    return std::array<std::string_view, size>{args...};
}

static constexpr auto QUICK_TEST_FILES = make_array(FILE_CAPITAL_A, FILE_BUNNY_SIMPLE, FILE_THINGY, FILE_AS);
static constexpr auto TEST_FILES = make_array(FILE_CAPITAL_A, FILE_BUNNY_SIMPLE, FILE_BUNNY_COMPLEX, FILE_THINGY, FILE_AS);

constexpr auto IS_EUCLIDEAN = false;
constexpr auto K_PARAM = 30;
constexpr auto GENUS = -1;
constexpr auto R_PARAM = 20;
constexpr auto THETA = 60;
constexpr auto N_PARAM = 50;

auto test_options()
{
    HMesh::RSR::RsROpts opts;
    opts.dist = IS_EUCLIDEAN ? Distance::EUCLIDEAN : Distance::NEIGHBORS;
    opts.k = K_PARAM;
    opts.genus = GENUS;
    opts.r = R_PARAM;
    opts.theta = THETA;
    opts.n = N_PARAM;
    return opts;
}

template <typename Collection>
auto indices_from(const Collection& collection) -> std::vector<size_t>
{
    const auto indices = [&collection] {
        std::vector<NodeID> temp(collection.size());
        std::iota(temp.begin(), temp.end(), 0);
        return temp;
    }();
    return indices;
}

template <typename T>
auto indexed_select(const std::vector<T>& vec, const std::vector<size_t>& indices) -> std::vector<T>
{
    std::vector<T> result;
    result.reserve(indices.size());
    for (auto idx : indices)
    {
        result.push_back(vec.at(idx));
    }
    return result;
}

inline
auto angle_cost(const Vec3d e1, const Vec3d e2) -> double
{
    return std::abs(2.0 * dot(e1, e2) / (e1.length() * e2.length() + 1e-18) - 1.0);
}

inline
auto edge_cost(const Vec3d e1, const Vec3d e2) -> double
{
    const Vec3d e3 = e1 - e2;
    return angle_cost(e1, e2) * e3.length();
}

inline
auto triangle_cost(const Vec3d p1, const Vec3d p2, const Vec3d p3) -> double
{
    return edge_cost((p2 - p1), (p3 - p1)) + edge_cost((p1 - p2), (p3 - p2)) + edge_cost((p1 - p3), (p2 - p3));
}

auto manifold_cost(const HMesh::Manifold& manifold) -> double
{
    double total_cost = 0.0;
    for (const auto f: manifold.faces()) {
        std::array<HMesh::VertexID, 3> arr;
        int idx = 0;
        HMesh::circulate_face_ccw(manifold, f, [&idx, &arr](const HMesh::VertexID v) {
            if (idx < 3)
                arr[idx++] = v;
        });
        if (idx == 3) {
            const auto p1 = manifold.positions[arr[0]];
            const auto p2 = manifold.positions[arr[1]];
            const auto p3 = manifold.positions[arr[2]];
            const auto cost = triangle_cost(p1, p2, p3);
            total_cost += cost;
        }
    }
    return total_cost;
}

auto manifold_average_edge_length(const HMesh::Manifold& manifold) -> double
{
    double total = 0.0;
    for (const auto halfedge : manifold.halfedges()) {
        HMesh::Walker walker = manifold.walker(halfedge);
        auto v1 = walker.vertex();
        auto v2 = walker.opp().vertex();
        GEL_ASSERT_NEQ(v1, v2);
        auto p1 = manifold.positions[v1];
        auto p2 = manifold.positions[v2];
        total += std::sqrt(dot((p2 - p1), (p2 - p1)));
    }
    return total / static_cast<double>(manifold.no_halfedges());
}

auto manifold_cost_normalized(const HMesh::Manifold& manifold) -> double
{
    const auto number_of_faces = manifold.no_faces();
    const auto total_cost = manifold_cost(manifold);

    return total_cost / static_cast<double>(number_of_faces) / manifold_average_edge_length(manifold);
}

auto manifold_is_triangular(const HMesh::Manifold& manifold) -> bool
{
    for (const auto f: manifold.faces()) {
        int idx = 0;
        HMesh::circulate_face_ccw(manifold, f, [&idx](const HMesh::VertexID v) {
            idx++;
        });
        if (idx != 3)
            return false;
    }
    return true;
}

auto manifold_face_vert_ratio(const HMesh::Manifold& manifold) -> double
{
    return static_cast<double>(manifold.no_faces()) / static_cast<double>(manifold.no_vertices());
}

/// Estimate quality of a Manifold based on dihedral angles.
/// @details The dihedral angle between each face is used to determine how expensive it would be to "cross" between
/// those two faces. The average of this for every half-edge is used to get the total cost. Boundary edges are considered
/// to be the most expensive to cross.
///
/// This gives a quality metric in terms of how "smooth" a Manifold is. In an infinite flat plane, this function is 0.
/// However, keep in mind that this is not a generalized quality measure. The Manifold might be intentionally bumpy which
/// will return a high cost.
auto manifold_dihedral_cost(const HMesh::Manifold& manifold) -> double
{
    double total_cost = 0.0;
    size_t bad_edges = 0;
    for (const auto half_edge: manifold.halfedges()) {
        const auto f1 = manifold.walker(half_edge).face();
        const auto f2 = manifold.walker(half_edge).opp().face();
        if (f1 != HMesh::InvalidFaceID && f2 != HMesh::InvalidFaceID) {
            const auto n1 = manifold.normal(f1);
            const auto n2 = manifold.normal(f2);
            const auto cost = std::abs(CGLA::dot(CGLA::normalize(n1), CGLA::normalize(n2)) - 1.0);
            if (!std::isnan(cost))
                total_cost += cost;
            else
                bad_edges++;
        } else {
            // boundary edges are considered to be the most expensive
            total_cost += 2.0;
        }
    }
    if (bad_edges > 0) {
        std::cerr << "dihedral cost found " << bad_edges << " bad edges\n";
    }
    return total_cost / static_cast<double>(manifold.no_halfedges());
}

struct ValencyHistogram {
    static constexpr size_t max_tracked_valency = 9;
    size_t total = 0;
    std::array<size_t, max_tracked_valency> common_valency{0};
    size_t high = 0;

    friend std::ostream& operator<<(std::ostream& os, const ValencyHistogram& hist);
};

std::ostream& operator<<(std::ostream& os, const ValencyHistogram& hist)
{
    os << "Total: " << hist.total << "\n";
    constexpr auto draw_factor = 40;
    for (size_t valency = 2; valency < hist.common_valency.size(); ++valency) {
        const auto ratio = static_cast<double>(hist.common_valency[valency]) / static_cast<double>(hist.total);
        const auto blocks = static_cast<int>(ratio * draw_factor);
        os << valency << ":  ";
        for (int i = 0; i < blocks; ++i) {
            os << "#";
        }
        for (int i = blocks; i < draw_factor; ++i) {
            os << " ";
        }
        os << " (" << (ratio * 100.0) << "%)\n";
    }
    const auto ratio = static_cast<double>(hist.high) / static_cast<double>(hist.total);
    const auto blocks = static_cast<int>(ratio * draw_factor);
    os  << hist.common_valency.size() << "+: ";
    for (int i = 0; i < blocks; ++i) {
        os << "#";
    }
    for (int i = blocks; i < draw_factor; ++i) {
        os << " ";
    }
    os << " (" << (ratio * 100.0) << "%)\n";
    return os;
}

auto manifold_valency_histogram(const HMesh::Manifold& manifold) -> ValencyHistogram
{
    ValencyHistogram hist;
    for (const auto v: manifold.vertices()) {
        if (const auto valency = manifold.valency(v); valency > 8) {
            hist.high++;
        } else {
            hist.common_valency[valency]++;
        }
    }
    hist.total = manifold.no_vertices();
    return hist;
}

auto manifold_valency_histogram_non_boundary(const HMesh::Manifold& manifold) -> ValencyHistogram
{
    ValencyHistogram hist;
    size_t non_boundary_vertices = 0;
    for (const auto v: manifold.vertices()) {
        if (!manifold.boundary(v)) {
            non_boundary_vertices++;
            if (const auto valency = manifold.valency(v); valency > 8) {
                hist.high++;
            } else {
                hist.common_valency[valency]++;
            }
        }
    }
    hist.total = non_boundary_vertices;
    return hist;
}

// Test functions begin

auto test_reconstruct_legacy(const std::string_view file_name, const HMesh::RSR::RsROpts& opts) -> std::optional<HMesh::Manifold>
{
    std::cout << "======================\n"
    << "Begin original function\n";
    auto input_maybe = Util::read_raw_obj(file_name);
    if (!input_maybe) {
        WARN_FALSE("File not found.");
        return std::nullopt;
    }
    auto input = *std::move(input_maybe);

    std::cout << "obj vertices: " << input.vertices.size() << "\n";
    std::cout << "obj normals: " << input.normals.size() << "\n";

    auto result = point_cloud_normal_estimate(input.vertices, {}, true);

    HMesh::Manifold output;
    input.normals = {};
    reconstruct_single(output, result.vertices, result.normals, opts.dist == Distance::EUCLIDEAN,  opts.genus, opts.k, opts.r, opts.theta, opts.n);
    // k: 70 is too large
    // r: needs isEuclidean false
    std::cout << output.positions.size() << "\n";

    return output;
}

auto test_reconstruct_collapse_reexpand(const std::string_view file_name, const HMesh::RSR::RsROpts& opts, const int max_iterations, const bool reexpand) -> std::optional<HMesh::Manifold>
{
    std::cout << "======================\n"
    << "Begin new function\n";
    auto input_maybe = Util::read_raw_obj(file_name);
    if (!input_maybe) {
        WARN_FALSE("File not found.");
        return std::nullopt;
    }
    const auto input = *std::move(input_maybe);

    std::cout << "obj vertices: " << input.vertices.size() << "\n";
    std::cout << "obj normals: " << input.normals.size() << "\n";

    HMesh::Manifold output = point_cloud_collapse_reexpand(input.vertices, {}, opts, max_iterations, reexpand);
    // k: 70 is too large
    // r: needs isEuclidean false
    std::cout << output.positions.size() << "\n";

    return output;
}

auto reconstruct_assertions(const HMesh::Manifold& manifold) -> void
{
    // Manifold properties:
    // Triangular
    // Approx. 2 face per vertex (if good reconstruction)
    // Each vertex around 6 edges
    CHECK_EQ(manifold_is_triangular(manifold), true);
    auto ratio = manifold_face_vert_ratio(manifold);
    CHECK_GT(ratio, 1.5);
    CHECK_LT(ratio, 2.25);
    auto faces = manifold.no_faces();
    auto cost = manifold_cost(manifold);
    auto average_edge = manifold_average_edge_length(manifold);
    auto normed_cost = cost / static_cast<double>(faces) / average_edge;
    auto valencies = manifold_valency_histogram(manifold);
    auto valencies_nb = manifold_valency_histogram_non_boundary(manifold);
    auto dihedral_cost = manifold_dihedral_cost(manifold);
    std::cout << "ratio: " << ratio << std::endl;
    std::cout << "total cost: " << cost << std::endl;
    std::cout << "num faces : " << faces << std::endl;
    std::cout << "average edge length: " << average_edge << std::endl;
    std::cout << "normed cost: " << normed_cost << std::endl;
    std::cout << "Valency histogram (all vertices): " << std::endl;
    std::cout << valencies;
    std::cout << "Valency histogram (non-boundary vertices): " << std::endl;
    std::cout << valencies_nb;
    std::cout << "Dihedral cost: " << dihedral_cost << std::endl;
}

TEST_CASE("reconstruct old")
{
    auto opts_euclidean = test_options();
    auto file = FILE_BUNNY_SIMPLE;
    opts_euclidean.dist = Distance::EUCLIDEAN;
    auto manifold = test_reconstruct_legacy(file, opts_euclidean);
    if (manifold.has_value()) {
        HMesh::obj_save("debug_obj_reexpand_euclidean.obj", *manifold);
        reconstruct_assertions(*manifold);
    }
    auto opts_neighbors = test_options();
    opts_neighbors.dist = Distance::NEIGHBORS;
    manifold = test_reconstruct_legacy(file, opts_neighbors);
    if (manifold.has_value()) {
        HMesh::obj_save("debug_obj_reexpand_neighbors.obj", *manifold);
        reconstruct_assertions(*manifold);
    }
}

TEST_CASE("reconstruct new")
{
    auto opts_euclidean = test_options();
    auto file = FILE_BUNNY_SIMPLE;
    auto max_iterations = 0;
    auto reexpand = false;
    opts_euclidean.dist = Distance::EUCLIDEAN;
    auto manifold = test_reconstruct_collapse_reexpand(file, opts_euclidean, max_iterations, reexpand);
    if (manifold.has_value()) {
        HMesh::obj_save("debug_obj_reexpand_euclidean.obj", *manifold);
        reconstruct_assertions(*manifold);
    }
    auto opts_neighbors = test_options();
    opts_neighbors.dist = Distance::NEIGHBORS;
    manifold = test_reconstruct_collapse_reexpand(file, opts_neighbors, max_iterations, reexpand);
    if (manifold.has_value()) {
        HMesh::obj_save("debug_obj_reexpand_neighbors.obj", *manifold);
        reconstruct_assertions(*manifold);
    }
}