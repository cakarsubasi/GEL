#ifndef GEL_HMesh_RsR2_hpp
#define GEL_HMesh_RsR2_hpp
#pragma once

// TODO: unnecessary imports will negatively affect compile times
#include <random>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <GEL/Geometry/Graph.h>
#include <GEL/HMesh/Manifold.h>
#include <GEL/Geometry/etf.h>
#include <GEL/Geometry/KDTree.h>
#include <GEL/Geometry/normal.h>

namespace GEL::HMesh::RSR
{

template<typename IndexTy = size_t, typename FloatTy = double>
class RSR {
public:
    RSR() = delete;

};

template <typename IndexTy>
struct Triangle {
    IndexTy v1, v2, v3;
    constexpr IndexTy operator[](const size_t idx) {
        assert(idx < 3);
        return std::bit_cast<IndexTy*>(this)[idx];
    }
};



using namespace CGLA;
using namespace Geometry;

using NodeID = AMGraph::NodeID;
using FaceType = std::vector<NodeID>;

using Vec3 = Vec3d;
using Point = Vec3;
using TEdge = std::pair<NodeID, NodeID>;

using Tree = Geometry::KDTree<Point, NodeID>;
using Record = Geometry::KDTreeRecord<Point, NodeID>;

struct Collapse {
    std::vector<NodeID> important_points;
    std::vector<std::vector<NodeID>> unimportant_points;
};

struct NeighborInfo {
    NodeID id;
    double distance; // the added precision doesn't justify the usage of doubles

    static constexpr NodeID invalid_id = -1;
    NeighborInfo() = delete;

    explicit NeighborInfo(const Record& record) noexcept : id(record.v), distance(std::sqrt(record.d))
    {}
};

using NeighborArray = std::vector<NeighborInfo>;
using NeighborMap = std::vector<NeighborArray>;

auto collapse_points(
    const std::vector<Point>& vertices,
    const std::vector<Vec3>& normals,
    const NeighborMap& neighbor_map,
    double collapse_factor
) -> Collapse;

double cal_radians_3d(const Vec3& branch, const Vec3& normal);

double cal_radians_3d(const Vec3& branch_vec, const Vec3& normal,
                      const Vec3& ref_vec);

///
/// TODO: documentation
struct RsROpts {
    int32_t genus = -1;
    int32_t k = 70;
    double r = 20;
    double theta = 60;
    int32_t n = 50;
    bool is_euclidean = false;

    /// Are normals included with the input?
    bool normals_included = true;
    bool isFaceNormal = true;
    bool isFaceLoop = true;
};

/// A trivial wrapper over bool to avoid the std::vector<bool> specialization
struct Boolean {
    bool inner;
};

/*Graph definition. The RsR graph here is integrated with the rotation system based on AMGraph*/
struct Vertex {
    NodeID id = 0;
    int normal_rep = -1;
    Vec3 coords = Vec3(0., 0., 0.);
    Vec3 normal = Vec3(0., 0., 0.);
    std::vector<Boolean> faceExist;
    float distance = 0.0f;

    // Hash implementation?
    struct Neighbor {
        double angle;
        uint v;
        mutable
        uint tree_id = 0;
        mutable
        bool faceExist = false;

        Neighbor(const Vertex& u, const Vertex& v, const uint id)
        {
            this->v = id;
            this->angle = cal_radians_3d(v.coords - u.coords, u.normal);
        }

        bool operator <(const Neighbor& rhs) const
        {
            return this->angle < rhs.angle || this->angle == rhs.angle && this->v != rhs.v;
        }

        bool operator<(Neighbor& rhs) const
        {
            return this->angle < rhs.angle || this->angle == rhs.angle && this->v != rhs.v;
        }
    };

    // TODO: mutable elements inside std::set is potentially unsound
    std::set<Neighbor> ordered_neighbors;
};

// TODO: 32 bytes per edge is quite heavy
struct Edge {
    NodeID source = -1;
    NodeID target = -1;
    double weight = 0.;
    // ?
    int ref_time = 0;
    // int count_weight = 1;
};

typedef Vertex::Neighbor Neighbor;

class SimpGraph : public AMGraph {
public:
    Util::AttribVec<AMGraph::EdgeID, Edge> m_edges;

    EdgeID connect_nodes(const NodeID source, const NodeID target, const float weight = 0.)
    {
        const EdgeID id = AMGraph::connect_nodes(source, target);
        m_edges[id].weight = weight;
        return id;
    }

    [[nodiscard]] double get_weight(const NodeID n1, const NodeID n2) const
    {
        return m_edges[find_edge(n1, n2)].weight;
    }

    /** Disconnect nodes. This operation removes the edge from the edge maps of the two formerly connected
         vertices, but the number of edges reported by the super class AMGraph is not decremented, so the edge is only
         invalidated. Call cleanup to finalize removal. */
    void disconnect_nodes(const NodeID n0, const NodeID n1)
    {
        if (valid_node_id(n0) && valid_node_id(n1)) {
            edge_map[n0].erase(n1);
            edge_map[n1].erase(n0);
        }
    }
};

// TODO: make this thread safe
class RSGraph : public AMGraph {
public:
    double total_edge_length = 0.;
    // int face_loop_id = 0;
    bool isEuclidean = false;
    bool isFinal = false;
    int exp_genus = -1;
    ETF etf;

    int current_no_edges = 0;

    Util::AttribVec<NodeID, Vertex> m_vertices;
    Util::AttribVec<AMGraph::EdgeID, Edge> m_edges;

    /// Compute sqr distance between two nodes - not necessarily connected.
    [[nodiscard]] double dist(const NodeID n0, const NodeID n1) const
    {
        if (valid_node_id(n0) && valid_node_id(n1))
            return (m_vertices[n0].coords - m_vertices[n1].coords).length();
        else
            return CGLA::CGLA_NAN;
    }

    /// Compute the average edge length
    [[nodiscard]] double cal_average_edge_length() const
    {
        const double sum_len = this->total_edge_length;
        return sum_len / current_no_edges;
    }

    void remove_edge(const NodeID source, const NodeID target)
    {
        this->total_edge_length -= this->m_edges[edge_map[source][target]].weight;
        if (valid_node_id(source) && valid_node_id(target)) {
            edge_map[source].erase(target);
            edge_map[target].erase(source);
        }
        current_no_edges--;
    }

    void remove_neighbor(const NodeID root, const NodeID neighbor)
    {
        auto& u = m_vertices[root];
        const auto& v = m_vertices[neighbor];
        u.ordered_neighbors.erase(Neighbor(u, v, neighbor));
    }

    void insert_neighbor(const NodeID root, const NodeID neighbor)
    {
        const auto& u = m_vertices[root];
        const auto& v = m_vertices[neighbor];
        m_vertices[root].ordered_neighbors.insert(Neighbor(u, v, neighbor));
    }

    EdgeID add_edge(const NodeID source, const NodeID target, const float weight = 0.)
    {
        const EdgeID id = this->connect_nodes(source, target);
        assert(id != InvalidEdgeID);
        // TODO: eliminating this check results in a crash, figure out where the UB actually is
        if (id != InvalidEdgeID) {
            current_no_edges++;
            m_edges[id].weight = weight;
            m_edges[id].source = source;
            m_edges[id].target = target;
            this->total_edge_length += weight;
            insert_neighbor(source, target);
            insert_neighbor(target, source);
        } else {
            /*std::cout << "weird" << std::endl;
            std::cout << valid_node_id(source) << std::endl;
            std::cout << valid_node_id(target) << std::endl;
            std::cout << source << std::endl;
            std::cout << target << std::endl;
            std::cout << (find_edge(target, source) == InvalidEdgeID)
                << std::endl;*/
        }

        return id;
    }

    NodeID add_node(const Vec3& p)
    {
        const NodeID n = AMGraph::add_node();
        Vertex v;
        v.id = n;
        v.coords = p;
        m_vertices[n] = std::move(v);
        return n;
    }

    NodeID add_node(const Vec3& p, const Vec3& in_normal)
    {
        const NodeID n = AMGraph::add_node();
        Vertex v;
        v.id = n;
        v.coords = p;
        v.normal = in_normal;
        m_vertices[n] = std::move(v);
        return n;
    }

    void init(const std::vector<Point>& vertices, const std::vector<Vec3>& normals)
    {
        for (size_t i = 0; i < vertices.size(); i++) {
            const NodeID id = this->add_node(vertices[i]);
            m_vertices[id].normal = normals[i];
        }
    }

    void init(const std::vector<Point>& vertices)
    {
        for (const auto& vertex : vertices) {
            this->add_node(vertex);
        }
    }

    void init(const int no_vertex)
    {
        for (int i = 0; i < no_vertex; i++) {
            AMGraph::add_node();
        }
    }

    void get_node_set(NodeSet& sets) const
    {
        for (NodeID i = 0; i < edge_map.size(); i++) {
            sets.insert(i);
        }
    }
};



void NN_search(const Point&, const Tree&, double,
               std::vector<NodeID>&, std::vector<double>&, bool isContain = true);

void init_graph(const std::vector<Point>& vertices, const std::vector<Point>& smoothed_v,
                const std::vector<Vec3>& normals, const Tree& kdTree, SimpGraph& dist_graph,
                std::vector<float>& max_length, std::vector<float>& pre_max_length, float cross_conn_thresh, int k,
                bool isEuclidean);

int find_shortest_path(const RSGraph& mst, NodeID start, NodeID target,
                       int threshold, std::vector<NodeID>& path);


void minimum_spanning_tree(
    const SimpGraph& g,
    NodeID root,
    RSGraph& gn,
    const std::vector<Vec3>& normals,
    const std::vector<Point>& vertices,
    bool isEuclidean
);

void minimum_spanning_tree(const SimpGraph& g, NodeID root, SimpGraph& gn);

void add_face(RSGraph& G, const std::vector<NodeID>& item,
              std::vector<std::vector<NodeID>>& faces);

void connect_handle(const std::vector<Point>& smoothed_v, Tree& KDTree,
                    RSGraph& mst, std::vector<NodeID>& connected_handle_root, int k,
                    int step_thresh, bool isEuclidean);

// Face Loop

void init_face_loop_label(RSGraph& g);

const Neighbor& successor(const RSGraph& g,
                          const NodeID& root,
                          const NodeID& branch);


const Neighbor& predecessor(const RSGraph& g,
                            const NodeID& root,
                            const NodeID& branch);

void maintain_face_loop(RSGraph& g, NodeID source, NodeID target);

const Neighbor& get_neighbor_info(const RSGraph& g, const NodeID& root, const NodeID& branch);

// Utils

void find_common_neighbor(NodeID neighbor, NodeID root, std::vector<NodeID>& share_neighbor, RSGraph& g);

// Algorithm

// bool geometry_check(const RSGraph& mst, const TEdge& candidate, const Tree& kdTree);

bool Vanilla_check(RSGraph& mst, const TEdge& candidate, const Tree& kdTree);

bool isIntersecting(const RSGraph& mst, NodeID v1, NodeID v2, NodeID v3, NodeID v4);

bool routine_check(const RSGraph& mst, const std::vector<NodeID>& triangle);

/// Standalone point cloud collapse, normals are also updated if empty
/// @param vertices vertices of the point cloud
/// @param normals normals of the point cloud or empty vector
/// @return point cloud collapse information
auto point_cloud_collapse(const std::vector<Point>& vertices, std::vector<Vec3>& normals) -> Collapse;

/// Convert point cloud to a Manifold
/// @param vertices vertices of the point cloud
/// @param normals normals of the point cloud or empty vector
/// @param opts collapse options
/// @return reconstructed manifold mesh
auto point_cloud_to_mesh(const std::vector<Point>& vertices, const std::vector<Vec3>& normals,
                         const RsROpts& opts) -> ::HMesh::Manifold;
} // namespace GEL::HMesh::RSR

#endif // GEL_HMesh_RsR2_hpp
