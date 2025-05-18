#include <GEL/HMesh/RsR2.h>
#include <GEL/HMesh/Timer.h>

#include "GEL/Util/ParallelAdapters.h"


namespace GEL::HMesh::RSR
{

struct NeighborInfo {
    NodeID id;
    double distance;
    NeighborInfo() = delete;
    explicit NeighborInfo(const Record& record) noexcept : id(record.v), distance(std::sqrt(record.d))
    {}
};

using NeighborArray = std::vector<NeighborInfo>;
using NeighborMap = std::vector<NeighborArray>;

struct m_cmp {
    bool operator()(const std::pair<std::vector<NodeID>, float>& left,
                    const std::pair<std::vector<NodeID>, float>& right) const
    {
        return (left.second) > (right.second);
    }
};

// TODO: two and three layers of indirection here. Probably good to look into this later
using m_priority_queue = std::priority_queue<
    std::pair<std::vector<NodeID>, float>,
    std::vector<std::pair<std::vector<NodeID>, float>>,
    m_cmp>;

// TODO: It might be better if some of these are explicit structs
using m_Edge_length = std::pair<float, int>;
using m_face_pair = std::pair<int, std::string>;
using m_neighbor_pair = std::pair<double, NodeID>;

using ::Util::AttribVec;

inline bool edge_comparator(const m_Edge_length& l, const m_Edge_length& r)
{
    return l.first < r.first;
}

inline bool face_comparator(const m_face_pair& l, const m_face_pair& r)
{
    return l.first > r.first;
}

inline bool neighbor_comparator(const m_neighbor_pair& l, const m_neighbor_pair& r)
{
    return l.first > r.first;
}

template <typename Heap>
void heap_sort(Heap& heap)
{
    for (size_t i = heap.size(); i > 0; --i) {
        std::pop_heap(heap.begin(), heap.begin() + i);
    }
}

template <typename Q>
constexpr typename Q::value_type queue_pop_front(Q& queue)
{
    if constexpr (std::is_same_v<Q, std::priority_queue<typename Q::value_type>>) {
        typename Q::value_type front = std::move(queue.top());
        queue.pop();
        return front;
    } else {
        typename Q::value_type front = std::move(queue.front());
        queue.pop();
        return front;
    }
}

/// Erases an element and moves the last element in the vector in its place
/// running in O(1) time
template <typename T>
T erase_swap(std::vector<T>& vec, typename std::vector<T>::size_type idx)
{
    T back = std::move(vec.back());
    vec.pop_back();
    T current = std::move(vec[idx]);
    vec[idx] = std::move(back);
    return current;
}

template <typename Float>
Float safe_div(Float lhs, Float rhs)
{
    if (rhs == 0) {
        return static_cast<Float>(0);
    } else {
        return lhs / rhs;
    }
}

template <typename Float>
Vec3 safe_div(Vec3 lhs, Float rhs)
{
    if (rhs == 0) {
        return Vec3(0);
    } else {
        return lhs / rhs;
    }
}

/**
    * @brief k nearest neighbor search
    *
    * @param query: the coordinate of the point to be queried
    * @param kdTree: kd-tree for knn query
    * @param num: number of nearest neighbors to be queried
    * @param neighbors: [OUT] indices of k nearest neighbors
    * @param neighbor_distance: [OUT] corresponding distance to the query point
    * @param isContain: whether the query point itself count as a neighbor
    *
    * @return None
    */
// TODO: We can cache search result for every point after smoothing
// TODO: weighted smoothing version doesn't use distances
template <typename DistType>
void knn_search(const Point& query, const Tree& kdTree,
                int num, std::vector<NodeID>& neighbors,
                std::vector<DistType>& neighbor_distance, const bool isContain)
{
    if (!isContain)
        num -= 1;
    std::vector<Record> records = kdTree.m_closest(num + 1, query, INFINITY);

    // Sort to the normal order from heap
    heap_sort(records);

    int idx = 0;
    for (auto& record : records) {
        if (idx == 0 && isContain) {
            idx++;
            continue;
        }
        neighbors.push_back(record.v);
        neighbor_distance.push_back(std::sqrt(record.d));
        idx++;
    }
}



void knn_search(const Point& query, const Tree& kdTree,
                int num, NeighborArray& neighbors, const bool isContain)
{
    // TODO: this seems like such a pointless difference.
    // It might be a better idea for the caller to handle this to reduce some clutter
    // TODO: maybe make NeighborInfo ABI compatible with Record to directly move the records in, or alternatively
    // TODO: change the implementation in KdTree to have an inout parameter
    const std::vector<Record> records = [num, isContain, &kdTree, &query]() -> std::vector<Record> {
        std::vector<Record> temp =
            (!isContain) ?
                kdTree.m_closest(num, query, INFINITY) :
                kdTree.m_closest(num + 1, query, INFINITY);
        heap_sort(temp);
        return temp;
    }();

    int idx = 0;
    for (auto record : records) {
        [[unlikely]]
        if (idx == 0 && isContain) {

        } else {
            neighbors.emplace_back(record);
        }
        idx++;
    }
}

auto calculate_neighbors(
    Util::ThreadPool& pool,
    const std::vector<Point>& vertices,
    const Tree& kdTree,
    const int k,
    NeighborMap&& neighbors_memoized = NeighborMap())
-> NeighborMap
{
    if (neighbors_memoized.empty()) {
        neighbors_memoized = NeighborMap(vertices.size());
        for (auto &neighbors : neighbors_memoized) {
            neighbors.reserve(k);
        }
    } else if (neighbors_memoized.at(0).capacity() < k) {;
        // TODO: this seems to lose performance because calls to reserve block yet the allocator is multithreaded
        // for (auto &neighbors : neighbors_memoized) {
        //     neighbors.reserve(k);
        // }
    }

    auto cache_kNN_search = [kdTree, k](auto vertex, auto& neighbor) {
        knn_search(vertex, kdTree, k, neighbor, true);
    };
    GEL::Util::parallel_foreach2(pool, cache_kNN_search, vertices, neighbors_memoized);
    return neighbors_memoized;
}

/// TODO: consider more efficient ways of culling duplicate vertices
void remove_duplicate_vertices(
    Util::ThreadPool& pool,
    std::vector<Point>& vertices,
    std::vector<Vec3>& normals,
    const Tree& kdTree,
    const int k)
{
    std::vector<Point> new_vertices;
    std::vector<Vec3> new_normals;
    const auto vertices_num = vertices.size();

    NeighborMap neighbors_memoized = calculate_neighbors(pool, vertices, kdTree, k);

    auto lambda = [&neighbors_memoized](const size_t this_idx, Point vertex, Vec3 _normal)
    -> std::optional<std::pair<Point, Vec3>> {

        for (const auto& neighbors = neighbors_memoized[this_idx]; const auto neighbor : neighbors) {
            const NodeID idx = neighbor.id;
            const double length = neighbor.distance;

            if (this_idx == idx)
                continue;

            // Remove duplicate vertices
            if (length < 1e-8 && this_idx != idx) {
                return std::nullopt;
            }
        }
        return std::make_optional(std::make_pair(vertex, _normal));
    };

    auto lambda2 = [lambda](size_t this_idx, Point vertex) -> std::optional<Point> {
        if (auto value = std::invoke(lambda, this_idx, vertex, Point(0,0,0)); value.has_value())
            return std::make_optional(std::get<0>(*value));
        else
            return std::nullopt;
    };

    if (!normals.empty())
        GEL::Util::parallel_enumerate_map2_filter2(pool, lambda, vertices, normals, new_vertices, new_normals);
    else
        GEL::Util::parallel_enumerate_map_filter(pool, lambda2, vertices, new_vertices);

    std::cout << vertices_num - vertices.size() << " duplicate vertices removed." << std::endl;

    vertices = std::move(new_vertices);
    normals = std::move(new_normals);
}

/**
    * @brief Calculate the reference vector for the rotation system
    *
    * @param normal: normal direction for the target vertex
    *
    * @return the reference vector
    */
Vec3 calculate_ref_vec(const Vec3& normal)
{
    constexpr float eps = 1e-6;
    double second = normal[1];
    if (second == 0.)
        second += eps;
    auto ref_vec = Vec3(0, -normal[2] / second, 1);
    if (normal[2] == 1.)
        ref_vec = Vec3(0., 1., 0.);
    ref_vec.normalize();
    return ref_vec;
}

/**
    * @brief Calculate the radian in the rotation system
    *
    * @param branch: vector of the out-going edge
    * @param normal: normal of the root vertex
    *
    * @return radian
    */
double cal_radians_3d(const Vec3& branch, const Vec3& normal)
{
    const Vec3 proj_vec = branch - dot(normal, branch) /
        normal.length() * normal;

    const auto ref_vec = calculate_ref_vec(normal);

    if (proj_vec.length() == 0.0)
        return 0.;

    const Vec3 proj_ref = ref_vec - dot(normal, ref_vec) /
        normal.length() * normal;
    const auto value = std::clamp<double>(dot(proj_vec, proj_ref) / proj_vec.length() /
                                          proj_ref.length(), -1, 1);
    double radian = std::acos(value);
    if (dot(cross(proj_vec, proj_ref), normal) > 0)
        radian = 2 * M_PI - radian;

    [[unlikely]]
    if (std::isnan(radian)) {
        std::cout << normal << std::endl;
        std::cout << ref_vec << std::endl;
        std::cout << "error" << std::endl;
    }
    return radian;
}

/**
    * @brief Calculate the radian given the reference vector
    *
    * @param branch_vec: vector of the out-going edge
    * @param normal: normal of the root vertex
    * @param ref_vec: the reference vector
    *
    * @return radian
    */
double cal_radians_3d(const Vec3& branch_vec, const Vec3& normal, const Vec3& ref_vec)
{
    const Vec3 proj_vec = branch_vec - dot(normal, branch_vec) /
        normal.length() * normal;
    if (std::abs(proj_vec.length()) < 1e-8)
        return 0.;

    const Vec3 proj_ref = ref_vec - dot(normal, ref_vec) /
        normal.length() * normal;
    const auto value = std::clamp<double>(
        dot(proj_vec, proj_ref) / proj_vec.length() /
        proj_ref.length(), -1, 1);
    double radian = std::acos(value);
    if (dot(CGLA::cross(proj_vec, proj_ref), normal) > 0)
        radian = 2 * M_PI - radian;
    return radian;
}

//TODO: adapt kdtree
Tree build_KDTree(const std::vector<Point>& vertices, const std::vector<NodeID>& indices)
{
    Tree kdTree;
    kdTree.reserve(vertices.size());
    // safety precondition
    assert(indices.size() >= vertices.size());
    int idx = 0;
    for (const auto& vertex : vertices) {
        kdTree.insert(vertex, indices[idx]);
        idx++;
    }
    kdTree.build();
    return kdTree;
}

/**
    * @brief neighbor search within a specific radius
    *
    * @param query: the coordinate of the point to be queried
    * @param kdTree: kd-tree for knn query
    * @param dist: the radius of the search ball
    * @param neighbors: [OUT] indices of k nearest neighbors
    * @param neighbor_dist: [OUT] corresponding distance to the query point
    * @param isContain: does the query point itself count as a neighbor
    *
    * @return None
    */
void NN_search(const Point& query, const Tree& kdTree,
               const double dist, std::vector<NodeID>& neighbors,
               std::vector<double>& neighbor_dist, const bool isContain)
{
    std::vector<Point> neighbor_coords;
    kdTree.in_sphere(query, dist, neighbor_coords, neighbors);

    std::vector<m_neighbor_pair> paired;
    int i = 0;
    for (auto& neighbor_coord : neighbor_coords) {
        neighbor_dist.push_back((neighbor_coord - query).length());
        paired.emplace_back(neighbor_dist[i], neighbors[i]);
        i++;
    }
    // TODO: given how slow std::ranges::transform was, I need to actually benchmark
    // TODO: std::ranges::sort before replacing std::sort with it
    std::sort(paired.begin(), paired.end(), neighbor_comparator);
    if (!isContain && !paired.empty()) {
        neighbor_dist[0] = paired[0].first;
        neighbors[0] = paired[0].second;
    }
    for (size_t idx = 1; idx < paired.size(); ++idx) {
        neighbor_dist[0] = paired[0].first;
        neighbors[0] = paired[0].second;
    }
}

/**
    * @brief Calculate projection distance
    *
    * @param edge: the edge to be considered
    * @param this_normal: normal of one vertex
    * @param neighbor_normal: normal of another vertex
    *
    * @return projection distance
    */
double cal_proj_dist(const Vec3& edge, const Vec3& this_normal, const Vec3& neighbor_normal)
{
    const double Euclidean_dist = edge.length();
    const double neighbor_normal_length = dot(edge, normalize(neighbor_normal));
    const double normal_length = dot(edge, normalize(this_normal));
    double projection_dist = sqrt((Euclidean_dist * Euclidean_dist) - (normal_length * normal_length));
    projection_dist += sqrt((Euclidean_dist * Euclidean_dist) -
        (neighbor_normal_length * neighbor_normal_length));
    projection_dist /= 2.;
    if (std::abs(dot(normalize(this_normal), normalize(neighbor_normal))) < std::cos(15. / 180. * M_PI))
        projection_dist = Euclidean_dist;
    return projection_dist;
}

/**
    * @brief initialize the graph and related information
    *
    * @param vertices: vertices of the componnet
    * @param smoothed_v: smoothed vertices of the component
    * @param normals: normal of the component vertices
    * @param kdTree: kd-tree for neighbor query
    * @param dist_graph: [OUT] a light-weight graph with essential connection for building MST
    * @param max_length: [OUT] the distance of the longest connection each vertex involved
    * @param pre_max_length: [OUT] the maximum length of connection before connecting handles (conservative connection)
    * @param cross_conn_thresh: angle threshold to avoid connecting vertices on different surface
    * @param k
    * @param isEuclidean
    *
    *
    * @return None
    */
void init_graph(
    const std::vector<Point>& vertices,
    const std::vector<Point>& smoothed_v,
    const std::vector<Vec3>& normals,
    const Tree& kdTree,
    SimpGraph& dist_graph,
    std::vector<float>& max_length,
    std::vector<float>& pre_max_length,
    const float cross_conn_thresh,
    const int k,
    const bool isEuclidean)
{
    for (int i = 0; i < vertices.size(); i++) {
        dist_graph.add_node();
    }

    NodeID i = 0;
    for (auto& vertex : vertices) {
        Vec3 this_normal = normals[i];

        std::vector<NodeID> neighbors;
        std::vector<float> dists;
        // TODO: move this outside
        knn_search(smoothed_v[i], kdTree, k, neighbors, dists, true);
        pre_max_length[i] = dists[static_cast<size_t>(k * 2. / 3.)];

        // Filter out the cross connection
        {
            // TODO: we can probably do this in place without using temp at all
            std::vector<NodeID> temp;
            for (const auto idx : neighbors) {
                Vec3 neighbor_normal = normals[idx];
                const double cos_theta = dot(this_normal, neighbor_normal) /
                    this_normal.length() /
                    neighbor_normal.length();
                double cos_thresh = std::cos(cross_conn_thresh / 180. * M_PI);
                if (isEuclidean)
                    cos_thresh = 0.;
                if (cos_theta >= cos_thresh) {
                    temp.push_back(idx);
                }
            }
            if (temp.empty())
                std::cerr << "Bad normal input" << std::endl;
            else {
                neighbors.clear();
                neighbors = std::move(temp);
            }
        }

        // TODO: once again, connect_nodes can be run in parallel with one protected variable
        for (const unsigned long idx : neighbors) {
            if (dist_graph.find_edge(i, idx) != AMGraph::InvalidEdgeID)
                continue;
            [[unlikely]]
            if (idx == i) {
                std::cerr << "Vertex " << idx << " connect back to its own." << std::endl;
                continue;
            }
            Vec3 neighbor_normal = normals[idx];
            Point neighbor_pos = vertices[idx];
            Vec3 edge = neighbor_pos - vertex;
            const double Euclidean_dist = edge.length();
            double weight = Euclidean_dist;
            if (!isEuclidean) {
                weight = cal_proj_dist(edge, this_normal, neighbor_normal);
            }
            if (weight > max_length[i])
                max_length[i] = weight;
            if (weight > max_length[idx])
                max_length[idx] = weight;
            [[unlikely]]
            if (weight < 1e-8)
                std::cerr << "error" << std::endl;

            dist_graph.connect_nodes(i, idx, weight);
        }
        i++;
    }
}

/**
    * @brief Find the shortest path from one vertex to another in the graph
    *
    * @param mst: the graph
    * @param start: the source vertex
    * @param target: the target vertex
    * @param threshold: the step threshold, if longer than this threshold, the algorithm early stop (never used?)
    * @param path: stores the indices of vertex in the shortest path (for visualization)
    *
    * @return the number of steps of the shortest path
    */
int find_shortest_path(const RSGraph& mst, NodeID start, NodeID target, int threshold, std::vector<NodeID>& path)
{
    std::queue<NodeID> q;
    std::vector<int> dist(mst.no_nodes(), -1); // Distance from start to each node
    std::vector<NodeID> pred(mst.no_nodes(), -1); // Predecessor array for path reconstruction
    std::set<NodeID> visited;

    dist[start] = 0;
    q.push(start);

    while (!q.empty()) {
        const int u = q.front();
        q.pop();

        // If the target node is reached, stop early
        if (u == target) {
            break;
        }

        // If the node has already been visited, skip it
        if (visited.contains(u)) {
            continue;
        }

        visited.insert(u);

        // Explore neighbors
        for (const auto& v : mst.neighbors(u)) {
            if (dist[v] == -1) {
                // If the node hasn't been visited
                dist[v] = dist[u] + 1; // Increment distance
                pred[v] = u; // Record predecessor
                q.push(v);
            }
        }
    }

    // Reconstruct path from start to target
    if (dist[target] != -1) {
        // Target is reachable
        for (NodeID v = target; v != -1; v = pred[v]) {
            path.push_back(v);
        }
        std::ranges::reverse(path); // Reverse to get the correct order
    }

    return dist[target];
}

/**
* @brief weighted smoothing method using defined neighborhood with tangential distance weighted
*
* @param pool thread pool
* @param vertices: vertices of the point cloud
* @param normals: normal of the point cloud
* @param neighbors_
* @param smoothed_v: [OUT] vertices after smoothing
*
* @return None
*/
void weighted_smooth(
    Util::ThreadPool& pool,
    const std::vector<Point>& vertices,
    const std::vector<Vec3>& normals,
    const NeighborMap& neighbors_,
    std::vector<Point>& smoothed_v)
{
    auto lambda = [&normals, &vertices](const size_t idx, const Point& vertex, const NeighborArray& neighbors) {
        const Vec3 normal = normals[idx];

        double weight_sum = 0.;
        double amp_sum = 0.;
        double max_dist = 0.;

        std::vector<double> vertical_length;
        std::vector<double> weights;
        const ssize_t limit = (neighbors.size() < 192) ? neighbors.size() : 192;
        vertical_length.reserve(limit);
        weights.reserve(limit);
        for (auto begin = neighbors.cbegin(); begin != neighbors.cbegin() + limit; ++begin) {
            const auto& neighbor = *begin;
        //for (const auto& neighbor : neighbors) {
            Point neighbor_pos = vertices[neighbor.id];
            Vec3 n2this = neighbor_pos - vertex;
            if (dot(normals[neighbor.id], normal) < std::cos(30. / 180. * M_PI)) {
                continue;
            }
            double vertical = dot(n2this, normal);
            const double n_dist = (neighbor_pos - vertex).length();

            const double tangential_square = n_dist * n_dist -
                vertical * vertical;
            double tangential_dist = 0.;
            if (tangential_square > 0.)
                tangential_dist = std::sqrt(tangential_square);

            [[unlikely]]
            if (!std::isfinite(tangential_dist)) {
                std::cerr << n_dist << " " << vertical << std::endl;
                std::cerr << "error" << std::endl;
            }

            const double weight = -tangential_dist;
            if (tangential_dist > max_dist)
                max_dist = tangential_dist;

            weights.push_back(weight);
            vertical_length.push_back(vertical);
        }
        for (int i = 0; i < vertical_length.size(); i++) {
            amp_sum += vertical_length[i] * (weights[i] + max_dist);
            weight_sum += weights[i] + max_dist;
        }

        weight_sum = (weight_sum == 0.) ? 1. : weight_sum;

        amp_sum /= weight_sum;
        if (!std::isfinite(amp_sum))
            std::cout << "error" << std::endl;
        const Vec3 move = amp_sum * normal;
        return vertex + move;
    };
    Util::parallel_enumerate_map2(pool, lambda, vertices, neighbors_, smoothed_v);
}

auto normalize_normals(std::vector<Vec3>& normals) -> void
{
    for (auto& normal : normals) {
        normal.normalize();
    }
}

void estimate_normal_no_normals_memoized(
    Util::ThreadPool& pool,
    const std::vector<Point>& vertices,
    const NeighborMap& neighbors,
    std::vector<Vec3>& normals)
{
    normals.clear();
    // Data type transfer & Cal diagonal size
    auto lambda = [&](const NeighborArray& neighbors_of_this) {
        // need id, distance and coords anyway

        std::vector<Point> neighbor_coords;
        for (const auto neighbor_id : neighbors_of_this) {
            neighbor_coords.push_back(vertices[neighbor_id.id]);
        }
        const Vec3 normal = estimateNormal(neighbor_coords);
        [[unlikely]]
        if (std::isnan(normal.length())) {
            std::cerr << neighbors_of_this.size() << std::endl;
            std::cerr << "error" << std::endl;
        }
        return normal;
    };
    Util::parallel_map(pool, lambda, neighbors, normals);
}

/**
    * @brief Calculate cos angle weight for correcting normal orientation
    *
    * @param this_normal: normal of current vertex
    * @param neighbor_normal: normal of its neighbor vertex
    *
    * @return angle weight calculated
    */
float cal_angle_based_weight(const Vec3& this_normal, const Vec3& neighbor_normal)
{
    float dot_pdt = std::abs(CGLA::dot(this_normal, neighbor_normal) / this_normal.length() / neighbor_normal.length());
    dot_pdt = std::clamp<float>(dot_pdt, 0., 1.0);
    if (1. - dot_pdt < 0)
        std::cout << "error" << std::endl;
    return 1. - dot_pdt;
}

void minimum_spanning_tree(
    const SimpGraph& g,
    NodeID root,
    RSGraph& gn,
    const std::vector<Vec3>& normals,
    const std::vector<Point>& vertices,
    const bool isEuclidean)
{
    using QElem = std::tuple<double, NodeID, NodeID>;
    if (root == AMGraph::InvalidNodeID)
        root = 0;

    for (const auto n : g.node_ids())
        gn.add_node(vertices[n], normals[n]);

    AttribVec<NodeID, unsigned char> in_tree(gn.no_nodes(), false);

    std::priority_queue<QElem> queue;
    for (auto n : g.neighbors(root)) {
        const auto d = CGLA::sqr_length(vertices[n] - vertices[root]);
        queue.emplace(-d, root, n);
    }

    while (!queue.empty()) {
        auto [d, n, m] = queue_pop_front(queue);

        if (!in_tree[m]) {
            in_tree[m] = true;

            Vec3 edge = gn.m_vertices[m].coords - gn.m_vertices[n].coords;
            const double Euclidean_dist = edge.length();
            const double projection_dist = cal_proj_dist(edge, gn.m_vertices[m].normal, gn.m_vertices[n].normal);
            [[unlikely]]
            if (std::isnan(projection_dist) || std::isnan(Euclidean_dist))
                std::cerr << "debug" << std::endl;

            if (isEuclidean)
                gn.add_edge(m, n, Euclidean_dist);
            else
                gn.add_edge(m, n, projection_dist);

            //gn.connect_nodes(n, m);
            for (auto nn : g.neighbors(m)) {
                const auto d_nn_m = CGLA::sqr_length(vertices[nn] - vertices[m]);
                queue.emplace(-d_nn_m, m, nn);
            }
        }
    }
}

void minimum_spanning_tree(const SimpGraph& g, NodeID root, SimpGraph& gn)
{
    using NodeID = NodeID;
    using QElem = std::tuple<double, NodeID, NodeID>;
    if (root == AMGraph::InvalidNodeID)
        root = 0;

    for (auto n : g.node_ids())
        gn.add_node();

    AttribVec<NodeID, unsigned char> in_tree(gn.no_nodes(), false);

    std::priority_queue<QElem> Q;
    for (auto n : g.neighbors(root)) {
        const auto d = g.get_weight(n, root);
        //auto d = CGLA::sqr_length(g.m_vertices[n].coords - g.m_vertices[root].coords);
        Q.emplace(-d, root, n);
    }

    while (!Q.empty()) {
        auto [d, n, m] = Q.top();
        Q.pop();

        if (!in_tree[m]) {
            in_tree[m] = true;
            gn.connect_nodes(n, m);
            for (auto nn : g.neighbors(m)) {
                const auto d_nn_m = g.get_weight(nn, m);
                Q.emplace(-d_nn_m, m, nn);
            }
        }
    }
}

/**
    * @brief Determine the normal orientation
    *
    * @param kdTree
    * \@param G_angle graph whose edges have angle-based weight
    * @param in_smoothed_v
    * @param normals: [OUT] normal of the point cloud with orientation corrected
    * @param k
    *
    * @return None
    */
void correct_normal_orientation(
    Util::ThreadPool& pool,
    const Tree& kdTree,
    const std::vector<Point>& in_smoothed_v,
    std::vector<Vec3>& normals,
    const int k)
{
    /// TODO: we rely on graph traversal here to flip normals as needed
    ///
    const auto [g_angle, sets] = [&]{
        SimpGraph g_angle_temp;
        AMGraph::NodeSet sets_temp;
        //sets.reserve(in_smoothed_v.size());

        for (int i = 0; i < in_smoothed_v.size(); i++) {
            sets_temp.insert(g_angle_temp.add_node());
        }

        //auto all_neighbors = calculate_neighbors(pool, vertices, kdTree, k);

        // Init angle based graph
        for (int i = 0; i < in_smoothed_v.size(); i++) {
            const Point vertex = in_smoothed_v[i];
            const Vec3 this_normal = normals[i];

            std::vector<NodeID> neighbors;
            std::vector<double> dists;
            knn_search(vertex, kdTree, k, neighbors, dists, false);

            for (const auto neighbor : neighbors) {
                if (g_angle_temp.find_edge(i, neighbor) != AMGraph::InvalidEdgeID)
                    continue;
                const Vec3 neighbor_normal = normals[neighbor];
                const float angle_weight = cal_angle_based_weight(this_normal, neighbor_normal);

                g_angle_temp.connect_nodes(i, neighbor, angle_weight);
            }
        }

        return std::make_tuple(g_angle_temp, sets_temp);
    }();


    const std::vector<AMGraph::NodeSet> components_vec = connected_components(g_angle, sets);

    // The number of components and their relative sizes is inconsistent, don't parallelize this
    for (const auto& i : components_vec) {
        SimpGraph mst_angle;
        NodeID root = *i.begin();
        minimum_spanning_tree(g_angle, root, mst_angle);

        auto visited_vertex = std::vector(g_angle.no_nodes(), Boolean{false});

        // This uses the MST to visit every node
        // Start from the root
        std::queue<NodeID> to_visit;
        to_visit.push(root);
        while (!to_visit.empty()) {
            const NodeID node_id = to_visit.front();
            to_visit.pop();

            visited_vertex[node_id].inner = true;
            const Vec3 this_normal = normals[node_id];
            for (auto neighbours = mst_angle.neighbors(node_id); NodeID vd : neighbours) {
                if (!visited_vertex[vd].inner) {
                    to_visit.push(vd);
                    const Vec3 neighbor_normal = normals[vd];
                    if (dot(this_normal, neighbor_normal) < 0) {
                        normals[vd] = -normals[vd];
                    }
                }
            }
        }
    }

}

/**
    * @brief Get the next neighbor information
    *
    * @param g: current graph
    * @param root: root vertex index
    * @param branch: current outgoing branch
    *
    * @return reference to the next neighbor struct
    */
const Neighbor& successor(const RSGraph& g, const NodeID& root, const NodeID& branch)
{
    const auto& u = g.m_vertices.at(root);
    const auto& v = g.m_vertices.at(branch);
    auto iter = u.ordered_neighbors.upper_bound(Neighbor(u, v, branch));
    if (iter == u.ordered_neighbors.end()) iter = u.ordered_neighbors.begin(); // Wrap around
    return (*iter); // This is honestly not good practice - ONLY modification of tree_id
}

/**
    * @brief Get last neighbor information
    *
    * @param g: current graph
    * @param root: root vertex index
    * @param branch: current outgoing branch
    *
    * @return reference to last neighbor struct
    */
const Neighbor& predecessor(const RSGraph& g, const NodeID& root, const NodeID& branch)
{
    const auto& u = g.m_vertices.at(root);
    const auto& v = g.m_vertices.at(branch);
    auto iter = u.ordered_neighbors.lower_bound({u, v, static_cast<uint>(branch)});
    if (iter == u.ordered_neighbors.begin()) iter = u.ordered_neighbors.end(); // Wrap around
    return (*(std::prev(iter)));
}

void init_face_loop_label(RSGraph& g)
{
    NodeID start_v = 0;
    NodeID last_vertex = start_v;
    int loop_step = 0;
    NodeID current_vertex = g.m_vertices[start_v].ordered_neighbors.begin()->v;
    std::vector<int> towers;
    do {
        auto& next_neighbor = predecessor(g, current_vertex, last_vertex);

        next_neighbor.tree_id = g.etf.accumulate();

        //std::cout << current_vertex << std::endl;
        last_vertex = current_vertex;
        current_vertex = next_neighbor.v;

        loop_step++;
    } while (current_vertex != g.m_vertices[start_v].ordered_neighbors.begin()->v || last_vertex != start_v);

    std::cout << "Loop step initialization finished after " + std::to_string(loop_step) + " steps." << std::endl;
}

/**
    * @brief Project a vector to a plane
    *
    * @param input: Vector to be projected
    * @param normal: normal to the plane
    *
    * @return projected Vector
    */
constexpr Vec3 projected_vector(const Vec3& input, const Vec3& normal)
{
    Vec3 normal_normed = normal;
    normal_normed.normalize();
    return input - CGLA::dot(input, normal_normed) * normal_normed;
}

/**
    * @brief Check if two segments are intersecting on both planes (defined by their normal) they belong
    *
    * @param mst: graph and vertex information
    * @param v1: 1st vertex of segment 1
    * @param v2: 2nd vertex of segment 1
    * @param v3: 1st vertex of segment 2
    * @param v4: 2nd vertex of segment 2
    *
    * @return if they are intersecting with each other
    */
bool isIntersecting(RSGraph& mst, NodeID v1, NodeID v2, NodeID v3, NodeID v4)
{
    Point p1 = mst.m_vertices[v1].coords;
    Point p2 = mst.m_vertices[v2].coords;
    Vec3 n1 = mst.m_vertices[v1].normal;
    Vec3 n2 = mst.m_vertices[v2].normal;
    Point midpoint_12 = p1 + (p2 - p1) / 2.;
    Vec3 normal_12 = (n1 + n2) / 2.;

    Point p3 = mst.m_vertices[v3].coords;
    Point p4 = mst.m_vertices[v4].coords;
    Vec3 n3 = mst.m_vertices[v3].normal;
    Vec3 n4 = mst.m_vertices[v4].normal;
    Point midpoint_34 = p3 + (p4 - p3) / 2.;
    Vec3 normal_34 = (n3 + n4) / 2.;

    // On the plane of edge 12
    {
        bool isIntersecting = true;
        Vec3 edge1 = p1 - midpoint_12;
        Vec3 edge2 = p3 - midpoint_12;
        Vec3 edge3 = p4 - midpoint_12;
        Vec3 proj_edge1 = projected_vector(edge1, normal_12);
        Vec3 proj_edge2 = projected_vector(edge2, normal_12);
        Vec3 proj_edge3 = projected_vector(edge3, normal_12);
        Vec3 pro1 = CGLA::cross(proj_edge2, proj_edge1);
        Vec3 pro2 = CGLA::cross(proj_edge3, proj_edge1);
        if (CGLA::dot(pro1, pro2) > 0)
            isIntersecting = false;
        if (isIntersecting) {
            edge1 = p3 - midpoint_34;
            edge2 = p1 - midpoint_34;
            edge3 = p2 - midpoint_34;
            proj_edge1 = projected_vector(edge1, normal_12);
            proj_edge2 = projected_vector(edge2, normal_12);
            proj_edge3 = projected_vector(edge3, normal_12);
            pro1 = CGLA::cross(proj_edge2, proj_edge1);
            pro2 = CGLA::cross(proj_edge3, proj_edge1);
            if (CGLA::dot(pro1, pro2) > 0)
                isIntersecting = false;
        }
        if (isIntersecting)
            return true;
    }

    // On the plane of edge 34
    if (true) {
        bool isIntersecting = true;
        Vec3 edge1 = p1 - midpoint_12;
        Vec3 edge2 = p3 - midpoint_12;
        Vec3 edge3 = p4 - midpoint_12;
        Vec3 proj_edge1 = projected_vector(edge1, normal_34);
        Vec3 proj_edge2 = projected_vector(edge2, normal_34);
        Vec3 proj_edge3 = projected_vector(edge3, normal_34);
        Vec3 pro1 = CGLA::cross(proj_edge2, proj_edge1);
        Vec3 pro2 = CGLA::cross(proj_edge3, proj_edge1);
        if (CGLA::dot(pro1, pro2) > 0)
            isIntersecting = false;
        if (isIntersecting) {
            edge1 = p3 - midpoint_34;
            edge2 = p1 - midpoint_34;
            edge3 = p2 - midpoint_34;
            proj_edge1 = projected_vector(edge1, normal_34);
            proj_edge2 = projected_vector(edge2, normal_34);
            proj_edge3 = projected_vector(edge3, normal_34);
            pro1 = CGLA::cross(proj_edge2, proj_edge1);
            pro2 = CGLA::cross(proj_edge3, proj_edge1);
            if (CGLA::dot(pro1, pro2) > 0)
                isIntersecting = false;
        }
        if (isIntersecting)
            return true;
    }
    return false;
}

/**
    * @brief Geometry check for connection
    *
    * @param mst: graph and vertex information
    * @param candidate: the edge to be examed
    * @param kdTree: kd-tree for knn query
    * @param tr_dist: distance container
    *
    * @return if the candidate pass the check
    */
bool geometry_check(RSGraph& mst, TEdge& candidate, Tree& kdTree)
{
    NodeID v1 = candidate.first;
    NodeID v2 = candidate.second;
    Point p1 = mst.m_vertices[v1].coords;
    Point p2 = mst.m_vertices[v2].coords;
    Vec3 n1 = mst.m_vertices[v1].normal;
    Vec3 n2 = mst.m_vertices[v2].normal;

    Vec3 mean_normal = (n1 + n2) / 2.;
    mean_normal.normalize();

    Point search_center = p1 + (p2 - p1) / 2.;
    float radius = (p2 - p1).length() / 2.;
    std::vector<NodeID> neighbors;
    std::vector<double> distance;
    NN_search(search_center, kdTree, float(radius * 3.), neighbors, distance, false);

    float query_radian1 = cal_radians_3d(p1 - search_center, mean_normal);
    float query_radian2 = cal_radians_3d(p2 - search_center, mean_normal);
    std::set<int> rejection_neighbor_set;
    for (int i = 0; i < neighbors.size(); i++) {
        if (neighbors[i] == v1 || neighbors[i] == v2)
            continue;
        if (CGLA::dot(mst.m_vertices[neighbors[i]].normal, mean_normal) > std::cos(60. / 180. * M_PI))
            rejection_neighbor_set.insert(neighbors[i]);
    }
    for (int i = 0; i < neighbors.size(); i++) {
        if (rejection_neighbor_set.find(neighbors[i]) ==
            rejection_neighbor_set.end())
            continue;
        //if (neighbors[i] == v1 || neighbors[i] == v2)
        //	continue;
        NodeID rejection_neighbor = neighbors[i];
        Point rej_neighbor_pos = mst.m_vertices[rejection_neighbor].coords;
        float min_radian, max_radian;

        for (auto& rej_neighbor_neighbor : mst.m_vertices[rejection_neighbor].ordered_neighbors) {
            if (rejection_neighbor_set.find(rej_neighbor_neighbor.v) ==
                rejection_neighbor_set.end())
                continue;

            if (false) {
                min_radian = cal_radians_3d(rej_neighbor_pos - search_center, mean_normal);

                Point rej_neighbor_neighbor_pos = mst.m_vertices[NodeID(rej_neighbor_neighbor.v)].coords;
                max_radian = cal_radians_3d(rej_neighbor_neighbor_pos - search_center,
                                            mean_normal);

                if (max_radian < min_radian) {
                    std::swap(max_radian, min_radian);
                }
                if (max_radian - min_radian > M_PI)
                    std::swap(max_radian, min_radian);

                bool is_in_between = false;
                if (max_radian < min_radian &&
                    (query_radian1 > min_radian || query_radian1 < max_radian))
                    is_in_between = true;
                if (max_radian > min_radian &&
                    (query_radian1 < max_radian && query_radian1 > min_radian))
                    is_in_between = true;
                if (max_radian < min_radian &&
                    (query_radian2 > min_radian || query_radian2 < max_radian))
                    is_in_between = true;
                if (max_radian > min_radian &&
                    (query_radian2 < max_radian && query_radian2 > min_radian))
                    is_in_between = true;

                if (is_in_between) {
                    Vec3 edge1 = p1 - rej_neighbor_pos;
                    Vec3 edge2 = p2 - rej_neighbor_pos;
                    Vec3 edge3 = rej_neighbor_neighbor_pos - rej_neighbor_pos;
                    Vec3 proj_edge1 = projected_vector(edge1, mean_normal);
                    Vec3 proj_edge2 = projected_vector(edge2, mean_normal);
                    Vec3 proj_edge3 = projected_vector(edge3, mean_normal);
                    Vec3 pro1 = CGLA::cross(proj_edge1, proj_edge3);
                    Vec3 pro2 = CGLA::cross(proj_edge2, proj_edge3);
                    if (CGLA::dot(pro1, pro2) <= 0)
                        return false;
                }
            } else {
                //if ((rejection_neighbor == 61 && rej_neighbor_neighbor.v == 345) ||
                //    (rejection_neighbor == 345 && rej_neighbor_neighbor.v == 61)) {
                //    std::cout << isIntersecting(mst, v1, v2, rejection_neighbor, rej_neighbor_neighbor.v) << std::endl;
                //    std::cout << isIntersecting(mst, v1, v2, rej_neighbor_neighbor.v, rejection_neighbor) << std::endl;
                //}
                bool result = isIntersecting(mst, v1, v2, rejection_neighbor, rej_neighbor_neighbor.v);
                if (result) {
                    //std::cout << rej_neighbor_neighbor.v << std::endl;
                    //std::cout << rejection_neighbor << std::endl;
                    return false;
                }
            }
        }
        rejection_neighbor_set.erase(neighbors[i]);
    }
    return true;
}

bool Vanilla_check(RSGraph& mst, TEdge& candidate, Tree& kdTree)
{
    const NodeID neighbor = candidate.second;
    const NodeID this_v = candidate.first;
    Vec3 this_normal = mst.m_vertices[this_v].normal;
    this_normal.normalize();
    Vec3 neighbor_normal = mst.m_vertices[neighbor].normal;
    neighbor_normal.normalize();

    // Topology check
    const auto this_v_tree = predecessor(mst, this_v, neighbor).tree_id;
    const auto neighbor_tree = predecessor(mst, neighbor, this_v).tree_id;

    if (!mst.etf.connected(this_v_tree, neighbor_tree)) {
        return false;
    }

    return geometry_check(mst, candidate, kdTree);
}

/**
    * @brief Find the common neighbors that two vertices are sharing
    *
    * @param neighbor: one vertex
    * @param root: the other vertex
    * @param share_neighbor: [OUT] common neighbors these two vertices share
    * @param g: tcurrent graph
    *
    * @return reference to last neighbor struct
    */
void find_common_neighbor(NodeID neighbor, NodeID root,
                          std::vector<NodeID>& share_neighbor, RSGraph& g)
{
    auto adj = g.neighbors(neighbor);
    std::set<NodeID> neighbor_neighbor(adj.begin(), adj.end());
    adj = g.neighbors(root);
    std::set<NodeID> vertex_neighbor(adj.begin(), adj.end());
    std::set_intersection(vertex_neighbor.begin(), vertex_neighbor.end(),
                          neighbor_neighbor.begin(), neighbor_neighbor.end(),
                          std::back_inserter(share_neighbor));
}

/**
    * @brief Get the neighbor information
    *
    * @param g: current graph
    * @param root: root vertex index
    * @param branch: the outgoing branch
    *
    * @return reference to the neighbor struct
    */
const Neighbor& get_neighbor_info(const RSGraph& g, const NodeID& root, const NodeID& branch)
{
    const auto& u = g.m_vertices.at(root);
    const auto& v = g.m_vertices.at(branch);
    auto iter = u.ordered_neighbors.lower_bound({u, v, static_cast<uint>(branch)});
    return (*iter);
}

void maintain_face_loop(RSGraph& g, const NodeID source, const NodeID target)
{
    auto this_v_tree = predecessor(g, source, target).tree_id;
    auto neighbor_tree = predecessor(g, target, source).tree_id;

    auto result = g.etf.insert(this_v_tree, neighbor_tree);
    auto u = result.first;
    auto v = result.second;
    get_neighbor_info(g, source, target).tree_id = u;
    get_neighbor_info(g, target, source).tree_id = v;
}

bool routine_check(const RSGraph& mst, const std::vector<NodeID>& triangle)
{
    const NodeID v1 = triangle[0];
    const NodeID v2 = triangle[1];
    const NodeID v3 = triangle[2];
    const Point p1 = mst.m_vertices[v1].coords;
    const Point p2 = mst.m_vertices[v2].coords;
    const Point p3 = mst.m_vertices[v3].coords;
    Vec3 n1 = mst.m_vertices[v1].normal;
    Vec3 n2 = mst.m_vertices[v2].normal;
    Vec3 n3 = mst.m_vertices[v3].normal;

    //Vector face_normal = normalize_vector(CGAL::cross_product(p2 - p1, p3 - p1));
    //bool isValid = (n1 * face_normal * (n2 * face_normal) < 0 ||
    //		n1 * face_normal * (n3 * face_normal) < 0);

    {
        const float len_ui = (p1 - p2).length();
        const float len_wi = (p3 - p2).length();
        const float len_uw = (p1 - p3).length();

        float max_value = std::acos(std::clamp<double>(
            dot((p3 - p2), (p1 - p2)) / len_ui /
            len_wi, -1, 1));
        float radian = std::acos(std::clamp<double>(
            dot((p2 - p1), (p3 - p1)) / len_ui /
            len_uw, -1, 1));
        if (radian > max_value)
            max_value = radian;
        radian = std::acos(std::clamp<double>(
            dot((p1 - p3), (p2 - p3)) / len_uw /
            len_wi, -1, 1));
        if (radian > max_value)
            max_value = radian;
        if (max_value > 175. / 180. * M_PI)
            return true;
    }


    if (mst.m_edges[mst.find_edge(v1, v3)].ref_time == 2 ||
        mst.m_edges[mst.find_edge(v2, v3)].ref_time == 2)
        return true;

    return false;
}


void add_face(RSGraph& G, const std::vector<NodeID>& item,
              std::vector<std::vector<NodeID>>& faces)
{
    NodeID v_i = item[0];
    NodeID v_u = item[1];
    NodeID v_w = item[2];

    // Maintain face exist for detecting holes
    {
        get_neighbor_info(G, v_i, v_u).faceExist = true;
        get_neighbor_info(G, v_u, v_w).faceExist = true;
        get_neighbor_info(G, v_w, v_i).faceExist = true;
    }

    G.m_edges[G.find_edge(v_u, v_w)].ref_time += 1;
    G.m_edges[G.find_edge(v_i, v_w)].ref_time += 1;
    G.m_edges[G.find_edge(v_u, v_i)].ref_time += 1;
    faces.push_back(item);
}

bool register_face(RSGraph& mst, NodeID v1, NodeID v2, std::vector<std::vector<NodeID>>& faces,
                   Tree& KDTree, float edge_length)
{
    const Point p1 = mst.m_vertices[v1].coords;
    const Point p2 = mst.m_vertices[v2].coords;

    if (mst.find_edge(v1, v2) != AMGraph::InvalidEdgeID)
        return false;

    std::vector<NodeID> share_neighbors;
    find_common_neighbor(v1, v2, share_neighbors, mst);
    if (share_neighbors.empty()) {
        mst.add_edge(v1, v2, edge_length);
        maintain_face_loop(mst, v1, v2);
        return true;
    }
    //if ((v1 == 30045 && v2 == 69461) || v1 == 69461 && v2 == 30045)
    //	std::cout << "debug here" << std::endl;

    const auto possible_root1 = predecessor(mst, v1, v2).v;
    float angle1 = cal_radians_3d(p1 - mst.m_vertices[possible_root1].coords,
                                  mst.m_vertices[possible_root1].normal, p2 - mst.m_vertices[possible_root1].coords);
    auto possible_root2 = predecessor(mst, v2, v1).v;
    float angle2 = cal_radians_3d(p2 - mst.m_vertices[possible_root2].coords,
                                  mst.m_vertices[possible_root2].normal, p1 - mst.m_vertices[possible_root2].coords);

    bool isValid = true;
    std::vector<std::vector<NodeID>> temp;
    for (auto v3 : share_neighbors) {
        std::vector<NodeID> triangle{v1, v2, v3};
        if (v3 == possible_root1 && angle1 < M_PI) {
            if (routine_check(mst, triangle)) {
                isValid = false;
                break;
            }
            if (successor(mst, v2, v1).v != v3) {
                isValid = false;
                break;
            }
            if (successor(mst, possible_root1, v2).v != v1) {
                isValid = false;
                break;
            }
            /*if (check_face_overlap(mst, triangle, KDTree, tr_dist)) {
                isValid = false;
                break;
            }*/
            temp.push_back(std::vector<NodeID>{v1, v3, v2});
        }

        if (v3 == possible_root2 && angle2 < M_PI) {
            if (routine_check(mst, triangle)) {
                isValid = false;
                break;
            }
            if (successor(mst, v1, v2).v != v3) {
                isValid = false;
                break;
            }
            if (successor(mst, possible_root2, v1).v != v2) {
                isValid = false;
                break;
            }
            /*if (check_face_overlap(mst, triangle, KDTree, tr_dist)) {
                isValid = false;
                break;
            }*/
            temp.push_back(std::vector<NodeID>{v1, v2, v3});
        }
    }

    if (temp.empty())
        isValid = false;

    if (isValid) {
        AMGraph::EdgeID added_edge = mst.add_edge(v1, v2,
                                                  edge_length);
        maintain_face_loop(mst, v1, v2);
        for (auto& face : temp) {
            add_face(mst, face, faces);
        }
    }

    return isValid;
}

/**
    * @brief Connect handle to raise the genus number
    *
    * @param smoothed_v: smoothed vertices of the point cloud
    * @param mst: graph and vertex information
    * @param kdTree: kd-tree for knn query
    * @param tr_dist: distance container
    * @param connected_handle_root: [OUT] log the connected handles
    * @param k: number of kNN search
    * @param isEuclidean: if to use Euclidean distance
    * @param step_thresh: step threshold for shortest distance path early stop
    *
    * @return None
    */
void connect_handle(
    const std::vector<Point>& smoothed_v, Tree& KDTree,
    RSGraph& mst, std::vector<NodeID>& connected_handle_root,
    int k,
    int step_thresh,
    bool isEuclidean)
{
    std::vector<NodeID> imp_node;
    int num = 0;
    int edge_num = 0;
    // Collect vertices w/ an open angle larger than pi
    {
        for (int i = 0; i < mst.no_nodes(); i++) {
            std::set<Neighbor>& neighbors = mst.m_vertices[i].ordered_neighbors;
            float last_angle = (--neighbors.end())->angle;
            float this_angle;

            for (auto& neighbor : neighbors) {
                this_angle = neighbor.angle;
                float angle_diff = this_angle - last_angle;
                if (angle_diff < 0)
                    angle_diff += 2 * M_PI;
                if (angle_diff > M_PI)
                    imp_node.push_back(i);
                last_angle = this_angle;
            }
        }
    }

    std::vector<NodeID> connect_p;
    std::vector<NodeID> to_connect_p;
    std::vector<uint> tree_id;
    std::vector<uint> to_tree_id;

    double last_dist = INFINITY;
    Point last_v(0., 0., 0.);
    for (const auto& this_v : imp_node) {
        // Point query = mst.m_vertices[this_v].coords;
        // Vec3 query_normal = mst.m_vertices[this_v].normal;
        std::vector<NodeID> neighbors;
        std::vector<double> dists;

        // Potential handle collection
        uint tree, to_tree;
        int validIdx = -1;

        last_dist += (smoothed_v[this_v] - last_v).length();
        knn_search(smoothed_v[this_v], KDTree, k, neighbors, dists, true);
        last_dist = dists[dists.size() - 1];
        last_v = smoothed_v[this_v];

        for (size_t i = 0; i < neighbors.size(); i++) {
            auto neighbor = neighbors[i];
            TEdge candidate(this_v, neighbor);
            if (mst.find_edge(this_v, neighbor) != AMGraph::InvalidEdgeID)
                continue;
            tree = mst.etf.representative((predecessor(mst, this_v, neighbor).tree_id));
            to_tree = mst.etf.representative(predecessor(mst, neighbor, this_v).tree_id);
            if (geometry_check(mst, candidate, KDTree) && tree != to_tree) {
                validIdx = i;
                break;
            }
        }
        // TODO: Check if any tree shares root, and return corresponding edges

        if (validIdx != -1) {
            connect_p.push_back(this_v);
            to_connect_p.push_back(neighbors[validIdx]);
            tree_id.push_back(tree);
            to_tree_id.push_back(to_tree);
        }
    }

    // Select one handle
    std::map<std::string, std::vector<int>> face_connections;
    for (int i = 0; i < connect_p.size(); i++) {
        uint tree = tree_id[i];
        uint to_tree = to_tree_id[i];
        if (to_tree > tree)
            std::swap(tree, to_tree);
        std::string key = std::to_string(tree) + "+" + std::to_string(to_tree);
        if (face_connections.find(key) == face_connections.end())
            face_connections[key] = std::vector<int>{i};
        else {
            face_connections[key].push_back(i);
        }
    }

    // Sort
    std::vector<m_face_pair> sorted_face;
    for (auto key = face_connections.begin(); key != face_connections.end(); ++key) {
        int length = face_connections[key->first].size();
        sorted_face.push_back(m_face_pair(length, key->first));
    }
    std::sort(sorted_face.begin(), sorted_face.end(), face_comparator);
    for (auto & i : sorted_face) {
        std::string key = i.second;
        std::vector<int> idx_vec = face_connections[key];
        if (idx_vec.size() <= 5)
            break;
        if (mst.exp_genus >= 0 && num >= mst.exp_genus)
            break;
        Point query;
        NodeID connected_neighbor, this_v;
        Edge added_edge;
        bool isFind = false;
        for (int idx : idx_vec) {
            this_v = connect_p[idx];
            query = mst.m_vertices[this_v].coords;
            connected_neighbor = to_connect_p[idx];
            std::vector<NodeID> path;
            int steps = find_shortest_path(mst, this_v, connected_neighbor, step_thresh, path);
            if (steps > step_thresh) {
                //if(steps >= 9){
                //std::cout << "This is connected" << std::endl;
                isFind = true;
                TEdge candidate(this_v, connected_neighbor);
                if (geometry_check(mst, candidate, KDTree)) {
                    Vec3 edge = query - mst.m_vertices[connected_neighbor].coords;
                    float Euclidean_dist = edge.length();
                    float projection_dist = cal_proj_dist(edge, mst.m_vertices[this_v].normal,
                                                          mst.m_vertices[connected_neighbor].normal);

                    if (isEuclidean) {
                        mst.add_edge(this_v, connected_neighbor, Euclidean_dist);
                        connected_handle_root.push_back(this_v);
                        connected_handle_root.push_back(connected_neighbor);
                    } else {
                        mst.add_edge(this_v, connected_neighbor, projection_dist);
                        connected_handle_root.push_back(this_v);
                        connected_handle_root.push_back(connected_neighbor);
                    }

                    edge_num++;
                    //fs::path out_edge_path("C:/Projects_output/letters/edge_" + std::to_string(edge_num) + ".obj");
                    //export_continuous_edges(mst, path, out_edge_path);
                }
            }
        }
        if (isFind) {
            num++;
        }
    }

    std::cout << "Handle Connection done :)" << std::endl;
    std::cout << std::to_string(num) << " pairs of faces are connected." << std::endl;
    std::cout << std::to_string(edge_num) << " edges are connected." << std::endl;
}

bool explore(RSGraph& G, int i, m_priority_queue& queue,
             std::unordered_set<std::string>& faces_in_queue, float avg_edge_length
             , std::vector<float>& length_thresh)
{
    const NodeID v_i = i;
    bool isFound = false;
    for (auto& neighbor : G.m_vertices[i].ordered_neighbors) {
        NodeID v_u = neighbor.v;
        NodeID v_w = successor(G, i, v_u).v;

        //if ((v_u == 1751375 && v_w == 1888150) ||
        //	(v_w == 1751375 && v_u == 1888150))
        //	std::cout << "debug here" << std::endl;

        //if ((v_u == 33020 && v_w == 30023) || (v_w == 33020 && v_u == 30023)) {
        //	std::cout << "debug here" << std::endl;
        //}

        Point w_pos = G.m_vertices[v_w].coords;
        Point u_pos = G.m_vertices[v_u].coords;
        Point i_pos = G.m_vertices[v_i].coords;
        Vec3 i_normal = G.m_vertices[v_i].normal;
        Vec3 u_normal = G.m_vertices[v_u].normal;
        Vec3 w_normal = G.m_vertices[v_w].normal;
        float angle = cal_radians_3d(w_pos - i_pos, i_normal,
                                     u_pos - i_pos);
        bool isLargerThanPi = angle < M_PI;
        std::vector<NodeID> face_vector{v_i, v_u, v_w};
        if (v_u != v_w && isLargerThanPi) {
            if (G.find_edge(v_u, v_w) == AMGraph::InvalidEdgeID) {
                float score = (G.m_vertices[v_u].coords - G.m_vertices[v_w].coords).length();
                if (!G.isEuclidean) {
                    score = cal_proj_dist(G.m_vertices[v_u].coords - G.m_vertices[v_w].coords,
                                          u_normal, w_normal);
                }
                if (score > length_thresh[v_u] || score > length_thresh[v_w])
                    continue;
                if (score >= 0) {
                    std::pair<std::vector<NodeID>, float> queue_item(face_vector, score);
                    queue.push(queue_item);
                    isFound = true;
                }
            }
        }
    }

    return isFound;
}

/**
    * @brief Calculate the local surface normal (averaged direction of normals of 3 vertices in the triangle)
    *
    * @return local surface normal
    */
Vec3 triangle_mean_normal(const Vec3& normal1, const Vec3& normal2, const Vec3& normal3)
{
    Vec3 normal1_norm = normal1;
    normal1_norm.normalize();
    Vec3 normal2_norm = normal2;
    normal2_norm.normalize();
    Vec3 normal3_norm = normal3;
    normal3_norm.normalize();
    Vec3 output = normal1_norm + normal2_norm + normal3_norm;
    output.normalize();
    return output;
}

bool check_branch_validity(RSGraph& G, NodeID root, NodeID branch1, NodeID branch2)
{
    Point pos_i = G.m_vertices[root].coords;
    Point pos_u = G.m_vertices[branch1].coords;
    Point pos_w = G.m_vertices[branch2].coords;
    Vec3 normal_i = G.m_vertices[root].normal;
    Vec3 normal_u = G.m_vertices[branch1].normal;
    Vec3 normal_w = G.m_vertices[branch2].normal;

    // Option 1
    std::vector<Point> triangle_pos{pos_i, pos_u, pos_w};
    //Vector face_normal = CGAL::cross_product(pos_i - pos_u, pos_w - pos_u);
    Vec3 face_normal = triangle_mean_normal(normal_i, normal_u, normal_w);

    float angle_thresh = 0. / 180. * M_PI;

    // Check u's RS validity
    float this_radian = cal_radians_3d(pos_w - pos_u, normal_u);
    auto former = predecessor(G, branch1, branch2);
    auto next = successor(G, branch1, branch2);
    if (G.isFinal) {
        bool isValid = false;
        if (next.v == root) {
            float diff = next.angle - this_radian;
            if (diff < 0)
                diff += 2 * M_PI;
            if (diff < M_PI)
                isValid = true;
        }
        if (former.v == root) {
            float diff = -former.angle + this_radian;
            if (diff < 0)
                diff += 2 * M_PI;
            if (diff < M_PI)
                isValid = true;
        }
        if (!isValid)
            return false;
    } else {
        float diff = next.angle - this_radian;
        if (diff < 0)
            diff += 2 * M_PI;
        if (next.v != root || diff > M_PI) {
            return false;
        }
    }

    // Thresh on angle
    {
        float diff_angle_thresh = this_radian - former.angle;
        if (diff_angle_thresh < 0)
            diff_angle_thresh += M_PI * 2.;
        if (diff_angle_thresh < angle_thresh)
            return false;
    }

    //Check w
    this_radian = cal_radians_3d(pos_u - pos_w, normal_w);
    former = predecessor(G, branch2, branch1);
    next = successor(G, branch2, branch1);
    if (G.isFinal) {
        bool isValid = false;
        if (next.v == root) {
            float diff = next.angle - this_radian;
            if (diff < 0)
                diff += 2 * M_PI;
            if (diff < M_PI)
                isValid = true;
        }
        if (former.v == root) {
            float diff = -former.angle + this_radian;
            if (diff < 0)
                diff += 2 * M_PI;
            if (diff < M_PI)
                isValid = true;
        }
        if (!isValid)
            return false;
    } else {
        float diff = -former.angle + this_radian;
        if (diff < 0)
            diff += 2 * M_PI;
        if (former.v != root || diff > M_PI) {
            return false;
        }
    }

    // Thresh on angle
    {
        float diff_angle_thresh = -this_radian + next.angle + 2. * M_PI;
        if (diff_angle_thresh < angle_thresh)
            return false;
    }
    return true;
}


// TODO: validity of what?
bool check_validity(RSGraph& G, std::pair<std::vector<NodeID>, float>& item,
                    const Tree& KDTree, bool isFaceloop, bool isFinalize)
{
    int i = item.first[0];
    int u = item.first[1];
    int w = item.first[2];
    NodeID v_i = i;
    NodeID v_u = u;
    NodeID v_w = w;
    Point pos_i = G.m_vertices[v_i].coords;
    Point pos_u = G.m_vertices[v_u].coords;
    Point pos_w = G.m_vertices[v_w].coords;
    Vec3 normal_i = G.m_vertices[v_i].normal;
    Vec3 normal_u = G.m_vertices[v_u].normal;
    Vec3 normal_w = G.m_vertices[v_w].normal;

    //if (i == 388237)
    //	std::cout << u << " " << w << std::endl;

    if (G.find_edge(v_u, v_w) != AMGraph::InvalidEdgeID)
        return false;

    // Non-manifold edge check
    if (G.m_edges[G.find_edge(v_i, v_u)].ref_time == 2 ||
        G.m_edges[G.find_edge(v_i, v_w)].ref_time == 2)
        return false;

    // Check this rotation system
    bool isValid = (successor(G, v_i, v_u).v == v_w);
    float angle = cal_radians_3d(pos_w - pos_i, normal_i, pos_u - pos_i);
    if (angle > M_PI)
        isValid = false;

    if (!isValid)
        return false;

    if (!isFinalize) {
        // Check the rotation system's validity of branch nodes
        if (!check_branch_validity(G, v_i, v_u, v_w)) {
            //std::cout << v_i << std::endl;
            return false;
        }
    }

    // Check face overlap
    // TODO: why is face overlap disabled
    /*if (check_face_overlap(G, item.first, KDTree, tr_dist))
        return false;*/

    return true;
}

void checkAndForce(NodeID v_u, NodeID v_w, RSGraph& G, m_priority_queue& queue,
                   std::vector<float>& length_thresh)
{
    std::vector<NodeID> check_v{v_u, v_w};
    for (int i = 0; i < check_v.size(); i++) {
        NodeID v_i = check_v[i];
        std::vector<int> connection_count(G.m_vertices[v_i].ordered_neighbors.size(), 0);
        std::vector<NodeID> neighbor_ids;
        std::set<NodeID> neighbor_ids_set;

        for (auto& neighbor : G.m_vertices[v_i].ordered_neighbors) {
            neighbor_ids_set.insert(neighbor.v);
            neighbor_ids.push_back(neighbor.v);
        }

        int idx = 0;
        for (auto& neighbor : G.m_vertices[v_i].ordered_neighbors) {
            for (auto& neighbor_neighbor : G.m_vertices[neighbor.v].ordered_neighbors) {
                if (neighbor_ids_set.find(neighbor_neighbor.v) != neighbor_ids_set.end()) {
                    connection_count[idx] += 1;
                }
            }
            idx++;
        }

        std::vector<NodeID> not_full;
        for (int i = 0; i < connection_count.size(); i++) {
            if (connection_count[i] < 2) {
                not_full.push_back(i);
            }
        }

        if (not_full.size() == 2) {
            NodeID v_u = neighbor_ids[not_full[0]];
            NodeID v_w = neighbor_ids[not_full[1]];
            Vec3 u_normal = G.m_vertices[v_u].normal;
            Vec3 w_normal = G.m_vertices[v_w].normal;
            float score = (G.m_vertices[v_u].coords - G.m_vertices[v_w].coords).length();
            if (!G.isEuclidean) {
                score = cal_proj_dist(G.m_vertices[v_u].coords - G.m_vertices[v_w].coords,
                                      u_normal, w_normal);
            }
            if (score > length_thresh[v_u] || score > length_thresh[v_w])
                return;
            std::vector<NodeID> face_vector{v_i, neighbor_ids[not_full[0]], neighbor_ids[not_full[1]]};
            std::pair<std::vector<NodeID>, float> queue_item(face_vector, -1);
            queue.push(queue_item);
            break;
        }
    }
}

void triangulate(
    std::vector<std::vector<NodeID>>& faces,
    RSGraph& G,
    const Tree& KDTree,
    bool isFaceLoop,
    bool isEuclidean,
    std::vector<float>& length_thresh,
    std::vector<NodeID>& connected_handle_root,
    bool isFinalize = false)
{
    std::unordered_set<std::string> faces_in_queue;
    std::unordered_set<NodeID> to_visit;

    m_priority_queue queue;

    float avg_edge_length = G.total_edge_length / G.no_edges();
    //std::cout << G.no_edges() << std::endl;

    // Init priority queue
    for (int i = 0; i < G.no_nodes(); i++) {
        to_visit.insert(i);
    }

    for (int i = 0; i < connected_handle_root.size(); i++) {
        bool result = explore(G, connected_handle_root[i], queue, faces_in_queue, avg_edge_length, length_thresh);
        to_visit.erase(connected_handle_root[i]);
    }

    std::cout << "Global init done :)" << std::endl;

    int loop_time = 0;

    while (to_visit.size() > 0) {
        while (!queue.empty()) {
            loop_time += 1;
            std::pair<std::vector<NodeID>, float> item = queue.top();
            queue.pop();

            //std::cout << item.first[1] << " " << item.first[2] << std::endl;
            //std::cout << item.second << std::endl;
            //if (item.first[1] == 345 && item.first[2] == 426)
            //    std::cout << item.second << std::endl;

            if (item.second >= 0) {
                // Validity check
                bool isValid = check_validity(G, item, KDTree, isFaceLoop, isFinalize);
                if (!isValid)
                    continue;
            }

            // Add the edge
            NodeID v_i = item.first[0];
            NodeID v_u = item.first[1];
            NodeID v_w = item.first[2];
            Point pos_i = G.m_vertices[v_i].coords;
            Point pos_u = G.m_vertices[v_u].coords;
            Point pos_w = G.m_vertices[v_w].coords;
            Vec3 normal_i = G.m_vertices[v_i].normal;
            Vec3 normal_u = G.m_vertices[v_u].normal;
            Vec3 normal_w = G.m_vertices[v_w].normal;

            float dist = (pos_u - pos_w).length();
            Vec3 edge = pos_u - pos_w;
            float Euclidean_dist = (edge).length();
            float projection_dist = cal_proj_dist(edge,
                                                  G.m_vertices[v_u].normal, G.m_vertices[v_w].normal);

            if (G.find_edge(v_u, v_w) == AMGraph::InvalidEdgeID) {
                Edge added_edge;
                if (isEuclidean)
                    G.add_edge(v_u, v_w, Euclidean_dist);
                else
                    G.add_edge(v_u, v_w, projection_dist);

                avg_edge_length = G.total_edge_length / G.no_edges();

                add_face(G, item.first, faces);
            } else
                continue;

            // Deal with incident triangles
            {
                std::vector<NodeID> share_neighbors;
                find_common_neighbor(v_u, v_w, share_neighbors, G);
                for (int idx = 0; idx < share_neighbors.size(); idx++) {
                    NodeID incident_root = share_neighbors[idx];
                    if (incident_root == v_i)
                        continue;
                    std::vector<NodeID> face{incident_root, v_w, v_u};

                    // Non-manifold edge check
                    int time1 = G.m_edges[G.find_edge(incident_root, v_u)].ref_time;
                    int time2 = G.m_edges[G.find_edge(incident_root, v_w)].ref_time;
                    int time3 = G.m_edges[G.find_edge(v_u, v_w)].ref_time;
                    if (time1 == 2 || time2 == 2 || time3 == 2)
                        continue;

                    add_face(G, face, faces);
                }
            }

            to_visit.erase(v_u);
            to_visit.erase(v_w);

            // Explore and sanity check
            bool isFound = false;
            bool result = explore(G, v_u, queue, faces_in_queue, avg_edge_length, length_thresh);
            isFound = isFound || result;
            result = explore(G, v_w, queue, faces_in_queue, avg_edge_length, length_thresh);
            isFound = isFound || result;

            if (isFinalize) {
                if ((!isFound)) {
                    checkAndForce(v_u, v_w, G, queue, length_thresh);
                }
            }
        }

        if (!to_visit.empty()) {
            NodeID pick = *to_visit.begin();
            to_visit.erase(pick);
            bool result = explore(G, pick, queue, faces_in_queue, avg_edge_length, length_thresh);
            // TODO: result is dropped?
        }
    }
}

/**
 * @brief Build minimum spanning tree (MST)
 *
 * @param out_mst: [OUT] constructed MST
 * @param g: connection information of the mst
 * @param root root node
 * @param isEuclidean: if to use Euclidean distance
 * @param vertices: coordinates of the point cloud
 * @param normals: normal of the point cloud
 *
 * @return None
 */
void build_mst(
    const SimpGraph& g,
    const NodeID root,
    RSGraph& out_mst,
    const std::vector<Vec3>& normals,
    const std::vector<Point>& vertices,
    const bool isEuclidean)
{
    //SimpGraph mst_temp;
    //for (int i = 0; i < mst.no_nodes(); i++) {
    //    mst_temp.add_node();
    //}
    RSGraph temp;
    minimum_spanning_tree(g, root, temp, normals, vertices, isEuclidean);

    // Fix strong ambiguous points
    if (!isEuclidean) {
        for (int i = 0; i < temp.m_edges.size(); i++) {
            NodeID source = temp.m_edges[i].source;
            NodeID target = temp.m_edges[i].target;
            Vec3 normal1 = temp.m_vertices[source].normal;
            normal1 /= normal1.length();
            Vec3 normal2 = temp.m_vertices[target].normal;
            normal2 /= normal2.length();
            Point pos1 = temp.m_vertices[source].coords;
            Point pos2 = temp.m_vertices[target].coords;
            if (temp.valence(source) >= 2 && temp.valence(target) >= 2)
                continue;
            Vec3 edge = pos2 - pos1;

            Vec3 normal_sum = normal1 + normal2;
            float cos_angle = std::abs(dot(edge, normal_sum / normal_sum.length() / edge.length()));
            if (cos_angle > std::cos(10. / 180. * M_PI)) {
                NodeID leaf, parent;
                if (temp.valence(source) == 1) {
                    temp.m_vertices[source].normal = temp.m_vertices[target].normal;
                    parent = target;
                    leaf = source;
                } else {
                    temp.m_vertices[target].normal = temp.m_vertices[source].normal;
                    parent = source;
                    leaf = target;
                }

                for (auto neighbors = temp.neighbors(parent); const auto neighbor : neighbors) {
                    if (temp.m_vertices[neighbor].normal_rep == -1) {
                        temp.m_vertices[neighbor].normal_rep = parent;
                    } else {
                        // Collision!
                        temp.m_vertices[neighbor].normal_rep = -2;
                    }
                }
            }
        }
        for (int i = 0; i < temp.m_vertices.size(); i++) {
            if (temp.m_vertices[i].normal_rep >= 0)
                temp.m_vertices[i].normal = temp.m_vertices[temp.m_vertices[i].normal_rep].normal;
        }
    }

    // Build corrected MST
    for (int i = 0; i < temp.m_vertices.size(); i++) {
        out_mst.add_node(temp.m_vertices[i].coords, temp.m_vertices[i].normal);
    }
    for (int i = 0; i < temp.m_edges.size(); i++) {
        const Vec3 edge = out_mst.m_vertices[temp.m_edges[i].source].coords -
            out_mst.m_vertices[temp.m_edges[i].target].coords;

        const double distance =
            (isEuclidean) ?
                edge.length() :
                cal_proj_dist(edge,out_mst.m_vertices[temp.m_edges[i].source].normal,
                                              out_mst.m_vertices[temp.m_edges[i].target].normal);

        [[unlikely]]
        if (std::isnan(distance))
            std::cerr << "debug" << std::endl;

        out_mst.add_edge(temp.m_edges[i].source,
                 temp.m_edges[i].target, distance);
    }
}

template <typename Collection>
auto indices_from(const Collection& collection) -> std::vector<NodeID>
{
    const auto indices = [&collection] {
        std::vector<NodeID> temp(collection.size());
        std::iota(temp.begin(), temp.end(), 0);
        return temp;
    }();
    return indices;
}

[[nodiscard]]
auto estimate_normals_and_smooth(
    Util::ThreadPool& pool,
    std::vector<Point>& org_vertices,
    std::vector<Vec3>& org_normals,
    const std::vector<NodeID>& indices,
    const RsROpts& opts) -> std::vector<Point>
{

    // Insert the number_of_data_points into the tree
    const auto tree_before_remove = build_KDTree(org_vertices, indices);
    remove_duplicate_vertices(pool, org_vertices, org_normals, tree_before_remove, opts.k);
    const int smoothing_size = std::max(static_cast<int>(static_cast<double>(org_vertices.size()) / 2000.), 192);

    const auto kdTree = build_KDTree(org_vertices, indices);
    auto neighbors =
        calculate_neighbors(pool, org_vertices, kdTree, smoothing_size);
    if (opts.normals_included) {
        normalize_normals(org_normals);
    } else {
        estimate_normal_no_normals_memoized(pool, org_vertices, neighbors, org_normals);
    }

    std::cout << "Start first round smoothing ..." << std::endl;
    std::vector<Point> in_smoothed_v;
    in_smoothed_v.reserve(org_vertices.size());
    if (!opts.is_euclidean)
        weighted_smooth(pool, org_vertices, org_normals, neighbors, in_smoothed_v);
    else
        // Note: this copies the entire vertices
        in_smoothed_v = org_vertices;


    if (!opts.normals_included) {
        const auto temp_tree1 = build_KDTree(in_smoothed_v, indices);
        neighbors = calculate_neighbors(pool, org_vertices, temp_tree1, smoothing_size);
        estimate_normal_no_normals_memoized(pool, in_smoothed_v, neighbors, org_normals);
    }

    // Another round of smoothing
    if (!opts.is_euclidean) {
        std::cout << "Start second round smoothing ..." << std::endl;

        std::vector<Point> temp;
        temp.reserve(in_smoothed_v.size());
        std::swap(temp, in_smoothed_v);
        in_smoothed_v.clear();
        weighted_smooth(pool, temp, org_normals, neighbors, in_smoothed_v);

        if (!opts.normals_included) {
            const Tree temp_tree2 = build_KDTree(in_smoothed_v, indices);
            neighbors = calculate_neighbors(pool, org_vertices, temp_tree2, smoothing_size, std::move(neighbors));
            estimate_normal_no_normals_memoized(pool, in_smoothed_v, neighbors, org_normals);
        }
    }
    assert(org_vertices.size() == org_normals.size());
    return in_smoothed_v;
}

struct Components {
    std::vector<std::vector<Point>> vertices;
    std::vector<std::vector<Point>> smoothed_v;
    std::vector<std::vector<Vec3>> normals;

    explicit Components(
        std::vector<std::vector<Point>>&& vertices,
        std::vector<std::vector<Point>>&& smoothed_v,
        std::vector<std::vector<Vec3>>&& normals) noexcept
    : vertices { std::move(vertices) }, smoothed_v { std::move(smoothed_v) }, normals { std::move(normals) } {}
    Components() = delete;
    Components(Components& other) = delete;
};

/**
    * @brief Find the number of connected components and separate them
    *
    * @param org_vertices: vertices of the point cloud
    * @param in_smoothed_v: smoothed vertices of the point cloud
    * @param org_normals: normal of the point cloud vertices
    * @param kdTree: kd-tree for neighbor query
    * @param cross_conn_thresh: angle threshold to avoid connecting vertices on different surface
    * @param outlier_thresh: threshold to remove outlier
    * @param k
    * @param isEuclidean
    *
    * @return None
    */
[[nodiscard]]
auto split_components(
    Util::ThreadPool& pool,
    const Tree& kdTree,
    std::vector<Point>&& org_vertices,
    std::vector<Vec3>&& org_normals,
    std::vector<Point>&& in_smoothed_v,
    const RsROpts& opts)
    -> Components
{
    assert(org_vertices.size() == org_normals.size());
    assert(org_vertices.size() == in_smoothed_v.size());
    std::vector<std::vector<Point>> component_vertices;
    std::vector<std::vector<Point>> component_smoothed_v;
    std::vector<std::vector<Vec3>> component_normals;

    // Identifies clusters of vertices which are reconstructed to disparate meshes
    const std::vector<Point>& vertices = org_vertices;
    const std::vector<Point>& smoothed_v = in_smoothed_v;
    const std::vector<Vec3>& normals = org_normals;
    double cross_conn_thresh = opts.theta;
    double outlier_thresh = opts.r;
    int k = opts.k;
    bool isEuclidean = opts.is_euclidean;
    assert(vertices.size() == smoothed_v.size());
    assert(vertices.size() == normals.size());
    double avg_edge_length = 0;
    // TODO: can't we cache this?
    AMGraph::NodeSet sets;
    SimpGraph components;
    for (int i = 0; i < vertices.size(); i++) {
        sets.insert(components.add_node());
    }
    NodeID this_idx = 0;
    // Construct graph
    for (auto& vertex : smoothed_v) {

        std::vector<NodeID> neighbors;
        std::vector<double> neighbor_distance;
        knn_search(vertex, kdTree, k, neighbors, neighbor_distance, true);

        // Filter out cross connection
        {
            std::vector<NodeID> temp;
            Vec3 this_normal = normals[this_idx];
            for (auto idx : neighbors) {
                Vec3 neighbor_normal = normals[idx];
                double cos_theta = dot(this_normal, neighbor_normal) /
                    this_normal.length() / neighbor_normal.length();
                double cos_thresh = cos(cross_conn_thresh / 180. * M_PI);
                if (isEuclidean)
                    cos_thresh = 0.;
                if (cos_theta >= cos_thresh) {
                    temp.push_back(idx);
                }
            }
            neighbors = std::move(temp);
        }

        for (int i = 0; i < neighbors.size(); i++) {
            NodeID idx = neighbors[i];
            double length = neighbor_distance[i];

            if (this_idx == idx)
                continue;

            avg_edge_length += length;

            for (int j = 0; j < k; j++) {
                if (components.find_edge(this_idx, idx) != AMGraph::InvalidEdgeID)
                    continue;
                components.connect_nodes(this_idx, idx);
            }
        }
        this_idx++;
    }
    const double thresh_r = avg_edge_length / static_cast<double>(components.no_edges()) * outlier_thresh;
    // Remove Edges Longer than the threshold
    std::vector<std::pair<NodeID, NodeID>> edge_rm_v;
    for (NodeID i = 0; i < components.no_nodes(); i++) {
        for (auto& edges = components.edges(i); const auto& pair : edges) {
            NodeID vertex1 = i;
            NodeID vertex2 = pair.first;
            double edge_length = (vertices[vertex1] - vertices[vertex2]).length();

            if (edge_length > thresh_r) {
                edge_rm_v.emplace_back(vertex1, vertex2);
            }
        }
    }
    for (auto& [fst, snd] : edge_rm_v) {
        components.disconnect_nodes(fst, snd);
    }
    // Find Components
    std::vector<AMGraph::NodeSet> components_vec;
    components_vec = connected_components(components, sets);
    const auto num = components_vec.size();
    std::cout << "The input contains " << num << " connected components." << std::endl;
    // Valid Components and create new vectors for components
    auto threshold = std::min<size_t>(vertices.size(), 100);
    for (auto& component : components_vec) {
        if (component.size() >= threshold) {
            std::vector<Vec3> this_normals;
            this_normals.reserve(component.size());
            std::vector<Point> this_vertices;
            this_vertices.reserve(component.size());
            std::vector<Point> this_smoothed_v;
            this_smoothed_v.reserve(component.size());
            for (const auto& element : component) {
                this_vertices.push_back(vertices[element]);
                this_smoothed_v.push_back(smoothed_v[element]);
                this_normals.push_back(normals[element]);
            }
            component_normals.emplace_back(std::move(this_normals));
            component_vertices.emplace_back(std::move(this_vertices));
            component_smoothed_v.emplace_back(std::move(this_smoothed_v));
        }
    }
    std::cout << component_vertices.size() << " of them will be reconstructed." << std::endl;

 return Components(
        std::move(component_vertices),
        std::move(component_smoothed_v),
        std::move(component_normals)
    );
}

auto component_to_manifold(
    Util::ThreadPool& pool,
    const RsROpts& opts,
    const std::vector<Point>&& vertices,
    const std::vector<Vec3>&& normals,
    const std::vector<Point>&& smoothed_v,
    const std::vector<NodeID>& indices) -> ::HMesh::Manifold
{
    std::vector<std::vector<NodeID>> faces;
    // Insert the number_of_data_points in the tree
    Tree kdTree = build_KDTree(smoothed_v, indices);

    std::cout << "Init mst" << std::endl;

    // Initial Structure
    RSGraph mst;
    /*mst.init(vertices.size());*/
    std::vector<TEdge> full_edges;
    std::vector<m_Edge_length> edge_length;
    std::vector<float> connection_max_length(vertices.size(), 0.);
    std::vector<float> pre_max_length(vertices.size(), 0.);
    mst.isEuclidean = opts.is_euclidean;
    mst.exp_genus = opts.genus;
    {
        SimpGraph g;
        init_graph(smoothed_v, smoothed_v, normals,
                   kdTree, g, connection_max_length,
                   pre_max_length, opts.theta, opts.k, opts.is_euclidean);

        // Generate MST
        build_mst(g, 0, mst, normals, smoothed_v, opts.is_euclidean);

        // Edge arrays and sort
        for (NodeID node : g.node_ids()) {
            for (NodeID node_neighbor : g.neighbors(node)) {
                if (node < node_neighbor) {
                    Vec3 edge = smoothed_v[node] - smoothed_v[node_neighbor];
                    double len = edge.length();

                    if (!opts.is_euclidean) {
                        len = cal_proj_dist(edge, normals[node], normals[node_neighbor]);
                    }

                    if (len > pre_max_length[node] ||
                        len > pre_max_length[node_neighbor])
                        continue;
                    edge_length.emplace_back(len, full_edges.size());
                    full_edges.emplace_back(node, node_neighbor);
                }
            }
        }
        std::ranges::sort(edge_length.begin(), edge_length.end(), edge_comparator);
    }

    // Initialize face loop label
    mst.etf.reserve(6 * vertices.size() - 11);
    init_face_loop_label(mst);


    // Vanilla MST imp
    // Edge connection
    // TODO: progress was tracked here
    for (auto & i : edge_length) {
        unsigned int edge_idx = i.second;
        TEdge this_edge = full_edges[edge_idx];

        if (mst.find_edge(this_edge.first, this_edge.second) != AMGraph::InvalidEdgeID)
            continue;

        // TODO: isAdded not checked?
        if (bool isValid = Vanilla_check(mst, this_edge, kdTree)) {
            bool isAdded = register_face(mst, this_edge.first, this_edge.second, faces, kdTree,
                                         i.first);
        }
    }
    std::cout << std::endl;

    // Create handles & Triangulation
    if (opts.genus != 0) {
        mst.isFinal = true;
        std::vector<NodeID> connected_handle_root;
        connect_handle(smoothed_v, kdTree, mst, connected_handle_root, opts.k, opts.n, opts.is_euclidean);
        bool isFaceLoop = false;
        triangulate(faces, mst, kdTree, isFaceLoop, opts.is_euclidean, connection_max_length, connected_handle_root);
    }

    ::HMesh::Manifold res;
    // Extract vertex position
    std::vector<double> pos;
    for (int i = 0; i < mst.no_nodes(); i++) {
        pos.push_back(vertices[i][0]);
        pos.push_back(vertices[i][1]);
        pos.push_back(vertices[i][2]);
    }

    std::vector<int> mesh_faces(faces.size(), 3);
    std::vector<int> flattened_face;
    for (auto& face : faces) {
        flattened_face.push_back(face[0]);
        flattened_face.push_back(face[1]);
        flattened_face.push_back(face[2]);
    }

    // TODO: This takes a raw pointer which is probably unwanted
    ::HMesh::build(res, mst.no_nodes(), &pos[0], faces.size(),
                   &mesh_faces[0], &flattened_face[0]);

    return std::move(res);
}

auto point_cloud_to_mesh(
    const std::vector<Point>& vertices,
    const std::vector<Vec3>& normals,
    const RsROpts& opts) -> ::HMesh::Manifold
{
    auto opts2 = opts;
    ::HMesh::Manifold output;
    Timer timer;

    timer.start("Whole process");

    auto vertices_copy = vertices;
    auto normals_copy = normals;
    Util::ThreadPool pool(15);
    if (normals.empty()) {
        opts2.normals_included = false;
    } else {
        assert(vertices.size() == normals.size());
    }

    // Estimate normals & orientation & weighted smoothing
    timer.start("Estimate and smooth normals");
    const auto indices = indices_from(vertices_copy);

    std::vector<Point> in_smoothed_v = estimate_normals_and_smooth(pool, vertices_copy, normals_copy, indices, opts2);
    timer.end("Estimate and smooth normals");

    timer.start("Correct normal orientation");
    const Tree kdTree = build_KDTree(in_smoothed_v, indices);
    std::cout << "correct normal orientation\n";

    if (!opts2.normals_included) {
        correct_normal_orientation(pool, kdTree, in_smoothed_v, normals_copy, opts.k);
    }
    timer.end("Correct normal orientation");

    // Find components
    timer.start("Split components");
    std::cout << "find components\n";
    auto [component_vertices,
        component_smoothed_v,
        component_normals] =
            split_components(pool, kdTree, std::move(vertices_copy), std::move(normals_copy), std::move(in_smoothed_v), opts2);
    timer.end("Split components");
    // There is no guarantee that there is more than one component, and components can
    // be highly non-uniform in terms of how many primitives they have. That means we cannot
    // rely on this loop for good parallelization opportunities.
    timer.start("Algorithm");
    for (size_t component_id = 0; component_id < component_vertices.size(); component_id++) {
        std::cout << "Reconstructing component " + std::to_string(component_id) + " ...\n";

        std::vector<Point> vertices_of_this = std::move(component_vertices[component_id]);
        std::vector<Vec3>  normals_of_this = std::move(component_normals[component_id]);
        std::vector<Point> smoothed_v_of_this = std::move(component_smoothed_v[component_id]);

        auto res = component_to_manifold(
            pool,
            opts2,
            std::move(vertices_of_this),
            std::move(normals_of_this),
            std::move(smoothed_v_of_this),
            indices);
        output.merge(res);
    }
    timer.end("Algorithm");
    timer.end("Whole process");

    const std::string line(40, '=');
    std::cout << line << "\n\n";

    timer.show();

    return output;
}

} // namespace GEL::HMesh::RsR
