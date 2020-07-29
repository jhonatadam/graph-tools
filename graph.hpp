#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <generator.hpp>

using namespace std;

struct None {};

template <class T_edata = void*,
          class T_vdata = void*,
          class T_vertex = long>
class Graph
{
public:
    struct Edge;

    using Vertex = T_vertex;
    using Adjacency = unordered_map<Vertex, shared_ptr<Edge>>;
    using VertexContent = pair<T_vdata, Adjacency>;
    using VertexMap = unordered_map<Vertex, VertexContent>;

    ////////////////////////////////////////////////////////////////
    /// \brief adj_map: adjacency map of the vertices.
    VertexMap vmap;

    ////////////////////////////////////////////////////////////////
    /// \brief size = number of edges (graph size).
    size_t num_edges;

    ////////////////////////////////////////////////////////////////
    /// \brief no_loops has value true if loops are allowed,
    ///        false otherwise.
    bool no_loops;

public:

    Graph(bool no_loops = true);

    // MODIFIERS
    void add_vertex(const Vertex &u, const T_vdata &data = T_vdata());
    void add_edge(const Vertex &u, const Vertex &v,
                  const T_edata &data = T_edata());

    void remove_vertex(const Vertex &u);
    void remove_edge(const Vertex &u, const Vertex &v);

    // ACCESS
    T_vdata const& vertex_data(const Vertex &u) const;
    T_edata const& edge_data(const Vertex &u, const Vertex &v) const;

    bool has_vertex(const Vertex &u) const;
    bool has_edge(const Vertex &u, const Vertex &v) const;

    auto vertices() const;
    auto edges() const;

    //    Graph * subgraph(vector<Vertex> &subgraph_vertices) const; // *
    //    Graph * subgraph(vector<pair<long, long>> &subgraph_edges) const; // *

    // PROPERTIES
    auto neighbors(const Vertex &u) const;
    size_t degree(const Vertex &u) const;
    size_t max_degree() const;

    size_t order() const;
    size_t size() const;
};


template<class T_edata, class T_vdata, class T_vertex>
struct Graph<T_edata, T_vdata, T_vertex>::Edge {
    Vertex u;
    Vertex v;
    T_edata data;
    Edge(const T_vertex &u, const T_vertex &v, const T_edata &data) :
        u(u), v(v), data(data) {}
};


////////////////////////////////////////////////////////////////
/// \brief Graph constructor (simple graph as default).
/// \param directed: true if digraph, false otherwise.
/// \param no_loops: true if loops aren't allowed, false
///        otherwise.
///
template<class T_edata, class T_vdata, class T_vertex>
Graph<T_edata, T_vdata, T_vertex>::Graph(bool no_loops)
{
    this->no_loops = no_loops;
    this->num_edges = 0;
}

////////////////////////////////////////////////////////////////
/// \brief add_vertex adds vertex u.
///        the graph.
/// \param u: vertex.
///
template<class T_edata, class T_vdata, class T_vertex>
void Graph<T_edata, T_vdata, T_vertex>::add_vertex(
        const Vertex &u, const T_vdata &data)
{
    if (!has_vertex(u)) {
        VertexContent &vc = vmap[u];
        vc.first = data;
    }
}

////////////////////////////////////////////////////////////////
/// \brief add_edge adds edge e.
/// \param e: edge.
///
template<class T_edata, class T_vdata, class T_vertex>
void Graph<T_edata, T_vdata, T_vertex>::add_edge(
        const Vertex &u, const Vertex &v, const T_edata &data)
{
    if (no_loops && u == v)
        return;

    if (!has_edge(u, v)) {
        add_vertex(u);
        add_vertex(v);

        shared_ptr<Edge> new_edge(new Edge(u, v, data));

        vmap[u].second[v] = new_edge;
        vmap[v].second[u] = new_edge;

        num_edges++;
    }
}

////////////////////////////////////////////////////////////////
/// \brief remove_vertex removes vertex u.
/// \param u_id: ID of vertex u.
///
template<class T_edata, class T_vdata, class T_vertex>
void Graph<T_edata, T_vdata, T_vertex>::remove_vertex(
        const Vertex &u)
{
    if (has_vertex(u)) {
        Adjacency &u_adj = vmap[u].second;
        num_edges = num_edges - u_adj.size();
        for (const auto &[v, uv_ptr] : u_adj) {
            vmap[v].second.erase(u);
        }
        vmap.erase(u);
    }
}

////////////////////////////////////////////////////////////////
/// \brief remove_edge removes the edge indentified by
///        (u_id, v_id).
/// \param u_id, v_id: IDs of the endpoints of the edge.
///
template<class T_edata, class T_vdata, class T_vertex>
void Graph<T_edata, T_vdata, T_vertex>::remove_edge(
        const Graph::Vertex &u, const Graph::Vertex &v)
{
    if (has_edge(u, v)) {
        vmap[u].second.erase(v);
        vmap[v].second.erase(u);
        num_edges--;
    }
}

////////////////////////////////////////////////////////////////
/// \brief vertex gives a constant access to the vertex with
///        ID u_id.
/// \param u_id: vertex ID.
/// \return a constant Vertex reference.
///
template<class T_edata, class T_vdata, class T_vertex>
const T_vdata & Graph<T_edata, T_vdata, T_vertex>::vertex_data(
        const Vertex &u) const
{
    return vmap.at(u).first;
}

////////////////////////////////////////////////////////////////
/// \brief vertex gives a constant access to the edge
/// identified by (u_id, v_id).
/// \param u_id, v_id: IDs of the endpoints of the edge.
/// \return a constant Edge reference.
///
template<class T_edata, class T_vdata, class T_vertex>
const T_edata & Graph<T_edata, T_vdata, T_vertex>::edge_data(
        const Vertex &u, const Vertex &v) const
{
    return vmap.at(u).second.at(v)->data;
}

////////////////////////////////////////////////////////////////
/// \brief has_vertex check if vertex u is in the graph
/// \param u_id: ID of the vertex u.
/// \return true if u in the graph, false otherwise.
///
template<class T_edata, class T_vdata, class T_vertex>
bool Graph<T_edata, T_vdata, T_vertex>::has_vertex(
        const Vertex &u) const
{
    return vmap.find(u) != vmap.end();
}

////////////////////////////////////////////////////////////////
/// \brief has_edge
/// \param u_id
/// \param v_id
/// \return
///
template<class T_edata, class T_vdata, class T_vertex>
bool Graph<T_edata, T_vdata, T_vertex>::has_edge(
        const Vertex &u, const Vertex &v) const
{
    typename VertexMap::const_iterator u_it = vmap.find(u);
    if (u_it != vmap.end()) {
        const Adjacency &u_adj = u_it->second.second;
        if (u_adj.find(v) != u_adj.end())
            return true;
    }
    return false;
}

////////////////////////////////////////////////////////////////
/// \brief neighbors
/// \param u
/// \return
///
template<class T_edata, class T_vdata, class T_vertex>
auto Graph<T_edata, T_vdata, T_vertex>::neighbors(const Vertex &u) const
{
    using AdjIt = typename Adjacency::const_iterator;
    const Adjacency &u_adj = vmap.at(u).second;

    return Generator<Vertex, AdjIt> (
        u_adj.begin(),
        u_adj.end(),
        [](AdjIt const& it)-> Vertex {return it->first;}
    );
}

////////////////////////////////////////////////////////////////
/// \brief degree
/// \param u
/// \return
///
template<class T_edata, class T_vdata, class T_vertex>
size_t Graph<T_edata, T_vdata, T_vertex>::degree(
        const Graph::Vertex &u) const
{
    return vmap.at(u).second.size();
}

template<class T_edata, class T_vdata, class T_vertex>
size_t Graph<T_edata, T_vdata, T_vertex>::max_degree() const
{
    size_t max = 0;
    for (auto const& u : vmap) {
        size_t u_degree = u.second.second.size();
        if (u_degree > max)
            max = u_degree;
    }
    return max;
}

////////////////////////////////////////////////////////////////
/// \brief has_edge
/// \param u_id
/// \param v_id
/// \return
///
template<class T_edata, class T_vdata, class T_vertex>
auto Graph<T_edata, T_vdata, T_vertex>::vertices() const
{
    using VMapIt = typename VertexMap::const_iterator;

    return Generator<Vertex, VMapIt> (
        vmap.begin(),
        vmap.end(),
        [](VMapIt it)->Vertex {return it->first;}
    );
}

////////////////////////////////////////////////////////////////
/// \brief edges
/// \return
///
template<class T_edata, class T_vdata, class T_vertex>
auto Graph<T_edata, T_vdata, T_vertex>::edges() const
{
    using VMapIt = typename VertexMap::const_iterator;
    using AdjIt = typename Adjacency::const_iterator;
    using State = tuple<VMapIt, AdjIt, unordered_set<Vertex>>;
    const Adjacency &begin_adj = vmap.begin()->second.second;

    return Generator<pair<Vertex, Vertex>, State> (
        make_tuple(vmap.begin(), begin_adj.begin(), unordered_set<Vertex>()),
        make_tuple(vmap.end(), begin_adj.end(), unordered_set<Vertex>()),
        [](State const &state)-> pair<Vertex, Vertex> {
            auto& [vmap_it, adj_it, visited] = state;
            return make_pair(vmap_it->first, adj_it->first);
        },
        [this](State &state)-> void {
            auto& [vmap_it, adj_it, seen] = state;

            do {
                if (adj_it != vmap_it->second.second.end()) {
                    adj_it++;
                } else {
                    seen.insert(vmap_it->first);
                    vmap_it++;
                    if (vmap_it != vmap.end())
                        adj_it = vmap_it->second.second.begin();
                }

            } while (vmap_it != vmap.end() &&
                     (adj_it == vmap_it->second.second.end() ||
                      seen.find(adj_it->first) != seen.end() ));
        },
        [](State const& s1, State const& s2)-> bool {
            return get<0>(s1) == get<0>(s2);
        }
    );
}

////////////////////////////////////////////////////////////////
/// \brief order
/// \return
///
template<class T_edata, class T_vdata, class T_vertex>
size_t Graph<T_edata, T_vdata, T_vertex>::order() const
{
    return vmap.size();
}

////////////////////////////////////////////////////////////////
/// \brief size
/// \return
///
template<class T_edata, class T_vdata, class T_vertex>
size_t Graph<T_edata, T_vdata, T_vertex>::size() const
{
    return num_edges;
}

#endif // GRAPH_HPP

