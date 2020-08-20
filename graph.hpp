#ifndef GRAPH_HPP
#define GRAPH_HPP


#include <memory>
#include <coroutine>
#include <unordered_map>
#include <unordered_set>

#include "generator.hpp"


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
    /// \brief vmap is a map where the vertices are the keys and,
    /// for each vertex (key), its content is the data and
    /// adjacency of the vertex.
    VertexMap vmap;

    ////////////////////////////////////////////////////////////////
    /// \brief size = number of edges (graph size).
    size_t num_edges;

    ////////////////////////////////////////////////////////////////
    /// \brief no_loops has value true if loops are allowed,
    /// false otherwise.
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

    generator<Vertex const&> vertices() const;
    generator<pair<Vertex, Vertex> const&> edges() const;

    Graph * subgraph(vector<Vertex> const& subgraph_vertices) const;
    Graph * subgraph(vector<pair<long, long>> const& subgraph_edges) const; // *

    // PROPERTIES
    generator<Vertex const&> neighbors(const Vertex &u) const;
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
generator<T_vertex const&>
Graph<T_edata, T_vdata, T_vertex>::neighbors(const Vertex &u) const
{
    const Adjacency &u_neigbors = vmap.at(u).second;
    for (auto const& [neighbor, _] : u_neigbors) {
        co_yield neighbor;
    }
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
generator<T_vertex const&>
Graph<T_edata, T_vdata, T_vertex>::vertices() const
{
    for (auto const& [vertex, _] : vmap)
        co_yield vertex;
}

////////////////////////////////////////////////////////////////
/// \brief edges
/// \return
///
template<class T_edata, class T_vdata, class T_vertex>
generator<pair<T_vertex, T_vertex> const&> Graph<T_edata, T_vdata, T_vertex>::edges() const
{
    unordered_set<Vertex> seen;
    for (auto const& [u, u_info] : vmap) {
        Adjacency const& neighbors = u_info.second;
        for (auto const& [v, _] : neighbors) {
            if (seen.find(v) == seen.end())
                co_yield {u, v};
        }
        seen.insert(u);
    }
}

////////////////////////////////////////////////////////////////
/// \brief subgraph: builds and returns a subgraph induced by
/// the intersection of this graph vertices and a given vector
/// of vertices.
/// \param subgraph_vertices: vector of vertices.
/// \return a pointer to a subgraph.
///
template<class T_edata, class T_vdata, class T_vertex>
Graph<T_edata, T_vdata, T_vertex> *Graph<T_edata, T_vdata, T_vertex>::subgraph(
    const vector<Graph::Vertex> &subgraph_vertices) const
{
    Graph<T_edata, T_vdata, T_vertex> *new_subgraph =
        new Graph<T_edata, T_vdata, T_vertex>(no_loops);

    for (Vertex const& u : subgraph_vertices) {
        typename VertexMap::const_iterator u_it = vmap.find(u);
        if (u_it != vmap.end()) {
            new_subgraph->add_vertex(u, u_it->second.first);

            const Adjacency & u_adj = u_it->second.second;
            for (auto const& [v, edge_ptr] : u_adj)
                if (new_subgraph->has_vertex(v))
                    new_subgraph->add_edge(u, v, edge_ptr->data);
        }
    }

    return new_subgraph;
}

////////////////////////////////////////////////////////////////
/// \brief subgraph: returns a subgraph induced by
/// the intersection of this graph edges and a given vector
/// of edges.
/// \param subgraph_vertices: vector of edges.
/// \return a pointer to the builded subgraph.
///
template<class T_edata, class T_vdata, class T_vertex>
Graph<T_edata, T_vdata, T_vertex> *Graph<T_edata, T_vdata, T_vertex>::subgraph(
        const vector<pair<long, long> > &subgraph_edges) const
{
    Graph<T_edata, T_vdata, T_vertex> *new_subgraph =
        new Graph<T_edata, T_vdata, T_vertex>(no_loops);

    for (auto const& [u, v] : subgraph_edges) {
        if (has_edge(u, v)) {
            new_subgraph->add_vertex(u, vertex_data(u));
            new_subgraph->add_vertex(v, vertex_data(v));

            new_subgraph->add_edge(u, v, edge_data(u, v));
        }
    }

    return new_subgraph;
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

