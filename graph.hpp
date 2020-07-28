#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <list>
#include <memory>
#include <unordered_set>
#include <unordered_map>
#include <functional>
#include <variant>

using namespace std;

struct None {};

template <class T_edata = void*,
          class T_vdata = void*,
          class T_vertex = long>
class Graph
{
public:
    struct Edge;
    class EdgeGenerator;
    class VertexGenerator;

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
    /// \brief directed has value true if digraph, false otherwise.
    bool directed;

    ////////////////////////////////////////////////////////////////
    /// \brief no_loops has value true if loops are allowed,
    ///        false otherwise.
    bool no_loops;

public:

    Graph(bool directed = false, bool no_loops = true);

    // MODIFIERS
    void add_vertex(const Vertex &u, const T_vdata &data = T_vdata());
    void add_edge(const Vertex &u, const Vertex &v,
                  const T_edata &data = T_edata());

    ////////////////////////////////////////////////////////////////
    /// \brief remove_vertex removes vertex u.
    /// \param u_id: ID of vertex u.
    ///
    void remove_vertex(const Vertex &u); // *

    ////////////////////////////////////////////////////////////////
    /// \brief remove_edge removes the edge indentified by
    ///        (u_id, v_id).
    /// \param u_id, v_id: IDs of the endpoints of the edge.
    ///
    void remove_edge(const Vertex &u, const Vertex &v); // *

    // ACCESS
    const T_vdata &vertex_data(const Vertex &u) const;
    const T_edata& edge_data(const Vertex &u, const Vertex &v) const;

    bool has_vertex(const Vertex &u) const;
    bool has_edge(const Vertex &u, const Vertex &v) const;

    VertexGenerator vertices() const; // *
    void edges() const; // *

    //    Graph * subgraph(vector<Vertex> &subgraph_vertices) const; // *
    //    Graph * subgraph(vector<pair<long, long>> &subgraph_edges) const; // *

    // PROPERTIES
    void adjacency(Vertex u) const; // *
    size_t degree(Vertex u) const; // *
    size_t max_degree() const; // *

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
Graph<T_edata, T_vdata, T_vertex>::Graph(
        bool directed, bool no_loops)
{
    this->directed = directed;
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
        if (!directed)
            vmap[v].second[u] = new_edge;

        num_edges++;
    }
}

template<class T_edata, class T_vdata, class T_vertex>
void Graph<T_edata, T_vdata, T_vertex>::remove_vertex(
        const Vertex &u)
{
    if (this->has_vertex(u)) {
        if (this->directed) {

        } else {
            for (pair<Vertex, Edge> &neigh : this->adj_map[u]) {

            }

        }

        this->adj_map.erase(u);
        this->vertex_map.erase(u);

        // // getting u iterator if it exists
        // auto u_it = this->adj_list.find(u);
        // // if vertex u is in the graph
        // if (u_it != this->adj_list.end()){
        //     // for each neighbor of u
        //     for (pair<Vertex, list<Edge*>>& uv_edges: u_it->second) {
        //         // implement remove_edge before this procedure
        //         // cout << x.first << ": " << x.second;
        //     }

        //     adj_list.erase(u);
        //     vtx_weights.erase(u);
        // }
    }
}

template<class T_edata, class T_vdata, class T_vertex>
void Graph<T_edata, T_vdata, T_vertex>::remove_edge(
        const Graph::Vertex &u, const Graph::Vertex &v)
{
    if (this->has_edge(u, v)) {
        this->adj_map[u].erase(v);
        if (!this->directed)
            this->adj_map[v].erase(u);
        this->num_edges--;
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

//template<class T_edata, class T_vdata, class T_vertex>
//typename Graph<T_edata, T_vdata, T_vertex>::VertexGenerator
//    Graph<T_edata, T_vdata, T_vertex>::vertices() const
//{
//    return VertexGenerator(adj_map);
//}

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
    return this->num_edges;
}

//template<class T_edata, class T_vdata, class T_vertex>
//class Graph<T_edata, T_vdata, T_vertex>::VertexGenerator
//{
//public:
//    enum GenType {
//        VERTEX_LIST,
//        IN_ADJACENCY,
//        OUT_ADJACENCY
//    };

//    struct iterator {
//        variant<typename VertexMap::const_iterator,
//                typename Adjacency::const_iterator> data_it;

//        GenType gentype;
//        Vertex u;

//        iterator& operator++()
//        {
//            if (gentype == VERTEX_LIST)
//                get<0>(data_it)++;
//            else if (gentype == OUT_ADJACENCY) {
//                get<1>(data_it)++;
//            } else if (IN_ADJACENCY) {

//            }

//            return *this;
//        }

//        iterator operator++(int) = delete;

//        bool operator==(iterator const& other) const
//        {
//            if (gentype == VERTEX_LIST)
//            return get<0>(data_it) == get<0>(other.data_it);
//        }

//        bool operator!=(iterator const& other) const
//        {
//            return !(*this == other);
//        }

//        const Vertex & operator*() const
//        {
//            if (gentype == VERTEX_LIST)
//                return get<0>(data_it)->first;

//            return 0;
//        }
//    };

//    VertexGenerator (const VertexMap &vmap, GenType gentype = VERTEX_LIST,
//                     const Vertex &u = Vertex())
//    {
//        first.u = u;
//        first.gentype = gentype;
//        if (gentype == VERTEX_LIST || gentype == IN_ADJACENCY) {
//            first.data_it = vmap.begin();
//            last.data_it = vmap.end();
//        } else if (gentype == OUT_ADJACENCY) {
//            Adjacency &u_adj = vmap[u].second;
//            first.data_it = u_adj.begin();
//            last.data_it = u_adj.end();
//        } else {
//            // exception!
//        }
//    }

//    iterator begin()
//    {
//        return first;
//    }

//    iterator end()
//    {
//        return last;
//    }

//private:
//    iterator first;
//    iterator last;
//};

#endif // GRAPH_HPP

