#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <list>
#include <memory>
#include <unordered_map>

using namespace std;

struct None {};

template <class T_edata = void*,
          class T_vdata = void*,
          class T_vertex = long>
class Graph
{
    typedef T_vertex Vertex;
    struct Edge {
        Vertex u;
        Vertex v;
        T_edata data;
        Edge(const Vertex &u, const Vertex &v, const T_edata &data) :
            u(u), v(v), data(data) {}
    };
    typedef unordered_map<Vertex,shared_ptr<Edge>> Adjacency;
    typedef unordered_map<Vertex,Adjacency> AdjacencyMap;
    typedef unordered_map<Vertex,T_vdata> VertexMap;
    
    ////////////////////////////////////////////////////////////////
    /// \brief adj_map: adjacency map of the vertices.
    ///
    AdjacencyMap adj_map;

    ////////////////////////////////////////////////////////////////
    /// \brief vertex_map
    ///
    VertexMap vertex_map;

    ////////////////////////////////////////////////////////////////
    /// \brief size = number of edges (graph size).
    ///
    size_t num_edges;

    ////////////////////////////////////////////////////////////////
    /// \brief directed has value true if digraph, false otherwise.
    ///
    bool directed;

    ////////////////////////////////////////////////////////////////
    /// \brief no_loops has value true if loops are allowed,
    ///        false otherwise.
    ///
    bool no_loops;

public:
    ////////////////////////////////////////////////////////////////
    /// \brief Graph constructor (simple graph as default).
    /// \param directed: true if digraph, false otherwise.
    /// \param no_loops: true if loops aren't allowed, false
    ///        otherwise.
    ///
    Graph(bool directed = false, bool no_loops = true);

    // MODIFIERS
    ////////////////////////////////////////////////////////////////
    /// \brief add_vertex adds vertex u.
    ///        the graph.
    /// \param u: vertex.
    ///
    void add_vertex(const Vertex &u, const T_vdata &data = T_vdata());

    ////////////////////////////////////////////////////////////////
    /// \brief add_edge adds edge e.
    /// \param e: edge.
    ///
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
    ////////////////////////////////////////////////////////////////
    /// \brief vertex gives a constant access to the vertex with
    ///        ID u_id.
    /// \param u_id: vertex ID.
    /// \return a constant Vertex reference.
    ///
    const T_vdata &vertex_data(const Vertex &u) const;

    ////////////////////////////////////////////////////////////////
    /// \brief vertex gives a constant access to the edge
    /// identified by (u_id, v_id).
    /// \param u_id, v_id: IDs of the endpoints of the edge.
    /// \return a constant Edge reference.
    ///
    const T_edata& edge_data(const Vertex &u, const Vertex &v) const;

    ////////////////////////////////////////////////////////////////
    /// \brief has_vertex check if vertex u is in the graph
    /// \param u_id: ID of the vertex u.
    /// \return true if u in the graph, false otherwise.
    ///
    bool has_vertex(const Vertex &u) const;

    ////////////////////////////////////////////////////////////////
    /// \brief has_edge
    /// \param u_id
    /// \param v_id
    /// \return
    ///
    bool has_edge(const Vertex &u, const Vertex &v) const;

    ////////////////////////////////////////////////////////////////
    /// \brief vertices
    /// \return
    ///
    void vertices() const; // *

    ////////////////////////////////////////////////////////////////
    /// \brief edges
    /// \return
    ///
    void edges() const; // *

    ////////////////////////////////////////////////////////////////
    /// \brief subgraph
    /// \param subgraph_vertices
    /// \return
    ///
//    Graph * subgraph(vector<Vertex> &subgraph_vertices) const; // *

    ////////////////////////////////////////////////////////////////
    /// \brief subgraph
    /// \param subgraph_edges
    /// \return
    ///
//    Graph * subgraph(vector<pair<long, long>> &subgraph_edges) const; // *

    // PROPERTIES
    ////////////////////////////////////////////////////////////////
    /// \brief adjacency
    /// \param u
    /// \return
    ///
    void adjacency(Vertex u) const; // *

    ////////////////////////////////////////////////////////////////
    /// \brief degree
    /// \param u
    /// \return
    ///
    size_t degree(Vertex u) const; // *

    ////////////////////////////////////////////////////////////////
    /// \brief max_degree
    /// \return
    ///
    size_t max_degree() const; // *

    ////////////////////////////////////////////////////////////////
    /// \brief order
    /// \return
    ///
    size_t order() const;

    ////////////////////////////////////////////////////////////////
    /// \brief size
    /// \return
    ///
    size_t size() const;

};

template<class T_edata, class T_vdata, class T_vertex>
Graph<T_edata, T_vdata, T_vertex>::Graph(
        bool directed, bool no_loops)
{
    this->directed = directed;
    this->no_loops = no_loops;
    this->num_edges = 0;
}

template<class T_edata, class T_vdata, class T_vertex>
void Graph<T_edata, T_vdata, T_vertex>::add_vertex(
        const Vertex &u, const T_vdata &data)
{
    if (!this->has_vertex(u)) {
        this->vertex_map[u] = data;
    }
}

template<class T_edata, class T_vdata, class T_vertex>
void Graph<T_edata, T_vdata, T_vertex>::add_edge(
        const Vertex &u, const Vertex &v, const T_edata &data)
{
    if (this->no_loops && u == v)
        return;

    if (!this->has_edge(u, v)) {
        this->add_vertex(u);
        this->add_vertex(v);

        shared_ptr<Edge> new_edge(new Edge(u, v, data));
        this->adj_map[u][v] = new_edge;
        if (!this->directed)
            this->adj_map[v][u] = new_edge;

        this->num_edges++;
    }
}

template<class T_edata, class T_vdata, class T_vertex>
void Graph<T_edata, T_vdata, T_vertex>::remove_vertex(
        const Vertex &u)
{
    if (this->has_vertex(u)) {
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

template<class T_edata, class T_vdata, class T_vertex>
const T_vdata & Graph<T_edata, T_vdata, T_vertex>::vertex_data(
        const Vertex &u) const
{
    return this->vertex_map.at(u);
}

template<class T_edata, class T_vdata, class T_vertex>
const T_edata & Graph<T_edata, T_vdata, T_vertex>::edge_data(
        const Vertex &u, const Vertex &v) const
{
    return this->adj_map.at(u).at(v)->data;
}

template<class T_edata, class T_vdata, class T_vertex>
bool Graph<T_edata, T_vdata, T_vertex>::has_vertex(
        const Vertex &u) const
{
    return this->vertex_map.find(u) != this->vertex_map.end();
}

template<class T_edata, class T_vdata, class T_vertex>
bool Graph<T_edata, T_vdata, T_vertex>::has_edge(
        const Vertex &u, const Vertex &v) const
{
    auto u_it = this->adj_map.find(u);
    if (u_it != this->adj_map.end()) {
        const Adjacency &u_adj = u_it->second;
        if (u_adj.find(v) != u_adj.end())
            return true;
    }
    return false;
}

template<class T_edata, class T_vdata, class T_vertex>
size_t Graph<T_edata, T_vdata, T_vertex>::order() const
{
    return vertex_map.size();
}

template<class T_edata, class T_vdata, class T_vertex>
size_t Graph<T_edata, T_vdata, T_vertex>::size() const
{
    return this->num_edges;
}

#endif // GRAPH_HPP
