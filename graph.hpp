#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <list>
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
    };
    typedef unordered_map<Vertex,Edge*> Adjacency;
    typedef unordered_map<Vertex,Adjacency> AdjacencyMap;
    typedef unordered_map<Vertex,T_vdata> VertexMap;
    typedef list<Edge> EdgeList;
    
    ////////////////////////////////////////////////////////////////
    /// \brief adj_map: adjacency map of the vertices.
    ///
    AdjacencyMap adj_map;

    ////////////////////////////////////////////////////////////////
    /// \brief vertex_map
    ///
    VertexMap vertex_map;
    
    ////////////////////////////////////////////////////////////////
    /// \brief edge_list:list of all edges in the graph.
    ///
    EdgeList edge_list;

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
                  const T_edata data = T_edata());

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
    void remove_edge(const long &u_id, const long &v_id); // *

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
        const Vertex &u, const Vertex &v, T_edata data)
{
    if (this->no_loops && u == v)
        return;

    if (!this->has_edge(u, v)) {
        this->add_vertex(u);
        this->add_vertex(v);
        this->edge_list.push_back({u, v, data});

        this->adj_map[u][v] = &edge_list.back();
        if (!this->directed)
            this->adj_map[v][u] = &edge_list.back();
    }
}

template<class T_edata, class T_vdata, class T_vertex>
void Graph<T_edata, T_vdata, T_vertex>::remove_vertex(
        const Vertex &u)
{
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
    return edge_list.size();
}

#endif // GRAPH_HPP
