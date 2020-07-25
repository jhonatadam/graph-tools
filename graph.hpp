#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <vector>
#include <list>
#include <unordered_map>

using namespace std;

struct Vertex 
{
    long id;
    Vertex(const long &id) : id(id) {}
};

struct Edge
{
    long u_id;
    long v_id;
    Edge(const long &u_id, const long &v_id) :
        u_id(u_id), v_id(v_id) {}
};

template <class Te = Edge, class Tv = Vertex>
class Graph
{
    typedef unordered_map<long,Edge*> Adjacency;
    typedef unordered_map<long,Adjacency> AdjacencyMap;
    typedef unordered_map<long,Vertex> VertexMap;
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
    void add_vertex(const Vertex &u);

    ////////////////////////////////////////////////////////////////
    /// \brief add_edge adds edge e.
    /// \param e: edge.
    ///
    void add_edge(const Edge &e);

    ////////////////////////////////////////////////////////////////
    /// \brief remove_vertex removes vertex u.
    /// \param u_id: ID of vertex u.
    ///
    void remove_vertex(const long &u_id); // *

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
    const Vertex& vertex(const long &u_id) const;

    ////////////////////////////////////////////////////////////////
    /// \brief vertex gives a constant access to the edge
    /// identified by (u_id, v_id).
    /// \param u_id, v_id: IDs of the endpoints of the edge.
    /// \return a constant Edge reference.
    ///
    const Edge& edge(const long &u_id, const long &v_id) const;

    ////////////////////////////////////////////////////////////////
    /// \brief has_vertex check if vertex u is in the graph
    /// \param u_id: ID of the vertex u.
    /// \return true if u in the graph, false otherwise.
    ///
    bool has_vertex(const long &u_id) const;

    ////////////////////////////////////////////////////////////////
    /// \brief has_edge
    /// \param u_id
    /// \param v_id
    /// \return
    ///
    bool has_edge(const long &u_id, const long &v_id) const;

    ////////////////////////////////////////////////////////////////
    /// \brief vertices
    /// \return
    ///
    vector<Vertex>& vertices() const; // *

    ////////////////////////////////////////////////////////////////
    /// \brief edges
    /// \return
    ///
    vector<pair<Vertex, Vertex>>& edges() const; // *

    ////////////////////////////////////////////////////////////////
    /// \brief subgraph
    /// \param subgraph_vertices
    /// \return
    ///
    Graph * subgraph(vector<long> &subgraph_vertices) const; // *

    ////////////////////////////////////////////////////////////////
    /// \brief subgraph
    /// \param subgraph_edges
    /// \return
    ///
    Graph * subgraph(vector<pair<long, long>> &subgraph_edges) const; // *

    // PROPERTIES
    ////////////////////////////////////////////////////////////////
    /// \brief adjacency
    /// \param u
    /// \return
    ///
    vector<Vertex>& adjacency(Vertex u) const; // *

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


template<class Te, class Tv>
Graph<Te, Tv>::Graph(bool directed, bool no_loops)
{
    this->directed = directed;
    this->no_loops = no_loops;
}

template<class Te, class Tv>
void Graph<Te, Tv>::add_vertex(const Vertex &u)
{
    if (!this->has_vertex(u.id)) {
        this->vertex_map.insert({u.id, u});
    }
}

template<class Te, class Tv>
void Graph<Te, Tv>::add_edge(const Edge &e)
{
    if (this->no_loops && e.u_id == e.v_id)
        return;

    if (!this->has_edge(e.u_id, e.v_id)) {
        this->add_vertex(Vertex(e.u_id));
        this->add_vertex(Vertex(e.v_id));
        this->edge_list.push_back(e);

        this->adj_map[e.u_id][e.v_id] = &edge_list.back();
        if (!this->directed)
            this->adj_map[e.v_id][e.u_id] = &edge_list.back();
    }
}

template<class Te, class Tv>
void Graph<Te, Tv>::remove_vertex(const long &u_id)
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

template<class Te, class Tv>
const Vertex & Graph<Te, Tv>::vertex(const long &u_id) const {
    return this->vertex_map.at(u_id);
}

template<class Te, class Tv>
const Edge & Graph<Te, Tv>::edge(const long &u_id, const long &v_id) const
{
    return *this->adj_map.at(u_id).at(v_id);
}

template<class Te, class Tv>
bool Graph<Te, Tv>::has_vertex(const long &u_id) const
{
    return this->vertex_map.find(u_id) != this->vertex_map.end();
}

template<class Te, class Tv>
bool Graph<Te, Tv>::has_edge(const long &u_id, const long &v_id) const
{
    AdjacencyMap::const_iterator u_it = this->adj_map.find(u_id);
    if (u_it != this->adj_map.end()) {
        const Adjacency &u_adj = u_it->second;
        if (u_adj.find(v_id) != u_adj.end())
            return true;
    }
    return false;
}

template<class Te, class Tv>
size_t Graph<Te, Tv>::order() const
{
    return vertex_map.size();
}

template<class Te, class Tv>
size_t Graph<Te, Tv>::size() const
{
    return edge_list.size();
}

#endif // GRAPH_HPP
