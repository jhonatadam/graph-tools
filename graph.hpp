#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <vector>
#include <list>
#include <unordered_map>
#include <unordered_set>

using namespace std;

enum GraphType
{
    GRAPH,
    DIGRAPH,
    MULTIGRAPH
};

template <class T1 = double, class T2 = long>
class Graph
{
public:
    typedef T1 Weight;
    typedef T2 Vertex;
    typedef long Label;

private:
    struct Edge
    {
        Vertex u;
        Vertex v;
        Label label;
        Weight weight;
        Edge(Vertex u, Vertex v, Label uv_label, Weight uv_weight) :
        	u(u), v(v), label(uv_label), weight(uv_weight){}
    };

    list<Edge> edge_list;
    unordered_map<Vertex, Weight> vtx_weights;
    unordered_map<Vertex, unordered_map<Vertex, list<Edge *>>> adj_list;

    GraphType graph_type;

public:
    // CONSTRUCTOR
    Graph(GraphType graph_type = GRAPH);

    // MODIFIERS
    void add_vertex(Vertex u);
    void add_vertex(Vertex u, Weight u_weight);
    void remove_vertex(Vertex u); // *
    void set_vertex_weight(Vertex u, Weight new_weight);

    void add_edge(Vertex u, Vertex v, Label uv_label = 0);
    void add_edge(Vertex u, Vertex v, Weight uv_weight, Label uv_label = 0);
    void remove_edge(Vertex u, Vertex v); // *
    void set_edge_weight(Vertex u, Vertex v, Weight new_weight); // *
    void set_edge_weight(Vertex u, Vertex v, Label uv_label, Weight new_weight); // *

    // ACCESS
    vector<Vertex>& vertices(); // *
    size_t order();
    Weight vertex_weight(Vertex u);
    vector<Vertex>& adjacency(Vertex u); // *
    size_t degree(Vertex u); // *
    size_t max_degree(); // *
    
    vector<pair<Vertex, Vertex>>& edges(); // *
    size_t size();
    Vertex edge_label(Vertex u, Vertex v); // *
    Weight edge_weight(Vertex u, Vertex v); // *

};


template<class T1, class T2>
Graph<T1, T2>::Graph(GraphType graph_type)
{
    this->graph_type = graph_type;
}

template<class T1, class T2>
void Graph<T1, T2>::add_vertex(Vertex u)
{
    if (this->adj_list.find(u) == this->adj_list.end())
        adj_list[u];
}

template<class T1, class T2>
void Graph<T1, T2>::add_vertex(Vertex u, Weight u_weight)
{
    if (this->adj_list.find(u) == this->adj_list.end()) {
        adj_list[u];
        vtx_weights[u] = u_weight;
    }
}

template<class T1, class T2>
void Graph<T1, T2>::remove_vertex(Vertex u)
{
    // getting u iterator if it exists
    auto u_it = this->adj_list.find(u);
    // if vertex u is in the graph
    if (u_it != this->adj_list.end()){
        // for each neighbor of u
        for (pair<Vertex, list<Edge*>>& uv_edges: u_it->second) {
            // implement remove_edge before this procedure
            // cout << x.first << ": " << x.second;
        }

        adj_list.erase(u);
        vtx_weights.erase(u);
    }
}

template<class T1, class T2>
void Graph<T1, T2>::set_vertex_weight(Vertex u, Weight new_weight)
{
    if (this->adj_list.find(u) != this->adj_list.end())
        vtx_weights[u] = new_weight;
}

template<class T1, class T2>
void Graph<T1, T2>::add_edge(Vertex u, Vertex v, Label uv_label)
{
    this->add_edge(u, v, Weight(0), uv_label);
}

template<class T1, class T2>
void Graph<T1, T2>::add_edge(Vertex u, Vertex v, Weight uv_weight, 
                             Label uv_label) 
{
    list<Edge*> &edges = adj_list[u][v];

	if (edges.empty() || this->graph_type == MULTIGRAPH) {
		edge_list.push_back(Edge(u, v, uv_label, uv_label));
		edges.push_back(&edge_list.back());

		if (this->graph_type != DIGRAPH)
			adj_list[v][u].push_back(&edge_list.back());
		else
			adj_list[v];
	}
}

template<class T1, class T2>
size_t Graph<T1, T2>::order() 
{
    return adj_list.size();
}

template <class T1, class T2>
T1 Graph<T1, T2>::vertex_weight(Vertex u) {
    auto it = this->vtx_weights.find(u);
    if (it != this->vtx_weights.end())
        return it->second;
    return Weight(0);
}

template<class T1, class T2>
size_t Graph<T1, T2>::size() 
{
    return edge_list.size();
}

#endif // GRAPH_HPP