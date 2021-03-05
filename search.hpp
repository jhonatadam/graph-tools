#ifndef SEARCH_HPP
#define SEARCH_HPP

#include <coroutine>
#include <queue>
#include <stack>

#include "graph.hpp"


template <class T_edata = void*,
          class T_vdata = void*,
          class T_vertex = long>
generator<T_vertex const&> depth_first_search(
        Graph<T_edata, T_vdata, T_vertex> const& graph,
        T_vertex const& source)
{
    stack<T_vertex> active;
    unordered_set<T_vertex> visited;

    active.push(source);
    visited.insert(source);
    while(!active.empty()) {
        T_vertex current = active.top();
        active.pop();
        co_yield current;

        for (T_vertex const& neighbor : graph.neighbors(current)) {
            if (visited.find(neighbor) == visited.end()) {
                active.push(neighbor);
                visited.insert(neighbor);
            }
        }
    }
}

template <class T_edata = void*,
          class T_vdata = void*,
          class T_vertex = long>
generator<T_vertex const&> breadth_first_search(
        Graph<T_edata, T_vdata, T_vertex> const& graph,
        T_vertex const& source)
{
    queue<T_vertex> active;
    unordered_set<T_vertex> visited;

    active.push(source);
    visited.insert(source);
    while(!active.empty()) {
        T_vertex current = active.front();
        active.pop();
        co_yield current;

        for (T_vertex const& neighbor : graph.neighbors(current)) {
            if (visited.find(neighbor) == visited.end()) {
                active.push(neighbor);
                visited.insert(neighbor);
            }
        }
    }
}

#endif // SEARCH_HPP
