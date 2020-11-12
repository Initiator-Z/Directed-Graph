#include <iostream>
#include <vector>
#include <queue>
#include <stack>
#include <unordered_set>
#include <unordered_map>
#include <set>
#include <array>
#include <list>
#include <forward_list>
#include <deque>
#include <map>
#include <cstddef>
#include <string>
#include <utility>
#include <algorithm>
#include <limits>
#include <optional>
#include <exception>
#include <stdexcept>

#include "directed_graph.h"

using namespace std;
template <typename T>
void test_topo(directed_graph<T>& g){
    vector<vertex<T>> result = topological_sort(g);
    for(vertex<T> vt : result){
        cout << vt.id << endl;
    }
}


/*
 * Computes the shortest distance from u to v in graph g.
 * The shortest path corresponds to a sequence of vertices starting from u and ends at v,
 * which has the smallest total weight of edges among all possible paths from u to v.
 */
template <typename T>
vector<vertex<T>> shortest_path(directed_graph<T>& g, int& u_id, int& v_id) {
    vector<vertex<T>> path;
    T inf = 99999;
    int size = g.max_id();

    // Store vertices with its predecessor in the path
    // predecessor of i is predecessor[i]
    int predecessor[size];
    // contain remaining vertices' id and edge weights
    vector<pair<T, int>> remain;

    // Initialise vertices' costs to reach
    for(int i = 0; i <= size; i++){
        if(i != u_id && g.contains(i)){
            if(g.adjacent(u_id, i)){
                remain.push_back(make_pair(g.get_edge(u_id, i), i));
                predecessor[i] = u_id;
            }
            else{
                remain.push_back(make_pair(inf, i));
            }
        }
    }

    int current = u_id;
    path.push_back(g.retrieve(current));

    // Continuously set each vertices as current to update costs until all vertices were set
    while(!remain.empty()){
        T current_weight;
        // find vertex in the remaining vector with smallest distance and select as current
        sort(remain.begin(), remain.end());

        current = remain[0].second;
        current_weight = remain[0].first;
        remain.erase(remain.begin());

        // Update weights
        for(int i = 0; i < remain.size(); i++){
            if(g.adjacent(current, remain[i].second)){
                if(remain[i].first > current_weight + g.get_edge(current, remain[i].second)){
                    remain[i].first = current_weight + g.get_edge(current, remain[i].second);
                    predecessor[remain[i].second] = current;
                }
            }
        }
    }

    // Starting from destination vertex, trace back to source vertex
    int last = v_id;
    while(last != u_id){
        path.emplace(path.begin()+1, g.retrieve(last));
        last = predecessor[last];
    }

    return path;
}

/*
 * Computes the strongly connected components of the graph.
 * A strongly connected component is a subset of the vertices
 * such that for every pair u, v of vertices in the subset,
 * v is reachable from u and u is reachable from v.
 */

// Stack keep track of visited vertices
stack<int> s;
// Attributes of vertices to trace visiting order and low values
int dfn[999] = {0}, low[999], dfn_cnt = 0;
// Check if a vertex is in stack
set<int> in_stack;
// Resulting SCCs
set<set<int>> all_scc;

template <typename T>
vector<vector<vertex<T>>> strongly_connected_components(directed_graph<T>& g) {
    vector<vector<vertex<T>>> results;
    vector<vertex<T>> vertices = g.get_vertices();

    // Perform Tarjans' algorithm on all vertices in the graph if not included in a SCC already
    for(vertex<T> vt : vertices){
        if(!dfn[vt.id]){
            tarjan(vt.id, g);
        }
    }

    for(set<int> st : all_scc){
        vector<vertex<T>> vector;
        for(int i : st){
            vector.push_back(g.retrieve(i));
        }
        results.push_back(vector);
    }

    return results;
}

template <typename T>
void tarjan(int u, directed_graph<T>& g){
    low[u] = dfn[u] = ++dfn_cnt;
    s.push(u); // Add u to stack
    in_stack.insert(u); // mark u as visited

    vector<vertex<T>> neighbours = g.get_neighbours(u);
    for(vertex<T> vt : neighbours){
        if(!dfn[vt.id]){// v not visited
            tarjan(vt.id, g);
            low[u] = min(low[u], low[vt.id]);
        }
        else if(in_stack.find(vt.id) != in_stack.end()){// v visited and in stack
            low[u] = min(low[u], dfn[vt.id]);
        }
    }

    if(dfn[u] == low[u]){// u is root of a scc
        set<int> scc;
        while(s.top() != u){// add all vertices above u in stack to a scc
            scc.insert(s.top());
            in_stack.erase(s.top());
            s.pop();
        }
        scc.insert(s.top());
        in_stack.erase(s.top());
        s.pop();
        all_scc.insert(scc);
    }
}

/*
 * Computes a topological ordering of the vertices.
 * For every vertex u in the order, and any of its
 * neighbours v, v appears later in the order than u.
 * You will be given a DAG as the argument.
 */

// Stack keep track of topological orders
stack<int> topo_order;
// Set marking visited vertices
set<int> marked;

template <typename T>
vector<vertex<T>> topological_sort(directed_graph<T>& g) {
    vector<vertex<T>> vertices = g.get_vertices();

    // Visit vertices if not marked already
    for(vertex<T> vt : vertices){
        while(marked.find(vt.id) == marked.end()){
            visit(vt.id, g);
        }
    }

    vector<vertex<T>> order;
    while(!topo_order.empty()){
        order.push_back(g.retrieve(topo_order.top()));
        topo_order.pop();
    }
    return order;
}
template <typename T>
void visit(int u, directed_graph<T>& g){
    // If u already marked, return
    if(marked.find(u) != marked.end()){
        return;
    }

    // If u not marked, recursively visit its neighbours
    vector<vertex<T>> neighbours = g.get_neighbours(u);
    for(vertex<T> vt : neighbours){
        visit(vt.id, g);
    }
    // Mark u and push to the stack
    marked.insert(u);
    topo_order.push(u);
}


/*
 * Computes the lowest cost-per-person for delivery over the graph.
 * u is the source vertex, which send deliveries to all other vertices.
 * vertices denote cities; vertex weights denote cities' population;
 * edge weights denote the fixed delivery cost between cities, which is irrelevant to
 * the amount of goods being delivered.
 */
template <typename T>
T low_cost_delivery(directed_graph<T>& g, int& u_id) {
    T population = 0;
    T cost = 0;
    vector<vertex<T>> vertices = g.get_vertices();

    // set containing vertices already delivered goods to.
    set<int> delivered;

    // Calculate total populations of cities to be delivered
    for(vertex<T> vt : vertices){
        if(vt.id != u_id && g.reachable(u_id, vt.id)){
            population += vt.weight;
        }
    }

    // set contains cities been considered in planning next delivery.
    // If this set is identical with 'delivered' set at certain step,
    // dynamic programming for planning deliveries is completed.
    set<int> considered;

    // Store potential cities to be delivered along with costs.
    vector<pair<T, int>> neighbour_weights;

    // Initialise dynamic programming by adding source city to delivered set.
    delivered.insert(u_id);

    // Continue planning for optimum delivery cost while not all
    // cities are considered as departure.
    while(considered.size() != delivered.size()){
        for(int ct_del : delivered){

            // Get neighbours of current departure city.
            // next lowest cost destination must by one of the neighbours.
            considered.insert(ct_del);
            vector<vertex<T>> neighbours = g.get_neighbours(ct_del);

            // If a neighbour is never delivered to, deliver to that city.
            for(vertex<T> vt : neighbours){
                if(delivered.find(vt.id) == delivered.end()){
                    neighbour_weights.push_back(make_pair(g.get_edge(ct_del, vt.id), vt.id));
                }
                // If a neighbour has been delivered to, check if the new route offers lower cost, and also
                // whether the neighbour comes earlier than current city in topological sequence(if yes cannot replace route).
                else{
                    for(int i = 0; i < neighbour_weights.size(); i++){
                        if(neighbour_weights[i].second == vt.id && neighbour_weights[i].first > g.get_edge(ct_del, vt.id) && !g.reachable(neighbour_weights[i].second, ct_del)){
                            cost -= neighbour_weights[i].first;
                            delivered.erase(neighbour_weights[i].second);
                            neighbour_weights.push_back(make_pair(g.get_edge(ct_del, vt.id), vt.id));
                        }
                    }
                }

            }

            // Sort all possible cities to be delivered and deliver to destination with lowest cost.
            sort(neighbour_weights.begin(), neighbour_weights.end());
            for(int i = 0; i < neighbour_weights.size(); i++){
                if(delivered.find(neighbour_weights[i].second) == delivered.end()){
                    delivered.insert(neighbour_weights[i].second);
                    cost += neighbour_weights[i].first;
                    break;
                }
            }
        }
    }
    return cost / population;
}