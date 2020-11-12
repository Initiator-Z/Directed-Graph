#ifndef DIRECTED_GRAPH_H
#define DIRECTED_GRAPH_H

#include<iostream>
#include<string>
#include<vector>
#include<queue>
#include<stack>
#include<utility>
#include<algorithm>

using namespace std;

template <typename T>
class vertex {

public:
    int id;
    T weight;
    T edge_weight;

    vertex(int x, T y) : id(x), weight(y), edge_weight() {}
    vertex() : id(), weight(), edge_weight() {}

};


template <typename T>
class directed_graph {

private:
    vector<vector<vertex<T>>> adj_list;

public:
    directed_graph(); //A constructor for directed_graph. The graph should start empty.
    ~directed_graph(); //A destructor. Depending on how you do things, this may not be necessary.

    bool contains(const int&) const; //Returns true if the graph contains the given vertex_id, false otherwise.
    bool adjacent(const int&, const int&) const; //Returns true if the first vertex is adjacent to the second, false otherwise.
    vertex<T> retrieve(const int&); //Return specific vertex according to the given vertex id.
    int max_id(); //Return the maximum vertex id in the graph.

    void add_vertex(const vertex<T>&); //Adds the passed in vertex to the graph (with no edges).
    void add_edge(const int&, const int&, const T&); //Adds a weighted edge from the first vertex to the second.
    T get_edge(const int&, const int&); // return the edge weight from first int to second int.

    void remove_vertex(const int&); //Removes the given vertex. Should also clear any incident edges.
    void remove_edge(const int&, const int&); //Removes the edge between the two vertices, if it exists.

    size_t in_degree(const int&) const; //Returns number of edges coming in to a vertex.
    size_t out_degree(const int&) const; //Returns the number of edges leaving a vertex.
    size_t degree(const int&) const; //Returns the degree of the vertex (both in edges and out edges).

    size_t num_vertices() const; //Returns the total number of vertices in the graph.
    size_t num_edges() const; //Returns the total number of edges in the graph.

    vector<vertex<T>> get_vertices(); //Returns a vector containing all the vertices.
    vector<vertex<T>> get_neighbours(const int&); //Returns a vector containing all the vertices reachable from the given vertex. The vertex is not considered a neighbour of itself.
    vector<vertex<T>> get_second_order_neighbours(const int&); // Returns a vector containing all the second_order_neighbours (i.e., neighbours of neighbours) of the given vertex.
    // A vector cannot be considered a second_order_neighbor of itself.
    bool reachable(const int&, const int&); //Returns true if the second vertex is reachable from the first (can you follow a path of out-edges to get from the first to the second?). Returns false otherwise.
    bool contain_cycles(); // Return true if the graph contains cycles (there is a path from any vertices directly/indirectly to itself), false otherwise.

    vector<vertex<T>> depth_first(const int&); //Returns the vertices of the graph in the order they are visited in by a depth-first traversal starting at the given vertex.
    vector<vertex<T>> breadth_first(const int&) ; //Returns the vertices of the graph in the order they are visisted in by a breadth-first traversal starting at the given vertex.

    directed_graph<T> out_tree(const int&); //Returns a spanning tree of the graph starting at the given vertex using the out-edges. This means every vertex in the tree is reachable from the root.

    vector<vertex<T>> pre_order_traversal(const int&, directed_graph<T>&); // returns the vertices in the visiting order of a pre-order traversal of the minimum spanning tree starting at the given vertex.
    vector<vertex<T>> in_order_traversal(const int&, directed_graph<T>&); // returns the vertices in the visiting order of an in-order traversal of the minimum spanning tree starting at the given vertex.
    vector<vertex<T>> post_order_traversal(const int&, directed_graph<T>&); // returns the vertices in ther visitig order of a post-order traversal of the minimum spanning tree starting at the given vertex.

    vector<vertex<T>> significance_sorting(); // Return a vector containing a sorted list of the vertices in descending order of their significance.

};

template <typename T>
directed_graph<T>::directed_graph() {}

template <typename T>
directed_graph<T>::~directed_graph() {}

template <typename T>
bool directed_graph<T>::contains(const int& u_id) const {
    for(int i = 0; i < adj_list.size(); i++){
        if(adj_list[i][0].id == u_id){
            return true;
        }
    }
    return false;
}

template <typename T>
bool directed_graph<T>::adjacent(const int& u_id, const int& v_id) const {
    for(int i = 0; i < adj_list.size(); i++){
        for(int j = 1; j < adj_list[i].size(); j++){
            if(adj_list[i][0].id == u_id && adj_list[i][j].id == v_id){
                return true;
            }
        }
    }
    return false;
}

template <typename T>
vertex<T> directed_graph<T>::retrieve(const int& u_id){
    for(int i = 0; i < adj_list.size(); i++){
        if(adj_list[i][0].id == u_id){
            return adj_list[i][0];
        }
    }
    return vertex<T>();
}

template <typename T>
int directed_graph<T>::max_id(){
    vector<int> id;
    for(int i = 0; i < adj_list.size(); i++){
        id.push_back(adj_list[i][0].id);
    }
    sort(id.begin(), id.end());
    int max_id = id.back();
    return max_id;
}

template <typename T>
void directed_graph<T>::add_vertex(const vertex<T>& u) {
    if(!contains(u.id)) {
        adj_list.resize(adj_list.size() + 1);
        adj_list[adj_list.size() - 1].push_back(u);
    }
    else {
        cout << "Vertex with same id already exist. Vertex id should be unique." << endl;
    }
}

template <typename T>
void directed_graph<T>::add_edge(const int& u_id, const int& v_id, const T& weight) {
    if(!adjacent(u_id, v_id)) {
        if (contains(u_id) && contains(v_id) && u_id != v_id) {
            for (int i = 0; i < adj_list.size(); i++) {
                if (adj_list[i][0].id == u_id) {
                    adj_list[i].push_back(retrieve(v_id));
                    adj_list[i].back().edge_weight = weight;
                }
            }
        }
        else if (u_id == v_id) {
            cout << "Self loops are not allowed." << endl;
        }
        else if (contains(u_id)) {
            cout << "vertex: " << v_id << " not exist" << endl;
        }
        else if (contains(v_id)) {
            cout << "vertex: " << u_id << " not exist" << endl;
        }
    }
    else{
        cout << "An edge already exist between: " << u_id << " and " << v_id << endl;
    }
}

template <typename T>
T directed_graph<T>::get_edge(const int& u_id, const int& v_id){
    for(int i = 0; i < adj_list.size(); i++){
        for(int j = 1; j < adj_list[i].size(); j++){
            if(adj_list[i][0].id == u_id && adj_list[i][j].id == v_id){
                return adj_list[i][j].edge_weight;
            }
        }
    }
    return 0;
}

template <typename T>
void directed_graph<T>::remove_vertex(const int& u_id) {
    if(contains(u_id)) {
        for(int i = 0; i < adj_list.size(); i++){
            for(int j = 1; j < adj_list[i].size(); j ++){
                if(adj_list[i][j].id == u_id){
                    remove_edge(adj_list[i][0].id, adj_list[i][j].id);
                }
            }
        }
        for (int i = 0; i < adj_list.size(); i++) {
            if (adj_list[i][0].id == u_id) {
                adj_list.erase(adj_list.begin() + i);
            }
        }
    }
    else{
        cout << "Vertex not exist, cannot remove." << endl;
    }
}

template <typename T>
void directed_graph<T>::remove_edge(const int& u_id, const int& v_id) {
    if(adjacent(u_id, v_id)){
        for(int i = 0; i < adj_list.size(); i++){
            for(int j = 1; j < adj_list[i].size(); j++){
                if(adj_list[i][0].id == u_id && adj_list[i][j].id == v_id){
                    adj_list[i].erase(adj_list[i].begin()+j);
                }
            }
        }
    }
    else{
        cout << "There's no edge between those vertices, cannot remove edge." << endl;
    }
}

template <typename T>
size_t directed_graph<T>::in_degree(const int& u_id) const {
    size_t count = 0;
    if(contains(u_id)) {
        for (int i = 0; i < adj_list.size(); i++) {
            for (int j = 1; j < adj_list[i].size(); j++) {
                if(adj_list[i][j].id == u_id){
                    count++;
                }
            }
        }
    }
    else{
        cout << "Vertex with given id not exist." << endl;
    }
    return count;
}

template <typename T>
size_t directed_graph<T>::out_degree(const int& u_id) const {
    size_t count = 0;
    if(contains(u_id)){
        for (int i = 0; i < adj_list.size(); i++) {
            if (adj_list[i][0].id == u_id) {
                count = adj_list[i].size()-1;
            }
        }
    }
    else{
        cout << "Vertex with given id not exist." << endl;
    }
    return count;
}

template <typename T>
size_t directed_graph<T>::degree(const int& u_id) const {
    size_t count;
    count  = in_degree(u_id) + out_degree(u_id);
    return count;
}

template <typename T>
size_t directed_graph<T>::num_vertices() const {
    return adj_list.size();
}

template <typename T>
size_t directed_graph<T>::num_edges() const {
    size_t count = 0;
    for(int i = 0; i < adj_list.size(); i++){
        count += out_degree(adj_list[i][0].id);
    }
    return count;
}

template <typename T>
vector<vertex<T>> directed_graph<T>::get_vertices() {
    vector<vertex<T>> vertices;
    for(int i = 0; i < adj_list.size(); i++){
        vertices.push_back(adj_list[i][0]);
    }
    return vertices;
}

template <typename T>
vector<vertex<T>> directed_graph<T>::get_neighbours(const int& u_id) {
    vector<vertex<T>> neighbours;

    //push all reachable vertices from given vertex.
    if(contains(u_id)){
        for (int i = 0; i < adj_list.size(); i++) {
            for (int j = 1; j < adj_list[i].size(); j++) {
                if (adj_list[i][0].id == u_id) {
                    neighbours.push_back(adj_list[i][j]);
                }
            }
        }
    }
    else{
        cout << "Given vertex not exist." << endl;
    }
    return neighbours;
}

template <typename T>
vector<vertex<T>> directed_graph<T>::get_second_order_neighbours(const int& u_id) {
    vector<vertex<T>> sec_neighbours;

    int size = max_id();
    bool added[size];
    for(int i = 0; i < size; i++){
        added[i] = false;
    }

    vector<vertex<T>> temp;
    vector<vertex<T>> temp2;

    //For every neighbour of given vertex, push all their neighbours.
    if(contains(u_id)){
        temp = get_neighbours(u_id);
        for(vertex<T> vt : temp){
            temp2 = get_neighbours(vt.id);
            for(vertex<T> vx : temp2){
                if(vx.id != u_id && !added[vx.id]){
                    sec_neighbours.push_back(vx);
                    added[vx.id] = true;
                }
            }
        }
    }
    else{
        cout << "Given vertex not exist." << endl;
    }

    return sec_neighbours;
}

template <typename T>
bool directed_graph<T>::reachable(const int& u_id, const int& v_id){
    if(contains(u_id) && contains(v_id)){

        // Employ breadth first search with reachable vertices from given vertex(u_id) only,
        // if target vertex(v_id) is found in the process, target is reachable(from u_id).
        int size = max_id();
        bool visited[size];
        for (int i = 0; i <= size; i++) {
            visited[i] = false;
        }

        queue<int> unprocessed;
        unprocessed.push(u_id);

        while (!unprocessed.empty()) {
            int n = unprocessed.front();
            unprocessed.pop();
            if (!visited[n]) {
                visited[n] = true;
                for (int i = 0; i < adj_list.size(); i++) {
                    for (int j = 1; j < adj_list[i].size(); j++) {
                        if (adj_list[i][0].id == n) {
                            unprocessed.push(adj_list[i][j].id);
                            if (adj_list[i][j].id == v_id) {
                                return true;
                            }
                        }
                    }
                }
            }
        }
        return false;
    }
    else{
        cout << "One or both vertices not exist." << endl;
        return false;
    }
}

template <typename T>
bool directed_graph<T>::contain_cycles() {
    for(int i = 0; i < adj_list.size(); i++){

        stack<int> vertices;
        vertices.push(adj_list[i][0].id);

        // Since self loops aren't considered in this assignment,
        // any cycles existed must involve at least two vertices that are reachable from each other.
        while(!vertices.empty()){
            int current = vertices.top();
            vertices.pop();
            vector<vertex<T>> neighbours = get_neighbours(current);

            for(vertex<T> vt : neighbours){
                if(reachable(vt.id, current)){
                    return true;
                }
                else{
                    vertices.push(vt.id);
                }
            }
        }
    }
    return false;
}

template <typename T>
vector<vertex<T>> directed_graph<T>::depth_first(const int& u_id) {
    // Iterative search via the use of a stack that keeps track of
    // vertices not yet finished visiting(has neighbours not yet visited).

    if(contains(u_id)){
        vector<vertex<T>> vertices;

        // bool arrays to keep track vertices that are visited or been on stack.
        int size = max_id();
        bool visited[size];
        for(int i = 0; i <= size; i++){
            visited[i] = false;
        }
        bool pushed[size];
        for(int i = 0; i <= size; i++){
            pushed[i] = false;
        }
        // bool variable to check if all vertices in the graph are visited.
        bool all_visited = true;

        stack<int> unprocessed;
        unprocessed.push(u_id);
        pushed[u_id] = true;

        while(!unprocessed.empty()){
            int n = unprocessed.top();
            bool n_pushed = true;
            bool n_visited = true;
            vector<vertex<T>> neighbours = get_neighbours(n);

            // check for current vertex, if all its neighbours are visited or been on stack.
            // n_visited = neighbours_visited, n_pushed = neighbours_pushed.
            for(vertex<T> vt : neighbours){
                if(!visited[vt.id]){
                    n_visited = false;
                    break;
                }
            }
            for(vertex<T> vt : neighbours){
                if(!pushed[vt.id]){
                    n_pushed = false;
                    break;
                }
            }

            // if all current vertex's neighbour has been visited, or been on stack,
            // or it has no neighbour, finished visiting current vertex, pop it off stack
            if(n_visited || neighbours.size() == 0 || n_pushed){
                unprocessed.pop();
            }

            // visiting current vertex, push all its neighbours onto stack for visiting later.
            if(!visited[n]){
                visited[n] = true;
                vertices.push_back(retrieve(n));
                for(vertex<T> vt : neighbours){
                    unprocessed.push(vt.id);
                    pushed[vt.id] = true;
                }
            }

            // Check if all vertices in the graph are visited.
            for(int i = 0; i < adj_list.size(); i++){
                if(visited[adj_list[i][0].id]){
                    all_visited = true;
                }
                else{
                    all_visited = false;
                    break;
                }
            }

            // If stack is empty and not all vertices in the graph are visited, there are vertices not reachable from the start vertex.
            // thus push them onto the stack to complete the traversal.
            if(unprocessed.empty() && !all_visited){
                for(int i = 0; i < adj_list.size(); i++){
                    if(!visited[adj_list[i][0].id]){
                        unprocessed.push(adj_list[i][0].id);
                    }
                }
            }
        }
        return vertices;
    }
    else{
        cout << "Given vertex not exist." << endl;
        return vector<vertex<T>>();
    }
}

template <typename T>
vector<vertex<T>> directed_graph<T>::breadth_first(const int& u_id){
    if(contains(u_id)) {
        vector<vertex<T>> vertices;

        // bool arrays to keep track vertices that are already visited.
        int size = max_id();
        bool visited[size];
        for (int i = 0; i <= size; i++) {
            visited[i] = false;
        }

        // bool value to check whether all vertices are already visited
        bool all_visited = true;
        queue<int> unprocessed;
        unprocessed.push(u_id);

        // Select the front element of queue, visit it and push all its neighbours onto the queue for visiting later.
        while (!unprocessed.empty()) {
            int n = unprocessed.front();
            unprocessed.pop();
            if (!visited[n]) {
                visited[n] = true;
                vertices.push_back(retrieve(n));
                for (int i = 0; i < adj_list.size(); i++) {
                    for (int j = 1; j < adj_list[i].size(); j++) {
                        if(adj_list[i][0].id == n) {
                            unprocessed.push(adj_list[i][j].id);
                        }
                    }
                }
            }

            // Check if all vertices in the graph are visited.
            for(int i = 0; i < adj_list.size(); i++){
                if(visited[adj_list[i][0].id]){
                    all_visited = true;
                }
                else{
                    all_visited = false;
                    break;
                }
            }

            // If queue is empty and not all vertices in the graph are visited, there are vertices not reachable from the start vertex.
            // thus push them onto the queue to complete the traversal.
            if(unprocessed.empty() && !all_visited){
                for(int i = 0; i < adj_list.size(); i++){
                    if(!visited[adj_list[i][0].id]){
                        unprocessed.push(adj_list[i][0].id);
                    }
                }
            }
        }
        return vertices;
    }
    else{
        cout << "Given vertex not exist." << endl;
        return vector<vertex<T>>();
    }
}

template <typename T>
directed_graph<T> directed_graph<T>::out_tree(const int& u_id) {
    // Implemented using Prim's algorithm.
    // Every step add a vertex that is reachable from tree and
    // has the minimum edge weight out of all reachable vertices from tree.
    directed_graph<T> tree;

    // Store reachable vertices along with the weight of edge connecting to that vertex.
    vector<pair<T, int>> weights;

    // Store vertices already added to tree.
    vector<vertex<T>> tree_nodes;

    // Keep track whether a vertex is added to the tree.
    int size = max_id();
    bool spanned[size];
    for(int i = 0; i < size; i++){
        spanned[i] = false;
    }

    // Keep track whether a potential vertex has already added to weights vector.
    bool added[size];
    for(int i = 0; i < size; i++){
        added[i] = false;
    }

    // Check if there are potential vertices in the graph not yet added to tree.
    bool potential = true;

    // 1. initialise tree with starting vertex
    int current = u_id;
    tree.add_vertex(retrieve(current));
    spanned[current] = true;
    tree_nodes.push_back(retrieve(current));

    // 2. find all reachable vertices from current tree
    while(tree.num_vertices() < adj_list.size() && potential){
        int tree_size = tree.num_vertices();
        for(vertex<T> vt : tree_nodes){

            // For each vertex in tree, find all reachable vertices in the graph.
            // add reachable vertices to weights vector if not already added.
            for(int i = 0; i < adj_list.size(); i++){
                for(int j = 1; j < adj_list[i].size(); j++){
                    if(vt.id == adj_list[i][0].id){

                        //if already added, but has less weight, update.
                        if(added[adj_list[i][j].id]){
                            for(int k = 0; k < weights.size(); k++){
                                if(weights[k].second == adj_list[i][j].id && weights[k].first > adj_list[i][j].edge_weight){
                                    weights.erase(weights.begin()+k);
                                    weights.push_back(make_pair(adj_list[i][j].edge_weight, adj_list[i][j].id));
                                }
                            }
                        }
                        else{
                            weights.push_back(make_pair(adj_list[i][j].edge_weight, adj_list[i][j].id));
                            added[adj_list[i][j].id] = true;
                        }
                    }
                }
            }
        }

        // 3. sort the weights vector with edge weights in ascending order.
        //    add vertex with the least edge weight that is not yet added to the tree.
        sort(weights.begin(), weights.end());

        // Find nodes with smallest edge weight not yet in tree
        for (int i = 0; i < weights.size(); i++) {
            if (!spanned[weights[i].second]) {
                for(int j = 0; j < adj_list.size(); j++){
                    for(int k = 1; k < adj_list[j].size(); k++){
                        if(adj_list[j][k].id == weights[i].second && adj_list[j][k].edge_weight == weights[i].first){
                            current = adj_list[j][0].id;
                        }
                    }
                }

                // add vertex to the tree.
                tree.add_vertex(retrieve(weights[i].second));
                tree.add_edge(current, weights[i].second, weights[i].first);
                tree_nodes.push_back(retrieve(weights[i].second));
                spanned[weights[i].second] = true;
                break;
            }
        }

        // Test if tree size has changed, if no, there are no more potential nodes to add to the tree, stop while loop.
        int new_size = tree.num_vertices();
        if(tree_size == new_size){
            potential = false;
        }
    }
    return tree;
}

template <typename T>
vector<vertex<T>> directed_graph<T>::pre_order_traversal(const int& u_id, directed_graph<T>& mst) {
    vector<vertex<T>> orders;
    stack<vertex<T>> nodes;

    int size = max_id();
    bool visited[size];
    for(int i = 0; i < size; i++){
        visited[i] = false;
    }

    // Visit parent node first before child nodes.
    // if a node have already been visited or is a leaf node, finished with it, pop.
    // else, push its child nodes onto stack for visiting later.
    nodes.push(mst.retrieve(u_id));
    while(!nodes.empty()){
        vertex<T> current = nodes.top();
        if(!visited[current.id]){
            orders.push_back(current);
        }
        vector<vertex<T>> child = mst.get_neighbours(current.id);
        if(child.size() > 0 && !visited[current.id]){
            visited[current.id] = true;
            for(int i = 0; i < child.size(); i++){
                nodes.push(child[i]);
            }
        }
        else{
            nodes.pop();
        }
    }
    return orders;
}

template <typename T>
vector<vertex<T>> directed_graph<T>::in_order_traversal(const int& u_id, directed_graph<T>& mst) {
    vector<vertex<T>> ordered;
    stack<vertex<T>> nodes;

    int size = max_id();
    bool visited[size];
    for(int i = 0; i < size; i++){
        visited[i] = false;
    }

    // Keep track of vertices that were encountered before and has at least one child
    // When a vertex marked encountered is encountered, it has one of its child visited already,
    // thus we can visit that vertex now.
    bool encountered[size];
    for(int i = 0; i < size; i++){
        encountered[i] = false;
    }

    // For a vertex, if it has at least one child node, visit a child node first,
    // then visit the that vertex(parent) before visiting its remaining child nodes.
    nodes.push(mst.retrieve(u_id));
    while(!nodes.empty()){
        vertex<T> current = nodes.top();
        vector<vertex<T>> child = mst.get_neighbours(current.id);

        // If its a leaf node, visit and pop.
        if(child.size() == 0){
            ordered.push_back(current);
            visited[current.id] = true;
            nodes.pop();
        }

            // If it's first encountered, mark it as encountered,
            // push one of its child onto stack for visit.
        else if(!encountered[current.id]){
            nodes.push(child[0]);
            encountered[current.id] = true;
        }

            // At this point, the vertex has one of its child visited,
            // we can move on visiting it and push the remaining of its child nodes
            // onto stack for visit.
        else{
            ordered.push_back(current);
            visited[current.id] = true;
            nodes.pop();
            for(int i = 0; i < child.size(); i++){
                if(!visited[child[i].id]){
                    nodes.push(child[i]);
                }
            }
        }
    }
    return ordered;
}

template <typename T>
vector<vertex<T>> directed_graph<T>::post_order_traversal(const int& u_id, directed_graph<T>& mst) {
    vector<vertex<T>> ordered;
    stack<vertex<T>> nodes;

    int size = max_id();
    bool visited[size];
    for(int i = 0; i < size; i++){
        visited[i] = false;
    }

    // Keep track whether a vertex have been encountered before
    // When a vertex marked encountered is encountered, it's been on hold and
    // its child nodes have been visited, we can finally visit that vertex.
    bool encountered[size];
    for(int i = 0; i < size; i++){
        encountered[i] = false;
    }

    // If encountered a vertex that has at least one child, holding on that vertex and
    // push its child nodes onto stack to visit first.
    nodes.push(mst.retrieve(u_id));
    while(!nodes.empty()){
        vertex<T> current = nodes.top();
        vector<vertex<T>> child = mst.get_neighbours(current.id);

        // Leaf node, simply visit and pop.
        if(child.size() == 0){
            ordered.push_back(current);
            visited[current.id] = true;
            nodes.pop();
        }

            // Not encountered before, hold on with the vertex and push its child nodes.
        else if(!encountered[current.id]){
            for(int i = 0; i < child.size(); i++){
                if(!visited[child[i].id]){
                    nodes.push(child[i]);
                }
            }
            encountered[current.id] = true;
        }

            // At this point, the vertex's child nodes have all been visited,
            // we can finally visit the vertex.
        else{
            visited[current.id] = true;
            ordered.push_back(current);
            nodes.pop();
        }
    }
    return ordered;
}

template <typename T>
vector<vertex<T>> directed_graph<T>::significance_sorting() {
    vector<vertex<T>> ordered;
    vector<pair<T, int>> vertices;

    // Sort vertices based on their weight.
    // Vertex with the greatest weight are considered most significant and pushed first.
    for(int i = 0; i < adj_list.size(); i++){
        vertices.push_back(make_pair(adj_list[i][0].weight, adj_list[i][0].id));
    }

    sort(vertices.begin(), vertices.end());
    reverse(vertices.begin(), vertices.end());

    for(int i = 0; i < vertices.size(); i++){
        ordered.push_back(retrieve(vertices[i].second));
    }
    return ordered;
}

#endif