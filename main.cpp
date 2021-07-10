/*
* A main function for you to build and run your
* own tests with.
* This file is not part of the marking, so you
* can do anything you want here.
*/
#include <iostream>

#include "directed_graph_algorithms.cpp"


int main() {
    directed_graph<int> g;

    g.add_vertex(vertex(1,800));
    g.add_vertex(vertex(2,300));
    g.add_vertex(vertex(3,400));
    g.add_vertex(vertex(4,710));
    g.add_vertex(vertex(5,221));
    g.add_vertex(vertex(6,400));
    g.add_vertex(vertex(7,710));
    g.add_vertex(vertex(8,221));

    g.add_edge(2,3,600);
    g.add_edge(2,4,900);
    g.add_edge(6,4,3000);
    g.add_edge(6,5,4000);
    g.add_edge(6,8,1);
    g.add_edge(1,3,700);
    g.add_edge(4,1,500);
    g.add_edge(5,1,500);
    g.add_edge(1,7,500);
    g.add_edge(8,7,500);

    //Returns the vertices of the graph in the order they are visited in by a depth-first traversal 
    //starting at the given vertex.
    g.depth_first(1); 
    
    //Returns a spanning tree of the graph starting at the given vertex using the out-edges. 
    //This means every vertex in the tree is reachable from the root.
    g.out_tree(2); 

    //Computes the shortest distance from 1 to 8 in the graph g
    shortest_path(g, 1, 8)
    
    //Test for strongly comment componnets in the given graph
    strongly_connected_components(g);

    //Computes a topological ordering of the vertices in the given graph. 
    //A Directed acyclic graph(DAG) is required.
    topological_sort(g);

}
