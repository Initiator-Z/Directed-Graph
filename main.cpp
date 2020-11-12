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


    /**
    g.add_vertex(vertex(1,800));
    g.add_vertex(vertex(2,3000));
    g.add_vertex(vertex(3,400));
    g.add_vertex(vertex(4,710));
    g.add_vertex(vertex(5,221));
    g.add_vertex(vertex(6,800));
    g.add_vertex(vertex(7,3000));
    g.add_vertex(vertex(8,400));
    g.add_vertex(vertex(9,710));

    g.add_edge(1,2,8);
    g.add_edge(1,3,3);
    g.add_edge(1,6,13);
    g.add_edge(2,3,2);
    g.add_edge(2,4,1);
    g.add_edge(3,2,3);
    g.add_edge(3,4,9);
    g.add_edge(3,5,2);


    g.add_edge(4,5,4);
    g.add_edge(4,8,2);
    g.add_edge(4,7,6);
    g.add_edge(5,4,6);
    g.add_edge(5,1,5);
    g.add_edge(5,6,5);
    g.add_edge(5,9,4);

    g.add_edge(6,9,7);
    g.add_edge(6,7,1);
    g.add_edge(7,5,3);
    g.add_edge(7,8,4);
    g.add_edge(8,9,1);
    g.add_edge(9,7,5);
    **/


    int a = 1;
    int b = 5;
    test_topo(g);


    //strongly_connected_components(g);

    //topological_sort(g);
    //cout << low_cost_delivery(g, a);

    //cout << g.get_edge(3,4) << endl;

}