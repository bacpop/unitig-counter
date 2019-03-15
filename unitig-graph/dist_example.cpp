/*
 * dist_example.cpp
 * Demonstrates used of Cdbg class
 *
 *  This code: John Lees
 *
 */

#include "node_dists.hpp"

// Simple example program
int main (int argc, char *argv[])
{
    cout << "Reading graph" << endl;
    Cdbg exampleGraph = Cdbg("test_data/graph");

    cout << "Calculating distances from node" << endl;
    // node id 22
    vector<int> dists = exampleGraph.node_distance("GTAATAAACAAAAAAAAAAGTTAAAAATTAAATAAAAATACTTGACTAAATAAAATATACCTGTTAGAATAAAAACAAGGAAAAAGAAAG");
    // dist to 23
    cout << dists[exampleGraph.get_vertex("AAAAAAAAAAGTTAAAAATTACATAAAAATA")] << endl;

    return 0;
}