/*
 * node_dists.cpp
 * Class to calculate distances between graph nodes
 *
 * Some code modified from DBGWAS by Jacob L., Jaillard M., Lima L.
 * Please cite:
 *  Magali Jaillard, Maud Tournoud, Leandro Lima, Vincent Lacroix,
 *  Jean-Baptiste Veyrieras and Laurent Jacob, "Representing Genetic
 *  Determinants in Bacterial GWAS with Compacted De Bruijn Graphs", 2017,
 *  Cold Spring Harbor Labs Journals, doi:10.1101/113563.
 *  (url: http://www.biorxiv.org/content/early/2017/03/03/113563)
 *
 *  This code: John Lees
 *
 */

#include "node_dists.hpp"

// Graph initialisation
Cdbg::Cdbg(const string& dbgPrefix)
   : Cdbg(dbgPrefix + ".nodes", dbgPrefix + ".edges.dbg")
{
}

Cdbg::Cdbg(const string& nodeFile, const string& edgeFile)
{
    int nbContigs = getNbLinesInFile(nodeFile);
    _dbgGraph = graph_t(nbContigs);

    // Create the nodes
    ifstream nodeIst(nodeFile.c_str());
    if (!nodeIst)
    {
        throw std::runtime_error("Could not open node file " + nodeFile + "\n");
    }

    int id;
    string sequence;
    while (nodeIst >> id >> sequence)
    {
        MyVertex vF = vertex(id, _dbgGraph);
        _dbgGraph[vF].name = sequence;
        _dbgGraph[vF].length = sequence.length();
        _dbgGraph[vF].id = id;
        _seqs[sequence] = id;
    }

    nodeIst.close();

    // Create the edges
    ifstream edgeIst(edgeFile.c_str());
    if (!edgeIst)
    {
        throw std::runtime_error("Could not open edge file " + edgeFile + "\n");
    }

    int from, to;
    int index = 0;
    string label;
    while (edgeIst >> from >> to >> label)
    {
        MyVertex fromVertex = vertex(from, _dbgGraph);
        MyVertex toVertex = vertex(to, _dbgGraph);

        // Add from -> to
        pair<MyEdge, bool> return_from_forward_edge = add_edge(fromVertex,
                                                           toVertex,
                                                           _dbgGraph);
        if (return_from_forward_edge.second) {
            _dbgGraph[return_from_forward_edge.first].id = index;
            _dbgGraph[return_from_forward_edge.first].weight = _dbgGraph[fromVertex].length;
            index++;
        }

        // Add to -> from
        pair<MyEdge, bool> return_from_to_edge = add_edge(toVertex,
                                                           fromVertex,
                                                           _dbgGraph);
        if (return_from_to_edge.second) {
            _dbgGraph[return_from_to_edge.first].id = index;
            _dbgGraph[return_from_to_edge.first].weight = _dbgGraph[toVertex].length;
            index++;
        }
    }
    edgeIst.close();
}

// Graph alogorithms and helper functions

vector<int> Cdbg::node_distance(const string& origin)
{
    vector<int> distances(num_vertices(_dbgGraph));
    fill(distances.begin(), distances.end(), 0);

    MyVertex originVertex = vertex(_seqs[origin], _dbgGraph);
    dijkstra_shortest_paths(_dbgGraph, originVertex, weight_map(boost::get(&EdgeInfo::weight, _dbgGraph)).
                            distance_map(boost::make_iterator_property_map(distances.begin(),
                                                                           boost::get(boost::vertex_index, _dbgGraph))));

    return(distances);
}

long int getNbLinesInFile(const string &filename) {
    ifstream ist(filename.c_str());
    if (!ist)
    {
        throw std::runtime_error("Could not open file " + filename + "\n");
    }

    // Number of lines in the file
    long int n = std::count(std::istreambuf_iterator<char>(ist), std::istreambuf_iterator<char>(), '\n');

    ist.close();

    return n;
}
