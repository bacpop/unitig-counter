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
 *  This code: John Lees 2019
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

// Graph algorithms

vector<int> Cdbg::node_distance(const int origin_id)
{
    vector<int> distances(num_vertices(_dbgGraph));
    fill(distances.begin(), distances.end(), 0);

    MyVertex originVertex = vertex(origin_id, _dbgGraph);
    dijkstra_shortest_paths(_dbgGraph, originVertex, weight_map(boost::get(&EdgeInfo::weight, _dbgGraph)).
                            distance_map(boost::make_iterator_property_map(distances.begin(),
                                                                           boost::get(boost::vertex_index, _dbgGraph))));

    return(distances);
}

vector<string> Cdbg::extend_hits(const int origin_id, const int length, const bool repeats)
{
    vector<string> pathSeqs;
    vector<vector<int>> uniquePaths;
    vector<vector<int>> paths = walk_enumeration(_dbgGraph, origin_id, length, repeats);

    // Paths from longest to shortest
    for (auto pathIt = paths.rbegin(); pathIt != paths.rend(); ++pathIt)
    {
        string pathSeq = "";
        vector<int> pathVisits;
        for (auto nodeIt = pathIt->begin(); nodeIt != pathIt->end(); ++nodeIt)
        {
            MyVertex vF = vertex(*nodeIt, _dbgGraph);
            pathVisits.push_back(_dbgGraph[vF].id);
            pathSeq = pathSeq + _dbgGraph[vF].name;
        }

        // Check if path already covered by another, longer path
        int covered = 0;
        for (auto uniqueIt = uniquePaths.begin(); uniqueIt != uniquePaths.end(); ++uniqueIt)
        {
            if (includes(uniqueIt->begin(), uniqueIt->end(), pathVisits.begin(), pathVisits.end()))
            {
                covered = 1;
                break;
            }
        }

        if (!covered)
        {
            pathSeqs.push_back(pathSeq);
            uniquePaths.push_back(pathVisits);
        }
    }
    return pathSeqs;
}

// Helper functions
vector<vector<int>> walk_enumeration(const graph_t& graph, const int start_node, const int length, const bool repeats)
{
    vector<vector<int>> path_list = {{start_node}};
    auto neighbours = boost::adjacent_vertices(start_node, graph);
    for (auto neighbour : boost::make_iterator_range(neighbours))
    {
        MyVertex vF = vertex(neighbour, graph);
        if (length - graph[vF].length >= 0)
        {
            auto paths = walk_enumeration(graph, graph[vF].id, length - graph[vF].length, repeats);
            for (auto path = paths.begin() ; path != paths.end() ; ++path)
            {
                auto search_it = find(path->begin(), path->end(), start_node);
                if (repeats || search_it == path->end())
                {
                    path->insert(path->begin(), start_node);
                    path_list.push_back(*path);
                }
                else
                {
                    break;
                }
            }
        }
    }
    return path_list;
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
