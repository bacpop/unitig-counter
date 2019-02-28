/*
## Copyright (C) <2017>  <bioMerieux, Universite Claude Bernard Lyon 1,
## Centre National de la Recherche Scientifique>

## 1. This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU Affero General Public License as published
## by the Free Software Foundation version 3 of the  License and under the
## terms of article 2 below.
## 2. This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
## or FITNESS FOR A PARTICULAR PURPOSE. See below the GNU Affero General
## Public License for more details.
## You should have received a copy of the GNU Affero General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.
## 3. Communication to the public by any means, in particular in the form of
## a scientific paper, a poster, a slideshow, an internet page, or a patent,
## of a result obtained directly or indirectly by running this program must
## cite the following paper :
##  Magali Jaillard, Maud Tournoud, Leandro Lima, Vincent Lacroix,
##  Jean-Baptiste Veyrieras and Laurent Jacob, "Representing Genetic
##  Determinants in Bacterial GWAS with Compacted De Bruijn Graphs", 2017,
##  Cold Spring Harbor Labs Journals, doi:10.1101/113563.
##  (url: http://www.biorxiv.org/content/early/2017/03/03/113563)
## -------------------------------------------------------------------------

## Authors (alphabetically): Jacob L., Jaillard M., Lima L.
## Modified by John Lees
*/

#ifndef KSGATB_GENERATE_OUPUT_H
#define KSGATB_GENERATE_OUPUT_H

#include <gatb/gatb_core.hpp>
#include <boost/config.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/labeled_graph.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/scope_exit.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/graph/subgraph.hpp>
#include <boost/graph/reverse_graph.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/subgraph.hpp>
#include <cstdlib>

#include "global.h"


using namespace std;

class UnitigStats {
private:
    long double qValue;
    long double weight;
    long double normalizedWeight;
    bool valid;
public:
    UnitigStats():valid(false){} //default constructor
    UnitigStats(const PatternFromStats *patternStat, int weightCorretion) {
        if (patternStat) { //valid stuff
            this->qValue = patternStat->qValue;
            this->weight = patternStat->weight * weightCorretion;
            this->normalizedWeight = (weightCorretion == 1 ? patternStat->normalizedWeight : (1 - patternStat->normalizedWeight));
            valid = true;
        }else{
            valid=false; //not valid
        }
    }


    string getQValueAsStr() const {
        stringstream ss;
        ss << scientific;
        if (valid) {
            ss << qValue;
            return ss.str();
        } else {
            return string("NA");
        }
    }

    string getWeightAsStr() const {
        stringstream ss;
        ss << scientific;
        if (valid) {
            ss << weight;
            return ss.str();
        } else {
            return string("NA");
        }
    }

};

//vertex informations
struct VertexInfo {
    string name; //probably represent this as (unitig id, pos)
    int id; //do not use unsigned values
    char strand;
    UnitigStats unitigStats;
};

//edge informations
struct EdgeInfo {
    int id; //do not use unsigned values
    int weight;
};

//some typedefs for making life easier
typedef boost::property<boost::vertex_index_t, std::size_t, VertexInfo> vertex_prop;
typedef boost::property<boost::edge_index_t, std::size_t, EdgeInfo> edge_prop;
typedef boost::adjacency_list <boost::setS, boost::vecS, boost::undirectedS, vertex_prop, edge_prop> adjlist_t;
typedef boost::subgraph< adjlist_t > graph_t;
typedef boost::graph_traits<adjlist_t>::vertex_descriptor MyVertex;
typedef boost::graph_traits<adjlist_t>::edge_descriptor MyEdge;



//class related to get the neighbourhood subgraph
struct TooDistant {
};

class StopWhenVeryDistantFromSourceDijkstraVisitor : public boost::dijkstra_visitor<> {
private:
    vector<int> &distances;
    int maxDistance;
    set <MyVertex> &verticesReachableByTheSourceWithinMaxDistance;
public:
    StopWhenVeryDistantFromSourceDijkstraVisitor(vector<int> &distances, int maxDistance,
                                                 set <MyVertex> &verticesReachableByTheSourceWithinMaxDistance) :
        distances(distances), maxDistance(maxDistance),
        verticesReachableByTheSourceWithinMaxDistance(verticesReachableByTheSourceWithinMaxDistance) { }

    template<class MyVertex, class Graph>
    void finish_vertex(MyVertex v, const Graph &g) {
        //check if we are already going over the max distance - since dijkstra is a greedy algorithm, then we can stop here!
        if (distances[boost::get(boost::vertex_index, g, v)] > maxDistance)
            throw TooDistant();
        else //this node is reachable!
            verticesReachableByTheSourceWithinMaxDistance.insert(v);
    }
};

class GraphWriter {
public:
    //some prints to help
    template<class Graph>
    static string toString(const MyVertex &v, const Graph &graph) {
        stringstream ss;
        ss << graph[v].id << "_" << graph[v].strand;
        return ss.str();
    }

    template<class Graph>
    static string toString(const MyEdge &e, const Graph &graph) {
        stringstream ss;
        ss << toString(source(e, graph), graph) << "\t" << toString(target(e, graph), graph) << " " <<
        graph[e].weight;
        return ss.str();
    }



    template<class T>
    static void writeGraphToFile(const T &graph, const string &fileName) {
        //output the graph
        ofstream graphOut;
        openFileForWriting(fileName, graphOut);

        typedef typename boost::graph_traits<T>::vertex_iterator vertex_iter;
        std::pair <vertex_iter, vertex_iter> vp;
        for (vp = vertices(graph); vp.first != vp.second; ++vp.first) {
            MyVertex v = *(vp.first);
            auto outEdges = out_edges(v, graph);
            for_each(outEdges.first, outEdges.second, [&](const MyEdge &e) {
                graphOut << toString(e, graph) << endl;
            });
        }
        graphOut.close();
    }
};


class generate_output : public Tool {
public:

    // Constructor
    generate_output();

    // Actual job done by the tool is here
    void execute();

    //overriding this in order to exit the tool when finding a problem with the arguments
    IProperties* run (int argc, char* argv[])
    {
        IProperties* toReturn = Tool::run(argc, argv);
        if (!toReturn)
            std::exit(1);
        return toReturn;
    }

};



#endif //KSGATB_GENERATE_OUPUT_H
