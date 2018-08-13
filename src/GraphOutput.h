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
*/

#include <unordered_map>
#include <functional>

#include <set>
#include <stdlib.h> // for exit()
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <boost/regex.hpp>
#include <assert.h> 
#include <gatb/gatb_core.hpp>


#ifndef _GRAPHOUTPUT_H
#define _GRAPHOUTPUT_H

using namespace std;

template<size_t span>
class GraphOutput {
 private:
    typedef typename Kmer<span>::ModelCanonical ModelCanonical;
    FILE *nodes_file,*edges_file;
    Graph *graph;
    string prefix;

    enum LeftOrRight { LEFT=0, RIGHT=1 };
    struct kMinus1_MerInfo {
        long node;
        Strand strand;
        LeftOrRight left_or_right;
        kMinus1_MerInfo(long node, Strand strand, LeftOrRight left_or_right) : node(node), strand(strand), left_or_right(left_or_right) {}
        bool operator<(const kMinus1_MerInfo &other) const {
            if (node != other.node)
                return (node < other.node);
            if (left_or_right != other.left_or_right)
                return left_or_right < other.left_or_right;
            return (strand < other.strand);
        }
    };


    //kMinus1_MerLinks maps the kmer to its unitig id, strand, and if it is the leftest or the rightest kmer of the unitig
    struct CompareKmerCanonical {
      bool operator() (const typename gatb::core::kmer::impl::Kmer<span>::KmerCanonical& lhs,
                       const typename gatb::core::kmer::impl::Kmer<span>::KmerCanonical& rhs) const {
        return lhs.value() < rhs.value();
      }
    };
    map<typename ModelCanonical::Kmer,set<kMinus1_MerInfo>,CompareKmerCanonical > kMinus1_MerLinks;

    //firstAndLastKmers stores the first and the last kmers of a unitig
    vector< pair<typename ModelCanonical::Kmer, typename ModelCanonical::Kmer> > firstAndLastKmers;


public:
    /************************************************************************************************************************/
    /*    Initialize first elements and files  (files are erasing)            */
    /*                              */
    /************************************************************************************************************************/
    GraphOutput(Graph *graph=NULL, const string &prefix="graph") : graph(graph), prefix(prefix){}

    void open(){
        string nodes_file_name=(prefix+".nodes");
        string edges_file_name=(prefix+".edges.dbg");
        nodes_file = fopen(nodes_file_name.c_str(),"w");
        edges_file = fopen(edges_file_name.c_str(),"w");
    }

    void close(){
        fclose(nodes_file);
        fclose(edges_file);
    }

    /************************************************************************************************************************/
    /*      output a single node to a file                                      */
    /*                                                          */
    /************************************************************************************************************************/
    void print_node(long index, char *ascii_node) // output a single node to a file
    {
        fprintf(nodes_file,"%ld\t%s\n",index,ascii_node);
    }


    /************************************************************************************************************************/
    /*      output a single edges to a file                                     */
    /*                                                          */
    /************************************************************************************************************************/
    void print_edge(long index, long id, long id2, string label)
    {
        fprintf(edges_file,"%ld\t%ld\t%s\n",id,id2,label.c_str());
    }

    /************************************************************************************************************************/
    /*      load nodes extremities                                          */
    /*                                                          */
    /************************************************************************************************************************/
    //function that goes through the unitig file and populates kMinus1_MerLinks
    //kMinus1_MerLinks stores as key a (k-1)-mer and as value the (k-1)-mer's unitig id, strand, and if it is the leftest or the rightest (k-1)-mer of the unitig (this for all unitigs. So, if a (k-1)-mer appears as leftest or rightest in n unitigs, we will have one entry and n values)
    //kMinus1_MerLinks stores all the leftest and rightest kmers of each unitig
    void load_nodes_extremities(const string &linear_seqs_name)
    {
        IBank *Nodes = Bank::open((char *)linear_seqs_name.c_str());
        LOCAL (Nodes);

        long nb_nodes = 0;

        // We loop over sequences.
        // for each node, output all the out-edges (in-edges will correspond to out-edges of neighbors)
        ProgressIterator<Sequence> it(*Nodes, "Loading endpoints of unitigs");
        ModelCanonical kMinus1_merModel(graph->getKmerSize()-1);
        ModelCanonical k_merModel(graph->getKmerSize());
        for (it.first(); !it.isDone(); it.next()) {
            string sequence = it.item().toString();

            //here we get the left and the right (k-1)-mer of the unitigs, reverse and forward, and code them
            typename ModelCanonical::Kmer leftkMinus1_mer, rightkMinus1_mer, leftk_mer, rightk_mer;

            //code kmers (the canonical (minimum between fw and rc will be saved))
            leftkMinus1_mer = kMinus1_merModel.codeSeed(sequence.c_str(), Data::ASCII, 0);
            leftk_mer = k_merModel.codeSeed(sequence.c_str(), Data::ASCII, 0);
            rightkMinus1_mer = kMinus1_merModel.codeSeed(sequence.c_str(), Data::ASCII, sequence.length()-(graph->getKmerSize()-1));
            rightk_mer = k_merModel.codeSeed(sequence.c_str(), Data::ASCII, sequence.length()-graph->getKmerSize());

            //kMinus1_MerLinks maps the kmer to its unitig id, strand, and if it is the leftest or the rightest kmer of the unitig
            kMinus1_MerLinks[leftkMinus1_mer].insert(kMinus1_MerInfo(nb_nodes, leftkMinus1_mer.strand(), LEFT));
            kMinus1_MerLinks[rightkMinus1_mer].insert(kMinus1_MerInfo(nb_nodes, rightkMinus1_mer.strand(), RIGHT));

            firstAndLastKmers.push_back(make_pair(leftk_mer, rightk_mer));

            nb_nodes++;
        }
    }


    /************************************************************************************************************************/
    /*      construct node file and edge file for graph file                            */
    /*                                                          */
    /************************************************************************************************************************/
    void construct_graph(string linear_seqs_name)
    {
        IBank *Nodes = Bank::open(linear_seqs_name);
        LOCAL (Nodes);
        long idNodes=0;
        long idEdges=0;

        // We loop over sequences.
        // for each node, output all the out-edges (in-edges will correspond to out-edges of neighbors)
        ProgressIterator<Sequence> it(*Nodes, "Building .nodes and .edges files");
        ModelCanonical kMinus1_merModel(graph->getKmerSize()-1);
        ModelCanonical k_merModel(graph->getKmerSize());
        int index=0;
        for (it.first(); !it.isDone(); it.next(), index++)
        {
            string sequence = it.item().toString();

            typename ModelCanonical::Kmer leftkMinus1_mer, rightkMinus1_mer;

            //code kmers (the canonical (minimum between fw and rc will be saved))
            leftkMinus1_mer = kMinus1_merModel.codeSeed(sequence.c_str(), Data::ASCII, 0);
            Node leftNode = graph->buildNode(k_merModel.toString(firstAndLastKmers[index].first.value()).c_str());

            if (firstAndLastKmers[index].first.strand() != leftNode.strand)
                leftNode = graph->reverse(leftNode);

            rightkMinus1_mer = kMinus1_merModel.codeSeed(sequence.c_str(), Data::ASCII, sequence.length()-(graph->getKmerSize()-1));
            Node rightNode = graph->buildNode(k_merModel.toString(firstAndLastKmers[index].second.value()).c_str());
            if (firstAndLastKmers[index].second.strand() != rightNode.strand)
                rightNode = graph->reverse(rightNode);

            // left edges (are revcomp extensions)
            // get the nodes that has a left kmer or right kmer identical to the kmer stored in leftkMinus1_mer
            typename set<kMinus1_MerInfo>::iterator it;
            for (it = kMinus1_MerLinks[leftkMinus1_mer].begin(); it != kMinus1_MerLinks[leftkMinus1_mer].end(); it++)
            {
                long cur_node = it->node;
                Strand cur_strand = it->strand;
                LeftOrRight cur_left_or_right = it->left_or_right;

                //build the node correctly
                Node cur_GATB_node;
                if (cur_left_or_right==LEFT) {
                    cur_GATB_node = graph->buildNode(k_merModel.toString(firstAndLastKmers[cur_node].first.value()).c_str());
                    if (firstAndLastKmers[cur_node].first.strand() != cur_GATB_node.strand)
                        cur_GATB_node = graph->reverse(cur_GATB_node);
                }else {
                    cur_GATB_node = graph->buildNode(k_merModel.toString(firstAndLastKmers[cur_node].second.value()).c_str());
                    if (firstAndLastKmers[cur_node].second.strand() != cur_GATB_node.strand)
                        cur_GATB_node = graph->reverse(cur_GATB_node);
                }


                /* TODO
                 * I do not understand this well...
                if (cur_node ==nb_els.node) // prevent self loops on same kmer
                     if (sequence.length() == sizeKmer)
                        continue;
                */
                //Replaced the preceding by:
                if (cur_node ==idNodes) // prevent self loops on same kmer
                    continue;
                
                string label = "R";

                if (cur_left_or_right == LEFT)
                {
                    if (cur_strand != leftkMinus1_mer.strand())
                        label+=(string)"F";
                    else
                        continue;
                }
                else
                {
                    if (cur_strand == leftkMinus1_mer.strand())
                        label+=(string)"R";
                    else
                        continue;
                }


                //fix the strands
                Node leftNodeCorrected=leftNode;
                if (label[0]=='R')
                    leftNodeCorrected = graph->reverse(leftNodeCorrected);
                if (label[1]=='R')
                    cur_GATB_node = graph->reverse(cur_GATB_node);

                //check if this edge exist in the graph (maybe it was removed in the sequencing-error removal step)
                vector<Node> neighbours;
                //cout << "Neighbours: " << endl;
                graph->successors(leftNodeCorrected).iterate([&](const Node &node) {
                    neighbours.push_back(node);
                    //cout << graph->toString(node) << endl;
                });
                if (find(neighbours.begin(), neighbours.end(), cur_GATB_node) != neighbours.end()) {
                    print_edge(idEdges, idNodes, cur_node, label);
                    idEdges++;
                }
                else {
                    //DEBUG
                    /*
                    string curGATBNodeAsStr=graph->toString(cur_GATB_node);
                    cout << "Edge skipped LEFT: [" <<
                            graph->toString(leftNodeCorrected) << " --" << (*(curGATBNodeAsStr.rbegin())) <<"--> " << curGATBNodeAsStr << "]" << endl;
                    */
                    //DEBUG
                }
                
            }

            // right edges
            for (it = kMinus1_MerLinks[rightkMinus1_mer].begin(); it != kMinus1_MerLinks[rightkMinus1_mer].end(); it++)
            {
                long cur_node = it->node;
                Strand cur_strand = it->strand;
                LeftOrRight cur_left_or_right = it->left_or_right;

                //build the node correctly
                Node cur_GATB_node;
                if (cur_left_or_right==LEFT) {
                    cur_GATB_node = graph->buildNode(k_merModel.toString(firstAndLastKmers[cur_node].first.value()).c_str());
                    if (firstAndLastKmers[cur_node].first.strand() != cur_GATB_node.strand)
                        cur_GATB_node = graph->reverse(cur_GATB_node);
                }else {
                    cur_GATB_node = graph->buildNode(k_merModel.toString(firstAndLastKmers[cur_node].second.value()).c_str());
                    if (firstAndLastKmers[cur_node].second.strand() != cur_GATB_node.strand)
                        cur_GATB_node = graph->reverse(cur_GATB_node);
                }

                /* TODO
                 * I do not understand this well...
                if (cur_node ==nb_els.node) // prevent self loops on same kmer
                     if (sequence.length() == sizeKmer)
                        continue;
                */
                //Replaced the preceding by:
                if (cur_node ==idNodes) // prevent self loops on same kmer
                    continue;


                string label = "F";
                if (cur_left_or_right == LEFT)
                {
                    if (cur_strand == rightkMinus1_mer.strand())
                        label+=(string)"F";
                    else
                        continue;
                }
                else
                {
                    if (cur_strand != rightkMinus1_mer.strand())
                        label+=(string)"R";
                    else
                        continue;
                }

                //fix the strands
                Node rightNodeCorrected=rightNode;
                if (label[0]=='R')
                    rightNodeCorrected = graph->reverse(rightNodeCorrected);
                if (label[1]=='R')
                    cur_GATB_node = graph->reverse(cur_GATB_node);

                //check if this edge exist in the graph (maybe it was removed in the sequencing-error removal step)
                vector<Node> neighbours;
                //cout << "Neighbours: " << endl;
                graph->successors(rightNodeCorrected).iterate([&](const Node &node) {
                    neighbours.push_back(node);
                    //cout << graph->toString(node) << endl;
                });
                if (find(neighbours.begin(), neighbours.end(), cur_GATB_node) != neighbours.end()) {
                    print_edge(idEdges, idNodes, cur_node, label);
                    idEdges++;
                }
                else {
                    //DEBUG
                    /*
                    string curGATBNodeAsStr=graph->toString(cur_GATB_node);
                    cout << "Edge skipped RIGHT: [" <<
                    graph->toString(rightNodeCorrected) << " --" << (*(curGATBNodeAsStr.rbegin())) <<"--> " << curGATBNodeAsStr << "]" << endl;
                     */
                    //DEBUG
                }
            }
            //print the node
            print_node(idNodes, const_cast<char*>(sequence.c_str()));

            //increase unitig id
            idNodes++;
        }
    }    
};
#endif //_GRAPHOUTPUT_H

