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

#include "build_dbg.hpp"
#include "global.h"
#include "GraphOutput.h"
#include "version.h"

using namespace std;

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
// The Tool constructor allows to give a name to our tool.
// This name appears when one gets help in the command line or in the final output
build_dbg::build_dbg ()  : Tool ("build_dbg") //give a name to our tool
{
    setVersion([](void* whatever) {
        cout << "unitig-counter v" << VERSION << endl;
    });
    populateParser(this);
}

void buildSequence (
        const Graph& graph,
        const Node& startingNode,
        size_t length,
        size_t nbContigs,
        const string& consensusRight,
        const string& consensusLeft,
        Sequence& seq
)
{
    /** Shortcuts. */
    Data&  data     = seq.getData();

    /** We set the sequence comment. */
    stringstream ss1;
    ss1 << nbContigs << "__len__" << length << " ";
    seq._comment = ss1.str();

    /** We set the data length. */
    seq.getData().resize (length);

    //We fill the data
    string finalSequence = consensusLeft + graph.toString (startingNode) + consensusRight;
    for (size_t i=0;i<finalSequence.size();i++)
        data[i] = finalSequence[i];
}

//Params:
//graph: the GATB graph
//node: a given node
//unitigPart: the exact substring of the unitig that refer to this node
//returns:
//'F' if the FORWARD node represented by the input node appears in the FORWARD sequence of the unitig
//'R' if the FORWARD node represented by the input node appears in the REVCOMPL sequence of the unitig
char getUnitigStrandTheForwardNodeMapsTo(const gatb::core::debruijn::impl::Graph& graph, const Node &node, const string &unitigPart) {
    //get the forward node of the node given as parameter
    Node forwardNode = node;
    if (node.strand == gatb::core::kmer::STRAND_REVCOMP) //if the node given is not the forward node
        forwardNode = graph.reverse(forwardNode); //reverse it
    if (unitigPart == graph.toString(forwardNode))
        return 'F';
    else if (reverse_complement(unitigPart) == graph.toString(forwardNode))
        return 'R';
    else
        throw new runtime_error("Fatal bug on build_dbg.cpp::getUnitigStrandTheForwardNodeMapsTo()");
}

void construct_linear_seqs (const gatb::core::debruijn::impl::Graph& graph, const string& linear_seqs_name,
                            vector< UnitigIdStrandPos >& nodeIdToUnitigId)
{
    using namespace gatb::core::debruijn::impl;
    using namespace gatb::core::tools::misc::impl;

    IBank* outputBank = new BankFasta (linear_seqs_name.c_str());
    LOCAL (outputBank);

    // We create a Terminator object - this will mark the nodes that are already built
    MPHFTerminator terminator (graph);

    // We create a BranchingTerminator object - this will mark the nodes where to stop the traversal
    BranchingTerminator branchingTerminator(graph);

    // We create a Traversal instance to traverse unitigs
    Traversal* traversal = Traversal::create (TRAVERSAL_UNITIG, graph, branchingTerminator);
    LOCAL (traversal);

    Path consensusRight;
    Path consensusLeft;
    Sequence seq (Data::ASCII);
    u_int64_t nbContigs=0;
    BankFasta::setDataLineSize(0);

    //We loop through the nodes and build the unitigs
    ProgressGraphIterator<Node, ProgressTimerAndSystem> it (graph.iterator(), "Graph: building unitigs");
    for (it.first(); !it.isDone(); it.next()) {
        auto &startingNode = it.item();

        if (terminator.is_marked(startingNode))
            continue;

        auto reversedNode = graph.reverse(startingNode);
        int lenRight = traversal->traverse (startingNode, DIR_OUTCOMING, consensusRight);
        int lenLeft = traversal->traverse (reversedNode, DIR_OUTCOMING, consensusLeft);
        int lenTotal = graph.getKmerSize() + lenRight + lenLeft;

        //mark the traversed nodes
        terminator.mark(startingNode);
        auto currentNode = startingNode;
        for_each(consensusRight.path.begin(), consensusRight.path.end(), [&](const Nucleotide &nucleotide) {
            currentNode = graph.successor(currentNode, nucleotide);
            terminator.mark(currentNode);
        });
        currentNode = reversedNode;
        for_each(consensusLeft.path.begin(), consensusLeft.path.end(), [&](const Nucleotide &nucleotide) {
            currentNode = graph.successor(currentNode, nucleotide);
            terminator.mark(currentNode);
        });

        // We get the unitig strings
        string consensusLeftStr;
        {
            stringstream ss;
            ss << consensusLeft;
            consensusLeftStr=reverse_complement(ss.str());
        }

        string consensusRightStr;
        {
            stringstream ss;
            ss << consensusRight;
            consensusRightStr=ss.str();
        }


        /** We create the contig sequence. */
        buildSequence(graph, startingNode, lenTotal, nbContigs, consensusRightStr, consensusLeftStr, seq);

        //associate the node id to its unitig id
        //Note: GATB kmer is any kmer... It is not the canonical one (i.e. smaller one). Maybe is the one that was added...
        //Anyway, for a given kmer, we have it as forward node and reverse node. The forward node is what it counts (and it is not necessarily the canonical kmer)
        char strand = getUnitigStrandTheForwardNodeMapsTo(graph, startingNode, seq.toString().substr(lenLeft, graph.getKmerSize()));
        nodeIdToUnitigId[graph.nodeMPHFIndex(startingNode)] = UnitigIdStrandPos(nbContigs,
                                                                                strand,
                                                                                (strand=='F' ? lenLeft : lenRight),
                                                                                lenTotal, graph.getKmerSize());
        currentNode = startingNode;
        int i=1;
        for_each(consensusRight.path.begin(), consensusRight.path.end(), [&](const Nucleotide &nucleotide) {
            currentNode = graph.successor(currentNode, nucleotide);
            strand = getUnitigStrandTheForwardNodeMapsTo(graph, currentNode, seq.toString().substr(lenLeft+i, graph.getKmerSize()));
            nodeIdToUnitigId[graph.nodeMPHFIndex(currentNode)] = UnitigIdStrandPos(nbContigs,
                                                                                   strand,
                                                                                   (strand=='F' ? lenLeft+i : lenRight-i),
                                                                                   lenTotal, graph.getKmerSize());
            i++;
        });

        currentNode = graph.reverse(startingNode);
        i=-1;
        for_each(consensusLeft.path.begin(), consensusLeft.path.end(), [&](const Nucleotide &nucleotide) {
            currentNode = graph.successor(currentNode, nucleotide);
            strand = getUnitigStrandTheForwardNodeMapsTo(graph, currentNode, seq.toString().substr(lenLeft+i, graph.getKmerSize()));
            nodeIdToUnitigId[graph.nodeMPHFIndex(currentNode)] = UnitigIdStrandPos(nbContigs,
                                                                                   strand,
                                                                                   (strand=='F' ? lenLeft+i : lenRight-i),
                                                                                   lenTotal, graph.getKmerSize());
            i--;
        });


        /** We add the sequence into the output bank. */
        outputBank->insert (seq);

        //increase the number of contigs
        nbContigs += 1;
    }

    outputBank->flush ();
}


class EdgeConstructionVisitor : public boost::static_visitor<>    {
private:
    const string& linear_seqs_name;

public:
    EdgeConstructionVisitor (const string &linear_seqs_name) : linear_seqs_name(linear_seqs_name) {}
    template<size_t span>
    void operator() (GraphOutput<span>& graphOutput) const
    {
        graphOutput.open();
        graphOutput.load_nodes_extremities(linear_seqs_name);
        graphOutput.construct_graph(linear_seqs_name);
        graphOutput.close();
    }
};


/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void build_dbg::execute ()
{
    cerr << "Step 1. Building DBG and mapping strains on the DBG..." << endl;
    checkParametersBuildDBG(this);
    if (run2) return;

    //get the parameters
    int kmerSize = getInput()->getInt(STR_KSKMER_SIZE);

    //TODO: seeveral questions are still unclear if we use the Freq count mode (how to run bugwas, the coloring, etc...). For now I am disabling this option
    //string countMode = getInput()->getStr(STR_COUNT_MODE);
    //presenceAbsenceCountMode = (countMode=="01");
    presenceAbsenceCountMode = true;
    //TODO: seeveral questions are still unclear if we use the Freq count mode (how to run bugwas, the coloring, etc...). For now I am disabling this option

    //create the step1 folder in the outputfolder
    string outputFolder = stripLastSlashIfExists(getInput()->getStr(STR_OUTPUT))+string("/unitigs");
    createFolder(outputFolder);

    //create the tmp folder of step1
    string tmpFolder = outputFolder+string("/tmp");
    createFolder(tmpFolder);

    int nbCores = getInput()->getInt(STR_NBCORES);

    //create the reads file
    string readsFile(tmpFolder+string("/readsFile"));
    Strain::createReadsFile(readsFile, strains);

    //Builds the DBG using GATB
    //TODO: by using create() and assigning to a Graph object, the copy constructor does a shallow or deep copy??
    graph = gatb::core::debruijn::impl::Graph::createAsPointer("-in %s -kmer-size %d -abundance-min 0 -out %s/graph -nb-cores %d",
                                                          readsFile.c_str(), kmerSize, outputFolder.c_str(), nbCores);


    // Finding the unitigs
    //nodeIdToUnitigId translates the nodes that are stored in the GATB graph to the id of the unitigs together with the unitig strand
    nodeIdToUnitigId = new vector< UnitigIdStrandPos >((size_t)graph->getInfo()["kmers_nb_solid"]->getInt()); //map nodeMPFHIndex() to unitigIds and strand
    string linear_seqs_name = outputFolder+"/graph.unitigs";
    construct_linear_seqs (*graph, linear_seqs_name, *nodeIdToUnitigId);

    //builds and outputs .nodes and .edges.dbg files
    typedef boost::variant <
        GraphOutput<KMER_SPAN(0)>,
        GraphOutput<KMER_SPAN(1)>,
        GraphOutput<KMER_SPAN(2)>,
        GraphOutput<KMER_SPAN(3)>
    >  GraphOutputVariant;

    GraphOutputVariant graphOutput;
    if (kmerSize < KMER_SPAN(0))  {  graphOutput = GraphOutput<KMER_SPAN(0)>(graph, outputFolder+string("/graph")); }
    else if (kmerSize < KMER_SPAN(1))  {  graphOutput = GraphOutput<KMER_SPAN(1)>(graph, outputFolder+string("/graph")); }
    else if (kmerSize < KMER_SPAN(2))  {  graphOutput = GraphOutput<KMER_SPAN(2)>(graph, outputFolder+string("/graph")); }
    else if (kmerSize < KMER_SPAN(3))  {  graphOutput = GraphOutput<KMER_SPAN(3)>(graph, outputFolder+string("/graph")); }
    else { throw gatb::core::system::Exception ("Graph failure because of unhandled kmer size %d", kmerSize); }
    boost::apply_visitor (EdgeConstructionVisitor(linear_seqs_name),  graphOutput);

    //save disk space
    remove(linear_seqs_name.c_str());

    //print some stats
    cout << "################################################################################" << endl;
    cout << "Stats: " << endl;
    cout << "Number of kmers: " << graph->getInfo()["kmers_nb_solid"]->getInt() << endl;
    cout << "Number of unitigs: " << getNbLinesInFile(outputFolder+string("/graph.nodes")) << endl;
    cout << "################################################################################" << endl;
}
