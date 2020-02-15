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

#include "global.h"
#include "map_reads.hpp"
#include "Utils.h"
#include <boost/algorithm/string/predicate.hpp>
#include <map>
#define NB_OF_READS_NOTIFICATION_MAP_AND_PHASE 10 //Nb of reads that the map and phase must process for notification
using namespace std;


void mapReadToTheGraphCore(const string &read, const Graph &graph, const vector< UnitigIdStrandPos > &nodeIdToUnitigId,
                           map<int,int> &unitigIdToCount) {
    int lastUnitig=-1;

    //goes through all nodes/kmers of the read
    if (read.size() >= graph.getKmerSize()) {
        for (int i = 0; i < read.size() - graph.getKmerSize() + 1; i++) {
            string LRKmer = string(read.c_str() + i, graph.getKmerSize());

            //TODO: From my tests, if you use buildNode with a Kmer containing an 'N', it will build a node with all Ns replaced by G
            //TODO: However, Kmers containing Ns are NOT included in any unitig (i.e. the graph builder does not convert an N to a G, and build the graph. It simply disregards kmers containing Ns)
            //TODO: this is a sanity check to also discard Kmers containing Ns
            //TODO: check this with GATB team

            //TODO: UPDATE
            //TODO: Magali had a dataset where we had the base 'K' in the fasta file
            //TODO: So I am just discarding all reads that are not composed by ACGT
            if (!boost::all(LRKmer, [](char c) -> bool {
                return c=='A' || c=='a' || c=='C' || c=='c' || c=='G' || c=='g' ||c=='T' || c=='t';
            })) {
                continue;
            }

            //build the node
            Node node = graph.buildNode(LRKmer.c_str());

            //get the unitig localization of this kmer
            u_int64_t index = graph.nodeMPHFIndex(node);
            UnitigIdStrandPos unitigIdStrandPos=nodeIdToUnitigId[index];

            if (lastUnitig != unitigIdStrandPos.unitigId) {
                if (unitigIdToCount.find(unitigIdStrandPos.unitigId) == unitigIdToCount.end() )
                    unitigIdToCount[unitigIdStrandPos.unitigId]=0;
                unitigIdToCount[unitigIdStrandPos.unitigId]++;
                lastUnitig = unitigIdStrandPos.unitigId;
            }
        }
    }
}

void mapReadToTheGraph(const string &read, int readfileIndex, unsigned long readIndex, const Graph &graph,
                       const vector< UnitigIdStrandPos > &nodeIdToUnitigId, map<int,int> &unitigIdToCount) {
    //map the read
    mapReadToTheGraphCore(read, graph, nodeIdToUnitigId, unitigIdToCount);
}


// We define a functor that will be cloned by the dispatcher
struct MapAndPhase
{
    const vector<string> &allReadFilesNames;
    const Graph& graph;
    const string &outputFolder;
    const string &tmpFolder;
    uint64_t &nbOfReadsProcessed;
    ISynchronizer* synchro;
    vector< UnitigIdStrandPos > &nodeIdToUnitigId;
    int nbContigs;

    struct MapAndPhaseIteratorListener : public IteratorListener {
        uint64_t &nbOfReadsProcessed;
        ISynchronizer* synchro;
        MapAndPhaseIteratorListener(uint64_t &nbOfReadsProcessed, ISynchronizer* synchro) :
            nbOfReadsProcessed(nbOfReadsProcessed), synchro(synchro){}

        virtual void inc (u_int64_t ntasks_done) {
            // We lock the synchronizer
            synchro->lock ();

            nbOfReadsProcessed+=NB_OF_READS_NOTIFICATION_MAP_AND_PHASE;
            cerr << '\r' << nbOfReadsProcessed << " reads processed.";
            cerr.flush();

            // We unlock the synchronizer
            synchro->unlock ();
        }
    };

    MapAndPhase (const vector<string> &allReadFilesNames, const Graph& graph,
                 const string &outputFolder, const string &tmpFolder, uint64_t &nbOfReadsProcessed, ISynchronizer* synchro,
                 vector< UnitigIdStrandPos > &nodeIdToUnitigId, int nbContigs) :
        allReadFilesNames(allReadFilesNames), graph(graph), outputFolder(outputFolder), tmpFolder(tmpFolder),
        nbOfReadsProcessed(nbOfReadsProcessed), synchro(synchro), nodeIdToUnitigId(nodeIdToUnitigId),
        nbContigs(nbContigs){}

    void operator()(int i) {
        // We declare an input Bank and use it locally
        IBank *inputBank = Bank::open(allReadFilesNames[i]);
        LOCAL(inputBank);

        // Create and use a progress iterator
        MapAndPhaseIteratorListener* mapAndPhaseIteratorListener = new MapAndPhaseIteratorListener(nbOfReadsProcessed, synchro);
        SubjectIterator <Sequence> it(inputBank->iterator(), NB_OF_READS_NOTIFICATION_MAP_AND_PHASE, mapAndPhaseIteratorListener);

        //XU_strain_i = how many times each unitig map to a strain
        ofstream mappingOutputFile;
        openFileForWriting(tmpFolder+string("/XU_strain_")+to_string(i), mappingOutputFile);

        // We loop over sequences.
        unsigned long readIndex = 0;
        map<int,int> unitigIdToCount; //TODO: change this to a vector
        for (it.first(); !it.isDone(); it.next()) {
            string read = (it.item()).toString();
            //transform the read to upper case
            for (int j=0;j<read.size();j++)
                read[j]=toupper(read[j]);

            //map this read to the graph
            mapReadToTheGraph(read, i, readIndex, graph, nodeIdToUnitigId, unitigIdToCount);

            readIndex++;
        }

        //output info for mapping - the number of times the unitig appear in the strain
        for (int i=0;i<nbContigs;i++) {
            if (unitigIdToCount.find(i) == unitigIdToCount.end() )
                mappingOutputFile << "0 ";
            else
                mappingOutputFile << "1 ";
        }

        mappingOutputFile.close();
    }
};

map_reads::map_reads ()  : Tool ("map_reads") //give a name to our tool
{
    populateParser(this);
}

//pattern is the unitig line
map< vector<int>, vector<int> > getUnitigsWithSamePattern (const vector< vector<int> > &XU, int nbContigs) {
    map< vector<int>, vector<int> > pattern2Unitigs;

    for (int i=0;i<XU.size();i++) { //goes through all unitigs
        if (pattern2Unitigs.count(XU[i])==0) { //pattern of unitig i is not in pattern2Unitigs
            //create a vector with unitig i
            vector<int> unitigs;
            unitigs.push_back(i);

            //insert this pattern and his new set to the map
            pattern2Unitigs.insert(make_pair(XU[i], unitigs));
        } else {
            //pattern of unitig i is already in pattern2Unitigs, just add
            pattern2Unitigs[XU[i]].push_back(i);
        }
    }

    return pattern2Unitigs;
}

void generate_XU(const string &filename, const string &nodesFile, const vector< vector<int> > &XU) {
    ofstream XUFile;
    openFileForWriting(filename, XUFile);

    //open file with list of node sequences
    ifstream nodesFileReader;
    openFileForReading(nodesFile, nodesFileReader);
    int id;
    string seq;

    for (int i=0;i<XU.size();i++) {
        // print the unitig sequence
        nodesFileReader >> id >> seq;
        XUFile << seq << " |";

        //print the strains present
        for (int j=0;j<XU[i].size();j++)
            if (XU[i][j] > 0) {
                XUFile << " " << (*strains)[j].id << ":" << XU[i][j];
            }
        XUFile << endl;
    }
    nodesFileReader.close();
    XUFile.close();
}

void generate_unique_id_to_original_ids(const string &filename,
                                        const map< vector<int>, vector<int> > &pattern2Unitigs) {
    ofstream uniqueIdToOriginalIdsFile;
    openFileForWriting(filename, uniqueIdToOriginalIdsFile);

    //for each pattern
    int i=0;
    auto it=pattern2Unitigs.begin();
    for (;it!=pattern2Unitigs.end();++it, ++i) {
        //print the id of this pattern
        uniqueIdToOriginalIdsFile << i << " = ";

        //and the unitigs in it
        for (auto id : it->second)
            uniqueIdToOriginalIdsFile << id << " ";

        uniqueIdToOriginalIdsFile << endl;
    }
    uniqueIdToOriginalIdsFile.close();
}

void generate_XU_unique(const string &filename, const vector< vector<int> > &XU,
                        const map< vector<int>, vector<int> > &pattern2Unitigs){
    ofstream XUUnique;
    openFileForWriting(filename, XUUnique);

    //print the header
    XUUnique << "pattern_id";
    for (const auto &strain : (*strains))
        XUUnique << " " << strain.id;
    XUUnique << endl;

    //for each pattern
    int i=0;
    auto it=pattern2Unitigs.begin();
    for (;it!=pattern2Unitigs.end();++it, ++i) {
        //print the id of this pattern
        XUUnique << i;

        //print the pattern
        for (const auto &v : it->first)
            XUUnique << " " << v;
        XUUnique << endl;
    }
    XUUnique.close();
}

//generate the pyseer input
void generatePyseerInput (const vector <string> &allReadFilesNames,
                          const string &outputFolder, const string &tmpFolder,
                          int nbContigs) {
    //Generate the XU (the pyseer input - the unitigs are rows with strains present)
    //XU_unique is XU is in matrix form (for Rtab input) with the duplicated rows removed
    cerr << endl << endl << "[Generating pyseer input]..." << endl;

    //Create XU
    vector< vector<int> > XU(nbContigs);
    for (auto & v : XU)
        v.resize(allReadFilesNames.size());

    //populate XU
    for (int j=0; j<allReadFilesNames.size(); j++) {
        ifstream inputFile;
        openFileForReading(tmpFolder+string("/XU_strain_")+to_string(j), inputFile);
        for (int i = 0; i < nbContigs; i++)
            inputFile >> XU[i][j];
        inputFile.close();
    }

    //create the files for pyseer
    {
        generate_XU(outputFolder+string("/unitigs.txt"), outputFolder+string("/graph.nodes"), XU);
        map< vector<int>, vector<int> > pattern2Unitigs = getUnitigsWithSamePattern(XU, nbContigs);
        generate_unique_id_to_original_ids(outputFolder+string("/unitigs.unique_rows_to_all_rows.txt"), pattern2Unitigs);
        generate_XU_unique(outputFolder+string("/unitigs.unique_rows.Rtab"), XU, pattern2Unitigs);
    }

    cerr << "[Generating pyseer input] - Done!" << endl;


}

void map_reads::execute ()
{
    //get the parameters
    string outputFolder = stripLastSlashIfExists(getInput()->getStr(STR_OUTPUT));
    string tmpFolder = outputFolder+string("/tmp");
    string longReadsFile = tmpFolder+string("/readsFile");
    int nbCores = getInput()->getInt(STR_NBCORES);

    //get the nbContigs
    int nbContigs = getNbLinesInFile(outputFolder+string("/graph.nodes"));

    //Do the Mapping
    //Maps all the reads back to the graph

    //get all the read files' name
    vector <string> allReadFilesNames = getVectorStringFromFile(longReadsFile);

    // We create an iterator over an integer range
    Range<int>::Iterator allReadFilesNamesIt(0, allReadFilesNames.size() - 1);

    //synchronizer object
    ISynchronizer *synchro = System::thread().newSynchronizer();

    // We create a dispatcher configured for 'nbCores' cores.
    Dispatcher dispatcher(nbCores, 1);

    cerr << "[Starting mapping process... ]" << endl;
    cerr << "Using " << nbCores << " cores to map " << allReadFilesNames.size() << " read files." << endl;

    // We iterate the range.  NOTE: we could also use lambda expression (easing the code readability)
    uint64_t nbOfReadsProcessed = 0;
    dispatcher.iterate(allReadFilesNamesIt,
                       MapAndPhase(allReadFilesNames, *graph, outputFolder, tmpFolder, nbOfReadsProcessed, synchro,
                                   *nodeIdToUnitigId, nbContigs));

    cerr << endl << "[Mapping process finished!]" << endl;

    //generate the pyseer input
    generatePyseerInput(allReadFilesNames, outputFolder, tmpFolder, nbContigs);

    cout << "Number of unique patterns: " << getNbLinesInFile(outputFolder+string("/unitigs.unique_rows.Rtab")) << endl;

    // Remove the global graph pointer, otherwise its destructor is called
    // twice by GATB (after main) giving a HDF5 error
    delete graph;

    //graph.~Graph();
    //delete nodeIdToUnitigId;

    //clean-up - saving some disk space
    //remove temp directory
    boost::filesystem::remove_all(tmpFolder);

    cerr.flush();
}
