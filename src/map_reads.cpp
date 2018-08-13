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

#include "global.h"
#include "map_reads.hpp"
#include "Utils.h"
#include "PhenoCounter.h"
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
                mappingOutputFile << (presenceAbsenceCountMode ? 1 : unitigIdToCount[i]) << " ";
        }

        mappingOutputFile.close();
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
// The Tool constructor allows to give a name to our tool.
// This name appears when one gets help in the command line or in the final output
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

void generate_XU(const string &filename, const vector< vector<int> > &XU) {
    ofstream XUFile;
    openFileForWriting(filename, XUFile);

    //print the header
    XUFile << "ps";
    for (const auto &strain : (*strains))
        XUFile << " " << strain.id;
    XUFile << endl;

    for (int i=0;i<XU.size();i++) {
        //print the unitig id
        XUFile << (i);

        //print the frequences/binaries
        for (int j=0;j<XU[i].size();j++)
            XUFile << " " << XU[i][j];
        XUFile << endl;
    }
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

void generate_unique_id_to_original_ids(const string &uniqueIdToOriginalIdsFilename,
                                        const string &gemmaPatternToNbUnitigsFilename,
                                        const string &gemmaUnitigToPatternFilename,
                                        const map< vector<int>, vector<int> > &pattern2Unitigs) {
    ofstream uniqueIdToOriginalIdsFile;
    openFileForWriting(uniqueIdToOriginalIdsFilename, uniqueIdToOriginalIdsFile);
    ofstream gemmaPatternToNbUnitigsFile;
    openFileForWriting(gemmaPatternToNbUnitigsFilename, gemmaPatternToNbUnitigsFile);
    ofstream gemmaUnitigToPatternFile;
    openFileForWriting(gemmaUnitigToPatternFilename, gemmaUnitigToPatternFile);

    //for each pattern
    int i=0;
    auto it=pattern2Unitigs.begin();
    for (;it!=pattern2Unitigs.end();++it, ++i) {
        //goes through the unitigos of this file
        for (auto id : it->second) {
            //uniqueIdToOriginalIdsFile
            //print all the unitigs of a pattern i in line i
            uniqueIdToOriginalIdsFile << id << " ";
            //uniqueIdToOriginalIdsFile

            //gemmaUnitigToPatternFile
            gemmaUnitigToPatternFile << id << " " << i << endl;
            //gemmaUnitigToPatternFile
        }
        uniqueIdToOriginalIdsFile << endl;
        //uniqueIdToOriginalIdsFile



        //gemmaPatternToNbUnitigsFile
        gemmaPatternToNbUnitigsFile << i << " " << (it->second).size() << endl;
        //gemmaPatternToNbUnitigsFile
    }

    uniqueIdToOriginalIdsFile.close();
    gemmaPatternToNbUnitigsFile.close();
    gemmaUnitigToPatternFile.close();
}

void generate_XU_unique(const string &filename, const vector< vector<int> > &XU,
                        const map< vector<int>, vector<int> > &pattern2Unitigs){
    ofstream XUUnique;
    openFileForWriting(filename, XUUnique);

    //print the header
    XUUnique << "ps";
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

//generate the bugwas input
void generateBugwasInput (const vector <string> &allReadFilesNames, const string &outputFolder, const string &tmpFolder, int nbContigs) {
    //Generate the XU (the bugwas input - the matrix where the unitigs are rows and the strains are columns)
    //XU_unique is XU with the duplicated rows removed
    cerr << endl << endl << "[Generating bugwas and gemma input]..." << endl;

    //create the ID and Phenotype file
    Strain::createIdPhenoFile(outputFolder+string("/bugwas_input.id_phenotype"), strains);

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

    //create a binary XU
    vector< vector<int> > XUbinary(XU);

    //creates also a file saying if the unitig was inverted (-1) or not (1)
    //this is a multiplicative factor that will correct the weight (estimated effect) from the statistical test
    ofstream weightCorrectionStream;
    openFileForWriting(outputFolder+string("/weight_correction"), weightCorrectionStream);

    for (int i=0;i<XUbinary.size();i++) {
        //1. Transform frequency to binary
        for (int j = 0; j < XUbinary[i].size(); j++)
            XUbinary[i][j] = ((int)((bool)(XUbinary[i][j])));

        //2. Transform to the enconding where 0 is the major allele and 1 is the minor one
        //count how many 0s and 1s we have
        int count0=0;
        int count1=0;
        for (int j=0;j<XUbinary[i].size();j++) {
            if (XUbinary[i][j]==0) count0++;
            if (XUbinary[i][j]==1) count1++;
        }

        //0 must be the major allele (we need to have more 0s than 1s). If it is not, 0 and 1 must be inverted
        int iMustInvert = ((int)(count0 < count1));

        //re-assign the values to XUbinary
        for (int j=0;j<XUbinary[i].size();j++)
            XUbinary[i][j] = ((XUbinary[i][j]+iMustInvert)%2);

        weightCorrectionStream << (iMustInvert ? -1 : 1) << endl;
    }
    weightCorrectionStream.close();

    //create the files for bugwas
    /*
     * TODO: THIS IS NOT CREATED RIGHT NOW BECAUSE WE CANNOT PROCESS IT - BUGWAS, FOR THE MOMENT, JUST ACCEPT THE 0/1 (BINARY) FILES
     * TODO: PUT THIS BACK WHEN WE ARE ABLE TO DO IT
     * TODO: IF THE COUNT MODE IS FREQ, ONLY THE BINARY VERSION IS USED FOR BUGWAS - FIX THIS
    {
        generate_XU(outputFolder+string("/bugwas_input.all_rows.frequency"), XU);
        map< vector<int>, vector<int> > pattern2Unitigs = getUnitigsWithSamePattern(XU, nbContigs);
        generate_unique_id_to_original_ids(outputFolder+string("/bugwas_input.unique_rows_to_all_rows.frequency"), pattern2Unitigs);
        generate_XU_unique(outputFolder+string("/bugwas_input.unique_rows.frequency"), XU, pattern2Unitigs);
    }
     */

    //create the files for bugwas - binary ones
    {
        generate_XU(outputFolder+string("/bugwas_input.all_rows.binary"), XUbinary);
        map< vector<int>, vector<int> > pattern2Unitigs = getUnitigsWithSamePattern(XUbinary, nbContigs);
        generate_unique_id_to_original_ids(outputFolder+string("/bugwas_input.unique_rows_to_all_rows.binary"),
                                           outputFolder+string("/gemma_input.pattern_to_nb_of_unitigs.binary"),
                                           outputFolder+string("/gemma_input.unitig_to_pattern.binary"),
                                           pattern2Unitigs);
        generate_XU_unique(outputFolder+string("/bugwas_input.unique_rows.binary"), XUbinary, pattern2Unitigs);
    }
    cerr << "[Generating bugwas and gemma input] - Done!" << endl;


    //create the file showing the overall frequencies of each unitig
    cerr << "[Generating the frequency files...]" << endl;
    vector<PhenoCounter > unitigs2PhenoCounter(nbContigs);
    for(int strainIndex=0;strainIndex<allReadFilesNames.size();strainIndex++) {
        ifstream unitigCountForStrain;
        openFileForReading(tmpFolder+string("/XU_strain_")+to_string(strainIndex), unitigCountForStrain);

        for (int unitigIndex=0; unitigIndex<nbContigs; unitigIndex++) {
            int count;
            unitigCountForStrain >> count;
            if ((*strains)[strainIndex].phenotype=="0")
                unitigs2PhenoCounter[unitigIndex].increasePheno0(count);
            else if ((*strains)[strainIndex].phenotype=="1")
                unitigs2PhenoCounter[unitigIndex].increasePheno1(count);
            else if ((*strains)[strainIndex].phenotype=="NA")
                unitigs2PhenoCounter[unitigIndex].increaseNA(count);
            else
                throw runtime_error("[FATAL ERROR] on map_reads::execute () [Generating the frequency files...]");
        }
        unitigCountForStrain.close();
    }

    ofstream frequencyFile;
    openFileForWriting(outputFolder+string("/frequency_unitig_to_total_pheno0_pheno1_NA_count"), frequencyFile);
    for (int unitigIndex=0; unitigIndex<nbContigs; unitigIndex++)
        frequencyFile << unitigs2PhenoCounter[unitigIndex].getTotal() << " " << unitigs2PhenoCounter[unitigIndex].getPheno0() << " "
                      << unitigs2PhenoCounter[unitigIndex].getPheno1() << " " << unitigs2PhenoCounter[unitigIndex].getNA() << endl;
    frequencyFile.close();


    //create a file with the total nb of Pheno0 and Pheno1 strains
    Strain::createFileWithAmountOfStrainsInEachPheno(outputFolder+string("/total_nb_of_strains_in_each_pheno"), strains);

    cerr << "[Generating the frequency files...] - Done!" << endl;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void map_reads::execute ()
{
    if (skip1) return;

    //get the parameters
    string outputFolder = stripLastSlashIfExists(getInput()->getStr(STR_OUTPUT))+string("/step1");
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

    //generate the bugwas input
    generateBugwasInput(allReadFilesNames, outputFolder, tmpFolder, nbContigs);

    //after the mapping, free some memory that will not be needed anymore
    delete graph;
    delete nodeIdToUnitigId;

    //clean-up - saving some disk space
    //remove temp directory
    boost::filesystem::remove_all(tmpFolder);
    //remove GATB's graph file
    remove((outputFolder+string("/graph.h5")).c_str());

    cerr << endl << "[Mapping process finished!]" << endl;
    cerr.flush();
}
