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

#include "generate_output.h"
#include "Utils.h"
#include "global.h"
#include <cstdlib>
#include <ctime>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/copy.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string.hpp>
#include "Blast.h"
#include <algorithm>
#define BOOST_NO_CXX11_SCOPED_ENUMS
#include <boost/filesystem.hpp>
#undef BOOST_NO_CXX11_SCOPED_ENUMS
#include "version.h"

using namespace std;

//variables regarding the size of the nodes
int minSize = 50;
int maxSize = 400;



// The Tool constructor allows to give a name to our tool.
// This name appears when one gets help in the command line or in the final output
generate_output::generate_output ()  : Tool ("generate_output") //give a name to our tool
{
  populateParser(this);
}

void generate_output::createIndexFile(int numberOfComponents, const string &visualisationsFolder, const string &step2OutputFolder, const vector<vector<MyVertex> > &nodesInComponent, graph_t& newGraph,
                     map<int, AnnotationRecord > &idComponent2Annotations, const vector<const PatternFromStats*> &unitigToPatternStats,
                     const vector<int> &selectedUnitigs, int nbCores) {
  cerr << "[Creating index file...]" << endl;

  //create the thumbnails in a multi-threaded way
  // We create an iterator over an integer range
  Range<int>::Iterator rangeOfComponentsIt(0, numberOfComponents - 1);

  // We create a dispatcher configured for 'nbCores' cores.
  Dispatcher dispatcher(nbCores, 1);

  // We iterate the range
  cerr << "[Rendering thumbnails...]" << endl;
  dispatcher.iterate(rangeOfComponentsIt, [&](int i) {
      string HTMLFile(boost::filesystem::canonical(visualisationsFolder+"/components/comp_"+std::to_string(i)+".html").string());
      string PNGFile = HTMLFile+".png";

      if (noPreview){
        boost::filesystem::copy_file(DBGWAS_lib+"/lib/resources/nopreview.png", PNGFile, boost::filesystem::copy_option::overwrite_if_exists);
      }
      else {
        stringstream commandSS;
        commandSS << phantomjsPath << " " << DBGWAS_lib << "/render_graph.js " << HTMLFile << " " << PNGFile;
        executeCommand(commandSS.str(), false);
      }

  });
  cerr << "[Rendering thumbnails...] - Done!" << endl;



  //create the index
  //create the object previews for each component
  vector<ObjectPreview> previews;

  //get the template preview
  string templatePreview="";
  {
    auto indexTableTemplateAsStringVector = getVectorStringFromFile(DBGWAS_lib + "/index_table_template.html");
    for (const auto &line : indexTableTemplateAsStringVector)
      templatePreview += line;
  }

  for (int i=0;i<numberOfComponents;i++) {
    string idString = std::to_string(i);
    string annotationsSQL=idComponent2Annotations[i].getSQLRepresentationForIndexPage();

    //compute the number of nodes in this component and the number of significant nodes
    int nbNodes = nodesInComponent[i].size();
    int nbOfSignificantNodes = count_if(nodesInComponent[i].begin(), nodesInComponent[i].end(), [&](const MyVertex &node) {
        return find(selectedUnitigs.begin(), selectedUnitigs.end(), newGraph[node].id) != selectedUnitigs.end();
    });


    //get the lowest qvalue of the nodes in the component
    long double lowestQValue = std::numeric_limits<long double>::max();
    for(const auto &node : nodesInComponent[i]) {
      int id = newGraph[node].id;

      //check if the patterns exists
      if (unitigToPatternStats[id])
        lowestQValue=min(lowestQValue, unitigToPatternStats[id]->qValue);
    }

    //add the true values to this preview
    string thisPreview(templatePreview);
    boost::replace_all(thisPreview, "<nb_unitigs>", to_string(nbNodes));
    boost::replace_all(thisPreview, "<nb_sig_unitigs>", to_string(nbOfSignificantNodes));
    boost::replace_all(thisPreview, "<id>", idString);

    string lowestQValueAsStr;
    {
      stringstream ss;
      ss << scientific;
      ss << lowestQValue;
      ss >> lowestQValueAsStr;
    }
    boost::replace_all(thisPreview, "<q-value>", lowestQValueAsStr);

    string annotationsForHOT=idComponent2Annotations[i].getAnnotationsForHOTForIndexPage(i);

    //add this preview to all previews
    previews.push_back(ObjectPreview(i, lowestQValue, annotationsSQL, thisPreview, annotationsForHOT));
  }


  //output the object previews
  stringstream ssPreview;
  ssPreview << "[";
  for (const auto &preview : previews)
    ssPreview << preview.toJSObject() << ", ";
  ssPreview << "]";

  //create the index file
  //read template file
  string templatePath = DBGWAS_lib + "/index_template.html";
  string indexOutput = readFileAsString(templatePath.c_str());


  //put the command-line in the index page
  string commandLineAsStr;
  {
    stringstream ss;

    class PrintCommandLine : public IPropertiesVisitor {
    private:
        stringstream &ss;
    public:
        PrintCommandLine(stringstream &ss):ss(ss) {}

        /** Called before the true visit of the IProperty instance. */
        virtual void visitBegin    () {}

        /** Visit of the IProperty instance.
         * \param[in] prop : the instance to be visited.
         */
        virtual void visitProperty (IProperty* prop)  {
          ss << prop->key;
          if (prop->value.size())
            ss << " = " << prop->value;
          ss << endl;
        }

        /** Called after the true visit of the IProperty instance. */
        virtual void visitEnd      () {}
    };
    PrintCommandLine printCommandLine(ss);
    getInput()->accept(&printCommandLine);

    commandLineAsStr = ss.str();
  }
  boost::replace_all(indexOutput, "<command_line>", commandLineAsStr);

  //put the info in the template file
  boost::replace_all(indexOutput, "<previews>", ssPreview.str());

  //populate the annotation dropdown filter
  {
    set<string> allTags;
    allTags.insert("No annotations found");
    for (const auto &idComponent2AnnotationsPair : idComponent2Annotations) {
      auto tagsOfThisComponent = idComponent2AnnotationsPair.second.getAnnotationIndexAsSet();
      allTags.insert(tagsOfThisComponent.begin(), tagsOfThisComponent.end());
    }

    stringstream ss;
    ss << "{";
    for (const auto &tag : allTags)
      ss << "'" << UNIQUE_SYMBOL_MARKER << tag << UNIQUE_SYMBOL_MARKER "' : '" << tag << "', ";
    ss << "}";
    boost::replace_all(indexOutput, "<all_tags_in_all_components>", ss.str());
  }

  //put the version on the index page
  boost::replace_all(indexOutput, "<version>", VERSION);

  //put the statistical figures in the output if the -newick parameter was given
  if (hasNewickFile) {
    boost::replace_all(indexOutput, "<stats_images_html>",
        "    p-values of tested unitigs sorted by principal component (each unitig is associated with the closest PC):<br/>\n"
        "    <img class=\"statImage\" src=\"components/stats/bugwas_out_SNPs_PC_manhattan.png\" /><br/><br/><br/>\n"
        "    p-value of each principal component, whose association with the phenotype is tested using a Bayesian Wald test:<br/>\n"
        "    <img class=\"statImage\" src=\"components/stats/bugwas_out_barplot_BayesianWald_PCs.png\" /><br/><br/><br/>\n"
        "    Phylogenetic tree annotated with the principal components which were found significantly associated with the phenotype using a Bayesian Wald test:<br/>\n"
        "    <img class=\"statImage\" src=\"components/stats/bugwas_out_tree_branchescolouredbyPC.png\" />");
    boost::filesystem::create_directories(visualisationsFolder + string("/components/stats/"));
    boost::filesystem::copy_file(step2OutputFolder + string("/bugwas_out_SNPs_PC_manhattan.png"), visualisationsFolder + string("/components/stats/bugwas_out_SNPs_PC_manhattan.png"));
    boost::filesystem::copy_file(step2OutputFolder + string("/bugwas_out_barplot_BayesianWald_PCs.png"), visualisationsFolder + string("/components/stats/bugwas_out_barplot_BayesianWald_PCs.png"));
    boost::filesystem::copy_file(step2OutputFolder + string("/bugwas_out_tree_branchescolouredbyPC.png"), visualisationsFolder + string("/components/stats/bugwas_out_tree_branchescolouredbyPC.png"));
  }else {
    boost::replace_all(indexOutput, "<stats_images_html>", "Re-run DBGWAS with a newick tree file (-newick parameter) to view figures on lineage effect.");
  }

  //output the file
  ofstream indexFile;
  openFileForWriting(visualisationsFolder+string("/index.html"), indexFile);
  indexFile << indexOutput;
  indexFile.close();
  cerr << "[Creating index file...] - Done!" << endl;
}


void generate_output::generateCytoscapeOutput(const graph_t &graph, const vector<MyVertex> &nodes, const string &typeOfGraph, int i,
                             const string &tmpFolder, const string &visualisationsFolder, const vector<int> &selectedUnitigs, int nbPheno0, int nbPheno1,
                             map<int, AnnotationRecord > &idComponent2Annotations,
                             int nbCores) {
  cerr << "Rendering " << typeOfGraph << "_" << i << "..." << endl;


  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //Annotation step
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //some indexes to help gathering info later
  AnnotationRecord annotationsOfThisComponent; //will store all the annotations in this component


  if (thereIsNucleotideDB || thereIsProteinDB) {
    cerr << "Annotating..." << endl;

    //create the input for blast
    string blastInputPath;
    {
      stringstream blastInputPathSS;
      blastInputPathSS << tmpFolder << "/nodes_comp_" << i << ".fasta";
      blastInputPath = blastInputPathSS.str();
    }

    //create a file containing all unitigs in this component
    {
      ofstream blastInputFile;
      openFileForWriting(blastInputPath, blastInputFile);
      for (const auto &node : nodes)
        blastInputFile << ">" << node << endl << graph[node].name << endl;
      blastInputFile.close();
    }

    //execute the blast
    vector <BlastRecord> records;
    if (thereIsNucleotideDB) {
      vector <BlastRecord> blastnRecords = Blast::blast("blastn", blastInputPath, nucleotideDBPath, nbCores);
      records.insert(records.end(), blastnRecords.begin(), blastnRecords.end());
    }
    if (thereIsProteinDB) {
      vector <BlastRecord> blastxRecords = Blast::blast("blastx", blastInputPath, proteinDBPath, nbCores);
      records.insert(records.end(), blastxRecords.begin(), blastxRecords.end());
    }

    //populate annotationsOfThisComponent
    for (const auto &record : records) {
      //here I have to use graph[v].id instead of simply record.nodeId because the original IDs of the node is used later
      MyVertex v = vertex(record.nodeId, graph);
      annotationsOfThisComponent.addAnnotation(record.DBGWAS_tags.at("specific"), graph[v].id, record.evalue, &record);
    }

    //populate idComponent2Annotations[i]
    for (const auto &record : records) {
      //here I have to use graph[v].id instead of simply record.nodeId because the original IDs of the node is used later
      MyVertex v = vertex(record.nodeId, graph);
      idComponent2Annotations[i].addAnnotation(record.DBGWAS_tags.at("general"), graph[v].id, record.evalue);
    }

    cerr << "Annotating... - Done!" << endl;
  }else {
    cerr << "Skipping annotation step - no DB provided" << endl;
  }
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //Annotation step
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //Cytoscape graph build step
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cerr << "Building Cytoscape graph..." << endl;

  //gets the maxCoverage of the nodes in this component - will be used to normalize the width and height of the nodes
  int maxCoverage=-1;
  for (const auto &node : nodes) {
    maxCoverage = max(maxCoverage, graph[node].phenoCounter.getTotal());
  }


  //declares the stringstream which will store the elements to be printed
  stringstream elementsSS;
  {
    elementsSS << scientific; //scientific notation for double

    //goes through the nodes and print them
    set<MyVertex> verticesInThisComponent; //keep track of the vertices in this component
    for (const auto &node : nodes) {
      verticesInThisComponent.insert(node);

      //print the node with the full data
      elementsSS << "{data: {id: 'n" << graph[node].id << "'" <<
      ", name: '" << graph[node].name << "'" <<
      ", sequenceLength: '" << graph[node].name.length() << "'" <<
      ", info: '" << graph[node].id << "'" <<
      ", total: '" << graph[node].phenoCounter.getTotal() << "'" <<
      ", annotations: " << annotationsOfThisComponent.getAllAnnotationsIDsFromANodeAsJSVector(graph[node].id) <<
      ", pheno0: '" << graph[node].phenoCounter.getPheno0() << "/" << nbPheno0 << "'" <<
      ", pheno1: '" << graph[node].phenoCounter.getPheno1() << "/" << nbPheno1 << "'" <<
      ", NA: '" << graph[node].phenoCounter.getNA() << "'" <<
      ", significant: '" << (find(selectedUnitigs.begin(), selectedUnitigs.end(), graph[node].id) == selectedUnitigs.end() ? "No" : "Yes") << "'" <<
      ", qValue: '" << graph[node].unitigStats.getQValueAsStr() << "'" <<
      ", weight: '" << graph[node].unitigStats.getWeightAsStr() << "'" <<
      ", waldStatistic: '" << graph[node].unitigStats.getWaldStatisticAsStr() << "'" <<
      ", background_color: rgbToHex(" << graph[node].unitigStats.getRGB() << ")" <<
      ", width: " << (minSize + (((double) graph[node].phenoCounter.getTotal()) / maxCoverage * (maxSize - minSize))) <<
      ", height: " << (minSize + (((double) graph[node].phenoCounter.getTotal()) / maxCoverage * (maxSize - minSize))) <<
      ", transparency: " <<
      (find(selectedUnitigs.begin(), selectedUnitigs.end(), graph[node].id) == selectedUnitigs.end() ? "76" : "255") <<

      //we print the style before in order to be able to export to Cytoscape Desktop
      "}, style: {'background-color': rgbToHex(" << graph[node].unitigStats.getRGB() << ")" <<
      ", 'width': " << (minSize + (((double) graph[node].phenoCounter.getTotal()) / maxCoverage * (maxSize - minSize))) <<
      ", 'height': " << (minSize + (((double) graph[node].phenoCounter.getTotal()) / maxCoverage * (maxSize - minSize))) <<
      //if it is not a selected unitig, then it becomes a little bit transparent
      (find(selectedUnitigs.begin(), selectedUnitigs.end(), graph[node].id) == selectedUnitigs.end()
       ? ", 'background-opacity': 0.3" : "") <<
      "}}, ";
    }

    //goes through the edges and print them
    for (auto ep = edges(graph); ep.first != ep.second; ++ep.first) {
      MyEdge e = *ep.first;
      //check if the edge is in this component
      if (verticesInThisComponent.find(source(e, graph)) != verticesInThisComponent.end() &&
          verticesInThisComponent.find(target(e, graph)) != verticesInThisComponent.end()) {
        //output edge
        elementsSS << "{data: {id: 'e" << graph[e].id << "', source: 'n" << graph[source(e, graph)].id
        << "', target: 'n" << graph[target(e, graph)].id << "'}}, ";
      }
    }
  }

  //this is what should replace <elementsTag>
  string elements = elementsSS.str();



  //create the graph file
  //read template file
  string templatePath = DBGWAS_lib + "/cytoscape_template.html";
  string cytoscapeOutput = readFileAsString(templatePath.c_str());

  //put the graph in the template file
  boost::replace_all(cytoscapeOutput, "<elementsTag>", elements);

  //put the annotation names in the template file
  boost::replace_all(cytoscapeOutput, "<allAnnotationsTag>", annotationsOfThisComponent.getAnnotationIndexAsJSVector());

  //put the annotation info into the template file
  boost::replace_all(cytoscapeOutput, "<componentAnnotationTag>", annotationsOfThisComponent.getJSRepresentationAnnotIdAnnotInfoGraphPage());

  //put the node2AnnotationEvalue info into the template file
  boost::replace_all(cytoscapeOutput, "<node2AnnotationEvalueTag>", annotationsOfThisComponent.getJSRepresentationNodeId2AnnotationsEvalueForGraphPage());

  //put the annotation2Nodes info into the template file
  boost::replace_all(cytoscapeOutput, "<annotation2NodesParTag>", annotationsOfThisComponent.getJSRepresentationAnnotation2NodesForGraphPage());

  //put the extraTags info into the template file
  boost::replace_all(cytoscapeOutput, "<extraTagsPar>", annotationsOfThisComponent.getExtraTagsAsJSVector());

  //output the file
  string outfilename;
  {
    stringstream ss;
    ss << visualisationsFolder << "/components/" << typeOfGraph << "_" << i << ".html";
    outfilename=ss.str();
  }
  ofstream outFile;
  openFileForWriting(outfilename, outFile);
  outFile << cytoscapeOutput;
  outFile.close();


  //copy the lib folder, if it is not already copied
  string fromLibPath = DBGWAS_lib + "/lib";
  string toLibPath = visualisationsFolder + "/components/lib";
  if (!boost::filesystem::exists(toLibPath))
    copyDirectoryRecursively(fromLibPath, toLibPath);
  cerr << "Building Cytoscape graph... - Done!" << endl;
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //Cytoscape graph build step
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  cerr << "Rendering " << typeOfGraph << "_" << i << "... - Done!" << endl;
}

void generate_output::execute () {
  checkParametersGenerateOutput(this);
  int neighbourhood = getInput()->getInt(STR_MAX_NEIGHBOURHOOD);
  string outputFolderParameter = stripLastSlashIfExists(getInput()->getStr(STR_OUTPUT));
  string outputFolder = outputFolderParameter+string("/step3");
  string tmpFolder = outputFolder+string("/tmp");
  string visualisationsFolder = outputFolderParameter+string("/visualisations");
  string step1OutputFolder = outputFolderParameter+string("/step1");
  string step2OutputFolder = outputFolderParameter+string("/step2");
  int nbCores = getInput()->getInt(STR_NBCORES);


  //get the nbContigs
  int nbContigs = getNbLinesInFile(step1OutputFolder+string("/graph.nodes"));

  cerr << "[Getting the significant unitigs from the patterns...]" << endl;

  //read the patterns
  auto sortedPatterns = PatternFromStats::readFile(step2OutputFolder + "/patterns.txt", true);

  //associates each unitig to its pattern stats
  vector<const PatternFromStats*> unitigToPatternStats(nbContigs, NULL);
  {
    //since not all patterns have stats, this maps the id of the pattern to its stats, so we know easily which pattern has stats
    map<int, const PatternFromStats *> patternsToStats;
    for (const auto &pattern : sortedPatterns)
      patternsToStats[pattern.pattern] = &pattern;

    //associates each unitig to its Pattern stats
    ifstream unitigToPatternFilestream;
    openFileForReading(step1OutputFolder + string("/gemma_input.unitig_to_pattern.binary"), unitigToPatternFilestream);

    int unitig, pattern;
    while (unitigToPatternFilestream >> unitig >> pattern) {
      if (patternsToStats.find(pattern) != patternsToStats.end()) {
        //the pattern has a stats
        unitigToPatternStats[unitig] = patternsToStats[pattern];
      }
    }
    unitigToPatternFilestream.close();
  }

  //read unitigToWeightCorrection
  vector<int> unitigToWeightCorrection(nbContigs);
  ifstream unitigToWeightCorrectionFilestream;
  openFileForReading(step1OutputFolder+string("/weight_correction"), unitigToWeightCorrectionFilestream);
  for (int i=0;i<unitigToWeightCorrection.size();i++)
    unitigToWeightCorrectionFilestream >> unitigToWeightCorrection[i];
  unitigToWeightCorrectionFilestream.close();

  //get the significant patterns according to SFF
  vector<PatternFromStats> significantPatterns;
  {
    boost::apply_visitor(GetSignificantPatterns(sortedPatterns, significantPatterns), SFF);
    cerr << "Selected " << significantPatterns.size() << " most significant patterns..." << endl;
  }

  //transform the significant patterns into significant unitigs
  {
    ofstream significantUnitigsFile;
    openFileForWriting(outputFolder + string("/significant_unitigs.txt"), significantUnitigsFile);
    ifstream pattern2UnitigsFile;
    openFileForReading(step1OutputFolder + string("/bugwas_input.unique_rows_to_all_rows.binary"), pattern2UnitigsFile);

    //read pattern2UnitigsFile
    string line;
    vector<vector<int> > pattern2Unitigs;
    while (getline(pattern2UnitigsFile, line)) {
      stringstream ss(line);
      int unitig;
      vector<int> unitigs;
      while (ss >> unitig)
        unitigs.push_back(unitig);
      if (unitigs.size() > 0)
        pattern2Unitigs.push_back(unitigs);
    }

    //process significantPatternsFile and produces significantUnitigsFile
    for (const auto &significantPattern : significantPatterns) {
      for (const auto &unitig : pattern2Unitigs[significantPattern.pattern]) {
        significantUnitigsFile << unitig << endl;
      }
    }

    pattern2UnitigsFile.close();
    significantUnitigsFile.close();
  }
  cerr << "[Getting significant unitigs from the patterns...] - Done!" << endl;

  cerr << "[Reading input and creating BOOST graph...]" << endl;

  //declare the graph
  graph_t graph(nbContigs);

  //get the unitigs selected for visualization
  vector<int> selectedUnitigs;
  {
    ifstream selectedUnitigsStream;
    openFileForReading(outputFolder + string("/significant_unitigs.txt"), selectedUnitigsStream);
    int selectedUnitig;
    while (selectedUnitigsStream >> selectedUnitig)
      selectedUnitigs.push_back(selectedUnitig);
    selectedUnitigsStream.close();
  }

  //read the frequency file
  vector< PhenoCounter > phenotypeCounters;
  {
    string frequencyFilename(step1OutputFolder+string("/frequency_unitig_to_total_pheno0_pheno1_NA_count"));
    ifstream frequencyFile;
    openFileForReading(frequencyFilename, frequencyFile);

    int temp1, temp2, temp3, temp4;
    while (frequencyFile >> temp1) {
      frequencyFile >> temp2 >> temp3 >> temp4;
      phenotypeCounters.push_back(PhenoCounter(temp2, temp3, temp4));
    }

    frequencyFile.close();
  }

  //read the total number of pheno0 and pheno1 strains
  int nbPheno0, nbPheno1;
  {
    ifstream totalNbOfStrainsInEachPheno;
    openFileForReading(step1OutputFolder + string("/total_nb_of_strains_in_each_pheno"), totalNbOfStrainsInEachPheno);
    totalNbOfStrainsInEachPheno >> nbPheno0 >> nbPheno1;
    totalNbOfStrainsInEachPheno.close();
  }




  //create the nodes of the graph
  //assume the graph is undirected, instead of directed
  string nodesFile = step1OutputFolder+string("/graph.nodes");
  {
    ifstream nodesFileReader;
    openFileForReading(nodesFile, nodesFileReader);
    int id;
    string seq;
    int index = 0;
    while (nodesFileReader >> id >> seq) {
      MyVertex vF = vertex(index, graph);
      //TODO: sequence should be stored as a 2-bit array
      graph[vF].name = seq;
      graph[vF].id = id;
      graph[vF].strand = 'F';
      graph[vF].phenoCounter = phenotypeCounters[id];
      graph[vF].unitigStats = UnitigStats(unitigToPatternStats[id], unitigToWeightCorrection[id]);
      index++;
    }
    nodesFileReader.close();
  }

  //create the edges of the graph
  string edgesFile = step1OutputFolder+string("/graph.edges.dbg");
  {
    ifstream edgesFileReader;
    openFileForReading(edgesFile, edgesFileReader);
    int from, to;
    string label;
    int index = 0;
    while (edgesFileReader >> from >> to >> label) {
        pair<MyEdge, bool> return_from_add_edge = add_edge(vertex(from, graph),
                                                           vertex(to, graph),
                                                           graph);
        if (return_from_add_edge.second) {
          graph[return_from_add_edge.first].id = index;
          graph[return_from_add_edge.first].weight = 1;
          index++;
        }
    }
    edgesFileReader.close();
  }
  cerr << "[Reading input and creating BOOST graph...] - Done!" << endl;


  cerr << "[Computing nodes' neighbourhoods...]" << endl;
  //compute the neighbourhood
  vector<int> distances(num_vertices(graph));
  set<MyVertex> verticesInTheNeighbourhood;
  map <int, set<MyVertex> > verticesInTheNeighbourhoodOfThisUnitig;

  for (auto unitig : selectedUnitigs) {
    MyVertex selectedVertex = vertex(unitig, graph);
    try {
      fill(distances.begin(), distances.end(), 0);
      dijkstra_shortest_paths(graph, selectedVertex, weight_map(get(&EdgeInfo::weight, graph)).
          distance_map(make_iterator_property_map(distances.begin(),
                                                  boost::get(boost::vertex_index, graph))).
          visitor(StopWhenVeryDistantFromSourceDijkstraVisitor(distances, neighbourhood,
                                                               verticesInTheNeighbourhoodOfThisUnitig[unitig])));
    } catch (TooDistant &e) { }

    //update verticesInTheNeighbourhood
    verticesInTheNeighbourhood.insert(verticesInTheNeighbourhoodOfThisUnitig[unitig].begin(), verticesInTheNeighbourhoodOfThisUnitig[unitig].end());
  }
  cerr << "[Computing nodes' neighbourhoods...] - Done!" << endl;



  cerr << "[Generating the visualisation files...]" << endl;
  //print one graph per component
  //create a subgraph containing only the nodes in the neighbourhoods
  graph_t& newGraph = graph.create_subgraph();
  vector<vector<MyVertex> > nodesInComponent; //care: the nodes in this variable are the nodes in newGraph, not in the normal graph
  map<int, AnnotationRecord > idComponent2Annotations;
  int numberOfComponents=0;
  {
    for (auto vp = vertices(graph); vp.first != vp.second; ++vp.first) {
      MyVertex v = *vp.first;
      if (verticesInTheNeighbourhood.find(v)!=verticesInTheNeighbourhood.end()) {
        add_vertex(v, newGraph);
      }
    }

    //get the components of this subgraphs
    vector<int> componentOfThisNode = vector<int>(num_vertices(newGraph));
    int num = connected_components(newGraph, &componentOfThisNode[0]);
    nodesInComponent = vector<vector<MyVertex> >(num);
    for (int i = 0; i != componentOfThisNode.size(); ++i)
      nodesInComponent[componentOfThisNode[i]].push_back(vertex(i, newGraph));

    for (int i = 0; i < nodesInComponent.size(); i++) {
      generateCytoscapeOutput(newGraph, nodesInComponent[i], "comp", i, tmpFolder, visualisationsFolder, selectedUnitigs, nbPheno0, nbPheno1,
                              idComponent2Annotations, nbCores);
    }
    numberOfComponents = nodesInComponent.size();
  }

  cerr << "[Generating the visualisation files...] - Done!" << endl;


  //create the index
  createIndexFile(numberOfComponents, visualisationsFolder, step2OutputFolder, nodesInComponent, newGraph,
                  idComponent2Annotations, unitigToPatternStats, selectedUnitigs, nbCores);

  //clean-up - saving some disk space
  //remove temp directory
  boost::filesystem::remove_all(tmpFolder);

  //tell we are done
  cout << endl << endl <<
      "******************************************************************************" << endl <<
      "We are done. The output can be found at " << visualisationsFolder << "/index.html" << endl <<
      "******************************************************************************" << endl << endl;
  cout.flush();
}