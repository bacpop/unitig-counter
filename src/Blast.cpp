//
// Created by Leandro Ishi Soares de Lima on 22/09/17.
//

#include "Blast.h"
#include "Utils.h"
#include "global.h"
#include <algorithm>
#include <boost/algorithm/string.hpp>

//parse a string and build a BlastRecord from it
BlastRecord BlastRecord::parseString (const string &str) {
  BlastRecord record;
  stringstream stream;
  stream << str;

  //read
  string header;
  stream >> record.nodeId >> header >> record.qcovs >> record.bitscore >> record.pident >> record.evalue;

  //escape '
  boost::replace_all(header, "'", "\\'");

  //parse all tags
  boost::regex expression("DBGWAS_(\\w+)_tag_*=_*([^;]+)_*;?");
  map<string, string> allTags = extractValuesWithRegex(header, expression);

  //check if general tag was specified
  if (allTags.count("general")==0)
    //add the header as DBGWAS_general_tag
    allTags["general"] = header;

  //check if general tag is not empty
  if (allTags["general"].size()==0) {
    cerr << "[WARNING] DBGWAS_general_tag of " << header << " is empty! Setting to <EMPTY>" << endl;
    allTags["general"] = "<EMPTY>";
  }

  //check if specific tag was specified
  if (allTags.count("specific")==0)
    //add the header as DBGWAS_specific_tag
    allTags["specific"] = header;

  //check if graph tag is not empty
  if (allTags["specific"].size()==0) {
    cerr << "[WARNING] DBGWAS_specific_tag of " << header << " is empty! Setting to <EMPTY>" << endl;
    allTags["specific"] = "<EMPTY>";
  }

  record.DBGWAS_tags=allTags;

  return record;
}

//parse the header and extract all the DBGWAS tags
//header is intentionally string and not const string &
map<string, string> BlastRecord::extractValuesWithRegex(string header, const boost::regex &expression) {
  map<string, string> extractedValues;
  boost::smatch matchResults;

  while(boost::regex_search(header, matchResults, expression))
  {
    string key = matchResults.str(1);
    boost::trim_if(key, [](char c) -> bool { return c=='_';});
    string value = matchResults.str(2);
    boost::trim_if(value, [](char c) -> bool { return c=='_';});
    extractedValues[key]=value;

    //go to the next field
    header = matchResults.suffix();
  }
  return extractedValues;
}


vector<BlastRecord> Blast::blast (const string &command, const string &queryPath, const string &dbPath, int nbCores) {
  string outFilePath = queryPath+"."+command+"Out";

  //build the command line
  stringstream ss;
  ss << blastPath << "/" << command << " -query " << queryPath << " -db " << dbPath << " -out " << outFilePath << " -num_threads " << nbCores << " -outfmt '6 qseqid sseqid qcovs bitscore pident evalue'";
  string commandLine=ss.str();
  executeCommand(commandLine, false);

  //read output
  vector<BlastRecord> records;
  {
    auto recordsAsStringVector = getVectorStringFromFile(outFilePath);
    for (const auto& recordString : recordsAsStringVector) {
      records.push_back(BlastRecord::parseString(recordString));
    }
  }

  return records;
}


string Blast::makeblastdb (const string &dbtype, const string &originalDBPath, const string &outputFolderPath) {
  string concatenatedDBPath;
  {
    stringstream ss;
    ss << outputFolderPath << "/" << dbtype << "_db";
    concatenatedDBPath = ss.str();
  }

  //concatenate all DBs into one
  {
    //TODO: I don't know why <(echo) is not allowed in the cat command, so we have to use this ugly thing...
    //create an empty file with only a newline
    string newlineFilepath = outputFolderPath + "/newline";
    executeCommand(string("echo > ") + newlineFilepath, false);
    //TODO: I don't know why <(echo) is not allowed in the cat command, so we have to use this ugly thing...

    string catCommand;
    stringstream ss;
    ss << "cat " << originalDBPath << " > " << concatenatedDBPath;
    catCommand = ss.str();
    boost::replace_all(catCommand, ",", string(" ") + newlineFilepath + " ");
    executeCommand(catCommand, false);
  }

  //replace spaces for underscores in the concatenated files, as this could create some problems...
  string fixedDBPath(concatenatedDBPath);
  fixedDBPath += "_fixed";
  {
    string commandLineFixSpaces = string("tr ' ' '_' <") + concatenatedDBPath + " >" + fixedDBPath;
    executeCommand(commandLineFixSpaces, false);
  }

  //create the DB using the fixed FASTA
  {
    string commandLineMakeblastdb = blastPath + "/makeblastdb -dbtype " + dbtype + " -in " + fixedDBPath;
    executeCommand(commandLineMakeblastdb);
  }

  //return the fixed db path
  return fixedDBPath;
}


void AnnotationRecord::AnnotationInfoGraphPage::addNode(int node, long double evalue, const BlastRecord* record) {
  nodes.insert(node);
  minEvalue = min(minEvalue, evalue);
  if (record)
    this->record=*record;
}

//get a representation of this annotation to be added to the SQL string in the index page
string AnnotationRecord::getSQLRepresentationForIndexPage() const {
  stringstream ss;
  if (annotations.size()==0)
    ss << UNIQUE_SYMBOL_MARKER << "No annotations found" << UNIQUE_SYMBOL_MARKER << " ";
  else {
    for (const auto & indexAndAnnotationInfoGraphPage : annotations)
      ss << UNIQUE_SYMBOL_MARKER << annotationIndex[indexAndAnnotationInfoGraphPage.first] << UNIQUE_SYMBOL_MARKER << " ";
  }
  return ss.str();
}

//get JS array representation of the annotation component for the index page
string AnnotationRecord::getAnnotationsForHOTForIndexPage(int componentId) const {
  stringstream ss;
  ss << "[";
  for (const auto & indexAndAnnotationInfoGraphPage : annotations)
    ss << "['" << annotationIndex[indexAndAnnotationInfoGraphPage.first] << "', " << indexAndAnnotationInfoGraphPage.second.getHTMLRepresentationForIndexPage() << "], ";
  ss << "]";

  return ss.str();
}


//transform to a javascript array
string AnnotationRecord::AnnotationInfoGraphPage::getHTMLRepresentationForGraphPage (const set<string>& allExtraTags) {
  stringstream ss;
  ss << scientific;
  ss << nodes.size() << ", " << minEvalue << ", ";
  for (const auto &extraTag : allExtraTags)
    ss << "'" << record.DBGWAS_tags[extraTag] << "', ";
  return ss.str();
}

//transform to a javascript array
string AnnotationRecord::AnnotationInfoGraphPage::getHTMLRepresentationForIndexPage () const {
  stringstream ss;
  ss << scientific;
  ss << nodes.size() << ", " << minEvalue;
  return ss.str();
}

set<string> AnnotationRecord::getAnnotationIndexAsSet() const {
  set<string> allAnnotationsNames;
  for (const auto & annotation : annotationIndex)
    allAnnotationsNames.insert(annotation);
  return allAnnotationsNames;
}


//get an HTML representation of the annotation component for the graph page, with the annotation index, and all other info like nb of nodes, evalue and extra tags
string AnnotationRecord::getJSRepresentationAnnotIdAnnotInfoGraphPage() {
  stringstream ss;
  ss << "[";
  for (auto & indexAndAnnotationInfoGraphPage : annotations)
    ss << "[" << indexAndAnnotationInfoGraphPage.first << ", " << indexAndAnnotationInfoGraphPage.second.getHTMLRepresentationForGraphPage(allExtraTags) << "], ";
  ss << "]";
  return ss.str();
}

//get the nodes as a javascript array
string AnnotationRecord::AnnotationInfoGraphPage::getNodesAsJSArray() const {
  stringstream ss;
  ss << "[";
  for (int node : nodes)
    ss << "'n" << node << "', ";
  ss << "]";
  return ss.str();
}

//get a JS array mapping annotation IDs to the nodes it maps to
string AnnotationRecord::getJSRepresentationAnnotation2NodesForGraphPage() const {
  stringstream ss;
  ss << "[";
  for (auto & indexAndAnnotationInfoGraphPage : annotations)
    ss << indexAndAnnotationInfoGraphPage.second.getNodesAsJSArray() << ", ";
  ss << "]";
  return ss.str();
}

//get all annotations names as a JS vector
string AnnotationRecord::getAnnotationIndexAsJSVector() const {
  stringstream ss;
  ss << "[";
  for (const auto & annotation : annotationIndex)
    ss << "'" << annotation << "', ";
  ss << "]";
  return ss.str();
}


//add an annotation to this set
void AnnotationRecord::addAnnotation(const string &tag, int node, long double evalue, const BlastRecord* record) {
  if (find(annotationIndex.begin(), annotationIndex.end(), tag) == annotationIndex.end())
    //did not find the tag, push it
    annotationIndex.push_back(tag);

  //find the index
  auto index = find(annotationIndex.begin(), annotationIndex.end(), tag)-annotationIndex.begin();

  //add to annotations
  annotations[index].addNode(node, evalue, record);
  nodeId2Annotations[node].addAnnotation(index, evalue);

  //add the extra tags if the record was given
  if(record) {
    for (const auto &pair : record->DBGWAS_tags) {
      if (pair.first!="general" && pair.first!="specific")
        allExtraTags.insert(pair.first);
    }
  }
}

void AnnotationRecord::AnnotationsAndEvalue::addAnnotation(int annotation, long double evalue) {
  if (annotation2Evalue.count(annotation)>0) //it is already present in the map
    annotation2Evalue[annotation] = min(annotation2Evalue[annotation], evalue);
  else
    annotation2Evalue[annotation] = evalue;
}


//get all the annotations IDs from a node as JS vector
string AnnotationRecord::getAllAnnotationsIDsFromANodeAsJSVector(int node) {
  stringstream ss;
  ss << "[";
  for (const auto &annotationAndEvalue : nodeId2Annotations[node].annotation2Evalue)
    ss << annotationAndEvalue.first << ", ";
  ss << "]";
  return ss.str();
}

//get a dictionary in JS where the key is the node id and the value is a pair annotation and evalue
string AnnotationRecord::getJSRepresentationNodeId2AnnotationsEvalueForGraphPage() const {
  stringstream ss;
  ss << scientific;

  ss << "{";

  for (const auto &nodeId2AnnotationsIt : nodeId2Annotations) {
    ss << "'n" << nodeId2AnnotationsIt.first << "': [";

    for (const auto &annotationEvalue : nodeId2AnnotationsIt.second.annotation2Evalue)
      ss << "[" << annotationEvalue.first << ", " << annotationEvalue.second << "], ";

    ss << "], ";
  }

  ss << "}";

  return ss.str();
};


//get the extra tags as a JS vector
string AnnotationRecord::getExtraTagsAsJSVector() const {
  stringstream ss;
  ss << "[";
  for (const string &extraTag : allExtraTags)
    ss << "'" << extraTag << "', ";
  ss << "]";
  return ss.str();
}