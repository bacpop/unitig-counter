//
// Created by Leandro Ishi Soares de Lima on 22/09/17.
//

#ifndef DBGWAS_BLAST_H
#define DBGWAS_BLAST_H

#include <sstream>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <boost/algorithm/string.hpp>
#include <boost/regex.hpp>

using namespace std;


struct ValueNotFound{};

class BlastRecord {
public:
    int nodeId; //the query id
    map<string, string> DBGWAS_tags;
    double qcovs, bitscore, pident; //the blast fields we care about
    long double evalue; //the blast fields we care about

    //constructors
    BlastRecord(){}
    BlastRecord(int nodeId, const map<string, string> &DBGWAS_tags,
                double qcovs, double bitscore, double pident, long double evalue):
        nodeId(nodeId), DBGWAS_tags(DBGWAS_tags), qcovs(qcovs),
        bitscore(bitscore), pident(pident), evalue(evalue) {}

    //parse a string and build a BlastRecord from it
    static BlastRecord parseString (const string &str);
private:

    static string extractValue (const string &header, const string &tag);

    //parse the header and extract all the DBGWAS tags
    //header is intentionally string and not const string &
    static map<string, string> extractValuesWithRegex(string header, const boost::regex &expression);
};


//For each set of annotations of a component, record the annotations' name, the nodes mapping to each and the lowest e-value of a hit to it
//This represents the annotations of a component
class AnnotationRecord {
private:
    vector<string> annotationIndex; //contains the annotations and their indexes

    //stores all the information of an annotation that is to be shown in the upper table in the graph page
    class AnnotationInfoGraphPage {
    private:
        set<int> nodes;
        long double minEvalue;
        BlastRecord record;

    public:
        AnnotationInfoGraphPage():nodes(), minEvalue(std::numeric_limits<long double>::max()), record(){}

        //add a new node to this set
        void addNode(int node, long double evalue, const BlastRecord* record);

        //get the nb of nodes mapping here
        int getNbOfNodes () const { return nodes.size(); }

        //get the evalue
        long double getMinEvalue () const { return minEvalue; }

        //transform to a javascript array
        string getHTMLRepresentationForGraphPage (const set<string>& allExtraTags);

        //transform to a javascript array
        string getHTMLRepresentationForIndexPage () const;

        //get the nodes as a javascript array
        string getNodesAsJSArray() const;
    };
    //maps an annotation (int) to its informations to be shown in the graph page
    map<int, AnnotationInfoGraphPage> annotations;
    set<string> allExtraTags; //records all extra tags in this component


    class AnnotationsAndEvalue {
    public:
        map<int, long double> annotation2Evalue;
        AnnotationsAndEvalue(){}

        //add annotation
        void addAnnotation(int annotation, long double evalue);
    };

    //maps a node (int) to the set of annotations and evalues it maps to
    //this is to be shown in the node handsontable in node table
    map<int, AnnotationsAndEvalue > nodeId2Annotations;
public:
    AnnotationRecord():annotations(){}

    //add an annotation to this set
    void addAnnotation(const string &tag, int node, long double evalue, const BlastRecord* record=NULL);

    //get a representation of this annotation to be added to the SQL string in the index page
    string getSQLRepresentationForIndexPage() const;

    //get JS array representation of the annotation component for the index page
    string getAnnotationsForHOTForIndexPage(int componentId) const;

    //get an HTML representation of the annotation component for the graph page, with the annotation index, and all other info like nb of nodes, evalue and extra tags
    string getJSRepresentationAnnotIdAnnotInfoGraphPage();

    //get a JS array mapping annotation IDs to the nodes it maps to
    string getJSRepresentationAnnotation2NodesForGraphPage() const;

    //gets the annotation index as a JS vector
    string getAnnotationIndexAsJSVector() const;

    //get all the annotations names as a set
    set<string> getAnnotationIndexAsSet() const;

    //get all the annotations IDs from a node as JS vector
    string getAllAnnotationsIDsFromANodeAsJSVector(int node);

    //get the extra tags as a JS vector
    string getExtraTagsAsJSVector() const;

    //get a dictionary in JS where the key is the node id and the value is a pair annotation and evalue
    string getJSRepresentationNodeId2AnnotationsEvalueForGraphPage() const;
};



class Blast {
public:
    //Blast (using command) the file in queryPath agains the db on dbPath and return the results as a vector<BlastRecord>
    static vector<BlastRecord> blast (const string &command, const string &queryPath, const string &dbPath, int nbCores);

    //make the blast DB
    //returns the path to the blast DB
    static string makeblastdb (const string &dbtype, const string &originalDBPath, const string &outputFolderPath);
};



#endif //DBGWAS_BLAST_H
