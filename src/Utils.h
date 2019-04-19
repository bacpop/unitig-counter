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

#ifndef KISSPLICE_UTILS_H
#define KISSPLICE_UTILS_H

#include <gatb/gatb_core.hpp>
#include <string>
#include <vector>
#include <fstream>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <limits>
#include <iterator>
#include <set>
#include <boost/algorithm/string/replace.hpp>
#define BOOST_NO_CXX11_SCOPED_ENUMS
#include <boost/filesystem.hpp>
#undef BOOST_NO_CXX11_SCOPED_ENUMS
#include <pstream.h>


using namespace std;
namespace fs = boost::filesystem;

char complement(char b);
string reverse_complement(const string &seq);
//Read all strings in the readsFile file and return them as a vector of strings
vector<string> getVectorStringFromFile(const string &readsFile);

//this function also populates strains if needed
void checkStrainsFile(const string &strainsFile);

string readFileAsString(const char* fileName);

//strips all last "/" if exists in the parameter
string stripLastSlashIfExists (string path);

void copyDirectoryRecursively(const fs::path& sourceDir, const fs::path& destinationDir);

int getNbLinesInFile(const string &filename);

void checkParametersBuildDBG(Tool *tool);
void fatalError (const string &message);
void executeCommand(const string &command, bool verbose=true, const string &messageIfItFails="");
void openFileForReading(const string &filePath, ifstream &stream);
void openFileForWriting(const string &filePath, ofstream &stream);
void createFolder(const string &path);



//global vars used by all programs
class UnitigIdStrandPos {
public:
    int unitigId;
    char strand; //TODO: change for bool
    int pos;
    int unitigSize;
    int kmerSize; //the size of the kmer (just for the reverseStrand() function) //TODO: static global var, sth like that

    UnitigIdStrandPos(int unitigId=0, char strand='?',int pos=0, int unitigSize=0, int kmerSize=0):
        unitigId(unitigId), strand(strand), pos(pos), unitigSize(unitigSize), kmerSize(kmerSize){}
    int getUnitigId () const {checkValidity(); return unitigId; }
    char getStrand () const {checkValidity(); return strand; }
    int getPos() const {checkValidity(); return pos; }
    string toString() const {
      stringstream ss;
      try {
        ss << getUnitigId() << " " << getStrand() << " " << getPos();
      }catch (const runtime_error &e) {
        ss << "-1 ? -1";
      }
      return ss.str();
    }
    void reverseStrand () {
      strand = (strand=='F' ? 'R' : 'F');
      pos = unitigSize-pos-kmerSize;
    }

    //check if the unitig is valid or not
    void checkValidity() const {
      if (strand=='?')
        throw runtime_error("Invalid unitig.");
    }
};

struct Strain {
    string id, path;
    Strain(const string &id, const string &path) : id(id) {
      //transfor to canonical path
      boost::filesystem::path boostPath(boost::filesystem::canonical(path));
      this->path = boostPath.string();
    }

    static void createReadsFile(const string &readsFile, vector< Strain >* strains) {
      ofstream fout;
      openFileForWriting(readsFile, fout);

      for (const auto &strain : (*strains))
        fout << strain.path << endl;

      fout.close();
    }

};

struct PatternFromStats {
    int pattern;
    long double qValue;
    long double weight;
    long double normalizedWeight;
    string waldStatistic;

    bool operator < (const PatternFromStats& other) const {
      return this->qValue < other.qValue;
    }

    static vector<PatternFromStats> readFile(const string &filename, bool header=false) {
      vector<PatternFromStats> patterns;

      //read
      {
        ifstream patternStream;
        openFileForReading(filename, patternStream);
        patternStream >> setprecision(std::numeric_limits<long double>::digits10 + 1);
        PatternFromStats pattern;

        //remove header
        if (header) {
          string tmp;
          getline(patternStream, tmp);
        }
        while (patternStream >> pattern.pattern >> pattern.qValue >> pattern.weight >> pattern.waldStatistic)
          patterns.push_back(pattern);
        patternStream.close();
      }

      //normalize
      vector<long double>  normalizedWeights;
      for (const auto &pattern : patterns)
        normalizedWeights.push_back(pattern.weight);

      //remove the minimum weight of everyone
      long double minWeight = *min_element(normalizedWeights.begin(), normalizedWeights.end());
      transform(normalizedWeights.begin(), normalizedWeights.end(), normalizedWeights.begin(), [&](long double weight) {
          return weight-minWeight;
      });

      //transform to [0,1]
      long double maxWeight = *max_element(normalizedWeights.begin(), normalizedWeights.end());
      transform(normalizedWeights.begin(), normalizedWeights.end(), normalizedWeights.begin(), [&](long double weight) {
          return weight/maxWeight;
      });

      //write back
      auto normalizedWeightsIt = normalizedWeights.begin();
      for (auto &pattern : patterns) {
        pattern.normalizedWeight = *normalizedWeightsIt;
        ++normalizedWeightsIt;
      }

      return patterns;
    }

    static void writeFile(const string &filename, const vector<PatternFromStats> &patterns) {
      ofstream patternSortedStream;
      openFileForWriting(filename, patternSortedStream);
      patternSortedStream << setprecision(std::numeric_limits<long double>::digits10 + 1);
      patternSortedStream << "pattern q-value weight wald_statistic" << endl;
      for (const auto &pattern : patterns)
        patternSortedStream << pattern.pattern << " " << pattern.qValue << " " << pattern.weight << " " << pattern.waldStatistic << endl;
      patternSortedStream.close();
    }
};


#endif //KISSPLICE_UTILS_H
