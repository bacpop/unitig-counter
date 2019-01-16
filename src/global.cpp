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
#include <string>

using namespace std;

const char* STR_STRAINS_FILE = "-strains";
const char* STR_MAX_NEIGHBOURHOOD = "-nh";
const char* STR_KSKMER_SIZE = "-k";
const char* STR_OUTPUT = "-output";
const char* STR_NBCORES = "-nb-cores";
const char* STR_RUN1 = "-run1";
const char* STR_RUN2 = "-run2";
const char* STR_SFF = "-SFF";
const char* STR_MAF_FILTER = "-maf";


//variables controlling where the executable is
string DBGWAS_lib = ""; //set on first command in main.cpp


bool run1 = true;
bool run2 = false;
bool presenceAbsenceCountMode = false;
boost::variant< int, double > SFF;

//global vars used by both programs
Graph* graph;
vector< UnitigIdStrandPos >* nodeIdToUnitigId;
vector< Strain >* strains = NULL;

void populateParser (Tool *tool) {
  // We add some custom arguments for command line interface
  tool->getParser()->push_front (new OptionOneParam (STR_MAF_FILTER, "Minor Allele Frequency Filter.",  false, "0.01"));
  tool->getParser()->push_front (new OptionOneParam (STR_MAX_NEIGHBOURHOOD, "Denotes the neighbourhood to be considered around the significant unitigs.",  false, "5"));
  tool->getParser()->push_front (new OptionOneParam (STR_SFF, "Denotes the Significant Features Filter - the features (or patterns) selected to create a visualisation around them. If it is a float number n, then only the features with q-value<=n are selected. If it is an integer n, then only the n first features are selected. Take a look at the output/step2/patterns.txt file to get a list of features ordered by q-value to better choose this parameter (re-run the tool with -skip2 in order to directly produce the visualisation of the features selected by your parameter).",  false, "100"));
  tool->getParser()->push_front (new OptionNoParam (STR_RUN2, "Run in mode 2: Find unitig neighbours",  false));
  tool->getParser()->push_front (new OptionNoParam (STR_RUN1, "Run in mode 1: Count unitigs",  false));
  tool->getParser()->push_front (new OptionOneParam (STR_OUTPUT, "Path to the folder where the final and temporary files will be stored.",  false, "output"));
  tool->getParser()->push_front (new OptionOneParam (STR_KSKMER_SIZE, "K-mer size.",  false, "31"));
  tool->getParser()->push_front (new OptionOneParam (STR_STRAINS_FILE, "A text file describing the strains containing 3 columns: 1) ID of the strain; 2) Phenotype (0/1/NA); 3) Path to a multi-fasta file containing the sequences of the strain. This file needs a header. Check the sample_example folder or https://gitlab.com/leoisl/dbgwas/raw/master/sample_example/strains for an example.",  true));

  //TODO: seeveral questions are still unclear if we use the Freq count mode (how to run bugwas, the coloring, etc...). For now I am disabling this option
  //tool->getParser()->push_front (new OptionOneParam (STR_COUNT_MODE, "The count mode. If \"01\", then the count mode is seen as presence/absence. If \"Freq\", then the count mode is seen as the frequency",  false, "01"));
  //TODO: seeveral questions are still unclear if we use the Freq count mode (how to run bugwas, the coloring, etc...). For now I am disabling this option
}
