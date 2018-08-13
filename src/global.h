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

#ifndef KSGATB_GLOBAL_H
#define KSGATB_GLOBAL_H
#include <gatb/gatb_core.hpp>
#include "Utils.h"

#define UNIQUE_SYMBOL_MARKER "#@#-322"

//global vars
extern Graph* graph;
extern vector< UnitigIdStrandPos >* nodeIdToUnitigId;
extern vector< Strain >* strains;
extern const char* STR_STRAINS_FILE;
extern const char* STR_KSKMER_SIZE;
extern const char* STR_OUTPUT;
extern const char* STR_NBCORES;
extern const char* STR_MAX_NEIGHBOURHOOD;
extern const char* STR_SKIP1;
extern const char* STR_SKIP2;
extern const char* STR_NEWICK_PATH;
extern const char* STR_SFF;
extern const char* STR_NUCLEOTIDE_DB;
extern const char* STR_PROTEIN_DB;
extern const char* STR_MAF_FILTER;
extern const char* STR_GEMMA_PATH;
extern const char* STR_BLAST_PATH;
extern const char* STR_PHANTOMJS_PATH;
extern const char* STR_RSCRIPT_PATH;
extern const char* STR_NO_PREVIEW;

//TODO: seeveral questions are still unclear if we use the Freq count mode (how to run bugwas, the coloring, etc...). For now I am disabling this option
//extern const char* STR_COUNT_MODE;
//TODO: seeveral questions are still unclear if we use the Freq count mode (how to run bugwas, the coloring, etc...). For now I am disabling this option
extern string DBGWAS_lib;
extern bool skip1;
extern bool skip2;
extern bool hasNewickFile;
extern bool presenceAbsenceCountMode;
extern bool thereIsNucleotideDB;
extern string nucleotideDBPath;
extern bool thereIsProteinDB;
extern string proteinDBPath;
extern boost::variant< int, double > SFF;
extern string gemmaPath;
extern string blastPath;
extern string phantomjsPath;
extern string RscriptPath;
extern bool noPreview;

void populateParser (Tool *tool);

#endif //KSGATB_GLOBAL_H
