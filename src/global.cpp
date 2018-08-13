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
const char* STR_SKIP1 = "-skip1";
const char* STR_SKIP2 = "-skip2";
const char* STR_NEWICK_PATH = "-newick";
const char* STR_SFF = "-SFF";
const char* STR_NUCLEOTIDE_DB = "-nc-db";
const char* STR_PROTEIN_DB = "-pt-db";
const char* STR_MAF_FILTER = "-maf";
const char* STR_GEMMA_PATH = "-GEMMA-path";
const char* STR_BLAST_PATH = "-Blast-path";
const char* STR_PHANTOMJS_PATH = "-phantomjs-path";
const char* STR_RSCRIPT_PATH = "-Rscript-path";
const char* STR_NO_PREVIEW = "-no-preview";

//TODO: seeveral questions are still unclear if we use the Freq count mode (how to run bugwas, the coloring, etc...). For now I am disabling this option
//const char* STR_COUNT_MODE = "-count";
//TODO: seeveral questions are still unclear if we use the Freq count mode (how to run bugwas, the coloring, etc...). For now I am disabling this option

//variables controlling where the executable is
string DBGWAS_lib = ""; //set on first command in main.cpp


bool skip1 = false;
bool skip2 = false;
bool hasNewickFile = false;
bool presenceAbsenceCountMode = false;
boost::variant< int, double > SFF;
bool thereIsNucleotideDB = false;
string nucleotideDBPath;
bool thereIsProteinDB = false;
string proteinDBPath;
string gemmaPath;
string blastPath;
string phantomjsPath;
string RscriptPath;
bool noPreview = false;

//global vars used by both programs
Graph* graph;
vector< UnitigIdStrandPos >* nodeIdToUnitigId;
vector< Strain >* strains = NULL;

void populateParser (Tool *tool) {
  // We add some custom arguments for command line interface
  tool->getParser()->push_front (new OptionOneParam (STR_RSCRIPT_PATH, "Path to Rscript.",  false, "Rscript"));
  tool->getParser()->push_front (new OptionNoParam (STR_NO_PREVIEW, "Do not produce the components preview in the summary output page.",  false));
  tool->getParser()->push_front (new OptionOneParam (STR_PHANTOMJS_PATH, "Path to phantomjs executable (DBGWAS was tested only with version 2.1.1).",  false, "<DBGWAS_lib>/phantomjs"));
  tool->getParser()->push_front (new OptionOneParam (STR_BLAST_PATH, "Path to the directory containing the Blast suite (should contain at least blastn, blastx, and makeblastdb).",  false, "<DBGWAS_lib>/"));
  tool->getParser()->push_front (new OptionOneParam (STR_GEMMA_PATH, "Path to the GEMMA executable.",  false, "<DBGWAS_lib>/gemma.0.93b"));
  tool->getParser()->push_front (new OptionOneParam (STR_MAF_FILTER, "Minor Allele Frequency Filter.",  false, "0.01"));
  tool->getParser()->push_front (new OptionOneParam (STR_MAX_NEIGHBOURHOOD, "Denotes the neighbourhood to be considered around the significant unitigs.",  false, "5"));
  tool->getParser()->push_front (new OptionOneParam (STR_SFF, "Denotes the Significant Features Filter - the features (or patterns) selected to create a visualisation around them. If it is a float number n, then only the features with q-value<=n are selected. If it is an integer n, then only the n first features are selected. Take a look at the output/step2/patterns.txt file to get a list of features ordered by q-value to better choose this parameter (re-run the tool with -skip2 in order to directly produce the visualisation of the features selected by your parameter).",  false, "100"));
  tool->getParser()->push_front (new OptionNoParam (STR_SKIP2, "Skips Steps 1 and 2, running only Step 3. Assumes that Steps 1 and 2 were correctly run and folders \"step1\" and \"step2\" are present in the output folder.",  false));
  tool->getParser()->push_front (new OptionNoParam (STR_SKIP1, "Skips Step 1, running only Steps 2 and 3. Assumes that Step 1 was correctly run and folder \"step1\" is present in the output folder.",  false));
  tool->getParser()->push_front (new OptionOneParam (STR_OUTPUT, "Path to the folder where the final and temporary files will be stored.",  false, "output"));
  tool->getParser()->push_front (new OptionOneParam (STR_PROTEIN_DB, "A list of Fasta files separated by comma containing annotations in a protein alphabet format (e.g.: -pt-db path/to/file_1.fa,path/to/file_2.fa,etc). You can customize these files to work better with DBGWAS (see https://gitlab.com/leoisl/dbgwas/tree/master#customizing-annotation-databases).",  false));
  tool->getParser()->push_front (new OptionOneParam (STR_NUCLEOTIDE_DB, "A list of Fasta files separated by comma containing annotations in a nucleotide alphabet format (e.g.: -nc-db path/to/file_1.fa,path/to/file_2.fa,etc). You can customize these files to work better with DBGWAS (see https://gitlab.com/leoisl/dbgwas/tree/master#customizing-annotation-databases).",  false));
  tool->getParser()->push_front (new OptionOneParam (STR_NEWICK_PATH, "Optional path to a newick tree file. If (and only if) a newick tree file is provided, the lineage effect analysis is computed and PCs figures are generated.",  false));
  tool->getParser()->push_front (new OptionOneParam (STR_KSKMER_SIZE, "K-mer size.",  false, "31"));
  tool->getParser()->push_front (new OptionOneParam (STR_STRAINS_FILE, "A text file describing the strains containing 3 columns: 1) ID of the strain; 2) Phenotype (0/1/NA); 3) Path to a multi-fasta file containing the sequences of the strain. This file needs a header. Check the sample_example folder or https://gitlab.com/leoisl/dbgwas/raw/master/sample_example/strains for an example.",  true));

  //TODO: seeveral questions are still unclear if we use the Freq count mode (how to run bugwas, the coloring, etc...). For now I am disabling this option
  //tool->getParser()->push_front (new OptionOneParam (STR_COUNT_MODE, "The count mode. If \"01\", then the count mode is seen as presence/absence. If \"Freq\", then the count mode is seen as the frequency",  false, "01"));
  //TODO: seeveral questions are still unclear if we use the Freq count mode (how to run bugwas, the coloring, etc...). For now I am disabling this option
}
