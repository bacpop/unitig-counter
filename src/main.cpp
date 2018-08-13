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

// We include the header file for the tool
#include "build_dbg.hpp"
#include "map_reads.hpp"
#include "statistical_test.h"
#include "generate_output.h"
#include "global.h"
#include "Utils.h"
#include <stdlib.h>
#include <time.h>
#define BOOST_NO_CXX11_SCOPED_ENUMS
#include <boost/filesystem.hpp>
#undef BOOST_NO_CXX11_SCOPED_ENUMS
#include <boost/regex.hpp>

using namespace std;
/********************************************************************************/

int main (int argc, char* argv[])
{
    // initialize random seed, in case we want to use rand(), we are already set
    srand (time(NULL));

    //get the path to the dir were the executable is
    DBGWAS_lib = getDirWhereDBGWASIsInstalled() + "/DBGWAS_lib/";

    try
    {
        //Build DBG
        build_dbg().run(argc, argv); //this call will set up graph and nodeIdToUnitigId
        map_reads().run(argc, argv);
        cerr << "Done!" << endl;

        //Run the statistical test
        cerr << "Step 2. Running statistical test (bugwas + gemma)..." << endl;
        statistical_test().run(argc, argv);
        cerr << "Done!" << endl;

        //Find the neighbourhood around significant unitigs...
        cerr << "Step 3. Building visualisation around significant unitigs..." << endl;
        generate_output().run(argc, argv);
        cerr << "Done!" << endl;
    }
    catch (Exception& e)
    {
        std::cerr << "EXCEPTION: " << e.getMessage() << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

