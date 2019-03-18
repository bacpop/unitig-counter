/*
 * dist_example.cpp
 * Demonstrates used of Cdbg class
 *
 *  This code: John Lees 2019
 *
 */

#include <boost/program_options.hpp>
#include "version.h"
#include "node_dists.hpp"

namespace po = boost::program_options; // Save some typing

// Use boost::program_options to parse command line input
int parseCommandLine (int argc, char *argv[], po::variables_map& vm)
{
   int failed = 0;

   po::positional_options_description mode;
   mode.add("mode", -1);

   po::options_description graph("Graph options");
   graph.add_options()
    ("graph", po::value<string>(), "Prefix of graph files")
    ("nodes", po::value<string>(), "Name of .node file")
    ("edges", po::value<string>(), "Name of .edges.dbg file");

   po::options_description dist("Distance options");
   dist.add_options()
    ("source", po::value<string>(), "Sequence of source node")
    ("target", po::value<string>(), "Sequence of target node");

   po::options_description extend("Extending options");
   extend.add_options()
    ("unitigs", po::value<string>(), "File containing unitigs to extend")
    ("length", po::value<int>()->default_value(100), "Maximum extension length")
    ("repeats", "Allow loops in extensions");

   po::options_description other("Other options");
   other.add_options()
    ("version", "prints version and exits")
    ("help,h", "full help message");

   po::options_description all;
   all.add(graph).add(dist).add(extend).add(other);

   try
   {
      po::store(
          po::command_line_parser(argc, argv).
          positional(mode).
          run(), vm);

      if (vm.count("help"))
      {
         cerr << "cdbg-ops dist" << endl;
         cerr << "Calculate distance between two nodes" << endl << endl;
         cerr << "cdbg-ops extend" << endl;
         cerr << "Extend sequence around a node by finding paths through it" << endl << endl;
         cerr << all << endl;
         failed = 1;
      }
      else if (vm.count("version"))
      {
         cout << VERSION << endl;
         failed = 1;
      }
      else
      {
         po::notify(vm);
         failed = 0;

         // Check input files exist, and can stat
         if (vm["mode"].as<string>() != "dist" && vm["mode"].as<string>() != "extend")
         {
            cerr << "Possible modes as 'dist' or 'extend" << endl;
            failed = 1;
         }
      }

   }
   catch (po::error& e)
   {
      // Report errors from boost library
      std::cerr << "Error in command line input: " << e.what() << "\n";
      std::cerr << "Run 'cdbg-ops --help' for full option listing\n\n";
      std::cerr << all << endl;

      failed = 1;
   }

   return failed;
}

// Simple example program
int main (int argc, char *argv[])
{
    // Do parsing and checking of command line params
    po::variables_map vm;
    if (parseCommandLine(argc, argv, vm))
    {
        return 1;
    }

    cerr << "Reading graph" << endl;

    string nodes, edges;
    if (vm.count("graph"))
    {
        nodes = vm["graph"].as<string>() + ".nodes";
        edges = vm["graph"].as<string>() + ".edges.dbg";
    }
    else if (vm.count("nodes") && vm.count("edges"))
    {
        nodes = vm["nodes"].as<string>();
        edges = vm["edges"].as<string>();
    }
    else
    {
        cerr << "Must give input graph with --graph or --nodes and --edges" << endl;
    }
    Cdbg graphIn(nodes, edges);

    if (vm["mode"].as<string>() == "dist")
    {
        cerr << "Calculating distances from node" << endl;
        vector<int> dists = graphIn.node_distance(vm["source"].as<string>());
        cout << dists[graphIn.get_vertex(vm["target"].as<string>())] << endl;
    }
    else if (vm["mode"].as<string>() == "extend")
    {
        cerr << "Getting paths through nodes" << endl;

        if (vm.count("unitigs"))
        {
            ifstream unitigsIst(vm["unitigs"].as<string>().c_str());
            if (!unitigsIst)
            {
                throw std::runtime_error("Could not open unitig file " + vm["unitigs"].as<string>() + "\n");
            }

            string sequence;
            while (unitigsIst >> sequence)
            {
                vector<string> paths = graphIn.extend_hits(sequence, vm["length"].as<int>(), vm.count["repeats"]);
                for (auto it = paths.begin(); it != paths.end(); ++it)
                {
                    cout << *it;
                    if (it != paths.end() - 1)
                    {
                        cout << ",";
                    }
                }
                cout << endl;;
            }

            unitigsIst.close();
        }
        else
        {
            cerr << "Must provide sequences in file with --unitigs" << endl;
            return 1;
        }

    }

    return 0;
}