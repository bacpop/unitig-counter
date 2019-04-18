/*
 * graph_ops.cpp
 * Command line interface to Cdbg class
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
   mode.add("mode", 1);

   po::options_description graph("Graph options");
   graph.add_options()
    ("graph", po::value<string>(), "Prefix of graph files")
    ("nodes", po::value<string>(), "Name of .node file")
    ("edges", po::value<string>(), "Name of .edges.dbg file");

   po::options_description dist("Distance options");
   dist.add_options()
    ("source", po::value<string>(), "Sequence of source node")
    ("source-list", po::value<string>(), "File containing sequences of source nodes")
    ("target", po::value<string>(), "Sequence of target node")
    ("all", "Generate distances to all other unitigs");

   po::options_description extend("Extending options");
   extend.add_options()
    ("unitigs", po::value<string>(), "File containing unitigs to extend")
    ("length", po::value<int>()->default_value(100), "Maximum extension length")
    ("repeats", "Allow loops in extensions");

   po::options_description other("Other options");
   other.add_options()
    ("mode", po::value<string>(), "Mode of operation")
    ("version", "prints version and exits")
    ("help,h", "full help message");

   po::options_description all;
   all.add(graph).add(dist).add(extend).add(other);

   try
   {
      po::store(
          po::command_line_parser(argc, argv).
          positional(mode).options(all).
          run(), vm);

      if (vm.count("help"))
      {
         cerr << "cdbg-ops dist: Calculate distance between two nodes" << endl;
         cerr << "cdbg-ops extend: Extend sequence around a node by finding paths through it" << endl;
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
         if (vm.count("mode") != 1 ||
              (vm["mode"].as<string>() != "dist" && vm["mode"].as<string>() != "extend"))
         {
            cerr << "Possible modes are 'dist' or 'extend" << endl;
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
    if (argc == 1)
    {
        cerr << "cdbg-ops dist --source AATCG --target TTGC" << endl;
        cerr << "cdbg-ops extend --unitigs significant_hits.txt" << endl;
        return 1;
    }
    else if (parseCommandLine(argc, argv, vm))
    {
        return 1;
    }

    cerr << "Reading graph" << endl;

    // Get nodes and edges files needed to create Cdbg object
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

    // Distance mode
    if (vm["mode"].as<string>() == "dist")
    {
        cerr << "Calculating distances from node" << endl;

        if (!vm.count("target") && !vm.count("all"))
        {
            throw std::runtime_error("Need to provide a target or all in dist mode");
        }

        vector<string> sources;
        if (vm.count("source-list"))
        {
            ifstream sourceIst(vm["source-list"].as<string>().c_str());
            if (!sourceIst)
            {
                throw std::runtime_error("Could not open node file " + vm["source-list"].as<string>() + "\n");
            }

            string sequence;
            while (sourceIst >> sequence)
            {
                sources.push_back(sequence);
            }

            sourceIst.close();
        }
        else if (vm.count("source"))
        {
            sources.push_back(vm["source"].as<string>());
        }
        else
        {
            throw std::runtime_error("Need to provide a source in dist mode");
        }

        cout << "Source\tTarget\tDistance" << endl;
        for (auto it = sources.begin(); it != sources.end(); it++)
        {
            vector<int> dists = graphIn.node_distance(*it);
            if (vm.count("all"))
            {
                for (size_t i = 0; i < dists.size(); i++)
                {
                    cout << *it << graphIn.node_seq(i) << "\t" << dists[i] << endl;
                }
            }
            else if (vm.count("target"))
            {
                cout << *it << dists[graphIn.get_vertex(vm["target"].as<string>())] << endl;
            }
        }

    }
    // Extend mode
    else if (vm["mode"].as<string>() == "extend")
    {
        cerr << "Getting paths through nodes" << endl;

        if (vm.count("unitigs"))
        {
            // Read in unitigs to extend
            ifstream unitigsIst(vm["unitigs"].as<string>().c_str());
            if (!unitigsIst)
            {
                throw std::runtime_error("Could not open unitig file " + vm["unitigs"].as<string>() + "\n");
            }

            string sequence;
            while (unitigsIst >> sequence)
            {
                // Get extensions and print
                vector<string> paths = graphIn.extend_hits(sequence, vm["length"].as<int>(), vm.count("repeats"));
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