#!/usr/bin/env python
# vim: set fileencoding=<utf-8> :
# by John Lees

import sys
import networkx as nx

__version__ = "1.0.0"

def get_options():

    import argparse

    parser = argparse.ArgumentParser(description='Extend unitig hits',
                                     prog='extend-hits')

    graphGroup = parser.add_argument_group('Graph input')
    graph = graphGroup.add_mutually_exclusive_group(required=True)
    graph.add_argument('--prefix',
            help='Prefix for DBG .nodes and .edges.dbg files')
    graph.add_argument('--load-graph',
            help='Load a previously saved graph')

    required = parser.add_argument_group('Required arguments')
    required.add_argument('--unitigs',
            help='Significant unitigs to extend',
            required=True)

    parser.add_argument('--save-graph',
            help='Prefix to save graph after loading',
            default=None)
    parser.add_argument('--length',
            help='Target length to extend to',
            default=100,
            type=int)

    parser.add_argument('--version', action='version',
                       version='%(prog)s '+__version__)

    return parser.parse_args()

def main():
    args = get_options()

    # read in graph as networkx
    sys.stderr.write("Loading graph\n")
    unitig_ids = {}
    if args.prefix:
        # Add nodes first
        node_list = []
        with open(args.prefix + ".nodes", 'r') as node_file:
            for node in node_file:
                (node_id, node_seq) = node.rstrip().split("\t")
                node_list.append((int(node_id), dict(seq=node_seq, seq_len=len(node_seq))))
                unitig_ids[node_seq] = node_id

        G = nx.DiGraph()
        G.add_nodes_from(node_list)

        # add edges
        edge_list = []
        with open(args.prefix + ".edges.dbg", 'r') as edge_file:
            for edge in edge_file:
                (start, end, label) = edge.rstrip().split("\t")

                edge_list.append((int(start), int(end), G.nodes[int(end)]['seq_len']))
                edge_list.append((int(end), int(start), G.nodes[int(start)]['seq_len']))

        G.add_weighted_edges_from(edge_list)

        if args.save_graph:
            nx.write_gpickle(G, args.save_graph + ".gpickle")
    else:
        # Load graph if already produced
        G = nx.read_gpickle(args.load_graph)
        for (node_id, node_data) in G.nodes(data=True):
            unitig_ids[node_data['seq']] = str(node_id)

    # extend each unitig
    sys.stderr.write("Extending unitigs\n")
    with open(args.unitigs, 'r') as unitig_file:
        for unitig in unitig_file:
            unitig_name = int(unitig_ids[unitig.rstrip()])
            paths = nx.single_source_dijkstra_path(G, unitig_name, args.length)
            extensions = []
            for path in paths.values():
                if len(path) > 1:
                    extension = [G.nodes[x]['seq'] for x in path]
                    extensions.append(("".join(extension)))
            print(",".join(extensions))

if __name__ == '__main__':
    main()

    sys.exit(0)