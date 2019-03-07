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

    parser.add_argument('nodes',
            help='DBG .nodes file')
    parser.add_argument('edges',
            help='DBG .edges.dbg file')
    parser.add_argument('unitigs',
            help='Significant unitigs to extend')

    parser.add_argument('--length',
            help='Target length to extend to',
            default=100,
            type=int)
    #TODO save graph and load graph

    parser.add_argument('--version', action='version',
                       version='%(prog)s '+__version__)

    return parser.parse_args()

def main():
    args = get_options()

    # read in graph as networkx. Add nodes first
    sys.stderr.write("Loading graph\n")
    node_list = []
    with open(args.nodes, 'r') as node_file:
        for node in node_file:
            (node_id, node_seq) = node.rstrip().split("\t")
            node_list.append((int(node_id), dict(seq=node_seq, seq_len=len(node_seq)))

    G = nx.DiGraph()
    G.add_nodes_from(node_list)

    # add edges
    edge_list = []
    with open(args.edges, 'r') as edge_file:
        for edge in edge_file:
            (start, end, label) = edge.rstrip().split("\t")
            start = int(start)
            end = int(start)

            edge1 = (node_list[start], node_list[end], len(node_list[end]))
            edge2 = (node_list[end], node_list[start], len(node_list[start]))
            edge_list.append([edge1, edge2])

    G.add_weighted_edges_from(edge_list)

    # extend each unitig
    sys.stderr.write("Extending unitigs\n")
    with open(args.unitigs, 'r') as unitig_file:
        for unitig in unitig_file:
            unitig_name = unitig.rstrip()
            paths = nx.single_source_dijkstra_path(G, unitig_name, args.length)
            extensions = []
            for extension in paths:
                extension.append(("".join(extension)))
            print(",".join(extensions) + "\n")

if __name__ == '__main__':
    main()

    sys.exit(0)