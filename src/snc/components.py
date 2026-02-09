import argparse
from Bio import SeqIO
import json
import sys
import math
import networkx as nx
import statistics
import tabulate

import signal
signal.signal(signal.SIGPIPE, signal.SIG_DFL)


def create_arg_parser():
    ap = argparse.ArgumentParser(description='Create clusters based on scores of NC type.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument('infile', help='Filename to data in 3-column format: "id1  id2  float". Each line represents a weighted edge in a graph.' )
    ap.add_argument('-c', '--count', action='store_true',
                    help='Count the number of clusters')
    ap.add_argument('-d', '--discard-size', type=int, default=2,
                    help='Discard components smaller (<) than this. This does not affect the --info or --count options. Default: %(default)s')
    ap.add_argument('-i', '--info', action='store_true',
                    help='Output information about the input to stderr')
    ap.add_argument('-j', '--json', metavar='FILENAME',
                    help='Save information and clusters (accessions only) in the named JSON file.')
    ap.add_argument('-p', '--prefix', default='c',
                    help='If using -s, use this string as a prefix for the output file names. Default is "c".')
    ap.add_argument('-s', '--seqfile', type=argparse.FileType('r'),
                    help='If used, output actual cluster FASTA files containing sequences from this file.')
    ap.add_argument('-t', '--threshold', type=float, default=0.7,
                    help='Threshold that defines cluster. Keep those edges whose weights are at least as high as the threshold.')
    ap.add_argument('-oh', '--output-help', action='store_true',
                    help='Output an explanation of the statistics output.')
    return ap


def output_help():
    """
    Describe the statistics that 'components' report.
    """
    info = 'The --info and --json output formats report a number of observations of the input data and output from clustering. Here is a summary.'
    variables = [
        ['n_nodes', 'Number of entities (probably sequences) for which there are at least one NC score reported.'],
        ['n_edges', 'Number of nodes pairs with a "good" NC score, i.e., NC larger than set by --threshold.'],
        ['discard_size', 'The smallest allowed component size. Smallest possible value is 2.'],
        ['n_components', 'Number of output clusters, i.e., number of connected components.'],
        ['n_discarded_nodes', 'Number of  entities in the input which are not in any component.'],
        ['median', 'Median size of components.'],
        ['largest', 'Size of the largest component.'],
        ['quantiles', 'The three sizes that split the distribution of component sizes into quartiles.'],
        ['graph_density', 'The number of edges in the graph, induced by the input (and --threshold), divided by the number of node pairs.'],
        ['mean_component_density', 'Reports the mean of component subgraph densities.'],
        ]
    table = tabulate.tabulate(variables, headers=['Variable', 'Explanation'])
    return info + '\n\n' + table


def get_stats(graph, components, n_edges, n_discarded_edges):
    stats = {}
    n_nodes = graph.number_of_nodes()
    stats['n_nodes'] = n_nodes
    sizes = list(map(len, components))
    stats['n_components'] = len(components)
    stats['n_discarded_nodes'] = n_nodes - sum(sizes) # nr seqs not in clusters
    stats['median'] = statistics.median(sizes)
    stats['largest'] = max(sizes)
    try:
        quantiles = statistics.quantiles(sizes)
    except statistics.StatisticsError:
        quantiles = 'N/A'
    stats['quantiles'] = quantiles

    vertex_pairs_in_components = map(lambda c: len(c)**2, components)
    component_densities = map(lambda c: nx.density(graph.subgraph(c)), components)
    stats['graph_density'] = nx.density(graph)
    stats['mean_component_density'] = statistics.mean(component_densities)

    return stats
    

def make_json_file(filename, components, stats):
    stats['components'] = list(map(lambda c: list(c), components))
    with open(filename, 'w') as jh:
        json.dump(stats, jh, indent=4)

    

def make_graph(infilename, threshold=0):
    graph = nx.Graph()
    n_edges = 0
    n_discarded_edges = 0

    with open(infilename) as h:
        for lineno, line in enumerate(h):
            try:
                id1, id2, val_s = line.split()
            except Exception as e:
                print(f'Error on line {lineno} in "{infilename}": {str(e)}', file=sys.stderr)
                sys.exit(1)
                
            val = float(val_s)
            if val > threshold:
                graph.add_edge(id1, id2)
                n_edges += 1
            else:
                graph.add_node(id1)
                graph.add_node(id2)
                n_discarded_edges += 1
    return graph, n_edges, n_discarded_edges


def get_components(graph, discard_size):
    components = nx.connected_components(graph)
    components = list(filter(lambda c: len(c) >= discard_size, components)) # Remove components deemed too small
    return components
              

def print_components(components):
    for component in components:
        print(' '.join(component))


def accession_lookup(table, acc_string):
    '''
    The accession string can contain several accessions, like
    "tr|G3XF63|G3XF63_9SPER", and we want to try all of them
    in the lookup.
    '''
    for acc in acc_string.split('|'):
        if acc in table:
            return table[acc]

    # accession not found in table
    return None


def create_fasta_cluster_files(components, seqfile, prefix):
    '''
    Create a Fasta file for each component. The function is slow, because
    the sequence file is traversed once (because assumed large) and
    clusterfiles are opened again and again. Saving RAM at the cost of
    speed.
    '''
    cluster_dict = dict()
    cluster_idx = 0
    assert len(components) > 0
    n_digits_needed = math.floor(math.log10(len(components))) + 1

    # First note for each sequence in which file to put it
    for component in components:
        cluster_name = f'{prefix}{cluster_idx:0{n_digits_needed}}.fa'
        cluster_idx += 1
        for v in component:
            cluster_dict[v] = cluster_name

    # Then go through all sequences and append them to the right cluster
    for record in SeqIO.parse(seqfile, 'fasta'):
        clusterfile = accession_lookup(cluster_dict, record.id)
        if clusterfile:
            with open(clusterfile, 'a') as h:
                print(f'>{record.id} {record.description}', file=h)
                print(record.seq, file=h)
        else:
            # This record was not among the clusters
            # print(f'Missing {record.id}', file=sys.stderr)
            continue


def components_main():
    ap = create_arg_parser()
    args = ap.parse_args()

    if args.output_help:
        print(output_help())
        sys.exit(0)

    graph, n_edges, n_discarded_edges = make_graph(args.infile, args.threshold)
    components = get_components(graph, args.discard_size)

    if args.count:
        print(len(components))
        return

    stats = get_stats(graph, components, n_edges, n_discarded_edges)
    stats['threshold'] = args.threshold
    stats['discard_size'] = args.discard_size
    stats['n_edges'] = n_edges
    stats['n_discarded_edges'] = n_discarded_edges # pairs we counted NC for but it is below threshold

    if args.info:
        print('# Graph statistics', file=sys.stderr)
        for key, val in stats.items():
            print(f'{key}: \t{val}', file=sys.stderr)

    if args.json:
        try:
            make_json_file(args.json, components, stats)
        except Exception as e:
            print(f'Error saving JSON data to "{args.json}": {str(e)}', file=sys.stderr)
    else:
        print_components(components)

    if args.seqfile:
        print('# Creating a FASTA file for each cluster', file=sys.stderr)
        create_fasta_cluster_files(components, args.seqfile, args.prefix)





if __name__=='__main__':
    components_main()
