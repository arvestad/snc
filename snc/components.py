import argparse
from Bio import SeqIO
import sys
import math
import networkx as nx
import statistics


def create_arg_parser():
    ap = argparse.ArgumentParser(description='Create clusters based on scores of NC type.')
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
    return ap


def get_component_stats(components, discard_size):
    stats = {}
    sizes = list(map(len, components))
    stats['n_singletons'] = sum(map(lambda n: n==1, sizes))
    stats['median'] = statistics.median(sizes)
    stats['largest'] = max(sizes)
    try:
        quantiles = statistics.quantiles(sizes)
    except statistics.StatisticsError:
        quantiles = 'N/A'
    stats['quantiles'] = quantiles

    stats['n_components'] = len(components)
    stats['large_enough'] = len(list(filter(lambda sz: sz>= discard_size, sizes)))
    stats['sizes'] = sizes
    return stats
    

def print_component_stats(components, discard_size):
    stats = get_component_stats(components, discard_size)

    print('\n# Component statistics', file=sys.stderr)
    print(f'Number of components: {stats["n_components"]}', file=sys.stderr)
    print(f'Number of singletons: {stats["n_singletons"]}', file=sys.stderr)
    print(f'Median component size: {stats["median"]}', file=sys.stderr)
    print(f'Largest component size: {stats["largest"]}', file=sys.stderr)
    print(f'Quantiles at: {stats["quantiles"]}', file=sys.stderr)


def make_json_file(filename, components, n_edges, n_non_edges, n_nodes, discard_size):
    stats = get_component_stats(components, discard_size)
    stats['n_edges'] = n_edges,
    stats['n_non_edges'] = n_non_edges,
    stats['n_nodes'] = n_nodes,
    stats['components'] = components,

    with open(filename, 'w') as jh:
        json.dump(stats, jh, indent=4)
    

def make_graph(infilename, threshold=0):
    graph = nx.Graph()
    n_edges = 0
    n_non_edges = 0

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
                n_non_edges += 1
    return graph, n_edges, n_non_edges
                

def print_components(components, sz_threshold):
    for component in components:
        if len(component) >= sz_threshold:
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


def create_fasta_cluster_files(components, seqfile, prefix, sz_threshold):
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
        if len(component) >= sz_threshold:
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

    graph, n_edges, n_non_edges = make_graph(args.infile, args.threshold)

    if args.info:
        print('# Graph statistics', file=sys.stderr)
        print(f'Number of edges recorded: {n_edges}', file=sys.stderr)
        print(f'Number of edges discarded (below threshold {args.threshold}): {n_non_edges}', file=sys.stderr)
        print(f'Number of vertices: {graph.number_of_nodes()}', file=sys.stderr)

    components = list(nx.connected_components(graph))
    if args.info:
        print_component_stats(components, args.discard_size)

    if args.json:
        try:
            make_json_file(args.json, components, args.discard_size)
        except Exception as e:
            print(f'Error saving JSON data to "{args.json}": {str(e)}', file=sys.stderr)

    if args.count:
        print(len(list(components)))
    else:
        print_components(components, args.discard_size)
        if args.seqfile:
            print('# Creating a FASTA file for each cluster', file=sys.stderr)
            create_fasta_cluster_files(components, args.seqfile, args.prefix, args.discard_size)




if __name__=='__main__':
    main()
