import argparse
import collections
import csv
import itertools
import logging
import math
import numpy as np
import scipy.stats as st
import scipy.sparse as sp
import sys

# Defaults
nc_thresh = 0.05
consideration_threshold = 30

def acc_dict_update(d, acc):
    if acc in d:
        return d[acc]
    else:
        n = len(d)
        d[acc] = n
        return n
    

# Assumptions:
#   * The query sequences are the ones we are interested in, and that we compute
#     NC for.
#   * Query sequences may or may not be in the database that has been compared
#     with.
#   * A sequence may (below) have two IDs: one as a query and one as a subject,
#     ie, part of DB. We don't care, because only query IDs, which will be used
#     as row IDs, will be used later. Column order does not matter.

def read_blast_tab(file_handles, transform=None):
    '''
    Read the file containing standard Blast tabular results (from BLAST or DIAMOND)
    and return the needed data.

    If transform is set, then the square root of bitscore is used.

    Returns:
      * A map (dict) ID to accession, for later presentation
      * A dict mapping pairs of IDs to bitscore
      * The number of queries
      * The number of subjects (#seqs in comparison DB)
    '''
    q_accession2id = dict()       # We need to map identifiers to integers
    s_accession2id = dict()       # ... more
    id2accession = dict()       # ... and back
    similar_pairs = dict()
    singletons = list()

    for fh in file_handles:
        reader = csv.reader(fh, delimiter='\t')
        for row in reader:
            query_id, subject_id, bit_score = row[0], row[1], row[11]
            if subject_id == '*':
                singletons.append(query_id)
            else:
                id1 = acc_dict_update(q_accession2id, query_id)
                id2 = acc_dict_update(s_accession2id, subject_id)
                id2accession[id1] = query_id
                if transform:
                    similar_pairs[(id1, id2)] = math.sqrt(float(bit_score))
                else:
                    similar_pairs[(id1, id2)] = float(bit_score)

    return id2accession, similar_pairs, max(q_accession2id.values()) + 1, max(s_accession2id.values()) + 1, singletons


def read_blast_3col_format(file_handles, transform=None):
    '''
    Same is read_blast_tab, but for three-column data like NC_standalone input.
    '''
    q_accession2id = dict()       # We need to map identifiers to integers
    s_accession2id = dict()       # ... more
    id2accession = dict()         # ... and back
    similar_pairs = dict()

    for fh in file_handles:
        reader = csv.reader(fh, delimiter=' ')
        for row in reader:
            queryid, subject_id, bit_score = row
            if transform:
                bit_score = math.sqrt(bit_score)
            id1 = acc_dict_update(q_accession2id, queryid)
            id2 = acc_dict_update(s_accession2id, subject_id)
            id2accession[id1] = queryid
            if transform:
                similar_pairs[(id1, id2)] = math.sqrt(float(bit_score))
            else:
                similar_pairs[(id1, id2)] = float(bit_score)

    return id2accession, similar_pairs, max(q_accession2id.values()) + 1, max(s_accession2id.values()) + 1


def scores_to_comparison_matrix(similar_pairs, nrows, ncols):
    '''
    Given a list of triples (idx1, idx2, bitscore), return 
    a sparse matrix containing the same information. The indices are integers pointing to the
    rows and columns where the score should be put. The nrows and ncols variables are
    needed to avoid computing maximum indices. 
    '''
    logging.info(f'Preparing matrix with similarity data ({nrows} by {ncols}, but a sparse matrix)')
    row_indices, col_indices = zip(*similar_pairs.keys())
    comparisons = sp.csr_array((list(similar_pairs.values()), (row_indices, col_indices)), shape=(nrows, ncols), dtype=np.float32)
    logging.info(f'There are {comparisons.getnnz()} (={100*comparisons.getnnz() / (nrows*ncols):.4}%) non-zero elements in the matrix')
    return comparisons


def pearson_correlation(comparison_matrix, i, j, cache={}):
    '''
    Compute the Pearson correlation for values of two rows in matrix.
    There are convenient functions for this in SciPy and elsewhere, 
    but it seems as if I have to roll my own for rows of sparse matrices.

    $\sum_i (x_i - \bar x)(y_i - \bar y)$ can be written as a vector 
    computation: $xy^T - x (1^T \bar y) - (\bar x1) y + n\bar x\bar y$. Input vectors
    $x$ and $y$ are sparse, but $x - 1\bar x$ would not be sparse. However,
    the multiplications in the vector computation are fast due to sparsity,
    end the final sum are on scalars. 
    '''
    row_i = comparison_matrix.getrow(i)
    row_j = comparison_matrix.getrow(j)
    n_cols = comparison_matrix.shape[1]

    if i in cache:
        m_i, s_i, root_variance_i = cache[i]
    else:
        m_i = row_i.mean()
        s_i = row_i.sum()
        root_variance_i = math.sqrt(row_i.dot(row_i.transpose())[0,0] - 2 * m_i*s_i + n_cols * m_i * m_i)
        cache[i] = (m_i, s_i, root_variance_i)
    if j in cache:
        m_j, s_j, root_variance_j = cache[j]
    else:
        m_j = row_j.mean()
        s_j = row_j.sum()
        root_variance_j = math.sqrt(row_j.dot(row_j.transpose())[0,0] - 2 * m_j*s_j + n_cols * m_j * m_j)
        cache[j] = (m_j, s_j, root_variance_j)

    numerator = row_i.dot(row_j.transpose())[0,0] - m_i * s_j - m_j * s_i + n_cols * m_i * m_j
    denominator = root_variance_i * root_variance_j

    correlation_coefficient = numerator / denominator
    return correlation_coefficient
    

def fst(x):
    '''Return first element of tuple/list x.'''
    return x[0]


def snd(x):
    '''Return second element of tuple/list x.'''
    return x[1]


def find_good_pairs(similarities, threshold):
    ref_hits = dict()
    for (a, b), score in similarities.items():
        if score >= threshold:
            if b in ref_hits:
                ref_hits[b].append(a)
            else:
                ref_hits[b] = [a]
            
    good_pairs = dict()
    for target, queries in ref_hits.items():
        for a, b in itertools.combinations(queries, 2):
            if (a,b) not in good_pairs:
                good_pairs[(a, b)] = True
    return good_pairs


def find_linked_pairs(comparison_matrix):
    '''We want to extract pairs of sequences that are either directly linked
    because the have reported similarity, or indirectly link through other
    similar sequences. In this latter case, it is assumed that the sequence
    similarity is so low that BLAST missed it, for example, but similarity
    should still be investigated.

    Returns an iterator with triples (idx1, idx2, grouplabel). The idx elements
    are row indices into the comparison matrix, and grouplabel is an integer
    indicating the subset the indexed sequences belong to.

    Implemented using connected components in a graph, so basically
    linear time to compute.

    I just realized that of two sequences are linked without having a
    hit to a common reference sequence, then their NC will be 0.
    '''
    n_rows, n_cols = comparison_matrix.shape
    # Add an n_rows x n_rows block on the left
    left_padding = sp.csr_array((n_rows, n_rows))
    G = sp.hstack( (left_padding, comparison_matrix) )
    
    # We need a square matrix, so will add some zeros as padding.
    lower_padding = sp.csr_array((n_cols, n_rows+n_cols))
    G = sp.vstack( (G, lower_padding) )

    # Now compute connected components
    n_components, labeling = sp.csgraph.connected_components(G, directed=False)
    del G                       # ... and forget this large matrix
    labeling = labeling[:n_rows]
    n_components = len(collections.Counter(labeling))
    logging.info(f'Identified {n_components} groups of sequences to compute NC for')
    
    # Use enumerate to get indices of group labels
    # Sort pairs based on group label
    # Then use groupby to collect the indices with the same group label
    groups = itertools.groupby(sorted(enumerate(labeling), key=snd), key=snd)
    for label, labeling in groups:
        component = map(fst, labeling) # extract indices for component
        for a, b in itertools.combinations(component, 2):
            yield a, b, label


def neighborhood_correlation(id2accession, similar_pairs, n_queries, n_ref_seqs, threshold=0):
    logging.info('Preparing NC data')
    comparison_matrix = scores_to_comparison_matrix(similar_pairs, n_queries, n_ref_seqs)

    good_pairs = find_good_pairs(similar_pairs, threshold)
    logging.info(f'Identified {len(good_pairs)} sequence pairs to compute NC for.')

    logging.info('Starting NC computations')
    cache = dict()          # Using a cache to avoid computing denominator factors more than once
    for a, b in good_pairs:
        cor = pearson_correlation(comparison_matrix, a, b, cache)
        yield id2accession[a], id2accession[b], cor
            


def snc_argparser():
    ap = argparse.ArgumentParser()
    ap.add_argument('infile', nargs='+', type=argparse.FileType('r'),
                    default=sys.stdin,
                    help='The infile is the path to a file containing BLAST or DIAMOND output in tabular format (using --outfmt 6 in both tools)')
    ap.add_argument('-a', '--all-vs-all', action='store_true',
                    help='If your comparisons are all-versus-all. Otherwise it is assumed that the sequences you are interested in have been compared with reference database.')
    ap.add_argument('-c', '--consider', type=float, metavar='TRESHOLD', default=consideration_threshold,
                    help=f'Consideration threshold. Only consider pairs of sequences linked by similarities (maybe in several steps) with this bitscores or higher. (Default:{consideration_threshold})')    
    ap.add_argument('-3', '--three-col', action='store_true',
                    help='Actually, assume the input file has three columns (acc1, acc2, and bitscore) separated by single blankspace, like NC_standalone.')
    ap.add_argument('-s', '--square-root', action='store_true',
                    help='Transform bitscores with square root')
    ap.add_argument('-t', '--nc-thresh', type=float, default=nc_thresh,
                    help=f'NC reporting threshold. Calculated values below this threshold will not be reported. Default: {nc_thresh}')
    ap.add_argument('-v', '--verbose', action='store_true', 
                    help='Output some progress information')
    return ap


def nc_main():
    ap = snc_argparser()
    args = ap.parse_args()

    if args.all_vs_all:
        sys.exit('Not implemented yet. Requires a symmetric scoring matrix. Currently, comparisonmatrix[i][i] refers to "query i" and "subject i", which is not the same sequence.')

    if args.verbose:
        logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%Y-%m-%d %H:%M:%S', level=logging.INFO)
    else:
        logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
        

    logging.info('Reading data')
    singletons = None
    if not args.three_col:
        id2accession, similarities, n_queries, n_ref_seqs, singletons = read_blast_tab(args.infile, args.square_root) # Note: args.infile is a list of filehandles
    else:
        id2accession, similarities, n_queries, n_ref_seqs = read_blast_3col_format(args.infile, args.square_root)

    if singletons:
        logging.info(f'Noted {len(singletons)} sequences without a hit in the reference data.')

    if similarities:
        counter = 0
        for acc1, acc2, nc in neighborhood_correlation(id2accession, similarities, n_queries, n_ref_seqs, args.consider):
            counter += 1
            if counter % 10000 == 0:
                logging.info(f'{counter} pairs analyzed')
            if nc >= args.nc_thresh:
                print(acc1, acc2, round(nc, 3))
        logging.info(f'{counter} pairs analyzed')
    logging.info(f'Done.')

if __name__ == '__main__':
    nc_main()
        
