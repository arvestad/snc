import argparse
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

def read_blast_tab(fh):
    '''
    Read the file containing standard Blast tabular results (from BLAST or DIAMOND)
    and return the needed data.

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

    reader = csv.reader(fh, delimiter='\t')
    for row in reader:
        query_id, subject_id, bit_score = row[0], row[1], row[11]
        
        id1 = acc_dict_update(q_accession2id, query_id)
        id2 = acc_dict_update(s_accession2id, subject_id)
        id2accession[id1] = query_id

        similar_pairs[(id1, id2)] = float(bit_score)

    return id2accession, similar_pairs, max(q_accession2id.values()) + 1, max(s_accession2id.values()) + 1


def read_blast_3col_format(fh):
    '''
    Same is read_blast_tab, but for three-column data like NC_standalone input.
    '''
    q_accession2id = dict()       # We need to map identifiers to integers
    s_accession2id = dict()       # ... more
    id2accession = dict()         # ... and back
    similar_pairs = dict()

    reader = csv.reader(fh, delimiter=' ')
    for row in reader:
        queryid, subject_id, bit_score = row
        id1 = acc_dict_update(q_accession2id, queryid)
        id2 = acc_dict_update(s_accession2id, subject_id)
        id2accession[id1] = queryid
        similar_pairs[(id1, id2)] = float(bit_score)

    return id2accession, similar_pairs, max(q_accession2id.values()) + 1, max(s_accession2id.values()) + 1


def scores_to_comparison_matrix(similar_pairs, nrows, ncols):
    logging.info(f'Preparing matrix with similarity data ({nrows} by {ncols}, but a sparse matrix)')
    row_indices, col_indices = zip(*similar_pairs.keys())
    comparisons = sp.csr_array((list(similar_pairs.values()), (row_indices, col_indices)), shape=(nrows, ncols), dtype=np.float32)
    return comparisons


#def pearson_correlation(row_i, row_j, cache={}):
def pearson_correlation(comparison_matrix, i, j, q_cache={}, s_cache={}):
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

    if i in q_cache:
        m_i, s_i, i_variance = q_cache[i]
    else:
        m_i = row_i.mean()
        s_i = row_i.sum()
        i_variance = row_i.dot(row_i.transpose())[0,0] - 2 * m_i*s_i + n_cols * m_i * m_i
        q_cache[i] = (m_i, s_i, i_variance)
    if j in s_cache:
        m_j, s_j, j_variance = s_cache[j]
    else:
        m_j = row_j.mean()
        s_j = row_j.sum()
        j_variance = row_j.dot(row_j.transpose())[0,0] - 2 * m_j*s_j + n_cols * m_j * m_j
        s_cache[j] = (m_j, s_j, j_variance)

    numerator = row_i.dot(row_j.transpose())[0,0] - m_i * s_j - m_j * s_i + n_cols * m_i * m_j
    denominator = math.sqrt(i_variance * j_variance)

    correlation_coefficient = numerator / denominator
    return correlation_coefficient
    

def fst(x):
    return x[0]


def snd(x):
    return x[1]


def find_good_pairs(comparison_matrix):
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
    '''
    n_rows, n_cols = comparison_matrix.shape
    # Add an n_rows x n_rows block on the left
    left_padding = sp.csr_array((n_rows, n_rows))
    G = sp.hstack( (left_padding, comparison_matrix) )
    
    # We need a square matrix, so will add some zeros as padding.
    lower_padding = sp.csr_array((n_cols, n_rows+n_cols))
    G = sp.vstack( (G, lower_padding) )

    # Now compute connected components
    n_components, labeling = sp.csgraph.connected_components(G)
    del G                       # ... and forget this large matrix
    logging.info(f'Identified {n_components} groups of sequences to compute NC for')
    for label, labeling in itertools.groupby(enumerate(labeling[:n_rows]), key=snd):
        component = map(fst, labeling) # extract indices for component
        for a, b in itertools.combinations(component, 2):
            yield a, b, label


def neighborhood_correlation(id2accession, similar_pairs, n_queries, n_ref_seqs, threshold=0):
    logging.info('Preparing NC data')
    comparison_matrix = scores_to_comparison_matrix(similar_pairs, n_queries, n_ref_seqs)

    logging.info('Identifying sequence pairs to compute NC for')
    good_pairs = find_good_pairs(comparison_matrix > threshold)

    logging.info('Starting NC computations')
    q_cache = dict()          # Using a cache to avoid computing denominator factors more than once
    s_cache = dict()
    for a, b, group in good_pairs:
        cor = pearson_correlation(comparison_matrix, a, b, q_cache, s_cache)
        yield id2accession[a], id2accession[b], cor, group
            
            
def nc_main():
    ap = argparse.ArgumentParser()
    ap.add_argument('infile',
                    help='The infile is the path to a file containing BLAST or DIAMOND output in tabular format (using --outfmt 6 in both tools)')
    ap.add_argument('-a', '--all-vs-all', action='store_true',
                    help='If your comparisons are all-versus-all. Otherwise it is assumed that the sequences you are interested in have been compared with reference database.')
    ap.add_argument('-c', '--consider', type=float, metavar='TRESHOLD', default=consideration_threshold,
                    help=f'Consideration threshold. Only consider pairs of sequences linked by similarities (maybe in several steps) with this bitscores or higher. (Default:{consideration_threshold})')    
    ap.add_argument('-3', '--three-col', action='store_true',
                    help='Actually, assume the input file has three columns (acc1, acc2, and bitscore) separated by single blankspace, like NC_standalone.')
    ap.add_argument('-t', '--nc-thresh', type=float, default=nc_thresh,
                    help=f'NC reporting threshold. Calculated values below this threshold will not be reported. Default: {nc_thresh}')
    ap.add_argument('-v', '--verbose', action='store_true', 
                    help='Output some progress information')

    args = ap.parse_args()

    if args.all_vs_all:
        sys.exit('Not implemented yet. Requires a symmetric scoring matrix. Currently, comparisonmatrix[i][i] refers to "query i" and "subject i", which is not the same sequence.')

    if args.verbose:
        logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%Y-%m-%d %H:%M:%S', level=logging.INFO)
    else:
        logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
        

    with open(args.infile) as f:
        logging.info('Reading data')
        if not args.three_col:
            id2accession, similar_pairs, n_queries, n_ref_seqs = read_blast_tab(f)
        else:
            id2accession, similar_pairs, n_queries, n_ref_seqs = read_blast_3col_format(f)
            
    if similar_pairs:
        for acc1, acc2, nc, group in neighborhood_correlation(id2accession, similar_pairs, n_queries, n_ref_seqs, args.consider):
            if nc >= args.nc_thresh:
                print(acc1, acc2, round(nc, 3), group)

    logging.info('Done')

if __name__ == '__main__':
    nc_main()
        
