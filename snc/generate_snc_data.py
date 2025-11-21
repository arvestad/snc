import argparse
import itertools
import random
import numpy.random as npr

description='''
Generate test data for snc (and similar tools.)

Sequence pairs that are supposed to be homologous are given a similarity results
with other sequences that are drawn from a normal distribution with mean `hom_sim`
and variance `hom_var`.

Non-homologous sequence pairs are given similarity data drawn from
a uniform distribution, [0, non_sim], where `non_sim` is given
by the option `--non_sim`.

`hom_sim` should be significantly higher than `non_sim` for easy testing.

Note that no domain overlap among sequences is simulated in this data.
'''

def setup_argparse():
    ap = argparse.ArgumentParser(description=description)
    ap.add_argument('n_clusters', type=int, default=5,
                    help='Choose the number of homologous sequence clusters.')
    ap.add_argument('cluster_size', type=int, default=5,
                    help='Decide how large every cluster is.')
    ap.add_argument('n_ref_seqs', type=int, default=20,
                    help='The number of simulated comparisons for each simulated homologous pair.')

    ap.add_argument('--hom_sim', type=float, default=100.0,
                    help='The average similarity for homologous sequence pairs.')
    ap.add_argument('--hom_var', type=float, default=10.0,
                    help='The variance on similarity for homologous sequence pairs.')

    ap.add_argument('--non_sim', type=float, default=50.0,
                    help='The maximum similarity for non-homologous sequence pairs.')

    ap.add_argument('--min_sim', type=float, default=25.0,
                    help='Simulated minimum detectable similarity. Pairs with lower similarity are not output.')
    
    return ap


def scoring_line(acc1, acc2, sc):
    print(f'{acc1}	{acc2}	{sc}')
    
   
counter = 0              # To give unique reference seq-ids
def correlated_pair(i, j, args):
    global counter
    for k in range(args.n_ref_seqs):
        ref_acc = f'r{counter}'
        counter += 1
        score = npr.default_rng().normal(args.hom_sim, args.hom_var, 2)
        scoring_line(i, ref_acc, score[0])
        scoring_line(j, ref_acc, score[1])

def uncorrelated_pair(i, j, args):
    global counter
    for k in range(args.n_ref_seqs):
        ref_acc = f'r{counter}'
        counter += 1
        score = random.uniform(0, args.non_sim)
        if score > args.min_sim:
            scoring_line(i, ref_acc, score)
        score = random.uniform(0, args.non_sim)
        if score > args.min_sim:
            scoring_line(j, ref_acc, score)
        

def main():
    ap = setup_argparse()
    args = ap.parse_args()

    sequences = {}              # Store info about our pretend sequences
    cluster_idx = 0
    for c in range(args.n_clusters):
        for i in range(args.cluster_size):
            accession = f'c{c}_s{i}' # Cluster c and sequence i in that cluster
            sequences[accession] = c # Map accession to cluster idx

    for a, b in itertools.pairwise(sequences): # Requires Python 3.10 or later.
        if sequences[a] == sequences[b]:       # Same cluster idx
            correlated_pair(a, b, args)
        else:
            uncorrelated_pair(a, b, args)


if __name__ == '__main__':
    main()
        
