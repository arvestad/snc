import argparse
import sys

def rbh_argparser():
    ap = argparse.ArgumentParser()
    ap.add_argument('infile', type=argparse.FileType('r'),
                    default=sys.stdin,
                    help='Infile on format: <acc1> <acc2> <nc-score>')
    return ap


def read_nc_file(datafile):
    hits = []

    for line in datafile:
        acc1, acc2, nc_score_s = line.split()
        nc_score = float(nc_score_s)
        hits.append( (acc1, acc2, nc_score) )
    return hits


def find_rbh(hitlist, epsilon=0):
    '''
    rbh == Reciprocal Best Hits

    If an accession has several "best hits" within epsilon, then all are output.
    '''
    rbh_list = []
    rbh_top = {}                # Save the best NC score for each accession
    sorted_hits = sorted(hitlist, key=lambda elem: -elem[2])

    for acc1, acc2, nc_score in sorted_hits:
        if acc1 in rbh_top:
            if rbh_top[acc1] > nc_score + epsilon:
                continue        # Not a RBH pair, acc1 has a strictly better hit already
        else:
            rbh_top[acc1] = nc_score    

        if acc2 in rbh_top:
            if rbh_top[acc2] > nc_score + epsilon:
                continue        # Not a RBH pair, acc1 has a strictly better hit already
        else:
            rbh_top[acc2] = nc_score

        # If we have come this far, both seqs have each other as "best hit"
        rbh_list.append( ( acc1, acc2, nc_score) )
    return rbh_list
    
    


def rbh_main():
    ap = rbh_argparser()
    args = ap.parse_args()

    data = read_nc_file(args.infile)
    rbh_result = find_rbh(data)
    for acc1, acc2, nc_score in rbh_result:
        print(f'{acc1}\t{acc2}\t{nc_score}')
    
      
