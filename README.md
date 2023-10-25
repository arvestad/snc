# SNC

An attempt at understanding and developing Dannie Durand's Neighborhood Correlation method
(see references below) for sequence comparisons.

There are several ideas tested here. 

1. Sparse matrices, from SciPy, are used for storing similarity scores. Since most comparisons do
   not have interesting similarity, we assume there are many zeros in a similarity matrix. This is
   deviating from an approach in the original Neighborhood Correlation papers of using as many
   pairwise scores as possible, including statistically insignificant scores, as input.
2. We are not assuming all-against-all comparisons. Rather, the input sequences we are interestedin,
   call them Q, are assumed to be compared to a reference database R, which may or may not overlap
   with Q. In fact, you can have R=Q, but we want to enable the use of a recurrently used reference
   database, and the goal is to reduce the number of sequence comparisons.
3. `NC_standalone` used a log-transform on E values instead of bitscores. The transform gave
   improved results. We are trying several other transforms, including log-transform on "bitscore +
   1" to avoid complications when the score/E value is close to zero.
   
The code for `snc` is way easier to understand and modify than the code for `NC_standalone` (well, I
am biased), but `snc` is also far slower.



## Install

For experimentation:

* Download this repo
* Ensure you have Python 3.8(?) and SciPy version 1.8.0 or later installed.
* Install the python module `flit` (e.g., `pip install flit`)
* Run `flit install`

You should then be able to run `snc similarities.tab` on a file produced by Diamond or BLAST
using the option `--outfmt 6`.

## Using snc

Notice the options!

* `-c x`: set the "consideration threshold" to _x_. If two seqeuences _a_ and _b_ have a similarity score being
  _x_ or higher, then NC will be calculated for _a_ and _b_. This is also transitive, so if similarity between _b_ and _c_ 
  is larger than _x_, but _a_ and _c_ is lower than _x_, then we will still compute NC for _a_ and _c_ thanks to the link 
  via _b_. This has a big effect on runtime! Default is set to 30, which is pretty low, in my opinion, but I have not
  benchmarked this at all.
* `-a`: Don't bother, it is not implemented yet. The idea was to mimic NC_standalone with this option, but the NC statistics 
  has not been adapted to that option.
* `-t`: This threshold is purely on output, so can be used as a means of reducing output space.


## References

* JM Joseph and D Durand, _Family classification without domain chaining_, Bioinformatics, 2009.
* N Song, R Sedgewick, D Duranda, _Domain architecture comparison for multidomain homology identification_, J Comp Biol, 2007.
