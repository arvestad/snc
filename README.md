# SNC

An attempt at understanding and developing Dannie Durand's Neighborhood Correlation method
for sequence comparisons.

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
