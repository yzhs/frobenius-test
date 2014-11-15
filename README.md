The Quadratic Frobenius Test
============================

This is an implementation of the probabilistic primality test described by Jon
Grantham in "A Probable Prime Test with High Confidence", published in the
"Journal of Number Theory" 72, pages 32--47 (1998).  The algorithm is
implemented in C using the GNU Multiple Precision Arithmetic Library (GMP).

In addition to the Quadratic Frobenius test from the above paper, there is
also an implementation of the Miller Rabin primality test for speed comparison.
Both tests are optimized to roughly the same extant (not very much) to provide
a reasonably fair comparison.

There are also implementations of both of these primality tests using machine
words (assuming a 64-bit architecture).  This implementation is significantly
faster for the numbers it can handle.  Therefore, it can be used to test a
large number of integers for compositeness, or checking a large number of the
possible parameters for a given integer.

The code in this repository was created as part of my bachelor's thesis in
mathematics.


Installation
------------

To build this project, you will need a C99 compatible C and a C++98 compatible
C++ compiler.  I have had success building it with both
[GCC 4.9.2](https://gcc.gnu.org/) and [clang 3.5.0](http://clang.llvm.org/).
You will also need a (not too old version of) GMP (I'm using GMP 6.0.0).

If you want to regenerate the plots, you will need a recent version of the Julia
language, as well as [gnuplot](http://www.gnuplot.info/).  The exact gnuplot
version probably won't matter because my code does not, as far as I know, use
any particularly new features.  It certainly works with version 4.6 patchlevel
6.

In case you want to run the included tests, you will also need `primes` from
djb's [primegen](http://cr.yp.to/primegen.html) project to generate a large
list of primes used during the test.  Alternatively you could just provide your
own file of prime numbers (one line per prime, no empty lines) for
`test/data/primelist.txt`.

To build the code all you have to do is type
```
make
```
assuming all the requirements are installed.  The tests can be built and run
using
```
make test
```

This project does not, at the moment, have any automatic installation method.
As there are only a few, quite specific programmes, it does not seem to be very
useful to install it.  You can obviously just run the programmes from this
directory.


Organization
------------

All files with the `_int` suffix contain the long long version while the
corresponding GMP version is in the file without the `_int` suffix.  Files with
an `_long` suffix contain the GMP implementation and the file without the
suffix contains the long long implementation.  In either case, the file without
the suffix contains some kind of a default version.  For example both versions of
`check_all_params` apply the QFT to the composites `n` with `3<n<10000` and try
all valid parameter pairs `(b,c)` modulo `n`.  As these numbers are all small
enough for the long long version of the QFT, the GMP version of this programme
exists only to make sure the long long implementation is correct and, less
importantly, how they compare from a performance standpoint.

The directory `test` contains tests for the different algorithms and versions.
The subdirectory `test/data` contains lists of primes and composites used by
those tests.  The `plots` directory contains plots to be included in a LaTeX
document (my thesis paper).  The `pictures` directory contains another set of
plots which is less suited for printing, but still pretty nice to look at.

A few of the more important programmes are

* `benchmark` runs the different tests (GMP version) on various inputs
  (primes and composites of the form `2^k+e` for a small `e`, and Mersenne
  numbers, both prime and composite).
* `check_all_params` tries all valid parameter pairs `(b,c)` for all non-square
  composites `n` between `3` and `10000`, counting the number of false
  positives.  There are both long long and GMP versions of this.
* `check_all_small_numbers` deterministically chooses a valid pair `(b,c)` for
  each non-square composite `n`.  It, too, counts the number of false positives.
  There are both long long and GMP versions of this.

There are also a few small helper programmes that compute parameters for the
benchmarks.

* `find_non_smooth_numbers` takes a file called `composites.txt` and reads from
  it a list of integers `k` (one per line).  For each `k` the programme searches
  for the smallest composite `n` greater than `2^k` which does not have a prime
  factor below a certain bound `B`.  Finally a list of those `k`s and `n`s
  (separated by a tab, one pair per line) is written to stdout.  This programme
  was used to generate the current content of `composites.txt`
* `nextprime` takes as input a number `k` and prints the smalles prime `p>2^k`.
  This was used to find the primes listed in `primes.txt`

Additionally, there are a few other programmes and scripts to analyse the
measurements produced by `benchmark` and generate plots from the results.

* `small_composites.sh` is a shell script that reads a log file produced by
  `check_all_params` and fills in the template in `small_composites.template` to
  generate a small snippet of LaTeX about how many numbers were tested and how
  many false positives were found.
* `process_timings.jl` is a script written in Julia which reads the raw data
  (run time for example) produced by `benchmark`.  It computes confidence
  intervalls for the run times, stores data concerning different aspects into
  different files for plotting with gnuplot and produces some plots of its own.
* `plot_algorithms.plt` is a gnuplot sript which plots run time grouped by the
  algorithm tested.
* `plot_sets.plt` plots the run time grouped by input set and
* `plot_multiplications.plt` plots the number of multiplications used by the QFT
  and the time per multiplication.

Each of the gnuplot scripts produces output using the `epslatex` terminal, which
contains the plot itself as an `.eps` file and the labels in a separate `.tex`
file.  In contrast to the plots produced by `process_timings.jl`, which are PNG
files that reside in `pictures/`, gnuplot produces these plots in the `plots`
directory which is linked into the directory for the thesis paper.

The inputs used by `benchmark` are stored in `composites.txt`, `primes.txt`
(both containing numbers of the form `2^k+e` for a given set of exponents `k`
and the smallest positives integers `e` such that `2^k+e` is prime, and a
composite without prime factors `p<B`, respectively), and
`mersenne_numbers.txt` containing Mersenne numbers (mostly but not exclusively
composite) and `mersenne_primes.txt` containing only Mersenne primes.
