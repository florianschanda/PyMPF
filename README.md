# PyMPF
This project contains two things:
  - An arbitrary precision IEEE-754 (2008) implementation in Python
  - A random floating point testcase generator for SMTLIB

The main motivation is the test-case generation, so a number of things
might seem a bit strange. This library has a bunch of limitations and
should never be used if you want things to be fast. Depending on what
you want the following libraries are much better:
  - SymFPU (if you want to reason about float in your SMT solver)
  - Z3 MPF (ditto, but has much more dependencies)
  - MPFR (if you just need arbitrary precision FP in your program)

Why "yet another" implementation?
  - This library supports RNA (MPFR does not)
  - This library supports subnormals and infinities (MPFR does, but only with tricks)
  - This library uses IEEE or SMTLIB terminology where possible (unlike MPFR)
  - This library is done completely in Python unlike gmpy, mpmath, etc.
  - This library cross-checks with MPFR when possible (rm is not RNA)
  - This library cross-checks with compiled code when possible (precision is single or double, and rm is not RNA)
  - This library is an independent implementation so can be used to check SMT solvers such as CVC4, Z3, MathSAT, SONOLAR, Colibri, etc.
  - This library uses a stupid but simple algorithm to do rounding

## The testcase generator
This has found real bugs in all solvers we have tried (CVC4, Z3,
MathSAT, Colibri, Alt-Ergo, goSAT, xsat, SONOLAR). I will use this to
submit benchmarks to the SMTLIB competition.

### Result tests
For all the basic operations we select 1, 2 or 3 semi-random inputs
and check if the output matches the expectation (which is computed by
PyMPF). We select not totally randomly, but from all combinations of
the following, for both positive and negative values:
  - Zero, Infinity
  - NaN (only one since SMTLIB doesn't do signalling NaNs)
  - The smallest and largest subnormal
  - The smallest and largest normal
  - A random subnormal
  - A random normal

This catches a surprising number of edge cases.

### Interval tests
Note that this is not supported yet.

For real -> float conversions we first select a float from the same
set as above, and then compute the real valued interval that would
yield the float after rounding.  Similar to above, we generate
multiple tests exporing all interesting points:
  - Something outside that interval, on both sides
  - Something just on the boundary
  - Something inside

# License and Copyright
Everything in this repository is licensed under the GNU GPL v3.  All
code is written by Florian Schanda and is :copyright: Altran UK
Limited and :copyright: Florian Schanda (as documented in each file).
