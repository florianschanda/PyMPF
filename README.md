# PyMPF
This is an arbitrary precision IEEE-754 floating point implementation.

The main motivation is correctness and test-case generation, so a
number of things might seem a bit strange. This library has a bunch of
limitations and should never be used if you want things to be
fast. MPFR, SymFPU or Z3's MPF library are what you want.

Why "yet another" implementation?
  - This library supports RNA (MPFR does not)
  - This library supports subnormals and infinities (MPFR does, but only with tricks)
  - This library uses IEEE or SMTLIB terminology where possible (unlike MPFR)
  - This library is done completely in Python unlike gmpy, mpmath, etc.
  - This library is an independent implementation so can be used to check Z3/CVC4
  - This library uses a stupid but simple algorithm to do rounding

The main use of this library is random test-case generation for
SMT-LIB. It has been used to validate the FP implementations of CVC4,
Z3, MathSAT, SONOLAR, Alt-Ergo, Colibri, goSAT, and xsat; and has
found bugs in all of them ;)

# License and Copyright
Everything in this repository is licensed under the GNU GPL v3.

All code is written by Florian Schanda and is :copyright: Altran UK
Limited and :copyright: Zenuity AB.
