# PyMPF
This is an arbitrary precision IEEE-754 floating point implementation.

The main motivation is correctness and test-case generation, so a
number of things might seem a bit strange. This library has a bunch of
limitations and should never be used if you want things to be
fast. MPFR, SymFPU or Z3's MPF library are what you want.

Why "yet another" implementation?
  - This library supports RNA (MPFR does not)
  - This library supports subnormals and infinities (MPFR does, but only with
    tricks)
  - This library uses IEEE or SMTLIB terminology where possible (MPFR tends to
    stick to more "maths" terminology)
  - This library is implemented completely in Python unlike gmpy, mpmath, etc.
  - This library is an independent implementation so can be used to check
    Z3/SymFPU
  - This library uses a stupid but simple algorithm to do rounding

The main use of this library is random test-case generation for
SMT-LIB. It has been used to validate the FP implementations of CVC4,
Z3, MathSAT, SONOLAR, Alt-Ergo, Colibri, goSAT, and xsat; and has
found bugs in all of them ;)

# SMT-LIB random testcase generator
This will move to a new project. For now it is still available in
the "python2" branch of PyMPF.

# Requirements
Python 3.5 or later.

There is also a Python2 port, but I don't plan to maintain it
(although do shout if you need it).

# Documentation
This is very much work in progress and entirely incomplete (I just
started adding this).

* [PyMPF API Documentation](https://florianschanda.github.io/PyMPF/)

# Installation
This package is available on PyPI. To install simply run:

```
$ pip3 install PyMPF
```

# License and Copyright
Everything in this repository is licensed under the GNU GPL v3.

Key copyright holders that contributed to this library are:
* Florian Schanda
* Altran UK Limited
* Zenuity AB
