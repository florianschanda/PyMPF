=====
PyMPF
=====

.. toctree::
   :maxdepth: 2

The PyMPF library contains two modules intended for most users:

* floats (arbitrary precision IEEE-754 floating point, see
  :class:`mpf.floats.MPF`)
* rationals (rational numbers, see :class:`mpf.rationals.Rational`)

It also contains the following modules indended for internal use and
the SMT-LIB testcase generator:

* bitvectors (very simple bitvector support, only literal printing currently)
* interval_q (rational intervals)
* bisect (binary search)

Fast tutorial
-------------

Import the relevant classes. Most likely you want to do this:

>>> from mpf.rationals import Rational
>>> from mpf.floats import *

You can now create a float like this (here we create a single
precision float):

>>> x = MPF(8, 24)

To quickly see what we have we can use the
:func:`mpf.floats.MPF.to_python_string` member function.

>>> x.to_python_string()
'0.0'

To set the float to a specific value, such as :math:`\frac{1}{3}` we
can use the :func:`mpf.floats.MPF.from_rational` member
function. Since we convert from rationals to floats we might need to
round.

>>> x.from_rational(RM_RNE, Rational(1, 3))
>>> x.to_python_string()
'0.3333333432674407958984375'

PyMPF supports all rounding modes defined by IEEE-754:

* RM_RNE (Round nearest even: to break ties, round to the nearest
  floating point number whos bit-pattern is even. This is the default
  on most desktop processors and programming languages.)

* RM_RNA (Round nearest away: to break ties, round to the floating
  point number furthest away from zero. Note that this is unsupported
  on most hardware (including i686 and amd64), other floating point
  libraries (e.g. MPFR).

* RM_RTZ (Round to zero: always round towards zero)

* RM_RTP (Round to positive: always round towards :math:`+\infty`)

* RM_RTN (Round to negative: always round towards :math:`-\infty`)

One of the main use-cases for this library is to generate test-cases
for SMT-LIB. To create an SMT-LIB literal you can use the
:func:`mpf.floats.MPF.smtlib_literal` function:

>>> x.smtlib_literal()
'(fp #b0 #b01111101 #b01010101010101010101011)'

The MPF class supports all floating-point comparisons:

>>> y = MPF(8, 24)
>>> x > y
True

>>> y.set_nan()
>>> x > y
False

Note that equality considers +0 and -0 to be equal. You can use the
:func:`mpf.floats.smtlib_eq` if you want bitwise equality:

>>> z = MPF(8, 24)
>>> z.set_zero(1)
>>> y.set_zero(0)
>>> y == z
True
>>> smtlib_eq(y, z)
False

To set values you can use the following functions:

* :func:`mpf.floats.MPF.from_rational` set to value closest to given rational

* :func:`mpf.floats.MPF.set_zero` set to +0 (if sign is 0) or -0 (if
  sign is 1)

* :func:`mpf.floats.MPF.set_infinite` set to :math:`+\infty` (if sign
  is 0) or :math:`-\infty` (if sign is 1)

* :func:`mpf.floats.MPF.set_nan` set to NaN. PyMPF does not support
  the distinction between signalling and non-signalling NaNs,
  similarly to SMT-LIB.

Finally, to do arithmetic, you can use the fp_* functions. Most take a
rounding mode and two parameters:

>>> x.from_rational(RM_RNE, Rational(1, 10))
>>> y.from_rational(RM_RNE, Rational(10))
>>> z = fp_mul(RM_RNE, x, y)
>>> z.to_python_string()
'1.0'

Here an example demonstrating accumulated rounding errors:

>>> y.set_zero(0)
>>> for i in range(10):
>>>     y = fp_add(RM_RNE, y, x)
>>>     print(y.to_python_string())
0.100000001490116119384765625
0.20000000298023223876953125
0.300000011920928955078125
0.4000000059604644775390625
0.5
0.60000002384185791015625
0.7000000476837158203125
0.80000007152557373046875
0.900000095367431640625
1.00000011920928955078125

=========
Rationals
=========

.. automodule:: mpf.rationals
   :members:
   :special-members:

===
MPF
===

.. automodule:: mpf.floats
   :members:
   :special-members:

=========
Changelog
=========

1.0
---

1.0.4
^^^^^
* We now run PyLint

1.0.3
^^^^^
* Add basic documentation (fast tutorial, and basic descriptions of
  MPF and Rational classes).

1.0.2
^^^^^
* First public release on PyPI.
