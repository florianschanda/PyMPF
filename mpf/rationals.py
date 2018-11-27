#!/usr/bin/env python
##############################################################################
##                                                                          ##
##                                PYMPF                                     ##
##                                                                          ##
##              Copyright (C) 2016-2017, Altran UK Limited                  ##
##                                                                          ##
##  This file is part of PyMPF.                                             ##
##                                                                          ##
##  PyMPF is free software: you can redistribute it and/or modify           ##
##  it under the terms of the GNU General Public License as published by    ##
##  the Free Software Foundation, either version 3 of the License, or       ##
##  (at your option) any later version.                                     ##
##                                                                          ##
##  PyMPF is distributed in the hope that it will be useful,                ##
##  but WITHOUT ANY WARRANTY; without even the implied warranty of          ##
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           ##
##  GNU General Public License for more details.                            ##
##                                                                          ##
##  You should have received a copy of the GNU General Public License       ##
##  along with PyMPF. If not, see <http://www.gnu.org/licenses/>.           ##
##                                                                          ##
##############################################################################

# Todo: Make this a subclass of Fraction instead of re-inventing the
# wheel.

import fractions
import floats

class Rational(object):
    def __init__(self, a=0, b=1):
        assert b != 0
        self.a = (a if b > 0 else -a)
        self.b = abs(b)
        d = fractions.gcd(self.a, self.b)
        assert d > 0
        self.a = self.a / d
        self.b = self.b / d

    def __mul__(self, other):
        return Rational(self.a*other.a,
                        self.b*other.b)

    def __div__(self, other):
        return Rational(self.a*other.b,
                        self.b*other.a)

    def __add__(self, other):
        return Rational(self.a*other.b + other.a*self.b,
                        self.b*other.b)

    def __sub__(self, other):
        return Rational(self.a*other.b - other.a*self.b,
                        self.b*other.b)

    def __abs__(self):
        return Rational(abs(self.a), self.b)

    def __neg__(self):
        return Rational(-self.a, self.b)

    def __repr__(self):
        if self.b == 1:
            return "Rational(%i)" % self.a
        else:
            return "Rational(%i, %u)" % (self.a, self.b)

    def __cmp__(self, other):
        return cmp(self.a*other.b, other.a*self.b)

    def isZero(self):
        return self.a == 0

    def isNegative(self):
        return self.a < 0

    def isIntegral(self):
        return (self.a % self.b) == 0

    def to_python_int(self):
        assert self.isIntegral()
        return self.a

    def to_python_float(self):
        try:
            return float(self.a) / float(self.b)
        except OverflowError:
            return None

    def to_smtlib(self):
        tmp = "%u.0" % abs(self.a)
        if self.b != 1:
            tmp = "(/ %s %u.0)" % (tmp, self.b)
        if self.a < 0:
            tmp = "(- %s)" % tmp
        return tmp

def q_pow2(n):
    if n >= 0:
        return Rational(2**n)
    else:
        return Rational(1, 2**(-n))

def q_round_to_nearest(n, tiebreak):
    lower = (n.a - (n.a % n.b)) /  n.b
    upper = lower + 1
    assert Rational(lower) <= n <= Rational(upper)

    if n - Rational(lower) < Rational(upper) - n:
        return Rational(lower)
    elif n - Rational(lower) > Rational(upper) - n:
        return Rational(upper)
    else:
        return Rational(tiebreak(lower, upper))

def q_round_rne(n):
    def tiebreak(lower, upper):
        if lower % 2 == 0:
            return lower
        else:
            assert upper % 2 == 0
            return upper
    return q_round_to_nearest(n, tiebreak)

def q_round_rna(n):
    def tiebreak(lower, upper):
        if abs(lower) > abs(upper):
            return lower
        else:
            return upper
    return q_round_to_nearest(n, tiebreak)

def q_round_rtz(n):
    i = Rational(abs(n.a) / n.b)
    if n.isNegative():
        return -i
    else:
        return i

def q_round_rtp(n):
    if n.isIntegral():
        return n
    i = abs(n.a) / n.b
    if n.isNegative():
        return Rational(-i)
    else:
        return Rational(i + 1)

def q_round_rtn(n):
    if n.isIntegral():
        return n
    i = abs(n.a) / n.b
    if n.isNegative():
        return Rational(-i - 1)
    else:
        return Rational(i)

def q_round(rm, n):
    assert rm in floats.MPF.ROUNDING_MODES
    rnd = {floats.RM_RNE : q_round_rne,
           floats.RM_RNA : q_round_rna,
           floats.RM_RTZ : q_round_rtz,
           floats.RM_RTP : q_round_rtp,
           floats.RM_RTN : q_round_rtn}[rm]
    return rnd(n)
