#!/usr/bin/env python3
##############################################################################
##                                                                          ##
##                                PYMPF                                     ##
##                                                                          ##
##              Copyright (C) 2016-2017, Altran UK Limited                  ##
##              Copyright (C) 2019,      Zenuity AB                         ##
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

class Rational(object):
    def __init__(self, a=0, b=1):
        assert type(a) is int
        assert type(b) is int

        assert b != 0
        self.a = (a if b > 0 else -a)
        self.b = abs(b)
        d = fractions.gcd(self.a, self.b)
        assert d > 0
        self.a = self.a // d
        self.b = self.b // d

        assert type(self.a) is int
        assert type(self.b) is int

    def __mul__(self, other):
        return Rational(self.a*other.a,
                        self.b*other.b)

    def __truediv__(self, other):
        return Rational(self.a*other.b,
                        self.b*other.a)

    def __add__(self, other):
        return Rational(self.a*other.b + other.a*self.b,
                        self.b*other.b)

    def __sub__(self, other):
        return Rational(self.a*other.b - other.a*self.b,
                        self.b*other.b)

    def __pow__(self, other):
        assert other.isIntegral()

        a = fractions.Fraction(self.a, self.b)
        p = a ** other.a

        return Rational(p.numerator, p.denominator)

    def __abs__(self):
        return Rational(abs(self.a), self.b)

    def __neg__(self):
        return Rational(-self.a, self.b)

    def __repr__(self):
        if self.b == 1:
            return "Rational(%i)" % self.a
        else:
            return "Rational(%i, %u)" % (self.a, self.b)

    def __lt__(self, other):
        return self.a*other.b < other.a*self.b

    def __le__(self, other):
        return self.a*other.b <= other.a*self.b

    def __eq__(self, other):
        return self.a*other.b == other.a*self.b

    def __ne__(self, other):
        return self.a*other.b != other.a*self.b

    def __gt__(self, other):
        return self.a*other.b > other.a*self.b

    def __ge__(self, other):
        return self.a*other.b >= other.a*self.b

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
        return float(fractions.Fraction(self.a, self.b))

    def to_smtlib(self):
        tmp = "%u.0" % abs(self.a)
        if self.b != 1:
            tmp = "(/ %s %u.0)" % (tmp, self.b)
        if self.a < 0:
            tmp = "(- %s)" % tmp
        return tmp

    def to_decimal_string(self):
        if self.b > 2:
            b = self.b
            while b > 1:
                if b % 2 == 0:
                    b = b // 2
                elif b % 5 == 0:
                    b = b // 5
                else:
                    raise Exception("decimal for %u / %u will not terminate" %
                                    (self.a, self.b))

        rv = str(abs(self.a) // self.b)

        a = abs(self.a) % self.b
        b = self.b

        if a > 0:
            rv += "."
            while a > 0:
                a *= 10
                rv += str(a // b)
                a = a % b
        else:
            rv += ".0"

        if self.isNegative():
            return "-" + rv
        else:
            return rv

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

def q_from_decimal_fragments(sign, integer_part, fraction_part, exp_part):
    """
    Build a rational from fragments of a decimal number.

    E.g. for "1.23E-1" we have fragments for
       sign          =   ""
       integer_part  =  "1"
       fraction_part = "23"
       exp_part      = "-1"
    """
    assert sign          is None or (type(sign) is str and
                                     len(sign) <= 1)
    assert integer_part  is None or type(integer_part)  is str
    assert fraction_part is None or type(fraction_part) is str
    assert exp_part      is None or type(exp_part)      is str

    if integer_part is None or len(integer_part) == 0:
        q = Rational(0)
    else:
        q = Rational(int(integer_part, 10))

    if fraction_part is not None and len(fraction_part) > 0:
        if fraction_part.startswith("."):
            fraction_part = fraction_part[1:]
        f = Rational(0)
        for digit in reversed(fraction_part):
            f += Rational(int(digit))
            f *= Rational(1, 10)
        q += f

    if exp_part is not None and len(exp_part) > 0:
        if exp_part.startswith("+"):
            exp_part = exp_part[1:]
        if exp_part.startswith("-"):
            exp_part = exp_part[1:]
            q *= Rational(1, 10**int(exp_part, 10))
        else:
            q *= Rational(10**int(exp_part, 10))

    if sign is not None and sign == "-":
        q = -q

    return q
