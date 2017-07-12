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

# This module implements IEEE floats using bitvectors as the in-memory
# format, but rationals (plus rounding and special case handling) for
# calculation. This allows us to directly implement the semantics
# described in IEEE-754 (2008) without having to consider any of the
# difficult normalisation problems.
#
# The main objective is that the implementation should be simple and
# simple to understand and as close to IEEE-754 and SMTLIB as possible
# (as opposed to fast) since the main use is the validation of other
# IEEE float implementations in SMT solvers (which have to be
# fast). It is also not a replacement for libraries such as MPFR
# (which is fast, but complicated and does not totally map onto IEEE
# floats).

# TODO: Implement RNA

import random

from rationals import *
from interval_q import Interval
from bisect import Bisect

##############################################################################
# IEEE Floats
##############################################################################

class Not_Implemented_Yet(Exception):
    pass

class Unspecified(Exception):
    pass

RM_RNE = "RNE"
RM_RNA = "RNA"
RM_RTP = "RTP"
RM_RTN = "RTN"
RM_RTZ = "RTZ"

class MPF(object):
    ROUNDING_MODES          = (RM_RNE, RM_RNA, RM_RTP, RM_RTN, RM_RTZ)
    ROUNDING_MODES_NEAREST  = (RM_RNE, RM_RNA)
    ROUNDING_MODES_DIRECTED = (RM_RTP, RM_RTN, RM_RTZ)

    def __init__(self, eb, sb, bitvec=0):
        # Naming chosen so that the interface is close to SMTLIB (eb,
        # sb) and the internals close to IEEE-754 (k, p, w, t, emax,
        # emin).
        assert eb >= 1
        assert sb >= 1
        assert 0 <= bitvec < 2 ** (eb + sb)

        self.k    = eb + sb

        self.p    = sb
        self.w    = self.k - self.p
        self.t    = self.p - 1
        self.emax = (2 ** (self.w - 1)) - 1
        self.emin = 1 - self.emax
        self.bias = self.emax

        self.bv   = bitvec

    ######################################################################
    # Constructors

    def new_mpf(self):
        return MPF(self.w, self.p, self.bv)

    ######################################################################
    # Internal utilities

    def compatible(self, other):
        return (self.w == other.w and
                self.p == other.p)

    def unpack(self):
        # TODO: remove ugly hack to convert via strings
        bits = "{0:b}".format(self.bv)
        assert len(bits) <= self.k
        bits = "0" * (self.k - len(bits)) + bits
        assert len(bits) == self.k

        S = int(bits[0], 2)
        E = int(bits[1:1 + self.w], 2)
        T = int(bits[1 + self.w:], 2)
        return (S, E, T)

    def pack(self, S, E, T):
        assert 0 <= S <= 1
        assert 0 <= E <= 2 ** self.w - 1
        assert 0 <= T <= 2 ** self.t - 1

        self.bv  = S << (self.w + self.t)
        self.bv |= E << self.t
        self.bv |= T

    def partial_order(self):
        # orders all non-NaN floats, with -oo .. {-0, 0} .. +oo
        # i.e. -0 and 0 evaluate to the same integer (0),
        #      -oo evaluates to an integer less than any other,
        #      +oo evaluates to one greater than any other.
        assert not self.isNaN()
        S, E, T = self.unpack()

        rv  = E << self.t
        rv |= T
        if S:
            return -rv
        else:
            return rv

    def from_rational(self, rm, q):
        assert rm in MPF.ROUNDING_MODES

        inf = q_pow2(self.emax) * \
              (Rational(2) - Rational(1, 2) * q_pow2(1 - self.p))

        target = abs(q)
        sign   = q.isNegative()

        # TODO: is this >= or > ???
        #
        # This is an area in the standard that is not 100% clear to
        # many people who have read this section. I chose to read it
        # as >= and so did Martin.
        if rm in MPF.ROUNDING_MODES_NEAREST and target >= inf:
            # infinite result (4.3.1)
            self.set_infinite(sign)
            return

        low  = 0
        high = ((2 ** self.w - 2) << self.t) | (2 ** self.t - 1)

        if target.isZero():
            # Converting 0 always gives +0
            self.set_zero(False)
            return

        # We need to check for infinity first
        self.bv = high
        if self.to_rational() < target:
            if rm == RM_RTP and not sign:
                # +oo
                self.set_infinite(sign)
            elif rm == RM_RTN and sign:
                # -oo
                self.set_infinite(sign)
            else:
                # maximum non-infinite value, with correct sign
                self.set_sign_bit(sign)
            return

        # We implement rounding through binary search as bitvectors
        # views of integers happen to align with the total order. This
        # is slow but easy to understand. We sign-correct at the end.
        for guess in Bisect(low, high):
            self.bv   = guess.value()
            v         = self.to_rational()
            low, high = guess.bounds()
            # print "%08X" % self.bv, v, v.to_python_float(), low, high

            if v > target:
                guess.too_high()
            elif v < target:
                guess.too_low()
            else:
                # Exact non-zero result, just need to fix up the sign
                # and we're good to go.
                self.set_sign_bit(sign)
                return
        assert low + 1 == high or low == high
        # print "low : %08X" % low
        # print "high: %08X" % high

        # Binary search has given us a lower and upper bound, now we
        # select based on rounding mode.
        self.bv = low
        v_low = self.to_rational()

        self.bv = high
        v_high = self.to_rational()

        assert v_low <= target <= v_high
        halfpoint = (v_low + v_high) * Rational(1, 2)
        assert v_low <= halfpoint <= v_high

        if rm in MPF.ROUNDING_MODES_NEAREST:
            # Implement rounding for 4.3.1 (nearest)
            if target == halfpoint:
                # exactly between, now we need to follow 4.3.1 or 4.3.2
                if rm == RM_RNE:
                    if high % 2 == 0:
                        self.bv = high
                    else:
                        self.bv = low
                    assert self.bv % 2 == 0
                else:
                    assert rm == RM_RNA
                    if high % 2 == 1:
                        self.bv = high
                    else:
                        self.bv = low
                    assert self.bv % 2 == 1
            elif target < halfpoint:
                self.bv = low
            else:
                assert target > halfpoint
                self.bv = high
        else:
            assert rm in MPF.ROUNDING_MODES_DIRECTED
            # Implement rounding for 4.3.2 (directed)

            if rm == RM_RTP:
                if sign:
                    self.bv = low
                else:
                    self.bv = high
            elif rm == RM_RTN:
                if sign:
                    self.bv = high
                else:
                    self.bv = low
            else:
                assert rm == RM_RTZ
                self.bv = low

        # And just need to fix up the sign
        self.set_sign_bit(sign)

    def to_rational(self):
        S, E, T = self.unpack()

        if E == 2 ** self.w - 1:
            # Infinity (T = 0) or NaN (T != 0)
            assert False
        elif 1 <= E or T != 0:
            v = q_pow2(1 - self.p) * Rational(T)
            if 1 <= E:
                # normal -1^S * S^(E-bias) * (1 + 2^(1-p) * T)
                assert E <= 2 ** self.w - 2
                v = q_pow2(E - self.bias) * (Rational(1) + v)
            else:
                assert E == 0 and T != 0
                # subnormal -1^S * 2^emin * (0 + 2^(1-p) * T)
                v = q_pow2(self.emin) * v
            if S == 1:
                v = -v
        else:
            assert E == 0 and T == 0
            # zero
            v = Rational(0)

        return v

    ######################################################################
    # Setters

    def set_zero(self, sign):
        assert 0 <= sign <= 1
        self.pack(sign, 0, 0)

    def set_infinite(self, sign):
        assert 0 <= sign <= 1
        self.pack(sign, 2 ** self.w - 1, 0)

    def set_nan(self):
        self.pack(1, 2 ** self.w - 1, 2 ** self.t - 1)

    def set_sign_bit(self, sign):
        assert type(sign) is bool

        S, E, T = self.unpack()
        S = (1 if sign else 0)
        self.pack(S, E, T)

    ######################################################################
    # Built-ins

    def __str__(self):
        rv = "Float%u(" % self.k

        S, E, T = self.unpack()

        if E == 2 ** self.w - 1:
            # Infinity (T = 0) or NaN (T != 0)
            if T == 0:
                if S:
                    rv += "-"
                else:
                    rv += "+"
                rv += "oo"
            else:
                rv += "NaN"
        elif 1 <= E or T != 0:
            rv += "0x%0*X" % (self.k / 4, self.bv)
            rv += " "
            rv += "[%s, %f]" % (self.to_rational(),
                                self.to_rational().to_python_float())
        else:
            if S:
                rv += "-"
            else:
                rv += "+"
            rv += "zero"

        rv += ")"
        return rv

    def __abs__(self):
        rv = self.new_mpf()
        S, E, T = rv.unpack()
        if not rv.isNaN():
            S = 0
            rv.pack(S, E, T)
        return rv

    def __neg__(self):
        rv = self.new_mpf()
        S, E, T = rv.unpack()
        if not rv.isNaN():
            S = 1 - S
            rv.pack(S, E, T)
        return rv

    def __le__(self, other):
        assert self.compatible(other)
        if self.isNaN() or other.isNaN():
            return False
        else:
            return self.partial_order() <= other.partial_order()

    def __lt__(self, other):
        assert self.compatible(other)
        if self.isNaN() or other.isNaN():
            return False
        else:
            return self.partial_order() < other.partial_order()

    def __ge__(self, other):
        assert self.compatible(other)
        if self.isNaN() or other.isNaN():
            return False
        else:
            return self.partial_order() >= other.partial_order()

    def __gt__(self, other):
        assert self.compatible(other)
        if self.isNaN() or other.isNaN():
            return False
        else:
            return self.partial_order() > other.partial_order()

    def __eq__(self, other):
        assert self.compatible(other)
        if self.isNaN() or other.isNaN():
            return False
        else:
            return self.partial_order() == other.partial_order()

    ######################################################################
    # Queries

    def isZero(self):
        S, E, T = self.unpack()
        return E == 0 and T == 0

    def isSubnormal(self):
        S, E, T = self.unpack()
        return E == 0 and T != 0

    def isNormal(self):
        S, E, T = self.unpack()
        return 1 <= E <= 2 ** self.w - 2

    def isNaN(self):
        S, E, T = self.unpack()
        return E == 2 ** self.w - 1 and T != 0

    def isInfinite(self):
        S, E, T = self.unpack()
        return E == 2 ** self.w - 1 and T == 0

    def isPositive(self):
        S, E, T = self.unpack()
        return not self.isNaN() and S == 0

    def isNegative(self):
        S, E, T = self.unpack()
        return not self.isNaN() and S == 1

    def isFinite(self):
        return self.isZero() or self.isSubnormal() or self.isNormal()

    def isIntegral(self):
        if self.isFinite():
            return self.to_rational().isIntegral()
        else:
            return False

    ######################################################################
    # SMTLIB support

    def smtlib_sort(self):
        if self.w == 5 and self.p == 11:
            return "Float16"
        elif self.w == 8 and self.p == 24:
            return "Float32"
        elif self.w == 11 and self.p == 53:
            return "Float64"
        elif self.w == 15 and self.p == 113:
            return "Float128"
        else:
            return "(_ FloatingPoint %u %u)" % (self.w, self.p)

    def smtlib_from_real(self):
        return "(_ to_fp %u %u)" % (self.w, self.p)

    def smtlib_from_ubv(self):
        return "(_ to_fp_unsigned %u %u)" % (self.w, self.p)

    def smtlib_from_sbv(self):
        return "(_ to_fp %u %u)" % (self.w, self.p)

    def smtlib_literals(self):
        choices = []

        # Binary interchange
        lit = "((_ to_fp %u %u) " % (self.w, self.p)
        if self.k % 4 == 0:
            # we can use a hex bitvector
            lit += "#x%0*X" % (self.k / 4, self.bv)
        else:
            # we need to use a binary bitvector
            lit += ("#b{0:0%ub}" % self.k).format(self.bv)
        lit += ")"
        choices.append(lit)

        # FP literal
        S, E, T = self.unpack()
        lit = ("(fp #b{0:0b} #b{1:0%ub} #b{2:0%ub})" % (self.w,
                                                        self.t)).format(S, E, T)
        choices.append(lit)

        # Special values
        lit = "(_ "
        if self.isZero():
            if self.isPositive():
                lit += "+"
            else:
                lit += "-"
            lit += "zero %u %u)" % (self.w, self.p)
            choices.append(lit)
        elif self.isInfinite():
            if self.isPositive():
                lit += "+"
            else:
                lit += "-"
            lit += "oo %u %u)" % (self.w, self.p)
            choices.append(lit)
        elif self.isNaN():
            lit += "NaN %u %u)" % (self.w, self.p)
            choices.append(lit)

        return choices

    def smtlib_literal(self):
        return self.smtlib_literals()[-1]

    def smtlib_random_literal(self):
        return random.choice(self.smtlib_literals())

def fp_add(rm, left, right):
    assert rm in MPF.ROUNDING_MODES
    assert left.compatible(right)
    rv = left.new_mpf() # rv == left
    if left.isNaN() or right.isNaN():
        rv.set_nan()
    elif left.isInfinite() and right.isInfinite():
        if left.isPositive() != right.isPositive():
            # -oo + +oo, +oo + -oo is NaN
            rv.set_nan()
    elif left.isInfinite():
        rv.bv = left.bv
    elif right.isInfinite():
        rv.bv = right.bv
    else:
        q = left.to_rational() + right.to_rational()
        if q.isZero():
            # This is implementing 6.3
            if left.isPositive() == right.isPositive():
                # result exactly zero, and same sign of operands; preserve sign
                rv.bv = left.bv
            elif rm == RM_RTN:
                # result is zero, signs differ, so -0 for RTN
                rv.set_zero(1)
            else:
                # or +0 otherwise
                rv.set_zero(0)
        else:
            # otherwise just round as normal
            rv.from_rational(rm, q)

    return rv

def fp_sub(rm, left, right):
    assert rm in MPF.ROUNDING_MODES
    assert left.compatible(right)
    rv = left.new_mpf() # rv == left
    if left.isNaN() or right.isNaN():
        rv.set_nan()
    elif left.isInfinite() and right.isInfinite():
        if left.isPositive() == right.isPositive():
            # -oo - -oo, +oo - +oo is NaN
            rv.set_nan()
    elif left.isInfinite():
        rv.bv = left.bv
    elif right.isInfinite():
        rv = -(right)
    else:
        q = left.to_rational() - right.to_rational()
        if q.isZero():
            # This is implementing 6.3
            if left.isPositive() != right.isPositive():
                # result exactly zero with different signs, preserve
                rv.bv = left.bv
            elif rm == RM_RTN:
                # result is zero, signs differ, so -0 for RTN
                rv.set_zero(1)
            else:
                # or +0 otherwise
                rv.set_zero(0)
        else:
            # otherwise just round as normal
            rv.from_rational(rm, q)

    return rv

def fp_mul(rm, left, right):
    assert rm in MPF.ROUNDING_MODES
    assert left.compatible(right)
    sign = (1 if left.isNegative() ^ right.isNegative() else 0)

    rv = left.new_mpf()
    if left.isNaN() or right.isNaN():
        rv.set_nan()
    elif left.isInfinite() or right.isInfinite():
        if left.isZero() or right.isZero():
            rv.set_nan()
        else:
            rv.set_infinite(sign)
    elif left.isZero() or right.isZero():
        rv.set_zero(sign)
    else:
        q = left.to_rational() * right.to_rational()
        rv.from_rational(rm, q)

    return rv

def fp_div(rm, left, right):
    assert rm in MPF.ROUNDING_MODES
    assert left.compatible(right)
    sign = (1 if left.isNegative() ^ right.isNegative() else 0)

    rv = left.new_mpf()
    if left.isNaN() or right.isNaN():
        rv.set_nan()
    elif left.isInfinite() and right.isInfinite():
        rv.set_nan()
    elif left.isZero() and right.isZero():
        rv.set_nan()
    elif left.isInfinite() or right.isZero():
        rv.set_infinite(sign)
    elif right.isInfinite():
        rv.set_zero(sign)
    else:
        q = left.to_rational() / right.to_rational()
        if q.isZero():
            rv.set_zero(sign)
        else:
            rv.from_rational(rm, q)

    return rv

def fp_fma(rm, x, y, z):
    # computes (x * y) + z
    assert rm in MPF.ROUNDING_MODES
    assert x.compatible(y)
    assert x.compatible(z)

    sign_xy = (1 if x.isNegative() ^ y.isNegative() else 0)
    sign_z  = (1 if z.isNegative() else 0)

    rv = x.new_mpf()
    if x.isNaN() or y.isNaN() or z.isNaN():
        rv.set_nan()
    elif x.isInfinite() or y.isInfinite():
        if x.isZero() or y.isZero():
            rv.set_nan()
        elif z.isInfinite() and sign_xy != z.isNegative():
            rv.set_nan()
        else:
            rv.set_infinite(sign_xy)
    elif z.isInfinite():
        rv.set_infinite(z.isNegative())
    else:
        q = (x.to_rational() * y.to_rational()) + z.to_rational()
        if q.isZero():
            # This is implementing 6.3
            if sign_xy == sign_z:
                # result exactly zero, and same sign of operands; preserve sign
                rv.set_zero(sign_xy)
            elif rm == RM_RTN:
                # result is zero, signs differ, so -0 for RTN
                rv.set_zero(1)
            else:
                # or +0 otherwise
                rv.set_zero(0)
        else:
            # otherwise just round as normal
            rv.from_rational(rm, q)

    return rv

def fp_sqrt(rm, op):
    assert rm in MPF.ROUNDING_MODES
    # We can get away with approximating the square root to sufficient
    # precision because of Theorem 19 (p168) of the Handbook of
    # Floating-Point Arithmetic: "In any radix, the square root of a
    # floating-point number cannot be the exact midpoint between two
    # consecutive floating-point numbers."

    root = op.new_mpf()
    if op.isNaN() or (op.isNegative() and not op.isZero()):
        root.set_nan()
    elif op.isInfinite() or op.isZero():
        pass # OK as is, preserve sign of zero
    else:
        # We now do a binary search to find the closest binary float. This is
        # much faster than trying to approximate in rationals as the numbers
        # can very quickly get very big.
        low  = op.new_mpf()
        high = op.new_mpf()

        low.bv  = 0
        high.bv = ((2 ** op.w - 2) << op.t) | (2 ** op.t - 1)

        target = op.to_rational()

        for guess in Bisect(low.bv, high.bv):
            root.bv         = guess.value()
            low.bv, high.bv = guess.bounds()
            q               = root.to_rational()
            sq              = q * q

            if sq > target:
                guess.too_high()
            elif sq < target:
                guess.too_low()
            else:
                # We've found an exact match, we should return that.
                return root

        if rm in MPF.ROUNDING_MODES_NEAREST:
            low_error = target - (low.to_rational() * low.to_rational())
            high_error = (high.to_rational() * high.to_rational()) - target
            if low_error < high_error:
                root.bv = low.bv
            else:
                assert high_error < low_error
                root.bv = high.bv
        elif rm in (RM_RTN, RM_RTZ):
            root.bv = low.bv
        else:
            root.bv = high.bv

    return root

def fp_rem(left, right):
    assert left.compatible(right)

    rv = left.new_mpf()
    if left.isNaN() or right.isNaN() or left.isInfinite() or right.isZero():
        rv.set_nan()
    elif right.isInfinite():
        # Result is left
        pass
    else:
        n = q_round_even(left.to_rational() / right.to_rational())
        r = left.to_rational() - (right.to_rational() * n)
        # Rounding mode is irrelevant here, r will be exact
        rv.from_rational(RM_RNE, r)
        if r.isZero():
            rv.set_sign_bit(left.isNegative())

    return rv

def fp_roundToIntegral(rm, op):
    assert rm in MPF.ROUNDING_MODES

    rv = op.new_mpf()
    if op.isZero() or op.isInfinite() or op.isNaN() or op.isIntegral():
        # Nothing to do here
        pass
    else:
        target = op.to_rational()
        assert abs(target) <= q_pow2(op.p)

        # Find FP integral just below and above
        high = 2 ** op.p
        low  = -high

        for guess in Bisect(low, high):
            v = guess.value()

            if Rational(v) > target:
                guess.too_high()
                high = v
            elif Rational(v) < target:
                guess.too_low()
                low = v
            else:
                assert False

        if rm in MPF.ROUNDING_MODES_NEAREST:
            if target - Rational(low) < Rational(high) - target:
                rv.from_rational(rm, Rational(low))
            elif target - Rational(low) > Rational(high) - target:
                rv.from_rational(rm, Rational(high))
            else:
                rv.from_rational(rm, Rational(low))
                if rm == RM_RNE:
                    if rv.bv % 2 == 1:
                        rv.from_rational(rm, Rational(high))
                    assert rv.bv % 2 == 0
                else:
                    assert rm == RM_RNA
                    if rv.bv % 2 == 0:
                        rv.from_rational(rm, Rational(high))
                    assert rv.bv % 2 == 1
        elif rm == RM_RTP:
            rv.from_rational(rm, Rational(high))
        elif rm == RM_RTN:
            rv.from_rational(rm, Rational(low))
        else:
            assert rm == RM_RTZ
            if abs(low) < abs(high):
                rv.from_rational(rm, Rational(low))
            else:
                rv.from_rational(rm, Rational(high))

        # We always preserve the sign in this operation (see 5.9)
        rv.set_sign_bit(target.isNegative())

    return rv

def fp_min(left, right):
    assert left.compatible(right)
    if (left.isZero() and right.isZero() and
        left.isPositive() != right.isPositive()):
        raise Unspecified

    if left.isNaN():
        return right.new_mpf()
    elif left > right:
        return right.new_mpf()
    else:
        return left.new_mpf()

def fp_max(left, right):
    assert left.compatible(right)
    if (left.isZero() and right.isZero() and
        left.isPositive() != right.isPositive()):
        raise Unspecified

    if left.isNaN():
        return right.new_mpf()
    elif left < right:
        return right.new_mpf()
    else:
        return left.new_mpf()

def smtlib_eq(left, right):
    assert left.compatible(right)
    return (left.isNaN() and right.isNaN()) or (left.bv == right.bv)

def fp_nextUp(op):
    rv = op.new_mpf()

    if op.isNaN():
        rv.set_nan()
    elif op.isInfinite():
        if op.isNegative():
            # Largest negative normal
            S = 1
            E = 2 ** rv.w - 2
            T = 2 ** rv.t - 1
            rv.pack(S, E, T)
        else:
            assert op.isPositive()
            rv.set_infinite(0)
    elif op.isPositive() or op.isZero():
        S, E, T = op.unpack()
        S = 0
        T += 1
        if T > 2 ** rv.t - 1:
            T = 0
            E += 1
        rv.pack(S, E, T)
    else:
        assert op.isNegative()
        S, E, T = op.unpack()
        T -= 1
        if T < 0:
            T = 2 ** rv.t - 1
            E -= 1
        rv.pack(S, E, T)

    return rv

def fp_nextDown(op):
    # This is how it is defined in 5.3.1
    return -fp_nextUp(-op)



def interval_nearest(rm, op):
    assert rm in MPF.ROUNDING_MODES_NEAREST
    assert not op.isNaN()

    interval = Interval()

    low  = fp_nextDown(op)
    high = fp_nextUp(op)

    inf = q_pow2(op.emax) * \
          (Rational(2) - Rational(1, 2) * q_pow2(1 - op.p))

    favoured_even = op.bv % 2 == 0
    # TODO: RNA should round to the largest magnitude

    #
    # -oo ... -inf ... -maxnormal ... -0 +0 ... +maxnormal ... +inf ... +oo
    #

    if op.isInfinite():
        # Eliminate the infinite case: values at least inf round to
        # infinity, so anything less will not.
        if op.isNegative():
            interval.set_high(-inf, inclusive=True)
        else:
            interval.set_low(inf, inclusive=True)
    elif low.isInfinite() or high.isInfinite():
        # We're just on the boundary, so we go from inf (exclusive) to the
        # relevant value below.
        if low.isInfinite():
            q = (op.to_rational() + high.to_rational()) * Rational(1, 2)
            interval.set_low(-inf, inclusive=False)
            interval.set_high(q, inclusive=favoured)
        else:
            q = (op.to_rational() + low.to_rational()) * Rational(1, 2)
            interval.set_low(q, inclusive=favoured)
            interval.set_high(inf, inclusive=False)
    elif op.isZero():
        # Zero is again a bit special as we need to carefully deal with
        # underflow.
        q = high.to_rational() * Rational(1, 2)
        if op.isNegative():
            interval.set_low(-q, inclusive=favoured)
            # 0 never rounds to -0
            interval.set_high(Rational(0), inclusive=False)
        else:
            # 0 always rounds to +0
            interval.set_low(Rational(0), inclusive=True)
            interval.set_high(q, inclusive=favoured)
    else:
        assert low.isNegative() == op.isNegative()
        assert op.isNegative() == high.isNegative()
        # We're either side of zero now, and we don't have to worry about
        # infinities.
        q = (op.to_rational() + low.to_rational()) * Rational(1, 2)
        interval.set_low(q, inclusive=favoured)

        q = (op.to_rational() + high.to_rational()) * Rational(1, 2)
        interval.set_high(q, inclusive=favoured)

    return interval

def interval_up(rm, op):
    assert rm in (RM_RTP, RM_RTZ)
    assert not op.isNaN()
    assert rm == RM_RTP or (op.isNaN() or op.isNegative())

    interval = Interval()

    low = fp_nextDown(op)

    if op.isZero():
        # Zero is a bit special of course.
        if op.isNegative():
            interval.set_low(low.to_rational(), inclusive=False)
            interval.set_high(Rational(0), inclusive=False)
        else:
            interval.set_low(Rational(0), inclusive=True)
            interval.set_high(Rational(0), inclusive=True)
    elif op.isInfinite():
        if op.isNegative():
            # RTP never rounds to -infinity
            return None
        else:
            interval.set_low(low.to_rational(), inclusive=False)
    else:
        if not low.isInfinite():
            interval.set_low(low.to_rational(), inclusive=False)
        interval.set_high(op.to_rational(), inclusive=True)

    return interval

def interval_down(rm, op):
    assert rm in (RM_RTN, RM_RTZ)
    assert not op.isNaN()
    assert rm == RM_RTN or (op.isNaN() or op.isPositive())

    interval = Interval()

    high = fp_nextUp(op)

    if op.isZero():
        # Zero is a bit special of course.
        if op.isNegative():
            return None
            # interval.set_low(Rational(0), inclusive=True)
            # interval.set_high(Rational(0), inclusive=True)
        else:
            interval.set_low(Rational(0), inclusive=True)
            interval.set_high(high.to_rational(), inclusive=False)

    elif op.isInfinite():
        if op.isNegative():
            interval.set_high(high.to_rational(), inclusive=False)
        else:
            # RTN never rounds to +infinity
            return None
    else:
        interval.set_low(op.to_rational(), inclusive=True)
        if not high.isInfinite():
            interval.set_high(high.to_rational(), inclusive=False)

    return interval

def fp_interval(rm, op):
    assert rm in MPF.ROUNDING_MODES
    assert not op.isNaN()

    return {
        RM_RNE: interval_nearest,
        RM_RNA: interval_nearest,
        RM_RTP: interval_up,
        RM_RTN: interval_down,
        RM_RTZ: (interval_up if op.isNegative() else interval_down),
    }[rm](rm, op)
