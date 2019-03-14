#!/usr/bin/env python3
##############################################################################
##                                                                          ##
##                                PYMPF                                     ##
##                                                                          ##
##              Copyright (C) 2016-2017, Altran UK Limited                  ##
##              Copyright (C) 2018,      Florian Schanda                    ##
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

"""
This module implements IEEE floats using bitvectors as the in-memory
format, but rationals (plus rounding and special case handling) for
calculation. This allows us to directly implement the semantics
described in IEEE-754 (2008) without having to consider any of the
difficult normalisation problems.

The main objective is that the implementation should be simple and
simple to understand and as close to IEEE-754 and SMTLIB as possible
(as opposed to fast) since the main use is the validation of other
IEEE float implementations in SMT solvers (which have to be
fast). It is also not a replacement for libraries such as MPFR
(which is fast, but complicated and does not totally map onto IEEE
floats).
"""

# TODO: Implement RNA in intervals

import random

from .rationals import *
from .interval_q import Interval
from .bitvector import BitVector
from .bisect import Bisect

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

class MPF:
    """Arbitrary precision IEEE-754 floating point number

    *eb* is the number of bits for the exponent.

    *sb* is the number of bits for the significand.

    This includes the "hidden bit", so for example to create a
    single-precision floating point number we use:

    >>> x = MPF(8, 24)

    By default the created float is +0, however *bitvec* can be used
    to set the initial value. The expected format is an integer
    :math:`0 \le bitvec < 2^{eb+sb}` corresponding to the binary
    interchange format. For example to create the single precision
    float representing +1:

    >>> x = MPF(8, 24, 0x3f800000)
    >>> x.to_python_string()
    '1.0'

    """

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

    def __repr__(self):
        return "MPF(%u, %u, 0x%x)" % (self.w, self.p, self.bv)

    ######################################################################
    # Constructors

    def new_mpf(self):
        """Copy constructor

        Returns a new MPF with the same precision and value.
        """
        return MPF(self.w, self.p, self.bv)

    ######################################################################
    # Internal utilities

    def compatible(self, other):
        """Test if another MPF has the same precision"""
        return (self.w == other.w and
                self.p == other.p)

    def unpack(self):
        """Unpack into sign, exponent, and significand

        Returns a tuple of integers (S, E, T) where *S* is the sign
        bit, *E* is the exponent, and *T* is the significand.

        This is the inverse of :func:`pack`.

        """

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
        """Pack sign, exponent, and significand

        Sets the value of the floating point number, such that

        The sign corresponds to *S*, the exponent is *E*, and the
        significand is *T*.

        This is the inverse of :func:`unpack`.

        """
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

    def inf_boundary(self):
        inf = q_pow2(self.emax) * \
              (Rational(2) - Rational(1, 2) * q_pow2(1 - self.p))
        return inf.to_python_int()

    def from_rational(self, rm, q):
        """Convert from rational to MPF

        Sets the value to the nearest representable floating-point
        value described by *q*, rounded according to *rm*.

        """
        assert rm in MPF.ROUNDING_MODES

        inf = Rational(self.inf_boundary())

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
            # print("%08X" % (self.bv, v, v.to_python_float(), low, high))

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
        # print("low : %08X" % low)
        # print("high: %08X" % high)

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
                    self.bv = high
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
        """Convert from MPF to :class:`.Rational`

        Raises AssertionError for infinities or NaN.
        """
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

    def to_int(self, rm):
        """Convert from MPF to Python int`

        Raises AssertionError for infinities and NaN. Returns the
        integer closest to the floating point, rounded according to
        *rm*.

        """
        assert rm in MPF.ROUNDING_MODES
        assert self.isFinite()

        q = q_round(rm, self.to_rational())

        assert q.isIntegral()
        return q.a

    def to_python_float(self):
        """Convert from MPF to Python float`

        Convert to python float. Do not rely on this to be precise for
        now.

        .. todo:: Use sys.float_info to first convert to the correct
          format and then directly interpret the bit-pattern.

        If you want something precise, then use
        :func:`to_python_string` instead.
        """
        if self.isNaN():
            return float("NaN")
        elif self.isInfinite():
            if self.isPositive():
                return float("Infinity")
            else:
                return float("-Infinity")
        else:
            return self.to_rational().to_python_float()

    def to_python_string(self):
        """Convert from MPF to Python string`

        Return a string describing the float. We return a decimal
        string (e.g. '0.25'), except for the following special cases:
        '-0', 'Infinity', '-Infinity', and 'NaN'.

        The decimal string might be very long (but precise), so if you
        want something compact then :func:`to_python_float` might be
        more appropriate.

        """
        if self.isNaN():
            return "NaN"
        elif self.isInfinite():
            if self.isPositive():
                return "Infinity"
            else:
                return "-Infinity"
        elif self.isZero() and self.isNegative():
            return "-0"
        else:
            return self.to_rational().to_decimal_string()

    ######################################################################
    # Setters

    def set_zero(self, sign):
        """Set value to zero with the given sign bit"""
        assert 0 <= sign <= 1
        self.pack(sign, 0, 0)

    def set_infinite(self, sign):
        """Set value to infinite with the given sign bit"""
        assert 0 <= sign <= 1
        self.pack(sign, 2 ** self.w - 1, 0)

    def set_nan(self):
        """Set value to NaN"""
        self.pack(1, 2 ** self.w - 1, 2 ** self.t - 1)

    def set_sign_bit(self, sign):
        """Set sign bit"""
        assert 0 <= sign <= 1
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
            rv += "0x%0*X" % (self.k // 4, self.bv)
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
        """Compute absolute value"""
        rv = self.new_mpf()
        S, E, T = rv.unpack()
        if not rv.isNaN():
            S = 0
            rv.pack(S, E, T)
        return rv

    def __neg__(self):
        """Compute inverse"""
        rv = self.new_mpf()
        S, E, T = rv.unpack()
        if not rv.isNaN():
            S = 1 - S
            rv.pack(S, E, T)
        return rv

    def __le__(self, other):
        """Test floating-point less than or equal"""
        assert self.compatible(other)
        if self.isNaN() or other.isNaN():
            return False
        else:
            return self.partial_order() <= other.partial_order()

    def __lt__(self, other):
        """Test floating-point less than"""
        assert self.compatible(other)
        if self.isNaN() or other.isNaN():
            return False
        else:
            return self.partial_order() < other.partial_order()

    def __ge__(self, other):
        """Test floating-point greater than or equal"""
        assert self.compatible(other)
        if self.isNaN() or other.isNaN():
            return False
        else:
            return self.partial_order() >= other.partial_order()

    def __gt__(self, other):
        """Test floating-point greater than"""
        assert self.compatible(other)
        if self.isNaN() or other.isNaN():
            return False
        else:
            return self.partial_order() > other.partial_order()

    def __eq__(self, other):
        """Test floating-point equality

        Note that this is floating-point equality, so comparisons to
        NaN always fail and -0 and +0 are considered equal.

        If you want true equality, use :func:`smtlib_eq`

        """
        assert self.compatible(other)
        if self.isNaN() or other.isNaN():
            return False
        else:
            return self.partial_order() == other.partial_order()

    ######################################################################
    # Queries

    def isZero(self):
        """Test if value is zero"""
        S, E, T = self.unpack()
        return E == 0 and T == 0

    def isSubnormal(self):
        """Test if value is subnormal"""
        S, E, T = self.unpack()
        return E == 0 and T != 0

    def isNormal(self):
        """Test if value is normal"""
        S, E, T = self.unpack()
        return 1 <= E <= 2 ** self.w - 2

    def isNaN(self):
        """Test if value is not a number"""
        S, E, T = self.unpack()
        return E == 2 ** self.w - 1 and T != 0

    def isInfinite(self):
        """Test if value is infinite"""
        S, E, T = self.unpack()
        return E == 2 ** self.w - 1 and T == 0

    def isPositive(self):
        """Test if value is positive

        Returns always false for NaN.
        """
        S, E, T = self.unpack()
        return not self.isNaN() and S == 0

    def isNegative(self):
        """Test if value is negative

        Returns always false for NaN.
        """
        S, E, T = self.unpack()
        return not self.isNaN() and S == 1

    def isFinite(self):
        """Test if value is finite

        Returns true for zero, subnormals, or normal numbers.

        Returns false for infinities, and not a number.
        """
        return self.isZero() or self.isSubnormal() or self.isNormal()

    def isIntegral(self):
        """Test if value is integral

        Returns true if the floating-point value is an integer, and
        false in all other cases (including infinities and not a
        number).
        """

        if self.isFinite():
            return self.to_rational().isIntegral()
        else:
            return False

    ######################################################################
    # SMTLIB support

    def smtlib_sort(self):
        """SMT-LIB sort

        Returns "Float16", "Float32", "Float64", or "Float128" if
        possible.

        Returns "(_ FloatingPoint *eb* *sb*)" otherwise.

        """
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

    def smtlib_from_binary_interchange(self):
        """SMT-LIB function for converting from bitvector"""
        return "(_ to_fp %u %u)" % (self.w, self.p)

    def smtlib_from_float(self):
        """SMT-LIB function for converting from FloatingPoint"""
        return "(_ to_fp %u %u)" % (self.w, self.p)

    def smtlib_from_real(self):
        """SMT-LIB function for converting from Real"""
        return "(_ to_fp %u %u)" % (self.w, self.p)

    def smtlib_from_int(self):
        """SMT-LIB function for converting from Int"""
        return "(_ to_fp %u %u)" % (self.w, self.p)

    def smtlib_from_ubv(self):
        """SMT-LIB function for converting from unsigned bitvector"""
        return "(_ to_fp_unsigned %u %u)" % (self.w, self.p)

    def smtlib_from_sbv(self):
        """SMT-LIB function for converting from signed bitvector"""
        return "(_ to_fp %u %u)" % (self.w, self.p)

    def smtlib_to_real(self):
        """SMT-LIB function for converting to Real"""
        return "fp.to_real"

    def smtlib_to_int(self):
        """SMT-LIB function for converting to Int"""
        return "fp.to_int"

    def smtlib_literals(self):
        """Return list of all possible SMTLIB literals

        Includes:

        * Binary interchange conversion, e.g. ((_ to_fp 8 24) #x00000000)
        * FP literal, e.g. (fp #b0 #b00 #b000)
        * Special value, e.g. (_ +zero 2 4)
        """
        choices = []

        # Binary interchange
        lit = "((_ to_fp %u %u) " % (self.w, self.p)
        if self.k % 4 == 0:
            # we can use a hex bitvector
            lit += "#x%0*X" % (self.k // 4, self.bv)
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
        """Return an SMT-LIB literal

        Favours special cases, such as (_ +zero 8 24), otherwise
        returns literals of the form (fp ...).

        """
        return self.smtlib_literals()[-1]

    def smtlib_random_literal(self):
        """Return an SMT-LIB literal (randomly chosen)

        Chooses randomly from :func:`smtlib_literals`.
        """
        return random.choice(self.smtlib_literals())

def q_round(rm, n):
    assert rm in MPF.ROUNDING_MODES
    rnd = {RM_RNE : q_round_rne,
           RM_RNA : q_round_rna,
           RM_RTZ : q_round_rtz,
           RM_RTP : q_round_rtp,
           RM_RTN : q_round_rtn}[rm]
    return rnd(n)

def fp_add(rm, left, right):
    """Floating-point addition

    Performs a correctly rounded :math:`left + right`. The following special
    cases apply:

    * NaN propagates
    * Adding opposite infinities results in NaN, otherwise infinities propagate
    * If the precise result is exactly 0 (i.e. 0 + 0 or a + -a) then we
      special case according to Section 6.3 in IEEE-754:

      * we preserve the sign if left and right have the same sign
      * otherwise, we return -0 if the rounding mode is RTN
      * otherwise, we return +0

    """
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
    """Floating-point substraction

    Performs a correctly rounded :math:`left - right`, which is
    equivalent to :math:`left + -right`.

    See :func:`fp_add` for special cases.

    """
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
    """Floating-point multiplication

    Performs a correctly rounded :math:`left * right`. The following
    special cases apply:

    * NaN propagates
    * Infinities propagate (for non-zero operands)
    * Multiplying by zero gives the following results

      * NaN if the other operand is infinite
      * 0 of the same sign as the 0 operand

    """
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
    """Floating-point division

    Performs a correctly rounded :math:`left \div right`. The
    following special cases apply:

    * NaN propagates
    * When dividing by infinity

      * NaN when dividing two infinities
      * 0 otherwise

    * When dividing by zero

      * NaN when dividing zero by zero
      * Infinity otherwise

    """
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
    """Floating-point fused multiply add

    Performs a correctly rounded :math:`x * y + z`. The special cases
    of :func:`fp_mul` and :func:`fp_add` apply, and in addition:

    * On a precise 0 result the sign of zero is:

      * Preserved if the sign of x * y is the same as z
      * -0 for RTN
      * +0 otherwise
    """
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
    """Floating-point square root"""
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
    """Floating-point remainder"""
    assert left.compatible(right)

    rv = left.new_mpf()
    if left.isNaN() or right.isNaN() or left.isInfinite() or right.isZero():
        rv.set_nan()
    elif right.isInfinite():
        # Result is left
        pass
    else:
        n = q_round_rne(left.to_rational() / right.to_rational())
        r = left.to_rational() - (right.to_rational() * n)
        # Rounding mode is irrelevant here, r will be exact
        rv.from_rational(RM_RNE, r)
        if r.isZero():
            rv.set_sign_bit(left.isNegative())

    return rv

def fp_roundToIntegral(rm, op):
    """Floating-point round to integer"""
    assert rm in MPF.ROUNDING_MODES

    rv = op.new_mpf()
    if op.isInfinite() or op.isNaN() or op.isIntegral():
        # Nothing to do here
        pass
    else:
        i = op.to_int(rm)
        rv.from_rational(rm, Rational(i))

        if i == 0 and op.isNegative():
            rv.set_sign_bit(True)

    return rv

def fp_min(left, right):
    """Floating-point minimum"""
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
    """Floating-point maximum"""
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
    """Bit-wise equality"""
    assert left.compatible(right)
    return (left.isNaN() and right.isNaN()) or (left.bv == right.bv)

def fp_nextUp(op):
    """Floating-point successor"""
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
    """Floating-point predecessor"""
    # This is how it is defined in 5.3.1
    return -fp_nextUp(-op)

def fp_from_ubv(eb, sb, rm, op):
    """Conversion from unsigned bitvector to MPF"""
    rv = MPF(eb, sb)
    rv.from_rational(rm, op.to_unsigned_int())
    return rv

def fp_to_ubv(op, rm, width):
    """Conversion from MPF to unsigned bitvector"""
    if op.isInfinite() or op.isNaN():
        raise Unspecified

    bv = BitVector(width)

    i = op.to_int(rm)
    if bv.min_unsigned <= i <= bv.max_unsigned:
        bv.from_unsigned_int(i)
        return bv
    else:
        raise Unspecified

def fp_from_sbv(eb, sb, rm, op):
    """Conversion from signed bitvector to MPF"""
    rv = MPF(eb, sb)
    rv.from_rational(rm, op.to_signed_int())
    return rv

def fp_to_sbv(op, rm, width):
    """Conversion from MPF to signed bitvector"""
    if op.isInfinite() or op.isNaN():
        raise Unspecified

    bv = BitVector(width)

    i = op.to_int(rm)
    if bv.min_signed <= i <= bv.max_signed:
        bv.from_signed_int(i)
        return bv
    else:
        raise Unspecified

# ((_ to_fp eb sb) rm op)
def fp_from_int(eb, sb, rm, op):
    """Conversion from Python integer to MPF"""
    rv = MPF(eb, sb)
    rv.from_rational(rm, Rational(op))
    return rv

# (fp.to_int rm op)
def fp_to_int(rm, op):
    """Conversion from MPF to Python integer"""
    if op.isInfinite() or op.isNaN():
        raise Unspecified

    n = op.to_rational()
    q = q_round(rm, n)
    return q.to_python_int()

# Convert op to (_ FloatingPoint eb sb) under rounding mode rm.
#
# IEEE-754 is a bit vague on what happens to zero in Section 4 (which
# is where you land when you read 5.4.2), but in 6.3 it says it
# doesn't change.
def fp_from_float(eb, sb, rm, op):
    """Conversion from MPF to MPF (of a different precision)"""
    rv = MPF(eb, sb)
    if op.isNaN():
        rv.set_nan()
    elif op.isInfinite():
        rv.set_infinite(int(op.isNegative()))
    elif op.isZero():
        rv.set_zero(int(op.isNegative()))
    else:
        rv.from_rational(rm, op.to_rational())
    return rv

##############################################################################
# Interval stuff
##############################################################################

def interval_nearest(rm, op):
    assert rm in MPF.ROUNDING_MODES_NEAREST
    assert not op.isNaN()

    DEBUG_INTERVAL = False

    if DEBUG_INTERVAL:
        print("Interval query: %s [%s]" % (op, rm))

    interval = Interval()

    op_is_even = (op.bv % 2) == 0
    if DEBUG_INTERVAL:
        if rm == RM_RNE:
            print("> even? : %s" % op_is_even)

    low  = fp_nextDown(op)
    high = fp_nextUp(op)
    if DEBUG_INTERVAL:
        print("> low  : %s" % low)
        print("> high : %s" % high)

    # Boundary for infinity, as described in IEEE 754 (Section 4.3.1)
    #
    # These is a question here as to what the standard really
    # means. It says "an infinitely precise result with magnitude at
    # least [...]"; now what does "at least" mean?
    #
    # I have chosen to interpret this as >=, instead of >.
    inf = Rational(op.inf_boundary())
    if DEBUG_INTERVAL:
        print("> inf  : %s" % inf)

    # Lets establish some basic bounds relevant to round-to-nearest
    #
    #      nextdown(op) ... | ... op ... | ... nextup(op)
    #                     q_low        q_high
    #
    # If the rounding mode is RNE, then q_low and q_high are in the
    # interval if op is even, and not in the interval if op is odd.
    #
    # If the rounding mode is RNA, then q_low is in the interval iff
    # op is positive, and q_high is in the interval iff op is
    # negative.
    #
    # We also need to special case 0 as it breaks intervals up:
    #
    #                   q_low           q_high
    #         nextdown(0) | ... -0 +0 ... | ... nextup(0)
    #
    # RNA -0        q_low ]      [ q_high = 0
    # RNE -0        q_low [      [ q_high = 0
    # RNA +0             q_low = 0 [      [ q_high
    # RNE +0             q_low = 0 [      ] q_high
    #
    # I.e. RNE is the same, except the non-zero low and high intervals
    # are inclusive since 0 is even; this means we do not need to do
    # anything special. The only special case is q_high for -0 and
    # q_low for +0.
    if op.isZero() and op.isNegative():
        q_high = Rational(0)
        high_inclusive = False
    elif not (op.isInfinite() or high.isInfinite()):
        q_high = (op.to_rational() + high.to_rational()) * Rational(1, 2)
        high_inclusive = {RM_RNE : op_is_even,
                          RM_RNA : op.isNegative()}[rm]
    elif not op.isInfinite():
        # We're on the boundary to infinity here
        assert high.isInfinite() and high.isPositive()
        q_high = inf
        high_inclusive = False
    elif op.isPositive():
        assert op.isInfinite()
        # We're +oo, so high bound does not exist
        q_high = None
        high_inclusive = True
    else:
        assert op.isNegative() and op.isInfinite()
        # We're -oo, so high bound is inf
        q_high = -inf
        high_inclusive = True

    if DEBUG_INTERVAL:
        print("> established high bound : %s %s" % (q_high,
                                                    "inclusive"
                                                    if high_inclusive
                                                    else ""))

    # Sanity check that the interval does or does not convert back
    if q_high is not None:
        tmp = op.new_mpf()
        tmp.from_rational(rm, q_high)
        assert smtlib_eq(tmp, op) == high_inclusive

    if op.isZero() and op.isPositive():
        q_low = Rational(0)
        low_inclusive = True
    elif not (op.isInfinite() or low.isInfinite()):
        q_low = (op.to_rational() + low.to_rational()) * Rational(1, 2)
        low_inclusive = {RM_RNE : op_is_even,
                         RM_RNA : op.isPositive()}[rm]
    elif not op.isInfinite():
        # We're on the boundary to infinity here
        assert low.isInfinite() and low.isNegative()
        q_low = -inf
        low_inclusive = False
    elif op.isNegative():
        assert op.isInfinite()
        # We're -oo, so low bound does not exist
        q_low = None
        low_inclusive = True
    else:
        assert op.isPositive() and op.isInfinite()
        # We're +oo, so low bound is inf
        q_low = inf
        low_inclusive = True

    if DEBUG_INTERVAL:
        print("> established low bound  : %s %s" % (q_low,
                                                    "inclusive"
                                                    if low_inclusive
                                                    else ""))

    # Sanity check that the interval does or does not convert back
    if q_low is not None:
        tmp = op.new_mpf()
        tmp.from_rational(rm, q_low)
        assert smtlib_eq(tmp, op) == low_inclusive

    if DEBUG_INTERVAL:
        tmp = ""
        if low_inclusive:
            tmp += "["
        else:
            tmp += "]"
        tmp += " "
        if q_low is None:
            tmp += "-oo"
        else:
            tmp += str(q_low)
        tmp += " .. "
        if q_high is None:
            tmp += "+oo"
        else:
            tmp += str(q_high)
        tmp += " "
        if high_inclusive:
            tmp += "]"
        else:
            tmp += "["
        print("> interval : %s" % tmp)

    if q_low is not None:
        interval.set_low(q_low, low_inclusive)
    if q_high is not None:
        interval.set_high(q_high, high_inclusive)
    if DEBUG_INTERVAL:
        print("> interval : %s" % interval)

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

##############################################################################
# SMTLIB Operations
##############################################################################

TYP_BOOL  = "boolean"
TYP_FLOAT = "float"
TYP_INT   = "int"
TYP_REAL  = "real"
TYP_BV    = "bitvector"

class Floating_Point_Operation(object):
    def __init__(self,
                 name,
                 arity,
                 result_type=TYP_FLOAT,
                 args_type=TYP_FLOAT,
                 rm_arg=True,
                 precision_arg=False):
        self.name          = name
        self.arity         = arity
        self.result_type   = result_type
        self.args_type     = args_type
        self.rm_arg        = rm_arg
        self.precision_arg = precision_arg

class Floating_Point_Predicate(Floating_Point_Operation):
    def __init__(self,
                 name,
                 arity,
                 result_type=TYP_BOOL,
                 args_type=TYP_FLOAT,
                 rm_arg=False,
                 precision_arg=False):
        super(Floating_Point_Predicate, self).__init__(name,
                                                       arity,
                                                       result_type,
                                                       args_type,
                                                       rm_arg,
                                                       precision_arg)

FP_OPS = {
    # Basic operations
    "fp.abs"             : Floating_Point_Operation("fp.abs",
                                                    1, rm_arg=False),
    "fp.neg"             : Floating_Point_Operation("fp.neg",
                                                    1, rm_arg=False),
    "fp.sqrt"            : Floating_Point_Operation("fp.sqrt", 1),
    "fp.roundToIntegral" : Floating_Point_Operation("fp.roundToIntegral", 1),
    "fp.add"             : Floating_Point_Operation("fp.add", 2),
    "fp.sub"             : Floating_Point_Operation("fp.sub", 2),
    "fp.mul"             : Floating_Point_Operation("fp.mul", 2),
    "fp.div"             : Floating_Point_Operation("fp.div", 2),
    "fp.rem"             : Floating_Point_Operation("fp.rem", 2, rm_arg=False),
    "fp.min"             : Floating_Point_Operation("fp.min", 2, rm_arg=False),
    "fp.max"             : Floating_Point_Operation("fp.max", 2, rm_arg=False),
    "fp.fma"             : Floating_Point_Operation("fp.fma", 3),

    # Predicates
    "fp.isNormal"        : Floating_Point_Predicate("fp.isNormal", 1),
    "fp.isSubnormal"     : Floating_Point_Predicate("fp.isSubnormal", 1),
    "fp.isZero"          : Floating_Point_Predicate("fp.isZero", 1),
    "fp.isInfinite"      : Floating_Point_Predicate("fp.isInfinite", 1),
    "fp.isNaN"           : Floating_Point_Predicate("fp.isNaN", 1),
    "fp.isPositive"      : Floating_Point_Predicate("fp.isPositive", 1),
    "fp.isNegative"      : Floating_Point_Predicate("fp.isNegative", 1),
    "fp.eq"              : Floating_Point_Predicate("fp.eq", 2),
    "fp.lt"              : Floating_Point_Predicate("fp.lt", 2),
    "fp.gt"              : Floating_Point_Predicate("fp.gt", 2),
    "fp.leq"             : Floating_Point_Predicate("fp.leq", 2),
    "fp.geq"             : Floating_Point_Predicate("fp.geq", 2),
    "smtlib.eq"          : Floating_Point_Predicate("=", 2),

    # Conversion to float
    "fp.from.real"       : Floating_Point_Operation("to_fp",
                                                    1,
                                                    args_type=TYP_REAL,
                                                    precision_arg=True),
    "fp.from.int"        : Floating_Point_Operation("to_fp",
                                                    1,
                                                    args_type=TYP_INT,
                                                    precision_arg=True),
    "fp.from.ubv"        : Floating_Point_Operation("to_fp_unsigned",
                                                    1,
                                                    args_type=TYP_BV,
                                                    precision_arg=True),
    "fp.from.sbv"        : Floating_Point_Operation("to_fp",
                                                    1,
                                                    args_type=TYP_BV,
                                                    precision_arg=True),
    "fp.from.binary"     : Floating_Point_Operation("to_fp",
                                                    1,
                                                    args_type=TYP_BV,
                                                    rm_arg=False,
                                                    precision_arg=True),

    # Conversion from float to float
    "fp.cast"            : Floating_Point_Operation("to_fp",
                                                    1,
                                                    precision_arg=True),

    # Conversion to other types
    "fp.to.real"         : Floating_Point_Operation("fp.to_real",
                                                    1,
                                                    result_type=TYP_REAL,
                                                    rm_arg=False),
    "fp.to.int"          : Floating_Point_Operation("fp.to_int",
                                                    1,
                                                    result_type=TYP_INT,
                                                    rm_arg=False),
    "fp.to.ubv"          : Floating_Point_Operation("fp.to_ubv",
                                                    1,
                                                    result_type=TYP_BV),
    "fp.to.sbv"          : Floating_Point_Operation("fp.to_sbv",
                                                    1,
                                                    result_type=TYP_BV),

    # Proposed extensions
    "fp.isFinite"        : Floating_Point_Predicate("fp.isFinite", 1),
    "fp.isIntegral"      : Floating_Point_Predicate("fp.isIntegral", 1),
    "fp.nextUp"          : Floating_Point_Operation("fp.nextUp",
                                                    1, rm_arg=False),
    "fp.nextDown"        : Floating_Point_Operation("fp.nextDown",
                                                    1, rm_arg=False),
}
