#!/usr/bin/env python3
##############################################################################
##                                                                          ##
##                                PYMPF                                     ##
##                                                                          ##
##              Copyright (C) 2017, Altran UK Limited                       ##
##              Copyright (C) 2018, Florian Schanda                         ##
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

import random

class BitVector(object):
    def __init__(self, width):
        self.width = width
        self.bv = [0] * width
        self.min_unsigned = 0
        self.max_unsigned = 2 ** width - 1
        self.min_signed = - (2 ** (width - 1))
        self.max_signed = 2 ** (width - 1) - 1

    def __str__(self):
        return self.smtlib_literal()

    ######################################################################
    # Setters

    def from_unsigned_int(self, value):
        assert self.min_unsigned <= value <= self.max_unsigned

        for i in xrange(self.width):
            self.bv[self.width - 1 - i] = 1 if (2 ** i) & value else 0

    def from_signed_int(self, value):
        assert self.min_signed <= value <= self.max_signed

        if value < 0:
            self.from_unsigned_int(self.max_unsigned + value + 1)
        else:
            self.from_unsigned_int(value)


    ######################################################################
    # Conversion

    def to_unsigned_int(self):
        as_ubv = 0
        for b in self.bv:
            as_ubv *= 2
            as_ubv |= b
        return as_ubv

    def to_signed_int(self):
        as_ubv = self.to_unsigned_int()
        if as_ubv <= self.max_signed:
            return as_ubv
        else:
            return -self.max_unsigned + as_ubv - 1

    ######################################################################
    # SMTLIB support

    def smtlib_sort(self):
        return "(_ BitVec %u)" % len(self.bv)

    def smtlib_literals(self):
        choices = []

        # Obvious binary
        choices.append("#b%s" % "".join(map(str, self.bv)))

        if len(self.bv) % 4 == 0:
            pass

        return choices

    def smtlib_literal(self):
        return self.smtlib_literals()[-1]

    def smtlib_random_literal(self):
        return random.choice(self.smtlib_literals())
