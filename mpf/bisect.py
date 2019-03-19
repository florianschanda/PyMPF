#!/usr/bin/env python3
##############################################################################
##                                                                          ##
##                                PYMPF                                     ##
##                                                                          ##
##              Copyright (C) 2016-2017, Altran UK Limited                  ##
##              Copyright (C) 2019,      Florian Schanda                    ##
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

# Helper class for binary integer search that can be used in a iterator.
#
# Intended use looks a bit like this:
#
#    for guess in Bisect(low, high):
#       # we use guess.value() to get the current guess
#       # we use guess.too_high() or too_low() to adjust

class Bisect:
    def __init__(self, low, high):
        assert isinstance(low, int)
        assert isinstance(high, int)
        assert low <= high

        self.search_min = low
        self.search_max = high
        self.low        = low
        self.high       = high
        self.v          = (self.low + self.high) // 2

    def too_low(self):
        if self.v == self.high:
            self.low = self.high
            self.v = None
            return
        elif self.low + 1 == self.high and self.high == self.search_max:
            self.v = self.high
            return
        elif self.low + 1 == self.high:
            self.v = None
            return
        self.low = self.v
        self.v   = (self.low + self.high) // 2

    def too_high(self):
        if self.low + 1 == self.high and self.low == self.search_min:
            self.high = self.low
            return
        if self.low + 1 >= self.high:
            self.v = None
            return
        self.high = self.v
        self.v   = (self.low + self.high) // 2

    def __iter__(self):
        class Guess:
            def __init__(self, bo):
                self.bo = bo

            def value(self):
                return self.bo.v

            def bounds(self):
                return (self.bo.low, self.bo.high)

            def too_low(self):
                self.bo.too_low()

            def too_high(self):
                self.bo.too_high()

        class It:
            def __init__(self, bo):
                self.bo = bo

            def __iter__(self):
                return self

            def __next__(self):
                if self.bo.v is None:
                    raise StopIteration
                else:
                    return Guess(self.bo)

        return It(self)
