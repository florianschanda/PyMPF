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

# Basic intervals for rationals (including infinity). This is good for
# the basic intervals for int/float conversions, however in the future
# we need more. For example float intervals should have special flags
# for including zero, infinity, nan in addition to a specific range.
#
# Representation uses the "German" method, i.e. "]2 3]" is the
# interval between 2 and 3, excluding 2 but including 3.

KIND_INFINITE  = "infinite"
KIND_INCLUSIVE = "inclusive"
KIND_EXCLUSIVE = "exclusive"

class Interval_Bound:
    INTERVAL_KINDS = (KIND_INFINITE, KIND_INCLUSIVE, KIND_EXCLUSIVE)

    def __init__(self, kind, q=None):
        assert kind in Interval_Bound.INTERVAL_KINDS
        assert not kind == KIND_INFINITE or q is None
        assert not kind in (KIND_INCLUSIVE, KIND_EXCLUSIVE) or q is not None
        self.kind  = kind
        self.value = q

class Interval:
    def __init__(self):
        self.low  = Interval_Bound(KIND_INFINITE)
        self.high = Interval_Bound(KIND_INFINITE)

    def __str__(self):
        rv = ""
        if self.low.kind == KIND_INFINITE:
            rv += "] -oo"
        else:
            if self.low.kind == KIND_INCLUSIVE:
                rv += "[ "
            else:
                rv += "] "
            rv += str(self.low.value)
        rv += " .. "
        if self.high.kind == KIND_INFINITE:
            rv += "+oo ["
        else:
            rv += str(self.high.value)
            if self.high.kind == KIND_INCLUSIVE:
                rv += " ]"
            else:
                rv += " ["
        return rv

    def set_low(self, q, inclusive):
        assert isinstance(inclusive, bool)
        if inclusive:
            self.low.kind = KIND_INCLUSIVE
        else:
            self.low.kind = KIND_EXCLUSIVE
        self.low.value = q

    def set_high(self, q, inclusive):
        assert isinstance(inclusive, bool)
        if inclusive:
            self.high.kind = KIND_INCLUSIVE
        else:
            self.high.kind = KIND_EXCLUSIVE
        self.high.value = q
