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

# Helper subprograms for writing SMTLIB files.

from out_common import *

##############################################################################
# SMTLIB output
##############################################################################

def smt_opsname(fp_ops):
    if fp_ops == "smtlib.eq":
        return "="
    else:
        return fp_ops

def smt_write_header(fd, status, comment=None, logic="QF_FP"):
    assert status in ("sat", "unsat")
    fd.write("(set-info :smt-lib-version 2.6)\n")
    fd.write("(set-logic %s)\n" % logic)
    if status == "sat":
        fd.write("(set-option :produce-models true)\n")
    fd.write("(set-info :source |%s|)\n" % SOURCE_STRING)
    fd.write("(set-info :license |%s|)\n" %
             "https://www.gnu.org/licenses/gpl-3.0.html")
    fd.write("(set-info :category random)\n")
    fd.write("(set-info :status %s)\n" % status)
    if comment is not None:
        smt_write_comment(fd, comment)

def smt_write_comment(fd, comment):
    fd.write(";; %s\n" % comment)

def smt_write_footer(fd):
    fd.write("(check-sat)\n")
    fd.write("(exit)\n")

def smt_write_var(fd, var_name, var_type, assertion = None, expectation = None):
    fd.write("(declare-const %s %s)\n" % (var_name, var_type))
    if assertion is not None:
        fd.write("(assert %s)\n" % assertion)
    if expectation is not None:
        smt_write_comment(fd, "%s should be %s\n" % (var_name, expectation))

def smt_write_vars(fd, test_vector):
    assert 1 <= len(test_vector["values"]) <= 3

    var_names = ["x", "y", "z"]
    for var, fp in zip(var_names, test_vector["values"]):
        smt_write_var(fd,
                      var,
                      fp.smtlib_sort(),
                      "(= %s %s)" % (var, fp.smtlib_random_literal()),
                      str(fp))

def smt_write_goal(fd, bool_expr, expectation, correct_answer = True):
    assert expectation in ("sat", "unsat")
    fd.write("(assert ")
    if expectation == "unsat":
        fd.write("(not ")
    if correct_answer == True:
        fd.write(bool_expr)
    else:
        fd.write("(not %s)" % bool_expr)
    if expectation == "unsat":
        fd.write(")")
    fd.write(")\n")
