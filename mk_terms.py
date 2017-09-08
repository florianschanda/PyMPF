#!/usr/bin/env python
##############################################################################
##                                                                          ##
##                                PYMPF                                     ##
##                                                                          ##
##              Copyright (C) 2017, Altran UK Limited                       ##
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

# This program will generate a large number of random terms.

import argparse
import itertools
import multiprocessing
import os
import resource
import shutil
import subprocess
import tempfile
from glob import glob
from copy import copy

from floats import *
from out_smtlib import *

REF_TYPE = MPF(8, 24)

def call(name, rounding, *args):
    return {
        "kind"      : "call",
        "name"      : name,
        "rounding"  : rounding,
        "args"      : list(args),
        "unordered" : False,     # +, *, or, etc.
        "distinct"  : False,     # true for <, <=, =
    }

def unordered(n):
    assert n["kind"] == "call"
    rv = n.copy()
    rv["unordered"] = True
    return rv

def distinct(n):
    assert n["kind"] == "call"
    rv = n.copy()
    rv["distinct"] = True
    return rv

def pp_tree(n, rounding="RNE"):
    if n["kind"] == "lit":
        return n["name"]
    elif n["kind"] == "var":
        return n["name"]
    elif n["kind"] == "call":
        if n["name"] == "fp.isFinite":
            return ("(not (or (fp.isInfinite %s) (fp.isNaN %s)))" %
                    (pp_tree(n["args"][0], rounding),
                     pp_tree(n["args"][0], rounding)))
        elif n["name"] == "fp.isIntegral":
            return ("(fp.eq (fp.roundToIntegral %s %s) %s)" %
                    (rounding,
                     pp_tree(n["args"][0], rounding),
                     pp_tree(n["args"][0], rounding)))
        else:
            rv = "(%s" % n["name"]
            if n["rounding"]:
                rv += " %s" % rounding
            for arg in n["args"]:
                rv += " %s" % pp_tree(arg, rounding)
            rv += ")"
            return rv
    else:
        assert False

def tree_contains(n, predicate):
    if predicate(n):
        return True
    if n["kind"] == "lit":
        pass
    elif n["kind"] == "var":
        pass
    elif n["kind"] == "call":
        for arg in n["args"]:
            if tree_contains(arg, predicate):
                return True
    else:
        assert False

def tree_count(n, predicate):
    c = 0
    if predicate(n):
        c += 1
    if n["kind"] == "lit":
        pass
    elif n["kind"] == "var":
        pass
    elif n["kind"] == "call":
        for arg in n["args"]:
            c += tree_count(arg, predicate)
    else:
        assert False
    return c

def tree_contains_var(n):
    return tree_contains(n, lambda x: x["kind"] == "var")

def tree_lt(n1, n2):
    return pp_tree(n1) < pp_tree(n2)

def tree_leq(n1, n2):
    return pp_tree(n1) <= pp_tree(n2)

def tree_eq(n1, n2):
    return pp_tree(n1) == pp_tree(n2)

def tree_sub(n, kind, var_asn):
    assert kind in ("lit", "var")
    todo = copy(var_asn)
    def rec(n):
        rv = n.copy()
        if rv["kind"] == kind:
            rv["name"] = todo.pop()
        elif rv["kind"] == "call":
            rv["args"] = map(rec, n["args"])
        return rv
    return rec(n)

def tree_check_constraint(n):
    if n["kind"] == "lit":
        return True
    elif n["kind"] == "var":
        return True
    elif n["kind"] == "call":
        for arg in n["args"]:
            if not tree_check_constraint(arg):
                return False

        if n["unordered"]:
            assert len(n["args"]) == 2
            if not tree_leq(n["args"][0], n["args"][1]):
                return False

        if n["distinct"]:
            assert len(n["args"]) == 2
            if tree_eq(n["args"][0], n["args"][1]):
                return False

        return True
    else:
        assert False

def tree_get_vars(n):
    if n["kind"] == "lit":
        return set()
    elif n["kind"] == "var":
        return set(n["name"])
    elif n["kind"] == "call":
        return reduce(lambda x, y: x | y,
                      map(tree_get_vars, n["args"]),
                      set())
    else:
        assert False

def gen_float_expressions(depth):
    if depth == 0:
        return

    yield {
        "kind" : "lit",
        "name" : "LIT",
    }
    yield {
        "kind" : "var",
        "name" : "VAR",
    }

    for child_1 in gen_float_expressions(depth - 1):
        if tree_contains_var(child_1):
            yield call("fp.abs", False, child_1)
            yield call("fp.neg", False, child_1)
            yield call("fp.sqrt", True, child_1)
            yield call("fp.roundToIntegral", True, child_1)

        for child_2 in gen_float_expressions(depth - 1):
            if tree_contains_var(child_1) or tree_contains_var(child_2):
                # yield call("fp.rem", False, child_1, child_2)
                yield call("fp.div", True, child_1, child_2)
                if tree_leq(child_1, child_2):
                    yield unordered(call("fp.add", True, child_1, child_2))
                    yield unordered(call("fp.mul", True, child_1, child_2))


def gen_bool_expressions(depth):
    if depth == 0:
        return

    for child_1 in gen_bool_expressions(depth - 1):
        yield call("not", False, child_1)

        for child_2 in gen_bool_expressions(depth - 1):
            if tree_leq(child_1, child_2):
                yield unordered(distinct(call("or", False, child_1, child_2)))

    for child_1 in gen_float_expressions(depth - 1):
        if tree_contains_var(child_1):
            yield call("fp.isZero",      False, child_1)
            yield call("fp.isSubnormal", False, child_1)
            yield call("fp.isNormal",    False, child_1)
            yield call("fp.isInfinite",  False, child_1)
            yield call("fp.isNaN",       False, child_1)
            yield call("fp.isPositive",  False, child_1)
            yield call("fp.isNegative",  False, child_1)
            yield call("fp.isFinite",    False, child_1)
            yield call("fp.isIntegral",  False, child_1)

        for child_2 in gen_float_expressions(depth - 1):
            if tree_contains_var(child_1) or tree_contains_var(child_2):
                yield distinct(call("fp.lt", False, child_1, child_2))
                yield distinct(call("fp.leq", False, child_1, child_2))
                if tree_leq(child_1, child_2):
                    yield unordered(distinct(call("=", False, child_1, child_2)))

def gen_vars(num_occ, num_vars):
    assert num_occ >= 1
    assert num_vars >= 1

    rv = [0] * num_occ
    yield rv

    while rv != [num_vars-1] * num_occ:
        for n in xrange(num_occ):
            rv[n] += 1
            if rv[n] == num_vars:
                rv[n] = 0
            else:
                break
        yield rv

def gen_expanded_bool_expressions(depth, num_vars, literals):
    def v_map(a):
        return [chr(ord('a') + x) for x in a]

    def l_map(a):
        return [literals[x] for x in a]

    for t in gen_bool_expressions(depth):
        var_occ = tree_count(t, lambda x: x["kind"] == "var")
        for v_asn in gen_vars(var_occ, num_vars):
            tmp_a = tree_sub(t, "var", v_map(v_asn))
            if tree_check_constraint(tmp_a):
                lit_occ = tree_count(tmp_a, lambda x: x["kind"] == "lit")
                if lit_occ == 0:
                    yield tmp_a
                else:
                    for l_asn in gen_vars(lit_occ, len(literals)):
                        tmp_b = tree_sub(tmp_a, "lit", l_map(l_asn))
                        if tree_check_constraint(tmp_b):
                            yield tmp_b

def check_status(timeout, solver, *terms):
    def set_limit():
        resource.setrlimit(resource.RLIMIT_CPU, (timeout, timeout))

    all_vars = reduce(lambda x, y: x | y,
                      map(tree_get_vars, terms),
                      set())

    need_tmp_file = "colibri" in solver[0]

    tmp = []
    tmp.append("(set-info :smt-lib-version 2.6)")
    tmp.append("(set-logic QF_FP)")
    for v in sorted(all_vars):
        tmp.append("(declare-const %s %s)" % (v, REF_TYPE.smtlib_sort()))
    for t in terms:
        tmp.append("(assert %s)" % pp_tree(t))
    tmp.append("(check-sat)")
    tmp.append("(exit)")

    if need_tmp_file:
        fd, fn = tempfile.mkstemp(".smt2", text=True)
        for line in tmp:
            os.write(fd, line + "\n")
        os.close(fd)
        cmd = copy(solver) + [fn]
    else:
        cmd = solver

    p = subprocess.Popen(cmd,
                         stdin=subprocess.PIPE,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE,
                         preexec_fn=set_limit)

    if need_tmp_file:
        stdout, stderr = p.communicate()
        os.unlink(fn)
    else:
        stdout, stderr = p.communicate("\n".join(tmp) + "\n")
    stdout = stdout.strip()
    stderr = stderr.strip()

    if stdout == "sat":
        return "sat"
    elif stdout == "unsat":
        return "unsat"
    elif stdout == "" and stderr == "":
        return "timeout"
    else:
        print "\n".join(tmp)
        raise Exception(stdout + "\n" + stderr.strip())

def process(task):
    timeout = 5

    rv = {"status" : None,
          "term"   : task["goal"]}

    status_a = check_status(timeout, task["solver"], task["goal"])
    if status_a == "timeout":
        rv["status"] = "hard"
        return rv

    not_goal = call("not", False, task["goal"])
    status_b = check_status(timeout, task["solver"], not_goal)

    if status_b == "timeout":
        rv["status"] = "hard"
        rv["term"]   = not_goal
        return rv

    if status_a == "unsat" and status_b == "sat":
        rv["status"] = "unsat"
    elif status_b == "unsat" and status_a == "sat":
        rv["status"] = "unsat"
        rv["term"]   = not_goal
    else:
        rv["status"] = "not interesting"

    return rv


def mk_theorems(depth, solver, single_thread):
    os.mkdir("theorem.inst")

    literals = []

    f = REF_TYPE.new_mpf()
    for sign in (0, 1):
        f.set_zero(sign)
        literals.append(f.smtlib_literal())

        f.set_infinite(sign)
        literals.append(f.smtlib_literal())

        f.from_rational(RM_RNE, Rational(1))
        f.set_sign_bit(sign == 0)
        literals.append(f.smtlib_literal())

    f.set_nan()
    literals.append(f.smtlib_literal())

    term_generator = gen_expanded_bool_expressions(depth,
                                                   2,
                                                   literals)
    tasks = [{"solver" : solver,
              "goal"   : x}
             for x in term_generator]


    if single_thread:
        iterator = itertools.imap(process, tasks)
    else:
        pool = multiprocessing.Pool()
        iterator = pool.imap(process, tasks, 10)

    n = 0
    for result in iterator:
        n += 1
        if result["status"] == "unsat":
            pass
            #print "Unsat:", pp_tree(result["term"])

        elif result["status"] == "hard":
            print "<%.1f%%> Hard: %s" % (float(n*100)/len(terms),
                                         pp_tree(result["term"]))

            with open("theorem.inst/hard_%08u.smt2" % n, "w") as fd:
                smt_write_header(fd, "unknown")
                all_vars = tree_get_vars(result["term"])
                for v in all_vars:
                    smt_write_var(fd, v, REF_TYPE.smtlib_sort())
                smt_write_goal(fd, pp_tree(result["term"]), "sat")
                smt_write_footer(fd)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--depth",
                    default=3)
    ap.add_argument("--solver",
                    default="cvc4",
                    choices=["cvc4", "z3", "mathsat", "colibri"])
    ap.add_argument("--single",
                    default=False,
                    action="store_true")
    options = ap.parse_args()

    for d in glob("theorem.*"):
        if os.path.isdir(d):
            shutil.rmtree(d)

    if options.solver == "cvc4":
        solver = ["../smtlib_schanda/cvc4_2017_08_18",
                  "--lang=smt2",
                  "--check-models"]
    elif options.solver == "z3":
        solver = ["../smtlib_schanda/z3_2017_08_02",
                  "-smt2",
                  "-in"]
    elif options.solver == "mathsat":
        solver = ["mathsat",
                  "-input=smt2"]
    elif options.solver == "colibri":
        solver = ["colibri"]
    else:
        assert False

    mk_theorems(options.depth, solver, options.single)

if __name__ == "__main__":
    main()

# 1. gen (assert term1) (assert term2) ... (assert termN) and check if unsat
# 2. gen (assert (distinct f_term f_term)) and check if unsat
