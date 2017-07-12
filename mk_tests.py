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

# This is the main testcase generator.

import os
import shutil
import random
import argparse
from glob import glob

from bitvector import *
from floats import *
from rationals import *
from interval_q import *
from evaluator import (fp_eval_predicate,
                       fp_eval_function,
                       is_rounding,
                       all_ops_where,
                       TYP_BOOL, TYP_BV, TYP_REAL, TYP_FLOAT)
from out_smtlib import *

##############################################################################
# Random floats
##############################################################################

def random_zero(sign=0):
    rv = MPF(8, 24)
    if sign > 0 or (sign == 0 and random.getrandbits(1)):
        rv.set_zero(0)
    else:
        rv.set_zero(1)
    return rv

def random_subnormal(sign=0):
    rv = MPF(8, 24)
    S = (0 if sign > 0 or (sign == 0 and random.getrandbits(1)) else 1)
    E = 0
    T = random.randrange(1, 2 ** rv.t)
    rv.pack(S, E, T)
    return rv

def random_normal(sign=0):
    rv = MPF(8, 24)
    S = (0 if sign > 0 or (sign == 0 and random.getrandbits(1)) else 1)
    E = random.randrange(1, 2 ** rv.w - 1)
    T = random.randrange(0, 2 ** rv.t)
    rv.pack(S, E, T)
    return rv

def random_infinite(sign=0):
    rv = MPF(8, 24)
    if sign > 0 or (sign == 0 and random.getrandbits(1)):
        rv.set_infinite(0)
    else:
        rv.set_infinite(1)
    return rv

def random_nan():
    rv = MPF(8, 24)
    S = random.getrandbits(1)
    E = 2 ** rv.w - 1
    T = random.randrange(1, 2 ** rv.t)
    rv.pack(S, E, T)
    return rv

def smallest_subnormal(sign=0):
    rv = MPF(8, 24)
    if sign > 0 or (sign == 0 and random.getrandbits(1)):
        S = 0
    else:
        S = 1
    E = 0
    T = 1
    rv.pack(S, E, T)
    return rv

def largest_subnormal(sign=0):
    rv = MPF(8, 24)
    if sign > 0 or (sign == 0 and random.getrandbits(1)):
        S = 0
    else:
        S = 1
    E = 0
    T = 2 ** rv.t - 1
    rv.pack(S, E, T)
    return rv

def smallest_normal(sign=0):
    rv = MPF(8, 24)
    if sign > 0 or (sign == 0 and random.getrandbits(1)):
        S = 0
    else:
        S = 1
    E = 1
    T = 0
    rv.pack(S, E, T)
    return rv

def largest_normal(sign=0):
    rv = MPF(8, 24)
    if sign > 0 or (sign == 0 and random.getrandbits(1)):
        S = 0
    else:
        S = 1
    E = 2 ** rv.w - 2
    T = 2 ** rv.t - 1
    rv.pack(S, E, T)
    return rv

def ubv_boundary(bv_width, sign=0):
    bv = BitVector(bv_width)
    if sign > 0 or (sign == 0 and random.getrandbits(1)):
        q = Rational(bv.max_unsigned)
    else:
        q = Rational(bv.min_unsigned)
    rv = MPF(8, 24)
    rv.from_rational(RM_RNE, q)
    return rv

def sbv_boundary(bv_width, sign=0):
    bv = BitVector(bv_width)
    if sign > 0 or (sign == 0 and random.getrandbits(1)):
        q = Rational(bv.max_signed)
    else:
        q = Rational(bv.min_signed)
    rv = MPF(8, 24)
    rv.from_rational(RM_RNE, q)
    return rv

##############################################################################
# Random rationals
##############################################################################

def random_rational(low, high):
    if low.kind == KIND_INFINITE:
        assert high.kind != KIND_INFINITE
        factor = min(10000, high.value.b ** 3)
        high_a = high.value.a * factor
        b      = high.value.b * factor

        if high.kind == KIND_EXCLUSIVE:
            high_a -= 1

        a = random.randint(-abs(high_a) - abs(factor), high_a)
        return Rational(a, b)

    elif high.kind == KIND_INFINITE:
        factor = min(10000, low.value.b ** 3)
        low_a = low.value.a * factor
        b     = low.value.b * factor

        if low.kind == KIND_EXCLUSIVE:
            low_a += 1

        a = random.randint(low_a, abs(low_a) + abs(factor))
        return Rational(a, b)

    else:
        if low.value < high.value:
            factor = 25
            low_a = low.value.a * high.value.b * factor
            high_a = high.value.a * low.value.b * factor
            b = low.value.b * high.value.b * factor

            if low.kind == KIND_EXCLUSIVE:
                low_a += 1
            if high.kind == KIND_EXCLUSIVE:
                high_a -= 1
            a = random.randint(low_a, high_a)
            return Rational(a, b)
        else:
            return low.value


##############################################################################
# Testvector generation
##############################################################################

def gen_rm(fp_ops):
    if is_rounding(fp_ops):
        for rm in MPF.ROUNDING_MODES:
            yield rm
    else:
        yield RM_RNE

def gen_vectors(fp_ops, n, test_dup):
    assert n >= 1
    assert test_dup >= 1

    CONSTRUCTORS = {
        "-0"         : lambda: random_zero(-1),
        "+0"         : lambda: random_zero(1),
        "-minsub"    : lambda: smallest_subnormal(-1),
        "+minsub"    : lambda: smallest_subnormal(1),
        "-subnormal" : lambda: random_subnormal(-1),
        "+subnormal" : lambda: random_subnormal(1),
        "-maxsub"    : lambda: largest_subnormal(-1),
        "+maxsub"    : lambda: largest_subnormal(1),
        "-minnormal" : lambda: smallest_normal(-1),
        "+minnormal" : lambda: smallest_normal(1),
        "-normal"    : lambda: random_normal(-1),
        "+normal"    : lambda: random_normal(1),
        "-maxnormal" : lambda: largest_normal(-1),
        "+maxnormal" : lambda: largest_normal(1),
        "-inf"       : lambda: random_infinite(-1),
        "+inf"       : lambda: random_infinite(1),
        "nan"        : lambda: random_nan(),
    }
    if fp_ops in ("fp.to.ubv", "fp.to.sbv"):
        for k in ("-minsub", "+minsub",
                  "-maxsub", "+maxsub"):
            del CONSTRUCTORS[k]
        CONSTRUCTORS["-sbv_8_bound"] = lambda: sbv_boundary(8, -1)
        CONSTRUCTORS["+sbv_8_bound"] = lambda: sbv_boundary(8, 1)
        CONSTRUCTORS["+ubv_8_bound"] = lambda: ubv_boundary(8, 1)

    TARGETS   = tuple(sorted(list(CONSTRUCTORS)))
    N_TARGETS = len(TARGETS)

    history = set()

    def mk_float(c):
        assert 0 <= c < N_TARGETS
        return CONSTRUCTORS[TARGETS[c]]()

    for rm in gen_rm(fp_ops):
        classes = [0] * n
        while True:
            text = map(lambda x: TARGETS[x], classes)
            for _ in xrange(test_dup):
                v_exp = random.choice(["sat", "unsat"])
                v_val = map(mk_float, classes)
                h     = tuple([v_exp] + [rm] + map(lambda x: x.bv, v_val))
                comment = "(" + fp_ops
                if is_rounding(fp_ops):
                    comment += " " + rm
                comment += " "
                comment += " ".join(text)
                comment += ")"
                if h not in history:
                    history.add(h)
                    yield {
                        "ops"         : fp_ops,
                        "rounding"    : rm,
                        "expectation" : v_exp,
                        "values"      : v_val,
                        "comment"     : comment,
                    }
            k = n - 1
            while (k >= 0):
                if classes[k] < (N_TARGETS - 1):
                    classes[k] += 1
                    break
                else:
                    assert classes[k] == (N_TARGETS - 1)
                    classes[k] = 0
                    k -= 1
            if k == -1:
                break

def gen_bv_vectors(fp_ops, n, test_dup):
    assert n == 1

    for bv_width in (8, 16, 32, 64, 128):
        bv = BitVector(bv_width)
        comment = "%s(BitVec %u)" % (fp_ops, bv_width)
        for bit_0 in (0, 1):
            bv.bv[0] = bit_0
            for bit_1 in (0, 1):
                bv.bv[1] = bit_1
                for bit_penultimate in (0, 1):
                    bv.bv[-2] = bit_penultimate
                    for bit_last in (0, 1):
                        bv.bv[-1] = bit_last
                        for rm in gen_rm(fp_ops):
                            # Zeros in middle
                            bv.bv[2:-2] = [0] * (bv_width - 4)
                            yield {
                                "ops"         : fp_ops,
                                "rounding"    : rm,
                                "expectation" : random.choice(["sat", "unsat"]),
                                "values"      : [bv],
                                "comment"     : comment,
                            }

                            # Ones in middle
                            bv.bv[2:-2] = [1] * (bv_width - 4)
                            yield {
                                "ops"         : fp_ops,
                                "rounding"    : rm,
                                "expectation" : random.choice(["sat", "unsat"]),
                                "values"      : [bv],
                                "comment"     : comment,
                            }

                            # Random in middle
                            for _ in xrange(test_dup):
                                for b in xrange(2, 2 + bv_width - 4):
                                    bv.bv[b] = random.randint(0, 1)
                                yield {
                                    "ops"         : fp_ops,
                                    "rounding"    : rm,
                                    "expectation" : random.choice(["sat",
                                                                   "unsat"]),
                                    "values"      : [bv],
                                    "comment"     : comment,
                                }


##############################################################################
# Test generation
##############################################################################

test_id = 0

def new_test(testvec):
    global test_id

    if not os.path.exists(testvec["ops"]):
        os.mkdir(testvec["ops"])

    test_id += 1

    test_name = testvec["ops"]
    if test_name.startswith("fp."):
        test_name = test_name[3:]
    if is_rounding(testvec["ops"]):
        test_name += "_" + testvec["rounding"].lower()
    test_name += "_%05u.smt2" % test_id

    print ">>> Generating test %u %s" % (test_id, testvec["comment"])
    return open(os.path.join(testvec["ops"], test_name), "w")


def mk_tests_for_classify(num_tests):
    for fp_ops in all_ops_where(arity=1, result=TYP_BOOL):
        # check if X is correctly classified
        for vec in gen_vectors(fp_ops, 1, num_tests):
            x      = vec["values"][0]
            result = fp_eval_predicate(fp_ops, x)

            with new_test(vec) as fd:
                smt_write_header(fd, vec["expectation"])
                smt_write_vars(fd, vec)
                smt_write_var(fd, "result", "Bool",
                              "(= result (%s x))" % smt_opsname(fp_ops))
                smt_write_goal(fd, "result", vec["expectation"], result)
                smt_write_footer(fd)

def mk_tests_for_unary(num_tests):
    for fp_ops in all_ops_where(arity=1, args=TYP_FLOAT, result=TYP_FLOAT):
        # pick x. compute op(x) = result. check result
        for vec in gen_vectors(fp_ops, 1, num_tests):
            x      = vec["values"][0]
            result = fp_eval_function(fp_ops, vec["rounding"], x)

            with new_test(vec) as fd:
                smt_write_header(fd, vec["expectation"])
                smt_write_vars(fd, vec)
                if is_rounding(fp_ops):
                    smt_write_var(fd, "result", result.smtlib_sort(),
                                  "(= result (%s %s x))" % (smt_opsname(fp_ops),
                                                            vec["rounding"]))
                else:
                    smt_write_var(fd, "result", result.smtlib_sort(),
                                  "(= result (%s x))" % smt_opsname(fp_ops))
                smt_write_goal(fd,
                               "(= result %s)" % result.smtlib_random_literal(),
                               vec["expectation"])
                smt_write_footer(fd)

        # pick result. compute op(x) = result. check x
        # TODO

def mk_tests_for_relations(num_tests):
    for fp_ops in all_ops_where(arity=2, result=TYP_BOOL):
        # pick x, y. compute result = x OP y. check result
        for vec in gen_vectors(fp_ops, 2, num_tests):
            x, y   = vec["values"]
            result = fp_eval_predicate(fp_ops, x, y)

            with new_test(vec) as fd:
                smt_write_header(fd, vec["expectation"])
                smt_write_vars(fd, vec)

                smt_write_var(fd, "result", "Bool",
                              "(= result (%s x y))" % smt_opsname(fp_ops))

                smt_write_goal(fd, "result", vec["expectation"], result)
                smt_write_footer(fd)

def mk_tests_for_binary(num_tests):
    for fp_ops in all_ops_where(arity=2, result=TYP_FLOAT):
        # pick x, y. compute result = x OP y. check result
        for vec in gen_vectors(fp_ops, 2, num_tests):
            x, y = vec["values"]
            try:
                result = fp_eval_function(fp_ops, vec["rounding"], x, y)
            except Unspecified:
                continue

            with new_test(vec) as fd:
                smt_write_header(fd, vec["expectation"])
                smt_write_vars(fd, vec)

                if is_rounding(fp_ops):
                    smt_write_var(fd, "result", result.smtlib_sort(),
                                  "(= result (%s %s x y))" %
                                  (smt_opsname(fp_ops), vec["rounding"]))
                else:
                    smt_write_var(fd, "result", result.smtlib_sort(),
                                  "(= result (%s x y))" % smt_opsname(fp_ops))
                smt_write_goal(fd,
                               "(= result %s)" %
                               result.smtlib_random_literal(),
                               vec["expectation"])
                smt_write_footer(fd)

def mk_tests_for_ternary(num_tests):
    fp_ops = "fp.fma"

    # pick x, y. compute result = x OP y. check result
    for vec in gen_vectors("fp.fma", 3, num_tests):
        x, y, z = vec["values"]
        result  = fp_eval_function(fp_ops, vec["rounding"], x, y, z)

        with new_test(vec) as fd:
            smt_write_header(fd, vec["expectation"])
            smt_write_vars(fd, vec)
            smt_write_var(fd, "result", result.smtlib_sort(),
                          "(= result (%s %s x y z))" % (smt_opsname(fp_ops),
                                                        vec["rounding"]))
            smt_write_goal(fd,
                           "(= result %s)" %
                           result.smtlib_random_literal(),
                           vec["expectation"])
            smt_write_footer(fd)

def mk_tests_from_bitvector(num_tests):
    # convert signed and unsigned bitvector to float
    def mk_test_from(vec):
        assert vec["ops"] in ("fp.from.ubv", "fp.from.sbv")

        interval = random.choice(["<=", "=", ">="])

        x = vec["values"][0]

        fmt = random.choice([
            (5, 11),
            (8, 24),
            (11, 53),
            #(15, 113),
            #(random.randint(3, 10), random.randint(3, 10)),
        ])
        f = MPF(fmt[0], fmt[1])

        if vec["ops"] == "fp.from_ubv":
            as_int = x.to_unsigned_int()
            bv_rel = {"<=" : "bvule",
                      "="  : "=",
                      ">=" : "bvuge"}[interval]
            fp_ops = f.smtlib_from_ubv()
        else:
            as_int = x.to_signed_int()
            bv_rel = {"<=" : "bvsle",
                      "="  : "=",
                      ">=" : "bvsge"}[interval]
            fp_ops = f.smtlib_from_sbv()

        q = Rational(as_int)
        f.from_rational(vec["rounding"], q)

        with new_test(vec) as fd:
            smt_write_header(fd,
                             status  = vec["expectation"],
                             comment = vec["comment"],
                             logic   = "QF_FPBV")
            smt_write_var(
                fd,
                var_name    = "x",
                var_type    = x.smtlib_sort(),
                assertion   = "(%s x %s)" % (bv_rel,
                                             x.smtlib_random_literal()),
                expectation = str(as_int))
            smt_write_var(
                fd,
                var_name    = "r",
                var_type    = f.smtlib_sort(),
                assertion   = "(= r (%s %s x))" % (fp_ops,
                                                   vec["rounding"]))
            smt_write_goal(
                fd,
                bool_expr   = "(%s r %s)" % ({"<=" : "fp.leq",
                                              "="  : "fp.eq",
                                              ">=" : "fp.geq"}[interval],
                                             f.smtlib_random_literal()),
                expectation = vec["expectation"])
            smt_write_footer(fd)

    # convert float (+/- some interesting values) to a bitvector
    def mk_test_to(vec, bv_width, shift=Rational(0)):
        assert vec["ops"] in ("fp.to.ubv", "fp.to.sbv")

        bv = BitVector(bv_width)
        for i in xrange(bv_width):
            bv.bv[i] = random.randint(0, 1)

        x = vec["values"][0]
        if not shift.isZero():
            if x.isInfinite() or x.isNaN():
                return
            fudge = x.new_mpf()
            fudge.from_rational(vec["rounding"], shift)
            x = fp_add(vec["rounding"], x, fudge)

        unspecified = x.isNaN() or x.isInfinite()

        if not unspecified:
            i = fp_roundToIntegral(vec["rounding"], x)
            q = i.to_rational()
            assert q.isIntegral()

            if vec["ops"] == "fp.to.ubv":
                unspecified = not (Rational(bv.min_unsigned) <=
                                   q <=
                                   Rational(bv.max_unsigned))
            else:
                unspecified = not (Rational(bv.min_signed) <=
                                   q <=
                                   Rational(bv.max_signed))

        if not unspecified:
            if vec["ops"] == "fp.to.ubv":
                bv.from_unsigned_int(q.a)
            else:
                bv.from_signed_int(q.a)
            expectation = vec["expectation"]
        else:
            expectation = "sat"

        with new_test(vec) as fd:
            smt_write_header(fd, expectation, vec["comment"], "QF_FPBV")
            smt_write_var(
                fd,
                var_name    = "x",
                var_type    = x.smtlib_sort(),
                assertion   = "(= x %s)" % x.smtlib_random_literal())
            smt_write_var(
                fd,
                var_name    = "y",
                var_type    = bv.smtlib_sort(),
                assertion   = "(= y ((_ %s %u) %s x))" % (
                    {"fp.to.ubv" : "fp.to_ubv",
                     "fp.to.sbv" : "fp.to_sbv"}[vec["ops"]],
                    bv.width,
                    vec["rounding"]),
                expectation = "unspecified" if unspecified else None)
            if not unspecified or random.randint(0, 1) == 0:
                z_assertion = "(= z %s)" % bv.smtlib_random_literal()
            else:
                z_assertion = None
            smt_write_var(
                fd,
                var_name    = "z",
                var_type    = bv.smtlib_sort(),
                assertion   = z_assertion)
            smt_write_goal(
                fd,
                bool_expr   = "(= y z)",
                expectation = expectation,
                correct_answer = random.choice([True, False]) if unspecified else True)
            smt_write_footer(fd)


    for vec in gen_bv_vectors("fp.from.ubv", 1, num_tests):
        mk_test_from(vec)

    for vec in gen_bv_vectors("fp.from.sbv", 1, num_tests):
        mk_test_from(vec)

    q_half = Rational(1, 2)
    for vec in gen_vectors("fp.to.ubv", 1, num_tests):
        for bv_width in (8, 32, 64):
            mk_test_to(vec, bv_width)
            mk_test_to(vec, bv_width, shift=q_half)
            mk_test_to(vec, bv_width, shift=-q_half)

    for vec in gen_vectors("fp.to.sbv", 1, num_tests):
        for bv_width in (8, 32, 64):
            mk_test_to(vec, bv_width)
            mk_test_to(vec, bv_width, shift=q_half)
            mk_test_to(vec, bv_width, shift=-q_half)

def mk_tests_for_real_to_float(num_tests):
    # We pick a random float; then determine the rational interval that would
    # round onto that float. we then select rationals from:
    #   - outisde that interval,
    #   - inside the interval,
    #   - just on the boundary
    # We then apply ((_ to_fp eb sb) rm RATIONAL) and make sure we get the
    # correct result.

    for vec in gen_vectors("fp.from.real", 1, num_tests):
        x = vec["values"][0]
        if x.isNaN():
            # This can't happen, so lets skip these
            continue
        if vec["rounding"] != RM_RNA:
            continue

        interval = fp_interval(vec["rounding"], x)
        if interval is None:
            assert vec["rounding"] in MPF.ROUNDING_MODES_DIRECTED
            # It is not possible to construct an interval in some cases
            continue

        def emit_test(q, result, expectation, comment):
            tmp = x.new_mpf()
            tmp.from_rational(vec["rounding"], q)
            #print x, q, result, expectation, comment
            #assert smtlib_eq(tmp, x) == result

            with new_test(vec) as fd:
                smt_write_header(fd, expectation, comment)
                smt_write_vars(fd, vec)
                smt_write_var(fd, "w", x.smtlib_sort(),
                              "(= w (%s %s %s))" % (x.smtlib_from_real(),
                                                    vec["rounding"],
                                                    q.to_smtlib()),
                              tmp.smtlib_literal())
                smt_write_goal(fd,
                               "(%s x w)" % ("="
                                             if result
                                             else "distinct"),
                               expectation)
                smt_write_footer(fd)

        def invert_expectation(expectation):
            if expectation == "unsat":
                return "sat"
            else:
                return "unsat"

        def check_interval_bound(bound, comment):
            if bound.kind == KIND_INCLUSIVE:
                emit_test(bound.value,
                          True,
                          vec["expectation"],
                          comment + " (inclusive)")
            elif bound.kind == KIND_EXCLUSIVE:
                emit_test(bound.value,
                          False,
                          invert_expectation(vec["expectation"]),
                          comment + " (exclusive)")

        # we pick a value somewhere inside the interval
        emit_test(random_rational(interval.low, interval.high),
                  True,
                  vec["expectation"],
                  "inside interval")

        # we then check the upper and lower bound
        check_interval_bound(interval.low, "on low bound")
        check_interval_bound(interval.high, "on high bound")

        # we also pick values outside the interval
        if interval.low.kind != KIND_INFINITE:
            l = Interval_Bound(KIND_INFINITE)
            h = Interval_Bound(KIND_EXCLUSIVE, interval.low.value)
            emit_test(random_rational(l, h), False, vec["expectation"], "below")
        if interval.high.kind != KIND_INFINITE:
            l = Interval_Bound(KIND_EXCLUSIVE, interval.high.value)
            h = Interval_Bound(KIND_INFINITE)
            emit_test(random_rational(l, h), False, vec["expectation"], "above")

def main():
    ap = argparse.ArgumentParser(
        description="Generate random SMTLIB testcases.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--test_classify", metavar="N", type=int,
                    default=0,
                    help="Generate tests for the classification operators.")
    ap.add_argument("--test_unary", metavar="N", type=int,
                    default=0,
                    help="Generate tests for all unary operators.")
    ap.add_argument("--test_binary", metavar="N", type=int,
                    default=0,
                    help="Generate tests for all binary operators.")
    ap.add_argument("--test_ternary", metavar="N", type=int,
                    default=0,
                    help="Generate tests for all ternary operators.")
    ap.add_argument("--test_relations", metavar="N", type=int,
                    default=0,
                    help="Generate tests for all binary relations.")
    ap.add_argument("--test_conversion", metavar="N", type=int,
                    default=0,
                    help="Generate tests for conversion to/from float.")
    options = ap.parse_args()

    for d in glob("fp.*") + glob("smtlib.*"):
        if os.path.isdir(d):
            shutil.rmtree(d)
    for f in glob("eval__*"):
        os.unlink(f)

    if options.test_classify >= 1:
        mk_tests_for_classify(options.test_classify)
    if options.test_relations >= 1:
        mk_tests_for_relations(options.test_relations)
    if options.test_unary >= 1:
        mk_tests_for_unary(options.test_unary)
    if options.test_binary >= 1:
        mk_tests_for_binary(options.test_binary)
    if options.test_ternary >= 1:
        mk_tests_for_ternary(options.test_ternary)
    if options.test_conversion >= 1:
        mk_tests_from_bitvector(options.test_conversion)
        # work in progress
        # mk_tests_for_real_to_float(options.test_to_float)

if __name__ == "__main__":
    main()
