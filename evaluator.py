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

# Here we evaluate floating point expressions using a single function,
# fp_eval_function(FUNCTION, ROUNDING, ARGS...).
#
# However, we also cross-check the answer using a compiled C program,
# and a compiled MPFR program.
#
# There are a number of big issues and omissions that we plan to lift:
#    - RNA is generally not supported
#    - We assume GCC and a modern processor
#    - We assume (and not even check via assertions) Float32

import os
import sys
import subprocess

from floats import *

GCC_FLAGS = " ".join(["-std=c99",
                      "-msse2",
                      "-mfpmath=sse",
                      "-frounding-math",
                      "-fsignaling-nans",
                      "-ffp-contract=off",
                      "-mfma",
                      "-mno-fma4",
                      "-pedantic",
                      "-Wall",
                      "-W",
                      "-Werror",
                  ])

VALID_MPFR_ROUNDING_MODES = (RM_RNE, RM_RTP, RM_RTN, RM_RTZ)

TYP_BOOL  = "boolean"
TYP_FLOAT = "float"
TYP_REAL  = "real"

FP_OPS = {
    "fp.abs"             : {"arity"  : 1,
                            "rnd"    : False,
                            "args"   : TYP_FLOAT,
                            "result" : TYP_FLOAT},
    "fp.neg"             : {"arity"  : 1,
                            "rnd"    : False,
                            "args"   : TYP_FLOAT,
                            "result" : TYP_FLOAT},
    "fp.sqrt"            : {"arity"  : 1,
                            "rnd"    : True,
                            "args"   : TYP_FLOAT,
                            "result" : TYP_FLOAT},
    "fp.roundToIntegral" : {"arity"  : 1,
                            "rnd"    : True,
                            "args"   : TYP_FLOAT,
                            "result" : TYP_FLOAT},
    "fp.isNormal"        : {"arity"  : 1,
                            "rnd"    : False,
                            "args"   : TYP_FLOAT,
                            "result" : TYP_BOOL},
    "fp.isSubnormal"     : {"arity"  : 1,
                            "rnd"    : False,
                            "args"   : TYP_FLOAT,
                            "result" : TYP_BOOL},
    "fp.isZero"          : {"arity"  : 1,
                            "rnd"    : False,
                            "args"   : TYP_FLOAT,
                            "result" : TYP_BOOL},
    "fp.isInfinite"      : {"arity"  : 1,
                            "rnd"    : False,
                            "args"   : TYP_FLOAT,
                            "result" : TYP_BOOL},
    "fp.isNaN"           : {"arity"  : 1,
                            "rnd"    : False,
                            "args"   : TYP_FLOAT,
                            "result" : TYP_BOOL},
    "fp.isNegative"      : {"arity"  : 1,
                            "rnd"    : False,
                            "args"   : TYP_FLOAT,
                            "result" : TYP_BOOL},
    "fp.isPositive"      : {"arity"  : 1,
                            "rnd"    : False,
                            "args"   : TYP_FLOAT,
                            "result" : TYP_BOOL},
    "fp.add"             : {"arity"  : 2,
                            "rnd"    : True,
                            "args"   : TYP_FLOAT,
                            "result" : TYP_FLOAT},
    "fp.sub"             : {"arity"  : 2,
                            "rnd"    : True,
                            "args"   : TYP_FLOAT,
                            "result" : TYP_FLOAT},
    "fp.mul"             : {"arity"  : 2,
                            "rnd"    : True,
                            "args"   : TYP_FLOAT,
                            "result" : TYP_FLOAT},
    "fp.div"             : {"arity"  : 2,
                            "rnd"    : True,
                            "args"   : TYP_FLOAT,
                            "result" : TYP_FLOAT},
    "fp.rem"             : {"arity"  : 2,
                            "rnd"    : False,
                            "args"   : TYP_FLOAT,
                            "result" : TYP_FLOAT},
    "fp.min"             : {"arity"  : 2,
                            "rnd"    : False,
                            "args"   : TYP_FLOAT,
                            "result" : TYP_FLOAT},
    "fp.max"             : {"arity"  : 2,
                            "rnd"    : False,
                            "args"   : TYP_FLOAT,
                            "result" : TYP_FLOAT},
    "fp.eq"              : {"arity"  : 2,
                            "rnd"    : False,
                            "args"   : TYP_FLOAT,
                            "result" : TYP_BOOL},
    "fp.lt"              : {"arity"  : 2,
                            "rnd"    : False,
                            "args"   : TYP_FLOAT,
                            "result" : TYP_BOOL},
    "fp.gt"              : {"arity"  : 2,
                            "rnd"    : False,
                            "args"   : TYP_FLOAT,
                            "result" : TYP_BOOL},
    "fp.leq"             : {"arity"  : 2,
                            "rnd"    : False,
                            "args"   : TYP_FLOAT,
                            "result" : TYP_BOOL},
    "fp.geq"             : {"arity"  : 2,
                            "rnd"    : False,
                            "result" : TYP_BOOL,
                            "args"   : TYP_FLOAT},
    "fp.fma"             : {"arity"  : 3,
                            "rnd"    : True,
                            "args"   : TYP_FLOAT,
                            "result" : TYP_FLOAT},
    # Pseudo ops
    "smtlib.eq"          : {"arity"  : 2,
                            "rnd"    : False,
                            "args"   : TYP_FLOAT,
                            "result" : TYP_BOOL},
    "fp.from.real"       : {"arity"  : 1,
                            "rnd"    : True,
                            "args"   : TYP_REAL,
                            "result" : TYP_FLOAT},
}

def all_ops_where(arity = None, args = None, result = None):
    rv = set()
    for op in FP_OPS:
        print op, FP_OPS[op], arity, args, result
        if arity is not None and arity != FP_OPS[op]["arity"]:
            continue
        if args is not None and args != FP_OPS[op]["args"]:
            continue
        if result is not None and result != FP_OPS[op]["result"]:
            continue
        rv.add(op)
    return rv

def is_rounding(fp_ops):
    return FP_OPS[fp_ops]["rnd"]

##############################################################################
# Writing C
##############################################################################

def c_write_header(fd):
    fd.write("#include <stdio.h>\n")
    fd.write("#include <stdint.h>\n")
    fd.write("#include <mpfr.h>\n")
    fd.write("#include <math.h>\n")
    fd.write("#include <assert.h>\n")
    fd.write("#include <fenv.h>\n")
    fd.write("\n")

    fd.write("int conforms(uint32_t a__i, uint32_t b__i)\n")
    fd.write("{\n")
    fd.write("  const float a = *((float*)&a__i);\n")
    fd.write("  const float b = *((float*)&b__i);\n")
    fd.write("  if (isnan(a) && isnan(b)) {\n")
    fd.write("    return 1;\n")
    fd.write("  } else {\n")
    fd.write("    return a__i == b__i;\n")
    fd.write("  }\n")
    fd.write("}\n")
    fd.write("\n")

    fd.write("int main(int argc, char **argv) {\n")
    fd.write("  int arg_id = 0;\n")
    fd.write("  mpfr_set_default_prec(24);\n")
    # Ok, now time for some special 'UP YOURS' brand MPFR exceptionalism.
    # "We don't need to use a sane value, nor do we need to document it
    # clearly because why wouldn't everyone just use the same."
    #
    # For reference the correct values for emin and emax are:
    #      Float32 Float64
    # emin    -148   -1073
    # emax     128    1024
    fd.write("  mpfr_set_emin(-148);\n")
    fd.write("  mpfr_set_emax(128);\n")

def c_write_footer(fd):
    fd.write("  return 0;\n")
    fd.write("}\n")

def c_write_input(fd, var):
    fd.write("\n")
    fd.write("  // input %s\n" % var)
    fd.write("  mpfr_t %s;\n" % var)
    fd.write("  uint32_t %s__i;\n" % var)
    fd.write("  arg_id += 1; assert (arg_id < argc);\n")
    fd.write('  sscanf(argv[arg_id], "%%x", &%s__i);\n' % var)
    fd.write("  const float %s__f = *((float*)&%s__i);\n" % (var, var))
    fd.write("  mpfr_init2(%s, 24);\n" % var)
    fd.write("  mpfr_set_flt(%s, %s__f, MPFR_RNDN);\n" % (var, var))

def c_write_output(fd, fp_var, fp_var_from_c = None):
    fd.write("  {\n")
    fd.write("    const float %s__f = mpfr_get_flt(%s, MPFR_RNDN);\n" %
             (fp_var, fp_var))
    fd.write("    const uint32_t %s__i = *((uint32_t*)&%s__f);\n" %
             (fp_var, fp_var))
    if fp_var_from_c is not None:
        fd.write("    if (!conforms(%s__i, %s)) {\n" % (fp_var, fp_var_from_c))
        fd.write('      fprintf(stderr, "from C   : %%08X [%%f]\\n", %s, *((float*)&%s));\n'
                 % (fp_var_from_c, fp_var_from_c))
        fd.write('      fprintf(stderr, "from MPFR: %%08X [%%f]\\n", %s__i, %s__f);\n'
                 % (fp_var, fp_var))
        fd.write("    }\n")
        fd.write("    assert(conforms(%s__i, %s));\n" % (fp_var, fp_var_from_c))
    fd.write('    printf("%%u\\n", %s__i);\n' % fp_var)
    fd.write("  }\n")

def c_write_mpfr_rm(rm):
    assert rm in VALID_MPFR_ROUNDING_MODES
    return {RM_RNE : "MPFR_RNDN",
            RM_RTZ : "MPFR_RNDZ",
            RM_RTP : "MPFR_RNDU",
            RM_RTN : "MPFR_RNDD"}[rm]

def c_write_c_rm(rm):
    assert rm in VALID_MPFR_ROUNDING_MODES
    return {RM_RNE : "FE_TONEAREST",
            RM_RTZ : "FE_TOWARDZERO",
            RM_RTP : "FE_UPWARD",
            RM_RTN : "FE_DOWNWARD"}[rm]

def c_write_mpfr_call(fd, function, rounding, result, *args):
    assert function.startswith("mpfr_")
    assert rounding in VALID_MPFR_ROUNDING_MODES or rounding is None
    assert len(args) >= 1

    fd.write("  {\n")
    if rounding is not None:
        fd.write("    int t;\n")
    fd.write("    mpfr_init2(%s, 24);\n" % result)
    fd.write("    mpfr_clear_flags();\n")
    fd.write("    ")
    if rounding is not None:
        fd.write("t = ")
    fd.write("%s(%s" % (function, result))
    for arg in args:
        fd.write(", %s" % arg)
    if rounding is not None:
        fd.write(", %s" % c_write_mpfr_rm(rounding))
    fd.write(");\n")
    if rounding is not None:
        fd.write("    t = mpfr_check_range(%s, t, %s);\n" %
                 (result, c_write_mpfr_rm(rounding)))
        fd.write("    mpfr_subnormalize(%s, t, %s);\n" %
                 (result, c_write_mpfr_rm(rounding)))
    fd.write("  }\n")

def c_write_sanity_check(fd, c_expr, rm):
    if rm not in VALID_MPFR_ROUNDING_MODES:
        return None

    result_name = "result_from_c"

    fd.write("\n")
    fd.write("  // sanity check C result is equal to MPFR result\n")
    fd.write("  uint32_t %s;\n" % result_name)
    fd.write("  {\n")
    if rm != RM_RNE:
        fd.write("    fesetround(%s);\n" % c_write_c_rm(rm))
    fd.write("    const float result__f = %s;\n" % c_expr)
    if rm != RM_RNE:
        fd.write("    fesetround(%s);\n" % c_write_c_rm(RM_RNE))
    fd.write("    %s = *((uint32_t*)&result__f);\n" % result_name)
    fd.write("  }\n")

    return result_name

def c_compile(c_file, x_file):
    assert os.path.isfile(c_file) and c_file.endswith(".c")
    # assert not os.path.exists(x_file)
    print "Compiling [%s]" % os.path.basename(x_file)
    rv = os.system("gcc -lmpfr -lm %s %s -o%s" % (GCC_FLAGS, c_file, x_file))
    if rv != 0:
        sys.exit(rv)

def c_run(x_file, *fp_args):
    assert os.path.isfile(x_file)
    assert len(fp_args) >= 1

    cmd = ["./" + x_file]
    for fp in fp_args:
        assert 0 <= fp.bv <= 0xFFFFFFFF
        cmd.append("0x%08X" % fp.bv)

    out = subprocess.check_output(cmd).strip().split()
    results = [MPF(8, 24, int(x)) for x in out]
    assert len(results) >= 1
    if len(results) == 1:
        return results[0]
    else:
        return results

##############################################################################
# FP evaluation and sanity checking
##############################################################################

def fp_eval_predicate(fp_ops, *args):
    assert fp_ops in FP_OPS
    op_arity = FP_OPS[fp_ops]["arity"]
    op_rnd   = FP_OPS[fp_ops]["rnd"]
    op_rt    = FP_OPS[fp_ops]["result"]
    assert op_rt == TYP_BOOL
    assert op_arity == len(args)
    assert not op_rnd

    mpf_fn = {
        "fp.isNormal"    : lambda x: x.isNormal(),
        "fp.isSubnormal" : lambda x: x.isSubnormal(),
        "fp.isZero"      : lambda x: x.isZero(),
        "fp.isInfinite"  : lambda x: x.isInfinite(),
        "fp.isNaN"       : lambda x: x.isNaN(),
        "fp.isNegative"  : lambda x: x.isNegative(),
        "fp.isPositive"  : lambda x: x.isPositive(),
        "fp.eq"          : lambda x, y: x == y,
        "fp.lt"          : lambda x, y: x < y,
        "fp.gt"          : lambda x, y: x > y,
        "fp.leq"         : lambda x, y: x <= y,
        "fp.geq"         : lambda x, y: x >= y,
        "smtlib.eq"      : smtlib_eq,
    }[fp_ops]

    return mpf_fn(*args)

FP_EVALUATION_PROGRAMS = {}

def fp_eval_function(fp_ops, rm, *args):
    assert fp_ops in FP_OPS
    op_arity = FP_OPS[fp_ops]["arity"]
    op_rnd   = FP_OPS[fp_ops]["rnd"]
    op_rt    = FP_OPS[fp_ops]["result"]
    assert op_rt == TYP_FLOAT
    assert op_arity == len(args)
    assert not op_rnd or rm in MPF.ROUNDING_MODES
    assert rm is None or rm in MPF.ROUNDING_MODES
    for arg in args[1:]:
        assert args[0].compatible(arg)

    results = {
        "MPF"    : None,
        "MPFR/C" : None,
    }

    precision_eb = args[0].w
    precision_sb = args[0].p

    # First we evaluate using our python floats
    mpf_fn = {
        "fp.abs"             : lambda x: abs(x),
        "fp.neg"             : lambda x: -x,
        "fp.sqrt"            : fp_sqrt,
        "fp.roundToIntegral" : fp_roundToIntegral,
        "fp.add"             : fp_add,
        "fp.sub"             : fp_sub,
        "fp.mul"             : fp_mul,
        "fp.div"             : fp_div,
        "fp.rem"             : fp_rem,
        "fp.min"             : fp_min,
        "fp.max"             : fp_max,
        "fp.fma"             : fp_fma,
    }[fp_ops]
    if op_rnd:
        results["MPF"] = mpf_fn(rm, *args)
    else:
        # Note we can throw the Unspecified exception here for min and max
        results["MPF"] = mpf_fn(*args)

    # Next we check C and MPFR results
    prog_name = "eval__" + fp_ops
    if op_rnd:
        prog_name += "__" + rm.lower()

    if prog_name not in FP_EVALUATION_PROGRAMS:
        mpfr_fn = {
            "fp.abs"             : "mpfr_abs",
            "fp.neg"             : "mpfr_neg",
            "fp.sqrt"            : "mpfr_sqrt",
            "fp.roundToIntegral" : ("mpfr_round"
                                    if rm == RM_RNA
                                    else "mpfr_rint"),
            "fp.add"             : "mpfr_add",
            "fp.sub"             : "mpfr_sub",
            "fp.mul"             : "mpfr_mul",
            "fp.div"             : "mpfr_div",
            "fp.rem"             : "mpfr_remainder",
            "fp.min"             : "mpfr_min",
            "fp.max"             : "mpfr_max",
            "fp.fma"             : "mpfr_fma",
        }[fp_ops]
        mpfr_rounding  = True
        mpfr_rm        = (rm if op_rnd else RM_RNE)
        mpfr_supported = True
        if rm == RM_RNA:
            if fp_ops == "fp.roundToIntegral":
                # We're using a special function here
                mpfr_rounding = False
                mpfr_rm       = None
            elif fp_ops == "fp.sqrt":
                # RNE and RNA do the same, see Theorem 19 in HoFPA
                mpfr_rm = RM_RNE
            else:
                # Otherwise, this operation is not supported by GNU/MPFR
                mpfr_supported = False

        if precision_eb == 8 and precision_sb == 24:
            c_type      = "float"
            c_supported = True
            c_indicator = "f"
        elif precision_eb == 11 and precision_sb == 53:
            c_type      = "double"
            c_supported = True
            c_indicator = ""
        else:
            c_type      = None
            c_supported = False
            c_indicator = ""
        def builtin(name, arity=1):
            args = ["%s__f"] * arity
            return "__builtin_%s%s(%s)" % (name, c_indicator, ", ".join(args))
        c_fn = {
            "fp.abs"             : builtin("fabs"),
            "fp.neg"             : "(-%s__f)",
            "fp.sqrt"            : builtin("sqrt"),
            "fp.roundToIntegral" : builtin("nearbyint"),
            "fp.add"             : "(%s__f + %s__f)",
            "fp.sub"             : "(%s__f - %s__f)",
            "fp.mul"             : "(%s__f * %s__f)",
            "fp.div"             : "(%s__f / %s__f)",
            "fp.rem"             : builtin("remainder", 2),
            "fp.min"             : builtin("fmin", 2),
            "fp.max"             : builtin("fmax", 2),
            "fp.fma"             : builtin("fma", 3),
        }[fp_ops]
        c_rounding = op_rnd
        c_rm       = rm
        if rm == RM_RNA:
            if fp_ops == "fp.sqrt":
                # RNE and RNA do the same, see Theorem 19 in HoFPA
                c_rm = RM_RNE
            else:
                # Otherwise, this operation is not supported by C
                c_supported = False

        # assert c_supported -> mpfr_supported
        assert not c_supported or mpfr_supported

        if mpfr_supported:
            # We do not yet have this program compiled...
            c_inputs = tuple(map(chr, map(lambda x: ord('x') + x,
                                          xrange(len(args)))))

            with open("tmp.c", "w") as fd:
                c_write_header(fd)
                for var in c_inputs:
                    c_write_input(fd, var)

                fd.write("\n")
                fd.write("  mpfr_t result;\n")
                c_write_mpfr_call(fd, mpfr_fn, mpfr_rm, "result", *c_inputs)

                if c_supported:
                    result_from_c = c_write_sanity_check(fd,
                                                         c_fn % c_inputs,
                                                         c_rm)
                else:
                    result_from_c = None

                c_write_output(fd, "result", result_from_c)
                c_write_footer(fd)

            c_compile("tmp.c", prog_name)
            FP_EVALUATION_PROGRAMS[prog_name] = prog_name

    if prog_name in FP_EVALUATION_PROGRAMS:
        results["MPFR/C"] = c_run(FP_EVALUATION_PROGRAMS[prog_name], *args)

        if not smtlib_eq(results["MPF"], results["MPFR/C"]):
            print "(%s %s" % (fp_ops, rm)
            for arg in args:
                print "  %s" % arg
            print ")"
            for k in results:
                print "%10s : %s" % (k, results[k])
            assert False

    return results["MPF"]

# a = MPF(8, 24)
# a.bv = 0x7F7FFFFF
# b = MPF(8, 24)
# b.bv = 0x007FFFFF
# c = MPF(8, 24)
# c.bv = 0
# fp_eval_function("fp.fma", RM_RTP, a, b, c)
