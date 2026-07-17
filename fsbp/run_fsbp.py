#!/usr/bin/env python3
# Validate the FSBP construction and compare optimal-node operators to standard
# Lobatto (LGL).
import math
from fsbp import (poly_space, exp_space, trig_space, fsbp_operator,
                  sbp_residuals, deriv_error, show)

def hdr(s): print("\n" + "="*74 + "\n  " + s + "\n" + "="*74)

# ── 1. Validation: polynomial space on [-1,1] must recover our LGL operators ──
hdr("1. POLYNOMIAL FSBP  (F = P3 on [-1,1])  ==  standard LGL / DGSEM")
op_lgl11 = fsbp_operator(poly_space(-1.0, 1.0, 3))
# reference LGL p=3 nodes & GLL weights used in the DGSEM
lgl_x = [-1.0, -1/math.sqrt(5), 1/math.sqrt(5), 1.0]
lgl_w = [1/6, 5/6, 5/6, 1/6]
print("  nodes  FSBP:", ["%+.10f" % v for v in op_lgl11["x"]])
print("  nodes  LGL :", ["%+.10f" % v for v in lgl_x])
print("  weights FSBP:", ["%.10f" % v for v in op_lgl11["w"]])
print("  weights GLL :", ["%.10f" % v for v in lgl_w])
print("  max|node diff|   = %.2e" % max(abs(op_lgl11["x"][i]-lgl_x[i]) for i in range(4)))
print("  max|weight diff| = %.2e" % max(abs(op_lgl11["w"][i]-lgl_w[i]) for i in range(4)))
c, s = sbp_residuals(op_lgl11)
print("  SBP checks:  |Df-f'| = %.2e   |Q+Q^T-B| = %.2e" % (c, s))

# ── 2. Validation: exponential space on [0,1]  ==  paper Example I ──
hdr("2. EXPONENTIAL FSBP  (F = span{1,x,e^x} on [0,1])  ==  paper Example I")
op_exp = fsbp_operator(exp_space())
paper_x = [0.0, 0.2956452974, 0.7423537958, 1.0]
paper_w = [0.0914828668, 0.4341375639, 0.3987262252, 0.0756533441]
paper_D = [[-5.465504277,  7.365125959, -2.802901094,  0.903279412],
           [-1.552003083,  0.0,          2.142484824, -0.590481741],
           [ 0.643091453, -2.332761387,  0.0,          1.689669934],
           [-1.092279411,  3.388486098, -8.905299859,  6.609093173]]
print("  nodes   FSBP :", ["%.10f" % v for v in op_exp["x"]])
print("  nodes   paper:", ["%.10f" % v for v in paper_x])
print("  weights FSBP :", ["%.10f" % v for v in op_exp["w"]])
print("  weights paper:", ["%.10f" % v for v in paper_w])
print("  max|node diff| = %.2e   max|weight diff| = %.2e" % (
      max(abs(op_exp["x"][i]-paper_x[i]) for i in range(4)),
      max(abs(op_exp["w"][i]-paper_w[i]) for i in range(4))))
print("  Differentiation matrix D (FSBP):"); show("D", op_exp["D"])
print("  max|D - D_paper| = %.2e" %
      max(abs(op_exp["D"][i][j]-paper_D[i][j]) for i in range(4) for j in range(4)))
c, s = sbp_residuals(op_exp)
print("  SBP checks:  |Df-f'| (f in F) = %.2e   |Q+Q^T-B| = %.2e" % (c, s))

# ── 3. Comparison: differentiation accuracy, exp-FSBP vs poly-LGL on [0,1] ──
hdr("3. DERIVATIVE ACCURACY  —  exp-FSBP  vs  standard Lobatto (both 4 nodes, [0,1])")
op_poly = fsbp_operator(poly_space(0.0, 1.0, 3))   # = LGL on [0,1]
tests = [
    ("e^x",        lambda x: math.exp(x),        lambda x: math.exp(x)),
    ("e^(1.7x)",   lambda x: math.exp(1.7*x),    lambda x: 1.7*math.exp(1.7*x)),
    ("cosh(2x-1)", lambda x: math.cosh(2*x-1),   lambda x: 2*math.sinh(2*x-1)),
    ("x^3",        lambda x: x**3,               lambda x: 3*x**2),
    ("x^5",        lambda x: x**5,               lambda x: 5*x**4),
    ("sin(3x)",    lambda x: math.sin(3*x),      lambda x: 3*math.cos(3*x)),
    ("1/(1+x)",    lambda x: 1/(1+x),            lambda x: -1/(1+x)**2),
]
print("  %-12s %16s %16s   %s" % ("f(x)", "poly-LGL err", "exp-FSBP err", "winner"))
print("  " + "-"*66)
for name, f, fp in tests:
    e_lgl = deriv_error(op_poly, f, fp)
    e_exp = deriv_error(op_exp,  f, fp)
    win = "exp-FSBP" if e_exp < e_lgl else "poly-LGL"
    print("  %-12s %16.3e %16.3e   %s" % (name, e_lgl, e_exp, win))
print("\n  (each operator is exact on functions in its own space:")
print("   poly-LGL exact for polynomials up to x^5; exp-FSBP exact for 1, x, e^x.)")

# ── 4. Trigonometric FSBP (5-node) vs standard Lobatto on oscillatory data ──
hdr("4. TRIG FSBP  (F = span{1,cos x,sin x,cos 2x,sin 2x} on [-1,1])  vs  LGL (P4)")
op_trig = fsbp_operator(trig_space(-1.0, 1.0))
op_p4   = fsbp_operator(poly_space(-1.0, 1.0, 4))     # standard 5-pt LGL
lgl5 = [-1.0, -math.sqrt(3.0/7.0), 0.0, math.sqrt(3.0/7.0), 1.0]
print("  nodes  trig-FSBP:", ["%+.7f" % v for v in op_trig["x"]])
print("  nodes  LGL (P4) :", ["%+.7f" % v for v in lgl5])
c, s = sbp_residuals(op_trig)
print("  trig SBP checks:  |Df-f'| (f in F) = %.2e   |Q+Q^T-B| = %.2e" % (c, s))
c, s = sbp_residuals(op_p4)
print("  LGL(P4) node err vs analytic = %.1e   SBP |Q+Q^T-B| = %.1e"
      % (max(abs(op_p4["x"][i]-lgl5[i]) for i in range(5)), s))
tests2 = [
    ("sin(2x)",   lambda x: math.sin(2*x),   lambda x: 2*math.cos(2*x)),
    ("cos(1.3x)", lambda x: math.cos(1.3*x), lambda x: -1.3*math.sin(1.3*x)),
    ("sin(3.5x)", lambda x: math.sin(3.5*x), lambda x: 3.5*math.cos(3.5*x)),
    ("e^x",       lambda x: math.exp(x),     lambda x: math.exp(x)),
    ("x^4",       lambda x: x**4,            lambda x: 4*x**3),
    ("x^7",       lambda x: x**7,            lambda x: 7*x**6),
]
print("\n  %-11s %16s %16s   %s" % ("f(x)", "LGL(P4) err", "trig-FSBP err", "winner"))
print("  " + "-"*63)
for name, f, fp in tests2:
    e_l = deriv_error(op_p4, f, fp)
    e_t = deriv_error(op_trig, f, fp)
    print("  %-11s %16.3e %16.3e   %s" % (name, e_l, e_t,
          "trig-FSBP" if e_t < e_l else "LGL(P4)"))
print("\n  => optimal nodes for a function space give an SBP operator that is")
print("     EXACT on that space; for polynomials it is precisely standard Lobatto.")
