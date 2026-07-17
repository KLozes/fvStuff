#!/usr/bin/env python3
# boundary_layer.py — Does a function-tailored (exponential) FSBP operator resolve
# a boundary layer that same-size polynomial-Lobatto cannot, with the layer
# thickness DISCOVERED from the solution (no a-priori delta)?
#
# Model (classic 1-D convection-diffusion outflow layer):
#     u' - eps*u'' = 0  on [0,1],  u(0)=0, u(1)=1
#     exact:  u(x) = (e^{x/eps}-1)/(e^{1/eps}-1),  layer width delta=eps at x=1.
#
# NOTE: BCs are imposed strongly (collocation with the FSBP differentiation
# matrix).  This isolates the *resolution/accuracy* benefit of the optimal nodes;
# SBP-SAT energy stability is a separate property, inherited from the FSBP
# framework (the operators satisfy Q+Q^T=B) and out of scope here.
# Keep eps >= ~0.01 so the stiff exponential stays well-conditioned (float64).

import math
from fsbp import fsbp_operator, poly_space, exp_layer_operator, matmul, solve

def u_exact(x, eps):
    return (math.exp(x/eps) - 1.0) / (math.exp(1.0/eps) - 1.0)

def solve_bvp(op, eps):
    """ -eps u'' + u' = 0, u(0)=0, u(1)=1  ->  L = D - eps*(D@D), strong Dirichlet."""
    x, D = op["x"], op["D"]; n = len(x)
    D2 = matmul(D, D)
    L = [[D[i][j] - eps*D2[i][j] for j in range(n)] for i in range(n)]
    b = [0.0]*n
    L[0]   = [1.0 if j == 0   else 0.0 for j in range(n)]; b[0]   = 0.0   # u(0)=0
    L[n-1] = [1.0 if j == n-1 else 0.0 for j in range(n)]; b[n-1] = 1.0   # u(1)=1
    return solve(L, b)

def nodal_error(op, U, eps):
    x = op["x"]
    return max(abs(U[i] - u_exact(x[i], eps)) for i in range(len(x)))

def adaptive_delta_solve(eps, polydeg=1, d0=0.3, tol=1e-6, itmax=30):
    """Discover the layer thickness from the solution: phi'=(1/delta)phi ->
       u/u' ~ delta at the layer.  Estimate delta = u(1)/u'(1), rebuild, repeat."""
    delta = d0
    hist = []
    for _ in range(itmax):
        op = exp_layer_operator(delta, 0.0, 1.0, polydeg)
        U  = solve_bvp(op, eps)
        x, D = op["x"], op["D"]
        n = len(x)
        DU1 = sum(D[n-1][j]*U[j] for j in range(n))          # u'(1)
        d_new = abs(U[n-1]/DU1) if abs(DU1) > 1e-30 else delta
        hist.append(d_new)
        if abs(d_new - delta) < tol:
            delta = d_new; break
        delta = d_new
    return op, U, delta, hist

def lagrange_interp(nodes, vals, x):
    s = 0.0
    for k in range(len(nodes)):
        t = vals[k]
        for m in range(len(nodes)):
            if m != k: t *= (x-nodes[m])/(nodes[k]-nodes[m])
        s += t
    return s

def write_svg(fname="bl_profile.svg", eps=0.05, polydeg=1, degF=3):
    """exact vs polynomial-LGL vs exp-FSBP (each a cubic through its 4 nodes)."""
    op_p = fsbp_operator(poly_space(0.0, 1.0, degF)); U_p = solve_bvp(op_p, eps)
    op_e, U_e, dfound, _ = adaptive_delta_solve(eps, polydeg)
    W, H, ml, mr, mt, mb = 720, 460, 64, 140, 34, 52
    def Xs(x): return ml + x*(W-ml-mr)
    def Ys(y): return H-mb - y*(H-mt-mb)
    fine = [i/500.0 for i in range(501)]
    def curve(fn, col, wd=2.5, dash=""):
        pts = " ".join("%.1f,%.1f" % (Xs(x), Ys(fn(x))) for x in fine)
        d = f' stroke-dasharray="{dash}"' if dash else ""
        return f'<polyline fill="none" stroke="{col}" stroke-width="{wd}"{d} points="{pts}"/>'
    def marks(op, col):
        return "".join('<circle cx="%.1f" cy="%.1f" r="4" fill="%s"/>' %
                       (Xs(x), Ys(u_exact(x, eps)), col) for x in op["x"])
    s = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" '
         f'font-family="sans-serif" font-size="13">',
         f'<rect width="{W}" height="{H}" fill="white"/>',
         f'<line x1="{ml}" y1="{Ys(0)}" x2="{Xs(1)}" y2="{Ys(0)}" stroke="#888"/>',
         f'<line x1="{ml}" y1="{Ys(0)}" x2="{ml}" y2="{mt}" stroke="#888"/>',
         f'<text x="{Xs(0.5)}" y="{H-16}" text-anchor="middle">x</text>',
         f'<text x="{Xs(0.5)}" y="20" text-anchor="middle">'
         f'convection-diffusion layer,  eps={eps}  (delta found = {dfound:.3f})</text>',
         curve(lambda x: u_exact(x, eps), "#111", 3),
         curve(lambda x: lagrange_interp(op_p["x"], U_p, x), "#d62728", 2),
         curve(lambda x: lagrange_interp(op_e["x"], U_e, x), "#1f77b4", 2),
         marks(op_p, "#d62728"), marks(op_e, "#1f77b4")]
    leg = [("exact", "#111"), (f"polynomial LGL ({degF+1} nodes)", "#d62728"),
           (f"exp-FSBP ({degF+1} nodes)", "#1f77b4")]
    for i, (lab, col) in enumerate(leg):
        yy = mt + 8 + i*22
        s.append(f'<line x1="{W-mr+6}" y1="{yy}" x2="{W-mr+30}" y2="{yy}" '
                 f'stroke="{col}" stroke-width="3"/>')
        s.append(f'<text x="{W-mr+34}" y="{yy+4}" font-size="11">{lab}</text>')
    s.append("</svg>")
    open(fname, "w").write("\n".join(s))
    print("  wrote %s  (eps=%.2f, %d nodes each)" % (fname, eps, degF+1))

def main():
    print("="*78)
    print("  Boundary layer  u' - eps u'' = 0,  u(0)=0, u(1)=1   (layer width = eps)")
    print("  max nodal error  |U_k - u_exact(x_k)|   (each method at its own nodes)")
    print("="*78)
    configs = [("4-node", 1, 3), ("5-node", 2, 4)]      # (label, exp polydeg, poly degF)
    for label, polydeg, degF in configs:
        print("\n  %s operators   (exp-FSBP polydeg=%d  vs  polynomial LGL P%d)" %
              (label, polydeg, degF))
        print("  %-6s %14s %16s %14s" % ("eps", "poly-LGL err", "exp-FSBP err",
                                         "delta found"))
        print("  " + "-"*54)
        for eps in [1.0, 0.1, 0.05, 0.02]:
            op_p = fsbp_operator(poly_space(0.0, 1.0, degF))
            e_p  = nodal_error(op_p, solve_bvp(op_p, eps), eps)
            op_e, U_e, dfound, _ = adaptive_delta_solve(eps, polydeg)
            e_e  = nodal_error(op_e, U_e, eps)
            print("  %-6.2f %14.2e %16.2e %10.4f  (eps=%.2f)" %
                  (eps, e_p, e_e, dfound, eps))
    print("\n  => exp-FSBP: adaptive delta -> ~eps, and it resolves the layer with")
    print("     the same node count where the polynomial Lobatto operator smears it.")
    write_svg("bl_profile.svg", eps=0.05, polydeg=1, degF=3)

if __name__ == "__main__":
    main()
