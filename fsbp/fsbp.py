#!/usr/bin/env python3
# fsbp.py — Function-space SBP (FSBP) operators with optimal nodes.
#
# Implements Hale, Harley, Nchupang & Nordstrom, "Summation-by-parts operators
# for general function spaces: optimal nodes" (arXiv:2604.23306v1).
#
# For a function space F, the optimal (minimal-dimension) diagonal-norm SBP
# operator D = P^{-1} Q uses the generalized Gauss-Lobatto quadrature (GGLQ)
# nodes for G = (FF)'.  For F = polynomials this recovers standard Lobatto (LGL);
# for other spaces (here: exponential) the optimal nodes shift.
#
# Pure Python (no numpy) — small dense matrices only.

import math

# ─────────────────────────── tiny linear algebra ───────────────────────────
def mat(r, c, v=0.0): return [[v]*c for _ in range(r)]
def matmul(A, B):
    r, k, c = len(A), len(B), len(B[0])
    C = mat(r, c)
    for i in range(r):
        for m in range(k):
            a = A[i][m]
            if a == 0.0: continue
            Bm = B[m]
            Ci = C[i]
            for j in range(c): Ci[j] += a*Bm[j]
    return C
def transpose(A): return [list(col) for col in zip(*A)]

def solve(A, b):
    """Solve A x = b (square) via Gaussian elimination w/ partial pivoting."""
    n = len(A)
    M = [row[:] + [b[i]] for i, row in enumerate(A)]
    for col in range(n):
        p = max(range(col, n), key=lambda r: abs(M[r][col]))
        if abs(M[p][col]) < 1e-300: raise ZeroDivisionError("singular")
        M[col], M[p] = M[p], M[col]
        piv = M[col][col]
        for r in range(n):
            if r == col: continue
            f = M[r][col]/piv
            if f == 0.0: continue
            for j in range(col, n+1): M[r][j] -= f*M[col][j]
    return [M[i][n]/M[i][i] for i in range(n)]

def lstsq(A, b):
    """Least-squares min ||A x - b|| via normal equations A^T A x = A^T b."""
    At = transpose(A)
    AtA = matmul(At, A)
    Atb = [sum(At[i][k]*b[k] for k in range(len(b))) for i in range(len(At))]
    return solve(AtA, Atb)

# ─────────────────────────── function-space defs ───────────────────────────
# A space is a dict with:
#   F   : list of (f, f')      basis of F              (dim_F functions)
#   G   : list of (g, g')      basis of G=(FF)'        (dim 2n functions)
#   m   : list of moments  ∫_a^b g_j dx                (len 2n)
#   ab  : (a, b)
#   n   : G has dimension 2n  ->  (n+1)-point closed GGLQ

E = math.e
def poly_space(a, b, degF):
    """F = polynomials of degree <= degF.  G = (FF)' = P_{2*degF-1}."""
    F = [((lambda x, k=k: x**k), (lambda x, k=k: (0.0 if k == 0 else k*x**(k-1))))
         for k in range(degF+1)]
    dG = 2*degF                                   # dim G = 2n
    G = [((lambda x, k=k: x**k), (lambda x, k=k: (0.0 if k == 0 else k*x**(k-1))))
         for k in range(dG)]
    m = [(b**(k+1)-a**(k+1))/(k+1) for k in range(dG)]
    return dict(F=F, G=G, m=m, ab=(a, b), n=degF)

def exp_space():
    """Paper Example I:  F = span{1, x, e^x} on [0,1].
       (FF)' = span{1, x, e^x, x e^x, e^{2x}}  (dim 5, odd) -> augment with x^2."""
    a, b = 0.0, 1.0
    F = [((lambda x: 1.0),      (lambda x: 0.0)),
         ((lambda x: x),        (lambda x: 1.0)),
         ((lambda x: math.exp(x)), (lambda x: math.exp(x)))]
    # G = {1, x, x^2, e^x, x e^x, e^{2x}}  (dim 6 = 2n, n=3)
    G = [((lambda x: 1.0),          (lambda x: 0.0)),
         ((lambda x: x),            (lambda x: 1.0)),
         ((lambda x: x*x),          (lambda x: 2*x)),
         ((lambda x: math.exp(x)),  (lambda x: math.exp(x))),
         ((lambda x: x*math.exp(x)),(lambda x: math.exp(x)*(1+x))),
         ((lambda x: math.exp(2*x)),(lambda x: 2*math.exp(2*x)))]
    m = [1.0, 0.5, 1.0/3.0, E-1.0, 1.0, (E*E-1.0)/2.0]
    return dict(F=F, G=G, m=m, ab=(a, b), n=3)

def _fact(n):
    r = 1
    for k in range(2, n+1): r *= k
    return r

def _int_xk_ecx(i, c, a, b):
    """∫_a^b x^i e^{c x} dx  (closed form, exact for any c -> robust for small δ)."""
    def antideriv(x):
        s = 0.0
        for k in range(i+1):
            s += ((-1)**k) * (_fact(i)//_fact(i-k)) * (x**(i-k)) / (c**(k+1))
        return math.exp(c*x)*s
    return antideriv(b) - antideriv(a)

def _int_xk_exp(i, c, x0, a, b):
    """∫_a^b x^i e^{c(x-x0)} dx."""
    return math.exp(-c*x0) * _int_xk_ecx(i, c, a, b)

def exp_layer_space(delta, a=0.0, b=1.0, polydeg=1):
    """Boundary-layer-enriched space  F = span{1, x, ..., x^polydeg, φ},
       φ(x) = e^{(x-b)/δ}  (anchored at the layer endpoint b so φ ∈ (0,1]).
       G = (FF)' = {1..x^{2p-1}} ∪ {φ, xφ,.., x^p φ} ∪ {e^{2(x-b)/δ}}, dim 3p+2;
       augment with x^{2p} when that is odd (as Example I augments with x²)."""
    p = polydeg
    c1, c2 = 1.0/delta, 2.0/delta
    ephi  = lambda x: math.exp((x-b)*c1)
    # F basis (span{1,x,..,x^p, φ})
    F = [((lambda x, k=k: x**k), (lambda x, k=k: (0.0 if k == 0 else k*x**(k-1))))
         for k in range(p+1)]
    F.append((ephi, (lambda x: c1*math.exp((x-b)*c1))))
    # G basis + moments
    dim = 3*p + 2
    pmax = 2*p - 1
    if dim % 2 == 1: pmax = 2*p; dim += 1
    G, m = [], []
    for k in range(pmax+1):                                   # x^0 .. x^pmax
        G.append(((lambda x, k=k: x**k),
                  (lambda x, k=k: (0.0 if k == 0 else k*x**(k-1)))))
        m.append((b**(k+1)-a**(k+1))/(k+1))
    for k in range(p+1):                                      # x^k φ
        G.append(((lambda x, k=k: (x**k)*math.exp((x-b)*c1)),
                  (lambda x, k=k: ((k*x**(k-1) if k > 0 else 0.0)
                                   + (x**k)*c1)*math.exp((x-b)*c1))))
        m.append(_int_xk_exp(k, c1, b, a, b))
    G.append(((lambda x: math.exp((x-b)*c2)),                 # e^{2(x-b)/δ}
              (lambda x: c2*math.exp((x-b)*c2))))
    m.append(_int_xk_exp(0, c2, b, a, b))
    return dict(F=F, G=G, m=m, ab=(a, b), n=dim//2)

def trig_space(a=-1.0, b=1.0):
    """F = span{1, cos x, sin x, cos 2x, sin 2x}  (first two Fourier modes).
       G = (FF)' = span{cos kx, sin kx : k=1..4}  (dim 8 = 2n, n=4) -> 5-pt."""
    import math as _m
    F = [((lambda x: 1.0),           (lambda x: 0.0))]
    for k in (1, 2):
        F.append(((lambda x, k=k: _m.cos(k*x)), (lambda x, k=k: -k*_m.sin(k*x))))
        F.append(((lambda x, k=k: _m.sin(k*x)), (lambda x, k=k:  k*_m.cos(k*x))))
    G, m = [], []
    for k in range(1, 5):                                    # k = 1..4
        G.append(((lambda x, k=k: _m.cos(k*x)), (lambda x, k=k: -k*_m.sin(k*x))))
        m.append((_m.sin(k*b)-_m.sin(k*a))/k)                # ∫cos(kx)
        G.append(((lambda x, k=k: _m.sin(k*x)), (lambda x, k=k:  k*_m.cos(k*x))))
        m.append((-_m.cos(k*b)+_m.cos(k*a))/k)               # ∫sin(kx)
    return dict(F=F, G=G, m=m, ab=(a, b), n=4)

# ─────────────────────────── GGLQ nodes (quasi-Newton) ─────────────────────
def gglq_nodes(space, tol=1e-14, itmax=100, x_init=None):
    """Closed (n+1)-point generalized Gauss-Lobatto rule for G (Sec 3.2).
       Fix x0=a, xn=b; iterate interior nodes:  x_k += ∫σ_k / ∫η_k.
       x_init: optional warm-start nodes (used by δ-continuation)."""
    a, b = space["ab"]; n = space["n"]; G = space["G"]; m = space["m"]
    npts = n+1
    if x_init is not None:
        x = list(x_init)
    else:
        mid, half = 0.5*(a+b), 0.5*(b-a)                        # Chebyshev-Lobatto init
        x = [mid - half*math.cos(math.pi*k/n) for k in range(npts)]
    x[0], x[n] = a, b
    for _ in range(itmax):
        # Hermite-Vandermonde V̂ : rows = [values@x0..xn, derivs@x1..x_{n-1}]
        Vh = mat(2*n, 2*n)
        for j in range(2*n):
            g, gp = G[j]
            for k in range(npts):      Vh[k][j]      = g(x[k])          # value rows
            for k in range(1, n):      Vh[n+k][j]    = gp(x[k])         # deriv rows
        # solve V̂^T z = m ;  z[value row k]=∫η_k , z[deriv row]=∫σ_k
        z = solve(transpose(Vh), m)
        dmax = 0.0
        for k in range(1, n):
            dx = z[n+k]/z[k]                      # ∫σ_k / ∫η_k
            x[k] += dx
            dmax = max(dmax, abs(dx))
        if dmax < tol: break
    # weights = ∫η_i at converged nodes
    Vh = mat(2*n, 2*n)
    for j in range(2*n):
        g, gp = G[j]
        for k in range(npts): Vh[k][j]   = g(x[k])
        for k in range(1, n): Vh[n+k][j] = gp(x[k])
    z = solve(transpose(Vh), m)
    w = [z[k] for k in range(npts)]
    return x, w

# ─────────────────────────── FSBP operator D = P^{-1}Q ─────────────────────
def exp_layer_operator(delta, a=0.0, b=1.0, polydeg=1, dsafe=0.5):
    """Build the boundary-layer FSBP operator with δ-continuation: ramp δ from a
       well-conditioned dsafe down to the target, warm-starting the GGLQ nodes.
       Needed because the stiff exponential makes a cold-start Vandermonde
       singular for small δ (the paper's continuation remedy)."""
    ds, d = [], max(dsafe, delta)
    while d > delta*1.0001:
        ds.append(d); d *= 0.75
    ds.append(delta)
    x_init = None
    for dk in ds:
        sp = exp_layer_space(dk, a, b, polydeg)
        x_init, _ = gglq_nodes(sp, x_init=x_init)
    return fsbp_operator(sp, x_init=x_init)

def fsbp_operator(space, x_init=None):
    x, w = gglq_nodes(space, x_init=x_init)
    npts = len(x); Fb = space["F"]; dF = len(Fb)
    Fm  = [[Fb[i][0](x[k]) for i in range(dF)] for k in range(npts)]   # f_i(x_k)
    Fx  = [[Fb[i][1](x[k]) for i in range(dF)] for k in range(npts)]   # f_i'(x_k)
    P   = [[(w[i] if i == j else 0.0) for j in range(npts)] for i in range(npts)]
    B   = mat(npts, npts); B[0][0] = -1.0; B[npts-1][npts-1] = 1.0
    # R = P Fx - 1/2 B F   ;  solve S F = R for skew-symmetric S
    PFx = matmul(P, Fx); BF = matmul(B, Fm)
    R = [[PFx[i][j] - 0.5*BF[i][j] for j in range(dF)] for i in range(npts)]
    # unknowns s_{ij}, i<j  (S[i][j]=s, S[j][i]=-s)
    idx = [(i, j) for i in range(npts) for j in range(i+1, npts)]
    # (S F)[r][c] = sum_m S[r][m] F[m][c]  linear in s
    A = []; rhs = []
    for r in range(npts):
        for c in range(dF):
            row = [0.0]*len(idx)
            for p, (i, j) in enumerate(idx):
                if r == i: row[p] += Fm[j][c]      # S[i][j]=+s
                if r == j: row[p] -= Fm[i][c]      # S[j][i]=-s
            A.append(row); rhs.append(R[r][c])
    s = lstsq(A, rhs)
    S = mat(npts, npts)
    for p, (i, j) in enumerate(idx): S[i][j] = s[p]; S[j][i] = -s[p]
    Q = [[0.5*B[i][j] + S[i][j] for j in range(npts)] for i in range(npts)]
    D = [[Q[i][j]/w[i] for j in range(npts)] for i in range(npts)]     # P^{-1}Q
    return dict(x=x, w=w, D=D, Q=Q, B=B, Fm=Fm, Fx=Fx)

# ─────────────────────────── checks / reporting ────────────────────────────
def sbp_residuals(op):
    npts = len(op["x"]); D = op["D"]; Q = op["Q"]; B = op["B"]
    # (i) consistency  D f = f'  for f in F
    DF = matmul(D, op["Fm"])
    cons = max(abs(DF[k][i]-op["Fx"][k][i]) for k in range(npts) for i in range(len(op["Fm"][0])))
    # (iii) SBP  Q + Q^T = B
    sbp = max(abs(Q[i][j]+Q[j][i]-B[i][j]) for i in range(npts) for j in range(npts))
    return cons, sbp

def deriv_error(op, f, fp):
    x, D = op["x"], op["D"]
    fx = [f(xi) for xi in x]
    Df = [sum(D[i][j]*fx[j] for j in range(len(x))) for i in range(len(x))]
    return max(abs(Df[i]-fp(x[i])) for i in range(len(x)))

def show(name, M, wid=13):
    print(f"  {name}:")
    for row in M:
        print("   " + " ".join(f"{v:>{wid}.7f}" for v in row))
