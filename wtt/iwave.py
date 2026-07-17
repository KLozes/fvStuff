#!/usr/bin/env python3
"""Donoho 4th-order interpolating wavelet transform (Deslauriers-Dubuc cubic).

Lifting form of the interpolating transform:
  * scaling coefficients = the even samples (interpolating => point values),
  * detail coefficients  = odd sample  -  cubic (4-point) prediction from evens.
Interior predict weights are the DD cubic [-1/16, 9/16, 9/16, -1/16]; near the
boundary a shifted 4-point Lagrange stencil keeps 4th order.  Perfect
reconstruction is exact (inverse just adds the same prediction back).
2-D is separable (rows then columns); multilevel recurses on the LL band, giving
the classic Mallat "squares" layout.
"""
import numpy as np


def _lagrange_weights(x, nodes):
    w = np.ones(len(nodes))
    for k in range(len(nodes)):
        for m in range(len(nodes)):
            if m != k:
                w[k] *= (x - nodes[m]) / (nodes[k] - nodes[m])
    return w


def predict_matrix(n_even):
    """P (n_even x n_even): predicts odd-sample values at positions 2j+1 from the
       even samples at positions 2k, using an order-4 (clamped) Lagrange stencil."""
    P = np.zeros((n_even, n_even))
    order = min(4, n_even)
    even_pos = 2 * np.arange(n_even)
    for j in range(n_even):                       # odd sample at position 2j+1
        lo = min(max(j - 1, 0), n_even - order)   # window of `order` consecutive evens
        idx = np.arange(lo, lo + order)
        P[j, idx] = _lagrange_weights(2 * j + 1, even_pos[idx])
    return P


_PCACHE = {}
def _P(n_even):
    if n_even not in _PCACHE:
        _PCACHE[n_even] = predict_matrix(n_even)
    return _PCACHE[n_even]


# ── 1-D one-level transform along an axis ──────────────────────────────────
def _fwd_axis(a, axis):
    a = np.moveaxis(a, axis, -1)
    even, odd = a[..., 0::2], a[..., 1::2]
    P = _P(even.shape[-1])
    detail = odd - even @ P.T                      # predict along last axis
    s = np.moveaxis(even,   -1, axis)
    d = np.moveaxis(detail, -1, axis)
    return s, d

def _inv_axis(s, d, axis):
    s = np.moveaxis(s, axis, -1); d = np.moveaxis(d, axis, -1)
    P = _P(s.shape[-1])
    odd = d + s @ P.T
    n = s.shape[-1] * 2
    out = np.empty(s.shape[:-1] + (n,))
    out[..., 0::2] = s; out[..., 1::2] = odd
    return np.moveaxis(out, -1, axis)


# ── 2-D multilevel transform (Mallat layout via (LL, [details])) ───────────
def iwt2(img, levels):
    """Return (LL, details) where details[l] = (LH, HL, HH) for level l (0=finest)."""
    LL = img.astype(float)
    details = []
    for _ in range(levels):
        L, H  = _fwd_axis(LL, axis=1)              # transform columns  -> low/high
        LL_, LH = _fwd_axis(L, axis=0)             # transform rows of L
        HL, HH  = _fwd_axis(H, axis=0)             # transform rows of H
        details.append((LH, HL, HH))
        LL = LL_
    return LL, details

def iiwt2(LL, details):
    for (LH, HL, HH) in reversed(details):
        L = _inv_axis(LL, LH, axis=0)
        H = _inv_axis(HL, HH, axis=0)
        LL = _inv_axis(L, H, axis=1)
    return LL


# ── 3-D multilevel transform ───────────────────────────────────────────────
def iwt3(vol, levels):
    """Return (LLL, details) with details[l] = 7 detail octants for level l."""
    LLL = vol.astype(float)
    details = []
    for _ in range(levels):
        L, H = _fwd_axis(LLL, 2)
        LL, LH = _fwd_axis(L, 1); HL, HH = _fwd_axis(H, 1)
        LLL_, LLH = _fwd_axis(LL, 0); LHL, LHH = _fwd_axis(LH, 0)
        HLL, HLH = _fwd_axis(HL, 0); HHL, HHH = _fwd_axis(HH, 0)
        details.append((LLH, LHL, LHH, HLL, HLH, HHL, HHH))
        LLL = LLL_
    return LLL, details

def iiwt3(LLL, details):
    for (LLH, LHL, LHH, HLL, HLH, HHL, HHH) in reversed(details):
        LL = _inv_axis(LLL, LLH, 0); LH = _inv_axis(LHL, LHH, 0)
        HL = _inv_axis(HLL, HLH, 0); HH = _inv_axis(HHL, HHH, 0)
        L = _inv_axis(LL, LH, 1); H = _inv_axis(HL, HH, 1)
        LLL = _inv_axis(L, H, 2)
    return LLL


def _selftest():
    rng = np.random.default_rng(0)
    N = 64
    img = rng.standard_normal((N, N))
    LL, det = iwt2(img, 4)
    rec = iiwt2(LL, det)
    pr = np.max(np.abs(rec - img))
    # order-4: details must vanish for a bicubic polynomial
    x = np.linspace(0, 1, N)
    X, Y = np.meshgrid(x, x)
    poly = 1 + 2*X - 3*Y + X*Y + X**2 - Y**3 + X**2*Y
    _, dpoly = iwt2(poly, 3)
    dmax = max(np.max(np.abs(b)) for lev in dpoly for b in lev)
    print(f"[iwave] 2D perfect-reconstruction max err = {pr:.2e}")
    print(f"[iwave] 2D bicubic detail magnitude (should ~0) = {dmax:.2e}")
    # 3D
    M = 32
    vol = rng.standard_normal((M, M, M))
    L3, d3 = iwt3(vol, 3)
    pr3 = np.max(np.abs(iiwt3(L3, d3) - vol))
    z = np.linspace(0, 1, M)
    Xx, Yy, Zz = np.meshgrid(z, z, z, indexing="ij")
    tri = 1 + Xx - 2*Yy + Zz + Xx*Yy - Zz**3 + Xx**2*Zz
    _, dt3 = iwt3(tri, 2)
    dmax3 = max(np.max(np.abs(b)) for lev in dt3 for b in lev)
    print(f"[iwave] 3D perfect-reconstruction max err = {pr3:.2e}")
    print(f"[iwave] 3D tricubic detail magnitude (should ~0) = {dmax3:.2e}")


if __name__ == "__main__":
    _selftest()
