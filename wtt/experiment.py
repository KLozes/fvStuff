#!/usr/bin/env python3
"""Wavelet x tensor-train compression experiment.

Pipeline under test (COMBO):
  1. Donoho 4th-order interpolating wavelet transform (Mallat squares).
  2. Hard-threshold small DETAIL coefficients to zero.
  3. QTT-decompose EACH thresholded detail sub-band (truncate small singular
     values).  Keep the coarse scaling (LL) band.
Compared against:
  * PLAIN WAVELET : threshold detail coeffs, keep the nonzeros (no TT).
  * PLAIN QTT     : QTT-decompose the whole image.

"Parameters" = number of stored float values (nonzero coeffs for wavelet,
core entries for TT).  Error = relative L2 (Frobenius).  For each method we
sweep its accuracy knob to trace a params-vs-error curve, then read off the
parameter count at matched error levels.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from iwave import iwt2, iiwt2
from tt import qtt_compress

LEVELS = 5


# ── test signals (N x N, N = 2^k) ──────────────────────────────────────────
def signals(N=256):
    x = np.linspace(0, 1, N, endpoint=False)
    X, Y = np.meshgrid(x, x)
    sig = {}
    # 1. smooth: sum of Gaussian bumps
    g = np.zeros((N, N))
    for (cx, cy, s, a) in [(.3, .3, .12, 1), (.7, .6, .08, .8), (.5, .8, .15, .6)]:
        g += a*np.exp(-((X-cx)**2+(Y-cy)**2)/(2*s*s))
    sig["smooth"] = g
    # 2. oscillatory: 2-D chirp (band-limited but high-frequency)
    sig["oscillatory"] = np.sin(2*np.pi*(6*X + 10*Y*Y)) * np.cos(2*np.pi*8*X*Y)
    # 3. piecewise: sharp geometric edges (disk + rectangle)
    pw = np.zeros((N, N))
    pw[((X-.4)**2 + (Y-.55)**2) < .22**2] = 1.0
    pw[(X > .55) & (X < .85) & (Y > .15) & (Y < .5)] = -0.7
    sig["piecewise"] = pw
    # 4. mixed: smooth + edges + texture
    sig["mixed"] = (0.6*g + pw
                    + 0.15*np.sin(2*np.pi*20*X)*np.sin(2*np.pi*18*Y))
    return sig


def rel_err(a, b):
    return float(np.linalg.norm(a-b) / (np.linalg.norm(b) + 1e-300))


# ── the three methods:  knob -> (params, error) ────────────────────────────
def plain_wavelet(img, tau):
    LL, details = iwt2(img, LEVELS)
    params = LL.size
    newdet = []
    for band in details:
        keep = []
        for b in band:
            bt = np.where(np.abs(b) < tau, 0.0, b)
            params += int(np.count_nonzero(bt))
            keep.append(bt)
        newdet.append(tuple(keep))
    rec = iiwt2(LL, newdet)
    return params, rel_err(rec, img)

def plain_qtt(img, eps):
    Ahat, p = qtt_compress(img, eps)
    return p, rel_err(Ahat, img)

def combo(img, tau, tt_eps, threshold=True):
    LL, details = iwt2(img, LEVELS)
    params = LL.size
    newdet = []
    for band in details:
        keep = []
        for b in band:
            bt = np.where(np.abs(b) < tau, 0.0, b) if threshold else b
            bhat, p = qtt_compress(bt, tt_eps)              # TT each detail band
            params += p
            keep.append(bhat)
        newdet.append(tuple(keep))
    rec = iiwt2(LL, newdet)
    return params, rel_err(rec, img)


# ── sweep helpers ──────────────────────────────────────────────────────────
def sweep(fn, knobs):
    pts = []
    for k in knobs:
        try:
            p, e = fn(k)
            if p > 0: pts.append((e, p))
        except Exception:
            pass
    pts.sort()
    return pts                                              # sorted by error

def params_at(pts, target):
    """Smallest param count achieving error <= target (from the swept curve)."""
    good = [p for (e, p) in pts if e <= target]
    return min(good) if good else None


def main():
    sig = signals(256)
    scale = {name: np.sqrt((im**2).mean()) for name, im in sig.items()}
    taus = np.concatenate([[0.0], np.logspace(-3, 0.2, 22)])
    epss = np.logspace(-3.2, -0.2, 22)
    targets = [1e-1, 3e-2, 1e-2, 3e-3, 1e-3]

    curves = {}
    print("="*86)
    print("  Wavelet(Donoho DD-4) x Tensor-Train compression   —  N=256, %d levels" % LEVELS)
    print("  parameters (stored values) to reach a target relative-L2 error")
    print("="*86)
    for name, im in sig.items():
        s = scale[name]
        cw = sweep(lambda t: plain_wavelet(im, t*s),           taus)
        cq = sweep(lambda e: plain_qtt(im, e),                 epss)
        cc = sweep(lambda e: combo(im, e*s, tt_eps=e),         epss)  # threshold+TT, one knob
        curves[name] = (cw, cq, cc)
        print("\n  signal: %-11s   (%d values total)" % (name, im.size))
        print("  %-10s %14s %14s %14s   %s" %
              ("err<=", "plain-wavelet", "plain-QTT", "wave+TT", "best"))
        print("  " + "-"*72)
        for tgt in targets:
            pw, pq, pc = (params_at(cw, tgt), params_at(cq, tgt), params_at(cc, tgt))
            vals = {"plain-wavelet": pw, "plain-QTT": pq, "wave+TT": pc}
            best = min((v, k) for k, v in vals.items() if v is not None)[1] \
                   if any(v is not None for v in vals.values()) else "-"
            fmt = lambda v: ("%14d" % v) if v is not None else ("%14s" % "-")
            print("  %-10.0e %s %s %s   %s" % (tgt, fmt(pw), fmt(pq), fmt(pc), best))

    # ── why: TT cost of a sparse thresholded band vs a smooth band ────────
    print("\n" + "="*86)
    print("  DIAGNOSTIC — TT is cheap for smooth/low-rank data, expensive for sparse")
    print("="*86)
    print("  %-24s %10s %12s %12s" % ("array (near-lossless QTT)", "size", "nonzeros", "TT params"))
    print("  " + "-"*62)
    for name in ("smooth", "piecewise"):
        im = sig[name]; s = scale[name]
        LL, details = iwt2(im, LEVELS)
        b = details[1][2]                         # a mid-level HH detail band
        bt = np.where(np.abs(b) < 3e-2*s, 0.0, b) # thresholded (sparse)
        _, p_sparse = qtt_compress(bt, 1e-6)
        _, p_full   = qtt_compress(b, 1e-6)       # un-thresholded (structured)
        print("  %-24s %10d %12d %12d" %
              (name+" detail (thresholded)", b.size, int(np.count_nonzero(bt)), p_sparse))
        print("  %-24s %10d %12d %12d" %
              (name+" detail (raw)",         b.size, int(np.count_nonzero(b)),  p_full))
    _, pfull = qtt_compress(sig["smooth"], 1e-2)
    print("  %-24s %10d %12s %12d   (QTT DOES compress the global smooth image)"
          % ("smooth FULL image", sig["smooth"].size, "-", pfull))
    print("  => wavelet thresholding already makes detail bands sparse; re-encoding")
    print("     sparse data as TT costs MORE than storing the few nonzeros.")

    # ── plot params-vs-error curves ───────────────────────────────────────
    fig, axs = plt.subplots(2, 2, figsize=(11, 8))
    for ax, name in zip(axs.ravel(), sig):
        cw, cq, cc = curves[name]
        for pts, lab, col in [(cw, "plain wavelet", "#d62728"),
                              (cq, "plain QTT", "#2ca02c"),
                              (cc, "wavelet + TT", "#1f77b4")]:
            if pts:
                e = [p[0] for p in pts]; p = [p[1] for p in pts]
                ax.loglog(p, e, "o-", ms=3, color=col, label=lab)
        ax.set_title(name); ax.set_xlabel("parameters"); ax.set_ylabel("rel. L2 error")
        ax.grid(True, which="both", alpha=.3); ax.legend(fontsize=8)
    fig.suptitle("Compression: Donoho-4 wavelet vs QTT vs wavelet+TT (N=256)")
    fig.tight_layout()
    fig.savefig("wtt_tradeoff.png", dpi=110)
    print("\n  saved wtt_tradeoff.png")


if __name__ == "__main__":
    main()
