#!/usr/bin/env python3
"""Wavelet vs tensor-train compression in 3-D.

Unlike 2-D (where a "TT" of a matrix is just an SVD), a 3-D array has a genuine
3-core tensor train.  We compare:
  * PLAIN WAVELET : 3-D Donoho-4 interpolating wavelet, threshold detail coeffs.
  * PLAIN TT      : 3-core tensor train of the N x N x N array (2 SVDs).
  * QTT           : quantized TT (bit-interleaved i,j,k, coarse-to-fine).
Parameters = stored float values; error = relative L2.  Knob swept per method.
"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from iwave import iwt3, iiwt3
from tt import tt3_compress, qtt_compress3

N = int(sys.argv[1]) if len(sys.argv) > 1 else 128
LEVELS = int(round(np.log2(N))) - 2          # coarsest LLL band = 4x4x4
MAXRANK = 256 if N >= 256 else (512 if N >= 128 else None)   # cap TT/QTT ranks


def signals(N=64):
    """Complex, higher-rank 3-D fields."""
    rng = np.random.default_rng(7)
    x = np.linspace(0, 1, N, endpoint=False)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    nrm = lambda a: a / np.abs(a).max()
    sig = {}
    # 1. turbulence: isotropic Gaussian random field with a k^{-5/3} energy spectrum
    k = np.fft.fftfreq(N) * N
    KX, KY, KZ = np.meshgrid(k, k, k, indexing="ij")
    kk = np.sqrt(KX**2 + KY**2 + KZ**2); kk[0, 0, 0] = 1.0
    F = np.fft.fftn(rng.standard_normal((N, N, N))) * kk**(-5.0/6.0)
    sig["turbulence"] = nrm(np.fft.ifftn(F).real)
    # 2. multiwave: sum of 40 random (non-separable) plane waves
    mw = np.zeros((N, N, N))
    for _ in range(40):
        kx, ky, kz = rng.integers(-8, 9, 3); ph = rng.uniform(0, 2*np.pi)
        mw += rng.uniform(.3, 1.)*np.cos(2*np.pi*(kx*X+ky*Y+kz*Z)+ph)
    sig["multiwave"] = nrm(mw)
    # 3. structured: multiscale Gaussian blobs + thin curved shells (surfaces)
    st = np.zeros((N, N, N))
    for _ in range(12):
        c = rng.uniform(.15, .85, 3); s = rng.uniform(.03, .15)
        st += rng.uniform(.4, 1.)*np.exp(-((X-c[0])**2+(Y-c[1])**2+(Z-c[2])**2)/(2*s*s))
    for _ in range(3):
        c = rng.uniform(.3, .7, 3); r = rng.uniform(.12, .25)
        rr = np.sqrt((X-c[0])**2+(Y-c[1])**2+(Z-c[2])**2)
        st += 0.6*(np.abs(rr-r) < 0.02)
    sig["structured"] = nrm(st)
    # 4. cusps: sum of 1/(r+eps) singular fields
    cu = np.zeros((N, N, N))
    for _ in range(6):
        c = rng.uniform(.2, .8, 3)
        cu += 1.0/(np.sqrt((X-c[0])**2+(Y-c[1])**2+(Z-c[2])**2) + 0.02)
    sig["cusps"] = nrm(cu)
    return sig


def rel_err(a, b): return float(np.linalg.norm(a-b)/(np.linalg.norm(b)+1e-300))

def plain_wavelet(vol, tau):
    LLL, details = iwt3(vol, LEVELS)
    params = LLL.size
    nd = []
    for band in details:
        keep = []
        for b in band:
            bt = np.where(np.abs(b) < tau, 0.0, b)
            params += int(np.count_nonzero(bt)); keep.append(bt)
        nd.append(tuple(keep))
    return params, rel_err(iiwt3(LLL, nd), vol)

def plain_tt(vol, eps):
    Ah, p = tt3_compress(vol, eps, MAXRANK);  return p, rel_err(Ah, vol)

def qtt(vol, eps):
    Ah, p = qtt_compress3(vol, eps, MAXRANK); return p, rel_err(Ah, vol)


def sweep(fn, knobs):
    pts = []
    for k in knobs:
        try:
            p, e = fn(k)
            if p > 0: pts.append((e, p))
        except Exception:
            pass
    pts.sort(); return pts

def params_at(pts, tgt):
    good = [p for (e, p) in pts if e <= tgt]
    return min(good) if good else None


def main():
    sig = signals(N)
    scale = {n: np.sqrt((im**2).mean()) for n, im in sig.items()}
    npts = 8 if N >= 256 else 12
    taus = np.concatenate([[0.0], np.logspace(-6, 0.3, 16)])   # fine: reach tight err
    epss = np.logspace(-3.2, -0.3, npts)
    targets = [1e-1, 3e-2, 1e-2, 3e-3, 1e-3]

    curves = {}
    print("="*84)
    print("  Wavelet vs Tensor-Train in 3-D    N=%d  (%d^3 voxels), %d levels" % (N, N, LEVELS))
    print("  parameters (stored values) to reach a target relative-L2 error")
    print("="*84)
    for name, im in sig.items():
        s = scale[name]
        cw = sweep(lambda t: plain_wavelet(im, t*s), taus)
        ct = sweep(lambda e: plain_tt(im, e),        epss)
        cq = [] if name == "turbulence" else sweep(lambda e: qtt(im, e), epss)  # skip incompressible QTT
        curves[name] = (cw, ct, cq)
        print("\n  signal: %-11s" % name)
        print("  %-10s %14s %14s %14s   %s" %
              ("err<=", "plain-wavelet", "plain-TT(3)", "QTT", "best"))
        print("  " + "-"*72)
        for tgt in targets:
            pw, pt, pq = params_at(cw, tgt), params_at(ct, tgt), params_at(cq, tgt)
            vals = {"plain-wavelet": pw, "plain-TT": pt, "QTT": pq}
            best = min((v, k) for k, v in vals.items() if v is not None)[1] \
                   if any(v is not None for v in vals.values()) else "-"
            fmt = lambda v: ("%14d" % v) if v is not None else ("%14s" % "-")
            print("  %-10.0e %s %s %s   %s" % (tgt, fmt(pw), fmt(pt), fmt(pq), best))

    fig, axs = plt.subplots(2, 2, figsize=(11, 8))
    for ax, name in zip(axs.ravel(), sig):
        cw, ct, cq = curves[name]
        for pts, lab, col in [(cw, "plain wavelet", "#d62728"),
                              (ct, "plain TT (3-core)", "#9467bd"),
                              (cq, "QTT", "#2ca02c")]:
            if pts:
                e = [p[0] for p in pts]; p = [p[1] for p in pts]
                ax.loglog(p, e, "o-", ms=3, color=col, label=lab)
        ax.set_title(name); ax.set_xlabel("parameters"); ax.set_ylabel("rel. L2 error")
        ax.grid(True, which="both", alpha=.3); ax.legend(fontsize=8)
    fig.suptitle("3-D compression: Donoho-4 wavelet vs 3-core TT vs QTT (N=%d)" % N)
    fig.tight_layout(); fig.savefig("wtt3d_tradeoff_N%d.png" % N, dpi=110)
    print("\n  saved wtt3d_tradeoff_N%d.png" % N)


if __name__ == "__main__":
    main()
