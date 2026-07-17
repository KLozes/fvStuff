#!/usr/bin/env python3
"""How does compression cost scale with 3-D resolution N, at fixed accuracy?
   Expectation:  wavelet ~ const (feature-count), 3-core TT ~ O(N) (size-N mode),
   QTT ~ O(log N).  Params to reach relative-L2 error ~1e-3, for smooth &
   oscillatory 3-D fields, across N = 32..256."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from iwave import iwt3, iiwt3
from tt import tt3_compress, qtt_compress3

TGT = 1e-3

def rel_err(a, b): return float(np.linalg.norm(a-b)/(np.linalg.norm(b)+1e-300))

def make(name, N):
    x = np.linspace(0, 1, N, endpoint=False)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    if name == "smooth":
        g = np.zeros((N, N, N))
        for (cx, cy, cz, s, a) in [(.3,.3,.4,.12,1), (.7,.6,.5,.09,.8), (.5,.7,.3,.14,.6)]:
            g += a*np.exp(-((X-cx)**2+(Y-cy)**2+(Z-cz)**2)/(2*s*s))
        return g
    return np.sin(2*np.pi*(4*X+5*Y+3*Z))*np.cos(2*np.pi*3*X*Y)

def wavelet_params(vol, levels):
    """Bisect the detail threshold to reach ~TGT relative error; return #params."""
    s = np.sqrt((vol**2).mean())
    LLL, details = iwt3(vol, levels)
    def eval_tau(tau):
        p = LLL.size; nd = []
        for band in details:
            keep = []
            for b in band:
                bt = np.where(np.abs(b) < tau, 0.0, b)
                p += int(np.count_nonzero(bt)); keep.append(bt)
            nd.append(tuple(keep))
        return p, rel_err(iiwt3(LLL, nd), vol)
    lo, hi = 1e-8*s, 2.0*s                       # small tau -> low err
    for _ in range(40):
        mid = np.sqrt(lo*hi)
        _, e = eval_tau(mid)
        if e > TGT: hi = mid
        else:       lo = mid
    return eval_tau(lo)[0]

def main():
    Ns = [32, 64, 128, 256]
    res = {}
    for name in ("smooth", "oscillatory"):
        res[name] = {"wavelet": [], "tt": [], "qtt": []}
        for N in Ns:
            vol = make(name, N); lv = int(round(np.log2(N))) - 2
            res[name]["wavelet"].append(wavelet_params(vol, lv))
            res[name]["tt"].append(tt3_compress(vol, TGT*0.7)[1])
            res[name]["qtt"].append(qtt_compress3(vol, TGT*0.7)[1])
            print("  %-11s N=%3d  wavelet=%8d  TT=%8d  QTT=%8d" %
                  (name, N, res[name]["wavelet"][-1], res[name]["tt"][-1],
                   res[name]["qtt"][-1]))
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.5))
    for ax, name in zip(axs, res):
        for key, lab, col in [("wavelet", "wavelet", "#d62728"),
                              ("tt", "3-core TT", "#9467bd"),
                              ("qtt", "QTT", "#2ca02c")]:
            ax.loglog(Ns, res[name][key], "o-", color=col, label=lab)
        ax.loglog(Ns, [res[name]["tt"][0]*(n/Ns[0]) for n in Ns], "k--",
                  lw=.8, alpha=.5, label="O(N) ref")
        ax.set_title(name); ax.set_xlabel("N (grid per axis)")
        ax.set_ylabel("params @ err~1e-3"); ax.grid(True, which="both", alpha=.3)
        ax.legend(fontsize=8)
    fig.suptitle("3-D compression cost vs resolution (fixed ~1e-3 error)")
    fig.tight_layout(); fig.savefig("wtt3d_scaling.png", dpi=120)
    print("  saved wtt3d_scaling.png")

if __name__ == "__main__":
    main()
