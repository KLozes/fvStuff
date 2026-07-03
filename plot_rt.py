#!/usr/bin/env python3
"""
plot_rt.py — visualise RT0/P0 DG simulation output

Usage:
  python3 plot_rt.py rt_vortex 5        # plot frame 5 of vortex run
  python3 plot_rt.py rt_sod    10       # plot frame 10 of Sod run
  python3 plot_rt.py rt_vortex 0 5      # plot frames 0..5 as a strip
  python3 plot_rt.py rt_vortex          # plot all available frames

Binary file format: 4-byte int N, then N*N float32 rho (row-major, j=0 bottom),
then N*N float32 p.
"""

import sys
import os
import glob
import numpy as np
import matplotlib
matplotlib.use("Agg")          # headless — saves PNG
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

# ── parse args ────────────────────────────────────────────────────────────────
prefix  = sys.argv[1] if len(sys.argv) > 1 else "rt_vortex"
f_start = int(sys.argv[2]) if len(sys.argv) > 2 else 0
f_end   = int(sys.argv[3]) if len(sys.argv) > 3 else f_start

# auto-detect if no explicit frame given
if len(sys.argv) <= 2:
    bins = sorted(glob.glob(f"{prefix}_*.bin"))
    if not bins:
        print(f"No files matching {prefix}_*.bin found.")
        sys.exit(1)
    f_start = 0
    f_end   = len(bins) - 1

frames = list(range(f_start, f_end + 1))

# ── load helper ───────────────────────────────────────────────────────────────
def load_frame(prefix, idx):
    fn = f"{prefix}_{idx:04d}.bin"
    if not os.path.exists(fn):
        raise FileNotFoundError(fn)
    with open(fn, "rb") as f:
        N = np.frombuffer(f.read(4), dtype=np.int32)[0]
        rho = np.frombuffer(f.read(N*N*4), dtype=np.float32).reshape(N, N)
        p   = np.frombuffer(f.read(N*N*4), dtype=np.float32).reshape(N, N)
    return N, rho, p

# ── single-frame detailed plot ────────────────────────────────────────────────
def plot_single(prefix, idx, save=True):
    N, rho, p = load_frame(prefix, idx)

    # derived: velocity magnitude if we only have rho and p, skip; just do rho & p
    # perturbation fields (background rho_inf=1, p_inf=1/gamma for vortex)
    is_vortex = "vortex" in prefix
    gamma = 1.4
    rho_bg = 1.0
    p_bg   = 1.0 / gamma if is_vortex else None

    fig, axes = plt.subplots(1, 2, figsize=(12, 5.5))
    fig.suptitle(f"{prefix}  frame {idx:04d}   N={N}×{N}", fontsize=13)

    # ── density ──────────────────────────────────────────────────────────────
    ax = axes[0]
    if is_vortex:
        drho = rho - rho_bg
        vmax = max(abs(drho.min()), abs(drho.max()), 1e-12)
        im = ax.imshow(drho, origin="lower", extent=[0,1,0,1],
                       cmap="RdBu_r", vmin=-vmax, vmax=vmax)
        ax.set_title(r"$\rho - \rho_\infty$")
    else:
        im = ax.imshow(rho, origin="lower", extent=[0,1,0,1],
                       cmap="plasma", vmin=rho.min(), vmax=rho.max())
        ax.set_title(r"$\rho$")
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    ax.set_xlabel("x"); ax.set_ylabel("y")

    # ── pressure ─────────────────────────────────────────────────────────────
    ax = axes[1]
    if is_vortex:
        dp = p - p_bg
        vmax = max(abs(dp.min()), abs(dp.max()), 1e-12)
        im = ax.imshow(dp, origin="lower", extent=[0,1,0,1],
                       cmap="RdBu_r", vmin=-vmax, vmax=vmax)
        ax.set_title(r"$p - p_\infty$")
    else:
        im = ax.imshow(p, origin="lower", extent=[0,1,0,1],
                       cmap="viridis", vmin=p.min(), vmax=p.max())
        ax.set_title(r"$p$")
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    ax.set_xlabel("x"); ax.set_ylabel("y")

    plt.tight_layout()
    outfn = f"{prefix}_{idx:04d}.png"
    if save:
        plt.savefig(outfn, dpi=150, bbox_inches="tight")
        print(f"  Saved {outfn}")
    return fig

# ── strip of multiple frames ──────────────────────────────────────────────────
def plot_strip(prefix, frames):
    nf = len(frames)
    if nf == 1:
        plot_single(prefix, frames[0])
        return

    fig, axes = plt.subplots(2, nf, figsize=(4*nf, 7))
    is_vortex = "vortex" in prefix
    gamma = 1.4
    rho_bg = 1.0
    p_bg   = 1.0 / gamma if is_vortex else None

    rhos, ps, Ns = [], [], []
    for idx in frames:
        try:
            N, rho, p = load_frame(prefix, idx)
            rhos.append(rho); ps.append(p); Ns.append(N)
        except FileNotFoundError:
            rhos.append(None); ps.append(None); Ns.append(0)

    # global symmetric range for perturbations
    if is_vortex:
        r_max = max((abs(r - rho_bg).max() for r in rhos if r is not None), default=1e-12)
        p_max = max((abs(q - p_bg).max()   for q in ps   if q is not None), default=1e-12)
    else:
        r_min = min((r.min() for r in rhos if r is not None), default=0)
        r_max = max((r.max() for r in rhos if r is not None), default=1)
        p_min = min((q.min() for q in ps   if q is not None), default=0)
        p_max = max((q.max() for q in ps   if q is not None), default=1)

    for col, idx in enumerate(frames):
        rho = rhos[col]; p = ps[col]; N = Ns[col]
        if rho is None:
            axes[0][col].axis("off"); axes[1][col].axis("off"); continue

        ext = [0, 1, 0, 1]

        ax = axes[0][col]
        if is_vortex:
            im = ax.imshow(rho - rho_bg, origin="lower", extent=ext,
                           cmap="RdBu_r", vmin=-r_max, vmax=r_max)
        else:
            im = ax.imshow(rho, origin="lower", extent=ext,
                           cmap="plasma", vmin=r_min, vmax=r_max)
        ax.set_title(f"t-frame {idx}", fontsize=9)
        if col == 0:
            ax.set_ylabel(r"$\rho - \rho_\infty$" if is_vortex else r"$\rho$", fontsize=9)
        ax.set_xticks([]); ax.set_yticks([])
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.02)

        ax = axes[1][col]
        if is_vortex:
            im = ax.imshow(p - p_bg, origin="lower", extent=ext,
                           cmap="RdBu_r", vmin=-p_max, vmax=p_max)
        else:
            im = ax.imshow(p, origin="lower", extent=ext,
                           cmap="viridis", vmin=p_min, vmax=p_max)
        if col == 0:
            ax.set_ylabel(r"$p - p_\infty$" if is_vortex else r"$p$", fontsize=9)
        ax.set_xticks([]); ax.set_yticks([])
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.02)

    N0 = next(n for n in Ns if n > 0)
    fig.suptitle(f"{prefix}   frames {frames[0]}–{frames[-1]}   N={N0}×{N0}", fontsize=12)
    plt.tight_layout()
    outfn = f"{prefix}_strip_{frames[0]:04d}_{frames[-1]:04d}.png"
    plt.savefig(outfn, dpi=150, bbox_inches="tight")
    print(f"  Saved {outfn}")

# ── main ──────────────────────────────────────────────────────────────────────
if len(frames) == 1:
    plot_single(prefix, frames[0])
else:
    plot_strip(prefix, frames)
