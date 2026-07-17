#!/usr/bin/env python3
"""Quantized tensor-train (QTT) decomposition of 2-D arrays.

A 2^L x 2^L array is reshaped to an order-2L tensor of mode size 2 by
interleaving the row/column bits coarse-to-fine (standard image QTT), then
TT-SVD (Oseledets) compresses it with a relative-Frobenius tolerance, keeping
singular values above the per-core threshold.  Parameter count = total size of
the TT cores.
"""
import numpy as np


# ── QTT reshape (bit-interleaved, coarse-to-fine) ──────────────────────────
def _perm(L):
    p = []
    for t in range(L):
        p += [t, L + t]                 # i_t (MSB-first) then j_t
    return p

def to_qtt(A):
    L = int(round(np.log2(A.shape[0])))
    assert A.shape == (2**L, 2**L)
    T = A.reshape((2,) * (2*L))         # axes: i_{L-1..0}, j_{L-1..0}
    return np.ascontiguousarray(np.transpose(T, _perm(L))), L

def from_qtt(T, L):
    T = T.reshape((2,) * (2*L))
    inv = np.argsort(_perm(L))
    return np.ascontiguousarray(np.transpose(T, inv)).reshape(2**L, 2**L)


# ── TT-SVD with truncation ─────────────────────────────────────────────────
def _trunc_rank(s, delta2):
    """Smallest r with sum of discarded squared singular values <= delta2."""
    tail = np.concatenate([np.cumsum((s**2)[::-1])[::-1], [0.0]])   # tail[r]=sum(s[r:]^2)
    r = int(np.argmax(tail <= delta2))
    return max(r, 1)

def tt_svd(T, eps, max_rank=None):
    """TT-SVD of tensor T (any shape) at relative tolerance eps.
       max_rank caps every core rank (bounds memory for high-rank data).
       Returns list of cores (each shape r_{k-1} x n_k x r_k)."""
    shape = T.shape
    d = len(shape)
    if d == 1:
        return [T.reshape(1, shape[0], 1)]
    nrm = np.linalg.norm(T)
    delta2 = (eps * nrm)**2 / (d - 1) if nrm > 0 else 0.0
    cores = []
    C = T.reshape(shape[0], -1)
    r_prev = 1
    for k in range(d - 1):
        n_k = shape[k]
        C = C.reshape(r_prev * n_k, -1)
        U, s, Vt = np.linalg.svd(C, full_matrices=False)
        r = _trunc_rank(s, delta2)
        if max_rank is not None: r = min(r, max_rank)
        U, s, Vt = U[:, :r], s[:r], Vt[:r, :]
        cores.append(U.reshape(r_prev, n_k, r))
        C = (s[:, None] * Vt)                       # (r x rest)
        r_prev = r
    cores.append(C.reshape(r_prev, shape[-1], 1))
    return cores

def tt_to_full(cores):
    out = cores[0]                                  # (1, n0, r0)
    for G in cores[1:]:
        out = np.tensordot(out, G, axes=([out.ndim - 1], [0]))
    return out.reshape(tuple(G.shape[1] for G in cores))

def tt_params(cores):
    return int(sum(G.size for G in cores))


# ── convenience: QTT-compress a 2-D array ──────────────────────────────────
def qtt_compress(A, eps):
    """Returns (reconstruction A_hat, parameter count)."""
    if A.size == 0 or np.max(np.abs(A)) == 0.0:
        return np.zeros_like(A), 0
    T, L = to_qtt(A)
    cores = tt_svd(T, eps)
    Ahat = from_qtt(tt_to_full(cores), L)
    return Ahat, tt_params(cores)


# ── 3-D QTT (bit-interleave i,j,k coarse-to-fine) + plain 3-core TT ────────
def _perm3(L):
    p = []
    for t in range(L):
        p += [t, L + t, 2*L + t]
    return p

def to_qtt3(A):
    L = int(round(np.log2(A.shape[0])))
    assert A.shape == (2**L, 2**L, 2**L)
    T = A.reshape((2,) * (3*L))
    return np.ascontiguousarray(np.transpose(T, _perm3(L))), L

def from_qtt3(T, L):
    T = T.reshape((2,) * (3*L))
    inv = np.argsort(_perm3(L))
    return np.ascontiguousarray(np.transpose(T, inv)).reshape(2**L, 2**L, 2**L)

def qtt_compress3(A, eps, max_rank=None):
    if A.size == 0 or np.max(np.abs(A)) == 0.0:
        return np.zeros_like(A), 0
    T, L = to_qtt3(A)
    cores = tt_svd(T, eps, max_rank)
    return from_qtt3(tt_to_full(cores), L), tt_params(cores)

def tt3_compress(A, eps, max_rank=None):
    """Plain 3-core tensor train of an N x N x N array (2 SVDs)."""
    if A.size == 0 or np.max(np.abs(A)) == 0.0:
        return np.zeros_like(A), 0
    cores = tt_svd(A, eps, max_rank)
    return tt_to_full(cores), tt_params(cores)


def _selftest():
    rng = np.random.default_rng(1)
    L = 6; N = 2**L
    # QTT reshape round-trip
    A = rng.standard_normal((N, N))
    T, LL = to_qtt(A)
    assert np.max(np.abs(from_qtt(T, LL) - A)) < 1e-12
    # full-rank TT-SVD reproduces exactly
    cores = tt_svd(T, 0.0)
    err_exact = np.max(np.abs(from_qtt(tt_to_full(cores), LL) - A))
    # low-rank structure: separable/smooth array compresses well
    x = np.linspace(0, 1, N)
    S = np.exp(-((x[:, None] - .5)**2 + (x[None, :] - .5)**2) / 0.05)
    Shat, p = qtt_compress(S, 1e-3)
    rel = np.linalg.norm(Shat - S) / np.linalg.norm(S)
    print(f"[tt] QTT reshape round-trip ok; full-rank TT-SVD err = {err_exact:.2e}")
    print(f"[tt] smooth {N}x{N} ({N*N} vals) -> QTT params={p}  rel.err={rel:.2e}")
    # 3D round-trips
    M = 16
    B = rng.standard_normal((M, M, M))
    T3, L3 = to_qtt3(B)
    assert np.max(np.abs(from_qtt3(T3, L3) - B)) < 1e-11
    assert np.max(np.abs(tt_to_full(tt_svd(T3, 0.0)).reshape(M,M,M)*0 +
                         from_qtt3(tt_to_full(tt_svd(T3, 0.0)), L3) - B)) < 1e-10
    zz = np.linspace(0, 1, M)
    G3 = np.exp(-((zz[:,None,None]-.5)**2 + (zz[None,:,None]-.5)**2 +
                  (zz[None,None,:]-.5)**2)/0.05)
    _, pq = qtt_compress3(G3, 1e-3); _, pt = tt3_compress(G3, 1e-3)
    print(f"[tt] 3D smooth {M}^3 ({M**3} vals) -> QTT params={pq}, 3-core TT params={pt}")


if __name__ == "__main__":
    _selftest()
