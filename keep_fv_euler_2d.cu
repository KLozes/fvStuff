/*
 * keep_fv_euler_2d.cu
 *
 * Collocated cell-centred finite-volume solver for the 2-D compressible
 * Euler equations, implementing the structure-preserving KEEP method of
 *
 *   K. Bahuguna, R. Kolluru, S.V. Raghurama Rao (2025),
 *   "Structure-preserving schemes conserving entropy and kinetic energy",
 *   arXiv:2505.13347   (KEEP.md in this repo).
 *
 * Discretisation
 * --------------
 *   State    : U = [ρ, ρu, ρv, E]  — ALL collocated P0 cell averages.
 *              No enriched / Raviart-Thomas momentum DOFs (plain FV).
 *   Recon.   : none — the ECKEP central flux is 2nd-order on cell averages
 *              (paper Sec 4.2); interface states are the cell averages.
 *   Flux     : ECKEP entropy-conserving + kinetic-energy-preserving central
 *              flux  F_EC  (Sec 3.2) PLUS MOVERS-RH entropy-stable scalar
 *              dissipation  F_RH  (Sec 5) — i.e. the ES scheme, eq. (41):
 *                     F = F_EC + F_RH .
 *   Time     : SSP-RK3.
 *   BC       : periodic.
 *
 * The flux is built in the rotated normal/tangential frame; the normal
 * velocity uₙ plays the role of "u" in the 1-D paper and the tangential
 * momentum ρuₜ is passively advected with a KEP convective form.
 *
 * Compile:
 *   nvcc -O3 -arch=native --expt-relaxed-constexpr -o keep_fv keep_fv_euler_2d.cu -lm
 *
 * Usage (same CLI as the rt_dg driver):
 *   ./keep_fv N            Circular Sod shock tube, [0,1]^2
 *   ./keep_fv N lmv [eps]  Low-Mach vortex (Barsukow et al. Sec 2), default eps=0.1
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <cuda_runtime.h>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

/* ── compile-time constants ─────────────────────────────────────────────── */
#define GAMMA_V  1.4
#define BLOCK1D  256
#define GS_NBLK  256
#define G_GHOST  2            /* ghost layers (two-point flux needs ±1)       */

#define KEEP_DELTA 1e-16      /* denominator floor in α₃        (paper)       */
#define KEEP_THETA 0.1        /* sonic-fix parameter Θ          (paper eq.40) */

/* HES hybrid (Sec 6):  F = F_EC + (1-φ)F_R + φ F_RH                          */
#define HES_Q    12.0         /* entropy-distance sensor scaling q ∈ [8,16]   */
#define HES_EPS  0.05         /* sensor clip threshold ε       ∈ [1e-2,1e-1]  */
#define HES_JST  (1.0/32.0)   /* JST 4th-order coefficient  α_R = λ̃/32       */

/* conserved-variable indices (collocated) */
#define NVAR  4
#define Q_RHO 0
#define Q_MX  1   /* ρu */
#define Q_MY  2   /* ρv */
#define Q_E   3

#define real double

/* ── CUDA error check ───────────────────────────────────────────────────── */
#define CK(x) do { \
    cudaError_t _e = (x); \
    if (_e != cudaSuccess) { \
        fprintf(stderr, "CUDA error at %s:%d: %s\n", \
                __FILE__, __LINE__, cudaGetErrorString(_e)); \
        exit(1); \
    } \
} while (0)

/* padded SoA index */
#define IDX_P(q,jp,ip,Np)  ((size_t)(q)*(Np)*(Np) + (size_t)(jp)*(Np) + (ip))

/* ═══════════════════════════════════════════════════════════════════════════
 * Device helpers
 * ═══════════════════════════════════════════════════════════════════════════ */
__device__ __forceinline__ real
sound_speed(real rho, real p)
{
    return sqrt((real)GAMMA_V * p / rho);
}

__device__ __forceinline__ real
pressure_cons(real rho, real mx, real my, real E)
{
    real u = mx / rho, v = my / rho;
    return ((real)GAMMA_V - 1.) * (E - 0.5 * rho * (u*u + v*v));
}

/* primitives W=[ρ,u,v,p] from conserved at padded cell (jp,ip) */
__device__ __forceinline__ void
prim_at(const real * __restrict__ Qp, int jp, int ip, int Np, real W[4])
{
    real rho = fmax(Qp[IDX_P(Q_RHO, jp, ip, Np)], 1e-14);
    real mx  = Qp[IDX_P(Q_MX,  jp, ip, Np)];
    real my  = Qp[IDX_P(Q_MY,  jp, ip, Np)];
    real E   = Qp[IDX_P(Q_E,   jp, ip, Np)];
    W[0] = rho;
    W[1] = mx / rho;
    W[2] = my / rho;
    W[3] = fmax(pressure_cons(rho, mx, my, E), 1e-14);
}

/* Cartesian conserved vector U=[ρ,ρu,ρv,E] at padded cell (jp,ip) */
__device__ __forceinline__ void
cons_at(const real * __restrict__ Qp, int jp, int ip, int Np, real U[4])
{
    U[0] = fmax(Qp[IDX_P(Q_RHO, jp, ip, Np)], 1e-14);
    U[1] = Qp[IDX_P(Q_MX, jp, ip, Np)];
    U[2] = Qp[IDX_P(Q_MY, jp, ip, Np)];
    U[3] = Qp[IDX_P(Q_E,  jp, ip, Np)];
}

/* Cartesian entropy variables V from primitives W=[ρ,u,v,p] (s=ln(p/ρ^γ)) */
__device__ __forceinline__ void
entropy_vars(const real W[4], real V[4])
{
    real rho=W[0], u=W[1], v=W[2], p=W[3];
    real s = log(p / pow(rho, (real)GAMMA_V));
    V[0] = (GAMMA_V - s)/(GAMMA_V-1.) - rho*(u*u + v*v)/(2.*p);
    V[1] = rho*u/p;
    V[2] = rho*v/p;
    V[3] = -rho/p;
}

/* ═══════════════════════════════════════════════════════════════════════════
 * KEEP central flux  F_EC  +  MOVERS-RH coefficient αS
 *
 *   Fec[4]  : ECKEP entropy-conserving + kinetic-energy-preserving central
 *             flux (Sec 3.2), returned in the Cartesian frame.
 *   *alphaS : scalar MOVERS-RH dissipation coefficient (Sec 5, eqs. 39–40).
 *             The RH flux  -½ αS (U_R-U_L)  is scalar×identity, hence frame
 *             independent, so the caller applies it to Cartesian jumps.
 *
 * WL,WR = primitive cell averages [ρ,u,v,p];  n=(nx,ny).
 * ═══════════════════════════════════════════════════════════════════════════ */
__device__ void
keep_flux(const real WL[4], const real WR[4],
          real nx, real ny, real Fec[4], real *alphaS)
{
    real rL=WL[0], uL=WL[1], vL=WL[2], pL=WL[3];
    real rR=WR[0], uR=WR[1], vR=WR[2], pR=WR[3];

    /* rotate velocities into normal (n) / tangential (t) frame */
    real unL = uL*nx + vL*ny,  utL = -uL*ny + vL*nx;
    real unR = uR*nx + vR*ny,  utR = -uR*ny + vR*nx;

    real EL = pL/(GAMMA_V-1.) + 0.5*rL*(unL*unL + utL*utL);
    real ER = pR/(GAMMA_V-1.) + 0.5*rR*(unR*unR + utR*utR);

    /* conserved & physical-flux states (normal frame) */
    real UL[4] = { rL, rL*unL, rL*utL, EL };
    real UR[4] = { rR, rR*unR, rR*utR, ER };
    real FnL[4] = { rL*unL, rL*unL*unL + pL, rL*unL*utL, (EL+pL)*unL };
    real FnR[4] = { rR*unR, rR*unR*unR + pR, rR*unR*utR, (ER+pR)*unR };

    /* entropy variables V (s = ln(p/ρ^γ)), normal frame */
    real sL = log(pL / pow(rL, (real)GAMMA_V));
    real sR = log(pR / pow(rR, (real)GAMMA_V));
    real VL[4] = { (GAMMA_V - sL)/(GAMMA_V-1.) - rL*(unL*unL+utL*utL)/(2.*pL),
                    rL*unL/pL, rL*utL/pL, -rL/pL };
    real VR[4] = { (GAMMA_V - sR)/(GAMMA_V-1.) - rR*(unR*unR+utR*utR)/(2.*pR),
                    rR*unR/pR, rR*utR/pR, -rR/pR };
    real dV[4] = { VR[0]-VL[0], VR[1]-VL[1], VR[2]-VL[2], VR[3]-VL[3] };

    /* entropy-flux potential ψ = ρuₙ  →  Δψ */
    real dpsi = rR*unR - rL*unL;

    /* ── ECKEP central flux  F_EC  (Sec 3.2), normal frame ──────────────── */
    real Fc0  = 0.5*(FnL[0] + FnR[0]);           /* mass = mean ρuₙ (KEP)    */
    real un_b = 0.5*(unL + unR);
    real ut_b = 0.5*(utL + utR);
    real p_b  = 0.5*(pL  + pR);
    real Fc1  = Fc0*un_b + p_b;                  /* normal mom: ṁ·ūₙ + p̄   */
    real Fc2  = Fc0*ut_b;                         /* tang.  mom: ṁ·ūₜ        */
    real Fe_b = 0.5*(FnL[3] + FnR[3]);
    real num  = dV[0]*Fc0 + dV[1]*Fc1 + dV[2]*Fc2 + dV[3]*Fe_b - dpsi;
    real Fc3  = Fe_b - (num / (dV[3]*dV[3] + KEEP_DELTA)) * dV[3];

    /* rotate F_EC back to Cartesian */
    Fec[0] = Fc0;
    Fec[1] = Fc1*nx - Fc2*ny;
    Fec[2] = Fc1*ny + Fc2*nx;
    Fec[3] = Fc3;

    /* ── MOVERS-RH dissipation coefficient αS (Sec 5, eqs. 39–40) ───────── */
    real r_av = 0.5*(rL + rR);
    real a_av = sound_speed(fmax(r_av,1e-14), fmax(p_b,1e-14));
    real l1 = fabs(un_b - a_av), l2 = fabs(un_b), l3 = fabs(un_b + a_av);
    real lam_min = fmin(l1, fmin(l2, l3));
    real lam_max = fmax(l1, fmax(l2, l3));

    real aS = lam_max;       /* αS = min over fields of sₖ = |ΔFₖ/ΔUₖ|     */
    for (int q = 0; q < 4; q++) {
        real dU = UR[q] - UL[q];
        real sk;
        if (fabs(dU) < 1e-12) {
            sk = lam_min;                     /* ΔU→0 : collapse to λmin     */
        } else {
            sk = fabs((FnR[q] - FnL[q]) / dU);
            if      (sk >= lam_max) sk = lam_max;
            else if (sk <= lam_min) sk = lam_min;
        }
        if (sk < aS) aS = sk;
    }
    /* sonic-point fix (eq. 40): lift small αS to ≥ Θ/2 so the RH branch is
     * never degenerate at stagnation faces.                                */
    if (aS < KEEP_THETA)
        aS = (aS*aS + KEEP_THETA*KEEP_THETA) / (2.*KEEP_THETA);
    *alphaS = aS;
}

/* ═══════════════════════════════════════════════════════════════════════════
 * Periodic boundary conditions (fill ghost layers)
 * ═══════════════════════════════════════════════════════════════════════════ */
__global__ void
apply_bc(real * __restrict__ Qp, int N, int Np)
{
    const int g   = G_GHOST;
    const int Np2 = Np * Np;
    for (int k = blockIdx.x * blockDim.x + threadIdx.x; k < Np2;
             k += blockDim.x * gridDim.x) {
        int jp = k / Np, ip = k % Np;
        bool gx = (ip < g || ip >= N + g);
        bool gy = (jp < g || jp >= N + g);
        if (!gx && !gy) continue;
        int ir = ip, jr = jp;
        if (ip < g)      ir = ip + N;
        if (ip >= N + g) ir = ip - N;
        if (jp < g)      jr = jp + N;
        if (jp >= N + g) jr = jp - N;
        for (int q = 0; q < NVAR; q++)
            Qp[IDX_P(q, jp, ip, Np)] = Qp[IDX_P(q, jr, ir, Np)];
    }
}

/* ═══════════════════════════════════════════════════════════════════════════
 * Entropy distance  ED = (U_R-U_C)·(V_R-V_C)  per cell (max over its R,T faces)
 * Used to normalise the shock sensor (Sec 6).  d_ed[k] = max(ED_right, ED_top).
 * ═══════════════════════════════════════════════════════════════════════════ */
__device__ __forceinline__ real
edist(const real Wc[4], const real Uc[4], const real Wn[4], const real Un[4])
{
    real Vc[4], Vn[4];
    entropy_vars(Wc, Vc);
    entropy_vars(Wn, Vn);
    real ed = 0.;
    for (int q = 0; q < 4; q++) ed += (Un[q]-Uc[q]) * (Vn[q]-Vc[q]);
    return fabs(ed);
}

__global__ void
compute_ed(const real * __restrict__ Qp, real * __restrict__ d_ed, int N, int Np)
{
    const int g  = G_GHOST;
    const int N2 = N * N;
    for (int k = blockIdx.x * blockDim.x + threadIdx.x; k < N2;
             k += blockDim.x * gridDim.x) {
        int j = k / N, i = k % N;
        int jp = j + g, ip = i + g;
        real Wc[4],Uc[4],Wxp[4],Uxp[4],Wyp[4],Uyp[4];
        prim_at(Qp, jp, ip,   Np, Wc);  cons_at(Qp, jp, ip,   Np, Uc);
        prim_at(Qp, jp, ip+1, Np, Wxp); cons_at(Qp, jp, ip+1, Np, Uxp);
        prim_at(Qp, jp+1, ip, Np, Wyp); cons_at(Qp, jp+1, ip, Np, Uyp);
        real er = edist(Wc, Uc, Wxp, Uxp);
        real et = edist(Wc, Uc, Wyp, Uyp);
        d_ed[k] = fmax(er, et);
    }
}

/* shock sensor φ at a face (eq. 46/48): 1 at shocks, 0 in smooth regions */
__device__ __forceinline__ real
shock_phi(real ed_face, real ed_max)
{
    real sed = ed_face / (ed_max + 1e-30);     /* scaled entropy distance     */
    real phi = 1. - exp(-HES_Q * sed);
    return (phi > HES_EPS) ? 1. : 0.;
}

/* ═══════════════════════════════════════════════════════════════════════════
 * RHS — collocated finite volume, HES hybrid (Sec 6, eq. 52):
 *     F = F_EC + (1-φ) F_R + φ F_RH
 *   F_EC : ECKEP entropy-conserving + KEP central flux (2nd-order, no recon).
 *   F_R  : JST 4th-order background dissipation, α_R = λ̃/32 (eqs. 49–51) —
 *          damps acoustic/odd-even modes in smooth regions (low-Mach stable).
 *   F_RH : MOVERS-RH scalar dissipation, switched on by the entropy-distance
 *          shock sensor φ at shocks only.
 * Each cell evaluates all four faces (computed twice, no atomics).
 * ═══════════════════════════════════════════════════════════════════════════ */
__global__ void
compute_rhs(const real * __restrict__ Qp,
            real       * __restrict__ RHS,
            real       * __restrict__ lam_out,
            int N, int Np, real h, real ed_max)
{
    const int g  = G_GHOST;
    const int N2 = N * N;

    for (int k = blockIdx.x * blockDim.x + threadIdx.x; k < N2;
             k += blockDim.x * gridDim.x) {
        int j = k / N, i = k % N;
        int jp = j + g, ip = i + g;

        /* primitives on the 5-point cross stencil (centre ±1 ±2) */
        real Wc[4],Wxm[4],Wxmm[4],Wxp[4],Wxpp[4],Wym[4],Wymm[4],Wyp[4],Wypp[4];
        prim_at(Qp, jp,   ip,   Np, Wc);
        prim_at(Qp, jp,   ip-1, Np, Wxm);   prim_at(Qp, jp,   ip-2, Np, Wxmm);
        prim_at(Qp, jp,   ip+1, Np, Wxp);   prim_at(Qp, jp,   ip+2, Np, Wxpp);
        prim_at(Qp, jp-1, ip,   Np, Wym);   prim_at(Qp, jp-2, ip,   Np, Wymm);
        prim_at(Qp, jp+1, ip,   Np, Wyp);   prim_at(Qp, jp+2, ip,   Np, Wypp);

        /* conserved on the same stencil */
        real Uc[4],Uxm[4],Uxmm[4],Uxp[4],Uxpp[4],Uym[4],Uymm[4],Uyp[4],Uypp[4];
        cons_at(Qp, jp,   ip,   Np, Uc);
        cons_at(Qp, jp,   ip-1, Np, Uxm);   cons_at(Qp, jp,   ip-2, Np, Uxmm);
        cons_at(Qp, jp,   ip+1, Np, Uxp);   cons_at(Qp, jp,   ip+2, Np, Uxpp);
        cons_at(Qp, jp-1, ip,   Np, Uym);   cons_at(Qp, jp-2, ip,   Np, Uymm);
        cons_at(Qp, jp+1, ip,   Np, Uyp);   cons_at(Qp, jp+2, ip,   Np, Uypp);

        real F[4][4];   /* face fluxes: 0=Right 1=Left 2=Top 3=Bottom         */

#define FACE(slot, WLs, WRs, ULs, URs, U_p2, U_m1, nxs, nys, edf) do {        \
        real fec[4], aS;                                                      \
        keep_flux((WLs), (WRs), (nxs), (nys), fec, &aS);                      \
        /* JST 4th-order:  α_R (U_{+2} - 3U_{+1} + 3U_C - U_{-1}) */          \
        real unL = (WLs)[1]*(nxs) + (WLs)[2]*(nys);                           \
        real unR = (WRs)[1]*(nxs) + (WRs)[2]*(nys);                           \
        real aR  = sound_speed(0.5*((WLs)[0]+(WRs)[0]),                       \
                               0.5*((WLs)[3]+(WRs)[3]));                      \
        real lamt = fmax(fabs(unL), fabs(unR)) + aR;       /* RICCA λ̃        */\
        real aJST = HES_JST * lamt;                                          \
        real phi  = shock_phi((edf), ed_max);                                \
        for (int q = 0; q < 4; q++) {                                         \
            real jst = aJST * ((U_p2)[q] - 3.*(URs)[q] + 3.*(ULs)[q] - (U_m1)[q]); \
            real rh  = -0.5 * aS * ((URs)[q] - (ULs)[q]);                     \
            F[slot][q] = fec[q] + (1.-phi)*jst + phi*rh;                      \
        } } while (0)

        /* RIGHT face (i|i+1): stencil U_{i-1},U_i,U_{i+1},U_{i+2}            */
        FACE(0, Wc,  Wxp, Uc,  Uxp, Uxpp, Uxm, 1., 0.,
             edist(Wc, Uc, Wxp, Uxp));
        /* LEFT face (i-1|i): stencil U_{i-2},U_{i-1},U_i,U_{i+1}            */
        FACE(1, Wxm, Wc,  Uxm, Uc,  Uxp,  Uxmm, 1., 0.,
             edist(Wxm, Uxm, Wc, Uc));
        /* TOP face (j|j+1)                                                   */
        FACE(2, Wc,  Wyp, Uc,  Uyp, Uypp, Uym, 0., 1.,
             edist(Wc, Uc, Wyp, Uyp));
        /* BOTTOM face (j-1|j)                                                */
        FACE(3, Wym, Wc,  Uym, Uc,  Uyp,  Uymm, 0., 1.,
             edist(Wym, Uym, Wc, Uc));
#undef FACE

        real inv_h = 1. / h;
        for (int q = 0; q < 4; q++)
            RHS[(size_t)q*N2 + k] = -(F[0][q] - F[1][q] + F[2][q] - F[3][q]) * inv_h;

        /* CFL spectral radius */
        real cs = sound_speed(Wc[0], Wc[3]);
        lam_out[k] = (fabs(Wc[1]) + cs) + (fabs(Wc[2]) + cs);
    }
}

/* ═══════════════════════════════════════════════════════════════════════════
 * SSP-RK3 stage kernels (operate on conserved DOFs)
 * ═══════════════════════════════════════════════════════════════════════════ */
__global__ void
rk_stage(real * __restrict__ Uout, const real * __restrict__ Ua,
         const real * __restrict__ Ub, const real * __restrict__ L,
         real ca, real cb, real cl, real dt, int N, int Np)
{
    const int g   = G_GHOST;
    const int N2  = N * N;
    const int Np2 = Np * Np;
    for (int k = blockIdx.x * blockDim.x + threadIdx.x; k < N2;
             k += blockDim.x * gridDim.x) {
        int j = k / N, i = k % N;
        int idxp = (j + g) * Np + (i + g);
        for (int q = 0; q < NVAR; q++)
            Uout[(size_t)q*Np2 + idxp] = ca * Ua[(size_t)q*Np2 + idxp]
                                       + cb * Ub[(size_t)q*Np2 + idxp]
                                       + cl * dt * L[(size_t)q*N2 + k];
        Uout[(size_t)Q_RHO*Np2 + idxp] = fmax(Uout[(size_t)Q_RHO*Np2 + idxp], 1e-14);
        Uout[(size_t)Q_E  *Np2 + idxp] = fmax(Uout[(size_t)Q_E  *Np2 + idxp], 1e-14);
    }
}

/* ═══════════════════════════════════════════════════════════════════════════
 * Reductions
 * ═══════════════════════════════════════════════════════════════════════════ */
__global__ void
reduce_max(const real * __restrict__ in, real * __restrict__ out, int n)
{
    extern __shared__ real sm[];
    int tid = threadIdx.x;
    real v = 0.;
    for (int k = blockIdx.x * blockDim.x + tid; k < n; k += blockDim.x * gridDim.x)
        v = fmax(v, in[k]);
    sm[tid] = v; __syncthreads();
    for (int s = BLOCK1D / 2; s > 0; s >>= 1) {
        if (tid < s) sm[tid] = fmax(sm[tid], sm[tid + s]);
        __syncthreads();
    }
    if (tid == 0) out[blockIdx.x] = sm[0];
}
__global__ void
reduce_sum(const real * __restrict__ in, real * __restrict__ out, int n)
{
    extern __shared__ real sm[];
    int tid = threadIdx.x;
    real v = 0.;
    for (int k = blockIdx.x * blockDim.x + tid; k < n; k += blockDim.x * gridDim.x)
        v += in[k];
    sm[tid] = v; __syncthreads();
    for (int s = BLOCK1D / 2; s > 0; s >>= 1) {
        if (tid < s) sm[tid] += sm[tid + s];
        __syncthreads();
    }
    if (tid == 0) out[blockIdx.x] = sm[0];
}
static real gpu_max(const real *d_in, real *d_tmp, int n) {
    reduce_max<<<GS_NBLK, BLOCK1D, BLOCK1D*sizeof(real)>>>(d_in, d_tmp, n);
    reduce_max<<<1, BLOCK1D, BLOCK1D*sizeof(real)>>>(d_tmp, d_tmp, GS_NBLK);
    real v; CK(cudaMemcpy(&v, d_tmp, sizeof(real), cudaMemcpyDeviceToHost));
    return v;
}
static real gpu_sum(const real *d_in, real *d_tmp, int n) {
    reduce_sum<<<GS_NBLK, BLOCK1D, BLOCK1D*sizeof(real)>>>(d_in, d_tmp, n);
    reduce_sum<<<1, BLOCK1D, BLOCK1D*sizeof(real)>>>(d_tmp, d_tmp, GS_NBLK);
    real v; CK(cudaMemcpy(&v, d_tmp, sizeof(real), cudaMemcpyDeviceToHost));
    return v;
}

/* ═══════════════════════════════════════════════════════════════════════════
 * Initial conditions
 * ═══════════════════════════════════════════════════════════════════════════ */

/* Circular Sod shock tube on [0,1]^2: high (ρ,p) disc inside r<0.25.
 * The interface is smoothed with a tanh profile of width ~2 cells so the
 * jump is resolved (no sub-cell discontinuity for the unlimited scheme).   */
__global__ void
ic_sod(real * __restrict__ Qp, int N, int Np, real h)
{
    const int g  = G_GHOST;
    const int N2 = N * N;
    for (int k = blockIdx.x * blockDim.x + threadIdx.x; k < N2;
             k += blockDim.x * gridDim.x) {
        int j = k / N, i = k % N;
        int jp = j + g, ip = i + g;
        real xc = (i + 0.5) * h, yc = (j + 0.5) * h;
        real r  = sqrt((xc-0.5)*(xc-0.5) + (yc-0.5)*(yc-0.5));
        /* phi → 1 inside the disc, 0 outside, smoothed over ~2h */
        real delta = 1.0 * h;
        real phi = 0.5 * (1. + tanh((0.25 - r) / delta));
        real rho = 0.125 + (1.0 - 0.125) * phi;   /* 1.0 inside, 0.125 outside */
        real p   = 0.1   + (1.0 - 0.1  ) * phi;   /* 1.0 inside, 0.1   outside */
        Qp[IDX_P(Q_RHO, jp, ip, Np)] = rho;
        Qp[IDX_P(Q_MX,  jp, ip, Np)] = 0.;
        Qp[IDX_P(Q_MY,  jp, ip, Np)] = 0.;
        Qp[IDX_P(Q_E,   jp, ip, Np)] = p / (GAMMA_V - 1.);
    }
}

/* Low-Mach vortex (Barsukow et al. 2025, Sec 2): stationary; p0 sets Ma=eps */
__device__ real lmv_vphi(real r)
{
    if (r < 0.2) return 5. * r;
    if (r < 0.4) return 2. - 5. * r;
    return 0.;
}
__device__ real lmv_pressure(real r, real p0)
{
    if (r < 0.2) return p0 + 12.5 * r * r;
    if (r < 0.4) return p0 + 4.*log(5.*r) + 4. - 20.*r + 12.5*r*r;
    return p0 + 4.*log(2.) - 2.;
}
__device__ void
lmv_exact(real x, real y, real p0, real *rho, real *u, real *v, real *p)
{
    real dx = x - 0.5, dy = y - 0.5;
    real r  = sqrt(dx*dx + dy*dy);
    real vp = lmv_vphi(r);
    real ir = (r > 1e-30) ? 1./r : 0.;
    *rho = 1.;
    *u   = -vp * dy * ir;
    *v   =  vp * dx * ir;
    *p   = lmv_pressure(r, p0);
}
__global__ void
ic_lmv(real * __restrict__ Qp, int N, int Np, real h, real eps)
{
    const real p0 = 1. / (GAMMA_V * eps * eps) - 0.5;
    const int g  = G_GHOST;
    const int N2 = N * N;
    for (int k = blockIdx.x * blockDim.x + threadIdx.x; k < N2;
             k += blockDim.x * gridDim.x) {
        int j = k / N, i = k % N;
        int jp = j + g, ip = i + g;
        real xc = (i + 0.5) * h, yc = (j + 0.5) * h;
        real rho, u, v, p; lmv_exact(xc, yc, p0, &rho, &u, &v, &p);
        Qp[IDX_P(Q_RHO, jp, ip, Np)] = rho;
        Qp[IDX_P(Q_MX,  jp, ip, Np)] = rho * u;
        Qp[IDX_P(Q_MY,  jp, ip, Np)] = rho * v;
        Qp[IDX_P(Q_E,   jp, ip, Np)] = p/(GAMMA_V-1.) + 0.5*rho*(u*u+v*v);
    }
}

/* L2 error vs. exact stationary low-Mach vortex (exact = IC) */
__global__ void
err_lmv(const real * __restrict__ Qp, real * __restrict__ erho,
        real * __restrict__ ep, int N, int Np, real h, real eps)
{
    const real p0 = 1. / (GAMMA_V * eps * eps) - 0.5;
    const int g  = G_GHOST;
    const int N2 = N * N;
    for (int k = blockIdx.x * blockDim.x + threadIdx.x; k < N2;
             k += blockDim.x * gridDim.x) {
        int j = k / N, i = k % N;
        int jp = j + g, ip = i + g;
        real xc = (i + 0.5) * h, yc = (j + 0.5) * h;
        real re, ue, ve, pe; lmv_exact(xc, yc, p0, &re, &ue, &ve, &pe);
        real rho = fmax(Qp[IDX_P(Q_RHO, jp, ip, Np)], 1e-14);
        real p   = pressure_cons(rho, Qp[IDX_P(Q_MX,jp,ip,Np)],
                                 Qp[IDX_P(Q_MY,jp,ip,Np)], Qp[IDX_P(Q_E,jp,ip,Np)]);
        real dr = rho - re, dp = p - pe;
        erho[k] = dr*dr*h*h;
        ep  [k] = dp*dp*h*h;
    }
}

/* ═══════════════════════════════════════════════════════════════════════════
 * Output: extract two host fields and write a 2-panel jet PNG
 * ═══════════════════════════════════════════════════════════════════════════ */
__global__ void
extract_fields(const real * __restrict__ Qp, real * __restrict__ left,
               real * __restrict__ right, int N, int Np, int lmv)
{
    const int g  = G_GHOST;
    const int N2 = N * N;
    for (int k = blockIdx.x * blockDim.x + threadIdx.x; k < N2;
             k += blockDim.x * gridDim.x) {
        int j = k / N, i = k % N;
        int jp = j + g, ip = i + g;
        real rho = fmax(Qp[IDX_P(Q_RHO, jp, ip, Np)], 1e-14);
        real mx  = Qp[IDX_P(Q_MX, jp, ip, Np)];
        real my  = Qp[IDX_P(Q_MY, jp, ip, Np)];
        real E   = Qp[IDX_P(Q_E,  jp, ip, Np)];
        real p   = pressure_cons(rho, mx, my, E);
        left [k] = lmv ? sqrt(mx*mx + my*my)/rho : rho;   /* |vel| or ρ      */
        right[k] = p;
    }
}

static void
write_png_2panel(const char *fn, const real *L, real Lmin, real Lmax,
                 const real *R, real Rmin, real Rmax, int N)
{
    const int W = 2*N, H = N;
    unsigned char *px = (unsigned char*)malloc((size_t)W*H*3);
    if (!px) { fprintf(stderr, "OOM png\n"); return; }
    for (int r = 0; r < H; r++) {
        int pj = N-1-r;
        for (int c = 0; c < W; c++) {
            const real *pan = (c < N) ? L : R;
            real vmin = (c < N) ? Lmin : Rmin, vmax = (c < N) ? Lmax : Rmax;
            int ci = (c < N) ? c : c - N;
            real t = (pan[pj*N+ci] - vmin) / (vmax - vmin + 1e-30);
            if (!(t > 0.)) t = 0.;  if (t > 1.) t = 1.;
            /* jet colormap */
            real cr = fmax(0., fmin(1., 1.5 - fabs(4.*t - 3.)));
            real cg = fmax(0., fmin(1., 1.5 - fabs(4.*t - 2.)));
            real cb = fmax(0., fmin(1., 1.5 - fabs(4.*t - 1.)));
            unsigned char *pix = px + (r*W + c)*3;
            pix[0] = (unsigned char)(cr*255.);
            pix[1] = (unsigned char)(cg*255.);
            pix[2] = (unsigned char)(cb*255.);
        }
    }
    stbi_write_png(fn, W, H, 3, px, W*3);
    free(px);
    printf("  Saved %s\n", fn);
}

/* ═══════════════════════════════════════════════════════════════════════════
 * main
 * ═══════════════════════════════════════════════════════════════════════════ */
int main(int argc, char **argv)
{
    int N = (argc > 1) ? atoi(argv[1]) : 200;
    const char *mode = (argc > 2) ? argv[2] : "";
    int  do_lmv = (strcmp(mode, "lmv") == 0);
    real lmv_eps = 0.1;
    if (do_lmv && argc > 3) {
        char *endp; real m = strtod(argv[3], &endp);
        if (endp != argv[3] && *endp == '\0' && m > 0.) lmv_eps = m;
    }

    const real L_domain = 1.0;
    const real CFL = 0.2;
    real t_end = do_lmv ? 2.*(real)M_PI*0.2 : 0.25;  /* lmv: one vortex period; sod: 0.25 */

    real h = L_domain / N;
    const int g = G_GHOST;
    int Np = N + 2*g;

    int dev; cudaDeviceProp prop;
    CK(cudaGetDevice(&dev));
    CK(cudaGetDeviceProperties(&prop, dev));
    printf("========================================================\n");
    printf("  Device : %s\n", prop.name);
    printf("  Collocated FV + KEEP (ECKEP+ES, eq.41)  --  %s\n",
           do_lmv ? "Low-Mach Vortex" : "Circular Sod");
    if (do_lmv) printf("  eps=%.4f  (p0=%.3f)\n", lmv_eps, 1./(GAMMA_V*lmv_eps*lmv_eps)-0.5);
    printf("  N=%dx%d  Np=%d  h=%.5f  CFL=%.2f  t_end=%.4f\n",
           N, N, Np, h, CFL, t_end);
    printf("  DOFs: rho, rho*u, rho*v, E  (collocated P0, no enriched DOFs)\n");
    printf("========================================================\n");

    size_t sz1  = (size_t)N*N*sizeof(real);
    size_t szc  = (size_t)NVAR*sz1;
    size_t szp  = (size_t)NVAR*Np*Np*sizeof(real);

    real *d_U,*d_U0,*d_U1,*d_U2,*d_RHS,*d_lam,*d_ed,*d_tmp,*d_left,*d_right,*d_erho,*d_ep;
    CK(cudaMalloc(&d_U,  szp));  CK(cudaMalloc(&d_U0, szp));
    CK(cudaMalloc(&d_U1, szp));  CK(cudaMalloc(&d_U2, szp));
    CK(cudaMalloc(&d_RHS, szc)); CK(cudaMalloc(&d_lam, sz1));
    CK(cudaMalloc(&d_ed, sz1));
    CK(cudaMalloc(&d_left, sz1)); CK(cudaMalloc(&d_right, sz1));
    CK(cudaMalloc(&d_erho, sz1)); CK(cudaMalloc(&d_ep, sz1));
    CK(cudaMalloc(&d_tmp, GS_NBLK*sizeof(real)));

    /* compute RHS(U) into d_RHS with the entropy-distance shock sensor       */
#define RHS_OF(U) do { \
        compute_ed<<<GS_NBLK,BLOCK1D>>>((U), d_ed, N, Np); \
        real _edm = gpu_max(d_ed, d_tmp, N*N); \
        CK(cudaMemset(d_RHS, 0, szc)); \
        compute_rhs<<<GS_NBLK,BLOCK1D>>>((U), d_RHS, d_lam, N, Np, h, _edm); \
    } while (0)
    real *h_left = (real*)malloc(sz1), *h_right = (real*)malloc(sz1);

    const char *prefix = do_lmv ? "keep_lmv" : "keep_sod";

    /* IC */
    if (do_lmv) ic_lmv<<<GS_NBLK,BLOCK1D>>>(d_U, N, Np, h, lmv_eps);
    else        ic_sod<<<GS_NBLK,BLOCK1D>>>(d_U, N, Np, h);
    apply_bc<<<GS_NBLK,BLOCK1D>>>(d_U, N, Np);
    CK(cudaDeviceSynchronize());

    /* color bounds from IC */
    extract_fields<<<GS_NBLK,BLOCK1D>>>(d_U, d_left, d_right, N, Np, do_lmv);
    CK(cudaDeviceSynchronize());
    real Lmin=0., Lmax=gpu_max(d_left, d_tmp, N*N);
    real Rmin=1e30, Rmax=gpu_max(d_right, d_tmp, N*N);
    { /* min via negate trick */
        reduce_max<<<GS_NBLK,BLOCK1D,BLOCK1D*sizeof(real)>>>(d_right, d_tmp, N*N);
        /* simple host min */
        CK(cudaMemcpy(h_right, d_right, sz1, cudaMemcpyDeviceToHost));
        for (int q=0;q<N*N;q++) Rmin = fmin(Rmin, h_right[q]);
    }
    if (Lmax <= 0.) Lmax = 1.;
    real rpad = 0.05*(Rmax-Rmin); if (rpad<1e-9) rpad=0.1;
    Rmin -= rpad; Rmax += rpad;

#define WRITE_FRAME(idx) do { \
    extract_fields<<<GS_NBLK,BLOCK1D>>>(d_U, d_left, d_right, N, Np, do_lmv); \
    CK(cudaDeviceSynchronize()); \
    real _Lhi = gpu_max(d_left, d_tmp, N*N); \
    CK(cudaMemcpy(h_left,  d_left,  sz1, cudaMemcpyDeviceToHost)); \
    CK(cudaMemcpy(h_right, d_right, sz1, cudaMemcpyDeviceToHost)); \
    real _Rlo=1e30,_Rhi=-1e30; \
    for (int _q=0;_q<N*N;_q++){_Rlo=fmin(_Rlo,h_right[_q]);_Rhi=fmax(_Rhi,h_right[_q]);} \
    char _fn[80]; sprintf(_fn,"figures/%s_%04d.png", prefix,(idx)); \
    write_png_2panel(_fn, h_left, 0., fmax(_Lhi,1e-12), h_right, _Rlo, _Rhi, N); \
    if (do_lmv) { \
        err_lmv<<<GS_NBLK,BLOCK1D>>>(d_U, d_erho, d_ep, N, Np, h, lmv_eps); \
        CK(cudaDeviceSynchronize()); \
        real _sr=gpu_sum(d_erho,d_tmp,N*N), _sp=gpu_sum(d_ep,d_tmp,N*N); \
        printf("    L2_rho=%.4e  L2_p=%.4e\n", sqrt(_sr), sqrt(_sp)); \
    } \
} while(0)

    WRITE_FRAME(0);

    RHS_OF(d_U);
    CK(cudaDeviceSynchronize());
    real lam_max = gpu_max(d_lam, d_tmp, N*N);

    const int N_FRAMES = 10;
    int frame = 1; real t_next = t_end/N_FRAMES;
    int step = 0; real t = 0.;
    struct timespec ts0, ts1; clock_gettime(CLOCK_MONOTONIC, &ts0);

    while (t < t_end) {
        if (!(lam_max > 0.)) { fprintf(stderr,"NaN lam at step %d\n",step); break; }
        real dt = CFL * h / lam_max;
        real ttar = (t_next < t_end) ? t_next : t_end;
        if (t + dt > ttar) dt = ttar - t;
        if (dt < 1e-14) { fprintf(stderr,"dt underflow\n"); break; }

        CK(cudaMemcpy(d_U0, d_U, szp, cudaMemcpyDeviceToDevice));

        /* stage 1:  U1 = U0 + dt L(U0) */
        rk_stage<<<GS_NBLK,BLOCK1D>>>(d_U1, d_U0, d_U0, d_RHS, 1.,0.,1., dt, N, Np);
        apply_bc<<<GS_NBLK,BLOCK1D>>>(d_U1, N, Np);

        /* stage 2:  U2 = 3/4 U0 + 1/4 U1 + 1/4 dt L(U1) */
        RHS_OF(d_U1);
        rk_stage<<<GS_NBLK,BLOCK1D>>>(d_U2, d_U0, d_U1, d_RHS, 0.75,0.25,0.25, dt, N, Np);
        apply_bc<<<GS_NBLK,BLOCK1D>>>(d_U2, N, Np);

        /* stage 3:  U = 1/3 U0 + 2/3 U2 + 2/3 dt L(U2) */
        RHS_OF(d_U2);
        rk_stage<<<GS_NBLK,BLOCK1D>>>(d_U, d_U0, d_U2, d_RHS, 1./3.,2./3.,2./3., dt, N, Np);
        apply_bc<<<GS_NBLK,BLOCK1D>>>(d_U, N, Np);

        RHS_OF(d_U);
        lam_max = gpu_max(d_lam, d_tmp, N*N);
        t += dt; step++;

        if (t >= t_next - 1e-12 && frame <= N_FRAMES) {
            clock_gettime(CLOCK_MONOTONIC, &ts1);
            real el = (ts1.tv_sec-ts0.tv_sec)+(ts1.tv_nsec-ts0.tv_nsec)*1e-9;
            printf("  frame %2d/%d  step %5d  t=%.5f  lam=%.3e  elapsed=%.1fs\n",
                   frame, N_FRAMES, step, t, lam_max, el);
            fflush(stdout);
            WRITE_FRAME(frame);
            frame++; t_next = frame * t_end / N_FRAMES;
        }
    }

    if (do_lmv) {
        err_lmv<<<GS_NBLK,BLOCK1D>>>(d_U, d_erho, d_ep, N, Np, h, lmv_eps);
        CK(cudaDeviceSynchronize());
        real sr=gpu_sum(d_erho,d_tmp,N*N), sp=gpu_sum(d_ep,d_tmp,N*N);
        printf("  ── FINAL L2 ERRORS (Low-Mach Vortex) ────────────────\n");
        printf("  eps=%.4f  N=%d  h=%.5f  t=%.4f\n", lmv_eps, N, h, t);
        printf("  L2(rho) = %.6e\n", sqrt(sr));
        printf("  L2(p)   = %.6e\n", sqrt(sp));
        printf("  ─────────────────────────────────────────────────────\n");
    }

    CK(cudaDeviceSynchronize());
    clock_gettime(CLOCK_MONOTONIC, &ts1);
    real wall = (ts1.tv_sec-ts0.tv_sec)+(ts1.tv_nsec-ts0.tv_nsec)*1e-9;
    printf("  Done: %d steps  t=%.6f  wall=%.2fs\n", step, t, wall);
    printf("  Output: figures/%s_0000.png .. figures/%s_%04d.png\n",
           prefix, prefix, N_FRAMES);
    printf("  Left panel: %s   Right panel: pressure  (jet)\n",
           do_lmv ? "|velocity|" : "density");

    free(h_left); free(h_right);
    cudaFree(d_U); cudaFree(d_U0); cudaFree(d_U1); cudaFree(d_U2);
    cudaFree(d_RHS); cudaFree(d_lam); cudaFree(d_ed); cudaFree(d_tmp);
    cudaFree(d_left); cudaFree(d_right); cudaFree(d_erho); cudaFree(d_ep);
    (void)Lmin; (void)Lmax; (void)Rmin; (void)Rmax;
    return 0;
}
