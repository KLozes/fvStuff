/*
 * p1_dg_euler_2d.cu
 * Mixed DG/FV — FULL-LINEAR (P1) momentum + P0 density/energy:
 *   density ρ   : P0 (piecewise constant, cell average)
 *   momentum m  : P1 vector — EACH component carries the full linear basis {1,ξ,η}
 *                   ρu = mxa + mxs·ξ + mxsy·η     (3 DOF)
 *                   ρv = mya + mysx·ξ + mys·η     (3 DOF)
 *                 Full linear space, NOT bilinear (no ξη term), and NOT the
 *                 div-conforming RT0 subspace (which drops mxsy, mysx).
 *   energy E    : P0 (piecewise constant, cell average)
 *   Riemann     : HLLC (default; override with -Driemann_n=roe_n|hll_n)
 *   Time        : SSP-RK3
 *   Quadrature  : 1-point (face midpoint) for all boundary integrals
 *   Volume term : analytic (P1 quadratic forms, incl. cross ρuv term)
 *   Precision   : DOUBLE  (real = reall = double) — needed to see whether the
 *                 low-Mach stationary state survives without float noise.
 *
 * Derived from rt_dg_euler_2d.cu.  Experiment: does the full P1 momentum space
 * still preserve low-Mach stationary states (Gresho / low-Mach vortex,
 * ./p1_dg N lmv [eps]) the way the RT0 space does?
 *
 * Compile:
 *   nvcc -O3 -arch=native --expt-relaxed-constexpr -o p1_dg p1_dg_euler_2d.cu -lm
 *
 * Test cases (Barsukow, Ciallella, Ricchiuto, Torlo 2025 — arXiv:2506.21700):
 *
 *   ./rt_dg N              Circular Sod shock tube (Sec 6.2.4), [0,1]^2
 *   ./rt_dg N Ma           Original scaled isentropic vortex, [0,1]^2
 *   ./rt_dg N pv           Paper isentropic vortex, STATIONARY (Sec 6.2.1)
 *                            domain [0,10]^2, h=10/N, tf=1
 *                            expect superconvergence: ~2nd order in L2
 *   ./rt_dg N pvmove       Paper isentropic vortex, MOVING  u0=v0=1
 *                            tf=10, expect 1st order (non-stationary)
 *   ./rt_dg N kh           Kelvin-Helmholtz, Ma=0.01 (Sec 6.2.5)
 *                            domain [0,2]^2 (y shifted by 1), tf=80
 *   ./rt_dg N dsl [Ma]    Doubly periodic shear layer (Bell-Colella-Glaz 1989)
 *                            domain [0,1]^2, default Ma=0.1
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
/* Circular-Sod IC strength (inside the disk); outside fixed at rho=0.125,p=0.1.
 * Defaults reproduce the standard case (ratio ~100); override for stronger blasts:
 *   nvcc ... -DSOD_RHOIN=1000 -DSOD_PIN=1000   → pressure ratio 1e4, etc. */
#ifndef SOD_RHOIN
#define SOD_RHOIN 10.0
#endif
#ifndef SOD_PIN
#define SOD_PIN 10.0
#endif
#define BLOCK1D  256
#define GS_NBLK  256

/* DOF indices in the padded state array (NVAR = 8, full-P1 momentum) */
#define NVAR  16  /* FULL BILINEAR Q1 for ALL fields: ρ,ρu,ρv,E each {1,ξ,η,ξη} */
#define Q_RHO     0   /* density   cell average       (Q1 mode 1)    */
#define Q_MXA     1   /* x-momentum cell average      (Q1 mode 1)    */
#define Q_MXS     2   /* x-momentum ξ-slope  (∂x)     (Q1 mode ξ)    */
#define Q_MXSY    3   /* x-momentum η-slope  (∂y)     (Q1 mode η)    */
#define Q_MYA     4   /* y-momentum cell average      (Q1 mode 1)    */
#define Q_MYSX    5   /* y-momentum ξ-slope  (∂x)     (Q1 mode ξ)    */
#define Q_MYS     6   /* y-momentum η-slope  (∂y)     (Q1 mode η)    */
#define Q_E       7   /* total energy cell average    (Q1 mode 1)    */
#define Q_RHO_SX  8   /* density   ξ-slope  (∂x)      (Q1 mode ξ)    */
#define Q_RHO_SY  9   /* density   η-slope  (∂y)      (Q1 mode η)    */
#define Q_E_SX   10   /* energy    ξ-slope  (∂x)      (Q1 mode ξ)    */
#define Q_E_SY   11   /* energy    η-slope  (∂y)      (Q1 mode η)    */
#define Q_RHO_SXY 12  /* density   ξη cross          (Q1 mode ξη)   [bilinear] */
#define Q_MXSXY  13   /* x-momentum ξη cross         (Q1 mode ξη)   [bilinear] */
#define Q_MYSXY  14   /* y-momentum ξη cross         (Q1 mode ξη)   [bilinear] */
#define Q_E_SXY  15   /* energy    ξη cross          (Q1 mode ξη)   [bilinear] */

/*
 * Modal ↔ nodal P1 relations  (ξ = 2(x-xc)/h ∈ [-1,1], η = 2(y-yc)/h ∈ [-1,1]):
 *   (ρu)(ξ,η) = mxa + mxs·ξ + mxsy·η
 *   (ρv)(ξ,η) = mya + mysx·ξ + mys·η
 *
 * Face-midpoint traces (the perpendicular slope vanishes at a face midpoint):
 *   right (ξ=+1,η=0): ρu = mxa+mxs, ρv = mya+mysx
 *   left  (ξ=-1,η=0): ρu = mxa-mxs, ρv = mya-mysx
 *   top   (ξ=0,η=+1): ρu = mxa+mxsy, ρv = mya+mys
 *   bot   (ξ=0,η=-1): ρu = mxa-mxsy, ρv = mya-mys
 *
 *   Projection from the four face-midpoint normal fluxes (mr,ml,mt,mb) and the
 *   two off-diagonal midpoint fluxes still needs the tangential midpoint values;
 *   the IC kernels build all six modes directly from face-midpoint momenta.
 */

/* g = 2 ghost-cell layers (needed for 3rd-order linear upwind stencil) */
#define G_GHOST 2

/* ── CUDA error check ───────────────────────────────────────────────────── */
#define CK(x) do { \
    cudaError_t _e = (x); \
    if (_e != cudaSuccess) { \
        fprintf(stderr, "CUDA error at %s:%d: %s\n", \
                __FILE__, __LINE__, cudaGetErrorString(_e)); \
        exit(1); \
    } \
} while (0)

/* ── Array indexing: SoA layout Q[q][jp][ip], padded size Np×Np ─────────── */
#define IDX_P(q,jp,ip,Np)  ((size_t)(q)*(Np)*(Np) + (size_t)(jp)*(Np) + (ip))

/* ═══════════════════════════════════════════════════════════════════════════
 * Device helpers
 * ═══════════════════════════════════════════════════════════════════════════ */

/* precision: override with  nvcc -Dreal=double -Dreall=double  */
#ifndef real
#define real float
#endif
#ifndef reall
#define reall float
#endif

/* Computed in reall (double): the (E − ½ρu²) subtraction cancels at low Mach. */
__device__ __forceinline__ real
sound_speed(real rho, real p)
{
    reall cs2 = (reall)GAMMA_V * (reall)p / (reall)rho;
    return (real)sqrt(cs2);
}

/* Pressure from P0 ρ, E and modal RT0 cell averages mxa, mya.
 * Internals in reall: (E − ½ρ|u|²) catastrophically cancels at low Mach,
 * so float here loses all pressure information. Result demoted to real
 * for storage — the cancellation has already been resolved.              */
__device__ __forceinline__ real
pressure_from_rt0(real rho, real mxa, real mya, real E)
{
    reall u_avg = (reall)mxa / (reall)rho;
    reall v_avg = (reall)mya / (reall)rho;
    reall p = ((reall)GAMMA_V - 1.) * ((reall)E - 0.5 * (reall)rho * (u_avg*u_avg + v_avg*v_avg));
    return (real)(p);
}

/* ── Positivity floors (shared by the limiter and the pointwise helpers) ──── */
#define PFLOOR    ((reall)1e-9)    /* pressure floor kept by the limiter      */
#define RHOFLOOR  ((real )1e-12)   /* density floor                           */

/* Conserved point state (ρ,ρu,ρv,E) → primitives W=[ρ,u,v,p], with floors.
 * Used for both face Riemann traces and volume flux quadrature now that ALL
 * fields are P1 (ρ,E vary linearly inside the cell). reall for low-Mach. */
__device__ __forceinline__ void
cons_pt_to_W(real rho, real mx, real my, real E, real W[4])
{
    real  r = fmax(rho, RHOFLOOR);
    reall u = (reall)mx / (reall)r;
    reall v = (reall)my / (reall)r;
    reall p = ((reall)GAMMA_V - 1.0) * ((reall)E - 0.5 * (reall)r * (u*u + v*v));
    W[0] = r;
    W[1] = (real)u;
    W[2] = (real)v;
    W[3] = fmax((real)p, (real)PFLOOR);
}

/* Conserved point state → x-flux F[4] and y-flux G[4] (Euler), floored p. */
__device__ __forceinline__ void
euler_FG_pt(real rho, real mx, real my, real E, real F[4], real G[4])
{
    real W[4]; cons_pt_to_W(rho, mx, my, E, W);
    real u = W[1], v = W[2], pp = W[3];
    F[0] = mx;             G[0] = my;
    F[1] = mx*u + pp;      G[1] = my*u;        /* ρu²+p ; ρvu   */
    F[2] = mx*v;           G[2] = my*v + pp;   /* ρuv   ; ρv²+p */
    F[3] = (E + pp)*u;     G[3] = (E + pp)*v;
}

/* ── Standard HLLC flux with normal n=(nx,ny) ───────────────────────────── */
/* WL,WR: primitives [rho, u, v, p] (u,v are full Cartesian velocities)     */
/* F[4] = [mass, x-mom, y-mom, energy] fluxes in direction n                */
__device__ void
hllc_n(const real WL[4], const real WR[4],
           real nx, real ny, real F[4])
{
    real rL=WL[0], uL=WL[1], vL=WL[2], pL=WL[3];
    real rR=WR[0], uR=WR[1], vR=WR[2], pR=WR[3];

    real unL = uL*nx + vL*ny;
    real unR = uR*nx + vR*ny;
    real utL =-uL*ny + vL*nx;
    real utR =-uR*ny + vR*nx;

    real cL = sound_speed(rL, pL), cR = sound_speed(rR, pR);
    real EL = pL/(GAMMA_V-1.) + 0.5*rL*(uL*uL + vL*vL);
    real ER = pR/(GAMMA_V-1.) + 0.5*rR*(uR*uR + vR*vR);

    /* Einfeldt wave speed estimates */
    real SL = fmin(unL - cL, unR - cR);
    real SR = fmax(unL + cL, unR + cR);

    /* Normal-frame fluxes */
    real FnL[4] = { rL*unL,
                     rL*unL*unL + pL,
                     rL*unL*utL,
                     (EL + pL)*unL };
    real FnR[4] = { rR*unR,
                     rR*unR*unR + pR,
                     rR*unR*utR,
                     (ER + pR)*unR };

    if (SL >= 0.) {
        F[0] = FnL[0];
        F[1] = FnL[1]*nx - FnL[2]*ny;
        F[2] = FnL[1]*ny + FnL[2]*nx;
        F[3] = FnL[3];
        return;
    }
    if (SR <= 0.) {
        F[0] = FnR[0];
        F[1] = FnR[1]*nx - FnR[2]*ny;
        F[2] = FnR[1]*ny + FnR[2]*nx;
        F[3] = FnR[3];
        return;
    }

    /* Contact wave speed (Batten et al. 1997) */
    real dL = rL*(SL - unL), dR = rR*(SR - unR);
    real Ss = (pR - pL + dL*unL - dR*unR) / (dL - dR);

    /* Intermediate star states */
    real factL  = rL*(SL - unL)/(SL - Ss);
    real UsL[4] = { factL,
                     factL * Ss,
                     factL * utL,
                     factL * (EL/rL + (Ss-unL)*(Ss + pL/(rL*(SL-unL)))) };

    real factR  = rR*(SR - unR)/(SR - Ss);
    real UsR[4] = { factR,
                     factR * Ss,
                     factR * utR,
                     factR * (ER/rR + (Ss-unR)*(Ss + pR/(rR*(SR-unR)))) };

    /* Conservative states in normal/tangential frame */
    real UL[4] = { rL, rL*unL, rL*utL, EL };
    real UR[4] = { rR, rR*unR, rR*utR, ER };

    /* Standard HLLC: F*L if Ss>=0, F*R if Ss<0 */
    real Fn[4];
    if (Ss >= 0.)
        for (int q = 0; q < 4; q++)
            Fn[q] = FnL[q] + SL*(UsL[q] - UL[q]);
    else
        for (int q = 0; q < 4; q++)
            Fn[q] = FnR[q] + SR*(UsR[q] - UR[q]);

    /* Rotate back from normal/tangential to Cartesian */
    F[0] = Fn[0];
    F[1] = Fn[1]*nx - Fn[2]*ny;
    F[2] = Fn[1]*ny + Fn[2]*nx;
    F[3] = Fn[3];
}

/* ── Roe flux with normal n=(nx,ny) ─────────────────────────────────────── */
/* WL,WR: primitives [rho, u, v, p] (u,v are full Cartesian velocities)     */
/* F[4] = [mass, x-mom, y-mom, energy] fluxes in direction n                */
/*                                                                           */
/* Standard Roe average with Harten-Hyman entropy fix:                      */
/*   δ = max(0, 2(λR-λL))  per wave; |λ| replaced by max(|λ|, δ/2)         */
__device__ void
roe_n(const real WL[4], const real WR[4],
      real nx, real ny, real F[4])
{
    real rL=WL[0], uL=WL[1], vL=WL[2], pL=WL[3];
    real rR=WR[0], uR=WR[1], vR=WR[2], pR=WR[3];

    real cL = sound_speed(rL, pL), cR = sound_speed(rR, pR);

    /* Roe averages */
    real sqrL = sqrt(rL), sqrR = sqrt(rR);
    real denom = sqrL + sqrR;

    real uRoe  = (sqrL*uL  + sqrR*uR)  / denom;
    real vRoe  = (sqrL*vL  + sqrR*vR)  / denom;
    real HL    = (pL/(GAMMA_V-1.) + 0.5*rL*(uL*uL+vL*vL) + pL) / rL;
    real HR    = (pR/(GAMMA_V-1.) + 0.5*rR*(uR*uR+vR*vR) + pR) / rR;
    real HRoe  = (sqrL*HL + sqrR*HR) / denom;
    real c2Roe = (GAMMA_V-1.) * (HRoe - 0.5*(uRoe*uRoe + vRoe*vRoe));
    real cRoe  = sqrt(c2Roe > 1e-14 ? c2Roe : 1e-14);

    /* Rotate to normal/tangential frame */
    real unL = uL*nx + vL*ny,  utL = -uL*ny + vL*nx;
    real unR = uR*nx + vR*ny,  utR = -uR*ny + vR*nx;
    real unRoe = uRoe*nx + vRoe*ny;
    real utRoe = -uRoe*ny + vRoe*nx;

    /* Physical fluxes (normal-frame) */
    real EL = pL/(GAMMA_V-1.) + 0.5*rL*(uL*uL+vL*vL);
    real ER = pR/(GAMMA_V-1.) + 0.5*rR*(uR*uR+vR*vR);

    real FnL[4] = { rL*unL,        rL*unL*unL + pL, rL*unL*utL, (EL+pL)*unL };
    real FnR[4] = { rR*unR,        rR*unR*unR + pR, rR*unR*utR, (ER+pR)*unR };

    /* Jump in conserved variables (normal frame) */
    real drho = rR - rL;
    real dun  = unR - unL;
    real dut  = utR - utL;
    real dp   = pR  - pL;

    /* Roe-average density: sqrt(rL)*sqrt(rR) */
    real rRoe = sqrL * sqrR;

    /* Wave strengths (eigenvector decomposition of the jump) */
    real b1 = (dp - rRoe*cRoe*dun) / (2.*c2Roe);  /* left-acoustic   */
    real b2 = drho - dp/c2Roe;                       /* entropy         */
    real b3 = rRoe * dut;                            /* transverse shear */
    real b4 = (dp + rRoe*cRoe*dun) / (2.*c2Roe);  /* right-acoustic   */

    /* Eigenvalues */
    real lam1 = unRoe - cRoe;
    real lam2 = unRoe;
    real lam3 = unRoe;
    real lam4 = unRoe + cRoe;

    /* Harten-Hyman entropy fix */
    real lam1L = unL - cL, lam1R = unR - cR;
    real lam4L = unL + cL, lam4R = unR + cR;
    real eps1 = fmax(0., 2.*(lam1R - lam1L));
    real eps4 = fmax(0., 2.*(lam4R - lam4L));
    real al1 = fabs(lam1); if (al1 < eps1*0.5) al1 = (lam1*lam1 + eps1*eps1*0.25) / eps1;
    real al2 = fabs(lam2);
    real al3 = fabs(lam3);
    real al4 = fabs(lam4); if (al4 < eps4*0.5) al4 = (lam4*lam4 + eps4*eps4*0.25) / eps4;

    /* Roe dissipation: sum_k al_k * b_k * r_k
     * Right eigenvectors in normal/tangential frame:
     *   r1 = [1, unRoe-cRoe, utRoe, HRoe-unRoe*cRoe]
     *   r2 = [1, unRoe,      utRoe, ½(unRoe²+utRoe²)]
     *   r3 = [0, 0,          1,     utRoe            ]
     *   r4 = [1, unRoe+cRoe, utRoe, HRoe+unRoe*cRoe  ]     */
    real d_rho = al1*b1*1.           + al2*b2*1.      + al3*b3*0.      + al4*b4*1.;
    real d_un  = al1*b1*(unRoe-cRoe)  + al2*b2*unRoe    + al3*b3*0.      + al4*b4*(unRoe+cRoe);
    real d_ut  = al1*b1*utRoe         + al2*b2*utRoe     + al3*b3*1.      + al4*b4*utRoe;
    real d_E   = al1*b1*(HRoe-unRoe*cRoe)
                + al2*b2*0.5*(unRoe*unRoe+utRoe*utRoe)
                + al3*b3*utRoe
                + al4*b4*(HRoe+unRoe*cRoe);

    real Fn[4];
    Fn[0] = 0.5*(FnL[0]+FnR[0]) - 0.5*d_rho;
    Fn[1] = 0.5*(FnL[1]+FnR[1]) - 0.5*d_un;
    Fn[2] = 0.5*(FnL[2]+FnR[2]) - 0.5*d_ut;
    Fn[3] = 0.5*(FnL[3]+FnR[3]) - 0.5*d_E;

    /* Rotate back to Cartesian */
    F[0] = Fn[0];
    F[1] = Fn[1]*nx - Fn[2]*ny;
    F[2] = Fn[1]*ny + Fn[2]*nx;
    F[3] = Fn[3];
}

/* ── HLL flux with normal n=(nx,ny) ─────────────────────────────────────── */
/* Two-wave HLL: F_HLL = (SR*FL - SL*FR + SL*SR*(UR-UL)) / (SR-SL)         */
__device__ void
hll_n(const real WL[4], const real WR[4],
      real nx, real ny, real F[4])
{
    real rL=WL[0], uL=WL[1], vL=WL[2], pL=WL[3];
    real rR=WR[0], uR=WR[1], vR=WR[2], pR=WR[3];

    real unL = uL*nx + vL*ny,  utL = -uL*ny + vL*nx;
    real unR = uR*nx + vR*ny,  utR = -uR*ny + vR*nx;

    real cL = sound_speed(rL, pL), cR = sound_speed(rR, pR);
    real EL = pL/(GAMMA_V-1.) + 0.5*rL*(uL*uL+vL*vL);
    real ER = pR/(GAMMA_V-1.) + 0.5*rR*(uR*uR+vR*vR);

    /* Einfeldt wave speed estimates */
    real SL = fmin(unL - cL, unR - cR);
    real SR = fmax(unL + cL, unR + cR);

    /* Physical fluxes in normal/tangential frame */
    real FnL[4] = { rL*unL, rL*unL*unL+pL, rL*unL*utL, (EL+pL)*unL };
    real FnR[4] = { rR*unR, rR*unR*unR+pR, rR*unR*utR, (ER+pR)*unR };

    if (SL >= 0.) {
        F[0] = FnL[0];
        F[1] = FnL[1]*nx - FnL[2]*ny;
        F[2] = FnL[1]*ny + FnL[2]*nx;
        F[3] = FnL[3];
        return;
    }
    if (SR <= 0.) {
        F[0] = FnR[0];
        F[1] = FnR[1]*nx - FnR[2]*ny;
        F[2] = FnR[1]*ny + FnR[2]*nx;
        F[3] = FnR[3];
        return;
    }

    /* Conservative states in normal/tangential frame */
    real UL[4] = { rL, rL*unL, rL*utL, EL };
    real UR[4] = { rR, rR*unR, rR*utR, ER };

    real denom = 1. / (SR - SL);
    real Fn[4];
    for (int q = 0; q < 4; q++)
        Fn[q] = (SR*FnL[q] - SL*FnR[q] + SL*SR*(UR[q] - UL[q])) * denom;

    F[0] = Fn[0];
    F[1] = Fn[1]*nx - Fn[2]*ny;
    F[2] = Fn[1]*ny + Fn[2]*nx;
    F[3] = Fn[3];
}

/* ── Select Riemann solver below (override: nvcc -Driemann_n=roe_n …) ────── */
/* hllc_n   : standard HLLC                                                   */
/* hll_n    : standard HLL                                                    */
/* roe_n    : standard Roe                                                    */
#ifndef riemann_n
#  define riemann_n hllc_n
#endif

/* ═══════════════════════════════════════════════════════════════════════════
 * Apply boundary conditions — periodic wrap
 * Ghost layer (g=1): left ghost ← right interior edge, etc.
 * No sign flips: DOF values copied verbatim (the ring wraps).
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
        if (!gx && !gy) continue;  /* interior — skip */

        /* Periodic wrap: map ghost index to the opposite interior edge */
        int ir = ip, jr = jp;
        if (ip < g)      ir = ip + N;   /* left  ghost  ← right interior */
        if (ip >= N + g) ir = ip - N;   /* right ghost  ← left  interior */
        if (jp < g)      jr = jp + N;   /* bottom ghost ← top   interior */
        if (jp >= N + g) jr = jp - N;   /* top    ghost ← bottom interior */

        for (int q = 0; q < NVAR; q++)
            Qp[IDX_P(q, jp, ip, Np)] = Qp[IDX_P(q, jr, ir, Np)];
    }
}

#ifdef NEDELEC
/* ═══════════════════════════════════════════════════════════════════════════
 * Nédélec projection — restrict the full P1 momentum to the lowest-order
 * (first-kind) Nédélec / curl-conforming subspace by zeroing the two
 * DIVERGENCE slopes.  Leaves:
 *     ρu = mxa + mxsy·η   (const in ξ, linear in η)
 *     ρv = mya + mysx·ξ   (linear in ξ, const in η)
 * i.e. the rotated RT0 element: div-free in-cell, carries vorticity
 *     ω = (2/h)(mysx − mxsy),  tangential component continuous.
 * Applied to the whole padded array after every BC fill so mxs≡mys≡0 holds
 * in interior AND ghost cells for the next RHS evaluation. ─────────────── */
__global__ void
project_nedelec(real * __restrict__ Qp, int Np)
{
    const int Np2 = Np * Np;
    for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < Np2;
             idx += blockDim.x * gridDim.x) {
        Qp[Q_MXS * Np2 + idx] = 0.;
        Qp[Q_MYS * Np2 + idx] = 0.;
    }
}
#define NEDELEC_PROJECT(Qarr) project_nedelec<<<GS_NBLK, BLOCK1D>>>((Qarr), Np)
#else
#define NEDELEC_PROJECT(Qarr) ((void)0)
#endif

/* ═══════════════════════════════════════════════════════════════════════════
 * RHS kernel — FULL-P1 momentum / P0 DG
 *
 * For each cell (i,j) [physical], padded index (ip,jp) = (i+g, j+g):
 *
 *   Momentum DOFs (modal): mxa,mxs,mxsy  (ρu = mxa+mxs·ξ+mxsy·η)
 *                          mya,mysx,mys  (ρv = mya+mysx·ξ+mys·η)
 *   Density/energy: P0 averages ρ, E.
 *
 *   FACE TRACES come from each cell's OWN P1 polynomial evaluated at the face
 *   midpoint (a genuine DG-P1 trace):
 *     Right face (ξ=+1,η=0): ρu = mxa+mxs (normal),  ρv = mya+mysx (tangential)
 *     Top   face (ξ=0,η=+1): ρv = mya+mys (normal),  ρu = mxa+mxsy (tangential)
 *   The perpendicular slope vanishes at the midpoint, so the NORMAL momentum
 *   trace is identical to RT0 — only the TANGENTIAL momentum now carries an
 *   in-cell slope (mysx on x-faces, mxsy on y-faces) instead of a reconstructed
 *   neighbour slope.  Density and pressure keep the 3rd-order neighbour
 *   reconstruction (they are P0).
 *
 *   HLLC at each face: F_R, G_T  [mass, x-mom, y-mom, energy], scattered to both
 *   adjacent cells so each face flux is built once.
 *
 *   Weak form → mass-matrix inverse per mode:
 *     d(avg)  = (1/h²)·RHS_1 ,  d(ξ-slope)=d(η-slope)=(3/h²)·RHS_ξ,η
 *   with RHS_a = ∫_K (F·∇φ_a) − ∮_∂K φ_a F_n.
 *
 *   Volume terms (analytic, from ∫∫ over [-1,1]²):
 *     vol_x  = (1/h)∫ (ρu²+p)  = h[(mxa²+(mxs²+mxsy²)/3)/ρ + p]        → mxs
 *     vol_y  = (1/h)∫ (ρv²+p)  = h[(mya²+(mysx²+mys²)/3)/ρ + p]        → mys
 *     vol_xy = (1/h)∫ (ρuv)    = (h/ρ)(mxa·mya+(mxs·mysx+mxsy·mys)/3)  → mxsy,mysx
 *   (the cross ρuv term feeds the off-diagonal slopes; pressure only enters the
 *    diagonal vol_x/vol_y — the correct momentum-flux tensor structure.)
 *
 *   Face scatter (midpoint quadrature; a perpendicular face contributes 0 to a
 *   slope whose basis is odd along that face):
 *     x-face flux F_R → mxa(±), mxs(−3), mysx(−3, tangential-slope);   mxsy,mys 0
 *     y-face flux G_T → mya(±), mys(−3), mxsy(−3, tangential-slope);   mxs,mysx 0
 * ═══════════════════════════════════════════════════════════════════════════ */
__global__ void
compute_rhs(const real * __restrict__ Qp,
            real       * __restrict__ RHS,
            real       * __restrict__ lam_out,
            int N, int Np, real h)
{
    const int g   = G_GHOST;
    const int N2  = N * N;

    for (int k = blockIdx.x * blockDim.x + threadIdx.x; k < N2;
             k += blockDim.x * gridDim.x) {
        int j = k / N, i = k % N;
        int jp = j + g, ip = i + g;
        int kR = j * N + ((i + 1 < N) ? (i + 1) : 0);
        int kT = ((j + 1 < N) ? (j + 1) : 0) * N + i;

        /* ── Load center + right + top neighbour FULL P1 states (12 DOFs) ─
         * All fields are P1 now, so every face trace comes from the cell's OWN
         * polynomial (the Hermite reconstruction) sampled at the face Gauss
         * points — no neighbour ρ,p reconstruction.  Each cell builds only its
         * RIGHT and TOP faces and scatters the flux to both adjacent cells. */
        real Uc[NVAR], Ur[NVAR], Ut[NVAR];
#pragma unroll
        for (int q = 0; q < NVAR; q++) {
            Uc[q] = Qp[IDX_P(q, jp,   ip,   Np)];
            Ur[q] = Qp[IDX_P(q, jp,   ip+1, Np)];
            Ut[q] = Qp[IDX_P(q, jp+1, ip,   Np)];
        }

        const real GP = (real)0.577350269189625764;   /* 1/√3 : 2-pt Gauss node */
        real inv_h2 = 1. / (h * h);

        /* Conserved DOF index map: var v∈{mass,xmom,ymom,energy} → {avg,ξ,η,ξη}. */
        const int A  [4] = { Q_RHO,     Q_MXA,   Q_MYA,   Q_E     };
        const int SX [4] = { Q_RHO_SX,  Q_MXS,   Q_MYSX,  Q_E_SX  };
        const int SY [4] = { Q_RHO_SY,  Q_MXSY,  Q_MYS,   Q_E_SY  };
        const int SXY[4] = { Q_RHO_SXY, Q_MXSXY, Q_MYSXY, Q_E_SXY };

        /* ── Hermite face reconstruction (compact, 3rd order) ─────────────
         * Bilinear Q1: along the normal at a fixed tangential node t the field is
         *   f(n) = (avg + tang·t) + (norm-slope + cross·t)·n ,
         * so the NORMAL slope picks up the ξη cross term (cross·t).  The trace is
         * the compact quadratic:
         *   own (downwind)  @ +1 :  5/6·Vc + 2/3·s  + 1/6·Vn
         *   neigh (upwind)  @ −1 :  5/6·Vn − 2/3·sN + 1/6·Vc
         * with Vc,Vn the cell/neighbour VALUES at t and s,sN the (t-dependent)
         * normal slopes.  Tangential variation sampled at the 2 Gauss nodes. */
        const real C0 = (real)(5./6.), C1 = (real)(2./3.), C2 = (real)(1./6.);

        real I0R[4]={0.,0.,0.,0.}, IeR[4]={0.,0.,0.,0.};  /* right face: 0th & η-moment */
        real I0T[4]={0.,0.,0.,0.}, JxT[4]={0.,0.,0.,0.};  /* top   face: 0th & ξ-moment */

        /* ── RIGHT face (normal ξ, tangential η): reconstruct in ξ ──────── */
#pragma unroll
        for (int m = 0; m < 2; m++) {
            real eta = (m == 0) ? -GP : GP;
            real UL[4], UR[4];
#pragma unroll
            for (int v = 0; v < 4; v++) {
                real Vc = Uc[A[v]]  + Uc[SY[v]]*eta;    /* cell value at η_g       */
                real Vn = Ur[A[v]]  + Ur[SY[v]]*eta;    /* right-neigh value       */
                real sc = Uc[SX[v]] + Uc[SXY[v]]*eta;   /* cell  normal slope @ η_g */
                real sn = Ur[SX[v]] + Ur[SXY[v]]*eta;   /* neigh normal slope @ η_g */
                UL[v] = C0*Vc + C1*sc + C2*Vn;          /* downwind quad @ ξ=+1 */
                UR[v] = C0*Vn - C1*sn + C2*Vc;          /* upwind   quad @ ξ=−1 */
            }
            real WL[4], WR[4];
            cons_pt_to_W(UL[0],UL[1],UL[2],UL[3], WL);
            cons_pt_to_W(UR[0],UR[1],UR[2],UR[3], WR);
            real Fg[4]; riemann_n(WL, WR, 1., 0., Fg);
            for (int q = 0; q < 4; q++) { I0R[q] += Fg[q]; IeR[q] += eta*Fg[q]; }
        }
        for (int q = 0; q < 4; q++) { I0R[q] *= 0.5*h; IeR[q] *= 0.5*h; }

        /* ── TOP face (normal η, tangential ξ): reconstruct in η ─────────── */
#pragma unroll
        for (int m = 0; m < 2; m++) {
            real xi = (m == 0) ? -GP : GP;
            real UB[4], UT[4];
#pragma unroll
            for (int v = 0; v < 4; v++) {
                real Vc = Uc[A[v]]  + Uc[SX[v]]*xi;     /* cell value at ξ_g       */
                real Vn = Ut[A[v]]  + Ut[SX[v]]*xi;     /* top-neigh value         */
                real sc = Uc[SY[v]] + Uc[SXY[v]]*xi;    /* cell  normal slope @ ξ_g */
                real sn = Ut[SY[v]] + Ut[SXY[v]]*xi;    /* neigh normal slope @ ξ_g */
                UB[v] = C0*Vc + C1*sc + C2*Vn;          /* downwind quad @ η=+1 */
                UT[v] = C0*Vn - C1*sn + C2*Vc;          /* upwind   quad @ η=−1 */
            }
            real WB[4], WT[4];
            cons_pt_to_W(UB[0],UB[1],UB[2],UB[3], WB);
            cons_pt_to_W(UT[0],UT[1],UT[2],UT[3], WT);
            real Gg[4]; riemann_n(WB, WT, 0., 1., Gg);
            for (int q = 0; q < 4; q++) { I0T[q] += Gg[q]; JxT[q] += xi*Gg[q]; }
        }
        for (int q = 0; q < 4; q++) { I0T[q] *= 0.5*h; JxT[q] *= 0.5*h; }

        /* ── VOLUME: 2×2 Gauss of the Euler flux over [-1,1]² ─────────────
         * Vol1=(h/2)ΣF, Vol2=(h/2)ΣG feed the ξ-,η-slope eqns; the ξη eqn needs
         * ∂ₓ(ξη)=(2/h)η, ∂_y(ξη)=(2/h)ξ  ⇒  Vol3=(h/2)Σ(η·F+ξ·G).  The state is
         * the FULL bilinear field (ξη term included). */
        real Vol1[4]={0.,0.,0.,0.}, Vol2[4]={0.,0.,0.,0.}, Vol3[4]={0.,0.,0.,0.};
#pragma unroll
        for (int m = 0; m < 4; m++) {
            real xi  = (m & 1) ? GP : -GP;
            real eta = (m & 2) ? GP : -GP;
            real xe  = xi * eta;
            real Fq[4], Gq[4];
            euler_FG_pt(Uc[Q_RHO] + Uc[Q_RHO_SX]*xi + Uc[Q_RHO_SY]*eta + Uc[Q_RHO_SXY]*xe,
                        Uc[Q_MXA] + Uc[Q_MXS ]  *xi + Uc[Q_MXSY ]*eta + Uc[Q_MXSXY ]*xe,
                        Uc[Q_MYA] + Uc[Q_MYSX]  *xi + Uc[Q_MYS  ]*eta + Uc[Q_MYSXY ]*xe,
                        Uc[Q_E  ] + Uc[Q_E_SX]  *xi + Uc[Q_E_SY ]*eta + Uc[Q_E_SXY ]*xe,
                        Fq, Gq);
            for (int q = 0; q < 4; q++) {
                Vol1[q] += Fq[q];
                Vol2[q] += Gq[q];
                Vol3[q] += eta*Fq[q] + xi*Gq[q];
            }
        }
        for (int q = 0; q < 4; q++) { Vol1[q]*=0.5*h; Vol2[q]*=0.5*h; Vol3[q]*=0.5*h; }

        /* ── Assemble dU/dt onto {avg, ξ, η, ξη} of each conserved var.
         * M⁻¹ folded in: avg×(1/h²), ξ,η-slope×(3/h²), ξη×(9/h²)
         * [M_avg=h², M_ξ=M_η=h²/3, M_ξη=h²/9].  On the right face φ_ξη=ξη→η
         * (moment IeR); on the top face φ_ξη→ξ (moment JxT). ──────────────── */
#pragma unroll
        for (int q = 0; q < 4; q++) {
            /* self cell k: owns right (ξ=+1) & top (η=+1) faces + volume */
            atomicAdd(&RHS[A  [q]*N2 + k], (-I0R[q] - I0T[q]) * inv_h2);
            atomicAdd(&RHS[SX [q]*N2 + k], (3.*Vol1[q] - 3.*I0R[q] - 3.*JxT[q]) * inv_h2);
            atomicAdd(&RHS[SY [q]*N2 + k], (3.*Vol2[q] - 3.*IeR[q] - 3.*I0T[q]) * inv_h2);
            atomicAdd(&RHS[SXY[q]*N2 + k], (9.*Vol3[q] - 9.*IeR[q] - 9.*JxT[q]) * inv_h2);
            /* right neighbour kR: shared face is its left (ξ=−1) face */
            atomicAdd(&RHS[A  [q]*N2 + kR],       I0R[q] * inv_h2);
            atomicAdd(&RHS[SX [q]*N2 + kR], -3.*I0R[q] * inv_h2);
            atomicAdd(&RHS[SY [q]*N2 + kR],  3.*IeR[q] * inv_h2);
            atomicAdd(&RHS[SXY[q]*N2 + kR], -9.*IeR[q] * inv_h2);
            /* top neighbour kT: shared face is its bottom (η=−1) face */
            atomicAdd(&RHS[A  [q]*N2 + kT],       I0T[q] * inv_h2);
            atomicAdd(&RHS[SX [q]*N2 + kT],  3.*JxT[q] * inv_h2);
            atomicAdd(&RHS[SY [q]*N2 + kT], -3.*I0T[q] * inv_h2);
            atomicAdd(&RHS[SXY[q]*N2 + kT], -9.*JxT[q] * inv_h2);
        }

        /* ── CFL spectral radius (cell-average state) ────────────────────── */
        real rho_a  = fmax(Uc[Q_RHO], (real)RHOFLOOR);
        real inv_ra = 1. / rho_a;
        real u_avg  = Uc[Q_MXA] * inv_ra;
        real v_avg  = Uc[Q_MYA] * inv_ra;
        real pa = pressure_from_rt0(rho_a, Uc[Q_MXA], Uc[Q_MYA], Uc[Q_E]);
        real cs = sound_speed(rho_a, fmax(pa, (real)PFLOOR));
        lam_out[k] = fabs(u_avg) + fabs(v_avg) + cs;
    }
}

/* ═══════════════════════════════════════════════════════════════════════════
 * RHS kernel — RECOVERY-TRACE variant (separate kernel, same 2-cell stencil)
 *
 * Identical to compute_rhs except for the face traces fed to the Riemann
 * solver.  The 4 per-face DOFs (Vc,sc,Vn,sn) determine a UNIQUE cubic (van
 * Leer recovery); its face value is
 *     u*   = (Vc+Vn)/2 − (sn−sc)/3                     (P3-exact, O(h⁴))
 * and the jump of the two one-sided Hermite quadratics is
 *     jmp  = UL_h − UR_h = (2/3)[(Vc−Vn) + (sc+sn)]    (O(h³) smooth)
 * The Riemann solver is fed the NARROWED pair
 *     UL = u* + (τ/2)·jmp ,   UR = u* − (τ/2)·jmp
 * τ=1 reproduces compute_rhs exactly (UL=UL_h, UR=UR_h);  τ→0 is the pure
 * central recovery flux: same 4th-order cell-average order, ~6× smaller
 * dispersion constant (θ⁵/1080 vs θ⁵/180), but NO face dissipation — the
 * spurious slope mode is then undamped, so keep τ ≥ ~0.05 (or blend with a
 * shock sensor) for nonlinear runs.
 *
 * Select at compile time:  nvcc ... -DRECOVERY_TAU=0.1
 * ═══════════════════════════════════════════════════════════════════════════ */
__global__ void
compute_rhs_recovery(const real * __restrict__ Qp,
                     real       * __restrict__ RHS,
                     real       * __restrict__ lam_out,
                     int N, int Np, real h, real tau)
{
    const int g   = G_GHOST;
    const int N2  = N * N;

    for (int k = blockIdx.x * blockDim.x + threadIdx.x; k < N2;
             k += blockDim.x * gridDim.x) {
        int j = k / N, i = k % N;
        int jp = j + g, ip = i + g;
        int kR = j * N + ((i + 1 < N) ? (i + 1) : 0);
        int kT = ((j + 1 < N) ? (j + 1) : 0) * N + i;

        real Uc[NVAR], Ur[NVAR], Ut[NVAR];
#pragma unroll
        for (int q = 0; q < NVAR; q++) {
            Uc[q] = Qp[IDX_P(q, jp,   ip,   Np)];
            Ur[q] = Qp[IDX_P(q, jp,   ip+1, Np)];
            Ut[q] = Qp[IDX_P(q, jp+1, ip,   Np)];
        }

        const real GP = (real)0.577350269189625764;   /* 1/√3 : 2-pt Gauss node */
        real inv_h2 = 1. / (h * h);

        const int A  [4] = { Q_RHO,     Q_MXA,   Q_MYA,   Q_E     };
        const int SX [4] = { Q_RHO_SX,  Q_MXS,   Q_MYSX,  Q_E_SX  };
        const int SY [4] = { Q_RHO_SY,  Q_MXSY,  Q_MYS,   Q_E_SY  };
        const int SXY[4] = { Q_RHO_SXY, Q_MXSXY, Q_MYSXY, Q_E_SXY };

        const real THIRD = (real)(1./3.), TT = (real)(2./3.);

        real I0R[4]={0.,0.,0.,0.}, IeR[4]={0.,0.,0.,0.};
        real I0T[4]={0.,0.,0.,0.}, JxT[4]={0.,0.,0.,0.};

        /* ── RIGHT face: recovery value ± (τ/2)·Hermite jump ─────────────── */
#pragma unroll
        for (int m = 0; m < 2; m++) {
            real eta = (m == 0) ? -GP : GP;
            real UL[4], UR[4];
#pragma unroll
            for (int v = 0; v < 4; v++) {
                real Vc = Uc[A[v]]  + Uc[SY[v]]*eta;
                real Vn = Ur[A[v]]  + Ur[SY[v]]*eta;
                real sc = Uc[SX[v]] + Uc[SXY[v]]*eta;
                real sn = Ur[SX[v]] + Ur[SXY[v]]*eta;
                real us  = (real)0.5*(Vc + Vn) - (sn - sc)*THIRD;  /* recovery */
                real jmp = TT*((Vc - Vn) + (sc + sn));             /* UL_h-UR_h */
                UL[v] = us + (real)0.5*tau*jmp;
                UR[v] = us - (real)0.5*tau*jmp;
            }
            real WL[4], WR[4];
            cons_pt_to_W(UL[0],UL[1],UL[2],UL[3], WL);
            cons_pt_to_W(UR[0],UR[1],UR[2],UR[3], WR);
            real Fg[4]; riemann_n(WL, WR, 1., 0., Fg);
            for (int q = 0; q < 4; q++) { I0R[q] += Fg[q]; IeR[q] += eta*Fg[q]; }
        }
        for (int q = 0; q < 4; q++) { I0R[q] *= 0.5*h; IeR[q] *= 0.5*h; }

        /* ── TOP face: recovery value ± (τ/2)·Hermite jump ───────────────── */
#pragma unroll
        for (int m = 0; m < 2; m++) {
            real xi = (m == 0) ? -GP : GP;
            real UB[4], UT[4];
#pragma unroll
            for (int v = 0; v < 4; v++) {
                real Vc = Uc[A[v]]  + Uc[SX[v]]*xi;
                real Vn = Ut[A[v]]  + Ut[SX[v]]*xi;
                real sc = Uc[SY[v]] + Uc[SXY[v]]*xi;
                real sn = Ut[SY[v]] + Ut[SXY[v]]*xi;
                real us  = (real)0.5*(Vc + Vn) - (sn - sc)*THIRD;
                real jmp = TT*((Vc - Vn) + (sc + sn));
                UB[v] = us + (real)0.5*tau*jmp;
                UT[v] = us - (real)0.5*tau*jmp;
            }
            real WB[4], WT[4];
            cons_pt_to_W(UB[0],UB[1],UB[2],UB[3], WB);
            cons_pt_to_W(UT[0],UT[1],UT[2],UT[3], WT);
            real Gg[4]; riemann_n(WB, WT, 0., 1., Gg);
            for (int q = 0; q < 4; q++) { I0T[q] += Gg[q]; JxT[q] += xi*Gg[q]; }
        }
        for (int q = 0; q < 4; q++) { I0T[q] *= 0.5*h; JxT[q] *= 0.5*h; }

        /* ── VOLUME: identical to compute_rhs ─────────────────────────────── */
        real Vol1[4]={0.,0.,0.,0.}, Vol2[4]={0.,0.,0.,0.}, Vol3[4]={0.,0.,0.,0.};
#pragma unroll
        for (int m = 0; m < 4; m++) {
            real xi  = (m & 1) ? GP : -GP;
            real eta = (m & 2) ? GP : -GP;
            real xe  = xi * eta;
            real Fq[4], Gq[4];
            euler_FG_pt(Uc[Q_RHO] + Uc[Q_RHO_SX]*xi + Uc[Q_RHO_SY]*eta + Uc[Q_RHO_SXY]*xe,
                        Uc[Q_MXA] + Uc[Q_MXS ]  *xi + Uc[Q_MXSY ]*eta + Uc[Q_MXSXY ]*xe,
                        Uc[Q_MYA] + Uc[Q_MYSX]  *xi + Uc[Q_MYS  ]*eta + Uc[Q_MYSXY ]*xe,
                        Uc[Q_E  ] + Uc[Q_E_SX]  *xi + Uc[Q_E_SY ]*eta + Uc[Q_E_SXY ]*xe,
                        Fq, Gq);
            for (int q = 0; q < 4; q++) {
                Vol1[q] += Fq[q];
                Vol2[q] += Gq[q];
                Vol3[q] += eta*Fq[q] + xi*Gq[q];
            }
        }
        for (int q = 0; q < 4; q++) { Vol1[q]*=0.5*h; Vol2[q]*=0.5*h; Vol3[q]*=0.5*h; }

        /* ── Assemble: identical to compute_rhs ───────────────────────────── */
#pragma unroll
        for (int q = 0; q < 4; q++) {
            atomicAdd(&RHS[A  [q]*N2 + k], (-I0R[q] - I0T[q]) * inv_h2);
            atomicAdd(&RHS[SX [q]*N2 + k], (3.*Vol1[q] - 3.*I0R[q] - 3.*JxT[q]) * inv_h2);
            atomicAdd(&RHS[SY [q]*N2 + k], (3.*Vol2[q] - 3.*IeR[q] - 3.*I0T[q]) * inv_h2);
            atomicAdd(&RHS[SXY[q]*N2 + k], (9.*Vol3[q] - 9.*IeR[q] - 9.*JxT[q]) * inv_h2);
            atomicAdd(&RHS[A  [q]*N2 + kR],       I0R[q] * inv_h2);
            atomicAdd(&RHS[SX [q]*N2 + kR], -3.*I0R[q] * inv_h2);
            atomicAdd(&RHS[SY [q]*N2 + kR],  3.*IeR[q] * inv_h2);
            atomicAdd(&RHS[SXY[q]*N2 + kR], -9.*IeR[q] * inv_h2);
            atomicAdd(&RHS[A  [q]*N2 + kT],       I0T[q] * inv_h2);
            atomicAdd(&RHS[SX [q]*N2 + kT],  3.*JxT[q] * inv_h2);
            atomicAdd(&RHS[SY [q]*N2 + kT], -3.*I0T[q] * inv_h2);
            atomicAdd(&RHS[SXY[q]*N2 + kT], -9.*JxT[q] * inv_h2);
        }

        real rho_a  = fmax(Uc[Q_RHO], (real)RHOFLOOR);
        real inv_ra = 1. / rho_a;
        real u_avg  = Uc[Q_MXA] * inv_ra;
        real v_avg  = Uc[Q_MYA] * inv_ra;
        real pa = pressure_from_rt0(rho_a, Uc[Q_MXA], Uc[Q_MYA], Uc[Q_E]);
        real cs = sound_speed(rho_a, fmax(pa, (real)PFLOOR));
        lam_out[k] = fabs(u_avg) + fabs(v_avg) + cs;
    }
}

/* ═══════════════════════════════════════════════════════════════════════════
 * RHS kernel — SENSOR-BLENDED recovery traces (separate kernel)
 *
 * Same recovery structure as compute_rhs_recovery, but τ is chosen PER FACE
 * GAUSS NODE from a jump sensor instead of being a global constant:
 *     sens = |jmp_ρ|/(|u*_ρ|+ε) + |jmp_E|/(|u*_E|+ε)
 *     τ    = τmax · min(1, sens/s0)
 * The Hermite-trace jump is O(h³) in smooth flow and O(1) at discontinuities,
 * so τ→0 (pure recovery: 4th order, ~3× smaller dispersion constant, zero
 * face dissipation) in smooth regions and τ→τmax at shocks/contacts.
 *
 * STABILITY: the fully discrete SSP-RK3 CFL limit depends on τ —
 *   ν_max(τ): 0.577 (τ=0), 0.774 (τ=0.5), 0.741 (τ=0.75), 0.610 (τ=1).
 * A blended run is bounded by the smallest limit over the τ range in use:
 * with τ ∈ [0, 0.75] that is ν=0.577 → run at CFL ≤ ~0.55 (e.g. -DCFL_DEFAULT=0.5).
 *
 * Select:  nvcc ... -DBLEND_TAU=0.75 [-DBLEND_S0=0.05] [-DCFL_DEFAULT=0.5]
 * ═══════════════════════════════════════════════════════════════════════════ */
__global__ void
compute_rhs_blend(const real * __restrict__ Qp,
                  real       * __restrict__ RHS,
                  real       * __restrict__ lam_out,
                  int N, int Np, real h, real tau_max, real s0)
{
    const int g   = G_GHOST;
    const int N2  = N * N;

    for (int k = blockIdx.x * blockDim.x + threadIdx.x; k < N2;
             k += blockDim.x * gridDim.x) {
        int j = k / N, i = k % N;
        int jp = j + g, ip = i + g;
        int kR = j * N + ((i + 1 < N) ? (i + 1) : 0);
        int kT = ((j + 1 < N) ? (j + 1) : 0) * N + i;

        real Uc[NVAR], Ur[NVAR], Ut[NVAR];
#pragma unroll
        for (int q = 0; q < NVAR; q++) {
            Uc[q] = Qp[IDX_P(q, jp,   ip,   Np)];
            Ur[q] = Qp[IDX_P(q, jp,   ip+1, Np)];
            Ut[q] = Qp[IDX_P(q, jp+1, ip,   Np)];
        }

        const real GP = (real)0.577350269189625764;
        real inv_h2 = 1. / (h * h);

        const int A  [4] = { Q_RHO,     Q_MXA,   Q_MYA,   Q_E     };
        const int SX [4] = { Q_RHO_SX,  Q_MXS,   Q_MYSX,  Q_E_SX  };
        const int SY [4] = { Q_RHO_SY,  Q_MXSY,  Q_MYS,   Q_E_SY  };
        const int SXY[4] = { Q_RHO_SXY, Q_MXSXY, Q_MYSXY, Q_E_SXY };

        const real THIRD = (real)(1./3.), TT = (real)(2./3.);
        const real SEPS  = (real)1e-12;

        real I0R[4]={0.,0.,0.,0.}, IeR[4]={0.,0.,0.,0.};
        real I0T[4]={0.,0.,0.,0.}, JxT[4]={0.,0.,0.,0.};

        /* ── RIGHT face: sensor-blended recovery traces ───────────────────── */
#pragma unroll
        for (int m = 0; m < 2; m++) {
            real eta = (m == 0) ? -GP : GP;
            real us[4], jmp[4], UL[4], UR[4];
#pragma unroll
            for (int v = 0; v < 4; v++) {
                real Vc = Uc[A[v]]  + Uc[SY[v]]*eta;
                real Vn = Ur[A[v]]  + Ur[SY[v]]*eta;
                real sc = Uc[SX[v]] + Uc[SXY[v]]*eta;
                real sn = Ur[SX[v]] + Ur[SXY[v]]*eta;
                us [v] = (real)0.5*(Vc + Vn) - (sn - sc)*THIRD;
                jmp[v] = TT*((Vc - Vn) + (sc + sn));
            }
            real sens = fabs(jmp[0])/(fabs(us[0])+SEPS)
                      + fabs(jmp[3])/(fabs(us[3])+SEPS);
            real tau  = tau_max * fmin((real)1., sens/s0);
#pragma unroll
            for (int v = 0; v < 4; v++) {
                UL[v] = us[v] + (real)0.5*tau*jmp[v];
                UR[v] = us[v] - (real)0.5*tau*jmp[v];
            }
            real WL[4], WR[4];
            cons_pt_to_W(UL[0],UL[1],UL[2],UL[3], WL);
            cons_pt_to_W(UR[0],UR[1],UR[2],UR[3], WR);
            real Fg[4]; riemann_n(WL, WR, 1., 0., Fg);
            for (int q = 0; q < 4; q++) { I0R[q] += Fg[q]; IeR[q] += eta*Fg[q]; }
        }
        for (int q = 0; q < 4; q++) { I0R[q] *= 0.5*h; IeR[q] *= 0.5*h; }

        /* ── TOP face: sensor-blended recovery traces ─────────────────────── */
#pragma unroll
        for (int m = 0; m < 2; m++) {
            real xi = (m == 0) ? -GP : GP;
            real us[4], jmp[4], UB[4], UT[4];
#pragma unroll
            for (int v = 0; v < 4; v++) {
                real Vc = Uc[A[v]]  + Uc[SX[v]]*xi;
                real Vn = Ut[A[v]]  + Ut[SX[v]]*xi;
                real sc = Uc[SY[v]] + Uc[SXY[v]]*xi;
                real sn = Ut[SY[v]] + Ut[SXY[v]]*xi;
                us [v] = (real)0.5*(Vc + Vn) - (sn - sc)*THIRD;
                jmp[v] = TT*((Vc - Vn) + (sc + sn));
            }
            real sens = fabs(jmp[0])/(fabs(us[0])+SEPS)
                      + fabs(jmp[3])/(fabs(us[3])+SEPS);
            real tau  = tau_max * fmin((real)1., sens/s0);
#pragma unroll
            for (int v = 0; v < 4; v++) {
                UB[v] = us[v] + (real)0.5*tau*jmp[v];
                UT[v] = us[v] - (real)0.5*tau*jmp[v];
            }
            real WB[4], WT[4];
            cons_pt_to_W(UB[0],UB[1],UB[2],UB[3], WB);
            cons_pt_to_W(UT[0],UT[1],UT[2],UT[3], WT);
            real Gg[4]; riemann_n(WB, WT, 0., 1., Gg);
            for (int q = 0; q < 4; q++) { I0T[q] += Gg[q]; JxT[q] += xi*Gg[q]; }
        }
        for (int q = 0; q < 4; q++) { I0T[q] *= 0.5*h; JxT[q] *= 0.5*h; }

        /* ── VOLUME + assembly: identical to compute_rhs ──────────────────── */
        real Vol1[4]={0.,0.,0.,0.}, Vol2[4]={0.,0.,0.,0.}, Vol3[4]={0.,0.,0.,0.};
#pragma unroll
        for (int m = 0; m < 4; m++) {
            real xi  = (m & 1) ? GP : -GP;
            real eta = (m & 2) ? GP : -GP;
            real xe  = xi * eta;
            real Fq[4], Gq[4];
            euler_FG_pt(Uc[Q_RHO] + Uc[Q_RHO_SX]*xi + Uc[Q_RHO_SY]*eta + Uc[Q_RHO_SXY]*xe,
                        Uc[Q_MXA] + Uc[Q_MXS ]  *xi + Uc[Q_MXSY ]*eta + Uc[Q_MXSXY ]*xe,
                        Uc[Q_MYA] + Uc[Q_MYSX]  *xi + Uc[Q_MYS  ]*eta + Uc[Q_MYSXY ]*xe,
                        Uc[Q_E  ] + Uc[Q_E_SX]  *xi + Uc[Q_E_SY ]*eta + Uc[Q_E_SXY ]*xe,
                        Fq, Gq);
            for (int q = 0; q < 4; q++) {
                Vol1[q] += Fq[q];
                Vol2[q] += Gq[q];
                Vol3[q] += eta*Fq[q] + xi*Gq[q];
            }
        }
        for (int q = 0; q < 4; q++) { Vol1[q]*=0.5*h; Vol2[q]*=0.5*h; Vol3[q]*=0.5*h; }

#pragma unroll
        for (int q = 0; q < 4; q++) {
            atomicAdd(&RHS[A  [q]*N2 + k], (-I0R[q] - I0T[q]) * inv_h2);
            atomicAdd(&RHS[SX [q]*N2 + k], (3.*Vol1[q] - 3.*I0R[q] - 3.*JxT[q]) * inv_h2);
            atomicAdd(&RHS[SY [q]*N2 + k], (3.*Vol2[q] - 3.*IeR[q] - 3.*I0T[q]) * inv_h2);
            atomicAdd(&RHS[SXY[q]*N2 + k], (9.*Vol3[q] - 9.*IeR[q] - 9.*JxT[q]) * inv_h2);
            atomicAdd(&RHS[A  [q]*N2 + kR],       I0R[q] * inv_h2);
            atomicAdd(&RHS[SX [q]*N2 + kR], -3.*I0R[q] * inv_h2);
            atomicAdd(&RHS[SY [q]*N2 + kR],  3.*IeR[q] * inv_h2);
            atomicAdd(&RHS[SXY[q]*N2 + kR], -9.*IeR[q] * inv_h2);
            atomicAdd(&RHS[A  [q]*N2 + kT],       I0T[q] * inv_h2);
            atomicAdd(&RHS[SX [q]*N2 + kT],  3.*JxT[q] * inv_h2);
            atomicAdd(&RHS[SY [q]*N2 + kT], -3.*I0T[q] * inv_h2);
            atomicAdd(&RHS[SXY[q]*N2 + kT], -9.*JxT[q] * inv_h2);
        }

        real rho_a  = fmax(Uc[Q_RHO], (real)RHOFLOOR);
        real inv_ra = 1. / rho_a;
        real u_avg  = Uc[Q_MXA] * inv_ra;
        real v_avg  = Uc[Q_MYA] * inv_ra;
        real pa = pressure_from_rt0(rho_a, Uc[Q_MXA], Uc[Q_MYA], Uc[Q_E]);
        real cs = sound_speed(rho_a, fmax(pa, (real)PFLOOR));
        lam_out[k] = fabs(u_avg) + fabs(v_avg) + cs;
    }
}

/* ═══════════════════════════════════════════════════════════════════════════
 * RHS kernel — GENUINELY MULTIDIMENSIONAL Osher-type corner flux
 * (ported from ../wavelet_cfd CompressibleSolver multiDRhsKernel;
 *  Gaburro, Ricchiuto & Dumbser, arXiv:2506.00207, Sec. 3-4)
 *
 * Face integral by SIMPSON assembly (Balsara-style):
 *     ∮_face F ≈ h[ (1/6)F̂_corner1 + (2/3)F̂_mid + (1/6)F̂_corner2 ]
 *
 *   • CORNER flux tensors F̂x,F̂y at each vertex from 4 QUADRANT STATES built
 *     by the TENSOR-PRODUCT RECOVERY operator: 1D recovery+τ-narrowing in x
 *     of the (A,SX) and (SY,SXY) mode pairs per cell row, then recovery+
 *     narrowing in y across the rows.  τ=0 collapses all four states to the
 *     unique P3⊗P3 recovery corner value; τ>0 opens O(h³) recovery jumps
 *     (NOT the O(h) P0 jumps of the FV original) so the Green-Gauss
 *     dissipation −(1/2)|A_d(Q̄)|(h∇Q)_d stays high-order in smooth flow.
 *     Each vertex is evaluated identically by all sharing cells → conservative.
 *   • MIDPOINT fluxes: the recovery-narrowed face traces at η=0 through the
 *     1D Osher-Solomon solver (3-pt Gauss path integral of |A_d| — a 1-pt
 *     rule degrades to Roe and admits expansion shocks).
 *   • Volume terms and DG weak-form slope/cross updates unchanged; the face
 *     η-moment under Simpson is  (h/6)(F̂_c2 − F̂_c1).
 *
 * Select:  nvcc ... -DMD_OSHER_TAU=0.5 [-DCFL_DEFAULT=0.5]
 * ═══════════════════════════════════════════════════════════════════════════ */

/* |A_n|·dQ for 2D Euler at mean primitive state (rb,ub,vb,pb), normal (nx,ny) */
__device__ void
md_absjac_dq(real rb, real ub, real vb, real pb,
             real nx, real ny, const real dQ[4], real out[4])
{
    real cb = sound_speed(rb, pb);
    real q2 = ub*ub + vb*vb;
    real Hb = (real)0.5*q2 + cb*cb/((real)GAMMA_V - 1.);

    real d0 = dQ[0];
    real du = (dQ[1] - ub*d0)/rb;
    real dv = (dQ[2] - vb*d0)/rb;
    real dp = ((real)GAMMA_V - 1.)*(dQ[3] - (real)0.5*q2*d0 - rb*(ub*du + vb*dv));

    real un  = ub*nx + vb*ny;
    real dun = du*nx + dv*ny;
    real am  = (dp - rb*cb*dun)/((real)2.*cb*cb);
    real ap  = (dp + rb*cb*dun)/((real)2.*cb*cb);
    real ae  = d0 - dp/(cb*cb);
    real dtx = du - dun*nx, dty = dv - dun*ny;

    real l1 = fabs(un - cb), l2 = fabs(un), l3 = fabs(un + cb);

    out[0] = l1*am              + l2*ae                 + l3*ap;
    out[1] = l1*am*(ub - cb*nx) + l2*(ae*ub + rb*dtx)   + l3*ap*(ub + cb*nx);
    out[2] = l1*am*(vb - cb*ny) + l2*(ae*vb + rb*dty)   + l3*ap*(vb + cb*ny);
    out[3] = l1*am*(Hb - cb*un) + l2*(ae*(real)0.5*q2 + rb*(ub*dtx + vb*dty))
                                + l3*ap*(Hb + cb*un);
}

/* 1D Osher-Solomon flux, direction d (0=x,1=y); conserved inputs UL,UR.
 * f̂ = ½(fL+fR) − ½[∫₀¹|A_d(ψ(ξ))|dξ](UR−UL), straight path, 3-pt Gauss. */
__device__ void
md_osher1d(const real UL[4], const real UR[4], int d, real F[4])
{
    real dQ[4];
    for (int n = 0; n < 4; n++) dQ[n] = UR[n] - UL[n];
    real WL[4], WR[4];
    cons_pt_to_W(UL[0],UL[1],UL[2],UL[3], WL);
    cons_pt_to_W(UR[0],UR[1],UR[2],UR[3], WR);
    real fL[4], gL[4], fR[4], gR[4];
    euler_FG_pt(UL[0],UL[1],UL[2],UL[3], fL, gL);
    euler_FG_pt(UR[0],UR[1],UR[2],UR[3], fR, gR);
    const real *fLd = (d == 0) ? fL : gL;
    const real *fRd = (d == 0) ? fR : gR;

    const real xg[3] = {(real)0.11270166537925831, (real)0.5, (real)0.88729833462074169};
    const real wg[3] = {(real)(5./18.), (real)(8./18.), (real)(5./18.)};
    real nx = (d == 0) ? (real)1. : (real)0., ny = (real)1. - nx;
    real acc[4] = {0.,0.,0.,0.};
    for (int g = 0; g < 3; g++) {
        real Qg[4], Wg[4];
        for (int n = 0; n < 4; n++) Qg[n] = UL[n] + xg[g]*dQ[n];
        cons_pt_to_W(Qg[0],Qg[1],Qg[2],Qg[3], Wg);
        real tmp[4];
        md_absjac_dq(Wg[0],Wg[1],Wg[2],Wg[3], nx, ny, dQ, tmp);
        for (int n = 0; n < 4; n++) acc[n] += wg[g]*tmp[n];
    }
    for (int n = 0; n < 4; n++) F[n] = (real)0.5*(fLd[n] + fRd[n]) - (real)0.5*acc[n];
}

/* Tensor-product recovery quadrant states at the vertex BETWEEN cells
 * (ipL,jpB) (ipL+1,jpB) (ipL,jpB+1) (ipL+1,jpB+1)  [padded indices].
 * Qm[q][v]: quadrant q = LB,RB,LT,RT (x fastest), conserved var v. */
__device__ void
md_corner_states(const real * __restrict__ Qp, int Np,
                 int jpB, int ipL, real taux, real tauy, real Qm[4][4])
{
    const int A_  [4] = { Q_RHO,     Q_MXA,   Q_MYA,   Q_E     };
    const int SX_ [4] = { Q_RHO_SX,  Q_MXS,   Q_MYSX,  Q_E_SX  };
    const int SY_ [4] = { Q_RHO_SY,  Q_MXSY,  Q_MYS,   Q_E_SY  };
    const int SXY_[4] = { Q_RHO_SXY, Q_MXSXY, Q_MYSXY, Q_E_SXY };
    const real THIRD = (real)(1./3.), TT = (real)(2./3.);

    for (int v = 0; v < 4; v++) {
        real ULa[2], UL1[2], URa[2], UR1[2];    /* per row: value & η-slope */
        for (int r = 0; r < 2; r++) {
            int jp = jpB + r;
            real aL  = Qp[IDX_P(A_  [v], jp, ipL,   Np)];
            real sxL = Qp[IDX_P(SX_ [v], jp, ipL,   Np)];
            real syL = Qp[IDX_P(SY_ [v], jp, ipL,   Np)];
            real xyL = Qp[IDX_P(SXY_[v], jp, ipL,   Np)];
            real aR  = Qp[IDX_P(A_  [v], jp, ipL+1, Np)];
            real sxR = Qp[IDX_P(SX_ [v], jp, ipL+1, Np)];
            real syR = Qp[IDX_P(SY_ [v], jp, ipL+1, Np)];
            real xyR = Qp[IDX_P(SXY_[v], jp, ipL+1, Np)];
            /* x-recovery of the value pair (A,SX) and η-slope pair (SY,SXY) */
            real v0 = (real)0.5*(aL + aR)   - (sxR - sxL)*THIRD;
            real j0 = TT*((aL - aR)   + (sxL + sxR));
            real v1 = (real)0.5*(syL + syR) - (xyR - xyL)*THIRD;
            real j1 = TT*((syL - syR) + (xyL + xyR));
            ULa[r] = v0 + (real)0.5*taux*j0;  UL1[r] = v1 + (real)0.5*taux*j1;
            URa[r] = v0 - (real)0.5*taux*j0;  UR1[r] = v1 - (real)0.5*taux*j1;
        }
        /* y-recovery across the rows (row 0 below vertex, row 1 above) */
        real usL = (real)0.5*(ULa[0] + ULa[1]) - (UL1[1] - UL1[0])*THIRD;
        real jmL = TT*((ULa[0] - ULa[1]) + (UL1[0] + UL1[1]));
        real usR = (real)0.5*(URa[0] + URa[1]) - (UR1[1] - UR1[0])*THIRD;
        real jmR = TT*((URa[0] - URa[1]) + (UR1[0] + UR1[1]));
        Qm[0][v] = usL + (real)0.5*tauy*jmL;   /* LB */
        Qm[1][v] = usR + (real)0.5*tauy*jmR;   /* RB */
        Qm[2][v] = usL - (real)0.5*tauy*jmL;   /* LT */
        Qm[3][v] = usR - (real)0.5*tauy*jmR;   /* RT */
    }
}

/* Corner flux tensors F̂x,F̂y at the vertex between the 2×2 cell block whose
 * lower-left padded cell is (jpB,ipL) — wavelet_cfd mdCornerFlux structure:
 * central ¼Σf(Q_m) − ½|A_d(Q̄)|·(Green-Gauss h∇Q)_d, with recovery states. */
__device__ void
md_corner_flux(const real * __restrict__ Qp, int Np,
               int jpB, int ipL, real tau, real Fx[4], real Fy[4])
{
    /* per-direction quadrant states (wavelet_cfd mdCornerFlux structure):
     * the direction-d tensor opens recovery jumps only along d; the
     * tangential direction stays pure (central) recovery.  The fully-coupled
     * variant closes a slope<->slope cross-advection loop and is weakly
     * unstable, exactly as documented in the FV original.                  */
    real Qmx[4][4], Qmy[4][4];
    md_corner_states(Qp, Np, jpB, ipL, tau, 0.,  Qmx);
    md_corner_states(Qp, Np, jpB, ipL, 0.,  tau, Qmy);

    real fx[4][4], fy[4][4];
    real hgx[4], hgy[4], Qbx[4], Qby[4];
    for (int m = 0; m < 4; m++) {
        real fdum[4];
        euler_FG_pt(Qmx[m][0],Qmx[m][1],Qmx[m][2],Qmx[m][3], fx[m], fdum);
        euler_FG_pt(Qmy[m][0],Qmy[m][1],Qmy[m][2],Qmy[m][3], fdum, fy[m]);
    }
    for (int n = 0; n < 4; n++) {
        hgx[n] = (real)0.5*((Qmx[1][n] - Qmx[0][n]) + (Qmx[3][n] - Qmx[2][n]));
        hgy[n] = (real)0.5*((Qmy[2][n] - Qmy[0][n]) + (Qmy[3][n] - Qmy[1][n]));
        Qbx[n] = (real)0.25*(Qmx[0][n] + Qmx[1][n] + Qmx[2][n] + Qmx[3][n]);
        Qby[n] = (real)0.25*(Qmy[0][n] + Qmy[1][n] + Qmy[2][n] + Qmy[3][n]);
    }
    real Wbx[4], Wby[4];
    cons_pt_to_W(Qbx[0],Qbx[1],Qbx[2],Qbx[3], Wbx);
    cons_pt_to_W(Qby[0],Qby[1],Qby[2],Qby[3], Wby);
    real ax[4], ay[4];
    md_absjac_dq(Wbx[0],Wbx[1],Wbx[2],Wbx[3], 1., 0., hgx, ax);
    md_absjac_dq(Wby[0],Wby[1],Wby[2],Wby[3], 0., 1., hgy, ay);

    /* ½ dissipation coefficient: calibrated so the 1D (tangential-uniform)
     * limit is exactly ½(fL+fR) − ½|A|ΔQ on the narrowed states.          */
    for (int n = 0; n < 4; n++) {
        real cx = 0., cy = 0.;
        for (int m = 0; m < 4; m++) { cx += fx[m][n]; cy += fy[m][n]; }
        Fx[n] = (real)0.25*cx - (real)0.5*ax[n];
        Fy[n] = (real)0.25*cy - (real)0.5*ay[n];
    }
}

__global__ void
compute_rhs_mdosher(const real * __restrict__ Qp,
                    real       * __restrict__ RHS,
                    real       * __restrict__ lam_out,
                    int N, int Np, real h, real tau)
{
    const int g   = G_GHOST;
    const int N2  = N * N;

    for (int k = blockIdx.x * blockDim.x + threadIdx.x; k < N2;
             k += blockDim.x * gridDim.x) {
        int j = k / N, i = k % N;
        int jp = j + g, ip = i + g;
        int kR = j * N + ((i + 1 < N) ? (i + 1) : 0);
        int kT = ((j + 1 < N) ? (j + 1) : 0) * N + i;

        const int A  [4] = { Q_RHO,     Q_MXA,   Q_MYA,   Q_E     };
        const int SX [4] = { Q_RHO_SX,  Q_MXS,   Q_MYSX,  Q_E_SX  };
        const int SY [4] = { Q_RHO_SY,  Q_MXSY,  Q_MYS,   Q_E_SY  };
        const int SXY[4] = { Q_RHO_SXY, Q_MXSXY, Q_MYSXY, Q_E_SXY };
        const real THIRD = (real)(1./3.), TT = (real)(2./3.);
        const real GP = (real)0.577350269189625764;
        real inv_h2 = 1. / (h * h);

        real Uc[NVAR], Ur[NVAR], Ut[NVAR];
#pragma unroll
        for (int q = 0; q < NVAR; q++) {
            Uc[q] = Qp[IDX_P(q, jp,   ip,   Np)];
            Ur[q] = Qp[IDX_P(q, jp,   ip+1, Np)];
            Ut[q] = Qp[IDX_P(q, jp+1, ip,   Np)];
        }

        /* ── corner flux tensors at the right-face and top-face vertices ──
         * right face endpoints: vertices (i+½, j∓½) → blocks (jp-1,ip),(jp,ip)
         * top   face endpoints: vertices (i∓½, j+½) → blocks (jp,ip-1),(jp,ip)
         * each evaluated identically by all sharing cells → conservative.  */
        real FxRB[4], FyRB[4], FxRT[4], FyRT[4];   /* right-face corners */
        real FxTL[4], FyTL[4];                     /* top-left vertex     */
        md_corner_flux(Qp, Np, jp-1, ip,   tau, FxRB, FyRB);
        md_corner_flux(Qp, Np, jp,   ip,   tau, FxRT, FyRT);   /* shared TR */
        md_corner_flux(Qp, Np, jp,   ip-1, tau, FxTL, FyTL);

        /* ── face MIDPOINT fluxes: recovery-narrowed traces at η(ξ)=0
         *    through the 1D Osher-Solomon solver ─────────────────────────── */
        real FmR[4], FmT[4];
        {
            real UL[4], UR[4];
            for (int v = 0; v < 4; v++) {
                real us  = (real)0.5*(Uc[A[v]] + Ur[A[v]]) - (Ur[SX[v]] - Uc[SX[v]])*THIRD;
                real jmp = TT*((Uc[A[v]] - Ur[A[v]]) + (Uc[SX[v]] + Ur[SX[v]]));
                UL[v] = us + (real)0.5*tau*jmp;
                UR[v] = us - (real)0.5*tau*jmp;
            }
            md_osher1d(UL, UR, 0, FmR);
            for (int v = 0; v < 4; v++) {
                real us  = (real)0.5*(Uc[A[v]] + Ut[A[v]]) - (Ut[SY[v]] - Uc[SY[v]])*THIRD;
                real jmp = TT*((Uc[A[v]] - Ut[A[v]]) + (Uc[SY[v]] + Ut[SY[v]]));
                UL[v] = us + (real)0.5*tau*jmp;
                UR[v] = us - (real)0.5*tau*jmp;
            }
            md_osher1d(UL, UR, 1, FmT);
        }

        /* ── Simpson face integrals: I0 = h[(1/6)(Fc1+Fc2) + (2/3)Fm],
         *    η-moment Ie = (h/6)(Fc2 − Fc1)  (corners at η=∓1) ───────────── */
        const real W6 = (real)(1./6.), W23 = (real)(2./3.);
        real I0R[4], IeR[4], I0T[4], JxT[4];
        for (int q = 0; q < 4; q++) {
            I0R[q] = h*(W6*(FxRB[q] + FxRT[q]) + W23*FmR[q]);
            IeR[q] = h*W6*(FxRT[q] - FxRB[q]);
            I0T[q] = h*(W6*(FyTL[q] + FyRT[q]) + W23*FmT[q]);
            JxT[q] = h*W6*(FyRT[q] - FyTL[q]);
        }

        /* ── VOLUME + assembly: identical to compute_rhs ──────────────────── */
        real Vol1[4]={0.,0.,0.,0.}, Vol2[4]={0.,0.,0.,0.}, Vol3[4]={0.,0.,0.,0.};
#pragma unroll
        for (int m = 0; m < 4; m++) {
            real xi  = (m & 1) ? GP : -GP;
            real eta = (m & 2) ? GP : -GP;
            real xe  = xi * eta;
            real Fq[4], Gq[4];
            euler_FG_pt(Uc[Q_RHO] + Uc[Q_RHO_SX]*xi + Uc[Q_RHO_SY]*eta + Uc[Q_RHO_SXY]*xe,
                        Uc[Q_MXA] + Uc[Q_MXS ]  *xi + Uc[Q_MXSY ]*eta + Uc[Q_MXSXY ]*xe,
                        Uc[Q_MYA] + Uc[Q_MYSX]  *xi + Uc[Q_MYS  ]*eta + Uc[Q_MYSXY ]*xe,
                        Uc[Q_E  ] + Uc[Q_E_SX]  *xi + Uc[Q_E_SY ]*eta + Uc[Q_E_SXY ]*xe,
                        Fq, Gq);
            for (int q = 0; q < 4; q++) {
                Vol1[q] += Fq[q];
                Vol2[q] += Gq[q];
                Vol3[q] += eta*Fq[q] + xi*Gq[q];
            }
        }
        for (int q = 0; q < 4; q++) { Vol1[q]*=0.5*h; Vol2[q]*=0.5*h; Vol3[q]*=0.5*h; }

#pragma unroll
        for (int q = 0; q < 4; q++) {
            atomicAdd(&RHS[A  [q]*N2 + k], (-I0R[q] - I0T[q]) * inv_h2);
            atomicAdd(&RHS[SX [q]*N2 + k], (3.*Vol1[q] - 3.*I0R[q] - 3.*JxT[q]) * inv_h2);
            atomicAdd(&RHS[SY [q]*N2 + k], (3.*Vol2[q] - 3.*IeR[q] - 3.*I0T[q]) * inv_h2);
            atomicAdd(&RHS[SXY[q]*N2 + k], (9.*Vol3[q] - 9.*IeR[q] - 9.*JxT[q]) * inv_h2);
            atomicAdd(&RHS[A  [q]*N2 + kR],       I0R[q] * inv_h2);
            atomicAdd(&RHS[SX [q]*N2 + kR], -3.*I0R[q] * inv_h2);
            atomicAdd(&RHS[SY [q]*N2 + kR],  3.*IeR[q] * inv_h2);
            atomicAdd(&RHS[SXY[q]*N2 + kR], -9.*IeR[q] * inv_h2);
            atomicAdd(&RHS[A  [q]*N2 + kT],       I0T[q] * inv_h2);
            atomicAdd(&RHS[SX [q]*N2 + kT],  3.*JxT[q] * inv_h2);
            atomicAdd(&RHS[SY [q]*N2 + kT], -3.*I0T[q] * inv_h2);
            atomicAdd(&RHS[SXY[q]*N2 + kT], -9.*JxT[q] * inv_h2);
        }

        real rho_a  = fmax(Uc[Q_RHO], (real)RHOFLOOR);
        real inv_ra = 1. / rho_a;
        real u_avg  = Uc[Q_MXA] * inv_ra;
        real v_avg  = Uc[Q_MYA] * inv_ra;
        real pa = pressure_from_rt0(rho_a, Uc[Q_MXA], Uc[Q_MYA], Uc[Q_E]);
        real cs = sound_speed(rho_a, fmax(pa, (real)PFLOOR));
        lam_out[k] = fabs(u_avg) + fabs(v_avg) + cs;
    }
}

/* ═══════════════════════════════════════════════════════════════════════════
 * RHS kernel — DUCROS-GATED recovery traces (separate kernel)
 *
 * Same recovery structure as compute_rhs_recovery, but τ per FACE comes from
 * a Ducros-like dilatation sensor at the two adjacent cell centres (from the
 * Q1 slopes):
 *     divU = ∂ₓu + ∂_y v ,   sens = divU² / (divU² + 0.1·c²/h²)
 *     τ_face = τmax · max(sens_L, sens_R)
 * Fires only where the compression rate |divU| rivals the acoustic rate c/h
 * (genuine shocks); vortical / low-Mach smooth flow has divU ≈ 0 → τ ≈ 0
 * (pure recovery flux).  UNION with the jump sensor of
 * compute_rhs_blend:  τ_face = τmax·max( min(1, s_jump/s0), sens_ducros ) —
 * the jump term catches t=0 discontinuities (divU has not developed yet) and
 * floors the marginal-mode noise; the dilatation term catches developed
 * shocks Mach-consistently and stays dark in vortical flow.
 *
 * Select:  nvcc ... -DDUCROS_TAU=0.75 [-DCFL_DEFAULT=0.5]
 * (τ→0 in smooth regions ⇒ the ν=0.577 imaginary-axis CFL limit applies)
 * ═══════════════════════════════════════════════════════════════════════════ */
#ifndef DUCROS_TAU_RE
#  define DUCROS_TAU_RE 0.5   /* fixed tau for rho & E in compute_rhs_ducros */
#endif

__device__ __forceinline__ real
ducros_sensor(const real U[NVAR], real h)
{
    real r  = fmax(U[Q_RHO], (real)RHOFLOOR);
    real iu = 1. / r;
    real u  = U[Q_MXA] * iu, v = U[Q_MYA] * iu;
    /* cell-centre velocity divergence from the modal slopes (chain rule) */
    real dxu = (2./h) * (U[Q_MXS] - u*U[Q_RHO_SX]) * iu;
    real dyv = (2./h) * (U[Q_MYS] - v*U[Q_RHO_SY]) * iu;
    real div = dxu + dyv;
    real p   = pressure_from_rt0(r, U[Q_MXA], U[Q_MYA], U[Q_E]);
    real c2  = (real)GAMMA_V * fmax(p, (real)PFLOOR) * iu;
    real d2  = div*div;
    return d2 / (d2 + (real)0.1*c2/(h*h) + (real)1e-30);
}


__global__ void
compute_rhs_ducros(const real * __restrict__ Qp,
                  real       * __restrict__ RHS,
                  real       * __restrict__ lam_out,
                  int N, int Np, real h, real tau_max, real s0)
{
    const int g   = G_GHOST;
    const int N2  = N * N;

    for (int k = blockIdx.x * blockDim.x + threadIdx.x; k < N2;
             k += blockDim.x * gridDim.x) {
        int j = k / N, i = k % N;
        int jp = j + g, ip = i + g;
        int kR = j * N + ((i + 1 < N) ? (i + 1) : 0);
        int kT = ((j + 1 < N) ? (j + 1) : 0) * N + i;

        real Uc[NVAR], Ur[NVAR], Ut[NVAR];
#pragma unroll
        for (int q = 0; q < NVAR; q++) {
            Uc[q] = Qp[IDX_P(q, jp,   ip,   Np)];
            Ur[q] = Qp[IDX_P(q, jp,   ip+1, Np)];
            Ut[q] = Qp[IDX_P(q, jp+1, ip,   Np)];
        }

        /* Ducros dilatation sensors of the two face-sharing cells */
        real dsc  = ducros_sensor(Uc, h);
        real ducR = fmax(dsc, ducros_sensor(Ur, h));
        real ducT = fmax(dsc, ducros_sensor(Ut, h));

        const real GP = (real)0.577350269189625764;
        real inv_h2 = 1. / (h * h);

        const int A  [4] = { Q_RHO,     Q_MXA,   Q_MYA,   Q_E     };
        const int SX [4] = { Q_RHO_SX,  Q_MXS,   Q_MYSX,  Q_E_SX  };
        const int SY [4] = { Q_RHO_SY,  Q_MXSY,  Q_MYS,   Q_E_SY  };
        const int SXY[4] = { Q_RHO_SXY, Q_MXSXY, Q_MYSXY, Q_E_SXY };

        const real THIRD = (real)(1./3.), TT = (real)(2./3.);
        const real SEPS  = (real)1e-12;

        real I0R[4]={0.,0.,0.,0.}, IeR[4]={0.,0.,0.,0.};
        real I0T[4]={0.,0.,0.,0.}, JxT[4]={0.,0.,0.,0.};

        /* ── RIGHT face: sensor-blended recovery traces ───────────────────── */
#pragma unroll
        for (int m = 0; m < 2; m++) {
            real eta = (m == 0) ? -GP : GP;
            real us[4], jmp[4], UL[4], UR[4];
#pragma unroll
            for (int v = 0; v < 4; v++) {
                real Vc = Uc[A[v]]  + Uc[SY[v]]*eta;
                real Vn = Ur[A[v]]  + Ur[SY[v]]*eta;
                real sc = Uc[SX[v]] + Uc[SXY[v]]*eta;
                real sn = Ur[SX[v]] + Ur[SXY[v]]*eta;
                us [v] = (real)0.5*(Vc + Vn) - (sn - sc)*THIRD;
                jmp[v] = TT*((Vc - Vn) + (sc + sn));
            }
            real sens = fabs(jmp[0])/(fabs(us[0])+SEPS)
                      + fabs(jmp[3])/(fabs(us[3])+SEPS);
            real tau  = tau_max * fmax(fmin((real)1., sens/s0), ducR);
#pragma unroll
            for (int v = 0; v < 4; v++) {
                /* sensor-gated tau on MOMENTUM only; rho & E keep a fixed
                 * tau floor (Mach-independent acoustic/entropy damping)   */
                real tv = (v == 1 || v == 2) ? tau : (real)DUCROS_TAU_RE;
                UL[v] = us[v] + (real)0.5*tv*jmp[v];
                UR[v] = us[v] - (real)0.5*tv*jmp[v];
            }
            real WL[4], WR[4];
            cons_pt_to_W(UL[0],UL[1],UL[2],UL[3], WL);
            cons_pt_to_W(UR[0],UR[1],UR[2],UR[3], WR);
            real Fg[4]; riemann_n(WL, WR, 1., 0., Fg);
            for (int q = 0; q < 4; q++) { I0R[q] += Fg[q]; IeR[q] += eta*Fg[q]; }
        }
        for (int q = 0; q < 4; q++) { I0R[q] *= 0.5*h; IeR[q] *= 0.5*h; }

        /* ── TOP face: sensor-blended recovery traces ─────────────────────── */
#pragma unroll
        for (int m = 0; m < 2; m++) {
            real xi = (m == 0) ? -GP : GP;
            real us[4], jmp[4], UB[4], UT[4];
#pragma unroll
            for (int v = 0; v < 4; v++) {
                real Vc = Uc[A[v]]  + Uc[SX[v]]*xi;
                real Vn = Ut[A[v]]  + Ut[SX[v]]*xi;
                real sc = Uc[SY[v]] + Uc[SXY[v]]*xi;
                real sn = Ut[SY[v]] + Ut[SXY[v]]*xi;
                us [v] = (real)0.5*(Vc + Vn) - (sn - sc)*THIRD;
                jmp[v] = TT*((Vc - Vn) + (sc + sn));
            }
            real sens = fabs(jmp[0])/(fabs(us[0])+SEPS)
                      + fabs(jmp[3])/(fabs(us[3])+SEPS);
            real tau  = tau_max * fmax(fmin((real)1., sens/s0), ducT);
#pragma unroll
            for (int v = 0; v < 4; v++) {
                /* sensor-gated tau on MOMENTUM only; rho & E keep a fixed
                 * tau floor (Mach-independent acoustic/entropy damping)   */
                real tv = (v == 1 || v == 2) ? tau : (real)DUCROS_TAU_RE;
                UB[v] = us[v] + (real)0.5*tv*jmp[v];
                UT[v] = us[v] - (real)0.5*tv*jmp[v];
            }
            real WB[4], WT[4];
            cons_pt_to_W(UB[0],UB[1],UB[2],UB[3], WB);
            cons_pt_to_W(UT[0],UT[1],UT[2],UT[3], WT);
            real Gg[4]; riemann_n(WB, WT, 0., 1., Gg);
            for (int q = 0; q < 4; q++) { I0T[q] += Gg[q]; JxT[q] += xi*Gg[q]; }
        }
        for (int q = 0; q < 4; q++) { I0T[q] *= 0.5*h; JxT[q] *= 0.5*h; }

        /* ── VOLUME + assembly: identical to compute_rhs ──────────────────── */
        real Vol1[4]={0.,0.,0.,0.}, Vol2[4]={0.,0.,0.,0.}, Vol3[4]={0.,0.,0.,0.};
#pragma unroll
        for (int m = 0; m < 4; m++) {
            real xi  = (m & 1) ? GP : -GP;
            real eta = (m & 2) ? GP : -GP;
            real xe  = xi * eta;
            real Fq[4], Gq[4];
            euler_FG_pt(Uc[Q_RHO] + Uc[Q_RHO_SX]*xi + Uc[Q_RHO_SY]*eta + Uc[Q_RHO_SXY]*xe,
                        Uc[Q_MXA] + Uc[Q_MXS ]  *xi + Uc[Q_MXSY ]*eta + Uc[Q_MXSXY ]*xe,
                        Uc[Q_MYA] + Uc[Q_MYSX]  *xi + Uc[Q_MYS  ]*eta + Uc[Q_MYSXY ]*xe,
                        Uc[Q_E  ] + Uc[Q_E_SX]  *xi + Uc[Q_E_SY ]*eta + Uc[Q_E_SXY ]*xe,
                        Fq, Gq);
            for (int q = 0; q < 4; q++) {
                Vol1[q] += Fq[q];
                Vol2[q] += Gq[q];
                Vol3[q] += eta*Fq[q] + xi*Gq[q];
            }
        }
        for (int q = 0; q < 4; q++) { Vol1[q]*=0.5*h; Vol2[q]*=0.5*h; Vol3[q]*=0.5*h; }

#pragma unroll
        for (int q = 0; q < 4; q++) {
            atomicAdd(&RHS[A  [q]*N2 + k], (-I0R[q] - I0T[q]) * inv_h2);
            atomicAdd(&RHS[SX [q]*N2 + k], (3.*Vol1[q] - 3.*I0R[q] - 3.*JxT[q]) * inv_h2);
            atomicAdd(&RHS[SY [q]*N2 + k], (3.*Vol2[q] - 3.*IeR[q] - 3.*I0T[q]) * inv_h2);
            atomicAdd(&RHS[SXY[q]*N2 + k], (9.*Vol3[q] - 9.*IeR[q] - 9.*JxT[q]) * inv_h2);
            atomicAdd(&RHS[A  [q]*N2 + kR],       I0R[q] * inv_h2);
            atomicAdd(&RHS[SX [q]*N2 + kR], -3.*I0R[q] * inv_h2);
            atomicAdd(&RHS[SY [q]*N2 + kR],  3.*IeR[q] * inv_h2);
            atomicAdd(&RHS[SXY[q]*N2 + kR], -9.*IeR[q] * inv_h2);
            atomicAdd(&RHS[A  [q]*N2 + kT],       I0T[q] * inv_h2);
            atomicAdd(&RHS[SX [q]*N2 + kT],  3.*JxT[q] * inv_h2);
            atomicAdd(&RHS[SY [q]*N2 + kT], -3.*I0T[q] * inv_h2);
            atomicAdd(&RHS[SXY[q]*N2 + kT], -9.*JxT[q] * inv_h2);
        }

        real rho_a  = fmax(Uc[Q_RHO], (real)RHOFLOOR);
        real inv_ra = 1. / rho_a;
        real u_avg  = Uc[Q_MXA] * inv_ra;
        real v_avg  = Uc[Q_MYA] * inv_ra;
        real pa = pressure_from_rt0(rho_a, Uc[Q_MXA], Uc[Q_MYA], Uc[Q_E]);
        real cs = sound_speed(rho_a, fmax(pa, (real)PFLOOR));
        lam_out[k] = fabs(u_avg) + fabs(v_avg) + cs;
    }
}

/* ── RHS dispatch:
 *   -DDUCROS_TAU=<τmax>                   → Ducros-gated recovery kernel
 *   -DMD_OSHER_TAU=<τ>                    → multidimensional Osher kernel
 *   -DBLEND_TAU=<τmax> [-DBLEND_S0=<s0>]  → sensor-blended kernel
 *   -DRECOVERY_TAU=<τ>                    → constant-τ recovery kernel
 *   (neither)                             → original compute_rhs            */
#ifndef BLEND_S0
#  define BLEND_S0 0.05
#endif
#if defined(DUCROS_TAU)
#  define COMPUTE_RHS(Q)  compute_rhs_ducros<<<GS_NBLK, BLOCK1D>>>((Q), d_RHS, d_lam, N, Np, h, (real)(DUCROS_TAU), (real)(BLEND_S0))
#elif defined(MD_OSHER_TAU)
#  define COMPUTE_RHS(Q)  compute_rhs_mdosher<<<GS_NBLK, BLOCK1D>>>((Q), d_RHS, d_lam, N, Np, h, (real)(MD_OSHER_TAU))
#elif defined(BLEND_TAU)
#  define COMPUTE_RHS(Q)  compute_rhs_blend<<<GS_NBLK, BLOCK1D>>>((Q), d_RHS, d_lam, N, Np, h, (real)(BLEND_TAU), (real)(BLEND_S0))
#elif defined(RECOVERY_TAU)
#  define COMPUTE_RHS(Q)  compute_rhs_recovery<<<GS_NBLK, BLOCK1D>>>((Q), d_RHS, d_lam, N, Np, h, (real)(RECOVERY_TAU))
#else
#  define COMPUTE_RHS(Q)  compute_rhs<<<GS_NBLK, BLOCK1D>>>((Q), d_RHS, d_lam, N, Np, h)
#endif

/* ═══════════════════════════════════════════════════════════════════════════
 * IC kernel  — circular Sod shock tube
 * ═══════════════════════════════════════════════════════════════════════════ */
__global__ void
ic_kernel(real * __restrict__ Qp, int N, int Np, real h)
{
    const int g   = G_GHOST;
    const int N2  = N * N;

    for (int k = blockIdx.x * blockDim.x + threadIdx.x; k < N2;
             k += blockDim.x * gridDim.x) {
        int j = k / N, i = k % N;
        int jp = j + g, ip = i + g;

        real xc = (i + 0.5) * h, yc = (j + 0.5) * h;
        real r  = sqrt((xc - 0.5)*(xc - 0.5) + (yc - 0.5)*(yc - 0.5));

        /* Smooth tanh interface (width ~ h) */
        real delta = 0.0 * h;
        real phi   = 0.5 * (1. + tanh((0.25 - r) / delta));

        real rho = 0.125 * (1. - phi) + (real)SOD_RHOIN * phi;  /* SOD_RHOIN inside, 0.125 outside */
        real p   = 0.1   * (1. - phi) + (real)SOD_PIN   * phi;  /* SOD_PIN   inside, 0.1   outside */
        real E   = p / (GAMMA_V - 1.);     /* zero velocity IC */

        /* P1 modal DOFs: zero velocity → all momentum modes = 0 */
        Qp[IDX_P(Q_RHO,  jp, ip, Np)] = rho;
        Qp[IDX_P(Q_MXA,  jp, ip, Np)] = 0.;
        Qp[IDX_P(Q_MXS,  jp, ip, Np)] = 0.;
        Qp[IDX_P(Q_MXSY, jp, ip, Np)] = 0.;
        Qp[IDX_P(Q_MYA,  jp, ip, Np)] = 0.;
        Qp[IDX_P(Q_MYSX, jp, ip, Np)] = 0.;
        Qp[IDX_P(Q_MYS,  jp, ip, Np)] = 0.;
        Qp[IDX_P(Q_E,    jp, ip, Np)] = E;
    }
}

/* ═══════════════════════════════════════════════════════════════════════════
 * Shock stabilisation for the Hermite-reconstructed BILINEAR-Q1 DG (all fields Q1).
 * Zhang-Shu positivity limiter: keep ρ≥RHOFLOOR and p≥PFLOOR at the four cell
 * corners (ρ linear, p concave ⇒ extrema at vertices), which makes every
 * interior/face quadrature point positive and kills the negative-pressure NaN
 * blow-up at shocks.  A NO-OP in smooth / low-Mach flow.
 * ═══════════════════════════════════════════════════════════════════════════ */
#ifndef SHOCK_LIMIT
#define SHOCK_LIMIT 1
#endif

/* Pressure of the P1 corner state at (ξ,η) from conserved averages+slopes. */
__device__ __forceinline__ reall
corner_pressure(reall ra, reall rsx, reall rsy, reall rsxy,
                reall mxa, reall mxsx, reall mxsy, reall mxsxy,
                reall mya, reall mysx, reall mys, reall mysxy,
                reall Ea, reall Esx, reall Esy, reall Esxy, reall xi, reall eta)
{
    reall xe = xi*eta;
    reall r  = ra  + rsx*xi  + rsy*eta  + rsxy *xe;
    reall mx = mxa + mxsx*xi + mxsy*eta + mxsxy*xe;
    reall my = mya + mysx*xi + mys *eta + mysxy*xe;
    reall E  = Ea  + Esx*xi  + Esy*eta  + Esxy *xe;
    if (r < (reall)RHOFLOOR) r = (reall)RHOFLOOR;
    return ((reall)GAMMA_V - 1.0) * (E - 0.5*(mx*mx + my*my)/r);
}

/* Zhang-Shu positivity limiter for the FULL bilinear Q1 state.
 * Bilinear ρ attains its cell min at a VERTEX (interior critical point is a
 * saddle; ρ is linear on each edge), so 4 corners suffice for density.  p is
 * concave in the conserved state but the state is bilinear, so p is NOT
 * guaranteed corner-minimised — check the 4 corners AND the 4 volume Gauss
 * points for pressure.  Two passes: (1) shrink ρ-slopes {∂x,∂y,ξη}, then
 * (2) shrink ALL slopes by one θ for pressure.  No-op in smooth flow. */
__device__ int g_pfix_count = 0;

__global__ void
limit_positivity(real * __restrict__ Qp, int N, int Np)
{
    const int g  = G_GHOST;
    const int N2 = N * N;
    const reall GPl = (reall)0.577350269189625764;
    /* point set: 4 corners then 4 Gauss (2x2) */
    const reall PX[8] = {-1.,  1., -1.,  1.,  -GPl,  GPl, -GPl,  GPl};
    const reall PY[8] = {-1., -1.,  1.,  1.,  -GPl, -GPl,  GPl,  GPl};
    for (int k = blockIdx.x * blockDim.x + threadIdx.x; k < N2;
             k += blockDim.x * gridDim.x) {
        int j = k / N, i = k % N;
        int jp = j + g, ip = i + g;

        reall ra   = fmax(Qp[IDX_P(Q_RHO,   jp, ip, Np)], (real)RHOFLOOR);
        reall rsx  = Qp[IDX_P(Q_RHO_SX, jp, ip, Np)];
        reall rsy  = Qp[IDX_P(Q_RHO_SY, jp, ip, Np)];
        reall rsxy = Qp[IDX_P(Q_RHO_SXY,jp, ip, Np)];
        reall mxa  = Qp[IDX_P(Q_MXA,  jp, ip, Np)];
        reall mxsx = Qp[IDX_P(Q_MXS,  jp, ip, Np)];
        reall mxsy = Qp[IDX_P(Q_MXSY, jp, ip, Np)];
        reall mxsxy= Qp[IDX_P(Q_MXSXY,jp, ip, Np)];
        reall mya  = Qp[IDX_P(Q_MYA,  jp, ip, Np)];
        reall mysx = Qp[IDX_P(Q_MYSX, jp, ip, Np)];
        reall mys  = Qp[IDX_P(Q_MYS,  jp, ip, Np)];
        reall mysxy= Qp[IDX_P(Q_MYSXY,jp, ip, Np)];
        reall Ea   = Qp[IDX_P(Q_E,    jp, ip, Np)];
        reall Esx  = Qp[IDX_P(Q_E_SX, jp, ip, Np)];
        reall Esy  = Qp[IDX_P(Q_E_SY, jp, ip, Np)];
        reall Esxy = Qp[IDX_P(Q_E_SXY,jp, ip, Np)];
        Qp[IDX_P(Q_RHO, jp, ip, Np)] = (real)ra;   /* commit density floor */

        /* ── 1. DENSITY: scale ρ-slopes so ρ(corner) ≥ RHOFLOOR ─────────── */
        reall rho_min = ra;
        for (int c = 0; c < 4; c++) {
            reall rc = ra + rsx*PX[c] + rsy*PY[c] + rsxy*PX[c]*PY[c];
            if (rc < rho_min) rho_min = rc;
        }
        reall rfloor = (reall)RHOFLOOR;
        if (rho_min < rfloor) {
            reall t1 = (ra - rfloor) / (ra - rho_min);   /* ra>rfloor guaranteed */
            if (t1 < 0.0) t1 = 0.0;
            rsx *= t1; rsy *= t1; rsxy *= t1;
            Qp[IDX_P(Q_RHO_SX, jp, ip, Np)] = (real)rsx;
            Qp[IDX_P(Q_RHO_SY, jp, ip, Np)] = (real)rsy;
            Qp[IDX_P(Q_RHO_SXY,jp, ip, Np)] = (real)rsxy;
        }

        /* ── 2. cell-average pressure floor: raise E if p̄ < PFLOOR ──────── */
        reall pbar = ((reall)GAMMA_V - 1.0) * (Ea - 0.5*(mxa*mxa + mya*mya)/ra);
        if (pbar < PFLOOR) {
            atomicAdd(&g_pfix_count, 1);
            Ea = 0.5*(mxa*mxa + mya*mya)/ra + PFLOOR / ((reall)GAMMA_V - 1.0);
            Qp[IDX_P(Q_E, jp, ip, Np)] = (real)Ea;
        }

        /* ── 3. PRESSURE: scale ALL slopes by θ so p(point) ≥ PFLOOR at the
         *   4 corners + 4 Gauss pts.  p(mean)≥PFLOOR and p concave along
         *   mean→point ⇒ bisect for the largest admissible t; θ = min. ──── */
        reall theta = 1.0;
        for (int c = 0; c < 8; c++) {
            reall pc = corner_pressure(ra,rsx,rsy,rsxy, mxa,mxsx,mxsy,mxsxy,
                                       mya,mysx,mys,mysxy, Ea,Esx,Esy,Esxy,
                                       PX[c], PY[c]);
            if (pc >= PFLOOR) continue;
            reall lo = 0.0, hi = 1.0;
            for (int it = 0; it < 40; it++) {
                reall t = 0.5*(lo + hi);
                reall pt = corner_pressure(ra, t*rsx, t*rsy, t*rsxy,
                                           mxa, t*mxsx, t*mxsy, t*mxsxy,
                                           mya, t*mysx, t*mys, t*mysxy,
                                           Ea, t*Esx, t*Esy, t*Esxy,
                                           PX[c], PY[c]);
                if (pt >= PFLOOR) lo = t; else hi = t;
            }
            if (lo < theta) theta = lo;
        }
        if (theta < 1.0) {
            Qp[IDX_P(Q_RHO_SX, jp, ip, Np)] = (real)(rsx  *theta);
            Qp[IDX_P(Q_RHO_SY, jp, ip, Np)] = (real)(rsy  *theta);
            Qp[IDX_P(Q_RHO_SXY,jp, ip, Np)] = (real)(rsxy *theta);
            Qp[IDX_P(Q_MXS,  jp, ip, Np)]   = (real)(mxsx *theta);
            Qp[IDX_P(Q_MXSY, jp, ip, Np)]   = (real)(mxsy *theta);
            Qp[IDX_P(Q_MXSXY,jp, ip, Np)]   = (real)(mxsxy*theta);
            Qp[IDX_P(Q_MYSX, jp, ip, Np)]   = (real)(mysx *theta);
            Qp[IDX_P(Q_MYS,  jp, ip, Np)]   = (real)(mys  *theta);
            Qp[IDX_P(Q_MYSXY,jp, ip, Np)]   = (real)(mysxy*theta);
            Qp[IDX_P(Q_E_SX, jp, ip, Np)]   = (real)(Esx  *theta);
            Qp[IDX_P(Q_E_SY, jp, ip, Np)]   = (real)(Esy  *theta);
            Qp[IDX_P(Q_E_SXY,jp, ip, Np)]   = (real)(Esxy *theta);
        }
    }
}

#if SHOCK_LIMIT
#define LIMIT_POS(Qarr) limit_positivity<<<GS_NBLK, BLOCK1D>>>((Qarr), N, Np)
#else
#define LIMIT_POS(Qarr) ((void)0)
#endif

/* ═══════════════════════════════════════════════════════════════════════════
 * Cockburn-Shu TVB minmod limiter on the Q1 slopes (componentwise, conserved
 * variables).  For each field: the ξ-slope is compared against the forward/
 * backward neighbour-average differences; TVB: slopes with |s| ≤ M·h² are
 * left untouched (protects smooth extrema, M=0 → strict TVD minmod).  If
 * either slope of a field is limited, its ξη cross mode is zeroed.
 * Needs FRESH ghost averages → call after apply_bc; run BEFORE the Zhang-Shu
 * positivity limiter (slope reduction preserves its corner-positivity).
 * Enable with -DSHU_LIMITER [-DSHU_TVB_M=<M>]; composes with any RHS kernel.
 * ═══════════════════════════════════════════════════════════════════════════ */
#ifndef SHU_TVB_M
#define SHU_TVB_M 0.0
#endif

__device__ __forceinline__ real
shu_minmod3(real a, real b, real c)
{
    if (a > 0. && b > 0. && c > 0.) return fmin(a, fmin(b, c));
    if (a < 0. && b < 0. && c < 0.) return fmax(a, fmax(b, c));
    return 0.;
}

__global__ void
shu_limit_kernel(real * __restrict__ Qp, int N, int Np, real h)
{
    const int g  = G_GHOST;
    const int N2 = N * N;
    const int A_  [4] = { Q_RHO,     Q_MXA,   Q_MYA,   Q_E     };
    const int SX_ [4] = { Q_RHO_SX,  Q_MXS,   Q_MYSX,  Q_E_SX  };
    const int SY_ [4] = { Q_RHO_SY,  Q_MXSY,  Q_MYS,   Q_E_SY  };
    const int SXY_[4] = { Q_RHO_SXY, Q_MXSXY, Q_MYSXY, Q_E_SXY };
    const real Mh2 = (real)SHU_TVB_M * h * h;

    for (int k = blockIdx.x * blockDim.x + threadIdx.x; k < N2;
             k += blockDim.x * gridDim.x) {
        int j = k / N, i = k % N;
        int jp = j + g, ip = i + g;

        for (int v = 0; v < 4; v++) {
            real aC = Qp[IDX_P(A_[v], jp, ip,   Np)];
            real aL = Qp[IDX_P(A_[v], jp, ip-1, Np)];
            real aR = Qp[IDX_P(A_[v], jp, ip+1, Np)];
            real aB = Qp[IDX_P(A_[v], jp-1, ip, Np)];
            real aT = Qp[IDX_P(A_[v], jp+1, ip, Np)];
            real sx = Qp[IDX_P(SX_[v], jp, ip, Np)];
            real sy = Qp[IDX_P(SY_[v], jp, ip, Np)];

            real sxn = (fabs(sx) <= Mh2) ? sx : shu_minmod3(sx, aR - aC, aC - aL);
            real syn = (fabs(sy) <= Mh2) ? sy : shu_minmod3(sy, aT - aC, aC - aB);

            if (sxn != sx || syn != sy) {
                Qp[IDX_P(SX_ [v], jp, ip, Np)] = sxn;
                Qp[IDX_P(SY_ [v], jp, ip, Np)] = syn;
                Qp[IDX_P(SXY_[v], jp, ip, Np)] = 0.;
            }
        }
    }
}

#ifdef SHU_LIMITER
#define SHU_LIMIT(Qarr) do { \
    apply_bc<<<GS_NBLK, BLOCK1D>>>((Qarr), N, Np); \
    shu_limit_kernel<<<GS_NBLK, BLOCK1D>>>((Qarr), N, Np, h); \
} while (0)
#else
#define SHU_LIMIT(Qarr) ((void)0)
#endif

/* ═══════════════════════════════════════════════════════════════════════════
 * SSP-RK3 stage kernels — operate directly on raw DOFs (no prim conversion)
 * ═══════════════════════════════════════════════════════════════════════════ */
__global__ void
rk3_s1(real * __restrict__ U1p,
       const real * __restrict__ U0p,
       const real * __restrict__ L,
       real dt, int N, int Np)
{
    const int g   = G_GHOST;
    const int N2  = N * N;
    const int Np2 = Np * Np;
    for (int k = blockIdx.x * blockDim.x + threadIdx.x; k < N2;
             k += blockDim.x * gridDim.x) {
        int j = k / N, i = k % N;
        int idxp = (j + g) * Np + (i + g);
        for (int q = 0; q < NVAR; q++)
            U1p[q * Np2 + idxp] = U0p[q * Np2 + idxp] + dt * L[q * N2 + k];
        /* floor rho and E */
        U1p[Q_RHO * Np2 + idxp] = fmax(U1p[Q_RHO * Np2 + idxp], 1e-14);
        U1p[Q_E   * Np2 + idxp] = fmax(U1p[Q_E   * Np2 + idxp], 1e-14);
    }
}

__global__ void
rk3_s2(real * __restrict__ U2p,
       const real * __restrict__ U0p,
       const real * __restrict__ U1p,
       const real * __restrict__ L,
       real dt, int N, int Np)
{
    const int g   = G_GHOST;
    const int N2  = N * N;
    const int Np2 = Np * Np;
    for (int k = blockIdx.x * blockDim.x + threadIdx.x; k < N2;
             k += blockDim.x * gridDim.x) {
        int j = k / N, i = k % N;
        int idxp = (j + g) * Np + (i + g);
        for (int q = 0; q < NVAR; q++)
            U2p[q * Np2 + idxp] = 0.75 * U0p[q * Np2 + idxp]
                                 + 0.25 * (U1p[q * Np2 + idxp] + dt * L[q * N2 + k]);
        U2p[Q_RHO * Np2 + idxp] = fmax(U2p[Q_RHO * Np2 + idxp], 1e-14);
        U2p[Q_E   * Np2 + idxp] = fmax(U2p[Q_E   * Np2 + idxp], 1e-14);
    }
}

__global__ void
rk3_s3(real * __restrict__ Up,
       const real * __restrict__ U0p,
       const real * __restrict__ U2p,
       const real * __restrict__ L,
       real dt, int N, int Np)
{
    const int g   = G_GHOST;
    const int N2  = N * N;
    const int Np2 = Np * Np;
    for (int k = blockIdx.x * blockDim.x + threadIdx.x; k < N2;
             k += blockDim.x * gridDim.x) {
        int j = k / N, i = k % N;
        int idxp = (j + g) * Np + (i + g);
        for (int q = 0; q < NVAR; q++)
            Up[q * Np2 + idxp] = (1./3.) * U0p[q * Np2 + idxp]
                                + (2./3.) * (U2p[q * Np2 + idxp] + dt * L[q * N2 + k]);
        Up[Q_RHO * Np2 + idxp] = fmax(Up[Q_RHO * Np2 + idxp], 1e-14);
        Up[Q_E   * Np2 + idxp] = fmax(Up[Q_E   * Np2 + idxp], 1e-14);
    }
}

/* ═══════════════════════════════════════════════════════════════════════════
 * max reduction
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

static real
gpu_max(const real *d_in, real *d_tmp, int n)
{
    reduce_max<<<GS_NBLK, BLOCK1D, BLOCK1D * sizeof(real)>>>(d_in, d_tmp, n);
    reduce_max<<<1,        BLOCK1D, BLOCK1D * sizeof(real)>>>(d_tmp, d_tmp, GS_NBLK);
    real v;   /* matches sizeof(d_tmp[0]) — promoted to real on return */
    CK(cudaMemcpy(&v, d_tmp, sizeof(real), cudaMemcpyDeviceToHost));
    return v;
}

/* ── min reduction (mirrors reduce_max) ──────────────────────────────────── */
__global__ void
reduce_min(const real * __restrict__ in, real * __restrict__ out, int n)
{
    extern __shared__ real sm[];
    int tid = threadIdx.x;
    real v = 1e38;
    for (int k = blockIdx.x * blockDim.x + tid; k < n; k += blockDim.x * gridDim.x)
        v = fmin(v, in[k]);
    sm[tid] = v; __syncthreads();
    for (int s = BLOCK1D / 2; s > 0; s >>= 1) {
        if (tid < s) sm[tid] = fmin(sm[tid], sm[tid + s]);
        __syncthreads();
    }
    if (tid == 0) out[blockIdx.x] = sm[0];
}

static real
gpu_min(const real *d_in, real *d_tmp, int n)
{
    reduce_min<<<GS_NBLK, BLOCK1D, BLOCK1D * sizeof(real)>>>(d_in, d_tmp, n);
    reduce_min<<<1,        BLOCK1D, BLOCK1D * sizeof(real)>>>(d_tmp, d_tmp, GS_NBLK);
    real v;   /* matches sizeof(d_tmp[0]) — promoted to real on return */
    CK(cudaMemcpy(&v, d_tmp, sizeof(real), cudaMemcpyDeviceToHost));
    return v;
}

/* ═══════════════════════════════════════════════════════════════════════════
 * Download helper — extract density from padded array for output
 * ═══════════════════════════════════════════════════════════════════════════ */
__global__ void
extract_rho_p(const real * __restrict__ Qp,
              real * __restrict__ out_rho,
              real * __restrict__ out_p,
              int N, int Np)
{
    const int g  = G_GHOST;
    const int N2 = N * N;
    for (int k = blockIdx.x * blockDim.x + threadIdx.x; k < N2;
             k += blockDim.x * gridDim.x) {
        int j = k / N, i = k % N;
        int jp = j + g, ip = i + g;
        real rho = Qp[IDX_P(Q_RHO, jp, ip, Np)];
        real mxa = Qp[IDX_P(Q_MXA, jp, ip, Np)];
        real mya = Qp[IDX_P(Q_MYA, jp, ip, Np)];
        real E   = Qp[IDX_P(Q_E,   jp, ip, Np)];
        rho = fmax(rho, 1e-14);
        out_rho[k] = rho;
        out_p  [k] = pressure_from_rt0(rho, mxa, mya, E);
    }
}

/* Extract density (left panel) and cell-centre v-velocity (right panel)    *
 * for Kelvin-Helmholtz visualisation.                                       *
 *   mya = cell-average y-momentum  →  v_avg = mya/rho                      */
__global__ void
extract_rho_v(const real * __restrict__ Qp,
              real * __restrict__ out_rho,
              real * __restrict__ out_v,
              int N, int Np)
{
    const int g  = G_GHOST;
    const int N2 = N * N;
    for (int k = blockIdx.x * blockDim.x + threadIdx.x; k < N2;
             k += blockDim.x * gridDim.x) {
        int j = k / N, i = k % N;
        int jp = j + g, ip = i + g;
        real rho = fmax(Qp[IDX_P(Q_RHO, jp, ip, Np)], 1e-14);
        real mya = Qp[IDX_P(Q_MYA, jp, ip, Np)];
        out_rho[k] = rho;
        out_v  [k] = mya / rho;  /* cell-average v */
    }
}

/* Extract cell-centre v-velocity (left panel) and vorticity (right panel)    *
 * for the doubly periodic shear layer.                                       *
 *   ω = ∂v/∂x − ∂u/∂y                                                        *
 * The cross-derivatives are NOT carried by the RT0 DOFs (which store only    *
 * ∂(ρu)/∂x and ∂(ρv)/∂y), so they are formed by central differences of the  *
 * cell-average velocities over the periodic ghost cells (filled by apply_bc).*/
__global__ void
extract_vort(const real * __restrict__ Qp,
             real * __restrict__ out_v,
             real * __restrict__ out_w,
             int N, int Np, real h)
{
    const int g  = G_GHOST;
    const int N2 = N * N;
    for (int k = blockIdx.x * blockDim.x + threadIdx.x; k < N2;
             k += blockDim.x * gridDim.x) {
        int j = k / N, i = k % N;
        int jp = j + g, ip = i + g;

        real rc = fmax(Qp[IDX_P(Q_RHO, jp,   ip,   Np)], 1e-14);
        real rR = fmax(Qp[IDX_P(Q_RHO, jp,   ip+1, Np)], 1e-14);
        real rL = fmax(Qp[IDX_P(Q_RHO, jp,   ip-1, Np)], 1e-14);
        real rT = fmax(Qp[IDX_P(Q_RHO, jp+1, ip,   Np)], 1e-14);
        real rB = fmax(Qp[IDX_P(Q_RHO, jp-1, ip,   Np)], 1e-14);

        real vc = Qp[IDX_P(Q_MYA, jp,   ip,   Np)] / rc;
        real vR = Qp[IDX_P(Q_MYA, jp,   ip+1, Np)] / rR;
        real vL = Qp[IDX_P(Q_MYA, jp,   ip-1, Np)] / rL;
        real uT = Qp[IDX_P(Q_MXA, jp+1, ip,   Np)] / rT;
        real uB = Qp[IDX_P(Q_MXA, jp-1, ip,   Np)] / rB;

        real dvdx = (vR - vL) / (2. * h);
        real dudy = (uT - uB) / (2. * h);

        out_v[k] = vc;
        out_w[k] = dvdx - dudy;
    }
}

/* Extract density (left panel) and vorticity (right panel) for the circular   *
 * Sod shock tube.                                                             *
 *   ω = ∂v/∂x − ∂u/∂y                                                         *
 * As in extract_vort, the cross-derivatives are not carried by the RT0 DOFs,  *
 * so they are formed by central differences of the cell-average velocities    *
 * over the ghost cells (filled by apply_bc).                                  */
__global__ void
extract_rho_vort(const real * __restrict__ Qp,
                 real * __restrict__ out_rho,
                 real * __restrict__ out_w,
                 int N, int Np, real h)
{
    const int g  = G_GHOST;
    const int N2 = N * N;
    for (int k = blockIdx.x * blockDim.x + threadIdx.x; k < N2;
             k += blockDim.x * gridDim.x) {
        int j = k / N, i = k % N;
        int jp = j + g, ip = i + g;

        real rc = fmax(Qp[IDX_P(Q_RHO, jp,   ip,   Np)], 1e-14);
        real rR = fmax(Qp[IDX_P(Q_RHO, jp,   ip+1, Np)], 1e-14);
        real rL = fmax(Qp[IDX_P(Q_RHO, jp,   ip-1, Np)], 1e-14);
        real rT = fmax(Qp[IDX_P(Q_RHO, jp+1, ip,   Np)], 1e-14);
        real rB = fmax(Qp[IDX_P(Q_RHO, jp-1, ip,   Np)], 1e-14);

        real vR = Qp[IDX_P(Q_MYA, jp,   ip+1, Np)] / rR;
        real vL = Qp[IDX_P(Q_MYA, jp,   ip-1, Np)] / rL;
        real uT = Qp[IDX_P(Q_MXA, jp+1, ip,   Np)] / rT;
        real uB = Qp[IDX_P(Q_MXA, jp-1, ip,   Np)] / rB;

        real dvdx = (vR - vL) / (2. * h);
        real dudy = (uT - uB) / (2. * h);

        out_rho[k] = rc;
        out_w[k]   = dvdx - dudy;
    }
}

/* ═══════════════════════════════════════════════════════════════════════════
 * Isentropic stationary vortex — exact solution
 *
 * Background: ρ∞=1, p∞=1/γ  →  c∞=1, T∞=1/γ
 * Vortex strength β = 2π·Ma  so peak |δu| = Ma·c∞ = Ma
 * Core radius R = 0.2 (in domain units [0,1])
 *
 * f(r) = β/(2π) · exp((1-r²)/2),   r² = [(x-½)²+(y-½)²]/R²
 * δu = -f·(y-½)/R,   δv = +f·(x-½)/R
 * δT = -(γ-1)β²/(8γπ²)·exp(1-r²)
 * T   = T∞ + δT
 * ρ   = ρ∞·(T/T∞)^{1/(γ-1)},   p = p∞·(T/T∞)^{γ/(γ-1)}
 * ═══════════════════════════════════════════════════════════════════════════ */
#define VORTEX_R  0.2

__device__ void
vortex_exact(real x, real y, real Ma,
             real *rho_out, real *u_out, real *v_out, real *p_out)
{
    const real g     = GAMMA_V;
    const real rhoI  = 1.;
    const real pI    = 1. / g;            /* c_inf = 1 */
    const real TI    = pI / rhoI;          /* = 1/g     */
    const real beta  = 2. * (real)M_PI * Ma;
    const real R     = VORTEX_R;

    real dx = x - 0.5, dy = y - 0.5;
    real r2 = (dx*dx + dy*dy) / (R*R);
    real f  = (beta / (2.*(real)M_PI)) * exp(0.5*(1. - r2));

    real u  = -f * dy / R;
    real v  = +f * dx / R;

    real dT  = -(g-1.) * beta*beta / (8.*g*(real)M_PI*(real)M_PI)
                * exp(1. - r2);
    real T   = TI + dT;
    if (T < 1e-6) T = 1e-6;
    real Tr  = T / TI;                     /* T / T_inf */
    real rho = rhoI * pow(Tr, 1./(g-1.));
    real p   = pI   * pow(Tr, g  /(g-1.));

    *rho_out = rho;
    *u_out   = u;
    *v_out   = v;
    *p_out   = p;
}

/* IC kernel — isentropic vortex projected onto RT0/P0 DOFs */
__global__ void
ic_vortex(real * __restrict__ Qp, int N, int Np, real h, real Ma)
{
    const int g  = G_GHOST;
    const int N2 = N * N;
    for (int k = blockIdx.x * blockDim.x + threadIdx.x; k < N2;
             k += blockDim.x * gridDim.x) {
        int j = k / N, i = k % N;
        int jp = j + g, ip = i + g;

        /* cell centre */
        real xc = (i + 0.5) * h, yc = (j + 0.5) * h;
        real rho, uc, vc, p;
        vortex_exact(xc, yc, Ma, &rho, &uc, &vc, &p);
        real E = p / (GAMMA_V - 1.) + 0.5 * rho * (uc*uc + vc*vc);

        /* RT0 DOFs: evaluate exact ρu at each face midpoint */
        /* Right face midpoint: x=(i+1)*h, y=yc */
        real rR, uR, vR_d, pR;
        vortex_exact((i+1)*h, yc, Ma, &rR, &uR, &vR_d, &pR);
        /* Left  face midpoint: x=i*h, y=yc */
        real rL, uL, vL_d, pL;
        vortex_exact(i*h,     yc, Ma, &rL, &uL, &vL_d, &pL);
        /* Top   face midpoint: x=xc, y=(j+1)*h */
        real rT, uT, vT, pT;
        vortex_exact(xc, (j+1)*h, Ma, &rT, &uT, &vT, &pT);
        /* Bottom face midpoint: x=xc, y=j*h */
        real rB, uB, vB, pB;
        vortex_exact(xc, j*h,     Ma, &rB, &uB, &vB, &pB);

        /* Full-P1 projection from face-midpoint momenta.
         * Normal:  q_x at right/left, q_y at top/bottom.
         * Tangential: q_x at top/bottom (→mxsy), q_y at right/left (→mysx). */
        real qxR = rR * uR,  qxL = rL * uL;     /* ρu at right/left  */
        real qxT = rT * uT,  qxB = rB * uB;     /* ρu at top/bottom  */
        real qyT = rT * vT,  qyB = rB * vB;     /* ρv at top/bottom  */
        real qyR = rR * vR_d, qyL = rL * vL_d;  /* ρv at right/left  */

        Qp[IDX_P(Q_RHO,  jp, ip, Np)] = rho;
        Qp[IDX_P(Q_MXA,  jp, ip, Np)] = 0.5 * (qxR + qxL);
        Qp[IDX_P(Q_MXS,  jp, ip, Np)] = 0.5 * (qxR - qxL);
        Qp[IDX_P(Q_MXSY, jp, ip, Np)] = 0.5 * (qxT - qxB);
        Qp[IDX_P(Q_MYA,  jp, ip, Np)] = 0.5 * (qyT + qyB);
        Qp[IDX_P(Q_MYSX, jp, ip, Np)] = 0.5 * (qyR - qyL);
        Qp[IDX_P(Q_MYS,  jp, ip, Np)] = 0.5 * (qyT - qyB);
        Qp[IDX_P(Q_E,    jp, ip, Np)] = E;
        (void)pL; (void)pT; (void)pB;
    }
}

/* ── L2 error vs. exact stationary vortex ──────────────────────────────── */
/* Writes per-cell squared error  (rho-rho_ex)²  and  (p-p_ex)²            */
__global__ void
compute_l2_err(const real * __restrict__ Qp,
               real * __restrict__ err_rho,
               real * __restrict__ err_p,
               int N, int Np, real h, real Ma)
{
    const int g  = G_GHOST;
    const int N2 = N * N;
    for (int k = blockIdx.x * blockDim.x + threadIdx.x; k < N2;
             k += blockDim.x * gridDim.x) {
        int j = k / N, i = k % N;
        int jp = j + g, ip = i + g;

        real rho = fmax(Qp[IDX_P(Q_RHO, jp, ip, Np)], 1e-14);
        real mxa = Qp[IDX_P(Q_MXA, jp, ip, Np)];
        real mya = Qp[IDX_P(Q_MYA, jp, ip, Np)];
        real E   = Qp[IDX_P(Q_E,   jp, ip, Np)];
        real p   = pressure_from_rt0(rho, mxa, mya, E);

        real xc = (i + 0.5) * h, yc = (j + 0.5) * h;
        real re, ue, ve, pe;
        vortex_exact(xc, yc, Ma, &re, &ue, &ve, &pe);

        real dr = rho - re, dp = p - pe;
        err_rho[k] = dr * dr * h * h;   /* volume-weighted */
        err_p  [k] = dp * dp * h * h;
    }
}

/* ── sum reduction (single precision) ──────────────────────────────────── */
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

static real
gpu_sum(const real *d_in, real *d_tmp, int n)
{
    reduce_sum<<<GS_NBLK, BLOCK1D, BLOCK1D * sizeof(real)>>>(d_in, d_tmp, n);
    reduce_sum<<<1,        BLOCK1D, BLOCK1D * sizeof(real)>>>(d_tmp, d_tmp, GS_NBLK);
    real v;
    CK(cudaMemcpy(&v, d_tmp, sizeof(real), cudaMemcpyDeviceToHost));
    return v;
}

/* ═══════════════════════════════════════════════════════════════════════════
 * Paper isentropic vortex  (Barsukow et al. 2025, Sec 6.2.1)
 *
 * Domain  : [0,10] x [0,10],  periodic
 * Centre  : (5,5),   strength ε = 5
 * Background: ρ∞=1, p∞=1   →  T∞=1,  c∞=sqrt(γ)
 *
 * δu = -ε/(2π) exp((1-r²)/2) (y-5)   (sign matches paper eq.)
 * δv = +ε/(2π) exp((1-r²)/2) (x-5)
 * δT = -(γ-1)ε²/(8γπ²) exp(1-r²)
 * T  = 1 + δT
 * ρ  = T^{1/(γ-1)},   p = T^{γ/(γ-1)}
 *
 * u0, v0 : optional background advection velocity
 * ═══════════════════════════════════════════════════════════════════════════ */

__device__ void
vortex_paper_exact(real x, real y, real u0, real v0,
                   real *rho_out, real *u_out, real *v_out, real *p_out)
{
    const real g  = GAMMA_V;
    const real eps = 5.;
    const real pi  = (real)M_PI;

    real dx = x - 5., dy = y - 5.;
    real r2 = dx*dx + dy*dy;
    real f  = eps / (2.*pi) * exp(0.5*(1. - r2));

    real u  = u0 - f * dy;
    real v  = v0 + f * dx;

    real dT = -(g-1.)*eps*eps / (8.*g*pi*pi) * exp(1. - r2);
    real T  = 1. + dT;
    if (T < 1e-6) T = 1e-6;

    *rho_out = pow(T, 1./(g-1.));
    *p_out   = pow(T, g  /(g-1.));
    *u_out   = u;
    *v_out   = v;
}

/* IC kernel — paper isentropic vortex on [0,10]^2, h = 10/N */
__global__ void
ic_vortex_paper(real * __restrict__ Qp, int N, int Np,
                real h, real u0, real v0)
{
    const int g  = G_GHOST;
    const int N2 = N * N;
    for (int k = blockIdx.x * blockDim.x + threadIdx.x; k < N2;
             k += blockDim.x * gridDim.x) {
        int j = k / N, i = k % N;
        int jp = j + g, ip = i + g;

        real xc = (i + 0.5) * h, yc = (j + 0.5) * h;

        /* ── L2 projection of ALL 16 Q1 modes by 2×2 Gauss quadrature ─────
         * Modal basis {1,ξ,η,ξη}:  avg = (1/4)Σq,  sx = (3/4)Σξ_g·q,
         * sy = (3/4)Ση_g·q,  sxy = (9/4)Σξ_gη_g·q .
         * Gauss cell-average error is O(h⁴) — required for the evolved
         * average to display any superconvergence past 2nd order (midpoint
         * sampling is only O(h²) and caps every measured rate at 2).      */
        const real GPq = (real)0.577350269189625764;
        real m[4][4] = {{0.,0.,0.,0.},{0.,0.,0.,0.},{0.,0.,0.,0.},{0.,0.,0.,0.}};
        for (int gq = 0; gq < 4; gq++) {
            real xi  = (gq & 1) ? GPq : -GPq;
            real eta = (gq & 2) ? GPq : -GPq;
            real re, ue, ve, pe;
            vortex_paper_exact(xc + (real)0.5*h*xi, yc + (real)0.5*h*eta,
                               u0, v0, &re, &ue, &ve, &pe);
            real q[4] = { re, re*ue, re*ve,
                          pe/(GAMMA_V-1.) + (real)0.5*re*(ue*ue + ve*ve) };
            for (int v = 0; v < 4; v++) {
                m[v][0] += (real)0.25 * q[v];
                m[v][1] += (real)0.75 * xi        * q[v];
                m[v][2] += (real)0.75 * eta       * q[v];
                m[v][3] += (real)2.25 * (xi*eta)  * q[v];
            }
        }
        const int A_  [4] = { Q_RHO,     Q_MXA,   Q_MYA,   Q_E     };
        const int SX_ [4] = { Q_RHO_SX,  Q_MXS,   Q_MYSX,  Q_E_SX  };
        const int SY_ [4] = { Q_RHO_SY,  Q_MXSY,  Q_MYS,   Q_E_SY  };
        const int SXY_[4] = { Q_RHO_SXY, Q_MXSXY, Q_MYSXY, Q_E_SXY };
        for (int v = 0; v < 4; v++) {
            Qp[IDX_P(A_  [v], jp, ip, Np)] = m[v][0];
            Qp[IDX_P(SX_ [v], jp, ip, Np)] = m[v][1];
            Qp[IDX_P(SY_ [v], jp, ip, Np)] = m[v][2];
            Qp[IDX_P(SXY_[v], jp, ip, Np)] = m[v][3];
        }
    }
}

/* L2 error vs. exact paper vortex
 * Computes errors in rho, rho*u, rho*v, rho*E  (matches Tables 2 & 3 of paper).
 * For the stationary case (u0=v0=0, tf=1) the exact solution equals the IC.
 * For the moving case (u0=v0=1, tf=10) the vortex traverses exactly one period
 * on the [0,10]^2 periodic domain, so the exact solution also equals the IC.
 */
__global__ void
compute_l2_err_paper(const real * __restrict__ Qp,
                     real * __restrict__ err_rho,
                     real * __restrict__ err_ru,
                     real * __restrict__ err_rv,
                     real * __restrict__ err_re,
                     int N, int Np, real h, real u0, real v0)
{
    const int g  = G_GHOST;
    const int N2 = N * N;
    for (int k = blockIdx.x * blockDim.x + threadIdx.x; k < N2;
             k += blockDim.x * gridDim.x) {
        int j = k / N, i = k % N;
        int jp = j + g, ip = i + g;

        real rho = fmax(Qp[IDX_P(Q_RHO, jp, ip, Np)], 1e-14);
        real mxa = Qp[IDX_P(Q_MXA, jp, ip, Np)];
        real mya = Qp[IDX_P(Q_MYA, jp, ip, Np)];
        real E   = Qp[IDX_P(Q_E,   jp, ip, Np)];
        /* cell-average momenta = mxa, mya directly (modal mode 0) */
        real ru  = mxa;
        real rv  = mya;
        /* total energy density stored as E = rho*E_spec */
        real rE  = E;

        real xc = (i + 0.5) * h, yc = (j + 0.5) * h;
        /* exact CELL AVERAGES by 2×2 Gauss (O(h⁴)) — midpoint values are
         * only O(h²)-accurate averages and would cap the measured rate at 2 */
        const real GPq = (real)0.577350269189625764;
        real re = 0., rue = 0., rve = 0., rEe = 0.;
        for (int gq = 0; gq < 4; gq++) {
            real xi  = (gq & 1) ? GPq : -GPq;
            real eta = (gq & 2) ? GPq : -GPq;
            real rg, ug, vg, pg;
            vortex_paper_exact(xc + (real)0.5*h*xi, yc + (real)0.5*h*eta,
                               u0, v0, &rg, &ug, &vg, &pg);
            re  += (real)0.25 * rg;
            rue += (real)0.25 * rg*ug;
            rve += (real)0.25 * rg*vg;
            rEe += (real)0.25 * (pg/(GAMMA_V-1.) + (real)0.5*rg*(ug*ug + vg*vg));
        }

        real vol = h * h;
        err_rho[k] = (rho - re)  * (rho - re)  * vol;
        err_ru [k] = (ru  - rue) * (ru  - rue) * vol;
        err_rv [k] = (rv  - rve) * (rv  - rve) * vol;
        err_re [k] = (rE  - rEe) * (rE  - rEe) * vol;
    }
}

/* ═══════════════════════════════════════════════════════════════════════════
 * Kelvin-Helmholtz instability  (Barsukow et al. 2025, Sec 6.2.5)
 *
 * Reference domain [0,2] x [-0.5,0.5],  periodic everywhere.
 * Implemented here on [0,2] x [0,2]  with y_shifted = y_phys - 1
 * so shear layers sit at y_phys = 0.75 and 1.25  (y_shifted = ±0.25).
 * h = 2/N  →  square cells. Use ./rt_dg N 0 kh
 *
 * M = Ma (cli),  r = Ma*0.1,  δ = 0.1,  ω = 1/16
 * ρ = γ + H(y)*r,  u = M*H(y),  v = δ*M*sin(2πx),  p = 1
 * ═══════════════════════════════════════════════════════════════════════════ */

__device__ real
kh_H(real y)   /* y = shifted coordinate ∈ [-0.5, 0.5] */
{
    const real omega = 1. / 16.;
    const real y1 = -0.25, y2 = 0.25;
    const real pi  = (real)M_PI;

    /* Transition at y1: +1 → -1 */
    if (y >= y1 - omega*0.5 && y < y1 + omega*0.5)
        return -sin(pi * (y - y1) / omega);
    /* Central region */
    if (y >= y1 + omega*0.5 && y < y2 - omega*0.5)
        return -1.;
    /* Transition at y2: -1 → +1 */
    if (y >= y2 - omega*0.5 && y < y2 + omega*0.5)
        return  sin(pi * (y - y2) / omega);
    /* Outer region */
    return 1.;
}

__global__ void
ic_kh(real * __restrict__ Qp, int N, int Np, real h, real Mkh)
{
    const int g  = G_GHOST;
    const int N2 = N * N;
    const real rkh  = Mkh * 0.1;  /* density jump ∝ Ma */
    const real del  = 0.1;
    const real pi   = (real)M_PI;

    for (int k = blockIdx.x * blockDim.x + threadIdx.x; k < N2;
             k += blockDim.x * gridDim.x) {
        int j = k / N, i = k % N;
        int jp = j + g, ip = i + g;

        real xc = (i + 0.5) * h;
        real yc = (j + 0.5) * h;
        real ys = yc - 1.0;  /* centre domain at y=1 → y_shifted ∈ [-1,1] */

        /* Cell-center primitives */
        real Hc  = kh_H(ys);
        real rho = GAMMA_V + Hc * rkh;
        real u   = Mkh * Hc;
        real v   = del * Mkh * sin(2. * pi * xc);
        real p   = 1.;
        real E   = p / (GAMMA_V - 1.) + 0.5 * rho * (u*u + v*v);

        /* Full-P1 momentum projection from face-midpoint momenta.
         * x-mom ρu = ρ(y)·Mkh·H(y)  varies with y  → mxs=0, mxsy from y-faces
         * y-mom ρv = ρ(y)·δMa·sin(2πx) varies with x → mys from y-faces, mysx from x-faces */
        real ys_t = (j+1)*h - 1.0;
        real ys_b = j*h     - 1.0;
        real Ht   = kh_H(ys_t),  Hb = kh_H(ys_b);
        real rhoT = GAMMA_V + Ht * rkh;
        real rhoB = GAMMA_V + Hb * rkh;

        /* x-momentum at the four face midpoints */
        real qxR = rho * u,          qxL = rho * u;                 /* y=yc → cell value */
        real qxT = rhoT * (Mkh*Ht),  qxB = rhoB * (Mkh*Hb);         /* y-faces carry shear */
        /* y-momentum at the four face midpoints (ρ at right/left = cell ρ(yc)) */
        real vR  = del * Mkh * sin(2.*pi*((i+1)*h));
        real vL  = del * Mkh * sin(2.*pi*( i   *h));
        real qyR = rho * vR,         qyL = rho * vL;                /* x-faces carry sin(2πx) */
        real qyT = rhoT * v,         qyB = rhoB * v;                /* v(xc) same top/bottom */

        Qp[IDX_P(Q_RHO,  jp, ip, Np)] = rho;
        Qp[IDX_P(Q_MXA,  jp, ip, Np)] = 0.5 * (qxR + qxL);
        Qp[IDX_P(Q_MXS,  jp, ip, Np)] = 0.5 * (qxR - qxL);   /* = 0 */
        Qp[IDX_P(Q_MXSY, jp, ip, Np)] = 0.5 * (qxT - qxB);
        Qp[IDX_P(Q_MYA,  jp, ip, Np)] = 0.5 * (qyT + qyB);
        Qp[IDX_P(Q_MYSX, jp, ip, Np)] = 0.5 * (qyR - qyL);
        Qp[IDX_P(Q_MYS,  jp, ip, Np)] = 0.5 * (qyT - qyB);
        Qp[IDX_P(Q_E,    jp, ip, Np)] = E;
    }
}

/* ═══════════════════════════════════════════════════════════════════════════
 * Doubly periodic shear layer  (Bell, Colella & Glaz 1989)
 *
 * Domain  : [0,1]^2, periodic everywhere ("doubly periodic").
 * Two shear layers at y=0.25 and y=0.75 give a y-periodic profile that rolls
 * up into a pair of counter-rotating billows — the classic incompressible
 * benchmark, run here as a low-Mach compressible flow.
 *
 * Background  : ρ = γ, p = 1  →  c∞ = 1, so velocity amplitudes = Mach number.
 *   u = Ma·tanh((y-0.25)/δ_s)   y ≤ 0.5
 *       Ma·tanh((0.75-y)/δ_s)   y > 0.5
 *   v = Ma·δ·sin(2πx)
 * with shear thickness δ_s = 1/30 (thin layer) and perturbation δ = 0.05.
 *
 * Usage: ./rt_dg N dsl [Ma]    default Ma = 0.1
 * ═══════════════════════════════════════════════════════════════════════════ */
#define DSL_THICK 30.0    /* 1/δ_s : shear-layer sharpness */
#define DSL_PERT  0.05    /* δ     : transverse perturbation amplitude */

__device__ real
dsl_u(real y)   /* streamwise velocity profile in units of Ma, y ∈ [0,1] */
{
    return (y <= 0.5) ? tanh(DSL_THICK * (y - 0.25))
                      : tanh(DSL_THICK * (0.75 - y));
}

__global__ void
ic_dsl(real * __restrict__ Qp, int N, int Np, real h, real Mdsl)
{
    const int g  = G_GHOST;
    const int N2 = N * N;
    const real pi = (real)M_PI;

    for (int k = blockIdx.x * blockDim.x + threadIdx.x; k < N2;
             k += blockDim.x * gridDim.x) {
        int j = k / N, i = k % N;
        int jp = j + g, ip = i + g;

        real xc = (i + 0.5) * h;
        real yc = (j + 0.5) * h;

        real rho = GAMMA_V;
        real u   = Mdsl * dsl_u(yc);
        real v   = Mdsl * DSL_PERT * sin(2. * pi * xc);
        real p   = 1.;
        real E   = p / (GAMMA_V - 1.) + 0.5 * rho * (u*u + v*v);

        /* Full-P1 momentum projection from face-midpoint momenta.
         * ρu = ρ·Ma·dsl_u(y)  varies with y  → mxs=0, mxsy from y-faces
         * ρv = ρ·Ma·δ·sin(2πx) varies with x → mys=0, mysx from x-faces */
        real qxR = rho * u,  qxL = rho * u;                            /* y=yc → cell value */
        real qxT = rho * (Mdsl * dsl_u((j+1)*h));                      /* ρu at top    */
        real qxB = rho * (Mdsl * dsl_u( j   *h));                      /* ρu at bottom */
        real qyT = rho * v,  qyB = rho * v;                            /* x=xc → cell value */
        real qyR = rho * (Mdsl * DSL_PERT * sin(2.*pi*((i+1)*h)));     /* ρv at right */
        real qyL = rho * (Mdsl * DSL_PERT * sin(2.*pi*( i   *h)));     /* ρv at left  */

        Qp[IDX_P(Q_RHO,  jp, ip, Np)] = rho;
        Qp[IDX_P(Q_MXA,  jp, ip, Np)] = 0.5 * (qxR + qxL);   /* = rho*u */
        Qp[IDX_P(Q_MXS,  jp, ip, Np)] = 0.5 * (qxR - qxL);   /* = 0     */
        Qp[IDX_P(Q_MXSY, jp, ip, Np)] = 0.5 * (qxT - qxB);
        Qp[IDX_P(Q_MYA,  jp, ip, Np)] = 0.5 * (qyT + qyB);   /* = rho*v */
        Qp[IDX_P(Q_MYSX, jp, ip, Np)] = 0.5 * (qyR - qyL);
        Qp[IDX_P(Q_MYS,  jp, ip, Np)] = 0.5 * (qyT - qyB);   /* = 0     */
        Qp[IDX_P(Q_E,    jp, ip, Np)] = E;
    }
}

/* ═══════════════════════════════════════════════════════════════════════════
 * Low-Mach vortex  (Barsukow et al. 2025, Sec 2 / "behaviour of DG")
 *
 * Domain  : [0,1]^2, periodic, centre (0.5, 0.5)
 * ρ = 1 everywhere.
 *
 * Tangential velocity (r = distance from centre):
 *   v_φ = 5r          r < 0.2
 *         2 - 5r  0.2 ≤ r < 0.4
 *         0           else
 *
 * Pressure (radial equilibrium, dp/dr = ρ v_φ²/r):
 *   p = p₀ + 12.5 r²                              r < 0.2
 *       p₀ + 4 ln(5r) + 4 - 20r + 12.5 r²   0.2 ≤ r < 0.4
 *       p₀ + 4 ln 2 - 2                           else
 *
 * p₀ = 1/(γ ε²) - 1/2  sets max local Ma = ε.
 *
 * Usage: ./rt_dg N lmv [eps]    default eps = 0.1
 * ═══════════════════════════════════════════════════════════════════════════ */

__device__ real
lmv_vphi(real r)
{
    if (r < 0.2) return 5. * r;
    if (r < 0.4) return 2. - 5. * r;
    return 0.;
}

__device__ real
lmv_pressure(real r, real p0)
{
    if (r < 0.2) return p0 + 12.5 * r * r;
    if (r < 0.4) return p0 + 4. * log(5. * r) + 4. - 20. * r + 12.5 * r * r;
    return p0 + 4. * log(2.) - 2.;
}

/* Helper: primitives at point (x,y) for the low-Mach vortex */
__device__ void
lmv_exact(real x, real y, real p0,
          real *rho_out, real *u_out, real *v_out, real *p_out)
{
    const real xvc = 0.5, yvc = 0.5;
    real dx = x - xvc, dy = y - yvc;
    real r  = sqrt(dx*dx + dy*dy);
    real vp = lmv_vphi(r);
    real ir = (r > 1e-30) ? 1. / r : 0.;
    *rho_out = 1.;
    *u_out   = -vp * dy * ir;
    *v_out   =  vp * dx * ir;
    *p_out   = lmv_pressure(r, p0);
}

__global__ void
ic_lmv(real * __restrict__ Qp, int N, int Np, real h, real eps)
{
    const real p0   = 1. / (GAMMA_V * eps * eps) - 0.5;
    const int g  = G_GHOST;
    const int N2 = N * N;
    for (int k = blockIdx.x * blockDim.x + threadIdx.x; k < N2;
             k += blockDim.x * gridDim.x) {
        int j = k / N, i = k % N;
        int jp = j + g, ip = i + g;

        real xc = (i + 0.5) * h, yc = (j + 0.5) * h;
        real rho, uc, vc, p;
        lmv_exact(xc, yc, p0, &rho, &uc, &vc, &p);
        real E = p / (GAMMA_V - 1.) + 0.5 * rho * (uc*uc + vc*vc);

        /* Face-midpoint momenta for the full-P1 projection (Gresho: ρ=1). */
        real rR, uR, vR, pR;  lmv_exact((i+1)*h, yc,     p0, &rR, &uR, &vR, &pR);
        real rL, uL, vL, pL;  lmv_exact( i   *h, yc,     p0, &rL, &uL, &vL, &pL);
        real rT, uT, vT, pT;  lmv_exact(xc,     (j+1)*h, p0, &rT, &uT, &vT, &pT);
        real rB, uB, vB, pB;  lmv_exact(xc,      j   *h, p0, &rB, &uB, &vB, &pB);

        real qxR = rR * uR, qxL = rL * uL;   /* ρu at right/left */
        real qxT = rT * uT, qxB = rB * uB;   /* ρu at top/bottom */
        real qyT = rT * vT, qyB = rB * vB;   /* ρv at top/bottom */
        real qyR = rR * vR, qyL = rL * vL;   /* ρv at right/left */

        Qp[IDX_P(Q_RHO,  jp, ip, Np)] = rho;
        Qp[IDX_P(Q_MXA,  jp, ip, Np)] = 0.5 * (qxR + qxL);
        Qp[IDX_P(Q_MXS,  jp, ip, Np)] = 0.5 * (qxR - qxL);
        Qp[IDX_P(Q_MXSY, jp, ip, Np)] = 0.5 * (qxT - qxB);
        Qp[IDX_P(Q_MYA,  jp, ip, Np)] = 0.5 * (qyT + qyB);
        Qp[IDX_P(Q_MYSX, jp, ip, Np)] = 0.5 * (qyR - qyL);
        Qp[IDX_P(Q_MYS,  jp, ip, Np)] = 0.5 * (qyT - qyB);
        Qp[IDX_P(Q_E,    jp, ip, Np)] = E;
        (void)pL; (void)pT; (void)pB; (void)pR;
    }
}

/* L2 error vs. exact stationary low-Mach vortex (exact = IC).
 * Reports ρ, p, and VELOCITY error — the velocity field is the primary
 * quantity a low-Mach vortex scheme must preserve. */
__global__ void
compute_l2_err_lmv(const real * __restrict__ Qp,
                   real * __restrict__ err_rho,
                   real * __restrict__ err_p,
                   real * __restrict__ err_vel,
                   int N, int Np, real h, real eps)
{
    const real p0 = 1. / (GAMMA_V * eps * eps) - 0.5;
    const int g  = G_GHOST;
    const int N2 = N * N;
    for (int k = blockIdx.x * blockDim.x + threadIdx.x; k < N2;
             k += blockDim.x * gridDim.x) {
        int j = k / N, i = k % N;
        int jp = j + g, ip = i + g;
        real xc = (i + 0.5) * h, yc = (j + 0.5) * h;

        real rho_ex, u_ex, v_ex, p_ex;
        lmv_exact(xc, yc, p0, &rho_ex, &u_ex, &v_ex, &p_ex);

        real rho = fmax(Qp[IDX_P(Q_RHO, jp, ip, Np)], 1e-14);
        real mxa = Qp[IDX_P(Q_MXA, jp, ip, Np)];
        real mya = Qp[IDX_P(Q_MYA, jp, ip, Np)];
        real E   = Qp[IDX_P(Q_E,   jp, ip, Np)];
        real p   = pressure_from_rt0(rho, mxa, mya, E);

        real u = mxa / rho, v = mya / rho;   /* cell-average velocity */
        real du = u - u_ex,  dv = v - v_ex;
        real dr = rho - rho_ex, dp = p - p_ex;
        err_rho[k] = dr * dr * h * h;
        err_p  [k] = dp * dp * h * h;
        err_vel[k] = (du*du + dv*dv) * h * h;
    }
}

/* ═══════════════════════════════════════════════════════════════════════════
 * Vortex-acoustic wave interaction test  (Coiffier 2025, Sec 6.6.3)
 *
 * A compact-support stationary vortex superposed with a right-traveling
 * low-Mach acoustic wave, on [0,1]^2 with periodic BCs.
 *
 * Background state: ρ∞=1, p∞=1, c∞=sqrt(γ)
 *
 * Vortex (centred at (0.5,0.5), radius Rv=0.2, C∞ compact support):
 *   r̄ = r/Rv
 *   For r̄<1: bump(r̄) = r̄·exp(1/(r̄²-1)) / fmax
 *               uθ = Mref·c∞·bump(r̄)  [max = Mref·c∞]
 *   Isentropic: c² = c∞² - (gamma-1)/2·uθ²,  rho=rho∞·(c/c∞)^{2/(gamma-1)}
 *
 * Acoustic wave (right-traveling, centred at xw=0.20, half-width sw=0.05):
 *   x̄ = (x-xw)/sw,  bump_w = exp(1/(x̄²-1))  for |x̄|<1
 *   δu = Mref·c∞·bump_w,  δρ = ρ∞·Mref·bump_w,  δp = γ·p∞·Mref·bump_w
 *
 * Usage:  ./rt_dg N va [Mref]   default Mref=0.1
 * ═══════════════════════════════════════════════════════════════════════════ */

/* Compact-support isentropic vortex perturbation at point (x,y). */
__device__ void
vacwav_vortex(real x, real y, real Mref,
              real *drho, real *du, real *dv, real *dp)
{
    const real gm  = GAMMA_V;
    const real c0  = sqrt(gm);   /* c_infty = sqrt(gamma*p0/rho0), p0=rho0=1 */
    const real Rv  = 0.2;
    const real xvc = 0.5, yvc = 0.5;

    real dx  = x - xvc, dy = y - yvc;
    real r   = sqrt(dx*dx + dy*dy);
    real rbar = r / Rv;

    if (rbar >= 1.) { *drho = 0.; *du = 0.; *dv = 0.; *dp = 0.; return; }

    /* C-inf bump: f(s) = s * exp(1/(s^2-1)),  max at s* = (sqrt(6)-sqrt(2))/2 */
    const real s_star = 0.5 * (sqrt(6.) - sqrt(2.));
    const real fmax   = s_star * exp(1. / (s_star*s_star - 1.));
    real bump = rbar * exp(1. / (rbar*rbar - 1.));
    real uth  = Mref * c0 * bump / fmax;

    /* Isentropic balance: c^2 = c0^2 - (gamma-1)/2 * uth^2  (Bernoulli) */
    real c2    = c0*c0 - 0.5*(gm - 1.)*uth*uth;
    if (c2 < 1e-6*c0*c0) c2 = 1e-6*c0*c0;
    real ratio = c2 / (c0*c0);
    real rho_v = pow(ratio, 1./(gm - 1.));   /* rho_infty = 1 */
    real p_v   = pow(ratio, gm  /(gm - 1.));  /* p_infty   = 1 */

    /* Cartesian velocity from tangential */
    real ir = 1. / (r + 1e-30);
    *du   = -uth * dy * ir;
    *dv   =  uth * dx * ir;
    *drho = rho_v - 1.;
    *dp   = p_v   - 1.;
}

/* Right-traveling acoustic wave perturbation at x. */
__device__ void
vacwav_wave(real x, real Mref, real *drho, real *du, real *dp)
{
    const real gm  = GAMMA_V;
    const real c0  = sqrt(gm);
    const real xw  = 0.10;   /* wave centre */
    const real sw  = 0.05;   /* half-width  */

    real xbar = (x - xw) / sw;
    if (fabs(xbar) >= 1.) { *drho = 0.; *du = 0.; *dp = 0.; return; }

    real bump = exp(1. / (xbar*xbar - 1.));
    *du   = Mref * c0 * bump;           /* acoustic velocity perturbation     */
    *drho = Mref * bump;                /* drho/rho_infty = Mref * bump        */
    *dp   = gm * Mref * bump;           /* dp/p_infty = gamma * drho/rho_infty */
}

/* IC kernel: compact vortex + right-traveling acoustic wave */
__global__ void
ic_vacwav(real * __restrict__ Qp, int N, int Np, real h, real Mref)
{
    const int g  = G_GHOST;
    const int N2 = N * N;
    for (int k = blockIdx.x * blockDim.x + threadIdx.x; k < N2;
             k += blockDim.x * gridDim.x) {
        int j = k / N, i = k % N;
        int jp = j + g, ip = i + g;

        real xc = (i + 0.5) * h, yc = (j + 0.5) * h;

        /* --- Cell-centre primitives (for energy) --- */
        real drho_v, du_v, dv_v, dp_v;
        vacwav_vortex(xc, yc, Mref, &drho_v, &du_v, &dv_v, &dp_v);
        real drho_w, du_w, dp_w;
        vacwav_wave(xc, Mref, &drho_w, &du_w, &dp_w);

        real rho = fmax(1. + drho_v + drho_w, 1e-14);
        real uc  = du_v + du_w;
        real vc  = dv_v;
        real p   = fmax(1. + dp_v  + dp_w,  1e-14);
        real E   = p / (GAMMA_V - 1.) + 0.5 * rho * (uc*uc + vc*vc);

        /* --- RT0 face-midpoint normal fluxes --- */
        /* Right face (x = (i+1)*h, y = yc) */
        real drv_r, duv_r, dvv_r, dpv_r;
        vacwav_vortex((i+1)*h, yc, Mref, &drv_r, &duv_r, &dvv_r, &dpv_r);
        real drw_r, duw_r, dpw_r;
        vacwav_wave((i+1)*h, Mref, &drw_r, &duw_r, &dpw_r);
        real mr = fmax(1.+drv_r+drw_r, 1e-14) * (duv_r + duw_r);

        /* Left face (x = i*h, y = yc) */
        real drv_l, duv_l, dvv_l, dpv_l;
        vacwav_vortex(i*h, yc, Mref, &drv_l, &duv_l, &dvv_l, &dpv_l);
        real drw_l, duw_l, dpw_l;
        vacwav_wave(i*h, Mref, &drw_l, &duw_l, &dpw_l);
        real ml = -(fmax(1.+drv_l+drw_l, 1e-14) * (duv_l + duw_l));

        /* Top face (x = xc, y = (j+1)*h): wave has no v-component */
        real drv_t, duv_t, dvv_t, dpv_t;
        vacwav_vortex(xc, (j+1)*h, Mref, &drv_t, &duv_t, &dvv_t, &dpv_t);
        real mt = fmax(1.+drv_t, 1e-14) * dvv_t;

        /* Bottom face (x = xc, y = j*h) */
        real drv_b, duv_b, dvv_b, dpv_b;
        vacwav_vortex(xc, j*h, Mref, &drv_b, &duv_b, &dvv_b, &dpv_b);
        real mb = -(fmax(1.+drv_b, 1e-14) * dvv_b);

        /* Extra P1 slopes: tangential-momentum variation across the cell.
         * The acoustic wave depends on x only, so at the y-faces (x=xc) it
         * equals the cell-centre wave (drho_w, du_w). */
        real qyR = fmax(1.+drv_r+drw_r, 1e-14) * dvv_r;   /* ρv at right */
        real qyL = fmax(1.+drv_l+drw_l, 1e-14) * dvv_l;   /* ρv at left  */
        real qxT = fmax(1.+drv_t+drho_w, 1e-14) * (duv_t + du_w);  /* ρu at top    */
        real qxB = fmax(1.+drv_b+drho_w, 1e-14) * (duv_b + du_w);  /* ρu at bottom */

        Qp[IDX_P(Q_RHO,  jp, ip, Np)] = rho;
        Qp[IDX_P(Q_MXA,  jp, ip, Np)] = 0.5 * (mr - ml);
        Qp[IDX_P(Q_MXS,  jp, ip, Np)] = 0.5 * (mr + ml);
        Qp[IDX_P(Q_MXSY, jp, ip, Np)] = 0.5 * (qxT - qxB);
        Qp[IDX_P(Q_MYA,  jp, ip, Np)] = 0.5 * (mt - mb);
        Qp[IDX_P(Q_MYSX, jp, ip, Np)] = 0.5 * (qyR - qyL);
        Qp[IDX_P(Q_MYS,  jp, ip, Np)] = 0.5 * (mt + mb);
        Qp[IDX_P(Q_E,    jp, ip, Np)] = E;
    }
}

/* Extract cell-centre local Mach number M = |u|/c */
__global__ void
extract_mach(const real * __restrict__ Qp,
             real       * __restrict__ out_mach,
             int N, int Np)
{
    const int g  = G_GHOST;
    const int N2 = N * N;
    for (int k = blockIdx.x * blockDim.x + threadIdx.x; k < N2;
             k += blockDim.x * gridDim.x) {
        int j = k / N, i = k % N;
        int jp = j + g, ip = i + g;
        real rho = fmax(Qp[IDX_P(Q_RHO, jp, ip, Np)], 1e-14);
        real mxa = Qp[IDX_P(Q_MXA, jp, ip, Np)];
        real mya = Qp[IDX_P(Q_MYA, jp, ip, Np)];
        real E   = Qp[IDX_P(Q_E,   jp, ip, Np)];
        real p   = pressure_from_rt0(rho, mxa, mya, E);
        real u   = mxa / rho, v = mya / rho;
        real c   = sound_speed(rho, p);
        out_mach[k] = sqrt(u*u + v*v) / c;
    }
}

/* Extract cell-centre |velocity| (left) and pressure (right) */
__global__ void
extract_magv_p(const real * __restrict__ Qp,
               real       * __restrict__ out_magv,
               real       * __restrict__ out_p,
               int N, int Np)
{
    const int g  = G_GHOST;
    const int N2 = N * N;
    for (int k = blockIdx.x * blockDim.x + threadIdx.x; k < N2;
             k += blockDim.x * gridDim.x) {
        int j = k / N, i = k % N;
        int jp = j + g, ip = i + g;
        real rho = fmax(Qp[IDX_P(Q_RHO, jp, ip, Np)], 1e-14);
        real mxa = Qp[IDX_P(Q_MXA, jp, ip, Np)];
        real mya = Qp[IDX_P(Q_MYA, jp, ip, Np)];
        real E   = Qp[IDX_P(Q_E,   jp, ip, Np)];
        real u   = mxa / rho, v = mya / rho;
        out_magv[k] = sqrt(u*u + v*v);
        out_p[k]    = pressure_from_rt0(rho, mxa, mya, E);
    }
}

/* ═══════════════════════════════════════════════════════════════════════════
 * PNG output
 * ═══════════════════════════════════════════════════════════════════════════ */
static const unsigned char MAGMA[256][3] = {
{0,0,3},{0,0,4},{0,0,6},{1,0,7},{1,1,9},{1,1,11},{2,2,13},{2,2,15},
{3,3,17},{4,3,19},{4,4,21},{5,4,23},{6,5,25},{7,5,27},{8,6,29},{9,7,31},
{10,7,34},{11,8,36},{12,9,38},{13,10,40},{14,10,42},{15,11,44},{16,12,47},{17,12,49},
{18,13,51},{20,13,53},{21,14,56},{22,14,58},{23,15,60},{24,15,63},{26,16,65},{27,16,68},
{28,16,70},{30,16,73},{31,17,75},{32,17,77},{34,17,80},{35,17,82},{37,17,85},{38,17,87},
{40,17,89},{42,17,92},{43,17,94},{45,16,96},{47,16,98},{48,16,101},{50,16,103},{52,16,104},
{53,15,106},{55,15,108},{57,15,110},{59,15,111},{60,15,113},{62,15,114},{64,15,115},{66,15,116},
{67,15,117},{69,15,118},{71,15,119},{72,16,120},{74,16,121},{75,16,121},{77,17,122},{79,17,123},
{80,18,123},{82,18,124},{83,19,124},{85,19,125},{87,20,125},{88,21,126},{90,21,126},{91,22,126},
{93,23,126},{94,23,127},{96,24,127},{97,24,127},{99,25,127},{101,26,128},{102,26,128},{104,27,128},
{105,28,128},{107,28,128},{108,29,128},{110,30,129},{111,30,129},{113,31,129},{115,31,129},{116,32,129},
{118,33,129},{119,33,129},{121,34,129},{122,34,129},{124,35,129},{126,36,129},{127,36,129},{129,37,129},
{130,37,129},{132,38,129},{133,38,129},{135,39,129},{137,40,129},{138,40,129},{140,41,128},{141,41,128},
{143,42,128},{145,42,128},{146,43,128},{148,43,128},{149,44,128},{151,44,127},{153,45,127},{154,45,127},
{156,46,127},{158,46,126},{159,47,126},{161,47,126},{163,48,126},{164,48,125},{166,49,125},{167,49,125},
{169,50,124},{171,51,124},{172,51,123},{174,52,123},{176,52,123},{177,53,122},{179,53,122},{181,54,121},
{182,54,121},{184,55,120},{185,55,120},{187,56,119},{189,57,119},{190,57,118},{192,58,117},{194,58,117},
{195,59,116},{197,60,116},{198,60,115},{200,61,114},{202,62,114},{203,62,113},{205,63,112},{206,64,112},
{208,65,111},{209,66,110},{211,66,109},{212,67,109},{214,68,108},{215,69,107},{217,70,106},{218,71,105},
{220,72,105},{221,73,104},{222,74,103},{224,75,102},{225,76,102},{226,77,101},{228,78,100},{229,80,99},
{230,81,98},{231,82,98},{232,84,97},{234,85,96},{235,86,96},{236,88,95},{237,89,95},{238,91,94},
{238,93,93},{239,94,93},{240,96,93},{241,97,92},{242,99,92},{243,101,92},{243,103,91},{244,104,91},
{245,106,91},{245,108,91},{246,110,91},{246,112,91},{247,113,91},{247,115,92},{248,117,92},{248,119,92},
{249,121,92},{249,123,93},{249,125,93},{250,127,94},{250,128,94},{250,130,95},{251,132,96},{251,134,96},
{251,136,97},{251,138,98},{252,140,99},{252,142,99},{252,144,100},{252,146,101},{252,147,102},{253,149,103},
{253,151,104},{253,153,105},{253,155,106},{253,157,107},{253,159,108},{253,161,110},{253,162,111},{253,164,112},
{254,166,113},{254,168,115},{254,170,116},{254,172,117},{254,174,118},{254,175,120},{254,177,121},{254,179,123},
{254,181,124},{254,183,125},{254,185,127},{254,187,128},{254,188,130},{254,190,131},{254,192,133},{254,194,134},
{254,196,136},{254,198,137},{254,199,139},{254,201,141},{254,203,142},{253,205,144},{253,207,146},{253,209,147},
{253,210,149},{253,212,151},{253,214,152},{253,216,154},{253,218,156},{253,220,157},{253,221,159},{253,223,161},
{253,225,163},{252,227,165},{252,229,166},{252,230,168},{252,232,170},{252,234,172},{252,236,174},{252,238,176},
{252,240,177},{252,241,179},{252,243,181},{252,245,183},{251,247,185},{251,249,187},{251,250,189},{251,252,191},
};

static const unsigned char VIRIDIS[256][3] = {
{68,1,84},{68,2,86},{69,4,87},{69,5,89},{70,7,90},{70,8,92},
{70,10,93},{70,11,94},{71,13,96},{71,14,97},{71,16,99},{71,17,100},
{71,19,101},{72,20,103},{72,22,104},{72,23,105},{72,25,107},{72,26,108},
{72,28,110},{72,29,111},{72,31,112},{72,32,113},{72,34,115},{72,35,116},
{72,37,117},{72,38,118},{72,40,119},{72,41,121},{72,43,122},{71,45,123},
{71,46,124},{71,48,125},{71,49,126},{71,51,127},{71,52,128},{70,54,129},
{70,55,130},{70,57,131},{70,58,132},{70,60,133},{70,61,134},{69,63,135},
{69,64,135},{69,66,136},{69,67,137},{69,69,138},{68,70,139},{68,72,140},
{68,73,140},{67,75,141},{67,76,142},{67,78,142},{67,79,143},{66,81,144},
{66,82,144},{65,84,145},{65,85,146},{65,87,146},{64,88,147},{64,90,148},
{63,91,148},{63,93,149},{63,94,149},{62,96,150},{62,97,150},{61,99,151},
{61,100,151},{60,102,152},{60,103,152},{59,105,153},{59,106,153},
{58,108,153},{58,109,154},{57,111,154},{57,112,155},{56,114,155},
{56,115,155},{55,117,156},{55,118,156},{54,120,157},{54,121,157},
{53,123,157},{53,124,157},{52,126,158},{52,127,158},{51,129,158},
{50,130,159},{50,132,159},{49,133,159},{49,135,160},{48,136,160},
{48,138,160},{47,140,161},{46,141,161},{46,143,161},{45,144,162},
{45,146,162},{44,147,162},{44,149,162},{43,151,163},{42,152,163},
{42,154,163},{41,155,163},{41,157,164},{40,159,164},{39,160,164},
{39,162,164},{38,164,165},{38,165,165},{37,167,165},{36,169,165},
{36,170,165},{35,172,166},{35,174,166},{34,175,166},{34,177,166},
{33,179,167},{32,180,167},{32,182,167},{31,184,167},{30,185,167},
{30,187,167},{29,189,168},{28,191,168},{28,192,168},{27,194,168},
{27,196,168},{26,197,168},{25,199,168},{24,201,169},{24,203,169},
{23,204,169},{22,206,169},{22,208,169},{21,210,169},{20,212,169},
{20,213,169},{19,215,169},{18,217,169},{17,219,169},{17,221,169},
{16,222,169},{15,224,169},{14,226,169},{14,228,169},{13,230,169},
{12,232,169},{11,234,169},{11,235,169},{10,237,169},{9,239,169},
{8,241,169},{8,243,168},{7,245,168},{6,247,168},{5,248,168},
{5,250,167},{4,252,167},{3,254,167},{2,254,166},{2,254,164},
{1,254,163},{1,254,162},{0,254,160},{0,254,159},{0,253,157},
{1,253,156},{1,253,154},{2,253,153},{3,253,152},{4,252,150},
{5,252,149},{7,252,147},{8,252,146},{9,251,144},{11,251,143},
{12,251,141},{13,250,140},{15,250,138},{16,249,137},{18,249,135},
{20,249,133},{21,248,132},{23,248,130},{24,247,129},{26,247,127},
{28,246,125},{30,246,124},{31,245,122},{33,244,120},{35,244,119},
{37,243,117},{38,243,115},{40,242,113},{42,241,112},{44,241,110},
{46,240,108},{48,239,106},{50,239,104},{52,238,103},{54,237,101},
{56,237,99},{58,236,97},{61,235,95},{63,234,93},{65,233,91},
{67,233,89},{70,232,87},{72,231,85},{74,230,83},{77,229,81},
{79,228,79},{82,228,77},{84,227,75},{87,226,73},{89,225,71},
{92,224,69},{94,223,67},{97,222,64},{99,221,62},{102,220,60},
{105,219,58},{107,218,56},{110,217,54},{113,216,52},{116,215,50},
{118,214,48},{121,213,45},{124,211,43},{127,210,41},{130,209,39},
{133,208,37},{136,206,35},{139,205,33},{142,204,31},{145,203,28},
{148,201,26},{151,200,24},{154,198,22},{157,197,20},{160,196,18},
{164,194,16},{167,193,14},{170,191,12},{173,190,10},{177,188,8},
{180,187,6},{183,185,3},{187,184,1},{190,182,2},{193,181,4},
{197,179,6},{200,178,8},{203,176,10},{207,175,13},{210,173,16}
};

/* RdYlBu (diverging): red (low) → yellow (mid) → blue (high).  Used for the
 * signed vorticity panel, where 0 maps to the pale yellow centre.            */
static const unsigned char RDYLBU[256][3] = {
{165,0,38},{167,2,38},{169,4,38},{171,6,38},{173,8,38},{175,9,38},
{177,11,38},{179,13,38},{181,15,38},{183,17,38},{185,19,38},{187,21,38},
{189,23,38},{190,24,39},{192,26,39},{194,28,39},{196,30,39},{198,32,39},
{200,34,39},{202,36,39},{204,38,39},{206,40,39},{208,41,39},{210,43,39},
{212,45,39},{214,47,39},{216,49,40},{217,52,41},{218,54,42},{219,56,43},
{220,59,44},{221,61,45},{222,64,46},{224,66,47},{225,68,48},{226,71,49},
{227,73,51},{228,76,52},{229,78,53},{230,80,54},{231,83,55},{233,85,56},
{234,87,57},{235,90,58},{236,92,59},{237,95,60},{238,97,62},{239,99,63},
{241,102,64},{242,104,65},{243,107,66},{244,109,67},{244,112,68},{245,114,69},
{245,117,71},{245,119,72},{246,122,73},{246,124,74},{246,127,75},{247,129,76},
{247,132,78},{248,134,79},{248,137,80},{248,140,81},{249,142,82},{249,145,83},
{249,147,85},{250,150,86},{250,152,87},{250,155,88},{251,157,89},{251,160,91},
{251,163,92},{252,165,93},{252,168,94},{252,170,95},{253,173,96},{253,175,98},
{253,177,100},{253,179,102},{253,181,103},{253,183,105},{253,185,107},{253,187,109},
{253,189,111},{253,191,113},{253,193,115},{253,195,116},{253,197,118},{253,199,120},
{254,200,122},{254,202,124},{254,204,126},{254,206,127},{254,208,129},{254,210,131},
{254,212,133},{254,214,135},{254,216,137},{254,218,138},{254,220,140},{254,222,142},
{254,224,144},{254,225,146},{254,226,148},{254,228,150},{254,229,151},{254,230,153},
{254,231,155},{254,233,157},{254,234,159},{254,235,161},{254,236,162},{254,237,164},
{254,239,166},{255,240,168},{255,241,170},{255,242,172},{255,243,173},{255,245,175},
{255,246,177},{255,247,179},{255,248,181},{255,250,183},{255,251,185},{255,252,186},
{255,253,188},{255,254,190},{254,255,192},{253,254,194},{252,254,197},{251,253,199},
{250,253,201},{248,252,203},{247,252,206},{246,251,208},{245,251,210},{243,251,212},
{242,250,214},{241,250,217},{240,249,219},{239,249,221},{237,248,223},{236,248,226},
{235,247,228},{234,247,230},{233,246,232},{231,246,235},{230,245,237},{229,245,239},
{228,244,241},{226,244,244},{225,243,246},{224,243,248},{222,242,247},{220,241,247},
{218,240,246},{216,239,246},{214,238,245},{212,237,244},{209,236,244},{207,235,243},
{205,234,243},{203,233,242},{201,232,242},{199,231,241},{197,230,240},{195,229,240},
{193,228,239},{191,227,239},{189,226,238},{187,225,237},{185,224,237},{182,223,236},
{180,222,236},{178,221,235},{176,220,234},{174,219,234},{172,218,233},{170,216,233},
{168,214,232},{166,213,231},{163,211,230},{161,209,229},{159,208,228},{157,206,227},
{155,204,226},{153,202,225},{151,201,224},{148,199,223},{146,197,222},{144,195,221},
{142,194,220},{140,192,219},{138,190,218},{135,189,217},{133,187,217},{131,185,216},
{129,183,215},{127,182,214},{125,180,213},{122,178,212},{120,176,211},{118,175,210},
{116,173,209},{114,171,208},{112,169,207},{110,166,206},{109,164,204},{107,162,203},
{105,160,202},{103,158,201},{101,155,200},{99,153,199},{98,151,198},{96,149,196},
{94,147,195},{92,144,194},{90,142,193},{88,140,192},{87,138,191},{85,136,190},
{83,133,189},{81,131,187},{79,129,186},{77,127,185},{75,125,184},{74,122,183},
{72,120,182},{70,118,181},{69,116,179},{68,113,178},{67,111,177},{66,108,176},
{65,106,175},{65,103,173},{64,101,172},{63,98,171},{62,96,170},{62,94,168},
{61,91,167},{60,89,166},{59,86,165},{58,84,164},{58,81,162},{57,79,161},
{56,76,160},{55,74,159},{54,71,158},{54,69,156},{53,66,155},{52,64,154},
{51,61,153},{51,59,151},{50,56,150},{49,54,149},
};

static void
write_png_2panel(const char *fname,
                 const real *tl, real tl_min, real tl_max,
                 const real *tr, real tr_min, real tr_max,
                 int N, int tr_rdylbu)
{
    const int W = 2 * N, H = N;
    unsigned char *pixels = (unsigned char *)malloc((size_t)W * H * 3);
    if (!pixels) { fprintf(stderr, "OOM in write_png_2panel\n"); return; }

    /* Clamp vmin so floor-clipped outliers (1e-14) don't compress the useful
     * data range into the near-black end of the colormap.  Only for
     * non-negative fields — signed fields (e.g. vorticity) keep their full
     * range so the negative half isn't collapsed into the dark end.         */
    real tl_vmin = tl_min, tr_vmin = tr_min;
    if (tl_vmin >= 0. && tl_vmin < tl_max * 1e-3) tl_vmin = tl_max * 1e-3;
    if (tr_vmin >= 0. && tr_vmin < tr_max * 1e-3) tr_vmin = tr_max * 1e-3;

    for (int r = 0; r < H; r++) {
        int phys_j = N - 1 - r;   /* flip y */
        for (int c = 0; c < W; c++) {
            const real *panel = (c < N) ? tl : tr;
            real vmin = (c < N) ? tl_vmin : tr_vmin;
            real vmax = (c < N) ? tl_max  : tr_max;
            int   ci   = (c < N) ? c : c - N;
            const unsigned char (*cmap)[3] =
                (c < N) ? MAGMA : (tr_rdylbu ? RDYLBU : VIRIDIS);

            real t = (panel[phys_j * N + ci] - vmin) / (vmax - vmin + 1e-30);
            if (!(t > 0.)) t = 0.;
            if (t > 1.)    t = 1.;
            int idx = (int)(t * 255.);
            if (idx < 0) idx = 0; if (idx > 255) idx = 255;
            unsigned char *pix = pixels + (r * W + c) * 3;
            pix[0] = cmap[idx][0]; pix[1] = cmap[idx][1]; pix[2] = cmap[idx][2];
        }
    }
    stbi_write_png(fname, W, H, 3, pixels, W * 3);
    free(pixels);
    printf("  Saved %s\n", fname);
}

/* Single-panel PNG with log-scale colormap — for Mach number visualisation.
 * vmin_log : floor value (e.g. 1e-6); values below are clamped to vmin_log.
 * vmax     : upper bound (e.g. peak Mach from IC).                          */
static void
write_png_1panel(const char *fname,
                 const real *data, real vmin, real vmax,
                 int N)
{
    unsigned char *pixels = (unsigned char *)malloc((size_t)N * N * 3);
    if (!pixels) { fprintf(stderr, "OOM in write_png_1panel\n"); return; }
    real range = vmax - vmin;
    for (int r = 0; r < N; r++) {
        int phys_j = N - 1 - r;
        for (int c = 0; c < N; c++) {
            real val = data[phys_j * N + c];
            real t = (val - vmin) / (range + 1e-30);
            if (!(t > 0.)) t = 0.;
            if (t > 1.)    t = 1.;
            /* Jet colormap: dark-blue -> cyan -> green -> yellow -> red */
            real cr = fmax(0., fmin(1., 1.5 - fabs(4.*t - 3.)));
            real cg = fmax(0., fmin(1., 1.5 - fabs(4.*t - 2.)));
            real cb = fmax(0., fmin(1., 1.5 - fabs(4.*t - 1.)));
            unsigned char *pix = pixels + (r * N + c) * 3;
            pix[0] = (unsigned char)(cr * 255.);
            pix[1] = (unsigned char)(cg * 255.);
            pix[2] = (unsigned char)(cb * 255.);
        }
    }
    stbi_write_png(fname, N, N, 3, pixels, N * 3);
    free(pixels);
    printf("  Saved %s\n", fname);
}

/* ═══════════════════════════════════════════════════════════════════════════
 * main
 *
 * Usage:
 *   ./rt_dg N               Circular Sod shock  (domain [0,1]^2)
 *   ./rt_dg N Ma            Isentropic vortex, scaled by Mach Ma
 *                            (original parameterization, domain [0,1]^2)
 *   ./rt_dg N pv            Paper isentropic vortex, stationary
 *                            (Barsukow et al. 2025 Sec 6.2.1, [0,10]^2)
 *   ./rt_dg N pvmove        Paper isentropic vortex, moving (u0=v0=1)
 *   ./rt_dg N kh [Ma]      Kelvin-Helmholtz instability
 *                            (Barsukow et al. 2025 Sec 6.2.5, default Ma=0.01)
 *   ./rt_dg N dsl [Ma]     Doubly periodic shear layer
 *                            (Bell, Colella & Glaz 1989, [0,1]^2, default Ma=0.1)
 *   ./rt_dg N va [Mref]    Vortex-acoustic wave interaction
 *                            (Coiffier 2025, Sec 6.6.3, default Mref=0.1)
 *   ./rt_dg N lmv [eps]    Low-Mach vortex (Barsukow et al. Sec 2, default eps=0.1)
 *                            Stationary vortex; Mach set via p0=1/(γε²)-1/2
 * ═══════════════════════════════════════════════════════════════════════════ */

static void q1_totals(const real *d_Uv, int N, int Np, double *mass, double *ener, real h)
{
    size_t sz = (size_t)NVAR*Np*Np*sizeof(real);
    real *hU = (real*)malloc(sz);
    cudaMemcpy(hU, d_Uv, sz, cudaMemcpyDeviceToHost);
    double m = 0., e = 0.;
    for (int j = 0; j < N; j++) for (int i = 0; i < N; i++) {
        int jp = j+G_GHOST, ip = i+G_GHOST;
        m += (double)hU[IDX_P(Q_RHO, jp, ip, Np)];
        e += (double)hU[IDX_P(Q_E,   jp, ip, Np)];
    }
    *mass = m*h*h; *ener = e*h*h; free(hU);
}

int main(int argc, char **argv)
{
    int   N    = (argc > 1) ? atoi(argv[1]) : 400;

    /* Determine test mode */
    const char *mode_str = (argc > 2) ? argv[2] : "";
    real Ma = 0.;
    int do_vortex       = 0;
    int do_paper_vortex = 0;  /* pv: stationary, pvmove: moving */
    int do_kh           = 0;
    int do_dsl          = 0;  /* dsl: doubly periodic shear layer (BCG 1989) */
    int do_vacwav       = 0;
    int do_lmv          = 0;  /* lmv: low-Mach vortex (Barsukow Sec 2) */
    int pv_moving       = 0;  /* 1 → pvmove: u0=v0=1, tf=10 */
    real kh_mach       = 0.01; /* KH convective Mach number */
    real dsl_mach      = 0.1;  /* shear-layer Mach number */
    real va_mach       = 0.1;  /* vortex-acoustic reference Mach number */
    real lmv_eps       = 0.1;  /* low-Mach vortex Mach number */

    /* Check string modes first, then fall back to numeric Ma */
    if (strcmp(mode_str, "pv") == 0) {
        do_paper_vortex = 1;
    } else if (strcmp(mode_str, "pvmove") == 0) {
        do_paper_vortex = 1;
        pv_moving       = 1;
    } else if (strcmp(mode_str, "kh") == 0) {
        do_kh = 1;
        if (argc > 3) {
            char *endp;
            real maybe = strtod(argv[3], &endp);
            if (endp != argv[3] && *endp == '\0' && maybe > 0.)
                kh_mach = maybe;
        }
    } else if (strcmp(mode_str, "dsl") == 0) {
        do_dsl = 1;
        if (argc > 3) {
            char *endp;
            real maybe = strtod(argv[3], &endp);
            if (endp != argv[3] && *endp == '\0' && maybe > 0.)
                dsl_mach = maybe;
        }
    } else if (strcmp(mode_str, "va") == 0) {
        do_vacwav = 1;
        if (argc > 3) {
            char *endp;
            real maybe = strtod(argv[3], &endp);
            if (endp != argv[3] && *endp == '\0' && maybe > 0.)
                va_mach = maybe;
        }
    } else if (strcmp(mode_str, "lmv") == 0) {
        do_lmv = 1;
        if (argc > 3) {
            char *endp;
            real maybe = strtod(argv[3], &endp);
            if (endp != argv[3] && *endp == '\0' && maybe > 0.)
                lmv_eps = maybe;
        }
    } else {
        /* Legacy: numeric Ma for original vortex */
        char *endp;
        real maybe_ma = strtod(mode_str, &endp);
        if (endp != mode_str && *endp == '\0' && maybe_ma > 0.) {
            Ma = maybe_ma;
            do_vortex = 1;
        }
    }

    /* Sod = the fallback mode (no other test selected) */
    int do_sod = !(do_vortex || do_paper_vortex || do_kh ||
                   do_dsl || do_vacwav || do_lmv);

    /* Domain size and run time per mode */
    real L_domain = 1.;           /* default: [0,1]^2 */
    if (do_paper_vortex) L_domain = 10.;
    if (do_kh)           L_domain =  2.;
    /* dsl, lmv, vacwav: [0,1]^2 — L_domain stays 1 */

#ifndef CFL_DEFAULT
#define CFL_DEFAULT 0.65
#endif
    real CFL   = CFL_DEFAULT;
    real t_end;
    if (do_vortex) {
        /* One full vortex orbital period: T = 2πR / |u|_max = 2πR / Ma
         * (peak velocity = Ma·c∞ = Ma since c∞=1; R = VORTEX_R = 0.2)
         * This gives the same number of vortex rotations regardless of Ma. */
        t_end = 2. * (real)M_PI * VORTEX_R / Ma;
    } else if (do_paper_vortex) t_end = pv_moving ? 10.0 : 1.0;
    /* t_end for KH: 8 convective times = 8 * L_domain / (2*kh_mach)  */
    else if (do_kh)      t_end = 8. * L_domain / (2. * kh_mach);
    /* t_end for DSL: dimensionless t*=1.2 (rollup), scaled by 1/Ma since |u|~Ma */
    else if (do_dsl)     t_end = 2.2 / dsl_mach;
    /* t_end for VA: run to 3.5 s (matches thesis Fig. 6.6) */
    else if (do_vacwav)  t_end = 3.5;
    /* t_end for LMV: one full vortex period T = 2π*0.2 / 1.0 ≈ 1.257 */
    else if (do_lmv)     t_end = 2. * (real)M_PI * 0.2;
    else                 t_end = 1.0;   /* Sod */

    real h    = L_domain / N;
    const int g = G_GHOST;
    int Np      = N + 2 * g;

    /* Background velocity for paper vortex */
    real u0 = pv_moving ? 1. : 0.;
    real v0 = pv_moving ? 1. : 0.;

    /* GPU info */
    int dev; cudaDeviceProp prop;
    CK(cudaGetDevice(&dev));
    CK(cudaGetDeviceProperties(&prop, dev));
    printf("========================================================\n");
    printf("  Device  : %s\n", prop.name);
    if (do_vortex)
        printf("  P1/P0 DG + HLLC  --  Isentropic Vortex  Ma=%.3f\n", Ma);
    else if (do_paper_vortex)
        printf("  P1/P0 DG + HLLC  --  Paper Isentropic Vortex  [0,10]^2  %s\n",
               pv_moving ? "MOVING (u0=v0=1)" : "STATIONARY");
    else if (do_kh)
        printf("  P1/P0 DG + HLLC  --  Kelvin-Helmholtz  Ma=0.01  [0,2]^2\n");
    else if (do_dsl)
        printf("  P1/P0 DG + HLLC  --  Doubly Periodic Shear Layer  Ma=%.3f  [0,1]^2\n", dsl_mach);
    else if (do_vacwav)
        printf("  P1/P0 DG + HLLC  --  Vortex-Acoustic Wave  Mref=%.3f  [0,1]^2\n", va_mach);
    else if (do_lmv)
        printf("  P1/P0 DG + HLLC  --  Low-Mach Vortex  eps=%.4f  [0,1]^2\n", lmv_eps);
    else
        printf("  BILINEAR-Q1 DG + HLLC  --  Circular Sod  (Zhang-Shu positivity limiter: %s)\n",
               SHOCK_LIMIT ? "on" : "off");
    printf("  N=%dx%d  Np=%d  h=%.5f  L=%.1f  CFL=%.2f  t_end=%.3f\n",
           N, N, Np, h, L_domain, CFL, t_end);
#ifdef NEDELEC
    printf("  DOFs: rho(P0) + [mxa,mxsy, mya,mysx](NEDELEC / curl-conforming momentum) + E(P0)  [DOUBLE]\n");
    printf("  >> Nedelec subspace: mxs=mys=0  (div-free in-cell, carries vorticity)\n");
#else
    printf("  DOFs: FULL BILINEAR Q1 for ALL fields — rho,rhou,rhov,E each {1,x,y,xy} = 16/cell  [DOUBLE precision]\n");
    printf("  faces: compact Hermite recon (value+slope+neighbour), HLLC;  volume: 2x2 Gauss;  limiter: Zhang-Shu positivity\n");
#endif
    printf("========================================================\n");

    size_t sz1  = (size_t)N  * N  * sizeof(real);
    size_t sz6  = (size_t)NVAR * sz1;
    size_t sz6p = (size_t)NVAR * Np * Np * sizeof(real);

    real *d_U, *d_U0, *d_U1, *d_U2, *d_RHS, *d_lam, *d_tmp;
    real *d_rho_out, *d_p_out, *d_err_rho, *d_err_p, *d_err_ru, *d_err_rv, *d_err_re;
    CK(cudaMalloc(&d_U,       sz6p));
    CK(cudaMalloc(&d_U0,      sz6p));
    CK(cudaMalloc(&d_U1,      sz6p));
    CK(cudaMalloc(&d_U2,      sz6p));
    CK(cudaMalloc(&d_RHS,     sz6));
    CK(cudaMalloc(&d_lam,     sz1));
    CK(cudaMalloc(&d_rho_out, sz1));
    CK(cudaMalloc(&d_p_out,   sz1));
    CK(cudaMalloc(&d_err_rho, sz1));
    CK(cudaMalloc(&d_err_p,   sz1));
    CK(cudaMalloc(&d_err_ru,  sz1));
    CK(cudaMalloc(&d_err_rv,  sz1));
    CK(cudaMalloc(&d_err_re,  sz1));
    CK(cudaMalloc(&d_tmp,     GS_NBLK * sizeof(real)));

    real *h_rho = (real*)malloc(sz1);
    real *h_p   = (real*)malloc(sz1);

    /* Color ranges — all cases use data-driven bounds computed from IC */
    real RHO_MIN = 0., RHO_MAX = 0., P_MIN = 0., P_MAX = 0.;
    static char _lmv_buf[32];
    snprintf(_lmv_buf, sizeof(_lmv_buf), "p1_lmv%04d", (int)round(lmv_eps*10000));
    static char _va_buf[32];
    snprintf(_va_buf, sizeof(_va_buf), "p1_va%04d", (int)round(va_mach*10000));
    static char _kh_buf[32];
    snprintf(_kh_buf, sizeof(_kh_buf), "p1_kh%03d", (int)round(kh_mach*1000));
    static char _dsl_buf[32];
    snprintf(_dsl_buf, sizeof(_dsl_buf), "p1_dsl%03d", (int)round(dsl_mach*1000));
    const char *prefix = do_vortex       ? "p1_vortex"  :
                         do_paper_vortex ? "p1_pvortex" :
                         do_vacwav       ? _va_buf :
                         do_kh           ? _kh_buf :
                         do_dsl          ? _dsl_buf :
                         do_lmv          ? _lmv_buf : "p1_sod";

    /* Helper: compute and print L2 errors (rho, rho*u, rho*v, rho*E — matches paper Tables 2 & 3) */
#define PRINT_L2_ERR() do { \
    if (do_vortex) { \
        compute_l2_err<<<GS_NBLK,BLOCK1D>>>(d_U, d_err_rho, d_err_p, N, Np, h, Ma); \
        CK(cudaDeviceSynchronize()); \
        real _sr = gpu_sum(d_err_rho, d_tmp, N*N); \
        real _sp = gpu_sum(d_err_p,   d_tmp, N*N); \
        printf("    L2_rho=%.4e  L2_p=%.4e\n", sqrt(_sr), sqrt(_sp)); \
    } else if (do_paper_vortex) { \
        /* exact at tf=1 (stationary) and tf=10 (moving) both equal the IC */ \
        compute_l2_err_paper<<<GS_NBLK,BLOCK1D>>>(d_U, d_err_rho, d_err_ru, \
                                                    d_err_rv, d_err_re, \
                                                    N, Np, h, u0, v0); \
        CK(cudaDeviceSynchronize()); \
        real _sr  = gpu_sum(d_err_rho, d_tmp, N*N); \
        real _sru = gpu_sum(d_err_ru,  d_tmp, N*N); \
        real _srv = gpu_sum(d_err_rv,  d_tmp, N*N); \
        real _sre = gpu_sum(d_err_re,  d_tmp, N*N); \
        printf("    L2_rho=%.4e  L2_ru=%.4e  L2_rv=%.4e  L2_rE=%.4e\n", \
               sqrt(_sr), sqrt(_sru), sqrt(_srv), sqrt(_sre)); \
    } else if (do_lmv) { \
        compute_l2_err_lmv<<<GS_NBLK,BLOCK1D>>>(d_U, d_err_rho, d_err_p, d_err_ru, N, Np, h, lmv_eps); \
        CK(cudaDeviceSynchronize()); \
        real _sr = gpu_sum(d_err_rho, d_tmp, N*N); \
        real _sp = gpu_sum(d_err_p,   d_tmp, N*N); \
        real _sv = gpu_sum(d_err_ru,  d_tmp, N*N); \
        printf("    L2_rho=%.4e  L2_p=%.4e  L2_vel=%.4e\n", sqrt(_sr), sqrt(_sp), sqrt(_sv)); \
    } \
} while(0)

#define WRITE_FRAME(idx) do { \
    if (do_vacwav) { \
        extract_mach<<<GS_NBLK,BLOCK1D>>>(d_U, d_rho_out, N, Np); \
        CK(cudaDeviceSynchronize()); \
        real _flo = gpu_min(d_rho_out, d_tmp, N*N); \
        real _fhi = gpu_max(d_rho_out, d_tmp, N*N); \
        CK(cudaMemcpy(h_rho, d_rho_out, sz1, cudaMemcpyDeviceToHost)); \
        char _fn[64]; sprintf(_fn, "figures/%s_%04d.png", prefix, (idx)); \
        write_png_1panel(_fn, h_rho, _flo, _fhi, N); \
    } else { \
    if (do_dsl) \
        extract_vort  <<<GS_NBLK,BLOCK1D>>>(d_U, d_rho_out, d_p_out, N, Np, h); \
    else if (do_kh) \
        extract_rho_v <<<GS_NBLK,BLOCK1D>>>(d_U, d_rho_out, d_p_out, N, Np); \
    else if (do_lmv) \
        extract_magv_p<<<GS_NBLK,BLOCK1D>>>(d_U, d_rho_out, d_p_out, N, Np); \
    else if (do_sod) \
        extract_rho_vort<<<GS_NBLK,BLOCK1D>>>(d_U, d_rho_out, d_p_out, N, Np, h); \
    else \
        extract_rho_p <<<GS_NBLK,BLOCK1D>>>(d_U, d_rho_out, d_p_out, N, Np); \
    CK(cudaDeviceSynchronize()); \
    real _rlo = gpu_min(d_rho_out, d_tmp, N*N); \
    real _rhi = gpu_max(d_rho_out, d_tmp, N*N); \
    real _plo = gpu_min(d_p_out,   d_tmp, N*N); \
    real _phi = gpu_max(d_p_out,   d_tmp, N*N); \
    if (do_dsl) { /* signed v (left) & vorticity (right): symmetric about 0 */ \
        real _wv = fmax(fabs(_rlo), fabs(_rhi)); _rlo = -_wv; _rhi = _wv; \
        real _ww = fmax(fabs(_plo), fabs(_phi)); _plo = -_ww; _phi = _ww; \
    } else if (do_sod) { /* vorticity (right) symmetric about 0 */ \
        real _ww = fmax(fabs(_plo), fabs(_phi)); _plo = -_ww; _phi = _ww; \
    } \
    CK(cudaMemcpy(h_rho, d_rho_out, sz1, cudaMemcpyDeviceToHost)); \
    CK(cudaMemcpy(h_p,   d_p_out,   sz1, cudaMemcpyDeviceToHost)); \
    char _fn[64]; sprintf(_fn, "figures/%s_%04d.png", prefix, (idx)); \
    write_png_2panel(_fn, h_rho, _rlo, _rhi, \
                          h_p,   _plo, _phi, N, do_sod); \
    } /* end non-va */ \
    PRINT_L2_ERR(); \
} while(0)

    /* --- IC and ghost cells --- */
    /* Zero the whole padded array first so the ρ,E slope DOFs (Q_RHO_SX/SY,
     * Q_E_SX/SY) start at 0 for every IC kernel (piecewise-const projection;
     * the Sod IC is cell-constant per region so this is exact). */
    CK(cudaMemset(d_U, 0, sz6p));
    if (do_vortex)
        ic_vortex      <<<GS_NBLK, BLOCK1D>>>(d_U, N, Np, h, Ma);
    else if (do_paper_vortex)
        ic_vortex_paper<<<GS_NBLK, BLOCK1D>>>(d_U, N, Np, h, u0, v0);
    else if (do_kh)
        ic_kh          <<<GS_NBLK, BLOCK1D>>>(d_U, N, Np, h, kh_mach);
    else if (do_dsl)
        ic_dsl         <<<GS_NBLK, BLOCK1D>>>(d_U, N, Np, h, dsl_mach);
    else if (do_vacwav)
        ic_vacwav      <<<GS_NBLK, BLOCK1D>>>(d_U, N, Np, h, va_mach);
    else if (do_lmv)
        ic_lmv         <<<GS_NBLK, BLOCK1D>>>(d_U, N, Np, h, lmv_eps);
    else
        ic_kernel      <<<GS_NBLK, BLOCK1D>>>(d_U, N, Np, h);
    SHU_LIMIT(d_U);
    LIMIT_POS(d_U);
    apply_bc  <<<GS_NBLK, BLOCK1D>>>(d_U, N, Np);
    NEDELEC_PROJECT(d_U);   /* project IC onto Nédélec subspace (no-op if !NEDELEC) */
    CK(cudaDeviceSynchronize());

    /* Compute data-driven color bounds from IC for all cases */
    {
        if (do_vacwav) {
            /* Linear Mach: find peak Mach from IC */
            extract_mach<<<GS_NBLK, BLOCK1D>>>(d_U, d_rho_out, N, Np);
            CK(cudaDeviceSynchronize());
            real mach_hi = gpu_max(d_rho_out, d_tmp, N*N);
            if (mach_hi < 1e-6) mach_hi = va_mach;
            RHO_MIN = 0.;
            RHO_MAX = mach_hi * 1.1;
            printf("  Mach range: [0, %.4f]\n", RHO_MAX);
        } else {
            if (do_dsl)
                extract_vort<<<GS_NBLK, BLOCK1D>>>(d_U, d_rho_out, d_p_out, N, Np, h);
            else if (do_kh)
                extract_rho_v<<<GS_NBLK, BLOCK1D>>>(d_U, d_rho_out, d_p_out, N, Np);
            else if (do_sod)
                extract_rho_vort<<<GS_NBLK, BLOCK1D>>>(d_U, d_rho_out, d_p_out, N, Np, h);
            else
                extract_rho_p<<<GS_NBLK, BLOCK1D>>>(d_U, d_rho_out, d_p_out, N, Np);
            CK(cudaDeviceSynchronize());
            RHO_MIN = gpu_min(d_rho_out, d_tmp, N*N);
            RHO_MAX = gpu_max(d_rho_out, d_tmp, N*N);
            P_MIN   = gpu_min(d_p_out,   d_tmp, N*N);
            P_MAX   = gpu_max(d_p_out,   d_tmp, N*N);
            /* 10% padding */
            real rho_pad = (RHO_MAX - RHO_MIN) * 0.1;
            real p_pad   = (P_MAX   - P_MIN  ) * 0.1;
            if (rho_pad < 1e-6) rho_pad = 0.1;
            if (p_pad   < 1e-6) p_pad   = 0.1;
            RHO_MIN -= rho_pad;  RHO_MAX += rho_pad;
            P_MIN   -= p_pad;    P_MAX   += p_pad;
            printf("  Color range  rho: [%.4f, %.4f]  p: [%.4f, %.4f]\n",
                   RHO_MIN, RHO_MAX, P_MIN, P_MAX);
        }
    }

    WRITE_FRAME(0);

    /* precompute RHS for first step */
    CK(cudaMemset(d_RHS, 0, sz6));
    COMPUTE_RHS(d_U);
    CK(cudaDeviceSynchronize());
    real lam_max = gpu_max(d_lam, d_tmp, N*N);

    const int  N_FRAMES = do_vacwav ? 14 : 10;  /* va: 0.25s spacing to match Fig. 6.6 */
    int   frame   = 1;
    real t_next  = t_end / N_FRAMES;
    int   step    = 0;
    real t       = 0.;
    struct timespec ts0, ts1;
    clock_gettime(CLOCK_MONOTONIC, &ts0);

    double m0_, e0_; q1_totals(d_U, N, Np, &m0_, &e0_, h);
    while (t < t_end) {
        if (!(lam_max > 0.)) {
            fprintf(stderr, "NaN/Inf lam_max at step %d t=%.5f\n", step, t); break;
        }
        real dt = CFL * h / lam_max;
        real t_target = (t_next < t_end) ? t_next : t_end;
        if (t + dt > t_target) dt = t_target - t;
        if (dt < 1e-14) { fprintf(stderr, "dt underflow\n"); WRITE_FRAME(99); break; }

        CK(cudaMemcpy(d_U0, d_U, sz6p, cudaMemcpyDeviceToDevice));

        /* stage 1 */
        rk3_s1<<<GS_NBLK, BLOCK1D>>>(d_U1, d_U0, d_RHS, dt, N, Np);
        SHU_LIMIT(d_U1);
        LIMIT_POS(d_U1);
        apply_bc<<<GS_NBLK, BLOCK1D>>>(d_U1, N, Np);
        NEDELEC_PROJECT(d_U1);

        /* stage 2 */
        CK(cudaMemset(d_RHS, 0, sz6));
        COMPUTE_RHS(d_U1);
        rk3_s2<<<GS_NBLK, BLOCK1D>>>(d_U2, d_U0, d_U1, d_RHS, dt, N, Np);
        SHU_LIMIT(d_U2);
        LIMIT_POS(d_U2);
        apply_bc<<<GS_NBLK, BLOCK1D>>>(d_U2, N, Np);
        NEDELEC_PROJECT(d_U2);

        /* stage 3 */
        CK(cudaMemset(d_RHS, 0, sz6));
        COMPUTE_RHS(d_U2);
        rk3_s3<<<GS_NBLK, BLOCK1D>>>(d_U, d_U0, d_U2, d_RHS, dt, N, Np);
        SHU_LIMIT(d_U);
        LIMIT_POS(d_U);
        apply_bc<<<GS_NBLK, BLOCK1D>>>(d_U, N, Np);
        NEDELEC_PROJECT(d_U);

        /* update RHS and lam for next step */
        CK(cudaMemset(d_RHS, 0, sz6));
        COMPUTE_RHS(d_U);
        lam_max = gpu_max(d_lam, d_tmp, N*N);
        if (!(lam_max > 0.)) { fprintf(stderr, "blow-up step %d\n", step); break; }

        t += dt; step++;

        if (t >= t_next - 1e-12 && frame <= N_FRAMES) {
            clock_gettime(CLOCK_MONOTONIC, &ts1);
            real el = (ts1.tv_sec-ts0.tv_sec) + (ts1.tv_nsec-ts0.tv_nsec)*1e-9;
            printf("  frame %2d/%d  step %5d  t=%.5f  lam=%.3e  elapsed=%.1fs\n",
                   frame, N_FRAMES, step, t, lam_max, el);
            fflush(stdout);
            WRITE_FRAME(frame);
            frame++;
            t_next = frame * t_end / N_FRAMES;
        } else if (step % 100 == 0) {
            clock_gettime(CLOCK_MONOTONIC, &ts1);
            real el = (ts1.tv_sec-ts0.tv_sec) + (ts1.tv_nsec-ts0.tv_nsec)*1e-9;
            printf("  step %5d  t=%.5f  dt=%.3e  lam=%.3e  elapsed=%.1fs\n",
                   step, t, dt, lam_max, el);
            fflush(stdout);
        }
    }
    /* Final L2 error summary */
    if (do_vortex) {
        compute_l2_err<<<GS_NBLK,BLOCK1D>>>(d_U, d_err_rho, d_err_p, N, Np, h, Ma);
        CK(cudaDeviceSynchronize());
        real sr = gpu_sum(d_err_rho, d_tmp, N*N);
        real sp = gpu_sum(d_err_p,   d_tmp, N*N);
        printf("  ── FINAL L2 ERRORS ──────────────────────────────────\n");
        printf("  Ma=%.3f  N=%d  h=%.5f  t=%.4f\n", Ma, N, h, t);
        printf("  L2(rho) = %.6e\n", sqrt(sr));
        printf("  L2(p)   = %.6e\n", sqrt(sp));
        printf("  ────────────────────────────────────────────────────\n");
    } else if (do_paper_vortex) {
        /* exact at tf=1 (stationary) and tf=10 (moving) both equal the IC:
         * the vortex traverses exactly one period on the [0,10]^2 periodic domain. */
        compute_l2_err_paper<<<GS_NBLK,BLOCK1D>>>(d_U, d_err_rho, d_err_ru,
                                                    d_err_rv, d_err_re,
                                                    N, Np, h, u0, v0);
        CK(cudaDeviceSynchronize());
        real sr  = gpu_sum(d_err_rho, d_tmp, N*N);
        real sru = gpu_sum(d_err_ru,  d_tmp, N*N);
        real srv = gpu_sum(d_err_rv,  d_tmp, N*N);
        real sre = gpu_sum(d_err_re,  d_tmp, N*N);
        printf("  ── FINAL L2 ERRORS  [Paper Vortex, Sec 6.2.1]  ─────\n");
        printf("  N=%d  h=%.5f  L=10  t=%.4f  %s\n", N, h, t,
               pv_moving ? "(moving u0=v0=1)" : "(stationary)");
        printf("  L2(rho) = %.6e%s\n", sqrt(sr),
               pv_moving ? "" : "   (superconvergence → ~2nd order expected)");
        printf("  L2(rho*u)= %.6e\n", sqrt(sru));
        printf("  L2(rho*v)= %.6e\n", sqrt(srv));
        printf("  L2(rho*E)= %.6e\n", sqrt(sre));
        printf("  ────────────────────────────────────────────────────\n");
    }

    CK(cudaDeviceSynchronize());
    clock_gettime(CLOCK_MONOTONIC, &ts1);
    real wall = (ts1.tv_sec-ts0.tv_sec) + (ts1.tv_nsec-ts0.tv_nsec)*1e-9;
    printf("  Done: %d steps  t=%.6f  wall=%.2fs\n", step, t, wall);
    { int pc_; CK(cudaMemcpyFromSymbol(&pc_, g_pfix_count, sizeof(int)));
      double m1_, e1_; q1_totals(d_U, N, Np, &m1_, &e1_, h);
      printf("  p-mean energy-fix activations: %d\n", pc_);
      printf("  conservation drift: mass %.3e  energy %.3e (relative)\n",
             fabs(m1_-m0_)/fabs(m0_), fabs(e1_-e0_)/fabs(e0_)); }
    printf("  Output: %s_0000.png .. %s_%04d.png\n", prefix, prefix, N_FRAMES);
    printf("  Left panel: %s (magma), Right panel: %s (viridis)\n",
           do_lmv ? "|velocity|" : do_dsl ? "v-velocity" : "density",
           (do_dsl || do_sod) ? "vorticity" : do_kh ? "v-velocity" : "pressure");

    free(h_rho); free(h_p);
    cudaFree(d_U); cudaFree(d_U0); cudaFree(d_U1); cudaFree(d_U2);
    cudaFree(d_RHS); cudaFree(d_lam); cudaFree(d_tmp);
    cudaFree(d_rho_out); cudaFree(d_p_out);
    cudaFree(d_err_rho); cudaFree(d_err_p); cudaFree(d_err_ru);
    cudaFree(d_err_rv); cudaFree(d_err_re);
    return 0;
}
