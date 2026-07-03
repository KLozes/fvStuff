/*
 * enriched_fv_euler_2d.cu
 * Enriched finite-volume scheme for low-Mach accurate compressible Euler
 * on Cartesian quadrangular meshes.
 *
 * Reference:
 *   J. Jung & V. Perrier,
 *   "A finite volume scheme accurate at low Mach number on quadrangular
 *   mesh by space velocity enrichment",
 *   ECCOMAS 2024, Lisbon, Portugal (hal-04836529).
 *
 * Method:
 *   density ρ   : P0 (piecewise constant, cell average)
 *   momentum m  : P0 cell average (mxa, mya) + divergence-free enrichment α
 *   energy E    : P0 (piecewise constant, cell average)
 *   Riemann     : HLLC (Toro 1994, standard)
 *   Time        : SSP-RK3
 *
 * Divergence-free enrichment (per cell):
 *   The velocity field is enriched with a single divergence-free mode α:
 *     u(x,y) = mxa/ρ + α*(y-yc)/(ρ·h²)   ← NOTE: sign defined so that
 *     v(x,y) = mya/ρ - α*(x-xc)/(ρ·h²)      div(δu,δv) = 0 identically
 *
 *   Face-normal velocities (midpoints) for reconstruction:
 *     Right  (x=xc+h/2, y=yc): u_n = mxa/ρ,  v_t = mya/ρ + α/(2ρ)
 *     Left   (x=xc-h/2, y=yc): u_n = mxa/ρ,  v_t = mya/ρ - α/(2ρ)  [outward -x]
 *     Top    (x=xc, y=yc+h/2): v_n = mya/ρ,  u_t = mxa/ρ - α/(2ρ)
 *     Bottom (x=xc, y=yc-h/2): v_n = mya/ρ,  u_t = mxa/ρ + α/(2ρ)  [outward -y]
 *
 *   The update for α is obtained by projecting the momentum RHS onto the
 *   divergence-free subspace:
 *     dα/dt = (d(mxs)/dt − d(mys)/dt) / 2,  where mxs=α, mys=−α
 *
 * DOFs per cell: NVAR=5 — (ρ, mxa, mya, E, α)
 *
 * Compile:
 *   nvcc -O3 -arch=native -o enriched_fv enriched_fv_euler_2d.cu -lm
 *
 * Usage:
 *   ./enriched_fv N              Circular Sod shock tube, [0,1]^2
 *   ./enriched_fv N Ma           Isentropic vortex (scaled), [0,1]^2
 *   ./enriched_fv N pv           Paper isentropic vortex, STATIONARY [0,10]^2
 *   ./enriched_fv N pvmove       Paper isentropic vortex, MOVING u0=v0=1
 *   ./enriched_fv N kh [Ma]      Kelvin-Helmholtz instability
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <cuda_runtime.h>

/* ── compile-time constants ─────────────────────────────────────────────── */
#define GAMMA_V  1.4f
#define BLOCK1D  256
#define GS_NBLK  256

/* DOF indices in the padded state array (NVAR = 5) */
#define NVAR  5
#define Q_RHO  0   /* density                      (P0)          */
#define Q_MXA  1   /* x-momentum cell average      (P0)          */
#define Q_MYA  2   /* y-momentum cell average      (P0)          */
#define Q_E    3   /* total energy                 (P0)          */
#define Q_ALP  4   /* divergence-free enrichment α (scalar mode) */

/*
 * Divergence-free enrichment relations:
 *   The divergence-free mode α contributes to face tangential velocities:
 *     Right face (x=xc+h/2): v_t = mya/ρ + α/(2ρ)
 *     Left  face (x=xc-h/2): v_t = mya/ρ - α/(2ρ)   [velocity, not outward]
 *     Top   face (y=yc+h/2): u_t = mxa/ρ - α/(2ρ)
 *     Bottom face(y=yc-h/2): u_t = mxa/ρ + α/(2ρ)
 *   Normal velocities are unaffected: u_n = mxa/ρ, v_n = mya/ρ
 *
 *   In terms of RT0 slope modes: mxs = +α, mys = -α  (div-free constraint)
 */

/* g = 1 ghost-cell layer (sufficient for 1st order stencil) */
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

__device__ __forceinline__ float
sound_speed(float rho, float p)
{
    float cs2 = GAMMA_V * p / rho;
    return sqrtf(cs2 > 1e-14f ? cs2 : 1e-14f);
}

/* Pressure from P0 ρ, E and cell-average momenta mxa, mya */
__device__ __forceinline__ float
pressure_from_cons(float rho, float mxa, float mya, float E)
{
    /* cell-average velocity = mxa/ρ, mya/ρ */
    float u_avg = mxa / rho;
    float v_avg = mya / rho;
    float p = (GAMMA_V - 1.f) * (E - 0.5f * rho * (u_avg*u_avg + v_avg*v_avg));
    return p > 1e-14f ? p : 1e-14f;
}

/* ── Standard HLLC flux with normal n=(nx,ny) ──────────────────────────── */
/* WL,WR: primitives [rho, u, v, p] (u,v are full Cartesian velocities)     */
/* F[4] = [mass, x-mom, y-mom, energy] fluxes in direction n                */
/*                                                                           */
/* Standard HLLC (Toro 1994) with Einfeldt wave speed estimates.            */
__device__ void
hllc_n(const float WL[4], const float WR[4],
       float nx, float ny, float F[4])
{
    float rL=WL[0], uL=WL[1], vL=WL[2], pL=WL[3];
    float rR=WR[0], uR=WR[1], vR=WR[2], pR=WR[3];

    float unL = uL*nx + vL*ny;
    float unR = uR*nx + vR*ny;
    float utL =-uL*ny + vL*nx;
    float utR =-uR*ny + vR*nx;

    float cL = sound_speed(rL, pL), cR = sound_speed(rR, pR);
    float EL = pL/(GAMMA_V-1.f) + 0.5f*rL*(uL*uL + vL*vL);
    float ER = pR/(GAMMA_V-1.f) + 0.5f*rR*(uR*uR + vR*vR);

    /* Einfeldt wave speed estimates */
    float SL = fminf(unL - cL, unR - cR);
    float SR = fmaxf(unL + cL, unR + cR);

    /* Normal-frame fluxes */
    float FnL[4] = { rL*unL,
                     rL*unL*unL + pL,
                     rL*unL*utL,
                     (EL + pL)*unL };
    float FnR[4] = { rR*unR,
                     rR*unR*unR + pR,
                     rR*unR*utR,
                     (ER + pR)*unR };

    if (SL >= 0.f) {
        /* rotate back from normal/tangential to Cartesian */
        F[0] = FnL[0];
        F[1] = FnL[1]*nx - FnL[2]*ny;
        F[2] = FnL[1]*ny + FnL[2]*nx;
        F[3] = FnL[3];
        return;
    }
    if (SR <= 0.f) {
        F[0] = FnR[0];
        F[1] = FnR[1]*nx - FnR[2]*ny;
        F[2] = FnR[1]*ny + FnR[2]*nx;
        F[3] = FnR[3];
        return;
    }

    /* Contact wave speed (Batten et al. 1997) — uses original SL, SR */
    float dL = rL*(SL - unL), dR = rR*(SR - unR);
    float Ss = (pR - pL + dL*unL - dR*unR) / (dL - dR);

    /* Intermediate star states — computed from original SL, SR */
    float factL  = rL*(SL - unL)/(SL - Ss);
    float UsL[4] = { factL,
                     factL * Ss,
                     factL * utL,
                     factL * (EL/rL + (Ss-unL)*(Ss + pL/(rL*(SL-unL)))) };

    float factR  = rR*(SR - unR)/(SR - Ss);
    float UsR[4] = { factR,
                     factR * Ss,
                     factR * utR,
                     factR * (ER/rR + (Ss-unR)*(Ss + pR/(rR*(SR-unR)))) };

    /* Conservative states in normal/tangential frame */
    float UL[4] = { rL, rL*unL, rL*utL, EL };
    float UR[4] = { rR, rR*unR, rR*utR, ER };

    /* Standard HLLC: upwind on either side of contact wave */
    float Fn[4];
    if (Ss >= 0.f)
        for (int q = 0; q < 4; q++)
            Fn[q] = FnL[q] + SL * (UsL[q] - UL[q]);
    else
        for (int q = 0; q < 4; q++)
            Fn[q] = FnR[q] + SR * (UsR[q] - UR[q]);

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
roe_n(const float WL[4], const float WR[4],
      float nx, float ny, float F[4])
{
    float rL=WL[0], uL=WL[1], vL=WL[2], pL=WL[3];
    float rR=WR[0], uR=WR[1], vR=WR[2], pR=WR[3];

    float cL = sound_speed(rL, pL), cR = sound_speed(rR, pR);

    /* Roe averages */
    float sqrL = sqrtf(rL), sqrR = sqrtf(rR);
    float denom = sqrL + sqrR;

    float uRoe  = (sqrL*uL  + sqrR*uR)  / denom;
    float vRoe  = (sqrL*vL  + sqrR*vR)  / denom;
    float HL    = (pL/(GAMMA_V-1.f) + 0.5f*rL*(uL*uL+vL*vL) + pL) / rL;
    float HR    = (pR/(GAMMA_V-1.f) + 0.5f*rR*(uR*uR+vR*vR) + pR) / rR;
    float HRoe  = (sqrL*HL + sqrR*HR) / denom;
    float c2Roe = (GAMMA_V-1.f) * (HRoe - 0.5f*(uRoe*uRoe + vRoe*vRoe));
    float cRoe  = sqrtf(c2Roe > 1e-14f ? c2Roe : 1e-14f);

    /* Rotate to normal/tangential frame */
    float unL = uL*nx + vL*ny,  utL = -uL*ny + vL*nx;
    float unR = uR*nx + vR*ny,  utR = -uR*ny + vR*nx;
    float unRoe = uRoe*nx + vRoe*ny;
    float utRoe = -uRoe*ny + vRoe*nx;

    /* Physical fluxes (normal-frame) */
    float EL = pL/(GAMMA_V-1.f) + 0.5f*rL*(uL*uL+vL*vL);
    float ER = pR/(GAMMA_V-1.f) + 0.5f*rR*(uR*uR+vR*vR);

    float FnL[4] = { rL*unL,        rL*unL*unL + pL, rL*unL*utL, (EL+pL)*unL };
    float FnR[4] = { rR*unR,        rR*unR*unR + pR, rR*unR*utR, (ER+pR)*unR };

    /* Jump in conserved variables (normal frame) */
    float drho = rR - rL;
    float dun  = unR - unL;
    float dut  = utR - utL;
    float dp   = pR  - pL;

    /* Roe-average density: sqrt(rL)*sqrt(rR) */
    float rRoe = sqrL * sqrR;

    /* Wave strengths (eigenvector decomposition of the jump) */
    float b1 = (dp - rRoe*cRoe*dun) / (2.f*c2Roe);  /* left-acoustic   */
    float b2 = drho - dp/c2Roe;                       /* entropy         */
    float b3 = rRoe * dut;                            /* transverse shear */
    float b4 = (dp + rRoe*cRoe*dun) / (2.f*c2Roe);  /* right-acoustic   */

    /* Eigenvalues */
    float lam1 = unRoe - cRoe;
    float lam2 = unRoe;
    float lam3 = unRoe;
    float lam4 = unRoe + cRoe;

    /* Harten-Hyman entropy fix */
    float lam1L = unL - cL, lam1R = unR - cR;
    float lam4L = unL + cL, lam4R = unR + cR;
    float eps1 = fmaxf(0.f, 2.f*(lam1R - lam1L));
    float eps4 = fmaxf(0.f, 2.f*(lam4R - lam4L));
    float al1 = fabsf(lam1); if (al1 < eps1*0.5f) al1 = (lam1*lam1 + eps1*eps1*0.25f) / eps1;
    float al2 = fabsf(lam2);
    float al3 = fabsf(lam3);
    float al4 = fabsf(lam4); if (al4 < eps4*0.5f) al4 = (lam4*lam4 + eps4*eps4*0.25f) / eps4;

    /* Roe dissipation: sum_k al_k * b_k * r_k
     * Right eigenvectors in normal/tangential frame:
     *   r1 = [1, unRoe-cRoe, utRoe, HRoe-unRoe*cRoe]
     *   r2 = [1, unRoe,      utRoe, ½(unRoe²+utRoe²)]
     *   r3 = [0, 0,          1,     utRoe            ]
     *   r4 = [1, unRoe+cRoe, utRoe, HRoe+unRoe*cRoe  ]     */
    float d_rho = al1*b1*1.f           + al2*b2*1.f      + al3*b3*0.f      + al4*b4*1.f;
    float d_un  = al1*b1*(unRoe-cRoe)  + al2*b2*unRoe    + al3*b3*0.f      + al4*b4*(unRoe+cRoe);
    float d_ut  = al1*b1*utRoe         + al2*b2*utRoe     + al3*b3*1.f      + al4*b4*utRoe;
    float d_E   = al1*b1*(HRoe-unRoe*cRoe)
                + al2*b2*0.5f*(unRoe*unRoe+utRoe*utRoe)
                + al3*b3*utRoe
                + al4*b4*(HRoe+unRoe*cRoe);

    float Fn[4];
    Fn[0] = 0.5f*(FnL[0]+FnR[0]) - 0.5f*d_rho;
    Fn[1] = 0.5f*(FnL[1]+FnR[1]) - 0.5f*d_un;
    Fn[2] = 0.5f*(FnL[2]+FnR[2]) - 0.5f*d_ut;
    Fn[3] = 0.5f*(FnL[3]+FnR[3]) - 0.5f*d_E;

    /* Rotate back to Cartesian */
    F[0] = Fn[0];
    F[1] = Fn[1]*nx - Fn[2]*ny;
    F[2] = Fn[1]*ny + Fn[2]*nx;
    F[3] = Fn[3];
}

/* ── Select Riemann solver ───────────────────────────────────────────────── */
/* Default: standard HLLC (Toro 1994).                                      */
/* Define RIEMANN_ROE at compile time (-DRIEMANN_ROE) to use Roe instead.   */
#ifndef RIEMANN_ROE
#  define riemann_n hllc_n
#else
#  define riemann_n roe_n
#endif

/* ═══════════════════════════════════════════════════════════════════════════
 * Apply boundary conditions — periodic wrap
 * Ghost layer (g=1): left ghost ← right interior edge, etc.
 * No sign flips: DOF values copied verbatim (the ring wraps).
 * Works for NVAR=5: (rho, mxa, mya, E, alpha).
 * ═══════════════════════════════════════════════════════════════════════════ */
__global__ void
apply_bc(float * __restrict__ Qp, int N, int Np)
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

/* ═══════════════════════════════════════════════════════════════════════════
 * RHS kernel — Enriched FV (Jung & Perrier, ECCOMAS 2024)
 *
 * DOFs per cell: rho (P0), mxa (P0), mya (P0), E (P0), alpha (div-free)
 *
 * Face velocity reconstruction (from cell (i,j)):
 *   Right face:   u_n = mxa/rho,  v_t = mya/rho + alpha/(2*rho)
 *   Left  face:   u_n = mxa/rho,  v_t = mya/rho - alpha/(2*rho)
 *   Top   face:   v_n = mya/rho,  u_t = mxa/rho - alpha/(2*rho)
 *   Bottom face:  v_n = mya/rho,  u_t = mxa/rho + alpha/(2*rho)
 *
 * P0 updates (rho, mxa, mya, E):
 *   dq/dt = -(1/h)[FR[q] - FL[q] + GT[q] - GB[q]]
 *
 * Div-free enrichment (alpha) via projection onto div-free subspace:
 *   vol_x = (h/rho)*(mxa^2 + alpha^2/3) + p*h
 *   vol_y = (h/rho)*(mya^2 + alpha^2/3) + p*h
 *   dmxs = (3/h^2)*[2*vol_x - (FR[1]+FL[1])*h]
 *   dmys = (3/h^2)*[2*vol_y - (GT[2]+GB[2])*h]
 *   dalpha = (dmxs - dmys) / 2
 * ═══════════════════════════════════════════════════════════════════════════ */
__global__ void
compute_rhs(const float * __restrict__ Qp,
            float       * __restrict__ RHS,
            float       * __restrict__ lam_out,
            int N, int Np, float h)
{
    const int g   = G_GHOST;
    const int N2  = N * N;

    for (int k = blockIdx.x * blockDim.x + threadIdx.x; k < N2;
             k += blockDim.x * gridDim.x) {
        int j = k / N, i = k % N;
        int jp = j + g, ip = i + g;

        /* ── Load center cell DOFs ──────────────────────────────────── */
        float rho   = fmaxf(Qp[IDX_P(Q_RHO, jp, ip, Np)], 1e-14f);
        float mxa   = Qp[IDX_P(Q_MXA, jp, ip, Np)];
        float mya   = Qp[IDX_P(Q_MYA, jp, ip, Np)];
        float E     = Qp[IDX_P(Q_E,   jp, ip, Np)];
        float alpha = Qp[IDX_P(Q_ALP, jp, ip, Np)];
        float p       = pressure_from_cons(rho, mxa, mya, E);
        float inv_rho = 1.f / rho;
        float u_c   = mxa * inv_rho;
        float v_c   = mya * inv_rho;
        float alp2_c = alpha * (0.5f * inv_rho);  /* alpha / (2*rho) */

        /* ── Load ±1 and ±2 stencil primitives (x and y) ───────────── */
        /* Naming: lowercase = negative offset, Uppercase = positive.   *
         * rx1 = cell (i-1,j), Rx1 = cell (i+1,j), etc.                */
#define LPRIM(r_,u_,v_,p_,jj_,ii_) {                                    \
    float _r = fmaxf(Qp[IDX_P(Q_RHO,jj_,ii_,Np)],1e-14f);              \
    float _mx= Qp[IDX_P(Q_MXA,jj_,ii_,Np)];                            \
    float _my= Qp[IDX_P(Q_MYA,jj_,ii_,Np)];                            \
    float _E = Qp[IDX_P(Q_E,  jj_,ii_,Np)];                            \
    (r_)=_r; (u_)=_mx/_r; (v_)=_my/_r;                                 \
    (p_)=pressure_from_cons(_r,_mx,_my,_E); }

        float rx2,ux2,vx2,px2; LPRIM(rx2,ux2,vx2,px2, jp,   ip-2)
        float rx1,ux1,vx1,px1; LPRIM(rx1,ux1,vx1,px1, jp,   ip-1)
        float Rx1,Ux1,Vx1,Px1; LPRIM(Rx1,Ux1,Vx1,Px1, jp,   ip+1)
        float Rx2,Ux2,Vx2,Px2; LPRIM(Rx2,Ux2,Vx2,Px2, jp,   ip+2)
        float ry2,uy2,vy2,py2; LPRIM(ry2,uy2,vy2,py2, jp-2, ip)
        float ry1,uy1,vy1,py1; LPRIM(ry1,uy1,vy1,py1, jp-1, ip)
        float Ry1,Uy1,Vy1,Py1; LPRIM(Ry1,Uy1,Vy1,Py1, jp+1, ip)
        float Ry2,Uy2,Vy2,Py2; LPRIM(Ry2,Uy2,Vy2,Py2, jp+2, ip)
        float alp_xm1 = Qp[IDX_P(Q_ALP, jp,   ip-1, Np)];
        float alp_xp1 = Qp[IDX_P(Q_ALP, jp,   ip+1, Np)];
        float alp_ym1 = Qp[IDX_P(Q_ALP, jp-1, ip,   Np)];
        float alp_yp1 = Qp[IDX_P(Q_ALP, jp+1, ip,   Np)];
#undef LPRIM

        /* 3rd-order linear upwind reconstruction (Lagrange quadratic):
         *   REC_L(a,b,c): left  state at b+½, stencil [b-1=a, b, b+1=c]
         *   REC_R(a,b,c): right state at a+½, stencil [a, a+1=b, a+2=c]
         *   q = (-a + 6b + 3c)/8   and   q = (3a + 6b - c)/8            */
#define REC_L(a,b,c) (0.125f*(-(a)+6.f*(b)+3.f*(c)))
#define REC_R(a,b,c) (0.125f*(3.f*(a)+6.f*(b)-(c)))

        /* ── RIGHT face (+x): left from (i,j), right from (i+1,j) ───── */
        float WL_R[4] = {
            fmaxf(REC_L(rx1,rho,Rx1), 1e-14f),
            REC_L(ux1, u_c, Ux1),
            REC_L(vx1, v_c, Vx1) + alp2_c,
            fmaxf(REC_L(px1,  p, Px1), 1e-14f) };
        float WR_R[4] = {
            fmaxf(REC_R(rho,Rx1,Rx2), 1e-14f),
            REC_R(u_c, Ux1, Ux2),
            REC_R(v_c, Vx1, Vx2) - alp_xp1 * (0.5f / Rx1),
            fmaxf(REC_R(  p, Px1, Px2), 1e-14f) };
        float FR[4]; riemann_n(WL_R, WR_R, 1.f, 0.f, FR);

        /* ── LEFT face (+x): left from (i-1,j), right from (i,j) ────── */
        float WL_L[4] = {
            fmaxf(REC_L(rx2,rx1,rho), 1e-14f),
            REC_L(ux2, ux1, u_c),
            REC_L(vx2, vx1, v_c) + alp_xm1 * (0.5f / rx1),
            fmaxf(REC_L(px2, px1,   p), 1e-14f) };
        float WR_L[4] = {
            fmaxf(REC_R(rx1,rho,Rx1), 1e-14f),
            REC_R(ux1, u_c, Ux1),
            REC_R(vx1, v_c, Vx1) - alp2_c,
            fmaxf(REC_R(px1,   p, Px1), 1e-14f) };
        float FL[4]; riemann_n(WL_L, WR_L, 1.f, 0.f, FL);

        /* ── TOP face (+y): left from (i,j), right from (i,j+1) ─────── */
        float WL_T[4] = {
            fmaxf(REC_L(ry1,rho,Ry1), 1e-14f),
            REC_L(uy1, u_c, Uy1) - alp2_c,
            REC_L(vy1, v_c, Vy1),
            fmaxf(REC_L(py1,  p, Py1), 1e-14f) };
        float WR_T[4] = {
            fmaxf(REC_R(rho,Ry1,Ry2), 1e-14f),
            REC_R(u_c, Uy1, Uy2) + alp_yp1 * (0.5f / Ry1),
            REC_R(v_c, Vy1, Vy2),
            fmaxf(REC_R(  p, Py1, Py2), 1e-14f) };
        float GT[4]; riemann_n(WL_T, WR_T, 0.f, 1.f, GT);

        /* ── BOTTOM face (+y): left from (i,j-1), right from (i,j) ──── */
        float WL_B[4] = {
            fmaxf(REC_L(ry2,ry1,rho), 1e-14f),
            REC_L(uy2, uy1, u_c) - alp_ym1 * (0.5f / ry1),
            REC_L(vy2, vy1, v_c),
            fmaxf(REC_L(py2, py1,   p), 1e-14f) };
        float WR_B[4] = {
            fmaxf(REC_R(ry1,rho,Ry1), 1e-14f),
            REC_R(uy1, u_c, Uy1) + alp2_c,
            REC_R(vy1, v_c, Vy1),
            fmaxf(REC_R(py1,   p, Py1), 1e-14f) };
        float GB[4]; riemann_n(WL_B, WR_B, 0.f, 1.f, GB);

#undef REC_L
#undef REC_R

        /* ── P0 RHS for rho, mxa, mya, E ─────────────────────────────── */
        float inv_h  = 1.f / h;
        float inv_h2 = inv_h * inv_h;
        float drho = -inv_h * (FR[0] - FL[0] + GT[0] - GB[0]);
        float dmxa = -inv_h * (FR[1] - FL[1] + GT[1] - GB[1]);
        float dmya = -inv_h * (FR[2] - FL[2] + GT[2] - GB[2]);
        float dE   = -inv_h * (FR[3] - FL[3] + GT[3] - GB[3]);

        /* ── Div-free enrichment RHS (alpha) ──────────────────────────── */
        /* Alpha is the coefficient of the enrichment mode phi_alpha =
         * ((x-xc)/h, -(y-yc)/h)^T in the velocity approximation space S0.
         * Mass matrix: M_alpha = h^2/6.
         * DG weak form:
         *   M_alpha * d(alpha)/dt = cell_vol_integral - face_integral
         *
         * Cell volume integral (4-pt Gauss, pressure cancels analytically):
         *   = (h/rho) * (mxa^2 - mya^2)
         *
         * Face integral (phi_alpha evaluated at face midpoints):
         *   = (h/2) * (FR[1] + FL[1] - GT[2] - GB[2])
         *
         * After dividing by M_alpha = h^2/6:
         *   d(alpha)/dt = (6/(rho*h^2))*(mxa^2-mya^2) - (3/h)*(FR[1]+FL[1]-GT[2]-GB[2])
         */
        float cell_alp = 3.f * inv_h * inv_rho * (mxa*mxa - mya*mya);
        float face_alp = 1.5f * inv_h * (FR[1] + FL[1] - GT[2] - GB[2]);
        float dalpha   = cell_alp - face_alp;
        (void)dalpha;
        dalpha = 0.f; /* DEBUG: disable alpha update */

        /* ── Write RHS ────────────────────────────────────────────────── */
        RHS[Q_RHO * N2 + k] = drho;
        RHS[Q_MXA * N2 + k] = dmxa;
        RHS[Q_MYA * N2 + k] = dmya;
        RHS[Q_E   * N2 + k] = dE;
        RHS[Q_ALP * N2 + k] = dalpha;

        /* ── CFL spectral radius ─────────────────────────────────────── */
        float cs = sound_speed(rho, p);
        lam_out[k] = fabsf(u_c) + fabsf(v_c) + cs;
    }
}

/* ═══════════════════════════════════════════════════════════════════════════
 * IC kernel  — circular Sod shock tube
 * ═══════════════════════════════════════════════════════════════════════════ */
__global__ void
ic_kernel(float * __restrict__ Qp, int N, int Np, float h)
{
    const int g   = G_GHOST;
    const int N2  = N * N;

    for (int k = blockIdx.x * blockDim.x + threadIdx.x; k < N2;
             k += blockDim.x * gridDim.x) {
        int j = k / N, i = k % N;
        int jp = j + g, ip = i + g;

        float xc = (i + 0.5f) * h, yc = (j + 0.5f) * h;
        float r  = sqrtf((xc - 0.5f)*(xc - 0.5f) + (yc - 0.5f)*(yc - 0.5f));

        /* Smooth tanh interface (width ~ h) */
        float delta = 1.5f * h;
        float phi   = 0.5f * (1.f + tanhf((0.25f - r) / delta));

        float rho = 0.125f + 9.875f * phi;   /* 1.0 inside, 0.125 outside */
        float p   = 0.1f   + 9.9f   * phi;   /* 1.0 inside, 0.1   outside */
        float E   = p / (GAMMA_V - 1.f);     /* zero velocity IC */

        /* Zero velocity → alpha = 0 */
        Qp[IDX_P(Q_RHO, jp, ip, Np)] = rho;
        Qp[IDX_P(Q_MXA, jp, ip, Np)] = 0.f;
        Qp[IDX_P(Q_MYA, jp, ip, Np)] = 0.f;
        Qp[IDX_P(Q_E,   jp, ip, Np)] = E;
        Qp[IDX_P(Q_ALP, jp, ip, Np)] = 0.f;
    }
}

/* ═══════════════════════════════════════════════════════════════════════════
 * SSP-RK3 stage kernels — operate directly on raw DOFs (no prim conversion)
 * ═══════════════════════════════════════════════════════════════════════════ */
__global__ void
rk3_s1(float * __restrict__ U1p,
       const float * __restrict__ U0p,
       const float * __restrict__ L,
       float dt, int N, int Np)
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
        U1p[Q_RHO * Np2 + idxp] = fmaxf(U1p[Q_RHO * Np2 + idxp], 1e-14f);
        U1p[Q_E   * Np2 + idxp] = fmaxf(U1p[Q_E   * Np2 + idxp], 1e-14f);
    }
}

__global__ void
rk3_s2(float * __restrict__ U2p,
       const float * __restrict__ U0p,
       const float * __restrict__ U1p,
       const float * __restrict__ L,
       float dt, int N, int Np)
{
    const int g   = G_GHOST;
    const int N2  = N * N;
    const int Np2 = Np * Np;
    for (int k = blockIdx.x * blockDim.x + threadIdx.x; k < N2;
             k += blockDim.x * gridDim.x) {
        int j = k / N, i = k % N;
        int idxp = (j + g) * Np + (i + g);
        for (int q = 0; q < NVAR; q++)
            U2p[q * Np2 + idxp] = 0.75f * U0p[q * Np2 + idxp]
                                 + 0.25f * (U1p[q * Np2 + idxp] + dt * L[q * N2 + k]);
        U2p[Q_RHO * Np2 + idxp] = fmaxf(U2p[Q_RHO * Np2 + idxp], 1e-14f);
        U2p[Q_E   * Np2 + idxp] = fmaxf(U2p[Q_E   * Np2 + idxp], 1e-14f);
    }
}

__global__ void
rk3_s3(float * __restrict__ Up,
       const float * __restrict__ U0p,
       const float * __restrict__ U2p,
       const float * __restrict__ L,
       float dt, int N, int Np)
{
    const int g   = G_GHOST;
    const int N2  = N * N;
    const int Np2 = Np * Np;
    for (int k = blockIdx.x * blockDim.x + threadIdx.x; k < N2;
             k += blockDim.x * gridDim.x) {
        int j = k / N, i = k % N;
        int idxp = (j + g) * Np + (i + g);
        for (int q = 0; q < NVAR; q++)
            Up[q * Np2 + idxp] = (1.f/3.f) * U0p[q * Np2 + idxp]
                                + (2.f/3.f) * (U2p[q * Np2 + idxp] + dt * L[q * N2 + k]);
        Up[Q_RHO * Np2 + idxp] = fmaxf(Up[Q_RHO * Np2 + idxp], 1e-14f);
        Up[Q_E   * Np2 + idxp] = fmaxf(Up[Q_E   * Np2 + idxp], 1e-14f);
    }
}

/* ═══════════════════════════════════════════════════════════════════════════
 * max reduction
 * ═══════════════════════════════════════════════════════════════════════════ */
__global__ void
reduce_max(const float * __restrict__ in, float * __restrict__ out, int n)
{
    extern __shared__ float sm[];
    int tid = threadIdx.x;
    float v = 0.f;
    for (int k = blockIdx.x * blockDim.x + tid; k < n; k += blockDim.x * gridDim.x)
        v = fmaxf(v, in[k]);
    sm[tid] = v; __syncthreads();
    for (int s = BLOCK1D / 2; s > 0; s >>= 1) {
        if (tid < s) sm[tid] = fmaxf(sm[tid], sm[tid + s]);
        __syncthreads();
    }
    if (tid == 0) out[blockIdx.x] = sm[0];
}

static float
gpu_max(const float *d_in, float *d_tmp, int n)
{
    reduce_max<<<GS_NBLK, BLOCK1D, BLOCK1D * sizeof(float)>>>(d_in, d_tmp, n);
    reduce_max<<<1,        BLOCK1D, BLOCK1D * sizeof(float)>>>(d_tmp, d_tmp, GS_NBLK);
    float v;
    CK(cudaMemcpy(&v, d_tmp, sizeof(float), cudaMemcpyDeviceToHost));
    return v;
}

/* ── min reduction (mirrors reduce_max) ──────────────────────────────────── */
__global__ void
reduce_min(const float * __restrict__ in, float * __restrict__ out, int n)
{
    extern __shared__ float sm[];
    int tid = threadIdx.x;
    float v = 1e38f;
    for (int k = blockIdx.x * blockDim.x + tid; k < n; k += blockDim.x * gridDim.x)
        v = fminf(v, in[k]);
    sm[tid] = v; __syncthreads();
    for (int s = BLOCK1D / 2; s > 0; s >>= 1) {
        if (tid < s) sm[tid] = fminf(sm[tid], sm[tid + s]);
        __syncthreads();
    }
    if (tid == 0) out[blockIdx.x] = sm[0];
}

static float
gpu_min(const float *d_in, float *d_tmp, int n)
{
    reduce_min<<<GS_NBLK, BLOCK1D, BLOCK1D * sizeof(float)>>>(d_in, d_tmp, n);
    reduce_min<<<1,        BLOCK1D, BLOCK1D * sizeof(float)>>>(d_tmp, d_tmp, GS_NBLK);
    float v;
    CK(cudaMemcpy(&v, d_tmp, sizeof(float), cudaMemcpyDeviceToHost));
    return v;
}

/* ═══════════════════════════════════════════════════════════════════════════
 * Download helper — extract density from padded array for output
 * ═══════════════════════════════════════════════════════════════════════════ */
__global__ void
extract_rho_p(const float * __restrict__ Qp,
              float * __restrict__ out_rho,
              float * __restrict__ out_p,
              int N, int Np)
{
    const int g  = G_GHOST;
    const int N2 = N * N;
    for (int k = blockIdx.x * blockDim.x + threadIdx.x; k < N2;
             k += blockDim.x * gridDim.x) {
        int j = k / N, i = k % N;
        int jp = j + g, ip = i + g;
        float rho = Qp[IDX_P(Q_RHO, jp, ip, Np)];
        float mxa = Qp[IDX_P(Q_MXA, jp, ip, Np)];
        float mya = Qp[IDX_P(Q_MYA, jp, ip, Np)];
        float E   = Qp[IDX_P(Q_E,   jp, ip, Np)];
        rho = fmaxf(rho, 1e-14f);
        out_rho[k] = rho;
        out_p  [k] = pressure_from_cons(rho, mxa, mya, E);
    }
}

/* Extract density (left panel) and cell-centre v-velocity (right panel)    *
 * for Kelvin-Helmholtz visualisation.                                       *
 *   mya = cell-average y-momentum  →  v_avg = mya/rho                      */
__global__ void
extract_rho_v(const float * __restrict__ Qp,
              float * __restrict__ out_rho,
              float * __restrict__ out_v,
              int N, int Np)
{
    const int g  = G_GHOST;
    const int N2 = N * N;
    for (int k = blockIdx.x * blockDim.x + threadIdx.x; k < N2;
             k += blockDim.x * gridDim.x) {
        int j = k / N, i = k % N;
        int jp = j + g, ip = i + g;
        float rho = fmaxf(Qp[IDX_P(Q_RHO, jp, ip, Np)], 1e-14f);
        float mya = Qp[IDX_P(Q_MYA, jp, ip, Np)];
        out_rho[k] = rho;
        out_v  [k] = mya / rho;  /* cell-average v */
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
#define VORTEX_R  0.2f

__device__ void
vortex_exact(float x, float y, float Ma,
             float *rho_out, float *u_out, float *v_out, float *p_out)
{
    const float g     = GAMMA_V;
    const float rhoI  = 1.f;
    const float pI    = 1.f / g;            /* c_inf = 1 */
    const float TI    = pI / rhoI;          /* = 1/g     */
    const float beta  = 2.f * (float)M_PI * Ma;
    const float R     = VORTEX_R;

    float dx = x - 0.5f, dy = y - 0.5f;
    float r2 = (dx*dx + dy*dy) / (R*R);
    float f  = (beta / (2.f*(float)M_PI)) * expf(0.5f*(1.f - r2));

    float u  = -f * dy / R;
    float v  = +f * dx / R;

    float dT  = -(g-1.f) * beta*beta / (8.f*g*(float)M_PI*(float)M_PI)
                * expf(1.f - r2);
    float T   = TI + dT;
    if (T < 1e-6f) T = 1e-6f;
    float Tr  = T / TI;                     /* T / T_inf */
    float rho = rhoI * powf(Tr, 1.f/(g-1.f));
    float p   = pI   * powf(Tr, g  /(g-1.f));

    *rho_out = rho;
    *u_out   = u;
    *v_out   = v;
    *p_out   = p;
}

/* IC kernel — isentropic vortex projected onto enriched FV DOFs */
__global__ void
ic_vortex(float * __restrict__ Qp, int N, int Np, float h, float Ma)
{
    const int g  = G_GHOST;
    const int N2 = N * N;
    for (int k = blockIdx.x * blockDim.x + threadIdx.x; k < N2;
             k += blockDim.x * gridDim.x) {
        int j = k / N, i = k % N;
        int jp = j + g, ip = i + g;

        /* cell centre */
        float xc = (i + 0.5f) * h, yc = (j + 0.5f) * h;
        float rho, uc, vc, p;
        vortex_exact(xc, yc, Ma, &rho, &uc, &vc, &p);
        float E = p / (GAMMA_V - 1.f) + 0.5f * rho * (uc*uc + vc*vc);

        /* Divergence-free enrichment alpha from div-free part of velocity:
         * The div-free mode gives tangential velocity offsets alpha/(2*rho) at each face.
         * Using the curl of rho*u: alpha ~ rho*h*(du/dy or -dv/dx) via 4 face evaluations:
         *   alpha = rho * h^2 * (d(rho*v)/dx - d(rho*u)/dy) / (2) -- stream function
         * Approximate: compute exact rho*u at all 4 face midpoints.
         */
        float rR, uR, vR_d, pR; vortex_exact((i+1)*h, yc,     Ma, &rR, &uR, &vR_d, &pR);
        float rL, uL, vL_d, pL; vortex_exact(i*h,     yc,     Ma, &rL, &uL, &vL_d, &pL);
        float rT, uT, vT,   pT; vortex_exact(xc,     (j+1)*h, Ma, &rT, &uT, &vT,   &pT);
        float rB, uB, vB,   pB; vortex_exact(xc,      j*h,    Ma, &rB, &uB, &vB,   &pB);

        /* P0 cell-average momenta: use cell-centre values */
        float mxa_c = rho * uc;
        float mya_c = rho * vc;

        /* Divergence-free alpha: project face tangential velocity offsets
         * alpha = rho*h^2 * curl(u) / 2  (stream function perturbation)
         * Discrete: alpha ~ rho * h * [(rT*uT/rT - rB*uB/rB) - (rR*vR/rR - rL*vL/rL)] / 4
         *         = (rho * h / 4) * [(uT - uB) - (vR - vL)] */
        float alpha = (rho * h * 0.25f) * ((uT - uB) - (vR_d - vL_d));

        Qp[IDX_P(Q_RHO, jp, ip, Np)] = rho;
        Qp[IDX_P(Q_MXA, jp, ip, Np)] = mxa_c;
        Qp[IDX_P(Q_MYA, jp, ip, Np)] = mya_c;
        Qp[IDX_P(Q_E,   jp, ip, Np)] = E;
        Qp[IDX_P(Q_ALP, jp, ip, Np)] = alpha;
        (void)pL; (void)pT; (void)pB; (void)pR;
    }
}

/* ── L2 error vs. exact stationary vortex ──────────────────────────────── */
/* Writes per-cell squared error  (rho-rho_ex)²  and  (p-p_ex)²            */
__global__ void
compute_l2_err(const float * __restrict__ Qp,
               float * __restrict__ err_rho,
               float * __restrict__ err_p,
               int N, int Np, float h, float Ma)
{
    const int g  = G_GHOST;
    const int N2 = N * N;
    for (int k = blockIdx.x * blockDim.x + threadIdx.x; k < N2;
             k += blockDim.x * gridDim.x) {
        int j = k / N, i = k % N;
        int jp = j + g, ip = i + g;

        float rho = fmaxf(Qp[IDX_P(Q_RHO, jp, ip, Np)], 1e-14f);
        float mxa = Qp[IDX_P(Q_MXA, jp, ip, Np)];
        float mya = Qp[IDX_P(Q_MYA, jp, ip, Np)];
        float E   = Qp[IDX_P(Q_E,   jp, ip, Np)];
        float p   = pressure_from_cons(rho, mxa, mya, E);

        float xc = (i + 0.5f) * h, yc = (j + 0.5f) * h;
        float re, ue, ve, pe;
        vortex_exact(xc, yc, Ma, &re, &ue, &ve, &pe);

        float dr = rho - re, dp = p - pe;
        err_rho[k] = dr * dr * h * h;   /* volume-weighted */
        err_p  [k] = dp * dp * h * h;
    }
}

/* ── sum reduction ─────────────────────────────────────────────────────── */
__global__ void
reduce_sum(const float * __restrict__ in, float * __restrict__ out, int n)
{
    extern __shared__ float sm[];
    int tid = threadIdx.x;
    float v = 0.f;
    for (int k = blockIdx.x * blockDim.x + tid; k < n; k += blockDim.x * gridDim.x)
        v += in[k];
    sm[tid] = v; __syncthreads();
    for (int s = BLOCK1D / 2; s > 0; s >>= 1) {
        if (tid < s) sm[tid] += sm[tid + s];
        __syncthreads();
    }
    if (tid == 0) out[blockIdx.x] = sm[0];
}

static float
gpu_sum(const float *d_in, float *d_tmp, int n)
{
    reduce_sum<<<GS_NBLK, BLOCK1D, BLOCK1D * sizeof(float)>>>(d_in, d_tmp, n);
    reduce_sum<<<1,        BLOCK1D, BLOCK1D * sizeof(float)>>>(d_tmp, d_tmp, GS_NBLK);
    float v;
    CK(cudaMemcpy(&v, d_tmp, sizeof(float), cudaMemcpyDeviceToHost));
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
vortex_paper_exact(float x, float y, float u0, float v0,
                   float *rho_out, float *u_out, float *v_out, float *p_out)
{
    const float g  = GAMMA_V;
    const float eps = 5.f;
    const float pi  = (float)M_PI;

    float dx = x - 5.f, dy = y - 5.f;
    float r2 = dx*dx + dy*dy;
    float f  = eps / (2.f*pi) * expf(0.5f*(1.f - r2));

    float u  = u0 - f * dy;
    float v  = v0 + f * dx;

    float dT = -(g-1.f)*eps*eps / (8.f*g*pi*pi) * expf(1.f - r2);
    float T  = 1.f + dT;
    if (T < 1e-6f) T = 1e-6f;

    *rho_out = powf(T, 1.f/(g-1.f));
    *p_out   = powf(T, g  /(g-1.f));
    *u_out   = u;
    *v_out   = v;
}

/* IC kernel — paper isentropic vortex on [0,10]^2, h = 10/N */
__global__ void
ic_vortex_paper(float * __restrict__ Qp, int N, int Np,
                float h, float u0, float v0)
{
    const int g  = G_GHOST;
    const int N2 = N * N;
    for (int k = blockIdx.x * blockDim.x + threadIdx.x; k < N2;
             k += blockDim.x * gridDim.x) {
        int j = k / N, i = k % N;
        int jp = j + g, ip = i + g;

        float xc = (i + 0.5f) * h, yc = (j + 0.5f) * h;
        float rho, uc, vc, p;
        vortex_paper_exact(xc, yc, u0, v0, &rho, &uc, &vc, &p);
        float E = p / (GAMMA_V - 1.f) + 0.5f * rho * (uc*uc + vc*vc);

        /* Face midpoint values for div-free enrichment projection */
        float rR, uR, vR_d, pR;  vortex_paper_exact((i+1)*h, yc,      u0, v0, &rR, &uR, &vR_d, &pR);
        float rL, uL, vL_d, pL;  vortex_paper_exact(i*h,     yc,      u0, v0, &rL, &uL, &vL_d, &pL);
        float rT, uT, vT,   pT;  vortex_paper_exact(xc,     (j+1)*h,  u0, v0, &rT, &uT, &vT,   &pT);
        float rB, uB, vB,   pB;  vortex_paper_exact(xc,      j*h,     u0, v0, &rB, &uB, &vB,   &pB);

        /* Div-free alpha: discrete curl of velocity field */
        float alpha = (rho * h * 0.25f) * ((uT - uB) - (vR_d - vL_d));

        Qp[IDX_P(Q_RHO, jp, ip, Np)] = rho;
        Qp[IDX_P(Q_MXA, jp, ip, Np)] = rho * uc;
        Qp[IDX_P(Q_MYA, jp, ip, Np)] = rho * vc;
        Qp[IDX_P(Q_E,   jp, ip, Np)] = E;
        Qp[IDX_P(Q_ALP, jp, ip, Np)] = alpha;
        (void)pL; (void)pT; (void)pB; (void)pR;
        (void)rR; (void)rL; (void)rT; (void)rB;
    }
}

/* L2 error vs. exact paper vortex
 * Computes errors in rho, rho*u, rho*v, rho*E  (matches Tables 2 & 3 of paper).
 * For the stationary case (u0=v0=0, tf=1) the exact solution equals the IC.
 * For the moving case (u0=v0=1, tf=10) the vortex traverses exactly one period
 * on the [0,10]^2 periodic domain, so the exact solution also equals the IC.
 */
__global__ void
compute_l2_err_paper(const float * __restrict__ Qp,
                     float * __restrict__ err_rho,
                     float * __restrict__ err_ru,
                     float * __restrict__ err_rv,
                     float * __restrict__ err_re,
                     int N, int Np, float h, float u0, float v0)
{
    const int g  = G_GHOST;
    const int N2 = N * N;
    for (int k = blockIdx.x * blockDim.x + threadIdx.x; k < N2;
             k += blockDim.x * gridDim.x) {
        int j = k / N, i = k % N;
        int jp = j + g, ip = i + g;

        float rho = fmaxf(Qp[IDX_P(Q_RHO, jp, ip, Np)], 1e-14f);
        float mxa = Qp[IDX_P(Q_MXA, jp, ip, Np)];
        float mya = Qp[IDX_P(Q_MYA, jp, ip, Np)];
        float E   = Qp[IDX_P(Q_E,   jp, ip, Np)];
        /* cell-average momenta = mxa, mya directly (modal mode 0) */
        float ru  = mxa;
        float rv  = mya;
        /* total energy density stored as E = rho*E_spec */
        float rE  = E;

        float xc = (i + 0.5f) * h, yc = (j + 0.5f) * h;
        float re, ue, ve, pe;
        vortex_paper_exact(xc, yc, u0, v0, &re, &ue, &ve, &pe);
        float rue = re * ue;
        float rve = re * ve;
        float rEe = pe / (GAMMA_V - 1.f) + 0.5f * re * (ue*ue + ve*ve);

        float vol = h * h;
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

__device__ float
kh_H(float y)   /* y = shifted coordinate ∈ [-0.5, 0.5] */
{
    const float omega = 1.f / 16.f;
    const float y1 = -0.25f, y2 = 0.25f;
    const float pi  = (float)M_PI;

    /* Transition at y1: +1 → -1 */
    if (y >= y1 - omega*0.5f && y < y1 + omega*0.5f)
        return -sinf(pi * (y - y1) / omega);
    /* Central region */
    if (y >= y1 + omega*0.5f && y < y2 - omega*0.5f)
        return -1.f;
    /* Transition at y2: -1 → +1 */
    if (y >= y2 - omega*0.5f && y < y2 + omega*0.5f)
        return  sinf(pi * (y - y2) / omega);
    /* Outer region */
    return 1.f;
}

__global__ void
ic_kh(float * __restrict__ Qp, int N, int Np, float h, float Mkh)
{
    const int g  = G_GHOST;
    const int N2 = N * N;
    const float rkh  = Mkh * 0.1f;  /* density jump ∝ Ma */
    const float del  = 0.1f;
    const float pi   = (float)M_PI;

    for (int k = blockIdx.x * blockDim.x + threadIdx.x; k < N2;
             k += blockDim.x * gridDim.x) {
        int j = k / N, i = k % N;
        int jp = j + g, ip = i + g;

        float xc = (i + 0.5f) * h;
        float yc = (j + 0.5f) * h;
        float ys = yc - 1.0f;  /* centre domain at y=1 → y_shifted ∈ [-1,1] */

        /* Cell-center primitives */
        float Hc  = kh_H(ys);
        float rho = GAMMA_V + Hc * rkh;
        float u   = Mkh * Hc;
        float v   = del * Mkh * sinf(2.f * pi * xc);
        float p   = 1.f;
        float E   = p / (GAMMA_V - 1.f) + 0.5f * rho * (u*u + v*v);

        /* y-faces: ρ varies with y */
        float ys_t = (j+1)*h - 1.0f;
        float ys_b = j*h     - 1.0f;
        float rhoT = GAMMA_V + kh_H(ys_t) * rkh;
        float rhoB = GAMMA_V + kh_H(ys_b) * rkh;
        float vT   = del * Mkh * sinf(2.f * pi * xc);
        float vB   = del * Mkh * sinf(2.f * pi * xc);

        /* For KH, initialize alpha=0 and let it evolve. */
        float alpha_kh = 0.f;

        Qp[IDX_P(Q_RHO, jp, ip, Np)] = rho;
        Qp[IDX_P(Q_MXA, jp, ip, Np)] = rho * u;
        Qp[IDX_P(Q_MYA, jp, ip, Np)] = rho * v;
        Qp[IDX_P(Q_E,   jp, ip, Np)] = E;
        Qp[IDX_P(Q_ALP, jp, ip, Np)] = alpha_kh;
        (void)rhoT; (void)rhoB; (void)vT; (void)vB;
    }
}

/* ═══════════════════════════════════════════════════════════════════════════
 * PPM output
 * ═══════════════════════════════════════════════════════════════════════════ */
static const unsigned char PLASMA[256][3] = {
{13,8,135},{16,7,136},{19,7,137},{22,7,138},{25,6,140},{27,6,141},
{29,6,142},{32,6,143},{34,5,144},{36,5,145},{38,5,145},{40,4,146},
{42,4,147},{44,4,148},{46,4,149},{47,3,150},{49,3,151},{51,3,152},
{53,3,153},{55,2,154},{57,2,155},{58,2,156},{60,2,157},{62,1,158},
{64,1,159},{65,1,160},{67,1,161},{69,0,162},{70,0,163},{72,0,164},
{74,0,165},{75,0,166},{77,0,167},{79,0,168},{80,0,168},{82,0,169},
{84,0,170},{85,0,171},{87,0,172},{89,0,173},{90,1,173},{92,1,174},
{94,1,175},{95,2,176},{97,2,176},{99,3,177},{100,4,178},{102,4,178},
{104,5,179},{105,6,180},{107,7,180},{109,7,181},{110,8,182},{112,9,182},
{114,10,183},{115,11,183},{117,12,184},{119,13,184},{120,14,185},
{122,15,185},{124,16,186},{125,17,186},{127,18,187},{129,19,187},
{130,20,188},{132,21,188},{134,22,189},{135,23,189},{137,25,190},
{139,26,190},{141,27,190},{142,28,191},{144,29,191},{146,30,191},
{147,32,192},{149,33,192},{151,34,192},{153,35,193},{154,36,193},
{156,38,193},{158,39,193},{159,40,194},{161,42,194},{163,43,194},
{165,44,194},{166,45,195},{168,47,195},{170,48,195},{171,49,195},
{173,51,195},{175,52,196},{176,54,196},{178,55,196},{180,56,196},
{182,58,196},{183,59,196},{185,61,196},{187,62,196},{188,64,196},
{190,66,196},{192,67,196},{194,69,196},{195,70,196},{197,72,196},
{199,74,196},{200,75,196},{202,77,196},{204,79,196},{205,80,196},
{207,82,196},{209,84,195},{210,85,195},{212,87,195},{214,89,195},
{215,91,195},{217,93,194},{219,94,194},{220,96,194},{222,98,193},
{224,100,193},{225,101,193},{227,103,192},{229,105,192},{230,107,191},
{232,109,191},{233,111,190},{235,113,190},{237,114,189},{238,116,189},
{240,118,188},{241,120,187},{243,122,187},{244,124,186},{246,126,185},
{247,128,185},{249,130,184},{250,132,183},{251,134,182},{253,136,181},
{254,138,181},{254,140,180},{254,142,179},{254,144,178},{254,146,177},
{254,148,175},{254,151,174},{254,153,173},{254,155,172},{254,157,171},
{254,159,169},{254,161,168},{254,163,167},{254,165,165},{254,167,164},
{254,169,163},{254,171,161},{254,173,160},{254,175,158},{254,177,157},
{254,179,155},{254,181,154},{254,183,152},{254,185,151},{254,187,149},
{254,189,148},{254,191,146},{254,193,144},{254,195,143},{254,197,141},
{254,199,139},{254,201,138},{254,203,136},{254,206,134},{254,208,132},
{254,210,131},{254,212,129},{254,214,127},{254,216,125},{254,218,124},
{254,220,122},{254,222,120},{254,224,118},{254,226,116},{254,228,115},
{254,230,113},{254,232,111},{254,234,109},{254,236,107},{254,238,105},
{254,240,103},{254,242,101},{254,244,99},{253,246,98},{253,248,96},
{253,250,94},{253,252,92},{252,254,90},{252,254,89},{252,254,87},
{252,254,85},{252,254,83},{251,254,81},{251,254,79},{251,254,77},
{250,254,75},{250,254,73},{249,254,71},{249,253,69},{248,253,68},
{248,253,66},{247,253,64},{247,253,62},{246,253,60},{246,252,58},
{245,252,56},{244,252,54},{244,252,52},{243,252,50},{242,252,48},
{242,251,46},{241,251,44},{240,251,42},{240,251,40},{239,250,38},
{238,250,36},{237,250,34},{237,249,32},{236,249,30},{235,249,28},
{234,248,26},{233,248,24},{232,247,22},{231,247,20},{231,246,18},
{230,246,16},{229,245,14},{228,245,12},{227,244,10},{226,244,8},
{225,243,6},{224,243,4},{223,242,2},{222,242,0}
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

static void
write_ppm_2panel(const char *fname,
                 const float *tl, float tl_min, float tl_max,
                 const float *tr, float tr_min, float tr_max,
                 int N)
{
    FILE *f = fopen(fname, "wb");
    if (!f) { fprintf(stderr, "Cannot open %s\n", fname); return; }
    fprintf(f, "P6\n%d %d\n255\n", 2*N, N);

    for (int r = 0; r < N; r++) {
        int phys_j = N - 1 - r;   /* flip y */
        for (int c = 0; c < 2*N; c++) {
            const float *panel = (c < N) ? tl : tr;
            float vmin = (c < N) ? tl_min : tr_min;
            float vmax = (c < N) ? tl_max : tr_max;
            int   ci   = (c < N) ? c : c - N;
            const unsigned char (*cmap)[3] = (c < N) ? PLASMA : VIRIDIS;

            float t = (panel[phys_j * N + ci] - vmin) / (vmax - vmin + 1e-30f);
            if (!(t > 0.f)) t = 0.f;
            if (t > 1.f)    t = 1.f;
            int idx = (int)(t * 255.f);
            if (idx < 0) idx = 0; if (idx > 255) idx = 255;
            unsigned char pix[3] = { cmap[idx][0], cmap[idx][1], cmap[idx][2] };
            fwrite(pix, 1, 3, f);
        }
    }
    fclose(f);
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
 * ═══════════════════════════════════════════════════════════════════════════ */
int main(int argc, char **argv)
{
    int   N    = (argc > 1) ? atoi(argv[1]) : 400;

    /* Determine test mode */
    const char *mode_str = (argc > 2) ? argv[2] : "";
    float Ma = 0.f;
    int do_vortex       = 0;
    int do_paper_vortex = 0;  /* pv: stationary, pvmove: moving */
    int do_kh           = 0;
    int pv_moving       = 0;  /* 1 → pvmove: u0=v0=1, tf=10 */
    float kh_mach       = 0.01f; /* KH convective Mach number */

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
            float maybe = strtof(argv[3], &endp);
            if (endp != argv[3] && *endp == '\0' && maybe > 0.f)
                kh_mach = maybe;
        }
    } else {
        /* Legacy: numeric Ma for original vortex */
        char *endp;
        float maybe_ma = strtof(mode_str, &endp);
        if (endp != mode_str && *endp == '\0' && maybe_ma > 0.f) {
            Ma = maybe_ma;
            do_vortex = 1;
        }
    }

    /* Domain size and run time per mode */
    float L_domain = 1.f;           /* default: [0,1]^2 */
    if (do_paper_vortex) L_domain = 10.f;
    if (do_kh)           L_domain =  2.f;

    float CFL   = 0.30f;
    float t_end;
    if (do_vortex) {
        /* One full vortex orbital period: T = 2πR / |u|_max = 2πR / Ma
         * (peak velocity = Ma·c∞ = Ma since c∞=1; R = VORTEX_R = 0.2)
         * This gives the same number of vortex rotations regardless of Ma. */
        t_end = 2.f * (float)M_PI * VORTEX_R / Ma;
    } else if (do_paper_vortex) t_end = pv_moving ? 10.0f : 1.0f;
    /* t_end for KH: 8 convective times = 8 * L_domain / (2*kh_mach)  */
    else if (do_kh)      t_end = 8.f * L_domain / (2.f * kh_mach);
    else                 t_end = 1.0f;   /* Sod */

    float h    = L_domain / N;
    const int g = G_GHOST;
    int Np      = N + 2 * g;

    /* Background velocity for paper vortex */
    float u0 = pv_moving ? 1.f : 0.f;
    float v0 = pv_moving ? 1.f : 0.f;

    /* GPU info */
    int dev; cudaDeviceProp prop;
    CK(cudaGetDevice(&dev));
    CK(cudaGetDeviceProperties(&prop, dev));
    printf("========================================================\n");
    printf("  Device  : %s\n", prop.name);
    if (do_vortex)
        printf("  Enriched FV + HLLC-LM  --  Isentropic Vortex  Ma=%.3f\n", Ma);
    else if (do_paper_vortex)
        printf("  Enriched FV + HLLC-LM  --  Paper Isentropic Vortex  [0,10]^2  %s\n",
               pv_moving ? "MOVING (u0=v0=1)" : "STATIONARY");
    else if (do_kh)
        printf("  Enriched FV + HLLC-LM  --  Kelvin-Helmholtz  Ma=0.01  [0,2]^2\n");
    else
        printf("  Enriched FV + HLLC-LM  --  Circular Sod  (no AV)\n");
    printf("  N=%dx%d  Np=%d  h=%.5f  L=%.1f  CFL=%.2f  t_end=%.3f\n",
           N, N, Np, h, L_domain, CFL, t_end);
    printf("  DOFs: rho(P0) + mxa(P0) + mya(P0) + E(P0) + alpha(div-free)  per cell\n");
    printf("========================================================\n");

    size_t sz1  = (size_t)N  * N  * sizeof(float);
    size_t sz5  = (size_t)NVAR * sz1;
    size_t sz6p = (size_t)NVAR * Np * Np * sizeof(float);

    float *d_U, *d_U0, *d_U1, *d_U2, *d_RHS, *d_lam, *d_tmp;
    float *d_rho_out, *d_p_out, *d_err_rho, *d_err_p, *d_err_ru, *d_err_rv, *d_err_re;
    CK(cudaMalloc(&d_U,       sz6p));
    CK(cudaMalloc(&d_U0,      sz6p));
    CK(cudaMalloc(&d_U1,      sz6p));
    CK(cudaMalloc(&d_U2,      sz6p));
    CK(cudaMalloc(&d_RHS,     sz5));
    CK(cudaMalloc(&d_lam,     sz1));
    CK(cudaMalloc(&d_rho_out, sz1));
    CK(cudaMalloc(&d_p_out,   sz1));
    CK(cudaMalloc(&d_err_rho, sz1));
    CK(cudaMalloc(&d_err_p,   sz1));
    CK(cudaMalloc(&d_err_ru,  sz1));
    CK(cudaMalloc(&d_err_rv,  sz1));
    CK(cudaMalloc(&d_err_re,  sz1));
    CK(cudaMalloc(&d_tmp,     GS_NBLK * sizeof(float)));

    float *h_rho = (float*)malloc(sz1);
    float *h_p   = (float*)malloc(sz1);

    /* Color ranges — vortex cases use data-driven symmetric bounds (set after IC) */
    float RHO_MIN, RHO_MAX, P_MIN, P_MAX;
    if (do_vortex || do_paper_vortex) {
        RHO_MIN = 0.f; RHO_MAX = 0.f;   /* sentinel: filled in after IC */
        P_MIN   = 0.f; P_MAX   = 0.f;
    } else if (do_kh) {
        /* color bounds scale with kh_mach; right panel = v-velocity */
        RHO_MIN = GAMMA_V - 2.f*kh_mach*0.1f; RHO_MAX = GAMMA_V + 2.f*kh_mach*0.1f;
        P_MIN   = -1.5f*kh_mach;              P_MAX   = 1.5f*kh_mach;
    } else {
        RHO_MIN = 0.125f; RHO_MAX = 10.0f;
        P_MIN   =  0.10f; P_MAX   = 10.0f;
    }
    const char *prefix = do_vortex       ? "rt_vortex"  :
                         do_paper_vortex ? "rt_pvortex" :
                         do_kh           ? ({static char _buf[32]; snprintf(_buf,sizeof(_buf),"rt_kh%03d",(int)roundf(kh_mach*1000)); _buf;}) : "rt_sod";

    /* Helper: compute and print L2 errors (rho, rho*u, rho*v, rho*E — matches paper Tables 2 & 3) */
#define PRINT_L2_ERR() do { \
    if (do_vortex) { \
        compute_l2_err<<<GS_NBLK,BLOCK1D>>>(d_U, d_err_rho, d_err_p, N, Np, h, Ma); \
        CK(cudaDeviceSynchronize()); \
        float _sr = gpu_sum(d_err_rho, d_tmp, N*N); \
        float _sp = gpu_sum(d_err_p,   d_tmp, N*N); \
        printf("    L2_rho=%.4e  L2_p=%.4e\n", sqrtf(_sr), sqrtf(_sp)); \
    } else if (do_paper_vortex) { \
        /* exact at tf=1 (stationary) and tf=10 (moving) both equal the IC */ \
        compute_l2_err_paper<<<GS_NBLK,BLOCK1D>>>(d_U, d_err_rho, d_err_ru, \
                                                    d_err_rv, d_err_re, \
                                                    N, Np, h, u0, v0); \
        CK(cudaDeviceSynchronize()); \
        float _sr  = gpu_sum(d_err_rho, d_tmp, N*N); \
        float _sru = gpu_sum(d_err_ru,  d_tmp, N*N); \
        float _srv = gpu_sum(d_err_rv,  d_tmp, N*N); \
        float _sre = gpu_sum(d_err_re,  d_tmp, N*N); \
        printf("    L2_rho=%.4e  L2_ru=%.4e  L2_rv=%.4e  L2_rE=%.4e\n", \
               sqrtf(_sr), sqrtf(_sru), sqrtf(_srv), sqrtf(_sre)); \
    } \
} while(0)

#define WRITE_FRAME(idx) do { \
    if (do_kh) \
        extract_rho_v<<<GS_NBLK,BLOCK1D>>>(d_U, d_rho_out, d_p_out, N, Np); \
    else \
        extract_rho_p<<<GS_NBLK,BLOCK1D>>>(d_U, d_rho_out, d_p_out, N, Np); \
    CK(cudaDeviceSynchronize()); \
    CK(cudaMemcpy(h_rho, d_rho_out, sz1, cudaMemcpyDeviceToHost)); \
    CK(cudaMemcpy(h_p,   d_p_out,   sz1, cudaMemcpyDeviceToHost)); \
    char _fn[64]; sprintf(_fn, "%s_%04d.ppm", prefix, (idx)); \
    write_ppm_2panel(_fn, h_rho, RHO_MIN, RHO_MAX, \
                          h_p,   P_MIN,   P_MAX,   N); \
    /* binary dump: 4-byte N, then N*N float32 rho, then N*N float32 p */ \
    { char _bf[64]; sprintf(_bf, "%s_%04d.bin", prefix, (idx)); \
      FILE *_fp = fopen(_bf, "wb"); \
      if (_fp) { fwrite(&N, sizeof(int), 1, _fp); \
                 fwrite(h_rho, sizeof(float), (size_t)N*N, _fp); \
                 fwrite(h_p,   sizeof(float), (size_t)N*N, _fp); \
                 fclose(_fp); } } \
    PRINT_L2_ERR(); \
} while(0)

    /* --- IC and ghost cells --- */
    if (do_vortex)
        ic_vortex      <<<GS_NBLK, BLOCK1D>>>(d_U, N, Np, h, Ma);
    else if (do_paper_vortex)
        ic_vortex_paper<<<GS_NBLK, BLOCK1D>>>(d_U, N, Np, h, u0, v0);
    else if (do_kh)
        ic_kh          <<<GS_NBLK, BLOCK1D>>>(d_U, N, Np, h, kh_mach);
    else
        ic_kernel      <<<GS_NBLK, BLOCK1D>>>(d_U, N, Np, h);
    apply_bc  <<<GS_NBLK, BLOCK1D>>>(d_U, N, Np);
    CK(cudaDeviceSynchronize());

    /* Compute data-driven symmetric color bounds for vortex cases from IC */
    if (do_vortex || do_paper_vortex) {
        extract_rho_p<<<GS_NBLK, BLOCK1D>>>(d_U, d_rho_out, d_p_out, N, Np);
        CK(cudaDeviceSynchronize());
        float rho_lo = gpu_min(d_rho_out, d_tmp, N*N);
        float rho_hi = gpu_max(d_rho_out, d_tmp, N*N);
        float p_lo   = gpu_min(d_p_out,   d_tmp, N*N);
        float p_hi   = gpu_max(d_p_out,   d_tmp, N*N);
        /* Center on the known background value and use the maximum one-sided
         * departure as the half-range (10% padding).
         * Paper vortex background: rho_bg=1, p_bg=1.
         * Original vortex background: rho_bg=1, p_bg=1/gamma. */
        float rho_bg = 1.f;
        float p_bg   = do_paper_vortex ? 1.f : (1.f / GAMMA_V);
        float rho_half = fmaxf(rho_bg - rho_lo, rho_hi - rho_bg) * 1.1f;
        float p_half   = fmaxf(p_bg   - p_lo,   p_hi   - p_bg)   * 1.1f;
        /* Guard against a flat field (e.g. after total diffusion) */
        if (rho_half < 1e-6f) rho_half = 0.1f;
        if (p_half   < 1e-6f) p_half   = 0.1f;
        RHO_MIN = rho_bg - rho_half;
        RHO_MAX = rho_bg + rho_half;
        P_MIN   = p_bg   - p_half;
        P_MAX   = p_bg   + p_half;
        printf("  Color range  rho: [%.4f, %.4f]  p: [%.4f, %.4f]\n",
               RHO_MIN, RHO_MAX, P_MIN, P_MAX);
    }

    WRITE_FRAME(0);

    /* precompute RHS for first step */
    compute_rhs<<<GS_NBLK, BLOCK1D>>>(d_U, d_RHS, d_lam, N, Np, h);
    CK(cudaDeviceSynchronize());
    float lam_max = gpu_max(d_lam, d_tmp, N*N);

    const int  N_FRAMES = 10;
    int   frame   = 1;
    float t_next  = t_end / N_FRAMES;
    int   step    = 0;
    float t       = 0.f;
    struct timespec ts0, ts1;
    clock_gettime(CLOCK_MONOTONIC, &ts0);

    while (t < t_end) {
        if (!(lam_max > 0.f)) {
            fprintf(stderr, "NaN/Inf lam_max at step %d t=%.5f\n", step, t); break;
        }
        float dt = CFL * h / lam_max;
        float t_target = (t_next < t_end) ? t_next : t_end;
        if (t + dt > t_target) dt = t_target - t;
        if (dt < 1e-14f) { fprintf(stderr, "dt underflow\n"); WRITE_FRAME(99); break; }

        CK(cudaMemcpy(d_U0, d_U, sz6p, cudaMemcpyDeviceToDevice));

        /* stage 1 */
        rk3_s1<<<GS_NBLK, BLOCK1D>>>(d_U1, d_U0, d_RHS, dt, N, Np);
        apply_bc<<<GS_NBLK, BLOCK1D>>>(d_U1, N, Np);

        /* stage 2 */
        compute_rhs<<<GS_NBLK, BLOCK1D>>>(d_U1, d_RHS, d_lam, N, Np, h);
        rk3_s2<<<GS_NBLK, BLOCK1D>>>(d_U2, d_U0, d_U1, d_RHS, dt, N, Np);
        apply_bc<<<GS_NBLK, BLOCK1D>>>(d_U2, N, Np);

        /* stage 3 */
        compute_rhs<<<GS_NBLK, BLOCK1D>>>(d_U2, d_RHS, d_lam, N, Np, h);
        rk3_s3<<<GS_NBLK, BLOCK1D>>>(d_U, d_U0, d_U2, d_RHS, dt, N, Np);
        apply_bc<<<GS_NBLK, BLOCK1D>>>(d_U, N, Np);

        /* update RHS and lam for next step */
        compute_rhs<<<GS_NBLK, BLOCK1D>>>(d_U, d_RHS, d_lam, N, Np, h);
        lam_max = gpu_max(d_lam, d_tmp, N*N);
        if (!(lam_max > 0.f)) { fprintf(stderr, "blow-up step %d\n", step); break; }

        t += dt; step++;

        if (t >= t_next - 1e-12f && frame <= N_FRAMES) {
            clock_gettime(CLOCK_MONOTONIC, &ts1);
            float el = (ts1.tv_sec-ts0.tv_sec) + (ts1.tv_nsec-ts0.tv_nsec)*1e-9f;
            printf("  frame %2d/%d  step %5d  t=%.5f  lam=%.3e  elapsed=%.1fs\n",
                   frame, N_FRAMES, step, t, lam_max, el);
            fflush(stdout);
            WRITE_FRAME(frame);
            frame++;
            t_next = frame * t_end / N_FRAMES;
        } else if (step % 100 == 0) {
            clock_gettime(CLOCK_MONOTONIC, &ts1);
            float el = (ts1.tv_sec-ts0.tv_sec) + (ts1.tv_nsec-ts0.tv_nsec)*1e-9f;
            printf("  step %5d  t=%.5f  dt=%.3e  lam=%.3e  elapsed=%.1fs\n",
                   step, t, dt, lam_max, el);
            fflush(stdout);
        }
    }
    /* Final L2 error summary */
    if (do_vortex) {
        compute_l2_err<<<GS_NBLK,BLOCK1D>>>(d_U, d_err_rho, d_err_p, N, Np, h, Ma);
        CK(cudaDeviceSynchronize());
        float sr = gpu_sum(d_err_rho, d_tmp, N*N);
        float sp = gpu_sum(d_err_p,   d_tmp, N*N);
        printf("  ── FINAL L2 ERRORS ──────────────────────────────────\n");
        printf("  Ma=%.3f  N=%d  h=%.5f  t=%.4f\n", Ma, N, h, t);
        printf("  L2(rho) = %.6e\n", sqrtf(sr));
        printf("  L2(p)   = %.6e\n", sqrtf(sp));
        printf("  ────────────────────────────────────────────────────\n");
    } else if (do_paper_vortex) {
        /* exact at tf=1 (stationary) and tf=10 (moving) both equal the IC:
         * the vortex traverses exactly one period on the [0,10]^2 periodic domain. */
        compute_l2_err_paper<<<GS_NBLK,BLOCK1D>>>(d_U, d_err_rho, d_err_ru,
                                                    d_err_rv, d_err_re,
                                                    N, Np, h, u0, v0);
        CK(cudaDeviceSynchronize());
        float sr  = gpu_sum(d_err_rho, d_tmp, N*N);
        float sru = gpu_sum(d_err_ru,  d_tmp, N*N);
        float srv = gpu_sum(d_err_rv,  d_tmp, N*N);
        float sre = gpu_sum(d_err_re,  d_tmp, N*N);
        printf("  ── FINAL L2 ERRORS  [Paper Vortex, Sec 6.2.1]  ─────\n");
        printf("  N=%d  h=%.5f  L=10  t=%.4f  %s\n", N, h, t,
               pv_moving ? "(moving u0=v0=1)" : "(stationary)");
        printf("  L2(rho) = %.6e%s\n", sqrtf(sr),
               pv_moving ? "" : "   (superconvergence → ~2nd order expected)");
        printf("  L2(rho*u)= %.6e\n", sqrtf(sru));
        printf("  L2(rho*v)= %.6e\n", sqrtf(srv));
        printf("  L2(rho*E)= %.6e\n", sqrtf(sre));
        printf("  ────────────────────────────────────────────────────\n");
    }

    CK(cudaDeviceSynchronize());
    clock_gettime(CLOCK_MONOTONIC, &ts1);
    float wall = (ts1.tv_sec-ts0.tv_sec) + (ts1.tv_nsec-ts0.tv_nsec)*1e-9f;
    printf("  Done: %d steps  t=%.6f  wall=%.2fs\n", step, t, wall);
    printf("  Output: %s_0000.ppm .. %s_%04d.ppm\n", prefix, prefix, N_FRAMES);
    printf("  Left panel: density (plasma), Right panel: %s (viridis)\n",
           do_kh ? "v-velocity" : "pressure");

    free(h_rho); free(h_p);
    cudaFree(d_U); cudaFree(d_U0); cudaFree(d_U1); cudaFree(d_U2);
    cudaFree(d_RHS); cudaFree(d_lam); cudaFree(d_tmp);
    cudaFree(d_rho_out); cudaFree(d_p_out);
    cudaFree(d_err_rho); cudaFree(d_err_p); cudaFree(d_err_ru);
    cudaFree(d_err_rv); cudaFree(d_err_re);
    return 0;
}
