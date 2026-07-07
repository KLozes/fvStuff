/*
 * carb_fv_2d.cu
 * Quirk odd-even decoupling / carbuncle test on a PLAIN collocated FV scheme.
 *
 *   Scheme  : first-order Godunov finite volume (no reconstruction, no limiter)
 *   Flux    : HLLC by default — the carbuncle-prone contact-preserving solver.
 *             Override with -Driemann_n=hll_n (carbuncle-free control) or roe_n.
 *   Time    : SSP-RK2, CFL 0.4
 *
 * Test (Quirk 1994, periodic-duct variant):
 *   Planar Mach-M shock (default M=6) moving +x through [0,1]^2, N x N cells.
 *   Pre-shock : rho1=1.4, p1=1, u=v=0  (c1=1, so shock speed = M)
 *   Post-shock: Rankine-Hugoniot state
 *   Trigger   : shock interface staggered by one cell on odd j-rows, plus a
 *               1e-3 odd-even density perturbation of the pre-shock gas.
 *   BCs       : x-left Dirichlet (post-shock), x-right zero-gradient,
 *               y periodic.
 *
 * Diagnostic: max |v| (transverse velocity). A stable planar shock keeps
 * max|v| small and decaying; odd-even decoupling / carbuncle grows it to
 * O(u2) and the front develops the classic sawtooth.
 *
 * Compile:  nvcc -O3 -arch=native --expt-relaxed-constexpr -o carb_fv carb_fv_2d.cu -lm
 * Usage  :  ./carb_fv [N] [Mach]        default N=400, Mach=6
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <cuda_runtime.h>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define GAMMA_V 1.4
#define BLOCK1D 256
#define GS_NBLK 256

#ifdef USE_DOUBLE
#define real double
#else
#define real float
#endif

#define CK(x) do { \
    cudaError_t _e = (x); \
    if (_e != cudaSuccess) { \
        fprintf(stderr, "CUDA error %s:%d: %s\n", __FILE__, __LINE__, \
                cudaGetErrorString(_e)); \
        exit(1); \
    } \
} while (0)

/* conservative state [rho, mx, my, E], SoA, 1 ghost layer */
#define NVAR 4
#define G_GHOST 1
#define IDX_P(q,jp,ip,Np) ((size_t)(q)*(Np)*(Np) + (size_t)(jp)*(Np) + (ip))

__device__ __forceinline__ real
sound_speed(real rho, real p) { return sqrt((real)GAMMA_V * p / rho); }

/* ── HLLC flux with normal n=(nx,ny) — verbatim from rt_dg_euler_2d.cu ──── */
__device__ void
hllc_n(const real WL[4], const real WR[4], real nx, real ny, real F[4])
{
    real rL=WL[0], uL=WL[1], vL=WL[2], pL=WL[3];
    real rR=WR[0], uR=WR[1], vR=WR[2], pR=WR[3];

    real unL = uL*nx + vL*ny,  utL = -uL*ny + vL*nx;
    real unR = uR*nx + vR*ny,  utR = -uR*ny + vR*nx;

    real cL = sound_speed(rL, pL), cR = sound_speed(rR, pR);
    real EL = pL/(GAMMA_V-1.) + 0.5*rL*(uL*uL + vL*vL);
    real ER = pR/(GAMMA_V-1.) + 0.5*rR*(uR*uR + vR*vR);

    real SL = fmin(unL - cL, unR - cR);
    real SR = fmax(unL + cL, unR + cR);

    real FnL[4] = { rL*unL, rL*unL*unL + pL, rL*unL*utL, (EL + pL)*unL };
    real FnR[4] = { rR*unR, rR*unR*unR + pR, rR*unR*utR, (ER + pR)*unR };

    if (SL >= 0.) {
        F[0]=FnL[0]; F[1]=FnL[1]*nx-FnL[2]*ny; F[2]=FnL[1]*ny+FnL[2]*nx; F[3]=FnL[3];
        return;
    }
    if (SR <= 0.) {
        F[0]=FnR[0]; F[1]=FnR[1]*nx-FnR[2]*ny; F[2]=FnR[1]*ny+FnR[2]*nx; F[3]=FnR[3];
        return;
    }

    real dL = rL*(SL - unL), dR = rR*(SR - unR);
    real Ss = (pR - pL + dL*unL - dR*unR) / (dL - dR);

    real factL  = rL*(SL - unL)/(SL - Ss);
    real UsL[4] = { factL, factL*Ss, factL*utL,
                    factL*(EL/rL + (Ss-unL)*(Ss + pL/(rL*(SL-unL)))) };
    real factR  = rR*(SR - unR)/(SR - Ss);
    real UsR[4] = { factR, factR*Ss, factR*utR,
                    factR*(ER/rR + (Ss-unR)*(Ss + pR/(rR*(SR-unR)))) };

    real UL[4] = { rL, rL*unL, rL*utL, EL };
    real UR[4] = { rR, rR*unR, rR*utR, ER };

    real Fn[4];
    if (Ss >= 0.)
        for (int q = 0; q < 4; q++) Fn[q] = FnL[q] + SL*(UsL[q] - UL[q]);
    else
        for (int q = 0; q < 4; q++) Fn[q] = FnR[q] + SR*(UsR[q] - UR[q]);

    F[0]=Fn[0]; F[1]=Fn[1]*nx-Fn[2]*ny; F[2]=Fn[1]*ny+Fn[2]*nx; F[3]=Fn[3];
}

/* ── HLL flux (2-wave, carbuncle-free control) ──────────────────────────── */
__device__ void
hll_n(const real WL[4], const real WR[4], real nx, real ny, real F[4])
{
    real rL=WL[0], uL=WL[1], vL=WL[2], pL=WL[3];
    real rR=WR[0], uR=WR[1], vR=WR[2], pR=WR[3];

    real unL = uL*nx + vL*ny,  utL = -uL*ny + vL*nx;
    real unR = uR*nx + vR*ny,  utR = -uR*ny + vR*nx;

    real cL = sound_speed(rL, pL), cR = sound_speed(rR, pR);
    real EL = pL/(GAMMA_V-1.) + 0.5*rL*(uL*uL + vL*vL);
    real ER = pR/(GAMMA_V-1.) + 0.5*rR*(uR*uR + vR*vR);

    real SL = fmin(unL - cL, unR - cR);
    real SR = fmax(unL + cL, unR + cR);

    real FnL[4] = { rL*unL, rL*unL*unL + pL, rL*unL*utL, (EL + pL)*unL };
    real FnR[4] = { rR*unR, rR*unR*unR + pR, rR*unR*utR, (ER + pR)*unR };

    if (SL >= 0.) {
        F[0]=FnL[0]; F[1]=FnL[1]*nx-FnL[2]*ny; F[2]=FnL[1]*ny+FnL[2]*nx; F[3]=FnL[3];
        return;
    }
    if (SR <= 0.) {
        F[0]=FnR[0]; F[1]=FnR[1]*nx-FnR[2]*ny; F[2]=FnR[1]*ny+FnR[2]*nx; F[3]=FnR[3];
        return;
    }

    real UL[4] = { rL, rL*unL, rL*utL, EL };
    real UR[4] = { rR, rR*unR, rR*utR, ER };

    real denom = 1. / (SR - SL);
    real Fn[4];
    for (int q = 0; q < 4; q++)
        Fn[q] = (SR*FnL[q] - SL*FnR[q] + SL*SR*(UR[q] - UL[q])) * denom;

    F[0]=Fn[0]; F[1]=Fn[1]*nx-Fn[2]*ny; F[2]=Fn[1]*ny+Fn[2]*nx; F[3]=Fn[3];
}

/* ── Roe flux with Harten-Hyman entropy fix (also carbuncle-prone) ──────── */
__device__ void
roe_n(const real WL[4], const real WR[4], real nx, real ny, real F[4])
{
    real rL=WL[0], uL=WL[1], vL=WL[2], pL=WL[3];
    real rR=WR[0], uR=WR[1], vR=WR[2], pR=WR[3];

    real cL = sound_speed(rL, pL), cR = sound_speed(rR, pR);
    real sqrL = sqrt(rL), sqrR = sqrt(rR);
    real denom = sqrL + sqrR;
    real uRoe = (sqrL*uL + sqrR*uR) / denom;
    real vRoe = (sqrL*vL + sqrR*vR) / denom;
    real HL = (pL/(GAMMA_V-1.) + 0.5*rL*(uL*uL+vL*vL) + pL) / rL;
    real HR = (pR/(GAMMA_V-1.) + 0.5*rR*(uR*uR+vR*vR) + pR) / rR;
    real HRoe = (sqrL*HL + sqrR*HR) / denom;
    real c2Roe = (GAMMA_V-1.) * (HRoe - 0.5*(uRoe*uRoe + vRoe*vRoe));
    real cRoe = sqrt(c2Roe > 1e-14 ? c2Roe : 1e-14);

    real unL = uL*nx + vL*ny,  utL = -uL*ny + vL*nx;
    real unR = uR*nx + vR*ny,  utR = -uR*ny + vR*nx;
    real unRoe = uRoe*nx + vRoe*ny;
    real utRoe = -uRoe*ny + vRoe*nx;

    real EL = pL/(GAMMA_V-1.) + 0.5*rL*(uL*uL+vL*vL);
    real ER = pR/(GAMMA_V-1.) + 0.5*rR*(uR*uR+vR*vR);
    real FnL[4] = { rL*unL, rL*unL*unL + pL, rL*unL*utL, (EL+pL)*unL };
    real FnR[4] = { rR*unR, rR*unR*unR + pR, rR*unR*utR, (ER+pR)*unR };

    real drho = rR - rL, dun = unR - unL, dut = utR - utL, dp = pR - pL;
    real rRoe = sqrL * sqrR;

    real b1 = (dp - rRoe*cRoe*dun) / (2.*c2Roe);
    real b2 = drho - dp/c2Roe;
    real b3 = rRoe * dut;
    real b4 = (dp + rRoe*cRoe*dun) / (2.*c2Roe);

    real lam1 = unRoe - cRoe, lam2 = unRoe, lam3 = unRoe, lam4 = unRoe + cRoe;
    real eps1 = fmax((real)0., 2.*((unR - cR) - (unL - cL)));
    real eps4 = fmax((real)0., 2.*((unR + cR) - (unL + cL)));
    real al1 = fabs(lam1); if (al1 < eps1*0.5) al1 = (lam1*lam1 + eps1*eps1*0.25)/eps1;
    real al2 = fabs(lam2);
    real al3 = fabs(lam3);
    real al4 = fabs(lam4); if (al4 < eps4*0.5) al4 = (lam4*lam4 + eps4*eps4*0.25)/eps4;

    real d_rho = al1*b1 + al2*b2 + al4*b4;
    real d_un  = al1*b1*(unRoe-cRoe) + al2*b2*unRoe + al4*b4*(unRoe+cRoe);
    real d_ut  = al1*b1*utRoe + al2*b2*utRoe + al3*b3 + al4*b4*utRoe;
    real d_E   = al1*b1*(HRoe-unRoe*cRoe) + al2*b2*0.5*(unRoe*unRoe+utRoe*utRoe)
               + al3*b3*utRoe + al4*b4*(HRoe+unRoe*cRoe);

    real Fn[4];
    Fn[0] = 0.5*(FnL[0]+FnR[0]) - 0.5*d_rho;
    Fn[1] = 0.5*(FnL[1]+FnR[1]) - 0.5*d_un;
    Fn[2] = 0.5*(FnL[2]+FnR[2]) - 0.5*d_ut;
    Fn[3] = 0.5*(FnL[3]+FnR[3]) - 0.5*d_E;

    F[0]=Fn[0]; F[1]=Fn[1]*nx-Fn[2]*ny; F[2]=Fn[1]*ny+Fn[2]*nx; F[3]=Fn[3];
}

/* flux selector (override: nvcc -Driemann_n=hll_n | roe_n) */
#ifndef riemann_n
#define riemann_n hllc_n
#endif
#define STR2(x) #x
#define STR(x)  STR2(x)

/* ── boundary conditions ────────────────────────────────────────────────── *
 * x-left : Dirichlet post-shock  |  x-right : zero-gradient  |  y : periodic */
__global__ void
apply_bc(real * __restrict__ Qp, int N, int Np,
         real r2, real m2, real E2)
{
    const int g = G_GHOST;
    const int Np2 = Np * Np;
    for (int k = blockIdx.x * blockDim.x + threadIdx.x; k < Np2;
             k += blockDim.x * gridDim.x) {
        int jp = k / Np, ip = k % Np;
        bool gx = (ip < g || ip >= N + g);
        bool gy = (jp < g || jp >= N + g);
        if (!gx && !gy) continue;

        /* y periodic first (also fills corners consistently) */
        int jr = jp;
        if (jp < g)      jr = jp + N;
        if (jp >= N + g) jr = jp - N;

        if (ip < g) {                       /* left: post-shock Dirichlet */
            Qp[IDX_P(0, jp, ip, Np)] = r2;
            Qp[IDX_P(1, jp, ip, Np)] = m2;
            Qp[IDX_P(2, jp, ip, Np)] = 0.;
            Qp[IDX_P(3, jp, ip, Np)] = E2;
        } else if (ip >= N + g) {           /* right: zero-gradient */
            for (int q = 0; q < NVAR; q++)
                Qp[IDX_P(q, jp, ip, Np)] = Qp[IDX_P(q, jr, N + g - 1, Np)];
        } else {                            /* pure y ghost */
            for (int q = 0; q < NVAR; q++)
                Qp[IDX_P(q, jp, ip, Np)] = Qp[IDX_P(q, jr, ip, Np)];
        }
    }
}

/* ── first-order Godunov RHS (gather form) ──────────────────────────────── */
__device__ __forceinline__ void
prim_at(const real * __restrict__ Qp, int jp, int ip, int Np, real W[4])
{
    real r = fmax(Qp[IDX_P(0, jp, ip, Np)], (real)1e-14);
    real u = Qp[IDX_P(1, jp, ip, Np)] / r;
    real v = Qp[IDX_P(2, jp, ip, Np)] / r;
    real E = Qp[IDX_P(3, jp, ip, Np)];
    W[0] = r; W[1] = u; W[2] = v;
    W[3] = fmax((real)((GAMMA_V-1.) * (E - 0.5*r*(u*u + v*v))), (real)1e-14);
}

__global__ void
compute_rhs(const real * __restrict__ Qp,
            real       * __restrict__ RHS,
            real       * __restrict__ lam_out,
            int N, int Np, real h)
{
    const int g = G_GHOST;
    const int N2 = N * N;
    for (int k = blockIdx.x * blockDim.x + threadIdx.x; k < N2;
             k += blockDim.x * gridDim.x) {
        int j = k / N, i = k % N;
        int jp = j + g, ip = i + g;

        real Wc[4], WL[4], WR[4], WB[4], WT[4];
        prim_at(Qp, jp, ip,   Np, Wc);
        prim_at(Qp, jp, ip-1, Np, WL);
        prim_at(Qp, jp, ip+1, Np, WR);
        prim_at(Qp, jp-1, ip, Np, WB);
        prim_at(Qp, jp+1, ip, Np, WT);

        real FR[4], FL[4], GT[4], GB[4];
        riemann_n(Wc, WR, 1., 0., FR);
        riemann_n(WL, Wc, 1., 0., FL);
        riemann_n(Wc, WT, 0., 1., GT);
        riemann_n(WB, Wc, 0., 1., GB);

        real inv_h = 1. / h;
        for (int q = 0; q < NVAR; q++)
            RHS[q * N2 + k] = -(FR[q] - FL[q] + GT[q] - GB[q]) * inv_h;

        lam_out[k] = fabs(Wc[1]) + fabs(Wc[2]) + sound_speed(Wc[0], Wc[3]);
    }
}

/* ── SSP-RK2 stages ─────────────────────────────────────────────────────── */
__global__ void
rk2_s1(real * __restrict__ U1, const real * __restrict__ U0,
       const real * __restrict__ L, real dt, int N, int Np)
{
    const int g = G_GHOST;
    const int N2 = N * N, Np2 = Np * Np;
    for (int k = blockIdx.x * blockDim.x + threadIdx.x; k < N2;
             k += blockDim.x * gridDim.x) {
        int j = k / N, i = k % N;
        int idxp = (j + g) * Np + (i + g);
        for (int q = 0; q < NVAR; q++)
            U1[q * Np2 + idxp] = U0[q * Np2 + idxp] + dt * L[q * N2 + k];
        U1[0 * Np2 + idxp] = fmax(U1[0 * Np2 + idxp], (real)1e-14);
        U1[3 * Np2 + idxp] = fmax(U1[3 * Np2 + idxp], (real)1e-14);
    }
}

__global__ void
rk2_s2(real * __restrict__ U, const real * __restrict__ U0,
       const real * __restrict__ U1, const real * __restrict__ L,
       real dt, int N, int Np)
{
    const int g = G_GHOST;
    const int N2 = N * N, Np2 = Np * Np;
    for (int k = blockIdx.x * blockDim.x + threadIdx.x; k < N2;
             k += blockDim.x * gridDim.x) {
        int j = k / N, i = k % N;
        int idxp = (j + g) * Np + (i + g);
        for (int q = 0; q < NVAR; q++)
            U[q * Np2 + idxp] = 0.5 * U0[q * Np2 + idxp]
                              + 0.5 * (U1[q * Np2 + idxp] + dt * L[q * N2 + k]);
        U[0 * Np2 + idxp] = fmax(U[0 * Np2 + idxp], (real)1e-14);
        U[3 * Np2 + idxp] = fmax(U[3 * Np2 + idxp], (real)1e-14);
    }
}

/* pseudo-random in [-0.5, 0.5) from cell indices (symmetry breaker) */
__device__ __forceinline__ real
cell_hash(int i, int j)
{
    unsigned int x = (unsigned int)(i * 73856093) ^ (unsigned int)(j * 19349663);
    x ^= x >> 13; x *= 0x85ebca6bu; x ^= x >> 16;
    return (real)(x & 0xFFFFu) / (real)65536. - (real)0.5;
}

/* ── IC: staggered planar shock + odd-even AND random pre-shock density
 *    perturbation.  The odd-even part seeds the decoupled mode; the random
 *    part breaks the exact period-2 symmetry so transverse velocity can
 *    couple (a perfectly alternating pattern has G_T ≡ G_B and generates
 *    NO v by symmetry — the decoupling itself!). ─────────────────────────── */
__global__ void
ic_quirk(real * __restrict__ Qp, int N, int Np, real h,
         real r1, real p1, real r2, real u2, real p2, real x0, real pert)
{
    const int g = G_GHOST;
    const int N2 = N * N;
    for (int k = blockIdx.x * blockDim.x + threadIdx.x; k < N2;
             k += blockDim.x * gridDim.x) {
        int j = k / N, i = k % N;
        int jp = j + g, ip = i + g;

        /* shock interface index, staggered one cell on odd rows */
        int is = (int)(x0 * N + 0.5) + (j & 1);

        real rho, u, p;
        if (i < is) { rho = r2; u = u2; p = p2; }             /* post-shock */
        else {                                                 /* pre-shock  */
            rho = r1 * (1. + pert * ((j & 1) ? 1. : -1.)
                           + pert * cell_hash(i, j));
            u = 0.; p = p1;
        }
        Qp[IDX_P(0, jp, ip, Np)] = rho;
        Qp[IDX_P(1, jp, ip, Np)] = rho * u;
        Qp[IDX_P(2, jp, ip, Np)] = 0.;
        Qp[IDX_P(3, jp, ip, Np)] = p / (GAMMA_V - 1.) + 0.5 * rho * u * u;
    }
}

/* ── extract rho and v for diagnostics/plots ────────────────────────────── */
__global__ void
extract_rho_v(const real * __restrict__ Qp,
              real * __restrict__ out_rho, real * __restrict__ out_v,
              real * __restrict__ out_absv, int N, int Np)
{
    const int g = G_GHOST;
    const int N2 = N * N;
    for (int k = blockIdx.x * blockDim.x + threadIdx.x; k < N2;
             k += blockDim.x * gridDim.x) {
        int j = k / N, i = k % N;
        int jp = j + g, ip = i + g;
        real r = fmax(Qp[IDX_P(0, jp, ip, Np)], (real)1e-14);
        real v = Qp[IDX_P(2, jp, ip, Np)] / r;
        out_rho[k] = r;
        out_v[k] = v;
        out_absv[k] = fabs(v);
    }
}

/* ── reductions ─────────────────────────────────────────────────────────── */
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
    real v;
    CK(cudaMemcpy(&v, d_tmp, sizeof(real), cudaMemcpyDeviceToHost));
    return v;
}

/* ── 2-panel PNG (jet colormap): rho | v ────────────────────────────────── */
static void
write_png_2panel(const char *fname,
                 const real *rho, real rlo, real rhi,
                 const real *vv,  real wv, int N)
{
    const int W = 2 * N, H = N;
    unsigned char *px = (unsigned char *)malloc((size_t)W * H * 3);
    if (!px) return;
    for (int r = 0; r < H; r++) {
        int pj = N - 1 - r;
        for (int c = 0; c < W; c++) {
            real t;
            if (c < N) t = (rho[pj * N + c] - rlo) / (rhi - rlo + (real)1e-30);
            else       t = (vv[pj * N + (c - N)] + wv) / (2. * wv + (real)1e-30);
            if (!(t > 0.)) t = 0.;
            if (t > 1.) t = 1.;
            real cr = fmax((real)0., fmin((real)1., (real)(1.5 - fabs(4.*t - 3.))));
            real cg = fmax((real)0., fmin((real)1., (real)(1.5 - fabs(4.*t - 2.))));
            real cb = fmax((real)0., fmin((real)1., (real)(1.5 - fabs(4.*t - 1.))));
            unsigned char *p = px + (r * W + c) * 3;
            p[0] = (unsigned char)(cr * 255.);
            p[1] = (unsigned char)(cg * 255.);
            p[2] = (unsigned char)(cb * 255.);
        }
    }
    stbi_write_png(fname, W, H, 3, px, W * 3);
    free(px);
    printf("  Saved %s\n", fname);
}

/* ═════════════════════════════════════════════════════════════════════════ */
int main(int argc, char **argv)
{
    int  N    = (argc > 1) ? atoi(argv[1]) : 400;
    real Mach = (argc > 2) ? atof(argv[2]) : 6.;

    const real g1 = GAMMA_V;
    /* pre-shock: c1 = 1 */
    const real r1 = 1.4, p1 = 1.0;
    /* Rankine-Hugoniot post-shock (shock speed s = Mach) */
    const real r2 = r1 * ((g1+1.)*Mach*Mach) / ((g1-1.)*Mach*Mach + 2.);
    const real p2 = p1 * (2.*g1*Mach*Mach - (g1-1.)) / (g1+1.);
    const real u2 = Mach * (1. - r1/r2);
    const real E2 = p2/(g1-1.) + 0.5*r2*u2*u2;
    const real x0 = 0.1;
    const real pert = 1e-3;          /* odd-even pre-shock density perturbation */

    real h = 1. / N;
    const int g = G_GHOST;
    int Np = N + 2 * g;
    real t_end = 0.75 / Mach;        /* shock travels x0=0.1 → ~0.85 */
    real CFL = 0.4;

    printf("========================================================\n");
    printf("  Quirk odd-even / carbuncle test — 1st-order FV + %s\n", STR(riemann_n));
    printf("  N=%dx%d  Mach=%.1f  shock x0=%.2f -> ~0.85  t_end=%.4f\n",
           N, N, (double)Mach, (double)x0, (double)t_end);
    printf("  post-shock: rho2=%.4f u2=%.4f p2=%.4f   pert=%.0e\n",
           (double)r2, (double)u2, (double)p2, (double)pert);
    printf("  diagnostic: max|v| (planar shock => should stay ~0)\n");
    printf("========================================================\n");

    size_t sz1  = (size_t)N * N * sizeof(real);
    size_t sz4  = (size_t)NVAR * sz1;
    size_t sz4p = (size_t)NVAR * Np * Np * sizeof(real);

    real *d_U, *d_U0, *d_U1, *d_RHS, *d_lam, *d_tmp, *d_rho, *d_v, *d_absv;
    CK(cudaMalloc(&d_U,   sz4p));
    CK(cudaMalloc(&d_U0,  sz4p));
    CK(cudaMalloc(&d_U1,  sz4p));
    CK(cudaMalloc(&d_RHS, sz4));
    CK(cudaMalloc(&d_lam, sz1));
    CK(cudaMalloc(&d_rho, sz1));
    CK(cudaMalloc(&d_v,   sz1));
    CK(cudaMalloc(&d_absv,sz1));
    CK(cudaMalloc(&d_tmp, GS_NBLK * sizeof(real)));

    real *h_rho = (real *)malloc(sz1);
    real *h_v   = (real *)malloc(sz1);

    ic_quirk<<<GS_NBLK, BLOCK1D>>>(d_U, N, Np, h, r1, p1, r2, u2, p2, x0, pert);
    apply_bc<<<GS_NBLK, BLOCK1D>>>(d_U, N, Np, r2, r2*u2, E2);
    CK(cudaDeviceSynchronize());

    const int N_FRAMES = 8;
    int frame = 1, step = 0;
    real t = 0., t_next = t_end / N_FRAMES;
    struct timespec ts0, ts1;
    clock_gettime(CLOCK_MONOTONIC, &ts0);

    /* odd-even sawtooth amplitude: max |rho_ij − ½(rho_i,j−1 + rho_i,j+1)|
     * (y-second-difference — isolates the period-2 decoupled mode)          */
    real saw = 0.;
#define DIAG_FRAME(idx) do { \
    extract_rho_v<<<GS_NBLK, BLOCK1D>>>(d_U, d_rho, d_v, d_absv, N, Np); \
    CK(cudaDeviceSynchronize()); \
    real vmax = gpu_max(d_absv, d_tmp, N*N); \
    real rmax = gpu_max(d_rho,  d_tmp, N*N); \
    CK(cudaMemcpy(h_rho, d_rho, sz1, cudaMemcpyDeviceToHost)); \
    CK(cudaMemcpy(h_v,   d_v,   sz1, cudaMemcpyDeviceToHost)); \
    saw = 0.; \
    for (int _j = 0; _j < N; _j++) { \
        int _jm = (_j == 0) ? N-1 : _j-1, _jp2 = (_j == N-1) ? 0 : _j+1; \
        for (int _i = 0; _i < N; _i++) { \
            real _s = fabs(h_rho[_j*N+_i] \
                     - 0.5*(h_rho[_jm*N+_i] + h_rho[_jp2*N+_i])); \
            if (_s > saw) saw = _s; } } \
    printf("  frame %2d/%d  step %6d  t=%.5f  max|v|=%.4e  sawtooth=%.4e  rho_max=%.4f\n", \
           (idx), N_FRAMES, step, t, vmax, saw, rmax); \
    real _rlo = h_rho[0], _rhi = h_rho[0]; \
    for (int _k = 1; _k < N*N; _k++) { \
        if (h_rho[_k] < _rlo) _rlo = h_rho[_k]; \
        if (h_rho[_k] > _rhi) _rhi = h_rho[_k]; } \
    char _fn[96]; \
    snprintf(_fn, sizeof(_fn), "figures/carb_%s_%02d.png", STR(riemann_n), (idx)); \
    write_png_2panel(_fn, h_rho, _rlo, _rhi, h_v, fmax(vmax, (real)1e-10), N); \
    fflush(stdout); \
} while (0)

    DIAG_FRAME(0);

    compute_rhs<<<GS_NBLK, BLOCK1D>>>(d_U, d_RHS, d_lam, N, Np, h);
    CK(cudaDeviceSynchronize());
    real lam_max = gpu_max(d_lam, d_tmp, N * N);

    while (t < t_end) {
        if (!(lam_max > 0.)) { fprintf(stderr, "blow-up step %d\n", step); break; }
        real dt = CFL * h / lam_max;
        real t_target = (t_next < t_end) ? t_next : t_end;
        if (t + dt > t_target) dt = t_target - t;
        if (dt < 1e-14) { fprintf(stderr, "dt underflow\n"); break; }

        CK(cudaMemcpy(d_U0, d_U, sz4p, cudaMemcpyDeviceToDevice));

        rk2_s1<<<GS_NBLK, BLOCK1D>>>(d_U1, d_U0, d_RHS, dt, N, Np);
        apply_bc<<<GS_NBLK, BLOCK1D>>>(d_U1, N, Np, r2, r2*u2, E2);
        compute_rhs<<<GS_NBLK, BLOCK1D>>>(d_U1, d_RHS, d_lam, N, Np, h);
        rk2_s2<<<GS_NBLK, BLOCK1D>>>(d_U, d_U0, d_U1, d_RHS, dt, N, Np);
        apply_bc<<<GS_NBLK, BLOCK1D>>>(d_U, N, Np, r2, r2*u2, E2);

        compute_rhs<<<GS_NBLK, BLOCK1D>>>(d_U, d_RHS, d_lam, N, Np, h);
        lam_max = gpu_max(d_lam, d_tmp, N * N);

        t += dt; step++;

        if (t >= t_next - 1e-12 && frame <= N_FRAMES) {
            DIAG_FRAME(frame);
            frame++;
            t_next = frame * t_end / N_FRAMES;
        }
    }

    extract_rho_v<<<GS_NBLK, BLOCK1D>>>(d_U, d_rho, d_v, d_absv, N, Np);
    CK(cudaDeviceSynchronize());
    real vmax = gpu_max(d_absv, d_tmp, N * N);
    clock_gettime(CLOCK_MONOTONIC, &ts1);
    double wall = (ts1.tv_sec - ts0.tv_sec) + (ts1.tv_nsec - ts0.tv_nsec) * 1e-9;

    printf("  ── VERDICT ─────────────────────────────────────────\n");
    printf("  flux=%s  Mach=%.1f  N=%d  steps=%d  wall=%.1fs\n",
           STR(riemann_n), (double)Mach, N, step, wall);
    printf("  final max|v| = %.4e  (max|v|/u2 = %.4e)\n", vmax, vmax / u2);
    printf("  final sawtooth = %.4e  (sawtooth/rho2 = %.4e)\n", saw, saw / r2);
    printf("  %s\n", (vmax / u2 > 0.05 || saw / r2 > 0.05)
           ? ">> ODD-EVEN DECOUPLING / CARBUNCLE"
           : ">> stable planar shock (no carbuncle)");
    printf("  ────────────────────────────────────────────────────\n");

    free(h_rho); free(h_v);
    cudaFree(d_U); cudaFree(d_U0); cudaFree(d_U1); cudaFree(d_RHS);
    cudaFree(d_lam); cudaFree(d_tmp); cudaFree(d_rho); cudaFree(d_v);
    cudaFree(d_absv);
    return 0;
}
