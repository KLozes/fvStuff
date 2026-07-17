/*
 * fv_rec3_2d.cu
 * Finite-volume Euler solver, UNLIMITED quadratic (rec3 / κ=1/3 parabolic)
 * face reconstruction + HLLC, stabilized ONLY by an a-posteriori Zhang-Shu
 * positivity limiter applied on non-overlapping 2×2 CELL BLOCKS:
 *
 *   • Reconstruction: 3rd-order upwind parabola in CONSERVED variables,
 *       UL_{i+1/2} = (−U_{i−1} + 5U_i + 2U_{i+1})/6
 *       UR_{i+1/2} = ( 2U_i + 5U_{i+1} − U_{i+2})/6
 *     dimension-by-dimension, face-midpoint quadrature.  NO slope/flux
 *     limiter of any kind.
 *   • Positivity (the q1_dg technique transplanted to FV): each 2×2 block
 *     plays the role of a Q1 element — its 4 cell averages are the "nodal
 *     DOFs".  After every RK stage, each block is scaled toward its BLOCK
 *     MEAN by a common θ ∈ [0,1] until every cell satisfies ρ ≥ ρ_floor and
 *     p ≥ p_floor.  Conservative within the block (block total unchanged),
 *     exactly like the elementwise Zhang-Shu limiter of the DG code.
 *   • Time: SSP-RK3.   Grid: periodic [0,L]², N even, g=2 ghosts.
 *
 * Compile:
 *   nvcc -O3 -arch=native --expt-relaxed-constexpr -o fv_rec3 fv_rec3_2d.cu -lm
 *   (precision: add -Dreal=double -Dreall=double ; CFL: -DCFL_DEFAULT=0.4)
 *
 * Run:
 *   ./fv_rec3 N            circular Sod blast, [0,1]², t=1
 *   ./fv_rec3 N pv         stationary paper vortex, [0,10]², t=1 (L2 errors)
 *   ./fv_rec3 N lmv [eps]  Gresho low-Mach vortex, [0,1]², t=1 (velocity err)
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <cuda_runtime.h>

#define GAMMA_V 1.4
#ifndef SOD_RHOIN
#define SOD_RHOIN 10.0
#endif
#ifndef SOD_PIN
#define SOD_PIN 10.0
#endif
#ifndef CFL_DEFAULT
#define CFL_DEFAULT 0.4
#endif
#define BLOCK1D 256
#define GS_NBLK 256
#define G_GHOST 2
#define NVAR 4          /* ρ, ρu, ρv, E — cell averages only */

#ifndef real
#define real float
#endif
#ifndef reall
#define reall float
#endif

#define PFLOOR   ((reall)1e-9)
#define RHOFLOOR ((real )1e-12)

#define CK(x) do { cudaError_t _e = (x); if (_e != cudaSuccess) { \
    fprintf(stderr,"CUDA error %s:%d: %s\n",__FILE__,__LINE__,cudaGetErrorString(_e)); exit(1);} } while(0)

#define IDX_P(q,jp,ip,Np) ((size_t)(q)*(Np)*(Np) + (size_t)(jp)*(Np) + (ip))

/* ── point helpers (same as q1_dg) ──────────────────────────────────────── */
__device__ __forceinline__ real
sound_speed(real rho, real p) { return (real)sqrt((reall)GAMMA_V*(reall)p/(reall)rho); }

__device__ __forceinline__ void
cons_pt_to_W(real rho, real mx, real my, real E, real W[4])
{
    real  r = fmax(rho, RHOFLOOR);
    reall u = (reall)mx / (reall)r;
    reall v = (reall)my / (reall)r;
    reall p = ((reall)GAMMA_V - 1.0) * ((reall)E - 0.5 * (reall)r * (u*u + v*v));
    W[0] = r; W[1] = (real)u; W[2] = (real)v; W[3] = fmax((real)p, (real)PFLOOR);
}

__device__ __forceinline__ reall
pressure_of(const real U[4])
{
    reall r = fmax(U[0], RHOFLOOR);
    reall u = (reall)U[1]/r, v = (reall)U[2]/r;
    return ((reall)GAMMA_V - 1.0) * ((reall)U[3] - 0.5*r*(u*u + v*v));
}

/* ── HLLC flux, normal (nx,ny) — identical to q1_dg_euler_2d.cu ─────────── */
__device__ void
hllc_n(const real WL[4], const real WR[4], real nx, real ny, real F[4])
{
    real rL=WL[0], uL=WL[1], vL=WL[2], pL=WL[3];
    real rR=WR[0], uR=WR[1], vR=WR[2], pR=WR[3];
    real unL = uL*nx + vL*ny, unR = uR*nx + vR*ny;
    real utL =-uL*ny + vL*nx, utR =-uR*ny + vR*nx;
    real cL = sound_speed(rL,pL), cR = sound_speed(rR,pR);
    real EL = pL/(GAMMA_V-1.) + 0.5*rL*(uL*uL+vL*vL);
    real ER = pR/(GAMMA_V-1.) + 0.5*rR*(uR*uR+vR*vR);
    real SL = fmin(unL-cL, unR-cR), SR = fmax(unL+cL, unR+cR);
    real FnL[4] = { rL*unL, rL*unL*unL+pL, rL*unL*utL, (EL+pL)*unL };
    real FnR[4] = { rR*unR, rR*unR*unR+pR, rR*unR*utR, (ER+pR)*unR };
    if (SL >= 0.) { F[0]=FnL[0]; F[1]=FnL[1]*nx-FnL[2]*ny; F[2]=FnL[1]*ny+FnL[2]*nx; F[3]=FnL[3]; return; }
    if (SR <= 0.) { F[0]=FnR[0]; F[1]=FnR[1]*nx-FnR[2]*ny; F[2]=FnR[1]*ny+FnR[2]*nx; F[3]=FnR[3]; return; }
    real dL = rL*(SL-unL), dR = rR*(SR-unR);
    real Ss = (pR - pL + dL*unL - dR*unR) / (dL - dR);
    real factL = rL*(SL-unL)/(SL-Ss);
    real UsL[4] = { factL, factL*Ss, factL*utL,
                    factL*(EL/rL + (Ss-unL)*(Ss + pL/(rL*(SL-unL)))) };
    real factR = rR*(SR-unR)/(SR-Ss);
    real UsR[4] = { factR, factR*Ss, factR*utR,
                    factR*(ER/rR + (Ss-unR)*(Ss + pR/(rR*(SR-unR)))) };
    real UL[4] = { rL, rL*unL, rL*utL, EL };
    real UR[4] = { rR, rR*unR, rR*utR, ER };
    real Fn[4];
    if (Ss >= 0.) for (int q=0;q<4;q++) Fn[q] = FnL[q] + SL*(UsL[q]-UL[q]);
    else          for (int q=0;q<4;q++) Fn[q] = FnR[q] + SR*(UsR[q]-UR[q]);
    F[0]=Fn[0]; F[1]=Fn[1]*nx-Fn[2]*ny; F[2]=Fn[1]*ny+Fn[2]*nx; F[3]=Fn[3];
}

/* ── Zhang-Shu scaling of ONE reconstructed state toward the owning cell's
 * (admissible) average until ρ ≥ floor and p ≥ floor.  This is the pointwise
 * half of the DG technique: positivity of the polynomial EVALUATIONS fed to
 * the Riemann solver.  Positivity-only — no TVD/flux limiting. ────────────── */
__device__ __forceinline__ void
scale_to_admissible(const real Ub[4], real U[4])
{
    /* ρ pass (linear) */
    if (U[0] < RHOFLOOR) {
        real d = Ub[0] - U[0];
        real th = (d > (real)0.) ? (Ub[0] - RHOFLOOR)/d : (real)0.;
        for (int q = 0; q < 4; q++) U[q] = Ub[q] + th*(U[q] - Ub[q]);
    }
    /* p pass (bisection; p concave along the segment toward the mean) */
    if (pressure_of(U) < PFLOOR) {
        real lo = 0., hi = 1.;
        for (int it = 0; it < 25; it++) {
            real mid = (real)0.5*(lo + hi), V[4];
            for (int q = 0; q < 4; q++) V[q] = Ub[q] + mid*(U[q] - Ub[q]);
            if (V[0] >= RHOFLOOR && pressure_of(V) >= PFLOOR) lo = mid; else hi = mid;
        }
        for (int q = 0; q < 4; q++) U[q] = Ub[q] + lo*(U[q] - Ub[q]);
    }
}

/* ── periodic BC on the padded array ────────────────────────────────────── */
__global__ void
apply_bc(real * __restrict__ Qp, int N, int Np)
{
    const int g = G_GHOST, Np2 = Np*Np;
    for (int k = blockIdx.x*blockDim.x + threadIdx.x; k < Np2; k += blockDim.x*gridDim.x) {
        int jp = k/Np, ip = k%Np;
        bool gx = (ip < g || ip >= N+g), gy = (jp < g || jp >= N+g);
        if (!gx && !gy) continue;
        int ir = ip, jr = jp;
        if (ip < g)     ir = ip + N;
        if (ip >= N+g)  ir = ip - N;
        if (jp < g)     jr = jp + N;
        if (jp >= N+g)  jr = jp - N;
        for (int q = 0; q < NVAR; q++)
            Qp[IDX_P(q,jp,ip,Np)] = Qp[IDX_P(q,jr,ir,Np)];
    }
}

/* ── RHS: unlimited rec3 faces + HLLC, scatter to both cells ────────────── */
__global__ void
compute_rhs(const real * __restrict__ Qp, real * __restrict__ RHS,
            real * __restrict__ lam_out, int N, int Np, real h)
{
    const int g = G_GHOST, N2 = N*N;
    for (int k = blockIdx.x*blockDim.x + threadIdx.x; k < N2; k += blockDim.x*gridDim.x) {
        int j = k/N, i = k%N;
        int jp = j+g, ip = i+g;
        int kR = j*N + ((i+1<N)?(i+1):0);
        int kT = ((j+1<N)?(j+1):0)*N + i;
        real inv_h = 1./h;

        /* right face (i+1/2, j): stencil i-1 .. i+2 */
        {
            real UL[4], UR[4];
            for (int q = 0; q < 4; q++) {
                real um = Qp[IDX_P(q,jp,ip-1,Np)], uc = Qp[IDX_P(q,jp,ip,Np)];
                real up = Qp[IDX_P(q,jp,ip+1,Np)], u2 = Qp[IDX_P(q,jp,ip+2,Np)];
                UL[q] = (-um + 5.*uc + 2.*up) * (real)(1./6.);
                UR[q] = ( 2.*uc + 5.*up - u2) * (real)(1./6.);
            }
            real UbL[4], UbR[4];
            for (int q = 0; q < 4; q++) {
                UbL[q] = Qp[IDX_P(q,jp,ip,  Np)];
                UbR[q] = Qp[IDX_P(q,jp,ip+1,Np)];
            }
#ifndef NO_FACE_LIMITER
            scale_to_admissible(UbL, UL);
            scale_to_admissible(UbR, UR);
#else
            (void)UbL; (void)UbR;
#endif
            real WL[4], WR[4], F[4];
            cons_pt_to_W(UL[0],UL[1],UL[2],UL[3], WL);
            cons_pt_to_W(UR[0],UR[1],UR[2],UR[3], WR);
            hllc_n(WL, WR, 1., 0., F);
            for (int q = 0; q < 4; q++) {
                atomicAdd(&RHS[q*N2 + k ], -F[q]*inv_h);
                atomicAdd(&RHS[q*N2 + kR],  F[q]*inv_h);
            }
        }
        /* top face (i, j+1/2): stencil j-1 .. j+2 */
        {
            real UB[4], UT[4];
            for (int q = 0; q < 4; q++) {
                real um = Qp[IDX_P(q,jp-1,ip,Np)], uc = Qp[IDX_P(q,jp,ip,Np)];
                real up = Qp[IDX_P(q,jp+1,ip,Np)], u2 = Qp[IDX_P(q,jp+2,ip,Np)];
                UB[q] = (-um + 5.*uc + 2.*up) * (real)(1./6.);
                UT[q] = ( 2.*uc + 5.*up - u2) * (real)(1./6.);
            }
            real UbB[4], UbT[4];
            for (int q = 0; q < 4; q++) {
                UbB[q] = Qp[IDX_P(q,jp,  ip,Np)];
                UbT[q] = Qp[IDX_P(q,jp+1,ip,Np)];
            }
#ifndef NO_FACE_LIMITER
            scale_to_admissible(UbB, UB);
            scale_to_admissible(UbT, UT);
#else
            (void)UbB; (void)UbT;
#endif
            real WB[4], WT[4], Gf[4];
            cons_pt_to_W(UB[0],UB[1],UB[2],UB[3], WB);
            cons_pt_to_W(UT[0],UT[1],UT[2],UT[3], WT);
            hllc_n(WB, WT, 0., 1., Gf);
            for (int q = 0; q < 4; q++) {
                atomicAdd(&RHS[q*N2 + k ], -Gf[q]*inv_h);
                atomicAdd(&RHS[q*N2 + kT],  Gf[q]*inv_h);
            }
        }
        /* spectral radius */
        real U0[4] = { Qp[IDX_P(0,jp,ip,Np)], Qp[IDX_P(1,jp,ip,Np)],
                       Qp[IDX_P(2,jp,ip,Np)], Qp[IDX_P(3,jp,ip,Np)] };
        real W[4]; cons_pt_to_W(U0[0],U0[1],U0[2],U0[3], W);
        lam_out[k] = fabs(W[1]) + fabs(W[2]) + sound_speed(W[0], W[3]);
    }
}

/* ── a-posteriori Zhang-Shu positivity on 2×2 CELL BLOCKS ────────────────
 * Non-overlapping tiling (N even).  The 4 cell averages of a block are the
 * "nodal DOFs" of a Q1 element; scale them toward the BLOCK MEAN by a common
 * θ until every cell has ρ ≥ ρ_floor and p ≥ p_floor.  Block total conserved.
 * ρ pass: θ from the linear formula.  p pass: bisection on the common θ
 * (p is concave along the segment toward the mean, so bisection is safe). */
__device__ int g_brutal_count = 0;

__global__ void
limit_positivity_blocks(real * __restrict__ Qp, int N, int Np)
{
    const int g = G_GHOST, Nb = N/2, NB2 = Nb*Nb;
    for (int kb = blockIdx.x*blockDim.x + threadIdx.x; kb < NB2; kb += blockDim.x*gridDim.x) {
        int jb = kb/Nb, ib = kb%Nb;
        int i0 = 2*ib + g, j0 = 2*jb + g;

        real U[4][4], Ub[4] = {0.,0.,0.,0.};
        for (int c = 0; c < 4; c++) {
            int ip = i0 + (c & 1), jp = j0 + (c >> 1);
            for (int q = 0; q < 4; q++) {
                U[c][q] = Qp[IDX_P(q,jp,ip,Np)];
                Ub[q]  += (real)0.25 * U[c][q];
            }
        }

        /* if even the block mean is inadmissible, floor it brutally (rare) */
        reall pb = pressure_of(Ub);
        if (Ub[0] < RHOFLOOR || pb < PFLOOR) {
            atomicAdd(&g_brutal_count, 1);
            real rb = fmax(Ub[0], (real)(10.*RHOFLOOR));
            real Efix = (real)(10.*PFLOOR)/((real)GAMMA_V - 1.)
                      + (real)0.5*(Ub[1]*Ub[1] + Ub[2]*Ub[2])/rb;
            Ub[0] = rb;  Ub[3] = fmax(Ub[3], Efix);
        }

        /* BILINEAR-EXTRAPOLATED BLOCK CORNER states: the 4 cell averages are
         * nodal values of a bilinear at the quadrant centres; extrapolate to
         * the block corners.  Any point of the bilinear is a convex combo of
         * the 4 corner states and p is concave, so corner admissibility ⇒
         * admissibility of the whole bilinear (the DG corner argument).
         * Corner values are LINEAR in the cell states, so θ-scaling the cells
         * toward the block mean scales the corner deviations identically.   */
        real C[4][4];
        for (int q = 0; q < 4; q++) {
            C[0][q] = ((real)9.*U[0][q] - (real)3.*U[1][q] - (real)3.*U[2][q] + U[3][q])*(real)0.25;
            C[1][q] = ((real)9.*U[1][q] - (real)3.*U[0][q] - (real)3.*U[3][q] + U[2][q])*(real)0.25;
            C[2][q] = ((real)9.*U[2][q] - (real)3.*U[3][q] - (real)3.*U[0][q] + U[1][q])*(real)0.25;
            C[3][q] = ((real)9.*U[3][q] - (real)3.*U[2][q] - (real)3.*U[1][q] + U[0][q])*(real)0.25;
        }

        /* ρ pass: common linear θ from the CORNER densities */
        real th = 1.;
        for (int c = 0; c < 4; c++)
            if (C[c][0] < RHOFLOOR) {
                real d = Ub[0] - C[c][0];
                if (d > (real)0.) th = fmin(th, (Ub[0] - RHOFLOOR)/d);
            }
        if (th < 1.)
            for (int c = 0; c < 4; c++)
                for (int q = 0; q < 4; q++) {
                    U[c][q] = Ub[q] + th*(U[c][q] - Ub[q]);
                    C[c][q] = Ub[q] + th*(C[c][q] - Ub[q]);
                }

        /* p pass: common θ by bisection on the CORNER states */
        bool bad = false;
        for (int c = 0; c < 4; c++)
            if (C[c][0] < RHOFLOOR || pressure_of(C[c]) < PFLOOR) bad = true;
        if (bad) {
            real lo = 0., hi = 1.;
            for (int it = 0; it < 30; it++) {
                real mid = (real)0.5*(lo + hi);
                bool ok = true;
                for (int c = 0; c < 4; c++) {
                    real V[4];
                    for (int q = 0; q < 4; q++) V[q] = Ub[q] + mid*(C[c][q] - Ub[q]);
                    if (V[0] < RHOFLOOR || pressure_of(V) < PFLOOR) { ok = false; break; }
                }
                if (ok) lo = mid; else hi = mid;
            }
            for (int c = 0; c < 4; c++)
                for (int q = 0; q < 4; q++)
                    U[c][q] = Ub[q] + lo*(U[c][q] - Ub[q]);
        }

        if (th < 1. || bad)
            for (int c = 0; c < 4; c++) {
                int ip = i0 + (c & 1), jp = j0 + (c >> 1);
                for (int q = 0; q < 4; q++)
                    Qp[IDX_P(q,jp,ip,Np)] = U[c][q];
            }
    }
}

/* ── SSP-RK3 stage kernels (interior cells of the padded array) ─────────── */
__global__ void
rk_stage(real * __restrict__ Qout, const real * __restrict__ Qa,
         const real * __restrict__ Qb, const real * __restrict__ RHS,
         real ca, real cb, real cdt, real dt, int N, int Np)
{   /* Qout = ca*Qa + cb*(Qb + dt*RHS)   (RHS indexed on the N×N grid) */
    const int g = G_GHOST, N2 = N*N;
    for (int k = blockIdx.x*blockDim.x + threadIdx.x; k < N2; k += blockDim.x*gridDim.x) {
        int j = k/N, i = k%N;
        int jp = j+g, ip = i+g;
        for (int q = 0; q < 4; q++)
            Qout[IDX_P(q,jp,ip,Np)] = ca*Qa[IDX_P(q,jp,ip,Np)]
                + cb*(Qb[IDX_P(q,jp,ip,Np)] + dt*RHS[q*N2 + k]);
        (void)cdt;
    }
}

/* ── ICs (cell averages by 4×4 Gauss where smooth matters) ──────────────── */
__device__ void
vortex_paper_exact(real x, real y, real *ro, real *uo, real *vo, real *po)
{
    const real gg = GAMMA_V, eps = 5., pi = (real)M_PI;
    real dx = x-5., dy = y-5., r2 = dx*dx+dy*dy;
    real f  = eps/(2.*pi)*exp(0.5*(1.-r2));
    real dT = -(gg-1.)*eps*eps/(8.*gg*pi*pi)*exp(1.-r2);
    real T  = 1.+dT; if (T < 1e-6) T = 1e-6;
    *ro = pow(T, 1./(gg-1.)); *po = pow(T, gg/(gg-1.));
    *uo = -f*dy; *vo = f*dx;
}

__device__ real lmv_vphi(real r)
{ if (r<0.2) return 5.*r; if (r<0.4) return 2.-5.*r; return 0.; }
__device__ real lmv_pressure(real r, real p0)
{ if (r<0.2) return p0+12.5*r*r;
  if (r<0.4) return p0+4.*log(5.*r)+4.-20.*r+12.5*r*r;
  return p0+4.*log(2.)-2.; }
__device__ void
lmv_exact(real x, real y, real p0, real *ro, real *uo, real *vo, real *po)
{
    real dx=x-0.5, dy=y-0.5, r=sqrt(dx*dx+dy*dy);
    real vp=lmv_vphi(r), ir=(r>1e-30)?1./r:0.;
    *ro=1.; *uo=-vp*dy*ir; *vo=vp*dx*ir; *po=lmv_pressure(r,p0);
}

/* case: 0 = sod, 1 = pv, 2 = lmv */
__global__ void
ic_kernel(real * __restrict__ Qp, int N, int Np, real h, int icase, real p0lmv)
{
    const int g = G_GHOST, N2 = N*N;
    /* 4-pt Gauss nodes/weights on [-1,1] */
    const real gx[4] = {(real)-0.861136311594053, (real)-0.339981043584856,
                        (real) 0.339981043584856, (real) 0.861136311594053};
    const real gw[4] = {(real)0.347854845137454, (real)0.652145154862546,
                        (real)0.652145154862546, (real)0.347854845137454};
    for (int k = blockIdx.x*blockDim.x + threadIdx.x; k < N2; k += blockDim.x*gridDim.x) {
        int j = k/N, i = k%N;
        int jp = j+g, ip = i+g;
        real xc = (i+0.5)*h, yc = (j+0.5)*h;
        real Q[4] = {0.,0.,0.,0.};
        if (icase == 0) {                        /* circular Sod: sharp, midpoint */
            real r = sqrt((xc-0.5)*(xc-0.5) + (yc-0.5)*(yc-0.5));
            real rho = (r < 0.25) ? (real)SOD_RHOIN : 0.125;
            real p   = (r < 0.25) ? (real)SOD_PIN   : 0.1;
            Q[0] = rho; Q[3] = p/(GAMMA_V-1.);
        } else {                                  /* smooth: 4×4 Gauss cell avg */
            for (int a = 0; a < 4; a++) for (int b = 0; b < 4; b++) {
                real x = xc + 0.5*h*gx[a], y = yc + 0.5*h*gx[b];
                real r,u,v,p;
                if (icase == 1) vortex_paper_exact(x,y,&r,&u,&v,&p);
                else            lmv_exact(x,y,p0lmv,&r,&u,&v,&p);
                real w = (real)0.25*gw[a]*gw[b];
                Q[0] += w*r; Q[1] += w*r*u; Q[2] += w*r*v;
                Q[3] += w*(p/(GAMMA_V-1.) + 0.5*r*(u*u+v*v));
            }
        }
        for (int q = 0; q < 4; q++) Qp[IDX_P(q,jp,ip,Np)] = Q[q];
    }
}

/* ── error norms vs exact cell averages (4×4 Gauss) ─────────────────────── */
__global__ void
l2_err_kernel(const real * __restrict__ Qp, real * __restrict__ e_rho,
              real * __restrict__ e_vel, int N, int Np, real h,
              int icase, real p0lmv)
{
    const int g = G_GHOST, N2 = N*N;
    const real gx[4] = {(real)-0.861136311594053, (real)-0.339981043584856,
                        (real) 0.339981043584856, (real) 0.861136311594053};
    const real gw[4] = {(real)0.347854845137454, (real)0.652145154862546,
                        (real)0.652145154862546, (real)0.347854845137454};
    for (int k = blockIdx.x*blockDim.x + threadIdx.x; k < N2; k += blockDim.x*gridDim.x) {
        int j = k/N, i = k%N;
        int jp = j+g, ip = i+g;
        real xc = (i+0.5)*h, yc = (j+0.5)*h;
        real re=0., rue=0., rve=0.;
        for (int a = 0; a < 4; a++) for (int b = 0; b < 4; b++) {
            real x = xc + 0.5*h*gx[a], y = yc + 0.5*h*gx[b];
            real r,u,v,p;
            if (icase == 1) vortex_paper_exact(x,y,&r,&u,&v,&p);
            else            lmv_exact(x,y,p0lmv,&r,&u,&v,&p);
            real w = (real)0.25*gw[a]*gw[b];
            re += w*r; rue += w*r*u; rve += w*r*v;
        }
        real rho = Qp[IDX_P(0,jp,ip,Np)];
        real du  = Qp[IDX_P(1,jp,ip,Np)]/rho - rue/re;
        real dv  = Qp[IDX_P(2,jp,ip,Np)]/rho - rve/re;
        e_rho[k] = (rho-re)*(rho-re) * h*h;
        e_vel[k] = (du*du + dv*dv) * h*h;
    }
}

/* ── reductions ─────────────────────────────────────────────────────────── */
__global__ void reduce_max(const real *in, real *out, int n)
{
    __shared__ real s[BLOCK1D];
    int tid = threadIdx.x;
    real m = -1e30;
    for (int k = blockIdx.x*blockDim.x + tid; k < n; k += blockDim.x*gridDim.x)
        m = fmax(m, in[k]);
    s[tid] = m; __syncthreads();
    for (int w = blockDim.x/2; w > 0; w >>= 1) {
        if (tid < w) s[tid] = fmax(s[tid], s[tid+w]);
        __syncthreads();
    }
    if (tid == 0) out[blockIdx.x] = s[0];
}
__global__ void reduce_sum(const real *in, real *out, int n)
{
    __shared__ real s[BLOCK1D];
    int tid = threadIdx.x;
    real m = 0.;
    for (int k = blockIdx.x*blockDim.x + tid; k < n; k += blockDim.x*gridDim.x)
        m += in[k];
    s[tid] = m; __syncthreads();
    for (int w = blockDim.x/2; w > 0; w >>= 1) {
        if (tid < w) s[tid] += s[tid+w];
        __syncthreads();
    }
    if (tid == 0) out[blockIdx.x] = s[0];
}
static real gpu_max(const real *d_in, real *d_tmp, int n)
{
    reduce_max<<<64,BLOCK1D>>>(d_in, d_tmp, n);
    reduce_max<<<1,BLOCK1D>>>(d_tmp, d_tmp+64, 64);
    real r; CK(cudaMemcpy(&r, d_tmp+64, sizeof(real), cudaMemcpyDeviceToHost));
    return r;
}
static real gpu_sum(const real *d_in, real *d_tmp, int n)
{
    reduce_sum<<<64,BLOCK1D>>>(d_in, d_tmp, n);
    reduce_sum<<<1,BLOCK1D>>>(d_tmp, d_tmp+64, 64);
    real r; CK(cudaMemcpy(&r, d_tmp+64, sizeof(real), cudaMemcpyDeviceToHost));
    return r;
}


static void totals(const real *d_Uv, int N, int Np, real h, double *mass, double *ener)
{
    size_t szP = (size_t)NVAR*Np*Np*sizeof(real);
    real *hU = (real*)malloc(szP);
    cudaMemcpy(hU, d_Uv, szP, cudaMemcpyDeviceToHost);
    double m = 0., e = 0.;
    for (int j = 0; j < N; j++) for (int i = 0; i < N; i++) {
        int jp = j+G_GHOST, ip = i+G_GHOST;
        m += (double)hU[IDX_P(0,jp,ip,Np)];
        e += (double)hU[IDX_P(3,jp,ip,Np)];
    }
    *mass = m*h*h; *ener = e*h*h; free(hU);
}

/* ═══════════════════════════════════════════════════════════════════════ */
int main(int argc, char **argv)
{
    int N = (argc > 1) ? atoi(argv[1]) : 256;
    if (N % 2) { fprintf(stderr, "N must be even (2x2 positivity blocks)\n"); return 1; }
    int icase = 0;
    real lmv_eps = 0.1;
    if (argc > 2) {
        if      (!strcmp(argv[2], "pv"))  icase = 1;
        else if (!strcmp(argv[2], "lmv")) icase = 2;
        if (icase == 2 && argc > 3) lmv_eps = atof(argv[3]);
    }
    real L = (icase == 1) ? 10. : 1.;
    real h = L / N, t_end = 1.0, CFL = CFL_DEFAULT;
    real p0lmv = 1./(GAMMA_V*lmv_eps*lmv_eps) - 0.5;

    int Np = N + 2*G_GHOST;
    size_t szP = (size_t)NVAR*Np*Np*sizeof(real);
    size_t szC = (size_t)NVAR*N*N*sizeof(real);
    real *d_U, *d_U0, *d_U1, *d_RHS, *d_lam, *d_tmp, *d_e1, *d_e2;
    CK(cudaMalloc(&d_U,  szP)); CK(cudaMalloc(&d_U0, szP)); CK(cudaMalloc(&d_U1, szP));
    CK(cudaMalloc(&d_RHS, szC));
    CK(cudaMalloc(&d_lam, (size_t)N*N*sizeof(real)));
    CK(cudaMalloc(&d_tmp, 65*sizeof(real)));
    CK(cudaMalloc(&d_e1,  (size_t)N*N*sizeof(real)));
    CK(cudaMalloc(&d_e2,  (size_t)N*N*sizeof(real)));

    const char *cname = icase==0 ? "circular Sod" : icase==1 ? "paper vortex (pv)" : "Gresho (lmv)";
    printf("FV rec3 (unlimited) + HLLC + 2x2-block Zhang-Shu positivity\n");
    printf("  N=%d  h=%.5f  L=%.1f  CFL=%.2f  case=%s", N, h, L, CFL, cname);
    if (icase == 2) printf("  eps=%.4f", lmv_eps);
    printf("  [%s]\n", sizeof(real) == 8 ? "double" : "float");

    ic_kernel<<<GS_NBLK,BLOCK1D>>>(d_U, N, Np, h, icase, p0lmv);
    limit_positivity_blocks<<<GS_NBLK,BLOCK1D>>>(d_U, N, Np);
    apply_bc<<<GS_NBLK,BLOCK1D>>>(d_U, N, Np);

    CK(cudaMemset(d_RHS, 0, szC));
    compute_rhs<<<GS_NBLK,BLOCK1D>>>(d_U, d_RHS, d_lam, N, Np, h);
    CK(cudaDeviceSynchronize());
    real lam_max = gpu_max(d_lam, d_tmp, N*N);

    double m0, e0; totals(d_U, N, Np, h, &m0, &e0);
    int step = 0; real t = 0.;
    struct timespec ts0, ts1; clock_gettime(CLOCK_MONOTONIC, &ts0);
    while (t < t_end) {
        if (!(lam_max > 0.)) { fprintf(stderr, "NaN/Inf lam_max step %d t=%.5f\n", step, t); break; }
        real dt = CFL * h / lam_max;
        if (t + dt > t_end) dt = t_end - t;
        if (dt < 1e-14) { fprintf(stderr, "dt underflow step %d t=%.5f\n", step, t); break; }

        CK(cudaMemcpy(d_U0, d_U, szP, cudaMemcpyDeviceToDevice));

        /* stage 1: U1 = U0 + dt L(U0)   (RHS already holds L(U0)) */
        rk_stage<<<GS_NBLK,BLOCK1D>>>(d_U1, d_U0, d_U0, d_RHS, 0., 1., 0., dt, N, Np);
        limit_positivity_blocks<<<GS_NBLK,BLOCK1D>>>(d_U1, N, Np);
        apply_bc<<<GS_NBLK,BLOCK1D>>>(d_U1, N, Np);

        /* stage 2: U1 <- 3/4 U0 + 1/4 (U1 + dt L(U1)) */
        CK(cudaMemset(d_RHS, 0, szC));
        compute_rhs<<<GS_NBLK,BLOCK1D>>>(d_U1, d_RHS, d_lam, N, Np, h);
        rk_stage<<<GS_NBLK,BLOCK1D>>>(d_U1, d_U0, d_U1, d_RHS, 0.75, 0.25, 0., dt, N, Np);
        limit_positivity_blocks<<<GS_NBLK,BLOCK1D>>>(d_U1, N, Np);
        apply_bc<<<GS_NBLK,BLOCK1D>>>(d_U1, N, Np);

        /* stage 3: U <- 1/3 U0 + 2/3 (U1 + dt L(U1)) */
        CK(cudaMemset(d_RHS, 0, szC));
        compute_rhs<<<GS_NBLK,BLOCK1D>>>(d_U1, d_RHS, d_lam, N, Np, h);
        rk_stage<<<GS_NBLK,BLOCK1D>>>(d_U, d_U0, d_U1, d_RHS, (real)(1./3.), (real)(2./3.), 0., dt, N, Np);
        limit_positivity_blocks<<<GS_NBLK,BLOCK1D>>>(d_U, N, Np);
        apply_bc<<<GS_NBLK,BLOCK1D>>>(d_U, N, Np);

        /* RHS + lam for next step */
        CK(cudaMemset(d_RHS, 0, szC));
        compute_rhs<<<GS_NBLK,BLOCK1D>>>(d_U, d_RHS, d_lam, N, Np, h);
        CK(cudaDeviceSynchronize());
        lam_max = gpu_max(d_lam, d_tmp, N*N);
        t += dt; step++;
    }
    clock_gettime(CLOCK_MONOTONIC, &ts1);
    double wall = (ts1.tv_sec-ts0.tv_sec) + (ts1.tv_nsec-ts0.tv_nsec)*1e-9;
    printf("  Done: %d steps  t=%.6f  wall=%.2fs\n", step, t, wall);
    { int bc_; CK(cudaMemcpyFromSymbol(&bc_, g_brutal_count, sizeof(int)));
      printf("  brutal mean-fix activations: %d\n", bc_); }
    { double m1, e1c; totals(d_U, N, Np, h, &m1, &e1c);
      printf("  conservation drift: mass %.3e  energy %.3e (relative)\n",
             fabs(m1-m0)/fabs(m0), fabs(e1c-e0)/fabs(e0)); }

    if (icase != 0) {
        l2_err_kernel<<<GS_NBLK,BLOCK1D>>>(d_U, d_e1, d_e2, N, Np, h, icase, p0lmv);
        CK(cudaDeviceSynchronize());
        printf("  L2(rho) = %.6e\n", sqrt((double)gpu_sum(d_e1, d_tmp, N*N)));
        printf("  L2(vel) = %.6e\n", sqrt((double)gpu_sum(d_e2, d_tmp, N*N)));
    } else {
        /* Sod diagnostics: final min rho / min p over cells */
        real *h_U = (real*)malloc(szP);
        CK(cudaMemcpy(h_U, d_U, szP, cudaMemcpyDeviceToHost));
        double rmin = 1e30, pmin = 1e30;
        for (int j = 0; j < N; j++) for (int i = 0; i < N; i++) {
            int jp = j+G_GHOST, ip = i+G_GHOST;
            double r = h_U[IDX_P(0,jp,ip,Np)], mx = h_U[IDX_P(1,jp,ip,Np)];
            double my = h_U[IDX_P(2,jp,ip,Np)], E = h_U[IDX_P(3,jp,ip,Np)];
            double p = (GAMMA_V-1.)*(E - 0.5*(mx*mx+my*my)/r);
            if (r < rmin) rmin = r;
            if (p < pmin) pmin = p;
        }
        printf("  final min rho = %.4e   min p = %.4e\n", rmin, pmin);
        free(h_U);
    }
    return 0;
}
