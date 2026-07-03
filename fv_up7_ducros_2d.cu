/*
 * fv_up7_ducros_2d.cu
 * 2D Circular Sod — FV: UP7 + HLLC + Ducros D2 AV, SSP-RK3
 *
 * Compile:
 *   nvcc -O3 -arch=native -o fv_up7 fv_up7_ducros_2d.cu -lm          # ORDER=7
 *   nvcc -O3 -arch=native -DORDER=5 -o fv_up5 fv_up7_ducros_2d.cu -lm
 *
 * Run:
 *   ./fv_up7          # N=400
 *   ./fv_up7 600      # N=600
 *
 * Output:
 *   sod_rho.ppm   — density field  (plasma colormap)
 *   sod_p.ppm     — pressure field (viridis colormap)
 *   sod_sd.ppm    — Ducros sensor  (hot colormap)
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <cuda_runtime.h>

/* ── compile-time constants ─────────────────────────────────────────────── */
#define GAMMA_V  1.4f
#define BLOCK1D  256         /* threads per block */
#define GS_NBLK  256         /* fixed grid size for all grid-stride kernels */
#define ORDER 7

/* Upwind reconstruction order: 3, 5, 7, or 9.                             */
/* Override at compile time:  nvcc -DORDER=5 ...                            */

/* Half-width g and stencil width sw derived from ORDER (compile-time).     */
#define ORDER_G   ((ORDER + 1) / 2)   /* 2,3,4,5  */
#define ORDER_SW  (2 * ORDER_G)       /* 4,6,8,10 */

/* ── CUDA error check ───────────────────────────────────────────────────── */
#define CK(x) do { \
    cudaError_t _e = (x); \
    if (_e != cudaSuccess) { \
        fprintf(stderr, "CUDA error at %s:%d: %s\n", \
                __FILE__, __LINE__, cudaGetErrorString(_e)); \
        exit(1); \
    } \
} while (0)

/* ── Array indexing: SoA layout U[q][j][i] ─────────────────────────────── */
/* q = conserved variable (0=rho,1=rhou,2=rhov,3=E), j=row(y), i=col(x)   */
#define IDX(q,j,i,N)   ((q)*(N)*(N) + (j)*(N) + (i))
/* Padded (ghost-cell) array indexing: Np = N + 2*g, jp/ip are padded coords */
#define IDX_P(q,j,i,Np) ((q)*(Np)*(Np) + (j)*(Np) + (i))

#define NVAR 4
enum { RHO_ID = 0, U_ID = 1, V_ID = 2, P_ID = 3 };

/* ═══════════════════════════════════════════════════════════════════════════
 * Device helpers
 * ═══════════════════════════════════════════════════════════════════════════ */

/* Sound speed from primitives */
__device__ __forceinline__ float
sound_speed(float p, float rho)
{
    float cs2 = GAMMA_V * p / rho;
    return sqrtf(cs2 > 1e-14f ? cs2 : 1e-14f);
}

/* Primitive [rho,u,v,p] → conserved [rho,rhou,rhov,E] */
__device__ __forceinline__ void
p2c(float r, float u, float v, float p, float cv[4])
{
    cv[0] = r;
    cv[1] = r * u;
    cv[2] = r * v;
    cv[3] = p / (GAMMA_V - 1.f) + 0.5f * r * (u*u + v*v);
}

/* Conserved [rho,rhou,rhov,E] → primitive [rho,u,v,p] */
__device__ __forceinline__ void
c2p(const float cv[4], float wp[4])
{
    wp[0] = cv[0];
    wp[1] = cv[1] / cv[0];
    wp[2] = cv[2] / cv[0];
    wp[3] = (GAMMA_V - 1.f) * (cv[3] - 0.5f*(cv[1]*cv[1] + cv[2]*cv[2]) / cv[0]);
}

/* 2-D spectral radius — primitives [rho, u, v, p] */
__device__ __forceinline__ float
lambda(float rho, float u, float v, float p)
{
    return fabsf(u) + fabsf(v) + sound_speed(p, rho);
}

/* Directional wave speeds — primitives */
__device__ __forceinline__ float
lambda_x(float rho, float u, float v_unused, float p)
{
    return fabsf(u) + sound_speed(p, rho);
}
__device__ __forceinline__ float
lambda_y(float rho, float u_unused, float v, float p)
{
    return fabsf(v) + sound_speed(p, rho);
}

/* ── Pressure-jump shock sensor at x-face (i+½, j) ────────────────────── */
/* P_p: padded pressure slice (Np×Np); jp,ip are padded coords of cell (j,i) */
__device__ __forceinline__
float sign(float x) { return (x > 0.f) - (x < 0.f); }

__device__ __forceinline__ float
ducros_x_p(const float * __restrict__ P_p, int jp, int ip, int Np)
{
    float pL = P_p[jp * Np + ip];
    float pR = P_p[jp * Np + ip + 1];
    float dp = pR - pL;
    float psum = pR + pL;
    return dp * dp / (dp * dp + 0.1f * psum * psum);
}

/* ── Pressure-jump shock sensor at y-face (j+½, i) ────────────────────── */
__device__ __forceinline__ float
ducros_y_p(const float * __restrict__ P_p, int jp, int ip, int Np)
{
    float pB = P_p[jp * Np + ip];
    float pT = P_p[(jp + 1) * Np + ip];
    float dp = pT - pB;
    float psum = pT + pB;
    return dp * dp / (dp * dp + 0.1f * psum * psum);
}

/* ── Cell-centre sensor (max of 4 surrounding faces) — visualisation only */
__device__ __forceinline__ float
ducros_sd_p(const float * __restrict__ P_p, int jp, int ip, int Np)
{
    return fmaxf(fmaxf(ducros_x_p(P_p, jp, ip,     Np),
                       ducros_x_p(P_p, jp, ip - 1, Np)),
                 fmaxf(ducros_y_p(P_p, jp, ip,     Np),
                       ducros_y_p(P_p, jp - 1, ip, Np)));
}

/* ── Upwind reconstructions at i+½ ─────────────────────────────────────── */
/*  All stencils are left-biased for vL (upwind) and right-biased for vR.  */

/* UP3  stencil[0..3] = [i-2, i-1, i, i+1]                                */
__device__ __forceinline__ void
up3(const float s[4], float *vL, float *vR)
{
    *vL = (-s[0] + 5.0*s[1] + 2.0*s[2]) / 6.0;
    *vR = ( 2.0*s[1] + 5.0*s[2] - s[3]) / 6.0;
}

/* UP5  stencil[0..5] = [i-3, i-2, i-1, i, i+1, i+2]                      */
__device__ __forceinline__ void
up5(const float s[6], float *vL, float *vR)
{
    *vL = (  2.0*s[0] - 13.0*s[1] + 47.0*s[2]
           + 27.0*s[3] -  3.0*s[4]) / 60.0;
    *vR = ( -3.0*s[1] + 27.0*s[2] + 47.0*s[3]
           - 13.0*s[4] +  2.0*s[5]) / 60.0;
}

/* UP7  stencil[0..7] = [i-4, i-3, i-2, i-1, i, i+1, i+2, i+3]            */
__device__ __forceinline__ void
up7(const float s[8], float *vL, float *vR)
{
    *vL = (-3.0*s[0] + 25.0*s[1] - 101.0*s[2] + 319.0*s[3]
           + 214.0*s[4] - 38.0*s[5] + 4.0*s[6]) / 420.0;
    *vR = (4.0*s[1] - 38.0*s[2] + 214.0*s[3] + 319.0*s[4]
           - 101.0*s[5] + 25.0*s[6] - 3.0*s[7]) / 420.0;
}

/* UP9  stencil[0..9] = [i-5, i-4, i-3, i-2, i-1, i, i+1, i+2, i+3, i+4] */
__device__ __forceinline__ void
up9(const float s[10], float *vL, float *vR)
{
    *vL = (  3.0*s[0] -  25.0*s[1] + 101.0*s[2] - 319.0*s[3]
           + 2765.0*s[4] + 2139.0*s[5] - 533.0*s[6]
           + 97.0*s[7] -  11.0*s[8]) / 2520.0;
    *vR = (-11.0*s[1] +  97.0*s[2] - 533.0*s[3] + 2139.0*s[4]
           + 2765.0*s[5] - 319.0*s[6] + 101.0*s[7]
           -  25.0*s[8] +   3.0*s[9]) / 2520.0;
}

/* ── Compile-time dispatch: select reconstruction based on ORDER ─────────── */
/* recon_x_p: var_p is a padded per-variable slice (Np×Np).
 *   jp = j+g (padded row), ib = padded base column for the stencil.
 *   Face i-½ → ib = i  (physical base i-g maps to padded i).
 *   Face i+½ → ib = i+1.
 */
__device__ __forceinline__ void
recon_x_p(const float * __restrict__ var_p, int jp, int ib, int Np,
          float *vL, float *vR)
{
    float s[ORDER_SW];
    for (int d = 0; d < ORDER_SW; d++)
        s[d] = var_p[jp * Np + ib + d];
#if   ORDER == 3
    up3(s, vL, vR);
#elif ORDER == 5
    up5(s, vL, vR);
#elif ORDER == 9
    up9(s, vL, vR);
#else
    up7(s, vL, vR);
#endif
}

/* recon_y_p: var_p is a padded per-variable slice (Np×Np).
 *   ip = i+g (padded col), jb = padded base row.
 *   Face j-½ → jb = j.
 *   Face j+½ → jb = j+1.
 */
__device__ __forceinline__ void
recon_y_p(const float * __restrict__ var_p, int jb, int ip, int Np,
          float *vL, float *vR)
{
    float s[ORDER_SW];
    for (int d = 0; d < ORDER_SW; d++)
        s[d] = var_p[(jb + d) * Np + ip];
#if   ORDER == 3
    up3(s, vL, vR);
#elif ORDER == 5
    up5(s, vL, vR);
#elif ORDER == 9
    up9(s, vL, vR);
#else
    up7(s, vL, vR);
#endif
}

/* ── HLLC flux along face normal n=(nx,ny) ─────────────────────────────── */
/* inputs are primitives [rho, u, v, p]                                     */
__device__ void
hllc_n(const float WL[4], const float WR[4],
    float nx, float ny,
    float F[4])
{
    float rL=WL[0], uL=WL[1], vL=WL[2], pL=WL[3];
    float rR=WR[0], uR=WR[1], vR=WR[2], pR=WR[3];

    float unL = uL*nx + vL*ny;
    float unR = uR*nx + vR*ny;
    float utL =-uL*ny + vL*nx;
    float utR =-uR*ny + vR*nx;

    float cL=sound_speed(pL,rL), cR=sound_speed(pR,rR);
    float EL=pL/(GAMMA_V-1.f)+0.5f*rL*(uL*uL+vL*vL);
    float ER=pR/(GAMMA_V-1.f)+0.5f*rR*(uR*uR+vR*vR);

    float FnL[4]={rL*unL, rL*unL*uL + pL*nx, rL*unL*vL + pL*ny, (EL+pL)*unL};
    float FnR[4]={rR*unR, rR*unR*uR + pR*nx, rR*unR*vR + pR*ny, (ER+pR)*unR};

    float SL = fminf(unL-cL, unR-cR);
    float SR = fmaxf(unL+cL, unR+cR);

    if (SL >= 0.0f) {
        for (int q = 0; q < 4; q++) F[q] = FnL[q];
    } else if (SR <= 0.0f) {
        for (int q = 0; q < 4; q++) F[q] = FnR[q];
    } else {
        float dL = rL*(SL-unL), dR = rR*(SR-unR);
        float Ss = (pR - pL + dL*unL - dR*unR) / (dL - dR);

        float rK = (Ss >= 0.0f) ? rL : rR;
        float unK = (Ss >= 0.0f) ? unL : unR;
        float utK = (Ss >= 0.0f) ? utL : utR;
        float pK = (Ss >= 0.0f) ? pL : pR;
        float EK = (Ss >= 0.0f) ? EL : ER;
        float SK = (Ss >= 0.0f) ? SL : SR;
        float fact = rK * (SK - unK) / (SK - Ss);

        float Us[4];
        Us[0] = fact;
        Us[1] = fact * Ss;
        Us[2] = fact * utK;
        Us[3] = fact * (EK/rK + (Ss - unK)*(Ss + pK/(rK*(SK - unK))));

        float FK[4] = { rK*unK, rK*unK*unK+pK, rK*unK*utK, (EK+pK)*unK };
        float UcK[4] = {rK, rK*unK, rK*utK, EK};
        float Floc[4];
        for (int q = 0; q < 4; q++) Floc[q] = FK[q] + SK*(Us[q] - UcK[q]);

        F[0] = Floc[0];
        F[1] = Floc[1]*nx - Floc[2]*ny;
        F[2] = Floc[1]*ny + Floc[2]*nx;
        F[3] = Floc[3];
    }
}

/* ── HLLE flux along face normal n=(nx,ny) ─────────────────────────────── */
/* inputs are primitives [rho, u, v, p]                                     */
__device__ void
hlle_n(const float WL[4], const float WR[4],
       float nx, float ny,
       float F[4])
{
    float rL=WL[0], uL=WL[1], vL=WL[2], pL=WL[3];
    float rR=WR[0], uR=WR[1], vR=WR[2], pR=WR[3];

    float unL = uL*nx + vL*ny;
    float unR = uR*nx + vR*ny;

    float cL = sound_speed(pL, rL), cR = sound_speed(pR, rR);
    float EL = pL/(GAMMA_V-1.f) + 0.5f*rL*(uL*uL + vL*vL);
    float ER = pR/(GAMMA_V-1.f) + 0.5f*rR*(uR*uR + vR*vR);

    float UcL[4]={rL, rL*uL, rL*vL, EL};
    float UcR[4]={rR, rR*uR, rR*vR, ER};
    float FnL[4]={rL*unL, rL*unL*uL + pL*nx, rL*unL*vL + pL*ny, (EL+pL)*unL};
    float FnR[4]={rR*unR, rR*unR*uR + pR*nx, rR*unR*vR + pR*ny, (ER+pR)*unR};

    float SL = fminf(unL - cL, unR - cR);
    float SR = fmaxf(unL + cL, unR + cR);

    if (SL >= 0.0f) {
        for (int q = 0; q < 4; q++) F[q] = FnL[q];
    } else if (SR <= 0.0f) {
        for (int q = 0; q < 4; q++) F[q] = FnR[q];
    } else {
        float dS = fmaxf(SR - SL, 1e-10f);
        for (int q = 0; q < 4; q++)
            F[q] = (SR*FnL[q] - SL*FnR[q] + SR*SL*(UcR[q] - UcL[q])) / dS;
    }
}

/* ═══════════════════════════════════════════════════════════════════════════
 * Ghost-cell wall-reflecting BC kernel  (grid-stride over ALL Np×Np cells)
 *
 * Physical cells: ip in [g, N+g), jp in [g, N+g).
 * Ghost cells are filled by reflecting the nearest interior cell and
 * flipping u (q=U_ID) for x-direction ghosts, v (q=V_ID) for y-direction.
 * Corners are doubly reflected (both u and v flipped).
 * ═══════════════════════════════════════════════════════════════════════════ */
__global__ void
apply_bc(float * __restrict__ Qp, int N, int Np)
{
    int g   = (Np - N) / 2;
    int Np2 = Np * Np;
    for (int k = blockIdx.x * blockDim.x + threadIdx.x; k < Np2;
             k += blockDim.x * gridDim.x) {
        int jp = k / Np, ip = k % Np;

        bool ghost_x = (ip < g || ip >= N + g);
        bool ghost_y = (jp < g || jp >= N + g);
        if (!ghost_x && !ghost_y) continue;   /* interior cell — skip */

        /* Reflect coordinates into interior */
        int ir = ip < g ? (2*g - 1 - ip) : (ip >= N+g ? (2*(N+g) - 1 - ip) : ip);
        int jr = jp < g ? (2*g - 1 - jp) : (jp >= N+g ? (2*(N+g) - 1 - jp) : jp);

        for (int q = 0; q < NVAR; q++) {
            float v = Qp[q * Np2 + jr * Np + ir];
            if (q == U_ID && ghost_x) v = -v;
            if (q == V_ID && ghost_y) v = -v;
            Qp[q * Np2 + jp * Np + ip] = v;
        }
    }
}

/* ═══════════════════════════════════════════════════════════════════════════
 * Main RHS kernel  (grid-stride)
 * ═══════════════════════════════════════════════════════════════════════════ */
/* Qp: padded primitive state [rho,u,v,p], layout Qp[q*Np*Np + jp*Np + ip].
 * RHS, lam_out: unpadded N×N output arrays.
 * Physical cell (j,i) → padded (j+g, i+g) where g = ORDER_G.              */
__global__ void
compute_rhs(const float * __restrict__ Qp,
            float       * __restrict__ RHS,
            float       * __restrict__ lam_out,
            int N, int Np, float h)
{
    const int g   = ORDER_G;
    int N2        = N * N;
    int Np2       = Np * Np;

    /* Per-variable padded slices */
    const float* Rho = Qp + RHO_ID * Np2;
    const float* Ux  = Qp + U_ID   * Np2;
    const float* Vy  = Qp + V_ID   * Np2;
    const float* P   = Qp + P_ID   * Np2;

    /* Per-variable unpadded RHS slices */
    float* ResRho = RHS + RHO_ID * N2;
    float* ResU   = RHS + U_ID   * N2;
    float* ResV   = RHS + V_ID   * N2;
    float* ResP   = RHS + P_ID   * N2;

    for (int k = blockIdx.x * blockDim.x + threadIdx.x; k < N2;
             k += blockDim.x * gridDim.x) {
        int j = k / N, i = k % N;
        /* Padded (jp,ip) coords of the center cell */
        int jp = j + g, ip = i + g;

        /* ── center-cell primitives ──────────────────────────────────── */
        float rho0 = Rho[jp * Np + ip];
        float u0   = Ux [jp * Np + ip];
        float v0   = Vy [jp * Np + ip];
        float p0   = P  [jp * Np + ip];
        float c0   = sound_speed(p0, rho0);

        /* ── Ducros sensor at each face (ghost cells already filled) ─── */
        /* ip is the padded column of cell (j,i):
         *   face i+½ lives between ip and ip+1  → pass ip
         *   face i-½ lives between ip-1 and ip  → pass ip-1                */
        float S_Dxp = ducros_x_p(P, jp, ip,     Np);   /* face i+½ */
        float S_Dxm = ducros_x_p(P, jp, ip - 1, Np);   /* face i-½ */
        float S_Dyp = ducros_y_p(P, jp,     ip, Np);   /* face j+½ */
        float S_Dym = ducros_y_p(P, jp - 1, ip, Np);   /* face j-½ */

        /* ── x-½ face ─────────────────────────────────────────────────── */
        /* High-order: UP7 reconstruction.  Base padded col = ip - g = i.  */
        float pvLm[4], pvRm[4];
        recon_x_p(Rho, jp, i,     Np, &pvLm[0], &pvRm[0]);
        recon_x_p(Ux,  jp, i,     Np, &pvLm[1], &pvRm[1]);
        recon_x_p(Vy,  jp, i,     Np, &pvLm[2], &pvRm[2]);
        recon_x_p(P,   jp, i,     Np, &pvLm[3], &pvRm[3]);
        pvLm[0] = pvLm[0] > 1e-14f ? pvLm[0] : 1e-14f;  pvRm[0] = pvRm[0] > 1e-14f ? pvRm[0] : 1e-14f;
        pvLm[3] = pvLm[3] > 1e-14f ? pvLm[3] : 1e-14f;  pvRm[3] = pvRm[3] > 1e-14f ? pvRm[3] : 1e-14f;
        float Fm_hi[4];  hllc_n(pvLm, pvRm, 1.0f, 0.0f, Fm_hi);
        /* Low-order: cell-centred HLLE */
        float WLm[4], WRm[4];
        WLm[0] = Rho[jp * Np + ip - 1];
        WLm[1] = Ux [jp * Np + ip - 1];
        WLm[2] = Vy [jp * Np + ip - 1];
        WLm[3] = P  [jp * Np + ip - 1];
        WRm[0] = rho0; WRm[1] = u0; WRm[2] = v0; WRm[3] = p0;
        float Fm_lo[4];  hlle_n(WLm, WRm, 1.0f, 0.0f, Fm_lo);
        float Fm[4];
        for (int q = 0; q < 4; q++) Fm[q] = (1.0f - S_Dxm)*Fm_hi[q] + S_Dxm*Fm_lo[q];

        /* ── x+½ face ─────────────────────────────────────────────────── */
        /* High-order: base padded col = ip - g + 1 = i + 1.              */
        float pvLp[4], pvRp[4];
        recon_x_p(Rho, jp, i + 1, Np, &pvLp[0], &pvRp[0]);
        recon_x_p(Ux,  jp, i + 1, Np, &pvLp[1], &pvRp[1]);
        recon_x_p(Vy,  jp, i + 1, Np, &pvLp[2], &pvRp[2]);
        recon_x_p(P,   jp, i + 1, Np, &pvLp[3], &pvRp[3]);
        pvLp[0] = pvLp[0] > 1e-14f ? pvLp[0] : 1e-14f;  pvRp[0] = pvRp[0] > 1e-14f ? pvRp[0] : 1e-14f;
        pvLp[3] = pvLp[3] > 1e-14f ? pvLp[3] : 1e-14f;  pvRp[3] = pvRp[3] > 1e-14f ? pvRp[3] : 1e-14f;
        float Fp_hi[4];  hllc_n(pvLp, pvRp, 1.0f, 0.0f, Fp_hi);
        /* Low-order */
        float WLp[4], WRp[4];
        WLp[0] = rho0; WLp[1] = u0; WLp[2] = v0; WLp[3] = p0;
        WRp[0] = Rho[jp * Np + ip + 1];
        WRp[1] = Ux [jp * Np + ip + 1];
        WRp[2] = Vy [jp * Np + ip + 1];
        WRp[3] = P  [jp * Np + ip + 1];
        float Fp_lo[4];  hlle_n(WLp, WRp, 1.0f, 0.0f, Fp_lo);
        float Fp[4];
        for (int q = 0; q < 4; q++) Fp[q] = (1.0f - S_Dxp)*Fp_hi[q] + S_Dxp*Fp_lo[q];

        /* ── y-½ face ─────────────────────────────────────────────────── */
        /* High-order: base padded row = jp - g = j.                       */
        float pvLm_y[4], pvRm_y[4];
        recon_y_p(Rho, j,     ip, Np, &pvLm_y[0], &pvRm_y[0]);
        recon_y_p(Ux,  j,     ip, Np, &pvLm_y[1], &pvRm_y[1]);
        recon_y_p(Vy,  j,     ip, Np, &pvLm_y[2], &pvRm_y[2]);
        recon_y_p(P,   j,     ip, Np, &pvLm_y[3], &pvRm_y[3]);
        pvLm_y[0] = pvLm_y[0] > 1e-14f ? pvLm_y[0] : 1e-14f;  pvRm_y[0] = pvRm_y[0] > 1e-14f ? pvRm_y[0] : 1e-14f;
        pvLm_y[3] = pvLm_y[3] > 1e-14f ? pvLm_y[3] : 1e-14f;  pvRm_y[3] = pvRm_y[3] > 1e-14f ? pvRm_y[3] : 1e-14f;
        float Gm_hi[4];  hllc_n(pvLm_y, pvRm_y, 0.0f, 1.0f, Gm_hi);
        float WLm_y[4], WRm_y[4];
        WLm_y[0] = Rho[(jp - 1) * Np + ip];
        WLm_y[1] = Ux [(jp - 1) * Np + ip];
        WLm_y[2] = Vy [(jp - 1) * Np + ip];
        WLm_y[3] = P  [(jp - 1) * Np + ip];
        WRm_y[0] = rho0; WRm_y[1] = u0; WRm_y[2] = v0; WRm_y[3] = p0;
        float Gm_lo[4];  hlle_n(WLm_y, WRm_y, 0.0f, 1.0f, Gm_lo);
        float Gm[4];
        for (int q = 0; q < 4; q++) Gm[q] = (1.0f - S_Dym)*Gm_hi[q] + S_Dym*Gm_lo[q];

        /* ── y+½ face ─────────────────────────────────────────────────── */
        /* High-order: base padded row = jp - g + 1 = j + 1.              */
        float pvLp_y[4], pvRp_y[4];
        recon_y_p(Rho, j + 1, ip, Np, &pvLp_y[0], &pvRp_y[0]);
        recon_y_p(Ux,  j + 1, ip, Np, &pvLp_y[1], &pvRp_y[1]);
        recon_y_p(Vy,  j + 1, ip, Np, &pvLp_y[2], &pvRp_y[2]);
        recon_y_p(P,   j + 1, ip, Np, &pvLp_y[3], &pvRp_y[3]);
        pvLp_y[0] = pvLp_y[0] > 1e-14f ? pvLp_y[0] : 1e-14f;  pvRp_y[0] = pvRp_y[0] > 1e-14f ? pvRp_y[0] : 1e-14f;
        pvLp_y[3] = pvLp_y[3] > 1e-14f ? pvLp_y[3] : 1e-14f;  pvRp_y[3] = pvRp_y[3] > 1e-14f ? pvRp_y[3] : 1e-14f;
        float Gp_hi[4];  hllc_n(pvLp_y, pvRp_y, 0.0f, 1.0f, Gp_hi);
        float WLp_y[4], WRp_y[4];
        WLp_y[0] = rho0; WLp_y[1] = u0; WLp_y[2] = v0; WLp_y[3] = p0;
        WRp_y[0] = Rho[(jp + 1) * Np + ip];
        WRp_y[1] = Ux [(jp + 1) * Np + ip];
        WRp_y[2] = Vy [(jp + 1) * Np + ip];
        WRp_y[3] = P  [(jp + 1) * Np + ip];
        float Gp_lo[4];  hlle_n(WLp_y, WRp_y, 0.0f, 1.0f, Gp_lo);
        float Gp[4];
        for (int q = 0; q < 4; q++) Gp[q] = (1.0f - S_Dyp)*Gp_hi[q] + S_Dyp*Gp_lo[q];

        /* ── local wave speed for CFL ───────────────────────────────── */
        float lam_c = fabsf(u0) + fabsf(v0) + c0;

        /* ── conservative RHS: dU/dt = -(dF/dx + dG/dy) ────────────── */
        ResRho[k] = -(Fp[0]-Fm[0])/h - (Gp[0]-Gm[0])/h;
        ResU  [k] = -(Fp[1]-Fm[1])/h - (Gp[1]-Gm[1])/h;
        ResV  [k] = -(Fp[2]-Fm[2])/h - (Gp[2]-Gm[2])/h;
        ResP  [k] = -(Fp[3]-Fm[3])/h - (Gp[3]-Gm[3])/h;

        lam_out[k] = lam_c;
    } /* grid-stride loop */
}

/* ── sensor-only kernel (post-processing, grid-stride) ────────────────── */
/* Qp: padded primitive state.  SD: unpadded N×N output. */
__global__ void
compute_sensor(const float * __restrict__ Qp,
               float       * __restrict__ SD,
               int N, int Np)
{
    const int g  = ORDER_G;
    int N2   = N * N;
    int Np2  = Np * Np;
    const float* P = Qp + P_ID * Np2;
    for (int k = blockIdx.x * blockDim.x + threadIdx.x; k < N2;
             k += blockDim.x * gridDim.x) {
        int j = k / N, i = k % N;
        int jp = j + g, ip = i + g;
        SD[k] = ducros_sd_p(P, jp, ip, Np);
    } /* grid-stride loop */
}

/* SSP-RK3 update kernels operating on interior cells of padded state arrays.
 * U*p = padded prim state [NVAR * Np * Np]; L = unpadded conservative RHS [NVAR * N*N].
 * Grid-stride loop over N*N interior cells; ghost cells filled by apply_bc after. */
__global__ void rk3_s1(float*U1p, const float*U0p, const float*L,
                       float dt, int N, int Np)
{
    const int g = ORDER_G;
    int N2 = N*N, Np2 = Np*Np;
    for (int k = blockIdx.x*blockDim.x+threadIdx.x; k < N2; k += blockDim.x*gridDim.x) {
        int j = k/N, i = k%N;
        int idxp = (j+g)*Np + (i+g);
        float w0[4] = {U0p[0*Np2+idxp], U0p[1*Np2+idxp], U0p[2*Np2+idxp], U0p[3*Np2+idxp]};
        float cv[4]; p2c(w0[0],w0[1],w0[2],w0[3],cv);
        for(int q=0;q<4;q++) cv[q] += dt*L[q*N2+k];
        float wp[4]; c2p(cv,wp);
        wp[0] = wp[0] > 1e-14f ? wp[0] : 1e-14f;  /* rho floor */
        wp[3] = wp[3] > 1e-14f ? wp[3] : 1e-14f;  /* p floor   */
        U1p[0*Np2+idxp]=wp[0]; U1p[1*Np2+idxp]=wp[1];
        U1p[2*Np2+idxp]=wp[2]; U1p[3*Np2+idxp]=wp[3];
    }
}
__global__ void rk3_s2(float*U2p, const float*U0p, const float*U1p,
                       const float*L, float dt, int N, int Np)
{
    const int g = ORDER_G;
    int N2 = N*N, Np2 = Np*Np;
    for (int k = blockIdx.x*blockDim.x+threadIdx.x; k < N2; k += blockDim.x*gridDim.x) {
        int j = k/N, i = k%N;
        int idxp = (j+g)*Np + (i+g);
        float w0[4]={U0p[0*Np2+idxp],U0p[1*Np2+idxp],U0p[2*Np2+idxp],U0p[3*Np2+idxp]};
        float w1[4]={U1p[0*Np2+idxp],U1p[1*Np2+idxp],U1p[2*Np2+idxp],U1p[3*Np2+idxp]};
        float cv0[4]; p2c(w0[0],w0[1],w0[2],w0[3],cv0);
        float cv1[4]; p2c(w1[0],w1[1],w1[2],w1[3],cv1);
        float cv[4];
        for(int q=0;q<4;q++) cv[q]=0.75f*cv0[q]+0.25f*(cv1[q]+dt*L[q*N2+k]);
        float wp[4]; c2p(cv,wp);
        wp[0] = wp[0] > 1e-14f ? wp[0] : 1e-14f;  /* rho floor */
        wp[3] = wp[3] > 1e-14f ? wp[3] : 1e-14f;  /* p floor   */
        U2p[0*Np2+idxp]=wp[0]; U2p[1*Np2+idxp]=wp[1];
        U2p[2*Np2+idxp]=wp[2]; U2p[3*Np2+idxp]=wp[3];
    }
}
__global__ void rk3_s3(float*Up, const float*U0p, const float*U2p,
                       const float*L, float dt, int N, int Np)
{
    const int g = ORDER_G;
    int N2 = N*N, Np2 = Np*Np;
    for (int k = blockIdx.x*blockDim.x+threadIdx.x; k < N2; k += blockDim.x*gridDim.x) {
        int j = k/N, i = k%N;
        int idxp = (j+g)*Np + (i+g);
        float w0[4]={U0p[0*Np2+idxp],U0p[1*Np2+idxp],U0p[2*Np2+idxp],U0p[3*Np2+idxp]};
        float w2[4]={U2p[0*Np2+idxp],U2p[1*Np2+idxp],U2p[2*Np2+idxp],U2p[3*Np2+idxp]};
        float cv0[4]; p2c(w0[0],w0[1],w0[2],w0[3],cv0);
        float cv2[4]; p2c(w2[0],w2[1],w2[2],w2[3],cv2);
        float cv[4];
        for(int q=0;q<4;q++) cv[q]=(1.f/3.f)*cv0[q]+(2.f/3.f)*(cv2[q]+dt*L[q*N2+k]);
        float wp[4]; c2p(cv,wp);
        wp[0] = wp[0] > 1e-14f ? wp[0] : 1e-14f;  /* rho floor */
        wp[3] = wp[3] > 1e-14f ? wp[3] : 1e-14f;  /* p floor   */
        Up[0*Np2+idxp]=wp[0]; Up[1*Np2+idxp]=wp[1];
        Up[2*Np2+idxp]=wp[2]; Up[3*Np2+idxp]=wp[3];
    }
}

/* ── max reduction (grid-stride two-pass) ──────────────────────────────── */
__global__ void reduce_max(const float*in, float*out, int n){
    extern __shared__ float sm[];
    int tid = threadIdx.x;
    float v = 0.;
    for (int k = blockIdx.x*blockDim.x + tid; k < n; k += blockDim.x*gridDim.x)
        v = fmaxf(v, in[k]);
    sm[tid] = v; __syncthreads();
    for (int s = BLOCK1D/2; s > 0; s >>= 1){
        if (tid < s) sm[tid] = fmaxf(sm[tid], sm[tid+s]);
        __syncthreads();
    }
    if (tid == 0) out[blockIdx.x] = sm[0];
}

static float gpu_max(const float*d_in, float*d_tmp, int n){
    reduce_max<<<GS_NBLK,BLOCK1D,BLOCK1D*sizeof(float)>>>(d_in,d_tmp,n);
    reduce_max<<<1,BLOCK1D,BLOCK1D*sizeof(float)>>>(d_tmp,d_tmp,GS_NBLK);
    float v; CK(cudaMemcpy(&v,d_tmp,sizeof(float),cudaMemcpyDeviceToHost));
    return v;
}

/* ═══════════════════════════════════════════════════════════════════════════
 * IC kernel  (grid-stride)
 * ═══════════════════════════════════════════════════════════════════════════ */
/* ─ IC kernel writes into interior cells of the padded state array ──────── */
__global__ void
ic_kernel(float* Up, int N, int Np, float h)
{
    const int g = ORDER_G;
    int N2  = N * N;
    int Np2 = Np * Np;
    for (int k = blockIdx.x * blockDim.x + threadIdx.x; k < N2;
             k += blockDim.x * gridDim.x) {
        int j = k / N, i = k % N;
        int jp = j + g, ip = i + g;

        float xc=(i+0.5f)*h, yc=(j+0.5f)*h;
        float r=sqrtf((xc-0.5f)*(xc-0.5f)+(yc-0.5f)*(yc-0.5f));
        float delta=1.0f*h;
        float phi=0.5f*(1.f+tanhf((0.25f-r)/delta));
        /* store primitives [rho, u, v, p] */
        Up[0*Np2 + jp*Np + ip] = 0.125f + 10.875f*phi;  /* rho */
        Up[1*Np2 + jp*Np + ip] = 0.f;                   /* u   */
        Up[2*Np2 + jp*Np + ip] = 0.f;                   /* v   */
        Up[3*Np2 + jp*Np + ip] = 0.1f + 9.9f*phi;      /* p   */
    } /* grid-stride loop */
}

/* ═══════════════════════════════════════════════════════════════════════════
 * PPM output
 * ═══════════════════════════════════════════════════════════════════════════ */

/* Plasma colormap (256 entries), sampled from matplotlib */
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

/* Viridis colormap (256 entries) */
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

/* Hot colormap (black->red->yellow->white) */
static void hot_rgb(float t, unsigned char *r, unsigned char *g, unsigned char *b)
{
    if      (t < 1./3.) { *r=(unsigned char)(t*3.*255); *g=0; *b=0; }
    else if (t < 2./3.) { *r=255; *g=(unsigned char)((t-1./3.)*3.*255); *b=0; }
    else                { *r=255; *g=255; *b=(unsigned char)((t-2./3.)*3.*255); }
}

static void
write_ppm(const char *fname, const float *data, int N,
          float vmin, float vmax, int cmap)
{
    /* cmap: 0=plasma, 1=viridis, 2=hot */
    FILE *f = fopen(fname, "wb");
    if (!f) { fprintf(stderr, "Cannot open %s\n", fname); return; }
    fprintf(f, "P6\n%d %d\n255\n", N, N);

    unsigned char pix[3];
    for (int jj = N-1; jj >= 0; jj--) {          /* flip y: j=0 at bottom */
        for (int ii = 0; ii < N; ii++) {
            float t = (data[jj*N + ii] - vmin) / (vmax - vmin + 1e-30f);
            if (!(t > 0.0f)) t = 0.0f;   /* also catches NaN */
            if (t > 1.0f) t = 1.0f;
            int idx = (int)(t * 255.0f);
            if (idx < 0) idx = 0; if (idx > 255) idx = 255;
            if (cmap == 0)      { pix[0]=PLASMA[idx][0]; pix[1]=PLASMA[idx][1]; pix[2]=PLASMA[idx][2]; }
            else if (cmap == 1) { pix[0]=VIRIDIS[idx][0];pix[1]=VIRIDIS[idx][1];pix[2]=VIRIDIS[idx][2];}
            else                { hot_rgb(t, &pix[0], &pix[1], &pix[2]); }
            fwrite(pix, 1, 3, f);
        }
    }
    fclose(f);
    printf("  Saved %s\n", fname);
}

/* Blue-white-red diverging colormap (for vorticity) */
static void bwr_rgb(float t, unsigned char *r, unsigned char *g, unsigned char *b)
{
    if (t < 0.5) {
        *r = (unsigned char)(t * 2.0 * 255);
        *g = (unsigned char)(t * 2.0 * 255);
        *b = 255;
    } else {
        *r = 255;
        *g = (unsigned char)((1.0 - t) * 2.0 * 255);
        *b = (unsigned char)((1.0 - t) * 2.0 * 255);
    }
}

/* Compute vorticity (dv/dx - du/dy) and velocity magnitude on GPU (grid-stride) */
/* Qp: padded state; g=ORDER_G ≥ 3 so the ±3 stencil is safely inside ghost region. */
__global__ void
compute_derived(const float * __restrict__ Qp,
                float       * __restrict__ vort,
                float       * __restrict__ vmag,
                int N, int Np, float h)
{
    const int g  = ORDER_G;     /* always ≥ 2; for ORDER=7 g=4 ≥ 3 ✓ */
    int N2   = N * N;
    int Np2  = Np * Np;
    const float* Ux = Qp + U_ID * Np2;
    const float* Vy = Qp + V_ID * Np2;

    for (int k = blockIdx.x * blockDim.x + threadIdx.x; k < N2;
             k += blockDim.x * gridDim.x) {
        int j = k / N, i = k % N;
        int jp = j + g, ip = i + g;

        /* 6th-order skew-symmetric stencil */
        const float coeff[6] = {-1.f, 9.f, -45.f, 45.f, -9.f, 1.f};
        const int    off[6]  = {-3, -2, -1,  1,  2,  3};

        float dv_dx = 0.f, du_dy = 0.f;
        for (int d = 0; d < 6; d++) {
            dv_dx += coeff[d] * Vy[jp * Np + (ip + off[d])];
            du_dy += coeff[d] * Ux[(jp + off[d]) * Np + ip];
        }
        vort[k] = (dv_dx - du_dy) / (60.0f * h);
        vmag[k] = sqrtf(Ux[jp*Np+ip]*Ux[jp*Np+ip] + Vy[jp*Np+ip]*Vy[jp*Np+ip]);
    } /* grid-stride loop */
}

/* Write a 2N×2N composite PPM.  Layout (physical coordinates):
 *   top-left  | top-right
 *   ----------+----------
 *   bot-left  | bot-right
 * cmap: 0=plasma, 1=viridis, 2=hot, 3=bwr
 */
static void
write_ppm_2x2(const char *fname,
              const float *tl, float tl_min, float tl_max, int tl_cmap,
              const float *tr, float tr_min, float tr_max, int tr_cmap,
              const float *bl, float bl_min, float bl_max, int bl_cmap,
              const float *br, float br_min, float br_max, int br_cmap,
              int N)
{
    FILE *f = fopen(fname, "wb");
    if (!f) { fprintf(stderr, "Cannot open %s\n", fname); return; }
    fprintf(f, "P6\n%d %d\n255\n", 2*N, 2*N);

    unsigned char pix[3];
    for (int r = 0; r < 2*N; r++) {        /* r=0 is the top of the image */
        int phys_j = 2*N - 1 - r;          /* y-flip: phys_j=0 → bottom  */
        for (int c = 0; c < 2*N; c++) {
            const float *panel;
            float vmin, vmax;
            int cmap, sub_j, sub_i;

            if (phys_j >= N) {              /* top row of panels */
                sub_j = phys_j - N;
                if (c < N) { panel=tl; vmin=tl_min; vmax=tl_max; cmap=tl_cmap; sub_i=c;   }
                else        { panel=tr; vmin=tr_min; vmax=tr_max; cmap=tr_cmap; sub_i=c-N; }
            } else {                        /* bottom row of panels */
                sub_j = phys_j;
                if (c < N) { panel=bl; vmin=bl_min; vmax=bl_max; cmap=bl_cmap; sub_i=c;   }
                else        { panel=br; vmin=br_min; vmax=br_max; cmap=br_cmap; sub_i=c-N; }
            }

            float t = (panel[sub_j*N + sub_i] - vmin) / (vmax - vmin + 1e-30f);
            if (!(t > 0.0f)) t = 0.0f;   /* also catches NaN */
            if (t > 1.0f) t = 1.0f;
            int idx = (int)(t * 255.0f);
            if (idx < 0) idx = 0; if (idx > 255) idx = 255;
            if      (cmap == 0) { pix[0]=PLASMA[idx][0];  pix[1]=PLASMA[idx][1];  pix[2]=PLASMA[idx][2];  }
            else if (cmap == 1) { pix[0]=VIRIDIS[idx][0]; pix[1]=VIRIDIS[idx][1]; pix[2]=VIRIDIS[idx][2]; }
            else if (cmap == 2) { hot_rgb(t, &pix[0], &pix[1], &pix[2]); }
            else                { bwr_rgb(t, &pix[0], &pix[1], &pix[2]); }
            fwrite(pix, 1, 3, f);
        }
    }
    fclose(f);
    printf("  Saved %s\n", fname);
}

/* ═══════════════════════════════════════════════════════════════════════════
 * main
 * ═══════════════════════════════════════════════════════════════════════════ */
int main(int argc, char **argv)
{
    int    N   = (argc > 1) ? atoi(argv[1]) : 400;
    float CFL  = 0.40f, t_end = 1.0f, h = 1.0f / N;
    const int g = ORDER_G;
    int Np      = N + 2 * g;   /* padded grid width (ghost cells each side) */

    /* GPU info */
    int dev; cudaDeviceProp prop;
    CK(cudaGetDevice(&dev));
    CK(cudaGetDeviceProperties(&prop, dev));
    printf("========================================================\n");
    printf("  Device  : %s\n", prop.name);
    printf("  UP%d+HLLC+Ducros FV 2D  --  Circular Sod\n", ORDER);
    printf("  n=%dx%d  Np=%d  h=%.5f  CFL=%.2f  order=%d\n",
           N, N, Np, h, CFL, ORDER);
    printf("========================================================\n");

    /* allocate */
    size_t sz1  = (size_t)N  * N  * sizeof(float);
    size_t sz4  = (size_t)4  * sz1;           /* unpadded RHS / scratch     */
    size_t sz4p = (size_t)NVAR * Np * Np * sizeof(float); /* padded state    */
    float *d_U, *d_U0, *d_U1, *d_U2, *d_RHS, *d_lam, *d_tmp, *d_SD, *d_vort, *d_vmag;
    CK(cudaMalloc(&d_U,    sz4p));
    CK(cudaMalloc(&d_U0,   sz4p));
    CK(cudaMalloc(&d_U1,   sz4p));
    CK(cudaMalloc(&d_U2,   sz4p));
    CK(cudaMalloc(&d_RHS,  sz4));    /* conservative RHS — unpadded N×N */
    CK(cudaMalloc(&d_lam,  sz1));
    CK(cudaMalloc(&d_SD,   sz1));
    CK(cudaMalloc(&d_vort, sz1));
    CK(cudaMalloc(&d_vmag, sz1));
    /* tmp for max reduction: one partial result per GS_NBLK block */
    CK(cudaMalloc(&d_tmp, GS_NBLK * sizeof(float)));

    float *h_rho  = (float*)malloc(sz1);
    float *h_p    = (float*)malloc(sz1);
    float *h_sd   = (float*)malloc(sz1);
    float *h_rhou = (float*)malloc(sz1);
    float *h_rhov = (float*)malloc(sz1);
    float *h_E    = (float*)malloc(sz1);
    float *h_vort = (float*)malloc(sz1);
    float *h_vmag = (float*)malloc(sz1);

    /* ── helper: download fields, compute pressure, write 3 PPMs ──────── */
    /* color scale fixed to IC extremes for temporal consistency           */
    const float RHO_MIN = 0.125, RHO_MAX = 10.1;

    /* local lambda via nested function trick — use a macro instead */
    /* Download N×N interior from padded array using cudaMemcpy2D */
#define DL_INTERIOR(h_dst, q) \
    CK(cudaMemcpy2D(h_dst, N*sizeof(float), \
                    d_U + (size_t)(q)*Np*Np + g*Np + g, Np*sizeof(float), \
                    N*sizeof(float), N, cudaMemcpyDeviceToHost))
#define WRITE_FRAME(frame_idx) do { \
    /* Download interior cells for each primitive variable */\
    DL_INTERIOR(h_rho,  RHO_ID); \
    DL_INTERIOR(h_rhou, U_ID); \
    DL_INTERIOR(h_rhov, V_ID); \
    DL_INTERIOR(h_p,    P_ID); \
    compute_sensor <<<GS_NBLK, BLOCK1D>>>(d_U, d_SD,   N, Np); \
    compute_derived<<<GS_NBLK, BLOCK1D>>>(d_U, d_vort, d_vmag, N, Np, h); \
    CK(cudaDeviceSynchronize()); \
    CK(cudaMemcpy(h_sd,   d_SD,   sz1, cudaMemcpyDeviceToHost)); \
    CK(cudaMemcpy(h_vort, d_vort, sz1, cudaMemcpyDeviceToHost)); \
    CK(cudaMemcpy(h_vmag, d_vmag, sz1, cudaMemcpyDeviceToHost)); \
    /* h_p already holds pressure; h_rhou/rhov hold u/v directly */ \
    /* symmetric vorticity color range */ \
    float _wmax = 1e-10f; \
    for (int _k = 0; _k < N*N; _k++) { \
        float _aw = h_vort[_k] < 0.0 ? -h_vort[_k] : h_vort[_k]; \
        if (_aw > _wmax) _wmax = _aw; \
    } \
    /* velocity magnitude range 0..max */ \
    float _vmax = 1e-10f; \
    for (int _k = 0; _k < N*N; _k++) \
        if (h_vmag[_k] > _vmax) _vmax = h_vmag[_k]; \
    char _fn[64]; \
    sprintf(_fn, "sod_%04d.ppm", (frame_idx)); \
    write_ppm_2x2(_fn, \
        h_rho,  RHO_MIN, RHO_MAX, 0, \
        h_vort, -_wmax,  _wmax,   3, \
        h_sd,   0.0,     1.0,     2, \
        h_vmag, 0.0,     _vmax,   1, \
        N); \
} while(0)

    /* IC — write interior, then fill ghost cells */
    ic_kernel <<<GS_NBLK, BLOCK1D>>>(d_U, N, Np, h);
    apply_bc  <<<GS_NBLK, BLOCK1D>>>(d_U, N, Np);
    CK(cudaDeviceSynchronize());

    /* frame 0 at t=0 */
    WRITE_FRAME(0);

    /* time loop */
    const int  N_FRAMES  = 10;
    int    frame     = 1;
    float t_next    = t_end / N_FRAMES;    /* next frame target time */
    int    step = 0;
    float t    = 0.0f;
    struct timespec ts0, ts1;
    clock_gettime(CLOCK_MONOTONIC, &ts0);

    /* precompute RHS and lam_max for the first step */
    compute_rhs<<<GS_NBLK, BLOCK1D>>>(d_U, d_RHS, d_lam, N, Np, h);
    CK(cudaDeviceSynchronize());
    float lam_max = gpu_max(d_lam, d_tmp, N*N);

    while (t < t_end) {
        /* lam_max is always computed on d_U at end of previous step (or IC) */
        if (!(lam_max > 0.0f)) { fprintf(stderr, "NaN/Inf lam_max at step %d t=%.5f\n", step, t); break; }
        float dt = CFL * h / lam_max;
        /* clamp to next frame boundary or t_end, whichever is closer */
        float t_target = (t_next < t_end) ? t_next : t_end;
        if (t + dt > t_target) dt = t_target - t;
        if (dt < 1e-14f) {
            fprintf(stderr, "  dt underflow at step %d t=%.6f lam_max=%.3e — writing blow-up frame\n",
                    step, t, lam_max);
            WRITE_FRAME(99);
            break;
        }

        /* SSP-RK3 */
        CK(cudaMemcpy(d_U0, d_U, sz4p, cudaMemcpyDeviceToDevice));

        /* stage 1: d_RHS already holds L(d_U) from previous iteration */
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

        /* compute RHS on updated d_U: health-check BEFORE frame write,
           and d_RHS / lam_max are ready for stage 1 of the next step     */
        compute_rhs<<<GS_NBLK, BLOCK1D>>>(d_U, d_RHS, d_lam, N, Np, h);
        lam_max = gpu_max(d_lam, d_tmp, N*N);
        if (!(lam_max > 0.0f)) { fprintf(stderr, "blow-up after step %d t=%.5f\n", step+1, t+dt); break; }

        t += dt; step++;

        /* check frame boundary */
        if (t >= t_next - 1e-12f && frame <= N_FRAMES) {
            clock_gettime(CLOCK_MONOTONIC, &ts1);
            float el = (ts1.tv_sec-ts0.tv_sec)+(ts1.tv_nsec-ts0.tv_nsec)*1e-9;
            printf("  frame %2d/%d  step %5d  t=%.5f  elapsed=%.1fs\n",
                   frame, N_FRAMES, step, t, el);
            fflush(stdout);
            WRITE_FRAME(frame);
            frame++;
            t_next = frame * t_end / N_FRAMES;
        } else if (step % 50 == 0) {
            clock_gettime(CLOCK_MONOTONIC, &ts1);
            float el = (ts1.tv_sec-ts0.tv_sec)+(ts1.tv_nsec-ts0.tv_nsec)*1e-9;
            float eta = (t > 0) ? (t_end-t)/t*el : 0;
            printf("  step %5d  t=%.5f  dt=%.3e  elapsed=%.1fs  ETA=%.0fs\n",
                   step, t, dt, el, eta);
            fflush(stdout);
        }
    }
    CK(cudaDeviceSynchronize());
    clock_gettime(CLOCK_MONOTONIC, &ts1);
    float wall = (ts1.tv_sec-ts0.tv_sec)+(ts1.tv_nsec-ts0.tv_nsec)*1e-9;
    printf("  Done: %d steps,  t=%.6f,  wall=%.2fs\n", step, t, wall);
    printf("  Frames written: sod_0000.ppm .. sod_%04d.ppm  (2x2 composite)\n", N_FRAMES);

    /* cleanup */
    free(h_rho); free(h_p); free(h_sd);
    free(h_rhou); free(h_rhov); free(h_E);
    free(h_vort); free(h_vmag);
    cudaFree(d_U); cudaFree(d_U0); cudaFree(d_U1); cudaFree(d_U2);
    cudaFree(d_RHS); cudaFree(d_lam); cudaFree(d_tmp); cudaFree(d_SD);
    cudaFree(d_vort); cudaFree(d_vmag);
    return 0;
}
