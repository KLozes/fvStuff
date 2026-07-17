/*
 * dgsem_lobatto_2d.cu
 * 2D Circular Sod — nodal Discontinuous Galerkin Spectral Element Method (DGSEM)
 * on Legendre–Gauss–Lobatto (LGL / "Lobatto") solution nodes.
 * Entropy-stable volume (Chandrashekar EC flux differencing) · HLLC interface
 * flux · Ducros shock-sensor artificial viscosity · SSP-RK3 · positivity floor.
 *
 * This is a GENUINE discontinuous-Galerkin scheme (not the spectral-difference
 * scheme in fr_lobatto_2d.cu).  The polynomial degree is p = ORDER (default 3),
 * so the default build is a P3 DG method.
 *
 * Compile (polynomial degree p=1..3):
 *   nvcc -O3 -arch=native --expt-relaxed-constexpr -DORDER=3 -o dgsem dgsem_lobatto_2d.cu -lm
 *
 * Run:
 *   ./dgsem          # N=50 elements/side
 *   ./dgsem 80       # N=80 elements/side
 *
 * ── Formulation ────────────────────────────────────────────────────────────
 * Collocation DGSEM on the reference element [-1,1] with LGL nodes {xi_i}
 * (which INCLUDE the endpoints ±1).  GLL quadrature makes the mass matrix
 * diagonal, M_ii = (h/2) w_i.  The strong form of  u_t + f_x + g_y = 0  is
 *
 *   du_ij/dt = -(2/h) [  Σ_k D_ik f(u_kj)            (x volume, differentiation matrix)
 *                      + Σ_k D_jk g(u_ik)            (y volume)
 *                      + (1/w_i)( δ_{i,p}(f*_R - f_pj) - δ_{i,0}(f*_L - f_0j) )
 *                      + (1/w_j)( δ_{j,p}(g*_T - g_ip) - δ_{j,0}(g*_B - g_i0) ) ]
 *
 *   D_ik = ℓ_k'(xi_i)         polynomial differentiation matrix on LGL nodes
 *   w_i                       GLL quadrature weights (Σ w = 2)
 *   f*, g*                    HLLC numerical flux, evaluated pointwise at each of
 *                             the (p+1) LGL nodes shared across the element face
 *   f_0j, f_pj, …             physical flux at the boundary solution nodes
 *
 * Because neighbouring elements share the same LGL nodes along a face, the
 * numerical flux couples matching face nodes one-to-one — the classic DGSEM
 * structure.  Boundaries are reflecting walls (as in fr_lobatto_2d.cu).
 *
 * ── Element-internal shock capturing (no limiter) ──────────────────────────
 * High-order DG needs shock capturing.  Instead of a slope limiter this solver
 * uses an isotropic artificial viscosity  ∂tW += ∇·(ν ∇U)  gated by a Ducros-
 * type compression sensor on the velocity divergence:
 *     θ = (∇·u)²/((∇·u)² + K c²/dx²),   dx = h/(2p+1),
 *     ν_e = C_av · θ_e · (h/(2p+1)) · λ   (element-wise, θ_e = max node).
 * θ→1 at a shock, θ→0 in smooth flow, so full P3 accuracy is retained away
 * from shocks.  A Zhang–Shu positivity floor is the final safeguard.
 * Knobs: -DAV_CAV, -DAV_KSENSOR, -DAV_ON=0.
 *
 * State layout (conservative):  Q[var * NEL * NN2 + el * NN2 + jn * NN + in]
 *   var ∈ {rho=0, rhou=1, rhov=2, E=3},  NEL = N*N,  NN = ORDER+1,  NN2 = NN*NN
 *   jn = node index in y,  in = node index in x,  both in [0, NN)
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

#ifndef ORDER
#  define ORDER 2          /* polynomial degree p: 1..3 */
#endif
#if ORDER > 3
#  error "ORDER > 3 is not supported: extend the LGL node / GLL weight tables."
#endif
#define NN      (ORDER + 1)   /* solution nodes per direction per element */
#define NN2     (NN * NN)     /* solution nodes per element               */

/* RHS kernel maps one thread per solution node and packs EPB elements per
 * block (≈ one warp), keeping each element's working set in shared memory. */
#ifndef EPB
#  define EPB    (32 / NN2)   /* elements per block (p3→2, p2→3, p1→8) */
#endif

/* ── Shock-sensor artificial viscosity (primary shock capturing) ───────────
 * Isotropic AV  ∂tW += ∇·(ν ∇U),  applied element-wise, gated by a Ducros-type
 * compression sensor keyed to the velocity divergence ∇·u:
 *     θ = (∇·u)² / ( (∇·u)² + K c²/dx² ),   dx = h/(2p+1),  K = AV_KSENSOR.
 * θ→1 where the compression rate rivals the acoustic rate c/dx (a shock),
 * θ→0 in smooth flow.  Element viscosity  ν_e = C_av · θ_e · (h/(2p+1)) · λ,
 * with θ_e = max over the element's nodes.  Zhang–Shu positivity floor is the
 * final safeguard.  Knobs: -DAV_CAV=, -DAV_KSENSOR=, -DAV_ON=0.              */
/* ── Entropy-stable volume operator ────────────────────────────────────────
 * EC_FLUX=1 replaces the collocation volume derivative Σ_k D_ik f_k with the
 * flux-differencing form  2 Σ_k D_ik F#(U_i,U_k)  using Chandrashekar's
 * entropy-conservative two-point flux.  This makes the volume operator
 * entropy-conservative (removes aliasing-driven blow-up); with the dissipative
 * HLLC interface the scheme is entropy-stable.  Set -DEC_FLUX=0 for the plain
 * collocation volume term.                                                    */
#ifndef EC_FLUX
#  define EC_FLUX 1
#endif


/* ── Circular-Sod blast strength ───────────────────────────────────────────
 * Ambient (outside): rho=0.125, p=0.1.  Inside: rho=RHO_HI, p=P_HI.
 * Baseline (S=1) is the 100:1 pressure / 88:1 density jump.  Scale both by a
 * strength multiplier to make a "bigger" shock: -DP_HI=<f> -DRHO_HI=<f>.      */
#ifndef P_HI
#  define P_HI   10.0f
#endif
#ifndef RHO_HI
#  define RHO_HI 11.0f
#endif

#ifndef AV_ON
#  define AV_ON 1
#endif
#ifndef AV_CAV
#  define AV_CAV .25f          /* artificial-viscosity strength */
#endif
#ifndef AV_KSENSOR
#  define AV_KSENSOR 0.1f      /* sensor threshold constant K   */
#endif

#define NVAR    4
enum { RHO_ID = 0, U_ID = 1, V_ID = 2, P_ID = 3 };

/* State index: var, element flat index el = ej*N+ei, node (jn, in) */
#define SIDX(var, el, jn, in, NEL)  ((var)*(NEL)*NN2 + (el)*NN2 + (jn)*NN + (in))

/* ── CUDA error check ───────────────────────────────────────────────────── */
#define CK(x) do { \
    cudaError_t _e = (x); \
    if (_e != cudaSuccess) { \
        fprintf(stderr, "CUDA error at %s:%d: %s\n", \
                __FILE__, __LINE__, cudaGetErrorString(_e)); \
        exit(1); \
    } \
} while(0)

/* ── DG operators in constant memory ───────────────────────────────────── */
__constant__ float c_D   [4][4];  /* differentiation matrix D_ik = ℓ_k'(xi_i) */
__constant__ float c_w   [4];     /* quadrature weights (sum = 2)             */
__constant__ float c_winv[4];     /* 1 / w_i  (surface lift)                  */
__constant__ float c_xi  [4];     /* solution node coordinates                */

/* ── LGL nodes on [-1,1] for p=1..3 (NN=2..4) ─────────────────────────── */
static const double lgl_xi[3][4] = {
    /* p=1 */ {-1.0, 1.0},
    /* p=2 */ {-1.0, 0.0, 1.0},
    /* p=3 */ {-1.0,-0.4472135954999579, 0.4472135954999579, 1.0},
};

/* ── GLL quadrature weights on [-1,1] for p=1..3 (sum = 2) ─────────────── */
static const double gll_w[3][4] = {
    /* p=1, n=2 */ {1.0, 1.0},
    /* p=2, n=3 */ {0.33333333333333333, 1.33333333333333333, 0.33333333333333333},
    /* p=3, n=4 */ {0.16666666666666667, 0.83333333333333333,
                    0.83333333333333333, 0.16666666666666667},
};


/* Barycentric differentiation matrix D_ik = ℓ_k'(x_i) on nodes x[0..n-1].   */
static void build_D(int n, const double *x, double D[4][4])
{
    double b[4];
    for (int k = 0; k < n; k++) {
        b[k] = 1.0;
        for (int m = 0; m < n; m++)
            if (m != k) b[k] /= (x[k] - x[m]);
    }
    for (int i = 0; i < n; i++) {
        double diag = 0.0;
        for (int k = 0; k < n; k++) {
            if (k != i) {
                D[i][k] = (b[k] / b[i]) / (x[i] - x[k]);
                diag   -= D[i][k];
            }
        }
        D[i][i] = diag;               /* negative-sum row property */
    }
}

/* ═══════════════════════════════════════════════════════════════════════════
 * Device helpers  (identical to the FV / FR solvers)
 * ═══════════════════════════════════════════════════════════════════════════ */
__device__ __forceinline__ float sound_speed(float p, float rho)
{
    float cs2 = GAMMA_V * p / rho;
    return sqrtf(cs2 > 1e-14f ? cs2 : 1e-14f);
}

__device__ __forceinline__ void p2c(float r, float u, float v, float p, float cv[4])
{
    cv[0] = r;
    cv[1] = r * u;
    cv[2] = r * v;
    cv[3] = p / (GAMMA_V - 1.f) + 0.5f * r * (u*u + v*v);
}

__device__ __forceinline__ void c2p(const float cv[4], float wp[4])
{
    wp[0] = cv[0];
    wp[1] = cv[1] / cv[0];
    wp[2] = cv[2] / cv[0];
    wp[3] = (GAMMA_V - 1.f) * (cv[3] - 0.5f*(cv[1]*cv[1] + cv[2]*cv[2]) / cv[0]);
}

__device__ __forceinline__ void sanitize_prim(float W[4]);

__device__ __forceinline__ float pressure_from_cons(const float U[4])
{
    float rho = fmaxf(U[0], 1e-14f);
    float ke  = 0.5f * (U[1]*U[1] + U[2]*U[2]) / rho;
    return (GAMMA_V - 1.f) * (U[3] - ke);
}

__device__ __forceinline__ void sanitize_cons(float U[4])
{
    U[0] = fmaxf(U[0], 1e-14f);
    float p = pressure_from_cons(U);
    if (p < 1e-14f)
        U[3] = 1e-14f/(GAMMA_V-1.f) + 0.5f*(U[1]*U[1] + U[2]*U[2]) / U[0];
}

__device__ __forceinline__ void cons_to_prim_sane(const float U[4], float W[4])
{
    float Us[4] = { U[0], U[1], U[2], U[3] };
    sanitize_cons(Us);
    c2p(Us, W);
    sanitize_prim(W);
}

/* ── Euler flux in direction (nx,ny) from primitives ────────────────────── */
__device__ __forceinline__ void euler_flux_n(float rho, float u, float v, float p,
                                             float nx, float ny, float F[4])
{
    float un = u*nx + v*ny;
    float E  = p/(GAMMA_V-1.f) + 0.5f*rho*(u*u+v*v);
    F[0] = rho * un;
    F[1] = rho * un * u + p * nx;
    F[2] = rho * un * v + p * ny;
    F[3] = (E + p) * un;
}

/* ── Stable logarithmic mean  (aL-aR)/(ln aL - ln aR) ────────────────────── */
__device__ __forceinline__ float logmean(float aL, float aR)
{
    float d  = aL/aR;
    float f  = (d-1.f)/(d+1.f);
    float u2 = f*f;
    float FF = (u2 < 1e-4f)
             ? (1.f + u2*(1.f/3.f + u2*(1.f/5.f + u2*(1.f/7.f))))
             : (logf(d)/(2.f*f));
    return (aL+aR)/(2.f*FF);
}

/* ── Chandrashekar (2013) entropy-conservative two-point flux ────────────── *
 * From primitive states (rho,u,v,p).  dir=0 → x-flux F#,  dir=1 → y-flux G#.
 * Symmetric & consistent: F#(W,W) = physical flux f(W).  Used in the
 * flux-differencing volume term to make the DG operator entropy-conservative. */
__device__ __forceinline__ void ec_flux(const float WL[4], const float WR[4],
                                        int dir, float F[4])
{
    float rL=WL[0], uL=WL[1], vL=WL[2], pL=WL[3];
    float rR=WR[0], uR=WR[1], vR=WR[2], pR=WR[3];
    float bL = 0.5f*rL/pL,  bR = 0.5f*rR/pR;         /* beta = rho/(2p)     */
    float r_ln = logmean(rL, rR);
    float b_ln = logmean(bL, bR);
    float r_av = 0.5f*(rL+rR);
    float b_av = 0.5f*(bL+bR);
    float u_av = 0.5f*(uL+uR);
    float v_av = 0.5f*(vL+vR);
    float p_hat = 0.5f*r_av/b_av;                    /* rho_avg/(2 beta_avg) */
    float vel2  = 0.5f*(uL*uL+vL*vL) + 0.5f*(uR*uR+vR*vR);
    float e_int = 0.5f*(1.f/((GAMMA_V-1.f)*b_ln) - vel2);
    if (dir==0) {
        float f1 = r_ln*u_av;
        float f2 = f1*u_av + p_hat;
        float f3 = f1*v_av;
        F[0]=f1; F[1]=f2; F[2]=f3;  F[3]=f1*e_int + f2*u_av + f3*v_av;
    } else {
        float g1 = r_ln*v_av;
        float g2 = g1*u_av;
        float g3 = g1*v_av + p_hat;
        F[0]=g1; F[1]=g2; F[2]=g3;  F[3]=g1*e_int + g2*u_av + g3*v_av;
    }
}

/* ── HLLC flux from primitives WL, WR along normal (nx,ny) ─────────────── */
__device__ void hllc_n(const float WL[4], const float WR[4],
                       float nx, float ny, float F[4])
{
    float rL=WL[0], uL=WL[1], vL=WL[2], pL=WL[3];
    float rR=WR[0], uR=WR[1], vR=WR[2], pR=WR[3];

    float unL = uL*nx + vL*ny,  unR = uR*nx + vR*ny;
    float utL =-uL*ny + vL*nx,  utR =-uR*ny + vR*nx;
    float cL = sound_speed(pL,rL), cR = sound_speed(pR,rR);
    float EL = pL/(GAMMA_V-1.f)+0.5f*rL*(uL*uL+vL*vL);
    float ER = pR/(GAMMA_V-1.f)+0.5f*rR*(uR*uR+vR*vR);

    float FnL[4] = {rL*unL, rL*unL*uL+pL*nx, rL*unL*vL+pL*ny, (EL+pL)*unL};
    float FnR[4] = {rR*unR, rR*unR*uR+pR*nx, rR*unR*vR+pR*ny, (ER+pR)*unR};

    float SL = fminf(unL-cL, unR-cR);
    float SR = fmaxf(unL+cL, unR+cR);

    if (SL >= 0.f) {
        for (int q=0;q<4;q++) F[q]=FnL[q];
    } else if (SR <= 0.f) {
        for (int q=0;q<4;q++) F[q]=FnR[q];
    } else {
        float dL = rL*(SL-unL), dR = rR*(SR-unR);
        float Ss = (pR-pL+dL*unL-dR*unR)/(dL-dR);

        float rK  = (Ss>=0.f)?rL:rR,  unK=(Ss>=0.f)?unL:unR;
        float utK = (Ss>=0.f)?utL:utR, pK=(Ss>=0.f)?pL:pR;
        float EK  = (Ss>=0.f)?EL:ER,   SK=(Ss>=0.f)?SL:SR;
        float fact = rK*(SK-unK)/(SK-Ss);

        float Us[4];
        Us[0] = fact;
        Us[1] = fact * Ss;
        Us[2] = fact * utK;
        Us[3] = fact*(EK/rK+(Ss-unK)*(Ss+pK/(rK*(SK-unK))));

        float FK[4]  = {rK*unK, rK*unK*unK+pK, rK*unK*utK, (EK+pK)*unK};
        float UcK[4] = {rK, rK*unK, rK*utK, EK};
        float Fl[4];
        for (int q=0;q<4;q++) Fl[q] = FK[q]+SK*(Us[q]-UcK[q]);

        F[0]=Fl[0];
        F[1]=Fl[1]*nx-Fl[2]*ny;
        F[2]=Fl[1]*ny+Fl[2]*nx;
        F[3]=Fl[3];
    }
}

/* ── Wall-reflecting ghost state from boundary-tangent primitives ───────── */
/* dir=0: x-boundary (flip u),  dir=1: y-boundary (flip v)                  */
__device__ __forceinline__ void wall_state(const float W[4], int dir, float Wg[4])
{
    Wg[0]=W[0]; Wg[1]=W[1]; Wg[2]=W[2]; Wg[3]=W[3];
    if (dir==0) Wg[1]=-W[1];
    else        Wg[2]=-W[2];
}

/* ── Sanitize primitives: floor density/pressure, clamp velocity ─────────── */
#define UMAX 1e4f
__device__ __forceinline__ void sanitize_prim(float W[4])
{
    W[0] = fmaxf(W[0], 1e-14f);
    W[3] = fmaxf(W[3], 1e-14f);
    if (W[0] < 1e-5f) { W[1] = 0.f; W[2] = 0.f; return; }
    W[1] = fmaxf(fminf(W[1], UMAX), -UMAX);
    W[2] = fmaxf(fminf(W[2], UMAX), -UMAX);
}

/* ═══════════════════════════════════════════════════════════════════════════
 * DGSEM RHS kernel  (strong form, collocated LGL nodes)
 * ═══════════════════════════════════════════════════════════════════════════ *
 * Q   : conservative state [NVAR * NEL * NN2]
 * RHS : conservative RHS [NVAR * NEL * NN2]
 * lam : max wave speed per element [NEL]
 * N   : elements per side,  h = 1/N (physical element width)
 */
__global__ void compute_dg_rhs(
    const float * __restrict__ Q,
    float       * __restrict__ RHS,
    float       * __restrict__ lam,
    int N, float h)
{
    const int NEL  = N * N;
    const int nn   = NN;
    const int nn2  = NN2;
    const float jac    = 2.f / h;                 /* dxi/dx = 2/h_el */
    const float len_p  = h / (2*ORDER+1);
    const float dx2inv = 1.f/(len_p*len_p);

    /* one thread = one node;  EPB elements share the block, working set in
     * shared memory (keeps registers low → high occupancy, coalesced loads). */
    __shared__ float sW [EPB][NVAR][NN2];         /* primitive               */
    __shared__ float sUc[EPB][NVAR][NN2];         /* conservative (sanitized) */
    __shared__ float sA [EPB][NVAR][NN2];         /* x flux, reused as AV g1  */
    __shared__ float sB [EPB][NVAR][NN2];         /* y flux, reused as AV g2  */
    __shared__ float sRed[EPB][2][NN2];           /* [.,0]=θ  [.,1]=λ reduce  */

    const int ell = threadIdx.x / nn2;            /* element-local index      */
    const int nd  = threadIdx.x % nn2;            /* node index               */
    const int in  = nd % nn, jn = nd / nn;

    for (int base = blockIdx.x*EPB; base < NEL; base += gridDim.x*EPB) {
        const int el     = base + ell;
        const int active = (el < NEL);
        const int ej = active ? el / N : 0;
        const int ei = active ? el % N : 0;

        /* ── Phase 1: load, derive primitive/conservative/fluxes → shared ── */
        if (active) {
            float U[4] = { Q[SIDX(0,el,jn,in,NEL)], Q[SIDX(1,el,jn,in,NEL)],
                           Q[SIDX(2,el,jn,in,NEL)], Q[SIDX(3,el,jn,in,NEL)] };
            float Wt[4]; cons_to_prim_sane(U, Wt);
            float Ut[4]; p2c(Wt[0],Wt[1],Wt[2],Wt[3], Ut);
            float F[4], G[4];
            euler_flux_n(Wt[0],Wt[1],Wt[2],Wt[3], 1.f,0.f, F);
            euler_flux_n(Wt[0],Wt[1],Wt[2],Wt[3], 0.f,1.f, G);
            for (int q=0;q<NVAR;q++){
                sW[ell][q][nd]=Wt[q]; sUc[ell][q][nd]=Ut[q];
                sA[ell][q][nd]=F[q];  sB[ell][q][nd]=G[q];
            }
        }
        __syncthreads();

        float R[NVAR] = {0.f,0.f,0.f,0.f};
        float lam_node = 0.f, theta_node = 0.f;

        if (active) {
            /* ── Phase 2: volume term (this node only) ────────────────── */
#if EC_FLUX
            /* Entropy-conservative flux differencing 2 Σ_k D_ik F#(U_i,U_k). */
            float Wi[4]={sW[ell][0][nd],sW[ell][1][nd],sW[ell][2][nd],sW[ell][3][nd]};
            float ax[4]={0.f,0.f,0.f,0.f}, ay[4]={0.f,0.f,0.f,0.f};
            for (int k=0;k<nn;k++){
                float Wx[4]={sW[ell][0][jn*nn+k],sW[ell][1][jn*nn+k],sW[ell][2][jn*nn+k],sW[ell][3][jn*nn+k]};
                float Wy[4]={sW[ell][0][k*nn+in],sW[ell][1][k*nn+in],sW[ell][2][k*nn+in],sW[ell][3][k*nn+in]};
                float Fs[4], Gs[4];
                ec_flux(Wi, Wx, 0, Fs);
                ec_flux(Wi, Wy, 1, Gs);
                float dix=c_D[in][k], djy=c_D[jn][k];
                for (int q=0;q<NVAR;q++){ ax[q]+=dix*Fs[q]; ay[q]+=djy*Gs[q]; }
            }
            for (int q=0;q<NVAR;q++) R[q] = -jac*(2.f*ax[q]+2.f*ay[q]);
#else
            for (int q=0;q<NVAR;q++){
                float dfx=0.f, dgy=0.f;
                for (int k=0;k<nn;k++){
                    dfx += c_D[in][k]*sA[ell][q][jn*nn+k];
                    dgy += c_D[jn][k]*sB[ell][q][k*nn+in];
                }
                R[q] = -jac*(dfx+dgy);
            }
#endif
            /* ── Phase 3: surface lift (boundary nodes only) ──────────── */
            if (in==0) {                                   /* x left face */
                float WR_in[4]={sW[ell][0][jn*nn+0],sW[ell][1][jn*nn+0],sW[ell][2][jn*nn+0],sW[ell][3][jn*nn+0]};
                float WL_nb[4];
                if (ei>0){ int nb=ej*N+(ei-1);
                    float U[4]={Q[SIDX(0,nb,jn,nn-1,NEL)],Q[SIDX(1,nb,jn,nn-1,NEL)],Q[SIDX(2,nb,jn,nn-1,NEL)],Q[SIDX(3,nb,jn,nn-1,NEL)]};
                    cons_to_prim_sane(U,WL_nb);
                } else wall_state(WR_in,0,WL_nb);
                float fL[4]; hllc_n(WL_nb,WR_in,1.f,0.f,fL);
                for(int q=0;q<NVAR;q++) R[q]+=jac*c_winv[0]*(fL[q]-sA[ell][q][jn*nn+0]);
            }
            if (in==nn-1) {                                /* x right face */
                float WL_in[4]={sW[ell][0][jn*nn+nn-1],sW[ell][1][jn*nn+nn-1],sW[ell][2][jn*nn+nn-1],sW[ell][3][jn*nn+nn-1]};
                float WR_nb[4];
                if (ei<N-1){ int nb=ej*N+(ei+1);
                    float U[4]={Q[SIDX(0,nb,jn,0,NEL)],Q[SIDX(1,nb,jn,0,NEL)],Q[SIDX(2,nb,jn,0,NEL)],Q[SIDX(3,nb,jn,0,NEL)]};
                    cons_to_prim_sane(U,WR_nb);
                } else wall_state(WL_in,0,WR_nb);
                float fR[4]; hllc_n(WL_in,WR_nb,1.f,0.f,fR);
                for(int q=0;q<NVAR;q++) R[q]-=jac*c_winv[nn-1]*(fR[q]-sA[ell][q][jn*nn+nn-1]);
            }
            if (jn==0) {                                   /* y bottom face */
                float WT_in[4]={sW[ell][0][in],sW[ell][1][in],sW[ell][2][in],sW[ell][3][in]};
                float WB_nb[4];
                if (ej>0){ int nb=(ej-1)*N+ei;
                    float U[4]={Q[SIDX(0,nb,nn-1,in,NEL)],Q[SIDX(1,nb,nn-1,in,NEL)],Q[SIDX(2,nb,nn-1,in,NEL)],Q[SIDX(3,nb,nn-1,in,NEL)]};
                    cons_to_prim_sane(U,WB_nb);
                } else wall_state(WT_in,1,WB_nb);
                float gB[4]; hllc_n(WB_nb,WT_in,0.f,1.f,gB);
                for(int q=0;q<NVAR;q++) R[q]+=jac*c_winv[0]*(gB[q]-sB[ell][q][in]);
            }
            if (jn==nn-1) {                                /* y top face */
                float WB_in[4]={sW[ell][0][(nn-1)*nn+in],sW[ell][1][(nn-1)*nn+in],sW[ell][2][(nn-1)*nn+in],sW[ell][3][(nn-1)*nn+in]};
                float WT_nb[4];
                if (ej<N-1){ int nb=(ej+1)*N+ei;
                    float U[4]={Q[SIDX(0,nb,0,in,NEL)],Q[SIDX(1,nb,0,in,NEL)],Q[SIDX(2,nb,0,in,NEL)],Q[SIDX(3,nb,0,in,NEL)]};
                    cons_to_prim_sane(U,WT_nb);
                } else wall_state(WB_in,1,WT_nb);
                float gT[4]; hllc_n(WB_in,WT_nb,0.f,1.f,gT);
                for(int q=0;q<NVAR;q++) R[q]-=jac*c_winv[nn-1]*(gT[q]-sB[ell][q][(nn-1)*nn+in]);
            }

            /* ── Phase 3.5: per-node wave speed + shock sensor ────────── */
            float rho=fmaxf(sW[ell][0][nd],1e-12f), c2=GAMMA_V*sW[ell][3][nd]/rho;
            float c=sqrtf(fmaxf(c2,1e-14f));
            lam_node = fabsf(sW[ell][1][nd])+fabsf(sW[ell][2][nd])+c;
#if AV_ON
            float du=0.f, dv=0.f;
            for (int k=0;k<nn;k++){ du+=c_D[in][k]*sW[ell][U_ID][jn*nn+k]; dv+=c_D[jn][k]*sW[ell][V_ID][k*nn+in]; }
            float divu=fminf(jac*(du+dv),0.f), du2=divu*divu;
            theta_node = du2/(du2 + (AV_KSENSOR)*c2*dx2inv + 1e-30f);
#endif
        }

        /* ── Per-element reduction: λ_e = max λ_node,  θ_e = max θ_node ── */
        sRed[ell][0][nd]=theta_node; sRed[ell][1][nd]=lam_node;
        __syncthreads();
        float lam_e=0.f, theta_e=0.f;
        for (int k=0;k<nn2;k++){ theta_e=fmaxf(theta_e,sRed[ell][0][k]); lam_e=fmaxf(lam_e,sRed[ell][1][k]); }

#if AV_ON
        /* ── Shock-sensor artificial viscosity  ∂tW += ∇·(ν ∇U) ───────── */
        float nu_e = (AV_CAV) * theta_e * len_p * lam_e;
        __syncthreads();                          /* fluxes done → reuse sA/sB as g */
        if (active) {                             /* Pass 1: g = ν ∇Uc */
            float gx[4]={0,0,0,0}, gy[4]={0,0,0,0};
            for (int k=0;k<nn;k++) for (int q=0;q<NVAR;q++){
                gx[q]+=c_D[in][k]*sUc[ell][q][jn*nn+k];
                gy[q]+=c_D[jn][k]*sUc[ell][q][k*nn+in];
            }
            for (int q=0;q<NVAR;q++){ sA[ell][q][nd]=nu_e*jac*gx[q]; sB[ell][q][nd]=nu_e*jac*gy[q]; }
        }
        __syncthreads();
        if (active) {                             /* Pass 2: weak divergence */
            for (int q=0;q<NVAR;q++){
                float sx=0.f, sy=0.f;
                for (int k=0;k<nn;k++){ sx+=c_w[k]*c_D[k][in]*sA[ell][q][jn*nn+k]; sy+=c_w[k]*c_D[k][jn]*sB[ell][q][k*nn+in]; }
                R[q]-=jac*(c_winv[in]*sx + c_winv[jn]*sy);
            }
        }
#endif /* AV_ON */

        /* ── Write RHS + element λ_max ─────────────────────────────────── */
        if (active) {
            for (int q=0;q<NVAR;q++) RHS[SIDX(q,el,jn,in,NEL)] = R[q];
            if (nd==0) lam[el] = lam_e;
        }
        __syncthreads();                          /* shared reused next iter */
    }
}

/* ── SSP-RK3 kernels on conservative state ─────────────────────────────── */
__global__ void rk3_s1(float*U1, const float*U0, const float*L, float dt, int NEL)
{
    int ntot = NEL * NN2;
    for (int c = blockIdx.x*blockDim.x+threadIdx.x; c < ntot; c += blockDim.x*gridDim.x) {
        int el = c / NN2, n = c % NN2;
        float Uc[4];
        for (int q=0;q<NVAR;q++)
            Uc[q] = U0[SIDX(q,el,n/NN,n%NN,NEL)] + dt * L[SIDX(q,el,n/NN,n%NN,NEL)];
        sanitize_cons(Uc);
        for (int q=0;q<NVAR;q++) U1[SIDX(q,el,n/NN,n%NN,NEL)] = Uc[q];
    }
}
__global__ void rk3_s2(float*U2, const float*U0, const float*U1, const float*L,
                       float dt, int NEL)
{
    int ntot = NEL * NN2;
    for (int c = blockIdx.x*blockDim.x+threadIdx.x; c < ntot; c += blockDim.x*gridDim.x) {
        int el = c / NN2, n = c % NN2;
        float Uc[4];
        for (int q=0;q<NVAR;q++) {
            float u0 = U0[SIDX(q,el,n/NN,n%NN,NEL)];
            float u1 = U1[SIDX(q,el,n/NN,n%NN,NEL)];
            Uc[q] = 0.75f*u0 + 0.25f*(u1 + dt*L[SIDX(q,el,n/NN,n%NN,NEL)]);
        }
        sanitize_cons(Uc);
        for (int q=0;q<NVAR;q++) U2[SIDX(q,el,n/NN,n%NN,NEL)] = Uc[q];
    }
}
__global__ void rk3_s3(float*U, const float*U0, const float*U2, const float*L,
                       float dt, int NEL)
{
    int ntot = NEL * NN2;
    for (int c = blockIdx.x*blockDim.x+threadIdx.x; c < ntot; c += blockDim.x*gridDim.x) {
        int el = c / NN2, n = c % NN2;
        float Uc[4];
        for (int q=0;q<NVAR;q++) {
            float u0 = U0[SIDX(q,el,n/NN,n%NN,NEL)];
            float u2 = U2[SIDX(q,el,n/NN,n%NN,NEL)];
            Uc[q] = (1.f/3.f)*u0 + (2.f/3.f)*(u2 + dt*L[SIDX(q,el,n/NN,n%NN,NEL)]);
        }
        sanitize_cons(Uc);
        for (int q=0;q<NVAR;q++) U[SIDX(q,el,n/NN,n%NN,NEL)] = Uc[q];
    }
}

/* ═══════════════════════════════════════════════════════════════════════════
 * Zhang–Shu positivity-preserving limiter.
 * Rescales the nodal polynomial toward the GLL cell mean by a factor θ∈[0,1]
 * so that all nodes satisfy rho ≥ eps_rho and p ≥ eps_p while the cell mean
 * (hence conservation) is unchanged.
 * ═══════════════════════════════════════════════════════════════════════════ */
__global__ void positivity_limiter(float *Q, int N, float eps_rho, float eps_p)
{
    const int NEL = N * N;
    const int nn  = NN;

    for (int el = blockIdx.x*blockDim.x+threadIdx.x; el < NEL;
             el += blockDim.x*gridDim.x) {

        /* GLL cell mean of the conservative state: Ubar = (1/4) Σ w_i w_j U_ij */
        float Ubar[NVAR] = {0.f,0.f,0.f,0.f};
        for (int jn=0; jn<nn; jn++) for (int in=0; in<nn; in++) {
            float wij = 0.25f * c_w[in] * c_w[jn];
            for (int q=0;q<NVAR;q++) Ubar[q] += wij * Q[SIDX(q,el,jn,in,NEL)];
        }
        /* Keep the mean admissible (should already hold; guard anyway). */
        Ubar[0] = fmaxf(Ubar[0], eps_rho);
        {
            float p = pressure_from_cons(Ubar);
            if (p < eps_p)
                Ubar[3] = eps_p/(GAMMA_V-1.f)+0.5f*(Ubar[1]*Ubar[1]+Ubar[2]*Ubar[2])/Ubar[0];
        }

        /* θ1: density.  θ = min(1, (rhobar-eps)/(rhobar-rho_min)). */
        float rho_min = 1e30f;
        for (int n=0; n<NN2; n++) rho_min = fminf(rho_min, Q[SIDX(0,el,n/nn,n%nn,NEL)]);
        float t_rho = 1.f;
        if (rho_min < eps_rho)
            t_rho = (Ubar[0]-eps_rho) / fmaxf(Ubar[0]-rho_min, 1e-30f);
        t_rho = fmaxf(0.f, fminf(1.f, t_rho));

        /* θ2: pressure.  Scale each node's deviation until p(node) ≥ eps_p.  */
        float theta = t_rho;
        for (int n=0; n<NN2; n++) {
            float U[4] = { Q[SIDX(0,el,n/nn,n%nn,NEL)], Q[SIDX(1,el,n/nn,n%nn,NEL)],
                           Q[SIDX(2,el,n/nn,n%nn,NEL)], Q[SIDX(3,el,n/nn,n%nn,NEL)] };
            /* Apply the density scaling first, then test pressure. */
            float Ud[4];
            for (int q=0;q<NVAR;q++) Ud[q] = Ubar[q] + t_rho*(U[q]-Ubar[q]);
            float p = pressure_from_cons(Ud);
            if (p < eps_p) {
                /* Bisection for the largest t∈[0,θ] giving p≥eps_p. */
                float lo=0.f, hi=t_rho;
                for (int it=0; it<20; it++) {
                    float tm=0.5f*(lo+hi);
                    float Um[4];
                    for (int q=0;q<NVAR;q++) Um[q]=Ubar[q]+tm*(U[q]-Ubar[q]);
                    if (pressure_from_cons(Um) >= eps_p) lo=tm; else hi=tm;
                }
                theta = fminf(theta, lo);
            }
        }

        /* Apply the single element-wide factor θ (keeps the mean exact). */
        for (int jn=0; jn<nn; jn++) for (int in=0; in<nn; in++)
            for (int q=0;q<NVAR;q++) {
                float u = Q[SIDX(q,el,jn,in,NEL)];
                Q[SIDX(q,el,jn,in,NEL)] = Ubar[q] + theta*(u-Ubar[q]);
            }
    }
}

/* ── max reduction ──────────────────────────────────────────────────────── */
__global__ void reduce_max(const float*in, float*out, int n){
    extern __shared__ float sm[];
    int tid=threadIdx.x; float v=0.f;
    for(int k=blockIdx.x*blockDim.x+tid; k<n; k+=blockDim.x*gridDim.x) v=fmaxf(v,in[k]);
    sm[tid]=v; __syncthreads();
    for(int s=BLOCK1D/2;s>0;s>>=1){if(tid<s)sm[tid]=fmaxf(sm[tid],sm[tid+s]);__syncthreads();}
    if(tid==0) out[blockIdx.x]=sm[0];
}
static float gpu_max(const float*d_in, float*d_tmp, int n){
    reduce_max<<<GS_NBLK,BLOCK1D,BLOCK1D*sizeof(float)>>>(d_in,d_tmp,n);
    reduce_max<<<1,BLOCK1D,BLOCK1D*sizeof(float)>>>(d_tmp,d_tmp,GS_NBLK);
    float v; CK(cudaMemcpy(&v,d_tmp,sizeof(float),cudaMemcpyDeviceToHost));
    return v;
}

/* ── Diagnostic: per-element shock sensor θ_e (max over nodes) ──────────── */
__global__ void compute_sensor(const float * __restrict__ Q,
                               float * __restrict__ sensor, int N, float h)
{
    const int NEL = N * N;
    const int nn  = NN;
    const float jac = 2.f / h;
    const float len_p = h / (2*ORDER+1);
    const float dx2inv = 1.f/(len_p*len_p);
    for (int el = blockIdx.x*blockDim.x+threadIdx.x; el < NEL;
             el += blockDim.x*gridDim.x) {
        float u[NN2], v[NN2], c2[NN2];
        for (int n=0; n<nn*nn; n++) {
            float U[4]={Q[SIDX(0,el,n/nn,n%nn,NEL)],Q[SIDX(1,el,n/nn,n%nn,NEL)],
                        Q[SIDX(2,el,n/nn,n%nn,NEL)],Q[SIDX(3,el,n/nn,n%nn,NEL)]};
            float Wt[4]; cons_to_prim_sane(U,Wt);
            u[n]=Wt[1]; v[n]=Wt[2]; c2[n]=GAMMA_V*Wt[3]/fmaxf(Wt[0],1e-12f);
        }
        float theta_e=0.f;
        for (int jn=0; jn<nn; jn++) for (int in=0; in<nn; in++) {
            float du=0.f, dv=0.f;
            for (int k=0;k<nn;k++){ du+=c_D[in][k]*u[jn*nn+k]; dv+=c_D[jn][k]*v[k*nn+in]; }
            float divu=jac*(du+dv), du2=divu*divu;
            float th = du2/(du2 + (AV_KSENSOR)*c2[jn*nn+in]*dx2inv + 1e-30f);
            theta_e=fmaxf(theta_e,th);
        }
        sensor[el]=theta_e;
    }
}

/* ── IC kernel: circular Sod on LGL nodes, stored in conservative form ─── */
__global__ void ic_kernel(float*U, int N, float h)
{
    const int NEL = N * N;
    int ntot = NEL * NN2;
    for (int c = blockIdx.x*blockDim.x+threadIdx.x; c < ntot; c += blockDim.x*gridDim.x) {
        int el = c / NN2, n = c % NN2;
        int ej = el / N, ei = el % N;
        int jn = n / NN, in_node = n % NN;

        /* Physical coordinates of this LGL node */
        float xl = ei * h, yb = ej * h;
        float xnode = xl + (c_xi[in_node] + 1.f) * 0.5f * h;
        float ynode = yb + (c_xi[jn]      + 1.f) * 0.5f * h;

        float r     = sqrtf((xnode-0.5f)*(xnode-0.5f)+(ynode-0.5f)*(ynode-0.5f));
        float delta = 0.5f * h;
        float phi   = 0.5f*(1.f+tanhf((0.25f-r)/delta));

        float C[4];
        p2c(0.125f + (RHO_HI-0.125f)*phi, 0.f, 0.f, 0.1f + (P_HI-0.1f)*phi, C);
        U[SIDX(RHO_ID, el, jn, in_node, NEL)] = C[0];
        U[SIDX(U_ID,   el, jn, in_node, NEL)] = C[1];
        U[SIDX(V_ID,   el, jn, in_node, NEL)] = C[2];
        U[SIDX(P_ID,   el, jn, in_node, NEL)] = C[3];
    }
}

/* ═══════════════════════════════════════════════════════════════════════════
 * PPM output helpers
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

/* Downsample conservative state to N×N cell-average rho/p images.          */
static void downsample_rho_p(const float*Q, float*rho_img, float*p_img, int N)
{
    const int NEL = N * N;
    for (int el=0; el<NEL; el++) {
        float rho_sum=0.f, p_sum=0.f;
        for (int n=0; n<NN2; n++) {
            float rho = Q[SIDX(RHO_ID, el, n/NN, n%NN, NEL)];
            float mx  = Q[SIDX(U_ID,   el, n/NN, n%NN, NEL)];
            float my  = Q[SIDX(V_ID,   el, n/NN, n%NN, NEL)];
            float E   = Q[SIDX(P_ID,   el, n/NN, n%NN, NEL)];
            if (rho < 1e-14f) rho = 1e-14f;
            float p = (GAMMA_V - 1.f) * (E - 0.5f*(mx*mx + my*my)/rho);
            rho_sum += rho;
            p_sum   += (p > 0.f) ? p : 0.f;
        }
        rho_img[el] = rho_sum / NN2;
        p_img  [el] = p_sum   / NN2;
    }
}

static void image_range(const float *img, int n, float *vmin, float *vmax)
{
    float lo = img[0], hi = img[0];
    for (int i = 1; i < n; i++) {
        if (img[i] < lo) lo = img[i];
        if (img[i] > hi) hi = img[i];
    }
    float pad = 0.05f * (hi - lo + 1e-12f);
    *vmin = lo - pad;
    *vmax = hi + pad;
    if (!(*vmax > *vmin)) { *vmin = lo - 1e-6f; *vmax = hi + 1e-6f; }
}

static void write_ppm_2x1(const char *fname,
                           const float *tl, float tl_min, float tl_max,
                           const float *tr, float tr_min, float tr_max,
                           int N)
{
    FILE *f = fopen(fname, "wb");
    if (!f) { fprintf(stderr, "Cannot open %s\n", fname); return; }
    fprintf(f, "P6\n%d %d\n255\n", 2*N, N);
    unsigned char pix[3];
    for (int r=0; r<N; r++) {
        int phys_j = N-1-r;
        for (int c=0; c<2*N; c++) {
            const float *panel = (c<N) ? tl : tr;
            float vmin = (c<N) ? tl_min : tr_min;
            float vmax = (c<N) ? tl_max : tr_max;
            int sub_i = (c<N) ? c : c-N;
            float t=(panel[phys_j*N+sub_i]-vmin)/(vmax-vmin+1e-30f);
            if(!(t>0.f))t=0.f; if(t>1.f)t=1.f;
            int idx=(int)(t*255.f);
            if(idx<0)idx=0; if(idx>255)idx=255;
            if(c<N){ pix[0]=PLASMA[idx][0];pix[1]=PLASMA[idx][1];pix[2]=PLASMA[idx][2]; }
            else   { pix[0]=VIRIDIS[idx][0];pix[1]=VIRIDIS[idx][1];pix[2]=VIRIDIS[idx][2]; }
            fwrite(pix,1,3,f);
        }
    }
    fclose(f);
    printf("  Saved %s\n", fname);
}

/* Single-panel colormap PPM (used for the shock-sensor field, range [0,1]). */
static void write_ppm_gray(const char *fname, const float *img,
                           float vmin, float vmax, int N)
{
    FILE *f = fopen(fname, "wb");
    if (!f) { fprintf(stderr, "Cannot open %s\n", fname); return; }
    fprintf(f, "P6\n%d %d\n255\n", N, N);
    unsigned char pix[3];
    for (int r=0; r<N; r++) {
        int phys_j = N-1-r;
        for (int cc=0; cc<N; cc++) {
            float t=(img[phys_j*N+cc]-vmin)/(vmax-vmin+1e-30f);
            if(!(t>0.f))t=0.f; if(t>1.f)t=1.f;
            int idx=(int)(t*255.f); if(idx<0)idx=0; if(idx>255)idx=255;
            pix[0]=PLASMA[idx][0]; pix[1]=PLASMA[idx][1]; pix[2]=PLASMA[idx][2];
            fwrite(pix,1,3,f);
        }
    }
    fclose(f);
}

/* ═══════════════════════════════════════════════════════════════════════════
 * main
 * ═══════════════════════════════════════════════════════════════════════════ */
int main(int argc, char **argv)
{
    int   N   = (argc > 1) ? atoi(argv[1]) : 50;
    /* dt = CFL * h / (lam * NN).  The old 0.1/ORDER was only ~4% of the
     * stability limit; CFL_COEF=0.4 (default) has ~20% margin (N=100 blows up
     * near 0.6) and is ~10x faster.  Override with -DCFL_COEF=<f>.            */
#ifndef CFL_COEF
#  define CFL_COEF 0.6f
#endif
    float CFL = CFL_COEF, t_end = 1.0f, h = 1.0f / N;

    /* ── Build DG operators on host ──────────────────────────────────── */
    double xs[4];
    for (int i=0;i<NN;i++) xs[i] = lgl_xi[ORDER-1][i];
    double D_h[4][4] = {};
    build_D(NN, xs, D_h);

    float D_f[4][4]={}, w_f[4]={}, winv_f[4]={}, xi_f[4]={};
    for (int i=0;i<NN;i++) {
        w_f[i]    = (float)gll_w[ORDER-1][i];
        winv_f[i] = (float)(1.0/gll_w[ORDER-1][i]);
        xi_f[i]   = (float)lgl_xi[ORDER-1][i];
        for (int k=0;k<NN;k++) D_f[i][k]=(float)D_h[i][k];
    }
    CK(cudaMemcpyToSymbol(c_D,    D_f,    sizeof(D_f)));
    CK(cudaMemcpyToSymbol(c_w,    w_f,    sizeof(w_f)));
    CK(cudaMemcpyToSymbol(c_winv, winv_f, sizeof(winv_f)));
    CK(cudaMemcpyToSymbol(c_xi,   xi_f,   sizeof(xi_f)));

    /* ── GPU info ────────────────────────────────────────────────────── */
    int dev; cudaDeviceProp prop;
    CK(cudaGetDevice(&dev));
    CK(cudaGetDeviceProperties(&prop,dev));
    printf("========================================================\n");
    printf("  Device  : %s\n", prop.name);
    printf("  DGSEM (nodal DG, LGL/Lobatto nodes) p=%d  --  Circular Sod\n", ORDER);
    printf("  N=%d  h=%.5f  NN=%d sol nodes/dir  CFL=%.3f\n", N, h, NN, CFL);
    printf("  volume=%s  faces=HLLC  shock-AV=%s(C=%.2f,K=%.3f)\n",
           EC_FLUX ? "EC-flux-diff" : "collocation",
           AV_ON ? "on" : "off", (float)AV_CAV, (float)AV_KSENSOR);
    printf("  Total DOF: %d\n", N*N*NN2*NVAR);
    printf("========================================================\n");

    /* ── Allocations ─────────────────────────────────────────────────── */
    int    NEL  = N * N;
    size_t szQ  = (size_t)NVAR * NEL * NN2 * sizeof(float);
    size_t sz1  = (size_t)NEL * sizeof(float);

    /* RHS kernel launch: one thread per node, EPB elements per block */
    const int RHS_THR = EPB * NN2;
    const int RHS_BLK = (NEL + EPB - 1) / EPB;

    float *d_Q, *d_Q0, *d_Q1, *d_Q2, *d_RHS, *d_lam, *d_tmp;
    CK(cudaMalloc(&d_Q,   szQ));
    CK(cudaMalloc(&d_Q0,  szQ));
    CK(cudaMalloc(&d_Q1,  szQ));
    CK(cudaMalloc(&d_Q2,  szQ));
    CK(cudaMalloc(&d_RHS, szQ));
    CK(cudaMalloc(&d_lam, sz1));
    CK(cudaMalloc(&d_tmp, GS_NBLK*sizeof(float)));

    float *d_sensor; CK(cudaMalloc(&d_sensor, sz1));
    float *h_Q   = (float*)malloc(szQ);
    float *h_rho = (float*)malloc(NEL*sizeof(float));
    float *h_p   = (float*)malloc(NEL*sizeof(float));
    float *h_sen = (float*)malloc(NEL*sizeof(float));

    /* colormap range: computed once from frame 0, reused for all frames */
    float g_rho_min=0.f, g_rho_max=1.f, g_p_min=0.f, g_p_max=1.f;
    int   g_range_set = 0;

#define WRITE_FRAME(fi) do { \
    CK(cudaDeviceSynchronize()); \
    CK(cudaMemcpy(h_Q, d_Q, szQ, cudaMemcpyDeviceToHost)); \
    downsample_rho_p(h_Q, h_rho, h_p, N); \
    { float _lo=h_rho[0],_hi=h_rho[0],_plo=h_p[0],_phi=h_p[0]; \
      for(int _i=1;_i<NEL;_i++){_lo=fminf(_lo,h_rho[_i]);_hi=fmaxf(_hi,h_rho[_i]); \
                                  _plo=fminf(_plo,h_p[_i]);_phi=fmaxf(_phi,h_p[_i]);} \
      printf("  [frame %d] rho=[%.4f,%.4f]  p=[%.4f,%.4f]\n",(fi),_lo,_hi,_plo,_phi); } \
    if (!g_range_set) { \
        image_range(h_rho, NEL, &g_rho_min, &g_rho_max); \
        image_range(h_p,   NEL, &g_p_min,   &g_p_max); \
        g_range_set = 1; \
    } \
    char _fn[64]; sprintf(_fn,"dgsem_sod_%04d.ppm",(fi)); \
    write_ppm_2x1(_fn, h_rho, g_rho_min, g_rho_max, h_p, g_p_min, g_p_max, N); \
    compute_sensor<<<GS_NBLK,BLOCK1D>>>(d_Q, d_sensor, N, h); \
    CK(cudaMemcpy(h_sen, d_sensor, NEL*sizeof(float), cudaMemcpyDeviceToHost)); \
    char _sf[64]; sprintf(_sf,"dgsem_sensor_%04d.ppm",(fi)); \
    write_ppm_gray(_sf, h_sen, 0.f, 1.f, N); \
} while(0)

    /* ── IC ─────────────────────────────────────────────────────────── */
    ic_kernel<<<GS_NBLK,BLOCK1D>>>(d_Q, N, h);
    positivity_limiter<<<GS_NBLK,BLOCK1D>>>(d_Q, N, 1e-6f, 1e-6f);
    CK(cudaDeviceSynchronize());
    WRITE_FRAME(0);

    /* ── Time loop ───────────────────────────────────────────────────── */
    const int N_FRAMES = 10;
    int   frame  = 1;
    float t_next = t_end / N_FRAMES;
    int   step   = 0;
    float t      = 0.f;
    struct timespec ts0, ts1;
    clock_gettime(CLOCK_MONOTONIC, &ts0);

    float lam_max = 0.f;

    while (t < t_end) {
        /* ── Compute RHS at current state, get lam for dt ────────────── */
        compute_dg_rhs<<<RHS_BLK,RHS_THR>>>(d_Q, d_RHS, d_lam, N, h);
        CK(cudaDeviceSynchronize());
        lam_max = gpu_max(d_lam, d_tmp, NEL);

        if (!(lam_max>0.f)) { fprintf(stderr,"NaN lam at step %d\n",step); break; }
        float dt = CFL * h / (lam_max * NN);
        float t_target = (t_next < t_end) ? t_next : t_end;
        if (t+dt > t_target) dt = t_target - t;
        if (dt < 1e-14f) {
            fprintf(stderr,"dt underflow at step %d t=%.6f lam=%.3e\n",step,t,lam_max);
            WRITE_FRAME(99); break;
        }

        CK(cudaMemcpy(d_Q0, d_Q, szQ, cudaMemcpyDeviceToDevice));

        /* RK3 stage 1: Q1 = Q0 + dt*L(Q0) */
        rk3_s1<<<GS_NBLK,BLOCK1D>>>(d_Q1, d_Q0, d_RHS, dt, NEL);
        positivity_limiter<<<GS_NBLK,BLOCK1D>>>(d_Q1, N, 1e-6f, 1e-6f);

        /* RK3 stage 2: Q2 = 3/4*Q0 + 1/4*(Q1 + dt*L(Q1)) */
        compute_dg_rhs<<<RHS_BLK,RHS_THR>>>(d_Q1, d_RHS, d_lam, N, h);
        rk3_s2<<<GS_NBLK,BLOCK1D>>>(d_Q2, d_Q0, d_Q1, d_RHS, dt, NEL);
        positivity_limiter<<<GS_NBLK,BLOCK1D>>>(d_Q2, N, 1e-6f, 1e-6f);

        /* RK3 stage 3: Q = 1/3*Q0 + 2/3*(Q2 + dt*L(Q2)) */
        compute_dg_rhs<<<RHS_BLK,RHS_THR>>>(d_Q2, d_RHS, d_lam, N, h);
        rk3_s3<<<GS_NBLK,BLOCK1D>>>(d_Q, d_Q0, d_Q2, d_RHS, dt, NEL);
        positivity_limiter<<<GS_NBLK,BLOCK1D>>>(d_Q, N, 1e-6f, 1e-6f);
        lam_max = gpu_max(d_lam, d_tmp, NEL);

        if (!(lam_max>0.f) || lam_max > 1e6f) {
            fprintf(stderr,"blow-up step %d t=%.6f lam=%.3e\n",step+1,t+dt,lam_max);
            WRITE_FRAME(99); break;
        }
        t += dt; step++;

        if (t >= t_next-1e-12f && frame <= N_FRAMES) {
            clock_gettime(CLOCK_MONOTONIC,&ts1);
            float el=(ts1.tv_sec-ts0.tv_sec)+(ts1.tv_nsec-ts0.tv_nsec)*1e-9f;
            printf("  frame %2d/%d  step %5d  t=%.5f  elapsed=%.1fs\n",
                   frame,N_FRAMES,step,t,el);
            fflush(stdout);
            WRITE_FRAME(frame);
            frame++; t_next=frame*t_end/N_FRAMES;
        } else if (step%50==0) {
            clock_gettime(CLOCK_MONOTONIC,&ts1);
            float el=(ts1.tv_sec-ts0.tv_sec)+(ts1.tv_nsec-ts0.tv_nsec)*1e-9f;
            float eta=(t>0)?(t_end-t)/t*el:0.f;
            printf("  step %5d  t=%.5f  dt=%.3e  elapsed=%.1fs  ETA=%.0fs\n",
                   step,t,dt,el,eta);
            fflush(stdout);
        }
    }
    CK(cudaDeviceSynchronize());
    clock_gettime(CLOCK_MONOTONIC,&ts1);
    float wall=(ts1.tv_sec-ts0.tv_sec)+(ts1.tv_nsec-ts0.tv_nsec)*1e-9f;
    printf("  Done: %d steps  t=%.6f  wall=%.2fs\n",step,t,wall);
    printf("  Frames: dgsem_sod_0000.ppm .. dgsem_sod_%04d.ppm\n",N_FRAMES);

    free(h_Q); free(h_rho); free(h_p);
    cudaFree(d_Q); cudaFree(d_Q0); cudaFree(d_Q1); cudaFree(d_Q2);
    cudaFree(d_RHS); cudaFree(d_lam); cudaFree(d_tmp);
    return 0;
}
