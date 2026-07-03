/*
 * fr_lobatto_2d.cu
 * 2D Circular Sod — Spectral Difference on LGL solution nodes
 * HLLC interface flux · no artificial viscosity · SSP-RK3
 *
 * Compile (polynomial degree p=1..3):
 *   nvcc -O3 -arch=native -DORDER=3 -o fr_p3 fr_lobatto_2d.cu -lm
 *
 * Run:
 *   ./fr_p3          # N=50 elements/side
 *   ./fr_p3 80       # N=80 elements/side
 *
 * State layout (conservative):  Q[var * NEL * NN2 + el * NN2 + jn * NN + in]
 *   var ∈ {rho=0, rhou=1, rhov=2, E=3},  NEL = N*N,  NN = ORDER+1,  NN2 = NN*NN
 *   jn = node index in y,  in = node index in x,  both in [0, NN)
 *
 * Spectral Difference flux divergence (x-direction at node (jn,in)):
 *   NF = NN+1 flux nodes per direction:  xf = { -1, GL_interior..., +1 }
 *   1. Interpolate Q from NN solution nodes to NF flux nodes via c_Isf
 *   2. Evaluate F at interior flux nodes (k=1..NF-2)
 *   3. F at k=0      = HLLC(neighbor_right_face, this_left_face)
 *      F at k=NF-1   = HLLC(this_right_face, neighbor_left_face)
 *   4. dF/dx|_in = (2/h) * Σ_k c_Dfs[in][k] * F_k
 *      where c_Dfs[i][k] = d(ℓ^f_k)/dξ at ξ = xs_i  (flux Lagrange basis)
 *   y-direction: symmetric
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
#  define ORDER 3          /* polynomial degree p: 1..3 */
#endif
#if ORDER > 3
#  error "ORDER > 3 is not supported: remove p=4..7 nodes/weights or extend tables."
#endif
#define NN      (ORDER + 1)   /* solution nodes per direction per element */
#define NF      (NN + 1)      /* flux nodes per direction per element     */
#define NN2     (NN * NN)     /* solution nodes per element               */

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

/* ── SD matrices in constant memory ────────────────────────────────────── */
/* Solution nodes: NN LGL nodes (include ±1).                               */
/* Flux nodes:     NF = NN+1 nodes:  {-1, GL_interior..., +1}               */
__constant__ float c_xf [5];      /* flux node positions in [-1,1]          */
__constant__ float c_Isf[5][4];   /* interp sol→flux:  c_Isf[k][j]=ℓ_j(xf_k) */
__constant__ float c_Dfs[4][5];   /* d(flux Lagrange)/dxi at sol nodes     */
__constant__ float c_xi[8];       /* LGL solution node coordinates          */
__constant__ float c_w [8];       /* GLL quadrature weights  (sum = 2)      */

/* ── LGL nodes on [-1,1] for p=1..3 (NN=2..4) ─────────────────────────── */
static const double lgl_xi[3][4] = {
    /* p=1 */ {-1.0, 1.0},
    /* p=2 */ {-1.0, 0.0, 1.0},
    /* p=3 */ {-1.0,-0.4472135954999579, 0.4472135954999579, 1.0},
};

/* ── Gauss-Legendre interior flux nodes for p=1..3 (NF-2 = NN-1 nodes) ── */
/* Combined with ±1 endpoints these give the NF = NN+1 flux nodes.          */
static const double gl_interior[3][3] = {
    /* p=1, 1 GL node  */ {0.0},
    /* p=2, 2 GL nodes */ {-0.5773502691896258, 0.5773502691896258},
    /* p=3, 3 GL nodes */ {-0.7745966692414834, 0.0, 0.7745966692414834},
};

/* ── GLL quadrature weights on [-1,1] for p=1..3 (sum = 2) ─────────────── */
static const double gll_w[3][4] = {
    /* p=1, n=2 */ {1.0, 1.0},
    /* p=2, n=3 */ {0.33333333333333333, 1.33333333333333333, 0.33333333333333333},
    /* p=3, n=4 */ {0.16666666666666667, 0.83333333333333333,
                    0.83333333333333333, 0.16666666666666667},
};

/* Evaluate Lagrange basis ℓ_j (defined on nodes xs[0..n-1]) at point x.  */
static double lagrange_val(int n, const double *xs, int j, double x)
{
    double v = 1.0;
    for (int k = 0; k < n; k++)
        if (k != j) v *= (x - xs[k]) / (xs[j] - xs[k]);
    return v;
}

/* Build interpolation matrix Isf[nf][ns]:  Isf[k][j] = ℓ_j(xf[k])        */
/* xs[0..ns-1] = solution nodes,  xf[0..nf-1] = flux nodes                 */
static void build_Isf(int ns, const double *xs, int nf, const double *xf,
                      double Isf[5][4])
{
    for (int k = 0; k < nf; k++)
        for (int j = 0; j < ns; j++)
            Isf[k][j] = lagrange_val(ns, xs, j, xf[k]);
}

/* Build flux-poly derivative matrix Dfs[ns][nf]:                           */
/* Dfs[i][k] = dℓ^f_k/dxi  at  xi = xs[i]                                 */
/* ℓ^f_k is the k-th Lagrange basis on the flux nodes xf[0..nf-1].         */
static void build_Dfs(int ns, const double *xs, int nf, const double *xf,
                      double Dfs[4][5])
{
    /* barycentric weights for flux nodes */
    double wf[5];
    for (int k = 0; k < nf; k++) {
        wf[k] = 1.0;
        for (int m = 0; m < nf; m++)
            if (m != k) wf[k] *= (xf[k] - xf[m]);
        wf[k] = 1.0 / wf[k];
    }
    for (int i = 0; i < ns; i++) {
        for (int k = 0; k < nf; k++) {
            /* check if xs[i] coincides with xf[k] */
            int hit = 0;
            for (int m = 0; m < nf; m++) {
                if (m == k) continue;
                if (fabs(xs[i] - xf[m]) < 1e-14) { hit = 1; break; }
            }
            if (fabs(xs[i] - xf[k]) < 1e-14) {
                /* xs[i] == xf[k]: use standard derivative formula */
                double d = 0.0;
                for (int m = 0; m < nf; m++)
                    if (m != k) d += (wf[m]/wf[k]) / (xf[k] - xf[m]);
                Dfs[i][k] = d;
            } else {
                /* xs[i] != xf[k]: d/dxi ℓ^f_k(xi) = (wf[k]/(sum_m wf[m]/(xi-xf[m]))) * ... */
                /* Direct formula: dL_k/dx at x = w_k/(x-x_k) * L(x) / sum_{m!=k} adapted  */
                /* Simpler: use product rule on Lagrange form */
                double num = 0.0;
                double den = 1.0;
                for (int m = 0; m < nf; m++)
                    if (m != k) den *= (xs[i] - xf[m]) / (xf[k] - xf[m]);
                /* den = L_k(xs[i]) without (xs[i]-xf[k]) factor... use direct sum */
                /* dL_k/dx|_{x=xs[i]} = L_k(xs[i]) * sum_{m!=k} 1/(xs[i]-xf[m]) */
                double Lk = lagrange_val(nf, xf, k, xs[i]);
                for (int m = 0; m < nf; m++)
                    if (m != k) num += 1.0 / (xs[i] - xf[m]);
                Dfs[i][k] = Lk * num;
                (void)den;
            }
        }
    }
}

/* ═══════════════════════════════════════════════════════════════════════════
 * Device helpers  (same as FV solver)
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

/* ── HLLE flux from primitives WL, WR along normal (nx,ny) ─────────────── */
__device__ void hlle_n(const float WL[4], const float WR[4],
                       float nx, float ny, float F[4])
{
    float rL=WL[0], uL=WL[1], vL=WL[2], pL=WL[3];
    float rR=WR[0], uR=WR[1], vR=WR[2], pR=WR[3];
    float unL=uL*nx+vL*ny, unR=uR*nx+vR*ny;
    float cL=sound_speed(pL,rL), cR=sound_speed(pR,rR);
    float EL=pL/(GAMMA_V-1.f)+0.5f*rL*(uL*uL+vL*vL);
    float ER=pR/(GAMMA_V-1.f)+0.5f*rR*(uR*uR+vR*vR);
    float SL=fminf(unL-cL,unR-cR), SR=fmaxf(unL+cL,unR+cR);
    float FL[4]={rL*unL,rL*unL*uL+pL*nx,rL*unL*vL+pL*ny,(EL+pL)*unL};
    float FR[4]={rR*unR,rR*unR*uR+pR*nx,rR*unR*vR+pR*ny,(ER+pR)*unR};
    float UL[4]={rL,rL*uL,rL*vL,EL};
    float UR[4]={rR,rR*uR,rR*vR,ER};
    if(SL>=0.f){ for(int q=0;q<4;q++) F[q]=FL[q]; }
    else if(SR<=0.f){ for(int q=0;q<4;q++) F[q]=FR[q]; }
    else{
        float dS=SR-SL+1e-30f;
        for(int q=0;q<4;q++) F[q]=(SR*FL[q]-SL*FR[q]+SL*SR*(UR[q]-UL[q]))/dS;
    }
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
 * SD RHS kernel
 * ═══════════════════════════════════════════════════════════════════════════ *
 * Q   : conservative state [NVAR * NEL * NN2]
 * RHS : conservative RHS [NVAR * NEL * NN2]
 * lam : max wave speed per element [NEL]
 * N   : elements per side,  h = 1/N (physical element width)
 */
__global__ void compute_sd_rhs(
    const float * __restrict__ Q,
    float       * __restrict__ RHS,
    float       * __restrict__ lam,
    int N, float h)
{
    const int NEL  = N * N;
    const int nn   = NN;             /* solution nodes per direction  */
    const int nf   = NF;             /* flux nodes per direction      */
    const int nn2  = NN2;
    const float jac = 2.f / h;      /* dxi/dx = 2/h_el               */

    for (int el = blockIdx.x*blockDim.x+threadIdx.x; el < NEL;
             el += blockDim.x*gridDim.x) {

        const int ej = el / N, ei = el % N;

        /* ── Load element conservative state and derive primitives ─── */
        float Qc[NVAR][NN2];          /* conservative [var][node]      */
        float Qe[NVAR][NN2];          /* primitive    [var][node]      */
        for (int q=0; q<NVAR; q++)
            for (int n=0; n<nn2; n++)
                Qc[q][n] = Q[SIDX(q, el, n/nn, n%nn, NEL)];
        for (int n=0; n<nn2; n++) {
            float Utmp[4] = { Qc[0][n], Qc[1][n], Qc[2][n], Qc[3][n] };
            float Wtmp[4];
            cons_to_prim_sane(Utmp, Wtmp);
            for (int q=0; q<NVAR; q++) Qe[q][n] = Wtmp[q];
        }

        /* ── Accumulate RHS ─────────────────────────────────────────── */
        float R[NVAR][NN2];
        for (int q=0; q<NVAR; q++)
            for (int n=0; n<nn2; n++) R[q][n] = 0.f;

        /* ── x-direction SD flux divergence ─────────────────────────── */
        for (int jn=0; jn<nn; jn++) {

            /* 1. Interpolate conservative row to NF flux nodes          */
            float Qf[NVAR][5];
            for (int k=0; k<nf; k++) {
                float Uk[4];
                for (int q=0; q<NVAR; q++) {
                    float v = 0.f;
                    for (int j=0; j<nn; j++) v += c_Isf[k][j]*Qc[q][jn*nn+j];
                    Uk[q] = v;
                }
                float Wk[4];
                cons_to_prim_sane(Uk, Wk);
                for (int q=0;q<NVAR;q++) Qf[q][k]=Wk[q];
            }

            /* 2. Euler flux at interior flux nodes k=1..nf-2            */
            float Ff[NVAR][5];
            for (int k=1; k<nf-1; k++) {
                float F4[4];
                euler_flux_n(Qf[0][k],Qf[1][k],Qf[2][k],Qf[3][k], 1.f,0.f, F4);
                for (int q=0; q<NVAR; q++) Ff[q][k] = F4[q];
            }

            /* 3. HLLC at left boundary (k=0) and right boundary (k=nf-1) */
            float WL_self[4] = {Qe[0][jn*nn+0],    Qe[1][jn*nn+0],
                                 Qe[2][jn*nn+0],    Qe[3][jn*nn+0]};
            float WR_self[4] = {Qe[0][jn*nn+nn-1], Qe[1][jn*nn+nn-1],
                                 Qe[2][jn*nn+nn-1], Qe[3][jn*nn+nn-1]};
            sanitize_prim(WL_self); sanitize_prim(WR_self);

            { /* left face */
                float Ftmp[4];
                if (ei > 0) {
                    int nb = ej*N+(ei-1);
                    float U_nb[4] = { Q[SIDX(0,nb,jn,nn-1,NEL)], Q[SIDX(1,nb,jn,nn-1,NEL)],
                                      Q[SIDX(2,nb,jn,nn-1,NEL)], Q[SIDX(3,nb,jn,nn-1,NEL)] };
                    float W_nb[4];
                    cons_to_prim_sane(U_nb, W_nb);
                    hllc_n(W_nb, WL_self, 1.f,0.f, Ftmp);
                } else {
                    float Wg[4]; wall_state(WL_self, 0, Wg);
                    hllc_n(Wg, WL_self, 1.f,0.f, Ftmp);
                }
                for (int q=0;q<NVAR;q++) Ff[q][0] = Ftmp[q];
            }
            { /* right face */
                float Ftmp[4];
                if (ei < N-1) {
                    int nb = ej*N+(ei+1);
                    float U_nb[4] = { Q[SIDX(0,nb,jn,0,NEL)], Q[SIDX(1,nb,jn,0,NEL)],
                                      Q[SIDX(2,nb,jn,0,NEL)], Q[SIDX(3,nb,jn,0,NEL)] };
                    float W_nb[4];
                    cons_to_prim_sane(U_nb, W_nb);
                    hllc_n(WR_self, W_nb, 1.f,0.f, Ftmp);
                } else {
                    float Wg[4]; wall_state(WR_self, 0, Wg);
                    hllc_n(WR_self, Wg, 1.f,0.f, Ftmp);
                }
                for (int q=0;q<NVAR;q++) Ff[q][nf-1] = Ftmp[q];
            }

            /* 4. dF/dx|_in = jac * Σ_k Dfs[in][k] * Ff[k]              */
            for (int in=0; in<nn; in++) {
                for (int q=0; q<NVAR; q++) {
                    float dF = 0.f;
                    for (int k=0; k<nf; k++) dF += c_Dfs[in][k]*Ff[q][k];
                    R[q][jn*nn+in] -= jac * dF;
                }
            }
        }

        /* ── y-direction SD flux divergence ─────────────────────────── */
        for (int in=0; in<nn; in++) {

            /* 1. Interpolate conservative column to NF flux nodes       */
            float Qf[NVAR][5];
            for (int k=0; k<nf; k++) {
                float Uk[4];
                for (int q=0; q<NVAR; q++) {
                    float v = 0.f;
                    for (int j=0; j<nn; j++) v += c_Isf[k][j]*Qc[q][j*nn+in];
                    Uk[q] = v;
                }
                float Wk[4];
                cons_to_prim_sane(Uk, Wk);
                for (int q=0;q<NVAR;q++) Qf[q][k]=Wk[q];
            }

            /* 2. Euler flux at interior flux nodes                      */
            float Gf[NVAR][5];
            for (int k=1; k<nf-1; k++) {
                float G4[4];
                euler_flux_n(Qf[0][k],Qf[1][k],Qf[2][k],Qf[3][k], 0.f,1.f, G4);
                for (int q=0; q<NVAR; q++) Gf[q][k] = G4[q];
            }

            /* 3. HLLC at bottom (k=0) and top (k=nf-1) faces           */
            float WB_self[4] = {Qe[0][0*nn+in],      Qe[1][0*nn+in],
                                 Qe[2][0*nn+in],      Qe[3][0*nn+in]};
            float WT_self[4] = {Qe[0][(nn-1)*nn+in], Qe[1][(nn-1)*nn+in],
                                 Qe[2][(nn-1)*nn+in], Qe[3][(nn-1)*nn+in]};
            sanitize_prim(WB_self); sanitize_prim(WT_self);

            { /* bottom face */
                float Gtmp[4];
                if (ej > 0) {
                    int nb = (ej-1)*N+ei;
                    float U_nb[4] = { Q[SIDX(0,nb,nn-1,in,NEL)], Q[SIDX(1,nb,nn-1,in,NEL)],
                                      Q[SIDX(2,nb,nn-1,in,NEL)], Q[SIDX(3,nb,nn-1,in,NEL)] };
                    float W_nb[4];
                    cons_to_prim_sane(U_nb, W_nb);
                    hllc_n(W_nb, WB_self, 0.f,1.f, Gtmp);
                } else {
                    float Wg[4]; wall_state(WB_self, 1, Wg);
                    hllc_n(Wg, WB_self, 0.f,1.f, Gtmp);
                }
                for (int q=0;q<NVAR;q++) Gf[q][0] = Gtmp[q];
            }
            { /* top face */
                float Gtmp[4];
                if (ej < N-1) {
                    int nb = (ej+1)*N+ei;
                    float U_nb[4] = { Q[SIDX(0,nb,0,in,NEL)], Q[SIDX(1,nb,0,in,NEL)],
                                      Q[SIDX(2,nb,0,in,NEL)], Q[SIDX(3,nb,0,in,NEL)] };
                    float W_nb[4];
                    cons_to_prim_sane(U_nb, W_nb);
                    hllc_n(WT_self, W_nb, 0.f,1.f, Gtmp);
                } else {
                    float Wg[4]; wall_state(WT_self, 1, Wg);
                    hllc_n(WT_self, Wg, 0.f,1.f, Gtmp);
                }
                for (int q=0;q<NVAR;q++) Gf[q][nf-1] = Gtmp[q];
            }

            /* 4. dG/dy|_jn = jac * Σ_k Dfs[jn][k] * Gf[k]             */
            for (int jn=0; jn<nn; jn++) {
                for (int q=0; q<NVAR; q++) {
                    float dG = 0.f;
                    for (int k=0; k<nf; k++) dG += c_Dfs[jn][k]*Gf[q][k];
                    R[q][jn*nn+in] -= jac * dG;
                }
            }
        }

        /* ── Write conservative RHS and element lambda_max ───────────── */
        float lam_el = 0.f;
        for (int n=0; n<nn2; n++) {
            float rho=Qe[RHO_ID][n], u=Qe[U_ID][n], v=Qe[V_ID][n], p=Qe[P_ID][n];
            float c=sound_speed(p,rho);
            lam_el=fmaxf(lam_el, fabsf(u)+fabsf(v)+c);
            for (int q=0;q<NVAR;q++)
                RHS[SIDX(q, el, n/nn, n%nn, NEL)] = R[q][n];
        }
        lam[el] = lam_el;
    }
}

/* ── Element-average density kernel (for PPM output) ───────────────────── */
__global__ void compute_rho_avg(
    const float * __restrict__ Q,
    float       * __restrict__ rho_avg,
    int N)
{
    const int NEL = N * N;
    for (int el = blockIdx.x*blockDim.x+threadIdx.x; el < NEL;
             el += blockDim.x*gridDim.x) {
        float s = 0.f;
        for (int n=0; n<NN2; n++) s += Q[SIDX(RHO_ID, el, n/NN, n%NN, NEL)];
        rho_avg[el] = s / NN2;
    }
}

/* ── SSP-RK3 kernels on conservative state ─────────────────────────────── */
__global__ void rk3_fr_s1(float*U1, const float*U0, const float*L, float dt,
                          int NEL)
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
__global__ void rk3_fr_s2(float*U2, const float*U0, const float*U1, const float*L,
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
__global__ void rk3_fr_s3(float*U, const float*U0, const float*U2, const float*L,
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

/* ── Nodal positivity repair on conservative state ─────────────────────── */
__global__ void positivity_limiter(float *Q, int N, float eps_rho, float eps_p)
{
    const int NEL = N * N;
    const int nn  = NN;

    for (int el = blockIdx.x*blockDim.x+threadIdx.x; el < NEL;
             el += blockDim.x*gridDim.x) {
        for (int jn=0; jn<nn; jn++) for (int in=0; in<nn; in++) {
            float Ufix[4] = {
                Q[SIDX(RHO_ID, el, jn, in, NEL)],
                Q[SIDX(U_ID,   el, jn, in, NEL)],
                Q[SIDX(V_ID,   el, jn, in, NEL)],
                Q[SIDX(P_ID,   el, jn, in, NEL)]
            };
            if (Ufix[0] < eps_rho) {
                Ufix[0] = eps_rho;
                Ufix[1] = 0.f;
                Ufix[2] = 0.f;
            }
            float p = pressure_from_cons(Ufix);
            if (p < eps_p)
                Ufix[3] = eps_p/(GAMMA_V-1.f) + 0.5f*(Ufix[1]*Ufix[1] + Ufix[2]*Ufix[2]) / Ufix[0];
            for (int q=0; q<NVAR; q++) Q[SIDX(q, el, jn, in, NEL)] = Ufix[q];
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
        float delta = 2.0f * h;
        float phi   = 0.5f*(1.f+tanhf((0.25f-r)/delta));

        float C[4];
        p2c(0.125f + 10.875f*phi, 0.f, 0.f, 0.1f + 9.9f*phi, C);
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
static void downsample_rho_p(const float*Q, float*rho_img, float*p_img,
                              float*sd_img, const float*sd_el, int N)
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
        sd_img [el] = sd_el ? sd_el[el] : 0.f;
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

/* ═══════════════════════════════════════════════════════════════════════════
 * main
 * ═══════════════════════════════════════════════════════════════════════════ */
int main(int argc, char **argv)
{
    int   N   = (argc > 1) ? atoi(argv[1]) : 50;
    float CFL = 0.1f / ORDER, t_end = 1.0f, h = 1.0f / N;

    /* ── Build flux nodes xf and SD matrices on host ─────────────────── */
    /* xf = { -1, GL_interior[0..NN-2], +1 }  (NF = NN+1 nodes)          */
    double xs[4], xf[5];
    for (int i=0;i<NN;i++) xs[i] = lgl_xi[ORDER-1][i];
    xf[0] = -1.0;
    for (int k=0; k<NN-1; k++) xf[1+k] = gl_interior[ORDER-1][k];
    xf[NF-1] = 1.0;

    double Isf_h[5][4] = {}, Dfs_h[4][5] = {};
    build_Isf(NN, xs, NF, xf, Isf_h);
    build_Dfs(NN, xs, NF, xf, Dfs_h);

    float xf_f[5]={}, Isf_f[5][4]={}, Dfs_f[4][5]={};
    for (int k=0;k<NF;k++) {
        xf_f[k]=(float)xf[k];
        for (int j=0;j<NN;j++) Isf_f[k][j]=(float)Isf_h[k][j];
    }
    for (int i=0;i<NN;i++) for (int k=0;k<NF;k++) Dfs_f[i][k]=(float)Dfs_h[i][k];
    CK(cudaMemcpyToSymbol(c_xf,  xf_f,  sizeof(xf_f)));
    CK(cudaMemcpyToSymbol(c_Isf, Isf_f, sizeof(Isf_f)));
    CK(cudaMemcpyToSymbol(c_Dfs, Dfs_f, sizeof(Dfs_f)));

    float xi_f[8] = {};
    for (int i=0;i<NN;i++) xi_f[i]=(float)lgl_xi[ORDER-1][i];
    CK(cudaMemcpyToSymbol(c_xi, xi_f, NN*sizeof(float)));
    float w_f[8] = {};
    for (int i=0;i<NN;i++) w_f[i]=(float)gll_w[ORDER-1][i];
    CK(cudaMemcpyToSymbol(c_w, w_f, NN*sizeof(float)));

    /* ── GPU info ────────────────────────────────────────────────────── */
    int dev; cudaDeviceProp prop;
    CK(cudaGetDevice(&dev));
    CK(cudaGetDeviceProperties(&prop,dev));
    printf("========================================================\n");
    printf("  Device  : %s\n", prop.name);
    printf("  Spectral Difference p=%d  --  Circular Sod\n", ORDER);
    printf("  N=%d  h=%.5f  NN=%d sol nodes  NF=%d flux nodes  CFL=%.3f\n",
           N, h, NN, NF, CFL);
    printf("  Total DOF: %d\n", N*N*NN2*NVAR);
    printf("========================================================\n");

    /* ── Allocations ─────────────────────────────────────────────────── */
    int    NEL  = N * N;
    size_t szQ  = (size_t)NVAR * NEL * NN2 * sizeof(float);
    size_t sz1  = (size_t)NEL * sizeof(float);

    float *d_Q, *d_Q0, *d_Q1, *d_Q2, *d_RHS, *d_lam, *d_tmp;
    CK(cudaMalloc(&d_Q,   szQ));
    CK(cudaMalloc(&d_Q0,  szQ));
    CK(cudaMalloc(&d_Q1,  szQ));
    CK(cudaMalloc(&d_Q2,  szQ));
    CK(cudaMalloc(&d_RHS, szQ));
    CK(cudaMalloc(&d_lam, sz1));
    CK(cudaMalloc(&d_tmp, GS_NBLK*sizeof(float)));

    float *h_Q   = (float*)malloc(szQ);
    float *h_rho = (float*)malloc(NEL*sizeof(float));
    float *h_p   = (float*)malloc(NEL*sizeof(float));
    float *h_sd  = (float*)malloc(NEL*sizeof(float));   /* unused, kept for write_ppm_2x1 */

    /* colormap range: computed once from frame 0, reused for all frames */
    float g_rho_min=0.f, g_rho_max=1.f, g_p_min=0.f, g_p_max=1.f;
    int   g_range_set = 0;

#define WRITE_FRAME(fi) do { \
    CK(cudaDeviceSynchronize()); \
    CK(cudaMemcpy(h_Q, d_Q, szQ, cudaMemcpyDeviceToHost)); \
    downsample_rho_p(h_Q, h_rho, h_p, h_sd, NULL, N); \
    { float _lo=h_rho[0],_hi=h_rho[0],_plo=h_p[0],_phi=h_p[0]; \
      for(int _i=1;_i<NEL;_i++){_lo=fminf(_lo,h_rho[_i]);_hi=fmaxf(_hi,h_rho[_i]); \
                                  _plo=fminf(_plo,h_p[_i]);_phi=fmaxf(_phi,h_p[_i]);} \
      printf("  [frame %d] rho=[%.4f,%.4f]  p=[%.4f,%.4f]\n",(fi),_lo,_hi,_plo,_phi); } \
    if (!g_range_set) { \
        image_range(h_rho, NEL, &g_rho_min, &g_rho_max); \
        image_range(h_p,   NEL, &g_p_min,   &g_p_max); \
        g_range_set = 1; \
    } \
    char _fn[64]; sprintf(_fn,"fr_sod_%04d.ppm",(fi)); \
    write_ppm_2x1(_fn, h_rho, g_rho_min, g_rho_max, h_p, g_p_min, g_p_max, N); \
} while(0)

    /* ── IC ─────────────────────────────────────────────────────────── */
    ic_kernel<<<GS_NBLK,BLOCK1D>>>(d_Q, N, h);
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
        compute_sd_rhs<<<GS_NBLK,BLOCK1D>>>(d_Q, d_RHS, d_lam, N, h);
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
        /* d_RHS = L(Q0) already computed above */
        rk3_fr_s1<<<GS_NBLK,BLOCK1D>>>(d_Q1, d_Q0, d_RHS, dt, NEL);
        positivity_limiter<<<GS_NBLK,BLOCK1D>>>(d_Q1, N, 1e-6f, 1e-6f);

        /* RK3 stage 2: Q2 = 3/4*Q0 + 1/4*(Q1 + dt*L(Q1)) */
        compute_sd_rhs<<<GS_NBLK,BLOCK1D>>>(d_Q1, d_RHS, d_lam, N, h);
        rk3_fr_s2<<<GS_NBLK,BLOCK1D>>>(d_Q2, d_Q0, d_Q1, d_RHS, dt, NEL);
        positivity_limiter<<<GS_NBLK,BLOCK1D>>>(d_Q2, N, 1e-6f, 1e-6f);

        /* RK3 stage 3: Q = 1/3*Q0 + 2/3*(Q2 + dt*L(Q2)) */
        compute_sd_rhs<<<GS_NBLK,BLOCK1D>>>(d_Q2, d_RHS, d_lam, N, h);
        rk3_fr_s3<<<GS_NBLK,BLOCK1D>>>(d_Q, d_Q0, d_Q2, d_RHS, dt, NEL);
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
    printf("  Frames: fr_sod_0000.ppm .. fr_sod_%04d.ppm\n",N_FRAMES);

    free(h_Q); free(h_rho); free(h_p); free(h_sd);
    cudaFree(d_Q); cudaFree(d_Q0); cudaFree(d_Q1); cudaFree(d_Q2);
    cudaFree(d_RHS); cudaFree(d_lam); cudaFree(d_tmp);
    return 0;
}
