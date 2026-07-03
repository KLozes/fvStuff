/*
 * rt_dg1_euler_2d.cu
 * ─────────────────────────────────────────────────────────────────────────
 * Degree-1 collocated Gauss-point DG for the 2D compressible Euler equations.
 * (Formerly a Raviart-Thomas mixed element; all DOFs now live on the 2x2 Gauss
 *  grid, so this is a standard Q1 Gauss-collocation DG with upwind faces.)
 *
 * Collocated nodal element on the reference cell [-1,1]^2  (g = 1/sqrt(3)):
 *
 *      +---------------+          *  blue  : all DOFs, interior 2x2 Gauss
 *      |   *       *   |
 *      |               |          ALL four conserved fields share the same
 *      |               |          2x2 Gauss nodal grid (x_i,y_j)=(+-g,+-g);
 *      |   *       *   |          face states are EXTRAPOLATED to xi/eta=+-1
 *      +---------------+          with the Gauss weights wGp/wGm.
 *
 *   rho   : Q1 nodal on 2x2 Gauss        (x_i,y_j) = (+-g,+-g)     4 DOFs
 *   E     : Q1 nodal on 2x2 Gauss                                  4 DOFs
 *   rho*u : Q1 nodal on 2x2 Gauss                                  4 DOFs
 *   rho*v : Q1 nodal on 2x2 Gauss                                  4 DOFs
 *                                                            total 16 DOFs/cell
 *
 *   Riemann   : compile-time choice (2-pt Gauss surface quad) — see RIEMANN macro
 *   Volume    : 2x2 Gauss quadrature of the physical fluxes
 *   Mass^-1   : exact, fully diagonal (Gauss nodal x Gauss quad) = (2/h)^2
 *   Time      : SSP-RK3
 *   Limiter   : Zhang-Shu positivity (density rescale + pressure bisection)
 *   BC        : periodic
 *
 * Compile (flux selected at build time, default HLLC):
 *   nvcc -O3 -arch=native --expt-relaxed-constexpr -o rt_dg1 rt_dg1_euler_2d.cu -lm
 *   nvcc ... -DRIEMANN=RIE_BARMAT  ...   low-Mach matrix-dissipation flux
 *   nvcc ... -DRIEMANN=RIE_CENTRAL ...   non-dissipative central flux
 *   nvcc ... -DRECON=REC_QUADRATIC ...   unlimited parabolic reconstruction
 *   nvcc ... -DRECON=REC_LINEAR    ...   unlimited centred-linear reconstruction
 *                                        (default REC_LIMITED, van Leer)
 * Run:
 *   ./rt_dg1 200                  circular Sod shock tube on [0,1]^2
 *   ./rt_dg1 200 lmv [eps]        low-Mach (Gresho) vortex, default eps=0.1
 *   ./rt_dg1 200 dsl [Ma]         doubly periodic shear layer, default Ma=0.1
 *   ./rt_dg1 200 va  [Mref]       vortex-acoustic wave interaction, default 0.1
 * ───────────────────────────────────────────────────────────────────────── */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <cuda_runtime.h>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

/* ── constants ──────────────────────────────────────────────────────────── */
#define GAMMA_V 1.4
#define BLOCK1D 256
#define GS_NBLK 256
#define G_GHOST 1                 /* one ghost layer (compact stencil) */
#define NDOF    16                /* DOFs per cell                     */
#define real    float
#define BARSUKOW_FLOOR 0.1        /* low-Mach dissipation floor f for barmat flux */

/* ── Riemann solver: compile-time selection ──────────────────────────────────
 * Build with -DRIEMANN=RIE_HLLC | RIE_BARMAT | RIE_CENTRAL (default HLLC).      */
#define RIE_HLLC    0
#define RIE_BARMAT  1
#define RIE_CENTRAL 2
#ifndef RIEMANN
#define RIEMANN RIE_HLLC
#endif
#if   RIEMANN == RIE_HLLC
#define RIEMANN_NAME "HLLC"
#elif RIEMANN == RIE_BARMAT
#define RIEMANN_NAME "barmat (low-Mach)"
#elif RIEMANN == RIE_CENTRAL
#define RIEMANN_NAME "central (non-dissipative)"
#else
#error "RIEMANN must be RIE_HLLC, RIE_BARMAT, or RIE_CENTRAL"
#endif

/* ── Face reconstruction: compile-time selection ─────────────────────────────
 * Build with -DRECON=REC_LIMITED | REC_QUADRATIC | REC_LINEAR (default LIMITED).
 *   REC_LIMITED   : flux-limited (van Leer) reconstruction  — TVD, default
 *   REC_QUADRATIC : unlimited 3-node parabolic reconstruction (no limiter)
 *   REC_LINEAR    : unlimited 2-node centred linear reconstruction
 * The two unlimited modes are *linear schemes* (no solution-dependent limiting);
 * they recover full design accuracy in smooth flow but are not monotone at
 * shocks.  REC_LIMITED is the robust choice for the Sod test.                 */
#define REC_LIMITED   0
#define REC_QUADRATIC 1
#define REC_LINEAR    2
#ifndef RECON
#define RECON REC_QUADRATIC
#endif
#if   RECON == REC_LIMITED
#define RECON_NAME "limited (van Leer)"
#elif RECON == REC_QUADRATIC
#define RECON_NAME "quadratic (unlimited)"
#elif RECON == REC_LINEAR
#define RECON_NAME "linear (unlimited)"
#else
#error "RECON must be REC_LIMITED, REC_QUADRATIC, or REC_LINEAR"
#endif

/* DOF channel layout in the padded SoA array Q[NDOF][Np][Np].
 * For all fields (RHO/E/MX/MY):  a = Gauss-x index (0,1),  b = Gauss-y index */
#define RHO(a,b)  (0  + ((b)*2 + (a)))
#define EN(a,b)   (4  + ((b)*2 + (a)))
#define MX(a,b)   (8  + ((b)*2 + (a)))
#define MY(a,b)   (12 + ((b)*2 + (a)))

#define IDX_P(q,jp,ip,Np) ((size_t)(q)*(Np)*(Np) + (size_t)(jp)*(Np) + (ip))

#define CK(x) do { cudaError_t _e=(x); if(_e!=cudaSuccess){ \
    fprintf(stderr,"CUDA error %s:%d: %s\n",__FILE__,__LINE__,\
    cudaGetErrorString(_e)); exit(1);} } while(0)

/* ── reference-cell nodal constants (g = 1/sqrt(3)) ─────────────────────────
 *   wGp[a]  = ell_a^Gauss(+1)   (extrapolate Gauss basis to the +1 face)
 *   wGm[a]  = ell_a^Gauss(-1)
 *   DGv[a]  = d/dxi of the two Gauss Lagrange basis fns (constant, = -+sqrt3/2) */
__constant__ real wGp[2] = {-0.3660254037844386, 1.3660254037844386};
__constant__ real wGm[2] = { 1.3660254037844386,-0.3660254037844386};
__constant__ real DGv[2] = {-0.8660254037844386, 0.8660254037844386};


/* ═══════════════════════════════════════════════════════════════════════════
 * Thermodynamics + HLLC Riemann solver
 * ═══════════════════════════════════════════════════════════════════════════ */
__device__ __forceinline__ real sound_speed(real rho, real p)
{ return sqrt(GAMMA_V * fmax(p,1e-30) / fmax(rho,1e-30)); }

/* primitives W=[rho,u,v,p] from conserved (rho,mx,my,E) */
__device__ __forceinline__ void
prim_of(real rho, real mx, real my, real E, real W[4])
{
    rho = fmax(rho, 1e-14);
    real u = mx/rho, v = my/rho;
    real p = (GAMMA_V-1.)*(E - 0.5*(mx*mx + my*my)/rho);
    W[0]=rho; W[1]=u; W[2]=v; W[3]=fmax(p,1e-14);
}

/* HLLC flux in direction n=(nx,ny); F=[mass,xmom,ymom,energy] (Cartesian) */
__device__ void
hllc_n(const real WL[4], const real WR[4], real nx, real ny, real F[4])
{
    real rL=WL[0],uL=WL[1],vL=WL[2],pL=WL[3];
    real rR=WR[0],uR=WR[1],vR=WR[2],pR=WR[3];
    real unL=uL*nx+vL*ny, unR=uR*nx+vR*ny;
    real utL=-uL*ny+vL*nx, utR=-uR*ny+vR*nx;
    real cL=sound_speed(rL,pL), cR=sound_speed(rR,pR);
    real EL=pL/(GAMMA_V-1.)+0.5*rL*(uL*uL+vL*vL);
    real ER=pR/(GAMMA_V-1.)+0.5*rR*(uR*uR+vR*vR);
    real SL=fmin(unL-cL,unR-cR), SR=fmax(unL+cL,unR+cR);
    real FnL[4]={rL*unL, rL*unL*unL+pL, rL*unL*utL, (EL+pL)*unL};
    real FnR[4]={rR*unR, rR*unR*unR+pR, rR*unR*utR, (ER+pR)*unR};
    real Fn[4];
    if (SL>=0.)      { Fn[0]=FnL[0];Fn[1]=FnL[1];Fn[2]=FnL[2];Fn[3]=FnL[3]; }
    else if (SR<=0.) { Fn[0]=FnR[0];Fn[1]=FnR[1];Fn[2]=FnR[2];Fn[3]=FnR[3]; }
    else {
        real dL=rL*(SL-unL), dR=rR*(SR-unR);
        real Ss=(pR-pL+dL*unL-dR*unR)/(dL-dR);
        real UL[4]={rL,rL*unL,rL*utL,EL};
        real UR[4]={rR,rR*unR,rR*utR,ER};
        if (Ss>=0.) {
            real f=rL*(SL-unL)/(SL-Ss);
            real Us[4]={f, f*Ss, f*utL, f*(EL/rL+(Ss-unL)*(Ss+pL/(rL*(SL-unL))))};
            for(int q=0;q<4;q++) Fn[q]=FnL[q]+SL*(Us[q]-UL[q]);
        } else {
            real f=rR*(SR-unR)/(SR-Ss);
            real Us[4]={f, f*Ss, f*utR, f*(ER/rR+(Ss-unR)*(Ss+pR/(rR*(SR-unR))))};
            for(int q=0;q<4;q++) Fn[q]=FnR[q]+SR*(Us[q]-UR[q]);
        }
    }
    F[0]=Fn[0];
    F[1]=Fn[1]*nx-Fn[2]*ny;   /* rotate normal/tangential -> Cartesian */
    F[2]=Fn[1]*ny+Fn[2]*nx;
    F[3]=Fn[3];
}

/* ── Barsukow matrix-dissipation flux (Barsukow 2018, eq. 7) — low-Mach ──────
 * Roe eigenvector decomposition on conservative jumps, but the entropy/shear
 * eigenvalues get a floor f (= BARSUKOW_FLOOR) instead of 0.  This makes the
 * momentum dissipation O(Ma^2) at low Mach (vs O(1) for Roe/HLLC), so it does
 * not over-damp slow vortical modes.  CFL scales with Δx.                     */
__device__ void
barmat_n(const real WL[4], const real WR[4], real nx, real ny, real F[4])
{
    real rL=WL[0],uL=WL[1],vL=WL[2],pL=WL[3];
    real rR=WR[0],uR=WR[1],vR=WR[2],pR=WR[3];
    real cL=sound_speed(rL,pL), cR=sound_speed(rR,pR);

    real sqrL=sqrt(rL), sqrR=sqrt(rR), denom=sqrL+sqrR;
    real uRoe=(sqrL*uL+sqrR*uR)/denom, vRoe=(sqrL*vL+sqrR*vR)/denom;
    real HL=(pL/(GAMMA_V-1.)+0.5*rL*(uL*uL+vL*vL)+pL)/rL;
    real HR=(pR/(GAMMA_V-1.)+0.5*rR*(uR*uR+vR*vR)+pR)/rR;
    real HRoe=(sqrL*HL+sqrR*HR)/denom;
    real c2Roe=(GAMMA_V-1.)*(HRoe-0.5*(uRoe*uRoe+vRoe*vRoe));
    real cRoe=sqrt(c2Roe>1e-14?c2Roe:1e-14);

    real unL=uL*nx+vL*ny, utL=-uL*ny+vL*nx;
    real unR=uR*nx+vR*ny, utR=-uR*ny+vR*nx;
    real unRoe=uRoe*nx+vRoe*ny, utRoe=-uRoe*ny+vRoe*nx;

    real EL=pL/(GAMMA_V-1.)+0.5*rL*(uL*uL+vL*vL);
    real ER=pR/(GAMMA_V-1.)+0.5*rR*(uR*uR+vR*vR);
    real FnL[4]={rL*unL, rL*unL*unL+pL, rL*unL*utL, (EL+pL)*unL};
    real FnR[4]={rR*unR, rR*unR*unR+pR, rR*unR*utR, (ER+pR)*unR};

    real drho=rR-rL, dun=unR-unL, dut=utR-utL, dp=pR-pL, rRoe=sqrL*sqrR;
    real b1=(dp-rRoe*cRoe*dun)/(2.*c2Roe);
    real b2=drho-dp/c2Roe;
    real b3=rRoe*dut;
    real b4=(dp+rRoe*cRoe*dun)/(2.*c2Roe);

    real lam1=unRoe-cRoe, lam2=unRoe, lam3=unRoe, lam4=unRoe+cRoe;
    real eps1=fmax(0.,2.*((unR-cR)-(unL-cL)));
    real eps4=fmax(0.,2.*((unR+cR)-(unL+cL)));
    real al1=fabs(lam1); if(al1<eps1*0.5) al1=(lam1*lam1+eps1*eps1*0.25)/eps1;
    real al2=fabs(lam2)+BARSUKOW_FLOOR;   /* entropy: |u_n| + f */
    real al3=fabs(lam3)+BARSUKOW_FLOOR;   /* shear:   |u_n| + f */
    real al4=fabs(lam4); if(al4<eps4*0.5) al4=(lam4*lam4+eps4*eps4*0.25)/eps4;

    real d_rho=al1*b1            + al2*b2          + al4*b4;
    real d_un =al1*b1*(unRoe-cRoe)+ al2*b2*unRoe    + al4*b4*(unRoe+cRoe);
    real d_ut =al1*b1*utRoe      + al2*b2*utRoe    + al3*b3 + al4*b4*utRoe;
    real d_E  =al1*b1*(HRoe-unRoe*cRoe)
              + al2*b2*0.5*(unRoe*unRoe+utRoe*utRoe)
              + al3*b3*utRoe + al4*b4*(HRoe+unRoe*cRoe);

    real Fn[4];
    Fn[0]=0.5*(FnL[0]+FnR[0])-0.5*d_rho;
    Fn[1]=0.5*(FnL[1]+FnR[1])-0.5*d_un;
    Fn[2]=0.5*(FnL[2]+FnR[2])-0.5*d_ut;
    Fn[3]=0.5*(FnL[3]+FnR[3])-0.5*d_E;
    F[0]=Fn[0];
    F[1]=Fn[1]*nx-Fn[2]*ny;
    F[2]=Fn[1]*ny+Fn[2]*nx;
    F[3]=Fn[3];
}

/* Pure central (non-dissipative) normal flux — for well-balancing diagnostics. */
__device__ void
central_n(const real WL[4], const real WR[4], real nx, real ny, real F[4])
{
    real rL=WL[0],uL=WL[1],vL=WL[2],pL=WL[3];
    real rR=WR[0],uR=WR[1],vR=WR[2],pR=WR[3];
    real unL=uL*nx+vL*ny, unR=uR*nx+vR*ny;
    real EL=pL/(GAMMA_V-1.)+0.5*rL*(uL*uL+vL*vL);
    real ER=pR/(GAMMA_V-1.)+0.5*rR*(uR*uR+vR*vR);
    real FL[4]={rL*unL, rL*unL*uL+pL*nx, rL*unL*vL+pL*ny, (EL+pL)*unL};
    real FR[4]={rR*unR, rR*unR*uR+pR*nx, rR*unR*vR+pR*ny, (ER+pR)*unR};
    for(int q=0;q<4;q++) F[q]=0.5*(FL[q]+FR[q]);
}

/* compile-time flux dispatch (see RIEMANN macro above) */
__device__ __forceinline__ void
riemann_n(const real WL[4], const real WR[4], real nx, real ny, real F[4])
{
#if   RIEMANN == RIE_CENTRAL
    central_n(WL,WR,nx,ny,F);
#elif RIEMANN == RIE_BARMAT
    barmat_n (WL,WR,nx,ny,F);
#else
    hllc_n   (WL,WR,nx,ny,F);
#endif
}

/* ═══════════════════════════════════════════════════════════════════════════
 * Boundary conditions — periodic, copies all 16 channels
 * ═══════════════════════════════════════════════════════════════════════════ */
__global__ void
apply_bc(real * __restrict__ Qp, int N, int Np)
{
    int Np2=Np*Np;
    for (int k=blockIdx.x*blockDim.x+threadIdx.x; k<Np2; k+=blockDim.x*gridDim.x){
        int jp=k/Np, ip=k%Np;
        bool gx=(ip<G_GHOST||ip>=N+G_GHOST), gy=(jp<G_GHOST||jp>=N+G_GHOST);
        if (!gx && !gy) continue;
        int ir=ip, jr=jp;
        if (ip<G_GHOST)      ir=ip+N;
        if (ip>=N+G_GHOST)   ir=ip-N;
        if (jp<G_GHOST)      jr=jp+N;
        if (jp>=N+G_GHOST)   jr=jp-N;
        for (int q=0;q<NDOF;q++) Qp[IDX_P(q,jp,ip,Np)]=Qp[IDX_P(q,jr,ir,Np)];
    }
}

/* Load the 16 DOFs of cell (jp,ip) into a local array. */
__device__ __forceinline__ void
load_cell(const real * __restrict__ Qp, int jp, int ip, int Np, real c[NDOF])
{
    for (int q=0;q<NDOF;q++) c[q]=Qp[IDX_P(q,jp,ip,Np)];
}

/* ── nodal face reconstruction (normal direction) ───────────────────────────
 * The 4 nearest normal nodes at a face are  far near | nnear nfar  (Gauss +-g),
 * located in the owner's reference coords at  -g, +g, 2-g  for far/near/nnear.
 * The face sits at xi=+1, midway between 'near' and 'nnear'.  Three build-time
 * reconstructions (RECON macro above) produce the owner-side trace u_face:
 *
 *   REC_LINEAR    : 2-node centred linear,  u_face = 1/2(near + nnear).  Exact
 *                   for linear data, non-dissipative, but a *linear* scheme.
 *   REC_QUADRATIC : the parabola through (far,near,nnear) sampled at xi=+1,
 *                   u_face = near + Q0*(near-far) + Q1*(nnear-near).  Third-order
 *                   in smooth flow; reduces to the linear value when r=1.  Still
 *                   an unlimited (linear) scheme — no monotonicity at shocks.
 *   REC_LIMITED   : MUSCL form u_face = near + 1/2 phi(r) (nnear-near) with van
 *                   Leer phi.  The smoothness ratio r=(sqrt3-1)(near-far)/(nnear
 *                   -near) is on GRADIENTS, so r=1 is a linear field: exact for
 *                   linear data, centred (non-dissipative) in smooth flow, and
 *                   first order at extrema/shocks (TVD).
 *
 * Q0,Q1 are the parabola's xi=+1 Lagrange weights on the (-g,+g,2-g) stencil:
 *   Q0 = sqrt(3)/3 - 1/2 = -ell_far(+1),   Q1 = (1+1/sqrt3)/4 = ell_nnear(+1). */
__device__ __forceinline__ real
recon_face(real u_far, real u_near, real n_near)
{
    real dp=n_near-u_near;
#if   RECON == REC_LINEAR
    (void)u_far;
    return u_near + 0.5*dp;                   /* = 1/2(u_near + n_near) */
#elif RECON == REC_QUADRATIC
    const real Q0=0.07735026918962576;        /* sqrt(3)/3 - 1/2        */
    const real Q1=0.39433756729740643;        /* (1+1/sqrt3)/4          */
    real dm=u_near-u_far;
    return u_near + Q0*dm + Q1*dp;
#else  /* REC_LIMITED */
    const real GR=0.7320508075688772;        /* sqrt(3)-1 : node-spacing ratio */
    real dm=u_near-u_far;
    real r  = GR*dm/(dp + copysign((real)1e-30, dp));
    real ar = fabs(r);
    real phi= (r+ar)/(1.0+ar);               /* van Leer limiter */
    return u_near + 0.5*phi*dp;
#endif
}

/* Convert a cell's 16 conserved DOFs to primitives, stored in the same layout:
 * RHO slot -> rho, MX -> u, MY -> v, EN -> p.  Face states are reconstructed
 * from these so the limiter acts on (rho,u,v,p) rather than (rho,mx,my,E). */
__device__ __forceinline__ void
to_prim_cell(const real c[NDOF], real cp[NDOF])
{
    for(int b=0;b<2;b++) for(int a=0;a<2;a++){
        real W[4]; prim_of(c[RHO(a,b)],c[MX(a,b)],c[MY(a,b)],c[EN(a,b)],W);
        cp[RHO(a,b)]=W[0]; cp[MX(a,b)]=W[1]; cp[MY(a,b)]=W[2]; cp[EN(a,b)]=W[3];
    }
}

/* ═══════════════════════════════════════════════════════════════════════════
 * RHS — degree-1 Gauss-collocation DG spatial operator (scatter form, atomics)
 *
 * Each interior face flux is evaluated exactly ONCE (by the cell that owns the
 * face on its + side) and atomically scattered to both adjacent cells: the owner
 * gets the -wGp face contribution, the neighbour the +wGm one.  This halves the
 * Riemann-solver work and makes the operator discretely conservative by
 * construction (the two cells share the identical flux value).
 *
 * Requires RHS to be zeroed before launch (the host memsets d_L each stage).
 * ═══════════════════════════════════════════════════════════════════════════ */
__global__ void
compute_rhs(const real * __restrict__ Qp, real * __restrict__ RHS,
            real * __restrict__ lam_out, int N, int Np, real h)
{
    int N2=N*N;
    const real hh  = 0.5*h;      /* both volume & surface carry a (h/2) factor */
    const real inv = 4.0/(h*h);  /* diagonal inverse mass (2/h)^2             */
    for (int k=blockIdx.x*blockDim.x+threadIdx.x; k<N2; k+=blockDim.x*gridDim.x){
        int j=k/N, i=k%N, jp=j+G_GHOST, ip=i+G_GHOST;

        real c[NDOF], e[NDOF], t[NDOF];
        load_cell(Qp, jp,   ip,   Np, c);   /* this cell           */
        load_cell(Qp, jp,   ip+1, Np, e);   /* east  (right)       */
        load_cell(Qp, jp+1, ip,   Np, t);   /* north (top)         */

        /* primitive (rho,u,v,p) node values — faces are reconstructed from these */
        real cp[NDOF],ep[NDOF],tp[NDOF];
        to_prim_cell(c,cp); to_prim_cell(e,ep); to_prim_cell(t,tp);

        /* self accumulator: volume + this cell's RIGHT/TOP face contributions.
         * The LEFT/BOTTOM face terms are scattered in by the west/south threads. */
        real Srho[2][2]={{0,0},{0,0}}, Sen[2][2]={{0,0},{0,0}};
        real Smx [2][2]={{0,0},{0,0}}, Smy[2][2]={{0,0},{0,0}};
        /* east neighbour accumulator (our RIGHT face = east's LEFT face)        */
        real Erho[2][2]={{0,0},{0,0}}, Een[2][2]={{0,0},{0,0}};
        real Emx [2][2]={{0,0},{0,0}}, Emy[2][2]={{0,0},{0,0}};
        /* north neighbour accumulator (our TOP face = north's BOTTOM face)      */
        real Nrho[2][2]={{0,0},{0,0}}, Nen[2][2]={{0,0},{0,0}};
        real Nmx [2][2]={{0,0},{0,0}}, Nmy[2][2]={{0,0},{0,0}};

        /* ── VOLUME term : 2x2 Gauss points (a=x-index, b=y-index) ─────────
         * All four fields are nodal on this same Gauss grid, so the state is
         * read directly and every field uses the Gauss test-function deriv. */
        for (int b=0;b<2;b++) for (int a=0;a<2;a++) {
            real rho=c[RHO(a,b)], En=c[EN(a,b)];
            real mx = c[MX(a,b)], my=c[MY(a,b)];
            real W[4]; prim_of(rho,mx,my,En,W);
            real u=W[1], v=W[2], p=W[3];
            real F[4]={mx, mx*u+p, mx*v,      (En+p)*u};
            real G[4]={my, my*u,   my*v+p,    (En+p)*v};

            /* +int(dx_psi F + dy_psi G), Gauss basis derivative in each dir */
            for(int ap=0;ap<2;ap++){
                Srho[ap][b] += hh*DGv[ap]*F[0];
                Smx [ap][b] += hh*DGv[ap]*F[1];
                Smy [ap][b] += hh*DGv[ap]*F[2];
                Sen [ap][b] += hh*DGv[ap]*F[3];
            }
            for(int bp=0;bp<2;bp++){
                Srho[a][bp] += hh*DGv[bp]*G[0];
                Smx [a][bp] += hh*DGv[bp]*G[1];
                Smy [a][bp] += hh*DGv[bp]*G[2];
                Sen [a][bp] += hh*DGv[bp]*G[3];
            }
        }

        /* ── RIGHT face (xi=+1, n=+x): this cell vs east — computed ONCE ──
         * Flux Fn is scattered to self (-wGp) and to the east cell (+wGm),
         * whose LEFT face sees exactly -Fn (conservation), i.e. +wGm*Fn.
         * Parabolic nodal traces: owner uses {c0,c1,e0}, east uses {c1,e0,e1}. */
        #define RECL(o0,o1,n0) recon_face((o0),(o1),(n0))
        #define RECR(o1,n0,n1) recon_face((n1),(n0),(o1))
        for (int b=0;b<2;b++) {
            real Fn[4];
            real WL[4]={ RECL(cp[RHO(0,b)],cp[RHO(1,b)],ep[RHO(0,b)]),
                         RECL(cp[MX (0,b)],cp[MX (1,b)],ep[MX (0,b)]),
                         RECL(cp[MY (0,b)],cp[MY (1,b)],ep[MY (0,b)]),
                         RECL(cp[EN (0,b)],cp[EN (1,b)],ep[EN (0,b)]) };
            real WR[4]={ RECR(cp[RHO(1,b)],ep[RHO(0,b)],ep[RHO(1,b)]),
                         RECR(cp[MX (1,b)],ep[MX (0,b)],ep[MX (1,b)]),
                         RECR(cp[MY (1,b)],ep[MY (0,b)],ep[MY (1,b)]),
                         RECR(cp[EN (1,b)],ep[EN (0,b)],ep[EN (1,b)]) };
            riemann_n(WL,WR,1.,0.,Fn);
            for(int ap=0;ap<2;ap++){
                Srho[ap][b]-=hh*wGp[ap]*Fn[0];  Smx[ap][b]-=hh*wGp[ap]*Fn[1];
                Smy [ap][b]-=hh*wGp[ap]*Fn[2];  Sen[ap][b]-=hh*wGp[ap]*Fn[3];
                Erho[ap][b]+=hh*wGm[ap]*Fn[0];  Emx[ap][b]+=hh*wGm[ap]*Fn[1];
                Emy [ap][b]+=hh*wGm[ap]*Fn[2];  Een[ap][b]+=hh*wGm[ap]*Fn[3];
            }
        }

        /* ── TOP face (eta=+1, n=+y): this cell vs north — computed ONCE ──
         * Scattered to self (-wGp) and to north's BOTTOM face (+wGm).
         * Parabolic nodal traces: owner uses {c0,c1,t0}, north uses {c1,t0,t1}. */
        for (int a=0;a<2;a++) {
            real Gn[4];
            real WL[4]={ RECL(cp[RHO(a,0)],cp[RHO(a,1)],tp[RHO(a,0)]),
                         RECL(cp[MX (a,0)],cp[MX (a,1)],tp[MX (a,0)]),
                         RECL(cp[MY (a,0)],cp[MY (a,1)],tp[MY (a,0)]),
                         RECL(cp[EN (a,0)],cp[EN (a,1)],tp[EN (a,0)]) };
            real WR[4]={ RECR(cp[RHO(a,1)],tp[RHO(a,0)],tp[RHO(a,1)]),
                         RECR(cp[MX (a,1)],tp[MX (a,0)],tp[MX (a,1)]),
                         RECR(cp[MY (a,1)],tp[MY (a,0)],tp[MY (a,1)]),
                         RECR(cp[EN (a,1)],tp[EN (a,0)],tp[EN (a,1)]) };
            riemann_n(WL,WR,0.,1.,Gn);
            for(int bp=0;bp<2;bp++){
                Srho[a][bp]-=hh*wGp[bp]*Gn[0];  Smx[a][bp]-=hh*wGp[bp]*Gn[1];
                Smy [a][bp]-=hh*wGp[bp]*Gn[2];  Sen[a][bp]-=hh*wGp[bp]*Gn[3];
                Nrho[a][bp]+=hh*wGm[bp]*Gn[0];  Nmx[a][bp]+=hh*wGm[bp]*Gn[1];
                Nmy [a][bp]+=hh*wGm[bp]*Gn[2];  Nen[a][bp]+=hh*wGm[bp]*Gn[3];
            }
        }
        #undef RECL
        #undef RECR

        /* ── apply inverse mass matrix and scatter dQ/dt with atomics ─────
         * Gauss-nodal basis + 2-pt Gauss quad => diagonal mass in both dirs.
         * Periodic wrap maps the +face neighbours back into the interior. */
        int ke = j*N + ((i+1==N)?0:i+1);   /* east  cell (periodic) */
        int kn = ((j+1==N)?0:j+1)*N + i;   /* north cell (periodic) */
        for (int b=0;b<2;b++) for (int a=0;a<2;a++) {
            atomicAdd(&RHS[RHO(a,b)*N2+k ], inv*Srho[a][b]);
            atomicAdd(&RHS[EN (a,b)*N2+k ], inv*Sen [a][b]);
            atomicAdd(&RHS[MX (a,b)*N2+k ], inv*Smx [a][b]);
            atomicAdd(&RHS[MY (a,b)*N2+k ], inv*Smy [a][b]);

            atomicAdd(&RHS[RHO(a,b)*N2+ke], inv*Erho[a][b]);
            atomicAdd(&RHS[EN (a,b)*N2+ke], inv*Een [a][b]);
            atomicAdd(&RHS[MX (a,b)*N2+ke], inv*Emx [a][b]);
            atomicAdd(&RHS[MY (a,b)*N2+ke], inv*Emy [a][b]);

            atomicAdd(&RHS[RHO(a,b)*N2+kn], inv*Nrho[a][b]);
            atomicAdd(&RHS[EN (a,b)*N2+kn], inv*Nen [a][b]);
            atomicAdd(&RHS[MX (a,b)*N2+kn], inv*Nmx [a][b]);
            atomicAdd(&RHS[MY (a,b)*N2+kn], inv*Nmy [a][b]);
        }

        /* ── CFL spectral radius from the cell-mean state ───────────────── */
        real rb=0.,mxb=0.,myb=0.,Eb=0.;
        for(int q=0;q<4;q++){ rb+=c[0+q]; Eb+=c[4+q]; mxb+=c[8+q]; myb+=c[12+q]; }
        rb*=0.25; mxb*=0.25; myb*=0.25; Eb*=0.25;
        real W[4]; prim_of(rb,mxb,myb,Eb,W);
        lam_out[k]=fabs(W[1])+fabs(W[2])+sound_speed(W[0],W[3]);
    }
}

/* ═══════════════════════════════════════════════════════════════════════════
 * Zhang-Shu positivity limiter (post-stage)
 *
 * The flux operator evaluates the solution not only at the 4 interior Gauss
 * points but also at the 8 face quadrature points, where ALL fields are
 * *extrapolated* (Gauss basis sampled at +-1, weights -0.366/1.366).  Positivity
 * must hold at all of these.  We therefore build the full 12-point set and:
 *   1. rescale rho   toward the mean so rho   >= eps_r at every point
 *   2. rescale (rho,m,E) toward the mean so p >= eps_p at every point
 * ═══════════════════════════════════════════════════════════════════════════ */
#define NPT 12
/* Fill the conserved state (rho,mx,my,E) at the 12 evaluation points. */
__device__ __forceinline__ void
eval_points(const real c[NDOF], real rp[NPT], real mxp[NPT],
            real myp[NPT], real Ep[NPT])
{
    int n=0;
    /* 4 interior Gauss points (a=x,b=y) — all fields are nodal here */
    for (int b=0;b<2;b++) for (int a=0;a<2;a++,n++) {
        rp[n]=c[RHO(a,b)]; Ep[n]=c[EN(a,b)];
        mxp[n]=c[MX(a,b)]; myp[n]=c[MY(a,b)];
    }
    /* right / left faces (eta = Gauss b) — extrapolate all fields in xi */
    for (int side=0;side<2;side++){
        const real *wx=(side==0)?wGp:wGm;   /* +1 then -1 */
        for (int b=0;b<2;b++,n++){
            rp[n]=wx[0]*c[RHO(0,b)]+wx[1]*c[RHO(1,b)];
            Ep[n]=wx[0]*c[EN(0,b)] +wx[1]*c[EN(1,b)];
            mxp[n]=wx[0]*c[MX(0,b)]+wx[1]*c[MX(1,b)];
            myp[n]=wx[0]*c[MY(0,b)]+wx[1]*c[MY(1,b)];
        }
    }
    /* top / bottom faces (xi = Gauss a) — extrapolate all fields in eta */
    for (int side=0;side<2;side++){
        const real *wy=(side==0)?wGp:wGm;   /* +1 then -1 */
        for (int a=0;a<2;a++,n++){
            rp[n]=wy[0]*c[RHO(a,0)]+wy[1]*c[RHO(a,1)];
            Ep[n]=wy[0]*c[EN(a,0)] +wy[1]*c[EN(a,1)];
            mxp[n]=wy[0]*c[MX(a,0)]+wy[1]*c[MX(a,1)];
            myp[n]=wy[0]*c[MY(a,0)]+wy[1]*c[MY(a,1)];
        }
    }
}

__global__ void
pos_limit(real * __restrict__ Qp, int N, int Np)
{
    const real eps_r=1e-12, eps_p=1e-12;
    int N2=N*N;
    for (int k=blockIdx.x*blockDim.x+threadIdx.x; k<N2; k+=blockDim.x*gridDim.x){
        int j=k/N, i=k%N, jp=j+G_GHOST, ip=i+G_GHOST;
        real c[NDOF]; load_cell(Qp,jp,ip,Np,c);

        real rb=0.,mxb=0.,myb=0.,Eb=0.;
        for(int q=0;q<4;q++){ rb+=c[0+q]; Eb+=c[4+q]; mxb+=c[8+q]; myb+=c[12+q]; }
        rb*=0.25; mxb*=0.25; myb*=0.25; Eb*=0.25;
        if (rb<eps_r) rb=eps_r;

        real rp[NPT],mxp[NPT],myp[NPT],Ep[NPT];
        eval_points(c,rp,mxp,myp,Ep);

        /* (1) density: largest theta keeping rho>=eps_r at all points */
        real thd=1.;
        for(int n=0;n<NPT;n++) if (rp[n]<eps_r){
            real th=(rb-eps_r)/(rb-rp[n]);
            thd=fmin(thd,fmin(fmax(th,0.),1.));
        }
        if (thd<1.) {
            for(int q=0;q<4;q++) c[0+q]=rb+thd*(c[0+q]-rb);
            for(int n=0;n<NPT;n++) rp[n]=rb+thd*(rp[n]-rb);
        }

        /* (2) pressure: bisection for theta keeping p>=eps_p at all points */
        #define MINP(TH) ({ real _pm=1e30; \
            for(int n=0;n<NPT;n++){ \
                real _r=rb +(TH)*(rp[n]-rb);  _r=fmax(_r,1e-14); \
                real _e=Eb +(TH)*(Ep[n]-Eb); \
                real _mx=mxb+(TH)*(mxp[n]-mxb); \
                real _my=myb+(TH)*(myp[n]-myb); \
                real _p=(GAMMA_V-1.)*(_e-0.5*(_mx*_mx+_my*_my)/_r); \
                _pm=fmin(_pm,_p);} _pm; })
        if (MINP(1.0) < eps_p) {
            real lo=0., hi=1.;
            for (int it=0; it<30; it++){ real mid=0.5*(lo+hi);
                if (MINP(mid) >= eps_p) lo=mid; else hi=mid; }
            real th=lo;
            for(int q=0;q<4;q++){
                c[0+q] =rb +th*(c[0+q]-rb);   c[4+q] =Eb +th*(c[4+q]-Eb);
                c[8+q] =mxb+th*(c[8+q]-mxb);  c[12+q]=myb+th*(c[12+q]-myb);
            }
        }
        #undef MINP
        for(int q=0;q<NDOF;q++) Qp[IDX_P(q,jp,ip,Np)]=c[q];
    }
}

/* ═══════════════════════════════════════════════════════════════════════════
 * Initial condition — circular Sod shock tube on [0,1]^2  (zero velocity)
 *   inside  r<0.25 : rho=1.0,   p=1.0
 *   outside        : rho=0.125, p=0.1
 * ═══════════════════════════════════════════════════════════════════════════ */
__device__ __forceinline__ void
sod_state(real x, real y, real *rho, real *p)
{
    real r=sqrt((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5));
    if (r<0.25) { *rho=10.0;   *p=10.0; }
    else        { *rho=0.125; *p=0.1; }
}

__global__ void
ic_kernel(real * __restrict__ Qp, int N, int Np, real h)
{
    /* 1D node coords on [-1,1]: Gauss +-g, Lobatto +-1, mapped to [0,1] */
    const real g=0.5773502691896257;
    real xG[2]={(1.-g)*0.5,(1.+g)*0.5};   /* momentum DOFs are 0 at IC */
    int N2=N*N;
    for (int k=blockIdx.x*blockDim.x+threadIdx.x; k<N2; k+=blockDim.x*gridDim.x){
        int j=k/N, i=k%N, jp=j+G_GHOST, ip=i+G_GHOST;
        real x0=i*h, y0=j*h, rho,p;
        /* rho,E on Gauss x Gauss */
        for(int b=0;b<2;b++) for(int a=0;a<2;a++){
            sod_state(x0+xG[a]*h, y0+xG[b]*h, &rho,&p);
            Qp[IDX_P(RHO(a,b),jp,ip,Np)]=rho;
            Qp[IDX_P(EN (a,b),jp,ip,Np)]=p/(GAMMA_V-1.);   /* zero velocity */
        }
        /* momentum DOFs : zero velocity */
        for(int b=0;b<2;b++) for(int a=0;a<2;a++){
            Qp[IDX_P(MX(a,b),jp,ip,Np)]=0.;
            Qp[IDX_P(MY(a,b),jp,ip,Np)]=0.;
        }
    }
}

/* ═══════════════════════════════════════════════════════════════════════════
 * Low-Mach vortex (lmv) — Barsukow et al., Sec 2.  Stationary Gresho-type
 * vortex on [0,1]^2, periodic.  rho=1; swirl v_phi(r); pressure in radial
 * equilibrium.  p0 = 1/(gamma eps^2) - 1/2  fixes the peak local Mach = eps.
 * The exact solution is steady, so the L2 error measures how well the scheme
 * preserves a low-Mach steady state.
 * ═══════════════════════════════════════════════════════════════════════════ */
__device__ __forceinline__ real lmv_vphi(real r){
    if (r<0.2) return 5.*r;
    if (r<0.4) return 2.-5.*r;
    return 0.;
}
__device__ __forceinline__ real lmv_pressure(real r, real p0){
    if (r<0.2) return p0 + 12.5*r*r;
    if (r<0.4) return p0 + 4.*log(5.*r) + 4. - 20.*r + 12.5*r*r;
    return p0 + 4.*log(2.) - 2.;
}
__device__ __forceinline__ void
lmv_exact(real x, real y, real p0, real *rho, real *u, real *v, real *p){
    real dx=x-0.5, dy=y-0.5, r=sqrt(dx*dx+dy*dy);
    real vp=lmv_vphi(r), ir=(r>1e-30)?1./r:0.;
    *rho=1.; *u=-vp*dy*ir; *v=vp*dx*ir; *p=lmv_pressure(r,p0);
}

/* Project the exact vortex onto the 16 nodal DOFs (interpolation at nodes). */
__global__ void
ic_lmv(real * __restrict__ Qp, int N, int Np, real h, real eps)
{
    const real p0=1./(GAMMA_V*eps*eps)-0.5;
    const real g=0.5773502691896257;
    real xG[2]={(1.-g)*0.5,(1.+g)*0.5};
    int N2=N*N;
    for (int k=blockIdx.x*blockDim.x+threadIdx.x; k<N2; k+=blockDim.x*gridDim.x){
        int j=k/N, i=k%N, jp=j+G_GHOST, ip=i+G_GHOST;
        real x0=i*h, y0=j*h, rho,u,v,p;
        for(int b=0;b<2;b++) for(int a=0;a<2;a++){
            /* all fields collocated on the Gauss x Gauss grid */
            lmv_exact(x0+xG[a]*h, y0+xG[b]*h, p0, &rho,&u,&v,&p);
            Qp[IDX_P(RHO(a,b),jp,ip,Np)]=rho;
            Qp[IDX_P(EN (a,b),jp,ip,Np)]=p/(GAMMA_V-1.)+0.5*rho*(u*u+v*v);
            Qp[IDX_P(MX (a,b),jp,ip,Np)]=rho*u;
            Qp[IDX_P(MY (a,b),jp,ip,Np)]=rho*v;
        }
    }
}

/* L2 error vs. exact (steady) vortex, by 2x2 Gauss quadrature per cell. */
__global__ void
compute_l2_err_lmv(const real * __restrict__ Qp, real * __restrict__ err_rho,
                   real * __restrict__ err_p, int N, int Np, real h, real eps)
{
    const real p0=1./(GAMMA_V*eps*eps)-0.5;
    const real g=0.5773502691896257;
    real xG[2]={(1.-g)*0.5,(1.+g)*0.5};
    int N2=N*N;
    for (int k=blockIdx.x*blockDim.x+threadIdx.x; k<N2; k+=blockDim.x*gridDim.x){
        int j=k/N, i=k%N, jp=j+G_GHOST, ip=i+G_GHOST;
        real c[NDOF]; load_cell(Qp,jp,ip,Np,c);
        real x0=i*h, y0=j*h, sr=0., sp=0.;
        for(int b=0;b<2;b++) for(int a=0;a<2;a++){
            real rho=c[RHO(a,b)], En=c[EN(a,b)];
            real mx=c[MX(a,b)], my=c[MY(a,b)];
            real ph=(GAMMA_V-1.)*(En-0.5*(mx*mx+my*my)/fmax(rho,1e-14));
            real re,ue,ve,pe;
            lmv_exact(x0+xG[a]*h, y0+xG[b]*h, p0, &re,&ue,&ve,&pe);
            sr += (rho-re)*(rho-re);
            sp += (ph-pe)*(ph-pe);
        }
        real w=0.25*h*h;   /* (h/2)^2 * (Gauss weight 1) */
        err_rho[k]=w*sr; err_p[k]=w*sp;
    }
}

/* ═══════════════════════════════════════════════════════════════════════════
 * Doubly periodic shear layer  (Bell, Colella & Glaz 1989)
 *
 * Domain [0,1]^2, periodic everywhere.  Two shear layers at y=0.25 and y=0.75
 * roll up into counter-rotating billows — the classic low-Mach benchmark.
 *   background : rho=gamma, p=1  ->  c=1, so |velocity| amplitude = Mach number
 *   u = Ma*tanh((y-0.25)/d_s)   y<=0.5 ;  Ma*tanh((0.75-y)/d_s)  y>0.5
 *   v = Ma*delta*sin(2 pi x)
 * shear thickness d_s = 1/30, transverse perturbation delta = 0.05.
 *
 * Usage: ./dg1 N dsl [Ma]    default Ma = 0.1
 * ═══════════════════════════════════════════════════════════════════════════ */
#define DSL_THICK 30.0    /* 1/d_s : shear-layer sharpness            */
#define DSL_PERT  0.05    /* delta : transverse perturbation amplitude */
__device__ __forceinline__ void
dsl_exact(real x, real y, real Ma, real *rho, real *u, real *v, real *p){
    real us = (y<=0.5) ? tanh(DSL_THICK*(y-0.25)) : tanh(DSL_THICK*(0.75-y));
    *rho=GAMMA_V; *u=Ma*us; *v=Ma*DSL_PERT*sin(2.*M_PI*x); *p=1.;
}
__global__ void
ic_dsl(real * __restrict__ Qp, int N, int Np, real h, real Ma){
    const real g=0.5773502691896257;
    real xG[2]={(1.-g)*0.5,(1.+g)*0.5};
    int N2=N*N;
    for (int k=blockIdx.x*blockDim.x+threadIdx.x; k<N2; k+=blockDim.x*gridDim.x){
        int j=k/N, i=k%N, jp=j+G_GHOST, ip=i+G_GHOST;
        real x0=i*h, y0=j*h, rho,u,v,p;
        for(int b=0;b<2;b++) for(int a=0;a<2;a++){
            dsl_exact(x0+xG[a]*h, y0+xG[b]*h, Ma, &rho,&u,&v,&p);
            Qp[IDX_P(RHO(a,b),jp,ip,Np)]=rho;
            Qp[IDX_P(EN (a,b),jp,ip,Np)]=p/(GAMMA_V-1.)+0.5*rho*(u*u+v*v);
            Qp[IDX_P(MX (a,b),jp,ip,Np)]=rho*u;
            Qp[IDX_P(MY (a,b),jp,ip,Np)]=rho*v;
        }
    }
}

/* ═══════════════════════════════════════════════════════════════════════════
 * Vortex-acoustic wave interaction  (Coiffier 2025, Sec 6.6.3)
 *
 * A compact-support stationary isentropic vortex superposed with a right-
 * traveling low-Mach acoustic pulse, on [0,1]^2 periodic.  Background
 * rho=1, p=1, c=sqrt(gamma).  Tests how the acoustic pulse passes through /
 * interacts with the slow vortical mode at low Mach.
 *
 * Usage: ./dg1 N va [Mref]    default Mref = 0.1
 * ═══════════════════════════════════════════════════════════════════════════ */
/* Compact-support isentropic vortex perturbation (centre (0.5,0.5), Rv=0.2). */
__device__ __forceinline__ void
vacwav_vortex(real x, real y, real Mref, real *drho, real *du, real *dv, real *dp){
    const real gm=GAMMA_V, c0=sqrt(gm), Rv=0.2;
    real dx=x-0.5, dy=y-0.5, r=sqrt(dx*dx+dy*dy), rbar=r/Rv;
    if (rbar>=1.){ *drho=0.; *du=0.; *dv=0.; *dp=0.; return; }
    /* C-inf bump f(s)=s*exp(1/(s^2-1)); normalise so peak u_theta = Mref*c0 */
    const real s_star=0.5*(sqrt(6.)-sqrt(2.));
    const real fmax=s_star*exp(1./(s_star*s_star-1.));
    real bump=rbar*exp(1./(rbar*rbar-1.));
    real uth=Mref*c0*bump/fmax;
    real c2=c0*c0-0.5*(gm-1.)*uth*uth;          /* isentropic Bernoulli balance */
    if (c2<1e-6*c0*c0) c2=1e-6*c0*c0;
    real ratio=c2/(c0*c0), ir=1./(r+1e-30);
    *du  =-uth*dy*ir;  *dv = uth*dx*ir;
    *drho=pow(ratio,1./(gm-1.))-1.;
    *dp  =pow(ratio,gm /(gm-1.))-1.;
}
/* Right-traveling acoustic pulse perturbation (centre xw=0.10, half-width 0.05). */
__device__ __forceinline__ void
vacwav_wave(real x, real Mref, real *drho, real *du, real *dp){
    const real gm=GAMMA_V, c0=sqrt(gm), xw=0.10, sw=0.05;
    real xbar=(x-xw)/sw;
    if (fabs(xbar)>=1.){ *drho=0.; *du=0.; *dp=0.; return; }
    real bump=exp(1./(xbar*xbar-1.));
    *du  =Mref*c0*bump;   *drho=Mref*bump;   *dp=gm*Mref*bump;
}
__global__ void
ic_vacwav(real * __restrict__ Qp, int N, int Np, real h, real Mref){
    const real g=0.5773502691896257;
    real xG[2]={(1.-g)*0.5,(1.+g)*0.5};
    int N2=N*N;
    for (int k=blockIdx.x*blockDim.x+threadIdx.x; k<N2; k+=blockDim.x*gridDim.x){
        int j=k/N, i=k%N, jp=j+G_GHOST, ip=i+G_GHOST;
        real x0=i*h, y0=j*h;
        for(int b=0;b<2;b++) for(int a=0;a<2;a++){
            real x=x0+xG[a]*h, y=y0+xG[b]*h;
            real drv,duv,dvv,dpv; vacwav_vortex(x,y,Mref,&drv,&duv,&dvv,&dpv);
            real drw,duw,dpw;     vacwav_wave  (x,  Mref,&drw,&duw,&dpw);
            real rho=fmax(1.+drv+drw,1e-14);
            real u=duv+duw, v=dvv, p=fmax(1.+dpv+dpw,1e-14);
            Qp[IDX_P(RHO(a,b),jp,ip,Np)]=rho;
            Qp[IDX_P(EN (a,b),jp,ip,Np)]=p/(GAMMA_V-1.)+0.5*rho*(u*u+v*v);
            Qp[IDX_P(MX (a,b),jp,ip,Np)]=rho*u;
            Qp[IDX_P(MY (a,b),jp,ip,Np)]=rho*v;
        }
    }
}

/* ═══════════════════════════════════════════════════════════════════════════
 * SSP-RK3 stage kernels (operate on all 16 DOFs of interior cells)
 * ═══════════════════════════════════════════════════════════════════════════ */
__global__ void rk3_s1(real *U1,const real *U0,const real *L,real dt,int N,int Np){
    int N2=N*N, Np2=Np*Np;
    for(int k=blockIdx.x*blockDim.x+threadIdx.x;k<N2;k+=blockDim.x*gridDim.x){
        int j=k/N,i=k%N,idp=(j+G_GHOST)*Np+(i+G_GHOST);
        for(int q=0;q<NDOF;q++) U1[q*Np2+idp]=U0[q*Np2+idp]+dt*L[q*N2+k];
    }
}
__global__ void rk3_s2(real *U2,const real *U0,const real *U1,const real *L,
                       real dt,int N,int Np){
    int N2=N*N, Np2=Np*Np;
    for(int k=blockIdx.x*blockDim.x+threadIdx.x;k<N2;k+=blockDim.x*gridDim.x){
        int j=k/N,i=k%N,idp=(j+G_GHOST)*Np+(i+G_GHOST);
        for(int q=0;q<NDOF;q++)
            U2[q*Np2+idp]=0.75*U0[q*Np2+idp]+0.25*(U1[q*Np2+idp]+dt*L[q*N2+k]);
    }
}
__global__ void rk3_s3(real *U,const real *U0,const real *U2,const real *L,
                       real dt,int N,int Np){
    int N2=N*N, Np2=Np*Np;
    for(int k=blockIdx.x*blockDim.x+threadIdx.x;k<N2;k+=blockDim.x*gridDim.x){
        int j=k/N,i=k%N,idp=(j+G_GHOST)*Np+(i+G_GHOST);
        for(int q=0;q<NDOF;q++)
            U[q*Np2+idp]=(1./3.)*U0[q*Np2+idp]+(2./3.)*(U2[q*Np2+idp]+dt*L[q*N2+k]);
    }
}

/* ═══════════════════════════════════════════════════════════════════════════
 * reductions
 * ═══════════════════════════════════════════════════════════════════════════ */
__global__ void reduce_max(const real *in,real *out,int n){
    extern __shared__ real sm[]; int tid=threadIdx.x; real v=0.;
    for(int k=blockIdx.x*blockDim.x+tid;k<n;k+=blockDim.x*gridDim.x) v=fmax(v,in[k]);
    sm[tid]=v; __syncthreads();
    for(int s=BLOCK1D/2;s>0;s>>=1){ if(tid<s) sm[tid]=fmax(sm[tid],sm[tid+s]); __syncthreads(); }
    if(tid==0) out[blockIdx.x]=sm[0];
}
__global__ void reduce_min(const real *in,real *out,int n){
    extern __shared__ real sm[]; int tid=threadIdx.x; real v=1e30;
    for(int k=blockIdx.x*blockDim.x+tid;k<n;k+=blockDim.x*gridDim.x) v=fmin(v,in[k]);
    sm[tid]=v; __syncthreads();
    for(int s=BLOCK1D/2;s>0;s>>=1){ if(tid<s) sm[tid]=fmin(sm[tid],sm[tid+s]); __syncthreads(); }
    if(tid==0) out[blockIdx.x]=sm[0];
}
static real gpu_max(const real *d,real *tmp,int n){
    reduce_max<<<GS_NBLK,BLOCK1D,BLOCK1D*sizeof(real)>>>(d,tmp,n);
    reduce_max<<<1,BLOCK1D,BLOCK1D*sizeof(real)>>>(tmp,tmp,GS_NBLK);
    real v; CK(cudaMemcpy(&v,tmp,sizeof(real),cudaMemcpyDeviceToHost)); return v;
}
static real gpu_min(const real *d,real *tmp,int n){
    reduce_min<<<GS_NBLK,BLOCK1D,BLOCK1D*sizeof(real)>>>(d,tmp,n);
    reduce_min<<<1,BLOCK1D,BLOCK1D*sizeof(real)>>>(tmp,tmp,GS_NBLK);
    real v; CK(cudaMemcpy(&v,tmp,sizeof(real),cudaMemcpyDeviceToHost)); return v;
}
__global__ void reduce_sum(const real *in,real *out,int n){
    extern __shared__ real sm[]; int tid=threadIdx.x; real v=0.;
    for(int k=blockIdx.x*blockDim.x+tid;k<n;k+=blockDim.x*gridDim.x) v+=in[k];
    sm[tid]=v; __syncthreads();
    for(int s=BLOCK1D/2;s>0;s>>=1){ if(tid<s) sm[tid]+=sm[tid+s]; __syncthreads(); }
    if(tid==0) out[blockIdx.x]=sm[0];
}
static real gpu_sum(const real *d,real *tmp,int n){
    reduce_sum<<<GS_NBLK,BLOCK1D,BLOCK1D*sizeof(real)>>>(d,tmp,n);
    reduce_sum<<<1,BLOCK1D,BLOCK1D*sizeof(real)>>>(tmp,tmp,GS_NBLK);
    real v; CK(cudaMemcpy(&v,tmp,sizeof(real),cudaMemcpyDeviceToHost)); return v;
}

/* ═══════════════════════════════════════════════════════════════════════════
 * Diagnostics: cell-mean density (left panel) and pressure (right panel)
 * ═══════════════════════════════════════════════════════════════════════════ */
__global__ void
extract_rho_p(const real * __restrict__ Qp, real * __restrict__ orho,
              real * __restrict__ op, int N, int Np)
{
    int N2=N*N;
    for(int k=blockIdx.x*blockDim.x+threadIdx.x;k<N2;k+=blockDim.x*gridDim.x){
        int j=k/N,i=k%N,jp=j+G_GHOST,ip=i+G_GHOST;
        real rb=0,mxb=0,myb=0,Eb=0;
        for(int q=0;q<4;q++){
            rb +=Qp[IDX_P(0+q,jp,ip,Np)]; Eb +=Qp[IDX_P(4+q,jp,ip,Np)];
            mxb+=Qp[IDX_P(8+q,jp,ip,Np)]; myb+=Qp[IDX_P(12+q,jp,ip,Np)];
        }
        rb*=0.25; mxb*=0.25; myb*=0.25; Eb*=0.25;
        rb=fmax(rb,1e-14);
        orho[k]=rb;
        op[k]=(GAMMA_V-1.)*(Eb-0.5*(mxb*mxb+myb*myb)/rb);
    }
}

/* |velocity| (left panel) and pressure (right panel), from the cell mean. */
__global__ void
extract_magv_p(const real * __restrict__ Qp, real * __restrict__ omv,
               real * __restrict__ op, int N, int Np)
{
    int N2=N*N;
    for(int k=blockIdx.x*blockDim.x+threadIdx.x;k<N2;k+=blockDim.x*gridDim.x){
        int j=k/N,i=k%N,jp=j+G_GHOST,ip=i+G_GHOST;
        real rb=0,mxb=0,myb=0,Eb=0;
        for(int q=0;q<4;q++){
            rb +=Qp[IDX_P(0+q,jp,ip,Np)]; Eb +=Qp[IDX_P(4+q,jp,ip,Np)];
            mxb+=Qp[IDX_P(8+q,jp,ip,Np)]; myb+=Qp[IDX_P(12+q,jp,ip,Np)];
        }
        rb*=0.25; mxb*=0.25; myb*=0.25; Eb*=0.25; rb=fmax(rb,1e-14);
        real u=mxb/rb, v=myb/rb;
        omv[k]=sqrt(u*u+v*v);
        op[k]=(GAMMA_V-1.)*(Eb-0.5*(mxb*mxb+myb*myb)/rb);
    }
}

/* |velocity| (left panel) and vorticity w = dv/dx - du/dy (right panel).
 * The cell velocity is the Q1 bilinear interpolant of the 4 Gauss-node values;
 * its derivatives use the Gauss-basis derivative DGv scaled by 2/h.  The value
 * stored is the cell mean of w (element-wise, so discontinuous across cells). */
__global__ void
extract_magv_vort(const real * __restrict__ Qp, real * __restrict__ omv,
                  real * __restrict__ ovort, int N, int Np, real h)
{
    int N2=N*N;
    for(int k=blockIdx.x*blockDim.x+threadIdx.x;k<N2;k+=blockDim.x*gridDim.x){
        int j=k/N,i=k%N,jp=j+G_GHOST,ip=i+G_GHOST;
        real c[NDOF]; load_cell(Qp,jp,ip,Np,c);
        real u[2][2], v[2][2];
        for(int b=0;b<2;b++) for(int a=0;a<2;a++){
            real r=fmax(c[RHO(a,b)],1e-14);
            u[a][b]=c[MX(a,b)]/r;  v[a][b]=c[MY(a,b)]/r;
        }
        real ub=0.25*(u[0][0]+u[1][0]+u[0][1]+u[1][1]);
        real vb=0.25*(v[0][0]+v[1][0]+v[0][1]+v[1][1]);
        omv[k]=sqrt(ub*ub+vb*vb);
        /* cell-mean dv/dx and du/dy of the bilinear field */
        real dvdx=(2./h)*(DGv[0]*0.5*(v[0][0]+v[0][1])+DGv[1]*0.5*(v[1][0]+v[1][1]));
        real dudy=(2./h)*(DGv[0]*0.5*(u[0][0]+u[1][0])+DGv[1]*0.5*(u[0][1]+u[1][1]));
        ovort[k]=dvdx-dudy;
    }
}

/* ── compact analytic colormap (viridis-like, 5 stops) ─────────────────────── */
static void cmap(double t, unsigned char rgb[3]){
    static const double s[6][4]={
        {0.00, 68,  1, 84},{0.25, 59, 82,139},{0.50, 33,145,140},
        {0.75, 94,201, 98},{1.00,253,231, 37},{2.00,253,231, 37}};
    if(t<0)t=0; if(t>1)t=1;
    int i=0; while(t>s[i+1][0]) i++;
    double f=(t-s[i][0])/(s[i+1][0]-s[i][0]+1e-30);
    for(int c=0;c<3;c++) rgb[c]=(unsigned char)(s[i][c+1]+f*(s[i+1][c+1]-s[i][c+1]));
}

static void
write_png_2panel(const char *fn,const real *tl,real tlmin,real tlmax,
                 const real *tr,real trmin,real trmax,int N)
{
    int W=2*N,H=N;
    unsigned char *px=(unsigned char*)malloc((size_t)W*H*3);
    if(!px){ fprintf(stderr,"OOM\n"); return; }
    for(int r=0;r<H;r++){ int pj=N-1-r;
        for(int col=0;col<W;col++){
            const real *pan=(col<N)?tl:tr;
            real vmin=(col<N)?tlmin:trmin, vmax=(col<N)?tlmax:trmax;
            int ci=(col<N)?col:col-N;
            real v=(pan[pj*N+ci]-vmin)/(vmax-vmin+1e-30);
            unsigned char c3[3]; cmap(v,c3);
            unsigned char *p=px+(r*W+col)*3; p[0]=c3[0];p[1]=c3[1];p[2]=c3[2];
        }
    }
    stbi_write_png(fn,W,H,3,px,W*3); free(px);
    printf("  saved %s\n",fn);
}

/* ═══════════════════════════════════════════════════════════════════════════
 * main
 * ═══════════════════════════════════════════════════════════════════════════ */
int main(int argc,char**argv)
{
    int  N    = (argc>1)?atoi(argv[1]):200;
    const char *mode = (argc>2)?argv[2]:"";
    int  do_lmv = (strcmp(mode,"lmv")==0);
    int  do_dsl = (strcmp(mode,"dsl")==0);   /* doubly periodic shear layer */
    int  do_va  = (strcmp(mode,"va" )==0);   /* vortex-acoustic wave        */
    /* third arg is a Mach-like parameter: lmv eps / dsl Ma / va Mref */
    real param = 0.1;
    if ((do_lmv||do_dsl||do_va) && argc>3) { char*e; real m=strtod(argv[3],&e);
        if (e!=argv[3] && *e=='\0' && m>0.) param=m; }
    real eps = param;                         /* alias used by the lmv paths */
    /* Riemann flux is a compile-time setting (RIEMANN macro); see top of file. */

    real L    = 1.0;
    real h    = L/N;
    int  Np   = N+2*G_GHOST;
    real CFL  = 0.20;            /* Q1 DG: ~1/(2k+1) of FV limit */
    real t_end;
    if      (do_lmv) t_end = 2.*M_PI*0.2;     /* one vortex rotation period   */
    else if (do_dsl) t_end = 2.2/param;       /* rollup time, scaled by 1/Ma  */
    else if (do_va)  t_end = 3.5;             /* matches Coiffier Fig. 6.6    */
    else             t_end = 1.00;            /* circular Sod                 */

    char prefix[32];
    if      (do_lmv) snprintf(prefix,sizeof(prefix),"rt_dg1_lmv%04d",(int)round(param*1e4));
    else if (do_dsl) snprintf(prefix,sizeof(prefix),"rt_dg1_dsl%04d",(int)round(param*1e4));
    else if (do_va)  snprintf(prefix,sizeof(prefix),"rt_dg1_va%04d", (int)round(param*1e4));
    else             snprintf(prefix,sizeof(prefix),"rt_dg1_sod");

    int dev; cudaDeviceProp prop; CK(cudaGetDevice(&dev));
    CK(cudaGetDeviceProperties(&prop,dev));
    printf("========================================================\n");
    printf("  Device : %s\n",prop.name);
    printf("  Degree-1 Gauss-collocation DG (16 DOF/cell) + %s + Zhang-Shu limiter\n",
           RIEMANN_NAME);
    printf("  Face reconstruction: %s\n", RECON_NAME);
    if      (do_lmv)
        printf("  Low-Mach vortex  eps=%.4f  N=%dx%d  h=%.5f  CFL=%.2f  tf=%.3f\n",
               param,N,N,h,CFL,t_end);
    else if (do_dsl)
        printf("  Doubly periodic shear layer  Ma=%.4f  N=%dx%d  h=%.5f  CFL=%.2f  tf=%.3f\n",
               param,N,N,h,CFL,t_end);
    else if (do_va)
        printf("  Vortex-acoustic wave  Mref=%.4f  N=%dx%d  h=%.5f  CFL=%.2f  tf=%.3f\n",
               param,N,N,h,CFL,t_end);
    else
        printf("  Circular Sod shock tube  N=%dx%d  h=%.5f  CFL=%.2f  tf=%.3f\n",
               N,N,h,CFL,t_end);
    printf("========================================================\n");

    size_t szc = (size_t)N*N*sizeof(real);
    size_t szp = (size_t)NDOF*Np*Np*sizeof(real);
    size_t szL = (size_t)NDOF*N*N*sizeof(real);

    real *d_U,*d_U0,*d_U1,*d_U2,*d_L,*d_lam,*d_tmp,*d_rho,*d_p;
    CK(cudaMalloc(&d_U ,szp)); CK(cudaMalloc(&d_U0,szp));
    CK(cudaMalloc(&d_U1,szp)); CK(cudaMalloc(&d_U2,szp));
    CK(cudaMalloc(&d_L ,szL)); CK(cudaMalloc(&d_lam,szc));
    CK(cudaMalloc(&d_tmp,GS_NBLK*sizeof(real)));
    CK(cudaMalloc(&d_rho,szc)); CK(cudaMalloc(&d_p,szc));
    CK(cudaMemset(d_U,0,szp));

    real *h_rho=(real*)malloc(szc), *h_p=(real*)malloc(szc);

    if      (do_lmv) ic_lmv   <<<GS_NBLK,BLOCK1D>>>(d_U,N,Np,h,param);
    else if (do_dsl) ic_dsl   <<<GS_NBLK,BLOCK1D>>>(d_U,N,Np,h,param);
    else if (do_va ) ic_vacwav<<<GS_NBLK,BLOCK1D>>>(d_U,N,Np,h,param);
    else             ic_kernel<<<GS_NBLK,BLOCK1D>>>(d_U,N,Np,h);
    pos_limit<<<GS_NBLK,BLOCK1D>>>(d_U,N,Np);
    apply_bc <<<GS_NBLK,BLOCK1D>>>(d_U,N,Np);
    CK(cudaDeviceSynchronize());

    /* velocity+pressure panels for the smooth low-Mach tests, density+pressure for Sod */
    int show_magv = (do_lmv||do_dsl||do_va);
#define WRITE_FRAME(idx) do { \
    if (show_magv) extract_magv_vort<<<GS_NBLK,BLOCK1D>>>(d_U,d_rho,d_p,N,Np,h); \
    else           extract_rho_p    <<<GS_NBLK,BLOCK1D>>>(d_U,d_rho,d_p,N,Np); \
    CK(cudaDeviceSynchronize()); \
    real _rl=gpu_min(d_rho,d_tmp,N*N), _rh=gpu_max(d_rho,d_tmp,N*N); \
    real _pl=gpu_min(d_p,  d_tmp,N*N), _ph=gpu_max(d_p,  d_tmp,N*N); \
    CK(cudaMemcpy(h_rho,d_rho,szc,cudaMemcpyDeviceToHost)); \
    CK(cudaMemcpy(h_p,  d_p,  szc,cudaMemcpyDeviceToHost)); \
    char _fn[80]; sprintf(_fn,"figures/%s_%04d.png",prefix,(idx)); \
    write_png_2panel(_fn,h_rho,_rl,_rh,h_p,_pl,_ph,N); \
} while(0)

    WRITE_FRAME(0);

    CK(cudaMemset(d_L,0,szL));
    compute_rhs<<<GS_NBLK,BLOCK1D>>>(d_U,d_L,d_lam,N,Np,h);
    CK(cudaDeviceSynchronize());
    real lam=gpu_max(d_lam,d_tmp,N*N);
    {   /* well-balancing check: residual of the operator on the exact IC */
        const char *nm[4]={"rho","E  ","mx ","my "};
        for(int grp=0;grp<4;grp++){ real mx=0;
            for(int q=0;q<4;q++){
                mx=fmax(mx,fabs(gpu_min(d_L+(grp*4+q)*N*N,d_tmp,N*N)));
                mx=fmax(mx,fabs(gpu_max(d_L+(grp*4+q)*N*N,d_tmp,N*N))); }
            printf("  [IC residual] max|d(%s)/dt| = %.4e\n",nm[grp],mx); }
    }

    const int NF=10; int frame=1; real tnext=t_end/NF, t=0.; int step=0;
    struct timespec ts0,ts1; clock_gettime(CLOCK_MONOTONIC,&ts0);

    while (t<t_end) {
        if (!(lam>0.)){ fprintf(stderr,"bad lam at step %d t=%.5f\n",step,t); break; }
        real dt=CFL*h/lam;
        real tt=(tnext<t_end)?tnext:t_end;
        if (t+dt>tt) dt=tt-t;
        if (dt<1e-14){ fprintf(stderr,"dt underflow\n"); break; }

        CK(cudaMemcpy(d_U0,d_U,szp,cudaMemcpyDeviceToDevice));

        rk3_s1<<<GS_NBLK,BLOCK1D>>>(d_U1,d_U0,d_L,dt,N,Np);
        pos_limit<<<GS_NBLK,BLOCK1D>>>(d_U1,N,Np);
        apply_bc <<<GS_NBLK,BLOCK1D>>>(d_U1,N,Np);
        CK(cudaMemset(d_L,0,szL));
        compute_rhs<<<GS_NBLK,BLOCK1D>>>(d_U1,d_L,d_lam,N,Np,h);

        rk3_s2<<<GS_NBLK,BLOCK1D>>>(d_U2,d_U0,d_U1,d_L,dt,N,Np);
        pos_limit<<<GS_NBLK,BLOCK1D>>>(d_U2,N,Np);
        apply_bc <<<GS_NBLK,BLOCK1D>>>(d_U2,N,Np);
        CK(cudaMemset(d_L,0,szL));
        compute_rhs<<<GS_NBLK,BLOCK1D>>>(d_U2,d_L,d_lam,N,Np,h);

        rk3_s3<<<GS_NBLK,BLOCK1D>>>(d_U,d_U0,d_U2,d_L,dt,N,Np);
        pos_limit<<<GS_NBLK,BLOCK1D>>>(d_U,N,Np);
        apply_bc <<<GS_NBLK,BLOCK1D>>>(d_U,N,Np);
        CK(cudaMemset(d_L,0,szL));
        compute_rhs<<<GS_NBLK,BLOCK1D>>>(d_U,d_L,d_lam,N,Np,h);
        lam=gpu_max(d_lam,d_tmp,N*N);
        t+=dt; step++;
        if (t>=tnext-1e-12 && frame<=NF) {
            clock_gettime(CLOCK_MONOTONIC,&ts1);
            real el=(ts1.tv_sec-ts0.tv_sec)+(ts1.tv_nsec-ts0.tv_nsec)*1e-9;
            printf("  frame %2d/%d  step %5d  t=%.5f  lam=%.3e  elapsed=%.1fs\n",
                   frame,NF,step,t,lam,el); fflush(stdout);
            WRITE_FRAME(frame); frame++; tnext=frame*t_end/NF;
        }
    }
    CK(cudaDeviceSynchronize());
    clock_gettime(CLOCK_MONOTONIC,&ts1);
    real wall=(ts1.tv_sec-ts0.tv_sec)+(ts1.tv_nsec-ts0.tv_nsec)*1e-9;
    printf("  done: %d steps  t=%.6f  wall=%.2fs\n",step,t,wall);
    printf("  output: figures/%s_0000.png .. _%04d.png\n",prefix,NF);
    if (do_lmv) {
        compute_l2_err_lmv<<<GS_NBLK,BLOCK1D>>>(d_U,d_rho,d_p,N,Np,h,eps);
        CK(cudaDeviceSynchronize());
        real er=sqrt(gpu_sum(d_rho,d_tmp,N*N));
        real ep=sqrt(gpu_sum(d_p,  d_tmp,N*N));
        printf("  ── FINAL L2 ERRORS  [low-Mach vortex, steady] ──────\n");
        printf("  eps=%.4f  N=%d  h=%.5f  t=%.4f\n",eps,N,h,t);
        printf("  L2(rho) = %.6e\n",er);
        printf("  L2(p)   = %.6e\n",ep);
        printf("  left panel: |velocity|   right panel: vorticity\n");
    } else if (do_dsl || do_va) {
        printf("  left panel: |velocity|   right panel: vorticity\n");
    } else {
        printf("  left panel: density   right panel: pressure\n");
    }

    free(h_rho); free(h_p);
    cudaFree(d_U);cudaFree(d_U0);cudaFree(d_U1);cudaFree(d_U2);
    cudaFree(d_L);cudaFree(d_lam);cudaFree(d_tmp);cudaFree(d_rho);cudaFree(d_p);
    return 0;
}
