Hybrid Entropy Stable HLL-Type Riemann Solvers for Hyperbolic
Conservation Laws
Birte Schmidtmanna,∗, Andrew R. Wintersb
aMathCCES, RWTH Aachen University, Schinkelstr. 2, 52062 Aachen
bMathematisches Institut, Universität zu Köln, Weyertal 86-90, 50931 Köln
Abstract
It is known that HLL-type schemes are more dissipative than schemes based on characteristic decompositions.
However, HLL-type methods offer greater flexibility to large systems of hyperbolic conservation laws because
the eigenstructure of the flux Jacobian is not needed. We demonstrate in the present work that several
HLL-type Riemann solvers are provably entropy stable. Further, we provide convex combinations of standard
dissipation terms to create hybrid HLL-type methods that have less dissipation while retaining entropy
stability. The decrease in dissipation is demonstrated for the ideal MHD equations with a numerical example.
Keywords: entropy stability, ideal magnetohydrodynamics, HLL, Riemann solver, discrete entropy
inequality
1. Introduction
We consider the numerical solution of systems of hyperbolic conservation laws of the form
∂q
+ f =0, (1.1)
∂t ∇·
on a domain Ω. For a one-dimensional approximation we divide Ω into K non-overlapping grid cells
C = [x , x ], i = 1,...,K which are not necessarily equidistant. In the context of finite volume
i i−1/2 i+1/2
schemes, hyperbolic equations, such as (1.1), require a numerical flux function which fully determines the
properties of the scheme [5]. The numerical flux function takes as input the left and right value of q at
the cell interface and solves a local Riemann problem. Smooth initial flows governed by (1.1) may develop
discontinuities (e.g. shocks) in finite time. Thus, solutions are sought in the weak sense [5]. However, weak
solutions are not unique and need to be supplemented with additional admissibility criteria. Following the
work of e.g. [8, 10], we use the concept of entropy to construct discretizations that agree with the second law
of thermodynamics. That is, the numerical flux function will possess entropy stability, cf. [8] and references
therein.
In particular, we prove entropy stability for the HLL scheme and present the construction of HLL-type
entropy stable numerical flux functions. It is known that HLL-type schemes are more dissipative than
upwind schemes. However, HLL-type methods need less information about the eigendecomposition of the
flux Jacobian. This is advantageous because the eigenstructure might be computationally expensive or no
analytical expression exists, especially for large systems. As such, we consider three standard dissipation
terms, namely Lax-Friedrichs (LF), HLL, and Lax-Wendroff (LW) and present two hybrid dissipation terms
introduced in [7]. We demonstrate that these five schemes are provably entropy stable.
The paper is organized as follows: Sec. 2 provides a brief background on entropy stable numerical fluxes.
In Sec. 3 we show entropy stability for the LF, HLL, and LW dissipation terms. The creation of two new
∗Correspondingauthor
Email address: schmidtmann@mathcces.rwth-aachen.de(BirteSchmidtmann)
Preprint submitted to Journal of Computational Physics October 24, 2016
6102
tcO
12
]AN.htam[
2v04260.7061:viXra

hybrid entropy stable dissipation operators is shown in Sec. 4. We demonstrate in Sec. 5 that the new hybrid
numericalflux reducethe overall dissipationin astandard finitevolumescheme. Ourconclusionsand outlook
| are | drawn   | in the | final     | section. |      |           |     |     |
| --- | ------- | ------ | --------- | -------- | ---- | --------- | --- | --- |
| 2.  | Entropy | stable | numerical |          | flux | functions |     |     |
A numerical method that recovers the local changes in entropy as predicted by the continuous entropy
conservation law is said to be entropy conservative. Entropy conservation is only valid for smooth flow
configurations. For discontinuous solutions, the entropy conservation law becomes the entropy inequality [8].
A numerical scheme is said to be entropy stable as long as the numerical approximation always obeys the
entropy inequality
|     | ∂S  | ∂F  |     |     |     |     |     |       |
| --- | --- | --- | --- | --- | --- | --- | --- | ----- |
|     |     | +   | 0,  |     |     |     |     | (2.1) |
|     | ∂t  | ∂x  | ≤   |     |     |     |     |       |
where we assume that the system of hyperbolic conservation laws is equipped with a strongly convex
mathematical entropy function, S, and a corresponding entropy flux, F, [8]. It is known that without
additional dissipation, entropy conservative numerical schemes produce high-frequency oscillations near
shocks, see e.g. [4, 10]. Thus, for the approximation to remain valid for general flow configurations we must
add a carefully designed dissipation term to ensure that (2.1) discretely holds.
To create an entropy stable (ES) numerical approximation we start with a baseline entropy conserving
(EC) numerical flux and then add a dissipation term. The resulting numerical flux at an arbitrary cell
1
| interface |     | i+  | takes the | form |     |     |     |     |
| --------- | --- | --- | --------- | ---- | --- | --- | --- | --- |
2
1
|     | f∗,ES | =f∗,EC |     | D q | ,   |     |     | (2.2) |
| --- | ----- | ------ | --- | --- | --- | --- | --- | ----- |
2
|     |     |     | −   | (cid:74) | (cid:75) |     |     |     |
| --- | --- | --- | --- | -------- | -------- | --- | --- | --- |
where q is the vector of conserved variables, D=D(q ,q ) is a suitable dissipation matrix evaluated at
i i+1
some mean state between the two cells, and =() () is the jump between the right and left cells.
i+1 i
|     |     |     |     |     |     | the(cid:74)i | · n(cid:75)dice · − · |     |
| --- | --- | --- | --- | --- | --- | ------------ | --------------------- | --- |
For simplicity of presentation we suppress s on t he numerical flux, the dissipation matrix, and any
jump terms. For an ES scheme, the baseline central flux from a classical Riemann solver is replaced by the
baseline EC flux. To guarantee entropy stability, the dissipation term in (2.2) must be carefully designed to
ensure that f∗,ES discretely satisfies the entropy inequality (2.1). To do so, we rewrite the dissipation term
[6]
|     | 1   |                           | 1   |                   |     |     |     |       |
| --- | --- | ------------------------- | --- | ----------------- | --- | --- | --- | ----- |
|     | D   | q                         | DH  | v ,               |     |     |     | (2.3) |
|     | 2   | (cid:74) (cid:75)(cid:39) | 2   | (cid:74) (cid:75) |     |     |     |       |
∂S ∂q
with the vector of entropy variables v = and the entropy Jacobian H= which relates the variables in
∂q ∂v
conserved and entropy space. Substituting (2.3) into (2.2) the entropy stable numerical flux becomes
1
|     | f∗,ES | =f∗,EC |     | DH  | v .               |     |     | (2.4) |
| --- | ----- | ------ | --- | --- | ----------------- | --- | --- | ----- |
|     |       |        | −   | 2   | (cid:74) (cid:75) |     |     |       |
The reformulation of the dissipation term, incorporating the jump in entropy variables (rather than the jump
in conservative variables) makes it possible to show entropy stability [1]. From the structure of the entropy
stable flux (2.4), we find a discrete version of the entropy inequality (2.1) in cell i to be [10]
|     | ∂S  | (cid:16) |     | (cid:17) | 1           |          | !                  |       |
| --- | --- | -------- | --- | -------- | ----------- | -------- | ------------------ | ----- |
|     |     | i + F    | F   |          |             | v T DH   | v 0.               | (2.5) |
|     |     |          | i+1 | i−1      |             |          |                    |       |
|     | ∂t  |          | 2 − | 2        | ≤−2(cid:74) | (cid:75) | (cid:74) (cid:75)≤ |       |
Thus, to guarantee discrete entropy stability, it is sufficient to show that DH is symmetric positive definite
(s.p.d).
| 3.  | Entropy | stable | classical |     | Riemann | solvers |     |     |
| --- | ------- | ------ | --------- | --- | ------- | ------- | --- | --- |
In thissection wedemonstrate entropy stabilityfor thenumerical fluxof theform (2.4) forthedissipation
matrix D of the LF, HLL, and LW scheme. To do so, we first assume that the flux Jacobian, A, or a suitable
| Roe | matrix,  | exists | with | the properties |     |     |     |       |
| --- | -------- | ------ | ---- | -------------- | --- | --- | --- | ----- |
|     | A=RΛR−1, |        |      |                |     | T   |     |       |
|     |          |        |      | H=(RZ)(RZ)     |     | ,   |     | (3.1) |
2

where R is the eigenvector matrix, Λ the diagonal corresponding eigenvalue matrix, and Z is a positive
diagonalscalingmatrixwhichcreatesasetofentropyscaledeigenvectorsRZ[1]. Weseethat,byconstruction
in (3.1), the matrix H is s.p.d. In the later proofs we only use the existence of the matrices R, Λ, and Z,
whereas, in practice, their explicit form does not need to be known. This is advantageous, because for large
systems of conservation laws, the eigendecomposition is expensive to compute or is not available.
To write the entropy stable LF scheme we substitute the dissipation matrix
∆x
D = I, (3.2)
LF
∆t
into the form (2.4). Now the complete dissipation term for (2.4) only depends of the known s.p.d matrix H,
so discrete entropy stability (2.5) for D follows immediately. We note that under the same assumption the
LF
local Lax-Friedrichs (LLF) and Roe-type dissipation terms satisfy the discrete entropy stability [10]. Next,
we consider the HLL flux. The numerical flux function of ES-HLL can be written in form (2.4), with the
dissipation matrix
λ λ λ λ λ λ
D =a I+a A, a = | L | R −| R | L , a = | R |−| L |. (3.3)
HLL 0 1 0 1
λ λ λ λ
R L R L
− −
Here, λ are the fastest signal velocities with λ <λ . >From assumption (3.1) it is straightforward to
L,R L R
show that the discrete entropy stability condition (2.5) is equivalent to showing that a +a λ 0 for all
0 1 i
≥
λ [λ ,λ ]. From the form of the coefficients a and a , keeping in mind that λ <λ and λ [λ ,λ ],
i L R 0 1 L R i L R
∈ ∈
we find
λ λ λ λ +(λ λ )λ 0 λ (λ λ )+ λ (λ λ ) 0. (3.4)
L R R L R L i L R i R i L
| | −| | | |−| | ≥ ⇔ | | − | | − ≥
Finally, we consider the dissipation matrix for LW
∆t
D = A2, (3.5)
LW
∆x
in the ES numerical flux (2.4) and find that the discrete entropy inequality (2.5) holds, i.e.,
∂S + (cid:16) F F (cid:17) 1 ∆t (cid:0) ZRT v (cid:1)T Λ2(cid:0) ZRT v (cid:1) 0. (3.6)
∂t i+ 2 1 − i−1 2 ≤−2∆x (cid:74) (cid:75) (cid:74) (cid:75) ≤
4. Hybrid entropy stable Riemann solvers
In this section, we consider dissipation matrices for hybrid Riemann solvers constructed in [7], motivated
by the work of Degond et al. [3]. The hybrid solvers were constructed using weighted combinations of the
dissipation matrices described in Sec. 3. The weighting is chosen to reduce dissipation, especially for signal
velocities close to zero. As such, the hybrid terms contain a parameter ω [0,1] which allows further control
∈
over the amount of dissipation added to the scheme.
4.1. HLLω
First we consider a generalization of the ES-HLL flux. The dissipation matrix, containing a parameter
ω [0,1], is a linear function of the eigenvalues of the associated flux Jacobian (or Roe matrix) A
∈
D =b (ω)I+b (ω)A,
HLLω 0 1
λ (ωλ2 +(1 ω)λ ) λ (ωλ2 +(1 ω)λ ) (1 ω)(λ λ )+ω(λ2 λ2)
b (ω)= R L − | L | − L R − | R | , b (ω)= − | R |−| L | R− L .
0 1
λ λ λ λ
R L R L
− −
(4.1)
We note that if ω =0 we recover the ES-HLL matrix (3.3). For D to fulfill the discrete entropy stability
HLLω
condition (2.5), it is sufficient to show that b (ω)+b (ω)λ 0 holds for all λ [λ ,λ ]. Inserting the
0 1 i i L R
≥ ∈
coefficients b (ω), b (ω), this condition rearranges to become
0 1
(cid:104) (cid:105) (cid:104) (cid:105)
ω λ2 (λ λ )+λ2 (λ λ ) +(1 ω) λ (λ λ )+ λ (λ λ ) 0, ω [0,1]. (4.2)
L R − i R i − L − | L | R − i | R | i − L ≥ ∀ ∈
3

Fig. 1. VisualizationofthedissipationmatrixDasafunctionoftheeigenvalueλ. Forthehybriddissipation
matrices we select ω =0.4.
4.2. HLLXω
Next, we treat a hybrid dissipation matrix, HLLXω, that includes the parameter ω and the quadratic
LW term. This term requires squaring the flux Jacobian (or applying the flux twice) but reduces the overall
magnitude of the dissipation, see Fig. 1. The dissipation matrix of HLLXω is a weighted combination of LF,
HLLω, and LW, see [7] for more details.
D =β (ω)D +β (ω)D +β (ω)D ,
HLLXω 0 LF 1 HLLω 2 LW
(cid:18) (cid:19)−1
(1 ω) λ λ 1 ω
L R
β (ω)=β(ω) − | | , β (ω)=1 β(ω) − +ω ,
0 (1 ω)+ω(λ + λ ) 1 − λ + λ (4.3)
L R L R
− | | | | | | | |
(cid:12) (cid:12)
λ R λ L (cid:12)λ R λ L (cid:12)
β (ω)=β(ω), β(ω)=ω+(1 ω) − − | |−| | ω+(1 ω)α.
2 − (λ λ )2 ≡ −
R L
−
The resulting numerical flux is not strictly monotone. However, we can guarantee entropy stability. We have
already shown discrete entropy stability for LF, HLLω and LW. Thus, to demonstrate ES for HLLXω, it
suffices to show the positivity of the coefficients β (ω), i=0,1,2. Note that for ω =1 we obtain LW, which
i
has already been shown to be discretely ES. Thus, let us consider ω [0,1). The claim holds since
∈
(cid:12) (cid:12)
β 2 (ω)=β(ω) 0 ω+(1 ω)α 0 α 0 λ R λ L (cid:12)λ R λ L (cid:12) 0. (4.4)
≥ ⇔ − ≥ ⇔ ≥ ⇔ − − | |−| | ≥
Hence, we directly see that β 0 because it is a combination of non-negative terms. In order to show that
0
≥
β (ω) 0 we consider the equivalent expression
1
≥
(cid:18) (cid:19) (cid:18) (cid:19)
1 ω 1 ω
β(ω) − +ω ω+(1 ω)α − +ω α(λ + λ ) 1. (4.5)
L R
≤ λ + λ ⇔ − ≤ λ + λ ⇔ | | | | ≤
L R L R
| | | | | | | |
(cid:12) (cid:12)
Wedistinguishtwocases. Case1: λ L andλ R arebothpositiveorbothnegative. Thenλ R λ L =(cid:12)λ R λ L (cid:12)
− | |−| |
and therefore α=0. Case 2: If λ and λ are of opposite sign, then λ λ = λ + λ , which yields
L R R L R L
− | | | |
α(λ
R
+ λ
L
)=1 (cid:12) (cid:12)λ
R
λ
L
(cid:12) (cid:12)(λ
R
+ λ
L
)/(λ
R
λ
L
) 2 1. (4.6)
| | | | − | |−| | | | | | − ≤
5. Numerical example - application to ideal MHD
WenowapplythenewhybridES-HLLXωRiemannsolvertotheequationsofidealmagnetohydrodynamics
(MHD). The ideal MHD equations are a hyperbolic system that decribes the flow of plasma assuming infinite
electric resistivity, see e.g. [10]. As a proof of concept we implement the hybrid ES-HLLXω solver into a first
orderfinitevolumeframework. Weusea1DshocktubeproblemfortheidealMHDequationstodemonstrate
the reduced dissipation of the new hybrid numerical flux. We consider the magnetic shock tube of Torrilhon
[9]
(cid:40)
(cid:2) (cid:3)T [1,0,0,0,1,1.5,0.5,0.6] T , if x 0,
(cid:37),(cid:37)u,(cid:37)v,(cid:37)w,p,B 1 ,B 2 ,B 3 = [1,0,0,0,1,1.5,1.6,0.2] T , if x ≤ >0, (5.1)
4

1.6
1.4
1.2
1.0
0.8
0.6
0.4
4 3 2 1 0 1 2 3 4
− − − −
x
B 2
0.9
ES-HLLXω,ω=0.925
ES-Roe
0.8
Exact
0.7
0.6
0.5
0.4
0.3
0.2
0.1
4 3 2 1 0 1 2 3 4
− − − −
x
B 3
Fig. 2. Comparison of the computed solution of B and B using ES-Roe (dashed) and ES-HLLXω (solid
2 3
with knots) with ω =0.25 for the magnetic shock tube problem (5.1) at T =1.0 on 300 regular grid cells.
on the domain Ω=[ 4,4] with Dirichlet boundary conditions and an adiabatic index γ = 5. The baseline
− 3
EC flux needed for the numerical flux ansatz (2.2) is chosen to be the entropy conserving and kinetic energy
preserving flux found in [10, App. B]. Briefly, we note that the non-linear stability of the scheme depends on
which underlying baseline EC flux is chosen to build the scheme. However, in our experience, the non-linear
stability properties of a scheme created with the available EC baseline fluxes is nearly identical for low Mach
number test cases (like the magnetic shock tube). We compare the ES-Roe flux of [10] and ES-HLLXω with
ω =0.925.
EachdissipationtermD=D(q ,q )intheES-HLLXω dissipationmatrix(4.3)iscreatedusingasimple
i i+1
arithmetic mean state of the primitive variables. We note that the value of ω could be chosen adaptively
using for example a pressure switch [2]. Fig. 2 presents the computed solution of B and B on 300 regular
2 3
grid cells against the exact solution of the Riemann problem at the final time T = 1.0. We see that the
entropy stable hybrid numerical flux has less dissipation than the entropy stable Roe-type scheme.
6. Conclusion
In this work we constructed two one-parameter families of hybrid entropy stable numerical fluxes. An
advantage of the new numerical flux functions is that they remain applicable even when the eigenstructure of
the flux Jacobian matrix is unknown. The derivations and proofs in this work are kept general such that the
hybrid entropy stable solvers can be applied to a broad range of non-linear hyperbolic conservation laws.
As an example, we applied the novel numerical fluxes to the ideal MHD equations and demonstrated the
decreased magnitude of dissipation for the hybrid solvers versus a standard solver. In the future we plan to
apply the hybrid entropy stable Riemann solvers to other complex hyperbolic systems, such as the two-layer
shallow water or the regularized 13-Moment Equations of Grad.
References
References
[1] TimothyJ.Barth.Numericalmethodsforgasdynamicsystemsonunstructuredmeshes.InDietmarKröner,MarioOhlberger,
andChristianRohde,editors,An Introduction to Recent Developments in Theory and Numerics for Conservation Laws,
volume5ofLecture Notes in Computational Science and Engineering,pages195–285.SpringerBerlinHeidelberg,1999.
5

[2] PraveenChandrashekar. KineticenergypreservingandentropystablefinitevolumeschemesforcompressibleEulerand
Navier-Stokesequations. Communications in Computational Physics,14:1252–1286,2013.
[3] P.Degond,P.-F.Peyrard,G.Russo,andP.Villedieu. Polynomialupwindschemesforhyperbolicsystems. ComptesRendus
| de l’Académie | des Sciences | – Series I – | Mathematics,328(6):479–483,1999. |
| ------------- | ------------ | ------------ | -------------------------------- |
[4] Ulrik S. Fjordholm, Siddhartha Mishra, and Eitan Tadmor. Arbitrarily high-order accurate entropy stable essentially
nonoscillatoryschemesforsystemsofconservationlaws. SIAM Journal on Numerical Analysis,50(2):544–573,2012.
[5] RandallJ.LeVeque. Finite volume methods for hyperbolic problems. CambridgeUniversityPress,firstedition,2002.
[6] PhilipL.Roe. Affordable,entropyconsistentfluxfunctions. InEleventhInternationalConferenceonHyperbolicProblems:
| Theory, Numerics | and Applications,Lyon,2006. |     |     |
| ---------------- | --------------------------- | --- | --- |
[7] BirteSchmidtmannandManuelTorrilhon. AhybridRiemannsolverforlargehyperbolicsystemsofconservationlaws.
| SIAM Journal | on Scientific | Computing | (submitted),arXiv:1607.05721,2016. |
| ------------ | ------------- | --------- | ---------------------------------- |
[8] EitanTadmor. Entropystabilitytheoryfordifferenceapproximationsofnonlinearconservationlawsandrelatedtime-
| dependentproblems. | Acta | Numerica,12:451–512,2003. |     |
| ------------------ | ---- | ------------------------- | --- |
[9] Manuel Torrilhon. Uniqueness conditions for Riemann problems of ideal magnetohydrodynamics. Journal of Plasma
Physics,69(3):253–276,2003.
[10] AndrewR.WintersandGregorJ.Gassner. Affordable,entropyconservingandentropystablefluxfunctionsfortheideal
| MHDequations. | Journal | of Computational | Physics,304:72–108,2016. |
| ------------- | ------- | ---------------- | ------------------------ |
6