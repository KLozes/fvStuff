|     |        | Genuinely |        | multi-dimensional |     |     |     | stationarity |            | preserving |      |     |
| --- | ------ | --------- | ------ | ----------------- | --- | --- | --- | ------------ | ---------- | ---------- | ---- | --- |
|     | Finite |           | Volume | formulation       |     |     | for | nonlinear    | hyperbolic |            | PDEs |     |
Wasilij Barsukowa, Mirco Ciallellab, Mario Ricchiutoc, Davide Torlod
aInstitut de Mathématiques de Bordeaux, Université de Bordeaux, CNRS UMR 5251, Talence, France
bLaboratoire
Jacques-Louis Lions, Université Paris Cité, CNRS UMR 7598, Paris, France
cCentre Inria de l’Université de Bordeaux, CNRS UMR 5251, Talence, France
|     |     |     | dDipartimento | di  | Matematica, | Università |     | di Roma | La Sapienza, | Rome, | Italy |     |
| --- | --- | --- | ------------- | --- | ----------- | ---------- | --- | ------- | ------------ | ----- | ----- | --- |
5202 ceD 51  ]AN.htam[  2v00712.6052:viXra
Abstract
Classical Finite Volume methods for multi-dimensional problems include stabilization (e.g.
via a Riemann solver), that is derived by considering several one-dimensional problems in
different directions. Such methods therefore ignore a possibly existing balance of contri-
butions coming from different directions, such as the one characterizing multi-dimensional
stationary states. Instead of being preserved, they are usually diffused away by such meth-
ods. Stationarity preserving methods use a better suited stabilization term that vanishes
at the stationary state, allowing the method to preserve it. This work presents a general
approach to stationarity preserving Finite Volume methods for nonlinear conservation/bal-
ance laws. It is based on a multi-dimensional stationarity preserving quadrature strategy
that allows to naturally introduce genuinely multi-dimensional numerical fluxes. The new
methods are shown to significantly outperform existing ones even if the latter are of higher
| order | of  | accuracy | and | even | on non-stationary |     | solutions. |     |     |     |     |     |
| ----- | --- | -------- | --- | ---- | ----------------- | --- | ---------- | --- | --- | --- | --- | --- |
Keywords: Stationarity preservation, Finite Volume, Multi-dimensional well-balancing,
| Hyperbolic |     | equations, |     | Global | flux, | residual | distribution |     |     |     |     |     |
| ---------- | --- | ---------- | --- | ------ | ----- | -------- | ------------ | --- | --- | --- | --- | --- |
1. Introduction
This paper focuses on the numerical solution of nonlinear hyperbolic systems of conser-
| vation | laws | in  | two | dimensions: |     |        |      |        |     |     |     |     |
| ------ | ---- | --- | --- | ----------- | --- | ------ | ---- | ------ | --- | --- | --- | --- |
|        |      |     |     |             |     | ∂ q +∂ | f +∂ | g = 0, |     |     |     | (1) |
|        |      |     |     |             |     | t      | x    | y      |     |     |     |     |
where q, f and g are the vectors of conservative variables and fluxes. Numerical methods for
hyperbolic partial differential equations (PDEs) need numerical diffusion to achieve entropy
stability and in order to deal with solutions characterized by strong gradients. The majority
of numerical methods for multi-dimensional problems, though, are developed following a
dimension-by-dimension approach, meaning that the numerical diffusion is usually derived
in a one-dimensional framework and that the diffusion term associated to an edge (or a face,
in 3D) usually involves only two states. Standard numerical methods with one-dimensional
| Riemann |     | solvers | typically | introduce |      | a diffusion |        | term of the | type |       |     |     |
| ------- | --- | ------- | --------- | --------- | ---- | ----------- | ------ | ----------- | ---- | ----- | --- | --- |
|         |     |         |           | ∂ q +∂    | f +∂ | g =         | ∆x∂ (ν | ∂ q)+∆y∂    | (ν   | ∂ q), |     | (2) |
|         |     |         |           | t         | x    | y           | x      | x x         | y y  | y     |     |     |

where ∆x and ∆y provide the size of the discretization, and ν and ν represent the diffusion
x y
coefficients, which are often chosen proportional to the spectral radius of the flux Jacobian.
This one-dimensional approach does not take into account possible multi-dimensional fea-
tures of the numerical solution, such as the stationary states characterized by a balance of
contributions coming from different directions [10]. For equation (1), stationary states are
governed by
∂ f +∂ g = 0. (3)
x y
Forclassicalmethods,withthetwo-dimensionaldiffusiontermdesignedfollowingone-dimensional
approaches, the solution will be completely diffused instead of being kept stationary [11, 13].
In contrast to Equation (3), the discrete stationary states are characterized by the much
more restrictive conditions ∂ f = 0 and ∂ g = 0. States where ∂ f ̸= 0 is balanced by
x y x
−∂ g is not a stationary state of the numerical method. This can be prevented by choosing
y
more sophisticated diffusion operators [46, 52, 34, 45]. Such methods are called stationarity
preserving [10].
Typically, solutions to (3) form a very large set; the equations might even be underde-
termined, as in the case, for example, in linear acoustics. In that case, they reduce merely
to the condition that the velocity field must be divergence-free. If the divergence is under-
stood weakly, this operator is not invertible even if all boundary conditions are prescribed.
Consequently, no numerical method can exactly preserve all stationary states whenever their
set is so rich. If we consider again the example of linear acoustics, given a finite set of
point values or averages, it is fundamentally impossible to establish whether they belong to
a divergence-free vector field or not.
Stationarity preserving methods guarantee that the stationary states of a numerical
method are a discretization of all the stationary states of the PDE. The adjective all is
very important here, and consistency of a numerical method is not enough to guarantee this.
As shown in [12], an important necessary condition for a method to be stationarity preserving
is that there exists a finite set of local approximations of (3) that characterize part of the
kernel of the numerical method. In other words, when these local operators vanish, then
the consistent part and the numerical dissipation of the numerical method vanish simultan-
eously. For this condition to be also sufficient, and for the stationary state uniquely defined,
the number of such approximations should be equal to the total number of unknowns. In
the context of linear problems, an alternative definition of a stationarity preserving method
(using the discrete Fourier transform) is given in [10].
These properties can hardly be obtained in the context of classical methods. As described
above, all too often the stationary states of a (consistent) numerical method are discretiz-
ations of ∂ f = 0 and ∂ g = 0, instead of (3). Consider, as another example, the trivial
x y
equation ∂ q = 0, which keeps everything stationary. A numerical method that adds diffu-
t
sion like ∂ q = ∆x∂2q will be consistent, but then only those q which are affine in space will
t x
be stationary states of the numerical method. In this paper, we aim to develop nonlinear
methods whose stationary states are discretizations of all the stationary states of the PDE,
while being stable under explicit time integration.
Recent examples of stationarity preserving diffusion operators were developed for geo-
strophic equilibria in the linear and nonlinear case [6, 7]. A connection to these equilibria in
2

the context of low Mach number limit of the Euler equations, which is related to the long-
time limit of linear acoustics, was provided in [36, 37] through the preservation of discrete
divergence with ad-hoc functional spaces. Early examples of stationarity preserving methods
for nonlinear conservation laws can be found in [9]. So far, however, no general theory for
the agnostic detection of stationary states of nonlinear multi-dimensional hyperbolic partial
| differential | equations | is available. |     |     |     |     |     |     |
| ------------ | --------- | ------------- | --- | --- | --- | --- | --- | --- |
The method presented in this work exploits an idea formulated in [13], which allows to
modifythequadraturestrategyinawaythatsystematicallyleadstoastationaritypreserving
method. Thismodificationisrelatedtoatechniqueoftenreferredtoastheglobalfluxmethod
[32, 21, 23] initially introduced for hyperbolic balance laws in one dimension,
|     |     |     |     |     | ∂ q +∂ | f = | s,  | (4) |
| --- | --- | --- | --- | --- | ------ | --- | --- | --- |
|     |     |     |     |     | t      | x   |     |     |
with the original goal of developing well-balanced methods [5, 17, 22], and the treatment of
sourcetermspresentinthemathematicalmodel. Theglobalfluxhasalreadybeensuccessfully
applied to different contexts and numerical methods [24, 25, 28, 41, 42, 38, 4] to preserve
one dimensional equilibria. The term “global” is unfortunately misleading. Indeed, quite
interestingly this approach leads to fully local discretizations (see e.g. [41, 42, 13, 4, 15]), and
indeed it is related to similar techniques used in the analysis of well balanced and asymptotic
preserving methods, and referred to as the spatial localization of source terms e.g. in [33]).
These methods also have very natural connections with the so-called residual distributions
| schemes | [3, 1]. |     |     |     |     |     |     |     |
| ------- | ------- | --- | --- | --- | --- | --- | --- | --- |
The underlying idea of these techniques is to recast the source term as a flux R:
(cid:90)
x
|     |     |     |     | ∂ R | = s ⇒ | R := | sdx. | (5) |
| --- | --- | --- | --- | --- | ----- | ---- | ---- | --- |
x
| In this framework, |     | equation | (4) | can be | recast | as     |      |     |
| ------------------ | --- | -------- | --- | ------ | ------ | ------ | ---- | --- |
|                    |     |          |     | ∂      | q +∂   | (f −R) | = 0, | (6) |
t x
| where discrete | steady | states | satisfy | the | relation |     |     |     |
| -------------- | ------ | ------ | ------- | --- | -------- | --- | --- | --- |
F
|     |     |     | ∂   | (f −R) | =   | 0 ⇔ f −R | ≡ , | (7) |
| --- | --- | --- | --- | ------ | --- | -------- | --- | --- |
|     |     |     |     | x      |     |          | 0   |     |
with F = f −R often referred to as global flux, as it embeds the entire evolution operator,
and F = F(x ) for a given x in the domain. In the same spirit, a similar approach that
| 0   |     | 0   | 0   |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
integrates the Coriolis term into an apparent bathymetry term was also developed in [19].
The concept of well-balancing is a particular case of the preservation of general stationary
solutions. Theoverarchingideaistodesignnumericalschemesinwhichtheartificialdiffusion
vanishesatrelevantequilibria. Thedevelopmentofwell-balancedschemesinone-dimensional
problems has reached high levels of maturity in the last decades, but the multi-dimensional
extensions are often tackled with trivial dimension-by-dimension approaches [24, 26, 44, 41],
whichonlyallowsthepreservationof1D-likeequilibria. In[13,15]someofthepresentauthors
propose, in the context of stabilized finite elements, a modified quadrature approach that
allows to guarantee the multi-dimensional stationarity preserving properties of the resulting
3

discretization. The underlying idea has some connections with the above mentioned global
flux methods. In two space dimensions, it is based on the idea of recasting the problem as
|     |     |     |     |     |     |     | (cid:18)(cid:90) | (cid:19) |     |
| --- | --- | --- | --- | --- | --- | --- | ---------------- | -------- | --- |
y
|     |     | ∂   | q +∂ | f +∂ | g = ∂ | q +∂ |     | ∂ f dy +g |     |
| --- | --- | --- | ---- | ---- | ----- | ---- | --- | --------- | --- |
|     |     |     | t x  |      | y t   |      | y   | x         |     |
(8)
|     |     |     |     |     |     |      | (cid:18) | (cid:90) x (cid:19) |     |
| --- | --- | --- | --- | --- | --- | ---- | -------- | ------------------- | --- |
|     |     |     |     |     | = ∂ | q +∂ | f +      | ∂ gdx = 0.          |     |
|     |     |     |     |     | t   |      | x        | y                   |     |
By symmetry, one combines the two above formulations, replacing the two-dimensional the
| conservation | law | (1) with |     |      |        |     |      |             |     |
| ------------ | --- | -------- | --- | ---- | ------ | --- | ---- | ----------- | --- |
|              |     |          | ∂ q | +∂ f | +∂ g = | ∂ q | +∂ ∂ | (F +G) = 0, | (9) |
|              |     |          | t   | x    | y      | t   | x    | y           |     |
by defining
|     |     |     | (cid:90) | y     |     |     |     | (cid:90) x |      |
| --- | --- | --- | -------- | ----- | --- | --- | --- | ---------- | ---- |
|     |     | F   | :=       | f dy, |     |     |     | G := gdx.  | (10) |
Thenewdivergenceoperator∂ ∂ (F+G)nowiseasytopreserveatthediscretelevel. Thanks
x y
to this formulation, it becomes also straightforward to consider multi-dimensional balance
| laws with     | source     | terms, |        |     |            |           |        |     |      |
| ------------- | ---------- | ------ | ------ | --- | ---------- | --------- | ------ | --- | ---- |
|               |            |        |        |     | ∂ q +∂     | f +∂      | g =    | s.  | (11) |
|               |            |        |        |     | t          | x         | y      |     |      |
| In this case, | the source |        | flux R | can | be defined | as        |        |     |      |
|               |            |        |        |     | (cid:90)   | x(cid:90) | y      |     |      |
|               |            |        |        |     | R :=       |           | sdxdy, |     | (12) |
and directly included in the modified divergence operator. Setting
F
|     |     |     |     |     | =   | F +G−R, |     |     | (13) |
| --- | --- | --- | --- | --- | --- | ------- | --- | --- | ---- |
the multi-dimensional stationarity preserving formulation of the original problem becomes
again
|     |     |     |     |     | ∂ q +∂ |     | F = 0. |     | (14) |
| --- | --- | --- | --- | --- | ------ | --- | ------ | --- | ---- |
|     |     |     |     |     | t      | xy  |        |     |      |
Inthiswork, wepresenthowthisideacanbeusedtodesignfirst-orderfinitevolumemeth-
ods preserving multi-dimensional steady states, not known a priori, for general nonlinear
hyperbolic PDEs. A thorough analysis of the method is presented for linear problems
showing a link with other stationarity preserving methods, as well as a discrete energy estim-
ate. The approach proposed here naturally leads to the introduction of nonlinear genuinely
multi-dimensional fluxes at cell corners, whichhavebeenshowntoprovidefundamental
enhancements in the numerical solutions, and enjoy many theoretical properties ([31, 12, 1]).
However, differently from previous works, the formulation proposed here naturally leads to
corner fluxes, without any hypotheses on the type of quadrature, or on the Riemann fluxes.
Thisapproachisastartingpointforthedevelopmentofnewfamiliesofstationaritypreserving
high-order methods based on high-degree polynomial reconstruction [27, 28], or discontinu-
ous Galerkin methods [41, 53]. Throughout the paper, we will refer to the new method as
4

stationarity preserving formulation or global flux quadrature as in [13, 15]. However, note
| that the | discretization |     | is fully | local, | as  | we will | show | in  | detail. |     |     |
| -------- | -------------- | --- | -------- | ------ | --- | ------- | ---- | --- | ------- | --- | --- |
The paper is organized as follows. In section 2, we present the examples of PDEs that will
be considered when assessing the performance of the method experimentally. In section 3,
we recall the global flux method in a one-dimensional framework for hyperbolic balance laws.
In section 4, we present the two-dimensional stationarity preserving/global flux quadrature
approach for hyperbolic PDEs. Here, we discuss the finite volume formulation, the stabiliza-
tion technique, boundary conditions, the treatment of source terms, as well as stability and
consistency of the method for a linear model. In section 5, we present the standard finite
volume method with piecewise constant and piecewise linear reconstructions used for com-
parison with the global flux method. Several numerical experiments are presented in section
6 to show the performance of the method. Finally, we draw some conclusions in section 7.
| 2. Mathematical |     |     | models |     |     |     |     |     |     |     |     |
| --------------- | --- | --- | ------ | --- | --- | --- | --- | --- | --- | --- | --- |
The numerical method presented in this work is rather general and, in order to show
its potential, we exemplify it on several mathematical models described by both linear and
nonlinear hyperbolic systems. In particular, herein we will focus on three systems: linear
acoustics, Euler equations for gas dynamics and the shallow water equations. For all of them,
| we focus    | on two-dimensional |     |        | problems. |     |     |     |     |     |     |     |
| ----------- | ------------------ | --- | ------ | --------- | --- | --- | --- | --- | --- | --- | --- |
| 2.1. Linear | acoustic           |     | system |           |     |     |     |     |     |     |     |
Thesystemoflinearacousticisasimplemodelthatdirectlyembedsnon-trivialdivergence-
free steady states. It can be written in the following 2D and vectorial forms as:

|     |     | ∂ u+∂ | p   |     | = 0, |     |     | (cid:40) |      |      |     |
| --- | --- | ----- | --- | --- | ---- | --- | --- | -------- | ---- | ---- | --- |
|     |     |   t | x   |     |      |     |     | ∂        | v+∇p | = 0, |     |
t
|     |     | ∂ v  | +∂ p |     | = 0,   |     |     |     |       |      | (15) |
| --- | --- | ---- | ---- | --- | ------ | --- | --- | --- | ----- | ---- | ---- |
|     |     | t    | y    |     |        |     |     |     |       |      |      |
|     |     |      |      |     |        |     |     | ∂   | p+∇·v | = 0, |      |
|     |     |  ∂ |      |     |        |     |     |     | t     |      |      |
|     |     | p+∂  | u+∂  |     | v = 0, |     |     |     |       |      |      |
|     |     | t    | x    | y   |        |     |     |     |       |      |      |
where p is the pressure and v = (u,v) is the velocity. The system can also be written in the
| compact    | form      | (1) with |       |          |           |     |      |     |          |        |      |
| ---------- | --------- | -------- | ----- | -------- | --------- | --- | ---- | --- | -------- | ------ | ---- |
|            |           |          |       |         |          |     |    |     |        |        |      |
|            |           |          |       |          | u         |     | p    |     | 0        |        |      |
|            |           |          |       | q = v, |           | f = | 0, |     | g = p. |        | (16) |
|            |           |          |       |          | p         |     | u    |     | v        |        |      |
| The steady | states    | of       | this  | system   | are given | by  |      |     |          |        |      |
|            |           |          | ∂ q ≡ | 0        | ⇔         | ∇·v | ≡ 0  | and | p ≡ p =  | const. | (17) |
|            |           |          | t     |          |           |     |      |     | 0        |        |      |
| 2.2. Euler | equations |          |       |          |           |     |      |     |          |        |      |
The Euler equations are a simplification of the full Navier-Stokes system that do not
include viscosity effects. Their use is widespread for the simulation of compressible gas
| dynamics. | The | system | can | be written |     | in vectorial |     | form | as: |     |     |
| --------- | --- | ------ | --- | ---------- | --- | ------------ | --- | ---- | --- | --- | --- |

|     |     |     |     | ∂   | ρ+∇·(ρv) |     |     |     | = 0, |     |     |
| --- | --- | --- | --- | --- | -------- | --- | --- | --- | ---- | --- | --- |
|     |     |     |     |    | t        |     |     |     |      |     |     |

|     |     |     |     | ∂   | (ρv)+∇·(ρv⊗v+pI) |     |     |     | = 0, |     | (18) |
| --- | --- | --- | --- | --- | ---------------- | --- | --- | --- | ---- | --- | ---- |
t

|     |     |     |     | ∂  | (ρE)+∇·(ρHv) |     |     |     | = 0, |     |     |
| --- | --- | --- | --- | --- | ------------ | --- | --- | --- | ---- | --- | --- |
t
5

e+∥v∥2/2
having denoted by ρ the density, by v the velocity field, by E = the specific total
energy, being e the specific internal energy and I is the identity matrix. Finally, the total
specific enthalpy is H = h + ∥v∥2/2, with h = e + p/ρ the specific enthalpy. To close the
system, we use the classical perfect gas equation of state p = (γ −1)ρe with γ the constant
| ratio | of specific | heats | (γ = | 1.4 for | air). |     |     |     |     |     |     |     |
| ----- | ----------- | ----- | ---- | ------- | ----- | --- | --- | --- | --- | --- | --- | --- |
The nonlinear system of Euler equations can also be recast in the compact form (1) with
|     |     |     |     |    |     |      |    |     |       |     |    |      |
| --- | --- | --- | ---- | --- | --- | ----- | --- | --- | ------ | --- | --- | ---- |
|     |     |     |      | ρ   |     | ρu    |     |     |        | ρv  |     |      |
|     |     |     | ρu |     |     | ρu2  | +p |     |       | ρuv |    |      |
|     |     |     | q =  | ,  | f = |       | ,  | g   | =      |     | .  | (19) |
|     |     |     |     |     |     |      |     |     |  ρv2 |     |     |      |
|     |     |     | ρv |     |     |  ρuv |    |     |        | +p |     |      |
|     |     |     |      | ρE  |     | ρHu   |     |     |        | ρHv |     |      |
Steady states of the Euler equations are more complex but, after some manipulations, the
smooth steady states can be characterized by the following relations:
|      |         | ∇·(ρv) |        | = 0, | (ρv·∇)v+∇p |     |     | = 0, | v·∇H |     | = 0. | (20) |
| ---- | ------- | ------ | ------ | ---- | ---------- | --- | --- | ---- | ---- | --- | ---- | ---- |
| 2.3. | Shallow | water  | system |      |            |     |     |      |      |     |      |      |
The Saint-Venant or shallow water equations describe the dynamics of hydrostatic free
surface waves influenced by gravity. This model is valid under the hypothesis of very large
wavelengths, or very shallow depths, and is applied in various engineering fields, including
river and estuarine hydrodynamics, urban flood management, and tsunami risk assessment.
In particular, when working with large scale problems, this simplified model becomes crucial
| to speed-up |        | the computational |            | time. |           |      |     |     |     |     |     |     |
| ----------- | ------ | ----------------- | ---------- | ----- | --------- | ---- | --- | --- | --- | --- | --- | --- |
| The         | system | can               | be written | in    | vectorial | form | as: |     |     |     |     |     |
(cid:40)
|     |     |     | ∂ h+∇·(hv) |     |         |     |       | =       | 0,     |     |     |      |
| --- | --- | --- | ---------- | --- | ------- | --- | ----- | ------- | ------ | --- | --- | ---- |
|     |     |     | t          |     |         |     |       |         |        |     |     | (21) |
|     |     |     |            |     | (cid:0) |     |       | (cid:1) |        |     |     |      |
|     |     |     | ∂ (hv)+∇·  |     | hv⊗v+   |     | 1gh2I | =       | −gh∇b, |     |     |      |
t
2
where h is the water height, v the velocity field, b is the bathymetry and g is the gravity
constant. The system can also be written in the classical compact notation (11) with
|     |    |    |    |     |    |     |    |     |    |     |   |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
|     | h   |     |     | hu  |     |     |     | hv  |     |     | 0   |     |
1gh2
| q   | = hu, |     | f = hu2 | +   | ,  | g   | =  | huv  | ,  | s   | = −gh∂ b. | (22) |
| --- | ------- | --- | -------- | --- | --- | --- | --- | ---- | --- | --- | ----------- | ---- |
|     |         |     |          |     | 2   |     |     |      |     |     | x           |      |
|     |         |     |          |     |     |     | hv2 | 1gh2 |     |     |             |      |
|     | hv      |     |          | huv |     |     |     | +    |     |     | −gh∂ b      |      |
|     |         |     |          |     |     |     |     | 2    |     |     | y           |      |
This system admits a large variety of equilibria depending on the interaction between the
flux and the source. The most studied equilibria in the context of well-balanced methods
are the so-called “lake at rest” states, which are characterized by a constant free surface
level η := h + b ≡ η and a zero velocity v ≡ 0. However, in the presence of a non-flat
0
bathymetry, the system can also admit non-trivial equilibria, which are characterized by a
| non-zero | velocity | and    | a non-flat |     | free surface |     | level:         |     |     |     |     |      |
| -------- | -------- | ------ | ---------- | --- | ------------ | --- | -------------- | --- | --- | --- | --- | ---- |
|          |          | ∇·(hv) |            | = 0 |              |     | (v·∇)v+g∇(h+b) |     |     |     | = 0 | (23) |
Several works have been devoted to the study of these equilibria in one dimension or in
a quasi-1D framework [43, 41, 26]. In this work, we are interested in truly multi-dimensional
well-balanced schemes that are capable of preserving all these equilibria at the discrete level.
6

| 3.  | Global | flux | for | 1D  | balance | laws |     |     |     |     |     |     |     |
| --- | ------ | ---- | --- | --- | ------- | ---- | --- | --- | --- | --- | --- | --- | --- |
In this section, we recall the main principle of the global flux method and its initial usage
F
in a 1D framework. Consider a general nonlinear balance law (4) and define a global flux
as
|     |     |     |     |     |     |     |       | (cid:90) | x    |     |     |     |      |
| --- | --- | --- | --- | --- | --- | --- | ----- | -------- | ---- | --- | --- | --- | ---- |
|     |     |     |     |     |     |     | F = f | −        | sdx, |     |     |     | (24) |
such that (4) can now be written in a pseudo-conservative form as
F
|     |     |     |     |     |     |     | ∂ q +∂ |     | = 0. |     |     |     | (25) |
| --- | --- | --- | --- | --- | --- | --- | ------ | --- | ---- | --- | --- | --- | ---- |
|     |     |     |     |     |     |     | t      | x   |      |     |     |     |      |
Steady states given by ∂ q = 0 are equivalently characterized by the condition ∂ F = 0 ⇔
|     |     |      |     | t   |     |     |     |     |     |     |     | x   |     |
| --- | --- | ---- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| F   | ≡ F | ∈ R. |     |     |     |     |     |     |     |     |     |     |     |
0
In a finite volume framework, the computational domain Ω is split into N cells and the
|          |     |        |            |     |           |      |          | (cid:104)  |          | (cid:105) |     |     |      |
| -------- | --- | ------ | ---------- | --- | --------- | ---- | -------- | ---------- | -------- | --------- | --- | --- | ---- |
| equation |     | (4) is | integrated |     | over each | cell | C        | = x        | ,x       | :         |     |     |      |
|          |     |        |            |     |           |      | i        | i−1        | i+1      |           |     |     |      |
|          |     |        |            |     |           |      |          |            | 2        | 2         |     |     |      |
|          |     |        |            |     |           |      | F        | −F         |          |           |     |     |      |
|          |     |        |            |     |           | d    | (cid:99) |            | (cid:99) |           |     |     |      |
|          |     |        |            |     |           |      | i+1      |            | i−1      |           |     |     |      |
|          |     |        |            |     |           | q¯   | +        | 2          | 2        | = 0,      |     |     | (26) |
|          |     |        |            |     |           | dt i |          | ∆x         |          |           |     |     |      |
| where    | the | cell   | average    | is  | defined   | as   |          |            |          |           |     |     |      |
|          |     |        |            |     |           |      | 1        | (cid:90) x | 1        |           |     |     |      |
i+2
|     |     |     |     |     |     | q¯  | :=  |     | qdx. |     |     |     | (27) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | ---- | --- | --- | --- | ---- |
i
∆x
x i−1
2
F
The numerical global flux (cid:99) is considered to be a function of the two values of the global
i+1
2
flux FL and FR reconstructed at both sides of interface x . For piecewise constant
|     | i+1 |     | i+1 |     |     |     |     |     |     |     | i+1 |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
2
|                 |     | 2   |        | 2      |      |     |        |     | FL  | F   | FR  | F       |     |
| --------------- | --- | --- | ------ | ------ | ---- | --- | ------ | --- | --- | --- | --- | ------- | --- |
| reconstructions |     |     | of the | global | flux | one | simply | has |     | =   | and | = .     |     |
|                 |     |     |        |        |      |     |        |     | i+1 | i   |     | i+1 i+1 |     |
|                 |     |     |        |        |      |     |        |     | 2   |     |     | 2       |     |
Remark 1 (Numerical global flux). It is important to underline that structure preservation
can only be achieved if the interface global flux F (cid:99) is constant ∀ i at steady state. That
i+1
2
implies that numerical dissipation should also vanish in correspondence of this state. Unless
this condition is imposed with some ad hoc modification (see e.g. [4]), this is only possible if
{F
the numerical dissipation depends only on global fluxes } in the cells, and not on the
|     |     |     |     |     |     |     |     |     |     | j   | j∈Z |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
values {q } of the conservative variables. This is due to the fact that, at equilibria, only
j j∈Z
| global | fluxes | are | constant |     | while conservative |     |     | variables |     | may vary. |     |     |     |
| ------ | ------ | --- | -------- | --- | ------------------ | --- | --- | --------- | --- | --------- | --- | --- | --- |
In our previous work [28], we employed the following upwind flux:
1±signΛ
|     |     |     | F        | 1+F |       | +1−F |      |       |     |        |     |     |      |
| --- | --- | --- | -------- | --- | ----- | ---- | ---- | ----- | --- | ------ | --- | --- | ---- |
|     |     |     | (cid:99) | =   | L     |      | R    | where | 1±  | := L−1 |     | L,  | (28) |
|     |     |     | i+       | 1   | i + 1 |      | i +1 |       |     |        |     | 2   |      |
|     |     |     |          | 2   | 2     |      | 2    |       |     |        |     |     |      |
where L is the matrix of left eigenvectors of the flux Jacobian ∂ f, and signΛ is the diagonal
q
matrix of the sign of the eigenvalues of the flux Jacobian, evaluated using any (average) state
at the interface1.
1±
1In principle, therefore, one should write to make clear that they differ from interface to interface;
i+1
2
| we do | not | make | this depends |     | explicit | to ensure | readability |     |     |     |     |     |     |
| ----- | --- | ---- | ------------ | --- | -------- | --------- | ----------- | --- | --- | --- | --- | --- | --- |
7

(cid:82)x
Forthedevelopmentoftheglobalfluxmethod, aconsistentapproximationofR := sdx
is necessary to define F = f − R. This integral can be computed in a recursive manner,
starting from the beginning of the domain, by integrating the source in each element. To
simplify the description of the method, we will assume that q and s are constant in each cell,
therefore the source integral can be computed as
R :=
(cid:90) xi−1
sdx+
(cid:90) x i−2 1
sdx+
(cid:90) xi
sdx = R +
∆x
s¯ +
∆x
s¯. (29)
i i−1 i−1 i
2 2
(cid:124) (cid:123)(cid:122) (cid:125) xi−1 x i−1 2
Ri−1
Hence, the global flux will now depend on both the conservative flux and the source term
∆x
F = f(q¯)−R = f(q¯)−R − (s¯ +s¯). (30)
i i i i i−1 i−1 i
2
Similarly, the recursive procedure gives us the following values for F and F :
i−1 i+1
F = f(q¯ )−R , (31)
i−1 i−1 i−1
(cid:18) (cid:19)
1 1
F = f(q¯ )−R −∆x s¯ +s¯ + s¯ . (32)
i+1 i+1 i−1 i−1 i i+1
2 2
It can be noticed that, when considering a simple numerical flux
F (cid:99) i+1 (F i ,F i+1 ) = F i (33)
2
equation (26) can be recast as
d F −F f(q¯)−f(q¯ ) s¯ +s¯
i i−1 i i−1 i i−1
q¯ = − = − + , (34)
i
dt ∆x ∆x 2
which shows already a difference with respect to the classical finite volume method, where
the source term would be treated in a centered way.
When the upwind numerical flux (28) is used, one has (having temporarily made the
dependence of 1± on the interface explicit)
1− f(q¯ )+(1+ +1− )f(q¯)−1+ f(q¯ )
d i+1 i+1 i+1 i−1 i i−1 i−1
q¯ = − 2 2 2 2 (35)
i
dt ∆x
s¯ +s¯ s¯ +s¯
+(1+ +1− −1− ) i i−1 +1− i i−1 +(1+ +1− −(1+ +1− ))R ,
i+1 i+1 i−1 2 i+1 2 i+1 i+1 i−1 i−1 i−1
2 2 2 2 2 2 2 2
(cid:124) (cid:123)(cid:122) (cid:125) (cid:124) (cid:123)(cid:122) (cid:125)
=1 =1
where the contribution from R cancels out even if 1± depend on the interface, since they
i−1
nevertheless add up to 1 on each of them.
Remark 2 (Compactness of global fluxes). It can be noticed that, although the global flux
in (30) is defined globally with R that depends on previous values, the time residual (34)
i−1
8

shows that the stencil is actually compact due to the cancellation of these terms. In particular,
simple manipulations show that the above upwind method can be also written as
d
F
(cid:99) i+1
−F
i
F
i
−F
(cid:99) i−1
q¯ = − 2 − 2
i
dt ∆x ∆x
1− (F −F ) 1+ (F −F )
i+1/2 i+1 i i−1/2 i i−1
= − −
∆x ∆x
(cid:18) (cid:19) (cid:18) (cid:19)
f(q¯ )−f(q¯) s¯ +s¯ f(q¯)−f(q¯ ) s¯ +s¯
= −1− i+1 i − i+1 i −1+ i i−1 − i i−1
i+1/2 ∆x 2 i−1/2 ∆x 2
which is the well known upwind splitting method dating back to the early works by P.L. Roe
[51], and later on Bermudez and Vazquez [16]. In the multidimensional case, a full analogy
with compact residual distribution methods on a dual cell is presented later in section 4.5.
4. Stationarity preserving formulation for multi-dimensional hyperbolic PDEs
4.1. Numerical method
When dealing with multi-dimensional conservation laws, non-trivial equilibria arise also
in absence of a source term in the equation. For steady states ∂ q = 0, it is no longer just
t
∂ f = 0 that follows, but instead the divergence ∂ f +∂ g = 0, which in general might have
x x y
many solutions.
To design stationarity preserving finite volumes schemes, we start by rewriting (1) as (9)
using the definitions in (10) to obtain
∂ q +∂ F = 0. (36)
t xy
where, by analogy with the one dimensional case, we refer to F := F +G as a global flux,
accounting for both contributions of the fluxes in the x and y direction (and later of the
source term).
(cid:104) (cid:105) (cid:104) (cid:105)
Integration of (14) over the cell C = x ,x × y ,y yields
i,j i−1 i+1 j−1 j+1
2 2 2 2
d
∆x∆y q¯ +F(t,x ,y )−F(t,x ,y )−F(t,x ,y )+F(t,x ,y ) = 0,
dt i,j i+1 2 j+1 2 i− 2 1 j+1 2 i+ 2 1 j−1 2 i−1 2 j−1 2
(37)
where the cell average is defined as
q¯ :=
1 (cid:90) x i+1
2
(cid:90) y j+1
2 qdxdy.
i,j
∆x∆y
x i−1
2
y j−1
2
From (37), we see that with this new formulation leads to the introduction of the numerical
corner fluxes F (cid:99) i±1,j±1 that then allow to write the evolution equation as
2 2
d
F
(cid:99) i+1,j+1
−F
(cid:99) i−1,j+1
−F
(cid:99) i+1,j−1
+F
(cid:99) i−1,j−1
q¯ + 2 2 2 2 2 2 2 2 = 0. (38)
i,j
dt ∆x∆y
9

|     |     |     |         |       |        |                  |     |            | C C           |         |     |
| --- | --- | --- | ------- | ----- | ------ | ---------------- | --- | ---------- | ------------- | ------- | --- |
|     |     |     |         |       |        |                  |     |            | i,j+1         | i+1,j+1 |     |
|     |     | C   |         | C     | C      |                  |     |            |               |         |     |
|     |     |     | i−1,j+1 | i,j+1 |        | i+1,j+1          |     |            |               |         |     |
|     |     |     | C       | C     |        | C                |     |            |               |         |     |
|     |     |     | i−1,j   | i,j   |        | i+1,j            |     |            |               |         |     |
|     |     |     |         |       |        |                  |     |            | C             | C       |     |
|     |     |     |         |       |        |                  |     |            | i,j           | i+1,j   |     |
|     |     | C   |         | C     | C      |                  |     |            |               |         |     |
|     |     |     | i−1,j−1 | i,j−1 |        | i+1,j−1          |     |            |               |         |     |
|     |     |     | (a)     | Main  | grid   |                  |     |            | (b) Dual grid |         |     |
|     |     |     |         |       | Figure | 1: Cell labeling |     | for the 2D | grid.         |         |     |
Recall that the global flux F accounts for contributions of the physical fluxes in all spatial
directions, as well as of the source, see equations (10), (12), and (13). In practice, the
transversally integrated fluxes F and G are computed via quadrature along, respectively, the
y and x directions in a 1D fashion. In particular, the value of F in the barycenter of a given
i
cell i can be computed recursively starting from the beginning of the domain, similarly to
| section | 3, as |     |          |      |          |        |       |     |           |     |      |
| ------- | ----- | --- | -------- | ---- | -------- | ------ | ----- | --- | --------- | --- | ---- |
|         |       |     | (cid:90) | yj−1 | (cid:90) | yj     |       | ∆y  |           |     |      |
|         |       | F   | =        | f dy | +        | f dy = | F     | +   | (f +f ),  | ∀i, | (39) |
|         |       | i,j |          |      |          |        | i,j−1 |     | i,j−1 i,j |     |      |
2
yj−1
where, for a first order method, trapezoidal rule is accurate enough. Similarly, G can be
i
| computed |     | as  |          |      |          |       |       |     |           |     |      |
| -------- | --- | --- | -------- | ---- | -------- | ----- | ----- | --- | --------- | --- | ---- |
|          |     |     | (cid:90) | xi−1 | (cid:90) | xi    |       | ∆x  |           |     |      |
|          |     | G   | =        | gdx+ |          | gdx = | G     | +   | (g +g ),  | ∀j. | (40) |
|          |     | i,j |          |      |          |       | i−1,j |     | i−1,j i,j |     |      |
2
xi−1
When dealing with hyperbolic PDEs with source terms, as for the shallow water equations
(cid:82)x(cid:82)y
in section 2.3, the integral of the source term R := sdxdy can also be embedded into
| the global |     | flux. | Similarly, | R can | be  | recursively | defined | as  |     |     |     |
| ---------- | --- | ----- | ---------- | ----- | --- | ----------- | ------- | --- | --- | --- | --- |
∆x∆y
|     | R   | =   | R     | +R    | −R  | +       |     | (s      | +s +s       | +s ), | (41) |
| --- | --- | --- | ----- | ----- | --- | ------- | --- | ------- | ----------- | ----- | ---- |
|     |     | i,j | i−1,j | i,j−1 |     | i−1,j−1 |     | i−1,j−1 | i−1,j i,j−1 | i,j   |      |
4
for every i,j. The treatment of source terms can be implemented through the compact
version of this formula, which involves only local source values, that does not use recursion,
and can be found in section 4.3. More details about the treatment of source terms for the
| shallow        | water | equations |        | 2.3, are | given | in section |     | 4.6. |     |     |     |
| -------------- | ----- | --------- | ------ | -------- | ----- | ---------- | --- | ---- | --- | --- | --- |
| 4.2. Numerical |       |           | corner | fluxes   |       |            |     |      |     |     |     |
The conservative formulation obtained using global fluxes (38) requires the definition
of numerical corner fluxes to update cell averages. This definition is better achieved by
10

considering the evolution of all the cells neighboring a given node. To this end it is preferable
to recast (38) as
d
F
(cid:99) i
(
+
i,j
1
)
,j+1
+F
(cid:99) i
(
−
i,j
1
)
,j+1
+F
(cid:99) i
(
+
i,j
1
)
,j−1
+F
(cid:99) i
(
−
i,j
1
)
,j−1
q¯ + 2 2 2 2 2 2 2 2 = 0. (42)
i,j
dt ∆x∆y
where we have added the superscript (i,j) to account for the fact that the signs used in the
formula involve the flux balance for q¯ . They thus are interpreted as depending on the
i,j
orientation of the corner normal for the four cells (i+ℓ,j +r) for ℓ,r ∈ {0,1} with respect
to the corner (i+ 1,j + 1), defined by
2 2
(cid:18) (cid:19)
(−1)ℓ+1
n(i+ℓ,j+r) := , i.e., (43)
i+1,j+1 (−1)r+1
2 2
(cid:18) (cid:19) (cid:18) (cid:19) (cid:18) (cid:19) (cid:18) (cid:19)
−1 1 −1 1
n(i,j) := ,n(i+1,j) := , n(i,j+1) := , n(i+1,j+1) := . (44)
i+1,j+1 −1 i+1,j+1 −1 i+1,j+1 1 i+1,j+1 1
2 2 2 2 2 2 2 2
Thus nc is a normal at corner p pointing into cell c. Corner normals and corner fluxes
p
alongside a modified concept of conservation associated to corners rather than edges are used
e.g. in [12] for general polygonal grids, where they are shown to be crucial for structure
preservation. There, corner normals are defined as the average of the two edge-normals
adjacent to the node, scaled with the respective edge lengths. We believe that the definition
above is sufficient in the context of Cartesian meshes.
We further define the scalar n(i+ℓ,j+r) as
i+1,j+1
2 2
n(i+ℓ,j+r) = n(n(i+ℓ,j+r)) = (−1)ℓ+1(−1)r+1 = (−1)ℓ+r, (45)
i+1,j+1 i+1,j+1
2 2 2 2
where n(n) := n n is the product of the two components.
x y
The reinterpreted definition of the corner fluxes above requires a different setting com-
paredtotheclassicalone, andappearsinthecontextofgenuinelymulti-dimensionalRiemann
solvers using more than two states as input (see e.g. [8, 30, 31, 1] and references therein).
Next, we aim at defining the notion of the numerical global flux in the multi-dimensional
context as generally as possible. Then, we elucidate the conditions imposed on the functional
form of the numerical flux by consistency, conservation, and preservation of steady states.
We define the numerical corner fluxes at the corner (i + 1,j + 1) with respect to the four
2 2
cells (i+ℓ,j +r) with ℓ,r ∈ {0,1} as
F (cid:99) (i+ℓ,j+r) = F (cid:99)(F ,F ,F ,F ;q¯ ,q¯ ,q¯ ,q¯ |n(i+ℓ,j+r)).
i+1,j+1 i,j i,j+1 i+1,j i+1,j+1 i,j i,j+1 i+1,j i+1,j+1 i+1,j+1
2 2 2 2
We can now formulate local consistency as
F (cid:99)(F,F,F,F;q,q,q,q|n) = F n(n). (46)
A stronger property is the steady state preservation requirement which can be expressed as
F (cid:99)(F,F,F,F;q¯ ,q¯ ,q¯ ,q¯ |n(i,j) ) = F n(i,j) . (47)
i,j i,j+1 i+1,j i+1,j+1 i+1,j+1 i+1,j+1
2 2 2 2
11

Remark 3 (Steady state subspace). Following [13], it will be shown below that steady state
preservation may be actually already proven if
F∗ = F +α +β (48)
i,j i j
for any two data distributions α and β, such that α only depends on i and β only on j. In
this case, a conservation property more general of (47) reads
F (cid:99) (i,j) (F∗ ,F∗ ,F∗ ,F∗ ;q¯ ,q¯ ,q¯ ,q¯ |n(i,j) ) = F n(i,j) .
i+1,j+1 i,j i,j+1 i+1,j i+1,j+1 i,j i,j+1 i+1,j i+1,j+1 i+1,j+1 i+1,j+1
2 2 2 2 2 2
(49)
This condition will be used also in the consistency analysis of Section 4.5.
Conservation cannot be expressed by face in this framework, as in standard finite volume
methods. It is instead formulated at corners as follows:
F (cid:99) (i,j) +F (cid:99) (i+1,j) +F (cid:99) (i,j+1) +F (cid:99) (i+1,j+1) = 0. (50)
i+1,j+1 i+1,j+1 i+1,j+1 i+1,j+1
2 2 2 2 2 2 2 2
Having established the necessary conditions that the numerical global fluxes have to
satisfy, we next propose a particular choice. As in [31], we define numerical global fluxes as
the sum of a consistent central flux plus a diffusion term D. In the present paper, corner
fluxes are obtained extending to quadrilaterals the multi-dimensional Osher-Solomon flux
proposed in [31], and combining it with the recent work [13] on stationarity preserving,
global flux quadrature SUPG stabilization. We define the numerical corner flux around the
corner (x ,y ), which involves the cells C , C , C and C (see figure 1).
i+1 j+1 i,j i+1,j i,j+1 i+1,j+1
2 2
The same principle is applied to the other corners. We set
F (cid:99) i ( + i,j 2 1 + ,j 1 + ) 1 2 = F i+1 2 ,j+1 2 n( i+ i,j 1 2 + ,j 1 + ) 1 2 +D i ( + i,j 1 2 + ,j 1 + ) 1 2 , F (cid:99) i ( + i+ 2 1 1 ,j ,j + + 1 2 1) = F i+1 2 ,j+1 2 n( i+ i+ 2 1 1 , , j j + + 1 2 1) +D i ( + i+ 1 2 1 , , j j + + 1 2 1),
F (cid:99) i ( + i,j 1 2 ) ,j+1 2 = F i+1 2 ,j+1 2 n( i+ i,j 1 2 ) ,j+1 2 +D i ( + i,j 1 2 ) ,j+1 2 , F (cid:99) i ( + i+ 2 1 1 ,j ,j + ) 1 2 = F i+1 2 ,j+ 2 1 n( i+ i+ 1 2 1 , , j j + ) 1 2 +D i ( + i+ 1 2 1 , , j j + ) 1 2 ,
(51)
where
1
F = (F +F +F +F )
i+1
2
,j+1
2 4
i,j i+1,j i+1,j+1 i,j+1
is the average of the global fluxes at the corner.
To define the numerical dissipation we consider corner dual cells, as depicted on the right
on figure 1. The conservation condition (50) requires that
D(i,j) +D(i+1,j) +D(i,j+1) +D(i+1,j+1) = 0.
i+1,j+1 i+1,j+1 i+1,j+1 i+1,j+1
2 2 2 2 2 2 2 2
To define the corner dissipation, we cannot proceed as in [31], since this would break the
stationarity preserving property. Instead, we take inspiration from the streamline upwind
stabilization (SUPG), studied in the global flux context in [13]. To this end, on the dual cell
C(cid:101)
i+1,j+1
we compute SUPG stabilizing terms
2 2
D i ( + i+ 1 2 ℓ , , j j + +r 1 2 ) :=D(F (cid:102) i+1 2 ,j+1 2 ,q¯ i+1 2 ,j+1 2 |n( i+ i+ 1 2 ℓ , , j j + +r 1 2 ))
(cid:90) (cid:18) 1 1 (cid:19) (52)
=α∆ ∆x Jx∂ ξ φ ℓ,r + ∆y Jy∂ η φ ℓ,r ∂ ξη F (cid:102) i+1 2 ,j+ 2 1 dξdη,
C(cid:101) i+2 1,j+1 2
12

whereF (cid:102) i+1,j+1 isabi-linearQ1 reconstructionoftheglobalfluxonthedualcellfromthefour
adjacent va 2 lues 2 F , F , F , F . Here, Jx and Jy are the Jacobians of the fluxes
i,j i+1,j i,j+1 i+1,j+1 √
f and g computed in the average value q¯ i+1 2 ,j+1 2 = qi,j+qi+1,j+qi 4 ,j+1+qi+1,j+1. ∆ = ∆x √ 2 2 +∆y2 is
the characteristic mesh size, and α = 1/λ with λ the maximal spectral radius of the flux
m m
Jacobians computed with the average state of the four reconstructed values at the corner.
The φ for ℓ,r ∈ {0,1} in the above definition are the standard bi-linear finite element basis
ℓ,r
functions on the quadrilateral C(cid:101) defined by
1
φ (ξ,η) = (1+(−1)ℓ+1ξ)(1+(−1)r+1η), (53)
ℓ,r
4
i.e.,
1 1
φ (ξ,η) = (1−ξ)(1−η), φ (ξ,η) = (1+ξ)(1−η)
0,0 1,0
4 4
(54)
1 1
φ (ξ,η) = (1−ξ)(1+η), φ (ξ,η) = (1+ξ)(1+η)
0,1 1,1
4 4
on the reference element ξ,η ∈ [−1,1]. With this, we can explicitly evaluate the streamline
upwind dissipation terms as
(cid:18) (cid:19)
α∆ n n
D(F (cid:102),q¯|n) = x Jx + y Jy Φ(F (cid:102)), (55)
4 ∆x ∆y
i.e.,
α∆ (cid:18) Jx Jy (cid:19) α∆ (cid:18) Jx Jy (cid:19)
D(i,j+1) = − + Φ , D(i+1,j+1) = + + Φ ,
i+ 2 1,j+1 2 4 ∆x ∆y i+1 2 ,j+1 2 i+1 2 ,j+ 2 1 4 ∆x ∆y i+1 2 ,j+1 2
α∆ (cid:18) Jx Jy (cid:19) α∆ (cid:18) Jx Jy (cid:19)
D(i,j) = − − Φ , D(i+1,j) = + − Φ ,
i+ 2 1,j+1 2 4 ∆x ∆y i+1 2 ,j+1 2 i+1 2 ,j+1 2 4 ∆x ∆y i+1 2 ,j+1 2
(56)
and where Φ is the dual cell residual
i+1,j+1
2 2
(cid:90)
Φ i+1,j+1 := Φ(F (cid:102) i+1,j+1 ) := F i+1,j+1 −F i,j+1 −F i+1,j +F i,j = ∂ xy F (cid:102) i+1,j+1 dxdy .
2 2 2 2 2 2
C(cid:101) i+1 2,j+1
2
(57)
The next sections are devoted to the analysis of some properties of the scheme obtained
with the above definitions, as well as some enhancements.
13

4.3. Compactness of the method
Taking into account only the central flux, without the diffusion, one obtains
d
F
(cid:99) i
(
+
i,j
1
)
,j+1
+F
(cid:99) i
(
−
i,j
1
)
,j+1
+F
(cid:99) i
(
+
i,j
1
)
,j−1
+F
(cid:99) i
(
−
i,j
1
)
,j−1
q¯ = − 2 2 2 2 2 2 2 2 (58a)
i,j
dt ∆x∆y
F n(i,j) +F n(i,j) +F n(i,j) +F n(i,j)
= − i+1 2 ,j+1 2 i+1 2 ,j+ 2 1 i−1 2 ,j+ 2 1 i−1 2 ,j+1 2 i+1 2 ,j−1 2 i+1 2 ,j−1 2 i−1 2 ,j− 2 1 i−1 2 ,j−1 2
∆x∆y
(58b)
F −F −F +F
i+1,j+1 i−1,j+1 i+1,j−1 i−1,j−1
= − 2 2 2 2 2 2 2 2. (58c)
∆x∆y
Define
1
F = (F +F +F +F ) (59)
i+1
2
,j+1
2 4
i,j i+1,j i,j+1 i+1,j+1
1
G = (G +G +G +G ) (60)
i+1
2
,j+1
2 4
i,j i+1,j i,j+1 i+1,j+1
1
R = (R +R +R +R ) (61)
i+1
2
,j+1
2 4
i,j i+1,j i,j+1 i+1,j+1
such that F = F + G + R , and use the recursions (39)–(40) to
i+1,j+1 i+1,j+1 i+1,j+1 i+1,j+1
2 2 2 2 2 2 2 2
obtain
∆y ∆y
F −F = (f +2f +f )+ (f +2f +f ). (62)
i+ 2 1,j+1 2 i+1 2 ,j−1 2 8 i+1,j+1 i+1,j i+1,j−1 8 i,j+1 i,j i,j−1
We find an analogous formula for F −F . The element source contribution is
i−1,j+1 i−1,j−1
2 2 2 2
given from recursion (41) and it reads:
R −R −R +R =
i+1,j+1 i−1,j+1 i+1,j−1 i−1,j−1
2 2 2 2 2 2 2 2
∆x∆y
(s +2s +s +2s +4s +2s +s +2s +s ),
i+1,j+1 i+1,j i+1,j−1 i,j+1 i,j i,j−1 i−1,j+1 i−1,j i−1,j−1
16
(63)
such that in the end
d 1 1
(cid:10) (cid:11)
q¯ = − ⟨[[f ]] ⟩ − [[⟨g ⟩ ]] + ⟨s ⟩ (64)
dt i,j ∆x ·,· i j ∆y ·,· i j ·,· i j
with the average and jump operators on the cells defined by
1
⟨a ⟩ := (a +2a +a ), (65)
i,· j i,j+1 i,j i,j−1
4
1
[[a ]] := (a −a ). (66)
i,· j i,j+1 i,j−1
2
One observes that all the global fluxes drop out and the central part of the method is local.
14

Next, we turn to the numerical stabilization. Defining
ΦF = F −F −F +F , (67)
i+1,j+1 i+1,j+1 i,j+1 i+1,j i,j
2 2
ΦG = G −G −G +G , (68)
i+1,j+1 i+1,j+1 i,j+1 i+1,j i,j
2 2
ΦR = R −R −R +R (69)
i+1,j+1 i+1,j+1 i,j+1 i+1,j i,j
2 2
suchthatΦ = ΦF +ΦG +ΦR , andusingagaintherecursions(39)–(40),
i+1,j+1 i+1,j+1 i+1,j+1 i+1,j+1
2 2 2 2 2 2 2 2
one obtains for every i,j
∆y
ΦF = (f +f −f −f ). (70)
i+1,j+1 2 i+1,j+1 i+1,j i,j+1 i,j
2 2
Analogously,
∆x
ΦG = (g +g −g −g ). (71)
i+1,j+1 2 i+1,j+1 i,j+1 i+1,j i,j
2 2
and
∆x∆y
ΦR = (s +s +s +s ). (72)
i+1,j+1 4 i+1,j i+1,j+1 i,j i,j+1
2 2
One observes again that all the global fluxes drop out and this shows that the method is
completely compact/local with a stencil 3×3.
The Jacobians involved in the update of q are evaluated at the four corners, such that
i,j
the method with only the numerical stabilization reads
D(i,j) +D(i,j) +D(i,j) +D(i,j)
d i+1,j+1 i−1,j+1 i+1,j−1 i−1,j−1
q¯ = − 2 2 2 2 2 2 2 2
i,j
dt ∆x∆y
(cid:34)(cid:32) Jx Jy (cid:33) (cid:32) Jx Jy (cid:33)
1 α∆ i+1,j+1 i+1,j+1 i−1,j+1 i−1,j+1
= − − 2 2 − 2 2 Φ + 2 2 − 2 2 Φ
∆x∆y 4 ∆x ∆y i+1 2 ,j+ 2 1 ∆x ∆y i−1 2 ,j+1 2
(cid:32) Jx Jy (cid:33) (cid:32) Jx Jy (cid:33) (cid:35)
i+1,j−1 i+1,j−1 i−1,j−1 i−1,j−1
+ − 2 2 + 2 2 Φ + 2 2 + 2 2 Φ .
∆x ∆y i+ 2 1,j−1 2 ∆x ∆y i−1 2 ,j−1 2
(73)
We can conclude that the scheme has a local character, as the central part is local after
the combination of the four corner fluxes, while the diffusion part is local at each corner
residual. This leads to a scheme with a compact 3×3 stencil.
4.3.1. Explicit characterization of the method in the linear case
To give the spirit of the method, assume for the moment that the Jacobians are evaluated
on same state (or that we deal with a linear problem). Then, the numerical stabilization
15

becomes
d 1 α∆ (cid:20) Jx (cid:16) (cid:17)
q¯ = − −Φ +Φ −Φ +Φ
dt i,j ∆x∆y 4 ∆x i+1 2 ,j+1 2 i−1 2 ,j+ 2 1 i+1 2 ,j−1 2 i−1 2 ,j−1 2
(74)
Jy (cid:16) (cid:17) (cid:21)
+ −Φ −Φ +Φ +Φ
∆y
i+1
2
,j+1
2
i−1
2
,j+1
2
i+1
2
,j−1
2
i−1
2
,j−1
2
α∆ (cid:20) Jx Jx Jx
= ⟨f −2f +f ⟩ + [[[[g ]] ]] + ⟨[[s ]] ⟩
4 ∆x2 i+1,· i,· i−1,· j ∆x∆y ·,· i j ∆x ·,· i j
(75)
Jy Jy Jy (cid:21)
+ [[[[f ]] ]] + ⟨g −2g +g ⟩ + [[⟨s ⟩ ]] .
∆x∆y ·,· i j ∆y2 ·,j+1 ·,j ·,j−1 i ∆y ·,· i j
Here, the average ⟨·⟩ and jump [[·]] operators introduced in (65) and (66) have been used
again.
Finally, the method can be expressed in classical flux form
ˆ ˆ
d f i+1,j −f i−1,j gˆ i,j+1 −gˆ i,j−1 (cid:10) (cid:11)
q + 2 2 + 2 2 = ⟨s ⟩ . (76)
dt i,j ∆x ∆y ·,· i j
This demonstrates that additionally to the notion (50), the method is also conservative in
the classical sense. The numerical flux through the edge (i+ 1,j) reads
2
Jx Φ +Jx Φ
f ˆ = 1 ⟨f +f ⟩ − α∆ i+1 2 ,j+ 2 1 i+1 2 ,j+1 2 i+1 2 ,j−1 2 i+1 2 ,j− 2 1 . (77)
i+1 2 ,j 2 i+1,· i,· j 2 2∆x∆y
In a quasi-1D situation, i.e. when nothing depends on j and when g = 0, the flux is
(cid:18) (cid:19)
1 α ∆x
f ˆ = (f +f )− Jx f −f − (s +s ) . (78)
i+1 2 2 i+1 i 2 i+1 2 i+1 i 2 i+1 i
With this, next we discuss the interplay between the numerical stabilization and stationarity
preservation.
4.4. Analysis of the method for the linear acoustic system
In this section, we focus on the analysis of the new numerical method by analyzing the
numerical diffusion that allows to achieve stationarity preservation, and obtaining an energy
estimate.
4.4.1. Numerical diffusion and stationarity preservation
We start by considering the linear acoustic system, but similar results can be shown for
nonlinear problems. A classical dimensionally split finite volume scheme with a local Lax-
Friedrichs numerical flux provides the following discretization of the linear acoustic system:
d λ ∆x λ ∆y
m m
p+D u+D v = D p+ D p,
x y xx yy
dt 2 2
d λ ∆x
u+D p = m D u, (79)
x xx
dt 2
d λ ∆y
m
v +D p = D v,
y yy
dt 2
16

j +1 −1 0 1 j +1 1 −1 1 j +1 −1 0 1
8 8 4 2 4 4 4
j −1 0 1 j 1 −1 1 j 0 0 0
4 4 2 2
j −1 −1 0 1 j −1 1 −1 1 j −1 1 0 −1
8 8 4 2 4 4 4
i−1 i i+1 i−1 i i+1 i−1 i i+1
(a) ∆xD¯ (b) ∆x2D¯ (c) ∆x∆yD
x xx xy
Figure 2: Finite difference-like stencils for global flux differential operators.
where D represent the discrete derivative operators given by
q −q q −2q +q
i+1,j i−1,j i+1,j i,j i−1,j
(D q) = , (D q) = , (80)
x i,j 2∆x xx i,j ∆x2
and similarly for D and D . As can be noticed, the stencil used in this discretization is a
y yy
simple 5-points stencil.
Contrary to this, the stencils involved in the new first order global flux method, equipped
with SUPG corner fluxes as described above, includes the cell itself and its eight neighbors
(9-points stencil):
d α∆
¯ ¯ ¯ ¯
p+D u+D v = (D p+D p),
x y xx yy
dt 2
d α∆
u+D ¯ p = (cid:0) D ¯ u+D v (cid:1) , (81)
x xx xy
dt 2
d α∆
¯ (cid:0) ¯ (cid:1)
v +D p = D u+D v ,
y xy yy
dt 2
¯ ¯
where the new finite difference operators are given by the stencils in figure 2. D and D are
x xx
standard discrete first and second order derivatives in x, but including a particular averaging
in y direction, introduced in (65) as ⟨a ⟩ := 1(a +2a +a ). These operators have
i,· j 4 i,j+1 i,j i,j−1
first appeared in [46] and then in virtually all subsequent works on stationarity and vorticity
preservation for linear acoustics on Cartesian grids, e.g. in [46, 52, 34, 45, 40, 10].
In [13] these finite difference operators appeared naturally as Kronecker products of uni-
¯
direction operators: D = D ⊗D I , D = D ⊗D I etc., with matrices I responsible
xy x y y xx xx y y y
for the integration and D I being the particular averaging matrix corresponding to ⟨·⟩.
y y
One observes that the diffusion operators for the velocity no longer depend on second
derivatives of individual components, but instead on the gradient of the divergence operator.
Although this characteristic of the scheme is more readily visible for a simplified model like
the linear acoustic system (81), similar considerations can be drawn also for more complex
nonlinear systems, as is obvious from (75).
Remark 4 (Numerical diffusion for nonlinear problems). For the shallow water equations,
the first order global flux method with SUPG corner fluxes leads to the following discrete
17

evolution equations:
d α∆ α∆
¯ ¯ ¯ ¯
h+D f +D g = (D f +D g )+ (D f +D g ),
x h y h xx hu xy hu xy hv yy hv
dt 2 2
d α∆ α∆
hu+D ¯ f +D ¯ g = (gh ¯ −u¯2)(D ¯ f +D g )− u¯v¯(D f +D ¯ g )+
x hu y hu xx h xy h xy h yy h
dt 2 2
α∆ α∆
+α∆u¯(D ¯ f +D g )+ v¯(D f +D ¯ g )+ u¯(D f +D ¯ g ), (82)
xx hu xy hu xy hu yy hu xy hv yy hv
2 2
d α∆ α∆
hv +D ¯ f +D ¯ g = (gh ¯ −v¯2)(D f +D ¯ g )− u¯v¯(D ¯ f +D g )+
x hv y hv xy h yy h xx h xy h
dt 2 2
α∆ α∆
¯ ¯ ¯
+α∆v¯(D f +D g )+ v¯(D f +D g )+ u¯(D f +D g ),
xy hv yy hv xx hu xy hu xx hv xy hv
2 2
¯ ¯ ¯
where we considered constant Jacobians defined in an average state q¯= (h,hu¯,hv¯) to regroup
the terms in a compact form. For notational convenience, we have introduced the terms f
hu
and g to denote the fluxes for the momentum equation in hu, and similarly for the other
hu
equations. Again, we obtain diffusion terms that depends on the gradient of the divergence
operator, which is essential for stationarity preservation.
4.4.2. Semi-discrete energy stability
In this section, we focus on the semi-discrete energy estimates for the linear acoustic
system. In particular, in the continuous setting, it can be easily proven that the conserved
energy of the system is E = u2+v2 + p2 , by multiplying (15) by qT, summing the three
2 2
equations and integrating over the whole domain Ω:
(cid:90) (cid:90)
[qT∂ q +qTJ∇q]dx = [u∂ u+u∂ p+v∂ v +v∂ p+p∂ p+p∂ u+p∂ v]dx
t t x t y t x y
Ω Ω
d (cid:90) (cid:20) u2 +v2 p2(cid:21) (cid:90)
= + dx+ pv·ndS (83)
dt 2 2
Ω ∂Ω
where the second term is zero for periodic boundary conditions. In our discrete framework,
we would like to prove that
d (cid:90) (cid:20) u2 +v2 p2(cid:21)
+ dx ≤ 0.
dt 2 2
Ω
To do that, we will write the new differential operators introduced above in a tensor
product form to split the contributions from the two dimensions, for more details see [13],
¯ ¯
D = (D M )⊗(M M ), D = (D D )⊗(M M ),
x + − + − xx + − + −
¯
D = (M M )⊗(D D ), D = (D M )⊗(D M )
yy + − + − xy + − + −
where the derivative, D, and average, M, operators are defined as
   
−1 1 0 ... ... 1 1 0 ... ...
2 2
 0 −1 1 ... ...  0 1 1 ... ...
   2 2 
D =  ... ... ... , M =  ... ... ... , (84)
+   +  
 ... ... 0 −1 1   ... ... 0 1 1 
   
2 2
1 ... ... 0 −1 1 ... ... 0 1
2 2
18

with periodic boundary conditions, and by
D = −DT =: D and M = MT =: M.
− + − +
Proposition 5 (Semi-discrete energy inequality). The following semi-discrete energy in-
equality holds,
d (cid:88)
E ≤ 0. (85)
i,j
dt
i,j
Proof. The central part of the method preserves the energy since, e.g. uTD ¯ p + pTD ¯ u =
x x
uTD ¯ p + uTD ¯Tp = 0 up to boundary terms. The evolution of the energy of the system is
x x
thus entirely given by the numerical stabilization as follows:
2 d (cid:88)
E = pTD ¯ p+pTD ¯ p+uTD ¯ u+uTD v +vTD u+vTD ¯ v, (86)
i,j xx yy xx xy xy yy
α∆dt
i,j
where the terms on the right-hand side can be recast as
pTD ¯ p = pT (D D )⊗(M M )p = −∥(D⊗M)p∥2 ≤ 0,
xx + − + −
pTD ¯ p = pT (M M )⊗(D D )p = −∥(M ⊗D)p∥2 ≤ 0,
yy + − + −
uTD ¯ u = uT (D D )⊗(M M )u = −[(D⊗M)u]T [(D⊗M)u],
xx + − + −
uTD v = uT (D M )⊗(M D )v = −[(D⊗M)u]T [(M ⊗D)v],
xy + + − −
vTD u = vT (M D )⊗(D M )u = −[(M ⊗D)v]T [(D⊗M)u],
xy + + − −
vTD ¯ v = vT (M M )⊗(D D )v = −[(M ⊗D)v]T [(M ⊗D)v],
yy + − + −
where the mixed operator was manipulated thanks to MD = DM. Hence, the semi-discrete
energy is found to decrease:
2 d (cid:88)
E ≤−[(D⊗M)u]T [(D⊗M)u+(M ⊗D)v]
i,j
α∆dt
i,j
−[(M ⊗D)v]T [(D⊗M)u+(M ⊗D)v]
=−∥(D⊗M)u+(M ⊗D)v∥2 ≤ 0.
(cid:90) (cid:90)
This is a discrete version of v·∇(∇·v)dx = − (∇·v)2dx+boundary terms.
Ω Ω
4.5. Analogy with Residual Distribution and discrete steady states
The recent work of [31] has provided a general analysis of the relations between multi-
dimensional finite volume methods with point fluxes and residual distribution schemes.
Earlier, it has been shown in [3] that residual distribution schemes can be reformulated
in terms of a global flux finite volume method. This section elaborates on these aspects for
the global flux finite volume approach proposed here. This allows us to give more details on
the discrete steady states of the method in a more general setting.
19

Following the last reference we start from the conservation condition (50) at each corner.
Consider, instead of (51), the following ansatz for the numerical global flux, given by the
trace of the cell (global) flux plus a fluctuation:
F (cid:99) (i+ℓ,j+m) = F n(i+ℓ,j+m) +Φ(i+ℓ,j+m), ℓ,m ∈ {0,1}. (87)
i+1,j+1 i,j i+1,j+1 i+1,j+1
2 2 2 2 2 2
All the properties of the numerical flux can be translated into requirements on the fluc-
tuations Φ(i+ℓ,j+m). The most interesting ones are related to conservation and stationarity
i+1,j+1
2 2
preservation. Corner conservation is written by using the above ansatz in (50), which leads
to the requirement
(cid:88)
Φ(i,j) +Φ(i+1,j) +Φ(i+1,j+1) +Φ(i,j+1) = − F n(i+ℓ,j+m) = −Φ ,
i+1 2 ,j+1 2 i+1 2 ,j+1 2 i+ 2 1,j+ 2 1 i+1 2 ,j+ 2 1 i+ℓ,j+m i+1 2 ,j+1 2 i+1 2 ,j+1 2
ℓ,m∈{0,1}
where Φ is the global flux integral on the corner dual cell as defined in (57). By virtue
i+1,j+1
2 2
of (87), defining a corner flux is thus equivalent to defining a residual distribution scheme
satisfying
(cid:88)
Φ(i+ℓ,j+m) = −Φ . (88)
i+1
2
,j+1
2
i+1
2
,j+
2
1
ℓ,m∈{0,1}
This analogy goes much further, and it is in fact a full equivalence. In particular, we can
prove the following facts.
Proposition 6 (Equivalence with RD). Consider the multi-dimensional global flux finite
volume method (42), with numerical fluxes written in terms of fluctuations (57). Then,
1. the multidimensional finite volume global flux method (42) with piecewise constant data
is equivalent to the Residual Distribution scheme
d
∆x∆y q¯ +Φ(i,j) +Φ(i,j) +Φ(i,j) +Φ(i,j) = 0, (89)
dt i,j i+1,j+1 i+1,j−1 i−1,j+1 i−1,j−1
2 2 2 2 2 2 2 2
with fluctuations Φ(i,j) verifying the conservation condition (88) at each corner (i±
i±1,j±1
2 2
1,j ± 1);
2 2
2. thefinitevolumeglobalfluxmethod (42)withaveragefluxF (cid:99) i ( + i+ 1 2 ℓ , , j j + + 1 2 m) = F i+1 2 ,j+1 2 n( i+ i+ 1 2 ℓ , , j j + +m 1 2 ),
∀ℓ,m ∈ {0,1} is equivalent to the residual distribution scheme with fluctuations
Φ(i+ℓ,j+m) = (F −F )n(i+ℓ,j+m)
i+1
2
,j+1
2
i+1
2
,j+1
2
i+ℓ,j+m i+1
2
,j+1
2
3. the finite volume method including the numerical dissipation in (51), defined by the
streamline upwind terms (56), is equivalent to the residual distribution scheme defined
by
1
Φ(i+ℓ,j+m) = (F −F )n(i+ℓ,j+m) − δ(i+ℓ,j+m)Φ , ℓ,m ∈ {0,1}
i+1 2 ,j+1 2 i+1 2 ,j+1 2 i+ℓ,j+m i+ 2 1,j+1 2 4 i+1 2 ,j+1 2 i+1 2 ,j+ 2 1
(cid:16) (cid:17)
with δ(i+ℓ,j+m) = α∆ (−1)ℓJx +(−1)mJy according to (56);
i+1,j+1 ∆x ∆y
2 2
20

4. both the centered and the stabilized method are steady state preserving with respect to
global fluxes of the type (48). In particular, they both admit discrete steady solutions
verifying
Φ = 0 ∀i,j ;
i+1,j+1
2 2
5. both the centered and the stabilized method are formally second order accurate at steady
state for smooth enough solutions.
Proof. Thefirstfactisaconsequenceoftheidentityn(i,j) +n(i,j) +n(i,j) +n(i,j) =
i+1,j+1 i−1,j+1 i+1,j−1 i−1,j−1
2 2 2 2 2 2 2 2
0, so when replacing (87) in (42) we obtain immediately (89): Using (89) or (42) to represent
the scheme is absolutely equivalent. The second and third properties are obtained by simply
subtracting from the average flux, and from (51) the cell contributions to obtain the corres-
ponding fluctuations.
Concerning the discrete kernel property, consider first the average flux without extra
dissipation, and compute explicitly the sum of the corner contributions:
F n(i,j) +F n(i,j) +F n(i,j) +F n(i,j)
i+1 2 ,j+ 2 1 i+1 2 ,j+ 2 1 i+ 2 1,j−1 2 i+ 2 1,j−1 2 i−1 2 ,j+1 2 i− 2 1,j+1 2 i−1 2 ,j−1 2 i− 2 1,j−1 2
F −F −F +F
i+1,j+1 i+1,j−1 i−1,j+1 i−1,j−1
= (90)
4
1 (cid:16) (cid:17)
= Φ +Φ +Φ +Φ . (91)
4 i+1 2 ,j+ 2 1 i−1 2 ,j+1 2 i+ 2 1,j−1 2 i−1 2 ,j− 2 1
This shows that Φ = 0 is in the kernel of the average flux method. The same is true
i+1,j+1
2 2
for the stabilized scheme for which we can write after assembly around cell i,j
(cid:20) (cid:21)
(cid:88) 1
F n(i,j) − δ(i,j) Φ
i+ 2 ℓ,j+m 2 i+ 2 ℓ,j+m 2 4 i+ 2 ℓ,j+m 2 i+ 2 ℓ,j+m 2
ℓ,m∈{−1,1}
(92)
(cid:88) 1 (cid:104)(cid:16) (cid:17) (cid:105)
= I−δ(i,j) Φ
4 i+ 2 ℓ,j+m 2 i+ 2 ℓ,j+m 2
ℓ,m∈{−1,1}
So to prove point 4., we just check that Φ = 0 is also a consequence of the steady
i+1,j+1
state preservation condition of Remark 3 in (4 2 9). 2 For F = F +α +β as in (48), in each
i,j i j
dual cell
Φ = F +α +β −F −α −β −F −α −β +F +α +β = 0
i+1,j+1 i+1 j+1 i j+1 i+1 j i j
2 2
and thus F (cid:99) (i+ℓ,j+m) = F. Finally, the second-order consistency at steady state is a con-
i+1,j+1
2 2
sequence of results for residual distribution schemes with bounded distribution coefficients
onlinearandbi-linearelements(see[29,48,2]andreferencestherein). Second-orderaccuracy
atsteadystatethusfollowsfromtheboundednessofthecoefficients1/4, and(I−δ(i+ℓ,j+m))/4
i+1,j+1
2 2
which appear in the equivalent forms of the scheme (91) and (92).
The switch from conservative finite volume to residual distribution methods contains
some nuances on which we would like to comment. In particular, the proof uses two different
21

writings of the average flux scheme. This may be confusing as to which is the proper form to
use. The confusion is originated from the fact that for this simple case, due to cancellation,
the same global assembly can be obtained from several different local contributions, not all
fitting into the same conservation framework. For example, the proof shows that the scheme
d Φ i+1,j+1 Φ i+1,j−1 Φ i−1,j+1 Φ i−1,j−1
∆x∆y q¯ + 2 2 + 2 2 + 2 2 + 2 2 = 0 (93)
i,j
dt 4 4 4 4
is equivalent to the average flux scheme. In light of (89), this may lead to the conclusion,
that the definition Φ(i,j) = Φ /4 is a viable one; which is, however, wrong. Such
i+1 2 ,j+ 2 1 i+1 2 ,j+1 2
a definition is not acceptable here, as it does not satisfy the local conservation constraint
(88), which has a minus on the right hand side. This change of sign is related to the differ-
ence between internal and exterior oriented normals. Scheme (93) can also be shown to be
conservative in the classical cell-vertex residual distribution framework, with the appropriate
conservation constraint. On Cartesian meshes, this mismatch cancels out and one ends up
with the same discretization after assembly. However, on general meshes the change in sign
must be carefully accounted for (see e.g. [31]), and the appropriate local conservation and
consistency conditions respected. In particular, the central RD (93) can indeed be written
in a flux form, but the numerical flux has a much more involved expression than the simple
average flux. The interested reader can refer to the last reference, and to [3] for more details.
4.6. Source modification to preserve solutions at rest
This subsection has the goal of providing a direct way to embed bathymetry source
terms typical of the shallow water system with bottom topography. For simplicity, this part
revolves around this PDE system, but the same approach can be also applied for other cases.
In particular, the goal is to achieve stationarity preservation for motionless equilibria, i.e.
lake at rest preservation, coming from the balance between the hydrodynamic pressure and
bottom topography, which is present in both 1D and 2D configurations:
(cid:40) (cid:40)
h(x)+b(x) ≡ η , h(x,y)+b(x,y) ≡ η ,
0 0
1D: 2D: (94)
u(x) ≡ 0, u(x,y) = v(x,y) ≡ 0.
As also shown in section 3, in a 1D global flux framework ([28]) it has been proposed to
integrate the source terms within the flux derivative, thus obtaining a quasi-conservative
formulation of the PDE:
(cid:40) (cid:40)
∂ h+∂ q = 0, ∂ h+∂ q = 0,
t x x t x x
(cid:16) (cid:17) =⇒ (cid:16) (cid:17)
∂ q +∂ q x 2 +gh2 = −gh∂ b, ∂ q +∂ q x 2 +gh2 + (cid:82)x gh∂ bdξ = 0.
t x x h 2 x t x x h 2 ξ
(95)
However, contrary to classical source terms, the bathymetry term shall be treated differently
given the presence of its derivative. To achieve consistency and well-balancedness for the
high order method, in [28] the integral of the bathymetry source is considered to jump at
each interface. In the current low order framework the approach amounts to
(cid:90) xi (cid:90) xi−1 (cid:90) xi h +h
Rx = gh∂ b dξ = gh∂ bdξ + gh∂ bdξ = Rx +g i i−1 (b −b ), (96)
i ξ ξ ξ i−1 2 i i−1
xi−1
22

where the second integral has been computed using a consistent approximation of ∂ b, while
ξ
| for the | h a simple |     | trapezoidal |     | rule has | been | used. |     |     |     |     |     |     |
| ------- | ---------- | --- | ----------- | --- | -------- | ---- | ----- | --- | --- | --- | --- | --- | --- |
The same approach can also be developed for 2D systems, including the the source term
containing ∂ in the x-flux, and the one containing ∂ in the y-flux. Starting from equation
|     |     | x   |     |     |     |     |     |     | y   |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
(21), we can integrate the source terms present in the momentum equation as follows:

|     |     |     | ∂ h+∂    | (hu)+∂ |          | (hv) | = 0       |     |          |          |     |      |      |
| --- | --- | --- | -------- | ------ | -------- | ---- | --------- | --- | -------- | -------- | --- | ---- | ---- |
|     |     |   | t        | x      | y        |      |           |     |          |          |     |      |      |
|     |     |    |          |        | (cid:16) |      | (cid:82)x |     | (cid:17) |          |     |      |      |
|     |     |     | ∂ (hu)+∂ |        | hu2 +gh  | 2    | +         | gh∂ | b dξ     | +∂ (huv) |     | = 0, |      |
|     |     |     | t        | x      |          |      |           | ξ   |          | y        |     |      | (97) |
2
|     |     |    |        |     |         |     | (cid:16) |      |           |       | (cid:17) |      |     |
| --- | --- | --- | ------ | --- | ------- | --- | -------- | ---- | --------- | ----- | -------- | ---- | --- |
|     |     |    |        |     |         |     | hv2      | +gh2 | (cid:82)y |       |          |      |     |
|     |     | ∂  | (hv)+∂ |     | (huv)+∂ |     |          |      | +         | gh∂ b | dη       | = 0. |     |
|     |     |     | t      | x   |         | y   |          | 2    |           | η     |          |      |     |
Proposition 7 (Lake at rest preservation). The 2D global flux scheme of the system (97)
with the source term quadrature provided in equation (96) is exactly well-balanced for the lake
| at rest | solution | (94). |     |     |     |     |     |     |     |     |     |     |     |
| ------- | -------- | ----- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
Proof. To prove that the family of equilibria (94) is exactly preserved when u = v ≡ 0
and η = h + b ≡ η , Since the mixed terms depend only on the velocity and thus vanish,
0
the two momentum equations can be treated separately for the x and y contribution. In
particular we want to show that given a zero velocity and constant free surface elevation,
f +Rx = f +Rx , ∀j. Without loss of generality, we will show this result only for
| i,j |     | i−1,j |       |     |     |     |     |     |     |     |     |     |     |
| --- | --- | ----- | ----- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
|     | i,j |       | i−1,j |     |     |     |     |     |     |     |     |     |     |
the x direction.
| By  | substitution |     | of the | relevant | quantities, |     |        | we obtain |     |       |     |      |      |
| --- | ------------ | --- | ------ | -------- | ----------- | --- | ------ | --------- | --- | ----- | --- | ---- | ---- |
|     |              |     |        |          |             |     | h2 −h2 |           |     | h +h  |     |      |      |
|     |              | f   | +Rx    | −f       | −Rx         | =   | g i    | i−1       | +g  | i i−1 | (b  | −b ) | (98) |
|     |              | i   |        | i−1      |             |     |        |           |     |       | i   | i−1  |      |
|     |              |     | i      |          | i−1         |     |        | 2         |     | 2     |     |      |      |
h +h
|     |     |     |     |     |     | =   | g i | i−1 | (η −η | )   | = 0, | ∀j, | (99) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | ----- | --- | ---- | --- | ---- |
|     |     |     |     |     |     |     |     |     | i     | i−1 |      |     |      |
2
where the last equality holds when η = η ≡ η , with η = h +b .
|     |     |     |     |     |     | i   | i−1 | 0   |     | i i | i   |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
Remark 8 (Compactness of bathymetry source term). Following the notation used in sec-
tion 4.3, we can derive a compact definition of the method also for a source term defined as
above. In particular, the central discretization will include in the update equation of hu the
i,j
| following | term |     |            |       |       |     |     |     |     |       |     |          |       |
| --------- | ---- | --- | ---------- | ----- | ----- | --- | --- | --- | --- | ----- | --- | -------- | ----- |
|           |      | 1   | (cid:28) h | +h    |       |     |     | h   | +h  |       |     | (cid:29) |       |
|           |      |     |            | i+1,· | i,·   |     |     |     | i,· | i−1,· |     |          |       |
|           |      |     | g          |       | (b    | −b  | )+g |     |     | (b    | −b  | ) ,      | (100) |
|           |      |     |            |       | i+1,· |     | i,· |     |     | i,·   |     | i−1,·    |       |
|           |      | ∆x  |            | 2     |       |     |     |     | 2   |       |     |          |       |
j
which is clearly compact. And similarly, for the equation of hv. On the other hand, in the
ΦR,hu
diffusion part, we will have the following term entering the definition of
i+1,j+1
|     |     |     |     |     |            |       |     |       |     |          |     | 2 2 |       |
| --- | --- | --- | --- | --- | ---------- | ----- | --- | ----- | --- | -------- | --- | --- | ----- |
|     |     |     |     |     | (cid:28) h | +h    |     |       |     | (cid:29) |     |     |       |
|     |     |     |     |     |            | i+1,· | i,· |       |     |          |     |     |       |
|     |     |     |     | ∆y  | g          |       | (b  |       | −b  | )        | ,   |     | (101) |
|     |     |     |     |     |            |       |     | i+1,· | i,· |          |     |     |       |
2
j+1
2
ΦR,hv
where ⟨z ⟩ = z + z and similarly for . This shows that the whole
|     | i+ 1,· | j+1 | i+1 | ,j+1 | i+1 ,j |     |     |     |     | i+1 ,j+1 |     |     |     |
| --- | ------ | --- | --- | ---- | ------ | --- | --- | --- | --- | -------- | --- | --- | --- |
|     | 2      | 2   | 2   |      | 2      |     |     |     |     | 2 2      |     |     |     |
method, also in the lake-at-rest well-balanced version, has a compact stencil of size 3×3.
23

|     |     |     |     |     | F       |     | F     | F       |     |     |     |     |
| --- | --- | --- | --- | --- | ------- | --- | ----- | ------- | --- | --- | --- | --- |
|     |     |     |     |     | N−1,j+1 |     | N,j+1 | N+1,j+1 |     |     |     |     |
i+1,j+1
|     |     |     |     |     |     |       |     | 2 2 |       |     |     |     |
| --- | --- | --- | --- | --- | --- | ----- | --- | --- | ----- | --- | --- | --- |
|     |     |     |     |     | F   |       | F   | F   |       |     |     |     |
|     |     |     |     |     |     | N−1,j | N,j |     | N+1,j |     |     |     |
i+1,j−1
|                 |     |          |     |            |         |            |       | 2 2     |                 |     |          |     |
| --------------- | --- | -------- | --- | ---------- | ------- | ---------- | ----- | ------- | --------------- | --- | -------- | --- |
|                 |     |          |     |            | F       |            | F     | F       |                 |     |          |     |
|                 |     |          |     |            | N−1,j−1 |            | N,j−1 | N+1,j−1 |                 |     |          |     |
|                 |     | Figure   | 3:  | Ghost      | cell    | labelling: | ghost | cells   | are highlighted |     | in blue. |     |
| 4.7. Compatible |     | boundary |     | conditions |         |            |       |         |                 |     |          |     |
Boundary conditions play an essential role in the practical application of numerical meth-
ods. Let us for example consider the usual ghost cell approach. In a classical dimension-by-
dimension finite volume method, homogeneous Neumann boundary conditions on the state
variables can be simply enforced by copying the state. This is somewhat consistent with
the internal treatment based on one dimensional Riemann fluxes using two states. However,
this approach cannot be steady state preserving since it is not based on the global flux. One
should consider corner fluxes on the boundaries, and choose the ghost states consistently with
the equilibrium condition. Here, we construct compatible Neumann boundaries based on the
discrete constraint Φ = 0, which we have shown to be the relevant multi-dimensional
corner
characterization of our discrete steady states. Compatible transmissive conditions require
this relation to be verified by the ghost cells at each boundary corner.
Let us consider the right boundary domain (see figure 3). We can impose F :=
N+1,j−1
F
for all j on the global fluxes, instead of computing the state variables as for classical
N,j−1
homogeneous Neumann conditions. Then, we can use use directly the global flux at the
| boundary | corners. |        | This will | lead      | to    | the following |     | relation |     |       |      |       |
| -------- | -------- | ------ | --------- | --------- | ----- | ------------- | --- | -------- | --- | ----- | ---- | ----- |
|          |          |        |           |           | F     | −F            |     | −F       | +F  |       |      |       |
|          |          |        | Φ         | =         |       |               |     |          |     |       | = 0, | (102) |
|          |          |        | N+1,j−1   |           | N+1,j |               | N,j | N+1,j−1  |     | N,j−1 |      |       |
|          |          |        | 2         | 2         |       |               |     |          |     |       |      |       |
| which    | gives a  | steady | state     | solution. |       |               |     |          |     |       |      |       |
On corner ghosts, similarly, one has to impose F := F , if also on the top
|                   |     |            |     |           |     |          |     | Nx+1,Ny+1 |         |       | Nx,Ny+1 |     |
| ----------------- | --- | ---------- | --- | --------- | --- | -------- | --- | --------- | ------- | ----- | ------- | --- |
| side transmissive |     | conditions |     | must      | be  | imposed, |     | this will | results | in    |         |     |
|                   |     |            | F   |           | =   | F        |     | = F       |         | = F   | .       |     |
|                   |     |            |     | Nx+1,Ny+1 |     | Nx,Ny+1  |     | Nx+1,Ny   |         | Nx,Ny |         |     |
Similarideascanbeusedalsoforotherboundarytypes, butinourexperiencetransmissive
conditions are the most critical to correctly maintain the internal structure, since they in-
volve no external data, which provide some link to the correct solution for other boundary
conditions.
Whenseekingtopreservesteadystates,itiscrucialtoensurethatthenumberofequations
imposed—either by boundary conditions or by the steady state equations—does not exceed
the number of unknowns, or that these equations are mutually compatible.
24

The steady state conditions enforced by Φ for i = 0,...,N and j = 0,...,N
|     |     |     |     |     |     |     |     | i+1,j+1 |     |     | x   |     | y   |
| --- | --- | --- | --- | --- | --- | --- | --- | ------- | --- | --- | --- | --- | --- |
|     |     |     |     |     |     |     |     | 2       | 2   |     |     |     |     |
introduce N (N +1)(N +1) linearly independent constraints on N N N unknowns q¯
|     | eq  | x   |     | y   |     |     |     |     |     |     | eq  | x y | i,j |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
(for i = 1,...,N , j = 1,...,N ). To satisfy these extra constraints, N (2N + 2N +
|     |     | x   |     |     | y   |     |     |     |     |     |     | eq x | y   |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | ---- | --- |
4) ghost cell values are introduced. This leaves N (N + N + 3) equations that can be
|     |     |     |     |     |     |     |     |     | eq x | y   |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | ---- | --- | --- | --- | --- |
specified at the boundaries, typically through Dirichlet conditions. Homogeneous Neumann
conditions, as previously discussed, are compatible with the internal constraints and do not
introduce additional equations. Therefore, the N (N +N +3) remaining degrees of freedom
|     |     |     |     |     |     |     |     | eq  | x   | y   |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
correspond to at most two sides where Dirichlet boundary conditions can be imposed, and if
| not consecutive |     | sides, | they | should |     | not have | all | corners | included. |     |     |     |     |
| --------------- | --- | ------ | ---- | ------ | --- | -------- | --- | ------- | --------- | --- | --- | --- | --- |
In summary, when seeking equilibria, the boundary conditions must be compatible with
the internal constraints. It is only possible, and necessary, to impose complete Dirichlet
conditions on (at most) two boundaries. In this respect, this work is changing the perspective
on this issue. For many years, schemes similar to the one obtained here have been considered
as flawed due to the existence of the steady states characterized by proposition 6. This is
due to the fact that spurious oscillating modes may also satisfy the condition Φ = 0
i+1/2,j+1/2
∀i, j (see e.g. [48, 2] and references therein). However, this is only true if one considers the
problem locally, which is a wrong way to define multidimensional steady states as they must
include boundary conditions. If these are imposed in a compatible manner, spurious modes
can be controlled. This work, as well as the work discussed in [13], contributes to rectifying
this notion.
| 5. Standard |     | finite | volume |     | scheme |     | used | for comparison |     |     |     |     |     |
| ----------- | --- | ------ | ------ | --- | ------ | --- | ---- | -------------- | --- | --- | --- | --- | --- |
In this section, we present the standard finite volume (FV) scheme used for comparison
with the novel global flux (GF) scheme in the numerical experiments presented in section 6.
The classical finite volume formulation for the 2D nonlinear hyperbolic problem (11) can be
| written | by integrating |     | it  | in the | cell | C   | :   |     |     |     |     |     |     |
| ------- | -------------- | --- | --- | ------ | ---- | --- | --- | --- | --- | --- | --- | --- | --- |
i,j
|     |     |     |     |     | ˆ     |     | ˆ     |       |     |       |        |     |       |
| --- | --- | --- | --- | --- | ----- | --- | ----- | ----- | --- | ----- | ------ | --- | ----- |
|     |     |     |     |     | f     | −f  |       | gˆ    | −gˆ |       |        |     |       |
|     |     |     |     | d   | i+1,j |     | i−1,j | i,j+1 |     | i,j−1 |        |     |       |
|     |     |     |     | q¯  | +     | 2   | 2     | +     | 2   | 2     | = s¯ , |     | (103) |
|     |     |     |     | i,j |       |     |       |       |     |       | i,j    |     |       |
|     |     |     |     | dt  |       | ∆x  |       |       | ∆y  |       |        |     |       |
where the numerical flux f ˆ is computed through the local Lax-Friedrichs (or Rusanov)
i+1,j
2
flux:
|     |     |     |     |     | 1 (cid:16) |     |     | (cid:17) | λ (cid:16) |        | (cid:17) |     |       |
| --- | --- | --- | --- | --- | ---------- | --- | --- | -------- | ---------- | ------ | -------- | --- | ----- |
|     |     |     | ˆ   |     |            |     |     |          | m          |        |          |     |       |
|     |     |     | f   | =   | fL         | +fR |     | −        | qR         | −qL    | ,        |     | (104) |
|     |     |     | i+1 | ,j  | i+1        | ,j  | i+1 | ,j       |            | i+ 1,j | i+ 1,j   |     |       |
|     |     |     |     | 2   | 2          | 2   | 2   |          | 2          | 2      | 2        |     |       |
where
|     |     |     |     | fL    | =   | f(qL  | ),  | fR    | =   | f(qR  | ),  |     | (105) |
| --- | --- | --- | --- | ----- | --- | ----- | --- | ----- | --- | ----- | --- | --- | ----- |
|     |     |     |     | i+1,j |     | i+1,j |     | i+1,j |     | i+1,j |     |     |       |
|     |     |     |     |       | 2   |       | 2   |       | 2   | 2     |     |     |       |
and similarly for the others. The source term is computed by integrating the source term
| over the | cell | C : |     |     |     |     |     |     |     |     |     |     |     |
| -------- | ---- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
i,j
|     |     |     |     |     |      | 1   | (cid:90) | x i+1 (cid:90) | y j+1  |     |     |     |       |
| --- | --- | --- | --- | --- | ---- | --- | -------- | -------------- | ------ | --- | --- | --- | ----- |
|     |     |     |     |     |      |     |          | 2              | 2      |     |     |     |       |
|     |     |     |     |     | s¯ = |     |          |                | sdxdy. |     |     |     | (106) |
i,j
∆x∆y
|     |     |     |     |     |     |     | x   | i−1 | y j−1 |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | ----- | --- | --- | --- | --- |
|     |     |     |     |     |     |     |     | 2   | 2     |     |     |     |     |
25

In this work, we are going to compare the first order global flux scheme against the standard
finitevolumewithbothpiece-wiseconstant(firstorderaccurate)andpiece-wiselinear(second
order accurate) reconstructions. Since the solution update can be performed in a dimension-
by-dimensionway, wecanfocusonlyonthex-directionforsimplicity. Forpiece-wiseconstant
reconstruction, the left and right states at interface i+ 1 simply are
2
|     |     |     |     | qL    |     |      | qR  |       |       |     |     |     |       |
| --- | --- | --- | --- | ----- | --- | ---- | --- | ----- | ----- | --- | --- | --- | ----- |
|     |     |     |     |       | =   | q¯ , |     | =     | q¯    | ,   | ∀j. |     | (107) |
|     |     |     |     | i+1,j |     | i,j  |     | i+1,j | i+1,j |     |     |     |       |
|     |     |     |     |       | 2   |      |     | 2     |       |     |     |     |       |
While in the second order case, the left and right states are computed through a piece-wise
| linear reconstruction |      |         | of    | the solution |     | as        |       |       |        |     |       |       |       |
| --------------------- | ---- | ------- | ----- | ------------ | --- | --------- | ----- | ----- | ------ | --- | ----- | ----- | ----- |
|                       |      | q˜(x,y) | =     | q¯ +(x−x     |     | )(∂       | q)    | +(y   | −y )(∂ | q)  | , x,y | ∈ C . | (108) |
|                       |      |         |       | i,j          |     | i         | x i,j |       | j      | y   | i,j   | i,j   |       |
| Hence, the            | left | and     | right | states       | at  | interface | i+    | 1 are |        |     |       |       |       |
2
|     |       |     |      | ∆x  |       |     |       |     |       | ∆x  |         |       |       |
| --- | ----- | --- | ---- | --- | ----- | --- | ----- | --- | ----- | --- | ------- | ----- | ----- |
|     | qL    |     |      |     |       |     | qR    |     |       |     |         |       |       |
|     |       | =   | q¯ + |     | (∂ q) | ,   |       | =   | q¯    | −   | (∂ q)   | , ∀j. | (109) |
|     | i+1,j |     | i,j  | 2   | x     | i,j | i+1,j |     | i+1,j | 2   | x i+1,j |       |       |
|     |       | 2   |      |     |       |     |       | 2   |       |     |         |       |       |
Here, the slopes (∂ q) are evaluated using the generalized minmod limiter [47]:
|     |     |       | x i,j |        |          |       |     |           |     |       |           |          |       |
| --- | --- | ----- | ----- | ------ | -------- | ----- | --- | --------- | --- | ----- | --------- | -------- | ----- |
|     |     |       |       |        | (cid:18) | q¯    | −q¯ | q¯        | −q¯ |       | q¯ −q¯    | (cid:19) |       |
|     |     |       |       |        |          | i+1,j |     | i,j i+1,j |     | i−1,j | i,j i−1,j |          |       |
|     |     | (∂ q) | =     | minmod |          | ϑ     |     | ,         |     |       | ,ϑ        | ,        | (110) |
|     |     | x     | i,j   |        |          |       |     |           |     |       |           |          |       |
|     |     |       |       |        |          |       | ∆x  |           | 2∆x |       | ∆x        |          |       |
whereϑisusedtocontroltheamountofdissipation. Inparticular, thelargerϑis, thesharper
and more oscillatory the reconstruction will be. For the simulations presented in section 6,
when not specified, we set the parameter ϑ = 1.3. The same approach has been used along
the y direction.
| The | classical | minmod |     | function |     | is defined |     | as  |     |     |     |     |     |
| --- | --------- | ------ | --- | -------- | --- | ---------- | --- | --- | --- | --- | --- | --- | --- |

|     |     |     |     |     |     |     | min(a,b,c), |     |     | if a,b,c | > 0, |     |     |
| --- | --- | --- | --- | --- | --- | --- | ----------- | --- | --- | -------- | ---- | --- | --- |


|     |     |     | minmod(a,b,c) |     |     | =   | max(a,b,c), |     |     | if a,b,c | < 0, |     | (111) |
| --- | --- | --- | ------------- | --- | --- | --- | ----------- | --- | --- | -------- | ---- | --- | ----- |

|              |     |             |     |     |     |     | 0, |     |     | otherwise. |     |     |     |
| ------------ | --- | ----------- | --- | --- | --- | --- | --- | --- | --- | ---------- | --- | --- | --- |
| 6. Numerical |     | experiments |     |     |     |     |     |     |     |            |     |     |     |
In this section, the goal is to show the performance of the new global flux scheme (GF)
when compared to the classical finite volume first order (FV-1) and second order (FV-2)
approaches presented in section 5. Several test cases are presented to study the impact
of the method on both linear and nonlinear hyperbolic problems: linear acoustics, shallow
water and the Euler equations. The numerical experiments are performed taking the gravity
g = 9.812, for the shallow water system, and the ratio of specific heats γ = 1.4 for the Euler
system. All convergence analyses are performed on a set of nested quadrilateral meshes with
N = N = 20,40,80,160,320. Time integration is performed through classical explicit Euler
x y
| and second | order | Runge-Kutta |     |     | methods. |     |     |     |     |     |     |     |     |
| ---------- | ----- | ----------- | --- | --- | -------- | --- | --- | --- | --- | --- | --- | --- | --- |
26

|      |     |     |     |     |     |     | 102.5 | FV  |     |     |     |
| ---- | --- | --- | --- | --- | --- | --- | ----- | --- | --- | --- | --- |
| 10−1 |     |     |     |     |     |     |       | FV2 |     |     |     |
GF
u
10−2
ni
y
N 102
rorre
| 10−3 |     |     |     |     |     | =   |     |     |     |     |     |
| ---- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
x
N
∞ 10−4
| L   |     | FV  |     |     |     |     |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
101.5
| 10−5 |     | FV2 |     |     |     |     |     |     |     |     |     |
| ---- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
GF
|     | 10−2 |               | 10−1 | 100 |      | 101 |     | 10−2 10−1     |     | 100 101  |     |
| --- | ---- | ------------- | ---- | --- | ---- | --- | --- | ------------- | --- | -------- | --- |
|     |      | Computational |      |     | time | [s] |     | Computational |     | time [s] |     |
Figure 4: Comparison of computational time vs L error in u and vs grid size N for FV, FV-2, and GF
|     |     |     |     |     |     | ∞   |     |     | x   |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
methods.
| 6.1. Linear |            | acoustic | system |     |     |     |     |     |     |     |     |
| ----------- | ---------- | -------- | ------ | --- | --- | --- | --- | --- | --- | --- | --- |
| 6.1.1.      | Stationary |          | vortex |     |     |     |     |     |     |     |     |
The first test case considered here concerns the simulation of the linear acoustic system.
The initial condition is a compactly supported vortex centered in (x ,y ) = (0.5,0.5) defined
0 0
on the square [0,1]×[0,1] with periodic boundary conditions, which is given by
|     |     |     |     |     | p(x,y) | = 1,                |     |     |     |     |     |
| --- | --- | --- | --- | --- | ------ | ------------------- | --- | --- | --- | --- | --- |
|     |     |     |     |     | u(x,y) | = (y −y )f(ρ(x,y)), |     |     |     |     |     |
0
|     |     |     |     |     | v(x,y) | = −(x−x | )f(ρ(x,y)), |     |     |     |     |
| --- | --- | --- | --- | --- | ------ | ------- | ----------- | --- | --- | --- | --- |
0
√
√
|            |     | (x−x0)2+(y−y0)2 |     |     |         |             |     | γ(1+cos(πρ))2 |      | 12π 0.981     |     |
| ---------- | --- | --------------- | --- | --- | ------- | ----------- | --- | ------------- | ---- | ------------- | --- |
| withρ(x,y) | =   |                 |     |     | ,wherer | = 0.45,f(ρ) | =   |               | andγ | = √           | .   |
|            |     |                 | r0  |     |         | 0           |     |               |      | r0 315π2−2048 |     |
This vortex is taken from the work [50], where its derivation is described. This initial condi-
| tion is | a steady | state | of  | the acoustic |     | system. |     |     |     |     |     |
| ------- | -------- | ----- | --- | ------------ | --- | ------- | --- | --- | --- | --- | --- |
In table 1, the errors computed with the L norm and convergence rates are shown. As
2
canbenoticed, theGFmethodoutperformsthestandardFV-1andFV-2methodsintermsof
discretizationerrors. AlthoughtheGFisinprinciplefirstorderaccurateduetothepiece-wise
constant reconstruction, a superconvergence behavior is experienced for stationary solutions
(compare Proposition 6). Hence, the GF method is not only able to preserve the vortex
structure, but does so at second order accuracy. In figure 4, we compare the computational
time with the L∞ error and with the mesh size of the three methods. The computational
cost of the GF-FV is slightly larger than the one of the classical FV-1 (on average less than a
factor 2), while it is comparable with the compuational costs of FV-2. Still, the error of the
GF is many order of magnitude better than classical schemes. In particular, this is evident
when increasing the final time of the simulation, as shown in figure 5. Classical methods
like FV-1 and FV-2 are not able to preserve the vortex structure for long times, due to their
| numerical | dissipation |     | that | spoils | the | final solution. |     |     |     |     |     |
| --------- | ----------- | --- | ---- | ------ | --- | --------------- | --- | --- | --- | --- | --- |
27

Table 1: Linear acoustic system: vortex (t =1). L error and order of accuracy n˜ for FV-1, FV-2 and GF.
f 2
p u v
N ,N L n˜ L n˜ L n˜
x y 2 2 2
FV-1
20 3.51E-05 – 6.51E-02 – 6.51E-02 –
40 4.58E-05 -0.38 5.42E-02 0.26 5.42E-02 0.26
80 2.71E-05 0.75 3.97E-02 0.44 3.97E-02 0.44
160 1.06E-05 1.35 2.54E-02 0.64 2.54E-02 0.64
320 3.34E-06 1.66 1.47E-02 0.79 1.47E-02 0.79
FV-2
20 4.31E-04 – 2.58E-02 – 2.58E-02 –
40 2.54E-04 0.76 6.12E-03 2.07 6.12E-03 2.07
80 6.37E-05 1.99 1.61E-03 1.93 1.61E-03 1.93
160 1.29E-05 2.30 4.46E-04 1.84 4.46E-04 1.84
320 2.73E-06 2.24 1.25E-04 1.83 1.25E-04 1.83
GF
20 4.72E-05 – 3.95E-04 – 3.95E-04 –
40 4.53E-05 0.05 9.17E-05 2.10 9.17E-05 2.10
80 1.95E-05 1.21 2.26E-05 2.01 2.26E-05 2.01
160 6.42E-06 1.60 5.58E-06 2.01 5.58E-06 2.01
320 1.85E-06 1.79 1.38E-06 2.01 1.38E-06 2.01
1.0
0.8 0.6 0.4
0.2
0.0
0.0 0.2 0.4 0.6 0.8 1.0
x
y
p
1.0
0.8 0.6 0.4
0.2
0.0
0.0 0.2 0.4 0.6 0.8 1.0
x
y
u
1.0
0.8 0.6 0.4
0.2
0.0
0.0 0.2 0.4 0.6 0.8 1.0
x
y
v
1.0
0.8 0.6 0.4
0.2
0.0
0.0 0.2 0.4 0.6 0.8 1.0
x
y
1e 13+1 1e 17 1e 15 norm vel 1e 15
0.921 1.0 1.6 3.00 3.00
0.682 2.25 0.8 1.2 2.64 0.442 0.8 1.50 2.28 0.202 0.6 0.4 0.75 1.92 0.039 0.0 0.00 1.56 0.279 0.4 0.4 0.75 1.20 1.50
0.520 0.2 0.8 0.84 2.25
0.759 1.2 0.48
0.0 3.00
0.999 0.12
0.0 0.2 0.4 0.6 0.8 1.0
x
(a) FV-1
y
p
1.0
0.8 0.6 0.4
0.2
0.0
0.0 0.2 0.4 0.6 0.8 1.0
x
y
u
1.0
0.8 0.6 0.4
0.2
0.0
0.0 0.2 0.4 0.6 0.8 1.0
x
y
v
1.0
0.8 0.6 0.4
0.2
0.0
0.0 0.2 0.4 0.6 0.8 1.0
x
y
1e 5+1 norm vel
0.92 1.0 0.0057 0.004 0.004
0.68 0.003 0.003 0.0051 0.8 0.0045 0.44 0.002 0.002 0.0039 0.20 0.6 0.001 0.001 0.0033 0.04 0.000 0.000 0.4 0.001 0.001 0.0027 0.28 0.0021 0.002 0.002
0.52 0.2 0.003 0.003 0.0015
0.76 0.004 0.004 0.0009
0.0
1.00 0.0003
0.0 0.2 0.4 0.6 0.8 1.0
x
(b) FV-2
y
p
1.0
0.8 0.6 0.4
0.2
0.0
0.0 0.2 0.4 0.6 0.8 1.0
x
y
u
1.0
0.8 0.6 0.4
0.2
0.0
0.0 0.2 0.4 0.6 0.8 1.0
x
y
v
1.0
0.8 0.6 0.4
0.2
0.0
0.0 0.2 0.4 0.6 0.8 1.0
x
y
1e 14+1 norm vel
0.4 0.24 0.24 0.27
0.24 0.2 0.18 0.18 0.21 0.0 0.12 0.12 0.18 0.2 0.06 0.06 0.15 0.00 0.00 0.4 0.06 0.06 0.12 0.6 0.12 0.12 0.09
0.8 0.18 0.18 0.06
1.0 0.24 0.24 0.03
1.2 0.00
(c) GF
Figure 5: Linear acoustic system: vortex. Isocontours of the velocity norm obtained with FV-1, FV-2 and
GF after a long time integration (t =200).
f
28

| 6.2. Euler | equations  |        |     |     |     |     |     |     |     |     |
| ---------- | ---------- | ------ | --- | --- | --- | --- | --- | --- | --- | --- |
| 6.2.1.     | Isentropic | vortex |     |     |     |     |     |     |     |     |
In this section, we test the proposed method on a smooth isentropic vortex [35]. The
initial condition is given in terms of primitive variables and it consists in superposition of a
| homogeneous |     | background |           | flow and | a perturbation: |        |     |      |        |     |
| ----------- | --- | ---------- | --------- | -------- | --------------- | ------ | --- | ---- | ------ | --- |
|             |     |            | (ρ,u,v,p) |          | = (1+δρ,        | u +δu, | v   | +δv, | 1+δp). |     |
|             |     |            |           |          |                 | 0      |     | 0    |        |     |
The test case is set up in a [0,10] × [0,10] domain with periodic boundary conditions and
(cid:112)
vortex radius r = (x−5)2 +(y −5)2. The vortex strength is ϵ = 5, and the entropy
perturbation is assumed to be zero. Given these hypotheses, the perturbations on velocity
| and temperature |          | can      | be written    |                      | as         |          |          |           |          |            |
| --------------- | -------- | -------- | ------------- | -------------------- | ---------- | -------- | -------- | --------- | -------- | ---------- |
|                 | (cid:20) | (cid:21) | (cid:18)      | 1−r2(cid:19)(cid:20) |            | (cid:21) |          |           |          |            |
|                 | δu       |          | ϵ             |                      | −(y        | −5)      |          |           | (γ −1)ϵ2 |            |
|                 |          | =        | exp           |                      |            | ,        | δT       | = −       |          | exp(1−r2). |
|                 | δv       |          | 2π            | 2                    | (x−5)      |          |          |           | 8γπ2     |            |
| It follows      | that     | the      | perturbations |                      | on density | and      | pressure | read      |          |            |
|                 |          |          |               |                      | 1          |          |          |           | γ        |            |
|                 |          |          | δρ =          | (1+δT)γ−1            | −1,        |          | δp =     | (1+δT)γ−1 | −1.      |            |
This test case is a stationary solution of the Euler equations. Note, that the maximum Mach
| number | for | this set | up is | about Ma= | 0.7. |     |     |     |     |     |
| ------ | --- | -------- | ----- | --------- | ---- | --- | --- | --- | --- | --- |
Intable2, theconvergenceanalysisfortheisentropicvortexispresentedbycomparingthe
FV-1, FV-2 and GF methods by running the simulation of a static vortex, i.e. u = v = 0,
0 0
until a final time t = 1. As observed above, the GF shows superconvergent behavior with
f
order 2. In terms of discretization errors, it outperforms not only the classical piecewise
constant finite volume method, but also the second-order approach equipped with a linear
reconstruction. In figure 6, we compare the computational costs of the method with respect
to the L∞ error of ρu and with the mesh size. We confirm also for the nonlinear case that the
GF is slightly more expensive of FV-1, but comparable to FV-2, while the error is increadibly
smaller. Even after very long simulations times (see figure 7) the new GF method is able
to maintain the vortex, while the first order FV-1 dissipates everything away, and FV-2
significantly distorts the vortex structure and still diffuses it more than GF. Observe that
the nonlinearity of the equations makes this test significantly more challenging than the
corresponding test for linear acoustics, where in particular advection is not present.
In table 4, the convergence analysis for a moving isentropic vortex is presented with
u = v = 1 and a final time t = 10. Here, we can directly observe that the GF method
| 0   | 0   |     |     | f   |     |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
is indeed first order accurate, as expected, since the reconstruction is piecewise constant.
No superconvergence is observed in this case, as the solution is not stationary. The results
show an relevant improvement of the GF method over the FV-1 in both discretization errors
and convergence rates, while the FV-2 method is, in this case, the best since it is able to
achieve second order accuracy. In figure 8, the solution at the final time is shown for the
three methods.
| 6.2.2. | Isentropic | vortex: | low | Mach | behaviour |     |     |     |     |     |
| ------ | ---------- | ------- | --- | ---- | --------- | --- | --- | --- | --- | --- |
As recalled in the introduction, the resolution of the long-time limit of linear acoustics
is tightly connected to the low Mach number limit of the Euler equations. This connection
29

Table 2: Euler equations: isentropic vortex with u =v =0 (t =1). L error and order of accuracy n˜ for
|     |     |     |     | 0 0 | f   | 2   |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
FV-1, FV-2 and GF methods.
|      |     | ρ   | ρu  |     | ρv  |     | ρE  |     |
| ---- | --- | --- | --- | --- | --- | --- | --- | --- |
| N ,N | L   | n˜  | L   | n˜  | L   | n˜  | L   | n˜  |
| x y  | 2   |     | 2   |     | 2   |     | 2   |     |
FV-1
| 20  | 3.58E-01 | –    | 6.77E-01 | –    | 6.77E-01 | –    | 1.16E+00 | –    |
| --- | -------- | ---- | -------- | ---- | -------- | ---- | -------- | ---- |
| 40  | 2.47E-01 | 0.53 | 4.40E-01 | 0.62 | 4.40E-01 | 0.62 | 8.29E-01 | 0.48 |
| 80  | 1.49E-01 | 0.72 | 2.59E-01 | 0.76 | 2.59E-01 | 0.76 | 5.15E-01 | 0.68 |
| 160 | 8.33E-02 | 0.84 | 1.43E-01 | 0.85 | 1.43E-01 | 0.85 | 2.91E-01 | 0.82 |
| 320 | 4.42E-02 | 0.91 | 7.56E-02 | 0.91 | 7.56E-02 | 0.91 | 1.56E-01 | 0.90 |
FV-2
| 20  | 1.06E-01 | –    | 2.05E-01 | –    | 2.00E-01 | –    | 4.32E-01 | –    |
| --- | -------- | ---- | -------- | ---- | -------- | ---- | -------- | ---- |
| 40  | 3.62E-02 | 1.55 | 6.74E-02 | 1.60 | 6.71E-02 | 1.57 | 1.20E-01 | 1.85 |
| 80  | 1.07E-02 | 1.76 | 1.93E-02 | 1.80 | 1.95E-02 | 1.78 | 2.91E-02 | 2.04 |
| 160 | 2.39E-03 | 2.16 | 5.58E-03 | 1.78 | 5.61E-03 | 1.79 | 7.04E-03 | 2.04 |
| 320 | 5.12E-04 | 2.22 | 1.39E-03 | 2.00 | 1.39E-03 | 2.01 | 1.56E-03 | 2.17 |
GF
| 20  | 1.52E-02 | –    | 3.67E-02 | –    | 3.67E-02 | –    | 4.59E-02 | –    |
| --- | -------- | ---- | -------- | ---- | -------- | ---- | -------- | ---- |
| 40  | 5.95E-03 | 1.35 | 1.15E-02 | 1.67 | 1.15E-02 | 1.67 | 1.54E-02 | 1.57 |
| 80  | 1.76E-03 | 1.76 | 3.06E-03 | 1.90 | 3.06E-03 | 1.90 | 4.35E-03 | 1.82 |
| 160 | 4.69E-04 | 1.90 | 7.87E-04 | 1.96 | 7.87E-04 | 1.96 | 1.16E-03 | 1.90 |
| 320 | 1.21E-04 | 1.95 | 2.00E-04 | 1.97 | 2.00E-04 | 1.97 | 3.02E-04 | 1.94 |
|     |          |      |          |      | 102.5    | FV   |          |      |
FV2
10−1
GF
uρ
ni y
N 102
10−2
rorre =
x
N
∞ 10−3
L FV
101.5
FV2
GF
10−4
| 10−2          | 10−1 | 100  | 101 |     |     | 10−2          | 10−1 | 100 101  |
| ------------- | ---- | ---- | --- | --- | --- | ------------- | ---- | -------- |
| Computational |      | time | [s] |     |     | Computational |      | time [s] |
Figure 6: Comparison of computational time vs L error in ρu and vs grid size N for FV, FV-2, and GF
|     |     |     |     | ∞   |     |     | x   |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
methods.
30

Table 3: Euler equations: isentropic vortex with u =v =1 (t =10). L error and order of accuracy n˜ for
0 0 f 2
FV-1, FV-2 and GF methods.
ρ ρu ρv ρE
N ,N L n˜ L n˜ L n˜ L n˜
x y 2 2 2 2
FV-1
20 6.50E-01 – 1.54E+00 – 1.54E+00 – 3.12E+00 –
40 6.21E-01 0.06 1.46E+00 0.07 1.46E+00 0.07 3.01E+00 0.05
80 5.82E-01 0.09 1.31E+00 0.15 1.31E+00 0.15 2.83E+00 0.09
160 5.13E-01 0.18 1.06E+00 0.30 1.07E+00 0.29 2.48E+00 0.19
320 4.01E-01 0.35 7.58E-01 0.49 7.63E-01 0.49 1.92E+00 0.36
FV-2
20 5.29E-01 – 1.07E+00 – 1.12E+00 – 2.48E+00 –
40 2.45E-01 1.10 4.42E-01 1.28 4.84E-01 1.21 1.17E+00 1.08
80 6.55E-02 1.90 1.26E-01 1.80 1.37E-01 1.81 2.82E-01 2.04
160 1.85E-02 1.82 3.22E-02 1.97 3.42E-02 2.00 6.03E-02 2.22
320 3.38E-03 2.45 7.49E-03 2.10 8.12E-03 2.07 1.28E-02 2.23
GF
20 5.46E-01 – 1.30E+00 – 1.13E+00 – 2.58E+00 –
40 4.44E-01 0.30 1.02E+00 0.34 7.71E-01 0.55 2.10E+00 0.29
80 3.18E-01 0.48 6.96E-01 0.56 4.92E-01 0.64 1.53E+00 0.45
160 1.99E-01 0.67 4.12E-01 0.75 2.96E-01 0.73 9.73E-01 0.65
320 1.12E-01 0.82 2.24E-01 0.87 1.68E-01 0.81 5.55E-01 0.80
10
8
6
4 2
0
0 2 4 6 8 10
x
y
rho 10
8
6
4 2
0
0 2 4 6 8 10
x
y
rhou 10
8
6
4 2
0
0 2 4 6 8 10
x
y
rhov 10
8
6
4 2
0
0 2 4 6 8 10
x
y
rhoE 10
8
6
4 2
0
0 2 4 6 8 10
x
y
+2.465 norm vel 0.9824440 10 0.00020 0.00020 0.000280 0.000276
0.9824365 8 0.00015 0.00015 0.000255 0.000246
0.00010 0.00010 0.000230 0.000216 0.9824290
0.9824215 6 0.00005 0.00005 0.000205 0.000186
0.00000 0.00000 0.000180 0.000156 0.9824140 4 0.00005 0.00005 0.000155 0.000126 0.9824065 0.000130 0.000096 0.00010 0.00010 0.9823990 2 0.000105 0.000066 0.00015 0.00015
0.9823915 0.000080 0.000036 0.9823840 0 0.00020 0.00020 0.000055 0.000006
0 2 4 6 8 10
x
(a) FV-1
y
rho 10
8
6
4 2
0
0 2 4 6 8 10
x
y
rhou 10
8
6
4 2
0
0 2 4 6 8 10
x
y
rhov 10
8
6
4 2
0
0 2 4 6 8 10
x
y
rhoE 10
8
6
4 2
0
0 2 4 6 8 10
x
y
norm vel 10 0.5 0.984 0.300 0.4 2.48 0.45
0.952 8 0.225 0.3 2.40 0.40
0.920 0.150 0.2 2.32 0.35
0.888 6 0.075 0.1 2.24 0.30
0.856 0.000 0.0 2.16 0.25 0.824 4 0.075 0.1 2.08 0.20 0.792 0.150 0.2 2.00 0.15 0.760 2 0.225 0.3 1.92 0.10
0.728 0.300 0.4 1.84 0.05 0.696 0 1.76 0.00
0 2 4 6 8 10
x
(b) FV-2
y
rho 10
8
6
4 2
0
0 2 4 6 8 10
x
y
rhou 10
8
6
4 2
0
0 2 4 6 8 10
x
y
rhov 10
8
6
4 2
0
0 2 4 6 8 10
x
y
rhoE 10
8
6
4 2
0
0 2 4 6 8 10
x
y
norm vel 0.97 0.60 0.60 2.52 0.675 0.600
0.45 0.45 2.32 0.91 0.525
0.30 0.30 2.12 0.85 0.450
0.15 0.15 1.92 0.79 0.375
0.73 0.00 0.00 1.72 0.300 0.67 0.15 0.15 1.52 0.225 0.61 0.30 0.30 1.32 0.150 0.45 0.45
0.55 1.12 0.075 0.60 0.60 0.49 0.92 0.000
(c) GF
Figure 7: Euler equations: isentropic vortex with u = v = 0. Isocontours of the velocity norm obtained
0 0
with FV-1, FV-2 and GF after a long time integration (t =200).
f
31

rho rhou rho rhov rhou rho rhoE rhov rhou norm vel rhoE rhov norm vel rhoE norm vel
| 10  | 1.0035 10 10 | 10 10 10 | 10 10 10 | 10 10 | 10  |     |     | 10 10 |     | 4.0 | 10  |     |     |
| --- | ------------ | -------- | -------- | ----- | --- | --- | --- | ----- | --- | --- | --- | --- | --- |
1.052 0.98 1.052 0.985 3.546 1.38 1.245 1.460 4.02 1.32 1.80 1.72
|     | 0.9960 |     | 1.38 |     |     |     |     |     |     | 3.8 |     |     |     |
| --- | ------ | --- | ---- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
1.028 0.93 1.028 0.955 3.510 1.23 1.155 1.435 3.72 1.22 1.65 1.62
| 8   | 8 8    | 8 8 8 | 1.23 8 8 8 | 8 8        | 8   |       |           | 8 8 |     |          | 8   |     |      |
| --- | ------ | ----- | ---------- | ---------- | --- | ----- | --------- | --- | --- | -------- | --- | --- | ---- |
|     | 0.9885 | 0.88  | 1.004      | 3.474      |     | 1.410 | 3.42 1.12 |     |     | 1.50 3.6 |     |     | 1.52 |
|     |        | 1.004 | 1.08 0.925 | 1.08 1.065 |     |       |           |     |     |          |     |     |      |
6 0.9810 6 6 0.83 6 6 6 6 6 6 3.438 6 6 6 3.12 1.02 6 6 1.35 3.4 6 1.42
|     |     | 0.980 | 0.980 0.93 0.895 | 0.93 0.975 |     | 1.385 |     |     |     |     |     |     |     |
| --- | --- | ----- | ---------------- | ---------- | --- | ----- | --- | --- | --- | --- | --- | --- | --- |
y 0.9735 y y y 0.78 y y y y y y y y 2.82 y 0.92 y 1.20 3.2 y 1.32
|     |        | 0.956      | 0.956 0.78 0.865 | 3.402 0.78 0.885 |     | 1.360 |           |     |     |      |     |     |      |
| --- | ------ | ---------- | ---------------- | ---------------- | --- | ----- | --------- | --- | --- | ---- | --- | --- | ---- |
| 4   | 4 4    | 0.73 4 4 4 | 4 4 4            | 4 4              | 4   |       | 2.52 0.82 | 4 4 |     | 1.05 | 4   |     | 1.22 |
|     | 0.9660 | 0.932      | 0.932 0.63 0.835 | 3.366 0.63 0.795 |     | 1.335 |           |     |     | 3.0  |     |     |      |
|     |        | 0.68       |                  |                  |     |       | 2.22 0.72 |     |     | 0.90 |     |     | 1.12 |
|     | 0.9585 | 0.908      | 0.908 0.805      | 3.330 0.48 0.705 |     | 1.310 |           |     |     | 2.8  |     |     |      |
| 2   | 2 2    | 0.63 2 2 2 | 0.48 2 2 2       | 2 2              | 2   |       | 1.92 0.62 | 2 2 |     | 0.75 | 2   |     | 1.02 |
|     | 0.9510 |            |                  |                  |     |       |           |     |     | 2.6  |     |     |      |
0.884 0.58 0.884 0.33 0.775 3.294 0.33 0.615 1.285 1.62 0.52 0.60 0.92
| 0   | 0 0 | 0 0 0 | 0 0 0 | 0 0 | 0   |     |     | 0 0 |     |     | 0   |     |     |
| --- | --- | ----- | ----- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
0.9435 0.860 0.53 0.860 0.18 0.745 3.258 0.18 0.525 1.260 1.32 0.42 0.45 2.4 0.82
0 2 4 x 6 8 10 0 2 0 4 2 x 6 4 x 8 6 10 8 10 0 2 0 4 2 0 x 6 4 2 x 8 6 4 10 x 8 6 10 8 10 0 2 0 4 2 0 x 6 4 2 x 8 6 4 10 x 8 6 10 8 10 0 2 0 4 2 0 x 6 4 2 x 8 6 4 10 x 8 6 10 8 10 0 2 0 4 2 x 6 4 x 8 6 10 8 10 0 2 4 x 6 8 10
|     |     |     |     |     | (a) FV-1 |     |     | (b) | FV-2 |     |     | (c) GF |     |
| --- | --- | --- | --- | --- | -------- | --- | --- | --- | ---- | --- | --- | ------ | --- |
Figure 8: Euler equations: isentropic vortex with u = v = 1. Isocontours of the velocity norm obtained
|     |     |     |     |            |          |         |      |     | 0 0 |     |     |     |     |
| --- | --- | --- | --- | ---------- | -------- | ------- | ---- | --- | --- | --- | --- | --- | --- |
|     |     |     |     | with FV-1, | FV-2 and | GF at t | =10. |     |     |     |     |     |     |
f
is studied in depth e.g. in [36, 37]. A connection between stationarity preservation and
preservation of asymptotic has not been rigorously proved, however previous studies have
shown that stationarity preserving methods are also well behaved for low Mach [14, 9].
Following the last references, in this section we investigate the low Mach behaviour of the
new stationarity preserving formulation. To this end we compute long time simulations of the
same vortex for Mach numbers Ma = 10−2, Ma = 10−4, and Ma = 10−6. The contours of the
velocity norm at time t = 200 on a 40×40 mesh are visualized on figure 9. The figure shows
f
clearly that the first order finite volume method without reconstruction is unable to provide a
reasonable approximation on this mesh level, giving essentially a constant state. The second
order method improves this, however considerably losing the circular symmetry of the vortex,
withsignificantdissipation. Thenewmethodprovidesaremarkableapproximation,withvery
little dissipation, and perfect circular symmetry despite of the coarse Cartesian mesh.
To provide a quantitative assessment of the enhancements brought by our approach, on
figure 10 we report the convergence of the error of the x momentum at t = 200. The error
f
is scaled by the exact maximum value (which is of the order of the Mach number itself), to
allow a comparison of the results for the three values of the Mach number. The results for the
FV method are similar to the unfiltered ones presented in [36] for the velocity: the method
barely converges on the coarser resolutions, and starts converging with a slope of about 1/2
only on the finest resolution (which is considerably finer that those used in the reference).
Despite of the poor qualitative results obtained on the coarsest meshes, the second order
FV method does provide the expected slope for all Mach numbers. However, our approach
provideserrorswhicharesystematicallytwoordersofmagnitudesmaller, forallgridsizesand
meshes. This supports the idea that stationarity preservation and asymptotic preservation
|     |     |     |     | are tightly | linked.      |        |            |        |     |     |     |     |     |
| --- | --- | --- | --- | ----------- | ------------ | ------ | ---------- | ------ | --- | --- | --- | --- | --- |
|     |     |     |     | 6.2.3.      | Perturbation | of the | isentropic | vortex |     |     |     |     |     |
In this section, we present a test case for the Euler equations that consists in a perturb-
ation of the isentropic vortex presented in the previous section. The initial conditions for
the three schemes FV-1, FV-2 and GF are taken as the final results q of the respective
eq
|     |     |     |     | simulations | run | until final | time t | = 50 with | a 80×80 | mesh. |     |     |     |
| --- | --- | --- | --- | ----------- | --- | ----------- | ------ | --------- | ------- | ----- | --- | --- | --- |
f
Then, we add to the initial conditions a density perturbation δρ centered in (4,4) of the
32

10
8
6
4
2
0
0 2 4 6 8 10
x
y
rho
10
8
6
4
2
0
0 2 4 6 8 10
x
y
rhou
10
8
6
4
2
0
0 2 4 6 8 10
x
y
rhov
10
8
6
4
2
0
0 2 4 6 8 10
x
y
rhoE
10
8
6
4
2
0
0 2 4 6 8 10
x
y
1e 11+9.999922747e 1 1e 8 1e 81.47 1e 10+2.499983777 norm vel 1e 8
9.978 10 1.23 1.17 6.523 1.96
9.621 8 0.93 0.87 6.024 1.76
9.264 0.63 0.57 5.525 1.56
8.907 6 0.33 0.27 5.026 1.36
8.550 0.03 0.03 4.527 1.16
8.194 4 0.27 0.33 4.028 0.96
7.837 0.57 0.63 3.529 0.76
7.480 2 0.87 0.93 3.031 0.56
7.123 1.17 1.23 2.532 0.36
0
6.766 1.47 2.033 0.16
0 2 4 6 8 10
x
(a) FV-1, Ma = 10−2
y
rho
10
8
6
4
2
0
0 2 4 6 8 10
x
y
rhou
10
8
6
4
2
0
0 2 4 6 8 10
x
y
rhov
10
8
6
4
2
0
0 2 4 6 8 10
x
y
rhoE
10
8
6
4
2
0
0 2 4 6 8 10
x
y
1e 13+9.999999992e 1 1e 10 1e 101.474 1e 13+2.499999998 norm vel 1e 10
10 5.352 1.233 1.173 7.496 1.96
5.145 8 0.932 0.872 7.243 1.76
4.937 0.632 0.571 6.986 1.56
4.728 6 0.331 0.271 6.732 1.36
4.521 0.030 0.030 6.479 1.16
4.313 4 0.271 0.331 6.222 0.96
4.106 0.571 0.632 5.969 0.76
3.897 2 0.872 0.932 5.715 0.56
3.689 1.173 1.233 5.458 0.36
0
3.482 1.474 5.205 0.16
0 2 4 6 8 10
x
(b) FV-1, Ma = 10−4
y
rho
10
8
6
4
2
0
0 2 4 6 8 10
x
y
rhou
10
8
6
4
2
0
0 2 4 6 8 10
x
y
rhov
10
8
6
4
2
0
0 2 4 6 8 10
x
y
rhoE
10
8
6
4
2
0
0 2 4 6 8 10
x
y
1e 13+1 1e 121.573 1e 121.573 1e 13+2.5 norm vel 1e 12
0.067 1.252 1.252 9.224 1.96
0.138 0.931 0.931 9.019 1.76
0.343 0.610 0.610 8.815 1.56
0.547 0.289 0.289 8.615 1.36
0.753 0.032 0.032 8.411 1.16
0.957 0.353 0.353 8.207 0.96
1.162 0.674 0.674 8.002 0.76
1.367 0.995 0.995 7.798 0.56
1.572 1.316 1.316 7.594 0.36
1.776 7.390 0.16
(c) FV-1, Ma = 10−6
10
8
6
4
2
0
0 2 4 6 8 10
x
y
rho
10
8
6
4
2
0
0 2 4 6 8 10
x
y
rhou
10
8
6
4
2
0
0 2 4 6 8 10
x
y
rhov
10
8
6
4
2
0
0 2 4 6 8 10
x
y
rhoE
10
8
6
4
2
0
0 2 4 6 8 10
x
y
1e 5+9.999e 1 +2.5 norm vel
9.868 10 0.003768 0.002935 0.0000275 0.0040
0.002189
9.422 0.003019 0.0000158 0.0035
8 0.001444
8.977 0.002270 0.000699 0.0000042 0.0030
8.531 6 0.001521 0.000047 0.0000074 0.0025
8.086 0.000772 0.0000191
0.000792 0.0020
7.640 4 0.000023 0.0000307
7.195 0.000726 0.001538 0.0000423 0.0015
0.002283
6.750 2 0.001475 0.0000540 0.0010
0.003028
6.304 0.002223 0.0000656 0.0005
5.859 0 0.002972 0.003774 0.0000772 0.0000
0 2 4 6 8 10
x
(d) FV-2, Ma = 10−2
y
rho
10
8
6
4
2
0
0 2 4 6 8 10
x
y
rhou
10
8
6
4
2
0
0 2 4 6 8 10
x
y
rhov
10
8
6
4
2
0
0 2 4 6 8 10
x
y
rhoE
10
8
6
4
2
0
0 2 4 6 8 10
x
y
1e 7+1 1e 5 1e 5 1e 7+2.5 norm vel 1e 5
0.889 10 3.780 2.958 3.149 4.0
2.210
0.671 3.031 2.387 3.5
8 1.461
0.454 2.282 0.712 1.624 3.0
0.236 6 1.534 0.037 0.862 2.5
0.018 0.785 0.099
0.785 2.0
0.200 4 0.036 0.663
0.418 0.713 1.534 1.426 1.5
2.283
0.636 2 1.461 2.189 1.0
3.032
0.854 2.210 2.951 0.5
1.071 0 2.959 3.780 3.714 0.0
0 2 4 6 8 10
x
(e) FV-2, Ma = 10−4
y
rho
10
8
6
4
2
0
0 2 4 6 8 10
x
y
rhou
10
8
6
4
2
0
0 2 4 6 8 10
x
y
rhov
10
8
6
4
2
0
0 2 4 6 8 10
x
y
rhoE
10
8
6
4
2
0
0 2 4 6 8 10
x
y
1e 9+1 1e 7 1e 7 1e 9+2.5 norm vel 1e 7
2.959 0.893 3.780 3.127 4.0
2.210
0.675 3.031 2.364 3.5
1.461
0.458 2.283 0.712 1.602 3.0
0.240 1.534 0.036 0.839 2.5
0.022 0.785 0.076
0.785 2.0
0.196 0.036 0.686
0.414 0.712 1.534 1.449 1.5
2.283
0.632 1.461 2.211 1.0
3.031
0.850 2.210 2.974 0.5
3.780 1.068 2.959 3.736 0.0
(f) FV-2, Ma = 10−6
10
8
6
4
2
0
0 2 4 6 8 10
x
y
rho
10
8
6
4
2
0
0 2 4 6 8 10
x
y
rhou
10
8
6
4
2
0
0 2 4 6 8 10
x
y
rhov
10
8
6
4
2
0
0 2 4 6 8 10
x
y
rhoE
10
8
6
4
2
0
0 2 4 6 8 10
x
y
+2.499 norm vel 0.016
10
0.9999807 0.01310 0.01310 0.0009441 0.014
0.9999564 0.00990 0.00990 0.0008588 8 0.012
0.9999320 0.00671 0.00671 0.0007735
0.99990776 0.00351 0.00351 0.0006881 0.010
0.9998833 0.00032 0.00032 0.0006028 0.008
0.99985904 0.00287 0.00287 0.0005175 0.006
0.9998346 0.00607 0.00607 0.0004321
0.99981032 0.00926 0.00926 0.0003468 0.004
0.9997859 0.01246 0.01246 0.0002615 0.002
0
0.9997616 0.01565 0.01565 0.0001761 0.000
0 2 4 6 8 10
x
(g) GF, Ma = 10−2
y
rho
10
8
6
4
2
0
0 2 4 6 8 10
x
y
rhou
10
8
6
4
2
0
0 2 4 6 8 10
x
y
rhov
10
8
6
4
2
0
0 2 4 6 8 10
x
y
rhoE
10
8
6
4
2
0
0 2 4 6 8 10
x
y
1e 8+9.999999e 1 1e 8+2.4999999 norm vel 0.00016
10
9.806 0.0001301 0.0001301 9.431 0.00014
9.563 0.0000983 0.0000983 8.579 8 0.00012
9.319 0.0000666 0.0000666 7.727
9.076 6 0.0000349 0.0000349 6.875 0.00010
8.833 0.0000032 0.0000032 6.023 0.00008
8.590 4 0.0000285 0.0000285 5.171 0.00006
8.347 0.0000603 0.0000603 4.319
8.103 2 0.0000920 0.0000920 3.467 0.00004
7.860 0.0001237 0.0001237 2.616 0.00002
0
7.617 0.0001554 0.0001554 1.764 0.00000
0 2 4 6 8 10
x
(h) GF, Ma = 10−4
y
rho
10
8
6
4
2
0
0 2 4 6 8 10
x
y
rhou
10
8
6
4
2
0
0 2 4 6 8 10
x
y
rhov
10
8
6
4
2
0
0 2 4 6 8 10
x
y
rhoE
10
8
6
4
2
0
0 2 4 6 8 10
x
y
1e 12+10e 1 1e 6 1e 61.554 1e 12+2.5 norm vel 1e 61.6
9.895 1.300 1.237 9.515 1.4
9.630 0.983 0.920 8.642 1.2
9.366 0.666 0.603 7.770
9.101 0.349 0.285 6.898 1.0
8.837 0.032 0.032 6.025 0.8
8.573 0.285 0.349 5.153 0.6
8.308 0.603 0.666 4.281
0.4 8.044 0.920 0.983 3.408
7.780 1.237 1.300 2.536 0.2
7.515 1.554 1.664 0.0
(i) GF, Ma = 10−6
Figure 9: Euler equations: Low Mach stationary isentropic vortex. Isocontours of the velocity norm at
t =200 obtained with FV-1 (top row), FV-2 (middle row), and GF (bottom row). Left: Ma=0.01. Center:
f
Ma = 0.0001. Right: Ma=0.000001.
33

| (a) FV-1 |     |     | (b) | FV-2 |     |     | (c) | GF  |
| -------- | --- | --- | --- | ---- | --- | --- | --- | --- |
Figure 10: Euler equations: Low Mach stationary isentropic vortex. Mesh convergence of the error on the
x momentum (scaled by the exact maximum value) at t = 200 and different Mach numbers for the FV-1
f
(left), FV-2 (center), and GF (right). On y-axis the log (∥ρu−(ρu) ∥ /∥(ρu) ∥ ) while on the
|     |     |     |     | 10  |     | exact L∞ | exact | L∞  |
| --- | --- | --- | --- | --- | --- | -------- | ----- | --- |
x-axis the log (∆x).
10
Table 4: Euler equations: isentropic vortex with u =v =1 (t =10). L error and order of accuracy n˜ for
|     |     |     |     | 0 0 | f   | 2   |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
FV-1, FV-2 and GF methods.
|      |     | ρ   | ρu  |     |     | ρv  | ρE  |     |
| ---- | --- | --- | --- | --- | --- | --- | --- | --- |
| N ,N | L   | n˜  | L   | n˜  | L   | n˜  | L   | n˜  |
| x y  | 2   |     | 2   |     | 2   |     | 2   |     |
FV-1
| 20  | 6.50E-01 | –    | 1.54E+00 | –    | 1.54E+00 | –    | 3.12E+00 | –    |
| --- | -------- | ---- | -------- | ---- | -------- | ---- | -------- | ---- |
| 40  | 6.21E-01 | 0.06 | 1.46E+00 | 0.07 | 1.46E+00 | 0.07 | 3.01E+00 | 0.05 |
| 80  | 5.82E-01 | 0.09 | 1.31E+00 | 0.15 | 1.31E+00 | 0.15 | 2.83E+00 | 0.09 |
| 160 | 5.13E-01 | 0.18 | 1.06E+00 | 0.30 | 1.07E+00 | 0.29 | 2.48E+00 | 0.19 |
| 320 | 4.01E-01 | 0.35 | 7.58E-01 | 0.49 | 7.63E-01 | 0.49 | 1.92E+00 | 0.36 |
FV-2
| 20  | 5.29E-01 | –    | 1.07E+00 | –    | 1.12E+00 | –    | 2.48E+00 | –    |
| --- | -------- | ---- | -------- | ---- | -------- | ---- | -------- | ---- |
| 40  | 2.45E-01 | 1.10 | 4.42E-01 | 1.28 | 4.84E-01 | 1.21 | 1.17E+00 | 1.08 |
| 80  | 6.55E-02 | 1.90 | 1.26E-01 | 1.80 | 1.37E-01 | 1.81 | 2.82E-01 | 2.04 |
| 160 | 1.85E-02 | 1.82 | 3.22E-02 | 1.97 | 3.42E-02 | 2.00 | 6.03E-02 | 2.22 |
| 320 | 3.38E-03 | 2.45 | 7.49E-03 | 2.10 | 8.12E-03 | 2.07 | 1.28E-02 | 2.23 |
GF
| 20  | 5.46E-01 | –    | 1.30E+00 | –    | 1.13E+00 | –    | 2.58E+00 | –    |
| --- | -------- | ---- | -------- | ---- | -------- | ---- | -------- | ---- |
| 40  | 4.44E-01 | 0.30 | 1.02E+00 | 0.34 | 7.71E-01 | 0.55 | 2.10E+00 | 0.29 |
| 80  | 3.18E-01 | 0.48 | 6.96E-01 | 0.56 | 4.92E-01 | 0.64 | 1.53E+00 | 0.45 |
| 160 | 1.99E-01 | 0.67 | 4.12E-01 | 0.75 | 2.96E-01 | 0.73 | 9.73E-01 | 0.65 |
| 320 | 1.12E-01 | 0.82 | 2.24E-01 | 0.87 | 1.68E-01 | 0.81 | 5.55E-01 | 0.80 |
34

rho rhou rho rhov rhou rho rhoE rhov rhou norm vel rhoE rhov norm vel rhoE norm vel
10 10 10 10 10 10 0.0100 10 0.0045 10 10 0.010 10 10 10 0.00276 0.00125 10 10 0.0135 10
|     |     | 0.00240 |     |     | 0.0024 0.005 | 0.0024 |               | 0.0100 | 0.0028 |         | 0.0125  |               | 0.00288 |
| --- | --- | ------- | --- | --- | ------------ | ------ | ------------- | ------ | ------ | ------- | ------- | ------------- | ------- |
|     |     |         |     |     |              |        | 0.0075 0.0040 | 0.008  |        | 0.00246 | 0.00100 | 0.0120 0.0008 |         |
|     |     |         |     |     | 0.0018 0.004 | 0.0018 |               | 0.0075 | 0.0024 |         | 0.0100  |               | 0.00252 |
8 0.00192 8 8 8 8 8 0.0050 8 0.0035 8 8 8 8 8 0.00216 0.00075 8 8 0.0105 0.0004 8
|     |     |         |     |     | 0.0012 0.003 | 0.0012 |               | 0.006 0.0050 | 0.0020 |               | 0.0075  |               | 0.00216 |
| --- | --- | ------- | --- | --- | ------------ | ------ | ------------- | ------------ | ------ | ------------- | ------- | ------------- | ------- |
|     |     | 0.00144 |     |     |              |        | 0.0025 0.0030 |              |        | 0 . 0 0 1 8 6 | 0.00050 | 0.0090 0.0000 |         |
6 6 6 0.0006 6 6 6 0.0006 6 6 6 0.004 0.0025 6 0.0016 6 6 0.0050 6 6 6 0.00180
|     |     | 0.00096 |     |     | 0.002 |     | 0.0000 0.0025 |     |     |     | 0.00025 | 0.0075 |     |
| --- | --- | ------- | --- | --- | ----- | --- | ------------- | --- | --- | --- | ------- | ------ | --- |
y y y 0.0000 y y y 0.0000 y y y 0.002 0.0000 y 0.0012 y y 0 . 0 0 1 5 6 0.0025 y y 0.0004 y
|     |     | 0.00048 |     |     | 0.001 |     | 0.0025 0.0020 |     |     |     | 0.00000 | 0.0060 | 0.00144 |
| --- | --- | ------- | --- | --- | ----- | --- | ------------- | --- | --- | --- | ------- | ------ | ------- |
4 4 4 0.0006 4 4 4 0.0006 4 4 4 0.000 0.0025 4 0.0008 4 4 0.00126 0.0000 4 4 0.0008 4
|     |     |     |     |     | 0.000 |     | 0.0050 0.0015 |     |     |     | 0.00025 |     | 0.00108 |
| --- | --- | --- | --- | --- | ----- | --- | ------------- | --- | --- | --- | ------- | --- | ------- |
0.00000 0.0012 0.0012 0.0050 0 . 0 0 0 4 0.00096 0.0025 0.0045 0.0012
|     |     |     |     |     |     |     | 0.0075 0.0010 | 0.002 | 0.0075 |     | 0.00050 |     | 0.00072 |
| --- | --- | --- | --- | --- | --- | --- | ------------- | ----- | ------ | --- | ------- | --- | ------- |
2 0.00048 2 2 0.0018 0.001 2 2 2 0.0018 2 2 2 2 0 . 0 0 0 0 2 2 0.00066 0.0050 2 2 0.0030 2
|     |     |     |     |     |     |     | 0.0100 0.0005 | 0.004 | 0.0100 |     | 0.00075 | 0.0016 |     |
| --- | --- | --- | --- | --- | --- | --- | ------------- | ----- | ------ | --- | ------- | ------ | --- |
0.00096 0.0024 0.002 0.0024 0.0004 0.00036 0.0075 0.0015 0.00036
0 0 0 0 0 0 0.0125 0 0.0000 0 0 0.006 0.0125 0 0 0 0.00100 0 0 0.0020 0
|     |     |     |     |     |     |     |     |     | 0.0008 | 0.00006 | 0.0100 | 0.0000 | 0.00000 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | ------ | ------- | ------ | ------ | ------- |
0 2 4 x 6 8 10 0 2 0 4 2 x 6 4 x 8 6 10 8 10 0 2 0 4 2 0 x 6 4 2 x 8 6 4 10 x 8 6 10 8 10 0 2 0 4 2 0 x 6 4 2 x 8 6 4 10 x 8 6 10 8 10 0 2 0 4 2 0 x 6 4 2 x 8 6 4 10 x 8 6 10 8 10 0 2 0 4 2 x 6 4 x 8 6 10 8 10 0 2 4 x 6 8 10
|     | (a) FV-1 |     |     | (b) FV-2 |     | (c) GF |     |     |     |     |     |     |     |
| --- | -------- | --- | --- | -------- | --- | ------ | --- | --- | --- | --- | --- | --- | --- |
Figure 11: Euler equations: perturbation of the isentropic vortex. Isocontours of the ρ−ρ norm obtained
eq
| with FV-1, | FV-2 and | GF at final | time t =2 | with a 80×80 | mesh. |     |     |     |     |     |     |     |     |
| ---------- | -------- | ----------- | --------- | ------------ | ----- | --- | --- | --- | --- | --- | --- | --- | --- |
f
form:
|     |     |     |     | Ae−(x−4)2 | + (y−4)2 |     |     |     |     |     |     |     |     |
| --- | --- | --- | --- | --------- | -------- | --- | --- | --- | --- | --- | --- | --- | --- |
δρ =
σ 2
5·10−3
where A = and σ = 0.8. The simulation is run until a final time t = 2 to compare
the effect of the numerical viscosity on the evolution of the perturbation.
In figure 11, we show the density contour plot at the final time for the three methods. The
GF method is able to capture the perturbation sharply, while the FV-1 and FV-2 methods
have discretization errors too large to capture it properly. By looking at the isocontours
scales, it is clear that the perturbation is completely dissipated for the FV-1 method, while
for the FV-2 method the perturbation is still visible but with a much larger error compared
| to the | expected solution. |         |     |     |     |     |     |     |     |     |     |     |     |
| ------ | ------------------ | ------- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 6.2.4. | Sod’s circular     | problem |     |     |     |     |     |     |     |     |     |     |     |
Here, we test the robustness of the global flux method on the Euler equations for the Sod
circular problem. This case is fully non-linear, and with Mach number of order one, allowing
to show that the new method has “all Mach” capabilities. The problem is a two-dimensional
extension of the classical shock tube problem. The simulation is performed on a domain
| [−1,1]×[−1,1] | and | the initial | condition | is given | by  |     |     |     |     |     |     |     |     |
| ------------- | --- | ----------- | --------- | -------- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
(cid:40)
|     |     |     |     | Q   | if r ≤ R, |     |     |     |     |     |     |     |     |
| --- | --- | --- | --- | --- | --------- | --- | --- | --- | --- | --- | --- | --- | --- |
i
Q(x,0) =
|     |     |     |     | Q   | if r > R, |     |     |     |     |     |     |     |     |
| --- | --- | --- | --- | --- | --------- | --- | --- | --- | --- | --- | --- | --- | --- |
e
(cid:112)
with r = x2 +y2. The circle of radius R = 0.5 is centered in the origin and separates
the inner state Q from the outer state Q , where Q = (ρ,u,v,p). The initial conditions are
|     | i   |     |     | e   |     |     |     |     |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
given by Q = (1,0,0,1) and Q = (0.125,0,0,0.1). For a reference solution of this problem,
|          | i        |     | e   |     |     |     |     |     |     |     |     |     |     |
| -------- | -------- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| we refer | to [18]. |     |     |     |     |     |     |     |     |     |     |     |     |
To have a smoother initial condition, the two states are connected by a smooth transition
| region | given by an | erfc function | defined | as         |              |     |     |     |     |     |     |     |     |
| ------ | ----------- | ------------- | ------- | ---------- | ------------ | --- | --- | --- | --- | --- | --- | --- | --- |
|        |             |               |         | 1 (cid:18) | r−R (cid:19) |     |     |     |     |     |     |     |     |
|        |             |               | ζ(r)    | = erfc     | ,            |     |     |     |     |     |     |     |     |
|        |             |               |         | 2          | δ            |     |     |     |     |     |     |     |     |
35

|     | 1.00 |     |     |     | 1.00 |     |     | 1.00 |     |     |     |
| --- | ---- | --- | --- | --- | ---- | --- | --- | ---- | --- | --- | --- |
|     | 0.75 |     |     |     | 0.75 |     |     | 0.75 |     |     |     |
|     | 0.50 |     |     |     | 0.50 |     |     | 0.50 |     |     |     |
|     | 0.25 |     |     |     | 0.25 |     |     | 0.25 |     |     |     |
|     | 0.00 |     |     |     | 0.00 |     |     | 0.00 |     |     |     |
|     | 0.25 |     |     |     | 0.25 |     |     | 0.25 |     |     |     |
|     | 0.50 |     |     |     | 0.50 |     |     | 0.50 |     |     |     |
|     | 0.75 |     |     |     | 0.75 |     |     | 0.75 |     |     |     |
|     | 1.00 |     |     |     | 1.00 |     |     | 1.00 |     |     |     |
1.00 0.75 0.50 0.25 0.00 0.25 0.50 0.75 1.00 1.00 0.75 0.50 0.25 0.00 0.25 0.50 0.75 1.00 1.00 0.75 0.50 0.25 0.00 0.25 0.50 0.75 1.00
|     | (a) | FV-1 | density |     |     | (b) FV-2 density |     |     | (c) GF | density |      |
| --- | --- | ---- | ------- | --- | --- | ---------------- | --- | --- | ------ | ------- | ---- |
|     | 1   |      | GF      |     | 1   | GF               |     |     | 6      |         | GF   |
|     |     |      | FV-1    |     |     | FV-1             |     |     |        |         | FV-1 |
|     |     |      | FV-2    |     |     | FV-2             |     |     |        |         | FV-2 |
4
| ρ   |     |     |     | v   | 0.5 |     |     | p   |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
0.5
2
0
|     | 0   | 0.5 | 1   |     | 0   | 0.5 | 1   |     | 0   | 0.5 | 1   |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
|     |     | y   |     |     |     | y   |     |     |     | y   |     |
(d) Slice of ρ at x=0 (e) Slice of v at x=0 (f) Slice of p at x=0
Figure 12: Euler equations: Sod’s circular problem. Numerical results obtained on a 400×400 mesh with
| FV-1, FV-2 | and GF | methods | run | until a | final time | t =0.2. |     |     |     |     |     |
| ---------- | ------ | ------- | --- | ------- | ---------- | ------- | --- | --- | --- | --- | --- |
f
36

where δ = 0.01. Therefore, we can define the smoothed initial condition as
|     |     |     | Q(r,0) | = ζ(r)Q | +(1−ζ(r))Q |     |     | .   |     |
| --- | --- | --- | ------ | ------- | ---------- | --- | --- | --- | --- |
|     |     |     |        |         | i          |     | e   |     |     |
The simulation is run until final time t = 0.2 before the shock waves reach the boundaries.
f
In figure 12, we present the numerical results obtained on a 400×400 mesh with the three
methods FV-1, FV-2 and GF. For all situations, we show the density contour plot at the final
time, alongwithasliceofthedensityandverticalvelocityonthex = 0axis. Itcanbenoticed
that, among all three simulations, the GF performs much better than the standard first order
scheme and it is clearly comparable to a second order one. It is able to sharply capture the
three waves, which are smoothed out by the classical FV-1. GF also avoids oscillations at
the beginning of the rarefaction, while the FV-2 shows small oscillations. On the foot of the
rarefaction, the GF show sharper results, while it is a little more diffusive on the contact
discontinuity with respect to the FV-2. On the shock, the GF does not oscillate, while the
FV-2 shows minimal oscillations and is a little more sharply representing the discontinuity.
This test case not only allows us to show the robustness of the method to deal with
unsteady shock propagation. It also provides interesting insights into its low dissipation even
though the method has not been designed to have any particular properties on unsteady
solutions.
| 6.2.5. Kelvin-Helmholtz |     | instability |     |     |     |     |     |     |     |
| ----------------------- | --- | ----------- | --- | --- | --- | --- | --- | --- | --- |
We consider a smooth Kelvin-Helmholtz instability for the Euler equations, introduced
in [39] to further confirm the ability of our approach to correctly cope with the low Mach
limit. numerical scheme to cope with low Mach number flow and to assess qualitatively the
numerical diffusion of the method. There is a large body of work available in the literature
concerning the shortcomings of classical Finite Volume methods in the subsonic regime (see
e.g. [11]). The effect of stabilizing diffusion becomes bigger as the Mach number decreases,
making it necessary to use highly resolved grids in order to capture the features of the flow.
Numerical methods that are not low Mach number compliant typically also stabilize Kelvin-
| Helmholtz | setups | in an artificial | way. |     |     |     |     |     |     |
| --------- | ------ | ---------------- | ---- | --- | --- | --- | --- | --- | --- |
The simulations are performed in the domain [0,2]×[−1/2,1/2] until a final time t = 80.
f
| The initial | condition | is given | by the | following | primitive |     | variables:   |     |        |
| ----------- | --------- | -------- | ------ | --------- | --------- | --- | ------------ | --- | ------ |
|             | ρ = γ     | +H(y)r,  | u      | = M H(y), |           | v = | δM sin(2πx), |     | p = 1, |
where the Mach number parameter is M = 10−2 and we use r = 10−3 and δ = 0.1. The
| function | H(y) is defined | as,  |         |                        |     |      |         |     |         |
| -------- | --------------- | ---- | ------- | ---------------------- | --- | ---- | ------- | --- | ------- |
|          |                 |     | (cid:0) | (cid:0) (cid:1)(cid:1) |     |      |         |     |         |
|          |                 | −sin | π       | y + 1                  | ,   | if − | 1 − ω ≤ | y < | −1 + ω, |
|          |                 |     | ω       | 4                      |     |      | 4 2     |     | 4 2     |

|     |      |  −1, |                 |                |     | if −   | 1 + ω ≤ | y <   | 1 − ω, |
| --- | ---- | ------ | --------------- | -------------- | --- | ------ | ------- | ----- | ------ |
|     | H(y) | =      |                 |                |     |        | 4 2     |       | 4 2    |
|     |      |        | (cid:0) (cid:0) | (cid:1)(cid:1) |     |        |         |       |        |
|     |      | sin    | π y             | − 1 ,          |     | if 1 − | ω ≤ y   | < 1 + | ω,     |
|     |      |      | ω               | 4              |     | 4      | 2       | 4     | 2      |

|     |     |  1 |     |     |     | else, |     |     |     |
| --- | --- | --- | --- | --- | --- | ----- | --- | --- | --- |
where ω = 1/16. Observe that the shear flow is smooth such that for short times, there exists
a solution to which numerical methods converge upon mesh refinement ([39]).
37

|     |     | FV first order |     |     | FV second order |     |     | Global Flux |
| --- | --- | -------------- | --- | --- | --------------- | --- | --- | ----------- |
| 0.4 |     |                |     | 0.4 |                 |     | 0.4 |             |
| 0.2 |     |                |     | 0.2 |                 |     | 0.2 |             |
| 0.0 |     |                |     | 0.0 |                 |     | 0.0 |             |
| 0.2 |     |                |     | 0.2 |                 |     | 0.2 |             |
| 0.4 |     |                |     | 0.4 |                 |     | 0.4 |             |
0.00 0.25 0.50 0.75 1.00 1.25 1.50 1.75 2.00 0.00 0.25 0.50 0.75 1.00 1.25 1.50 1.75 2.00 0.00 0.25 0.50 0.75 1.00 1.25 1.50 1.75 2.00
|     | (a) | FV-1 (64×32)   |     | (b) | FV-2            | (64×32) | (c) | GF (64×32)  |
| --- | --- | -------------- | --- | --- | --------------- | ------- | --- | ----------- |
|     |     | FV first order |     |     | FV second order |         |     | Global Flux |
| 0.4 |     |                |     | 0.4 |                 |         | 0.4 |             |
| 0.2 |     |                |     | 0.2 |                 |         | 0.2 |             |
| 0.0 |     |                |     | 0.0 |                 |         | 0.0 |             |
| 0.2 |     |                |     | 0.2 |                 |         | 0.2 |             |
| 0.4 |     |                |     | 0.4 |                 |         | 0.4 |             |
0.00 0.25 0.50 0.75 1.00 1.25 1.50 1.75 2.00 0.00 0.25 0.50 0.75 1.00 1.25 1.50 1.75 2.00 0.00 0.25 0.50 0.75 1.00 1.25 1.50 1.75 2.00
|     | (d) | FV-1 (128×64)  |     | (e) | FV-2            | (128×64) | (f) | GF (128×64) |
| --- | --- | -------------- | --- | --- | --------------- | -------- | --- | ----------- |
|     |     | FV first order |     |     | FV second order |          |     | Global Flux |
| 0.4 |     |                |     | 0.4 |                 |          | 0.4 |             |
| 0.2 |     |                |     | 0.2 |                 |          | 0.2 |             |
| 0.0 |     |                |     | 0.0 |                 |          | 0.0 |             |
| 0.2 |     |                |     | 0.2 |                 |          | 0.2 |             |
| 0.4 |     |                |     | 0.4 |                 |          | 0.4 |             |
0.00 0.25 0.50 0.75 1.00 1.25 1.50 1.75 2.00 0.00 0.25 0.50 0.75 1.00 1.25 1.50 1.75 2.00 0.00 0.25 0.50 0.75 1.00 1.25 1.50 1.75 2.00
|     | (g) | FV-1 (256×128) |     | (h) | FV-2 (256×128)  |     | (i) | GF (256×128) |
| --- | --- | -------------- | --- | --- | --------------- | --- | --- | ------------ |
|     |     | FV first order |     |     | FV second order |     |     | Global Flux  |
| 0.4 |     |                |     | 0.4 |                 |     | 0.4 |              |
| 0.2 |     |                |     | 0.2 |                 |     | 0.2 |              |
| 0.0 |     |                |     | 0.0 |                 |     | 0.0 |              |
| 0.2 |     |                |     | 0.2 |                 |     | 0.2 |              |
| 0.4 |     |                |     | 0.4 |                 |     | 0.4 |              |
0.00 0.25 0.50 0.75 1.00 1.25 1.50 1.75 2.00 0.00 0.25 0.50 0.75 1.00 1.25 1.50 1.75 2.00 0.00 0.25 0.50 0.75 1.00 1.25 1.50 1.75 2.00
|     | (j) | FV-1 (512×256) |     | (k) | FV-2 (512×256) |     | (l) | GF (512×256) |
| --- | --- | -------------- | --- | --- | -------------- | --- | --- | ------------ |
Figure 13: Euler equations: Kelvin-Helmholtz instability. Density isocontours are presented for a set of
| nested | meshes | to compare | FV-1 | (top), FV-2 | (middle) | and GF | (bottom). |     |
| ------ | ------ | ---------- | ---- | ----------- | -------- | ------ | --------- | --- |
38

In figure 13, we present the numerical results obtained with FV-1, FV-2 and GF for
the Kelvin-Helmholtz instability arising from the aforementioned initial conditions. The
simulations are performed on a set of 4 nested grids from a 64×32, the coarsest, to 512×256,
the finest.
The FV-1 scheme is not able to capture any of the features arising from the instability.
| No vortices | form, since | FV-1 is not | low Mach | number | compliant. |
| ----------- | ----------- | ----------- | -------- | ------ | ---------- |
Much improved results are obtained using the FV-2 with a linear reconstruction of the
conservative variables. Here, the higher order of accuracy helps to overcome excessive diffu-
sion at this Mach number and for this simulation time. However, the structures still appear
diffused and would need even more resolution for the vortex details to be captured.
Very differently from these methods, the GF method is able to capture all details of the
flow very accurately. Already on the coarsest mesh, the fluid structures start to appear and
develop. Here, some spurious vortices are visible, which are a known artefact of virtually any
numerical method (see e.g. [20]). When increasing the resolution, the fluid features converge
to the expected solution found in other references [39]. Comparison to the results obtained
with low Mach compliant methods studied in [39] shows that the GF method is at least as
good.
It has been suggested in [10] that numerical methods for the Euler equations whose
linearization (= method for linear acoustics) is stationarity preserving, are low Mach number
compliant. A nonlinear stationarity preserving method naturally has this property, and some
experimental examples of this behavior can also be found in [9]. Thus, even though we set
out to improve the performance of the numerical method at stationary state, here we observe
that this property is beneficial even for solutions far away from it.
| 6.3. Shallow     | water system |     |     |     |     |
| ---------------- | ------------ | --- | --- | --- | --- |
| 6.3.1. Potential | flow         |     |     |     |     |
The first test case implemented for the shallow water equations isan equilibrium(see [49])
characterized by a known exact solution, for which it is possible to perform a convergence
analysis. The initial condition is a potential flow defined on the square [−1,1]×[−1,1] with
| Dirichlet | boundary conditions, | which  | is given | by     |      |
| --------- | -------------------- | ------ | -------- | ------ | ---- |
|           |                      | h(x,y) | = (x−x   | )(y −y | )+C, |
|           |                      |        |          | 0      | 0    |
|           |                      | u(x,y) | = (x−x   | ),     |      |
0
|     |     | v(x,y) | = −(y | −y ) |     |
| --- | --- | ------ | ----- | ---- | --- |
0
where C = 3/2 and (x ,y ) = (0,0). The 2D equilibrium is achieved thanks to a special
0 0
| bathymetry | given by |        |          |             |         |
| ---------- | -------- | ------ | -------- | ----------- | ------- |
|            |          |        | (cid:18) | +y2(cid:19) |         |
|            |          |        | 1        | x2          |         |
|            |          | b(x,y) | = 30−    |             | −xy −C. |
|            |          |        | g        | 2           |         |
The solution of this potential flow is shown in figure 14, and the convergence rates com-
puted at final time t = 1 are presented in table 5, demonstrating the improvement brought
f
about by the global flux formulation. Again, since the setup is stationary, superconvergence
is observed. Moreover, the new method is even able to outperform FV-2.
39

|      | h   |      |      | hu  |     |      | hv  |     |      | norm vel |     |
| ---- | --- | ---- | ---- | --- | --- | ---- | --- | --- | ---- | -------- | --- |
| 1.00 |     |      | 1.00 |     |     | 1.00 |     |     | 1.00 |          | 3.6 |
|      |     | 2.40 |      |     | 2.4 |      |     | 2.4 |      |          |     |
3.2
| 0.75 |     | 2.16 | 0.75 |     | 1.8 | 0.75 |     | 1.8 | 0.75 |     |     |
| ---- | --- | ---- | ---- | --- | --- | ---- | --- | --- | ---- | --- | --- |
2.8
| 0.50   |     |      | 0.50   |     | 1.2 | 0.50     |     | 1.2 | 0.50   |     |     |
| ------ | --- | ---- | ------ | --- | --- | -------- | --- | --- | ------ | --- | --- |
|        |     | 1.92 |        |     |     |          |     |     |        |     | 2.4 |
| 0.25   |     |      | 0.25   |     | 0.6 | 0.25     |     | 0.6 | 0.25   |     |     |
|        |     | 1.68 |        |     |     |          |     |     |        |     | 2.0 |
| y 0.00 |     |      | y 0.00 |     | 0.0 | y 0.00   |     | 0.0 | y 0.00 |     |     |
|        |     | 1.44 |        |     |     |          |     |     |        |     | 1.6 |
|        |     |      |        |     |     | 0.6      |     | 0.6 |        |     |     |
| 0.25   |     | 1.20 | 0.25   |     |     | 0.25     |     |     | 0.25   |     | 1.2 |
|        |     |      |        |     |     | 1.2      |     | 1.2 |        |     |     |
| 0.50   |     | 0.96 | 0.50   |     |     | 0.50     |     |     | 0.50   |     | 0.8 |
| 0.75   |     | 0.72 | 0.75   |     |     | 1.8 0.75 |     | 1.8 | 0.75   |     | 0.4 |
| 1.00   |     | 0.48 | 1.00   |     |     | 2.4 1.00 |     | 2.4 | 1.00   |     | 0.0 |
1.0 0.5 0.0 0.5 1.0 1.0 0.5 0.0 0.5 1.0 1.0 0.5 0.0 0.5 1.0 1.0 0.5 0.0 0.5 1.0
|     | x   |     |     | x   |     |     | x   |     |     | x   |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
Figure 14: Shallow water system: potential flow. Reference solution of the conservative variables.
Table 5: Shallow water system: potential flow (t = 1). L error and order of accuracy n˜ for FV-1, FV-2
|     |     |     |     | f   | 2   |     |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
and GF.
|     |      |     | h   | hu  |     | hv  |     |     |     |     |     |
| --- | ---- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
|     | N ,N | L   | n˜  | L   | n˜  | L   | n˜  |     |     |     |     |
|     | x y  |     | 2   | 2   |     | 2   |     |     |     |     |     |
FV-1
|     | 20  | 1.54E-02 | –    | 1.57E-01 | –    | 1.71E-01 | –    |     |     |     |     |
| --- | --- | -------- | ---- | -------- | ---- | -------- | ---- | --- | --- | --- | --- |
|     | 40  | 8.25E-03 | 0.89 | 1.08E-01 | 0.53 | 1.11E-01 | 0.62 |     |     |     |     |
|     | 80  | 4.30E-03 | 0.94 | 6.60E-02 | 0.71 | 6.43E-02 | 0.78 |     |     |     |     |
|     | 160 | 2.18E-03 | 0.97 | 3.68E-02 | 0.84 | 3.47E-02 | 0.88 |     |     |     |     |
|     | 320 | 1.10E-03 | 0.99 | 1.95E-02 | 0.91 | 1.80E-02 | 0.94 |     |     |     |     |
FV-2
|     | 20  | 2.49E-04 | –    | 1.06E-03 | –    | 1.47E-03 | –    |     |     |     |     |
| --- | --- | -------- | ---- | -------- | ---- | -------- | ---- | --- | --- | --- | --- |
|     | 40  | 5.26E-05 | 2.24 | 2.61E-04 | 2.02 | 3.25E-04 | 2.17 |     |     |     |     |
|     | 80  | 1.09E-05 | 2.27 | 7.11E-05 | 1.87 | 8.17E-05 | 1.99 |     |     |     |     |
|     | 160 | 2.24E-06 | 2.28 | 1.86E-05 | 1.93 | 2.06E-05 | 1.98 |     |     |     |     |
|     | 320 | 4.81E-07 | 2.21 | 4.69E-06 | 1.98 | 5.17E-06 | 1.99 |     |     |     |     |
GF
|     | 20  | 1.15E-04 | –    | 4.29E-04 | –    | 1.11E-03 | –    |     |     |     |     |
| --- | --- | -------- | ---- | -------- | ---- | -------- | ---- | --- | --- | --- | --- |
|     | 40  | 2.69E-05 | 2.09 | 1.01E-04 | 2.08 | 2.39E-04 | 2.21 |     |     |     |     |
|     | 80  | 6.49E-06 | 2.05 | 2.46E-05 | 2.03 | 5.50E-05 | 2.11 |     |     |     |     |
|     | 160 | 1.59E-06 | 2.02 | 6.08E-06 | 2.01 | 1.32E-05 | 2.05 |     |     |     |     |
|     | 320 | 3.95E-07 | 2.01 | 1.51E-06 | 2.00 | 3.24E-06 | 2.02 |     |     |     |     |
40

Table 6: Shallow water system: lake at rest (t =0.1). L error and order of accuracy n˜ for FV-1 and FV-2
|         |                |          |     | f   | 2   |     |     |     |
| ------- | -------------- | -------- | --- | --- | --- | --- | --- | --- |
| schemes | with the novel | GF.      |     |     |     |     |     |     |
|         |                |          |     | h   | hu  |     | hv  |     |
|         |                | N x ,N y | L 2 | n˜  | L 2 | n˜  | L 2 | n˜  |
FV-1
|     |     | 20  | 7.13E-03 | –    | 4.67E-02 | –    | 3.61E-02 | –    |
| --- | --- | --- | -------- | ---- | -------- | ---- | -------- | ---- |
|     |     | 40  | 2.79E-03 | 1.35 | 2.42E-02 | 0.95 | 2.04E-02 | 0.82 |
|     |     | 80  | 1.20E-03 | 1.22 | 1.21E-02 | 1.00 | 1.10E-02 | 0.89 |
|     |     | 160 | 5.46E-04 | 1.12 | 5.99E-03 | 1.00 | 5.71E-03 | 0.94 |
|     |     | 320 | 2.60E-04 | 1.06 | 2.98E-03 | 1.00 | 2.91E-03 | 0.97 |
FV-2
|     |     | 20  | 2.92E-03 | –    | 2.17E-03 | –    | 2.26E-03 | –    |
| --- | --- | --- | -------- | ---- | -------- | ---- | -------- | ---- |
|     |     | 40  | 6.76E-04 | 2.10 | 3.65E-04 | 2.57 | 3.94E-04 | 2.51 |
|     |     | 80  | 1.60E-04 | 2.07 | 7.58E-05 | 2.26 | 8.10E-05 | 2.28 |
|     |     | 160 | 3.87E-05 | 2.04 | 1.70E-05 | 2.15 | 1.78E-05 | 2.18 |
|     |     | 320 | 9.50E-06 | 2.02 | 4.02E-06 | 2.08 | 4.13E-06 | 2.11 |
GF
|        |              | 20  | 4.50E-16 | –   | 6.24E-15 | –   | 6.36E-15 | –   |
| ------ | ------------ | --- | -------- | --- | -------- | --- | -------- | --- |
|        |              | 40  | 9.71E-16 | –   | 1.36E-14 | –   | 1.31E-14 | –   |
|        |              | 80  | 1.58E-15 | –   | 3.03E-14 | –   | 3.29E-14 | –   |
|        |              | 160 | 3.38E-15 | –   | 8.62E-14 | –   | 8.67E-14 | –   |
|        |              | 320 | 6.47E-15 | –   | 2.28E-13 | –   | 2.29E-13 | –   |
| 6.3.2. | Lake at rest |     |          |     |          |     |          |     |
In this section, we test the well-balanced property, proven in section 4.6, of the global
flux method for lake at rest solutions of the shallow water system. The problem is set in a
rectangular domain [0,1] × [0,1] with periodic boundary conditions. The initial and exact
| solution | is given       | by         |             |                        |        |     |        |      |
| -------- | -------------- | ---------- | ----------- | ---------------------- | ------ | --- | ------ | ---- |
|          |                | h(x,y)     | = 1−b(x,y), |                        | u(x,y) | =   | v(x,y) | ≡ 0, |
| where    | the bathymetry | is defined | as          |                        |        |     |        |      |
|          |                |            | b(x,y)      | = 0.1sin(2πx)cos(2πy). |        |     |        |      |
In table 6, a convergence study is presented at final time t = 0.1. As expected, thanks to the
f
well-balanced property of the global flux method, the GF is able to achieve machine precision
errors. The standard FV-1 and FV-2 methods show only the classical first and second order
convergence slopes. In figure 15, we present the comparison between the well-balanced GF
| method | and the          | non-well-balanced |       | FV-1 and   | FV-2 | methods. |     |     |
| ------ | ---------------- | ----------------- | ----- | ---------- | ---- | -------- | --- | --- |
| 6.3.3. | 2D supercritical | moving            | water | equilibria |      |          |     |     |
We consider two fully multi-dimensional moving water steady states of the shallow water
system, characterized by constant momentum in supercritical regimes. However, contrary
to the one-dimensional version of such equilibria [28, 26], no exact solution is known for
the simulations presented in this section. For this reason, the numerical results obtained
through the three schemes will be compared qualitatively. The problems are simulated on a
41

|     |     | h   |     |     |     | hu  |      |     | hv  |      |     | norm vel |       |
| --- | --- | --- | --- | --- | --- | --- | ---- | --- | --- | ---- | --- | -------- | ----- |
|     | 1.0 |     |     | 1.0 |     |     | 0.04 | 1.0 |     | 0.04 | 1.0 |          | 0.040 |
1.084
|     |     |     |     |     |     |     | 0.03 |     |     | 0.03 |     |     | 0.036 |
| --- | --- | --- | --- | --- | --- | --- | ---- | --- | --- | ---- | --- | --- | ----- |
1.064
|     | 0.8 |     |     | 0.8   |     |     |      | 0.8 |     |      | 0.8 |     | 0.032 |
| --- | --- | --- | --- | ----- | --- | --- | ---- | --- | --- | ---- | --- | --- | ----- |
|     |     |     |     | 1.044 |     |     | 0.02 |     |     | 0.02 |     |     |       |
0.028
|     | 0.6 |     |     | 1.024 0.6 |     |     | 0.01 | 0.6 |     | 0.01 | 0.6 |     |     |
| --- | --- | --- | --- | --------- | --- | --- | ---- | --- | --- | ---- | --- | --- | --- |
0.024
|     | y   |     |     | 1.004 y |     |     | 0.00 | y   |     | 0.00 | y   |     |     |
| --- | --- | --- | --- | ------- | --- | --- | ---- | --- | --- | ---- | --- | --- | --- |
0.020
|     | 0.4 |     |     | 0.984 0.4 |     |     | 0.01 | 0.4 |     | 0.01 | 0.4 |     |     |
| --- | --- | --- | --- | --------- | --- | --- | ---- | --- | --- | ---- | --- | --- | --- |
0.016
0.964
|     | 0.2 |     |     | 0.2 |     |     | 0.02 | 0.2 |     | 0.02 | 0.2 |     | 0.012 |
| --- | --- | --- | --- | --- | --- | --- | ---- | --- | --- | ---- | --- | --- | ----- |
0.944
|     |     |     |     |     |     |     | 0.03 |     |     | 0.03 |     |     | 0.008 |
| --- | --- | --- | --- | --- | --- | --- | ---- | --- | --- | ---- | --- | --- | ----- |
0.924
|     | 0.0 |     |     | 0.0 |     |     | 0.04 | 0.0 |     | 0.04 | 0.0 |     |     |
| --- | --- | --- | --- | --- | --- | --- | ---- | --- | --- | ---- | --- | --- | --- |
0.0 0.2 0.4 0.6 0.8 1.0 0.904 0.0 0.2 0.4 0.6 0.8 1.0 0.0 0.2 0.4 0.6 0.8 1.0 0.0 0.2 0.4 0.6 0.8 1.0 0.004
|     |     | x   |     |     |     | x   |     |     | x   |     |     | x   |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
(a) FV-1
|     |     | h   |     |           |     | hu  |        |     | hv  |        |     | norm vel |        |
| --- | --- | --- | --- | --------- | --- | --- | ------ | --- | --- | ------ | --- | -------- | ------ |
|     | 1.0 |     |     | 1.092 1.0 |     |     |        | 1.0 |     |        | 1.0 |          | 0.0027 |
|     |     |     |     |           |     |     | 0.0016 |     |     | 0.0016 |     |          |        |
0.0024
|     |     |     |     | 1.068 |     |     | 0.0012 |     |     | 0.0012 |     |     |        |
| --- | --- | --- | --- | ----- | --- | --- | ------ | --- | --- | ------ | --- | --- | ------ |
|     | 0.8 |     |     | 0.8   |     |     |        | 0.8 |     |        | 0.8 |     | 0.0021 |
1.044
|     |     |     |     |           |     |     | 0.0008 |     |     | 0.0008 |     |     | 0.0018 |
| --- | --- | --- | --- | --------- | --- | --- | ------ | --- | --- | ------ | --- | --- | ------ |
|     | 0.6 |     |     | 1.020 0.6 |     |     |        | 0.6 |     |        | 0.6 |     |        |
|     |     |     |     |           |     |     | 0.0004 |     |     | 0.0004 |     |     | 0.0015 |
|     | y   |     |     | y         |     |     |        | y   |     |        | y   |     |        |
|     |     |     |     | 0.996     |     |     | 0.0000 |     |     | 0.0000 |     |     | 0.0012 |
|     | 0.4 |     |     | 0.4       |     |     |        | 0.4 |     |        | 0.4 |     |        |
|     |     |     |     | 0.972     |     |     | 0.0004 |     |     | 0.0004 |     |     | 0.0009 |
|     | 0.2 |     |     | 0.948 0.2 |     |     | 0.0008 | 0.2 |     | 0.0008 | 0.2 |     | 0.0006 |
|     |     |     |     | 0.924     |     |     | 0.0012 |     |     | 0.0012 |     |     |        |
0.0003
|     | 0.0 |     |     | 0.0 |     |     | 0.0016 | 0.0 |     | 0.0016 | 0.0 |     |     |
| --- | --- | --- | --- | --- | --- | --- | ------ | --- | --- | ------ | --- | --- | --- |
0.0 0.2 0.4 0.6 0.8 1.0 0.900 0.0 0.2 0.4 0.6 0.8 1.0 0.0 0.2 0.4 0.6 0.8 1.0 0.0 0.2 0.4 0.6 0.8 1.0 0.0000
|     |     | x   |     |     |     | x   |     |     | x   |     |     | x   |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
(b) FV-2
|     |     | h   |     |           |     | hu  | 1e 145.00 |     | hv  | 1e 144.8 |     | norm vel | 1e 14 |
| --- | --- | --- | --- | --------- | --- | --- | --------- | --- | --- | -------- | --- | -------- | ----- |
|     | 1.0 |     |     | 1.092 1.0 |     |     |           | 1.0 |     |          | 1.0 |          |       |
6.75
|     |     |     |     |           |     |     | 3.75 |     |     | 3.6 |     |     |      |
| --- | --- | --- | --- | --------- | --- | --- | ---- | --- | --- | --- | --- | --- | ---- |
|     |     |     |     | 1.068     |     |     |      |     |     |     |     |     | 6.00 |
|     | 0.8 |     |     | 0.8       |     |     | 2.50 | 0.8 |     | 2.4 | 0.8 |     |      |
|     |     |     |     | 1.044     |     |     |      |     |     |     |     |     | 5.25 |
|     |     |     |     |           |     |     | 1.25 |     |     | 1.2 |     |     |      |
|     | 0.6 |     |     | 1.020 0.6 |     |     |      | 0.6 |     |     | 0.6 |     | 4.50 |
0.00
|     | y   |     |     | y     |     |     |      | y   |     | 0.0 | y       |     | 3.75    |
| --- | --- | --- | --- | ----- | --- | --- | ---- | --- | --- | --- | ------- | --- | ------- |
|     |     |     |     | 0.996 |     |     | 1.25 |     |     |     |         |     |         |
|     | 0.4 |     |     | 0.4   |     |     |      | 0.4 |     |     | 1.2 0.4 |     | 3.00    |
|     |     |     |     | 0.972 |     |     | 2.50 |     |     |     |         |     |         |
|     |     |     |     |       |     |     |      |     |     |     | 2.4     |     | 2 . 2 5 |
|     |     |     |     | 0.948 |     |     | 3.75 |     |     |     |         |     |         |
|     | 0.2 |     |     | 0.2   |     |     |      | 0.2 |     |     | 3.6 0.2 |     | 1 . 5 0 |
|     |     |     |     | 0.924 |     |     | 5.00 |     |     |     |         |     |         |
0.75
|     | 0.0 |     |     | 0.0 |     |     | 6.25 | 0.0 |     |     | 4.8 0.0 |     |     |
| --- | --- | --- | --- | --- | --- | --- | ---- | --- | --- | --- | ------- | --- | --- |
0.0 0.2 0.4 0.6 0.8 1.0 0.900 0.0 0.2 0.4 0.6 0.8 1.0 0.0 0.2 0.4 0.6 0.8 1.0 0.0 0.2 0.4 0.6 0.8 1.0 0.00
|     |     | x   |     |     |     | x   |     |     | x   |     |     | x   |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
(c) GF
Figure 15: Shallow water system: lake at rest. Numerical results for the lake at rest solution on the coarse
| mesh | 40×40 at | final time | t =0.1 | obtained | with | FV-1, FV-2 | and GF. |     |     |     |     |     |     |
| ---- | -------- | ---------- | ------ | -------- | ---- | ---------- | ------- | --- | --- | --- | --- | --- | --- |
f
42

|     |     | h, FV-1 |     |     |     |     |     |     | hv, FV-1 |     |
| --- | --- | ------- | --- | --- | --- | --- | --- | --- | -------- | --- |
|     |     | h, FV-2 |     |     |     |     |     |     | hv, FV-2 |     |
|     |     | h,      | GF  |     |     |     |     |     | hv, GF   |     |
Figure 16: Shallow water system: 2D supercritical equilibria. Numerical results obtained with FV-1, FV-2
| and GF | to steady | state | for N x =N | y =450. |     |     |     |     |     |     |
| ------ | --------- | ----- | ---------- | ------- | --- | --- | --- | --- | --- | --- |
rectangular domain [0,25]×[0,8], and are made fully multi-dimensional by employing a 2D
| bathymetry |     | that is | a function | of  | both x    | and       | y, given | by    |          |     |
| ---------- | --- | ------- | ---------- | --- | --------- | --------- | -------- | ----- | -------- | --- |
|            |     |         |            |     | (cid:18)  |           | (cid:19) |       |          |     |
|            |     |         |            |    | (cid:16)  | (cid:17)2 |          |       |          |     |
|            |     |         |            | 1  |           | r(x,y)    |          |       |          |     |
|            |     |         |            |     | 1−        |           | ,        | where | r(x,y) < | 2   |
|            |     |         |            | 5   |           | 2         |          |       |          |     |
|            |     |         | b(x,y) =   |     |           |           |          |       |          |     |
|            |     |         |            | 0, | elsewhere |           |          |       |          |     |
(cid:112)
with r(x,y) = (x−x )2 +(y −y )2 and (x ,y ) = (10,4). The initial conditions of the
|               |     |           | 0   |           | 0   |           | 0 0 |       |           |      |
| ------------- | --- | --------- | --- | --------- | --- | --------- | --- | ----- | --------- | ---- |
| first problem |     | are given | by  |           |     |           |     |       |           |      |
|               |     | h(x,y,0)  | =   | 2−b(x,y), |     | q (x,y,0) |     | = 24, | q (x,y,0) | = 0. |
|               |     |           |     |           |     | x         |     |       | y         |      |
Inlet boundary conditions (equal to the initial conditions) are imposed on the left boundary
of the domain, and outlet (homogeneous Neumann) on the right. Top and bottom of the
domain are periodic boundaries. In figure 16, we present the numerical solutions for the
conservative variables when the numerical steady state is reached (time residual close to
machine precision). All simulations are performed on a mesh of 450 × 450 elements. GF
is able to capture and resolve sharply the many shocks appearing behind the bathymetry
bump. Although the mesh resolution for this case is quite fine, FV-1 still presents a highly
diffused result. An improvement is experienced when using the linear reconstruction for
FV-2, although all the waves still appear as smooth transitions. However, the results of the
classical methods still remain significantly inferior to those of GF, which captures all waves
| sharply | in much | fewer | cells. |     |     |     |     |     |     |     |
| ------- | ------- | ----- | ------ | --- | --- | --- | --- | --- | --- | --- |
In figure 17, we show a perturbation of the numerical equilibrium obtained above. We
| add a small |     | drop of | water shaped |         | as  |              |     |            |     |     |
| ----------- | --- | ------- | ------------ | ------- | --- | ------------ | --- | ---------- | --- | --- |
|             |     |         |              |         |     | 10−4e−(x−16) |     | 2 + (y−3)2 |     |     |
|             |     |         |              | δh(x,y) |     | =            |     |            |     |     |
0 .8 2
43

h−h , FV-1 N = 150,N = 50
eq x y
8
7
6
5 4 3
2
1
0 0 5 10 x 15 20 25
y
h 8
7
6
5 4 3
2
1
0 0 5 10 x 15 20 25
y
hu 8
7
6
5 4 3
2
1
0 0 5 10 x 15 20 25
y
hv 8
7
6
5 4 3
2
1
0 0 5 10 x 15 20 25
y
norm vel 0.0001110 0.0001110 0.0001110 0.000288
0.0000886 0.0000886 0.0000886 0.000252
0.0000661 0.0000661 0.0000661 0.000216
0.0000437 0.0000437 0.0000437 0.0000212 0.0000212 0.0000212 0.000180 0.0000012 0.0000012 0.0000012 0.000144 0.0000237 0.0000237 0.0000237 0.000108
0.0000461 0.0000461 0.0000461 0.000072
0.0000686 0.0000686 0.0000686 0.000036
0.0000910 0.0000910 0.0000910 0.000000 h−h , FV-2 N = 150,N = 50 eq x y
8 7
6 5 4 3
2
1
0 0 5 10 x 15 20 25
y
h 8 7
6 5 4 3
2
1
0 0 5 10 x 15 20 25
y
hu 8 7
6 5 4 3
2
1
0 0 5 10 x 15 20 25
y
hv 8 7
6 5 4 3
2
1
0 0 5 10 x 15 20 25
y
h−h , GF N = 150,N = 50
eq x y
8
7
6
5 4 3
2
1
0 0 5 10 x 15 20 25
norm vel 0.0001110 0.0001110 0.0001110 0.0016 0.0000886 0.0000886 0.0000886 0.0014
0.0000661 0.0000661 0.0000661 0.0012 0.0000437 0.0000437 0.0000437 0.0010 0.0000212 0.0000212 0.0000212 0.0000012 0.0000012 0.0000012 0.0008 0.0000237 0.0000237 0.0000237 0.0006
0.0000461 0.0000461 0.0000461 0.0004
0.0000686 0.0000686 0.0000686 0.0002
0.0000910 0.0000910 0.0000910 0.0000
y
h 8
7
6
5 4 3
2
1
0 0 5 10 x 15 20 25
y
hu 8
7
6
5 4 3
2
1
0 0 5 10 x 15 20 25
y
hv 8
7
6
5 4 3
2
1
0 0 5 10 x 15 20 25
y
norm vel 0.0001110 0.0001110 0.0001110 0.000675
0.0000886 0.0000886 0.0000886 0.000600
0.0000661 0.0000661 0.0000661 0.000525
0.0000437 0.0000437 0.0000437 0.000450 0.0000212 0.0000212 0.0000212 0.000375 0.0000012 0.0000012 0.0000012 0.000300 0.0000237 0.0000237 0.0000237 0.000225
0.0000461 0.0000461 0.0000461 0.000150
0.0000686 0.0000686 0.0000686 0.000075
0.0000910 0.0000910 0.0000910 0.000000 h−h , FV-1 N = 450,N = 450 eq x y
8 7
6 5 4 3
2
1
0 0 5 10 x 15 20 25
y
h 8 7
6 5 4 3
2
1
0 0 5 10 x 15 20 25
y
hu 8 7
6 5 4 3
2
1
0 0 5 10 x 15 20 25
y
hv 8 7
6 5 4 3
2
1
0 0 5 10 x 15 20 25
y
8
7 6
5
4
norm vel 3 0.0001110 0.0001110 0.0001110 0.00048 0.0000886 0.0000886 0.0000886
0.0000661 0.0000661 0.0000661 0.00042 2 0.0000437 0.0000437 0.0000437 0.00036 0.0000212 0.0000212 0.0000212 0.00030 0.0000012 0.0000012 0.0000012 0.00024 1 0.0000237 0.0000237 0.0000237 0.00018
0.0000461 0.0000461 0.0000461 0.00012
0 0.0000686 0.0000686 0.0000686 0.00006
0 5 10 15 20 25 0.0000910 0.0000910 0.0000910 0.00000 x
y
h 8
7 6
5
4
3
2 1
0
0 5 10 15 20 25 x
y
hu 8
7 6
5
4
3
2 1
0
0 5 10 15 20 25 x
y
hv 8
7 6
5
4
3
2 1
0
0 5 10 15 20 25 x
y
norm vel
0.0001110 0.0001110 0.0001110
0.00048
0.0000886 0.0000886 0.0000886 0.00042 0.0000661 0.0000661 0.0000661
0.00036 0.0000437 0.0000437 0.0000437
0.0000212 0.0000212 0.0000212 0.00030
0.0000012 0.0000012 0.0000012 0.00024 0.0000237 0.0000237 0.0000237 0.00018
0.0000461 0.0000461 0.0000461 0.00012 0.0000686 0.0000686 0.0000686 0.00006
0.0000910 0.0000910 0.0000910
0.00000
Figure17: Shallowwatersystem: perturbationof2Dsupercriticalequilibria. Differencebetweenthesolution
and the equilibrium for different methods on the mesh N = 150, N = 50 and reference solution with FV
x y
on mesh N =450, N =450
x y
on the h variable and we continue the simulation until t = 0.4. We plot the difference
f
between each numerical equilibrium and solutions for all the methods on a coarse grid N =
x
150, N = 50, whilewerunafinertestfortheFV-1testcaseonthegridN = 450, N = 450.
y x y
Looking at the simulation of FV-1 on the coarse mesh, we barely see the perturbation as
the method is too dissipative. For the FV-2 simulations, the perturbation is visible, but
small oscillations are present due to the fact that the method does not reach proper time
convergence (the time residual stays around 10−6). On the other hand, the GF nicely shows
the perturbation moving towards the right without any spurious oscillations.
As a second problem, and in order to test the robustness of the method, we also consider
a crooked supercritical equilibrium with the same bathymetry but different initial conditions:
h(x,y,0) = 2−b(x,y), q (x,y,0) = 24, q (x,y,0) = 4π.
x y
In this case, left and bottom boundaries are inlet boundaries, while right and top are outlets.
It should be noticed, that this test case is even more challenging than the one shown before
since no part of the fluid is aligned with the background Cartesian mesh. In figure 18, the
resultsobtainedfortheconservativevariablesarepresented, wherethesameconclusionabout
the quality of the result of the GF method can be drawn as for the previous test. All the
physical features of the equilibria are well captured, while they are significantly more diffused
by the FV-1 and FV-2 methods.
7. Conclusions and perspectives
In this work, we have presented a new way to derive finite volume methods for nonlinear
multi-dimensional hyperbolic systems, which is based on the global flux approach (9), in-
troduced in [13]. It is a general way to obtain stationarity preserving schemes for nonlinear
problems. Besides its generality, the method is also able to achieve super-convergence on
steady problems, with error reductions of one or two orders of magnitude compared to a
standard second order finite volume approach. Despite a focus on stationary states during
its design, we observe remarkably high resolution for unsteady multi-dimensional problems,
44

h, FV-1 hv, FV-1
h, FV-2 hv, FV-2
h, GF hv, GF
Figure18: Shallowwatersystem: 2Dcrookedsupercriticalequilibria. NumericalresultsobtainedwithFV-1,
FV-2 and GF to steady state for N =N =450.
x y
outperforming standard first and even second-order finite volume methods, for a large span
of Mach/Froude numbers.
This work opens the way to several future developments. In particular, the extension
of the finite volume formulation to high order methods by using high-degree polynomial re-
construction techniques like WENO [35] is a natural next step, following the work on the
1D global flux WENO approach introduced in [28]. Moreover, the first order finite volume
method can also be seen as the starting point to develop a new family of multi-dimensional
high order discontinuous Galerkin methods based on the global flux formulation. More
investigations will also be dedicated to the extension of the method to deal with mathem-
atical models characterized by curl-free solutions, like the Maxwell equations. Extending,
for instance, the observation that stationarity preserving methods are also low Mach number
compliant, theoreticalworkwillincludefurtheranalysisofthemethodinunsteadysituations.
References
[1] R. Abgrall, P.-H. Maire, and M. Ricchiuto. Embedding general conservation constraints
in discretizations of hyperbolic systems on arbitrary meshes: A multidimensional frame-
work, 2025.
[2] R.AbgrallandM.Ricchiuto. HighordermethodsforCFD. InRenedeBorstErwinStein
and T. J.R. Hughes, editors, Encyclopedia of Computational Mechanics, Second Edition.
John Wiley and Sons, 2017.
[3] R. Abgrall and M. Ricchiuto. Hyperbolic balance laws: Residual distribution, local and
global fluxes. In Dia Zeidan, Jochen Merker, Eric Goncalves Da Silva, and Lucy T.
Zhang, editors, Numerical Fluid Dynamics: Methods and Computations, pages 177–222.
Springer Nature Singapore, Singapore, 2022.
45

[4] Remi Abgrall, Yongle Liu, and Mario Ricchiuto. Positivity-preserving well-balanced
pampa schemes with global flux quadrature for one-dimensional shallow water models,
2025.
[5] E. Audusse, François Bouchut, Marie-Odile Bristeau, Rupert Klein, and Benoıt Per-
thame. A fast and stable well-balanced scheme with hydrostatic reconstruction for
shallow water flows. SIAM Journal on Scientific Computing, 25(6):2050–2065, 2004.
[6] E. Audusse, Minh Hieu Do, Pascal Omnes, and Yohan Penel. Analysis of modified
Godunovtypeschemesforthetwo-dimensionallinearwaveequationwithCoriolissource
term on cartesian meshes. Journal of Computational Physics, 373:91–129, 2018.
[7] E. Audusse, V. Dubos, N. Gaveau, and Y. Penel. Energy-stable and linearly well-
balanced numerical schemes for the nonlinear shallow water equations with the coriolis
force. SIAM Journal on Scientific Computing, 47(1):A1–A23, 2025.
[8] D.S. Balsara. Multidimensional HLLE Riemann solver: Application to Euler and mag-
| netohydrodynamic |     | flows. | J.  | Comput. | Phys., | 229:1970–1993, |     | 2010. |
| ---------------- | --- | ------ | --- | ------- | ------ | -------------- | --- | ----- |
[9] W. Barsukow. Low Mach number finite volume methods for the acoustic and Euler
| equations. | Doctoral | thesis, |     | University | of  | Wuerzburg, | 2018. |     |
| ---------- | -------- | ------- | --- | ---------- | --- | ---------- | ----- | --- |
[10] W. Barsukow. Stationarity preserving schemes for multi-dimensional linear systems.
| Mathematics |     | of Computation, |     | 88(318):1621–1645, |     |     | 2019. |     |
| ----------- | --- | --------------- | --- | ------------------ | --- | --- | ----- | --- |
[11] W. Barsukow. Truly multi-dimensional all-speed schemes for the euler equations on
cartesian grids. Journal of Computational Physics, 435:110216, 2021.
[12] W. Barsukow, Raphael Loubère, and Pierre-Henri Maire. A node-conservative vorticity
preservingfinitevolumemethodforlinearacousticsonunstructuredgrids. Math. Comp.,
| 94:2299–2343, |     | 2025. |     |     |     |     |     |     |
| ------------- | --- | ----- | --- | --- | --- | --- | --- | --- |
[13] W. Barsukow, M. Ricchiuto, and D. Torlo. Structure preserving nodal continuous fi-
nite elements via global flux quadrature. Numerical Methods for Partial Differential
| Equations, | 41(1):e23167, |     | 2025. |     |     |     |     |     |
| ---------- | ------------- | --- | ----- | --- | --- | --- | --- | --- |
[14] Wasilij Barsukow. Stationarity preservation and the low mach number behaviour of the
| discontinuous |     | galerkin | method | on  | cartesian | grids, | 2025. |     |
| ------------- | --- | -------- | ------ | --- | --------- | ------ | ----- | --- |
[15] Wasilij Barsukow, Mario Ricchiuto, and Davide Torlo. Stationarity preserving nodal
finite element methods for multi-dimensional linear hyperbolic balance laws via a global
| flux quadrature |     | formulation, |     | 2025. |     |     |     |     |
| --------------- | --- | ------------ | --- | ----- | --- | --- | --- | --- |
[16] A. Bermudez and M.E. Vazquez. Upwind methods for hyperbolic conservation laws with
| source | terms. | Computers | &   | Fluids, | 23(8):1049 |     | – 1071, 1994. |     |
| ------ | ------ | --------- | --- | ------- | ---------- | --- | ------------- | --- |
[17] C. Berthon and C. Chalons. A fully well-balanced, positive and entropy-satisfying
Godunov-type method for the shallow-water equations. Mathematics of Computation,
| 85(299):1281–1307, |     | 2016. |     |     |     |     |     |     |
| ------------------ | --- | ----- | --- | --- | --- | --- | --- | --- |
46

[18] W. Boscheri and M. Dumbser. Arbitrary-Lagrangian-Eulerian one-step WENO finite
volume schemes on unstructured triangular meshes. Communications in Computational
| Physics, | 14(5):1174–1206, | 2013. |     |     |
| -------- | ---------------- | ----- | --- | --- |
[19] F. Bouchut, J. Le Sommer, and V. Zeitlin. Frontal geostrophic adjustment and nonlin-
ear wave phenomena in one-dimensional rotating shallow water. Part 2. high-resolution
numerical simulations. Journal of Fluid Mechanics, 514:35–63, 2004.
[20] David L Brown. Performance of under-resolved two-dimensional incompressible flow
simulations. Journal of Computational Physics, 122(1):165–183, 1995.
[21] V. Caselles, R. Donat, and G. Haro. Flux-gradient and source-term balancing for certain
high resolution shock-capturing schemes. Computers & Fluids, 38(1):16–36, 2009.
[22] M.J. Castro and C. Parés. Well-balanced high-order finite volume methods for systems
| of balance | laws. Journal | of Scientific | Computing, | 82(2):48, 2020. |
| ---------- | ------------- | ------------- | ---------- | --------------- |
[23] Y. Cheng, A. Chertock, M. Herty, A. Kurganov, and T. Wu. A new approach for
designing moving-water equilibria preserving schemes for the shallow water equations.
| Journal | of Scientific | Computing, | 80(1):538–554, | 2019. |
| ------- | ------------- | ---------- | -------------- | ----- |
[24] A. Chertock, S. Cui, A. Kurganov, Ş.N. Özcan, and E. Tadmor. Well-balanced schemes
for the Euler equations with gravitation: Conservative formulation using global fluxes.
| Journal | of Computational | Physics, | 358:36–52, | 2018. |
| ------- | ---------------- | -------- | ---------- | ----- |
[25] A. Chertock, A. Kurganov, Xin Liu, Yongle Liu, and T. Wu. Well-balancing via flux
globalization: Applications to shallow water equations with wet/dry fronts. Journal of
| Scientific | Computing, | 90(1):1–21, | 2022. |     |
| ---------- | ---------- | ----------- | ----- | --- |
[26] M. Ciallella, L. Micalizzi, V. Michel-Dansac, P. Öffner, and D. Torlo. A high-order, fully
well-balanced, unconditionally positivity-preserving finite volume framework for flood
simulations. GEM-International Journal on Geomathematics, 16(1):1–33, 2025.
[27] M. Ciallella, L. Micalizzi, P. Öffner, and D. Torlo. An arbitrary high order and positivity
preserving method for the shallow water equations. Computers & Fluids, 247:105630,
2022.
[28] M. Ciallella, D. Torlo, and M. Ricchiuto. Arbitrary high order WENO finite volume
scheme with flux globalization for moving equilibria preservation. Journal of Scientific
| Computing, | 96(2):53, | 2023. |     |     |
| ---------- | --------- | ----- | --- | --- |
[29] H. Deconinck and M. Ricchiuto. Residual distribution schemes: Foundations and ana-
lysis. In Encyclopedia of Computational Mechanics Second Edition, pages 1–53. John
| Wiley | & Sons, Ltd, | 2017. |     |     |
| ----- | ------------ | ----- | --- | --- |
[30] M. Dumbser D.S. Balsara and R. Abgrall. Multidimensional HLLC Riemann Solver for
Unstructured Meshes - With Application to Euler and MHD Flows. J. Comput. Phys.,
| 261:172–208, | 2014. |     |     |     |
| ------------ | ----- | --- | --- | --- |
47

[31] E. Gaburro, M. Ricchiuto, and M. Dumbser. On general and complete multidimensional
riemann solvers for nonlinear systems of hyperbolic conservation laws, 2025. arXiv
| eprint, | 2506.00207, math.na. |     |     |     |     |
| ------- | -------------------- | --- | --- | --- | --- |
[32] L.L. Gascón and J.M. Corberán. Construction of second-order TVD schemes for
nonhomogeneous hyperbolic conservation laws. Journal of Computational Physics,
| 172(1):261–297, | 2001. |     |     |     |     |
| --------------- | ----- | --- | --- | --- | --- |
[33] Laurent Gosse and Giuseppe Toscani. Space localization and well-balanced schemes
for discrete kinetic models in diffusive regimes. SIAM Journal on Numerical Analysis,
| 41(2):641–658, | 2003. |     |     |     |     |
| -------------- | ----- | --- | --- | --- | --- |
[34] R. Jeltsch and M. Torrilhon. On curl-preserving finite volume discretizations for shallow
| water equations. | BIT | Numerical | Mathematics, | 46:35–53, | 2006. |
| ---------------- | --- | --------- | ------------ | --------- | ----- |
[35] G.-S.JiangandC.-C.Wu. Ahigh-orderWENOfinitedifferenceschemefortheequations
of ideal magnetohydrodynamics. Journal of Computational Physics, 150(2):561–594,
1999.
[36] J. Jung and V. Perrier. Steady low mach number flows: identification of the spurious
| mode and | filtering method. | J.  | Comput. Phys., | 468:111462, | 2022. |
| -------- | ----------------- | --- | -------------- | ----------- | ----- |
[37] J. Jung and V. Perrier. Behavior of the discontinuous galerkin method for compress-
ible flows at low mach number on triangles and tetrahedrons. SIAM J. Sci. Comput.,
| 46(1):A452–A482, | 2024. |     |     |     |     |
| ---------------- | ----- | --- | --- | --- | --- |
[38] M. Kazolea, C. ParésMadroñal, and M. Ricchiuto. Approximate well-balanced WENO
finite difference schemes using a global-flux quadrature method with multi-step ode
| integrator | weights. arXiv | preprint | arXiv:2501.06155, |     | 2025. |
| ---------- | -------------- | -------- | ----------------- | --- | ----- |
[39] Giovanni L., R. Andrassy, W. Barsukow, J. Higl, P.V.F. Edelmann, and F.K. Röpke.
Performance of high-order godunov-type methods in simulations of astrophysical low
| mach number | flows. Astronomy |     | & Astrophysics, | 686:A34, | 2024. |
| ----------- | ---------------- | --- | --------------- | -------- | ----- |
[40] TB Lung and PL Roe. Toward a reduction of mesh imprinting. International Journal
| for Numerical | Methods | in Fluids, | 76(7):450–470, | 2014. |     |
| ------------- | ------- | ---------- | -------------- | ----- | --- |
[41] Y. Mantri, P. Öffner, and M. Ricchiuto. Fully well-balanced entropy controlled discon-
tinuous galerkin spectral element method for shallow water flows: global flux quadrature
and cell entropy correction. Journal of Computational Physics, 498:112673, 2024.
[42] L. Micalizzi, M. Ricchiuto, and R. Abgrall. Novel well-balanced continuous interior
penalty stabilizations. Journal of Scientific Computing, 100(1):14, 2024.
[43] V. Michel-Dansac, C. Berthon, S. Clain, and F. Foucher. A well-balanced scheme for the
shallow-water equations with topography or manning friction. Journal of Computational
| Physics, | 335:115–154, | 2017. |     |     |     |
| -------- | ------------ | ----- | --- | --- | --- |
48

[44] V. Michel-Dansac, C. Berthon, Stéphane Clain, and Françoise Foucher. A two-
dimensional high-order well-balanced scheme for the shallow water equations with topo-
graphy and manning friction. Computers & Fluids, 230:105152, 2021.
[45] S. Mishra and E. Tadmor. Constraint preserving schemes using potential-based fluxes
I. multidimensional transport equations. Communications in Computational Physics,
9(3):688–710, 2011.
[46] KeithWilliamMortonandPhilipLRoe. Vorticity-preservinglax–wendroff-typeschemes
for the system wave equation. SIAM Journal on Scientific Computing, 23(1):170–192,
2001.
[47] H. Nessyahu and E. Tadmor. Non-oscillatory central differencing for hyperbolic conser-
vation laws. Journal of computational physics, 87(2):408–463, 1990.
[48] M. Ricchiuto. Contributions to the development of residual discretizations for hyperbolic
conservation laws with application to shallow water flows. Habilitation à diriger des
recherches, Université Sciences et Technologies - Bordeaux I, December 2011. [pdf].
[49] M. Ricchiuto and A. Bollermann. Stabilized residual distribution for shallow water
simulations. Journal of Computational Physics, 228(4):1071–1115, 2009.
[50] M.RicchiutoandD.Torlo. Analyticaltravellingvortexsolutionsofhyperbolicequations
for validating very high order schemes. arXiv preprint arXiv:2109.10183, 2021.
[51] P.L. Roe. Upwind differencing schemes for hyperbolic conservation laws with source
terms. In Claude Carasso, Denis Serre, and Pierre-Arnaud Raviart, editors, Nonlinear
Hyperbolic Problems, pages 41–51, Berlin, Heidelberg, 1987. Springer Berlin Heidelberg.
[52] D. Sidilkover. Factorizable schemes for the equations of fluid flow. Applied numerical
mathematics, 41(3):423–436, 2002.
[53] J. Zhang, Y. Xia, and Y. Xu. Well-balanced discontinuous galerkin method with flux
globalization for rotating shallow water equations. Journal of Computational Physics,
page 114094, 2025.
49