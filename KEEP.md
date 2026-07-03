Structure-preserving schemes conserving entropy and kinetic energy
Kunal Bahugunaa,1,∗, Ramesh Kollurua,2, S.V. Raghurama Raoa
aDepartment of Aerospace Engineering, Indian Institute of Science, Bangalore-560012, India
Abstract
This paper presents a novel structure-preserving scheme for Euler equations, focusing on the numerical con-
servationofentropyandkineticenergy. Explicitfluxfunctionsengineeredtoconserveentropyareintroduced
within the finite-volume framework. Further, discrete kinetic energy conservation too is introduced. A
systematic inquiry is presented, commencing with an overview of numerical entropy conservation and formu-
lation of entropy-conserving and kinetic energy-preserving fluxes, followed by the study of their properties
and efficacy. A novelty introduced is to associate numerical entropy conservation to the discretization of
the energy conservation equation. Furthermore, an entropy-stable shock-capturing diffusion method and a
hybridapproachutilizingtheentropydistancetomanagesmoothregionseffectivelyarealsointroduced. The
addition of artificial viscosity in apprproiate regions ensures entropy generation sufficient to prevent numer-
ical instabilities. Various test cases, showcasing the efficacy and stability of the proposed methodology, are
presented.
1. Introduction
Thecurrentthrustinalgorithmdevelopmentisdirectedtowardsdevisingsuitablenumericalmethodsthat
adheretotheadditionalconstraintsinherentinEulerequations,suchaspreservingentropyandkineticenergy.
Owingtotheirnon-linearhyperbolicnature,Eulerequationspermitdiscontinuoussolutionslikeshockwaves
andcontactdiscontinuities. Thesesolutionslackuniqueness,promptingthepursuitofsolutionsthatadhereto
theEulerequationsinaweaksense. ThisisincontrasttothatofmodellingwiththeNavier-Stokesequations,
where physical viscosity guarantees adherence to the second law of thermodynamics, yielding physically
meaningful solutions. Numerical methods devised for Euler equations must also ensure the conservation or
productionofnumericalentropy, therebyconvergingtowardsthecorrectsolution. Furthermore, thestability
of numerical methods for non-linear partial differential equations hinges on entropy stability, rendering it a
highly desirable property.
Tadmor [1; 2] introduced the mathematical conditions required for numerical entropy conservation and
stabilityforsystemsofhyperbolicequations,withintheframeworkofthefinite-volumemethod. Theresulting
conservative fluxes of entropy did not have an explicit form and required an expensive numerical quadrature
for calculations [3]. Roe [4] introduced affordable and explicit entropy conservative flux functions for Euler
equations. Ranocha [5] introduced various parameterization of entropy variables in different ways. All of
these fluxes are found to be second-order accurate spatially. LeFloch and Rohde [6] introduced arbitrary
higher-order entropy conservative fluxes.
Ensuring accurate preservation of discrete kinetic energy is another crucial objective for numerical flux
functions, particularly in flows dominated by kinetic energy [7; 8]. Jameson [9] established conditions on
∗Correspondingauthor
Email addresses: kunalb@iisc.ac.in(KunalBahuguna),kollurur@alum.iisc.ac.in(RameshKolluru),
raghu@iisc.ac.in(S.V.RaghuramaRao)
1ResearchScholar
2Presentaffilication: SeniorCFDConsultant, BosonQPSI(BQP),IndiaOffice, 2743, 15thcross, 27thmainroad, HSR1st
sector,Bangalore-560102
Preprint submitted to Elsevier May 20, 2025
5202
yaM
91
]AN.htam[
1v47331.5052:viXra

the numerical momentum flux to maintain kinetic energy preservation, and a slightly modified version was
presentedbyRanocha[10]. Combiningthiswithentropyconservation,Chandrashekarintroducedanentropy-
conserving and kinetic energy-preserving flux [11], further addressing kinetic energy stability to dissipate
discretekineticenergyakintotheeffectofphysicalviscosityinNavier-Stokesequations. Entropyconservative
| schemes | can also be | obtained | through | an optimization | problem | [12]. |     |
| ------- | ----------- | -------- | ------- | --------------- | ------- | ----- | --- |
In the numerical simulation of viscous flows too, additional diffusion becomes necessary due to the in-
adequacy of physical viscosity in yielding oscillation-free solutions on coarse grids. However, any artificial
diffusion introduced must invariably support entropy generation to prevent nonphysical solutions from en-
tropy violations. Notably, improper entropy generation, as observed in Roe’s approximate Riemann solver
[13], leads to issues such as expansion shocks and shock instabilities such as carbuncle phenomenon [14].
To circumvent such anomalies, various schemes, especially the family of Riemann solvers, need to explicitly
ensure the satisfaction of entropy inequalities, thereby averting numerical aberrations. Ismail and Roe [15]
proposed a stabilizing matrix diffusion method to achieve a ‘consistent’ entropy across shocks. Nevertheless,
whilenumericaldiffusiontailoredforshockwavesproveseffective,itmightnotbesuitableforsmoothregions.
Consequently, many schemes adopt a hybrid approach in which numerical diffusion is adapted based on flow
gradients [16; 17]. Furthermore, achieving higher-order entropy stability is feasible through scaled entropy
| variable | reconstruction | or sign-preserving |     | ENO | reconstruction | [18]. |     |
| -------- | -------------- | ------------------ | --- | --- | -------------- | ----- | --- |
Inthiswork,anovelmethodologyispresented,combiningentropyconservationandkineticenergypreser-
vation strategies. a special feature of Euler equations, that the energy equation actually contains an entropy
conservation equation and a kinetic energy conservation equation, is exploited by associating numerical en-
tropy conservative fluxes with the energy equation. Further, entropy stability is ensured by appropriately
addingnumericaldiffusionasrequired. Thestructureofthispaperisasfollows: InSection2,aconciseintro-
duction is given to numerical entropy conservation. In Section 3, a diffusion-based strategy for conservation,
detailing the formulation of entropy-conserving and kinetic-energy-preserving numerical fluxes is outlined.
The numerical investigation of these fluxes, focusing on aspects like numerical entropy conservation, com-
putational efficiency and accurate contact discontinuity capturing, is presented in Section 4, including the
demonstration of the experimental order of convergence. Then, an entropy-stable shock-capturing numerical
diffusion method based on Rankine-Hugoniot conditions, highlighting its ability to generate appropriate en-
tropy, is introduced in Section 5. Addressing the insufficiency of this numerical diffusion in smooth regions,
a hybrid scheme based on entropy distance is proposed in Section 6.. In Section 7, the numerical results for
one-dimensional and two-dimensional benchmark test cases for Euler equations are presented. Finally, the
| conclusions | drawn from | the | work are | presented | in Section | 8.  |     |
| ----------- | ---------- | --- | -------- | --------- | ---------- | --- | --- |
2. Structures within Euler Equations: Entropy and Kinetic Energy
Governing equations for compressible inviscid flow are the unsteady Euler equations, with the closure
based on the assumption of a perfect gas. Focusing on the one-dimensional Euler equations with conserved
variable vector U=[ρ,ρu,ρE]T and the flux vector F(U)=[ρu,ρu2+p,ρuE+pu]T, they can be expressed
| as conserving | mass, | momentum | and | energy, as | in (1). |              |      |
| ------------- | ----- | -------- | --- | ---------- | ------- | ------------ | ---- |
|               |       |          |     | ∂U ∂F(U)   |         |              |      |
|               |       |          |     | +          | =0, or, | individually | (1a) |
|               |       |          |     | ∂t ∂x      |         |              |      |
|               |       |          |     |            | ∂ρ      | ∂            |      |
|               |       |          |     |            |         | + (ρu)=0     | (1b) |
|               |       |          |     |            | ∂t      | ∂x           |      |
|               |       |          |     |            | ∂ ∂     |              |      |
(ρu2+p)=0
|     |     |     |     |     | (ρu)+            |     | (1c) |
| --- | --- | --- | --- | --- | ---------------- | --- | ---- |
|     |     |     |     |     | ∂t ∂x            |     |      |
|     |     |     |     | ∂   | ∂                |     |      |
|     |     |     |     |     | (ρE)+ (ρuE+pu)=0 |     | (1d) |
|     |     |     |     | ∂t  | ∂x               |     |      |
=e+1u2
Here ρ is the density, u is the velocity, p is the thermodynamic pressure and E is the total specific
2
(internal + kinetic) energy. The system is closed using the ideal gas law p = (γ − 1)ρe where e is the
2

specific internal energy, and γ is the ratio of specific heats of the gas. This system of equations is hyperbolic
in nature, with the flux Jacobian matrix A˜ = ∂F(U)/∂U having real and distinct eigenvalues. Since the
system is non-linear, Euler equations may generate discontinuous solutions even when the initial conditions
are smooth. These solutions are also not unique, and it is important for the numerical schemes to converge
to the right physically admissible solutions (i.e., the vanishing viscosity solutions).
2.1. Entropy inequality, entropy variables, symmetric form and entropy potential functions
SincethesystemofEulerequationsishyperbolic,itadmitsanadditionalconservationlaw[19],represent-
ing the conservation of entropy in smooth regions. Consider a convex entropy function η(U). Multiplying
| the Euler | equations | by its | derivative, | we  | obtain |          |          |     |     |
| --------- | --------- | ------ | ----------- | --- | ------ | -------- | -------- | --- | --- |
|           |           |        |             |     |        | (cid:20) | (cid:21) |     |     |
|           |           |        |             |     | ∂η(U)  | ∂U ∂F(U) |          |     |     |
|           |           |        |             |     |        | +        | =0       |     | (2) |
|           |           |        |             |     | ∂U     | ∂t       | ∂x       |     |     |
or
|           |           |               |     |           | ∂η(U)   |         | ∂U     |     |     |
| --------- | --------- | ------------- | --- | --------- | ------- | ------- | ------ | --- | --- |
|           |           |               |     |           |         | ′ ′     |        |     |     |
|           |           |               |     |           |         | +η (U)F | (U) =0 |     |     |
|           |           |               |     |           | ∂t      |         | ∂x     |     |     |
| If we now | introduce | a consistency |     | condition | as      |         |        |     |     |
|           |           |               |     |           | ′       | ′ (U)TF | ′      |     |     |
|           |           |               |     |           | ζ (U)=η |         | (U)    |     | (3) |
we obtain
|     |     |     |     |     | ∂η(U) |          | ∂U  |     |     |
| --- | --- | --- | --- | --- | ----- | -------- | --- | --- | --- |
|     |     |     |     |     |       | +ζ ′ (U) | =0  |     |     |
|     |     |     |     |     | ∂t    |          | ∂x  |     |     |
or
|     |     |     |     |     | ∂η(U) | ∂ζ(U) |     |     |     |
| --- | --- | --- | --- | --- | ----- | ----- | --- | --- | --- |
|     |     |     |     |     |       | +     | =0  |     | (4) |
|     |     |     |     |     |       | ∂t ∂x |     |     |     |
which is the entropy conservation equation, with ζ(U) being the entropy flux function. The above entropy
conservation equation is valid only in smooth regions of the flow. If we include the regions of discontinuities
| (shock waves), | then       | we have | the         | inequality |            |          |          |     |     |
| -------------- | ---------- | ------- | ----------- | ---------- | ---------- | -------- | -------- | --- | --- |
|                |            |         |             |            | ∂η(U)      | ∂ζ(U)    |          |     |     |
|                |            |         |             |            |            | +        | ≤0       |     | (5) |
|                |            |         |             |            |            | ∂t ∂x    |          |     |     |
| For 1-D Euler  | equations, |         | the entropy | pair       | fulfilling | (3) is   | given by |     |     |
|                |            |         |             |            |            | ρs       | ρus      |     |     |
|                |            |         |             | η(U)=−     |            | , ζ(U)=− |          | ,   |     |
|                |            |         |             |            |            | γ−1      | γ−1      |     |     |
with s = ln(p/ργ), representing the specific entropy. Interestingly, this entropy pair is valid not only for
Euler equations but also for Navier-Stokes equations. Mock [20] showed that the Euler equations with an
entropy equation implies the existence of entropy variable vector V with a one-to-one mapping to U, and is
given by
|     |     |     |     |     |       | (cid:20) |        | (cid:21)T |     |
| --- | --- | --- | --- | --- | ----- | -------- | ------ | --------- | --- |
|     |     |     |     |     | dη(U) | γ−s      | ρu2 ρu | ρ         |     |
|     |     |     |     | V=  |       | = −      | , ,−   |           |     |
|     |     |     |     |     | dU    | γ−1      | 2p p   | p         |     |
Entropyvariablesymmetrizesthehyperbolicsystem(1),i.e.,changingthevariablesfromUtoVintheEuler
| equations | gives the | following | symmetric |     | system |       |     |     |     |
| --------- | --------- | --------- | --------- | --- | ------ | ----- | --- | --- | --- |
|           |           |           |           |     |        | ∂V ∂V |     |     |     |
A˜ +B˜
|     |     |     |     |     |     |       | =0  |     | (6) |
| --- | --- | --- | --- | --- | --- | ----- | --- | --- | --- |
|     |     |     |     |     |     | ∂t ∂x |     |     |     |
where A˜ is a symmetric positive definite matrix and B˜ is a symmetric matrix. Further, entropy variable
facilitates the introduction of entropy potential function ϕ and entropy flux potential function ψ (together
| called as entropy |     | pair) as |           |     |     |     |                     |     |     |
| ----------------- | --- | -------- | --------- | --- | --- | --- | ------------------- | --- | --- |
|                   |     | ϕ=VT     | ·U−η(U)=ρ |     |     | and | ψ =VT ·F(U)−ζ(U)=ρu |     | (7) |
Thisformulationofentropyinequalitytogetherwiththehyperbolicsystem, andthustheconceptofentropic
solutions, extends to arbitrary systems of conservation laws in n dimensions, as elucidated by Harten [21].
3

2.2. Structures in the total energy equation: entropy conservation and kinetic energy preservation
ThetwospecificstructureswithintheEulerequationsworthpreservinginanumericalschemeare: entropy
conservationandkineticenergypreservation. Bothofthemarepartsofthetotalenergyconservationequation
| and this | fact is exploited | in  | our | work. |     |     |     |     |     |     |     |     |     |
| -------- | ----------------- | --- | --- | ----- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
Letusfirstobtainthekineticenergyequationseparately; wemultiply(1c)byu,(1b)by−u2/2,andthen
| sum them | to obtain: |     |     |     |                      |     |                      |     |     |     |     |     |     |
| -------- | ---------- | --- | --- | --- | -------------------- | --- | -------------------- | --- | --- | --- | --- | --- | --- |
|          |            |     |     |     | (cid:18) ρu2(cid:19) |     | (cid:18) ρu3(cid:19) |     |     |     |     |     |     |
|          |            |     |     | ∂   |                      | ∂   |                      |     | ∂p  |     |     |     |     |
|          |            |     |     |     |                      | +   |                      | +u  | =0  |     |     |     | (8) |
|          |            |     |     | ∂t  | 2                    | ∂x  | 2                    |     | ∂x  |     |     |     |     |
The kinetic energy equation is thus implicit in the mass and momentum conservation equations, and the
numerical discretizations of mass and momentum may not necessarily enforce kinetic energy preserving.
| Thus, developing |     | kinetic energy |              | preserving | schemes  |     | is a separate |     | endeavor. |     |     |     |     |
| ---------------- | --- | -------------- | ------------ | ---------- | -------- | --- | ------------- | --- | --------- | --- | --- | --- | --- |
| Now consider     |     | the energy     | conservation |            | equation |     | (1d).         |     |           |     |     |     |     |
|                  |     |                |              |            | ∂        | ∂   |               |     |           |     |     |     |     |
|                  |     |                |              |            | (ρE)+    |     | (ρuE+pu)=0    |     |           |     |     |     |     |
|                  |     |                |              |            | ∂t       | ∂x  |               |     |           |     |     |     |     |
u2
Using the relation E = e+ in the energy equation, the internal energy and kinetic energy parts can be
2
separated.
|     |     | ∂     |     | ∂       |     | ∂u  | ∂ (cid:18) ρu2(cid:19) |     | ∂ (cid:18) ρu3(cid:19) |     | ∂p  |     |     |
| --- | --- | ----- | --- | ------- | --- | --- | ---------------------- | --- | ---------------------- | --- | --- | --- | --- |
|     |     | (ρe)+ |     | (ρue)+p |     | +   |                        | +   |                        | +u  | =0  |     |     |
|     |     | ∂t    |     | ∂x      |     | ∂x  | ∂t 2                   |     | ∂x                     | 2   | ∂x  |     |     |
or
|     |     | ∂               | ∂             |                    |          | ∂u   | ∂                   | (cid:18) ρu2(cid:19) | ∂            | (cid:18) ρu3(cid:19) | ∂p     |           |     |
| --- | --- | --------------- | ------------- | ------------------ | -------- | ---- | ------------------- | -------------------- | ------------ | -------------------- | ------ | --------- | --- |
|     |     | (ρe)+           | (ρue)+ρe(γ−1) |                    |          |      | +                   |                      | +            |                      | +u     | =0        | (9) |
|     |     | ∂t              | ∂x            |                    |          | ∂x   | ∂t                  | 2                    | ∂x           | 2                    | ∂x     |           |     |
|     |     | (cid:124)       |               | (cid:123)(cid:122) |          |      | (cid:125) (cid:124) |                      |              | (cid:123)(cid:122)   |        | (cid:125) |     |
|     |     | internal energy | consevation   |                    | + source | term |                     |                      |              |                      |        |           |     |
|     |     |                 |               |                    |          |      | kinetic             | energy               | conservation | +                    | source | term      |     |
Now, subtracting (8) from the above equation, one can focus only on the internal energy equation.
|     |     |     |     |     | ∂     | ∂   |         | ∂u  |     |     |     |     |      |
| --- | --- | --- | --- | --- | ----- | --- | ------- | --- | --- | --- | --- | --- | ---- |
|     |     |     |     |     | (ρe)+ |     | (ρue)+p |     | =0  |     |     |     | (10) |
|     |     |     |     |     | ∂t    | ∂x  |         | ∂x  |     |     |     |     |      |
p
| Utilizing | e=  | , we obtain |     |     |     |     |     |     |     |     |     |     |     |
| --------- | --- | ----------- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
ρ(γ−1)
|     |     |     |     |     | ∂p  | ∂(up) |         | ∂u  |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | ----- | ------- | --- | --- | --- | --- | --- | --- |
|     |     |     |     |     | +   |       | +(γ−1)p |     | =0  |     |     |     |     |
|     |     |     |     |     | ∂t  | ∂x    |         | ∂x  |     |     |     |     |     |
or
|     |     |     |     |     | ∂p  | ∂p  |     | ∂u  |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
|     |     |     |     |     |     | +u  | +γp | =0  |     |     |     |     |     |
|     |     |     |     |     | ∂t  | ∂x  |     | ∂x  |     |     |     |     |     |
dp γdρ,
| Using ds= | −   | we obtain |     | from | the above | equation |     |     |     |     |     |     |     |
| --------- | --- | --------- | --- | ---- | --------- | -------- | --- | --- | --- | --- | --- | --- | --- |
p ρ
|     |     |     |     | ∂s  | γp∂ρ | ∂s  | uγp∂ρ |      | ∂u  |     |     |     |     |
| --- | --- | --- | --- | --- | ---- | --- | ----- | ---- | --- | --- | --- | --- | --- |
|     |     |     |     | p + |      | +up | +     |      | +γp | =0  |     |     |     |
|     |     |     |     | ∂t  | ρ ∂t | ∂x  |       | ρ ∂x | ∂x  |     |     |     |     |
or
|           |           |                |     | (cid:20)  |     |     | (cid:21) | (cid:20) |     | (cid:21) |     |     |      |
| --------- | --------- | -------------- | --- | --------- | --- | --- | -------- | -------- | --- | -------- | --- | --- | ---- |
|           |           |                |     | γp ∂ρ     | ∂ρ  |     | ∂u       | ∂s       | ∂s  |          |     |     |      |
|           |           |                |     |           | +u  | +ρ  | +p       |          | +u  | =0       |     |     | (11) |
|           |           |                |     | ρ ∂t      | ∂x  |     | ∂x       | ∂t       | ∂x  |          |     |     |      |
| or (after | utilizing | the continuity |     | equation) |     |     |          |          |     |          |     |     |      |
|           |           |                |     |           |     | ∂s  | ∂s       |          |     |          |     |     |      |
|           |           |                |     |           |     | +u  |          | =0       |     |          |     |     | (12) |
|           |           |                |     |           |     | ∂t  | ∂x       |          |     |          |     |     |      |
Thus, the entropy convection equation is a part of the internal energy conservation equation. Together with
the mass conservation equation, we can obtain the entropy conservation equation as
|     |     |     |     |     | ∂   |       | ∂        |     |     |     |     |     |      |
| --- | --- | --- | --- | --- | --- | ----- | -------- | --- | --- | --- | --- | --- | ---- |
|     |     |     |     |     |     | (ρs)+ | (ρus)=0. |     |     |     |     |     | (13) |
|     |     |     |     |     | ∂t  |       | ∂x       |     |     |     |     |     |      |
4

Thus, the entropy conservation equation is implicit in the internal energy conservation equation. Again,
like the kinetic energy preservation, numerical discretizations of Euler equations do not necessarily satisfy
the discrete entropy conservation. Numerical schemes satisfying these additional conservation equations are
often termed as structure preserving schemes and are discussed in the next subsections. The fact that the
total energy equation contains both the entropy conservation and the kinetic energy conservation is used in
section 2 to construct entropy conservative flux by only modifying the energy flux. Note that the equations
(8) and (13) are only valid for smooth regions and not valid across discontinuities like shocks.
2.3. Numerical Entropy Conservation
In semi-discrete form, the update formula for a cell (with cell-centre at x and cell-interfaces x ) can
j j±1
2
be expressed as
dU 1 (cid:16) (cid:17)
j =− F −F (14)
dt ∆x j j+1 2 j− 2 1
When multiplied by the entropy variables at cell j, denoted as V , it yields the semi-discrete entropy con-
j
servation equation:
dη(U) 1 (cid:16) (cid:17)
j =− V · F −F (15)
dt ∆x j j j+1 2 j− 2 1
The right-hand side of the above equation is not in conservation form. Therefore, we modify it as
(cid:18) (cid:19) (cid:18) (cid:19)
dη(U) 1 1
∆x j =− V − ∆V ·F + V + ∆V ·F (16)
j dt j+1 2 2 j+1 2 j+ 2 1 j−1 2 2 j− 2 1 j−1 2
by introducing V = (V +V )/2 and ∆V = V −V . The numerical entropy flux consistent
j+1 j+1 j j+1 j+1 j
2 2
with a given potential function ψ at an interface is defined as ζ = V ·F −ψ , where ψ =
j+1 j+1 j+1 j+1 j+1
2 2 2 2 2
(ψ +ψ )/2. Thus, the above equation simplifies to
j+1 j
dη(U) 1 1 1 1
∆x j +ζ −ζ = ∆V ·F − ∆ψ + ∆V ·F − ∆ψ (17)
j dt j+1 2 j−1 2 2 j+1 2 j+1 2 2 j+ 2 1 2 j−1 2 j−1 2 2 j− 2 1
To ensure the numerical conservation of entropy, the right hand side must vanish. Therefore, the entropy
conservative numerical fluxes Fc should satisfy the condition [1; 2]:
(V −V )·Fc =ψ −ψ (18)
j+1 j j+1 j+1 j
2
2.4. Numerical Kinetic Energy Preservation
Asemi-discreteformofthekineticenergyequationcanbeobtainedbymultiplyingthesemi-discretemass
and momentum equations by −u2/2 and u , respectively, and adding them.
j j
(cid:32) (cid:33)
d u2 u2 (cid:16) (cid:17) (cid:16) (cid:17) (cid:16) (cid:17)
∆x ρ j = j (ρu) −(ρu) −u (ρu2) −(ρu2) −u p −p (19)
jdt j 2 2 j+1 2 j−1 2 j j+ 2 1 j− 2 1 j j+1 2 j− 2 1
Jameson [9] introduced the following form of convective momentum.
(ρu2) =(ρu) (u +u )/2 (20)
j+1 j+1 j+1 j
2 2
We get the local kinetic energy preservation equation using this flux in (19).
∆x d (cid:32) ρ j u2 j (cid:33) + (cid:18) ρu3(cid:19) − (cid:18) ρu3(cid:19) =−u (p −p ) (21)
jdt 2 2 j+1 2 j−1 j j+1 2 j−1 2
2 2
5

|     |     |     |     |     |     | (cid:16) ρu3 (cid:17) |     | (cid:0)uj+1uj | (cid:1) |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --------------------- | --- | ------------- | ------- | --- | --- | --- | --- |
Here the kinetic energy flux is defined as = (ρu) which is of the form F(U ,U ) =
|     |     |     |     |     |     | 2   |     | 2   | j+1 |     |     |     | L R |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
|     |     |     |     |     |     |     | j+1 |     |     | 2   |     |     |     |
2
−F(U R ,U L ) (given mass flux satisfies the same property). Summing equation (19) on all j = 1 to n,
using the form of convective momentum flux as in the above equation, and taking the boundary mass
fluxes as (ρu) = ρ u and (ρu) = ρ u and momentum fluxes as (ρu2) = ρ u2 + p and
|     |     | 1/2 | L L |     | n+1/2 | R   | R   |     |     |     | 1/2 | L   | L L |
| --- | --- | --- | --- | --- | ----- | --- | --- | --- | --- | --- | --- | --- | --- |
(ρu2) =ρ u2 +p , we recover the global kinetic energy preservation equation.
| n+1/2 |          | R R | R           |          |          |      |          |          |        |            |          |     |        |
| ----- | -------- | --- | ----------- | -------- | -------- | ---- | -------- | -------- | ------ | ---------- | -------- | --- | ------ |
|       | n        |     | (cid:32) u2 | (cid:33) | (cid:18) |      | (cid:19) | (cid:18) |        | (cid:19) n |          |     |        |
|       | (cid:88) |     | d ρ j       |          | ρ        | u2   |          | ρ        | u2     | (cid:88)   |          |     |        |
|       |          | ∆x  |             | j +u     | R        | R +p | −u       |          | L L +p | =          | p (u     | −u  | ) (22) |
|       |          | jdt |             |          | R        |      | R        | L        | L      |            | j+ 1 j+1 | j   |        |
|       |          |     | 2           |          | 2        |      |          |          | 2      |            | 2        |     |        |
|       | j=1      |     |             |          |          |      |          |          |        | j=1        |          |     |        |
Note that we are still to choose the forms of mass flux (ρu) j+1 , pressure flux p j+1 and the energy flux
|     |     |     |     |     |     |     |     |     | 2   |     | 2   |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
(ρuE +pu) . Gassner et al. [22] observed that the arithmetic pressure average preserves kinetic energy
j+1
2
better than the schemes using different pressure averages. Ranocha [10] showed that an arithmetic pressure
average is well suited as pressure flux for Euler equations at the incompressible limit. This flux formulation
stillleavesmassflux(ρu) ,andenergyflux(ρuE+pu) undefinedandallowsustofixtheseusingother
|              |     |              | j+1       |             |      |        |            | j+1 |      |     |     |     |     |
| ------------ | --- | ------------ | --------- | ----------- | ---- | ------ | ---------- | --- | ---- | --- | --- | --- | --- |
| criteria.    |     |              | 2         |             |      |        |            | 2   |      |     |     |     |     |
| 3. Entropy   |     | Conservative |           | and Kinetic |      | Energy | Preserving |     | Flux |     |     |     |     |
| 3.1. Entropy |     | Conservative | Numerical |             | Flux |        |            |     |      |     |     |     |     |
Entropy conservative fluxes can be constructed by assuming Fc to be of the following form (average
j+1
2
| flux + | diffusive | flux) | as  |     |     |     |     |     |     |     |     |     |     |
| ------ | --------- | ----- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
|        |           | 1     |     | 1   |     |     |     |     | 1   |     |     | 1   |     |
Fc
= (F j+1 +F j )− α(V j+1 −V j )=F j+1/2 − α(V j+1 −V j )=F j+1/2 − α∆V (23)
|     | j+  | 1 2 |     | 2   |     |     |     |     | 2   |     |     | 2   |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
2
Here, α is a scalar diffusion coefficient, adding diffusion to satisfy semi-discrete entropy conservation. It can
be obtained by substituting the flux (23) into (18). Thus, ∆V·Fc =∆ψ gives
j+1
2
|     |     |     |     |     |     | (cid:18) |     |     | (cid:19) |     |     |     |     |
| --- | --- | --- | --- | --- | --- | -------- | --- | --- | -------- | --- | --- | --- | --- |
1
|     |     |     |     |     | ∆V· | F   | −   | α∆V | =∆ψ |     |     |     | (24) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | ---- |
j+1/2
2
| This gives | one | form of | diffusion | coefficient, |     | termed | as α  | 1 here, | as    |        |     |     |      |
| ---------- | --- | ------- | --------- | ------------ | --- | ------ | ----- | ------- | ----- | ------ | --- | --- | ---- |
|            |     |         | α         | F·∆V−∆ψ      |     |        | F1∆V  | +F2∆V   | +F3∆V | −∆ψ    |     |     |      |
|            |     |         | 1         |              |     |        |       | 1       | 2     | 3      |     |     |      |
|            |     |         |           | =            |     | =      |       |         |       |        |     |     | (25) |
|            |     |         | 2         | ∆V·∆V        |     |        | ∆V ∆V | +∆V     | ∆V    | +∆V ∆V |     |     |      |
|            |     |         |           |              |     |        | 1     | 1       | 2 2   | 3      | 3   |     |      |
where F1, F2 and F3 are the average mass, momentum and energy fluxes of jth and (j+1)th cells. ∆(X)
of any quantity X is defined as ∆(X)=X j+1 −X j . This flux, termed EC1, can be defined by
1
|     |     |     |     |     |     | F1  | =F1− | α (∆V | )   |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | ---- | ----- | --- | --- | --- | --- | --- |
|     |     |     |     |     |     | EC1 |      | 2 1   | 1   |     |     |     |     |
1
F2
|     |     |     |     |     |     |     | =F2− | α (∆V | )   |     |     |     | (26) |
| --- | --- | --- | --- | --- | --- | --- | ---- | ----- | --- | --- | --- | --- | ---- |
|     |     |     |     |     |     | EC1 |      | 2 1   | 2   |     |     |     |      |
1
F3
|     |     |     |     |     |     |     | =F3− | α 1 (∆V | 3 ) |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | ---- | ------- | --- | --- | --- | --- | --- |
|     |     |     |     |     |     | EC1 |      | 2       |     |     |     |     |     |
It can be noted that α 1 is a scalar diffusion that acts on mass, momentum, and energy equations, resulting
in an entropy conservative flux. Note that this flux can also be obtained through an optimization procedure
(see appendix A) as given by Abgrall [23]. This flux, however, is not kinetic energy preserving. It was shown
in section 2.2 that the entropy equation is a part of the energy equation, and just adding a corresponding
6

]T
diffusion to the energy flux can ensure entropy conservation. This can be achieved by taking α =[0,0,α
2 2
| in equation | (23) | and | we get | the following |     | form of | α   |     |     |     |
| ----------- | ---- | --- | ------ | ------------- | --- | ------- | --- | --- | --- | --- |
2
|            |        |      | α   | 2 F·∆V−∆ψ  |       |     | F1∆V | 1 +F2∆V | 2 +F3∆V 3 −∆ψ |      |
| ---------- | ------ | ---- | --- | ---------- | ----- | --- | ---- | ------- | ------------- | ---- |
|            |        |      |     | =          |       | =   |      |         |               | (27) |
|            |        |      |     | 2          | ∆V ∆V |     |      |         | ∆V ∆V         |      |
|            |        |      |     |            | 3     | 3   |      |         | 3 3           |      |
| This flux, | termed | EC2, | can | be defined | by    |     |      |         |               |      |
F1
=F1
EC2
|     |     |     |     |     |     | F2  | =F2 |     |     |      |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | ---- |
|     |     |     |     |     |     | EC2 |     |     |     | (28) |
1
|     |     |     |     |     |     | F3  | =F3− | α (∆V | )   |     |
| --- | --- | --- | --- | --- | --- | --- | ---- | ----- | --- | --- |
|     |     |     |     |     |     | EC2 |      | 2     | 3   |     |
2
The entropy conservation of a numerical flux thus can theoretically be achieved by adding diffusion only to
theenergyequation. Thismakesmassandmomentumfluxavailableforothermodifications, asshowninthe
| following    | subsection. |              |     |             |        |            |     |           |      |     |
| ------------ | ----------- | ------------ | --- | ----------- | ------ | ---------- | --- | --------- | ---- | --- |
| 3.2. Entropy |             | Conservative |     | and Kinetic | Energy | Preserving |     | Numerical | Flux |     |
Jameson[9]givesthefollowingnumericalconditiononmassandmomentumfluxesforthediscretekinetic
| energy preservation |     |     | condition | to be | satisfied. |     |        |     |     |      |
| ------------------- | --- | --- | --------- | ----- | ---------- | --- | ------ | --- | --- | ---- |
|                     |     |     |           |       |            | F2  |        | F1  |     |      |
|                     |     |     |           |       |            | j+1 | =u j+1 | j+1 | +p  | (29) |
|                     |     |     |           |       |            | 2   |        | 2 2 |     |      |
1(u
Here u j+1 = j+1 +u j ) and p is any consistent pressure flux. An entropy conservative and kinetic energy
2 2
preservingfluxcanbeobtainedbytakingthecorrect formofmomentumfluxgivenby (29)inequation(23).
Like EC2 flux, the diffusion term α 3 is only applied to the energy equation. Thus, the entropy conservative
and kinetic energy preserving numerical flux (termed ECKEP here) is given as
|     |     |     |     |     |     | F1  | =F1 |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
ECKEP
|     |     |     |     |     |     | F2    | =F1u+p |     |     |      |
| --- | --- | --- | --- | --- | --- | ----- | ------ | --- | --- | ---- |
|     |     |     |     |     |     | ECKEP |        |     |     | (30) |
1
|         |     |     |     |     |     | F3    | =F3− | α   | (∆V ) |     |
| ------- | --- | --- | --- | --- | --- | ----- | ---- | --- | ----- | --- |
|         |     |     |     |     |     | ECKEP |      | 2   | 3 3   |     |
| where α | is  |     |     |     |     |       |      |     |       |     |
3
|     |     |     |     | α   | F1∆V | +(F1u+p)∆V |     |     | +F3∆V −∆ψ |      |
| --- | --- | --- | --- | --- | ---- | ---------- | --- | --- | --------- | ---- |
|     |     |     |     | 3   |      | 1          |     | 2   | 3         |      |
|     |     |     |     |     | =    |            |     |     |           | (31) |
|     |     |     |     | 2   |      |            | ∆V  | ∆V  |           |      |
3 3
The convective momentum flux in (30) is taken such that it results in the preservation of kinetic energy, and
| the pressure | flux | is taken | as  | an arithmetic |     | average. |     |     |     |     |
| ------------ | ---- | -------- | --- | ------------- | --- | -------- | --- | --- | --- | --- |
Note that to prevent denominators of diffusion terms from becoming zero in fluxes (25), (27) and (31), a
small parameter δ is added to the denominators (δ =10−16 in our computations).
4. Numerical investigation of the entropy conservative and kinetic energy preserving fluxes
| 4.1. Exact | Capturing |     | of Steady | Contact |     | Discontinuity |     |     |     |     |
| ---------- | --------- | --- | --------- | ------- | --- | ------------- | --- | --- | --- | --- |
Across a stationary contact discontinuity in one dimension, we have primitive variable vectors W on the
| left and | right | sides given | by  |     |        |       |     |     |            |     |
| -------- | ----- | ----------- | --- | --- | ------ | ----- | --- | --- | ---------- | --- |
|          |       |             |     |     |        |  ρ  |     |     |  ρ       |     |
|          |       |             |     |     |        | L     |     |     | R          |     |
|          |       |             |     |     | W L = | 0    | and |     | W R = 0  |     |
|          |       |             |     |     |        | p     |     |     | p          |     |
7

Resolving contact disctontinuities accurately is important for good resolution of shear and boundary layers
inmultidimensionalflow [24]. Wecaneasilycheckwhetheranumericalschemewillcaptureasteadycontact
discontinuity exactly. Flux at the interface across a steady contact disctontinuity is given by
|     |     |     |     |    |     |   |    |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
|     |     |     |     |     | ρu  | 0   |     |     |     |
ρu2+p
|     |     |     | F j+1/2 | =  |        |  =p |     |     |     |
| --- | --- | --- | ------- | --- | ------ | ------ | --- | --- | --- |
|     |     |     |         |     | ρuE+pu | 0      |     |     |     |
j+1
2
=[1.4,0,1]T
TheaboveconditionissatisfiedbyEC1,EC2,andECKEPfluxes. AReimannproblemwithW L
and W = [1,0,1]T is considered and the solution is shown in figure 1. The steady contact discontinuities
R
| are preserved | exactly for | all the | three numerical |     | fluxes | presented. |     |     |     |
| ------------- | ----------- | ------- | --------------- | --- | ------ | ---------- | --- | --- | --- |
1.4
|     |     |     |     | Exact |     | 1   |     | Exact |     |
| --- | --- | --- | --- | ----- | --- | --- | --- | ----- | --- |
EC1
|     |     |     |     | EC2   |     |     |     | EC1   |     |
| --- | --- | --- | --- | ----- | --- | --- | --- | ----- | --- |
|     | 1.3 |     |     |       |     |     |     | EC2   |     |
|     |     |     |     | ECKEP |     |     |     | ECKEP |     |
0.5
1.2
|     |     |     |     |     | U   | 0   |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
1.1
-0.5
1
|     | 0   | 0.2 0.4 | 0.6 | 0.8 | 1   | 0 0.2 0.4 | 0.6 | 0.8 1 |     |
| --- | --- | ------- | --- | --- | --- | --------- | --- | ----- | --- |
|     |     |         | x   |     |     |           | x   |       |     |
Figure1: DensityandvelocityplotsforastationarycontactdiscontinuityatT=2sbythethreeentropyconservativefluxes
| 4.2. Experimental | order | of convergence |     | (EOC) |     |     |     |     |     |
| ----------------- | ----- | -------------- | --- | ----- | --- | --- | --- | --- | --- |
Entropyconservativefluxes(26)and(28)andtheentropyconservativeandkineticenergypreservingflux
(30) are all second-order accurate in space. This can be demonstrated by taking a sinusoidal test case with
| the following | initial conditions. |     |     |                      |     |     |     |     |       |
| ------------- | ------------------- | --- | --- | -------------------- | --- | --- | --- | --- | ----- |
|               |                     |     |     | ρ(x,0)=1+0.2sin(2πx) |     |     |     |     | (32a) |
|               |                     |     |     | u(x,0)=0.1           |     |     |     |     | (32b) |
|               |                     |     |     | p(x,0)=1             |     |     |     |     | (32c) |
EOCwascomputedonadomainof[0,1]. Aperiodicboundarywasappliedatbothends. Theexactsolution
to the above problem is known. The velocity and pressure remain constant and the density wave advances
| with time | given by the following |     | equation.                   |     |     |     |     |     |      |
| --------- | ---------------------- | --- | --------------------------- | --- | --- | --- | --- | --- | ---- |
|           |                        |     | ρ(x,t)=1+0.2sin(2π(x−0.1t)) |     |     |     |     |     | (33) |
The computational domain was divided into N=40,80,160,320,640 and 1280 points, and L and L errors in
|     |     |     |     |     |     |     |     | 1 2 |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
density were computed at t=10. The SSPRK-3 [25] method was used for time integration. For a kth order
accurate scheme, the error on N points (dx=1/N) and 2N points (dx=1/2N) are given by
∥=C(∆x)k+O(∆x)k+1
|           |                      |     | ∥ϵ     | N          |                       |                         |     |     | (34a) |
| --------- | -------------------- | --- | ------ | ---------- | --------------------- | ----------------------- | --- | --- | ----- |
|           |                      |     |        |            | (cid:18) ∆x (cid:19)k | (cid:18) ∆x (cid:19)k+1 |     |     |       |
|           |                      |     | ∥ϵ 2N  | ∥=C        |                       | +O                      |     |     | (34b) |
|           |                      |     |        |            | 2                     | 2                       |     |     |       |
| Thus, the | order of convergence | k   | can be | calculated | using                 |                         |     |     |       |
|           |                      |     |        |            |                       | (cid:18) (cid:19)       |     |     |       |
∥ϵ ∥
|     |     |     |     | k =log |     | N   |     |     | (35) |
| --- | --- | --- | --- | ------ | --- | --- | --- | --- | ---- |
2
∥ϵ ∥
2N
8

L and L errors with computed order of accuracy for schemes EC1, EC2 and ECKEP are given in tables 1,
1 2
2 and 3 respectively. As shown, all fluxes were found to be of the order O(∆x)2.
| N L 1 Error        | EOC        | L 2 Error        | EOC        |
| ------------------ | ---------- | ---------------- | ---------- |
| 40 0.00328453      | –          | 0.00364681       | –          |
| 80 0.00082177      | 1.99888941 | 0.00091259       | 1.99859961 |
| 160 0.00020546     | 1.99988108 | 0.00022820       | 1.99965607 |
| 320 0.00005137     | 1.99998144 | 0.00005705       | 1.99991456 |
| 640 0.00001284     | 1.99999608 | 0.00001426       | 1.99997867 |
| 1280 0.00000321    | 1.99999922 | 0.00000357       | 1.99999468 |
| Table1: EOCusingL1 | andL2      | errorsforEC1flux |            |
| N L Error          | EOC        | L Error          | EOC        |
1 2
| 40 0.00325342      | –          | 0.00361227       | –          |
| ------------------ | ---------- | ---------------- | ---------- |
| 80 0.00081392      | 1.99899572 | 0.00090405       | 1.99843686 |
| 160 0.00020351     | 1.99979790 | 0.00022607       | 1.99962551 |
| 320 0.00005088     | 1.99995364 | 0.00005652       | 1.99990901 |
| 640 0.00001272     | 1.99998905 | 0.00001413       | 1.99997750 |
| 1280 0.00000318    | 1.99999723 | 0.00000353       | 1.99999439 |
| Table2: EOCusingL1 | andL2      | errorsforEC2flux |            |
| N L Error          | EOC        | L Error          | EOC        |
1 2
| 40 0.00325342      | –          | 0.00361227         | –          |
| ------------------ | ---------- | ------------------ | ---------- |
| 80 0.00081392      | 1.99899571 | 0.00090405         | 1.99843686 |
| 160 0.00020351     | 1.99979789 | 0.00022607         | 1.99962551 |
| 320 0.00005088     | 1.99995362 | 0.00005652         | 1.99990901 |
| 640 0.00001272     | 1.99998898 | 0.00001413         | 1.99997750 |
| 1280 0.00000318    | 1.99999768 | 0.00000353         | 1.99999441 |
| Table3: EOCusingL1 | andL2      | errorsforECKEPflux |            |
9

10-1
10-2
10-3
10-4
10-5
10-6
10-3 10-2
Grid size
rorre
1L
10-1
EC1 Flux
Slope=2
10-2
10-3
10-4
10-5
10-6
10-3 10-2
Grid size
rorre
1L
10-1
EC2 Flux
Slope=2
10-2
10-3
10-4
10-5
10-6
10-3 10-2
Grid size
rorre
1L
ECKEP Flux
Slope=2
10-1
10-2
10-3
10-4
10-5
10-6
10-3 10-2
Grid size
rorre
2L
10-1
EC1 Flux
Slope=2
10-2
10-3
10-4
10-5
10-6
10-3 10-2
Grid size
rorre
2L
10-1
EC2 Flux
Slope=2
10-2
10-3
10-4
10-5
10-6
10-3 10-2
Grid size
rorre
2L
ECKEP Flux
Slope=2
Figure2: Log-logplotoferrorswithgridsizeforEC1,EC2andECKEPfluxesrespectively
4.3. Conservation properties of the new flux functions
4.3.1. Taylor-Green Vortex
Taylor-Green vortex problem is a shock-free test case originally introduced for low-speed compressible
turbulent flows. It has also been used to test the inviscid portion of Navier-Stokes equations and shows
conservationpropertiesofnumericalschemes[26]. A[2π×2π]computationaldomainistakenwith100×100
gridsizewithperiodicboundaryconditions. Theproblemwasrunforanextendedtimeof20seconds. Initial
conditions are given by
ρ=1,
u=sin(x)×cos(y),
v =−cos(x)×sin(y), (36)
100 cos(2x)+cos(2y)
p= +
γ 4
Conservation quantities are only spatially conserved, so it is desirable to have a high order time-stepping
with a low CFL number. Thus, third-order SSPRK [25] was implemented for temporal integration with a
CFL of 0.1. The objective was to compare numerically different conservative fluxes. Comparison is made of
EC1, EC2, ECKEP, and the entropy conservative and kinetic energy preserving flux of Chandrasekar [11]
(cid:80)
(referred to as PC ECKEP). It can be seen that all fluxes preserve the total entropy in domain, η , with
j
∀j
sufficient accuracy for the given time, as shown in figure 3. Total kinetic energy, (cid:80) ρ u2/2, keeps reducing
j j
∀j
forEC1fromtheverystartofthesolutionanddivergesforEC2atalatertime(approximately25sec)which
is expected as these schemes are not designed to be necessarily kinetic energy preserving. Schemes ECKEP
and PC ECKEP are kinetic energy preserving and no change is observed as evident from figure 4. Only the
semi-discrete cases are considered for preserving the structures in this study. Extending these ideas to the
fully discrete case is more involved, is some times known to generate oscillations (for example, see [27]) and
is not attempted in this work.
10

Figure3: TotalEntropyforTaylor-Greensvortex Figure4: TotalKineticEnergyforTaylor-Greensvortex
4.3.2. Isentropic vortex convection
Isentropic vortex is another shock-free test case with the following initial data, W [28].
o
γ−1
ρ(x,y)=(1− ω2)γ− 1 1,
2
(y−y )ω
u(x,y)=Mcos(θ)− c ,
R
(37)
(x−x )ω
v(x,y)=Msin(θ)+ c ,
R
1
p(x,y)= ργ
γ
Here ω is a Gaussian function of the form ω =βexp(f). The strength of the vortex is determined by β, and
the perturbation function f is given by
1 (cid:18)(cid:16)x(cid:17)2 (cid:16)y(cid:17)2 (cid:19)
f(x,y)=− + .
2σ2 R R
(x ,y ) is initially the centre of the vortex. For this study, a square domain of 200x200 cells was taken
c c
with periodic boundaries on all sides. β = 1, R = 1 and σ = 1 were taken in (37). Mach number is 2/γ
and the flow angle (θ) is 45 degrees. The exact solution of the above problem is the above initial data
translated with Mcos(θ) and Msin(θ), i.e., W(x,y) = W (x−Mcos(θ)t,y −Msin(θ)t). As shown in
o
Figure5: Totalentropyandtotalkineticenergyv/sTimeforIsentropicvortextestcasefordifferententropyconservativefluxes
figure 5, entropy remains constant throughout the simulations for PC ECKEP, EC1, and ECKEP schemes
11

and for EC2 scheme up to t=25 s. Divergence of EC2 scheme is unexpected but can be attributed to the
accumulation of errors in temporal discretization (as the entropy conservation formulation is forced only for
the semi-discrete case). Kinetic energy preserving schemes, PC ECKEP and ECKEP perform as expected
and kinetic energy remains constant as shown in figure 5. Kinetic energy diverges with time for EC1 and
EC2 schemes since these schemes do not necessarily preserve kinetic energy numerically.
| 4.4. Computational |     | Cost |     |     |     |     |     |     |     |     |
| ------------------ | --- | ---- | --- | --- | --- | --- | --- | --- | --- | --- |
The entropy conservative fluxes presented do not need to evaluate computationally expensive averages
(unlike in Riemann solvers) and are expected to be more efficient. To test this, we ran the smooth test case
problem given in (32) for 1000 computational cells and with SSPRK-2 discretization of temporal derivatives
for a final time of 10 seconds. This resulted in 47,477,430 flux function calls for all the schemes, and
the average time per call is given in table 4. Out of all entropy conservative schemes, the EC2 scheme is
consecutively29%,11%,and5%moreefficientthanROEEC(referringtoRoe’sentropyconservativescheme
[15]), ECKEP, and EC1 schemes, respectively. Also, the new ECKEP scheme is computationally 8% more
| efficient than | the PC | ECKEP | scheme. |     |         |          |         |                |     |     |
| -------------- | ------ | ----- | ------- | --- | ------- | -------- | ------- | -------------- | --- | --- |
|                |        | Flux  |         |     | Average | time per | call    | (microseconds) |     |     |
|                |        | ROE   | EC      |     |         |          | 0.12947 |                |     |     |
|                |        | PC    | ECKEP   |     |         |          | 0.10401 |                |     |     |
|                |        | EC1   |         |     |         |          | 0.09648 |                |     |     |
|                |        | EC2   |         |     |         |          | 0.09215 |                |     |     |
|                |        | ECKEP |         |     |         |          | 0.09532 |                |     |     |
Table4: Wallclocktimesforvariousentropyconservativeschemes
| 5. Rankine–Hugoniot |     | condition |     | satisfying |     | Entropy | Stable | Scheme |     |     |
| ------------------- | --- | --------- | --- | ---------- | --- | ------- | ------ | ------ | --- | --- |
Entropy conservative flux is only valid at smooth regions, and additional diffusion must be added at
discontinuities such as shocks. Additionally, the diffusion must generate entropy such that inequality in (5)
is satisfied. Ismail and Roe [15] gave one possible entropy stable diffusion, which linearizes Euler equations
aboutanaveragestate. Morerecently,Chandrashekar[11]introducedadiffusionwhichsatisfiesbothentropy
stability and kinetic energy preservation. However, the above types of diffusion depend on the underlying
eigen-structure of Euler equations. With the motivation to avoid complicated Riemann solvers and their
strong dependence on the eigen-structure, the simple central scheme (which still catpures steady discontinu-
itires exactly) introduced by Jaisankar and Raghurama Rao [29] (Method of optimal viscosity for enhanced
resolution of shocks, MOVERS) is utilized as a basic framework for the additional diffusion proposed here.
The above scheme satisfies the Rankine-Hugonoit jump conditions across an interface in determining numer-
ical diffusion. Further advantage of this approach is that this strategy can be extended to any system of
hyperbolic conservation laws with an R-H-like condition. The general form of diffusive flux for this central
| scheme is | given by |     |     |     |     |     |     |     |     |     |
| --------- | -------- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
1 D˜
|     |     |     |     | Fd  | =−  | (U    |     | −U ) |     | (38) |
| --- | --- | --- | --- | --- | --- | ----- | --- | ---- | --- | ---- |
|     |     |     |     |     | j+1 | 2 j+1 | j+1 | j    |     |      |
|     |     |     |     |     | 2   | 2     |     |      |     |      |
For MOVERS-n the diffusion is given by D˜ = diag(s1,s2,s3) where wave-speeds sk are computed
|             |      |     |          |     | j+1 |     |     | j+1 |     |     |
| ----------- | ---- | --- | -------- | --- | --- | --- | --- | --- | --- | --- |
|             |      |     |          |     | 2   |     |     |     | 2   |     |
| using ∆F=D˜ | ·∆U. | D˜  | is given | by  |     |     |     |     |     |     |
j+1
2
|     |     |     |     |     | (cid:12) (cid:12)F 1 −F | 1(cid:12)  |     |     |    |     |
| --- | --- | --- | --- | --- | ------------------------ | ---------- | --- | --- | --- | --- |
|     |     |     |     |     | j +1                     | j (cid:12) | 0   |     | 0   |     |
(cid:12)U1 −U1(cid:12)
|     |     |     |     |     | j+1 | j                    |              |             |             |     |
| --- | --- | --- | --- | --- | --- | -------------------- | ------------ | ----------- | ----------- | --- |
|     |     |     |     |     |    | (cid:12) (cid:12)F 2 | −F 2(cid:12) |             |            |     |
|     |     |     | D˜  | =   |  0 | j +1                 | j (cid:12)   |             | 0          |     |
|     |     |     |     | j+1 |    | (cid:12)U 2          | −U 2(cid:12) |             |            |     |
|     |     |     |     | 2   |    | j +1                 | j            | (cid:12)    | 3(cid:12)  |     |
|     |     |     |     |     |     |                      |              | (cid:12)F 3 | −F          |     |
|     |     |     |     |     | 0   |                      | 0            | j +1        | j (cid:12)  |     |
|     |     |     |     |     |     |                      |              | (cid:12)U3  | −U3(cid:12) |     |
|     |     |     |     |     |     |                      |              | j+1         | j           |     |
12

To prevent unphysical values of wave speeds when ∆Uk is very small, a wave-speed correction is introduced
such that
sk =sign(sk )λ if |sk |≥λ
j+1 j+1 max j+1 max
2 2 2
sk =sign(sk )λ if |sk |≤λ
j+1 j+1 min j+1 min
2 2 2
where λ = min(|u−a|,|u|,|u+a|) and λ = max(|u−a|,|u|,|u+a|). When Uk → Uk wave speed
min max j+1 j
is taken as λ . We propose an accurate discontinuity capturing version of MOVERS flux (RH) by taking
min
the following form of diffusion flux.
1 1
FRH =− min(s ,s ,s ) I˜(U −U )=− (αS) I˜(U −U ) (39)
j+1 2 2 1 2 3 j+ 2 1 j+1 j 2 j+1/2 j+1 j
A fix is applied to the diffusion to enable a smooth transition at the sonic point and prevent a sonic glitch.
(αS)2+Θ2
αS = where Θ=0.1 (40)
2Θ
The above fix is only used if αS is not zero. Across any isolated contact discontinuity moving with a
speed of v in one dimension, all the wave speeds computed using ∆F/∆U collapse to the speed of the
d
contact discontinuity, i.e., s = s = s = v . This allows us to preserve exactly steady (and grid-aligned
1 2 3 d
in multi-dimensions) contact discontinuities since d v = 0 for them and thus no diffusion is added. For
d
moving contact discontinuities, the diffuion is proportional to the speed of discontinuity. This gives a better
wavespeedestimatethantheLocal-LaxFriedrichs(Rusnaov)schemeandadditionallypreservessteady(and
grid-aligned) discontinuities. An entropy stable flux at interface x is now given by
j+1
2
F =FEC +FRH (41)
j+1 j+1 j+1
2 2 2
Theabovecombinationwillpreservesteadycontactdiscontinuitiesexactlysinceboththeentropyconservative
fluxes (as discussed in section 4.1) and the entropy stable dissipation (based on R-H conditions) preserves
steady contact discontinuities exactly. Even though the entropy stable dissipation also preserves steady
shocks, this property is not satisfied by the entropy conservative fluxes. We can show that this leads to an
entropy stable scheme by taking the dot product of the semi-discrete conservation law by entropy variable
V .
j
dη(U)
∆x j =−V ·(F −F )
dt j j+1 2 j−1 2
(cid:18) (cid:19) (cid:18) (cid:19)
dη(U) 1 1 1 1
∆x j =− (V +V )− (V −V ) ·F + (V +V )+ (V −V ) F
dt 2 j+1 j 2 j+1 j j+1 2 2 j j−1 2 j j−1 j−1 2
(cid:20)
dη(U) 1 1 1
∆x j =− (V +V )·FEC − (V −V )·FEC − (V +V )·(αS) (U −U )
dt 2 j+1 j j+1 2 2 j+1 j j+1 2 2 j+1 j j+ 2 1 j+1 j
1 1 1
− (V +V )·FEC + (V −V )·FEC + (V +V )(αS) (U −U )
2 j j−1 j−1 2 2 j j−1 j−1 2 2 j j−1 j−1 2 j j−1
(cid:21)
1 1
+ (V −V )(αS) (U −U )+ (V −V )(αS) (U −U )
4 j+1 j j+ 2 1 j+1 j 4 j j−1 j− 2 1 j j−1
Using the entropy conservation relation ∆V ·FEC =∆ψ we get
(cid:20)
dη(U) 1
∆x j +ζ −ζ =− + (V −V )(αS) (U −U )
dt j+1 2 j−1 2 4 j+1 j j+ 2 1 j+1 j
(42)
(cid:21)
1
+ (V −V )(αS) (U −U )
4 j j−1 j− 2 1 j j−1
13

where the numerical entropy flux ζ is given by
j+1
2
1 1 1
ζ = (V +V )·FEC − (ψ +ψ )− (V +V )·(αS) (U −U ) (43)
j+1 2 2 j+1 j j+ 2 1 2 j+1 j 2 j+1 j j+ 2 1 j+1 j
In equation 42, the RHS is always negative given that αS is always positive and since ∆V · ∆U >= 0
(See appendix B). Thus, the scheme always leads to reduction in mathematical entropy and is thus entropy
stable. Scheme (given by (41)) can also be shown to be kinetic energy stable, i.e., it prevents kinetic energy
from growing spuriously. Taking ECKEP flux as the central flux in (41), we get the discrete kinetic energy
preservation equation as
∆x d (cid:32) ρ j u2 j (cid:33) + (cid:18) ρu3(cid:19) − (cid:18) ρu3(cid:19) +u (p −p )=
jdt 2 2 j+1 2 j−1 j j+1 2 j−1 2
2 2
1(cid:104) (cid:105)
− (αS) ρ (u −u )2+(αS) ρ (u −u )2
4 j+ 2 1 j+1 j+1 j j−1 2 j−1 j j−1
(44)
Here, the numerical kinetic energy flux can be written as follows.
(cid:18) ρu3(cid:19)
u u 1
=(ρu) j+1 j − (αS) (ρ u2 −ρ u2) (45)
2 j+1 j+1 2 2 4 j+ 2 1 j+1 j+1 j j
2
The above flux is consistent and conservative. The terms on the RHS of (44) are always negative (given
αS ≥0) and thus lead to kinetic energy stability. In smooth regions where ∆U→0, the diffusion coefficient
for MOVERS type flux tends to λ because of wave speed correction. This results in unnecessary high
max
diffusioninsmoothregionslikeexpansionfans. Tomitigatethisproblem, inthenextsection, weuseashock
sensor to apply RH condition-based diffusion only at shocks. This allows low diffusion in smooth regions
while maintaining accurate shock capturing.
6. Second order hybrid entropy stable scheme
This section proposes a hybrid scheme based on a higher-order diffusion term with entropy conservative
flux, which keeps diffusion low in smooth regions. At shocks, diffusion is switched to RH condition based on
entropy distance, as given by (39). Entropy distance is defined as the dot product between the change in
the entropy variable vector and the change in the conserved variable vector between two points. This senses
gradients of flow variables and thus can be used to identify regions of high gradients and shocks. It has been
used as a means for mesh movement and applying entropy fixes [30] [31] [32]. The entropy distance for a
given interface can be given as
ED =(U −Uj)·(V −Vj)
j+1 j+1 j+1
2
Here, we propose an exponentially scaled entropy distance-based sensor to identify shocks, which is given as
(cid:12) (cid:12)
ϕ =1−(cid:12)exp(−qSED )(cid:12) (46)
j+1 (cid:12) j+1 (cid:12)
2 2
Here SED is the scaled entropy distance given as ED /max(ED ) and q is the scaling factor. Since
j+1 j+1
2 ∀j 2
exponent is computationally expensive to compute, we can use quadratic approximation given as follows
(cid:12) (cid:12)
ϕ =1−(cid:12)1−qSED +(qSED )2(cid:12) (47)
j+1 (cid:12) j+1 j+1 (cid:12)
2 2 2
ϕ is then clipped and flattened as follows
(cid:40)
1 if ϕ >ϵ
j+1
ϕ = 2 (48)
j+1
2 0 otherwise
14

Formostproblems,qistakenbetween8to16andϵistakenbetween10−1 to10−2.
Theabovefunctionmight
notdetectcompleteshockprofilescorrectlyforsolutionsoncoarsegridswhereshockprofilesgetsmearedover
multiple cells. Thus a ”flattening” of was done i.e, ϕ j+1 = max(ϕ j−1 ,ϕ j+1 ,ϕ j+3 ). For two-dimensional
|     |     |     |     |     |     |     |     | 2   |     | 2 2 | 2   |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
problems, all the direct neighbours of a cell are considered while calculating ϕ. Entropy distance and the
shock sensor value are plotted for a sod tube test case, and as seen in figure 6a, it matches with the shock
location. Figure 6b shows the highlighted cells where the shock is present for an oblique shock reflection test
case.
Through numerical experiments, it was found that the entropy conservative scheme does not add sufficient
1.2
1
0.8
0.6
0.4
0.2
00
|     | 0.1 | 0.2 0.3 | 0.4 0x.5 | 0.6 0.7 | 0.8 0.9 1 |     |     |     |     |     |     |     |
| --- | --- | ------- | -------- | ------- | --------- | --- | --- | --- | --- | --- | --- | --- |
|     |     |         |          | (a)     |           |     |     |     | (b) |     |     |     |
Figure6: (a): ShocksensorvalueforSODtubetestcase. (b): Shocksensorvalueandcontourofdensityvariationforoblique
shockreflectiontestcase(cellswithshocksmarkedwithyellow).
stabilizingdiffusioninsmoothregionsoncoarsemeshes. Additionalfourth-ordernumericaldiffusionisadded
based on the JST scheme [17] to prevent oscillations. The diffusion is given as
1
|     |           |             |     |          | FR = | α (U  | −3U | +3U | −U  | )   |     | (49) |
| --- | --------- | ----------- | --- | -------- | ---- | ----- | --- | --- | --- | --- | --- | ---- |
|     |           |             |     |          | 1    | R j+2 |     | j+1 | j   | j−1 |     |      |
|     |           |             |     |          | j+ 2 | 2     |     |     |     |     |     |      |
| The | diffusion | coefficient |     | is given | by   |       |     |     |     |     |     |      |
1
λ˜
|     |     |     |     |     |     | α R | =   | j+1 |     |     |     | (50) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | ---- |
|     |     |     |     |     |     |     | 32  | 2   |     |     |     |      |
andλ˜ isthecoefficientofnumericaldiffusionbasedonthatoftheRiemanninvariantsbasedexactcontact
j+1
2
discontinuity capturing scheme (RICCA) [33]. The following equation defines it.
|u
|     |     |     |     |       | |+|u | |     |     |     |         |     |        |      |
| --- | --- | --- | --- | ----- | ------ | --- | --- | --- | ------- | --- | ------ | ---- |
|     |     |     |     |  j+1 | j      | if  | |F  | −F  | |≤δ and | |U  | −U |≤δ |      |
|     |     | λ˜  |     |       |        |     | j+1 |     | j       |     | j+1 j  |      |
|     |     | j+1 | =   |       | 2      |     |     |     |         |     |        | (51) |
2
|     |     |     |     |  max(|u | |,|u  | |)+sign(|∆p|)a |     |     | otherwise |     |     |     |
| --- | --- | --- | --- | -------- | ----- | -------------- | --- | --- | --------- | --- | --- | --- |
|     |     |     |     |          | j+1 j |                |     | j+1 |           |     |     |     |
2
Here δ is small parameter taken as 10−16. The above diffusion and the entropy conservative scheme result in
a numerical scheme that can capture steady contact discontinuities exactly. The final scheme can be written
as
|     |     |     |     |     | F=FEC | +(1−ϕ)FR+ϕFRH |     |     |     |     |     | (52) |
| --- | --- | --- | --- | --- | ----- | ------------- | --- | --- | --- | --- | --- | ---- |
Theabovehybridschemeisspatiallysecond-orderaccurateinsmoothregionsasthecentralentropyconserva-
tivefluxissecond-orderaccurate,asshowninsection4.2. Inthefollowingsections,theschemesgivenby(41)
and (52) will be referred to as ES (Entropy stable) and HES (Hybrid entropy stable) schemes, respectively.
The following sections test them for various one-dimensional and two-dimensional test cases.
15

1.2
|     |     |     |     |          |     | 1.5 |     |     |     | 3.5      |     |     |
| --- | --- | --- | --- | -------- | --- | --- | --- | --- | --- | -------- | --- | --- |
|     |     |     |     | Exact    |     |     |     |     |     | Exact    |     |     |
| 1   |     |     |     | Roe EC   |     |     |     |     |     |          |     |     |
|     |     |     |     | PC ECKEP |     |     |     |     |     | Roe EC   |     |     |
|     |     |     |     | ECKEP    |     |     |     |     |     | PC ECKEP |     |     |
ECKEP
| 0.8 |     |     |     |     |     | 1   |     |     |     | 3   |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
Exact
| 0.6 |     |     |     |     |     |     | Roe EC   |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | -------- | --- | --- | --- | --- | --- |
|     |     |     |     |     | U   |     | PC ECKEP |     | e   |     |     |     |
ECKEP
| 0.4 |     |     |     |     |     | 0.5 |     |     |     | 2.5 |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
0.2
| 0   |     |     |     |     |     | 0     |         |     |     | 2     |             |     |
| --- | --- | --- | --- | --- | --- | ----- | ------- | --- | --- | ----- | ----------- | --- |
| 0   | 0.2 | 0.4 | 0.6 | 0.8 | 1   | 0 0.2 | 0.4 0.6 | 0.8 | 1   | 0 0.2 | 0.4 0.6 0.8 | 1   |
|     |     |     | x   |     |     |       | x       |     |     |       | x           |     |
Figure7: Plotsofdensity,velocityandinternalenergyforvariousschemesforthefirsttestcase
| 7. Numerical |                 | Results |         |     |     |     |     |     |     |     |     |     |
| ------------ | --------------- | ------- | ------- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 7.1.         | One Dimensional |         | Results |     |     |     |     |     |     |     |     |     |
Numerous1-DEulertestcaseshavebeensolvedtodeterminetherobustnessandaccuracyoftheentropy
stable scheme given in (52). Note that EC3 flux, kinetic energy preserving, was taken as the entropy con-
servative flux in ES (41) and HES (52) schemes. For the one-dimensional problems, 100 points with a CFL
number of 0.1 is taken. Neumann boundary conditions are applied applied on both ends. Time-step is com-
puted using ∆t=CFL× ∆x . The initial conditions, final time and location of initial discontinuity
max(|uj|+aj)
∀j
| for all | test cases | are | given  | in table | 5.   |         |           |          |           |                |       |     |
| ------- | ---------- | --- | ------ | -------- | ---- | ------- | --------- | -------- | --------- | -------------- | ----- | --- |
|         | No.        | x   | ρ      |          | u    | p       | ρ         | u        |           | p              | t     |     |
|         |            | o   | L      |          | L    | L       | R         | R        |           | R              | f     |     |
|         | 1          | 0.3 | 1      |          | 0.75 | 1       | 0.125     | 0        |           | 0.1            | 0.2   |     |
|         | 2          | 0.5 | 1      |          | 0    | 1000    | 1         | 0        |           | 0.01           | 0.012 |     |
|         | 3          | 0.4 | 5.9924 | 19.5975  |      | 460.894 | 5.9924    | -6.19633 |           | 46.0950        | 0.035 |     |
|         |            |     |        |          |      |         | γ +1 PR+1 |          | (cid:113) |                |       |     |
|         | 4          | 0.5 | 1      |          | 1    | 1       | γ − 1     | PR       |           | γ(2+(γ−1)M2)uR | 5     |     |
|         |            |     |        |          |      | γM2     | γ +1      | γM2      |           | (2γM2+(1−γ))ρR |       |     |
γ − 1 +PR
|     | 5   | 0.5 | 1.4  |     | 0       | 1                               | 1   | 0     |     | 1   | 2   |     |
| --- | --- | --- | ---- | --- | ------- | ------------------------------- | --- | ----- | --- | --- | --- | --- |
|     | 6   | 0.1 | 3.86 |     | -0.81   | 10.33                           | 1   | -3.44 |     | 1   | 4   |     |
|     | 7   | .5  | 1.4  |     | 0.1     | 1.0                             | 1.0 | 0.1   |     | 1.0 | 1   |     |
|     |     |     |      |     | Table5: | Initialconditionsfor1Dtestcases |     |       |     |     |     |     |
Test case 1 is a sod shock tube problem with a right-moving shock, right-moving contact discontinuity and
a left-moving rarefaction with a sonic point. The ES scheme is compared with the Roe EC and PC ECKEP
schemes. ECKEP works well with this test case as shock and contact discontinuity are captured with less
diffusion, and no expansion shock or sonic glitch is observed in the expansion fan region. It is comparable
to the existing entropy stable and kinetic preserving schemes. Test cases 2 and 3 involve shocks with high-
pressuregradientsandareusedtotesttherobustnessofthenumericalschemes. Theschemesareabletogive
non-oscillatory results for these test cases. Test case 4 depicts a stationary shock, captured with one interior
pointbytheESscheme. Testcase5representsastationarycontactdiscontinuitywhichthenewschemescan
capture exactly. Test case 6 deals with a slowly moving shock wave moving to the right. Numerical schemes
withlownumericaldiffusionsufferfrompost-shockoscillationsforthistest. However,thenewentropystable
schemes capture the shock without oscillations. Test case 7 represents a slowly moving contact discontinuity
moving to the right. The entropy stable schemes are seen to capture the contact discontinuity with minimal
| diffusion | and | without | oscillations. |     |     |     |     |     |     |     |     |     |
| --------- | --- | ------- | ------------- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
16

| 1   |     |     |                | 1.4 |     |     | 3.5            |     |     |     |
| --- | --- | --- | -------------- | --- | --- | --- | -------------- | --- | --- | --- |
|     |     |     | Exact Solution |     |     |     | Exact Solution |     |     |     |
|     |     |     | HES Scheme     |     |     |     | HES Scheme     |     |     |     |
| 0.9 |     |     |                | 1.2 |     |     |                |     |     |     |
0.8
1
| 0.7 |     |     |     |     |     |     | 3   |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0.6 |     |     |     | 0.8 |     |     |     |     |     |     |
Exact Solution
|     |     |     |     | U   |     | HES Scheme | e   |     |     |     |
| --- | --- | --- | --- | --- | --- | ---------- | --- | --- | --- | --- |
| 0.5 |     |     |     | 0.6 |     |            |     |     |     |     |
| 0.4 |     |     |     |     |     |            | 2.5 |     |     |     |
0.4
0.3
| 0.2 |     |     |     | 0.2 |     |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
0.1 0 0.2 0.4 0.6 0.8 1 0 0 0.2 0.4 0.6 0.8 1 2 0 0.2 0.4 0.6 0.8 1
|     |     | x   |     |     | x   |     |     | x   |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
Figure8: Plotsofdensity,velocityandinternalenergyusingHESschemeforthefirsttestcase
| 6   |     |     |     | 20  |     |     | 2500 |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | ---- | --- | --- | --- |
Exact Solution
|     |     |     |     | 18  |     |     |     |     | ES Scheme |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --------- | --- |
5
|     |     |     | Exact Solution | 16  |     |     | 2000 |     |     |     |
| --- | --- | --- | -------------- | --- | --- | --- | ---- | --- | --- | --- |
ES Scheme
14
4
|     |     |     |     | 12   |     |                | 1500 |     |     |     |
| --- | --- | --- | --- | ---- | --- | -------------- | ---- | --- | --- | --- |
| 3   |     |     |     | U 10 |     |                | e    |     |     |     |
|     |     |     |     | 8    |     | Exact Solution | 1000 |     |     |     |
| 2   |     |     |     |      |     | ES Scheme      |      |     |     |     |
6
|     |     |     |     | 4   |     |     | 500 |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
1
2
| 0   |                |     |         | 0       |     |         | 0       |     |                |     |
| --- | -------------- | --- | ------- | ------- | --- | ------- | ------- | --- | -------------- | --- |
| 0   | 0.2            | 0.4 | 0.6 0.8 | 1 0 0.2 | 0.4 | 0.6 0.8 | 1 0 0.2 | 0.4 | 0.6 0.8        | 1   |
|     |                | x   |         |         | x   |         |         | x   |                |     |
| 6   |                |     |         | 20      |     |         | 2500    |     |                |     |
|     | Exact Solution |     |         |         |     |         |         |     | Exact Solution |     |
|     | HES Scheme     |     |         | 18      |     |         |         |     | HES Scheme     |     |
5
|     |     |     |     | 16  |     |     | 2000 |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | ---- | --- | --- | --- |
14
| 4   |     |     |     |      |     | Exact Solution |      |     |     |     |
| --- | --- | --- | --- | ---- | --- | -------------- | ---- | --- | --- | --- |
|     |     |     |     | 12   |     | HES Scheme     | 1500 |     |     |     |
| 3   |     |     |     | U 10 |     |                | e    |     |     |     |
|     |     |     |     | 8    |     |                | 1000 |     |     |     |
2
6
|     |     |     |     | 4   |     |     | 500 |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
1
2
| 0   |     |     |         | 0       |     |         | 0       |     |         |     |
| --- | --- | --- | ------- | ------- | --- | ------- | ------- | --- | ------- | --- |
| 0   | 0.2 | 0.4 | 0.6 0.8 | 1 0 0.2 | 0.4 | 0.6 0.8 | 1 0 0.2 | 0.4 | 0.6 0.8 | 1   |
|     |     | x   |         |         | x   |         |         | x   |         |     |
Figure9: Plotsofdensity,velocityandinternalenergyusingES(top)andHES(bottom)schemesforthethirdtestcase
17

| 35  |     | 20  |     |     | 300 |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- |
Exact Solution
ES Scheme
| 30  |     | 15  |     |     | 250 |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 25  |     | 10  |     |     | 200 |     |     |
Exact Solution
| 20  |     | U 5 |     |     | e 150 ES Scheme |     |     |
| --- | --- | --- | --- | --- | --------------- | --- | --- |
| 15  |     | 0   |     |     | 100             |     |     |
Exact Solution
ES Scheme
| 10    |             | -5      |     |         | 50      |             |     |
| ----- | ----------- | ------- | --- | ------- | ------- | ----------- | --- |
| 5     |             | -10     |     |         | 0       |             |     |
| 0 0.2 | 0.4 0.6 0.8 | 1 0 0.2 | 0.4 | 0.6 0.8 | 1 0 0.2 | 0.4 0.6 0.8 | 1   |
|       | x           |         | x   |         |         | x           |     |
| 35    |             | 20      |     |         | 300     |             |     |
Exact Solution
HES Scheme
| 30    |                | 15      |     |         | 250     |                |     |
| ----- | -------------- | ------- | --- | ------- | ------- | -------------- | --- |
| 25    |                | 10      |     |         | 200     |                |     |
| 20    |                | U 5     |     |         | e 150   |                |     |
| 15    |                | 0       |     |         | 100     |                |     |
|       | Exact Solution |         |     |         |         | Exact Solution |     |
|       | HES Scheme     |         |     |         |         | HES Scheme     |     |
| 10    |                | -5      |     |         | 50      |                |     |
| 5     |                | -10     |     |         | 0       |                |     |
| 0 0.2 | 0.4 0.6 0.8    | 1 0 0.2 | 0.4 | 0.6 0.8 | 1 0 0.2 | 0.4 0.6 0.8    | 1   |
|       | x              |         | x   |         |         | x              |     |
Figure10: Plotsofdensity,velocityandinternalenergyusingES(top)andHES(bottom)schemesforthefourthtestcase
| 2.8 |     | 1   |     |     | 0.8 |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- |
2.6
|     |     | 0.9 |     |     | 0.75 |     |     |
| --- | --- | --- | --- | --- | ---- | --- | --- |
2.4
0.7
| 2.2 |                | 0.8 |     |     |      |                |     |
| --- | -------------- | --- | --- | --- | ---- | -------------- | --- |
|     | Exact Solution |     |     |     |      | Exact Solution |     |
|     | ES Scheme      |     |     |     | 0.65 | ES Scheme      |     |
| 2   |                | 0.7 |     |     |      |                |     |
E x a ct  S o lu tion
| 1.8 |     | U   |     | E S  S c he m e | e 0.6 |     |     |
| --- | --- | --- | --- | --------------- | ----- | --- | --- |
0.6
| 1.6 |     |     |     |     | 0.55 |     |     |
| --- | --- | --- | --- | --- | ---- | --- | --- |
| 1.4 |     | 0.5 |     |     |      |     |     |
0.5
1.2
|     |     | 0.4 |     |     | 0.45 |     |     |
| --- | --- | --- | --- | --- | ---- | --- | --- |
1
| 0.8   |             | 0.3     |     |         | 0.4     |             |     |
| ----- | ----------- | ------- | --- | ------- | ------- | ----------- | --- |
| 0 0.2 | 0.4 0.6 0.8 | 1 0 0.2 | 0.4 | 0.6 0.8 | 1 0 0.2 | 0.4 0.6 0.8 | 1   |
| 2.8   | x           | 1.1     | x   |         | 0.8     | x           |     |
| 2.6   |             | 1       |     |         | 0.75    |             |     |
Exact Solution
| 2.4 |                |     |     | HES Scheme |      |                |     |
| --- | -------------- | --- | --- | ---------- | ---- | -------------- | --- |
|     |                | 0.9 |     |            | 0.7  |                |     |
| 2.2 | Exact Solution |     |     |            |      | Exact Solution |     |
|     | HES Scheme     | 0.8 |     |            | 0.65 | HES Scheme     |     |
2
| 1.8 |     | U 0.7 |     |     | e 0.6 |     |     |
| --- | --- | ----- | --- | --- | ----- | --- | --- |
| 1.6 |     | 0.6   |     |     | 0.55  |     |     |
1.4
|     |     | 0.5 |     |     | 0.5 |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- |
1.2
|     |     | 0.4 |     |     | 0.45 |     |     |
| --- | --- | --- | --- | --- | ---- | --- | --- |
1
| 0.8   |             | 0.3     |     |         | 0.4     |             |     |
| ----- | ----------- | ------- | --- | ------- | ------- | ----------- | --- |
| 0 0.2 | 0.4 0.6 0.8 | 1 0 0.2 | 0.4 | 0.6 0.8 | 1 0 0.2 | 0.4 0.6 0.8 | 1   |
|       | x           |         | x   |         |         | x           |     |
Figure 11: Plots of density, velocity and internal energy using ES (top) and HES (bottom) schemes for the stationary shock
wave.
18

| 1.4  |                | 1   |     |                | 2.5 |                |     |
| ---- | -------------- | --- | --- | -------------- | --- | -------------- | --- |
|      | Exact Solution |     |     | Exact Solution |     | Exact Solution |     |
|      | ES Scheme      | 0.8 |     | ES Scheme      |     | ES Scheme      |     |
| 1.35 |                |     |     |                | 2.4 |                |     |
0.6
| 1.3 |     |     |     |     | 2.3 |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- |
0.4
| 1.25 |     | 0.2 |     |     | 2.2 |     |     |
| ---- | --- | --- | --- | --- | --- | --- | --- |
| 1.2  |     | U 0 |     |     | 2.1 |     |     |
e
-0.2
| 1.15 |     |     |     |     | 2   |     |     |
| ---- | --- | --- | --- | --- | --- | --- | --- |
-0.4
| 1.1 |     |     |     |     | 1.9 |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- |
-0.6
| 1.05  |             | -0.8    |     |         | 1.8     |             |     |
| ----- | ----------- | ------- | --- | ------- | ------- | ----------- | --- |
| 1     |             | -1      |     |         | 1.7     |             |     |
| 0 0.2 | 0.4 0.6 0.8 | 1 0 0.2 | 0.4 | 0.6 0.8 | 1 0 0.2 | 0.4 0.6 0.8 | 1   |
|       | x           |         | x   |         |         | x           |     |
Figure12: Plotsofdensity,velocityandinternalenergyusingESschemeforstationarycontactwave.
| 4   |     | -0.5 |     |     | 7   |     |     |
| --- | --- | ---- | --- | --- | --- | --- | --- |
6.5
| 3.5 |     | -1  |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- |
6
| 3   |                |      |     |                | 5.5 |     |     |
| --- | -------------- | ---- | --- | -------------- | --- | --- | --- |
|     | Exact Solution | -1.5 |     | Exact Solution |     |     |     |
|     | ES Scheme      |      |     | ES Scheme      | 5   |     |     |
2.5
|     |     | U -2 |     |     | e 4.5 |     |     |
| --- | --- | ---- | --- | --- | ----- | --- | --- |
| 2   |     |      |     |     | 4     |     |     |
Exact Solution
|     |     | -2.5 |     |     | 3.5 | ES Scheme |     |
| --- | --- | ---- | --- | --- | --- | --------- | --- |
1.5
3
| 1   |     | -3  |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- |
2.5
| 0.5   |             | -3.5    |     |         | 2       |             |     |
| ----- | ----------- | ------- | --- | ------- | ------- | ----------- | --- |
| 0 0.2 | 0.4 0.6 0.8 | 1 0 0.2 | 0.4 | 0.6 0.8 | 1 0 0.2 | 0.4 0.6 0.8 | 1   |
|       | x           |         | x   |         |         | x           |     |
| 4     |             | -0.5    |     |         | 7       |             |     |
6.5
| 3.5 |                | -1  |     |     |     |     |     |
| --- | -------------- | --- | --- | --- | --- | --- | --- |
|     | Exact Solution |     |     |     | 6   |     |     |
HES Scheme
| 3   |     | -1.5 |     |     | 5.5 |     |     |
| --- | --- | ---- | --- | --- | --- | --- | --- |
5
| 2.5 |     | -2   |     |                          |     |                |     |
| --- | --- | ---- | --- | ------------------------ | --- | -------------- | --- |
|     |     | U    |     | E x a c t  S o l u tio n | 4.5 | Exact Solution |     |
|     |     |      |     | H E S   S c h e m e      | e   | HES Scheme     |     |
| 2   |     | -2.5 |     |                          | 4   |                |     |
3.5
| 1.5 |     | -3  |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- |
3
| 1     |             | -3.5    |     |         | 2.5     |             |     |
| ----- | ----------- | ------- | --- | ------- | ------- | ----------- | --- |
| 0.5   |             | -4      |     |         | 2       |             |     |
| 0 0.2 | 0.4 0.6 0.8 | 1 0 0.2 | 0.4 | 0.6 0.8 | 1 0 0.2 | 0.4 0.6 0.8 | 1   |
|       | x           |         | x   |         |         | x           |     |
Figure13: Plotsofdensity,velocityandinternalenergyusingES(top)andHES(bottom)schemesslowlymovingshockwave.
19

|     | 1.45 |     |     |     | 0.2  |     |     |                |     | 2.5 |                |     |     |     |
| --- | ---- | --- | --- | --- | ---- | --- | --- | -------------- | --- | --- | -------------- | --- | --- | --- |
|     |      |     |     |     |      |     |     | Exact Solution |     |     | Exact Solution |     |     |     |
|     | 1.4  |     |     |     | 0.18 |     |     | ES Scheme      |     |     | ES Scheme      |     |     |     |
2.4
|     | 1.35 Exact Solution |     |     |     | 0.16 |     |     |     |     |     |     |     |     |     |
| --- | ------------------- | --- | --- | --- | ---- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
|     | ES Scheme           |     |     |     |      |     |     |     |     | 2.3 |     |     |     |     |
|     | 1.3                 |     |     |     | 0.14 |     |     |     |     |     |     |     |     |     |
|     | 1.25                |     |     |     | 0.12 |     |     |     |     | 2.2 |     |     |     |     |
|     | 1.2                 |     |     | U   | 0.1  |     |     |     |     | 2.1 |     |     |     |     |
e
|     | 1.15 |     |     |     | 0.08 |     |     |     |     |     |     |     |     |     |
| --- | ---- | --- | --- | --- | ---- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
2
|     | 1.1 |     |     |     | 0.06 |     |     |     |     |     |     |     |     |     |
| --- | --- | --- | --- | --- | ---- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
1.9
|     | 1.05  |     |         |     | 0.04 |     |     |                |     |     |                |         |     |     |
| --- | ----- | --- | ------- | --- | ---- | --- | --- | -------------- | --- | --- | -------------- | ------- | --- | --- |
|     | 1     |     |         |     | 0.02 |     |     |                |     | 1.8 |                |         |     |     |
|     | 0.95  |     |         |     | 0    |     |     |                |     | 1.7 |                |         |     |     |
|     | 0 0.2 | 0.4 | 0.6 0.8 | 1   | 0    | 0.2 | 0.4 | 0.6 0.8        | 1   | 0   | 0.2            | 0.4 0.6 | 0.8 | 1   |
|     |       | x   |         |     |      |     | x   |                |     |     |                | x       |     |     |
|     | 1.45  |     |         |     | 0.6  |     |     |                |     | 2.6 |                |         |     |     |
|     |       |     |         |     |      |     |     | Exact Solution |     |     | Exact Solution |         |     |     |
|     | 1.4   |     |         |     | 0.5  |     |     | HES Scheme     |     | 2.5 | HES Scheme     |         |     |     |
|     | 1.35  |     |         |     | 0.4  |     |     |                |     |     |                |         |     |     |
2.4
|     | 1.3 |     |     |     | 0.3 |     |     |     |     |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
2.3
|     | 1.25 |     |     |     | 0.2 |     |     |     |     |     |     |     |     |     |
| --- | ---- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
2.2
|     | 1.2  |     | E x a c t  S o l u | tio n U | 0.1 |     |     |     |     | e   |     |     |     |     |
| --- | ---- | --- | ------------------ | ------- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
|     |      |     | H E S   S c h e m  | e       |     |     |     |     |     | 2.1 |     |     |     |     |
|     | 1.15 |     |                    |         | 0   |     |     |     |     |     |     |     |     |     |
2
|     | 1.1   |     |         |     | -0.1 |     |     |         |     |     |     |         |     |     |
| --- | ----- | --- | ------- | --- | ---- | --- | --- | ------- | --- | --- | --- | ------- | --- | --- |
|     | 1.05  |     |         |     | -0.2 |     |     |         |     | 1.9 |     |         |     |     |
|     | 1     |     |         |     | -0.3 |     |     |         |     | 1.8 |     |         |     |     |
|     | 0.95  |     |         |     | -0.4 |     |     |         |     | 1.7 |     |         |     |     |
|     | 0 0.2 | 0.4 | 0.6 0.8 | 1   | 0    | 0.2 | 0.4 | 0.6 0.8 | 1   | 0   | 0.2 | 0.4 0.6 | 0.8 | 1   |
|     |       | x   |         |     |      |     | x   |         |     |     |     | x       |     |     |
Figure14: Plotsofdensity,velocityandinternalenergyusingES(top)andHES(bottom)schemesslowlymovingcontactwave.
| 7.2. | Two Dimensional |     | Results |     |     |     |     |     |     |     |     |     |     |     |
| ---- | --------------- | --- | ------- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
This section presents solutions to several two-dimensional problems using the ES scheme (41) and HES
scheme(52). Thediffusionoperator(49)iscomputedutilizingtheformulationforunstructuredgridsoutlined
in [34]. These test cases serve to validate the solutions for both accuracy and robustness. A CFL number of
0.1 is utilized for all instances, and the time-step is determined using the equation:
CFL×Area
|     |     |     | ∆t= | min         |     |        |     |        |        |     |     |     |     | (53) |
| --- | --- | --- | --- | ----------- | --- | ------ | --- | ------ | ------ | --- | --- | --- | --- | ---- |
|     |     |     |     | ∀Cells(|V·n |     | |+c)∆s |     | +(|V·n | |+c)∆s |     |     |     |     |      |
|     |     |     |     |             |     | η      |     | η      | ζ      |     | ζ   |     |     |      |
Here,ηandζ representthegridcoordinatedirections,andthecorrespondingnormalsn andn arecomputed
|     |     |     |     |     |     |     |     |     |     |     |     | ζ η |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
by averaging the normals from opposite sides of the control volume. For more details, see chapter 6 of [35].
All steady test cases are run until all residuals (density, x-momentum, y-momentum, and total energy) fell
below 10−8. The value of q in equations 47 is taken as 10, and ϵ in 48 is taken between 0.1 and 0.5.
| 7.2.1. | Oblique | Shock Reflection |     |     |     |     |     |     |     |     |     |     |     |     |
| ------ | ------- | ---------------- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
As described by Yee [36], this test case involves an oblique shock encountering a wall, resulting in a
reflected shock wave. The free stream Mach number (M ) is 2.9, and the incident shock angle is 29o.
∞
The domain is discretized with Cartesian mesh cells in the [0,3]×[0,1] region. Supersonic inflow boundary
conditions are applied on the left boundary with inflow values (W = [1.0,2.9,0,1.0/1.4]) and post-shock
values are prescribed at the top boundary. A flow tangency boundary condition is used at the bottom
boundary (solid wall), while a supersonic outflow boundary condition is applied on the right boundary. The
flow field is initialized with the inflow conditions. Simulations are conducted with grid sizes of 240x80 and
480x160, and density contours are depicted in Figure 15 for both schemes. As illustrated, the incident and
reflected shocks are accurately resolved, with the shocks being less diffusive for higher-order schemes. Semi-
logplotsofdensity,momentum,andenergyresidualsversusthenumberofiterationsaredisplayedinFigures
16 and 17, demonstrating good convergence to steady-state solutions for both schemes.
| 7.2.2. | Supersonic | flow | over a | compression | ramp |     |     |     |     |     |     |     |     |     |
| ------ | ---------- | ---- | ------ | ----------- | ---- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
This test case, as studied by Levy et al. [37], involves a supersonic flow of Mach 2 over a compression
rampof150. Thedomainforthisproblemis[0,3]×[0,1]. Supersonicinflowboundaryconditionsareimposed
20

on the left boundary with primitive variables W = [1.0,2,0,1.0/1.4]. Flow tangency boundary conditions
are applied on the top and bottom boundaries, while a supersonic outflow boundary condition is imposed on
therightboundary. Theentiredomainisinitializedwiththeinflowconditions. Thissteadyproblemexhibits
an initial shock emanating from a ramp, further getting reflected from the top and bottom boundaries. An
expansion fan forms at the end of the ramp and further interacts with the reflected shock. Figure 18 shows
pressurecontours,capturingtheincidentshock,reflectedshock,anexpansionfanandwaveinteractionswith
both ES and HES schemes.
7.2.3. Hypersonic flow over a half-cylinder
This test case involves hypersonic flow at Mach 20 over a half-cylindrical body, resulting in a pronounced
bowshock. Anundesirablephenomenonknownasacarbuncleshock[38]isobservedinlow-diffusivenumerical
schemes like Roe scheme, with an unphysical perturbation disrupting the bow shock on the stagnation line.
The problem becomes more severe when the mesh is aligned with the bow shock. This undesirable feature,
which often occurs with Riemann solvers, has attracted considerable research and is counted together with
several shock instabilities that plague the Riemann solvers [7]. For generating the mesh, the mesh curvature
is determined using the formulation by Huang et al. [39], given by the equations:
(cid:112)
2x − 4x2−4(1+tan(θ)2)(x2−r2))
x= c c c i
2(1+tan(θ)2) (54)
y =−tan(θ)x
where θ = (j −1)∗5π/(6M)−5π/12, x = 1.8(N −i+1)/N, and r = 1+2.4(N −i+1)/N. Here,
c i
N and M are chosen as 40 and 320, respectively. The initial and inflow conditions (on the left boundary)
are W = [1.0,20,0,1.0/1.4]. A flow tangency boundary condition is imposed on the cylinder surface, while
a supersonic outflow boundary condition is applied on the remaining periphery. Figure 19 depicts density
contours obtained using different schemes. As illustrated in Figure 19a, the Roe scheme exhibits a carbuncle
shock, disrupting the shock structure, while the ES and HES schemes, shown in Figures 19b and 19c, avoid
this artifact.
7.2.4. Forward-facing step in supersonic flow
This test case, introduced by Emery [40], comprises a step facing a supersonic flow of Mach 3. Domain
is taken as [0,3]×[0,1] and inflow and initial conditions are: W = [1.0,3.0,0.0,1.0/1.4]. The height of the
step is 0.2 units and is located at a distance of 0.6 from the left boundary. Top and bottom boundaries
are solid walls, with flow tangency boundary conditions being applied. At the left boundary supersonic
inflow boundary condition is applied, and at the right boundary supersonic outflow boundary condition is
enforced. This is an unsteady test case; the final solution is taken at t = 4s, at which point the shocks
move very slowly. Flow features include a bow shock and its reflection from top and bottom boundaries.
Also, a powerful expansion fan is formed at the tip of the step, which interacts with the reflected shocks.
A slipstream that travels downstream is also formed at the λ-shock (at the triple-point), approximately at
y =0.8. Thecornerisasingularityintheflowfield,anditinducesanumericalboundarylayerthatruinsthe
solution downstream. ”Corner fix” given in [41] is used to fix this issue. Results are shown in figure 20, and
both the schemes can capture the discontinuities and wave interactions with good accuracy, with the HES
being better than the ES scheme. No expansion shock is seens with either of the schemes in the region near
the corner. For the ES scheme, the Mach stem is observed on the lower wall in the final solution. However,
no Mach stem is observed for the HES scheme, and the Mach stem at the top wall is also captured at the
correct location (x = 0.6). Slipstream at y = 0.8 just after the Mach stem on the top surface is also clearly
captured, which is a feature associated with low diffusion schemes.
7.2.5. Shock diffraction problem
ThistestcaseinvolvesaMach5.09flowdiffractingovera90o cornerofabackward-facingstep,forminga
strong expansion wave. Density and pressure after the expansion are minimal, and thus, schemes that fail to
21

preserve positivity fail in this test case, such as Roe’s approximate Riemann solver [7]. Domain for this test
caseistakenas[0,1]×[0,1]withthecorneroflength0.05unitslocatedaty =0.6. Thedomainisinitialized
withW=[1.4,0,0,1.0]forx>0.05andthepost-shockconditionsforx<0.05. Supersonicinflowboundary
condition is imposed on the left side from y =0.6 to 1.0 and flow tangecy boundary condition is imposed on
the surface of the step. Symmetry boundary condition is imposed on all the other boundaries. The solution
is computed at a final time of T = 0.1561s. Numerous schemes that suffer from shock instabilities do not
preservetheshockstructureinthistestcase,andanunphysicalexpansionshockoftenappearsatthecorner.
Also, the rarefaction from the corner creates near-zero density and pressures, and schemes fail to preserve
theirpositivity. ResultsfromtheESandHESschemescanbeseeninfigure21fordifferentgrids. Theplanar
shock is preserved, and the entropy stable schemes do not produce expansion shocks or any other anomalies.
| 7.2.6. Odd-even | decoupling |     |     |     |     |     |
| --------------- | ---------- | --- | --- | --- | --- | --- |
Another simulation anomaly is the grid-aligned planar shock structure deterioration caused by odd-even
decoupling. Given in [7], in this test case, a planar shock travels in a long rectangular tube. The centre-line
| of the grid | is perturbed | slightly, given | as  |     |     |     |
| ----------- | ------------ | --------------- | --- | --- | --- | --- |
(cid:40)
|     |     |     | y i,j +0.1 | for | i even, |      |
| --- | --- | --- | ---------- | --- | ------- | ---- |
|     |     |     | y =        |     |         | (55) |
i,j
|     |     |     | y i,j −0.1 | for | i odd. |     |
| --- | --- | --- | ---------- | --- | ------ | --- |
Domain is taken as [0,2400] × [0,20] and is partly shown in figure 22a. The domain is initialized with
W=[1.0,0,0,1.0] for x>5, and post-shock conditions corresponding to Mach 20 are specified for x<5. It
refers to a stronger shock than that given in [7]. Flow tangency boundary conditions are applied at the top
andbottomboundaries,andasupersonicoutflowboundaryconditionwasappliedattherightboundary. The
odd-evendecouplingphenomenoncausestheshockstructuretobreakdown,andmanyschemes(likeRiemann
solvers) fail before reaching the final time of T = 330s. This shock instability worsens with increasing the
shock strength and magnitude of perturbations [42]. Results with ES and HES schemes are shown in figures
22b and 22c; with both the schemes, shock structure remains intact. This test case and the hypersonic flow
over the half-cylinder (section 7.2.3) indicate that the schemes are free of any numerical shock instabilities.
| 7.2.7. Double | Mach | reflection |     |     |     |     |
| ------------- | ---- | ---------- | --- | --- | --- | --- |
Thisisanotherunsteadytestcase[41]whereashockofMach10meetsawallatanangleof30o,resulting
in a complicated shock structure including a reflected shock and a slipstream. Domain for this test case is
| taken to | be [0,4]×[0,1]. | The initial | conditions for | the problem | are |     |
| -------- | --------------- | ----------- | -------------- | ----------- | --- | --- |
√
(cid:40)
|     |     | [1.4,0,0,1.0] |                   |     | if y ≤ 3(x−1/6) |      |
| --- | --- | ------------- | ----------------- | --- | --------------- | ---- |
|     |     | W=            | √                 |     |                 | (56) |
|     |     | [8.0,33       | 3/8,−4.125,116.5] |     | otherwise.      |      |
This is an unsteady benchmark test case introduced by Woodward & Colella [41] to evaluate shock-
60°,
capturing schemes. A Mach 10 shock, inclined at travels through air (γ = 1.4) and reflects off a solid
wallthatbeginsatx=1/6alongthebottomboundary. Theshockisintroducednotbyasimplediscontinuity,
but through boundary-driven conditions: post-shock values are enforced on the bottom-left segment and the
topboundaryisdynamicallyupdatedtomatchtheshock’sangleandspeed. Thissetupproducesacomplex,
self-similar flow involving two triple points — each forming a Mach stem, a reflected shock, and a contact
discontinuity (slipstream). One of the key challenges is resolving the high-density jet formed near the wall,
which resembles shaped charge behavior and requires high-order methods like PPM to accurately capture.
Lower-order methods typically fail to resolve this jet and the secondary reflected shock, highlighting the
| importance | of accurate | shock capturing | and contact | discontinuity | resolution. |     |
| ---------- | ----------- | --------------- | ----------- | ------------- | ----------- | --- |
Density contours of the solutions obtained using ES and HES schemes are given in figure 23. Prominent
flow features like the Mach reflections and stems are captured well by both schemes. However, only the
higher order scheme is able to capture the complicated flow features such as the secondary reflected shock
| and jet | with sufficient | detail as shown | in the density | map in | figure 24. |     |
| ------- | --------------- | --------------- | -------------- | ------ | ---------- | --- |
22

7.2.8. Shock - vortex filament interaction
Insection7.2.3, thenewschemesaredemonstratednotproducinganyunphysicalcarbuncleshock, which
is ultimately a numerical artifact. Recent studies have demonstrated that carbuncles can be triggered under
certain physical conditions. Elling [43] shows that a vortex filament interacting with a steady shock of high
strengths will trigger a carbuncle structure that grows with time. This is true for many schemes, such as the
Gudunov, Lax-Friedrichs, and Osher schemes. Kemm [44] shows that schemes such as HLLE have very high
numerical shear viscosity, preventing the formation of the physical carbuncles induced by vortex filament.
Sincethistestcaseisamodelfortheshock-boundarylayerinteractionproblem(vortexfilamentrepresentinga
numericalboundarylayer),correctcapturingoftheinducedcarbuncleisessentialandmustnotbesuppressed.
Thedomainforthisproblemistakentobe[0,200]x[0,100]withasteadyshockofM=20atx=100. Attheleft
boundary, supersonic inflow conditions are imposed corresponding to W = 1.0,20.0,0,1.0/1.4 everywhere
except at a single cell in the centre where W=1.0,0,0,1.0/1.4 is imposed. The top and bottom boundaries
are taken as solid walls, and the right boundary is considered to represent supersonic outflow. Simulation is
run for a time of T = 20s. As shown in figure 25, the physical carbuncle is captured for both schemes, and
the flow features are much more pronounced for the HES scheme than for the ES scheme. Compared with
the LLF scheme in figure 25 (a), the carbuncle is captured more precisely with ES and HES schemes.
7.2.9. Flows NACA 0012 airfoil
The new schemes are tested for flows over a symmetric NACA 0012 airfoil, in trans-sonic and supersonic
regimes. The domain is taken as a circle of radius of 10 units centered at the trailing edge of an airfoil with
chord of unit length. A structured mesh of 200x300 elements with 200 points on the airfoil surface is used
for the simulation (see fig 26). Grid is stretched in radial direction by a 2.5%. Far-field boundary conditions
are imposed on the outer boundary, and flow tangency wall boundary conditions are imposed on the airfoil
surface. Flow field is initialized with W = 1.4,Mcos(θ),Msin(θ),1.0/1.4, where M is the mach number of
flowandθ istheangleofattack. Flowiscomputedforthefollowingconditions: M=0.85,θ =2o andM=1.2,
θ =0. Thesecomputationsareonlyperformedforthehybridentropystableflux(52). Pressureandpressure
coefficient C = (p−p )/(1ρ u2 ) for both test cases are shown in figure 28. For the transonic case, the
p ∞ 2 ∞ ∞
shocks on upper and lower surfaces can be seen to be captured accurately in these plots. For the supersonic
csse, crisply captured bow shock and fish-tail shocks can be observed. Pressure contours for both schemes
are shown in figure 27 and the C plots are shown in figure 28.
P
23

(a)Densitycontoursfor(30in0.9-2.7)for240x80gridwithES(left)andHES(right)Scheme
(b)Densitycontoursfor(30in0.9-2.7)for480x160gridwithES(left)andHES(right)Scheme
|                | Figure15:           | Obliqueshockreflectiontestcase |                |                     |     |
| -------------- | ------------------- | ------------------------------ | -------------- | ------------------- | --- |
| 10-4           |                     |                                | 10-4           |                     |     |
|                | Mass Residual       |                                |                | Mass Residual       |     |
|                | X-Momentum Residual |                                |                | X-Momentum Residual |     |
| 10-6           | Y - M o m           | e n tu m   R e sidual          | 10-6           | Y-Momentum Residual |     |
|                | E n er g            | y  R e s id u a l              |                | Energy Residual     |     |
| slaudiseR 10-8 |                     |                                | slaudiseR 10-8 |                     |     |
| 10-10          |                     |                                | 10-10          |                     |     |
| 10-12          |                     |                                | 10-12          |                     |     |
| 10-14          |                     |                                | 10-14          |                     |     |
| 0 1            | 2 3                 | 4 5                            | 0 1            | 2 3                 | 4 5 |
|                | Time                |                                |                | Time                |     |
Figure16: Semi-logplotofresidualsforESschemefor240x80grid(left)and480x160grid(right)
| 10-4 |                     |     | 10-4 |                     |     |
| ---- | ------------------- | --- | ---- | ------------------- | --- |
|      | Mass Residual       |     |      | Mass Residual       |     |
|      | X-Momentum Residual |     |      | X-Momentum Residual |     |
|      | Y-Momentum Residual |     |      | Y-Momentum Residual |     |
|      | Energy Residual     |     |      | Energy Residual     |     |
10-5
| slaudiseR |     |     | slaudiseR |     |     |
| --------- | --- | --- | --------- | --- | --- |
10-5
10-6
10-6
10-7
| 0 1 | 2 3  | 4 5 | 0 1 | 2 3  | 4 5 |
| --- | ---- | --- | --- | ---- | --- |
|     | Time |     |     | Time |     |
Figure17: Semi-logplotofresidualsforHESschemefor240x160grid(left)and480x160grid(right)
24

Figure18: Pressurecontours(50in0.6-2.4)for15o compressionrampwithES(left)andHES(right)schemes
| (a) | (b) |     | (c) |
| --- | --- | --- | --- |
Figure19: Densitycontours(30in1.3-9.0)forhypersonicflowovercylinderon320x40gridwith(a)Roescheme,(b)ESscheme
and(c)HESscheme
Figure20: Densitycontours(30in0.3-6.8)forforwardstepon480x160gridusingES(left)andHES(right)scheme
| (a) | (b) | (c) | (d) |
| --- | --- | --- | --- |
Figure 21: Density contours (30 in 0.05-7.1) for backwards-facing step. (a) and (b) using the ES scheme on 400x400 and
1200x1200grids,respectively. (c)and(d)usingtheHESschemeon400x400and1200x1200grids,respectively
25

(a) (b) (c)
Figure22: (a): Partofthegrid(x=0to35)usedforodd-evendecouplingwithcenterlineperturbations. (b): Contoursofdensity
(30 in 0.9-44) for odd-even decoupling using the ES scheme. (c): Contours of density (30 in 0.9-44) for odd-even decoupling
usingHESscheme
Figure 23: Density contours (30 in 1.3-22) for double mach reflection test case on a grid of 960x240 using ES (left) and HES
(right)schemes.
Figure24: Colorplotofdensityforthedoublemachreflectiontestcaseonagridof960x240usingtheHESscheme.
26

| (a) | (b) | (c) |
| --- | --- | --- |
Figure25: entropy(s=ln(p/ργ))contours(15in-0.6-4.2)forshockvortexfilamentinteractionwith(a)LLF(b)ES(c)HES
schemes
Figure26: StructuredmeshforNACA0012airfoil
27

(a)TransoniccasewithM=0.85and2o angleofattack. (b)TransoniccasewithM=1.2and0o angleofattack. Pressure
Pressurecontours(30in0.4-1.6)usingHESScheme. contours(30in0.5-2.4)usingHESScheme.
Figure27: NACA0012Airfoiltestcase
| -0.6 |     | 1.6      |                | -0.4   |     | 2.4          |                |
| ---- | --- | -------- | -------------- | ------ | --- | ------------ | -------------- |
|      |     |          | Top Surface    |        |     |              | Top Surface    |
|      |     |          | Bottom Surface | -0.2   |     |              | Bottom Surface |
| -0.4 |     | 1.4      |                |        |     |              |                |
|      |     |          |                | 0      |     | 2            |                |
| -0.2 |     | 1.2      |                | 0.2    |     |              |                |
|      |     | erusserP |                |        |     | erusserP 1.6 |                |
| pC   |     |          |                | pC 0.4 |     |              |                |
| 0    |     | 1        |                | 0.6    |     |              |                |
| 0.2  |     | 0.8      |                | 0.8    |     | 1.2          |                |
1
| 0.4 |                | 0.6 |     |     |                |     |       |
| --- | -------------- | --- | --- | --- | -------------- | --- | ----- |
|     | Top Surface    |     |     | 1.2 | Top Surface    | 0.8 |       |
|     | Bottom Surface |     |     |     | Bottom Surface |     |       |
| 0.6 |                | 0.4 |     |     |                |     |       |
| 0   | 0.5            | 1 0 | 0.5 | 1 0 | 0.5            | 1 0 | 0.5 1 |
|     | x/c            |     | x/c |     | x/c            |     | x/c   |
|     | (a)            |     | (b) |     | (c)            |     | (d)   |
Figure28: PressureandCp plotsalongtopandbottomsurfacesofairfoilusingHESscheme: (a)and(b)Transoniccase. (c)
and(d)Supersoniccase.
28

| 8. Summary |     | and Conclusions |     |     |     |     |     |     |
| ---------- | --- | --------------- | --- | --- | --- | --- | --- | --- |
A family of structure-preserving numerical fluxes for the Euler equations that conserve entropy and
preservekineticenergyinasemi-discretesenseareintroducedinthiswork. Beginningwithaderivationbased
ontheinherentstructureswithintheenergyequation,newentropy-conservativeandkineticenergy-preserving
fluxes that do not require logarithmic averages are proposed, thereby improving computational efficiency.
These fluxes maintain second-order accuracy and are validated on representative smooth problems such as
the Taylor-Green vortex and isentropic vortex convection, demonstrating their effectiveness in conserving
| physical | invariants | in  | the absence | of shocks. |     |     |     |     |
| -------- | ---------- | --- | ----------- | ---------- | --- | --- | --- | --- |
To extend applicability to problems involving discontinuities, a scalar numerical diffusion mechanism
based on the Rankine–Hugoniot conditions is introduced. This diffusion strategy ensures entropy stability
and kinetic energy damping while preserving steady contact discontinuities exactly. Furthermore, to reduce
numerical diffusion in smooth regions, a hybrid entropy-stable (HES) scheme is developed using an entropy-
distance-basedshocksensortoselectivelyswitchbetweenlowandhighdiffusionmodes. Theresultingscheme
balances accuracy and stability across both smooth and discontinuous regimes.
Extensiveone-andtwo-dimensionalnumericalexperimentsconfirmtherobustness,accuracy,andstability
of the proposed methods. Notably, the schemes successfully handle classical and challenging test cases such
as odd-even decoupling, double Mach reflection, shock-vortex interaction, and supersonic flow over airfoils,
without exhibiting common numerical artifacts like expansion shocks or carbuncle-like instabilities.
| CRediT | author | statement |     |     |     |     |     |     |
| ------ | ------ | --------- | --- | --- | --- | --- | --- | --- |
KunalBahuguna: Conceptualization,Methodology,Software,Validation,Formalanalysis,Writing—orig-
| inal draft, | Visualization. |     |     |     |     |     |     |     |
| ----------- | -------------- | --- | --- | --- | --- | --- | --- | --- |
Ramesh Kolluru: Conceptualization, Validation, Investigation, Writing — review & editing, Supervision.
S. V. Raghurama Rao: Conceptualization,Validation,Investigation,Resources,Writing—review&edit-
ing, Supervision.
| Declaration | of  | competing |     | interest |     |     |     |     |
| ----------- | --- | --------- | --- | -------- | --- | --- | --- | --- |
The authors declare that they have no known financial interests or personal relationships with any other
| people or | organisations |     | that could | influence |     | the work presented | here. |     |
| --------- | ------------- | --- | ---------- | --------- | --- | ------------------ | ----- | --- |
9. Appendices
| 9.1. Appendix |     | A: Entropy | conservative |     | flux using | an optimisation | procedure |     |
| ------------- | --- | ---------- | ------------ | --- | ---------- | --------------- | --------- | --- |
Entropy-conservative fluxes given in section 3 can also be obtained from an optimisation problem defined
as follows.
|     |     |     |     | min||Fc−F|| |     |           | ∆V·Fc |       |
| --- | --- | --- | --- | ----------- | --- | --------- | ----- | ----- |
|     |     |     |     |             |     | such that | =∆ψ   | (A-1) |
The above quadratic optimisation with an equality constraint can be analytically solved using the method
| of Lagrange | multipliers. |     | The Lagrange |                             | function | can be given | as  |       |
| ----------- | ------------ | --- | ------------ | --------------------------- | -------- | ------------ | --- | ----- |
|             |              |     |              | L(Fc,λ)=||Fc−F||+λ(V·Fc−∆ψ) |          |              |     | (A-2) |
whereλistheLagrangemultipliercorrespondingtotheentropyconservativeconstraint. Thesolutiontothis
| optimisation | problem |     | can be found | using | equation |     |     |     |
| ------------ | ------- | --- | ------------ | ----- | -------- | --- | --- | --- |
∇·L(Fc,λ)=0
(A-3)
which gives
λ
Fc
|     |     |     |     |     |     | =F− ∆V |     | (A-4a) |
| --- | --- | --- | --- | --- | --- | ------ | --- | ------ |
2
λ ∆V·F−∆ψ
|     |     |     |     |     |     | =       |     | (A-4b) |
| --- | --- | --- | --- | --- | --- | ------- | --- | ------ |
|     |     |     |     |     |     | 2 ∆V·∆V |     |        |
29

Note that the λ/2 here is exactly the dissipation term in EC1 flux (26). A similar procedure can be followed
to obtain an entropy conservative and kinetic energy preserving flux. The optimisation problem can be
| formulated | with | two constraints |     | as follows. |           |       |         |     |        |       |
| ---------- | ---- | --------------- | --- | ----------- | --------- | ----- | ------- | --- | ------ | ----- |
|            |      | min||Fc−F||     |     |             |           | ∆V·Fc |         | Fc  | =Fcu+p |       |
|            |      |                 |     |             | such that |       | =∆ψ and |     |        | (A-5) |
2 1
This, however, results in a complicated flux. To obtain a simple flux we can modify the above problem as
|     |     |     |     | min||Fc−F |     | || such | that ∆V·Fc | =∆ψ |     | (A-6) |
| --- | --- | --- | --- | --------- | --- | ------- | ---------- | --- | --- | ----- |
|     |     |     |     |           | 3 3 |         |            |     |     |       |
where we assume Fc =F and Fc =Fcu+p. Lagrange function for this can be given as
|          |           | 1        | 1        | 2   | 1           |              |               |     |     |       |
| -------- | --------- | -------- | -------- | --- | ----------- | ------------ | ------------- | --- | --- | ----- |
|          |           |          |          | L(F | c,λ)=||Fc−F |              | ||+λ(V·Fc−∆ψ) |     |     | (A-7) |
|          |           |          |          | 3   |             | 3            | 3             |     |     |       |
| with the | following | solution | obtained | by  | taking      | ∇·L(Fc,λ)=0. |               |     |     |       |
3
λ
Fc
|     |     |     |     |     |     | =F − | ∆V  |     |     | (A-8a) |
| --- | --- | --- | --- | --- | --- | ---- | --- | --- | --- | ------ |
|     |     |     |     |     |     | 3 3  | 2 3 |     |     |        |
λ ∆V·F−∆ψ
|                  |       |            |          |               |               | =     |            |          |           | (A-8b) |
| ---------------- | ----- | ---------- | -------- | ------------- | ------------- | ----- | ---------- | -------- | --------- | ------ |
|                  |       |            |          |               |               | 2 ∆V  | ·∆V        |          |           |        |
|                  |       |            |          |               |               |       | 3 3        |          |           |        |
| This corresponds |       | to the     | ECKEP    | flux          | (30) obtained | in    | section 3. |          |           |        |
| 9.2. Appendix    |       | B: Entropy | Distance | Positivity    |               |       |            |          |           |        |
| For              | small | changes    | we can   | write entropy | distance      | ∆V·∆U | as         |          |           |        |
|                  |       |            |          |               | (cid:18)      |       | (cid:19)T  | (cid:18) | (cid:19)T |        |
|                  |       |            |          |               |               | dV    |            | dV       |           |        |
|                  |       |            | ED       | =dVTdU=       |               | ·dU   | dU=dUT     |          | dU        | (B-1)  |
|                  |       |            |          |               |               | dU    |            | dU       |           |        |
d2η(U)
We know that V = dη(U) and dV = which is a positive definite matrix because of convexity of η(U)
|           |     | dU     |            | dU  | dU2     |     |     |     |     |     |
| --------- | --- | ------ | ---------- | --- | ------- | --- | --- | --- | --- | --- |
| and hence | ED  | ≥0 for | any states | U 1 | and U 2 | .   |     |     |     |     |
References
[1] E. Tadmor, The Numerical Viscosity of Entropy Stable Schemes for Systems of Conservation Laws. I,
Mathematics of Computation 49 (179) (1987) 91–103. doi:https://doi.org/10.2307/2008251.
[2] E.Tadmor,Entropystabilitytheoryfordifferenceapproximationsofnonlinearconservationlawsandre-
lated time-dependent problems, Acta Numerica 12 (2003) 451–512. doi:10.1017/S0962492902000156.
URL https://www.cambridge.org/core/product/identifier/S0962492902000156/type/journal_
article
[3] T. J. Barth, Numerical Methods for Gasdynamic Systems on Unstructured Meshes, in: M. Griebel,
D. E. Keyes, R. M. Nieminen, D. Roose, T. Schlick, D. Kr¨oner, M. Ohlberger, C. Rohde (Eds.), An
Introduction to Recent Developments in Theory and Numerics for Conservation Laws, Vol. 5, Springer
Berlin Heidelberg, Berlin, Heidelberg, 1999, pp. 195–285, series Title: Lecture Notes in Computational
| Science | and                                                  | Engineering. |     | doi:10.1007/978-3-642-58535-7_5. |     |     |     |     |     |     |
| ------- | ---------------------------------------------------- | ------------ | --- | -------------------------------- | --- | --- | --- | --- | --- | --- |
| URL     | http://link.springer.com/10.1007/978-3-642-58535-7_5 |              |     |                                  |     |     |     |     |     |     |
[4] Roe P L, Affordable, entropy consistent flux functions, in: talk presented at the Eleventh International
Conference on Hyperbolic Problems: Theory, Numerics, Applications, Lyon, 2006.
URL https://scholar.google.co.in/citations?view_op=view_citation&hl=en&user=
4fNzp4oAAAAJ&cstart=20&pagesize=80&citation_for_view=4fNzp4oAAAAJ:Tiz5es2fbqcC
30

[5] H. Ranocha, Comparison of some Entropy Conservative Numerical Fluxes for the Euler Equations,
Journal of Scientific Computing 76 (1) (2018) 216–242, arXiv:1701.02264 [math]. doi:10.1007/
s10915-017-0618-1.
URL http://arxiv.org/abs/1701.02264
[6] P. G. LeFloch, C. Rohde, High-Order Schemes, Entropy Inequalities, and Nonclassical Shocks, SIAM
Journal on Numerical Analysis 37 (6) (2000) 2023–2060. doi:10.1137/S0036142998345256.
URL http://epubs.siam.org/doi/10.1137/S0036142998345256
[7] J. J. Quirk, A contribution to the great Riemann solver debate, International Journal for Numerical
Methods in Fluids 18 (6) (1994) 555–574. doi:10.1002/fld.1650180603.
URL https://onlinelibrary.wiley.com/doi/10.1002/fld.1650180603
[8] S.Pirozzoli,NumericalMethodsforHigh-SpeedFlows,AnnualReviewofFluidMechanics43(1)(2011)
163–194. doi:10.1146/annurev-fluid-122109-160718.
URL https://www.annualreviews.org/doi/10.1146/annurev-fluid-122109-160718
[9] A. Jameson, Formulation of Kinetic Energy Preserving Conservative Schemes for Gas Dynamics and
Direct Numerical Simulation of One-Dimensional Viscous Compressible Flow in a Shock Tube Using
Entropy and Kinetic Energy Preserving Schemes, Journal of Scientific Computing 34 (2) (2008) 188–
208. doi:10.1007/s10915-007-9172-6.
URL http://link.springer.com/10.1007/s10915-007-9172-6
[10] H. Ranocha, Entropy Conserving and Kinetic Energy Preserving Numerical Methods for the Euler
Equations Using Summation-by-Parts Operators, in: S. J. Sherwin, D. Moxey, J. Peir´o, P. E. Vincent,
C. Schwab (Eds.), Spectral and High Order Methods for Partial Differential Equations ICOSAHOM
2018, Springer International Publishing, Cham, 2020, pp. 525–535.
[11] P. Chandrashekar, Kinetic energy preserving and entropy stable finite volume schemes for compressible
EulerandNavier-Stokesequations,CommunicationsinComputationalPhysics14(5)(2013)1252–1286,
arXiv:1209.4994 [cs, math]. doi:10.4208/cicp.170712.010313a.
URL http://arxiv.org/abs/1209.4994
[12] R. Abgrall, P. Offner, H. Ranocha, Reinterpretation and extension of entropy correction terms
for residual distribution and discontinuous Galerkin schemes: Application to structure pre-
serving discretization, Journal of Computational Physics 453 (110955) (2022). doi:https:
//doi.org/10.1016/j.jcp.2022.110955.
URL https://www.sciencedirect.com/science/article/abs/pii/S0021999122000171?via%
3Dihub
[13] P. L. Roe, Approximate Riemann solvers, parameter vectors, and difference schemes, Journal of Com-
putational Physics 43 (2) (1981) 357–372. doi:10.1016/0021-9991(81)90128-5.
URL https://www.sciencedirect.com/science/article/pii/0021999181901285
[14] F. Ismail, P. L. Roe, H. Nishikawa, A Proposed Cure to the Carbuncle Phenomenon, in: H. Deconinck,
E. Dick (Eds.), Computational Fluid Dynamics 2006, Springer Berlin Heidelberg, Berlin, Heidelberg,
2009, pp. 149–154. doi:10.1007/978-3-540-92779-2_21.
URL http://link.springer.com/10.1007/978-3-540-92779-2_21
[15] F.Ismail,P.L.Roe,Affordable,entropy-consistentEulerfluxfunctionsII:Entropyproductionatshocks,
Journal of Computational Physics 228 (15) (2009) 5410–5436. doi:10.1016/j.jcp.2009.04.021.
URL https://linkinghub.elsevier.com/retrieve/pii/S0021999109002113
31

[16] A. Harten, The Artificial Compression Method for Computation of Shocks and Contact Discontinuities:
III. Self-Adjusting Hybrid Schemes, Mathematics of Computation 32 (142) (1978) 363–389, publisher:
American Mathematical Society. doi:10.2307/2006149.
URL http://www.jstor.org/stable/2006149
[17] A.Jameson,W.Schmidt,E.Turkel,NumericalsolutionoftheEulerequationsbyfinitevolumemethods
using Runge Kutta time stepping schemes, in: 14th Fluid and Plasma Dynamics Conference, American
Institute of Aeronautics and Astronautics, Palo Alto,CA,U.S.A., 1981, AIAA Paper No. AIAA-81-1259.
doi:10.2514/6.1981-1259.
URL https://arc.aiaa.org/doi/10.2514/6.1981-1259
[18] U. S. Fjordholm, S. Mishra, E. Tadmor, Arbitrarily High-order Accurate Entropy Stable Essentially
Nonoscillatory Schemes for Systems of Conservation Laws, SIAM Journal on Numerical Analysis 50 (2)
(2012) 544–573. doi:10.1137/110836961.
URL http://epubs.siam.org/doi/10.1137/110836961
[19] K. O. Friedrichs, P. D. Lax, Systems of Conservation Equations with a Convex Extension, Proceedings
oftheNationalAcademyofSciencesoftheUnitedStatesofAmerica68(8)(1971)1686–1688,publisher:
National Academy of Sciences.
URL https://www.jstor.org/stable/61263
[20] M. Mock, Systems of conservation laws of mixed type, Journal of Differential Equations 37 (1) (1980)
70–88. doi:10.1016/0022-0396(80)90089-3.
URL https://linkinghub.elsevier.com/retrieve/pii/0022039680900893
[21] A. Harten, On the symmetric form of systems of conservation laws with entropy, Journal of Computa-
tional Physics 49 (1) (1983) 151–164. doi:10.1016/0021-9991(83)90118-3.
URL https://linkinghub.elsevier.com/retrieve/pii/0021999183901183
[22] G. J. Gassner, A. R. Winters, D. A. Kopriva, Split Form Nodal Discontinuous Galerkin Schemes with
Summation-By-PartsPropertyfortheCompressibleEulerEquations,JournalofComputationalPhysics
327 (2016) 39–66, arXiv:1604.06618 [math]. doi:10.1016/j.jcp.2016.09.013.
URL http://arxiv.org/abs/1604.06618
[23] R. Abgrall, A general framework to construct schemes satisfying additional conservation relations. Ap-
plication to entropy conservative and entropy dissipative schemes, Journal of Computational Physics
372 (2018) 640–666. doi:10.1016/j.jcp.2018.06.031.
URL https://linkinghub.elsevier.com/retrieve/pii/S0021999118304091
[24] B. Leer, Flux-Vector Splitting for the 1990s, 1990, NASA Lewis Research Center Computational Fluid
Dynamics Symposium on AeropropulsionAt: Cleaveland, OH, Volume: N91-21062 13-02.
URL https://www.researchgate.net/publication/241365933_Flux-vector_Splitting_for_the_
1990
[25] C.-W. Shu, S. Osher, Efficient implementation of essentially non-oscillatory shock-capturing schemes,
Journal of Computational Physics 77 (2) (1988) 439–471. doi:10.1016/0021-9991(88)90177-5.
URL https://linkinghub.elsevier.com/retrieve/pii/0021999188901775
[26] H. C. Yee, B. Sj¨ogreen, On Entropy Conservation and Kinetic Energy Preservation Methods, Journal of
Physics: Conference Series 1623 (1) (2020) 012020. doi:10.1088/1742-6596/1623/1/012020.
URL https://iopscience.iop.org/article/10.1088/1742-6596/1623/1/012020
[27] A. Gouasmi, S. Murman, K. Duraisamy, Entropy conservative schemes and the receding flow
problem, Journal of Scientific Computing 78 (2019) 971–994. doi:https://doi.org/10.1007/
s10915-018-0793-8.
URL https://link.springer.com/article/10.1007/s10915-018-0793-8#citeas
32

[28] S. C. Spiegel, H. Huynh, J. R. DeBonis, A Survey of the Isentropic Euler Vortex Problem using High-
Order Methods, in: 22nd AIAA Computational Fluid Dynamics Conference, American Institute of
Aeronautics and Astronautics, Dallas, TX, 2015. doi:10.2514/6.2015-2444.
URL https://arc.aiaa.org/doi/10.2514/6.2015-2444
[29] S.Jaisankar,S.V.RaghuramaRao,AcentralRankine–Hugoniotsolverforhyperbolicconservationlaws,
Journal of Computational Physics 228 (3) (2009) 770–798. doi:10.1016/j.jcp.2008.10.002.
URL https://linkinghub.elsevier.com/retrieve/pii/S0021999108005123
[30] D. Zaide, P. Roe, Entropy-Based Mesh Refinement, II: A New Approach to Mesh Movement, in: 19th
AIAA Computational Fluid Dynamics, American Institute of Aeronautics and Astronautics, San Anto-
nio, Texas, 2009. doi:10.2514/6.2009-3791.
URL https://arc.aiaa.org/doi/10.2514/6.2009-3791
[31] K. Shrinath, N. Maruthi, S.V. Raghurama Rao, V. Vasudeva Rao, A Kinetic Flux Difference Splitting
method for compressible flows, Computers & Fluids 250 (2023) 105702. doi:10.1016/j.compfluid.
2022.105702.
URL https://linkinghub.elsevier.com/retrieve/pii/S004579302200295X
[32] S. S. Roy, S.V. Raghurama Rao, A kinetic scheme with variable velocities and relative entropy, Com-
puters & Fluids 265 (2023) 106016. doi:10.1016/j.compfluid.2023.106016.
URL https://linkinghub.elsevier.com/retrieve/pii/S0045793023002414
[33] R. Kolluru, N. V. Raghavendra, S.V. Raghurama Rao, G. N. Sekhar, Novel, simple and robust contact-
discontinuity capturing schemes for high speed compressible flows, Applied Mathematics and Computa-
tion 414 (2022) 126660, arXiv:2003.10695 [cs, math]. doi:10.1016/j.amc.2021.126660.
URL http://arxiv.org/abs/2003.10695
[34] A. Jameson, Origins and Further Development of the Jameson-Schmidt-Turkel Scheme (Invited), in:
33rd AIAA Applied Aerodynamics Conference, American Institute of Aeronautics and Astronautics,
Dallas, TX, 2015. doi:10.2514/6.2015-2718.
URL https://arc.aiaa.org/doi/10.2514/6.2015-2718
[35] J. Blazek, Computational fluid and solid mechanics, Elsevier, Amsterdam, [Netherlands], 2001.
[36] H. C. Yee, R. F. Warming, A. Harten, A high-resolution numerical technique for inviscid gas-dynamic
problems with weak solutions, in: E. Krause (Ed.), Eighth International Conference on Numerical
Methods in Fluid Dynamics, Vol. 170, Springer Berlin Heidelberg, Berlin, Heidelberg, 1982, pp. 546–
552, series Title: Lecture Notes in Physics. doi:10.1007/3-540-11948-5_72.
URL http://link.springer.com/10.1007/3-540-11948-5_72
[37] D. W. Levy, K. G. Powell, B. van Leer, Use of a rotated Riemann solver for the two-dimensional Euler
equations, Journal of Computational Physics 106 (1993) 201–214, aDS Bibcode: 1993JCoPh.106..201L.
doi:10.1006/jcph.1993.1103.
URL https://ui.adsabs.harvard.edu/abs/1993JCoPh.106..201L
[38] K. Peery, S. Imlay, Blunt-body flow simulations, in: 24th Joint Propulsion Conference, American Insti-
tute of Aeronautics and Astronautics, Boston,MA,U.S.A., 1988. doi:10.2514/6.1988-2904.
URL https://arc.aiaa.org/doi/10.2514/6.1988-2904
[39] K. Huang, H. Wu, H. Yu, D. Yan, Cures for numerical shock instability in HLLC solver, International
Journal for Numerical Methods in Fluids 65 (9) (2011) 1026–1038. doi:10.1002/fld.2217.
URL https://onlinelibrary.wiley.com/doi/10.1002/fld.2217
33

[40] A. F. Emery, An evaluation of several differencing methods for inviscid fluid flow problems, Journal of
Computational Physics 2 (3) (1968) 306–331. doi:10.1016/0021-9991(68)90060-0.
URL https://www.sciencedirect.com/science/article/pii/0021999168900600
[41] P. Woodward, P. Colella, The numerical simulation of two-dimensional fluid flow with strong shocks,
Journal of Computational Physics 54 (1) (1984) 115–173. doi:10.1016/0021-9991(84)90142-6.
URL https://linkinghub.elsevier.com/retrieve/pii/0021999184901426
[42] N. Fleischmann, S. Adami, X. Y. Hu, N. A. Adams, A low dissipation method to cure the grid-aligned
shock instability, Journal of Computational Physics 401 (2020) 109004. doi:10.1016/j.jcp.2019.
109004.
URL https://linkinghub.elsevier.com/retrieve/pii/S0021999119307090
[43] V.Elling, Thecarbunclephenomenonisincurable, ActaMathematicaScientia29(6)(2009)1647–1656.
doi:10.1016/S0252-9602(10)60007-0.
URL https://linkinghub.elsevier.com/retrieve/pii/S0252960210600070
[44] F.Kemm,Heuristicalandnumericalconsiderationsforthecarbunclephenomenon,AppliedMathematics
and Computation 320 (2018) 596–613. doi:10.1016/j.amc.2017.09.014.
URL https://linkinghub.elsevier.com/retrieve/pii/S0096300317306355
34