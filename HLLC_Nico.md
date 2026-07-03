Journal of Computational Physics 423 (2020) 109762

Contents lists available at ScienceDirect

Journal  of  Computational  Physics

www.elsevier.com/locate/jcp

A  shock-stable  modiﬁcation  of  the  HLLC  Riemann  solver  with
reduced  numerical  dissipation

Nico Fleischmann

∗

,  Stefan Adami,  Nikolaus  A. Adams

Technical University of Munich, Department of Mechanical Engineering, Chair of Aerodynamics and Fluid Mechanics, Boltzmannstraße 15,
85748 Garching, Germany

a  r  t  i  c  l  e

i  n  f  o

a  b  s  t  r  a  c  t

Article history:
Received 1 February 2020
Received in revised form 31 July 2020
Accepted 4 August 2020
Available online 26 August 2020

Keywords:
Shock instability
Carbuncle phenomenon
HLLC
High-order schemes
WENO
Low-dissipation schemes

The  purpose  of  this  paper  is  twofold.  First,  the  application  of  high-order  methods
in  combination  with  the  popular  HLLC  Riemann  solver  demonstrates  that  the  grid-
aligned  shock  instability  can  strongly  affect  simulation  results  when  the  grid  resolution
is  increased.  Beyond  the  well-documented  two-dimensional  behavior,  the  problem  is
particularly  troublesome  with  three-dimensional  simulations.  Hence,  there  is  a  need  for
shock-stable modiﬁcations of HLLC-type solvers for high-speed ﬂow simulations.
Second,  the  paper  provides  a  stabilization  of  the  popular  HLLC  ﬂux  based  on  a  recently
proposed  mechanism  for  grid  aligned-shock  instabilities  Fleischmann  et al.  (2020)  [8].
The  instability  was  found  to  be  triggered  by  an  inappropriate  scaling  of  acoustic  and
advection dissipation for local low Mach numbers. These low Mach numbers occur during
the calculation of ﬂuxes in transverse direction of the shock propagation, where the local
velocity  component  vanishes.  A centralized  formulation  of  the  HLLC  ﬂux  is  provided  for
this purpose, which allows for a simple reduction of nonlinear signal speeds. In contrast to
other shock-stable versions of the HLLC ﬂux, the resulting HLLC-LM ﬂux reduces the inherent
numerical dissipation of the scheme.
The  robustness  of  the  proposed  scheme  is  tested  for  a  comprehensive  range  of  cases
involving strong shock waves. Three-dimensional single- and multi-component simulations
are performed with high-order methods to demonstrate that the HLLC-LM ﬂux also copes
with latest challenges of compressible high-speed computational ﬂuid dynamics.

© 2020 The Author(s). Published by Elsevier Inc. This is an open access article under the
CC BY-NC-ND license (http://creativecommons.org/licenses/by-nc-nd/4.0/).

1.  Introduction

Approximate Riemann solvers in combination with shock-capturing Godunov schemes [1] dominate modern computation
of  phenomena  that  involve  complex  ﬂow  interactions  across  scales  such  as  shock  interaction  with  multi-phase  interfaces
and turbulent scales. The application of high-order discretizations allows for an accurate prediction of many of such ﬂows.
However,  over  the  last  decades  the  grid-aligned  shock  instability  has  presented  a  barrier  for  robust  computation  of  high
Mach number ﬂows using high-order discretizations with state-of-the-art low-dissipation Riemann solvers such as Roe [2]
or HLLC [3,4]. Since the ﬁrst description of the problem by Peery and Imlay [5] and Quirk [6] an extensive research on the
topic resulted in a large number of scientiﬁc publications addressing various aspects. A summary of major developments to

* Corresponding author.

E-mail addresses: nico.ﬂeischmann@tum.de (N. Fleischmann), stefan.adami@tum.de (S. Adami), nikolaus.adams@tum.de (N.A. Adams).

https://doi.org/10.1016/j.jcp.2020.109762
0021-9991/© 2020 The Author(s). Published by Elsevier Inc. This is an open access article under the CC BY-NC-ND license
(http://creativecommons.org/licenses/by-nc-nd/4.0/).

N. Fleischmann, S. Adami and N.A. Adams

Journal of Computational Physics 423 (2020) 109762

the present day can be found in [7,8]. Even though most of the research has focused on the Roe solver, also the HLLC solver
is aﬄicted by the instability.

HLL-type  solvers  were  originally  developed  by  Harten,  Lax  and  van  Leer  [3].  In  combination  with  the  nonlinear  signal
speed  estimates  of  Einfeldt  [9] and  the  restoration  of  the  contact  wave  proposed  by  Toro  et al.  [4],  the  resulting  HLLC
Riemann approximation became one of the most successful and widespread Riemann solvers for hyperbolic systems [10,11].
An  accurate  estimation  of  the  contact  wave  speed  was  communicated  by  Batten  et al.  [12].  Due  to  the  explicit  modeling
of  each  wave  of  the  governing  Euler  equations,  HLLC  is  a  complete  Riemann  solver  with  signiﬁcantly  reduced  dissipation
near contact discontinuities compared to the HLL scheme. The design of the HLLC ﬂux allows for straightforward extensions
to other types of hyperbolic equations, e.g. for magneto-hydrodynamics [13–15], by introduction of additional wave types.
Moreover,  the  HLLC  ﬂux  has  been  applied  successfully  to  multi-component  ﬂows  [16,17],  and  capillary  forces  have  been
introduced to simulate surface tension effects at liquid/gas interfaces [18]. Further recent applications are reviewed in [11].
While  the  HLLC  ﬂux  is  known  to  suffer  from  the  shock  instability,  the  stable  behavior  of  the  HLL  ﬂux  was  described
already  by  Quirk  [6].  He  suggested  to  apply  the  HLL  scheme  near  strong  shocks  in  combination  with  lower-dissipation
schemes, such as HLLC, in the remaining domain. These hybrid schemes lead to stable, but nevertheless contact preserving
results.  The  switching  procedure  was  improved  by  Kim  et al.  [19],  where  the  dissipative  HLL  ﬂux  only  is  applied  for  the
ﬂuxes in transverse direction of the shock propagation. Another modiﬁcation of the hybrid scheme was suggested in [20],
where the dissipative HLL ﬂux only is applied for two components of the ﬂux. However, hybrid schemes may still signif-
icantly  increase  dissipation,  and  a  switching  procedure  has  to  be  provided.  Additionally,  the  authors  of  [20] successfully
tested the shock stability of the rotated Riemann solver method [21] applied to the HLLC ﬂux, but they found that the latter
approach is computationally rather expensive. The ﬁrst pure HLLC-type ﬂux with shock-stable properties, called HLLCM, was
developed by Shen et al. [22] via smearing of the shear velocities on both sides of the contact line. This procedure introduces
shear viscosity and stabilizes the calculation of strong shocks. However, the introduced amount of dissipation limits the ac-
curacy of boundary layer calculations and therefore the authors again suggested to apply a hybrid HLLC-HLLCM version for
complex ﬂows. Recently, Xie et al. [23] proposed an HLLC-type Riemann solver with an additional pressure-dissipation term
that is activated near shocks and damps spurious pressure perturbations. Simon and Mandal [24,25] proposed two different
approaches to avoid the shock instability. They separated the HLLC ﬂux into the inherent HLL part and an antidiffusive part.
In their ﬁrst approach [24], the activation of the antidiffusive term is controlled by a pressure-ratio-based multi-dimensional
shock sensor. The resulting solver called HLLC-ADC restores the shock stability of the HLL ﬂux. The second approach [25] is
to apply a selective wave modiﬁcation that increases the inherent dissipative HLL part in the vicinity of a shock wave. The
antidiffusive term of the resulting HLLC-SWM ﬂux remains identical to that of the original HLLC.

In comparison to the large number of proposed modiﬁcations of the Roe ﬂux, the grid-aligned shock instability of the
popular HLLC solver has found much less consideration in literature. The reason is probably, that the solution of most two-
dimensional  simulations  remains  bounded,  and  therefore  the  effect  of  the  introduced  disturbances  is  not  as  catastrophic
as  with  the  Roe  ﬂux.  However,  with  increased  resolution,  high-order  discretizations,  and  extension  to  three-dimensional
simulations, the application of the HLLC ﬂux is prone to develop severe carbuncles, similarly to that obtained with the Roe
ﬂux, as is shown in this paper.

In [8], the authors proposed a new possible mechanism of the grid-aligned shock instability. A wrong scaling behavior of
numerical dissipation due to the local low Mach number in transverse direction of the shock front propagation was found to
cause the numerical shock instability. A modiﬁcation for the popular Roe ﬂux and the local componentwise Lax-Friedrichs
ﬂux was proposed that proved to be shock stable. The present paper proposes a new shock-stable modiﬁcation of the HLLC
ﬂux called HLLC-LM that is based on these ﬁndings. As a straightforward reduction of nonlinear wave speeds is not sensible
for the classical HLLC formulation, a new centralized reformulation of the HLLC ﬂux is derived. This alternative formulation
allows for an analogous reduction of acoustic dissipation as with the modiﬁed Roe scheme without introducing additional
diﬃculties. Most of the present shock-stabilizing variants of the HLLC ﬂux restore the shock stability by adding additional
dissipation in one way or the other, as motivated by the stability of the stable, but highly dissipative HLL scheme. In contrast,
the proposed HLLC-LM ﬂux with less numerical dissipation than the classical HLLC ﬂux represents a fundamentally different
approach  in  comparison  to  the  earlier  HLLC-HLL  combination  models.  Moreover,  in  this  paper  the  shock  stability  of  both
HLLC  and  HLLC-LM  is  studied  using  high-order  methods  in  space  and  time,  unlike  the  low-order  examples  presented  in
most of the aforementioned publications. We also investigate the grid-aligned shock instability for the HLLC solver in three
dimensions and reveal that carbuncles are more likely to occur than in two dimensions.

The  paper  is  organized  as  follows.  In  Section 2,  the  governing  equations  and  the  general  framework  of  Godunov-type
methods  are  reviewed  together  with  the  classical  HLLC  ﬂux  formulation.  A  centralized  formulation  of  the  HLLC  ﬂux  is
derived in the ﬁrst part of Section 3, followed by the low Mach number adapted wave speed formulations resulting in the
newly  proposed  HLLC-LM  scheme.  In  Section 4,  a  comprehensive  set  of  test  cases  is  studied  to  verify  the  accuracy  and
shock  stability  of  the  new  scheme.  Results  are  also  provided  with  high  resolution  including  a  study  of  three-dimensional
effects. Finally in Section 5, calculations of complex ﬂow phenomena that take advantage of the applied high-order schemes,
such  as  a  ﬂow  around  a  diamond  and  multi-component  ﬂows  with  nontrivial  shock-interface  interactions,  are  studied  to
further demonstrate both the stability and the reduced numerical dissipation of the HLLC-LM ﬂux. Conclusions are drawn
in Section 6.

2

N. Fleischmann, S. Adami and N.A. Adams

Journal of Computational Physics 423 (2020) 109762

2.  Governing equations and numerical approach

We consider an inviscid compressible ﬂow that evolves according to the three-dimensional Euler equations

Ut + F (U)x

+ G (U) y

+ H (U)z

= 0,

(1)

where U is the density of the conserved quantities mass ρ, momentum ρv ≡ (ρu, ρ v, ρ w) and total energy E = ρe + 1
with e being the internal energy per unit mass. The ﬂuxes F , G and H are deﬁned as

2 ρv2,

⎛

F =

⎜
⎜
⎜
⎜
⎝

ρu
ρu2 + p
ρuv
ρu w
u (E + p)

⎞

⎛

⎟
⎟
⎟
⎟
⎠

, G =

⎜
⎜
⎜
⎜
⎝

ρ v
ρuv
ρ v 2 + p
ρ v w
v (E + p)

⎞

⎛

⎟
⎟
⎟
⎟
⎠

, H =

⎜
⎜
⎜
⎜
⎝

ρ w
ρu w
ρ v w
ρ w 2 + p
w (E + p)

⎞

⎟
⎟
⎟
⎟
⎠

.

(2)

The set of equations is closed by the ideal-gas equation of state, where the pressure  p is given by  p = (γ − 1) ρe with a
constant ratio of speciﬁc heats γ .

2.1.  Finite volume approach

Our  numerical  framework  is  identical  to  the  one  described  in  [8],  where  Godunov’s  approach  [1] for  ﬁnite  volumes  is
applied to solve the given set of equations. The time evolution of the vector of cell-averaged conservative states  ¯U is given
by

d

dt

¯Ui = 1
(cid:4)x

(Fi− 1

2 , j,k

− Fi+ 1

2 , j,k

+ Gi, j− 1
2 ,k

− Gi, j+ 1
2 ,k

+ Hi, j,k− 1

2

− Hi, j,k+ 1

2

),

(3)

where  (cid:4)x is  the  cell  size  of  a  uniform  Cartesian  grid  and  F,  G and  H approximate  the  cell-face  ﬂuxes  in  x-,  y- and
z-direction,  respectively.  These  ﬂuxes  are  determined  dimension-by-dimension  from  a  Riemann  solver  combined  with  a
high-order  WENO  spatial  reconstruction  scheme  [26].  Additional  volume  source  terms,  such  as  gravitational  acceleration,
are  omitted  here  for  simplicity.  The  resulting  system  of  ODE  (3) is  integrated  in  time  using  a  high-order  strong  stability-
preserving (SSP) Runge-Kutta scheme [27].

2.2.  The HLLC Riemann solver

In  order  to  avoid  computationally  expensive  iterative  solution  of  the  Riemann  problem,  approximate  Riemann  solvers
are  commonly  employed.  In  this  paper,  we  focus  on  one  speciﬁc  approximation,  the  HLLC  solver,  which  is  one  of  the
most  popular  and  versatile  Riemann  solvers.  It  has  been  extended  to  a  broad  range  of  applications,  also  beyond  classical
computational ﬂuid dynamics [11].

Toro et al. [4] deﬁne the HLLC ﬂux as

FH LLC =

⎧
⎪⎪⎪⎨
FL
F∗L = FL + S L · (U∗L − UL)
⎪⎪⎪⎩
F∗R = FR + S R · (U∗R − UR )
FR

if S L ≥ 0,
if S L < 0 ∩ S∗ ≥ 0,
if S R > 0 ∩ S∗ ≤ 0,
if S R ≤ 0,

where two intermediate states, U∗L and U∗R , are separated by the contact wave and are determined from

⎛

⎜
⎜
⎜
⎜
⎝

U ∗K = S K − u K
S K − S∗

ρK
ρK S∗
ρK v K
ρK w K
(cid:12)

⎞

⎟
⎟
⎟
⎟
⎠

(cid:13)

E K + (S∗ − u K )

ρK S∗ + p K

S K −u K

with  K = L, R, and UL , UR being the reconstructed left and right face states, respectively.

Following Einfeldt [9], the maximum left and right nonlinear signal speed estimates are obtained from

S L = min(u L − cL, ˆu − ˆc),

S R = max(u R + c R , ˆu + ˆc),

where  ˆu and ˆc are determined from the Roe average

ˆu = u L ·

√
√

√

ρR

ρL + u R ·
ρL +
ρR

√

3

(4)

(5)

(6)

(7)

N. Fleischmann, S. Adami and N.A. Adams

Journal of Computational Physics 423 (2020) 109762

Fig. 1. Schematic illustration of the vanishing velocity component v in transverse direction of the shock front propagation.

and

2 ·

ˆc2 = cL

√
√

ρL + c2
·
R
√
ρL +
ρR

√

ρR

+ 1
2

√

√

ρL
ρL +

ρR
√
ρR

(cid:14)√

(u R − u L)2 .

(cid:15)

2

The contact wave speed is estimated according to Batten et al. [12] from

S∗ = p R − p L + ρL u L (S L − u L) − ρR u R (S R − u R )

ρL (S L − u L) − ρR (S R − u R )

.

(8)

(9)

High-order approximations for the left and right face states, UL and UR , are obtained upon characteristic decomposition in
combination with a high-order WENO scheme as described in detail in [28].

3.  A shock-stable HLLC type solver with low Mach number modiﬁcation

An inaccurate scaling behavior of the acoustic and advection contribution to the numerical dissipation in the low Mach
number limit has been found to be the driving mechanism of the numerical grid-aligned shock instability [8]. The connec-
tion is motivated by the observation, that shock instabilities only occur when a high Mach number shock wave propagates
almost perfectly aligned with the computational grid.

When the shock wave moves in x-direction as shown in Fig. 1, the velocity components of the local transverse direction
v,  respectively  v and  w for  the  three-dimensional  case,  have  a  vanishing  magnitude.  Consequently,  the  local  directional
Mach number will also vanish during the computation of the ﬂuxes in transverse directions of the shock wave propagation.
Note that a perfect alignment with zero Mach number in transverse directions leads to a one-dimensional situation where
no instability occurs. A small deﬂection is always required to trigger the instability. There is a thorough documentation of the
shortcomings of Riemann solvers in the low Mach regime [29–31] which dates back to the ﬁndings of Guillard et al. [32,33].
In [32], the authors showed that a wrong scaling behavior of the numerical dissipation leads to pressure ﬂuctuations that
may ruin the prediction of low Mach number ﬂows using Godunov’s approach. This ﬂaw is now considered as the driving
mechanism of the grid-aligned shock instability. In their recent publication, Chen et at. [34] performed a stability analysis
to investigate the shock instability mechanism for simpliﬁed systems. Their results support the given argumentation as the
authors  also  detect  an  inaccurate  pressure  dissipation  of  the  Riemann  solver  at  the  vertical  transverse  face  of  the  shock
to  be  the  driving  mechanism  for  the  instability.  A  minor  modiﬁcation  that  reduces  the  acoustic  dissipation  of  the  Roe
Riemann  solver  in  the  low  Mach  number  limit  proved  to  be  effective  in  suppressing  the  instability  [8].  The  reduction  of
acoustic dissipation can be achieved by reduction of the nonlinear eigenvalues of the Roe dissipation matrix for small Mach
numbers. This procedure stabilizes simulations of supersonic ﬂows. For dealing with global low Mach number ﬂows near
the incompressible limit, there are other methods available in literature [29–31].

A straightforward modiﬁcation of the nonlinear signal speeds of the HLLC solver following [8] turns out to be ineffective
in suppressing grid-aligned shock instabilities. The reason for the ineffectiveness can be found when the limit solution of
the Roe-M ﬂux and a modiﬁed HLLC ﬂux with similarly reduced nonlinear signal speed are compared for vanishing Mach
numbers. While the Roe-M approximation [8] reduces in the low Mach number limit to the central ﬂux term

FRoe−M Ma→0−−−−→ 1
2

(FL + FR) ,

a modiﬁed HLLC approximation with identical reduction of nonlinear signal speeds results in

(10)

4

N. Fleischmann, S. Adami and N.A. Adams

Journal of Computational Physics 423 (2020) 109762

FH LLC−R E DU C E D Ma→0−−−−→

(cid:16)

FL
FR

if S∗ ≥ 0,
if S∗ ≤ 0.

(11)

Thus, differently from the limit solution of the Roe-M scheme, a straightforward modiﬁcation of the HLLC ﬂux leads to a
pure classical upwind scheme. Upwinding is not required in the absence of shocks and, moreover, introduces an undesir-
able  amount  of  numerical  dissipation,  which  counteracts  the  objective  of  reducing  dissipation.  Thus,  the  goal  is  to  ﬁnd  a
formulation of the HLLC ﬂux that continuously approaches the central ﬂux term in the limit of low Mach numbers.

3.1.  Central formulation of the HLLC ﬂux

In a ﬁrst step, the classical HLLC ﬂux will be reformulated motivated by the derivation of the central Roe ﬂux formulation.
The intermediate ﬂux  F ∗L can be determined using two alternative approaches

F∗L = FL + S L (U∗L − UL)

and

F∗L = FR + S R (U∗R − UR ) + S∗ (U∗L − U∗R ) .

(12)

(13)

While  the  traditional  derivation  of  Eq.  (12) applies  the  Rankine-Hugoniot  condition  only  once  starting  from  the  left  side,
alternatively,  the  Rankine-Hugoniot  condition  can  also  be  applied  twice  starting  from  the  right  side,  Eq.  (13).  A  central
formulation of F∗L can be established by averaging both formulations and is given by

F∗L = 1
2

(FL + FR ) + 1
2

[S L (U∗L − UL) + S∗ (U∗L − U∗R ) + S R (U∗R − UR )] .

Analogously, the right intermediate ﬂux can be determined by

F∗R = FR + S R (U∗R − UR )

and

F∗R = FL + S L (U∗L − UL) + S∗ (U∗R − U∗L)

resulting in

F∗R = 1
2

(FL + FR ) + 1
2

[S L (U∗L − UL) − S∗ (U∗L − U∗R ) + S R (U∗R − UR )] .

(14)

(15)

(16)

(17)

By  comparing  Eq.  (14) and  Eq.  (17) we  note  that  only  the  sign  of  the  third  term,  which  is  related  to  the  contact  wave,
differs for both expressions. Finally, considering the requirement that F∗L is applied if  S∗ ≥ 0 and F∗R is applied if  S∗ ≤ 0, a
central formulation of the HLLC ﬂux is obtained by

FH LLC =

⎧
⎪⎨

⎪⎩

FL
FR
F∗

if S L ≥ 0,
if S R ≤ 0,
else

with

F∗ = 1
2

(FL + FR )

1

2

[S L (U∗L − UL) + |S∗| (U∗L − U∗R ) + S R (U∗R − UR )] .

3.2.  On the numerical dissipation of HLL(C)-type solvers

(18)

(19)

Using the centralized formulation derived in Section 3.1 both the HLL and the HLLC ﬂux in the subsonic regime can be

written as

FH LL = 1
2
FH LLC = 1
2

(FL + FR ) − 1
2
(FL + FR ) − 1
2

[|S L| (U∗ − UL) + |S R | (UR − U∗)]

[|S L| (U∗L − UL) + |S∗| (U∗R − U∗L) + |S R | (UR − U∗R )] .

(20)

A connection to the Lax-Friedrichs ﬂux can be established, when |S L| = |S∗| = |S R | = |λ| is introduced into FH LL or FH LLC
resulting in

Fλ
H LL(C)

= 1
2

(FL + FR ) − 1
2

|λ| (UR − UL) .

5

(21)

N. Fleischmann, S. Adami and N.A. Adams

Journal of Computational Physics 423 (2020) 109762

Fig. 2. Dependence of the activation function φ on the local Mach number Malocal with Malimit = 0.1.

Now,  the  HLL(C)  ﬂux  can  be  seen  as  a  Lax-Friedrichs  ﬂux  where  the  dissipation  has  been  split  into  two  (HLL),  or  three
(HLLC) differently weighted contributions representing the general wave system of the underlying Riemann problem.

In the original formulation of the HLLC approximation (4), advection and acoustic contributions to the numerical dissipa-
tion are diﬃcult to separate. However, the proposed central formulation of the HLLC solver allows for a separation of both
contributions in analogy with the Roe ﬂux, which is given by

FRoe = 1
2

(FL + FR ) − 1
2

R |(cid:2)| R

−1 (UR − UL) .

(22)

The  ﬁrst  part  both  in  Eq.  (19) and  Eq.  (22) is  the  central  ﬂux  term,  and  the  second  term  is  the  dissipation  ﬂux  term,
which is characteristic for each solver. The advection dissipation of the Roe ﬂux is proportional to the eigenvalue |u|, and
the acoustic dissipation of the Roe ﬂux is proportional to the eigenvalues |u ± c|. Analogously, the acoustic dissipation of
the  HLLC  ﬂux  is  related  to  the  ﬁrst  and  third  term  of  the  dissipation  ﬂux  term  as  both  terms  are  proportional  to  the
acoustic signal speed  S L , respectively  S R . The advection dissipation is related to the center term, which is proportional to
the contact signal speed  S∗. Note that the situation for the HLLC ﬂux is more complex than for the Roe ﬂux since  S L and
S R also contribute to  S∗, U∗L and U∗R . However, the results of this paper indicate that the main contributions of advection
and acoustic dissipation can be distinguished as discussed.

3.3.  HLLC-LM ﬂux with low Mach number correction

The main goal of the proposed modiﬁcation is to balance the vanishing advective and dominant acoustic dissipation in
the  low  Mach  number  limit  by  a  reduction  of  overall  dissipation.  The  central  formulation  of  the  HLLC  ﬂux  given  by  Eq.
(18) and Eq. (19) enables a straightforward application of the Mach number dependent reduction of nonlinear signal speeds
according to

S H LLC−LM

L

= φ · S L,

S H LLC−LM

R

= φ · S R

with

and

(cid:17)

(cid:17)

φ = sin

min

1,

(cid:18)

(cid:18)

· π
2

Malocal
Malimit

Malocal = max

(cid:17)(cid:19)
(cid:19)
(cid:19)
(cid:19)

u L
cL

(cid:19)
(cid:19)
(cid:19)
(cid:19) ,

(cid:19)
(cid:19)
(cid:19)
(cid:19)

u R
c R

(cid:19)
(cid:18)
(cid:19)
(cid:19)
(cid:19)

.

(23)

(24)

(25)

and  S H LLC−LM
u denotes the velocity component dependent on the direction of the cell-face Riemann problem.  S H LLC−LM
are only applied for the ﬁnal ﬂux evaluation in Eq. (19). All previous procedures, especially the calculation of  S∗, U∗L and
U∗R , are performed using the original values for  S L and  S R .

The application of the sine function in Eq. (24) causes a smooth decay of the acoustic dissipation as depicted in Fig. 2.
The reference parameter Malimit is set to 0.1 for all calculations presented in this paper. This selection ensures that the mod-
iﬁcation will only be active if the local ﬂow speed component is less than ten percent of the local sound speed. Otherwise,
the classical HLLC formulation is fully recovered. The new scheme, denoted as HLLC-LM in the following, fully preserves the
favorable  low  dissipation  of  HLLC  at  the  contact  line  as  the  acoustic  dissipation  of  HLLC-LM  is  reduced  proportionally  to
the level of local velocities instead of the speed of sound for low Mach numbers while the advection dissipation remains
unchanged.

R

L

6

N. Fleischmann, S. Adami and N.A. Adams

Journal of Computational Physics 423 (2020) 109762

Fig. 3. Comparison of the classical HLLC formulation with the central HLLC formulation for the corner diffraction of a Mach 5.09 shock wave: logarithmic
gradients of density from 1 to 1, 000 at t = 0.157.

4.  Central aspects of the grid-aligned shock instability with HLLC-type solvers demonstrated for classical test cases

The calculations in this section serve to study the evolution of the numerical shock instability when HLLC-type solvers
are  applied  in  combination  with  high-order  schemes.  Moreover,  the  stability  of  the  HLLC-LM  scheme  with  respect  to  the
grid-aligned shock problem is demonstrated for a comprehensive set of cases with strong moving shocks that are prone to
exhibiting this instability. If not mentioned otherwise, all calculations were performed using the classical ﬁfth-order WENO
scheme [26] for spatial discretization combined with a third-order strong-stability-preserving Runge-Kutta time integration
[27] and the approximate Riemann solvers as described in the previous sections. The effective range of the shock-transverse
Mach  number  modiﬁcation  in  the  HLLC-LM  solver  is  always  limited  to  local  Mach  numbers  lower  than  0.1.  The  ﬂuid  is
modeled as  ideal  gas  with  γ = 1.4.  The  CFL  number  is  set  0.6 for  single-phase  cases  and  0.4 for  cases  with  interfaces
employing the level-set approach. The combination of a multiresolution procedure [35] and an adaptive local time stepping
[36] enables eﬃcient computation with high effective resolutions. In the following, the given resolution information deﬁnes
the  ﬁnest  level.  Shocks  are  discretized  with  the  highest  resolution  in  all  presented  cases  due  to  the  applied  reﬁnement
criteria, whereas material interfaces are by deﬁnition on the highest level.

4.1.  Corner ﬂow problem I: veriﬁcation of centralized HLLC formulation

As a ﬁrst step, the proposed centralized HLLC formulation given in Eq. (18) and Eq. (19) is veriﬁed against the classical
HLLC procedure for the diffraction of a shock wave around a sharp corner. This is a well-established test case, where the
instability  of  the  HLLC  ﬂux  becomes  apparent.  This  case  was  already  selected  by  Quirk  [6] to  demonstrate  the  failure  of
low-dissipation Riemann approximations. Additionally, the problem yields complex ﬂow patterns. Thus, it is well suited to
compare results of different solvers and to verify our reformulations.

We use a domain of size [0, 1] × [0, 1], that is uniformly initialized with (ρ, u, v, p) = (1, 0, 0, 1/1.4) and discretized by
1280 × 1280 cells. Reﬂecting-wall boundary conditions are set everywhere, except for the upper left boundary at x = 0 from
y = 0.5 to  y = 1. Here, the post-shock condition of a Mach 5.09 shock wave is prescribed. The ﬁnal time is set to 0.8/Ma.
Even  though  the  ﬁrst-order  HLLC  approximation  is  known  to  be  positivity  preserving,  this  property  is  not  guaranteed  for
high-order extensions [37]. We encountered instabilities in the vicinity of the corner point of the backward facing step at
the inﬂow for all tested variants of the HLLC ﬂux when combined with a ﬁfth-order WENO scheme. Therefore, simulations
were performed using a third-order WENO scheme [26].

Fig. 3 and Fig. 4 show schlieren images of the density gradients at the ﬁnal time of the simulation. The results shown
in the left frame of Fig. 3 are obtained applying the original HLLC formulation given in Eq. (4), whereas results shown in
the right frame of Fig. 3 are obtained applying the centralized HLLC formulation given in Eq. (18) and Eq. (19). As expected,
there are no distinguishable differences for both formulations. Moreover, all other test cases presented in this paper have
been investigated without encountering any differences exceeding the ﬂoating-point roundoff error. We therefore conclude
that Eq. (18) with Eq. (19) is a valid alternative representation of the HLLC ﬂux.

4.2.  Corner ﬂow problem II: stability of HLLC-LM formulation

As  a  second  step,  the  stability  of  the  HLLC-LM  scheme  is  demonstrated.  The  aforementioned  corner  ﬂow  simulations
show  severe  disturbances  in  the  backﬂow  of  the  leading  shock  front  similar  to  results  obtained  with  the  Roe  Riemann

7

N. Fleischmann, S. Adami and N.A. Adams

Journal of Computational Physics 423 (2020) 109762

Fig. 4. Corner diffraction of a Mach 5.09 shock wave: logarithmic gradients of density from 1 to 1, 000 at t = 0.157.

solver [6]. The left frame of Fig. 4 shows that the HLLC-LM is able to capture all details of the ﬂow while preventing any
disturbances of the shock wave. Additionally, a high-resolved simulation is performed using 16 times smaller cells. Typically,
the instability is enhanced by higher resolutions, however, the results presented in the right frame of Fig. 4 are still free of
any instability.

4.3.  Rayleigh Taylor instability: numerical dissipation at contact lines

The inherent numerical dissipation of the original HLLC ﬂux and HLLC-LM ﬂux is compared by investigating a classical
Rayleigh-Taylor instability. Two initial gas layers with different densities are exposed to gravity with unity magnitude, where
the resulting acceleration is directed towards the lighter ﬂuid. A small disturbance of the contact line triggers the instabil-
ity.  The  computational  domain  is  given  by  [0, 0.25] × [0, 1] and  the  interface  initially  is  placed  at  y = 0.5.  Initial  states
= (1, 0, −0.025c · cos(8π x), y + 1.5),
are given by (ρ, u, v, p) y≤0.5
γ p
where the speed of sound is c =
ρ with γ = 5
3 . Top and bottom boundary states are ﬁxed to (1, 0, 0, 2.5) and (2, 0, 0, 1),
respectively. Symmetry boundary conditions are imposed at the left and right boundary.

= (2, 0, −0.025c · cos(8π x), 2 y + 1) and (ρ, u, v, p) y>0.5

The ﬁnal density evolution for both solvers is shown in Fig. 5 for a resolution of 128 × 512. Results indicate a signiﬁcant

(cid:20)

reduction of dissipation at the contact line when the HLLC-LM ﬂux is applied instead of the original HLLC ﬂux.

4.4.  Quirk’s odd-even decoupling test: quantitative evaluation of the shock instability

The results of the Section 4.2 indicate the effectiveness of the proposed method qualitatively, however, a detailed quanti-
tative study is diﬃcult to perform for the corner ﬂow problem. For this purpose, the simple plane shock propagation along
a  rectangular  duct  with  a  deﬁned  disturbance  level  is  studied.  This  test  case  was  also  proposed  by  Quirk  [6] due  to  its
simple setup. Nevertheless, it provides an effective and reliable way to trigger the odd-even decoupling near strong shocks,
which is related to the grid-aligned shock instability. Moreover, it allows for a simple quantitative study of the rise of the
instability.

The domain is set to [0, 2400] × [0, 20], and discretized with 2400 × 20 cells. Inﬂow and outﬂow conditions are applied
at the left and at the right boundary, respectively. Reﬂecting wall conditions, which are equivalent to symmetry boundary
conditions for inviscid ﬂows, are enforced both at the top and at the bottom boundary of the domain. Pre-shock density and
pressure are set to unity, and all velocity components are set to zero. Artiﬁcial numerical noise is introduced to all primitive
variables in the initial state to trigger the instability [38,8]. We have performed simulations with the original Mach 6 setup
and with a Mach 20 setup with initial conditions given by

(ρ, u, v, p) =

⎧
⎪⎪⎨
⎪⎪⎩

(1, 0, 0, 1)
(cid:12)
√
35
36 , 0, 251
216
41 , 35
√
27 , 133
160

(cid:12)

6

8

1.4, 0, 466.5

(cid:13)

(cid:13)

if x > 5,
else (for Ma = 6 case),
else (for Ma = 20 case),

(26)

where the shock front is initially placed at x = 5. Both simulations are performed up to a late point in time till the shock
front approaches the end of the domain. The ﬁnal time is set to 330 for the low Mach number simulation and to 100 for
the high Mach number simulation, respectively.

8

N. Fleischmann, S. Adami and N.A. Adams

Journal of Computational Physics 423 (2020) 109762

Fig. 5. Rayleigh-Taylor instability t = 1.95: density contours from 0.85 (blue) to 2.25 (red). (For interpretation of the colors in the ﬁgure(s), the reader is
referred to the web version of this article.)

Fig. 6. Instability progress in Quirk’s test case for Mach 6.

The maximum magnitude of the y-velocity component  v in the domain provides a reasonable measure of the deviation
from  the  one-dimensional  solution,  and  therefore,  it  is well  suited  to  monitor  the  growth  rate  of  the  disturbance  quanti-
tatively over time. Fig. 6 and Fig. 8 show the evolution of the velocity deviation for the Mach 6 and the Mach 20 case for
different ﬂux approximations when all initial primitive variables are superposed by uniform random perturbations ranging
−3.  In  addition,  the  ﬁnal  density  distributions  are  presented  in  Fig. 7 and  Fig. 9.  Besides  the
from  −0.5 · 10
results for the discussed HLLC and HLLC-LM ﬂuxes, the results for the more dissipative HLL ﬂux [3] and the shock-stable
Roe-M ﬂux [8] are provided for comparison.

−3 to  0.5 · 10

Simulations with the classical HLLC solver show an exponential instability where instabilities saturate at O(1) at around
t = 20 for the low Mach number case and at around t = 5 for the high Mach number case. Unlike the Roe approximation,

9

N. Fleischmann, S. Adami and N.A. Adams

Journal of Computational Physics 423 (2020) 109762

Fig. 7. Quirk’s test case for Mach 6: color map of density from blue = 1.0 to red = 6.8 at t = 330.

Fig. 8. Instability progress in Quirk’s test case for Mach 20.

Fig. 9. Quirk’s test case for Mach 20: color map of density from blue = 1.0 to red = 8.0 at t = 100.

the  HLLC  ﬂux  forms  no  distinct  carbuncles,  and  density  disturbances  remain  bounded.  However,  the  instability  disturbs
the  shock  front  signiﬁcantly  as  shown  in  the  left  frame,  respectively  top  frame,  of  Fig. 7 and  Fig. 9.  Moreover,  when  the
ﬁnal  position  of  the  shock  front  is  compared  to  the  analytically  predicted  position,  an  incorrect  wave  speed  is  obtained.
This  effect  is  even  more  dominant  for  the  high  Mach  number  case.  With  the  modiﬁed  HLLC-LM  scheme,  the  stable  and
analytically predicted result is obtained as depicted in the middle frames of Fig. 7 and Fig. 9. The magnitude of disturbances

10

N. Fleischmann, S. Adami and N.A. Adams

Journal of Computational Physics 423 (2020) 109762

Fig. 10. Double Mach reﬂection of a Mach 10 shock wave: 40 density contours from 1.88783 to 20.9144 at t = 0.2.

is similar to the one obtained with the Roe-M formulation [8] and slightly higher than the one obtained with the HLL ﬂux
for both cases. The lower disturbances of the HLL ﬂux can be explained by a signiﬁcantly higher level of inherent dissipation
of the scheme. However, no major differences can be observed in the qualitative density results for HLL and HLLC-LM, e.g.
middle and bottom frame of Fig. 9.

4.5.  Double Mach reﬂection problem: effect of resolution

Several numerical schemes encounter diﬃculties when simulating a double Mach reﬂection as proposed by Woodward
and  Colella  [39].  The  leading  Mach  stem  may  be  kinked  in  consequence  of  the  numerical  shock  instability  [6,7].  The  test
case represents a Mach 10 shock wave hitting a solid ramp with an angle of 30 degrees. The initial shock wave is set up
with

(cid:16)

(ρ, u, v, p) =

(1.4, 0, 0, 1)
(cid:12)

√
8 , −4.125, 116.5

3

8, 33

√

3 (x − 1/6) ,

(cid:13)

if y <

else.

(27)

A  Neumann  boundary  condition  with  zero  gradients  for  all  variables  is  applied  at  the  left,  right  and  upper  boundary.
Along the bottom boundary, at  y = 0, the region from x = 0 to x = 1/6 is always assigned post-shock conditions, whereas
reﬂecting-wall conditions are imposed from  x = 1/6 to  x = 4. The domain size of [0, 4] × [0, 6.67] is chosen large enough
to avoid any disturbances entering the domain at the upper boundary. The domain is discretized with 960 × 1600 cells and
the ﬁnal time is set to t = 0.2. Besides the large vertical domain size, this setup is commonly chosen in literature [40].

The ﬁnal density contours for both HLLC and HLLC-LM are shown in Fig. 10. Both schemes deliver almost identical results
with no visible deﬂection at the leading Mach stem. However, if the resolution is increased to  1920 × 3200 cells and the
ﬁnal  time  is  set  to  t = 0.28 the  results  for  both  schemes  differ  signiﬁcantly  as  shown  in  Fig. 11.  A  kinked  Mach  stem,
together  with  a  severe  disturbance  of  the  wall  jet  can  be  observed  for  the  original  HLLC  scheme,  whereas  the  HLLC-LM
scheme is free of any instability.

11

N. Fleischmann, S. Adami and N.A. Adams

Journal of Computational Physics 423 (2020) 109762

Fig. 11. Double Mach reﬂection of a Mach 10 shock wave: 40 density contours from 1.88783 to 20.9144 at t = 0.28.

4.6.  Supersonic ﬂow around cylinder: steady shock position

The next case predicts the bow shock resulting from a supersonic ﬂow around a stationary cylinder. This case was ﬁrst
described by Peery and Imlay [5] to suffer from the carbuncle phenomenon. Unlike the other cases in this paper, the relevant
shock  wave  is  not  moving,  which  renders  the  case  particularly  challenging  for  high-order  shock-capturing  schemes  with
explicit time integration. Following the argumentation in [8], we do not change the Cartesian grid nor the time integration,
which likely results in a small resolved level of ﬂuctuations around the steady shock due to the high order of the applied
scheme. We include this case for the sake of completeness even though the application of high-order schemes here is not
expected to reveal additional information for such conﬁgurations compared to low-order schemes.

The circular reﬂecting-wall condition representing the cylinder is approximated using a level-set approach [41]. At the
left  and  the  remaining  right  boundary  inﬂow  and  outﬂow  conditions  are  applied,  respectively.  Top  and  bottom  boundary
conditions are set to Neumann boundary conditions with zero gradient for all variables. Two different Mach numbers, Ma =
(cid:12)
3 and Ma = 20, are studied with initial states (ρ, u, v, p) =
. The domain size is set to [0, 0.3] × [0, 0.8]
1,
for  the  lower  Mach  number,  and  [0, 0.3] × [0, 0.6] for  the  higher  Mach  number.  Final  times  are  chosen  large  enough  to
reach a fully developed bow shock. The cylinder with a diameter  D = 0.2 is placed at the center of the right boundary and
resolved by 160 cells per diameter.

(cid:13)
1.4 · Ma, 0, 1

√

Besides the HLLC and the HLLC-LM schemes, the more dissipative HLL scheme is also applied. Fig. 12 and Fig. 13 show
the resulting pressure distributions and Mach contour lines that are chosen identical to [7] for both Mach number ﬂows. All
three schemes show comparable results for both Mach numbers. Note that also the HLL scheme reveals some disturbances
in the backﬂow of the steady shock. These disturbances of the HLL scheme in combination with high-order methods have
been reported in literature [42]. None of the schemes suffers from the carbuncle phenomenon with the described Cartesian
setup.  Moreover,  the  HLLC-LM  scheme  has  been  tested  for  a  signiﬁcantly  increased  resolution  of  640 cells  per  diameter,
where it still delivers stable results as shown in the right frames of Fig. 12 and Fig. 13.

4.7.  The Sedov blast wave: comparison of shock instability in two and three dimensions

The  next  case  of  this  section  is  the  classical  Sedov  blast  wave  [43,37,7].  Due  to  its  symmetry,  the  Sedov  blast  wave
simulation is suitable to demonstrate the effect of the grid alignment on the numerical shock instability [8]. The test case
consists of a high pressure area covering only few cells that is initiated at the center of the domain. The rest of the domain
is set to a near vacuum state. The whole domain is initially at rest. The initial states are given by

(ρ, u, v, p) =

(cid:16)(cid:14)
(cid:14)

(cid:15)

,

1, 0, 0, 3.5 · 105
(cid:15)
1, 0, 0, 10
,

−10

(cid:21)

x2 + y2 < 0.005,

if
otherwise.

12

(28)

N. Fleischmann, S. Adami and N.A. Adams

Journal of Computational Physics 423 (2020) 109762

Fig. 12. Supersonic ﬂow around cylinder Ma = 3 at t = 1.5: color pressure map (blue = 1.0 to red = 12.1) is overlaid by 25 Mach contours (0.1 to 2.5).

Fig. 13. Supersonic ﬂow around cylinder Ma = 20 at t = 0.5: color pressure map (blue = 1.0 to red = 550) is overlaid by 25 Mach contours (0.1 to 2.5).

Reﬂecting-wall conditions are applied at all boundaries. The domain size is set to [−1.2, 1.2] × [−1.2, 1.2], and it is resolved
by 960 × 960 cells. The ﬁnal time is set to 0.1.

The schlieren image for logarithmic density gradients is given in Fig. 14 when using the HLLC ﬂux and the HLLC-LM ﬂux.
At  locations  where  the  shock  front  propagates  aligned  with  the  computational  grid,  disturbances  behind  the  shock  wave
can be observed. The magnitude of disturbances is smaller than for the Roe ﬂux [8] and no carbuncles occur. The results
obtained with the HLLC-LM ﬂux are free of any disturbance.

Finally, we extend the problem to three dimensions in a straightforward way. The domain size is set to [−2, 2] ×[−2, 2] ×
[−2, 2] and is resolved by a resolution of 640 × 640 × 640 cells. In order to save computational cost, only one-eighth of the

13

N. Fleischmann, S. Adami and N.A. Adams

Journal of Computational Physics 423 (2020) 109762

Fig. 14. Two-dimensional Sedov blast wave: logarithmic gradients of density from 1 to 500 at t = 0.1.

Fig. 15. Three-dimensional Sedov blast wave: logarithmic gradients of density from 1 to 500 at t = 0.1.

given  domain  is  simulated  and  appropriate  symmetry  boundary  conditions  are  applied.  The  setup  is  chosen  according  to
Tasker et al. [43], where initial states are

(ρ, u, v, w, p) =

(cid:16)(cid:14)
(cid:14)

1, 0, 0, 0, 23.757239 · 106
1, 0, 0, 0, 1 · 10

−10

(cid:15)

,

(cid:21)

(cid:15)

,

if

x2 + y2 + z2 < 0.0875,

otherwise

(29)

with γ = 5

3 . The ﬁnal time again is set to 0.1.

Results for the three-dimensional Sedov blast wave are shown in Fig. 15 for both HLLC and HLLC-LM. Differently from
the two dimensional case, the three-dimensional simulation reveals an increased level of disturbances and the occurrence
of signiﬁcant carbuncles for the HLLC ﬂux. This indicates that the instability is enhanced for three-dimensional simulations.
Following the argumentation of Section 3, this behavior can be explained as follows. A three-dimensional shock wave that
propagates along one coordinate axis suffers from an excessive acoustic dissipation that is now introduced from two sides
as the ﬂuxes in both other directions have a vanishing Mach number. As expected, the reduction of acoustic dissipation in
the HLLC-LM scheme also helps to prevent the grid-aligned shock instability in three-dimensional simulations.

4.8.  Subsonic ﬂow around cylinder: low Mach number ﬂow

In addition to the shock-dominated ﬂow problems presented before, the performance of the proposed HLLC-LM scheme
also is tested in the global low Mach number regime using the well-known test case of a subsonic ﬂow around a cylinder.

14

N. Fleischmann, S. Adami and N.A. Adams

Journal of Computational Physics 423 (2020) 109762

Fig. 16. Flow around a cylinder at Ma = 0.01 using ﬁrst order (top) and WENO5 (bottom): color density map (blue = 0.99993 to red = 1.00007) is overlaid
by 21 contour lines for normalized pressure ﬂuctuations from −7 · 10

−5 to 7 · 10

−5.

This ﬂow conﬁguration is troublesome for Godunov schemes in combination with Riemann solvers as comprehensively dis-
cussed in literature, e.g. [29]. Different modiﬁcations to Riemann solvers and preconditioning techniques have been proposed
to increase the simulation accuracy of low Mach number ﬂows [32,30,31].

The domain of size [0, 80D] × [0, 80D] is set large enough to avoid any interaction of reﬂected waves, which is crucial
for the high-order simulation. The cylinder is placed in the center of the domain with a diameter  D = 1. Initial density and
pressure are set to unity in the entire domain. The initial velocity of u = 0.01 ·
1.4 results in a free-stream Mach number of
0.01. At all boundaries we apply Neumann boundary conditions with zero gradient for all variables. The effective resolution
is set to 128 cells per diameter. The ﬁnal time t = 30 is large enough to approach a steady state before disturbances due
to reﬂections at the domain boundaries affect the region of interest around the cylinder. Note that the application of high-
order  schemes  in  combination  with  explicit  time  integration  for  the  fully  compressible  evolution  equations  renders  low
Mach number simulations particularly expensive.

√

Fig. 16 shows the density distribution in the relevant region around the cylinder and 21 isocontours for pressure ﬂuctu-
−5 similarly to [29] for HLLC and HLLC-LM using both a ﬁrst-order and a

ations δ p = p − p0 between ±γ Ma2/2 = ±7 · 10
WENO5-JS spatial discretization.

The fully symmetric ﬂow ﬁeld obtained with WENO5 shows excellent agreement with the expected result. In either case,
the  HLLC-LM  solver  shows  similar  or  better  performance  than  the  original  HLLC.  Nevertheless,  it  should  be  pointed  out
that  the  HLLC-LM  solver  is  primarily  designed  for  applications  in  the  high  Mach  number  regime  that  suffer  from  shock
instabilities.  Due  to  the  decreasing  numerical  dissipation  in  the  low  Mach  number  limit  we  expect  the  occurrence  of
pressure-velocity decoupling when the Mach number is further reduced.

15

N. Fleischmann, S. Adami and N.A. Adams

Journal of Computational Physics 423 (2020) 109762

Fig. 17. Supersonic ﬂow around a diamond-shaped obstacle with Ma = 2.85: logarithmic gradients of density from 1 to 1000 at t = 0.5.

5.  Application to complex ﬂow situations

The main motivation for the application of high-order low-dissipation schemes is an accurate prediction of highly com-
plex ﬂow situations. Therefore, we studied three additional types of test cases that involve interaction of shock waves with
nontrivial structures and recent examples of multi-component ﬂow simulations using the level-set approach [41].

5.1.  Supersonic ﬂow around diamond-shaped obstacle

The  ﬁrst  example  of  a  highly  complex  ﬂow  evolution  is  the  supersonic  ﬂow  around  a  diamond-shaped  obstacle.  The
Mach  number  of  2.85  is  chosen  to  be  high  enough  to  form  a  double  Mach  reﬂection  during  and  after  the  shock  wave
propagates over the diamond [44]. The sharp geometry changes result in extremely complex ﬂow patterns in the wake of
the diamond. In addition to the double Mach reﬂection, this case also involves a bow shock in front of the obstacle and the
classical odd-even decoupling situation near the leading shock wave. This makes the case particularly interesting to study
in the context of this paper.

The shock wave is initialized with

(cid:16)

(ρ, u, v, p) =

(3.714, 2.464, 0, 9.310)
(1, 0, 0, 1)

if x < 0.375

else,

(30)

and the center of the diamond is placed at x = 0.7 and  y = 1.6 with a distance  D = 0.6 from corner to corner. The domain
size is set to [0, 2.2] × [0, 3.2] and it is discretized with 7040 × 10240 cells. The ﬁnal time is set to 0.5. Neumann boundary
conditions with zero gradients for all variables are applied at the lower and upper boundary. Inﬂow and outﬂow conditions
are imposed at the left and right boundary. The reﬂecting-wall condition representing the diamond is again approximated
using a level-set approach [41].

Fig. 17 shows  the  ﬁnal  schlieren  images  of  density  gradients  using  both  HLLC  and  HLLC-LM.  An  obvious  disturbance
behind  the  leading  shock  wave  develops  when  the  classical  HLLC  approximation  is  applied.  This  is  caused  by  an  odd-
even  decoupling  effect,  similarly  to  the  corner  ﬂow  presented  in  Section 4.  Again,  the  HLLC-LM  ﬂux  fully  removes  the
disturbance. Note, that the complex ﬂow evolution is not affected by the low Mach number correction. Further details can
be observed within the double Mach reﬂection zone as shown in the zoomed region given in Fig. 18. The proposed HLLC-LM
scheme results in a stable and disturbance-free ﬂow ﬁeld behind the leading Mach stem. Moreover, the decreased numerical
dissipation of the HLLC-LM ﬂux becomes apparent when the resolution of the wave patterns is compared.

16

N. Fleischmann, S. Adami and N.A. Adams

Journal of Computational Physics 423 (2020) 109762

Fig. 18. Zoom on double Mach reﬂection in supersonic ﬂow around a diamond-shaped obstacle with Ma = 2.85: logarithmic gradients of density from 1 to
1000 at t = 0.5.

5.2.  Shock interface interaction: helium bubble in air

Another important application of high-order methods is the prediction of multi-component ﬂows. In [8], it was shown
that the grid-aligned shock instability limits the numerical investigation of shock-interface interaction problems. The same
case  of  the  interaction  of  a  Mach  6  shock  wave  in  air  (γ = 1.4)  with  a  helium  bubble  (γ = 1.66)  is  now  studied  with
HLLC-type solvers.

Initial states are given by

(cid:12)

(ρ, u, v, p) =

⎧
⎪⎨

⎪⎩

√
36 , 0, 251

35

6

216
41 , 35
(1, 0, 0, 1)
(0.138, 0, 0, 1)

(cid:13)

air post-shock,

air pre-shock,
helium,

(31)

where the shock is placed initially at x = 0.05. A helium bubble with initial diameter D = 0.05 is placed at x = 0.1,  y = 0.15
within  in  a  domain  of  size  [0, 0.4] × [0, 0.3].  Inﬂow  and  outﬂow  conditions  are  applied  at  the  left  and  right  boundary,
respectively. Neumann boundaries with zero gradient for all quantities are set at the remaining boundaries. The resolution
is set to 1280 × 960, which resolves the helium bubble with 160 cells per diameter. The ﬁnal time of the simulation is set
to 0.035.

Fig. 19 shows the ﬁnal density results for both HLLC and HLLC-LM. The numerical instabilities at the shock front induced
by the HLLC approximation are not as dominant as for the Roe approximation [8]. Especially, no carbuncles can be observed.
Instead, an odd-even decoupling develops in the backﬂow of the shock wave similar to the one observed for the previous
example. As before, the HLLC-LM scheme produces a clean shock front without any disturbances. Moreover, the stability of
the proposed scheme is tested for an extreme resolution of  1280 cells per diameter. The results shown in Fig. 20 still do
not indicate any instability.

5.3.  Shock interface interaction: air bubble in water in two and three dimensions

Finally, the challenging simulation of a strong 1.6 GPa shock wave in water interacting with an embedded air bubble was
studied. The strong transmitted shock wave in air may suffer from the grid aligned-shock instability. First, the simulations
in [8] were repeated in two dimensions with the HLLC-type solvers. Afterwards, a new fully three-dimensional simulation
of the problem is presented.

17

N. Fleischmann, S. Adami and N.A. Adams

Journal of Computational Physics 423 (2020) 109762

Fig. 19. Shock interface interaction of a helium bubble in air I: density contours from blue = 0.138 to red = 7.5 at t = 0.035.

Fig. 20. Shock interface interaction of a helium bubble in air II: density contours from blue = 0.138 to red = 7.5 at t = 0.035.

The setup is chosen similar to [45] with initial states

(ρ, u, v, p) =

⎧
⎪⎨

⎪⎩

(cid:14)
(cid:14)
(cid:14)

(cid:15)

1323.65, 661.81, 0, 1.6 · 109
(cid:15)
1000, 0, 0, 105
1, 0, 0, 105

(cid:15)

water post-shock
water pre-shock,
air,

(32)

where  water  is  modeled with  a  stiffened  equation  of  state  (γ = 4.4,  P inf = 6 · 108)  and  air  as  ideal  gas  (γ = 1.4).  The
domain size is set to [0, 0.024] × [0, 0.024], where an air bubble with diameter D = 0.006 is placed in the center. The shock
front is initially placed at x = 0.008. Inﬂow and outﬂow conditions are applied at the left and right boundary, respectively.

18

N. Fleischmann, S. Adami and N.A. Adams

Journal of Computational Physics 423 (2020) 109762

Fig. 21. Shock interface interaction of a air bubble in water at t = 3 · 10

−6: velocity magnitude contours from blue = 0 to red = 2850.

Fig. 22. 3D shock interface interaction of a air bubble in water: velocity magnitude within the air bubble from blue = 0 to red = 3500 at t = 2.6 · 10

−6.

Neumann boundary condition with zero gradient for all quantities is set at the remaining boundaries. The bubble initially is
resolved by 160 cells per diameter and the ﬁnal time is set to 3 · 10

−6.

Velocity magnitude results for the HLLC and HLLC-LM solver are shown in Fig. 21. Similarly to the previous case with
helium,  the  HLLC  approximation  does  not  create  any  carbuncles.  However,  the  ﬂow  behind  the  shock  wave  in  air  is  sig-
niﬁcantly  disturbed.  The  HLLC-LM  solver  enables  a  stable  prediction  of  the  ﬂow  ﬁeld.  Again,  the  stability  of  HLLC-LM  is
further demonstrated by an extremely increased resolution of 1280 cells per diameter. The result of this simulation still is
disturbance-free as shown in the right frame of Fig. 21.

We studied the same setup also in three dimensions with a straightforward extension of the domain in z-direction to
[0, 0.024] × [0, 0.024] × [0, 0.024]. The resolution is chosen identically to the original two-dimensional case with 160 cells
−6.
per diameter. Since the air bubble collapses faster in three dimensions, the ﬁnal simulation time was reduced to = 2.6 · 10
The  results  given  in  Fig. 22 demonstrate  that  the  numerical  instability  is  signiﬁcantly  stronger  for  three-dimensional
simulations  when  the  original  HLLC  solver  is  applied.  Similar  to  the  three-dimensional  results  for  the  Sedov  blast  wave
in  Section 4,  now,  small  carbuncles  can  be  observed,  which  never  occurred  in  any  of  our  two-dimensional  simulations.

19

N. Fleischmann, S. Adami and N.A. Adams

Journal of Computational Physics 423 (2020) 109762

Nevertheless,  the  low  Mach  number  modiﬁcation  in  the  HLLC-LM  scheme  leads  to  stable  and  carbuncle-free  results  as
shown in the right frame of Fig. 22.

6.  Conclusion

In this paper, the general idea that the low Mach number in transverse direction of the shock wave propagation is the
reason for the grid-aligned shock instability has been exploited to design a shock-stable version of the popular HLLC approx-
imate Riemann solver. A simple reduction of non-linear wave speeds as done for the Roe ﬂux would lead to pure upwinding
due to the one-sided deﬁnition of the HLLC ﬂux. Therefore, a centralized formulation of the HLLC ﬂux is proposed. Applying
this centralized formulation does not only avoid the switching, but also allows for a straightforward reduction of nonlinear
eigenvalues.  A  smooth  reduction  of  acoustic  dissipation  is  guaranteed  using  a  sine  function.  The  proposed  version  of  the
HLLC scheme with modiﬁed low Mach number behavior is denoted HLLC-LM. The modiﬁed ﬂux reduces the dissipation dur-
ing the ﬂux calculation in case of low directional Mach number, and fully recovers the original HLLC ﬂux otherwise. Thus,
shock stability is retained by a further reduction of the dissipation of the HLLC approximation.

Results obtained with the centralized formulation have been thoroughly compared to the ones obtained with the classical
formulation  and  found  to  be  identical  with  respect  to  ﬂoating-point  differences  for  all  studied  cases.  The  stability  and
accuracy  of  the  HLLC-LM  ﬂux  has  been  demonstrated  for  a  comprehensive  series  of  test  cases  commonly  related  to  the
grid-aligned  shock  instability.  However,  the  prime  goal  of  the  high-order  methods  as  applied  throughout  the  paper  is  to
simulate  more  complex  ﬂow  situations  than  the  classical  carbuncle  cases.  The  advantages  of  the  HLLC-LM  when  applied
to  supersonic  multi-component  ﬂows  have  been  presented  in  detail.  Stability  can  be  maintained  also  for  extremely  high-
resolved simulations. Although the HLLC ﬂux might still be considered as suitable for most two-dimensional situations due
to the fact that the occurring disturbances are commonly bounded and they rarely lead to large deviations unless resolution
is  drastically  increased,  this  is  not  valid  in  three  dimensions.  The  three-dimensional  simulations  presented  in  this  paper
demonstrate that the HLLC ﬂux is likely to produce severe carbuncles similar to the Roe scheme. The HLLC-LM ﬂux revealed
excellent results also for three-dimensional simulations. Hence, the combination of HLLC-LM with state-of-the-art high-order
methods allows for a robust and accurate simulation of current challenges in high-speed ﬂuid dynamics.

CRediT authorship contribution statement

Nico Fleischmann: Conceptualization, Data curation, Investigation, Methodology, Software, Validation, Visualization, Writ-
ing  - original  draft. Stefan Adami: Project  administration,  Supervision,  Writing  - review  &  editing. Nikolaus A. Adams:
Funding acquisition, Resources, Supervision, Writing - review & editing.

Declaration of competing interest

The  authors  declare  that  they  have  no  known  competing  ﬁnancial  interests  or  personal  relationships  that  could  have

appeared to inﬂuence the work reported in this paper.

Acknowledgements

This project has received funding from the European Research Council (ERC) under the European Union’s Horizon 2020

research and innovation programme (grant agreement No. 667483).

The  authors  gratefully  acknowledge  the  Gauss  Centre  for  Supercomputing  e.V.  (www.gauss -centre .eu)  for  funding  this
project by providing computing time on the GCS Supercomputer SuperMUC at Leibniz Supercomputing Centre (www.lrz .de).

References

[1] S.K. Godunov, A difference method for numerical calculation of discontinuous solutions of the equations of hydrodynamics, Mat. Sb. 89 (3) (1959)

271–306.

[2] P.L. Roe, Approximate Riemann solvers, parameter vectors, and difference schemes, J. Comput. Phys. 43 (2) (1981) 357–372.
[3] A. Harten, P.D. Lax, B.v. Leer, On upstream differencing and Godunov-type schemes for hyperbolic conservation laws, SIAM Rev. 25 (1) (1983) 35–61.
[4] E.F. Toro, M. Spruce, W. Speares, Restoration of the contact surface in the HLL-Riemann solver, Shock Waves 4 (1) (1994) 25–34.
[5] K. Peery, S. Imlay, Blunt-body ﬂow simulations, in: 24th Joint Propulsion Conference, 1988, p. 2904.
[6] J.J. Quirk, A contribution to the great Riemann solver debate, in: Upwind and High-Resolution Schemes, Springer, 1997, pp. 550–569.
[7] A.V. Rodionov, Artiﬁcial viscosity in Godunov-type schemes to cure the carbuncle phenomenon, J. Comput. Phys. 345 (2017) 308–329.
[8] N. Fleischmann, S. Adami, X.Y. Hu, N.A. Adams, A low dissipation method to cure the grid-aligned shock instability, J. Comput. Phys. 401 (2020) 109004.
[9] B. Einfeldt, On Godunov-type methods for gas dynamics, SIAM J. Numer. Anal. 25 (2) (1988) 294–318.
[10] E.F. Toro, Riemann Solvers and Numerical Methods for Fluid Dynamics: A Practical Introduction, Springer Science & Business Media, 2013.
[11] E.F. Toro, The HLLC Riemann solver, Shock Waves (2019) 1–18.
[12] P. Batten, N. Clarke, C. Lambert, D.M. Causon, On the choice of wavespeeds for the HLLC Riemann solver, SIAM J. Sci. Comput. 18 (6) (1997) 1553–1570.
[13] D.S. Balsara, M. Dumbser, R. Abgrall, Multidimensional HLLC Riemann solver for unstructured meshes – with application to Euler and MHD ﬂows, J.

Comput. Phys. 261 (2014) 172–208.

[14] S. Li, An HLLC Riemann solver for magneto-hydrodynamics, J. Comput. Phys. 203 (1) (2005) 344–357.
[15] K.F. Gurski, An HLLC-type approximate Riemann solver for ideal magnetohydrodynamics, SIAM J. Sci. Comput. 25 (6) (2004) 2165–2187.

20

N. Fleischmann, S. Adami and N.A. Adams

Journal of Computational Physics 423 (2020) 109762

[16] X. Hu, N.A. Adams, G. Iaccarino, On the HLLC Riemann solver for interface interaction in compressible multi-ﬂuid ﬂow, J. Comput. Phys. 228 (17) (2009)

6572–6589.

(2017) 46–67.

(2009) 7634–7642.

[17] E. Johnsen, T. Colonius, Implementation of WENO schemes in compressible multicomponent ﬂow problems, J. Comput. Phys. 219 (2) (2006) 715–732.
[18] D.P. Garrick, M. Owkes, J.D. Regele, A ﬁnite-volume HLLC-based scheme for compressible interfacial ﬂows with surface tension, J. Comput. Phys. 339

[19] S.D. Kim, B.J. Lee, H.J. Lee, I.-S. Jeung, Robust HLLC Riemann solver with weighted average ﬂux scheme for strong shock, J. Comput. Phys. 228 (20)

[20] K. Huang, H. Wu, H. Yu, D. Yan, Cures for numerical shock instability in HLLC solver, Int. J. Numer. Methods Fluids 65 (9) (2011) 1026–1038.
[21] Y.-X. Ren, A robust shock-capturing scheme based on rotated Riemann solvers, Comput. Fluids 32 (10) (2003) 1379–1403.
[22] Z. Shen, W. Yan, G. Yuan, A robust HLLC-type Riemann solver for strong shock, J. Comput. Phys. 309 (2016) 185–206.
[23] W. Xie, R. Zhang, J. Lai, H. Li, An accurate and robust HLLC-type Riemann solver for the compressible Euler system at various Mach numbers, Int. J.

Numer. Methods Fluids 89 (10) (2019) 430–463.

[24] S. Simon, J. Mandal, A cure for numerical shock instability in HLLC Riemann solver using antidiffusion control, Comput. Fluids 174 (2018) 144–166.
[25] S. Simon, J. Mandal, A simple cure for numerical shock instability in the HLLC Riemann solver, J. Comput. Phys. 378 (2019) 477–496.
[26] G.-S. Jiang, C.-W. Shu, Eﬃcient implementation of weighted ENO schemes, J. Comput. Phys. 126 (1) (1996) 202–228.
[27] S. Gottlieb, C.-W. Shu, E. Tadmor, Strong stability-preserving high-order time discretization methods, SIAM Rev. 43 (1) (2001) 89–112.
[28] N. Fleischmann, S. Adami, N.A. Adams, Numerical symmetry-preserving techniques for low-dissipation shock-capturing schemes, Comput. Fluids 189

(2019) 94–107.

Elsevier, 2017, pp. 203–231.

[29] H. Guillard, B. Nkonga, On the behaviour of upwind schemes in the low Mach number limit: a review, in: Handbook of Numerical Analysis, vol. 18,

[30] F. Rieper, A low-Mach number ﬁx for Roe’s approximate Riemann solver, J. Comput. Phys. 230 (13) (2011) 5263–5287.
[31] X.-s.  Li,  C.-w.  Gu,  An  all-speed  Roe-type  scheme  and  its  asymptotic  analysis  of  low  Mach  number  behaviour,  J.  Comput.  Phys.  227 (10)  (2008)

5144–5159.

[32] H. Guillard, C. Viozat, On the behaviour of upwind schemes in the low Mach number limit, Comput. Fluids 28 (1) (1999) 63–86.
[33] H. Guillard, A. Murrone, On the behavior of upwind schemes in the low Mach number limit: II. Godunov type schemes, Comput. Fluids 33 (4) (2004)

[34] Z. Chen, X. Huang, Y.-X. Ren, Z. Xie, M. Zhou, Mechanism study of shock instability in Riemann-solver-based shock-capturing scheme, AIAA J. 56 (9)

[35] A. Harten, Adaptive multiresolution schemes for shock computations, J. Comput. Phys. 115 (2) (1994) 319–338.
[36] J.W. Kaiser, N. Hoppe, S. Adami, N.A. Adams, An adaptive local time-stepping scheme for multiresolution simulations of hyperbolic conservation laws,

[37] X.Y. Hu, N.A. Adams, C.-W. Shu, Positivity-preserving method for high-order conservative schemes solving compressible Euler equations, J. Comput.

655–675.

(2018) 3636–3651.

J. Comput. Phys. X 4 (2019) 100038.

Phys. 242 (2013) 169–180.

[38] F. Kemm, Heuristical and numerical considerations for the carbuncle phenomenon, Appl. Math. Comput. 320 (2018) 596–613.
[39] P. Woodward, P. Colella, The numerical simulation of two-dimensional ﬂuid ﬂow with strong shocks, J. Comput. Phys. 54 (1) (1984) 115–173.
[40] F. Kemm, On the proper setup of the double Mach reﬂection as a test case for the resolution of gas dynamics codes, Comput. Fluids 132 (2016) 72–75.
[41] X.Y. Hu, B. Khoo, N.A. Adams, F. Huang, A conservative interface method for compressible ﬂows, J. Comput. Phys. 219 (2) (2006) 553–578.
[42] L. Fu, A low-dissipation ﬁnite-volume method based on a new TENO shock-capturing scheme, Comput. Phys. Commun. 235 (2019) 25–39.
[43] E.J. Tasker, R. Brunino, N.L. Mitchell, D. Michielsen, S. Hopton, F.R. Pearce, G.L. Bryan, T. Theuns, A test suite for quantitative comparison of hydrodynamic

codes in astrophysics, Mon. Not. R. Astron. Soc. 390 (3) (2008) 1267–1281.

[44] I. Glass, J. Kaca, D. Zhang, H. Glaz, J. Bell, J. Trangenstein, J. Collins, Diffraction of planar shock waves over half-diamond and semicircular cylinders: an

experimental and numerical comparison, in: AIP Conference Proceedings, vol. 208, AIP, 1990, pp. 246–251.

[45] O. Haimovich, S.H. Frankel, Numerical simulations of compressible multicomponent and multiphase ﬂow using a high-order targeted ENO (TENO)

ﬁnite-volume method, Comput. Fluids 146 (2017) 105–116.

21

