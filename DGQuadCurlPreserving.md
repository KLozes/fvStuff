A curl preserving ﬁnite volume scheme by space velocity
enrichment. Application to the low Mach number accuracy
problem.

Jonathan Jung∗ and Vincent Perrier†

Abstract

In this article, we address the problem of accuracy of ﬁnite volume schemes in the low
Mach number limit. It has been known for years that collocated ﬁnite volume schemes are
naturally correctly behaving in this limit on triangular meshes [21, 22, 16], but fail in general
on other types of mesh. We are ﬁrst interested in the general problem of the conservation of
vorticity for the wave system. By enriching the approximation space for vectors, we prove that
the Hodge-Helmholtz context developed for triangular meshes in [16] can be recovered in the
quadrangular mesh case. This leads to a numerical scheme for the wave system that naturally
preserves the vorticity under mild assumption on the numerical ﬂux. The new approximation
space is then used with the barotropic Euler system. Numerical tests show that the new
numerical scheme is accurate for both steady and acoustic problems at low Mach number.

Contents

1 Introduction

2 Cartesian and periodic case

2.1 A new ﬁnite element space for vectorial ﬁeld . . . . . . . . . . . . . . . . . . . . . .
2.2 Approximation space on Cartesian mesh . . . . . . . . . . . . . . . . . . . . . . . .
2.3 Discrete Hodge-Helmholtz decomposition . . . . . . . . . . . . . . . . . . . . . . .
2.4 First order wave system . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
2.4.1 Discretization . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
2.4.2 Godunov stabilization and stationary solution . . . . . . . . . . . . . . . . .
2.4.3 Conservation of the adjoint curl with the Godunov’ ﬂux . . . . . . . . . . .

3 Quadrangular case with boundary conditions

. . . . . . . . . . . . . . . . . . . . . . . .
3.1 Discretization and boundary conditions
3.2 Discrete Hodge-Helmholtz decomposition . . . . . . . . . . . . . . . . . . . . . . .
3.2.1 The two Piola transformations and the choice of the approximation space .
3.2.2 Determination of uϕ
. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
h
3.2.3 Existence and uniqueness of the decomposition . . . . . . . . . . . . . . . .
3.3 Structure preserved and long time behaviour with Godunov’ ﬂux . . . . . . . . . .

4 Discretization for the barotropic Euler system

5 Numerical results

5.1 Order of accuracy on the adjoint curl . . . . . . . . . . . . . . . . . . . . . . . . . .
5.2 Wave equation . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
5.2.1 Periodic vortex . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
5.2.2 Cylinder scattering . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
5.3 Euler equation . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
5.3.1 Cylinder scattering . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
5.3.2 Propagation of a low Mach number acoustic wave over a steady vortex . . .

2

3
4
4
5
8
8
9
10

11
11
12
12
13
14
17

19

19
20
20
20
21
23
23
25

6 Conclusion

26
∗LMAP, UPPA, Pau, France and Cagire team, Inria Bordeaux Sud-Ouest, France jonathan.jung@univ-pau.fr.
†Cagire team, Inria Bordeaux Sud-Ouest, France and LMAP, UPPA, Pau, France, vincent.perrier@inria.fr.

1

1

Introduction

In this article, we are interested in two types of problems that are strongly interconnected: the
accuracy problem at low Mach number of upwind ﬁnite volume numerical schemes for the compress-
ible barotropic Euler system, and the conservation of the vorticity in the ﬁrst order formulation of
the wave system. The strong connection between these two problems was for example discussed
in [27], and we mainly focus this state of the art section on the long time accuracy problem of the
wave system:

1
ρ0

∂τ p +
∂τ u + κ0∇




divxu = 0,

xp = 0,

(1)

where u is the velocity and p the pressure. System (1) depends on two strictly non-negative
parameters, κ0 and ρ0. The wave velocity is c0, linked with the parameters of the system by
u) = 0,
0 = κ0/ρ0. By taking the curl of the velocity equation of (1), we formally ﬁnd ∂τ (
c2
u is conserved. Conservation of this quantity is a necessary
which means that the vorticity
condition for having long time accuracy.



∇

∇

×

×

x

x

The problem of conservation of the vorticity for (1) has been tackled by several strategies,
some of them being linked with the conservation of the divergence of a vector v when v ensures
the following type of equation

These diﬀerent strategies may be gathered into the three following families

∂τ v +

x

∇

×

g = 0.

(2)

1. Projection method. Projection methods have been widely used in the context of preser-
vation of a divergence for incompressible ﬂows [6] but also for magnetohydrodynamic system
[10]. This method was rather designed for (2), and consists in computing a candidate update
vn+1 which is then corrected. For ensuring preservation of the divergence, a potential ϕ is
computed such that

vn+1

vn

∆ϕ = divx
and the update is projected as vn+1 = vn+1
xϕ, which leads to the conservation of the
(cid:0)
− ∇
divergence of v. As far as we know, it has never been proposed for the conservation of the curl
for (1), but could be easily adapted. This method raises several problems: ﬁrst, this requires
to solve the elliptic problem (3) at each time step, which is costly, and second, (3) should be
equipped with appropriate boundary conditions, which is not always straightforward.

−

(cid:1)

,

(3)

2. Discretizations based on staggered data. This type of discretization consists in having
data located on the faces or edges of the cells, whereas other data are located in the center of
the cells. For example, the MAC scheme [30] consists in discretizing the pressure at the center
of cells whereas the velocity components are normal to each face. The MAC scheme has been
thoroughly analyzed in [35, 36], and extended to compressible ﬂows in [18, 20, 19]. Similar
approach were designed in [24, 15, 23, 32], in which general deﬁnitions of discrete gradient,
divergence and curl were provided with a staggered design. Application of staggered schemes
to hyperbolic problems can be found for example in [8]. A high order discontinuous Galerkin
version of staggered schemes was proposed in [41], and high order ﬁnite volume schemes were
proposed in [2]. Still linked with staggered ideas, the scheme [7] was proposed, for which the
data are collocated, but the divergence is cleaned on a staggered grid. The main drawback
of staggered data is its diﬃculty of implementation on unstructured meshes at high order.

3. Generalized Lagrange multipliers. The general Lagrange multiplier method was de-
signed originally for the conservation of the divergence for the Maxwell system [34] and for
the MHD system [14]. This method consists in adding an extra variable replacing the di-
vergence to be preserved, and to deﬁne an extra equation on this additional variable, which
includes a relaxation between the divergence of the vector and the additional variable. The
extension to the curl constraint preservation was proposed in [17]. The main drawbacks of
this method are the increase of the number of unknowns to solve, and the addition of numer-
ical parameters, that should be tuned (propagation speed of the additional variable, stiﬀness
of the relaxation).

Apart from these types of schemes that are especially designed for conserving divergence or curl
constraints, some cell-centered ﬁnite volume schemes seem to be able to naturally conserve these
constraint. Even if it was not expressed in these terms, [16, 21], which was then extended to high

2

order in [29] are examples of ﬁnite volumes and discontinuous Galerkin methods that are long time
accurate for (1), and so must preserve the curl in a sense that was not yet deﬁned. Note that these
articles deal with triangular or tetrahedral meshes, and are known to fail if other type of meshes
are used (see e.g.
[29], in which the high order discontinuous Galerkin method on quadrangular
and triangular meshes are compared). A second family of collocated schemes seem to be able to
preserve correctly divergence or curl constraints: the node-based solvers, such as [9] for divergence
preservation or [3, 4, 5] for vorticity preservation and application to the low Mach number ﬂows.
These schemes may be seen as an extension of [25, 44, 43] to non Cartesian meshes.

In this article, we wish to ﬁnd an extension of the schemes studied in [16, 21, 29] to quadrangular
meshes. Note that as we are dealing with the two-dimensional case, the deﬁnition of the curl may
be tricky. We therefore deﬁne two operators:

• the rotational operator, denoted as rotx which maps vector ﬁelds to scalar ﬁelds as

rotx u :=

∂uy
∂x −

∂ux
∂y

,

• the curl operator, denoted as curlx, which maps scalar ﬁelds to vector ﬁelds as

curlx f =

∂f
∂y

,

∂f
∂x

−

T

.

(cid:18)
These two operators are such that rotx is the opposite of the adjoint of curlx and inversely. Also,
the quantity that is preserved for (1) is rotx u. The theoretical study of ﬁrst order schemes [16, 21]
can rely on a discrete Hodge-Helmholtz decomposition that reads on periodic domains as

(cid:19)

dPPP0 = R

2

curlx P1 ⊕ ∇

⊕

xCR,

(4)

where dPPP0 is the set of piecewise constant vectors, P1 is the space of continuous ﬁnite element of
degree 1, and CR is the Crouzeix-Raviart ﬁnite element space [13]. The focus of this article is
on ﬁnding a similar discrete Hodge-Helmholtz decomposition as (4), but on quadrangular meshes.
For this, we propose to enrich the approximation space of velocity.

This article is organized as follows. In section 2, the new approximation space is introduced
for Cartesian meshes, and the Cartesian meshes counterpart of (4) is proven. A ﬁnite volume
numerical scheme for (1) is proposed with this new approximation space, and we prove that this
numerical scheme preserves a curl. An originality of this curl is that it is deﬁned in an adjoint sense.
Then the scheme is extended in section 3 to the quadrangular case with non periodic boundary
conditions in the spirit of [27]. In section 4, the numerical scheme is extended to the barotropic
Euler system. Because the Euler system is nonlinear, and because the new basis for vectors is not
piecewise constant, the numerical scheme is no more a purely ﬁnite volume scheme, but includes
a cell integral as in the discontinuous Galerkin method. The section 5 is dedicated to numerical
results where tests are performed on the wave system with periodic Cartesian meshes and with
general quadrangular meshes. Tests are also performed for the accuracy problem at low Mach
number for stationary problems and for propagation of acoustic waves in a low Mach number ﬂow.
This article ﬁnishes with the section 6, a conclusion.

2 Cartesian and periodic case

On triangular mesh with periodic boundary conditions, the discrete Hodge-Helmholtz decomposi-
tion (4) holds. Moreover, it is adapted to the wave system (1) in the sense that its divergence free
component is preserved over time using a Godunov numerical scheme [16]. On Cartesian mesh,
this decomposition does not exist anymore for piecewise constant functions. In this section, a new
ﬁnite element space that is richer than piecewise constant functions is deﬁned and allows to recover
a discrete Hodge-Helmholtz decomposition. Using this vectorial approximation space, a numeri-
cal discretization of the wave system (1) is proposed. Then, it is proved that this discretization
coupled with a Godunov numerical ﬂux preserves the divergence free component of the discrete
Hodge-Helmholtz decomposition but also the adjoint curl over time.

3

(cid:98)
(5)

(6)

(7)

(8)

2.1 A new ﬁnite element space for vectorial ﬁeld

K = [0; 1]2 the unit square element. The ﬁnite dimensional space vector ﬁelds

Vh on

K

Let
considered matches with the space

S0(

K) deﬁned in [1] by

(cid:98)
S0

K

= Q0

K

2

+ span

(cid:16)

(cid:16)
(cid:17)
(cid:98)
We denote by Σ =

(cid:98)

(cid:17)
(cid:98)
σ1, σ2, σ3}
{

1/2
x
(cid:98)
(cid:98)
−
1/2)
y
(
−
−

(cid:98)
(cid:98)

= span

1
0

,

,

0
1

(cid:18)
(cid:18)
the set of linear forms deﬁned on L2(

(cid:18)(cid:18)

(cid:19)

(cid:19)

(cid:18)
(cid:19)
K)2 by

(cid:98)

.

(cid:19)(cid:19)

x
1/2
−
y + 1/2

−
(cid:98)
(cid:98)

1
0

0
1

σ1(

u) =

u(

x)

(cid:90) (cid:98)K

·

(cid:18)

σ2(

(cid:98)
u) =

(cid:98)
u(

(cid:98)
x)

(cid:98)

(cid:19)

d

x,

d

(cid:98)
x,

(cid:90) (cid:98)K
(cid:98)
u) = 12

σ3(

·

(cid:18)
x)

(cid:98)
u(

(cid:98)
(cid:90) (cid:98)K
is a ﬁnite element. The ﬁnite element

x.
d

(cid:19)

(cid:18)

·

(cid:19)
1/2
x
(cid:98)
−
y + 1/2

−
(cid:98)
(cid:98)

We directly obtain that
(cid:98)
(cid:98)
S0(
then strictly included in the the Raviart-Thomas RT0 ﬁnite element of degree zero [40, 1].
(cid:98)

K, Σ,

K, Σ,

(cid:98)
S0(

K)

(cid:17)

(cid:16)

(cid:16)

(cid:98)

(cid:98)

(cid:98)

(cid:98)

(cid:98)

K)

is

(cid:17)

(cid:98)

2.2 Approximation space on Cartesian mesh

In this section, we are interested in the Cartesian mesh case. For the sake of simplicity, the domain
l, divided into N cells in each direction. We denote by (x0, y0)
considered is a square of size l
the bottom left point of the domain. Then the domain is Ω = [x0; x0 + l]
[y0; y0 + l], the nodes
of the mesh are given by the couples (xi, yj) where

×

×

xi = x0 + ih,
yj = y0 + jh,

i = 0, 1, ..., N

i = 0, 1, ..., N

with h =
linear application

l
N

. All cell ci,j = [xi; xi+1]

[yj; yj+1] is the image of

K = [0; 1]2 under the following

×

Fci,j :

x
y

xi +
yj +

xh
yh

(cid:98)

.

(cid:55)→
We denote by C the set of cells of the mesh and by S the set of sides. For a given vectorial
(cid:98)
V on the reference square, the approximation space of vector ﬁelds on Ω is
approximation space
(cid:98)
deﬁned as

(cid:98)
(cid:98)

(cid:18)

(cid:19)

(cid:18)

(cid:19)

Vh =
(cid:98)
The average and jump of uh
deﬁned by

(cid:110)
∈

L2(Ω)2

C ,
uh
Vh on an interior side S

c
| ∀

∈

∈

V

Fc

|c ◦

uh
S separating two cells cL and cR are

∈

(cid:111)

.

(cid:98)

[[ uh ]]= uh
uh

uh

=

uh
|cL −
|cL + uh
2

|cR ,
|cR
,

[[ uh

uh

·
nS

uh
uh
|cL −
|cL + uh
uh
(cid:0)
2

|cR
|cR

(cid:1)
·

=

nS,

·
nS.

·
where nS is the normal vector to side S going from cell cL to cell cR. We denote by C (S) the set
of the two adjacent cells to a side S

S and by S (c) the set of the four sides of a cell c

C .

(cid:8)(cid:8)

(cid:9)(cid:9)

{{

}}

It was proven in [1] that the best approximation of u in Vh is of order 1 for the L2(Ω) norm if

∈

∈

and only if

V contains

S0. This leads to consider the following approximation space

(cid:98)

dSh(Ω) :=
(cid:98)

c
| ∀
Proposition 1 (Properties dSh(Ω)). dSh(Ω) is a ﬁnite element space approximating L2(Ω)2 at
order one, and its dimension is

|c ◦

uh

uh

Fc

S0

(cid:110)

∈

∈

∈

(cid:111)

(cid:98)

.

L2(Ω)2

C ,

(9)

dim (dSh(Ω)) = 3#C = 3N 2.

∈
nS ]]=

Also, the following properties hold for all uh

dSh

∈

•

•

C ,

∈

divx ( uh

c
∀
c
∀
nS ]] is constant along each side.

|c) = 0
uh

S (c),

|c ·

C ,

∈

∈

S

∀

nS is constant along S. This also means that the jump [[ uh

·

Proof. All these properties are direct consequences of straightforward computations.

4

y

(cid:98)

(0, 1)

4

L

1

T

B

(1, 1)

3

R

2

(0, 0)

(1, 0)

x

(cid:98)

Figure 1: Labels of the vertices (1, 2, 3, 4) and the sides of the references square (B,R,T,L for
Bottom, Right, Top, Left).

2.3 Discrete Hodge-Helmholtz decomposition
In this section, a discrete Hodge-Helmholtz decomposition of the new approximation space dSh(Ω)
is proposed. This result is the Cartesian version of the decomposition (4) that holds on triangular
mesh. We denote by Q1(Ω) the set of continuous, piecewise bilinear Lagrange ﬁnite element

Q1(Ω) =

ψ

0(Ω)

∈ C

c
| ∀

(cid:110)
K) is spanned by the following basis

C ,

ψc := ψ

|c ◦

Fc

∈

∈

(cid:98)

Q1(

K)

,

(cid:111)

(cid:98)

where Q1(

(cid:98)

y
1)(
x
−
−
1),
y
x(
−
y,
x
(cid:98)
(cid:98)
y,
1)
x
(
(cid:98)
−
(cid:98)
(cid:98)
see Figure 1 for the labelling of the vertices. On a Cartesian mesh with periodic boundary condi-
tions, the dimension of Q1(Ω) is the number of vertices, namely

y) = (
y) =
y) =
(cid:98)
y) =
(cid:98)
(cid:98)
(cid:98)

ψ1(
ψ2(
(cid:98)
ψ3(
(cid:98)
ψ4(
(cid:98)
(cid:98)

x,
x,
x,
(cid:98)
x,
(cid:98)
(cid:98)
(cid:98)





(10)

(cid:98)
−

1),

−

(cid:98)

(cid:98)

We denote by RaTu(Ω) the Rannacher-Turek [39] ﬁnite element space

dim Q1(Ω) = N 2.

(11)

RaTu(Ω) =

ϕ

(cid:26)
The space RaTu(

L2(Ω)

C ,

c
| ∀

∈

∈

ϕc := ϕ

|c ◦

Fc

∈

RaTu(

K) and

S ,

S

∀

∈

[[ ϕ ]] = 0

.

(cid:90)S

(cid:27)

(cid:98)
K) is spanned by the following basis functions

(cid:98)

(cid:98)

ϕB(

x,

y) =

(cid:98)
ϕR(

(cid:98)
x,

(cid:98)
y) =

(cid:98)
ϕT (

(cid:98)
x,

(cid:98)
y) =

(cid:98)
ϕL(

(cid:98)
x,

(cid:98)
y) =






3
4

1
4

1
4

−

−

y

−

(cid:98)

+

x

(cid:98)

+

y

3
4 −

x

(cid:98)

−

+

−

+

3
2 (cid:32)(cid:18)
3
2 (cid:32)(cid:18)
3
2 (cid:32)(cid:18)
3
2 (cid:32)(cid:18)

x

(cid:98)
x

(cid:98)
x

(cid:98)
x

−

−

−

−

1
2

1
2

1
2

1
2

2

(cid:19)
2

(cid:19)
2

(cid:19)
2

−

(cid:18)

−

(cid:18)

−

(cid:18)

−

(cid:18)

(cid:19)

y

(cid:98)
y

(cid:98)
y

(cid:98)
y

−

−

−

−

1
2

1
2

1
2

1
2

2

(cid:33)

(cid:19)

2

(cid:33)

(cid:19)

2

(cid:33)

(cid:19)

2

(cid:33)

(cid:19)

,

,

,

,

(12)

(cid:98)
and its degrees of freedom are the integrals on each side. On a Cartesian mesh with periodic
boundary conditions, the dimension of RaTu(Ω) is the number of sides, namely

(cid:98)

(cid:98)

(cid:98)

(cid:98)

(cid:98)

dim RaTu(Ω) = 2N 2.

(13)

5

Proposition 2 (Exact discrete Hodge-Helmholtz decomposition of dSh(Ω) in the Cartesian peri-
odic case). On a uniform Cartesian mesh

h of Ω, the following decomposition holds

dSh(Ω) = R

2

⊕ ∇

curlx (Q1(Ω)) .

⊕

M
x (RaTu(Ω))

More precisely, a velocity ﬁled uh

∈

dSh(Ω) can be uniquely written as

uh =

α
β

+

(cid:19)

(cid:18)

ϕ
h [uh] +

P

ψ
h [uh]

P

R2 and

where (α, β)

∈

•

•

ϕ
h [uh] is the gradient of a Rannacher-Turek scalar potential, i.e.
x is applied cellwise:

RaTu(Ω) where

P
ϕh

∈

∇

ϕ
h [uh] =

P

∇

xϕh with

C

c
∀

∈

xϕh)

(
∇

|c :=

∇

x ( ϕh

|c)

ψ
h [uh] is the curl of a continuous Q1 scalar function, i.e.
∂yψ, ∂xψ)T is applied cellwise:

P
Q1(Ω) where curlx(ψ) = (

ψ
h [uh] = curlx(ψh) with ψh

P

∈

|c) .
Moreover, this decomposition is orthogonal for the L2(Ω) scalar product

|c = curlx ( ψh

(curlx ψh)

∈

−
c
∀

C

Proof. We ﬁrst prove that if ψ
where
Q1(

K), we have

ψc

∈

(cid:98)

(cid:98)

u, v
(cid:104)

(cid:105)

=

v.

u

·

(cid:90)Ω

Q1(Ω), curlx ψ

∈

∈

dSh. Since for all cell c, we have ψ

curlx ( ψ

|c) =

1
h

curl

ψc.

(cid:98)x

ψc

|c =

◦

1

F−
c

(cid:98)

(14)

Then it is suﬃcient to prove that the curl of each basis function of Q1(
the curl of each of the functions deﬁned in (10), we get

(cid:98)

K) is included in

S0. Taking

(cid:98)

(cid:98)

S0,

curl

(cid:98)x(

ψ1) =

curl

(cid:98)x(

ψ2) =
(cid:98)

curl

(cid:98)x(

ψ3) =
(cid:98)

curl

(cid:98)x(

ψ4) =
(cid:98)

(cid:18)

(cid:18)

(cid:18)

(cid:18)






x + 1
−
1
y

∈

(cid:98)

∈

∈

(cid:19)

S0,
(cid:98)

(cid:19)
S0,

−
x
(cid:98)
y + 1
(cid:98)
−
x
(cid:98)
−
y
(cid:98)
x
(cid:98)
−
y
(cid:98)
−
(cid:98)
dSh. Since for all cell c, ϕ
(cid:98)
∈

S0.

(cid:19)
1

(cid:98)
∈

(cid:19)

(cid:98)

We now prove that if ϕ
∈
K), we have
ϕc

RaTu(

∈

RaTu(Ω), we have

(cid:98)

xϕ

∇

(cid:98)

|c) =
(cid:98)
Then, it is also suﬃcient to prove that the gradient of each basis function of the Rannacher-
. Taking the derivative of the basis functions
Turek ﬁnite element on

x ( ϕ

ϕc.

∇

(cid:98)

(cid:98)

(15)

1
h ∇(cid:98)x

(

ϕB,

ϕR,

ϕT ,

K is included in
ϕL) deﬁned in (12) leads to
(cid:98)

S0

K

|c =

1
F−
c

ϕc

◦

where

(cid:98)

(cid:98)

(cid:98)

(cid:98)

ϕB =

∇(cid:98)x

ϕR =
(cid:98)

∇(cid:98)x

ϕT =
(cid:98)

∇(cid:98)x

ϕL =
(cid:98)

∇(cid:98)x

(cid:18)

(cid:18)

(cid:18)

(cid:18)

(cid:98)






(cid:16)

(cid:17)

(cid:98)

(cid:98)

x + 3/2
3
−
5/2
y
3
−
1/2
x
3
(cid:98)
−
y + 3/2
3
(cid:98)
−
x + 3/2
3
(cid:98)
−
1/2
y
3
(cid:98)
−
5/2
x
3
(cid:98)
−
y + 3/2
3
(cid:98)
(cid:98)
(cid:98)

−

6

∈

∈

∈

∈

(cid:19)

(cid:19)

(cid:19)

(cid:19)

,

,

,

.

S0

K

(cid:17)

(cid:17)

(cid:17)

(cid:17)

(cid:16)

K
(cid:98)

(cid:16)

K
(cid:98)

(cid:16)

K
(cid:98)

(cid:16)

(cid:98)

S0
(cid:98)

S0
(cid:98)

S0
(cid:98)

(cid:98)

We now prove that the sum is orthogonal. For ψh

Q1(Ω) and ϕh

∈

∈

RaTu(Ω), we have

xϕh)

(
∇

·

(cid:90)Ω

(curlx ψh) =

=

=

(
∇

xϕh)

|c ·

(cid:90)c

(ϕh)

S

S (c)(cid:90)S
∈
(cid:80)

C
c
∈
(cid:80)
C
c
∈
(cid:80)

−
C
c
∈
(cid:80)

(cid:90)c
[[ ϕh curlx (ψh)

·

S
S
∈
(cid:80)

(cid:90)S

(curlx ψh)

nS

|c
|c curlx ( (ψh)
(ϕh)

|c)
|c divx (curlx ( (ψh)
nS ]],

·

|c))

where the Gauss’s theorem was used for the second computing line. As ψh
nS is constant along each side. This leads to
each side, so that curlx ψh

Q1, it is linear along

∈

·

xϕh)

(
∇

·

(cid:90)Ω

(curlx ψh) =

curlx (ψh)

[[ ϕh ]].

nS

·

(cid:90)S

S
S
∈
(cid:80)

Last, by deﬁnition of the degrees of freedom of RaTu(Ω), we have

which proves that

S

S

ϕh

RaTu(Ω)

∀

∈

∀
xϕh and curlx ψh are orthogonal.

∈

[[ ϕh ]] = 0,

(cid:90)S

∇

of

Following exactly the same lines, we can prove that any (α, β)
xRaTu(Ω).
∇
It remains to deal with the orthogonality of (α, β)
∈
consider ψh an element of Q1(Ω), and compute as follows

R2 is orthogonal to any element

∈

R2 and any element of curlx Q1(Ω). We

α
β

·

(cid:19)

(cid:90)Ω (cid:18)

curlx ψh =

=

α
β

·

(cid:19)

(cid:90)Ω (cid:18)

curlx ψh =

curlx ψh

·

α
β

(cid:90)c (cid:18)

C
c
∈
(cid:80)
−

C
c
∈
(cid:80)

(cid:90)c (cid:18)

divx

(cid:19)
α
β

⊥

(cid:19)

(cid:32)(cid:18)
ψh

=

−

=

−

C
c
∈
(cid:80)

(cid:90)c

S

C
c
∈
(cid:80)

−

S
S
∈
(cid:80)

S (c)(cid:90)S
∈
(cid:80)
[[ ψh

(cid:90)S

(cid:18)

(cid:18)

α
β

⊥

·

(cid:19)

xψh

⊥

ψh

· ∇
α
β

(cid:19)
α
β

nS

(cid:33)

⊥

·

(cid:19)
nS ]].

⊥

α
β

In this last equality,

nS is constant along the side, and so can be put out of both the jump
and the integral. But then it remains the jump of ψh, which vanishes because ψh is continuous.
We have thus proven that curlx Q1(Ω) is orthogonal to any (α, β)T
R2, and this ends the proof
of orthogonality.

(cid:18)

(cid:19)

∈

·

For the moment, we have proven that

2

R

⊕ ∇

x (RaTu(Ω))

curlx (Q1(Ω))

dSh(Ω).

⊂

⊕

It remains to prove the equality, and for this, we rely on the dimension of the diﬀerent approxima-
tion spaces.

First, we consider

xϕh = 0, then ϕh is constant
∇
on each cell. Still, as RaTu(Ω) is such that the average of the jump is zero along each face, ϕh
is actually constant on the whole domain. This means that dim (ker
x (RaTu(Ω))) = 1, so that
using (13):

RaTu(Ω) is such that

x (RaTu(Ω)). If ϕh

∇

∈

∇
1 = 2N 2

1.

−

∇
Doing the same for curlx Q1(Ω), and using (11), we ﬁnd

dim (

x (RaTu(Ω))) = dim RaTu(Ω)

−

This leads to

dim (curlx (Q1(Ω))) = dim Q1(Ω)

1 = N 2

1.

−

−

dim Sh(Ω)

dim (

∇

−

x (RaTu(Ω)))

−

dim (curlx (Q1(Ω)))

dim R

2 = 0,

−

and proves the equality of dimensions, which ends the proof.

7

Corollary 2.1. For all uh

dSh(Ω), we have

∈

ϕ
h [uh] = 0

∈
h [uh] = 0, from Proposition 2 we can write

⇔

P

∀

ϕ

S

S ,

nS ]] = 0.

[[ uh

·

Proof. If

P

uh =

α
β

(cid:18)

(cid:19)

+ curlx ψ

R2 and ψ

for some (α, β)
for all S

∈
S , [[ curlx ψ
We now assume that uh

∈
we know that (α, β)

∈

nS ]] = 0 and then [[ uh

nS ]] = 0.
dSh(Ω) is such that for all S
RaTu(Ω) and ψ

·
∈
R2, ϕ

·

∈

Q1(Ω) exist such that

·

∈

∈

∈

uh =

α
β

(cid:18)

+

∇

(cid:19)

xϕ + curlx ψ.

As [[ curlx ψ

·

n ]] = [[ (α, β)T

n ]] = 0, [[

·

xϕ

·

∇

n ]] = [[ uh

·

n ]] = 0. Then

Q1(Ω). Moreover, as recalled in the proof of Proposition 2, we have

S , [[ uh

nS ]] = 0. From Proposition 2,

2 =

xϕ
(cid:107)

(cid:90)Ω (cid:107)∇

=

=

2 =

xϕ
(cid:107)

(cid:90)Ω (cid:107)∇

(cid:90)c ∇

C
c
∈
(cid:80)
C (cid:18)(cid:90)c
c
∈
(cid:80)
C
c
∈
(cid:80)
S
S
∈
(cid:80)

(cid:90)S

(cid:90)

xϕ

xϕ

· ∇

divx (ϕ

∇

xϕ)

ϕdivx (

∇

−

(cid:90)c

xϕ)

(cid:19)
xϕ)

ϕdivx (

∇

ϕdivx (

xϕ) .

∇

(cid:90)c

(cid:90)c

C
c
∈
(cid:80)
C
c
∈
(cid:80)

xϕ

ϕ

∇

·

n

−

S (c)

[[ ϕ

n ]]

xϕ

·

−

∇

The second sum is zero, because as stated in Proposition 1, the divergence of all the basis functions
n is constant on each side,
is zero. Concerning the ﬁrst sum, still following Proposition 1,
which gives

xϕ

∇

·

(cid:90)Ω∇
and by hypothesis, this last integral vanishes on all sides, so that

S
∈
(cid:80)

S ∇

· ∇

(cid:90)S

·

xϕ

xϕ =

xϕ

n

[[ ϕ ]],

2 = 0,

xϕ
(cid:107)

(cid:90)Ω (cid:107)∇

which ends the proof.

2.4 First order wave system

In this section, a numerical discretization of the ﬁrst order waves system (1), based on the new
approximation space dSh, is proposed. Then, it is proved that this numerical scheme, coupled
with a Godunov numerical ﬂux, preserves the divergence free component of the Hodge-Helmholtz
decomposition and preserves the adjoint curl.

2.4.1 Discretization

The discontinuous Galerkin discretization of system (1) on a Cartesian periodic mesh consists in
ﬁnding (ph, uh) in some approximation space Vh = Vh
Vh and all
w

Vh such that for all v

Vh, we have

×

∈

∈

v∂τ ph

−

∂τ uh

w

·

(cid:90)c

(cid:90)c

C
c

∈
(cid:80)

C
c
∈
(cid:80)

where

(cid:90)c

C
c
∈
(cid:80)
−

C
c
∈
(cid:80)

1
ρ0

(cid:90)c

uh

v +

· ∇

S
S
∈
(cid:80)
κ0phdivxw +

(cid:90)S

[[ v ]]

1
ρ0

(cid:18)
[[ w ]]

S
S
∈
(cid:80)

(cid:90)S

uh

nS

·
(cid:9)(cid:9)
phnS

+ d11 [[ ph ]]

= 0,

(cid:19)

+ [[ D22(nS)uh ]]

(cid:8)(cid:8)
κ0

·

(cid:0)

(cid:8)(cid:8)

(cid:9)(cid:9)

is the stabilization matrix. We recall that Godunov stabilization is given by

= 0

(cid:1)

(16)

(17)

(18)

D(n) =

(cid:18)

0

d11
0 D22(n)

(cid:19)

DGodunov(n) =

1
0
0 nnT

(cid:19)

c0
2

8

(cid:18)

while Rusanov stabilization is given by

DRusanov(n) =

c0
2

(cid:18)

1
0

0
I2 (cid:19)

.

(19)

As a consequence of the Proposition 1, using the approximation space Vh(Ω) = dP0(Ω)
×
dSh(Ω), the cell integrals in (16) vanish and discretization (16) can be seen as a ﬁnite volume
scheme that consists in ﬁnding (ph, uh)
dP0(Ω) and all
w

dSh(Ω) such that for all v

dSh(Ω), we have

dP0(Ω)

×

∈

∈

∈

v∂τ ph +

[[ v ]]

(cid:90)c

w

S
S
∈
(cid:80)
∂τ uh +

(cid:90)S

·

1
ρ0

(cid:18)
[[ w ]]

(cid:8)(cid:8)
κ0

·

uh

nS

·
(cid:9)(cid:9)
phnS

+ d11 [[ ph ]]

= 0,

(cid:19)

+ [[ D22(nS)uh ]]

S
S
∈
(cid:80)
Note that, in the case of a Cartesian mesh, the mass matrix is diagonal.

(cid:90)S

(cid:90)c

(cid:9)(cid:9)

(cid:8)(cid:8)

(cid:0)

C
c
∈
(cid:80)
C
c
∈
(cid:80)






(20)

= 0.

(cid:1)

2.4.2 Godunov stabilization and stationary solution

Proposition 3. Using the discretization (20) of the wave system with Godunov stabilization (18),
dSh is a stationary solution of (20) if and only if ph is uniform and for
a state (ph, uh)
all side S

nS ]] = 0 (divergence free velocity ﬁeld).

∈
S , [[ uh

dP0 ×
·

∈

The following proof is an adaptation of the proof performed in [21] for triangular mesh. It can

also be found in [22].

Proof. If (ph, uh)
satisﬁes for all v

dP0 ×
∈
dP0 and all w
∈

dSh,

∈

dSh is a stationary solution of (20) with Godunov stabilization (18), it




(cid:90)S

S
S
∈
(cid:80)
S
S
(cid:90)S
∈
(cid:80)
dP0 ×
C ,

∈


Moreover, for (ph, uh)
∈

c
∀

and

[[ v ]]

(cid:18)

[[ w ]]

·

(cid:16)
dSh, v

1
ρ0

(cid:8)(cid:8)

κ0

uh

nS

·
phnS

+

(cid:9)(cid:9)

+

c0
2
c0
2

(cid:8)(cid:8)
(cid:9)(cid:9)
dP0, and w

∈

∈

[[ ph ]]

= 0,

[[ uh

(cid:19)
nS ]]nS

·

= 0.

(cid:17)

dSh, we have

0 =

1
ρ0

(cid:90)c

divx(uhv) =

S (c) (cid:90)S

(cid:88)S
∈

1
ρ0

uh

·

nSv

(21)

C ,

0 =

κ0divx(phw) =

c
∀

∈

S (c) (cid:90)S
Using test functions vanishing everywhere except on one cell c
equalities from (21) gives

(cid:88)S
∈

(cid:90)c

κ0phw

nS.

·

C , and subtracting the above

∈

C ,

c
∀

∈

w

∀

∈

S0(c),

(cid:98)

S

S (c)(cid:90)S
∈
(cid:80)

S

S (c)(cid:90)S
∈
(cid:80)





1
ρ0

[[ uh

·

1
2
w

−
(cid:18)
nS
·
2

(cid:0)

nS ]] + c0[[ ph ]]

= 0,

κ0[[ ph ]] + c0[[ uh

−

(cid:19)
nS ]]

·

= 0.

(22)

(cid:1)

, and multiplying the second

Last, multiplying the ﬁrst line of (22) by ρ0c0 and using κ0 = ρ0c2
0
line of (22) by

1 lead to

−

1
2
nS
(cid:0)
w
·
2

κ0[[ ph ]]

c0[[ uh

·

−

nS ]]

= 0,

κ0[[ ph ]]

c0[[ uh

−

(cid:1)

·

nS ]]

= 0.

(cid:1)

S , and so

(23)

From Proposition 1, [[ uh

·

C ,

c
∀

∈

w

∀

∈

S0(c),

S

S (c)(cid:90)S
∈
(cid:80)





S

(cid:98)

S (c)(cid:90)S
∈
(cid:80)
nS ]] is constant along each side S

(cid:0)

∈
nS ]],

γS := κ0[[ ph ]]

c0[[ uh

·

−

9

is also constant on each side S
R4
the unknown on the right, top, left and bottom faces of the cell c, (23) leads to the following linear
system to solve

C , denoting by (γR, γT , γL, γB)

S . For a square cell c

∈

∈

∈

1
1
0
1
Since the matrix is invertible, we get γR = γT = γL = γB = 0 and then

1
1
−
0
1

γR
γT
γL
γB

1
0
1
1

1
0
1
1

= 0.

−
−

















−









C ,

c
∀

∈

S

∀

∈

S (c),

κ0[[ ph ]]

c0[[ uh

·

−

nS ]] = 0.

S is a side, denoting by c(cid:48) the neighboring cell to c with respect to the normal nS and

If S
performing the same study for cell c(cid:48) we obtain

∈

Then, for all side S
is obvious.

∈

S , we obtain [[ ph ]] = 0 = [[ uh

κ0[[ ph ]] + c0[[ uh

nS ]] = 0.

nS ]] and the result is proven. The reciprocal

·

·

2.4.3 Conservation of the adjoint curl with the Godunov’ ﬂux

In Proposition 2, the following map is considered:

In this subsection, we consider the adjoint of this operator deﬁned as

curlx : Q1
ψ

dSh(Ω)
curlx ψ.

(cid:55)−→
(cid:55)−→

u

∀

∈

dSh(Ω)

ψ

∀

∈

Q1

ψ (curlx)(cid:63) u =

u

curlx ψ.

(24)

·

(cid:90)Ω
From a practical point of view, (curlx)(cid:63) can be computed by inverting the Q1 mass matrix. As the
operator rotx is the opposite of the adjoint of curlx, (curlx)(cid:63) can also be seen as the approximation
of the opposite of rotx.

(cid:90)Ω

Proposition 4 (Conservation of the adjoint curl). Consider the numerical scheme (20) of the
wave system (1) with Godunov’ stabilization (18). Then

ψ

∀

∈

Q1

curlx ψ

·

∂τ uh = 0,

(cid:90)Ω

which means that ∂τ

(curlx)(cid:63) uh

= 0.

Proof. Taking the equation on the velocity evolution of (20) gives
(cid:1)

(cid:0)

∂τ uh +

w

·

[[ w ]]

·

S
S
∈
(cid:80)

(cid:90)S

C
c
∈
(cid:80)

(cid:90)c

(cid:0)

(cid:8)(cid:8)

κ0

phnS

+ [[ D22(nS)uh ]]

= 0.

(cid:9)(cid:9)
c0
2

nS

nS

(cid:1)

T , which gives

With the Godunov’ numerical ﬂux, we have D22(nS) =

∂τ uh +

w

·

(cid:90)c

C
c
∈
(cid:80)

which can be rewritten

[[ w ]]

·

(cid:16)

S
S
∈
(cid:80)

(cid:90)S

κ0

phnS

+

(cid:1)
nS

(cid:0)
c0
[[
2

(cid:8)(cid:8)

(cid:9)(cid:9)

(cid:0)

(cid:1)

T

uh ]]nS

= 0,

(cid:17)

w

∂τ uh +

·

S
C
S
c
∈
∈
(cid:80)
(cid:80)
Suppose now that w = curlx ψ for ψ
in the side direction, which means that [[ curlx ψ

(cid:90)S

(cid:90)c

∈

·

[[ w

nS ]]

·

κ0 {{

ph

}}

+

c0
2

[[ uh

·

(cid:16)

nS ]]

= 0.

(cid:17)

Q1. As ψ is linear on each side, its gradient is continuous

nS ]] = 0 on each side, which ends the proof.

Note that from an implementation point of view, the quadrature formula used for the mass

matrix should match with the one used for computing the discrete operator (curlx)(cid:63) from (24).

10

3 Quadrangular case with boundary conditions

In this section, we still consider the wave system (1), but the domain Ω is no more periodic, and
we consider two types of boundary conditions:

• Inlet/Outlet boundary conditions, in which a pressure pb and a velocity ub are weakly imposed

• Wall boundary conditions, in which u

ﬂux is 0.

n = 0 is weakly imposed, and the imposed pressure

·

We suppose that pb is uniform, and that

long time limit of the wave system.

ub

·

(cid:90)∂Ωinlet/outlet

n = 0, which is necessary for ensuring a

When dealing with the bounded case, it is important to understand that (1) does not always
preserve the curl of the velocity, because of boundary conditions. The structure that is preserved is
a bit more complicated, and relies on a special Hodge-Helmholtz decomposition that was proposed
in [27], and is recalled here

Proposition 5 (Hodge-Helmholtz decomposition adapted to the wave system, see [27, Prop. 1
& 2] ). Consider the problem (1), with inlet/oulet boundary conditions (pb, ub) or wall boundary
H 1(Ω) a
conditions, and with an initial pressure p0 and velocity u0. For any u, consider ϕ
solution of the variational problem

∈

H 1(Ω)

g
∀

∈

(cid:90)Ω
and which is unique up to a constant, and deﬁne uϕ =
Helmholtz decomposition

∇

(cid:90)Ω∇

· ∇

·

· ∇

−

(cid:90)Ω
ϕ, and uψ = u

g

ϕ =

u

g

gub

n,

uϕ. Then the Hodge-

−

u = uϕ + uψ,

is adapted to the wave system in the sense that the component uψ (divergence free component) is
constant in time. Moreover, the long time limit of the wave system exists and is given by (pb, uψ(0))
where uψ(0) corresponds to the divergence free component of the initial velocity ﬁeld u0.

Note that in Proposition 5, due to boundary conditions, the set of all the uϕ and of all the uψ
are not vectorial spaces, but aﬃne spaces, which means that the decomposition of 0 provided by
Proposition 5 is not 0 + 0.

The aim of this section is to develop a discrete counterpart of the structure depicted in Propo-

sition 5 on quadrangular meshes.

3.1 Discretization and boundary conditions

M

M

h, by Si the set of interior faces of

h a conformal quadrangular mesh, by P the set of points of

h, C the set of cells
We denote by
of
h, and by Sb the set of boundary faces. Boundary faces
are supposed to be oriented, such that the normal is outgoing. Each interior side S is arbitrarily
oriented, and denoting by nS its normal, we will call the left cell the one from which nS is outgoing,
and the right cell the one in which nS is ingoing. The discontinuous Galerkin discretization of the
wave system (1) consists in ﬁnding (ph, uh)
dP0(Ω) and
all w

dSh(Ω) such that for all v

dSh(Ω), we have

dP0(Ω)

M

M

×

∈

∈

∈

v∂τ ph +

C
c
∈
(cid:80)

(cid:90)c

[[ v ]]

1
ρ0

(cid:18)

nS

uh

·

(cid:8)(cid:8)

(cid:9)(cid:9)

S

Si(cid:90)S
∈
(cid:80)

+ d11 [[ ph ]]

(cid:19)






∂τ uh +

w

·

[[ w ]]

·

κ0

phnS

(cid:0)

(cid:8)(cid:8)

(cid:9)(cid:9)

S

Si(cid:90)S
∈
(cid:80)

C
c
∈
(cid:80)

(cid:90)c

where the boundary ﬂuxes in (25) are given by

+

Sb(cid:90)S
S
∈
(cid:80)
+ [[ D22(nS)uh ]]

vL

1
ρ0

(cid:20)

nS

uh

·

(cid:21)b

= 0,

(25)

+

S

Sb(cid:90)S
∈
(cid:80)

(cid:1)
wL

·

κ0phnS

(cid:2)

b = 0
(cid:3)

1
nS
uh
ρ0
·
κ0phnS 


wall





=

κ0ph

·

(cid:18)

0
nS + c0(uh

nS)nS

·

(cid:19)

(26)

11

for wall boundary condition and

1
nS
uh
ρ0
·
κ0phnS 






SW

= 




uh

1
ρ0
nS +

κ0

ph + pb
2

nS + ub
2
(uh

nS

·
c0
2

·

nS

·

ub

·

−

nS)nS



(27)




for inlet/outlet boundary condition. Note that in the case of inlet/outlet boundary condition, a full
state (pb, ub) is enforced. Still, this enforcement is performed weakly through Riemann problems
between the interior state and the imposed boundary state. Details on the derivation of expressions
(27) and (26) can be found in [29, Equation (2.9)] and [27, Appendix A.].

3.2 Discrete Hodge-Helmholtz decomposition

For triangular meshes, a discrete counterpart of the decomposition of Proposition 5 was proven
in [27]. It is based on the deﬁnition of ϕ as follows

Find ϕ

CR

∈

g
∀

∈

CR

xϕ

· ∇

xg =

(cid:90)c ∇

C
c
∈
(cid:80)

S

Si(cid:90)S
∈
(cid:80)

nS ]] +

g[[ u

·

S

Sb(cid:90)S
∈
(cid:80)

g

u

nS

·

ub

·

−

nS

,

(cid:0)

(cid:1)
(28)

where CR is the Crouzeix-Raviart ﬁnite element space. Then taking uϕ
is computed cellwise) and uψ
meshes.

xϕ (where the gradient
gives a discrete counterpart of Proposition 5 for triangular

h = u

h =

uϕ
h

∇

−

For Cartesian meshes with periodic boundary conditions, from Proposition 2 and Theorem 2.1,

we can write all uh

dSh, as

∈

where uϕ

x (RaTu(Ω)) and uψ

h

ensures for all S

h ∈ ∇

We wish to extend this discrete decomposition for general quads and for bounded domain with
boundary conditions of type (27) and/or (26), which ﬁrst raises the problem of the approximation
space to choose for the velocity space.

uh = uϕ

h + uψ

h
Si, [[ uψ

∈

nS ]] = 0.

·

3.2.1 The two Piola transformations and the choice of the approximation space

The contravariant Piola transformation takes

u :

K

(cid:55)→

R2 to u : c

(cid:55)→

R2 deﬁned by

Pdiv
Fc (

u) :=

1
(cid:98)
(cid:98)
(cid:98)xFc)
det (D

(D

(cid:98)xFc)

u

F−
c

1

.

◦

(29)

Using this deﬁnition gives

(cid:98)

(cid:98)

curlx ψ = Pdiv
Fc

curl

ψ

(cid:98)x

On the other hand, the covariant Piola transformation takes
(cid:98)

(cid:16)

(cid:17)
u :

K

(cid:55)→

R2 to u : c

(cid:55)→

R2 deﬁned by

Pcurl
Fc (

u) := (D

(cid:98)xFc)−

T

u

F−
(cid:98)
c

1

.

(cid:98)

◦

(30)

Using this deﬁnition gives

xϕ = Pcurl
Fc (

(cid:98)
∇
are no
As we are considering the general quadrangular case, the transformations Pdiv
Fc
more linear inside the elements, and the vector approximation space generated by the parametric
ﬁnite elements curlx Q1 does not match any more with the parametric
xRaTu. Based on the
proofs of the properties derived in section 2, we consider that the properties ensured by curlx Q1
are important to preserve, which leads us to deﬁne the approximation space dSh of vector ﬁelds
on Ω as

and Pcurl
Fc

(cid:98)
ϕ) .

∇(cid:98)x

∇

(cid:98)

dSh(Ω) =

uh

L2(Ω)2

∈

c
| ∀

∈

C ,

uh

|c ∈

Pdiv
Fc

S0

.

(31)

The contravariant Piola transformation gives for all uh
c

C ,

∈

(cid:110)

dSh(Ω), for all p
(cid:98)

(cid:16)

(cid:17)(cid:111)
∈

dP0(Ω), and all cell

∈

(cid:90)c
where n and

p divx(u) =

p div(cid:98)x(

u)

and

S (c),

S

∀

∈

(cid:90) (cid:98)K

n denote the unit outward normals on S and
(cid:98)

(cid:98)

n =

p u

·

(cid:90)S

(cid:90) (cid:98)S

p

u

n

·

(cid:98)

(cid:98)

(cid:98)

(32)

(cid:98)

12

S.

(cid:98)

Proposition 6 (Finite element for vectors on quadrangular mesh). dSh deﬁned by (31) is a ﬁnite
element space approximating L2(Ω)2 at order one with dim (dSh) = 3#C . Moreover, all uh
dSh
deﬁned by (31) satisﬁes for all cell c

∈

C

∈

divx ( uh

and

|c) = 0
∈
nS ]] is also constant along an interior face S

|c ·

S

∀

uh

S (c),

nS is constant along S.

Si.

∈

Then, the jump [[ uh

·

Proof. This is a direct consequence of Proposition 1 and (32).

3.2.2 Determination of uϕ
h
The determination of the uϕ
component is slightly more complicated than in the Cartesian or in the
h
triangular case. This comes from our choice of approximation space (31), which does not contain
the gradient of the parametric Rannacher-Turek elements. Note that the parametric Rannacher-
Turek ﬁnite element has a lot of drawbacks (for example, it is not optimal order on general quads
[39]), and the deﬁnition of nonconforming ﬁnite elements on general quads is still an active research
topic [12, 33, 26, 31].

In this article we propose to deﬁne directly a subset of the vectorial space (31), based on the
properties on the sides that we wish, instead of relying on a variant of the Rannacher-Turek ﬁnite
element. For this we ﬁrst remark that the solution of (28) can actually be explicitly computed.
Indeed, denoting by ϕS the Crouzeix-Raviart basis function associated to the side S, such that

and denote by uS :=

ϕS. Then changing

∇

ϕ by uϕ
h

, g by ϕS and

∇

∇

g by uS in (28) gives:

ϕS = 1,

(cid:90)S

uϕ
h ·

uS = 


c

C (S)(cid:90)c
∈
(cid:80)

(cid:90)S

(cid:90)S

(cid:0)



[[ u

u

·

·

nS ]]

nS

ub

·

−

nS

if S

if S

(cid:1)

In the quadrangular case, based on (33), for each side S, we deﬁne uS
C (S) (the cells adjacent to the side S), such that

Si

Sb.

(33)

dSh, which support is in

∈

∈

∈

Si

Sb

S

∀

S

∀

∈

∈

u

∀

∈

u

∀

∈

dSh(Ω)

dSh(Ω)

C (S)

(cid:90)

u

u

C (S)

(cid:90)

·

·

uS =

[[ u

uS =

(cid:90)S

(cid:90)S

u

·

nS ]],

·
nS.

(34)

Such a function can be built by just solving the mass-matrix system based on the equalities of (34)
on each cell in C (S) for each basis member of dSh of C (S). Then equation (34) holds for each
u

dSh by linearity.

∈

Deﬁnition 3.1 (Deﬁnition of the space in which uϕ
h

will be searched).

dSϕ

h = Span

uS, S

S

.

∈

). If a vector α

(cid:8)

R#S is such that
(cid:9)

Proposition 7 (Dimension of dSϕ
h

then

This leads to

∈
αSuS = 0,

S
S
∈
(cid:80)

α1 = α2 =

= α#S .

· · ·

h = #S
−
Proof. Consider a linear combination of uS that is zero

dim dSϕ

1.

Consider now a cell c and an element uc of dSh(Ω), which support is reduced to c, such that its
1 along the other side,
trace is 0 along two of the sides, and 1 along one side denoted Si, and

−

αSuS = 0.

S
S
∈
(cid:80)

13

denoted Sj. Such an element uc exists because of Proposition 6: the three degrees of freedom in
the cell c allow to deﬁne a uc such that it is 0 on two sides and 1 on a third side, and as uc is
divergence free, its trace is

1 on the last side. Then

−

uc

·

αSuS

=

αS

uc

uS = αSi

[[ uc

nSi ]] + αSj

(cid:90)Ω

S
S
(cid:18)
∈
(cid:80)
which leads to αSi = αSj

(cid:19)

S
S
(cid:90)Ω
∈
(cid:80)
, and gives dim dSϕ

·

·

(cid:90)Si
#S

(cid:90)Sj
1. We now consider the vector

[[ uc

nSj ]] = αSi −

·

αSj ,

h ≤

−

s :=

uS.

S
S
∈
(cid:80)

Then for all uh

dSh(Ω)

∈

uh

·

(cid:90)Ω

nS ]] +

[[ uh

·

nS

uh

·

S

Sb(cid:90)S
∈
(cid:80)

divxuh

s =

S

=

Si(cid:90)S
∈
(cid:80)
C
c
∈
= 0,
(cid:80)

(cid:90)c

because all the elements of dSh(Ω) have a zero divergence on each cell. This ends the proof.

Now, based on the variational formulation (28), we propose to deﬁne uϕ
h

as follows

Proposition 8 (Deﬁnition of uϕ
h

). We denote by uh and element of dSh(Ω). Suppose that

n = 0. Then the problem to ﬁnd β

R#S /(1, 1 . . . 1) such that vh =

∈

βSuS, ensuring

S
(cid:80)

ub

·

(cid:90)∂Ω

#S

α

∀

∈

R

/(1, 1 . . . 1)

uS =

v

·

αS

S
(cid:80)

(cid:90)Ω

S

Si
∈
(cid:80)

αS

[[ uh

(cid:90)S

nS ]] +

·

αS

(cid:90)S

S

Sb
∈
(cid:80)

(uh

ub)

·

−

nS,

(35)

has a unique solution βϕ, and we denote by

uϕ

h :=

βϕ

SuS.

S
(cid:80)
Proof. The proof relies on the Lax-Milgram theorem. The right hand side is clearly a linear form.
R#S ,
The left hand side is clearly a positive bilinear form. If the problem (35) is considered for α
n = 0, the
then the problem is singular because of Proposition 7. Note however that as
R#S has an inﬁnite number of solutions. Last, the problem (35)

problem (35) considered for α
is coercive because we are working on R#S /(1, 1 . . . 1), on which uϕ

(cid:90)∂Ω
h = 0 if and only if β = 0.

ub

∈

∈

·

3.2.3 Existence and uniqueness of the decomposition

We are now able to deﬁne a Hodge-Helmholtz decomposition on general quadrangular mesh with
boundary condition.

Proposition 9 (Discrete Hodge-Helmholtz decomposition). Let Ω a connected bounded domain
N hole(s) and boundary conditions of type (26) or/and (27). Suppose that
with possibly r
n = 0. Then a velocity

n is known on the boundary such that

a boundary value ub

ub

∈

uh

∈

dSh(Ω) on a general quadrangular mesh

h can be uniquely written as

·

(cid:90)∂Ω

·

M
uh = uϕ

h + uψ

h

where

• uϕ

h ∈

dSϕ
h

• uψ

h ∈

dSψ

h such that

dSψ

h :=

vh

dSh

∈

S

∀

∈

Si, [[ vh

·

nS ]] = 0 and

S

∀

∈

Sb, vh

·

nS = ub

(cid:8)

(cid:12)
(cid:12)

14

nS

.

·

(cid:9)

Moreover

dim(dSψ

h ) = 3#C

#S + 1.

−

We note that unlike the periodic case of Proposition 2, due to the boundary conditions, the

decomposition of Proposition 9 is no more orthogonal.

Proof. We ﬁrst address the existence of the decomposition. Given uh, we deﬁne uϕ
h
tion 8, and deﬁne uψ
h

as

as in Proposi-

uψ

h := uh

uϕ
h .

−

It remains to check that uψ
except for one, which we denote as Sk, for which αSk = 1. Then

dSψ
h

h ∈

. For this, the formula (35) is tested with αS = 0 for all S,

• If Sk

∈

Si then the left hand side is

where (34) was used. The right hand side is

uϕ
h ·

uSk =

(cid:90)Ω

(cid:90)Sk

[[ uϕ
h ·

nSk ]],

nSk ]].

[[ uh

·

(cid:90)Sk

This leads to

[[ uψ
h ·

nSk ]] =

(cid:90)Sk

(cid:90)Sk

[[ (uh

uϕ
h )

·

−

nSk ]] = 0,

and as all the elements of dSh have a constant trace on each side, this gives uψ
h ·

nSk = 0.

• If Sk

∈

Sb then the left hand side is still

uϕ
h ·

uSk =

(cid:90)Ω

(cid:90)Sk

[[ uϕ
h ·

nSk ]],

[[ (uh

ub)

·

−

nSk ]].

(cid:90)Sk

nSk ]] =

ub

·

(cid:17)

(cid:90)Sk

[[ (uh

uϕ
h )

·

−

nSk ]] = 0,

whereas the right hand side is

This leads to

[[

and so uψ
h ·

nSk = ub

·

uψ

h −

(cid:90)Sk

(cid:16)
nSk .

, and this proves the existence of the decomposition.

We now would like to prove uniqueness of the decomposition. Suppose that a given uh can be

Therefore, uψ

h ∈
decomposed as

dSψ
h

Then we deﬁne dh such that

uh = uϕ

h + uψ

h = vϕ

h + vψ
h .

dh := uϕ

h = vψ
vϕ

h −

h −

uψ
h ,
h , uψ

h ∈

As [[ dh ]] = 0 on all interior and boundary sides (because vψ

dSψ
h

), this gives

αS

uS = 0.

dh

·

(cid:90)Ω

S
(cid:80)

, and by uniqueness of the solution of (35), we have dh = 0, which proves

As dh = uϕ
uniqueness.

h −

vϕ
h ∈

dSϕ
h

As dim dSh = 3#C and dim dSϕ

h = #S

gives

dim dSψ

h = dim dSh

which ends the proof.

Proposition 10 (Characterization of dSψ
h
null vector given by Proposition 9. Then

1, the existence and uniqueness of the decomposition

dim dSϕ

h = 3#C

#S + 1,

−

−

−

). We denote by uϕ and uψ the decomposition of the

15

• if ψp denotes an element of the canonical Q1 basis function for an interior point p

Pi,

∈

then curlx ψp + uψ

dSψ
h

∈

• if ψ∂k denotes the element of Q1 that is 1 on the whole kth connected component of ∂Ω and

0 elsewhere, then curlx ψ∂k + uψ

dSψ
h .

∈

h can be written as a linear combination of uψ, curlx ψp for p

Reciprocally, any element of dSψ
and curlx ψ∂k .
are such that the jump is 0 along all boundary and interior
Proof. All the curlx ψp and curlx ψ∂k
sides. Then it is clear that by adding uψ, they all belong to dSψ
. Also, we remark that all these
h
functions are not linearly independent, as the function of Q1 equal to 1 on all points belong to this
space. This gives

Pi

∈

dim span (curlx ψp, curlx ψ∂k ) = #Pi + r.

It remains to prove that

For this, we ﬁrst use the Euler relation in the mesh

dim dSψ

h = #Pi + r.

We then remark that the sum

#P

−

#S + #C = 1

r.

−

(36)

S (c)
∈
(cid:80)
can be computed in two ways, and gives 4#C on one hand, and 2#Si + #Sb on the other hand:

C
c
∈
(cid:80)

S

1,

Using this last formula for eliminating #S = #Si + #Sb from (36) leads to

4#C = 2#Si + #Sb.

#P

−

(4#C

−

#Si) + #C = 1

r,

−

which provides the equality

3#C = #P + #Si + r

Last, remarking that #Pb = #Sb, we get

1.

−

3#C

−

#S + 1 = #Pi + r,

which ends the proof.

Proposition 11 (Characterization of (dSϕ

h )⊥). We have

span (curlx ψp, curlx ψ∂k ) = (dSϕ

h )⊥

Proof. If p

∈

Pi, we have for all α

R#S ,

∈

curlx ψp

·

(cid:90)Ω

αSuS

(cid:18)

S
(cid:80)

(cid:19)

=

=

C
c
∈
(cid:80)

(cid:90)c
αS

curlx ψp

αSuS

·

S
(cid:18)
(cid:80)
curlx ψp

(cid:19)
uS

·

S
(cid:80)

c

C (S)(cid:90)c
∈
(cid:80)

=

S

Si
∈
(cid:80)
= 0.

αS

(cid:90)S

[[ curlx ψp

nS ]] +

·

S

Sb
∈
(cid:80)

αS

curlx ψp

(cid:90)S

nS

·

Similarly, if ∂k is a connected component of ∂Ω, we have for all α

R#S ,

∈

(cid:90)Ω

curlx ψ∂k ·

S
(cid:18)
(cid:80)
h )⊥. Since

(dSϕ

αSuS

= 0,

(cid:19)

so that span (curlx ψp, curlx ψ∂k )

⊂
dim span (curlx ψp, curlx ψ∂k ) = #Pi + r = 3#C

#S + 1 = dim (dSϕ

h )⊥ ,

−

it ends the proof.

16

3.3 Structure preserved and long time behaviour with Godunov’ ﬂux
Discretization (25) with Godunov stabilization can be written as ﬁnding (ph, uh)
such that for all v

dSh, we have

dP0(Ω) and all w

dP0(Ω)

×

∈

dSh(Ω)

(37a)

pb)

= 0,

(cid:19)

−

(cid:17)

∈

∈
1
ρ0

1
ρ0

(cid:18)
nS ]]

·

v∂τ ph +

[[ v ]]

S

Si(cid:90)S
∈
(cid:80)

C
c
∈
(cid:80)

(cid:90)c

(cid:18)

C
c
∈
(cid:80)

(cid:90)c

+

∂τ uh +

w

·

+

+

S

S

S

S

v

SSW (cid:90)S
∈
(cid:80)
[[ w

Si(cid:90)S
∈
(cid:80)

w

nS

nS

·

·

(cid:0)

(cid:18)

∈

SWall(cid:90)S
(cid:80)

w

SSW (cid:90)S
∈
(cid:80)






nS

uh

·

(cid:8)(cid:8)
uh

(cid:9)(cid:9)
nS + ub
2

·

·

+

c0
2
nS

[[ ph ]]

(cid:19)

(ph

+

c0
2

κ0 {{

ph

}}

+

(cid:16)
κ0ph + c0uh

c0
2

nS ]]

[[ uh

·

nS

·

κ0

ph + pb
2

+

c0
2

(cid:1)
uh

nS

·

ub

·

−

nS

(cid:0)

= 0.

(37b)

(cid:19)

(cid:1)

The following propositions are the discrete version of Proposition 5.

Proposition 12 (Discrete Hodge-Helmholtz decomposition adapted to the wave system). Denote
dSh(Ω) the solution of the wave system discretization (25) with Godunov
(ph, uh)
stabilization (18). The Hodge-Helmholtz decomposition uh = uϕ
h of Proposition 9 is adapted
to the wave system in the sense that uψ

h + uψ

dP0(Ω)

×

∈

Proof. For all p
·
·
(37b) with forward or backward Euler method for the time integration gives

Pi, since [[ curlx ψp

Si and curlx ψp

∈

∈

nS = 0 for all S

h is constant in time.
nS ]] = 0 for all S

so that

Pi

p

∀

∈

curlx ψp

·

∂tuh = 0

C
c
∈
(cid:80)

(cid:90)c

∈
Using Proposition 11, it gives

∀

p

Pi

curlx ψp

·

(cid:16)

C
c
∈
(cid:80)

(cid:90)c

∂tuϕ

h + ∂tuψ

h

= 0.

(cid:17)

In the same way, we obtain for all ∂k connected component of ∂Ω that

Pi

p

∀

∈

curlx ψp

·

(cid:90)c

C
c
∈
(cid:80)

∂tuψ
h

= 0.

(cid:16)

(cid:17)

(cid:90)c
Then, using the charaterization of dSψ
h

C
c
∈
(cid:80)

which ends the proof.

curlx ψ∂k ·

∂tuψ
h

= 0.

(cid:16)

(cid:17)

of Proposition 10, we obtain with (38) and (39) that

∂tuψ

h = 0

Proposition 13 (Long time limit). We consider the discretization (25) of the wave system (1)
with Godunov stabilization (18). The initial conditions are denoted by (p0
h). Suppose that the
boundary conditions ensure

h, u0

• p = pb where pb is uniform.

• u = ub where

n = 0.

ub

·

(cid:90)∂Ω

Then the long time limit exists and is denoted by (p∞h , u∞h ), it satisﬁes p∞h is uniform (equals to pb)
and u∞h = uψ
h (0) corresponds to the divergence free component of the initial velocity
ﬁeld u0

h obtained with Proposition 9. Moreover, deﬁning the relative energy as

h (0) where uψ

the volumic average of the relative energy is a Liapunov functional.

EUψ

h (0)(Uh) :=

ρ2
0c2
0
2

(ph

−

pb)2 +

1
2

uψ

uh

−

2

,

h (0)
(cid:13)
(cid:13)
(cid:13)

(cid:13)
(cid:13)
(cid:13)

17

Sb,

∈

(38)

(39)

Proof. Since for all v

dP0(Ω) and all w

dSh, we have

∈

∈

divx(vw)

0 =

=

=

C
c
∈
(cid:80)

(cid:90)c

S

Si(cid:90)S
∈
(cid:80)

S

Si(cid:90)S
∈
(cid:80)

vw

nS

·

[[ vw

·

nS ]] +

S

[[ v ]]

nS

w

·

Sb(cid:90)S
∈
(cid:80)
v
+

{{

(cid:8)(cid:8)

(cid:9)(cid:9)

[[ w

·

}}

nS ]] +

nS,

vw

·

S

Sb(cid:90)S
∈
(cid:80)

so that (37a) can be rewritten as

C
c
∈
(cid:80)

(cid:90)c

v∂τ ph +

S

−

S

+

S

Si(cid:90)S (cid:18)
∈
(cid:80)

∈

SWall(cid:90)S
(cid:80)
v

SSW (cid:90)S
∈
(cid:80)
dSψ
h

×

[[ uh

v

}}

·

nS ]] +

c0
2

[[ v ]] [[ ph ]]

(cid:19)

1
ρ0 {{
1
ρ0

uh

−

v

nS

·

1
ρ0

uh

−

nS + ub
2

·

·

nS

+

c0
2

(ph

−

(cid:18)

pb)

= 0.

(cid:19)

Moreover, (pb, uψ

h (0))

∈

P0(Ω)

satisﬁes for all v

dP0(Ω)

∈

1
ρ0 {{

v

−

}}

[[ uψ

h (0)

·

nS ]] +

c0
2

[[ v ]] [[ pb ]]

= 0,

(cid:19)

S

Si(cid:90)S (cid:18)
∈
(cid:80)

and for all w

dSh,

∈

(40)

(41)

nS ]]

κ0pb +

[[ w

·

c0
2

[[ uψ

h (0)

·

nS ]]

+

S

Si(cid:90)S
∈
(cid:80)

=
S

Si(cid:90)S
∈
(cid:80)
=0.

(cid:16)
w

nS

·

(cid:9)(cid:9)

κ0

−

(cid:16)

(cid:8)(cid:8)

[[ pb ]] +

c0
2

[[ w

·

S

SWall(cid:90)S
(cid:17)
∈
(cid:80)
nS ]] [[ uψ
h (0)

nS ]]

·

(cid:17)

nSκ0pb +

w

·

nSκ0pb

(42)

w

·

S

SSW (cid:90)S
∈
(cid:80)

Denoting by pRel
we obtain for all v

h = ph

pb and uRel
−
dP0(Ω) and all w

h = uh

uψ
h (0) and substracting (41) to (40) and (42) to (37b),
−
dSh,

∈

∈

1
ρ0 {{
1
ρ0

uh

−

v

nS

·

[[ uRel
h

v

}}

·

nS ]] +

c0
2

[[ v ]] [[ pRel

h ]]

(cid:19)

v∂τ ph +

C
c
∈
(cid:80)

(cid:90)c

S

Si(cid:90)S (cid:18)
∈
(cid:80)

−

S

+

C
c
∈
(cid:80)

(cid:90)c

∂τ uh +

w

·

+

S

S

S

∈

SWall(cid:90)S
(cid:80)
v

SSW (cid:90)S
∈
(cid:80)
[[ w

Si(cid:90)S
∈
(cid:80)

1
ρ0

uh

−

nS + ub
2

·

·

nS

+

c0
2

(ph

−

κ0

pRel
h

+

c0
2

[[ uRel
h

nS ]]

·

(cid:18)
nS ]]

·

(cid:16)
(cid:8)(cid:8)
κ0(ph

nS

w

(cid:9)(cid:9)

pb) + c0uh

−

nS

·

(43a)

pb)

= 0

(cid:19)

(cid:17)

(ph

−

pb) +

c0
2

uh

·

(cid:1)
nS

−

ub

·

nS

= 0

(43b)

∈

SWall(cid:90)S
(cid:80)

w

+

·

·

nS

κ0
2

(cid:0)

(cid:16)

SSW (cid:90)S
∈
(cid:80)
Multiplying (43a) with v = pRel

S

by ρ2

0c2
0

h

and adding to (43b) with w = uRel

(cid:0)

(cid:1)(cid:17)
h , we obtain

C
c
∈
(cid:80)

(cid:90)c

so that

∂τ EUψ

h (0)(Uh)+

1
2

0c3
ρ2

0[[ pRel

h ]]2 + c0[[ uRel

h

nS ]]2

·

+

+

S

S

S

∈

Si(cid:90)S
∈
(cid:80)
SWall(cid:90)S
(cid:80)
SSW (cid:90)S
∈
(cid:80)

c0

uh

·

2

nS

0c3
ρ2
(cid:0)
0
2

(ph

−

(cid:1)
pb)2 +

c0
2

nS

uh

·

ub

·

−

2

nS

= 0.

(cid:0)

(cid:1)

(44)

∂τ EUψ

h (0)(Uh)

0.

≤

C
c
∈
(cid:80)

(cid:90)c

18






dSh(Ω) the long time limit, (44) gives

Denoting by (p∞h , u∞h )

∈




×

dP0(Ω)
[[ p∞h ]] = 0 and [[ u∞h ·
u∞h ·
p∞h = pb

nS = 0,

and u∞h ·

nS ]] = 0,

nS = ub

nS,

·

Si,
Swall
SSW,

S
S
S

∀
∀
∀

∈
∈
∈

so that p∞h = pb and u∞h = uψ



h (0) by uniqueness of the Hodge decomposition of Proposition 9.

4 Discretization for the barotropic Euler system

The isentropic Euler system is shortly noted as

∂tW +

f (W) = 0

(45)

∇ ·
where W = (ρ, ρu)T is the vector of conservative variables, ρ is the density and u the velocity. f
is the ﬂux deﬁned as

f (W) =

ρu
u + pI

.

ρu
p is the pressure, and is linked with the density through the equation of state: p = p(ρ), which is
p(cid:48)(ρ), and the
supposed to be convex and strictly increasing. The sound velocity c is deﬁned as
Mach number is deﬁned as

⊗

(cid:18)

(cid:19)

The numerical solutions of (45) will be searched in the ﬁnite element space

(cid:112)

/c.

u
|
|

We consider the discontinuous Galerkin method for (45):

Vh = dP0 ×

dSh.

(46)

Find Wh

Vh

∈

ϕ

∀

∈

Vh

(ϕ

∂tWh

f (Wh)

−

·
(cid:90)c
f (Wh, nS) +

C
c
∈
(cid:80)

+

[[ ϕ ]]

·

ϕ)

· ∇
f b(Wh)

ϕLef t

·

nS = 0.

·

S

S

Si(cid:90)S
∈
(cid:80)

Sb(cid:90)S
∈
(cid:80)
The numerical ﬂux
f (Wh, nS) may be any of known numerical ﬂux [42]: Roe, Lax-Friedrich,
HLL, HLLC, exact Godunov... Note that in the case of the wave system (16) and (25), the cell
integration term vanishes, because of the linearity of the ﬂux and the fact that the ﬂux of velocity
of (1) includes only a gradient, which leads to a weak formulation involving only the divergence
of the test functions which is 0 for all test functions. This is no longer the case for the barotropic
Euler system, and the cell integral must be taken into account.

(cid:101)

(cid:101)

The link between the behaviour of (45) at low Mach number and the wave system was exten-
sively discussed in [27, 29]. The main result, which may be obtained by a two-scale asymptotic
expansion of the numerical scheme (46) [29, Section 3.2], is the classiﬁcation of the numerical ﬂux
for Euler with respect to the possible numerical ﬂux for the wave system:

• A ﬁrst family of ﬂux is asymptotically consistent with the Godunov’ ﬂux for the wave system.
This family includes for example the HLLC scheme, the Roe scheme, and the Osher scheme.
For this family of numerical ﬂux, Proposition 12 ensures the preservation of the divergence
free component of the discrete Hodge-Helmholtz decomposition and Proposition 13 ensures
to reach the right long time limit (according to continuous case of Proposition 5) for the wave
system. Then, the numerical scheme (46) is expected to be low Mach number accurate.

• A second family of ﬂux is asymptotically consistent with the Rusanov scheme for the wave
system. This family includes for example the Rusanov, HLL or Lax-Friedrich ﬂux for the
barotropic Euler system. In this case, the numerical scheme (46) is expected to be inaccurate
at low Mach number.

5 Numerical results

In this section, numerical tests are performed to illustrate the results of section 2 and section 3
concerning the discretization (20) of the wave system (1) obtained the new approximation space
dSh. Finally, the low Mach number behavior of the Euler discretization (46) using the new ap-
proximation space is studied.

19

From a numerical point of view, cells integrals of (20) and (46) are performed using a 2

2 = 4
points Gauss’ quadrature formula while side integrals of (20) and (46) are performed using midpoint
approximation. Concerning time integration, all the results presented below were obtained with
a forward Euler scheme. Last, all the initializations are performed with a L2 projection on the
ﬁnite element space. Note that none of the propositions (and especially the main ones, namely
Proposition 4, Proposition 12, Proposition 13) depend on the way the initial condition is done.

×

5.1 Order of accuracy on the adjoint curl
The aim of this test case is to assess the order of approximation given by the adjoint curl, (curlx)(cid:63),
f (r). We
deﬁned in (24). For this, we ﬁrst consider a function f deﬁned in polar coordinate r
then consider the vector ﬁeld deﬁned as

(cid:55)→

We have then

u = curlx f =

∂f
∂r

eθ.

1
r
For assessing the order, we choose a function that is regular and has a compact support in [0, 1]2:

rotx u =

∂f
∂r

∂
∂r

−

(cid:18)

(cid:19)

r

.

f : r

(cid:55)−→

f (r) =

α
1−( r

r0 )2

−

r0e
0

(cid:40)

if
r0
r
≤
otherwise,

where r is the L2 distance deﬁned with respect to a point [xc, yc]
circle of center [xc, yc] and of radius r0 is strictly inside [0, 1]2. Then, denoting by r = r
r0

[0, 1]2, and r0 is such that the

∈

,

(47)

(48)

u =









and

2αy
r0

2αx
(cid:0)
r0

−

−α
1−r2

e

1

−
e

2

2

r2
−α
1−r2
(cid:1)
r2

1

−

(cid:0)
αr2 + r4

4α

(cid:1)
1
4
(cid:1)

−
r2

,









α
1−r2

e−

.

rotx u = −

(cid:1)

(cid:0)

1

∈

×

×

×

−

40, 80

10, 20

20, 40

80, and 160

r2
0
The numerical test is led as follows. We ﬁrst project u of (47) on the ﬁnite element space dSh,
(cid:0)
for ﬁnding uh. Then (curlx)(cid:63) uh
Q1 is computed using (24). Last, the error with respect to the
exact solution (48) is computed. The computation is led on four regular Cartesian meshes with
160 cells. The data used are α = 4 (for avoiding a too
10
sharp behaviour near r0), r0 = 0.15, and xc = yc = 0.5. The error obtained is presented in Table 1,
and shows a second order of accuracy of the adjoint curl. This second order of accuracy of the
adjoint curl ﬁrst surprised us, because the projection of the vector ﬁeld u is done on a piecewise
constant approximation space, for which the order of accuracy is 1, and we could expect that
the error on the adjoint curl be bounded by the initial projection. Still, the approximation space
of (curlx)(cid:63) u is piecewise linear, and so a second order of approximation can be obtained. The
observed order of accuracy probably comes from a commutation property between the projection
operators and the continuous and discrete (curlx)(cid:63) u. The study of this commutation property is
out of the scope of this paper.

×

×

5.2 Wave equation

5.2.1 Periodic vortex

The aim of this section is to assess the result of Proposition 4. For this, we consider the following
vortex initial solution on the quad [0, 1]2

y
r0
−
x
r0

r2

e−

,

r2 

e−



p0 = 0,

u0 =





20

h
0.1
0.05
0.025
0.0125
0.00625
0.003125

Error
0.111664
0.0588788
0.0159502
0.00379819
0.000855902
0.000208047

Local order

−
0.923
1.884
2.070
2.149
2.040

Table 1: Order of accuracy obtained after a projection on dSh and a computation of the adjoint
curl, (curlx)(cid:63).

Figure 2: Quadrangular meshes used for the conservation of the adjoint of the curl test described
in subsubsection 5.2.1: Cartesian on the left and quadrangular unstructured on the right.

with the same notations as in the previous subsection, and with periodic boundary conditions. The
vortex is centered around (xc, yc) = (0.5, 0.5), and we take r0 = 0.15. The test is run until time
t = 1 on a Cartesian mesh of 10
10 cells, and on an unstructured quadrangular mesh with similar
resolution represented in Figure 2. At each time step, the diﬀerence between the adjoint curl of the
velocity, (curlx)(cid:63) u and its initial value (curlx)(cid:63) u0 is computed. The L2 norm of this diﬀerence is
dSh, so that the (curlx)(cid:63) deﬁned in (24) is also deﬁned in
plotted in Figure 3. Note that (dQ0)2
(dQ0)2. As proven in Proposition 4, (curlx)(cid:63) u is conserved when the approximation space is dSh
and the numerical ﬂux is the Godunov’ ﬂux. In the other cases, namely either when the velocity
approximation space is (dQ0)2 or when the Lax-Friedrich numerical ﬂux is used, the adjoint curl
is not preserved.

×

⊂

5.2.2 Cylinder scattering
The purpose of this test case is to illustrate Proposition 13. The domain is an annulus [r0, r1]
×
[0; 2π[ where r0 = 0.5 and r1 = 5.5. Wall boundary condition (26) is applied in r = r0 and

Figure 3: Plot of L2 norm of the diﬀerence between the adjoint curl of the velocity (curlx)(cid:63) u
and the initial adjoint curl of the velocity (curlx)(cid:63) u0 for the periodic test case described in
subsubsection 5.2.1. The left plot matches with a Cartesian mesh and the right plot with a
quadrangular unstructured mesh (the two meshes are represented in Figure 2).

21

0.00.20.40.60.81.0time0.00.51.01.52.02.5k(curlx)?u−(curlx)?u0k2dSh,Lax-Friedrich(dQ0)2,GodunovdSh,Godunov(dQ0)2,Lax-Friedrich0.00.20.40.60.81.0time0.00.51.01.52.02.5k(curlx)?u−(curlx)?u0k2(dQ0)2,GodunovdSh,Lax-FriedrichdSh,Godunov(dQ0)2,Lax-FriedrichFigure 4: Wave cylinder scattering - Energy residual as a function of the computational time.

inlet/outlet boundary condition (27) is applied in r = r1 with pb = 0 and ub = (1, 0)T . Initial data
are given by p0 = 0 and u0 = 0.

In Figure 4, we plot the relative energy residual deﬁned in (44) by

1
2

0c3
ρ2

0[[ ph ]]2 + c0[[ uh

nS ]]2

+

·

c0

uh

2

nS

·

S

S

S

∈

(cid:0)

(cid:1)

×

+

−

uh

nS

(ph

c0
2

0c3
ρ2
(cid:0)
0
2

(cid:1)
pb)2 +

Si(cid:90)S
∈
(cid:80)

SWall(cid:90)S
(cid:80)
SSW (cid:90)S
∈
(cid:80)
as a function of the computational time for the numerical solution (ph, uh). A mesh containing
20 cells, where the ﬁrst number corresponds to radial discretization and the second one
10
to orthoradial discretization, is used. Using (dQ0)2 or the new approximation space dSh, the
numerical schemes converge to a long-time limit (p∞h , u∞h ). This result was proved in [28] for (dQ0)2
with Rusanov and Godunov stabilization. Concerning dSh, no numerical diﬃculty was observed.
In agreement with the result of Proposition 13, the relative energy residual of the long time limit
obtained with dSh and Godunov stabilization is zero such that the long time limit satisﬁes p∞h = pb
and u∞h = uψ
h (0) by uniqueness of the Hodge-Helmholtz decomposition of Proposition 9. Using
(dQ0)2 or dSh with Rusanov stabilization, the energy residual of the long time limit is not zero.
For this test case, the exact long time wave solution can be computed analytically and is given

nS

ub

−

(cid:1)

(cid:0)

2

·

,

·

by

p∞ex(r, θ) = pb = 0

1

r2
0
r2 cos(2θ)
−
r2
0
r2 sin(2θ)

.





r2
0

r2
1

r2
1 −

u∞ex(r, θ) = uψ(0) =






In Figure 5, a mesh convergence study is performed. The computation is performed on ﬁve meshes
160 cells. Using the new
containing respectively 6
approximation space dSh with Godunov stabilization, the long time pressure is exactly zero and
the long time velocity converges to the exact solution with a rate of one. This means that we have
a commutation between the mesh convergence and the long time limit convergence. In the two
other approximation type, namely with (dQ0)2 with the Godunov’ scheme, and with dSh and the
Rusanov scheme, we know that for each mesh, we have a long time limit as proven in [28], but the
long time limit does not converge to the correct solution when the mesh is reﬁned. This means
that in these cases, the limits in long time and the mesh convergence do not commute.

80 and 80

12, 10

20, 20

40, 40




(49)

×

−

×

×

×

×

The results found in Figure 5 deserves additional comments. In [29, Section 5.1.1], the numerical
results obtained on a diﬀerent test case show that the order of accuracy obtained on quads is
generally one order lower than the optimal one. This means that a ﬁrst order of accuracy should
be obtained with the classical DG1 numerical method, still at the price of having eight degrees
of freedom per cell for the velocity (four by velocity component), whereas the numerical method
proposed here has only three degrees of freedom per cell for the velocity. The enrichment procedure
for the velocity space was further discussed in the two submitted article [38, 37] for the higher order
case.

22

050100150200250300350400Time10−2610−2210−1810−1410−1010−610−2102RelativeenergyresidualGodunov,(dQ0)2Godunov,dShRusanov,dShFigure 5: Wave cylinder scattering - L2 norm of the error between the exact long time solution
and the numerical long time solution. A log-log plot is used.

5.3 Euler equation

5.3.1 Cylinder scattering

∈

dP0 ×

We illustrate the good behavior at low Mach numbers of the numerical scheme obtained with the
dSh and a Roe numerical ﬂux. Note that similar results
new approximation space (ρ, ρu)h
are obtained with any other numerical ﬂux that degenerates to a Godunov stabilization for the
asymptotic wave system. The domain and meshes used are the same as in subsubsection 5.2.2.
Wall boundary condition is applied in r = r0 while Steger-Warming boundary condition is applied
in r = r1 with a state characterized by its density at inﬁnity ρb and its velocity at inﬁnity ub =
(ub, 0)T . For all the computations, ρb is set to ρb = 2 and ub is deduced from the Mach number
Mb by ub = Mb
p(cid:48)(ρb). Initial data are uniform and set equal to ρ0 = ρb and u0 = 0. The exact
incompressible solution is given by
(cid:112)

uIncomp

ex

(r, θ) = ubu∞ex(r, θ),

ρIncomp
ex

(r, θ) = ρb + ρb




is given by (49).


where u∞ex

2

r2
1

r2
1 −

r0
0 (cid:19)

(cid:18)

(cid:18)

r2
0
r2 cos(2θ)

−

r4
0
2r4

M 2
b ,

(cid:19)

In Figure 6, the iso-contours of the steady velocity ﬁeld are plotted for Mb = 10−

4 with the
mesh containing 20
40 cells. The solution obtained with the new approximation space dSh and
Roe’ numerical ﬂux seems to be in agreement with the exact incompressible solution. However,
the study of isolines is not suﬃcient to conclude on low Mach number behavior. A more detailed
study needs to be carried out not only on the mesh convergence of the compressible low Mach
number solution to the incompressible solution, but also on the order of magnitude of the density
gradient and the momentum divergence with respect to the Mach number.

×

In Figure 7, a mesh convergence at Mach number Mb = 10−

4 is performed. The steady solution
obtained with the new approximation space dSh and Roe’ numerical ﬂux converges towards the
incompressible exact solution with a rate close to one.

In Figure 8, we study the order of the dimensionless density gradient and the dimensionless

momentum divergence with respect to the Mach number. Then, we plot the semi-norms

˜ρh

|

|H 1

23

10−1100h10−1310−1110−910−710−510−310−1kp∞h−0k2slope=1.0Godunov,(dQ0)2Godunov,dShRusanov,dSh10−1100h10−1100k(u∞x)h−(u∞x)exk2slope=1.0Godunov,(dQ0)2Godunov,dShRusanov,dSh10−1100h10−1k(cid:0)u∞y(cid:1)h−(cid:0)u∞y(cid:1)exk2slope=1.0Godunov,(dQ0)2Godunov,dShRusanov,dShExact incompressible

Roe, (dQ0)2

Rusanov, dSh

Roe, dSh

Figure 6: Euler cylinder scattering - Iso-contours of the norm of the velocity obtained at Mach
4 are plotted.
number Mb = 10−

4. Twenty equally reparted contours between 8

6 and 3

10−

10−

×

×

24

Figure 7: Euler cylinder scattering - L2 norm of the error between the exact incompressible solution
and the numerical solution. Results are shown for a Mach number of Mb = 10−
4. A log-log plot is
used.

and

˜ρ˜uh

|

|H div , deﬁned by

ph
|

2
H 1 =
|

uh
|

2
H div =
|

C
c
∈
(cid:80)

ph

2 +
(cid:107)

(cid:90)c (cid:107)∇

S

Si(cid:90)S
∈
(cid:80)
(divxuh)2 +

[[ ph ]]2

[[ uh

nS ]]2

C
c
∈
(cid:80)

(cid:90)c

S

Si(cid:90)S
∈
(cid:80)

as a function of the Mach number for a ﬁxed mesh containing 10
imation space dSh and Roe’ numerical ﬂux, we observe that

(M ) and so

O

as expected.

˜ρ =

∇

O

M 2

and

(cid:0)

(cid:1)

divx(˜ρ˜u) =

O

·
20 cells. Using the new approx-
×
and
|H div =
|H 1 =
˜ρh
|
(M )

˜ρ˜uh
|

M 2

O

(cid:1)

(cid:0)

5.3.2 Propagation of a low Mach number acoustic wave over a steady vortex

In this test case, we evaluate the ability of the numerical scheme to handle both incompressible
phenomena and acoustic wave propagation at low Mach numbers. This test case was proposed
in [11], because it had been remarked that some low Mach number ﬁxes spoil the propagation of
acoustic waves. The domain is the rectangle [
[0; 1]. Periodic conditions are applied to
top and bottom boundaries, and input-output conditions are used for left and right boundaries.
The vortex is centered around (xc, yc) = (0.5, 0.5), characterized by its reference Mach number
Mref, and is given by

0.1; 1.1]

−

×

M 2
ref
λ2

max

2α2
1−r2

e−

,

ρ(x) = ρ0

1

(cid:18)

u(x) = u0




−
y
x

−

(cid:18)

λmaxr0

(cid:19)

(cid:19)
2α
1

−

α2
1−r2 ,

e−

r2


where ρ0 = 2, u0 = Mref
(x
λmax is chosen such that the maximal velocity norm is u0, i.e. λmax =
(cid:112)
(cid:112)

p(cid:48)(ρ0), α = 2, r0 = 0.45, r = r/r0, r =

(cid:0)

(cid:1)

25

yc)2 and
−
α2
max with
1−r2

xc)2 + (y

−
2αrmax
r2
1
max

−

e−

10−1100h10−910−810−710−610−510−4kρh−ρincompexk2slope=1.0slope=0.8Roe,(dQ0)2Rusanov,dShRoe,dSh10−1100h10−4k(ρux)h−(ρux)incompexk2slope=1.0Roe,(dQ0)2Rusanov,dShRoe,dSh10−1100h10−510−4k(ρuy)h−(ρuy)incompexk2slope=0.9Roe,(dQ0)2Rusanov,dShRoe,dShFigure 8: Euler cylinder scattering - Semi-norms
number for Mb = 10−

|H1 and
˜ρ
|
7. A log-log plot is used.

1 to Mb = 10−

˜ρ˜u
|

|Hdiv with respect to the Mach

α2 + √1 + α4. The low Mach number acoustic wave is centered in x = 0 and given for

rmax =
[
x
−

∈

−

0.05; 0.05] by
(cid:112)

ρ(x) = ρ0

(cid:16)
u(x) = u0 +

1 + Mref e1
−
2

1
1−x2

,

(cid:17)

p(cid:48) (ρ(x))




γ

1

−

p(cid:48) (ρ0)

,

−



where x = x/0.05.

(cid:16)(cid:112)
In Figure 9, the Mach number is plotted at diﬀerent time with a mesh containing 480

400
cells for a reference Mach number Mref = 10−
2. Using the new approximation space dSh with
Roe’ numerical ﬂux, the acoustic wave propagates correctly and the vortex is preserved over time.
Using dSh with Rusanov’ numerical ﬂux or (dQ0)2, the acoustic wave propagates correctly but the
vortex is quickly diﬀused over time.

(cid:112)

×

(cid:17)

6 Conclusion

In this article, the problem of low Mach number accuracy was addressed through the problem
of discrete curl preservation. We proposed to use a special approximation space for the velocity,
which is not the tensor product of the approximation space for scalars. This idea is partly inspired
by the Raviart-Thomas or Nédéléc ﬁnite element basis, which are vector ﬁnite element spaces,
but is diﬀerent because all the degrees of freedom are located within the cells. With this new
approximation space for velocities, we were able to design a curl, deﬁned in the adjoint sense,
that is preserved by the numerical scheme under mild assumption on the numerical ﬂux. Going
further, in the general quadrangular case, with boundary conditions, we were also able to prove the
existence and uniqueness of a discrete Hodge-Helmholtz decomposition, and also to prove that this
Hodge-Helmholtz decomposition is preserved by ﬁnite volume numerical scheme under the same
assumption on the numerical ﬂux. All these properties were conﬁrmed with numerical experiments
on the wave system.

The new approximation space for vectors was then used with the barotropic Euler system for
the momentum variable, and was thoroughly tested in diﬀerent low Mach number conﬁgurations.
All numerical tests show that this single change of approximation space provides a strong beneﬁt
for solving the low Mach number accuracy on quadrangular meshes.

This paper deals with the case of dimension two. The new approximation space deﬁned by
(1) for Cartesian mesh and by (31) for general quadrangular mesh can be extended to dimension
enriched by the two basis
S0 deﬁned by (5) must be replaced by Q3
three. For dimension three,
z + 1/2)T . All the results of the paper can be
y + 1/2, 0)T and (
functions (
x
−
extended to dimension three.
(cid:98)

−
Current investigations are focused on the higher order extension of this method, on extensions to
the problem of conservation of divergence constraints, and also on the three dimensional extension
of the scheme.

1/2, 0,

1/2,

−

−

x

(cid:98)

(cid:98)

(cid:98)

(cid:98)

0

26

10−710−610−510−410−310−210−1Machnumber10−1310−1110−910−710−510−3|˜ρh|H1slope=1.0slope=2.0Roe,(dQ0)2Rusanov,dShRoe,dSh10−710−610−510−410−310−210−1Machnumber10−710−610−510−410−310−210−1100|(˜ρ˜u)h|Hdivslope=1.0Roe,(dQ0)2Rusanov,dShRoe,dShInitial condition

Roe, (dQ0)2 , t = 0.25 s

Roe, (dQ0)2, t = 0.5 s

Roe, (dQ0)2, t = 3.5 s

Rusanov, dSh, t = 0.25 s

Rusanov, dSh, t = 0.5 s

Rusanov, dSh, t = 3.5 s

Roe, dSh, t = 0.25 s

Roe, dSh, t = 0.5 s

Roe, dSh, t = 3.5 s

Figure 9: Euler equations - Propagation of a low Mach acoutic wave over a steady vortex - Mach
number obtained at diﬀerent times for a reference Mach number Mref = 10−

2.

27

References

[1] Douglas N. Arnold, Daniele Boﬃ, and Richard S. Falk. Quadrilateral H (div) ﬁnite elements.

SIAM Journal on Numerical Analysis, 42(6):2429–2451, 2005.

[2] Dinshaw S Balsara, Roger Käppeli, Walter Boscheri, and Michael Dumbser. Curl constraint-
preserving reconstruction and the guidance it gives for mimetic scheme design. Communica-
tions on Applied Mathematics and Computation, pages 1–60, 2021.

[3] Wasilij Barsukow. Stationarity preserving schemes for multi-dimensional linear systems. Math-

ematics of Computation, 88(318):1621–1645, 2019.

[4] Wasilij Barsukow. Truly multi-dimensional all-speed schemes for the Euler equations on Carte-

sian grids. Journal of Computational Physics, 435:110216, 2021.

[5] Wasilij Barsukow, Pierre-Henri Maire, and Raphaël Loubère. A node-conservative vorticity-
preserving ﬁnite volume method for linear acoustics on unstructured grids. Math. of Comp.,
2023. Submitted.

[6] John B. Bell, Phillip Colella, and Harland M. Glaz. A second-order projection method for
the incompressible Navier-Stokes equations. Journal of computational physics, 85(2):257–283,
1989.

[7] Walter Boscheri, Giacomo Dimarco, and Lorenzo Pareschi. Locally structure-preserving div-
curl operators for high order discontinuous Galerkin schemes. Journal of Computational
Physics, 486:112130, 2023.

[8] Walter Boscheri, Michael Dumbser, Matteo Ioriatti, Ilya Peshkov, and Evgeniy Romenski. A
structure-preserving staggered semi-implicit ﬁnite volume scheme for continuum mechanics.
Journal of Computational Physics, 424:109866, 2021.

[9] Walter Boscheri, Raphaël Loubère, and Pierre-Henri Maire. An unconventional divergence
preserving ﬁnite-volume discretization of Lagrangian ideal MHD. Communications on Applied
Mathematics and Computation, pages 1–55, 2023.

[10] Jeremiah U. Brackbill and Daniel C. Barnes. The eﬀect of nonzero

B on the numerical so-
lution of the magnetohydrodynamic equations. Journal of Computational Physics, 35(3):426–
430, 1980.

∇ ·

[11] Pascal Bruel, Simon Delmas, Jonathan Jung, and Vincent Perrier. A low Mach correction
able to deal with low Mach acoustics. Journal of Computational Physics, 378:723–759, 2019.

[12] Zhiqiang Cai, Jim Douglas, and Xiu Ye. A stable nonconforming quadrilateral ﬁnite element
method for the stationary Stokes and Navier–Stokes equations. Calcolo, 36:215–232, 1999.

[13] Michel Crouzeix and P-A Raviart. Conforming and nonconforming ﬁnite element methods for
solving the stationary Stokes equations. Revue francaise d’automatique informatique recherche
opérationnelle. Mathématique, 7(R3):33–75, 1973.

[14] Andreas Dedner, Friedemann Kemm, Dietmar Kröner, C-D Munz, Thomas Schnitzer, and
Matthias Wesenberg. Hyperbolic divergence cleaning for the MHD equations. Journal of
Computational Physics, 175(2):645–673, 2002.

[15] Sarah Delcourte, Komla Domelevo, and Pascal Omnes. A discrete duality ﬁnite volume ap-
proach to hodge decomposition and div-curl problems on almost arbitrary two-dimensional
meshes. SIAM Journal on Numerical Analysis, 45(3):1142–1174, 2007.

[16] Stéphane Dellacherie, Pascal Omnes, and Felix Rieper. The inﬂuence of cell geometry on
the Godunov scheme applied to the linear wave equation. Journal of Computational Physics,
229(14):5315–5338, 2010.

[17] Michael Dumbser, Francesco Fambri, Elena Gaburro, and Anne Reinarz. On GLM curl clean-
ing for a ﬁrst order reduction of the CCZ4 formulation of the Einstein ﬁeld equations. Journal
of Computational Physics, 404:109088, 2020.

[18] Robert Eymard, Thierry Gallouët, Raphaele Herbin, and Jean-Claude Latché. Convergence of
the MAC scheme for the compressible Stokes equations. SIAM Journal on Numerical Analysis,
48(6):2218–2246, 2010.

28

[19] Robert Eymard, Thierry Gallouët, Raphaele Herbin, and Jean-Claude Latché. A convergent
ﬁnite element-ﬁnite volume scheme for the compressible Stokes problem. part II: the isentropic
case. Mathematics of Computation, 79(270):649–675, 2010.

[20] Thierry Gallouët, Raphaele Herbin, and Jean-Claude Latché. A convergent ﬁnite element-
ﬁnite volume scheme for the compressible Stokes problem. part I: The isothermal case. Math-
ematics of Computation, 78(267):1333–1352, 2009.

[21] Hervé Guillard. On the behavior of upwind schemes in the low Mach number limit. IV: P0
approximation on triangular and tetrahedral cells. Computers & Fluids, 38(10):1969–1972,
2009.

[22] Hervé Guillard and Boniface Nkonga. On the behaviour of upwind schemes in the low Mach

number limit: A review. Handbook of Numerical Analysis, 18:203–231, 2017.

[23] Francois Hermeline, Siham Layouni, and Pascal Omnes. A ﬁnite volume method for the
approximation of maxwell’s equations in two space dimensions on arbitrary meshes. Journal
of Computational Physics, 227(22):9365–9388, 2008.

[24] James M Hyman and Mikhail Shashkov. Natural discretizations for the divergence, gradient,
and curl on logically rectangular grids. Computers & Mathematics with Applications, 33(4):81–
104, 1997.

[25] Rolf Jeltsch and Manuel Torrilhon. On curl-preserving ﬁnite volume discretizations for shallow

water equations. BIT Numerical Mathematics, 46:35–53, 2006.

[26] Youngmok Jeon, Hyun Nam, Dongwoo Sheen, and Kwangshin Shim. A class of nonparametric
DSSY nonconforming quadrilateral elements. ESAIM: Mathematical Modelling and Numerical
Analysis, 47(6):1783–1796, 2013.

[27] Jonathan Jung and Vincent Perrier. Steady low Mach number ﬂows:

identiﬁcation of the
spurious mode and ﬁltering method. Journal of Computational Physics, page 111462, 2022.

[28] Jonathan Jung and Vincent Perrier. Long time behavior of ﬁnite volume discretization of
symmetrizable linear hyperbolic systems. IMA Journal of Numerical Analysis, 43(1):326–356,
2023.

[29] Jonathan Jung and Vincent Perrier. Behavior of the discontinuous Galerkin method for
compressible ﬂows at low Mach number on triangles and tetrahedrons. SIAM Journal on
Scientiﬁc Computing, 46(1):A452–A482, 2024.

[30] Vyacheslav Ivanovich Lebedev. Diﬀerence analogues of orthogonal decompositions, basic dif-
ferential operators and some boundary problems of mathematical physics. I. USSR Compu-
tational Mathematics and Mathematical Physics, 4(3):69–92, 1964.

[31] Youai Li. A new family of nonconforming ﬁnite elements on quadrilaterals. Computers &

Mathematics with Applications, 70(4):637–647, 2015.

[32] Konstantin Lipnikov, Gianmarco Manzini, and Mikhail Shashkov. Mimetic ﬁnite diﬀerence

method. Journal of Computational Physics, 257:1163–1227, 2014.

[33] Zhaoliang Meng, Jintao Cui, and Zhongxuan Luo. A new rotated nonconforming quadrilateral

element. Journal of Scientiﬁc Computing, 74:324–335, 2018.

[34] Claus-Dieter Munz, Pascal Omnes, Rudolf Schneider, Éric Sonnendrücker, and Ursula Voss.
Divergence correction techniques for Maxwell solvers based on a hyperbolic model. Journal
of Computational Physics, 161(2):484–511, 2000.

[35] Roy A Nicolaides. Analysis and convergence of the MAC scheme. I. the linear problem. SIAM

Journal on Numerical Analysis, 29(6):1579–1591, 1992.

[36] Roy A. Nicolaides and X. Wu. Analysis and convergence of the MAC scheme. II. Navier-Stokes

equations. Mathematics of Computation, 65(213):29–44, 1996.

[37] Vincent Perrier. Development of discontinuous Galerkin methods for hyperbolic systems that

preserve a curl or a divergence constraint. Submitted, May 2024.

29

[38] Vincent Perrier. discrete de-Rham complex involving a discontinuous ﬁnite element space for
velocities: the case of periodic straight triangular and Cartesian meshes. Submitted, April
2024.

[39] Rolf Rannacher and Stefan Turek. Simple nonconforming quadrilateral Stokes element. Nu-

merical Methods for Partial Diﬀerential Equations, 8(2):97–111, 1992.

[40] Pierre-Arnaud Raviart and Jean-Marie Thomas. A mixed ﬁnite element method for 2-nd
order elliptic problems. In Mathematical Aspects of Finite Element Methods: Proceedings of
the Conference Held in Rome, December 10–12, 1975, pages 292–315. Springer, 2006.

[41] Maurizio Tavelli and Michael Dumbser. A pressure-based semi-implicit space–time discontin-
uous Galerkin method on staggered unstructured meshes for the solution of the compressible
Navier–Stokes equations at all Mach numbers. Journal of Computational Physics, 341:341–
376, 2017.

[42] Eleuterio F. Toro. Riemann solvers and numerical methods for ﬂuid dynamics. Springer-

Verlag, Berlin, third edition, 2009. A practical introduction.

[43] Manuel Torrilhon. Locally divergence-preserving upwind ﬁnite volume schemes for magneto-
hydrodynamic equations. SIAM Journal on Scientiﬁc Computing, 26(4):1166–1191, 2005.

[44] Manuel Torrilhon and Michael Fey. Constraint-preserving upwind methods for multidimen-
sional advection equations. SIAM journal on numerical analysis, 42(4):1694–1728, 2004.

30

