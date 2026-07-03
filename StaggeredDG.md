| Numerical | analysis | and simulation |        | of staggered | schemes | for |
| --------- | -------- | -------------- | ------ | ------------ | ------- | --- |
|           |          | low Mach       | number | flows        |         |     |
Esteban Coiffier
| To cite this | version: |     |     |     |     |     |
| ------------ | -------- | --- | --- | --- | --- | --- |
EstebanCoiffier. NumericalanalysisandsimulationofstaggeredschemesforlowMachnumberflows. Nu-
mericalAnalysis[math.NA].UniversitédePauetdesPaysdel’Adour,2025. English. ⟨NNT:2025PAUU3031⟩.
⟨tel-05528307⟩
|     |     | HAL Id: | tel-05528307 |     |     |     |
| --- | --- | ------- | ------------ | --- | --- | --- |
https://theses.hal.science/tel-05528307v1
Submittedon26Feb2026
HAL is a multi-disciplinary open access archive L’archiveouvertepluridisciplinaireHAL,estdes-
for the deposit and dissemination of scientific re- tinée au dépôt et à la diffusion de documents scien-
searchdocuments,whethertheyarepublishedornot. tifiquesdeniveaurecherche,publiésounon,émanant
Thedocumentsmaycomefromteachingandresearch des établissements d’enseignement et de recherche
institutionsinFranceorabroad,orfrompublicorpri- français ou étrangers, des laboratoires publics ou
| vateresearchcenters. |     |     | privés. |     |     |     |
| -------------------- | --- | --- | ------- | --- | --- | --- |
HALAuthorization

|           |     | Universite´ |          | de               | Pau et       | des        | Pays     |              | de  | l’Adour      |     |
| --------- | --- | ----------- | -------- | ---------------- | ------------ | ---------- | -------- | ------------ | --- | ------------ | --- |
|           |     |             |          | the`se           | de           | doctorat   |          |              |     |              |     |
|           |     |             | ED       | 211 - Sciences   | exactes      | et         | leurs    | applications |     |              |     |
|           |     |             |          | Th`ese           | pour obtenir |            | le grade | de           |     |              |     |
|           |     | Docteur     |          | de l’Universit´e | de           | Pau        | et des   | Pays         | de  | l’Adour      |     |
|           |     |             |          | Math´ematiques   | appliqu´ees  | –          | analyse  | num´erique   |     |              |     |
| Numerical |     |             | analysis |                  | and          | simulation |          |              |     | of staggered |     |
|           |     | schemes     |          | for              | low Mach     |            |          | number       |     | flows        |     |
Analyse et simulation num´eriques de sch´emas d´ecal´es pour les ´ecoulements `a
|             |             |     |      | bas | nombre | de  | Mach |              |     |                      |      |
| ----------- | ----------- | --- | ---- | --- | ------ | --- | ---- | ------------ | --- | -------------------- | ---- |
|             |             |     |      |     |        |     |      |              |     | Vincent PERRIER      |      |
| Pr´esent´ee | et soutenue |     | le : |     |        |     |      | Directeur    | de  | th`ese :             |      |
| 8 d´ecembre | 2025        |     |      |     |        |     |      | Co-directeur |     | de th`ese : Jonathan | JUNG |
Par : Esteban COIFFIER Encadrant/Examinateur : Michael NDJINGA
|     |     |     |     |     | Jury | :   |     |     |     |     |     |
| --- | --- | --- | --- | --- | ---- | --- | --- | --- | --- | --- | --- |
Christophe CHALONS Pr, Universit´e Versailles Saint-Quentin-en-Yvelines Rapporteur
Nicolas SEGUIN Dr INRIA, Universit´e de Montpellier Rapporteur
Rapha`ele HERBIN Pr, Universit´e d’Aix-Marseille Examinatrice
| Alexiane |         | PLESSIER |     | CEA Saclay |     |     |     |     |     | Examinatrice |     |
| -------- | ------- | -------- | --- | ---------- | --- | --- | --- | --- | --- | ------------ | --- |
| Michael  | NDJINGA |          |     |            |     |     |     |     |     | Examinateur  |     |
|          |         |          |     | CEA Saclay |     |     |     |     |     |              |     |
Jonathan JUNG Mcf, Universit´e de Pau et des Pays de l’Adour Co-directeur
Vincent PERRIER Dr INRIA, Universit´e de Pau et des Pays de l’Adour Directeur
`
Rapha¨el LOUBERE Dr CNRS, Universit´e de Bordeaux Pr´esident du Jury
Service de Thermohydraulique et de M´ecanique de Fluides, CEA Saclay/ISAS/DM2S
E´quipe-projet
Laboratoire de Math´ematiques et de leurs applications de Pau UMR CNRS 5142 inria Cagire

2
Abstract.
This thesis focuses on the simulation of flows in nuclear reactors with the aim of improving
safety conditions. The objective of a numerical simulation is, by definition, to faithfully preserve
the properties of the underlying physical model. Since, under nominal operating conditions,
the flows in such reactors are considered to be at low Mach number, it is essential to have a
numerical tool capable of accurately capturing this flow regime. On the other hand, the equations
of fluid mechanics, which govern these flows, fundamentally derive from conservation principles.
Replicating these conservation properties at the discrete level is possible, notably through the
widely used finite volume methods. However, it is well known that such methods generally lack
accuracy in the low-Mach-number regime. There is, nevertheless, a particular case in which this
limitation disappears: on simplicial meshes, these methods can recover this desirable accuracy
without any correction. A possible interpretation of this unexpected property is that, in this
specific case, a Hodge–Helmholtz decomposition exists on the discrete velocity space. Furthermore,
there exist so-called staggered methods, which appear to be ideal candidates for preserving this
structure on more general meshes, although the notion of conservation is not well defined for them.
Thus, our approach is twofold: on the one hand, we aim to extend the conservation properties of
classical schemes to the staggered framework; on the other hand, we introduce a methodology for
constructing numerical schemes that are accurate at low Mach numbers, which we then apply to a
staggered discretization.
In practice, a well-established link connects the low-Mach-number limit with the long-time limit of
a wave system, through an asymptotic analysis with respect to the Mach number. This connection
allows us to reformulate the problem of low-Mach-number accuracy as the simpler problem of
capturing the long-time limit of a linear system. At the continuous level, this limit is characterized
by a Hodge–Helmholtz decomposition, itself arising from a more general framework: the de Rham
complexes. The first step therefore consists in introducing a staggered discretization of the wave
system based on a discrete de Rham complex. Although this structure may seem, at first glance,
unrelated to standard numerical analysis issues, it enables the preservation, at the discrete level, of
fundamental properties of partial differential equations. We show in particular that it leads to the
identification of a long-time limit corresponding to a specific discrete Hodge–Helmholtz decomposi-
tion. This discretization is then extended to the barotropic Euler system. Conservation is achieved
by working directly on this system, while low-Mach-number accuracy follows from the long-time
consistency of the proposed staggered scheme on the wave system. Finally, implementations on the
C++ platform SolverLab illustrate numerically the validity of these results for both the wave
system and the barotropic Euler system.
Keywords: Compressible Euler system ; Linear acoustic wave system ; Staggered
schemes ; Low Mach number flows ; Discrete de Rham complexes ; Hodge-Helmholtz
Decomposition .

3
R´esum´e.
Cette th`ese s’int´eresse `a la simulation d’´ecoulements dans des r´eacteurs nucl´eaires dans un but
d’am´eliorer les conditions de suˆret´e. L’objectif d’une simulation est, par d´efinition, de pr´eserver
fid`element les propri´et´es du mod`ele sous-jacent. De fait, puisqu’en r´egime nominal de ces r´eacteurs
les ´ecoulements sont consid´er´es comme ´etant `a bas nombre de Mach, il est essentiel d’avoir un
outil num´erique capable de capturer convenablement ce type d’´ecoulements. D’un autre cˆot´e, les
´equations de la m´ecanique des fluides, qui r´egissent ces ´ecoulements, d´erivent fondamentalement de
principes de conservation. R´epliquer cette propri´et´e de conservation num´eriquement est possible,
notamment grˆace aux tr`es populaires m´ethodes de volumes finis. Cependant, il est bien connu
que ces derni`eres ne sont g´en´eralement pas pr´ecises `a bas nombre de Mach. Il existe toutefois un
cas particulier ou` cette limitation disparaˆıt : sur des maillages simpliciaux, ces derni`eres peuvent
r´ecup´erer sans correction cette pr´ecieuse pr´ecision. Une interpr´etation possible de cette propri´et´e
inattendue est qu’il existe dans ce cas particulier une d´ecomposition de Hodge-Helmholtz sur
l’espace des vitesses discr`etes. Par ailleurs, il existe des m´ethodes, dites d´ecal´ees, qui semblent
ˆetre les candidates id´eales pour pr´eserver sur des maillages plus vari´es cette structure, mais pour
lesquelles la notion de conservation est mal d´efinie.
Ainsi, notre proposition est double : d’une part nous cherchons `a ´etendre les propri´et´es de
conservation de sch´emas classiques au formalisme d´ecal´e, d’autre part, nous introduisons une
m´ethodologie permettant d’obtenir des sch´emas num´eriques pr´ecis `a bas nombre de Mach que nous
appliquons au cas d’une discr´etisation d´ecal´ee.
En pratique, un lien bien ´etabli relie la limite `a bas nombre de Mach et la limite en temps long
d’un syst`eme d’ondes, via une analyse asymptotique en nombre de Mach. Ce lien permet de
ramener le probl`eme de la pr´ecision `a bas nombre de Mach `a celui, plus simple, de la capture de
la limite en temps long d’un syst`eme lin´eaire. Cette limite, au niveau continu, est caract´eris´ee par
une d´ecomposition de Hodge–Helmholtz, elle-mˆeme issue d’un cadre plus g´en´eral : les complexes
de de Rham. La premi`ere ´etape consiste `a introduire une discr´etisation d´ecal´ee du syst`eme
d’ondes fond´ee sur un complexe de de Rham discret. Cette structure, bien qu’a priori ´eloign´ee
des probl´ematiques classiques d’analyse num´erique, permet de pr´eserver au niveau discret des
propri´et´es fondamentales des ´equations aux d´eriv´ees partielles. Nous montrons qu’elle conduit `a
l’identification d’une limite en temps long correspondant `a une d´ecomposition de Hodge–Helmholtz
discr`ete particuli`ere. Cette discr´etisation est ensuite ´etendue au syst`eme d’Euler barotrope. La
conservation est obtenue en travaillant directement sur ce syst`eme, tandis que la pr´ecision `a
bas nombre de Mach r´esulte de la consistance en temps long du sch´ema d´ecal´e sur le syst`eme
d’ondes. Enfin, des impl´ementations sur la plateforme C++ SolverLab permettent d’illustrer
num´eriquementlavalidit´edespropositionspourlesyst`emed’ondesetlesyst`emed’Eulerbarotrope.
Mots cl´es : Syst`eme d’Euler compressible ; Syst`eme des ondes lin´eaire ; Sch´emas
d´ecal´es ; E´coulements `a bas nombre de Mach ; Complexes de de Rham discrets ;
D´ecomposition de Hodge-Helmholtz.

5
Remerciements.
Unmanuscritdeth`esen’enestunqueparcequ’ila´et´e’valid´eparlespairs’,etc’estpourquoi
je tiens a` remercier en premier lieu les membres du jury, et en particulier Christophe Chalons
| et Nicolas | Seguin, | qui | ont accept´e | de  | rapporter | cette | th`ese. |     |
| ---------- | ------- | --- | ------------ | --- | --------- | ----- | ------- | --- |
Je remercie ´egalement Rapha`ele Herbin, Rapha¨el Loub`ere et Alexiane Plessier, pour leurs
questions, leurs remarques pertinentes et l’int´erˆet qu’ils ont port´e `a ce travail. C’est un
| honneur | de vous | avoir | tous | compt´es | dans mon | jury. |     |     |
| ------- | ------- | ----- | ---- | -------- | -------- | ----- | --- | --- |
Maintenant, une th`ese ne se conduit pas seul et je tiens donc ´evidemment `a exprimer ma
profonde gratitude `a Vincent, Jonathan et Michael pour leur encadrement tout au long de
ces trois ann´ees. Jonathan, merci encore pour ta disponibilit´e, ton ouverture et ta sagacit´e.
Vincent, je te remercie pour ta rigueur et ton souci du d´etail. Vous avez port´e ce manuscrit `a
un niveau sup´erieur, j’ai beaucoup appris de vous. Michael, merci de m’avoir toujours pouss´e
`a d´ecouvrir plus. Je vous suis sinc`erement et absolument reconnaissant de cette opportunit´e.
Pour toutes les conversations brillantes (et parfois compl`etement stupides) je remercie
Corentin, Th´eotime, Matthieu (que j’accuserai d’ˆetre les principaux instigateurs de ce second
point) mais aussi Quentin, Damien, Lidija, Mathis, Pierre-Lo¨ıc, Ayoub, Thomas, et tout
les membres du LMEC de Didier (puis de Sandrine). Merci pour tous ces bons moments,
| dommage | que | nous ayons | eu  | a travaill´e | entre | ces pauses | caf´es. |     |
| ------- | --- | ---------- | --- | ------------ | ----- | ---------- | ------- | --- |
Je remercie aussi chaleureusement les membres du LDEL de Julie pour leur accueil.
E´douard, E´tienne, E´lie,
J’adresse une pens´ee particuli`ere `a Capucines, Andrew, Cl´ement,
Gr´eg, Antonin avec qui nous avons r´eussi `a ne jamais parler boulot (vraiment ?) mais aussi
tout les membres du LDEL que j’ai pu cˆotoyer. Vous avez r´eussi `a rendre ma fin de th`ese
| (p´eriode | notoirement | d´esagr´eable) |     | moins | p´enible. | Merci | pour | ¸ca. |
| --------- | ----------- | -------------- | --- | ----- | --------- | ----- | ---- | ---- |
Enfin, je remercie mes parents pour leur soutien constant tout au long de ces ann´ees, sans
lequel cette th`ese et bien d’autres choses (ma naissance) n’auraient pas ´et´e possibles. Merci `a
| mon fr`ere | (d´esol´e | pour | l’avion, | c’´etait | moi). | Merci `a | Chlo´e. |     |
| ---------- | --------- | ---- | -------- | -------- | ----- | -------- | ------- | --- |
P.S : Ces quelques lignes m’ont pris autant de temps `a r´ediger que le reste du manuscrit,
| soyez indulgents. |     | Ceux | que | j’ai oubli´es, | je pense | `a vous | ici. |     |
| ----------------- | --- | ---- | --- | -------------- | -------- | ------- | ---- | --- |

Contents
1 Introduction 10
1.1 Context . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 11
1.2 Flow modelling: the Euler equations . . . . . . . . . . . . . . . . . . . . . . . . 13
1.2.1 The Euler equations as a system of conservation laws . . . . . . . . . . 14
1.2.2 The dimensionless Euler equations . . . . . . . . . . . . . . . . . . . . . 16
1.3 Low Mach number flows simulation: an overview . . . . . . . . . . . . . . . . . 17
1.3.1 Collocated schemes . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 19
1.3.2 Staggered schemes . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 22
1.4 Plan of the manuscript . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 25
2 An approach based on asymptotic expansions to simplify the study of the
low Mach number limit 27
2.1 Introduction . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 27
2.2 Link between the low Mach number limit and the long time limit of a wave
system through low Mach number asymptotic expansion . . . . . . . . . . . . . 28
2.3 Long time limit of the wave system . . . . . . . . . . . . . . . . . . . . . . . . . 32
2.4 Conclusion . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 34
3 Study of the stability of the staggered schemes for the one dimensional wave
system 36
3.1 Introduction . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 36
3.2 The staggered schemes on the one dimensional wave system . . . . . . . . . . . 39
3.2.1 Staggered schemes . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 39
3.2.2 Finite volume interpretation of the staggered scheme . . . . . . . . . . . 40
3.3 von-Neumann analysis and energy dissipation . . . . . . . . . . . . . . . . . . . 42
3.3.1 von Neumann stability analysis . . . . . . . . . . . . . . . . . . . . . . . 45
3.3.2 Energy dissipation . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 56
3.4 l -stability on the characteristic variables . . . . . . . . . . . . . . . . . . . . . 63
∞
3.4.1 Characteristic variables defined on half cells . . . . . . . . . . . . . . . . 63
3.4.2 Characteristic variables defined on primal cells . . . . . . . . . . . . . . 65
3.5 Discussion on some preexisting staggered schemes through low Mach number
asymptotics . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 66
3.6 Numerical results . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 67
7

CONTENTS 8
3.6.1 Numerical study of the amplification matrices . . . . . . . . . . . . . . . 68
3.6.2 Tests on the numerical schemes . . . . . . . . . . . . . . . . . . . . . . . 69
3.7 Conclusion . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 70
4 Hodge-Helmholtz decomposition and de Rham complexes: continuous and
discrete aspects 79
4.1 Introduction . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 80
4.2 Continuous de Rham complexes and Hodge-Helmholtz decomposition . . . . . 80
4.2.1 A discussion on the case of the 2D de Rham complex. . . . . . . . . . . 80
4.2.2 The de Rham formalism and harmonic forms . . . . . . . . . . . . . . . 83
4.3 Discrete de Rham complexes . . . . . . . . . . . . . . . . . . . . . . . . . . . . 85
4.4 The N´ed´elec-Raviart-Thomas de Rham staggered approximation space . . . . . 87
4.4.1 Mesh and notation . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 87
4.4.2 The Raviart-Thomas space . . . . . . . . . . . . . . . . . . . . . . . . . 87
4.4.3 Properties of the discrete complex . . . . . . . . . . . . . . . . . . . . . 91
4.5 Conclusion . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 102
5 Development of a class of long time consistent staggered schemes on the
wave system 104
5.1 Introduction . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 105
5.2 The Raviart-Thomas staggered scheme for the multi-dimensional wave system . 107
5.2.1 Deriving the ’centred’ finite volume scheme . . . . . . . . . . . . . . . . 107
5.2.2 Adding the appropriate, curl-preserving, diffusion . . . . . . . . . . . . . 112
5.3 L2 Stability analysis . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 115
−
5.3.1 Tools for stability. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 115
5.3.2 Stability results . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 119
5.3.3 Proof of Proposition 5.3.1 . . . . . . . . . . . . . . . . . . . . . . . . . . 121
5.3.4 Proof of Proposition 5.3.3 . . . . . . . . . . . . . . . . . . . . . . . . . . 126
5.4 Discrete long time behaviour . . . . . . . . . . . . . . . . . . . . . . . . . . . . 127
5.5 Discussion on some preexisting staggered schemes through low Mach number
asymptotics . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 136
5.6 Numerical Results . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 138
5.6.1 Numerical long time . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 139
5.6.2 On the necessity of a stationary preserving diffusion . . . . . . . . . . . 142
5.7 Conclusion . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 147
6 Extension of the Raviart-Thomas staggered scheme to compressible baro-
tropic flows 150
6.1 Introduction . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 150
6.2 The Raviart-Thomas staggered scheme for the two-dimensional Euler barotropic
equations . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 152
6.2.1 Deriving the ’centred’ scheme . . . . . . . . . . . . . . . . . . . . . . . . 152
6.2.2 Adding the appropriate diffusion . . . . . . . . . . . . . . . . . . . . . . 155

CONTENTS 9
6.2.3 Numerical treatment of the boundary conditions . . . . . . . . . . . . . 158
6.3 Conservation . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 160
6.4 Low Mach number analysis . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 162
6.5 Discussion on some preexisting staggered schemes . . . . . . . . . . . . . . . . . 165
6.5.1 Low Mach number behaviour . . . . . . . . . . . . . . . . . . . . . . . . 165
6.5.2 Discrete entropy dissipation . . . . . . . . . . . . . . . . . . . . . . . . 166
6.5.3 Other computational questions . . . . . . . . . . . . . . . . . . . . . . . 167
6.6 Numerical results . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 168
6.6.1 1D Riemann problems . . . . . . . . . . . . . . . . . . . . . . . . . . . . 168
6.6.2 Cylinder scattering . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 169
6.6.3 Propagation of a low Mach number acoustic wave through a stationary
low Mach number vortex . . . . . . . . . . . . . . . . . . . . . . . . . . 170
6.7 Conclusion . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 172
7 Conclusion 179
A Proof of Energy Dissipation in ImEx 183
B Computation of the lumped rotated gradient 192
C Low Mach number asymptotic analysis of some staggered schemes 196
C.1 First example: the explicit scheme of Duran .A, Vila. J-P and Baraille . R . . 196
C.2 Second example: the implicit scheme of Herbin .R, Kheriji. W and Latch´e .J-C 202
D Entropy for Explicit Crouzeix-Raviart staggered discretizations using ∇div
stabilization 207
D.1 The numerical scheme . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 209
D.2 Properties of the momentum equation . . . . . . . . . . . . . . . . . . . . . . . 212
D.3 Properties of the mass equation . . . . . . . . . . . . . . . . . . . . . . . . . . . 215
E Solution of the 1d Riemann Problem in the barotropic case 224
E.1 The Riemann problem in one dimension . . . . . . . . . . . . . . . . . . . . . . 225
E.1.1 Shocks . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 226
E.1.2 Rarefaction waves . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 227
E.1.3 The Newton solver . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 227
F Computation of the convection term for the Raviart-Thomas staggered
scheme 231

Chapter 1
Introduction
Contents
1.1 Context . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 11
1.2 Flow modelling: the Euler equations . . . . . . . . . . . . . . . . . . 13
1.2.1 The Euler equations as a system of conservation laws . . . . . . . . . 14
1.2.2 The dimensionless Euler equations . . . . . . . . . . . . . . . . . . . . 16
1.3 Low Mach number flows simulation: an overview . . . . . . . . . . 17
1.3.1 Collocated schemes . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 19
1.3.2 Staggered schemes . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 22
1.4 Plan of the manuscript . . . . . . . . . . . . . . . . . . . . . . . . . . 25
10

CHAPTER 1. INTRODUCTION 11
1.1 Context
This thesis has been conducted in the French Atomic and Alternative Energies Commission
(denoted CEA in French acronym), which is an institute at the intersection between academic
research and industrial questions. At CEA Saclay, the Thermohydraulics and Fluid Mechanics
Unit is focused on the development of simulation tools for fluid dynamics, with the aim to
applythesetoolstotheconceptionandthecontinuoussearchforsafetyimprovementofnuclear
reactors.
Asoftenwhendealingwithnuclearenergyproduction,thegoalofanuclearpowerplantisto
extract thermal energy in order to produce electricity. Using fissile matter, generally enriched
uranium isotopes, a nuclear chain reaction is provoked and controlled in order to free heat
with fissions, or ’breaking heavy atoms’; the heat will be then transported by the water of the
primary circuit (see Figure 1.1). For contamination reasons, the primary circuit is separated
from a secondary circuit in which a steam generator uses the heat of the primary circuit to
generate steam. In return, this steam impulses the rotation of an alternator: through this,
mechanical energy from a turbine is transformed into electricity. Actually, 80 % of the French
Figure 1.1: Schematic representation of PWR nuclear plant (Source: radioactivity.eu.com)
nuclear stock is composed of Pressurized Water Reactors (PWRs) and one of the particularities
of this technology is, as the name suggests it, that the pressure in the primary circuit is kept
purposefully around 150 times the atmospheric pressure. In such conditions, water neighboring
320°C is still liquid. Fittingly, under liquid form, water has a much better calorific capacity
thanundergaseousform, sosincetheheatgeneratedbythefissionwillbringthewateratthese

CHAPTER 1. INTRODUCTION 12
temperatures, it is particularly interesting to maintain water in liquid state in order to increase
energy efficiency.
Then,becauseoftheaforementionedphysicalconditions,inthecaseofnominalfunctioning,
the flow in such reactors is considered to be in a low Mach number regime, in the sense
that the sound propagation speed is of multiple orders larger than the fluid velocity.
In parallel, when looking at the simulation of flows in the primary circuit, it is industrially
convenient to model the uranium rods and circuits with a porosity that will account for the
geometry. In this case, when the flow will encounter, for example, a rod, a porosity jump
will represent this obstacle in the model. This kind of modelling, if not treated carefully
can introduce spurious oscillations and nonphysical mass and pressure losses. Industrially,
thesepurelynumerical phenomenaareprohibitivebecausetheyleadtoaconservationproperty
breakdown. The deterioration of this conservation property is antinomic with the continuous
models,inwhichitisakeyattribute. Hence,itisprimordialtodefineproperlyconservation
at the level of our numerical approximation.
In this context, some remarks are of particular interest :
• In hydrodynamics simulation, it is very natural to use finite volumes methods which
mimic trivially (at the numerical scale) the conservative nature of the models. As a
natural consequence, they are good candidates to treat conservation related issues.
• These methods however, are infamously inefficient and/or imprecise in this low Mach
number/subsonicregimeonquadrangularmeshes,whicharewidelypopularinthereactor
design field.
• By opposition, staggered methods are known for their accurate approximation of the
acoustic related operators that are sensitive in the subsonic regime.
• Staggered schemes are originally coming from the fully incompressible community and
necessitate adaptations to weakly compressible flows, in particular with respect to mim-
icking conservation properties.
Generally, a numerical scheme that is able to tackle low Mach number flows in a conservative
manner is still considered as a holy Grail in computational fluid mechanics. Since experiment-
ation at the same scale as a reactor is very costly, simulation is still a fundamental tool in
the conception of nuclear vessel; in the special case of Pressurized Water Reactors, it is thus
essential to be able to mimic in the simulations the conditions in which the vessel runs usu-
ally. As a consequence a numerical method that is both low-Mach-number accurate and
mathematically conservative, is needed.
Weaimatdevelopingaclassofstaggeredschemesthatanswerstothisproblematicsituation.
For this purpose, we will present the model of interest in section 1.2 then we will introduce
multiple approaches for the approximation of this model in the low Mach number regime in
section 1.3. This section will enable us to distill ingredients that we interpret as important for
low Mach number accuracy; so that in section 1.4 we will present our methodology to obtain
low-Mach-number accurate and conservative staggered schemes.

| CHAPTER  | 1. INTRODUCTION |     |     |           |           | 13  |
| -------- | --------------- | --- | --- | --------- | --------- | --- |
| 1.2 Flow | modelling:      |     |     | the Euler | equations |     |
The primary circuit of the vessel contains water that overwhelms the combustible (uranium)
rods (see, Figure 1.2). It is said to act as a moderator and a coolant: on the one hand it helps
controlling the nuclear reaction, and on the other hand the water transports/diffuses the heat
generated by the nuclear fissions to the secondary circuit, where the heat is extracted from.
A nuclear reactor and particularly the core, is a highly multi-physical object combined with
| rapidly retro-active |     | phenomena:   |     |           |                            |     |
| -------------------- | --- | ------------ | --- | --------- | -------------------------- | --- |
| • fluid-structure    |     | interactions |     | occurring | at very high temperatures, |     |
•
neutrons induce fluid density variations which in return impact the neutronic’spectrum
| of         | absorption,  | fission | (called | Doppler        | effect) etc, |     |
| ---------- | ------------ | ------- | ------- | -------------- | ------------ | --- |
| • chemical | interactions |         | caused  | by neutronical | effects,     |     |
•
| and | the list | goes on... |     |     |     |     |
| --- | -------- | ---------- | --- | --- | --- | --- |
Figure 1.2: Schematic representation of PWR core (Source: energyencyclopedia.com)
Yet, even when ignoring neutronical effects and phenomena linked to the mechanical structure,
a complete modelling of the flow in the primary circuit of a nuclear vessel would need to take
into account, in itself, at least (see [1] for an extensive review of the thermohydraulics of a
nuclear reactor):

| CHAPTER | 1. INTRODUCTION |     |          |     |     |     |     |     |     | 14  |
| ------- | --------------- | --- | -------- | --- | --- | --- | --- | --- | --- | --- |
| i) weak | compressibility |     | effects, |     |     |     |     |     |     |     |
ii) turbulence induced by thermal and viscous effects near the combustible and the walls,
iii) biphasic phenomena, as small thermodynamical changes near the combustible rods will
| lead         | to bubbles | nucleation, |     |              |     |     |     |     |     |     |
| ------------ | ---------- | ----------- | --- | ------------ | --- | --- | --- | --- | --- | --- |
| iv) external | forces     | such        | as  | the gravity. |     |     |     |     |     |     |
We are far from numerical tools that are able to treat the whole coupled chain of computations,
we will thus have to humbly restrict to the conception of a numerical tool dedicated to the fluid
mechanics only. ii) and iii) are out of the scope of this study and we will neglect the gravity
in order to focus on i) which defines the behaviour of the flow at low Mach number.
With this in mind we can restrict our analysis to the equations for monophasic inviscid
flows without source terms given by the seminal Euler system. Suppose that thermodynamical
Rd,
equilibriumisachieved,andthattheentropyisconstant,thenaflowinanopendomainΩ
⊂
d 1,2,3
|     | is fully | described |     | with: |     |     |     |     |     |     |
| --- | -------- | --------- | --- | ----- | --- | --- | --- | --- | --- | --- |
| ∈ { | }        |           |     |       |     |     |     |     |     |     |
• a density ρ : Ω R, accounting for information on the mass and matter,
→
| •   |     |     | Rd, |     |     |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
a velocity u : Ω carrying information on the kinetic variation of the fluids.
→
The Euler isentropic/barotropic system is therefore derived using first physical principles such
as: the conservation of mass (1st equation), Newton’s 2nd law which leads to the equation of
| conservation | of the | momentum |             | (2nd | equation). |     |               |     |          |     |
| ------------ | ------ | -------- | ----------- | ---- | ---------- | --- | ------------- | --- | -------- | --- |
|              |        |          | ∂ ρ+div(ρu) |      | = 0        |     | (conservation |     | of mass) |     |
t
(1.1)

|     |     | ∂ (ρu)+div(ρu |     |     | u)+∇p | = 0 | (Newton’s |     | 2nd law). |     |
| --- | --- | ------------- | --- | --- | ----- | --- | --------- | --- | --------- | --- |
 t
⊗
isclosed
The system with an equation of state p := p(ρ), the pressure is a function of the
density.
The system (1.1) can be rewritten under the form of a system of conservation laws.
| 1.2.1 | The Euler | equations |     |     | as a system | of  | conservation |     | laws |     |
| ----- | --------- | --------- | --- | --- | ----------- | --- | ------------ | --- | ---- | --- |
A conservation law is a physical principle stating that an isolated system will preserve in time
a set of measurable quantities. In the case of the Euler system, because the viscosity of the
fluid and the friction effects on walls are neglected, the system is isolated and only interacting
with itself.
Thiscanbetranslatedinmathematicaltermswiththefollowingformulation; letU = (ρ,ρu)t a
new vectorial unknown, then we can rewrite the Euler system under the form of a conservation
law :
|     |     |     |           |     |      |       | Rd   |     | Rk, |       |
| --- | --- | --- | --------- | --- | ---- | ----- | ---- | --- | --- | ----- |
|     |     | ∂   | U+divF(U) |     | = 0, | U : Ω | [0,+ | [   |     | (1.2) |
|     |     |     | t         |     |      | ⊂     | ×    | ∞   | −→  |       |

| CHAPTER   | 1. INTRODUCTION |     |     |     |     |     | 15  |
| --------- | --------------- | --- | --- | --- | --- | --- | --- |
| with here | k = d+1 and     |     |     |     |     |     |     |
ρut
|     |     | F(U) | =   | .     |     |     |     |
| --- | --- | ---- | --- | ----- | --- | --- | --- |
|     |     |      | pI  | +ρu u |     |     |     |
d
|     |     |     | (cid:18) | ⊗ (cid:19) |     |     |     |
| --- | --- | --- | -------- | ---------- | --- | --- | --- |
The divergence operator div acts row-wise and is responsible for the conservation law formu-
lation;
div(ρu)
|     |     | divF(U) | =   |     | Rd+1. |     |     |
| --- | --- | ------- | --- | --- | ----- | --- | --- |
u+pI
|     |     |     | div(ρu     | )          | ∈   |     |     |
| --- | --- | --- | ---------- | ---------- | --- | --- | --- |
|     |     |     | (cid:18) ⊗ | d (cid:19) |     |     |     |
Under the assumption that p(ρ) > 0, this system carries information with finite-speed waves.
(cid:48)
This characteristic makes it a hyperbolic system, which is formalized in the following defin-
ition
Definition 1.2.1 (Hyperbolicity, see [2, Introduction p 2] ). A system of conservation laws of
| the form | is said to | be hyperbolic | if  |     |     |     |     |
| -------- | ---------- | ------------- | --- | --- | --- | --- | --- |
(1.2)
n Rd, U Rk, the matrix ∇(F(U)n) is diagonalizable with real eigenvalues
| ∀   | ∈ ∀ ∈ |     |     |     |     |     |     |
| --- | ----- | --- | --- | --- | --- | --- | --- |
So a general hyperbolic system will carry information at (real) finite velocities, and these
velocities will be given by the eigenvalues of the Jacobian matrix of the flux. In the case of the
|     |     |     |     | n   | Rd, |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- |
Euler system these velocities, in any normalized direction are given by :
∈
|     |     | u   | n, u n+c, | u n c, |     |     | (1.3) |
| --- | --- | --- | --------- | ------ | --- | --- | ----- |
|     |     |     | · ·       | · −    |     |     |       |
| u   |     |     | c p(ρ)    |        |     |     |       |
where is the velocity of the fluid and := is the sound velocity. For very simple case
(cid:48)
of conservation laws (for example see [2, Introduction p.20] on Burgers’ equation), one can
(cid:112)
exhibit multiple solutions for given initial conditions. Hence, a supplementary constraint needs
to be verified in order to ensure uniqueness of the solution. This constraint is imposed with
| the mathematical | entropy. |     |     |     |     |     |     |
| ---------------- | -------- | --- | --- | --- | --- | --- | --- |
Definition 1.2.2 52]). An entropy for a system of conservation laws is
|     | (Entropy | [3, p |     |     |     |     | (1.2) |
| --- | -------- | ----- | --- | --- | --- | --- | ----- |
a twice differentiable and strictly convex function η : Rk R. Let 1 f d Rk the
j
|     |     |     |     | −→  |     | ≤ ≤ | ⊂   |
| --- | --- | --- | --- | --- | --- | --- | --- |
columns of the flux F(U) Rk d. The associated entropy flux is the function ξ : Ω Rk
×
|     |     | ∈   |     |     |     |     | −→  |
| --- | --- | --- | --- | --- | --- | --- | --- |
that verifies
t
|     |     | [∇η(U)] | f j (U) = ξ j | (U) for 1 | j d, |     |     |
| --- | --- | ------- | ------------- | --------- | ---- | --- | --- |
|     |     |         |               | ≤         | ≤    |     |     |
and a weak solution of (1.2) will verify the following entropy inequality (in the distributional
sense)
|     |     |     | ∂ η(U)+div(ξ(U)) | 0.  |     |     |     |
| --- | --- | --- | ---------------- | --- | --- | --- | --- |
t
≤
The mathematical entropy defined here is an abstract notion that coincides in fact, in the
case of the Euler system, with the well-known physical notion of entropy. Indeed, the second
principle of thermodynamics states that the physical entropy of a system must either stay the

| CHAPTER |     | 1.  | INTRODUCTION |     |     |     |     |     |     |     | 16  |
| ------- | --- | --- | ------------ | --- | --- | --- | --- | --- | --- | --- | --- |
same, or increase. In the case of the Euler system, the second principle of thermodynamics is
| written | for | the | specific | entropy | s as | :               |     |     |     |     |     |
| ------- | --- | --- | -------- | ------- | ---- | --------------- | --- | --- | --- | --- | --- |
|         |     |     |          |         |      | ∂ (ρs)+div(ρsu) |     | 0.  |     |     |     |
t
≥
Leading in this case to a mathematical entropy (Definition 1.2.2) equal to
|       |     |     |               | η   | := ρs | and the   | entropy flux | ξ := | ρsu. |     |     |
| ----- | --- | --- | ------------- | --- | ----- | --------- | ------------ | ---- | ---- | --- | --- |
|       |     |     |               |     | −     |           |              |      | −    |     |     |
| 1.2.2 |     | The | dimensionless |     | Euler | equations |              |      |      |     |     |
WhentalkingaboutthelowMachnumberregime, itisveryconvenienttorewritetheequations
intheirdimensionlessform. ItwillenabletomakeappeartheMachnumber,whichencapsulates
the weakly or strongly compressible nature of the flow in a dimensionless parameter. Given
the domain Ω in which the flow will evolve and initial conditions, it is generally reasonable
to suppose that we know characteristic quantities such as: a time scale t , a space scale x , a
|     |     |     |     |     |     |     |     |     |     | 0   | 0   |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
density scale ρ 0 . With this we can define dimensionless corresponding variables such as
|     |     |     |     | t   |      | x       | ρ     | p   |     | u   |     |
| --- | --- | --- | --- | --- | ---- | ------- | ----- | --- | --- | --- | --- |
|     |     |     |     | t˜= | , x˜ | = , ρ˜= | , p˜= |     | ,u˜ | = , |     |
|     |     |     |     | t   |      | x       | ρ     | p   |     | u   |     |
|     |     |     |     | 0   |      | 0       | 0     | 0   |     | 0   |     |
x
with u = 0 , p = p(ρ ). Then we have c2 = p(ρ ) and the reference Mach number
|     | 0   |     | 0   | 0   |     | 0   | (cid:48) 0 |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | ---------- | --- | --- | --- | --- |
t 0
u
0.
M :=
c
0
Plugging these variables in system (1.1), we obtain the following system
|     |     |     |     |     |     | ∂ ρ˜+div | (ρ˜u˜) = | 0,  |     |     |       |
| --- | --- | --- | --- | --- | --- | -------- | -------- | --- | --- | --- | ----- |
|     |     |     |     |     |     | t˜       | x˜       |     |     |     |       |
|     |     |     |     |    |     |          |          |     |     | ,   | (1.4) |
1

|     |     |     |     |   ∂ | t˜ (ρ˜u˜)+div | x˜ (ρ˜u˜ | u˜)+  | ∇   | x˜ p˜= | 0,  |     |
| --- | --- | --- | --- | ----- | ------------- | -------- | ----- | --- | ------ | --- | --- |
|     |     |     |     |       |               |          | ⊗ γM2 |     |        |     |     |


|     |     | c2  | ρ0. |    |     |     |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
where γ := 0 As we can see it under the pressure gradient, the system obviously exhibits
p0
a singular limit when the Mach number tends to 0. This singular limit has been extensively
studied on the continuous Euler and Navier-Stokes barotropic system [4, 5, 6, 7] and strong
convergencetotheincompressiblesystem(onperiodicdomain)hasbeenidentifiedtooccuronly
when constraints on the initial conditions are verified. Initial conditions that verify these hy-
pothesis are said to be well-prepared initial conditions. On a bounded domain, these conditions
| must | be  | enriched | with | well-prepared |     | boundary | conditions. |     |     |     |     |
| ---- | --- | -------- | ---- | ------------- | --- | -------- | ----------- | --- | --- | --- | --- |
Definition 1.2.3 (Well-prepared initial and boundary conditions). We say that the initial

| CHAPTER    |     | 1.     | INTRODUCTION |            |           |          |                   |      |       |            | 17    |
| ---------- | --- | ------ | ------------ | ---------- | --------- | -------- | ----------------- | ---- | ----- | ---------- | ----- |
| conditions |     | of the | Euler        | barotropic | equations |          | are well-prepared |      | if    |            |       |
|            |     |        | ρ˜(x˜,t˜=    |            |           | ρ˜(0)+   | (M2)              |      | ρ˜(0) | R+,        |       |
|            |     |        |              | 0,M)       | =         |          |                   | with |       |            |       |
|            |     |        |              |            |           | O        |                   |      |       | ∈          |       |
|            |     |        |             |            |           |          |                   |      |       |            | (1.5) |
|            |     |        |  u˜(x˜,t˜=  |            |           | ( 0)     |                   |      |       | ( 0)       |       |
|            |     |        |             | 0,M)       | =         | u˜ (x˜)+ | (M),              | with | div   | x˜ u˜ = 0, |       |
|            |     |        |              |            |           | 0        | O                 |      |       | 0          |       |

| as  | in [4] | and that | the boundary |           | conditions | are     | well-prepared |     | if [8] | :   |     |
| --- | ------ | --------- | ------------ | --------- | ---------- | ------- | ------------- | --- | ------ | --- | --- |
|     |        |           | ρ˜           | (x˜,t˜) = | ρ˜ (0)     | + (M2), |               |     |        |     |     |
|     |        |           | b            |           | b          |         |               |     |        |     |     |
O
|     |     |     |              |     | (0)      |     |      |     | (0) |          |     |
| --- | --- | --- | ------------ | --- | -------- | --- | ---- | --- | --- | -------- | --- |
|     |     |     |  u˜ (x˜,t˜) | =   | u˜ (x˜)+ | (M) | with |     | u˜  | ndΓ = 0. |     |
|     |     |     | b            |     | b        |     |      |     | b   |          |     |
|     |     |     |             |     |          | O   |      | ∂Ω  | ·   |          |     |
(cid:90)
(cid:101)
|       | (0) |      |            |     |        | t˜. |     |     |     |     |     |
| ----- | --- | ---- | ----------- | --- | ------ | --- | --- | --- | --- | --- | --- |
| where | ρ˜  | does | no t depend | on  | x˜ nor |     |     |     |     |     |     |
b
In other words, these particular constraints on the initial and boundary conditions will
encapsulate in mathematical terms the physical intuition that a fluid is near the incompressible
state: the initial velocity is ’close’ to a divergence-free velocity and the density is ’close’ to a
constant.
|     | Sum | up  |     |     |     |     |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
Despite the physical phenomena in a nuclear reactor being very rich and complex, we
motivate here the restriction to the isentropic/barotropic Euler equations; a system that
describesaninviscidflowwithoutgravity. Itisfullydescribedbyadensityandavelocity.
|     | While | far from    | applications, |              | this | model |     |     |     |     |     |
| --- | ----- | ----------- | ------------- | ------------ | ---- | ----- | --- | --- | --- | --- | --- |
|     | 1)    | is a system | of            | conservation |      | laws, |     |     |     |     |     |
2) already contains the problematic of low Mach number approximation,
As a consequence, its approximation forces, on the one side the numerical emulation
of conservation and, on the other side, its careful approximation in the low Mach
number regime. The two challenges we are trying to tackle concomitantly.
Let’s now dive deeper in the problematic and history of the low Mach number flows ap-
proximation.
| 1.3 | Low |     | Mach | number |     | flows | simulation: |     | an  | overview |     |
| --- | --- | --- | ---- | ------ | --- | ----- | ----------- | --- | --- | -------- | --- |
A plethora of numerical methods for the approximation and simulation of fluid dynamics exists
and their relevance depends on the application and underlying model. Also, the difficulties of
the numerical approximation of the Euler equations at low Mach number have been known for
decades and have been extensively studied; depending on the numerical schemes, two types
of symptoms have been highlighted: a computational efficiency problem mainly due to time
integration strategies and an accuracy problem, linked to the spatial discretization of the

CHAPTER 1. INTRODUCTION 18
acoustic operator.
The efficiency problem is related to the disparity between the waves velocity (1.3), the
convective term evolves with the fluid velocity, which is by definition of the low Mach number
regime, smaller than the acoustic speed. In fact, for stability necessities, explicit schemes
will require the time step to verify a condition dependent of the mesh size and the inverse of
the fastest wave: the well-known Courant-Friedrichs-Lewy condition (CFL); thus leading to
inacceptably small time step (this is illustrated in Figure 1.3).
Similarly, implicit time integration leads to stiff systems whose fan of eigenvalues is
particularly large in the complex plane. This is a problem since iterative solvers for implicit
systems are generally preferred to direct solvers because of the latter’s lack of memory scaling
potential. However, as shown for example in [9, Proposition 6.15 p 194], iteratively solving
a system is more difficult when the eigenvalues are scattered. Relying on semi-time implicit
t
1
u
| |
∆t( u )
| |
1
u c
| ± |
∆t( u c )
| ± |
x
∆x
Figure 1.3: Space-time representation of the time stepping constraint when u << c.
integrations will help to free oneself from the acoustic waves CFL restriction [10, 11].
Theaccuracyproblemisrelatedtothenon-convergenceofthediscretecompressiblesolution
totheincompressibleonewhentheMachnumberdropsto0. Itisstronglyrelatedtothespatial

| CHAPTER       | 1. INTRODUCTION |               |                    |     | 19  |
| ------------- | --------------- | ------------- | ------------------ | --- | --- |
| approximation | of the          | kernel of the | acoustic operator: |     |     |
|               |                 |               | 0                  | div |     |
.
|     |     |     | ∇        | 0        |     |
| --- | --- | --- | -------- | -------- | --- |
|     |     |     | (cid:18) | (cid:19) |     |
In this thesis, we will mostly focus on this spatial discretization problem, therefore, we will
mention two families of spatial methods that are of interest for us: collocated and staggered
| finite volume | methods. |     |     |     |     |
| ------------- | -------- | --- | --- | --- | --- |
Ontheonehand,thecollocated finite volume methodarisesnaturallyinthesimulation
of compressible flows and its efficiency for transonic and supersonic flows is known [12, 13, 14].
On the other hand, staggered schemes are more naturally formulated in the context of
| incompressible | flows      | [15, 16]. Let’s | see in detail | these two categories. |     |
| -------------- | ---------- | --------------- | ------------- | --------------------- | --- |
| 1.3.1          | Collocated | schemes         |               |                       |     |
ρ
u
|     |     |     | (cid:18) | (cid:19) |     |
| --- | --- | --- | -------- | -------- | --- |
Figure1.4: RepresentationofthelocalisationofunknownsforCollocatedschemesonCartesian
grids
Collocated finite volumes are based on the formulation of equations under a system of
| conservation | laws (1.2): | in one space | dimension | we have   |     |
| ------------ | ----------- | ------------ | --------- | --------- | --- |
|              |             |              | ∂ U+∂     | F(U) = 0. |     |
t x
(x
Now let i+1/2 ) i 1,...,N the nodes of a mesh defined on the domain of definition of the equa-
tions: the main idea ∈{ is } to integrate on a cell [x ,x ] the previous equation (this will
i 1/2 i+1/2
−
correspond to a cell of the mesh). Setting ∆x := [x ,x ] the Lebesgue measure of the
|     |     |     |     | i i 1/2 i+1/2 |     |
| --- | --- | --- | --- | ------------- | --- |
| − |
cell and
1
|     |     | U   | :=  | Udx, |     |
| --- | --- | --- | --- | ---- | --- |
i
∆x i
|     |     |     | (cid:90)[x | i−1/2 ,x i+1/2] |     |
| --- | --- | --- | ---------- | --------------- | --- |

| CHAPTER |       | 1.            | INTRODUCTION |     |         |       |      | 20  |
| ------- | ----- | ------------- | ------------ | --- | ------- | ----- | ---- | --- |
| this    | gives | the following | formulation  |     |         |       |      |     |
|         |       |               |              |     | F i+1,i | F i,i |      |     |
|         |       |               |              | ∂   | U +     | − 1   | = 0, |     |
|         |       |               |              | t   | i       | −     |      |     |
|         |       |               |              |     |         | ∆x i  |      |     |
where F and F are arbitrary approximations of the physical fluxes at the interfaces
|     |     | i+1,i | i,i 1 |     |     |     |     |     |
| --- | --- | ----- | ----- | --- | --- | --- | --- | --- |
−
x and x . As a consequence, under this formulation, it remains to determine relevant
| i+1/2 |     | i 1/2 |     |     |     |     |     |     |
| ----- | --- | ----- | --- | --- | --- | --- | --- | --- |
−
approximationofthefluxesF ,F inordertoensureconsistencywiththeinitialproblem,
|     |     |     |     | i+1,i | i,i 1 |     |     |     |
| --- | --- | --- | --- | ----- | ----- | --- | --- | --- |
−
stability, entropy dissipation etc. They benefit from a well-developed theory [2, 14, 17, 18] and
| also | from | simple, | convenient | properties: |     |     |     |     |
| ---- | ---- | ------- | ---------- | ----------- | --- | --- | --- | --- |
1) These schemes are simple to generalize in multiple space dimensions.
2) Similarly, they are convenient when it comes to discretizing different conservation laws
|     | because | it  | all boils | down to modifying |     | the flux . |     |     |
| --- | ------- | --- | --------- | ----------------- | --- | ---------- | --- | --- |
3) They are easy to implement in a code, which makes them operationally efficient.
4) Finallyandthisisgoingthemostimportantremarkinourcontext: bynaturallyimposing
that the flux (at an interface) going out of a cell is the opposite of the flux arriving in a
|     | cell | (with            | this common | interface): |     |       |       |     |
| --- | ---- | ---------------- | ----------- | ----------- | --- | ----- | ----- | --- |
|     |      |                  |             |             | i   | F =   | F ,   |     |
|     |      |                  |             |             |     | i+1,i | i,i+1 |     |
|     |      |                  |             |             | ∀   | −     |       |     |
|     | they | are conservative |             | .           |     |       |       |     |
As it is, it almost seems like these schemes are the perfect first order schemes. However it is
widely known that they are not accurate at low Mach number on quadrangular (2 dimensions),
hexahedral (3 dimensions) and Cartesian grids [19, 20, 21, 22, 23, 24]. Indeed, this can be
| highly | problematic |     | for conception | in  | an industrial | context. |     |     |
| ------ | ----------- | --- | -------------- | --- | ------------- | -------- | --- | --- |
It is notoriously difficult to simulate numerically the Euler equations at low Mach number
with collocated schemes. A huge body of literature is available on the subject: in [20, 23], the
behaviouroftheRoeschemeinthelowMachnumberregimeisstudiedwithaonescaleasymp-
totic analysis in Mach number: the authors recall that a solution arising from well-prepared
(Definition 1.2.3) initial conditions should preserve the (M2) behaviour on the density. They
O
show that even if the initial data is well-prepared in the sense of Definition 1.2.3, classical
Riemann solvers will introduce a spurious acoustic mode of order (M) in the numerical dens-
O
ity. Similarly, [25] proposes to study the low Mach number behaviour by restricting to the
study of the wave system kernel; they show the loss of preservation of continuous stationary
state. More generally, the excessive diffusion intrinsic to these numerical schemes leads to the
significant loss of structure (acoustic kernel) preservation. In other words, even for initial data
that are close to the incompressible regime, so by definition, with almost negligible acoustic
information, these schemes inoculate a spurious (acoustic) mode in the approximated solution.
A plethora of solutions has been proposed to reduce the numerical diffusion, some will focus on
preconditioning techniques of the diffusion matrix [26, 27] while others will affect directly the

| CHAPTER | 1. INTRODUCTION |     |     |     |     | 21  |
| ------- | --------------- | --- | --- | --- | --- | --- |
diffusionoperator[20,25,28,29,30,31,32]. Mostofthemboildown, insomeway, tocentering
the pressure gradient. In a similar spirit semi-implicit schemes in time have been developed
[33, 34, 35, 36, 37], in which only the (centered) pressure gradient is implicited, improving at
once the CFL condition and the behaviour with respect to the Mach number.
Yet, pressure centered schemes exhibit at least one big flaw: they create what has been
designated in the literature as ’checkerboard modes’ [38]. Until the recent work of [39], these
modes were badly characterized. Crudely speaking, they can be understood as high frequen-
cies modes generated by the velocity even-odd mode. For classical, typically highly diffusive,
schemes, the high frequencies are smoothed rapidly thus hiding such modes [38, Theorem 3.4,
3)]. By contrast, in the pressure-centered case, the remaining numerical diffusion will enforce
that stationary pressures are indeed constant whereas the lack of velocity diffusion will harm
the discrete stationary velocity space, thus disrupting the properties of the numerical scheme
[40]. Moreover, [39, Chapter 6 , Lemma 6.6, p 118] shows that these spurious ’checkerboard’
modes will live in the velocity kernel and their actual number is huge: the dimension of the
spurious space is of the order of the number of interior nodes of the mesh.
Note that in the case of triangular (2 dimensions) and tetrahedral (3 dimensions) meshes,
the low Mach number accuracy is recovered with the Roe scheme without any modification
of the diffusion [22, 41, 42], furthermore this particularity is kept when increasing the order
of approximation with the Discontinuous Galerkin method [8]. This special case puts in light
a very interesting behaviour of numerical schemes with respect to mesh topology, degrees of
freedom localization AND above all structure preservation in the matter of low Mach number
precision.
Structure preservation: the recipe for low Mach number precision ? The particu-
larities of simplicial meshes, and what distinguishes them from the quadrilateral case, are well
illustrated in [43]. In fact, thanks to a two scale asymptotic analysis in Mach number an equi-
valence between the low Mach number limit and the long-time limit of a wave system coupling
the first order pressure and the zeroth order velocity is established. From there, it is concluded
that a Godunov-type scheme is low-Mach-number accurate on Euler equations if and only if its
matching discretization is consistent in long time on a wave system of the type:
1
|     |     | p        | 0            | div | p                 |       |
| --- | --- | -------- | ------------ | --- | ----------------- | ----- |
|     |     | ∂        | +            | ρ   | = 0.              | (1.6) |
|     |     | τ        |              | 0   |                   |       |
|     |     | u        |             |     |  u               |       |
|     |     | (cid:18) | (cid:19) κ ∇ | 0   | (cid:18) (cid:19) |       |
0
|     |     |     |    |     |    |     |
| --- | --- | --- | --- | --- | --- | --- |
Through this type of analysis we understand that the core of the problem is a double limit
| problem: | do we have | for the discrete | solution     | U ∆x of | (1.6)          |     |
| -------- | ---------- | ---------------- | ------------ | ------- | -------------- | --- |
|          |            | lim              | lim U ∆x (τ) | = lim   | lim U ∆x (τ) ? |     |
|          |            | τ                |              | 0τ      |                |     |
|          |            | + ∆x             | 0            | ∆x      | +              |     |
|          |            | −→ ∞             | →            | →       | −→ ∞           |     |

| CHAPTER | 1. INTRODUCTION |     |     |     |     | 22  |
| ------- | --------------- | --- | --- | --- | --- | --- |
In particular, with this approach, the low Mach number approximation problem boils down to
capturing numerically the accurate space of long time stationary states of this wave system.
In fact, ironically, while the case of a periodic domain is not evident to study as far as
the existence of a long-time limit is concerned, it is actually very convenient when introducing
the matter of structure preservation. Indeed, the long-time limit (1.6) is identified through
decompositions that are similar to the following, given on a 2D torus by:
R2
|     |     | u = ∇ϕ+∇ | ⊥ β +ω | ω ϕ,β | fields, |     |
| --- | --- | -------- | ------ | ----- | ------- | --- |
∈
which is called a Hodge-Helmholtz decomposition (HHD). The case of collocated schemes on
triangular meshes and quadrangular meshes can be separated here: on triangular meshes, this
| type of decomposition |     | exists | at the discrete | scale: |     |     |
| --------------------- | --- | ------ | --------------- | ------ | --- | --- |
Theorem 1.3.1 (HHD for simplicial meshes in collocated framework, particular case of [43,
Proposition 5 ] ). Let Ω = T2 R2 the two dimensional torus. We have the following HHD
⊂
| on simplicial | meshes; |         |          |           |      |     |
| ------------- | ------- | ------- | -------- | --------- | ---- | --- |
|               |         |         |          | L2        | L2   |     |
|               |         | dP0(Ω)2 | ∇ CR1(Ω) | ⊥ ∇ P1(Ω) | ⊥ R2 |     |
= ⊥
|     |     |     |     | ⊕   | ⊕   |     |
| --- | --- | --- | --- | --- | --- | --- |
where dP0(Ω)2 is the space of cellw(cid:2)ise const(cid:3)ant vecto(cid:2)r fields(cid:3)of R2, P1(Ω) is the space of
continuous finite element of degree 1, CR1(Ω) is the space of first-order Crouzeix-Raviart [44].
dQ0(Ω)d
Whereas, this decomposition is lost in quadrangles/hexahedrals for (fields that
are constant by cell on quad/hexa meshes) due to a potential insufficient amount of degrees
of freedom. Confirming the importance of decompositions of the type Theorem 1.3.1 in the
theoretical analysis of low Mach number flows, it is shown in [45] that low Mach number
precisioncanberecoveredinquadrangles, byanastuteenrichmentofthespacevelocityvectors
that leads to the existence of a similar HHD on this new velocity space.
All in all, the low Mach number approximation fits in the larger class of structure preser-
vation problems, andtheparticularstructureofinteresthereisaHHD.Letusnowintroduce
| staggered | schemes   | and how they | are linked | to this. |     |     |
| --------- | --------- | ------------ | ---------- | -------- | --- | --- |
| 1.3.2     | Staggered | schemes      |            |          |     |     |
Contrasting with collocated schemes, the staggered schemes will not place all of the unknowns
at the center of the cell: scalar unknowns are indeed placed in the center of the cell whereas the
normal component of vectorial unknowns, such as the velocity, are placed on the faces of the
meshFigure1.5. Naturally, thisstaggeredschemewasintroducedintheincompressiblesetting,
with the MAC scheme [15, 16], with the aim of bettering the preservation of the divergence-free
velocity constraint. As a matter of fact, the motivation behind the staggering of the velocities
degree of freedom is simply found in the Green-Ostrogradsky formula; integrating the velocity
divergence operator on a volume is equal to the boundary integral of the normal component of
the velocity:

| CHAPTER | 1. INTRODUCTION |     |          |     |        | 23  |
| ------- | --------------- | --- | -------- | --- | ------ | --- |
|         |                 | 1   |          |     | 1      |     |
|         | div(u)          |     | div(u)dx | =   | u ndΓ. |     |
K
|     |     | | ≈ K |            |                   | K ·             |     |
| --- | --- | ----- | ---------- | ----------------- | --------------- | --- |
|     |     | | |   | (cid:90) K | Green-Ostogradsky | | | (cid:90) ∂K |     |
where K is a cell , K its Lebesgue measure, n (cid:124)(cid:123)(cid:122)R(cid:125)d the exterior normal. As a result, it
|     | |   | |   |     | ∈   |     |     |
| --- | --- | --- | --- | --- | --- | --- |
is easily computed with the given degrees of freedom of velocities (Figure 1.5). Moreover,
u
y
ρ
u
x
Figure 1.5: Representation of the localisation of unknowns for the original Staggered schemes
| (MAC) on | Cartesian grids |     |     |     |     |     |
| -------- | --------------- | --- | --- | --- | --- | --- |
in this incompressible setting, the convergence of the MAC scheme for incompressible flows
has been thoroughly studied in [46, 47]. Two commonly known limitations of the MAC
scheme is that it is only defined on Cartesian grids, and solely for incompressible flows. For
the latter, the scheme was in fact extended to compressible flows in [48, 49, 50, 51]. As for
the former, mixes of finite volumes and finite elements approaches have been introduced in
order to address this problem. Indeed, more generally, multiple links between finite volume
schemes and finite element discretizations have been established in various setups; by, for
example, using mass-lumping strategies: [52] for nodal finite elements on Euler equations,
[53, 54] for mixed formulation of the diffusion equations, and mostly, for our interest, vectorial
(to be understood on each component of the velocity field) Crouzeix-Raviart elements for
Navier-Stokes equations [55, 56, 57]. The latter is the result of a long line of research started
on the scalar linear advection diffusion equation, connecting a finite volume discretization
and the scalar Crouzeix-Raviart finite elements through sub-integration and mass-lumping
techniques [58]. It has also been extended to non-linear scalar problems [59, 60] and to
the incompressible Navier-Stokes equations with, this time, vectorial Crouzeix-Raviart [44]
elements in [61, 62]. Similar extension have been made on quadrangular/hexahedral meshes
| with Rannacher-Turek[63] |     | finite | elements | [64]. |     |     |
| ------------------------ | --- | ------ | -------- | ----- | --- | --- |
The authors of [55, 56, 57] use this formalism to develop staggered finite volumes schemes
for which the convection term is stable while with a similar approach [65, 66, 67, 68] work
on staggered discretizations that are entropic or entropic up to a remainder dependent of the

| CHAPTER | 1.  | INTRODUCTION |     |     |     |     |     | 24  |
| ------- | --- | ------------ | --- | --- | --- | --- | --- | --- |
mesh size and time step. Using the same core arguments and combining them with particular
diffusion terms [69] obtained energy-dissipative schemes with Euler explicit time integration
and semi-implicit for the shallow water equations. Similarly, [70] obtained, with a mass flux
| inspired | by kinetic | theory, | entropy | dissipation |     | properties. |     |     |
| -------- | ---------- | ------- | ------- | ----------- | --- | ----------- | --- | --- |
Structure preservation: Staggered schemes as natural structure preserving discret-
| izations | ?   |     |     |     |     |     |     |     |
| -------- | --- | --- | --- | --- | --- | --- | --- | --- |
The intuition leads to believe that staggered schemes do behave better in the low
Mach number regime. In the aim to conceive a low-Mach-number accurate scheme we wish to
adapt the work of [43], (which we have introduced in the previous section treating low-Mach-
number accurate collocated schemes) and reverse-engineer their methodology founded on the
existence of a HHD for the velocity approximation space. In parallel, our interpretation, en-
couraged by [71, Theorem 4.1 and Theorem 4.2], is that normal-component based staggered
| schemes | are perfect | candidates |     | to have | a discrete | HHD. |     |     |
| ------- | ----------- | ---------- | --- | ------- | ---------- | ---- | --- | --- |
Based on these premises, we notice that even if we restrict to simplicial meshes, we are not
yet able to exhibit such decomposition for vectorial Crouzeix-Raviart elements.
| Sum | up  |     |     |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
In short,
1) Collocated schemes are convenient and conceptually simple schemes and above all
they are naturally conservative. Despite these desirable qualities, they are not low
|     | Mach accurate |     | on quadrangular |     | or Cartesian | Grids. |     |     |
| --- | ------------- | --- | --------------- | --- | ------------ | ------ | --- | --- |
2) CollocatedschemesthatarenaturallypreciseatlowMachnumber(soonsimplicial
meshes) are based on approximation spaces that hide a discrete Hodge-Helmholtz
Decomposition. The preservation of this structure at the discrete level seem im-
|     | portant | in order | to prove | and | ensure | low Mach number | precision. |     |
| --- | ------- | -------- | -------- | --- | ------ | --------------- | ---------- | --- |
3) Staggered schemes, by their very nature, seems to be good candidates to be low-
Mach-number accurate because of their natural approximation of the divergence
operator. Alineofresearchhasextendedtocompressibleflows,aclassofstaggered
schemes using the full velocity vector at the faces with multiple properties; conser-
vative form of the convection term and entropy dissipation types theorems.
Nevertheless, we are not aware of the existence of a Hodge-Helmholtz Decomposition for
the Crouzeix-Raviart/Rannacher-Turek velocity staggering. We thus propose a different
| staggered | discretization |     | that | preserves | a HHD | at the discrete | scale. |     |
| --------- | -------------- | --- | ---- | --------- | ----- | --------------- | ------ | --- |

CHAPTER 1. INTRODUCTION 25
1.4 Plan of the manuscript
This study consists in developing a class of staggered conservative low-Mach-
number accurate schemes for fluid mechanics. To this aim, we propose to restrict the
analysis to Euler barotropic equations; it is divided as follow:
1) In chapter 2, an equivalence between the low Mach number limit and the long-time
limit of a wave system is recalled. This key simplification thereby reduces the problem
of approximating low Mach number flows to approaching the ω limit set of a known
−
linear system. In this simplified framework, we present the core arguments that lead to
convergenceinlongtime. Inparticular,theidentificationofthelong-timelimitisachieved
through a fundamental mathematical object: a Hodge-Helmholtz decomposition.
2) In chapter 3 we will select the staggered schemes potentially stable on the multi-
dimensional wave system by first looking at an extensive one dimensional von Neumann
and L2 stability and energy dissipation study. In practice, in one space dimension,
staggered schemes (even in the nodal sense) are similar in the eyes of the approxima-
tion space. In general, staggered schemes are centered, but in this chapter, the stability
and oscillation properties of staggered schemes without and with numerical diffusion are
studied. The study will be carried out with different choices of time integration. These
staggered schemes have been implemented in a one dimensional Python code in order to
numerically illustrate the theoretical results.
3) In chapter 4 is introduced the concept of de Rham complexes and their discrete coun-
terparts. We demonstrate why abstract structures such as complexes are relevant in the
context of numerical analysis by highlighting their natural link with Hodge-Helmholtz
type decompositions. Finally, we introduce the discrete de Rham complex that will lead
to our de Rham based staggered discretization .
4) Then, chapter 5 focuses on the study of the long-time limit on the multi-dimensional
wave system of the de Rham numerical scheme. In this multi-dimensional context, we
introduce in a detailed manner our ’de Rham staggered scheme’ for the spatial discretiz-
ation while we investigate different time steppings: an Euler explicit time integration, an
ImEx pressure centered time stepping and a fully Implicit time integration. Using the de
Rham formalism we show the existence of a long-time limit and finally, we prove that the
selected discretizations converge to a consistent long-time limit. Numerical results will
illustrate the obtained theorems based on an implementation of the scheme in the C++
code SolverLab.
5) Finally, chapter 6 focuses on the extension of the de Rham staggered scheme on Euler
barotropic equations. The goal of this chapter is to tackle the difficulty of defining a
conservative convection term. In parallel, the extension is designed in such way that the
low Mach number behaviour will be consistent with our previous findings. The resulting
scheme is implemented in the C++ code SolverLab for different time steppings and is
tested numerically.

| Chapter    |          | 2   |      |       |          |     |            |     |
| ---------- | -------- | --- | ---- | ----- | -------- | --- | ---------- | --- |
| An         | approach |     |      | based |          | on  | asymptotic |     |
| expansions |          |     | to   |       | simplify |     | the study  | of  |
| the        | low      |     | Mach |       | number   |     | limit      |     |
Contents
2.1 Introduction . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 27
2.2 Link between the low Mach number limit and the long time limit
of a wave system through low Mach number asymptotic expansion 28
2.3 Long time limit of the wave system . . . . . . . . . . . . . . . . . . 32
2.4 Conclusion . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 34
2.1 Introduction
In this chapter we investigate a formal simplification that links the low Mach number limit of
| the | dimensionless | barotropic | Euler | equations, |          |        |     |     |
| --- | ------------- | ---------- | ----- | ---------- | -------- | ------ | --- | --- |
|     |               |            |       |            | ∂ ρ˜+div | (ρ˜u˜) | 0,  |     |
|     |               |            |       |            | t˜ x˜    | =      |     |     |
, (2.1)

1
|     |     |     |   ∂ (ρ˜u˜)+div |     | (ρ˜u˜ | u˜)+ | ∇ p˜= 0, |     |
| --- | --- | --- | ---------------- | --- | ----- | ---- | -------- | --- |
|     |     |     |  t˜             |     | x˜    |      | x˜       |     |
|     |     |     |                  |     | ⊗     | γM2  |          |     |
 
|     | c2  | ρ0  |    |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
where γ := 0 to the lo ng time limit of a wave system. Th chapter is presented in three parts:
p0
1) First, in section 2.2 it is quickly shown through a low Mach number asymptotic analysis
how a wave system emerges and, in parallel, in which sense should be defined the time
|     | scale linked | to  | the low Mach | number | limit. |     |     |     |
| --- | ------------ | --- | ------------ | ------ | ------ | --- | --- | --- |
2) Then,thequestionoflongtimeconvergenceofthiswavesystemisdevelopedinsection2.3.
Theprocess ofidentifying thelimit isseparatedmethodically insuch way thatitbecomes
|     | clear what | should | be preserved |     | at the continuous |     | scale. |     |
| --- | ---------- | ------ | ------------ | --- | ----------------- | --- | ------ | --- |
27

| CHAPTER  |       | 2.        | AN APPROACH |        | BASED           |      | ON ASYMPTOTIC |               | EXPANSIONS |         | TO   |        |
| -------- | ----- | --------- | ----------- | ------ | --------------- | ---- | ------------- | ------------- | ---------- | ------- | ---- | ------ |
| SIMPLIFY |       | THE       | STUDY       | OF     | THE             | LOW  | MACH          | NUMBER        | LIMIT      |         |      | 28     |
|          | 3) In | section   | 2.4, we     | gather | the conclusions |      | of            | this chapter. |            |         |      |        |
| 2.2      | Link  | between   |             | the    | low             | Mach | number        |               | limit      | and the | long | time   |
|          | limit |           | of a        | wave   | system          |      | through       | low           | Mach       | number  |      | asymp- |
|          | totic | expansion |             |        |                 |      |               |               |            |         |      |        |
In the seminal paper [5], a relation between the low Mach number limit of the solution of
Euler equations and a wave system is established. However this link is obtained through
Fourier analysis, so naturally, on unbounded domains. Consequently, the extension of the
Rd
result to bounded domains of is not straightforward. In fact, this rigorous result is initially
motivated by an asymptotic analysis in Mach number, and in [43] this theorem is recovered
formally following the two time scales asymptotic expansion given in [72]. The basic principle
is:
t˜/M
1) first, define an acoustic time scale τ = which accounts for the fast acoustic phenom-
ena,
2) then, plug a two-time-scale asymptotic expansion in Mach number M:
N
|     |     |     |     | ϕ˜(x˜,t˜,M) |     | =   | Mnϕ˜(n)(x˜,t˜,τ)+ |     | (MN+1), |     |     | (2.2) |
| --- | --- | --- | --- | ----------- | --- | --- | ----------------- | --- | ------- | --- | --- | ----- |
O
n=0
(cid:88)
|     | in  | (2.1). Taking |     | the time | derivative |     | of (2.2) | we obtain | :   |     |     |     |
| --- | --- | ------------- | --- | -------- | ---------- | --- | -------- | --------- | --- | --- | --- | --- |
N
1
∂ ϕ˜(x˜,t˜,M) = Mn ∂ ϕ˜(n)(x˜,t˜,τ)+ ∂ ϕ˜(n)(x˜,t˜,τ) + (MN+1).
|     |     | t˜  |     |     |     | t˜       |     |     | τ   |          |     |     |
| --- | --- | --- | --- | --- | --- | -------- | --- | --- | --- | -------- | --- | --- |
|     |     |     |     |     |     |          |     | M   |     |          | O   |     |
|     |     |     |     | n=0 |     | (cid:18) |     |     |     | (cid:19) |     |     |
(cid:88)
Plugging the expansion ϕ ρ,ux,uy,... (2.2) in (2.1) will enable to obtain the following
|     |     |     |     | ∈   | {   |     | }   |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
result:
Proposition 2.2.1. Assume that the initial conditions are well-prepared in the sense of
| Definition |     | 1.2.3. | Then, |     |      |      |            |              |     |     |     |       |
| ---------- | --- | ------ | ----- | --- | ---- | ---- | ---------- | ------------ | --- | --- | --- | ----- |
|            |     |        |       |     | ρ(0) | does | not depend | on x˜,t˜,τ˜, |     |     |     | (2.3) |
and the first order pressure p˜(1) is coupled with the zeroth order momentum (ρu)(0) through the
| following |     | wave system |     | (in acoustic | time | scale): |     |     |     |     |     |     |
| --------- | --- | ----------- | --- | ------------ | ---- | ------- | --- | --- | --- | --- | --- | --- |
2
|     |     |     |     |     | ∂ p˜(1)+γ |         | c˜(0) div | (ρ˜(0)u˜(0)) | = 0, |     |     |       |
| --- | --- | --- | --- | --- | --------- | ------- | --------- | ------------ | ---- | --- | --- | ----- |
|     |     |     |     |     | τ         |         |           | x˜           |      |     |     |       |
|     |     |     |     |    |           | (cid:0) | (cid:1)   |              |      |     |     | (2.4) |
1
|     |     |     |     |   | ∂   | (ρ˜(0)u˜(0))+ |     | ∇ p˜(1) | 0,  |     |     |     |
| --- | --- | --- | --- | --- | --- | ------------- | --- | ------- | --- | --- | --- | --- |
|     |     |     |     |    | τ   |               |     | x˜      | =   |     |     |     |
γ

 

CHAPTER 2. AN APPROACH BASED ON ASYMPTOTIC EXPANSIONS TO
SIMPLIFY THE STUDY OF THE LOW MACH NUMBER LIMIT 29
with initial conditions,
p˜(1)(x˜,τ = 0) = 0,
(2.5)
(cid:40) u˜(0)(x˜,τ = 0) = u˜ ( 0 0) (x˜), with div x˜ (u˜ ( 0 0) ) = 0
and boundary conditions:
p˜ (1) (x˜,τ) = 0,
b
(2.6)
 u˜ (0) (x˜,τ) = u˜ (0) (x˜) with u˜ (0) ndΓ˜ = 0.
 b b
(cid:90)
∂Ω b ·
For Proposition 2.2.1 to be proven we have to notice first that under regularity assumption
on the equation of state, the following stands
Lemma 2.2.1 (Asymptotic expansion of the pressure). Suppose that the equation of state
ρ˜ p˜(ρ˜) is 2 and that the asymptotic expansion of ρ˜ is given by
−→ C
ρ˜:= ρ˜(0)+Mρ˜(1)+M2ρ˜(2)+ (M3).
O
Then, the asymptotic expansion in Mach number of the dimensionless pressure is equal to
1
p˜= p˜(ρ˜(0))+M p˜(ρ˜(0))ρ(1) +M2 ρ(2)p˜(ρ˜(0))+ (ρ(1))2p˜ (ρ˜(0)) + (M3). (2.7)
(cid:48) (cid:48) (cid:48)(cid:48)
2 O
(cid:35)
(cid:20)
(cid:104) (cid:105)
c
and if we denote c˜:= , then
c
0
γ(c˜(0))2 = p˜(ρ˜(0)) (2.8)
(cid:48)
Proof of Lemma 2.2.1. A simple Taylor expansion yields:
p˜ ρ˜(0)+Mρ˜(1)+M2ρ˜(2)+ (M3) =p˜(ρ˜(0))+p˜(ρ˜(0)) Mρ˜(1)+M2ρ˜(2)+ (M3)
(cid:48)
O O
(cid:18) (cid:19) (cid:18) (cid:19)
2
1
+ p˜ (ρ˜(0)) Mρ˜(1)+M2ρ˜(2)+ (M3) + (M3).
(cid:48)(cid:48)
2 O O
(cid:18) (cid:19)
(2.9)
But
2
Mρ˜(1)+M2ρ˜(2)+ (M3) = M2(ρ˜(1))2+ (M3). (2.10)
O O
(cid:18) (cid:19)
Inserting (2.10) in (2.9) concludes for (2.7). Similarly a Taylor expansion yields
γc˜2 = p˜ ρ˜(0)+Mρ˜(1)+M2ρ˜(2)+ (M3) = p˜(ρ˜(0))+ (M)
(cid:48) (cid:48)
O O
(cid:16) (cid:17)

| CHAPTER      | 2.  | AN APPROACH    |     | BASED  |          | ON ASYMPTOTIC |        | EXPANSIONS |     | TO  |     |
| ------------ | --- | -------------- | --- | ------ | -------- | ------------- | ------ | ---------- | --- | --- | --- |
| SIMPLIFY     | THE | STUDY          | OF  | THE    | LOW MACH |               | NUMBER | LIMIT      |     |     | 30  |
| which yields | by  | identification |     | (2.8). |          |               |        |            |     |     |     |
Proof of Proposition 2.2.1 . Plugging the asymptotic expansion (2.2) on each variable in (2.1)
yields:
1
| At order |     |     | the momentum |     | equation |     | gives; |     |     |     |     |
| -------- | --- | --- | ------------ | --- | -------- | --- | ------ | --- | --- | --- | --- |
O M2
(cid:18) (cid:19)
|     |     |     |     |     | ∇   | p˜(0) | = 0. |     |     |     | (2.11) |
| --- | --- | --- | --- | --- | --- | ----- | ---- | --- | --- | --- | ------ |
x˜
1
At order
O M
(cid:18) (cid:19)
•
| The | density | equation | gives; |     |     |         |      |     |     |     |        |
| --- | ------- | -------- | ------ | --- | --- | ------- | ---- | --- | --- | --- | ------ |
|     |         |          |        |     |     | ∂ ρ˜(0) | = 0. |     |     |     | (2.12) |
τ
| • The | momentum |     | equation | gives; |     |     |     |     |     |     |     |
| ----- | -------- | --- | -------- | ------ | --- | --- | --- | --- | --- | --- | --- |
1
|     |     |     |     |     | ∂ (ρ˜(0)u˜(0))+ |     | ∇   | p˜(1) = 0. |     |     | (2.13) |
| --- | --- | --- | --- | --- | --------------- | --- | --- | ---------- | --- | --- | ------ |
|     |     |     |     |     | τ               |     | x˜  |            |     |     |        |
γ
| At order | (1) |     |     |     |     |     |     |     |     |     |     |
| -------- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
O
•
| The   | density  | equation | gives;   |       |         |           |              |      |     |     |        |
| ----- | -------- | -------- | -------- | ----- | ------- | --------- | ------------ | ---- | --- | --- | ------ |
|       |          |          |          |       | ρ˜(1)+∂ | ρ˜(0)+div | (ρ˜(0)u˜(0)) |      |     |     |        |
|       |          |          |          | ∂ τ   | t˜      |           | x˜           | = 0. |     |     | (2.14) |
| • The | momentum |          | equation | gives | ;       |           |              |      |     |     |        |
1
|                   |     | ∂ (ρ˜u˜)(1)+∂ |                                    | (ρ˜u˜)(0)+div |     |     | (ρ˜(0)u˜(0) | u(0))+                      | ∇ p˜(2) | = 0. |     |
| ----------------- | --- | ------------- | ---------------------------------- | ------------- | --- | --- | ----------- | --------------------------- | ------- | ---- | --- |
|                   |     | τ             |                                    | t˜            |     | x˜  |             |                             | x˜      |      |     |
|                   |     |               |                                    |               |     |     |             | ⊗                           | γ       |      |     |
| ByLemma2.2.1p˜(0) |     |               | p˜(ρ˜(0))sothat,theequationatorder |               |     |     |             | 1                           |         |      |     |
|                   |     | =             |                                    |               |     |     |             | (2.11)andthedensityequation |         |      |     |
M2
| at order | 1 (2.12) | give | together |     |     |     |     |     |     |     |     |
| -------- | -------- | ---- | -------- | --- | --- | --- | --- | --- | --- | --- | --- |
M
|            |             |         | p˜(0)(x˜,t˜,τ) |     | p˜(0)(t˜), |     | ρ˜(0)(x˜,t˜,τ) | ρ˜(0)(t˜). |     |     |     |
| ---------- | ----------- | ------- | -------------- | --- | ---------- | --- | -------------- | ---------- | --- | --- | --- |
|            |             |         |                |     | =          |     |                | =          |     |     |     |
| It implies | that (2.14) | becomes |                |     |            |     |                |            |     |     |     |
d
|     |     |     |     | ∂ ρ˜(1)+div |     | (ρ˜(0)u˜(0)) | =   | ρ˜(0). |     |     | (2.15) |
| --- | --- | --- | --- | ----------- | --- | ------------ | --- | ------ | --- | --- | ------ |
|     |     |     |     | τ           | x˜  |              |     | −dt˜   |     |     |        |
So the zeroth order density is equal to the zeroth order boundary density, which is, by Defin-
ition 1.2.3, time-independent, hence giving (2.3). As a consequence, the right hand side of
|     |     |     |     |     |     | p˜(ρ˜(0)) |     | (c˜(0))2 |     |     |     |
| --- | --- | --- | --- | --- | --- | --------- | --- | -------- | --- | --- | --- |
(2.15) vanishes. Now, multiplying (2.15) by (cid:48) (= by definition) and using again

| CHAPTER  | 2.    | AN APPROACH |       |                   | BASED | ON   | ASYMPTOTIC |     |       | EXPANSIONS |     | TO  |     |
| -------- | ----- | ----------- | ----- | ----------------- | ----- | ---- | ---------- | --- | ----- | ---------- | --- | --- | --- |
| SIMPLIFY | THE   | STUDY       |       | OF THE            | LOW   | MACH | NUMBER     |     | LIMIT |            |     |     | 31  |
| Lemma    | 2.2.1 | to identify | p˜(1) | = ρ˜(1)p˜(ρ˜(0)), |       | we   | have       |     |       |            |     |     |     |
(cid:48)
|     |     |     |     | ∂ p˜(1)+p˜(ρ˜(0))div |     |          | (ρ˜(0)u˜(0)) |     | =   | 0.  |     |     | (2.16) |
| --- | --- | --- | --- | -------------------- | --- | -------- | ------------ | --- | --- | --- | --- | --- | ------ |
|     |     |     |     | τ                    |     | (cid:48) | x˜           |     |     |     |     |     |        |
Gathering (2.13)and (2.16) gives the result for the wave system (2.4). Now the asymptotic
| expansion | on  | the initial | density | is      |      |     |       |     |       |     |     |     |     |
| --------- | --- | ----------- | ------- | ------- | ---- | --- | ----- | --- | ----- | --- | --- | --- | --- |
|           |     |             |         |         | (0)  | (1) |       | (2) |       |     |     |     |     |
|           |     |             |         | ρ˜ = ρ˜ | +Mρ˜ |     | +M2ρ˜ | +   | (M3). |     |     |     |     |
|           |     |             |         | 0       | 0    | 0   |       | 0   |       |     |     |     |     |
O
Since we supposed that the initial conditions are well-prepared Definition 1.2.3, we have by
identification:
(1)
|     |     |     |     |     |     | ρ˜  | = 0. |     |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | ---- | --- | --- | --- | --- | --- | --- |
0
|     |     |     |     | p˜ (1) p˜(ρ˜ | (0) )ρ(1) |     |     |     |     |     |     |     |     |
| --- | --- | --- | --- | ------------ | --------- | --- | --- | --- | --- | --- | --- | --- | --- |
By Lemma 2.2.1 we have = (cid:48) = 0 . Similarly by identification on the asymptotic
|     |     |     |     | 0   | 0   | 0   |     |     |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
expansion of the velocity and the boundary conditions we obtain (2.5) and (2.6).
t˜
|     |     |     |     | When | M   | 0,  | then τ | :=  |     | 0,  |     |     |     |
| --- | --- | --- | --- | ---- | --- | --- | ------ | --- | --- | --- | --- | --- | --- |
M
|     |     |     |     |     | −→  |     |     |     | −→  |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
so this issue of low Mach number limit can be viewed as the long time limit of a wave system.
In particular, (2.4),(2.5),(2.6) pertains to a wider range of first order wave systems:
1
|             |              |     |          | p                 |            | 0      |     | div  | p        |          |     |     |        |
| ----------- | ------------ | --- | -------- | ----------------- | ---------- | ------ | --- | ---- | -------- | -------- | --- | --- | ------ |
|             |              |     | ∂        |                   | +          |        | ρ   |      |          | = 0,     |     |     | (2.17) |
|             |              |     |          | τ                 |            |        | 0   |      |          |          |     |     |        |
|             |              |     |          | u                 |           |        |     |     | u        |          |     |     |        |
|             |              |     |          | (cid:18) (cid:19) |            | κ 0 ∇p | 0   |      | (cid:18) | (cid:19) |     |     |        |
|             |              |     |          |                   |           |        |     |     |          |          |     |     |        |
| with weakly | inlet/outlet |     | boundary |                   | conditions |        |     |      |          |          |     |     |        |
|             |              |     |          |                   |            |        | 1 u | n+u  | n        | c        |     |     |        |
|             |              | 1   |          |                   |            |        |     |      | b        | 0        |     |     |        |
|             |              |     | u n      |                   |            |        | ·   |      | · +      | (p       | p ) |     |        |
|             |              |     |          |                   |            | ρ      |     | 2    |          | 2 −      | b   |     |        |
|             |              | ρ   | ·        |                   | =          |        | 0   |      |          |          |     | ,   | (2.18) |
|             |              | 0   |          |                   |            |       | p+p |      | c        |          |     |    |        |
|             |              |  κ | pn       |                  |            | κ      |     | b n+ | 0 (u     | n u      | n)n |     |        |
|             |              |     | 0        |                   |            | 0      |     |      |          | b        |     |     |        |
|             |              |     |          | inlet/outlet      |            |        | 2   |      | 2 ·      | −        | ·   |     |        |
|             |              |     |          |                   |            |       |     |      |          |          |     |    |        |
|             |              |    |          |                  |            |       |     |      |          |          |     |    |        |
with
|     |     |     |     | with | p   | = cst, |     | u ndΓ | = 0. |     |     |     |     |
| --- | --- | --- | --- | ---- | --- | ------ | --- | ----- | ---- | --- | --- | --- | --- |
|     |     |     |     |      | b   |        |     | b     |      |     |     |     |     |
·
(cid:90) ∂Ω
| and with | weakly | wall | boundary | conditions |     |     |     |     |     |     |     |     |     |
| -------- | ------ | ---- | -------- | ---------- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
1
|     |     |     |     | u   | n   |     |     |     |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
0
|     |     |     |     | ρ      | ·   | =   |          |     |      | .        |     |     | (2.19) |
| --- | --- | --- | --- | ------ | --- | --- | -------- | --- | ---- | -------- | --- | --- | ------ |
|     |     |     |     | 0      |     |     | κ pn+c   |     | u nn |          |     |     |        |
|     |     |     |     |  κ pn |    |     | 0        |     | 0    |          |     |     |        |
|     |     |     |     | 0      |     |     | (cid:20) |     | ·    | (cid:21) |     |     |        |
wall
|     |     |     |     |    |    |     |     |     |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
Remark 2.2.1. Such analysis can be also made on the full Euler system (with total energy
conservation equation) [39, 73] and results in a very similar structure: a wave system dependent

| CHAPTER  |     | 2.  | AN APPROACH |     | BASED |     | ON   | ASYMPTOTIC |     |       | EXPANSIONS | TO  |     |
| -------- | --- | --- | ----------- | --- | ----- | --- | ---- | ---------- | --- | ----- | ---------- | --- | --- |
| SIMPLIFY |     | THE | STUDY       | OF  | THE   | LOW | MACH | NUMBER     |     | LIMIT |            |     | 32  |
on the fast acoustic scale will also arise, though, in this case the coefficients will be non-constant
| because |      | ρ˜(0) will | depend | on x. |     |     |      |        |     |     |     |     |     |
| ------- | ---- | ---------- | ------ | ----- | --- | --- | ---- | ------ | --- | --- | --- | --- | --- |
| 2.3     | Long |            | time   | limit | of  | the | wave | system |     |     |     |     |     |
While it is a simple linear set up, the wave system has a very rich structure. In this regard, it
is known that on periodic domains (in 2D for the sake of conciseness) the system yields for the
| velocity |     | the following |     | preservation |         | of a peculiar |      | operator: |        |     |     |     |     |
| -------- | --- | ------------- | --- | ------------ | ------- | ------------- | ---- | --------- | ------ | --- | --- | --- | --- |
|          |     |               |     |              |         |               |      |           |        | uy  | ux. |     |     |
|          |     |               |     | ∂ τ          | (curlu) | = 0           | with | curlu     | := ∂ x | ∂   | y   |     |     |
−
Inthecaseofaboundeddomain,thestructurepreservedisalittlebitmoresubtle. Tointroduce
| it, | we need | the | following | theorem, |     | that we | do  | not prove: |     |     |     |     |     |
| --- | ------- | --- | --------- | -------- | --- | ------- | --- | ---------- | --- | --- | --- | --- | --- |
Theorem 2.3.1 (Hodge-Helmholtz decomposition adapted to the wave system, [43, 74] ).
Let d 1,2,3 , Rd a Lipschitz open bounded set, u L2(Ω)d, and u is such that
|     |     |          | Ω   |              |     |               |     |               |     |           |     | b   |     |
| --- | --- | -------- | --- | ------------ | --- | ------------- | --- | ------------- | --- | --------- | --- | --- | --- |
|     | ∈   | {        | }   | ⊂            |     |               |     |               |     | ∈         |     |     |     |
| u   | n   | H1/2(∂Ω) |     | and respects |     | the following |     | compatibility |     | condition |     |     |     |
b ∂Ω
|     | · | | ∈   |     |     |     |     |     |       |     |     |     |     |        |
| --- | --- | --- | --- | --- | --- | --- | --- | ----- | --- | --- | --- | --- | ------ |
|     |     |     |     |     |     |     | u   | ndΓ = | 0,  |     |     |     | (2.20) |
b
|     |     |     |     |     |     | ∂Ω  | ·   |     |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
(cid:90)
| then, | there | exists | a unique | decomposition |     |     |        |     |     |     |     |     |     |
| ----- | ----- | ------ | -------- | ------------- | --- | --- | ------ | --- | --- | --- | --- | --- | --- |
|       |       |        |          |               |     | u   | = ∇ϕ+u |     | ,   |     |     |     |     |
Ψ
| with | ϕ   | H1(Ω)/R | and |     |     |     |     |     |     |     |     |     |     |
| ---- | --- | ------- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
∈
|     |     |     |     | u   | v L2(Ω)d,div(v) |     |     | = 0; | v n | = u | n ,   |     |     |
| --- | --- | --- | --- | --- | --------------- | --- | --- | ---- | --- | --- | ----- | --- | --- |
|     |     |     |     | Ψ   |                 |     |     |      |     | ∂Ω  | b ∂Ω  |     |     |
|     |     |     |     | ∈ { | ∈               |     |     |      | · | |     | · | } |     |     |
This decomposition is useful to show that in the case of the wave system on a bounded
| domain, |     | the following |     | invariance | property |     | stands: |     |     |     |     |     |     |
| ------- | --- | ------------- | --- | ---------- | -------- | --- | ------- | --- | --- | --- | --- | --- | --- |
Lemma 2.3.1 (Invariance of divergence-free part, [43, Lemma 2, Appendix B.2] ). Suppose
p
| (2.20) | stands | and | let | U = |     | the solution |     | of the | wave | system. | Let |     |     |
| ------ | ------ | --- | --- | --- | --- | ------------ | --- | ------ | ---- | ------- | --- | --- | --- |
u
|     |     |     |     | (cid:18) | (cid:19) |            |     |        |     |     |         |     |     |
| --- | --- | --- | --- | -------- | -------- | ---------- | --- | ------ | --- | --- | ------- | --- | --- |
|     |     |     | u   | = u ϕ    | +u       | with div(u |     | ) = 0, | u   | n = | u b n , |     |     |
|     |     |     |     |          | Ψ        |            | Ψ   |        | Ψ · | ∂Ω  | · ∂Ω    |     |     |
|     |     |     |     |          |          |            |     |        |     | |   | |       |     |     |
the decomposition of the velocity field (solution of the wave system) as in Theorem 2.3.1. Then,
| the | divergence-free |     | part | is invariant |     | in time |     |      |     |     |     |     |        |
| --- | --------------- | --- | ---- | ------------ | --- | ------- | --- | ---- | --- | --- | --- | --- | ------ |
|     |                 |     |      |              |     |         | ∂ u | = 0. |     |     |     |     | (2.21) |
τ Ψ

| CHAPTER  |     | 2. AN     | APPROACH |        | BASED | ON   | ASYMPTOTIC |        |       | EXPANSIONS | TO  |     |
| -------- | --- | --------- | -------- | ------ | ----- | ---- | ---------- | ------ | ----- | ---------- | --- | --- |
| SIMPLIFY |     | THE STUDY |          | OF THE | LOW   | MACH |            | NUMBER | LIMIT |            |     | 33  |
Thepreservationofthesestructuresareinfactdeeplyconnectedtothelongtimebehaviour:
a reasonable long time limit (p ,u ) is, if it exists, a steady physical state, which is of the
|     |     |     |     |     | ∞ ∞ |     |     |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
form:
|     |     |     |     |     | p   | = cst, | divu | =   | 0.  |     |     |     |
| --- | --- | --- | --- | --- | --- | ------ | ---- | --- | --- | --- | --- | --- |
|     |     |     |     |     | ∞   |        |      | ∞   |     |     |     |     |
Now, steadiness is obtained if the boundary conditions do not inject energy, so that a sensible
candidate should also verify the boundary conditions on the domain’s boundary:
|     |     |     |     | p   | = p  | u   | n   | = u | n . |     |     |     |
| --- | --- | --- | --- | --- | ---- | --- | --- | --- | --- | --- | --- | --- |
|     |     |     |     |     | ∞ ∂Ω | b ∞ | ∂Ω  | b   | ∂Ω  |     |     |     |
|     |     |     |     |     | |    |     | · | |     | · | |     |     |     |
Naturally, now, a candidate for the limit pressure is readily identified as p = p , but it also
|     |     |     |     |     |     |     |     |     |     | ∞   | b   |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
possible to exhibit the existence of the relevant velocity state thanks to Theorem 2.3.1. In fact,
we can show:
Lemma 2.3.2 ( Relative energy dissipation, [43, Lemma 1, Appendix B.2 ] ). Let U = (p,u)t
the solution of the wave system (2.17) with p ,u as initial conditions and with weakly imposed
0 0
boundary conditions (2.18) and wall conditions (2.19) such that (2.20)stands. Define (u 0 ) Ψ the
divergence free part of the initial condition from the decomposition Theorem 2.3.1. Then, the
relative energy of the solution with respect to this state is a Lyapunov functional
d
|     |      |             |      |                   | ρ κ (p | p )2+ | u   | (u  | ) 2dx    | 0.  |     |     |
| --- | ---- | ----------- | ---- | ----------------- | ------ | ----- | --- | --- | -------- | --- | --- | --- |
|     |      |             |      |                   | 0 0    | b     |     | 0   | Ψ        |     |     |     |
|     |      |             | dτ   |                   |        | −     | |   | −   | |        | ≤   |     |     |
|     |      |             |      | (cid:20)(cid:90)Ω |        |       |     |     | (cid:21) |     |     |     |
|     | This | shows that, | with |                   |        |       |     |     |          |     |     |     |
p
|     |     |     |     |     | U   | =   |     | b , |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
∞
|     |     |     |     |     |     |     | (u 0     | ) Ψ      |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | -------- | -------- | --- | --- | --- | --- |
|     |     |     |     |     |     |     | (cid:18) | (cid:19) |     |     |     |     |
Lemma 2.3.2 leads to relative energy dissipation of the solution with respect to this state: it is
then clear that at the limit p p is uniform with 0 as boundary trace (thus trivially, in the
|     |     |     |     | ∞   | b   |     |     |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
−
limit p = p ) and that u (u ) is divergence free. But by Lemma 2.3.1 we know that u
|     | ∞   | b   |     |     | 0 Ψ |     |     |     |     |     |     | Ψ   |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
−
is constant in time, so that in the limit u (u ) is equal to u . However it is noted that
|     |     |     |     |     |     | ∞   | 0   | Ψ   |     | ∞ϕ  |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
−
u (u ) = u is divergence-free and in parallel u n = u n , yielding u n = 0.
|     | 0   | Ψ ∞ϕ |     |     |     |     |     | ∞   | ∂Ω  | b ∂Ω | ∞ϕ ∂Ω |     |
| --- | --- | ---- | --- | --- | --- | --- | --- | --- | --- | ---- | ----- | --- |
| −   |     |      |     |     |     |     |     | · | |     | · |  | · |   |     |
By uniqueness of the decomposition it thus equal to 0. This yields the following fundamental
theorem:
Theorem 2.3.2 (Long time limit of the wave system, [43] ). Let d 2,3 , Ω Rd an open
|     |     |     |     |     |     |     |     |     |     | ∈ { } | ⊂   |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | ----- | --- | --- |
p
| bounded |     | set and U | =   | the | solution | of  | the wave | system | with: |     |     |     |
| ------- | --- | --------- | --- | --- | -------- | --- | -------- | ------ | ----- | --- | --- | --- |
u
(cid:18) (cid:19)
|     | • initial | conditions | p(τ | = 0,x) | :=  | p (x) and | u(τ | = 0,x) | :=  | u (x) |     |     |
| --- | --------- | ---------- | --- | ------ | --- | --------- | --- | ------ | --- | ----- | --- | --- |
|     |           |            |     |        |     | 0         |     |        |     | 0     |     |     |
•
wall boundary conditions (2.19) and weakly imposed boundary conditions (2.18), where

| CHAPTER  |     | 2.  | AN APPROACH |     | BASED   | ON   | ASYMPTOTIC | EXPANSIONS | TO  |     |
| -------- | --- | --- | ----------- | --- | ------- | ---- | ---------- | ---------- | --- | --- |
| SIMPLIFY |     | THE | STUDY       | OF  | THE LOW | MACH | NUMBER     | LIMIT      |     | 34  |
p
|     | U   |     | b are | time-independent |     | and | such that |     |     |     |
| --- | --- | --- | ----- | ---------------- | --- | --- | --------- | --- | --- | --- |
|     | b   | =   |       |                  |     |     |           |     |     |     |
u
|     |     | (cid:18) | b (cid:19) |            |       |          |                   |           |        |     |
| --- | --- | -------- | ---------- | ---------- | ----- | -------- | ----------------- | --------- | ------ | --- |
|     |     |          | p          | is uniform | and u | verifies | the compatibility | condition | (2.20) |     |
|     |     |          | b          |            |       | b        |                   |           |        |     |
p
| Then | the | long | time limit | U   | = ∞ | is  | such that |     |     |     |
| ---- | --- | ---- | ---------- | --- | --- | --- | --------- | --- | --- | --- |
∞
u ∞
|     |     |            |       |     | (cid:18) | (cid:19) |     |     |     |     |
| --- | --- | ---------- | ----- | --- | -------- | -------- | --- | --- | --- | --- |
|     | • p | is uniform | equal | to  | p        |          |     |     |     |     |
|     | ∞   |            |       |     | b        |          |     |     |     |     |
•
u ∞ is equal to (u 0 ) Ψ , the divergence-free part of the initial velocity u 0 extracted with
|     | (2.3.1), |     | which matches |     | u on the | boundary. |     |     |     |     |
| --- | -------- | --- | ------------- | --- | -------- | --------- | --- | --- | --- | --- |
b
2.4 Conclusion
A few insights can be extracted from this chapter, of which, most importantly, the driving
argument of our methodology: we can study the behaviour at low Mach number of the solution
of Euler barotropic equations by examining a wave system. It is key to understand the long
time limit of this wave system at the continuous scale. In particular we have shown, firstly that
this limit, if it exists, is identified thanks to a peculiar Hodge-Helmhlotz decomposition (HHD)
with boundary conditions, secondly that this decomposition yields an invariance property on
the HHD, and finally, that the dissipation of a relative energy with respect to a state identified
by this HHD stands. Hence, we find that the low Mach number flows numerical approximation
boils down to conceiving a discretization that respects the structure of the wave system: the
chosen discretization
•
should enable the identification of the limit with a discrete Hodge-Helmholtz decomposi-
tion,
• should also preserve the fundamental structures of the system, especially the invariance
|     | property |     | linked | to the | HHD, |     |     |     |     |     |
| --- | -------- | --- | ------ | ------ | ---- | --- | --- | --- | --- | --- |
•
|     | and | finally, | should | dissipate | the | relative | energy. |     |     |     |
| --- | --- | -------- | ------ | --------- | --- | -------- | ------- | --- | --- | --- |
Our starting point in this methodology consists in be to investigating energy dissipation
properties. This is the goal of the next chapter, in which energy dissipation of one space
| dimension |     | staggered | discretizations |     | of  | the wave | system is | studied. |     |     |
| --------- | --- | --------- | --------------- | --- | --- | -------- | --------- | -------- | --- | --- |

Chapter 3
Study of the stability of the
staggered schemes for the one
dimensional wave system
Contents
3.1 Introduction . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 36
3.2 The staggered schemes on the one dimensional wave system . . . 39
3.2.1 Staggered schemes . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 39
3.2.2 Finite volume interpretation of the staggered scheme . . . . . . . . . . 40
3.3 von-Neumann analysis and energy dissipation . . . . . . . . . . . . 42
3.3.1 von Neumann stability analysis . . . . . . . . . . . . . . . . . . . . . . 45
3.3.2 Energy dissipation . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 56
3.4 l -stability on the characteristic variables . . . . . . . . . . . . . . 63
∞
3.4.1 Characteristic variables defined on half cells . . . . . . . . . . . . . . . 63
3.4.2 Characteristic variables defined on primal cells . . . . . . . . . . . . . 65
3.5 Discussion on some preexisting staggered schemes through low
Mach number asymptotics . . . . . . . . . . . . . . . . . . . . . . . . 66
3.6 Numerical results . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 67
3.6.1 Numerical study of the amplification matrices . . . . . . . . . . . . . . 68
3.6.2 Tests on the numerical schemes . . . . . . . . . . . . . . . . . . . . . . 69
3.7 Conclusion . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 70
3.1 Introduction
We’ve recalled in chapter 2 that performing a two time scales asymptotic expansion in Mach
number on the Euler system gives a first order wave system (3.2) coupling the first-order
pressure and the zeroth-order momentum (see for example [72, 23, 75, 43]). Then, to propagate
36

CHAPTER 3. STUDY OF THE STABILITY OF THE STAGGERED SCHEMES FOR
| THE | ONE | DIMENSIONAL |     | WAVE | SYSTEM |     |     |     |     | 37  |
| --- | --- | ----------- | --- | ---- | ------ | --- | --- | --- | --- | --- |
acoustic waves in a low Mach number flow, a numerical scheme for the Euler system should
be asymptotically consistent with a discretization of the wave system [75]. Energy dissipation
is a key ingredient when looking at the long time convergence of the continuous wave system
[43] (equipped with boundary conditions). As a consequence, the first step of the process is
ensuring such property, and in particular, stability; in one space dimension it can easily be
studied through a von Neumann analysis of the numerical scheme. Thus, we consider the one
| dimensional, |     | first order | linear | wave | system |     |     |     |     |     |
| ------------ | --- | ----------- | ------ | ---- | ------ | --- | --- | --- | --- | --- |
1
|     |     |     |     |     |     | ∂ p+ | ∂ u = 0, |     |     | (3.1a) |
| --- | --- | --- | --- | --- | --- | ---- | -------- | --- | --- | ------ |
|     |     |     |     |     |     | t    | x        |     |     |        |
ρ
0

|     |     |     |     |     | ∂  | u+κ | ∂ p = 0, |     |     | (3.1b) |
| --- | --- | --- | --- | --- | --- | --- | -------- | --- | --- | ------ |
|     |     |     |     |     |     | t   | 0 x      |     |     |        |
wherepisthepressure,uisthevelocityandκ
|     |     |     |     |     |     |     | andρ aretwostrictlynon-negativeparameters. |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | ------------------------------------------ | --- | --- | --- |
|     |     |     |     |     |     | 0   | 0                                          |     |     |     |
(p,u)T,
The wave velocity c is linked with the two parameters by c2 = κ /ρ . Noting U = the
|        |     |              | 0       |         |     |     |     | 0 0 | 0   |     |
| ------ | --- | ------------ | ------- | ------- | --- | --- | --- | --- | --- | --- |
| system |     | (3.1) can be | briefly | written | as  |     |     |     |     |     |
1
0
|     |     |     | ∂ U+B∂ |     | U = 0 | where | B = | ρ   | .   | (3.2) |
| --- | --- | --- | ------ | --- | ----- | ----- | --- | --- | --- | ----- |
|     |     |     | t      |     | x     |       |     | 0   |     |       |
|     |     |     |        |     |       |       |    |     |    |       |
κ 0
0
|     |     |     |     |     |     |     |    |     |    |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
The continuous wave system (3.2) verifies two fundamental properties; a supplementary con-
servation equation on the energy, and, due to its hyperbolic linear nature, transport equations
| on  | characteristic | variables: |     |     |     |     |     |     |     |     |
| --- | -------------- | ---------- | --- | --- | --- | --- | --- | --- | --- | --- |
• Concerning the energy; multiplying the system (3.2) by (SU)t with S a symmetrizer of
|     | the | wave system |     |     |     |      |     |     |     |     |
| --- | --- | ----------- | --- | --- | --- | ---- | --- | --- | --- | --- |
|     |     |             |     |     |     |      | 1 0 |     |     |     |
|     |     |             |     |     |     | S := | 1   | ,   |     |     |
|     |     |             |     |     |     |     | 0  |     |     |     |
ρ2c2
0 0
|     |       |               |     |            |              |      |         |      |     |     |
| --- | ----- | ------------- | --- | ---------- | ------------ | ----- | -------- | ---- | --- | --- |
|     | gives | the following |     | additional | conservation |       | law      |      |     |     |
|     |       |               |     |            | p2           | u2    | 1        |      |     |     |
|     |       |               |     |            | ∂            | +     | + ∂ (pu) | = 0. |     |     |
|     |       |               |     |            | t 2          | 2ρ2c2 | ρ x      |      |     |     |
0
|     |     |     |     |     | (cid:18) | 0   | 0(cid:19) |     |     |     |
| --- | --- | --- | --- | --- | -------- | --- | --------- | --- | --- | --- |
Assuming that the domain is periodic and corresponds to [0;1[, the energy E defined by
1 1
UTSUdx,
|     |     |     |     |     | E(U) | =   |     |     |     | (3.3) |
| --- | --- | --- | --- | --- | ---- | --- | --- | --- | --- | ----- |
2
(cid:90)0
|     | verifies | as a consequence |     |     |     |        |      |     |     |       |
| --- | -------- | ---------------- | --- | --- | --- | ------ | ---- | --- | --- | ----- |
|     |          |                  |     |     |     | ∂ E(U) | = 0. |     |     | (3.4) |
t
• Concerningthecharacteristic variables weobtainthroughdiagonalizationofthesystema
L -stabilityproperty: indeed, itisknownthatthematrixB from(3.2)hascharacteristic
∞

CHAPTER 3. STUDY OF THE STABILITY OF THE STAGGERED SCHEMES FOR
THE ONE DIMENSIONAL WAVE SYSTEM 38
polynomial
κ
p(X) := det(B XI ) = X2 0 = X2 c2,
− d − ρ − 0
0
so it is diagonalizable with
1
1
B = P − c 0 0 P 1 with P = 1 1 and P 1 = 1 −ρ 0 c 0 .
0 c − ρ c ρ c − 2  1 
(cid:18) 0 (cid:19) (cid:18) − 0 0 0 0 (cid:19) 1
ρ c
 0 0 
 
Now (3.2) is equivalent to
c 0 c 0
∂ t U+P − 0 P − 1∂ x U = 0 ∂ t U+P − 0 ∂ x (P − 1U) = 0.
0 c ⇐⇒ 0 c
0 0
(cid:18) (cid:19) (cid:18) (cid:19)
Applying P 1 to last equality we get the following two linear advection equations on the
−
characteristic variables p u/(ρ c ) and p+u/(ρ c ):
0 0 0 0
−
u u
∂ p +c ∂ p = 0,
t 0 x
− ρ c − ρ c
 (cid:18) 0 0(cid:19) (cid:18) 0 0(cid:19)
 u u

 ∂
t
p+ c
0
∂
x
p+ = 0.
ρ c − ρ c
(cid:18) 0 0(cid:19) (cid:18) 0 0(cid:19)



In the whole chapter, weassume that the domain is periodic and equals to [0;1[. It leads
then to the following L stability on the characteristic variables
∞
−
u(t, ) u(t = 0, )
t > 0, x [0;1[, p(t, ) · = p(t = 0, ) · .
∀ ∀ ∈ · ± ρ c · ± ρ c
(cid:13) 0 0 (cid:13)L∞ (cid:13) 0 0 (cid:13)L∞
(cid:13) (cid:13) (cid:13) (cid:13)
(cid:13) (cid:13) (cid:13) (cid:13)
In this chapter, we propose to s(cid:13)tudy the stabilit(cid:13)y of the(cid:13)discretizations of schem(cid:13)es for which
some data are located on the faces or edges of the cells in the spirit of the MAC scheme
introduced in [15] for incompressible flows. The stability will be studied in terms of energy
dissipation and L -stability on the characteristic variables: two properties that are satisfied at
∞
the continuous level. Energy dissipation ensures that the numerical scheme satisfies a discrete
entropy inequality, while L -stability ensures that the scheme does not introduce oscillations.
∞
The study of L -stability allows us to glimpse the difficulty of defining variables depending on
∞
variables located at the faces and at the cells. As we will see, definitions of such variables will
have an impact on the properties obtained.
The chapter is organized as follow:
1) In section 3.2 the staggered discretization of interest is presented
2) In section 3.3, stability is addressed in the following senses: first, the von Neumann
necessary condition is considered; then, the discrete energy dissipation is investigated.

CHAPTER 3. STUDY OF THE STABILITY OF THE STAGGERED SCHEMES FOR
| THE | ONE | DIMENSIONAL |     |     | WAVE | SYSTEM |     |     |     |     |     |     | 39  |
| --- | --- | ----------- | --- | --- | ---- | ------ | --- | --- | --- | --- | --- | --- | --- |
3) Finally, in section 3.4 the L stability on the characteristic variables is studied.
∞
4) In section 3.5, the schemes are put in perspective with staggered discretizations arising
from the low Mach number asymptotic analysis of pre-existing schemes.
5) The theoretical results obtained are illustrated numerically in section 3.6.
6) These chapter’s findings will lay the ground for the extension in multiple dimension of
the discretizations on this very same system; we will thus sum them up in section 3.7.
| 3.2   |     | The       | staggered |         | schemes |     | on  | the | one | dimensional |     | wave | system |
| ----- | --- | --------- | --------- | ------- | ------- | --- | --- | --- | --- | ----------- | --- | ---- | ------ |
| 3.2.1 |     | Staggered |           | schemes |         |     |     |     |     |             |     |      |        |
In this section, we present the staggered discretizations of interest. The low Mach number
precision problem being intrinsic to multi-dimensions set ups [25, 24] , we will postpone until
the following chapter the introduction of the de Rham staggered setting. In more than one
space dimension our formalism with Raviart-Thomas will become relevant, here it is sufficient
to study the naive one dimensional staggered scheme: as in the case of the MAC scheme of [76],
the pressure unknowns are located at the cells, while the velocity unknowns are located at the
faces. We denote by p the value of p in cell C and by u the value of u at face x .
|       |     |        |     | i         |       | h    |     | i   |            | i+1/2 |       | h   | i+1/2 |
| ----- | --- | ------ | --- | --------- | ----- | ---- | --- | --- | ---------- | ----- | ----- | --- | ----- |
| Since | the | domain | is  | periodic, | we    | have |     |     |            |       |       |     |       |
|       |     |        |     | p         | p ,   |      |     |     | p          | p     | ,     |     |       |
|       |     |        |     | N         | = 0   |      |     |     |            | 1 =   | N 1   |     |       |
|       |     |        |     |           |       |      |     | and | −          |       | −     |     |       |
|       |     |        |     | u         | =     | u ,  |     |     | u          | =     | u     | .   |       |
|       |     |        |     |           | N+1/2 | 1/2  |     |     |            | 1/2   | N 1/2 |     |       |
|       |     |        |     | (cid:26)  |       |      |     |     | (cid:26) − |       | −     |     |       |
The general expression of the semi-discrete staggered scheme is given by
1
|     |     |     | ∆x∂    | p +     | u         |     | u     | = d c       | (p    | 2p  | +p    | ),      | (3.5a) |
| --- | --- | --- | ------ | ------- | --------- | --- | ----- | ----------- | ----- | --- | ----- | ------- | ------ |
|     |     |     |        | t i     | i+1/2     |     | i 1/2 | p           | 0 i+1 |     | i i   | 1       |        |
|     |     |     |        |         | ρ         | −   | −     |             |       | −   |       | −       |        |
|     |     |     |       |         | 0         |     |       |             |       |     |       |         |        |
|     |     |     |  ∆x∂ | u       | +(cid:0)κ | (p  | p )   | =(cid:1)d c | u     | 2u  | +u    | ,       | (3.5b) |
|     |     |     |        | t i+1/2 | 0         | i+1 | i     | u 0         | i+3/2 |     | i+1/2 | i 1/2   |        |
|     |     |     |        |         |           |     | −     |             |       | −   |       | −       |        |
|     |     |     |        |         |           |     |       | (cid:0)     |       |     |       | (cid:1) |        |
where d ,d  0 are stabilization terms. We not e that the original MAC scheme is fully-
|     |     | p u |     |     |     |     |     |     |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
≥
centered (FC) namely d = d = 0 but for a more general study, we propose to analyze also
|     |           |          |     | p   | u   |     |     |     |     |     |     |     |     |
| --- | --------- | -------- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| the | following | schemes: |     |     |     |     |     |     |     |     |     |     |     |
•
|     | fully-upwind        |     | (FU) | scheme: |         | d p = | d u > | 0;    |      |     |     |     |     |
| --- | ------------------- | --- | ---- | ------- | ------- | ----- | ----- | ----- | ---- | --- | --- | --- | --- |
|     | • pressure-centered |     |      | (PC)    | scheme: | d     | > 0   | and d | = 0; |     |     |     |     |
|     |                     |     |      |         |         |       | p     | u     |      |     |     |     |     |
•
|     | velocity-centered |     |     | (VC) | scheme | :   | d = 0 | and d | > 0. |     |     |     |     |
| --- | ----------------- | --- | --- | ---- | ------ | --- | ----- | ----- | ---- | --- | --- | --- | --- |
|     |                   |     |     |      |        |     | p     | u     |      |     |     |     |     |

CHAPTER 3. STUDY OF THE STABILITY OF THE STAGGERED SCHEMES FOR
THE ONE DIMENSIONAL WAVE SYSTEM 40
Three types of time integration are considered: explicit time integration
pn+1 pn 1
∆x i − i + un un = d c pn 2pn+pn , (3.6a)

∆t ρ
0
i+1/2− i
−
1/2 p 0 i+1− i i
−
1
(cid:16) (cid:17)


un+1 un (cid:0) (cid:1)
 ∆x i+1/2− i+1/2 +κ pn pn = d c un 2un +un , (3.6b)
∆t 0 i+1− i u 0 i+3/2− i+1/2 i 1/2
−
  (cid:0) (cid:1) (cid:16) (cid:17)

implicittime integration
pn+1 pn 1
∆x i − i + un+1 un+1 = d c pn+1 2pn+1+pn+1 ,
∆t ρ i+1/2− i 1/2 p 0 i+1 − i i 1
 0 − − (3.7)
   ∆x
un
i+
+
1
1
/2−
un
i+1/2 +
(cid:16)
κ pn+1 pn+
(cid:17)
1 = d c (cid:0) un+1 2un+1 + (cid:1) un+1 ,
∆t 0 i+1 − i u 0 i+3/2− i+1/2 i 1/2
−
  (cid:0) (cid:1) (cid:16) (cid:17)
and implicit-explicit time integration in the case of the fully-centered and pressure-centered
schemes (d = 0)
u
pn+1 pn 1
∆x i − i + un un = d c pn 2pn+pn ,
∆t ρ i+1/2− i 1/2 p 0 i+1− i i 1
 0 − − (3.8)
   ∆x
un
i+
+
1
1
/2−
un
i+1/2 +
(cid:16)
κ pn+1 pn+
(cid:17)
1 = 0. (cid:0) (cid:1)
∆t 0 i+1 − i
  (cid:0) (cid:1)
We note thatscheme (3.8) gives a fully explicit scheme in the sense that it is sufficient to start
by updating the pressure and then to update the velocity. This scheme (3.8) with d = 0 was
p
firstly proposed in [77] and we refer to [78] for a stability analysis.
3.2.2 Finite volume interpretation of the staggered scheme
Inonespacedimensionitispossibletolayabridgebetweenstaggereddiscretizationspresented
in the previous section and Godunov-type (collocated) methods with particular values of d
u
and d . Numerical scheme (3.5) can be seen as a finite volume scheme where the pressure p
p
is approximated by a constant per cell function p while the velocity u is approximated by
h
a constant per dual cell function u . This approach provides a natural upwinding for (3.5).
h
Denoting by p the value of p in cell C and by u the value of u in the dual cell C =
i h i i+1/2 h i+1/2
]x ;x [, the staggered scheme is obtained by
i i+1
1) First, writing a Godunov finite volume scheme on the half-cells x ;x and
i 1/2 i
−
x ;x : piecewise constant pressure and velocity are obtained on the half-cells. We
i i+1/2 (cid:3) (cid:2)
denote by (p ,u ) (resp. (p ,u )) the values of (p ,u ) in the half cell
(cid:3) (cid:2) i,L i 1/2,R i,R i+1/2,L h h
−
x ;x (resp. x ;x ).
i 1/2 i i i+1/2
−
2) (cid:3)Then, aver(cid:2)aging th(cid:3)e quantiti(cid:2)es on the half-cells to get a constant per cell pressure and a
constant per dual cell velocity.

CHAPTER 3. STUDY OF THE STABILITY OF THE STAGGERED SCHEMES FOR
THE ONE DIMENSIONAL WAVE SYSTEM 41
Integrating (3.1) on x ;x and using an explicit time integration, we obtain
i 1/2 i
−
(cid:3) (cid:2)
∆xpn
i,
+
L
1
−
pn
i + 1 u(cid:63) u(cid:63) = 0,
2 ∆t ρ i − i 1/2
 0 − (3.9)
   ∆ 2
xun
i −
+
1
1
/2,R ∆ − t
un
i − 1/2
(cid:16)
+κ 0 p(cid:63) i −
(cid:17)
p(cid:63) i 1/2 = 0,
−
 (cid:16) (cid:17)

where (p(cid:63),u(cid:63)) is the solution of the Riemann problem at face x between (pn,un ) and
i i i i i 1/2
(pn,un ) and (p(cid:63) ,u(cid:63) ) is the solution of the Riemann problem at face x −between
i i+1/2 i 1/2 i 1/2 i 1/2
(pn ,un ) and (−pn,un − ). −
i 1 i 1/2 i i 1/2
− − −
   u p(cid:63) i (cid:63) i = = p u n i n i − − 1/ 2 2 ρ + 2 1 0 c u 0 n i ( + u 1 n i / + 2 1 , /2− un i − 1/2 ), and   p u (cid:63) i (cid:63) i − − 1 1 / / 2 2 = = u pn i n i − − 1 1 2 / + 2− pn i ρ , 0 2 c 0 (pn i − pn i − 1 ),
 
sothat (3.9) gives
∆xpn
i,
+
L
1
−
pn
i + 1 un un = c 0 pn pn , (3.10a)
 2 ∆t 2ρ 0 i+1/2− i − 1/2 − 2 i − i − 1
 
   ∆ 2
xun
i −
+
1
1
/2,R ∆ − t
un
i − 1/2
(cid:16)
+ κ 2 0 pn i − pn i − 1
(cid:17)
= c 2 0 u
(cid:0)
n i+1/2− un i −
(cid:1)
1/2 . (3.10b)
   (cid:0) (cid:1) (cid:16) (cid:17)

Integrating now(3.1) on x ;x , we similarly obtain
i i+1/2
(cid:3) (cid:2)
∆xpn
i,
+
R
1
−
pn
i + 1 un un = c 0 pn pn , (3.11a)
 2 ∆t 2ρ i+1/2− i 1/2 2 i+1− i
0 −
 
  
∆xun
i+
+
1
1
/2,L−
un
i+1/2
(cid:16)
+ κ 0 pn pn
(cid:17)
= c 0
(cid:0)
un
(cid:1)
un . (3.11b)
2 ∆t 2 i+1− i − 2 i+1/2− i 1/2
−
   (cid:0) (cid:1) (cid:16) (cid:17)

Adding (3.10a) to (3.11a), we get
∆x
pn+1+pn+1
1 c
i,L i,R pn + un un = 0 pn 2pn+pn ,
∆t (cid:32) 2 − i (cid:33) ρ 0 i+1/2− i − 1/2 2 i+1− i i − 1
(cid:16) (cid:17)
(cid:0) (cid:1)
pn+1+pn+1
that corresponds to (3.6a) with d = 1/2 and pn+1 = i,L i,R . Moreover, adding the
p i 2
equation obtained by replacing i by i+1 in (3.10b) to (3.11b), we get
∆x
un
i+
+
1
1
/2,L
+un
i+
+
1
1
/2,R un +κ pn pn =
c
0 un 2un +un
∆t 2 − i+1/2 0 i+1− i 2 i+3/2− i+1/2 i 1/2
(cid:32) (cid:33) −
(cid:16) (cid:17)
(cid:0) (cid:1)
un+1 +un+1
that corresponds to (3.6b) with d = 1/2 and un+1 = i+1/2,L i+1/2,R .
u i+1/2 2

CHAPTER 3. STUDY OF THE STABILITY OF THE STAGGERED SCHEMES FOR
| THE ONE | DIMENSIONAL |     | WAVE     | SYSTEM |     |        |             | 42  |
| ------- | ----------- | --- | -------- | ------ | --- | ------ | ----------- | --- |
| 3.3     | von-Neumann |     | analysis |        | and | energy | dissipation |     |
The main objective of this chapter consists in studying L2 dissipation property of the staggered
schemes section 3.2 in order to guide our choice of numerical integrations and numerical dif-
fusions when extending the schemes in the multi-dimensional wave system. Since we aim at
dissipating the energy for long time convergence, we focus on mimicking this property in one
space dimension:
Definition 3.3.1 (Energy dissipation). A numerical scheme dissipates energy if for all Un, we
h
have
|     |     |     |     | E   | Un+1 | E(Un). |     |     |
| --- | --- | --- | --- | --- | ---- | ------ | --- | --- |
|     |     |     |     |     | h    | h      |     |     |
≤
| where the | energy | is defined | by  | (3.3). | (cid:0) (cid:1) |     |     |     |
| --------- | ------ | ---------- | --- | ------ | --------------- | --- | --- | --- |
von-NeumannstabilityanalysisisbasedonFourieranalysistostudyL2-stabilityofanumer-
ical scheme, which necessitates the strong assumption that the domain is periodic. It is inter-
esting to point out the periodicity since, one could argue that in such context the energy is con-
served (3.4) and dissipating it at the discrete level would be inconsistent with the model. How-
ever, we are more interested in studying eventually a wave system with boundary conditions
because of its relevance for low Mach number behaviour. The boundary conditions will lead to
adissipatingmechanism,itisthuscompletelysensibletoensurenumericalenergydissipationin
this (periodic) case. In general, numerical time integration that preserves at the discrete scale
the energy conservation (in the linear setting) would lead to a symplectic integration which is
| a complete | field | of research | and | is out | of the scope | of this | work. |     |
| ---------- | ----- | ----------- | --- | ------ | ------------ | ------- | ----- | --- |
Now, coming back to the numerical scheme, since p is constant by cell, as for classical
h
| collocated | scheme, | the kth | Fourier | coefficient | is  | given , | for k = 0, by |     |
| ---------- | ------- | ------- | ------- | ----------- | --- | ------- | ------------- | --- |
(cid:54)
|     |     |     |     | sin(kπ∆x)N |     | 1     |         |     |
| --- | --- | --- | --- | ---------- | --- | ----- | ------- | --- |
|     |     |     |     | p =        |     | − p e | 2jkπxi. |     |
|     |     |     |     | k          |     | i     | −       |     |
kπ
i=0
(cid:80)
n(cid:98)umber
where j denotes the imaginary such that j2 = 1. For u , we use the natural finite
|     |     |     |     |     |     |     | − h |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
kth
volume interpretation of subsection 3.2.2 on dual cell and then, the Fourier coefficient is
k
| given, | for = | 0, by |     |     |     |     |     |     |
| ------ | ----- | ----- | --- | --- | --- | --- | --- | --- |
(cid:54)
1
2jkπxdx
|     |     |     | u k | = u | h (x)e − |     |     |     |
| --- | --- | --- | --- | --- | -------- | --- | --- | --- |
(cid:90)0
|     |     |     |          | N 1 | xi+1  |         |     |     |
| --- | --- | --- | -------- | --- | ----- | ------- | --- | --- |
|     |     |     |          | −   |       | 2jkπxdx |     |     |
|     |     |     | (cid:98) | =   | u     | e       |     |     |
|     |     |     |          |     | i+1/2 | −       |     |     |
|     |     |     |          | i=0 | xi    |         |     |     |
(cid:90)
(cid:80)
|     |     |     |     | s in(kπ∆x)N |     | 1     |              |        |
| --- | --- | --- | --- | ----------- | --- | ----- | ------------ | ------ |
|     |     |     |     | =           | −   | u e   | 2jkπx i+1/2. | (3.12) |
|     |     |     |     |             |     | i+1/2 | −            |        |
kπ
i=0
(cid:80)
Then,usinganexplicitoranimplicittimeintegration,thenumericalschemescanbewritten
| in term | of Fourier | coefficients | as  |     |                |     |     |     |
| ------- | ---------- | ------------ | --- | --- | -------------- | --- | --- | --- |
|         |            |              |     |     | Un+1 = A(k)Un, |     |     |     |
k
k
|     |     |     |     |     | (cid:98) | (cid:98) |     |     |
| --- | --- | --- | --- | --- | -------- | -------- | --- | --- |

CHAPTER 3. STUDY OF THE STABILITY OF THE STAGGERED SCHEMES FOR
| THE | ONE | DIMENSIONAL |     | WAVE |     | SYSTEM |     |     |     |     | 43  |
| --- | --- | ----------- | --- | ---- | --- | ------ | --- | --- | --- | --- | --- |
pn
|     | A(k) |     |     |     |     |     | Un  | k   |     |     |     |
| --- | ---- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
where is the amplification matrix and = . In order to prove the discrete
|        |     |             |        |     |           |              | k   | un         |          |     |     |
| ------ | --- | ----------- | ------ | --- | --------- | ------------ | --- | ---------- | -------- | --- | --- |
|        |     |             |        |     |           |              |     | (cid:18) k | (cid:19) |     |     |
| energy |     | dissipation | we use | the | following | propositions |     | (cid:98)   |          |     |     |
(cid:98)
(cid:98)
Proposition 3.3.1 (energy dissipation). Denoting by A(k) the amplification matrix of the kth
Fourier coefficient and assuming that for all k A( k) = A(k), then the following properties are
−
equivalent
|     | 1. the | scheme | is  | dissipative; |     |     |     |     |     |     |     |
| --- | ------ | ------ | --- | ------------ | --- | --- | --- | --- | --- | --- | --- |
energy
2.
|     |     |     |     |     | k   |     | S1/2A(k)S | 1/2 |     | 1,  | (3.13) |
| --- | --- | --- | --- | --- | --- | --- | --------- | --- | --- | --- | ------ |
−
|     |        |          |            |          | ∀     |                          |     |                   | 2 ≤                      |          |     |
| --- | ------ | -------- | ---------- | -------- | ----- | ------------------------ | --- | ----------------- | ------------------------ | -------- | --- |
|     |        |          |            |          |       | (cid:12)(cid:12)(cid:12) |     |                   | (cid:12)(cid:12)(cid:12) |          |     |
|     | where  |          | is the     | operator | norm; | (cid:12)(cid:12)(cid:12) |     |                   | (cid:12)(cid:12)(cid:12) |          |     |
|     |        | |||·|||2 |            |          |       | (cid:12)(cid:12)(cid:12) |     |                   | (cid:12)(cid:12)(cid:12) |          |     |
|     | 3. for | all k    | the matrix |          |       |                          |     |                   |                          |          |     |
|     |        |          |            |          | I     | 1                        | 1   | 1                 |                          | 1        |     |
|     |        |          |            |          |       | S2 A(k)S                 | 2   | ∗ S 2A(k)S        |                          | 2        |     |
|     |        |          |            |          | 2 −   |                          | −   |                   | −                        |          |     |
|     |        |          |            |          |       | (cid:16)                 |     | (cid:17) (cid:16) |                          | (cid:17) |     |
T
is positive, where the operator corresponds to the conjugate transpose, i.e. M = M .
|     |     |     |     |     |     | ∗   |     |     |     |     | ∗   |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
·
Proof. We firstly prove that point 2 implies point 1. From Parseval identity, we get that the
| energy |     | (3.3) is | given by |     |       |            |               |                |               |     |     |
| ------ | --- | -------- | -------- | --- | ----- | ---------- | ------------- | -------------- | ------------- | --- | --- |
|        |     |          |          |     |       | 1          |               | 1              |               |     |     |
|        |     |          |          |     | E(U ) | =          | p 2 +         | u              | 2             |     |     |
|        |     |          |          |     | h     |            | h (cid:107)L2 |                | h (cid:107)L2 |     |     |
|        |     |          |          |     |       | 2(cid:107) |               | 2ρ2c2(cid:107) |               |     |     |
|        |     |          |          |     |       |            |               | 0 0            |               |     |     |
|        |     |          |          |     |       | 1          |               | 1              |               |     |     |
|        |     |          |          |     |       | =          | p 2+          |                | u 2           |     |     |
|        |     |          |          |     |       |            | k             | 2ρ2c2          | k             |     |     |
|        |     |          |          |     |       | 2          | | |           |                | | |           |     |     |
|        |     |          |          |     |       |            | k             | 0 0 k          |               |     |     |
|        |     |          |          |     |       | (cid:80)   |               | (cid:80)       |               |     |     |
|        |     |          |          |     |       |            | 1             | 2              |               |     |     |
|        |     |          |          |     |       | =          | S(cid:98)2U   | ,              | (cid:98)      |     |     |
k
|     |     |     |     |     |     | k        |          | 2        |     |     |     |
| --- | --- | --- | --- | --- | --- | -------- | -------- | -------- | --- | --- | --- |
|     |     |     |     |     |     | (cid:80) | (cid:13) | (cid:13) |     |     |     |
where we recall that S is the symmetriz er o(cid:13)f th(cid:98)e w(cid:13)ave system. Then, we have
|     |     |     |     |        |     |     | (cid:13) | (cid:13) |     |     |     |
| --- | --- | --- | --- | ------ | --- | --- | -------- | -------- | --- | --- | --- |
|     |     |     |     |        |     |     | 1        | 2        |     |     |     |
|     |     |     |     | E Un+1 | =   |     | S 2Un+1  |          |     |     |     |
|     |     |     |     |        | h   |     | k        |          |     |     |     |
2
k
|     |     |     |     | (cid:0) | (cid:1) | (cid:13)         |                | (cid:13)   |     |     |     |
| --- | --- | --- | --- | ------- | ------- | ---------------- | -------------- | ---------- | --- | --- | --- |
|     |     |     |     |         |         | (cid:80)(cid:13) | 1              | (cid:13) 2 |     |     |     |
|     |     |     |     |         | =       | (cid:13)S        | 2A(cid:98)(k)U | (cid:13)n  |     |     |     |
k
2
k
|     |     |     |     |     |     | (cid:80)(cid:13) (cid:13) |               | (cid:13)                 |                |          |     |
| --- | --- | --- | --- | --- | --- | ------------------------- | ------------- | ------------------------ | -------------- | -------- | --- |
|     |     |     |     |     |     |                           | 1             | 1(cid:13) 1              | 2              |          |     |
|     |     |     |     |     | =   | (cid:13)S2                | A(k)S(cid:98) | (cid:13)S 2Un            |                |          |     |
|     |     |     |     |     |     |                           |               | −2                       | k              |          |     |
|     |     |     |     |     |     | k                         |               |                          | 2              |          |     |
|     |     |     |     |     |     | (cid:80) (cid:13)         |               |                          | (cid:13)       |          |     |
|     |     |     |     |     |     | (cid:13)                  | 1             | 1 (cid:98) 2             | (cid:13) 1     | 2        |     |
|     |     |     |     |     |     | (cid:13)                  | S2 A(k)S      | 2                        | S (cid:13) 2Un | ,        |     |
|     |     |     |     |     |     |                           |               | −                        |                | k        |     |
|     |     |     |     |     | ≤   | k                         |               | 2                        |                | 2        |     |
|     |     |     |     |     |     | (cid:12)(cid:12)(cid:12)  |               | (cid:12)(cid:12)(cid:12) | (cid:13)       | (cid:13) |     |
(cid:80)(cid:12)(cid:12)(cid:12)
|      |       |     |     |     |     |                          |     | (cid:12)(cid:12)(cid:12) | (cid:13) (cid:98) | (cid:13) |     |
| ---- | ----- | --- | --- | --- | --- | ------------------------ | --- | ------------------------ | ----------------- | -------- | --- |
| that | gives |     |     |     |     | (cid:12)(cid:12)(cid:12) |     | (cid:12)(cid:12)(cid:12) | (cid:13)          | (cid:13) |     |
2
|     |           |      |            |     | Un +1   |                          | 1 Un     | E(Un     |     |     |     |
| --- | --------- | ---- | ---------- | --- | ------- | ------------------------ | -------- | -------- | --- | --- | --- |
|     |           |      |            | E   |         |                          | S2       | =        | ),  |     |     |
|     |           |      |            |     | h       | ≤                        |          | k        | h   |     |     |
|     |           |      |            |     |         |                          | k        | 2        |     |     |     |
|     |           |      |            |     |         |                          | (cid:13) | (cid:13) |     |     |     |
|     |           |      |            |     | (cid:0) | (cid:1) (cid:80)(cid:13) |          |          |     |     |     |
| if  | the point | 2 is | satisfied. |     |         |                          | (cid:98) | (cid:13) |     |     |     |
|     |           |      |            |     |         |                          | (cid:13) | (cid:13) |     |     |     |

CHAPTER 3. STUDY OF THE STABILITY OF THE STAGGERED SCHEMES FOR
THE ONE DIMENSIONAL WAVE SYSTEM 44
We now prove that point 1 implies point 3. From the previous computation, we have
E Un h +1 − E(Un h ) = k S2 1 A(k)S −2 1 S 1 2Un k 2 2− k S 1 2Un k 2 2
(cid:13) (cid:13) (cid:13) (cid:13)
(cid:0) (cid:1) = (cid:80)(cid:13) (cid:13)S2 1 Un k ∗ S 1 2A(cid:98)(k (cid:13) (cid:13))S −2 1(cid:80) ∗ (cid:13) (cid:13)S 1 2 (cid:98)A(k (cid:13) (cid:13))S − 1 2 − I 2 S 1 2Un k .
k
(cid:16) (cid:17) (cid:16)(cid:16) (cid:17) (cid:16) (cid:17) (cid:17)(cid:16) (cid:17)
(cid:80)
(cid:98) (cid:98)
Let k and Y C2. Taking for Un the vector Un(k ) defined on all cell C by
0 ∈ h h 0 i
Un
i
(k
0
) = 2cos(2k
0
πx
i
)
sin
k
(
0
k
π∆
π∆
x
x)
S
−
1 2Y,
0
the computation of the kth Fourier coefficient gives
S 1 2Un
k
= Y
0,
, i
e
f
ls
k
ew
=
h
−
er
k
e
0
,
or k = k 0 ,
(cid:26)
(cid:98)
so that E Un+1 E(Un) 0 gives
h − h ≤
(cid:0) (cid:1)
Y ∗ S2 1 A( k 0 )S − 1 2 ∗ S 1 2A( k 0 )S − 1 2 I 2 Y
− − −
(cid:16)(cid:16) +Y ∗ S2 1(cid:17) A( (cid:16) k 0 )S − 1 2 ∗ S 1 2A (cid:17) (k 0 )S − (cid:17)1 2 I 2 Y 0.
− ≤
(cid:16)(cid:16) (cid:17) (cid:16) (cid:17) (cid:17)
Since A( k ) = A(k ), we have
0 0
−
Y ∗ S 1 2A( k 0 )S − 1 2 ∗ S 1 2A( k 0 )S − 1 2 Y = Y ∗ S 1 2A(k 0 )S −2 1 ∗ S 1 2A(k 0 )S −2 1 Y
− −
(cid:16) (cid:17) (cid:16) (cid:17) (cid:16) (cid:17) (cid:16) (cid:17)
= Y ∗ S 1 2A(k 0 )S − 1 2 ∗ S 1 2A(k 0 )S − 1 2 Y
(cid:16) (cid:17) (cid:16) (cid:17)
= Y ∗ S 1 2A(k 0 )S − 1 2 ∗ S 1 2A(k 0 )S −2 1 Y
= Y ∗ (cid:16) S 1 2A(k 0 )S − 1 2 (cid:17) ∗ (cid:16) S 1 2A(k 0 )S −2 1(cid:17) Y,
(cid:16) (cid:17) (cid:16) (cid:17)
so that
Y ∗ S2 1 A(k 0 )S − 1 2 ∗ S 1 2A(k 0 )S − 1 2 I 2 Y 0,
− ≤
(cid:16)(cid:16) (cid:17) (cid:16) (cid:17) (cid:17)
for all Y C2. Then, the matrix I 2 S 1 2A(k)S −2 1 ∗ S 1 2A(k)S −2 1 is positive.
∈ −
We now prove that point 3 implies(cid:16)point 2. We(cid:17)rec(cid:16)all that for a(cid:17)ll matrix M, for ρ(M) its
spectral radius, we have
M 2 = ρ(M M).
||| |||2 ∗
Since the matrix S2 1 A(k)S − 1 2 ∗ S 1 2A(k)S − 1 2 is Hermitian, it is diagonalizable with real
(cid:16) (cid:17) (cid:16) (cid:17)

CHAPTER 3. STUDY OF THE STABILITY OF THE STAGGERED SCHEMES FOR
| THE          | ONE | DIMENSIONAL |     |      | WAVE     | SYSTEM |      |          |       |          |     | 45  |
| ------------ | --- | ----------- | --- | ---- | -------- | ------ | ---- | -------- | ----- | -------- | --- | --- |
| eigenvalues. |     | Then,       | we  | have |          |        |      |          |       |          |     |     |
|              |     |             |     | I    | 1        |        | 1    | 1        | 1     |          |     |     |
|              |     |             |     |      | S 2A(k)S |        | −2 ∗ | S 2A(k)S | −2 is | positive |     |     |
2 −
|     |     |     |     |                  | (cid:16)1 |     | (cid:17) (cid:16)1 |        | (cid:17)         |     |     |     |
| --- | --- | --- | --- | ---------------- | --------- | --- | ------------------ | ------ | ---------------- | --- | --- | --- |
|     |     |     |     |                  |           | 1   | ∗                  |        | 1                |     |     |     |
|     |     |     |     | ρ                | S 2A(k)S  | − 2 | S                  | 2A(k)S | − 2              | 1   |     |     |
|     |     |     |     | ⇒                |           |     |                    |        | ≤                |     |     |     |
|     |     |     |     | (cid:16)(cid:16) |           |     | (cid:17)2 (cid:16) |        | (cid:17)(cid:17) |     |     |     |
|     |     |     |     |                  | 1         | 1   |                    |        |                  |     |     |     |
|     |     |     |     | S2A(k)S          |           | −2  | 1                  |        |                  |     |     |     |
|     |     |     |     | ⇒                |           |     | ≤                  |        |                  |     |     |     |
2
|       |           |     |                                | (cid:12)(cid:12)(cid:12) |     | (cid:12)(cid:12)(cid:12) |     |     |     |     |     |     |
| ----- | --------- | --- | ------------------------------ | ------------------------ | --- | ------------------------ | --- | --- | --- | --- | --- | --- |
| which | concludes | the | proof(cid:12).(cid:12)(cid:12) | (cid:12)(cid:12)(cid:12) |     | (cid:12)(cid:12)(cid:12) |     |     |     |     |     |     |
(cid:12)(cid:12)(cid:12)
|     | Then, | following | [79, | 80] |     |     |     |     |     |     |     |     |
| --- | ----- | --------- | ---- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
Definition 3.3.2 (von–Neumann necessary stability condition). If A(k) denotes the amplifica-
| tion | matrix       | of the | kth Fourier |     | coefficient, |     | then    | the         |     |                     |           |        |
| ---- | ------------ | ------ | ----------- | --- | ------------ | --- | ------- | ----------- | --- | ------------------- | --------- | ------ |
|      |              |        |             |     |              |     |         | von–Neumann |     | necessary stability | condition |        |
| for  | L2-stability | is     | defined     | as  |              |     |         |             |     |                     |           |        |
|      |              |        |             |     |              | max | ρ(A(k)) |             | 1   |                     |           | (3.14) |
≤
k
| where | ρ(A) | is the | spectral | radius | of  | the matrix |     | A.  |     |     |     |     |
| ----- | ---- | ------ | -------- | ------ | --- | ---------- | --- | --- | --- | --- | --- | --- |
and,
Proposition 3.3.2. If a numerical scheme ensures (3.13), then it ensures Definition 3.3.2.
Proof.
|     | The | property                 | follows   | from | the                      | inequality |           |     |          |                |     |     |
| --- | --- | ------------------------ | --------- | ---- | ------------------------ | ---------- | --------- | --- | -------- | -------------- | --- | --- |
|     |     |                          | S1/2A(k)S |      | 1/2                      |            | S1/2A(k)S |     | 1/2      |                |     |     |
|     |     |                          |           |      |                          | max        | ρ         |     |          | = max ρ(A(k)). |     |     |
|     |     |                          |           | −    |                          |            |           |     | −        |                |     |     |
|     |     |                          |           |      | 2                        | ≥ k        |           |     |          | k              |     |     |
|     |     | (cid:12)(cid:12)(cid:12) |           |      | (cid:12)(cid:12)(cid:12) |            | (cid:16)  |     | (cid:17) |                |     |     |
|     |     | (cid:12)(cid:12)(cid:12) |           |      | (cid:12)(cid:12)(cid:12) |            |           |     |          |                |     |     |
|     |     | (cid:12)(cid:12)(cid:12) |           |      | (cid:12)(cid:12)(cid:12) |            |           |     |          |                |     |     |
Proposition 3.3.3. If a numerical scheme ensures the von–Neumann necessary stability con-
dition (3.3.2) and for all k, the matrix S1/2A(k)S 1/2 is normal then it is energy dissipative.
−
|        |       |         | S1/2A(k)S |     |     | 1/2        |     |         |     |     |     |     |
| ------ | ----- | ------- | --------- | --- | --- | ---------- | --- | ------- | --- | --- | --- | --- |
| Proof. | Since | for all | k,        |     |     | is normal, |     | we have |     |     |     |     |
−
|     |     |     |     |     | S1/2A(k)S |     | 1/2 | =   | ρ(A(k)), |     |     |     |
| --- | --- | --- | --- | --- | --------- | --- | --- | --- | -------- | --- | --- | --- |
−
2
|         |     |          |     |     | (cid:12)(cid:12)(cid:12) |           |     | (cid:12)(cid:12)(cid:12) |     |     |     |     |
| ------- | --- | -------- | --- | --- | ------------------------ | --------- | --- | ------------------------ | --- | --- | --- | --- |
|         |     |          |     |     | (cid:12)(cid:12)(cid:12) |           |     | (cid:12)(cid:12)(cid:12) |     |     |     |     |
| so that | by  | (3.3.2), |     |     | (cid:12)(cid:12)(cid:12) |           |     | (cid:12)(cid:12)(cid:12) |     |     |     |     |
|         |     |          |     |     |                          | S1/2A(k)S |     | 1/2                      | 1,  |     |     |     |
−
2 ≤
|       |         |         |             |           |                                       | (cid:12)(cid:12)(cid:12) |     | (cid:12)(cid:12)(cid:12)       |     |     |     |     |
| ----- | ------- | ------- | ----------- | --------- | ------------------------------------- | ------------------------ | --- | ------------------------------ | --- | --- | --- | --- |
| which | implies | energy  | dissipation |           | by(cid:12)(cid:12)(cid:12)Proposition |                          |     | 3.3.(cid:12)1(cid:12)(cid:12). |     |     |     |     |
|       |         |         |             |           |                                       | (cid:12)(cid:12)(cid:12) |     | (cid:12)(cid:12)(cid:12)       |     |     |     |     |
| 3.3.1 | von     | Neumann |             | stability |                                       | analysis                 |     |                                |     |     |     |     |
For each time integration, we study the behaviour of the Fourier modes.

CHAPTER 3. STUDY OF THE STABILITY OF THE STAGGERED SCHEMES FOR
THE ONE DIMENSIONAL WAVE SYSTEM 46
Explicit time integration
The staggered explicit scheme (3.6) can be written in terms of Fourier coefficients as
n+1 n
p p
k = A(k) k ,
u u
k k
(cid:18) (cid:19) (cid:18) (cid:19)
(cid:98) (cid:98)
where the amplification matrix A(k) is given by
(cid:98) (cid:98)
c ∆t ∆t
1 4d 0 sin2(kπ∆x) j2 sin(kπ∆x)
p
AStag,Exp(k) = − ∆x − ρ 0 ∆x . (3.15)
 κ ∆t c ∆t 
j2 0 sin(kπ∆x) 1 4d 0 sin2(kπ∆x)
u
 − ∆x − ∆x 
 
Proposition 3.3.4 (staggered explicit schemes and von–Neumann necessary condition). The
staggered explicit scheme (3.6) satisfies the von–Neumann necessary condition (3.14) under the
CFL condition
d +d
p u
, if d d 1,
c 0 ∆t 1+4d p d u | p − u | ≤
0 
1
≤ ∆x ≤  , elsewhere.

 d +d + (d d )2 1
p u p u
− −

Then, we obtain for   (cid:112)
• the staggered explicit fully-upwind (d = d > 0)
p u
c ∆t 2d
0 p
0 (3.16)
≤ ∆x ≤ 1+4d2
p
and so the CFL condition is optimal for d = d = 1/2.
p u
• the staggered explicit pressure-centered (d > 0 and d = 0)
p u
d , if d 1,
p p
c 0 ∆t 1 ≤
0 ≤ ∆x ≤  , elsewhere (3.17)
 d + d2 1
 p p
−
(cid:113)

and so the CFL condition is optimal for d = 1.
p
• the staggered explicit velocity-centered (d = 0 and d > 0)
p u
d , if d 1,
0 c 0 ∆t u 1 u ≤ (3.18)
≤ ∆x ≤  , elsewhere
 d u + d2 u 1
−
and so the CFL condition is optimal for d (cid:112) = 1.
u

CHAPTER 3. STUDY OF THE STABILITY OF THE STAGGERED SCHEMES FOR
THE ONE DIMENSIONAL WAVE SYSTEM 47
Proof. Let us study the eigenvalues of the following matrix (3.15):
c ∆t ∆t
1 4d 0 sin2(kπ∆x) j2 sin(kπ∆x)
p
AStag,Exp(k) = − ∆x − ρ 0 ∆x .
 κ ∆t c ∆t 
j2 0 sin(kπ∆x) 1 4d 0 sin2(kπ∆x)
u
 − ∆x − ∆x 
 
The trace of the matrix is
c ∆t
Tr(AStag,Exp(k)) = 2 1 2(d +d ) 0 sin2(kπ∆x) .
p u
− ∆x
(cid:18) (cid:19)
Its determinant is
c ∆t c ∆t
det(AStag,Exp(k))= 1 4d 0 sin2(kπ∆x) 1 4d 0 sin2(kπ∆x)
p u
− ∆x − ∆x
(cid:18) (cid:19)(cid:18) (cid:19)
c2(∆t)2
+4 0 sin2(kπ∆x)
(∆x)2
c ∆t 16d d c2(∆t)2
=1 4(d +d ) 0 sin2(kπ∆x)+ p u 0 sin4(kπ∆x)
− p u ∆x (∆x)2
c2(∆t)2
+4 0 sin2(kπ∆x).
(∆x)2
The reduced determinant of the characteristic polynomial is
c ∆t 2 c ∆t
∆ = 1 2(d +d ) 0 sin2(kπ∆x) 1+4(d +d ) 0 sin2(kπ∆x)
(cid:48) p u p u
− ∆x − ∆x
(cid:18) (cid:19)
16d d c2(∆t)2 c2(∆t)2
p u 0 sin4(kπ∆x) 4 0 sin2(kπ∆x)
− (∆x)2 − (∆x)2
c ∆t c2(∆t)2
=1 4(d +d ) 0 sin2(kπ∆x)+4(d +d )2 0 sin4(kπ∆x)
− p u ∆x p u (∆x)2
c ∆t
1+4(d +d ) 0 sin2(kπ∆x)
p u
− ∆x
16d d c2(∆t)2 c2(∆t)2
p u 0 sin4(kπ∆x) 4 0 sin2(kπ∆x)
− (∆x)2 − (∆x)2
c2(∆t)2
=4 0 sin2(kπ∆x) (d d )2sin2(kπ∆x) 1 .
(∆x)2 p − u −
(cid:16) (cid:17)
Case d d 1 : If d d 1, then the reduced discriminant is always negative, and
p u p u
| − | ≤ | − | ≤
the eigenvalues are
c ∆t
λStag,Exp =1 2(d +d ) 0 sin2(kπ∆x)
p u
± − ∆x
c ∆t
2j 0 sin(kπ∆x) 1 (d d )2sin2(kπ∆x).
p u
± ∆x | | − −
(cid:113)

CHAPTER 3. STUDY OF THE STABILITY OF THE STAGGERED SCHEMES FOR
THE ONE DIMENSIONAL WAVE SYSTEM 48
The module of these eigenvalues is
2 c ∆t 2
λStag,Exp = 1 2(d +d ) 0 sin2(kπ∆x)
p u
± − ∆x
(cid:12) (cid:12) (cid:18) c2(∆t)2 (cid:19)
(cid:12) (cid:12) (cid:12) (cid:12) +4 ( 0 ∆x)2 sin2(kπ∆x) 1 − (d p − d u )2sin2(kπ∆x)
c ∆t
(cid:16) c2(∆t)2(cid:17)
=1 4(d +d ) 0 sin2(kπ∆x)+4(d +d )2 0 sin4(kπ∆x)
− p u ∆x p u (∆x)2
c2(∆t)2
+4 0 sin2(kπ∆x) 1 (d d )2sin2(kπ∆x) .
(∆x)2 − p − u
(cid:16) (cid:17)
2
Then
λStag,Exp
1 0 if and only if
± − ≤
(cid:12) (cid:12)
(cid:12) (cid:12)
(cid:12) (cid:12) c ∆t c2(∆t)2
4(d +d ) 0 sin2(kπ∆x)+4(d +d )2 0 sin4(kπ∆x)
− p u ∆x p u (∆x)2
c2(∆t)2
+4 0 sin2(kπ∆x) 1 (d d )2sin2(kπ∆x) 0.
(∆x)2 − p − u ≤
(cid:16) (cid:17)
4c ∆t
In this last expression, we simplify by 0 sin2(kπ∆x) for obtaining
∆x
c ∆t c ∆t
(d +d )+(d +d )2 0 sin2(kπ∆x)+ 0 1 (d d )2sin2(kπ∆x) 0.
p u p u p u
− ∆x ∆x − − ≤
(cid:16) (cid:17)
This inequality is equivalent to
c ∆t
0 1+4d d sin2(kπ∆x) d +d ,
p u p u
∆x ≤
(cid:0) (cid:1)
which gives the following CFL condition
c ∆t d +d
0 p u .
∆x ≤ 1+4d d
p u
Case d d 1: In this case, the reduced discriminant is negative for small k, whereas
p u
| − | ≥
it is positive for the k such that sin2(kπ∆x) is close to 1.
When the discriminant is negative: Forsmallk,thestudyledintheprevioussectioncanbe
done, except for the upper bound on sin2(kπ∆x), which is no more 1, but rather 1/(d d )2,
p u
−
which gives the following first stability condition
c ∆t d +d (d +d )(d d )2 (d d )2
0 p u = p u p − u = p − u .
∆x ≤ 1+ (d 4 p dp d d u u )2 (d p +d u )2 d p +d u
−
When the discriminant is positive: Fork suchthatsin2(kπ∆x)isgreaterthan1/(d d )2,
p u
−

CHAPTER 3. STUDY OF THE STABILITY OF THE STAGGERED SCHEMES FOR
THE ONE DIMENSIONAL WAVE SYSTEM 49
the reduced discriminant is positive, and the eigenvalues are real. The eigenvalues are
c ∆t
λStag,Exp =1 2(d +d ) 0 sin2(kπ∆x)
p u
± − ∆x
c ∆t
2 0 sin(kπ∆x) (d d )2sin2(kπ∆x) 1.
p u
± ∆x | | − −
(cid:113)
• Investigation on λStag,Exp 1; we first investigate conditions for ensuring λStag,Exp
− ≥ − − ≥
1, which read
−
c ∆t c ∆t
2 2(d +d ) 0 sin2(kπ∆x) 2 0 sin(kπ∆x) (d d )2sin2(kπ∆x) 1 0,
p u p u
− ∆x − ∆x | | − − ≥
(cid:113)
which can be rewritten
c ∆t
1 0 (d +d ) sin2(kπ∆x)+ sin(kπ∆x) (d d )2sin2(kπ∆x) 1 0,
p u p u
− ∆x | | − − ≥
(cid:18) (cid:113) (cid:19)
or
c ∆t 1
0
∆x ≤ (d +d ) sin2(kπ∆x)+ sin(kπ∆x) (d d )2sin2(kπ∆x) 1
p u p u
| | − −
(cid:113)
because the denominator is positive. We therefore want to study the maximum of the
function
X (d +d )X +√X (d d )2X 1
p u p u
(cid:55)→ − −
(cid:113)
1
for X ;1 . As the function is increasing, the maximum is reached in 1,
∈ (d d )2
p u
(cid:20) − (cid:21)
which gives the stability condition
c ∆t 1
0 .
∆x ≤ (d +d )+ (d d )2 1
p u p u
− −
(cid:113)
• Investigation on λStag,Exp 1 ; we now investigate conditions for ensuring λStag,Exp 1,
+ ≤ + ≤
which read
c ∆t c ∆t
2(d +d ) 0 sin2(kπ∆x)+2 0 sin(kπ∆x) (d d )2sin2(kπ∆x) 1 0,
p u p u
− ∆x ∆x | | − − ≤
(cid:113)
it can be simplified as
sin(kπ∆x) (d d )2sin2(kπ∆x) 1 (d +d )sin2(kπ∆x).
p u p u
| | − − ≤
(cid:113)
Taking the square and simplifying by sin2(kπ∆x) gives
(d d )2sin2(kπ∆x) 1 (d +d )2sin2(kπ∆x),
p u p u
− − ≤

CHAPTER 3. STUDY OF THE STABILITY OF THE STAGGERED SCHEMES FOR
| THE ONE | DIMENSIONAL |              |     | WAVE | SYSTEM |     |            |     |     |     |     | 50  |
| ------- | ----------- | ------------ | --- | ---- | ------ | --- | ---------- | --- | --- | --- | --- | --- |
| which   | may         | be rewritten |     |      |        |     |            |     |     |     |     |     |
|         |             |              |     |      | 1+4d   | d   | sin2(kπ∆x) |     | 0,  |     |     |     |
p u
≥
which is always ensured. This means that we always have λStag,Exp 1.
+
≤
Conclusion: Ensuring the necessary von-Neumann condition, requires
|     |     | c ∆t |     |          | (d   | d )2 |     |      |       |           |     |     |
| --- | --- | ---- | --- | -------- | ---- | ---- | --- | ---- | ----- | --------- | --- | --- |
|     |     | 0    |     |          | p    | u    |     | 1    |       |           |     |     |
|     |     |      |     | min      | −    | ,    |     |      |       | .         |     |     |
|     |     | ∆x   |     |          | d +d |      |     |      |       |           |     |     |
|     |     |      | ≤   | (cid:32) | p    | u d  | +d  | + (d | d )2  | 1(cid:33) |     |     |
|     |     |      |     |          |      |      | p u |      | p − u | −         |     |     |
(cid:112)
| The minimum | between |     | the | two is | quite | easy to | determine. |     | Indeed, |     |     |     |
| ----------- | ------- | --- | --- | ------ | ----- | ------- | ---------- | --- | ------- | --- | --- | --- |
|             |         |     |     | d +d   | d     | +d      | (d         | d   | )2 1,   |     |     |     |
|             |         |     |     | p      | u     | p u     | +          | p u |         |     |     |     |
|             |         |     |     |        | ≤     |         |            | −   | −       |     |     |     |
(cid:113)
which gives
|     |     |     |     | 1   |     |     | 1   |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
.
|             |     |     |       | d p +d | u ≥ d | +d  | + (d      | d     | )2 1 |     |     |     |
| ----------- | --- | --- | ----- | ------ | ----- | --- | --------- | ----- | ---- | --- | --- | --- |
|             |     |     |       |        |       | p u |           | p − u | −    |     |     |     |
| We multiply | by  | (d  | d )2, | which  | gives |     | (cid:112) |       |      |     |     |     |
|             |     | p   | u     |        |       |     |           |       |      |     |     |     |
−
|     |     |     | (d  | d    | )2  |      | (d  | d )2 |      |     |     |     |
| --- | --- | --- | --- | ---- | --- | ---- | --- | ---- | ---- | --- | --- | --- |
|     |     |     |     | p    | u   |      | p   | u    |      | ,   |     |     |
|     |     |     |     | −    |     |      | −   |      |      |     |     |     |
|     |     |     |     | d +d | ≥   | d +d | +   | (d d | )2 1 |     |     |     |
|     |     |     |     | p    | u   | p    | u   | p    | u    |     |     |     |
|     |     |     |     |      |     |      |     | −    | −    |     |     |     |
(cid:112)
| and as | we are in | the case | in  | which | (d  | d )2 | 1, we | find |     |     |     |     |
| ------ | --------- | -------- | --- | ----- | --- | ---- | ----- | ---- | --- | --- | --- | --- |
|        |           |          |     |       | p   | u    |       |      |     |     |     |     |
|        |           |          |     |       | −   | ≥    |       |      |     |     |     |     |
|        |           |          | (d  | d     | )2  |      |       | 1    |     |     |     |     |
p u
|            |      |                 |     | −    |           |      |                     |         |      | .   |     |     |
| ---------- | ---- | --------------- | --- | ---- | --------- | ---- | ------------------- | ------- | ---- | --- | --- | --- |
|            |      |                 |     | d +d | ≥         | d +d | +                   | (d d    | )2 1 |     |     |     |
|            |      |                 |     | p    | u         | p    | u                   | p       | u    |     |     |     |
|            |      |                 |     |      |           |      |                     | −       | −    |     |     |     |
| This means | that | the von-Neumann |     |      | necessary |      | condition (cid:112) | in this | case | is  |     |     |
c ∆t
|     |     |     |     | 0   |     |     | 1   |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
.
∆x
|     |     |     |     |     | ≤ d | +d + | (d        | d )2 | 1   |     |     |     |
| --- | --- | --- | --- | --- | --- | ---- | --------- | ---- | --- | --- | --- | --- |
|     |     |     |     |     | p   | u    | p         | − u  | −   |     |     |     |
|     |     |     |     | d   |     | d    | (cid:112) |      |     | d   | d   |     |
This expression is symmetric in p and u . Also , we observe that for p 0 and u 0, the
|        |         |     |        |         |      |     |          |     |     | ≥   | ≥   |     |
| ------ | ------- | --- | ------ | ------- | ---- | --- | -------- | --- | --- | --- | --- | --- |
| region | (d d )2 | 1   | can be | divided | into | two | regions: |     |     |     |     |     |
p u
− ≥
|               |      |     |          | d      |      |        |         | d          |         |     |     |     |
| ------------- | ---- | --- | -------- | ------ | ---- | ------ | ------- | ---------- | ------- | --- | --- | --- |
|               |      |     |          | p      | 0    |        |         | u          | 0       |     |     |     |
|               |      |     |          | ≥      |      |        | or      | ≥          |         |     |     |     |
|               |      |     |          | d      | d +1 |        |         | d          | d 1     |     |     |     |
|               |      |     |          | u      | p    |        |         | u          | p       |     |     |     |
|               |      |     | (cid:26) | ≥      |      |        |         | (cid:26) ≤ | −       |     |     |     |
| If we suppose | that | d   | d        | +1 and | d    | 0 then | we have | d          | 1. Then |     |     |     |
|               |      | u   | p        |        | p    |        |         | u          |         |     |     |     |
|               |      |     | ≥        |        | ≥    |        |         |            | ≥       |     |     |     |
|               |      |     |          | d      | +d + | (d     | d )2    | 1          | 1,      |     |     |     |
|               |      |     |          | p      | u    |        | p u     |            |         |     |     |     |
|               |      |     |          |        |      |        | −       | − ≥        |         |     |     |     |
(cid:113)
| which | directly gives | a   | CFL | number | lower | than | 1.  |     |     |     |     |     |
| ----- | -------------- | --- | --- | ------ | ----- | ---- | --- | --- | --- | --- | --- | --- |

CHAPTER 3. STUDY OF THE STABILITY OF THE STAGGERED SCHEMES FOR
| THE      | ONE DIMENSIONAL |             |     | WAVE | SYSTEM |     |     |     |     |     |     |     | 51  |
| -------- | --------------- | ----------- | --- | ---- | ------ | --- | --- | --- | --- | --- | --- | --- | --- |
| Implicit | time            | integration |     |      |        |     |     |     |     |     |     |     |     |
The amplification matrix of the staggered implicit scheme (3.7) is given by
1
|     |     |     |      | c   | ∆t           |     |     | ∆t  |           |     |     | −   |     |
| --- | --- | --- | ---- | --- | ------------ | --- | --- | --- | --------- | --- | --- | --- | --- |
|     |     |     | 1+4d |     | 0 sin2(kπ∆x) |     |     | j2  | sin(kπ∆x) |     |     |     |     |
p
|     | AStag,Imp(k) |     |     |        | ∆x        |     |      | ρ   | ∆x         |     |     | ,   |        |
| --- | ------------ | --- | --- | ------ | --------- | --- | ---- | --- | ---------- | --- | --- | --- | ------ |
|     |              | =   |    |        |           |     |      | 0   |            |     |    |     | (3.19) |
|     |              |     |     | κ 0 ∆t |           |     |      |     | c 0 ∆t     |     |     |     |        |
|     |              |     |     | j2     | sin(kπ∆x) |     | 1+4d |     | sin2(kπ∆x) |     |     |     |        |
u
|     |     |     |      | ∆x   |      |            |     |      | ∆x  |           |            |     |     |
| --- | --- | --- | ----- | ---- | ---- | ---------- | --- | ---- | --- | --------- | ----------- | --- | --- |
|     |     |     |      |      | c    | ∆t         |     |      |     | ∆t        |            |     |     |
|     |     |     |       | 1+4d | 0    | sin2(kπ∆x) |     |      | j2  | sin(kπ∆x) |             |     |     |
|     |     |     | 1     |      | u    |            |     |      |     |           |             |     |     |
|     |     |     |       |      | ∆x   |            |     |      | − ρ | ∆x        |             |     |     |
|     |     | =   |       |      |      |            |     |      |     | 0         |             |     |     |
|     |     |     | µStag |     | κ ∆t |            |     |      |     | c ∆t      |             |    |     |
|     |     |     |       |      | j2 0 | sin(kπ∆x)  |     | 1+4d |     | 0         | sin2(kπ∆x). |     |     |
p
|     |     |     |     |  − | ∆x  |     |     |     |     | ∆x  |     |    |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
|     |     |     |     |    |     |     |     |     |     |     |     |    |     |
where
|       |          |     | c ∆t         |     |                  |      | c ∆t |            |     |          | c ∆t     | 2           |     |
| ----- | -------- | --- | ------------ | --- | ---------------- | ---- | ---- | ---------- | --- | -------- | -------- | ----------- | --- |
| µStag |          |     | 0 sin2(kπ∆x) |     |                  |      | 0    | sin2(kπ∆x) |     |          | 0        | sin2(kπ∆x). |     |
|       | = 1+4d   |     |              |     |                  | 1+4d |      |            |     | +4       |          |             |     |
|       |          | p   | ∆x           |     |                  |      | u ∆x |            |     |          | ∆x       |             |     |
|       | (cid:18) |     |              |     | (cid:19)(cid:18) |      |      |            |     | (cid:19) | (cid:18) | (cid:19)    |     |
Proposition 3.3.5 (staggered implicit schemes and von–Neumann necessary condition). The
staggered implicit schemes with d 0 and d 0 satisfy the von–Neumann necessary condi-
|      |            |     |         |     | p ≥ |     | u ≥ |     |     |     |     |     |     |
| ---- | ---------- | --- | ------- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| tion | (3.14) for | all | ∆t > 0. |     |     |     |     |     |     |     |     |     |     |
Proof. We are interested in the eigenvalues of the following matrix (3.19):
|       |              |     |        |         | c ∆t       |             |      |        | ∆t         |              |     | 1   |     |
| ----- | ------------ | --- | ------ | ------- | ---------- | ----------- | ---- | ------ | ---------- | ------------ | --- | --- | --- |
|       |              |     |        |         | 0          | sin2(kπ∆x)  |      |        |            |              |     | −   |     |
|       |              |     |        | 1+4d    | p          |             |      |        | j2         | sin(kπ∆x)    |     |     |     |
|       |              |     |        |         | ∆x         |             |      |        | ρ ∆x       |              |     |     |     |
|       | AStag,Imp(k) |     | =      |         |            |             |      |        | 0          |              |     |     | ,   |
|       |              |     |       |         | κ ∆t       |             |      |        | c          | ∆t           |     |    |     |
|       |              |     |        |         | 0          |             |      |        |            | 0 sin2(kπ∆x) |     |     |     |
|       |              |     |        | j2      |            | sin(kπ∆x)   |      | 1+4d   | u          |              |     |     |     |
|       |              |     |        |         | ∆x         |             |      |        |            | ∆x           |     |     |     |
|       |              |     |       |         |            |             |      |        |            |              |     |    |     |
|       |              |     |       |         |            |             |      |        |            |              |     |    |     |
| which | correspond   |     | to the | inverse | of the     | eigenvalues |      | of the | matrix     |              |     |     |     |
|       |              |     |        | c 0 ∆t  |            |             |      | ∆t     |            |              |     |     |     |
|       |              |     | 1+4d   |         | sin2(kπ∆x) |             |      | j2     | sin(kπ∆x)  |              |     |     |     |
|       |              |     |        | p ∆x    |            |             |      | ρ ∆x   |            |              |     |     |     |
|       |              |     |        |         |            |             |      | 0      |            |              | ,   |     |     |
|       |              |     |       | κ ∆t    |            |             |      | c      | ∆t         |              |    |     |     |
|       |              |     |        | 0       |            |             |      | 0      |            |              |     |     |     |
|       |              |     | j2     |         | sin(kπ∆x)  |             | 1+4d |        | sin2(kπ∆x) |              |     |     |     |
|       |              |     |        | ∆x      |            |             |      | u ∆x   |            |              |     |     |     |
|       |              |     |       |         |            |             |      |        |            |              |    |     |     |
|       |              |     |       |         |            |             |      |        |            |              |    |     |     |
that is exactly the matrix of the explicit case (3.15) where we have replaced d p by d p , d u by
−
| d   | and k by | k.  | The trace | of  | this matrix |     | is  |     |     |     |     |     |     |
| --- | -------- | --- | --------- | --- | ----------- | --- | --- | --- | --- | --- | --- | --- | --- |
u
| −   |     | −   |     |     |     |     |     |     |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
c ∆t
|     |     |     |     |     | 1+2(d | +d  | 0   | sin2(kπ∆x) |     | .   |     |     |     |
| --- | --- | --- | --- | --- | ----- | --- | --- | ---------- | --- | --- | --- | --- | --- |
|     |     |     |     | 2   |       | p   | u ) |            |     |     |     |     |     |
∆x
|       |                 |     |               | (cid:18) |     |       |         |              |     | (cid:19) |         |             |     |
| ----- | --------------- | --- | ------------- | -------- | --- | ----- | ------- | ------------ | --- | -------- | ------- | ----------- | --- |
| and   | its determinant |     | is            |          |     |       |         |              |     |          |         |             |     |
|       |                 |     | c ∆t          |          |     | 16d d | c2(∆t)2 |              |     |          | c2(∆t)2 |             |     |
|       |                 |     | 0 sin2(kπ∆x)+ |          |     | p     | u 0     | sin4(kπ∆x)+4 |     |          | 0       | sin2(kπ∆x), |     |
| 1+4(d | p +d            | u ) |               |          |     |       |         |              |     |          |         |             |     |
|       |                 |     | ∆x            |          |     | (∆x)2 |         |              |     |          | (∆x)2   |             |     |

CHAPTER 3. STUDY OF THE STABILITY OF THE STAGGERED SCHEMES FOR
| THE ONE | DIMENSIONAL |     |     | WAVE |     | SYSTEM |     |     |     |     |     |     | 52  |
| ------- | ----------- | --- | --- | ---- | --- | ------ | --- | --- | --- | --- | --- | --- | --- |
so that the reduced discriminant of the characteristic polynomial is given by
c2(∆t)2
|     |     |     |          | 0     | sin2(kπ∆x) |     |          | )2sin2(kπ∆x) |     |     |          |     |     |
| --- | --- | --- | -------- | ----- | ---------- | --- | -------- | ------------ | --- | --- | -------- | --- | --- |
|     |     |     | ∆ =4     |       |            |     | (d       | d            |     | 1   | .        |     |     |
|     |     |     | (cid:48) | (∆x)2 |            |     | p        | − u          |     | −   |          |     |     |
|     |     |     |          |       |            |     | (cid:16) |              |     |     | (cid:17) |     |     |
Case d d 1 : If d d 1, then the reduced discriminant is always negative, and
|                 | p   | u         |        | p   | u    |           |            |     |      |               |     |     |     |
| --------------- | --- | --------- | ------ | --- | ---- | --------- | ---------- | --- | ---- | ------------- | --- | --- | --- |
| |               | −   | | ≤       |        | | − | | ≤  |           |            |     |      |               |     |     |     |
| the eigenvalues |     | satisfy   |        |     |      |           |            |     |      |               |     |     |     |
|                 |     | 1         |        |     |      | c ∆t      |            |     |      |               |     |     |     |
|                 |     |           |        |     |      | 0         | sin2(kπ∆x) |     |      |               |     |     |     |
|                 |     |           | =1+2(d |     | p +d | u )       |            |     |      |               |     |     |     |
|                 |     | λStag,Imp |        |     |      | ∆x        |            |     |      |               |     |     |     |
|                 |     | ±         |        |     | c    | ∆t        |            |     |      |               |     |     |     |
|                 |     |           |        |     | 2j 0 | sin(kπ∆x) |            | 1   | (d d | )2sin2(kπ∆x). |     |     |     |
|                 |     |           |        |     |      |           |            |     | p u  |               |     |     |     |
|                 |     |           |        | ±   | ∆x   | |         |            | | − | −    |               |     |     |     |
(cid:113)
| The module | of  | these | eigenvalues |     | satisfies |     |     |     |     |     |     |     |     |
| ---------- | --- | ----- | ----------- | --- | --------- | --- | --- | --- | --- | --- | --- | --- | --- |
2
|     |           | 1   |     |          |     | c    | ∆t           |     |          |     |     |     |     |
| --- | --------- | --- | --- | -------- | --- | ---- | ------------ | --- | -------- | --- | --- | --- | --- |
|     |           |     | =   | 1+2(d    |     | +d ) | 0 sin2(kπ∆x) |     |          |     |     |     |     |
|     |           |     |     |          | p   | u    |              |     |          |     |     |     |     |
|     | λStag,Imp |     | 2   |          |     |      | ∆x           |     |          |     |     |     |     |
|     |           |     |     | (cid:18) |     |      |              |     | (cid:19) |     |     |     |     |
±
|     | (cid:12) |     | (cid:12) |     | c2(∆t)2 |     |     |     |     |     |     |     |     |
| --- | -------- | --- | -------- | --- | ------- | --- | --- | --- | --- | --- | --- | --- | --- |
)2sin2(kπ∆x)
|     | (cid:12) |     | (cid:12) | +4  | 0     | sin2(kπ∆x) |     | 1        | (d d |     |     |          |     |
| --- | -------- | --- | -------- | --- | ----- | ---------- | --- | -------- | ---- | --- | --- | -------- | --- |
|     | (cid:12) |     | (cid:12) |     | (∆x)2 |            |     |          | p    | u   |     |          |     |
|     |          |     |          |     |       |            |     | −        | −    |     |     |          |     |
|     |          |     |          |     |       |            |     | (cid:16) |      |     |     | (cid:17) |     |
1,
≥
2
|         | λStag,Exp |     |       |     | ∆t > |     |     |     |     |     |     |     |     |
| ------- | --------- | --- | ----- | --- | ---- | --- | --- | --- | --- | --- | --- | --- | --- |
| so that |           |     | 1 for | all |      | 0.  |     |     |     |     |     |     |     |
≤
±
|     | (cid:12) | (cid:12) |     |     |     |     |     |     |     |     |     |     |     |
| --- | -------- | -------- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
|     | (cid:12) | (cid:12) |     |     |     |     |     |     |     |     |     |     |     |
Case d(cid:12) d (cid:12)1 : In this case, the reduced discriminant is negative for small k, whereas
| |   | p − | u | ≥ |     |     |     |     |     |     |     |     |     |     |     |
| --- | --- | ----- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
sin2(kπ∆x)
| it is positive |     | for the | k such | that |     |     | is close | to 1. |     |     |     |     |     |
| -------------- | --- | ------- | ------ | ---- | --- | --- | -------- | ----- | --- | --- | --- | --- | --- |
k,
When the discriminant is negative For small the study led in the previous section can be
2
|          |     |        |        |             |     |           | λStag,Imp |     |       | ∆t  | > 0. |     |     |
| -------- | --- | ------ | ------ | ----------- | --- | --------- | --------- | --- | ----- | --- | ---- | --- | --- |
| done and | the | module | of the | eigenvalues |     | satisfies |           |     | 1 for | all |      |     |     |
≤
±
When the discriminant is positive For k such(cid:12)that sin2((cid:12)kπ∆x) is greater than 1/(d d )2,
|     |     |     |     |     |     |     | (cid:12) |     | (cid:12) |     |     | p   | u   |
| --- | --- | --- | --- | --- | --- | --- | -------- | --- | -------- | --- | --- | --- | --- |
−
the reduced discriminant is positive, and the eige(cid:12)nvalues ar(cid:12)e real. The eigenvalues satisfy
|     |     | 1         |        |     |      | c         | ∆t         |     |                |     |     |     |     |
| --- | --- | --------- | ------ | --- | ---- | --------- | ---------- | --- | -------------- | --- | --- | --- | --- |
|     |     |           |        |     |      | 0         | sin2(kπ∆x) |     |                |     |     |     |     |
|     |     |           | =1+2(d |     | p +d | u )       |            |     |                |     |     |     |     |
|     |     | λStag,Imp |        |     |      | ∆x        |            |     |                |     |     |     |     |
|     |     | ±         |        |     | c ∆t |           |            |     |                |     |     |     |     |
|     |     |           |        |     | 2 0  | sin(kπ∆x) |            | (d  | d )2sin2(kπ∆x) |     |     | 1   |     |
p u
|     |     |     |     |     | ± ∆x | |   |     | |   | −   |     | −   |     |     |
| --- | --- | --- | --- | --- | ---- | --- | --- | --- | --- | --- | --- | --- | --- |
(cid:113)
so that
|     |     |     |     |     |     | 1   |     | 1   |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
.
|     |     |     |     |     | λStag,Imp |     | ≥ λStag,Imp |     |     |     |     |     |     |
| --- | --- | --- | --- | --- | --------- | --- | ----------- | --- | --- | --- | --- | --- | --- |
+
−
We then investigate conditions for ensuring 1/λStag,Imp 1, which reads
|     |     |     |            |     |     |      | −         | ≥   |      |              |     |      |     |
| --- | --- | --- | ---------- | --- | --- | ---- | --------- | --- | ---- | ------------ | --- | ---- | --- |
|     |     | c   | ∆t         |     |     | c ∆t |           |     |      |              |     |      |     |
| 2(d | +d  | ) 0 | sin2(kπ∆x) |     |     | 2 0  | sin(kπ∆x) |     | (d d | )2sin2(kπ∆x) |     | 1 0, |     |
|     | p   | u   |            |     |     |      |           |     | p u  |              |     |      |     |
|     |     | ∆x  |            |     | −   | ∆x   | |         | |   | −    |              |     | − ≥  |     |
(cid:113)

CHAPTER 3. STUDY OF THE STABILITY OF THE STAGGERED SCHEMES FOR
THE ONE DIMENSIONAL WAVE SYSTEM 53
which can be simplified as
sin(kπ∆x) (d d )2sin2(kπ∆x) 1 (d +d )sin2(kπ∆x).
p u p u
| | − − ≤
(cid:113)
Taking the square and simplifying by sin2(kπ∆x) gives
(d d )2sin2(kπ∆x) 1 (d +d )2sin2(kπ∆x),
p u p u
− − ≤
which may be rewritten
1+4d d sin2(kπ∆x) 0,
p u
≥
which is always ensured. This means that we always have
λStag,Imp
1.
± ≤
Conclusion No matter the case, λStag,Imp is always bounded by 1.
±
Implicit-explicit time integration
The amplification matrix of the staggered implicit-explicit pressure-centered scheme (3.8) is
given by
AStag,PC,Imp-Exp(k) =
1 c ∆t ∆t
1 0 − 1 4d 0 sin2(kπ∆x) j2 sin(kπ∆x)
p
κ ∆t − ∆x − ρ ∆x
 j2 0 sin(kπ∆x) 1   0 
∆x 0 1
   
c ∆t ∆t
1 4d 0 sin2(kπ∆x) j2 sin(kπ∆x)
p
− ∆x − ρ ∆x
=  κ ∆t c ∆t c 0 ∆t 2 .
j2 0 sin(kπ∆x) 1 4d 0 sin2(kπ∆x) 1 4 0 sin2(kπ∆x)
p
 − ∆x − ∆x − ∆x 
 (cid:18) (cid:19) (cid:18) (cid:19) 
 
(3.20)
Proposition 3.3.6 (staggered implicit-explicit schemes and von–Neumann necessary condi-
tion). The staggered implicit-explicit fully-centered (d = d = 0) and pressure-centered (d > 0
p u p
and d = 0) schemes satisfy the von–Neumann necessary condition (3.14) under the CFL con-
u
dition
c ∆t
0 0 d + 1+d2 (3.21)
≤ ∆x ≤ − p p
(cid:113)
so that the CFL condition is optimal for d = 0.
p
Proof. We are interested in the eigenvalues of the following matrix (3.20):
AStag,PC,Imp-Exp(k) =

CHAPTER 3. STUDY OF THE STABILITY OF THE STAGGERED SCHEMES FOR
THE ONE DIMENSIONAL WAVE SYSTEM 54
c ∆t ∆t
1 4d 0 sin2(kπ∆x) j2 sin(kπ∆x)
p
− ∆x − ρ ∆x
 κ ∆t c ∆t c 0 ∆t 2 .
j2 0 sin(kπ∆x) 1 4d 0 sin2(kπ∆x) 1 4 0 sin2(kπ∆x)
p
 − ∆x − ∆x − ∆x 
 (cid:18) (cid:19) (cid:18) (cid:19) 
 
The trace of the matrix is
c ∆t c ∆t
Tr AStag,PC,Imp-Exp(k) = 2 4 0 sin2(kπ∆x) d + 0 ,
p
− ∆x ∆x
(cid:18) (cid:19)
(cid:0) (cid:1)
and its determinant is
c ∆t
det AStag,PC,Imp-Exp(k) = 1 4d 0 sin2(kπ∆x).
p
− ∆x
(cid:0) (cid:1)
The reduced discriminant of the characteristic polynomial is
c ∆t 2 c ∆t 2
∆ = 4 0 sin2(kπ∆x) 1 sin2(kπ∆x) d + 0 .
(cid:48) p
− ∆x − ∆x
(cid:32) (cid:33)
(cid:18) (cid:19) (cid:18) (cid:19)
c ∆t 2
When the discriminant is negative : If 1 sin2(kπ∆x) d + 0 0, the discrim-
p
− ∆x ≥
(cid:18) (cid:19)
inant is negative and the eigenvalues of the amplification matrix (3.20) are given by
c ∆t c ∆t
λStag,PC,Imp-Exp =1 2 0 sin2(kπ∆x) d + 0
p
± − ∆x ∆x
(cid:18) (cid:19)
c ∆t c ∆t 2
j2 0 sin(kπ∆x) 1 sin2(kπ∆x) d + 0 ,
p
± ∆x | |(cid:115) − ∆x
(cid:18) (cid:19)
so that
2 c ∆t
λStag,PC,Imp-Exp = 1 4d 0 sin2(kπ∆x),
p
± − ∆x
(cid:12) (cid:12)
and
λStag,PC,Imp-Exp (cid:12)
(cid:12) 1 is always satis
(cid:12)
(cid:12)fied.
± ≤
(cid:12) (cid:12)
(cid:12) (cid:12)
(cid:12) (cid:12) c ∆t 2
When the discriminant is positive : If 1 sin2(kπ∆x) d + 0 0, the discrim-
p
− ∆x ≤
(cid:18) (cid:19)
inant is positive and the eigenvalues are real and given by
c ∆t c ∆t
λStag,PC,Imp-Exp =1 2 0 sin2(kπ∆x) d + 0
p
± − ∆x ∆x
(cid:18) (cid:19)
c ∆t c ∆t 2
2 0 sin(kπ∆x) 1+sin2(kπ∆x) d + 0
p
± ∆x | |(cid:115)− ∆x
(cid:18) (cid:19)
so that
λStag,PC,Imp-Exp λStag,PC,Imp-Exp.
− ≤ +

CHAPTER 3. STUDY OF THE STABILITY OF THE STAGGERED SCHEMES FOR
THE ONE DIMENSIONAL WAVE SYSTEM 55
Investigation on
λStag,PC,Imp-Exp
1. This condition reads
+ ≤
c ∆t c ∆t c ∆t c ∆t 2
2 0 sin2(kπ∆x) d + 0 +2 0 sin(kπ∆x) 1+sin2(kπ∆x) d + 0 0,
p p
− ∆x ∆x ∆x | |(cid:115)− ∆x ≤
(cid:18) (cid:19) (cid:18) (cid:19)
which can be simplified as
c ∆t 2 c ∆t
sin(kπ∆x) 1+sin2(kπ∆x) d + 0 sin2(kπ∆x) d + 0 .
p p
| |(cid:115)− ∆x ≤ ∆x
(cid:18) (cid:19) (cid:18) (cid:19)
Taking the squared and simplifying by sin2(kπ∆x), we obtain
c ∆t 2 c ∆t 2
1+sin2(kπ∆x) d + 0 sin2(kπ∆x) d + 0 ,
p p
− ∆x ≤ ∆x
(cid:18) (cid:19) (cid:18) (cid:19)
which is always ensured. Then, we always have λStag,PC,Imp-Exp 1.
+ ≤
Investigation on
λStag,PC,Imp-Exp
1. We have
λStag,PC,Imp-Exp
1 that reads
− ≥ − − ≥ −
c ∆t c ∆t c ∆t c ∆t 2
0 sin2(kπ∆x) d + 0 + 0 sin(kπ∆x) 1+sin2(kπ∆x) d + 0 1.
p p
∆x ∆x ∆x | |(cid:115)− ∆x ≤
(cid:18) (cid:19) (cid:18) (cid:19)
Since the maximal value is reached for sin2(kπ∆x) = 1, we investigate conditions such that
c ∆t c ∆t c ∆t c ∆t 2
0 d + 0 + 0 1+ d + 0 1
p p
∆x ∆x ∆x (cid:115)− ∆x ≤
(cid:18) (cid:19) (cid:18) (cid:19)
c ∆t c ∆t c ∆t 2
0 d + 0 + 1+ d + 0 1.
p p
⇔ ∆x  ∆x (cid:115)− ∆x  ≤
(cid:18) (cid:19)
 

CHAPTER 3. STUDY OF THE STABILITY OF THE STAGGERED SCHEMES FOR
THE ONE DIMENSIONAL WAVE SYSTEM 56
c ∆t c ∆t 2
Multiplying by d + 0 1+ d + 0 > 0, we obtain
p p
∆x −(cid:115)− ∆x
(cid:18) (cid:19)
c ∆t c ∆t 2 c ∆t 2 c ∆t c ∆t 2
0 d + 0 1+ d + 0 d + 0 1+ d + 0
p p p p
∆x
(cid:32)
∆x −
(cid:32)
− ∆x
(cid:33)(cid:33)
≤ ∆x −(cid:115)− ∆x
(cid:18) (cid:19) (cid:18) (cid:19) (cid:18) (cid:19)
c ∆t c ∆t c ∆t 2
0 d + 0 1+ d + 0
p p
⇔ ∆x ≤ ∆x −(cid:115)− ∆x
(cid:18) (cid:19)
c ∆t 2
0 d 1+ d + 0
p p
⇔ ≤ −(cid:115)− ∆x
(cid:18) (cid:19)
c ∆t 2
1+ d + 0 d2
⇔− p ∆x ≤ p
(cid:18) (cid:19)
c ∆t 2 c ∆t
0 +2d 0 1 0
p
⇔ ∆x ∆x − ≤
(cid:18) (cid:19)
c ∆t
d 1+d2 0 d + 1+d2,
⇔− p − p ≤ ∆x ≤ − p p
(cid:113) (cid:113)
which gives the stability condition (3.21).
For the fully-centered scheme, injecting d = 0 on the expression of the eigenvalues ob-
p
tained for the pressure-centered scheme, under CFL condition (3.21) with d = 0, we have
p
λStag,FC,Imp-Exp
= 1 that concludes the proof.
±
(cid:12) (cid:12)
(cid:12) (cid:12)
(cid:12) (cid:12)
3.3.2 Energy dissipation
The von-Neumann condition only gives a necessary condition for stability. Indeed, some
schemes lose the normal structure of the amplification matrix that enables to conclude on
the dissipation of their energy. As detailed in Proposition 3.3.1, we should thereby study the
matrix I
2
(S 1 2A(k)S
−
1 2)
∗
(S 1 2A(k)S
−
1 2) to draw general conclusions.
−
Explicit time integration
Proposition 3.3.7 (energy dissipation of the staggered explicit fully-upwind scheme). The
staggered explicit fully-upwind (d = d > 0) scheme (3.6) dissipates energy under the CFL
p u
condition (3.16).
1 1
Proof. For the explicit fully-upwind scheme the matrix S2A(k)S −2 is given by
c ∆t 2c ∆t
1 4d 0 sin2(kπ∆x) j 0 sin(kπ∆x)
S 1 2AStag,Up,Exp(k)S − 1 2 =  − 2 p c ∆ ∆ x t − ∆ c x ∆t  ,
j 0 sin(kπ∆x) 1 4d 0 sin2(kπ∆x)
p
− ∆x − ∆x
 
 

CHAPTER 3. STUDY OF THE STABILITY OF THE STAGGERED SCHEMES FOR
THE ONE DIMENSIONAL WAVE SYSTEM 57
and then, it is normal because
S 1 2AStag,Up,Exp(k)S − 1 2 ∗ S 1 2AStag,Up,Exp(k)S − 1 2
(cid:16) c ∆t (cid:17) (cid:16) 2 2c ∆t 2 (cid:17) 1 0
= 1 4d 0 sin2(kπ∆x) + 0 sin2(kπ∆x)
p
− ∆x ∆x 0 1
(cid:32) (cid:33)
(cid:18) (cid:19) (cid:18) (cid:19) (cid:18) (cid:19)
= S2 1 AStag,Up,Exp(k)S − 1 2 S 1 2AStag,Up,Exp(k)S − 1 2 ∗.
(cid:16) (cid:17)(cid:16) (cid:17)
The result is obtained using and Proposition 3.3.4, Proposition 3.3.3.
Concerning the explicit pressure-centered scheme, the matrix S2 1 AStag,PC,Exp(k)S − 1 2 is
given by
c ∆t 2c ∆t
1 4d 0 sin2(kπ∆x) j 0 sin(kπ∆x)
S2 1 AStag,PC,Exp(k)S − 1 2 =  − 2 p c ∆ ∆ x t − ∆x  ,
j 0 sin(kπ∆x) 1
− ∆x
 
 
and so is not normal because
S 1 2AStag,PC,Exp(k)S − 1 2 ∗ S 1 2AStag,PC,Exp(k)S − 1 2
(cid:16) c ∆t (cid:17) (cid:16) 2 2c ∆t 2 (cid:17) 2c ∆t 2
1 4d 0 sin2(kπ∆x) + 0 sin2(kπ∆x) j2d 0 sin3(kπ∆x)
p p
− ∆x ∆x ∆x
= (cid:18) (cid:19) (cid:18) (cid:19) (cid:18) (cid:19) 
2c ∆t 2 2c ∆t 2
j2d 0 sin3(kπ∆x) 1+ 0 sin2(kπ∆x)
 − p ∆x ∆x 
 (cid:18) (cid:19) (cid:18) (cid:19) 
 
= S 1 2AStag,PC,Exp(k)S − 1 2 S 1 2AStag,PC,Exp(k)S −2 1 ∗.
(cid:16) (cid:17)(cid:16) (cid:17)
We can prove that the staggered pressure-centered scheme does not dissipate energy. More
precisely, we have
Proposition 3.3.8 (energy dissipation of the staggered explicit pressure-centered scheme).
The staggered explicit pressure-centered scheme (3.6) does not dissipate energy. Indeed, for all
∆t > 0 satisfying (3.17), we have
Un R2 N , E Un+1 > E(Un).
∃ h ∈ h h
More precisely, for all frequency k (cid:0) = 0 (cid:1) , we can (cid:0) build a (cid:1) vector Un(k) of frequency k such that
(cid:54) h
E Un+1(k) > E(Un(k)).
h h
Pr(cid:0)oof. From(cid:1)the proof of Proposition 3.3.1, we have
E Un h +1 − E(Un h ) = S 1 2Un k ∗ S2 1 A(k)S2 1 ∗ S 1 2A(k)S 1 2 − I 2 S2 1 Un k .
k
(cid:16) (cid:17) (cid:16)(cid:16) (cid:17) (cid:16) (cid:17) (cid:17)(cid:16) (cid:17)
(cid:0) (cid:1) (cid:80)
(cid:98) (cid:98)

CHAPTER 3. STUDY OF THE STABILITY OF THE STAGGERED SCHEMES FOR
THE ONE DIMENSIONAL WAVE SYSTEM 58
For the explicit pressure-centered scheme, we have
I 2 S 1 2AStag,PC,Exp(k)S 1 2 ∗ S 1 2AStag,PC,Exp(k)S 1 2
−
(cid:16) (cid:17)c ∆(cid:16)t (cid:17) c ∆t
=4 c 0 ∆t sin2(kπ∆x) 2d p − ∆ 0 x 1+4d2 p sin2(kπ∆x) − j2d p ∆ 0 x sin(kπ∆x) ,
∆x  j2d c 0(cid:0) ∆t sin(kπ∆x) (cid:1) c 0 ∆t 
p
∆x − ∆x
 
 
M
with (cid:124) (cid:123)(cid:122) (cid:125)
c ∆t
Tr(M) = 2 d 0 1+2d2sin2(kπ∆x) ,
p − ∆x p
(cid:18) (cid:19)
(cid:0) (cid:1)
and
c ∆t c ∆t
det(M) = 0 2d 0 ,
p
− ∆x − ∆x
(cid:18) (cid:19)
so that the reduced discriminant ∆ is given by
(cid:48)
c ∆t c ∆t
∆ = d2 1 4 0 sin2(kπ∆x) d 0 1+d2sin2(kπ∆x) .
(cid:48) p − ∆x p − ∆x p
(cid:18) (cid:18) (cid:19)(cid:19)
(cid:0) (cid:1)
Then, the matrix I
2
S2 1 AStag,PC,Exp(k)S2 1 ∗ S 1 2AStag,PC,Exp(k)S2 1 has two different ei-
−
genvalues λ and λ+ g(cid:16)iven by (cid:17) (cid:16) (cid:17)
−k k
4c ∆t c ∆t
λ = 0 sin2(kπ∆x) d 0 1+2d2sin2(kπ∆x)
±k ∆x p − ∆x p
(cid:18)
(cid:0) (cid:1)
c ∆t c ∆t
d 1 4 0 sin2(kπ∆x) d 0 1+d2sin2(kπ∆x) .
± p (cid:115) − ∆x p − ∆x p (cid:33)
(cid:18) (cid:19)
(cid:0) (cid:1)
c ∆t
Since ∆t satisfies (3.17), we have 0 d < 2d so that det(M) < 0 and λ < 0 < λ+ if
∆x ≤ p p −k k
sin(kπ∆x)
(cid:54)
= 0 which means that k
(cid:54)
= 0. Taking for S 1 2Un
k
an eigenvector associated to the
eigenvalue λ given by
−k
S2 1 Un
k
=
jδ
1
(k)
, (cid:98)
(cid:18) − (cid:19)
where (cid:98)
∆x c ∆t
δ(k) = 1 2d 0 sin2(kπ∆x)
p
2c ∆tsin(kπ∆x) − ∆x
0 (cid:18)
c ∆t c ∆t
+ 1 4 0 sin2(kπ∆x) d 0 1+d2sin2(kπ∆x) ,
(cid:115) − ∆x p − ∆x p (cid:33)
(cid:18) (cid:19)
(cid:0) (cid:1)

CHAPTER 3. STUDY OF THE STABILITY OF THE STAGGERED SCHEMES FOR
THE ONE DIMENSIONAL WAVE SYSTEM 59
we have
S 1 2Un k ∗ S2 1 AStag,PC,Exp(k)S − 1 2 ∗S 1 2AStag,PC,Exp(k)S − 1 2 I 2 S2 1 Un k
−
= (cid:16)
−
λ −k (cid:98) S (cid:17) 2 1 U (cid:16)(cid:16) n k 2
2
(cid:17) (cid:17)(cid:16) (cid:98) (cid:17)
(cid:13) (cid:13)
>0. (cid:13) (cid:13)
(cid:98)
(cid:13) (cid:13)
Then, the pressure pn(k) defined on all cell C by
h i
1 1
pn(k) = √2 pn e 2jkπxi + pne2jkπxi
i √2
−
k − √2 k
(cid:18) (cid:19)
= 2√2cos(2πkx ), (3.22)
(cid:98) i (cid:98)
and the velocity un(k) defined on all cell C by
h i+1/2
1 1
un i+1/2 (k) = √2ρ 0 c 0 (cid:18) √2ρ 0 c 0 un − k e − 2jkπx i+1/2 + √2ρ 0 c 0 un k e2jkπx i+1/2 (cid:19)
= 2√2ρ c δ(k)sin(2πkx ), (3.23)
0 0 (cid:98) i+1/2 (cid:98)
−
gives E Un+1(k) > E(Un(k)).
h h
(cid:0) (cid:1)
Implicit time integration
Proposition 3.3.9 (staggeredimplicitschemesandenergydissipation). The staggered implicit
fully-upwind (d = d > 0), pressure-centered (d > 0 and d = 0), velocity-centered (d = 0
p u p u p
and d > 0) and fully-centered (d = d = 0) schemes dissipate energy for all ∆t > 0.
u p u
1 1
Proof. For the implicit staggered schemes (3.7), the matrix S2A(k)S −2 is given by
c ∆t 2c ∆t
S 1 2AStag,Imp(k)S − 1 2 = µS 1 tag  1+4 j d 2 u c 0 ∆ ∆ 0 x t si s n in ( 2 kπ (k ∆ π x ∆ ) x) 1+ − j 4d ∆ 0 c x 0 ∆t si s n in ( 2 k ( π k ∆ π x ∆ ) x)  ,
p
− ∆x ∆x
 
 
so that
S2 1 AStag,Imp(k)S − 1 2 ∗ S 1 2AStag,Imp(k)S − 1 2
(cid:16) (cid:17) (cid:16) (cid:17) c ∆t 3
α j8(d d ) 0 sin3(kπ∆x)
1 p − u ∆x
=  (cid:18) (cid:19) ,
|
µStag
|
2
 − j8(d p − d u )
c
∆ 0
∆
x
t 3
sin3(kπ∆x) β 
 (cid:18) (cid:19) 
 
= S 1 2AStag,Imp(k)S − 1 2 S2 1 AStag,Imp(k)S − 1 2 ∗,
(cid:16) (cid:17)(cid:16) (cid:17)

CHAPTER 3. STUDY OF THE STABILITY OF THE STAGGERED SCHEMES FOR
THE ONE DIMENSIONAL WAVE SYSTEM 60
where
c ∆t 2 c ∆t 3
α = 1+4d 0 sin2(kπ∆x) +4 0 sin3(kπ∆x),
u
∆x ∆x
(cid:18) (cid:19) (cid:18) (cid:19)
c ∆t 2 c ∆t 3
β = 1+4d 0 sin2(kπ∆x) +4 0 sin3(kπ∆x).
p
∆x ∆x
(cid:18) (cid:19) (cid:18) (cid:19)
Then, if d p = d u , the matrix S 1 2AStag,Imp(k)S −2 1 is normal and the result concerning the
staggeredimplicitfully-upwindandfully-centeredschemesareobtainedusingProposition3.3.3
and Proposition 3.3.5.
For the pressure-centered scheme (d p > 0 and d u = 0), since the matrix S2 1 AStag,Imp(k)S − 1 2
is no more normal, Proposition 3.3.1 is used to prove the energy dissipation. Since
c ∆t c ∆t 2
µStag,PC = 1+4d 0 sin2(kπ∆x)+4 0 sin2(kπ∆x),
p
∆x ∆x
(cid:18) (cid:19)
the matrix I 2 S2 1 A(k)S −2 1 ∗ S 1 2A(k)S − 1 2 is given by
−
(cid:16) (cid:17) (cid:16) (cid:17)
I 2 S 1 2AStag,PC,Imp(k)S −2 1 ∗ S 1 2AStag,PC,Imp(k)S − 1 2
−
(cid:16) (cid:17) (cid:16) (cid:17) c ∆t 2
γ j2d 0 sin(kπ∆x)
= 4c 0 ∆tsin2(kπ∆x)  − p (cid:18) ∆x (cid:19) ,
∆x(µStag,PC) 2
j2d
c
0
∆t 2
sin(kπ∆x) δ
 p ∆x 
 (cid:18) (cid:19) 
 
M
where (cid:124) (cid:123)(cid:122) (cid:125)
c ∆t c ∆t 2 c ∆t 3
γ = 2d + 0 1+4d2sin2(kπ∆x) +8d 0 sin2(kπ∆x)+4 0 sin2(kπ∆x),
p ∆x p p ∆x ∆x
(cid:18) (cid:19) (cid:18) (cid:19)
c ∆t (cid:0) c ∆t (cid:1) c ∆t 2
δ = 0 1+8d 0 sin2(kπ∆x)+4 0 sin2(kπ∆x) .
p
∆x ∆x ∆x
(cid:32) (cid:33)
(cid:18) (cid:19)
Then, for all ∆t > 0, we have
Tr(M) = γ +δ 0,
≥
c ∆t 4
det(M) = γδ 4d2 0 sin4(kπ∆x)
− p ∆x
(cid:18) (cid:19)
c ∆t 4 c ∆t 4
64d2 0 sin4(kπ∆x) 4d2 0 sin4(kπ∆x)
≥ p ∆x − p ∆x
(cid:18) (cid:19) (cid:18) (cid:19)
0,
≥
that concludes the proof.

CHAPTER 3. STUDY OF THE STABILITY OF THE STAGGERED SCHEMES FOR
THE ONE DIMENSIONAL WAVE SYSTEM 61
Implicit-explicit time integration
Proposition 3.3.10 (energy dissipation of the staggered implicit-explicit pressure-centered
scheme). The staggered implicit-explicit pressure-centered (d > 0,d = 0) scheme (3.8) dissip-
p u
ates energy under the CFL condition
c ∆t 2d
0 0 p . (3.24)
≤ ∆x ≤ 1+4d2
p
and so the CFL condition is optimal for d = 1/2.
p
Proof. For the implicit-explicit pressure-centered scheme, the matrix
S 1 2AStag,PC,Imp-Exp(k)S − 1 2 is given by
S2 1 AStag,PC,Imp-Exp(k)S − 1 2
c ∆t 2c ∆t
1 4d 0 sin2(kπ∆x) j 0 sin(kπ∆x)
p
− ∆x − ∆x
=  2c ∆t c ∆t c ∆t 2  .
j 0 sin(kπ∆x) 1 4d 0 sin2(kπ∆x) 1 4 0 sin2(kπ∆x)
p
− ∆x − ∆x − ∆x
 (cid:18) (cid:19) (cid:18) (cid:19) 
 
Then, we have
S2 1 AStag,PC,Imp-Exp(k)S − 1 2 ∗ S 1 2AStag,PC,Imp-Exp(k)S −2 1
(cid:16) c ∆t 2 (cid:17) (cid:16) c ∆t 3 (cid:17)
α2 1+4 0 sin2(kπ∆x) j8 0 sin3(kπ∆x)α
∆x − ∆x
= (cid:32) (cid:18) (cid:19) (cid:33) (cid:18) (cid:19) ,
c ∆t 3 c ∆t 2
 j8 0 sin3(kπ∆x)α β2+4 0 sin2(kπ∆x) 
 ∆x ∆x 
 (cid:18) (cid:19) (cid:18) (cid:19) 
 
and
S 1 2AStag,PC,Imp-Exp(k)S − 1 2 S 1 2AStag,PC,Imp-Exp(k)S − 1 2 ∗
(cid:16) c ∆t 2 (cid:17)(cid:16) c ∆t (cid:17)
α2+4 0 sin2(kπ∆x) j2 0 sin(kπ∆x) α2 β
∆x ∆x −
= (cid:18) (cid:19) ,
c ∆t c ∆t 2 (cid:0) (cid:1)
j2 0 sin(kπ∆x) α2 β β2+4α2 0 sin2(kπ∆x)
 − ∆x − ∆x 
 (cid:18) (cid:19) 
 (cid:0) (cid:1) 
where
c ∆t
α = 1 4d 0 sin2(kπ∆x),
p
− ∆x
c ∆t 2
β = 1 4 0 sin2(kπ∆x),
− ∆x
(cid:18) (cid:19)
so that thematrix S 1 2AStag,PC,Imp-Exp(k)S − 1 2 is not normal. From Proposition3.3.1, it issuffi-
cienttoprovethatthematrixI 2 S 1 2AStag,PC,Imp-Exp(k)S −2 1 ∗ S 1 2AStag,PC,Imp-Exp(k)S − 1 2
−
(cid:16) (cid:17) (cid:16) (cid:17)

CHAPTER 3. STUDY OF THE STABILITY OF THE STAGGERED SCHEMES FOR
THE ONE DIMENSIONAL WAVE SYSTEM 62
is positive to prove the energy dissipation. We get that
I 2 S 1 2AStag,PC,Imp-Exp(k)S − 1 2 ∗ S 1 2AStag,PC,Imp-Exp(k)S −2 1
−
(cid:16) (cid:17) (cid:16) c ∆ (cid:17) t 2
γ j2 0 sin(kπ∆x)α
c ∆t ∆x
=4 0 sin2(kπ∆x) (cid:18) (cid:19) ,
∆x c ∆t 2 c ∆t c ∆t 2
 j2 0 sin(kπ∆x)α 0 1 4 0 sin2(kπ∆x) 
 − ∆x ∆x (cid:32) − ∆x (cid:33) 
 (cid:18) (cid:19) (cid:18) (cid:19) 
 
=M
where (cid:124) (cid:123)(cid:122) (cid:125)
c ∆t c ∆t c ∆t 2
γ = 2d 0 1+4d2sin2(kπ∆x) 8d 0 sin2(kπ∆x)+16d2 0 sin4(kπ∆x) .
p − ∆x p − p ∆x p ∆x
(cid:32) (cid:33)
(cid:18) (cid:19)
Then, we have
c ∆t c ∆t
det(M) = 0 2d 0 1+4d2sin2(kπ∆x)
∆x p − ∆x p
(cid:18) (cid:19)
c ∆t c ∆t (cid:0) (cid:1)
0 2d 0 1+4d2 ,
≥ ∆x p − ∆x p
(cid:18) (cid:19)
(cid:0) (cid:1)
and det(M) 0 under CFL condition (3.24). Moreover, we get
≥
c ∆t c ∆t 2 c ∆t 2
Tr(M) = 2d 4 0 sin2(kπ∆x) d 0 +4d2 0 sin2(kπ∆x)
p − ∆x p − ∆x p ∆x
(cid:32) (cid:33)
(cid:18) (cid:19) (cid:18) (cid:19)
c ∆t c ∆t 2 c ∆t 2
2d 4 0 d 0 +4d2 0
≥ p − ∆x p − ∆x p ∆x
(cid:32) (cid:33)
(cid:18) (cid:19) (cid:18) (cid:19)
c ∆t c ∆t c ∆t
= 2d 4 0 d2 0 2d 0 1+4d2 ,
p − ∆x p − ∆x p − ∆x p
(cid:18) (cid:18) (cid:19)(cid:19)
(cid:0) (cid:1)
so that, under the CFL condition (3.24), we obtain
c ∆t c ∆t
Tr(M) 2d 4d2 0 = 2d 1 2d 0 .
≥ p − p ∆x p − p ∆x
(cid:18) (cid:19)
Then, if
c ∆t 2d 1 2d
0 0 min p , = p ,
≤ ∆x ≤ 1+4d2 2d 1+4d2
(cid:18) p p (cid:19) p
we have Tr(M) 0 and det(M) 0 so that the matrix
≥ ≥
I 2 S 1 2AStag,PC,Imp-Exp(k)S − 1 2 ∗ S 1 2AStag,PC,Imp-Exp(k)S − 1 2 ,
−
(cid:16) (cid:17) (cid:16) (cid:17)
2x
is positive. Concerning optimality, it is sufficient to remark that function x has its
(cid:55)→ 1+4x2

CHAPTER 3. STUDY OF THE STABILITY OF THE STAGGERED SCHEMES FOR
| THE ONE | DIMENSIONAL |      | WAVE |     | SYSTEM |     |     |     |     |     | 63  |
| ------- | ----------- | ---- | ---- | --- | ------ | --- | --- | --- | --- | --- | --- |
| maximum | in x =      | 1/2. |      |     |        |     |     |     |     |     |     |
Remark 3.3.1. Theimplicit-explicitpressure-centeredschemedissipatesenergyundertheCFL
condition (3.24) which is more restrictive than the CFL condition (3.21) to satisfy the von–
| Neumann | necessary | condition |     | because |     |     |     |     |     |     |     |
| ------- | --------- | --------- | --- | ------- | --- | --- | --- | --- | --- | --- | --- |
2d
|     |     |     |     |      |       | p   |     | 1+d2. |     |     |     |
| --- | --- | --- | --- | ---- | ----- | --- | --- | ----- | --- | --- | --- |
|     |     |     | d   | p 0, |       |     | d p | +     |     |     |     |
|     |     |     | ∀   | ≥    | 1+4d2 |     | ≤ − |       | p   |     |     |
p
(cid:113)
Proposition 3.3.11. The staggered implicit-explicit fully-centered scheme (3.8) does not dis-
sipate energy.
|     |     | I   | 1   |     |     |     | 1   | 1   |     | 1   |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
Proof. The matrix S2 AStag,FC,Imp-Exp(k)S 2 ∗ S 2AStag,FC,Imp-Exp(k)S 2 of the
|     |     | 2   |     |     |     |     | −   |     |     | −   |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
−
fully-centered scheme is o(cid:16)btained by replacing d by(cid:17)zer(cid:16)o in the expression of the m(cid:17)atrix of
p
the pressure-centered scheme in the proof of Proposition 3.3.10. Then, we get that det(M) =
2
c ∆t
0 < 0 and so the matrix M has a negative eigenvalue. A counter-example can then
− ∆x
| (cid:18) | (cid:19)  |       |                |     |        |     |     |     |     |     |     |
| -------- | --------- | ----- | -------------- | --- | ------ | --- | --- | --- | --- | --- | --- |
| be build | as in the | proof | of Proposition |     | 3.3.8. |     |     |     |     |     |     |
Remark 3.3.2. As noted in section 3.2, some ImEx schemes can be seen as fully explicit in
the sense that it is sufficient to plug the pressure equation in the pressure gradient to obtain an
explicit formulation (in other words, no matrix inversion is needed). In this case, the previous
proposition Proposition 3.3.11 could be intuited from Proposition 3.3.8.
| 3.4 | l -stability |     | on  | the | characteristic |     |     | variables |     |     |     |
| --- | ------------ | --- | --- | --- | -------------- | --- | --- | --------- | --- | --- | --- |
∞
Since for staggered schemes, p h is piecewise constant per primal cell C i = x ;x
i 1/2 i+1/2
−
andu h ispiecewiseconstantperdualcellC i+1/2 = ]x i ;x i+1 [,thedefinitionofthecharacteristic
|     |     |     |     |     |     |     |     |     |     | (cid:3) | (cid:2) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | ------- | ------- |
variables is not obvious. In this section, we study l -stability of staggered schemes based on
∞
two different definitions of the characteristic variables. Each definition of the characteristic
| variables | will lead      | to different |           | properties | on      | the | staggered | schemes. |     |     |     |
| --------- | -------------- | ------------ | --------- | ---------- | ------- | --- | --------- | -------- | --- | --- | --- |
| 3.4.1     | Characteristic |              | variables |            | defined |     | on half   | cells    |     |     |     |
The characteristic variables are defined such that they are piecewise constant on half cells.
Definition 3.4.1 (Characteristic variables on half cells). The characteristic variables are
Ch±
| defined | on primal | cells | C = | x     | ;x    | as  |     |     |     |     |     |
| ------- | --------- | ----- | --- | ----- | ----- | --- | --- | --- | --- | --- | --- |
|         |           |       | i   | i 1/2 | i+1/2 |     |     |     |     |     |     |
−
|     |     |     | (cid:3)  |       |     | (cid:2) |      |         |          |     |     |
| --- | --- | --- | -------- | ----- | --- | ------- | ---- | ------- | -------- | --- | --- |
|     |     |     |          |       | p u | /(ρ     | c ), | on      | x ;x ,   |     |     |
|     |     |     |          | ± =   | i i | 1 / 2   | 0 0  |         | i 1/ 2 i |     |     |
|     |     |     | = C      | i , L | ±   |         |      |         |          |     |     |
|     |     | Ch± |          | =     | p u | − /(ρ   | c ), | on      | x − ;x . |     |     |
|     |     |     | (cid:40) | ±     | i i | + 1 / 2 | 0 0  |         | i i+ 1/2 |     |     |
|     |     |     | C        | i , R | ±   |         |      | (cid:3) | (cid:2)  |     |     |
|     |     |     |          |       |     |         |      | (cid:3) | (cid:2)  |     |     |

CHAPTER 3. STUDY OF THE STABILITY OF THE STAGGERED SCHEMES FOR
THE ONE DIMENSIONAL WAVE SYSTEM 64
Using the finite volume interpretation of subsection 3.2.2, we can prove this partial l -
∞
stability.
Proposition 3.4.1 (partial l -stability). Under the CFL condition (3.16), the staggered ex-
∞
plicit Godunov (d = d = 1/2) scheme (3.6) satisfies
p u
un+1 +un+1
n n
0 m i a N x 1(cid:12) pn i +1 ± i − 1/2,R 2ρ 0 c 0 i+1/2,L (cid:12) ≤ 0 m i a N x 1 max Ci±,L , Ci±,R ,
≤≤ − (cid:12) (cid:12) ≤≤ − (cid:16) (cid:16)(cid:12)(cid:16) (cid:17) (cid:12) (cid:12)(cid:16) (cid:17) (cid:12)(cid:17)(cid:17)
(cid:12) (cid:12) (cid:12) (cid:12) (cid:12) (cid:12)
(cid:12) (cid:12) (cid:12) (cid:12) (cid:12) (cid:12)
where Ci±,L and
C
(cid:12) i±,R are defined in Definition(cid:12)3.4.1. This stability is only partial because the left
term can be different from n+1 . Indeed, we recall that un+1 = un+1 +un+1 /2.
Ci± i 1/2 i 1/2,L i 1/2,R
− − −
(cid:16) (cid:17)
(cid:0) (cid:1)
Proof. Adding 1/(ρ c ) times (3.10b) to (3.10a), the explicit staggered Godunov scheme (d =
0 0 p
d = 1/2) gives
u
2 ∆ ∆ x t (cid:32) (cid:18) p i,L + u i ρ − 0 1 c / 0 2,R (cid:19) n+1 − (cid:18) p i + u ρ i − 0 c 1 0 /2 (cid:19) n (cid:33) +c 0 pn i − pn i − 1 = 0, (3.25)
(cid:0) (cid:1)
u u
Using that pn pn = pn+ i − 1/2 pn + i − 1/2 , we obtain
i − i 1 i ρ c − i 1 ρ c
− 0 0 (cid:18) − 0 0 (cid:19)
p i,L +
u
i − 1/2,R
n+1
= 1
2c
0
∆t
p i +
u
i − 1/2
n
+
2c
0
∆t
p i 1 +
u
i − 1/2
n
. (3.26)
(cid:18) ρ 0 c 0 (cid:19) (cid:18) − ∆x (cid:19)(cid:18) ρ 0 c 0 (cid:19) ∆x (cid:18) − ρ 0 c 0 (cid:19)
Moreover, adding 1/(ρ c ) times (3.11b) to (3.11a), we get
0 0
p i,R +
u
i+1/2,L
n+1
= 1
2c
0
∆t
p i +
u
i+1/2
n
+
2c
0
∆t
p i +
u
i − 1/2
n
. (3.27)
ρ c − ∆x ρ c ∆x ρ c
(cid:18) 0 0 (cid:19) (cid:18) (cid:19)(cid:18) 0 0 (cid:19) (cid:18) 0 0 (cid:19)
Then, multiplying (3.26) and (3.27) by 1/2 and summing both equations, we have
p i,L +p i,R + u i − 1/2,R +u i+1/2,L n+1 = c 0 ∆t p i 1 + u i − 1/2 n + 1 p i + u i − 1/2 n
(cid:18) 2 ρ 0 c 0 (cid:19) ∆x (cid:18) − ρ 0 c 0 (cid:19) 2 (cid:18) ρ 0 c 0 (cid:19)
+ 1 c 0 ∆t p + u i+1/2 n ,
i
2 − ∆x ρ c
(cid:18) (cid:19)(cid:18) 0 0 (cid:19)
that gives the result under CFL condition (3.16). Concerning the characteristic variable p
−
u/(ρ c ), subtracting 1/(ρ c ) times (3.10b) to (3.10a), the explicit staggered Godunov scheme
0 0 0 0
(d = d = 1/2) gives
p u
p i,L
u
i − 1/2,R
n+1
= 1
2c
0
∆t
p i
u
i − 1/2
n
+
2c
0
∆t
p i
u
i+1/2
n
. (3.28)
− ρ c − ∆x − ρ c ∆x − ρ c
(cid:18) 0 0 (cid:19) (cid:18) (cid:19)(cid:18) 0 0 (cid:19) (cid:18) 0 0 (cid:19)

CHAPTER 3. STUDY OF THE STABILITY OF THE STAGGERED SCHEMES FOR
THE ONE DIMENSIONAL WAVE SYSTEM 65
Moreover, subtracting 1/(ρ c ) times (3.11b) to (3.11a), we get
0 0
p
u
i+1/2,L
n+1
= 1
2c
0
∆t
p
u
i+1/2
n
+
2c
0
∆t
p
u
i+1/2
n
. (3.29)
i,R i i+1
− ρ c − ∆x − ρ c ∆x − ρ c
(cid:18) 0 0 (cid:19) (cid:18) (cid:19)(cid:18) 0 0 (cid:19) (cid:18) 0 0 (cid:19)
Then, multiplying (3.28) and (3.29) by 1/2 and summing both equations, we have
p i,L +p i,R u i − 1/2,R +u i+1/2,L n+1 = 1 c 0 ∆t p i u i − 1/2 n + 1 p i u i+1/2 n
2 − ρ c 2 − ∆x − ρ c 2 − ρ c
(cid:18) 0 0 (cid:19) (cid:18) (cid:19)(cid:18) 0 0 (cid:19) (cid:18) 0 0 (cid:19)
+
c
0
∆t
p
u
i+1/2
n
,
i+1
∆x − ρ c
(cid:18) 0 0 (cid:19)
that gives the result under CFL condition (3.16).
Remark 3.4.1. Using a finite volume interpretation on half-cells, Proposition 3.4.1 theoretic-
ally guarantees a certain stability of the solution. From a practical point of view, it has little
interest. Indeed, the left-hand quantity in the inequality of Proposition 3.4.1 can not be com-
puted because the staggered scheme (3.6) gives only the average pressure p on cell C and the
i i
average velocity u on cell C but not the half-cell values p , p u and u .
i+1/2 i+1/2 i,L i,R i+1/2,L i+1/2,R
3.4.2 Characteristic variables defined on primal cells
In this section, we use another definition of the characteristic variables for staggered grids such
that they are piecewise constant on primal cells taking into account an upwinding.
Definition 3.4.2 (Characteristic variables on primal cells). The characteristic variables
Ch±
are defined on primal cells C = x ;x as
i i 1/2 i+1/2
− (cid:98)
(cid:3)u (cid:2) u
Ci− = p i − ρ i − c 1/2 and Ci + = p i + ρ i+ c 1/2 .
0 0 0 0
Proposition 3.4.2. (cid:98)• The explicit staggered sche(cid:98)me (3.6) provides the following explicit
upwinding schemes on the characteristic variables defined in Definition 3.4.2
Ch±
∆x + n+1 + n +c (cid:98)+ n + n = 0,
∆t Ci − Ci 0 Ci − Ci 1
 (cid:18) (cid:19) −
∆x (cid:16) (cid:17)n+1 (cid:16) (cid:17)n (cid:16)(cid:16) (cid:17) n (cid:16) (cid:17)n(cid:17)
   ∆t C (cid:98) i− − C (cid:98) i− − c 0 C (cid:98) i−1 − (cid:98) Ci− = 0,
−
(cid:18) (cid:19)
(cid:16) (cid:17) (cid:16) (cid:17) (cid:16)(cid:16) (cid:17) (cid:16) (cid:17) (cid:17)

if and only if  d = 1 and (cid:98) d = 0. (cid:98) (cid:98) (cid:98)
p u
• Theimplicitstaggeredscheme (3.7)providesimplicitupwindschemesonthecharacteristic
variables if and only if d = 1 and d = 0.
Ch± p u
(cid:98)

CHAPTER 3. STUDY OF THE STABILITY OF THE STAGGERED SCHEMES FOR
| THE | ONE | DIMENSIONAL |     | WAVE | SYSTEM |     |     |     |     |     |     |     | 66  |
| --- | --- | ----------- | --- | ---- | ------ | --- | --- | --- | --- | --- | --- | --- | --- |
Proof. We firstly consider the explicit staggered scheme (3.6). Adding 1/(ρ c ) times (3.6b)
|     |            |       |     |     |       |     |          |     |       |     |     | 0 0   |     |
| --- | ---------- | ----- | --- | --- | ----- | --- | -------- | --- | ----- | --- | --- | ----- | --- |
| to  | (3.6a), we | get   |     |     |       |     |          |     |       |     |     |       |     |
|     |            | u     | n+1 |     | u     | n   |          |     | u     | n   |     | u     | n   |
|     | ∆x         | i+1/2 |     |     | i+1/2 |     |          |     | i+1/2 |     |     | i 1/2 |     |
|     | p          | +     |     |     | p +   |     | +c       | p   | +     |     | p   | + −   |     |
|     | ∆t         | i ρ c |     |     | i ρ   | c   |          | 0   | i+1 ρ | c   |     | i ρ c |     |
|     | (cid:32)   | 0     | 0   | −   | 0     | 0   | (cid:33) |     | 0     | 0   | −   | 0 0   |     |
(cid:18) (cid:19) (cid:18) (cid:19) (cid:18)(cid:18) (cid:19) (cid:18) (cid:19) (cid:19)
d
|     |     |     |     |     |            |     |         | n     | u         |     |          |      | n       |
| --- | --- | --- | --- | --- | ---------- | --- | ------- | ----- | --------- | --- | -------- | ---- | ------- |
|     |     |     |     | =   | c 0 d p (p | i+1 | 2p i +p | i 1 ) | + u i+3/2 |     | 2u i+1/2 | +u i | 1/2 ,   |
|     |     |     |     |     |            | −   |         |       | ρ         | −   |          |      |         |
|     |     |     |     |     |            |     |         | −     | 0         |     |          | −    |         |
|     |     |     |     |     |            |     |         |       | (cid:0)   |     |          |      | (cid:1) |
so that
|     | ∆x       | u     | n+1 |     | u     | n   |          |     | u     | n   |     | u     | n   |
| --- | -------- | ----- | --- | --- | ----- | --- | -------- | --- | ----- | --- | --- | ----- | --- |
|     |          | i+1/2 |     |     | i+1/2 |     |          |     | i+1/2 |     |     | i 1/2 |     |
|     | p        | +     |     |     | p +   |     | +c       | p   | +     |     | p   | + −   |     |
|     | ∆t       | i ρ c |     | −   | i ρ   | c   |          | 0   | i ρ c | −   | i   | 1 ρ c |     |
|     | (cid:32) | 0     | 0   |     | 0     | 0   | (cid:33) |     | 0 0   |     | −   | 0 0   |     |
(cid:18) (cid:19) (cid:18) (cid:19) (cid:18)(cid:18) (cid:19) (cid:18) (cid:19) (cid:19)
d
|     |     |     | =   | c (d | 1)(p |     | 2p +p | ) n | + u u |     | 2u    | +u  | n , |
| --- | --- | --- | --- | ---- | ---- | --- | ----- | --- | ----- | --- | ----- | --- | --- |
|     |     |     |     | 0    | p    | i+1 | i     | i 1 | i+3/2 |     | i+1/2 | i   | 1/2 |
|     |     |     |     |      | −    | −   |       | −   | ρ     | −   |       | −   |     |
0
|     |     |     |     |     |     |     |     |     | (cid:0) |     |     |     | (cid:1) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | ------- | --- | --- | --- | ------- |
that gives the result if and only if d = 1 and d = 0. Subtracting 1/(ρ c ) times (3.6b)
|     |             |       |       |     | p   |     |     | u   |     |     |     | 0 0 |     |
| --- | ----------- | ----- | ----- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| to  | (3.6a) with | i = i | 1, we | get |     |     |     |     |     |     |     |     |     |
−
|     |     | u     | n+1 |     | u   | n   |     |     | u     | n   |     | u     | n   |
| --- | --- | ----- | --- | --- | --- | --- | --- | --- | ----- | --- | --- | ----- | --- |
|     | ∆x  | i 1/2 |     |     | i   | 1/2 |     |     | i+1/2 |     |     | i 1/2 |     |
|     | p   | −     |     |     | p − |     | c   | p   |       |     | p   | −     |     |
|     |     | i     |     |     | i   |     |     | 0   | i+1   |     |     | i     |     |
∆t (cid:32) − ρ 0 c 0 − − ρ 0 c 0 (cid:33) − − ρ 0 c 0 − − ρ 0 c 0
(cid:18) (cid:19) (cid:18) (cid:19) (cid:18)(cid:18) (cid:19) (cid:18) (cid:19) (cid:19)
d
|     |     |     |     |        |        |     |         | n     | u         |     |          |      | n       |
| --- | --- | --- | --- | ------ | ------ | --- | ------- | ----- | --------- | --- | -------- | ---- | ------- |
|     |     |     | =   | c 0 (d | p 1)(p | i+1 | 2p i +p | i 1 ) | + u i+3/2 |     | 2u i+1/2 | +u i | 1/2 ,   |
|     |     |     |     |        | −      | −   |         |       | ρ         | −   |          |      |         |
|     |     |     |     |        |        |     |         | −     | 0         |     |          | −    |         |
|     |     |     |     |        |        |     |         |       | (cid:0)   |     |          |      | (cid:1) |
|     |     |     |     |        | d      |     | d       |       |           |     |          |      |         |
that gives the result if and only if p = 1 and u = 0. Same results are obtained for the implicit
| staggered | scheme | by  | considering |     | (3.7) | instead | of (3.6). |     |     |     |     |     |     |
| --------- | ------ | --- | ----------- | --- | ----- | ------- | --------- | --- | --- | --- | --- | --- | --- |
Corollary 3.4.2. The staggered explicit (resp. implicit) pressure-centered scheme with d = 1
p
and d = 0 is l -stable and TVD on the characteristic variables under CFL condition (3.17)
|        | u       | ∞    |      |     |     |     |     |     | Ch± |     |     |     |     |
| ------ | ------- | ---- | ---- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| (resp. | for all | ∆t > | 0.). |     |     |     |     |     |     |     |     |     |     |
(cid:98)
Remark 3.4.3. By contrast with usual properties of collocated schemes, getting discrete trans-
port equations on the characteristic variables does not imply energy dissipation. Indeed, the
staggered explicit pressure-centered scheme with d = 1 gives upwind schemes on the character-
p
istic variables (see Proposition 3.4.2) but is not energy dissipative (see Proposition 3.3.8).
±
C
(cid:98)
| 3.5 | Discussion |      | on     | some | preexisting |     |     | staggered |     | schemes |     | through |     |
| --- | ---------- | ---- | ------ | ---- | ----------- | --- | --- | --------- | --- | ------- | --- | ------- | --- |
|     | low        | Mach | number |      | asymptotics |     |     |           |     |         |     |         |     |
Performingthetwotimescalesasymptoticexpansionintroducedinchapter2,itcanbeobtained
that most of the staggered schemes on Euler system are asymptotically consistent with a fully-
centered (d = d = 0) staggered discretization (3.5) of the wave system. The fully discrete low
|     |     | p u |     |     |     |     |     |     |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |

CHAPTER 3. STUDY OF THE STABILITY OF THE STAGGERED SCHEMES FOR
| THE | ONE | DIMENSIONAL |     | WAVE | SYSTEM |     |     |     |     |     | 67  |
| --- | --- | ----------- | --- | ---- | ------ | --- | --- | --- | --- | --- | --- |
MachnumberasymptoticanalysisontwoexamplesareconductedinAppendixC;theotherswill
follow equivalently. The particular case of one space dimension masks some differences between
schemes, and as a consequence they can mainly be sorted into three families, depending only
| on  | the time   | integration; | the              | ones | asymptotically  | consistent |        | with:      |              |         |     |
| --- | ---------- | ------------ | ---------------- | ---- | --------------- | ---------- | ------ | ---------- | ------------ | ------- | --- |
|     | • Implicit |              | (fully-centered) |      | discretizations |            | (3.7): |            |              |         |     |
|     |            |              |                  |      |                 |            |        | we mention | the implicit | schemes |     |
of [49, 66, 81]. In [66], discrete entropy inequalities are obtained for the barotropic and
full Euler system, while in [81] the authors study the convergence of the scheme at low
|     | Mach | number | limit. |     |     |     |     |     |     |     |     |
| --- | ---- | ------ | ------ | --- | --- | --- | --- | --- | --- | --- | --- |
•
ImEx (fully-centered) discretizations (3.8): for example, the scheme of [65] and the
semi-implicit scheme of [82, 69]. In [65], the scheme is described as an explicit scheme
because the two equations of (3.8) can be decoupled and resolved sequentially. In [65],
the stability and the consistency of the scheme for the barotropic Euler system is studied.
In [82, 69], the energy stability is studied for shallow water equations. Generally, the
|     | results    | are | in line with    | the | ones obtained | here.     |     |                    |     |     |       |
| --- | ---------- | --- | --------------- | --- | ------------- | --------- | --- | ------------------ | --- | --- | ----- |
|     | • Explicit |     | discretizations |     | with          | diffusion |     | on both equations: |     |     |       |
|     |            |     |                 |     | (3.6)         |           |     |                    | To  | our | know- |
ledge, the only reference that proposes an explicit scheme is [69]. The authors introduce
stabilizationcoefficientsγ andαthatarelinkedtothecoefficientsd andd ofthepresent
p u
|     | chapter | with | the following |     | relation |     |     |         |     |     |     |
| --- | ------- | ---- | ------------- | --- | -------- | --- | --- | ------- | --- | --- | --- |
|     |         |      |               |     | c ∆t     |     |     | c ∆t    |     |     |     |
|     |         |      |               |     | 0        |     |     | 0       |     |     |     |
|     |         |      |               |     | d = γ    | and |     | d = α . |     |     |     |
|     |         |      |               |     | p ∆x     |     |     | u ∆x    |     |     |     |
In [69], the stabilization terms are factor of c 0 ∆t/∆x which does not seem to be classical,
∆t
because it corresponds to add a dependency of the numerical flux on the time step or
the CFL number. It also means that the stabilization coefficients depend on the mesh.
Thestabilizationcoefficientsαandγ arethencalibratedtoensureaglobalenergydecrease
in the non-linear case, which is completed with the study of the von-Neumann necessary
condition in the linearized case. In the present chapter, the study provides stabilizations
terms d p and d u , independent of the mesh, that guarantee energy dissipation under some
|     | CFL | condition. |         |     |     |     |     |     |     |     |     |
| --- | --- | ---------- | ------- | --- | --- | --- | --- | --- | --- | --- | --- |
| 3.6 |     | Numerical  | results |     |     |     |     |     |     |     |     |
In this section, numerical results are performed. Firstly, amplification matrices are studied
numerically to illustrate the stability results of section 3.2. Finally, numerical tests are per-
formed to analyze the behavior of the different schemes. From a numerical point of view, the
parameters c 0 , ρ 0 and κ 0 of the linear wave system (3.2) are set equal to 1. We recall that the
computational domain is [0,1] equipped with periodic boundary conditions.

CHAPTER 3. STUDY OF THE STABILITY OF THE STAGGERED SCHEMES FOR
| THE ONE | DIMENSIONAL |       | WAVE | SYSTEM            |     |     |          | 68  |
| ------- | ----------- | ----- | ---- | ----------------- | --- | --- | -------- | --- |
| 3.6.1   | Numerical   | study | of   | the amplification |     |     | matrices |     |
The aim of this section is to illustrate numerically the results of section 3.2 obtained on the
von-Neumann necessary condition and the energy dissipation. The domain is meshed with 200
| uniform  | cells.  |             |     |           |     |           |     |     |
| -------- | ------- | ----------- | --- | --------- | --- | --------- | --- | --- |
| Spectral | radius: | von-Neumann |     | necessary |     | condition |     |     |
Inthissection,weproposetotestnumericallythevon-Neumannnecessaryconditionbyplotting
the spectral radius of the amplification matrix with respect to the CFL number for different
values of the frequency values k. Following Definition 3.3.2, the spectral radius should be lower
| than 1 | for ensuring | the von-Neumann |     |     | necessary | condition. |     |     |
| ------ | ------------ | --------------- | --- | --- | --------- | ---------- | --- | --- |
In Figure 3.1, results obtained with staggered schemes are shown. The two top pictures
give the results for the explicit, fully upwind scheme (d p = d u = 1/2), and explicit, pressure
d
centered scheme ( p = 1). In agreement with Proposition 3.3.4, these schemes are respectively
ensuring the von-Neumann necessary condition for CFL 1/2 and 1. The two middle positioned
pictures show the results for implicit time integration. In agreement with Proposition 3.3.5,
these schemes ensure the von-Neumann necessary condition unconditionally. In the bottom
pictures of Figure 3.1, results for the implicit-explicit pressure-centered scheme with d = 1/2
p
and the fully-centered scheme are shown. In agreement with Proposition 3.3.6, the first one
ensuresthevon-NeumannnecessaryconditionforCFLnumberslowerthan(√5
1)/2 0.5618
− ≈
while for the fully-centered scheme, maximal CFL number ensuring the condition is 1.
| Norm | of the amplification |     | matrices: |     | energy |     | dissipation |     |
| ---- | -------------------- | --- | --------- | --- | ------ | --- | ----------- | --- |
In this section, we would like to assess the energy dissipation property of the staggered schemes
for different time integrations by plotting the norm S 1/2A(k)S1/2 with respect to the
−
2
| CFL number, | and | for different | k.  |     |     |     |     |     |
| ----------- | --- | ------------- | --- | --- | --- | --- | --- | --- |
(cid:12)(cid:12)(cid:12) (cid:12)(cid:12)(cid:12)
In Figure 3.2, the norm of S 1/2A(k)S1/2 i(cid:12)s(cid:12)(cid:12)plotted with di(cid:12)ff(cid:12)(cid:12)erent time integration
|     |     |     |     | −   |     |     | 2   |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
types and different frequency numbers k. In the top left picture, the norm for the staggered
|     |     |     | (cid:12)(cid:12)(cid:12) |     |     | (cid:12)(cid:12)(cid:12) |     |     |
| --- | --- | --- | ------------------------ | --- | --- | ------------------------ | --- | --- |
explicit fully-upwind scheme wit(cid:12)(cid:12)h(cid:12)d = d = 1/2,(cid:12)(cid:12)w(cid:12) ith explicit time stepping is plotted, which
|     |     |     |     | p   | u   |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
shows in agreement with Proposition 3.3.1, that the norm of the matrix remains lower than 1
1/2,
for CFL lower than which means that the scheme is energy dissipative, in agreement with
Proposition 3.3.7. In the same manner, the bottom-left picture (implicit-explicit, pressure-
centered with d = 1/2) shows that the norm of the matrix is lower than 1 for CFL lower than
p
1/2, in agreement with Proposition 3.3.10. Still in Figure 3.2, the top-right picture (explicit,
pressure centered with d = 1), and the bottom-right picture (implicit-explicit, fully-centered)
p
1/2A(k)S1/2
show that S is always greater than 1, which means that these schemes are
|     | −   |     | 2   |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
never energy dissipative, whatever the CFL number. These pictures are in agreement with
|     | (cid:12)(cid:12)(cid:12) |     | (cid:12)(cid:12)(cid:12) |     |     |     |     |     |
| --- | ------------------------ | --- | ------------------------ | --- | --- | --- | --- | --- |
Proposition(cid:12)(cid:12)3(cid:12).3.8andPropos(cid:12)(cid:12)it(cid:12)ion3.3.11. Last,thetwomiddlepositionedpicturesofFigure3.2,
whichmatcheswiththeimplicittimeintegrationofthepressure-centeredandthefully-centered
staggered schemes show a norm of the matrix S 1/2A(k)S1/2 always lower than 1, which means
−
that these schemes are unconditionally stables. This is in agreement with Proposition 3.3.9.

CHAPTER 3. STUDY OF THE STABILITY OF THE STAGGERED SCHEMES FOR
| THE   | ONE DIMENSIONAL |     |               | WAVE | SYSTEM |         |     |     |     |     | 69  |
| ----- | --------------- | --- | ------------- | ---- | ------ | ------- | --- | --- | --- | --- | --- |
| 3.6.2 | Tests           | on  | the numerical |      |        | schemes |     |     |     |     |     |
In this section, we test the different numerical schemes on two different test cases. The results
presented below are obtained with a uniform mesh of 200 cells and with a CFL number equal
to 0.45. In each test case the pressure is also initialized using a midpoint formula while the
| velocity | is initialized |             | by computing |     | the | value | at the | face. |     |     |     |
| -------- | -------------- | ----------- | ------------ | --- | --- | ----- | ------ | ----- | --- | --- | --- |
| On       | energy         | dissipation |              |     |     |       |        |       |     |     |     |
The aim of this test case is to illustrate the property of energy dissipation for the different
staggered schemes. For a given initial condition, the energy is plotted with respect to time.
The initial solution that is used is given by (3.22)-(3.23), namely, this matches with the mode
inducinganincreaseoftheenergyforthestaggered,pressure-centeredschemewithexplicittime
stepping. The parameter k is set to 5. The results are shown in the right picture of Figure 3.3.
Forthestaggeredexplicitfully-upwindscheme,thestaggeredimplicit-explicitpressure-centered
schemeandallthestaggeredimplicitschemes,weseethattheenergyisalwaysdecreasing,which
is in agreement with Proposition 3.3.7, Proposition 3.3.9 and Proposition 3.3.10. Concerning
the explicit, pressure centred scheme, we observe that the energy increases (until nearly time
0.02), then decreases until approximately time 0.08, and then restarts increasing. This is
in agreement with Proposition 3.3.8. Last, for the staggered implicit-explicit fully-centered
scheme, the energy starts decreasing, which is not surprising, because the initial condition does
not match a priori with an increasing mode. Still, at around time 0.02, the energy starts
increasing, which matches with what was proven in Proposition 3.3.11.
| On  | l -stability |     | and TVD | property: |     | oscillating |     | schemes |     |     |     |
| --- | ------------ | --- | ------- | --------- | --- | ----------- | --- | ------- | --- | --- | --- |
∞
After having extensively studied the energy dissipation property, we wish to study whether the
schemes are oscillating or not. For that purpose, the initial condition considered is a Riemann
| problem | given    | by  |              |          |             |      |       |       |          |           |        |
| ------- | -------- | --- | ------------ | -------- | ----------- | ---- | ----- | ----- | -------- | --------- | ------ |
|         |          |     |              | p ,      | x           | 1/2  |       |       | u ,      | x 1/2     |        |
|         |          |     | p0(x)        | L        |             |      | u0(x) |       | L        |           |        |
|         |          |     | =            |          | ≤           |      |       | =     |          | ≤         | (3.30) |
|         |          |     |              | p        | , otherwise |      |       |       | u ,      | otherwise |        |
|         |          |     |              | R        |             |      |       |       | R        |           |        |
|         |          |     |              | (cid:26) |             |      |       |       | (cid:26) |           |        |
| where   | the left | and | right states | are      | given       | by   |       |       |          |           |        |
|         |          |     |              |          | p = u       | = 1, | p     | = u = | 1.       |           |        |
|         |          |     |              |          | L           | L    | R     | R     |          |           |        |
−
Since p u /(ρ c ) = p u /(ρ c ) = 0, this test case corresponds to a simple advection
|     | L                  | L   | 0 0      | R R    | 0   | 0    |     |     |     |     |     |
| --- | ------------------ | --- | -------- | ------ | --- | ---- | --- | --- | --- | --- | --- |
|     | −                  |     |          | −      |     |      |     |     |     |     |     |
| on  | the characteristic |     | variable | p+u/(ρ |     | c ). |     |     |     |     |     |
0 0
In Figure 3.4, the pressure p, velocity u and the characteristic variables (as defined in
Ch±
Definition 3.4.2) obtained with the explicit, implicit and implicit-explicit staggered schemes are
oscillatio(cid:98)ns
plotted. Numerical solutions are consistent with the exact one but appear for all
the staggered schemes. While the implicit-explicit fully-centered scheme appears to be highly
oscillating, for the other schemes the number of oscillations appears to be limited. Indeed, for
these different schemes, only two oscillations seem to appear following the numerical solution

CHAPTER 3. STUDY OF THE STABILITY OF THE STAGGERED SCHEMES FOR
THE ONE DIMENSIONAL WAVE SYSTEM 70
of the Riemann problem, one localized in x = 0.4 and the other in x = 0.9. In Figure 3.5,
the l -norm and the total variation of the characteristic variables and + are plotted with
∞ Ch− Ch
respecttotime. Inagreementwith3.4.2,theexplicitandtheimplicitpressure-centeredschemes
with d = 1 are l -stable and TVD on the characteristic variables(cid:98) and(cid:98) +. For the other
p ∞ Ch− Ch
schemes, itisnotthecase. However, exceptfortheimplicit-explicitfully-centeredschemethese
quantities seem to be controlled over time. Based on the numerical(cid:98)results,(cid:98)we infer that the
energy dissipative schemes have a lower amplitude of the oscillations than the ones that are not
energydissipative. InFigure3.6, theenergyofthenumericalsolutionisplottedasafunctionof
the time. Energy dissipation is not a sufficient condition to guarantee a non-oscillating scheme.
Indeed, for this test case, the energy decreases for all the schemes except the implicit-explicit
fully-centered staggered scheme but all these schemes are oscillating.
3.7 Conclusion
In this chapter, several kinds of stability were investigated: the von–Neumann necessary con-
dition, the energy dissipation and the l -stability on the characteristic variables of staggered
∞
schemeswithdifferentupwindinganddifferenttimeintegrations(explicit,implicitandimplicit-
explicit).
It is known that the necessary von-Neumann condition becomes a sufficient condition only
when the amplification matrix is normal. In this chapter, several examples show that when
the amplification matrix is not normal, the von-Neumann necessary condition may be ensured
eventhoughschemeisnotenergydissipative. Forexample, staggeredexplicitpressure-centered
schemes satisfy the necessary von–Neumann stability condition under some CFL condition but
are not energy dissipative. In this case we were able to constuct examples for which the energy
increases on the first iteration in Proposition 3.3.8. Switching to semi-implicit or implicit time
integration generally result in an energy-dissipative scheme. These results are in line with the
literature about discrete entropy inequality [66, 81]. We draw the following broad conclusions
on the energy dissipation
• for explicit time stepping, staggered schemes require numerical diffusion to ensure energy
dissipation; diffusion is needed in the form of both pressure terms and velocity terms
whether it be as spatial diffusion or numerical diffusion arising from the time integration
(typically ImEx, Implicit),
• whileforimplicittimestepping, theenergydissipationmaybeensuredwithoutnumerical
diffusion.
However, ensuring the energy dissipation is not sufficient to guarantee that the scheme will
not be oscillatory, so for a complete analysis, l -stability and TVD property on the character-
∞
istic variables have been be studied. Regarding this analysis we observe that:
• The most popular time integrations in the staggered scheme community give, in the low
Mach number asymptotics, either the fully-centred Implicit or the fully-centred ImEx

CHAPTER 3. STUDY OF THE STABILITY OF THE STAGGERED SCHEMES FOR
THE ONE DIMENSIONAL WAVE SYSTEM 71
schemes: if the implicit one is energy dissipative, the implicit-explicit one is not and both
schemes are oscillating.
• Thel -stabilityandTVDpropertydependonthedefinitionofthecharacteristicvariables
∞
that is not obvious to define for staggered schemes because the pressure and the velocity
are not defined in the same location. A definition of the characteristic variables on
primal cells for staggered schemes is proposed but other choice could be made. Still,
as the numerical results on the original staggered variables are oscillating, it is likely
that the schemes cannot be TVD nor l stable, whatever the choice of reconstruction of
∞
characteristic variables.
• Different upwinding and time integrators have been proposed and tested. If all staggered
schemes seem to be a bit oscillating, adding numerical dissipation allows to reduce the
amplitude of these oscillations.
To conclude, we proposed an explicit fully-upwind staggered scheme that is energy dissipative
and provides numerical results comparable to these obtained with implicit upwind schemes
and the implicit-explicit pressure-centered scheme. This scheme can simply be interpreted as
a classical Godunov scheme written in half-cells. All in all, these three schemes will be our
starting point for the extension on the multi-dimensional wave system.

CHAPTER 3. STUDY OF THE STABILITY OF THE STAGGERED SCHEMES FOR
| THE ONE | DIMENSIONAL | WAVE SYSTEM |      |     |     | 72  |     |
| ------- | ----------- | ----------- | ---- | --- | --- | --- | --- |
| k=1     |             |             | 1.75 |     |     |     |     |
3.0 k=11
k=21
1.50
k=31
k=41
| 2.5 |     |     | 1.25 |     |     |     |     |
| --- | --- | --- | ---- | --- | --- | --- | --- |
k=51
| ))k(A(ρ |     |     | ))k(A(ρ |     |     |     | k=1 |
| ------- | --- | --- | ------- | --- | --- | --- | --- |
k=61
| k=71 |     |     |      |     |     |     | k=11 |
| ---- | --- | --- | ---- | --- | --- | --- | ---- |
| 2.0  |     |     | 1.00 |     |     |     |      |
| k=81 |     |     |      |     |     |     | k=21 |
k=31
k=91
|     |     |     | 0.75 |     |     |     | k=41 |
| --- | --- | --- | ---- | --- | --- | --- | ---- |
1.5
k=51
|     |     |     | 0.50 |     |     |     | k=61 |
| --- | --- | --- | ---- | --- | --- | --- | ---- |
k=71
1.0
|     |     |     | 0.25 |     |     |     | k=81 |
| --- | --- | --- | ---- | --- | --- | --- | ---- |
k=91
| 0.2 | 0.4 0.6 0.8 | 1.0 1.2 | 1.4 | 0.2 0.4 | 0.6 0.8 | 1.0 1.2 | 1.4 |
| --- | ----------- | ------- | --- | ------- | ------- | ------- | --- |
|     | CFL         |         |     |         | CFL     |         |     |
Explicit, fully-upwind with d = d = 1/2 Explicit, pressure-centered with d = 1
|     |     | p u |     |     |     | p   |     |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 1.0 |     |     | 1.0 |     |     |     |     |
| 0.9 |     |     | 0.9 |     |     |     |     |
0.8
0.8
| k=1         |     |     |         | k=1  |     |     |     |
| ----------- | --- | --- | ------- | ---- | --- | --- | --- |
| ))k(A(ρ 0.7 |     |     | ))k(A(ρ |      |     |     |     |
| k=11        |     |     | 0.7     | k=11 |     |     |     |
| k=21        |     |     |         | k=21 |     |     |     |
0.6
| k=31 |     |     | 0.6 | k=31 |     |     |     |
| ---- | --- | --- | --- | ---- | --- | --- | --- |
| k=41 |     |     |     | k=41 |     |     |     |
0.5
| k=51 |     |     |     | k=51 |     |     |     |
| ---- | --- | --- | --- | ---- | --- | --- | --- |
| k=61 |     |     | 0.5 | k=61 |     |     |     |
0.4
| k=71 |     |     |     | k=71 |     |     |     |
| ---- | --- | --- | --- | ---- | --- | --- | --- |
| k=81 |     |     | 0.4 | k=81 |     |     |     |
0.3
| k=91 |             |         |     | k=91    |         |         |     |
| ---- | ----------- | ------- | --- | ------- | ------- | ------- | --- |
| 0.2  | 0.4 0.6 0.8 | 1.0 1.2 | 1.4 | 0.2 0.4 | 0.6 0.8 | 1.0 1.2 | 1.4 |
|      | CFL         |         |     |         | CFL     |         |     |
Implicit, pressure-centered with d = 1 Implicit, fully-centered (d = d = 0)
|        |     | p   |     |      |     | p u |     |
| ------ | --- | --- | --- | ---- | --- | --- | --- |
| k=1    |     |     |     | k=1  |     |     |     |
| 8 k=11 |     |     |     | k=11 |     |     |     |
5
| k=21    |     |     |         | k=21 |     |     |     |
| ------- | --- | --- | ------- | ---- | --- | --- | --- |
| k=31    |     |     |         | k=31 |     |     |     |
| k=41    |     |     |         | k=41 |     |     |     |
| 6       |     |     | 4       |      |     |     |     |
| k=51    |     |     |         | k=51 |     |     |     |
| ))k(A(ρ |     |     | ))k(A(ρ |      |     |     |     |
| k=61    |     |     |         | k=61 |     |     |     |
| k=71    |     |     |         | k=71 |     |     |     |
| 4 k=81  |     |     |         | k=81 |     |     |     |
3
| k=91 |             |         |     | k=91    |         |         |     |
| ---- | ----------- | ------- | --- | ------- | ------- | ------- | --- |
| 2    |             |         | 2   |         |         |         |     |
| 0    |             |         | 1   |         |         |         |     |
| 0.2  | 0.4 0.6 0.8 | 1.0 1.2 | 1.4 | 0.2 0.4 | 0.6 0.8 | 1.0 1.2 | 1.4 |
|      | CFL         |         |     |         | CFL     |         |     |
Implicit-explicit, pressure-centered with d p = 1/2 Implicit-explicit, fully-centered
Figure 3.1: Spectral radius of the amplification matrices with explicit, implicit and implicit-
| explicit | time integrations. |     |     |     |     |     |     |
| -------- | ------------------ | --- | --- | --- | --- | --- | --- |

CHAPTER 3. STUDY OF THE STABILITY OF THE STAGGERED SCHEMES FOR
| THE ONE | DIMENSIONAL | WAVE SYSTEM |     |     |     | 73  |     |
| ------- | ----------- | ----------- | --- | --- | --- | --- | --- |
6
| k=1                                                      |     |     |                              | k=1                           |     |     |     |
| -------------------------------------------------------- | --- | --- | ---------------------------- | ----------------------------- | --- | --- | --- |
| 3.0 k=11                                                 |     |     |                              | k=11                          |     |     |     |
| k=21                                                     |     |     |                              | k=21                          |     |     |     |
| 2 (cid:12)(cid:12)(cid:12) (cid:12)(cid:12)(cid:12) k=31 |     |     | 2 5 (cid:12)(cid:12)(cid:12) | (cid:12)(cid:12)(cid:12) k=31 |     |     |     |
| k=41                                                     |     |     |                              | k=41                          |     |     |     |
| 2/1S)k(A2/1 2.5                                          |     |     | 2/1S)k(A2/1                  |                               |     |     |     |
| k=51                                                     |     |     |                              | k=51                          |     |     |     |
4
| k=61 |     |     |     | k=61 |     |     |     |
| ---- | --- | --- | --- | ---- | --- | --- | --- |
| k=71 |     |     |     | k=71 |     |     |     |
2.0
| k=81                                              |     |     |                            | k=81                     |     |     |     |
| ------------------------------------------------- | --- | --- | -------------------------- | ------------------------ | --- | --- | --- |
| k=91                                              |     |     | 3                          | k=91                     |     |     |     |
| −S 1.5                                            |     |     | −S                         |                          |     |     |     |
| (cid:12)(cid:12)(cid:12) (cid:12)(cid:12)(cid:12) |     |     | 2 (cid:12)(cid:12)(cid:12) | (cid:12)(cid:12)(cid:12) |     |     |     |
1.0
1
| 0.2 | 0.4 0.6 0.8 | 1.0 1.2 | 1.4 | 0.2 0.4 | 0.6 0.8 | 1.0 1.2 | 1.4 |
| --- | ----------- | ------- | --- | ------- | ------- | ------- | --- |
|     | CFL         |         |     |         | CFL     |         |     |
Explicit, fully-upwind with d = d = 1/2 Explicit, pressure-centered with d = 1
|     |     | p u |     |     |     | p   |     |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 1.0 |     |     | 1.0 |     |     |     |     |
0.9
| 2                                                     |     |     | 2                        |                          |     |     |     |
| ----------------------------------------------------- | --- | --- | ------------------------ | ------------------------ | --- | --- | --- |
| 0.9 (cid:12)(cid:12)(cid:12) (cid:12)(cid:12)(cid:12) |     |     | (cid:12)(cid:12)(cid:12) | (cid:12)(cid:12)(cid:12) |     |     |     |
| 2/1S)k(A2/1                                           |     |     | 2/1S)k(A2/1 0.8          |                          |     |     |     |
| k=1                                                   |     |     |                          | k=1                      |     |     |     |
0.8
| k=11 |     |     | 0.7 | k=11 |     |     |     |
| ---- | --- | --- | --- | ---- | --- | --- | --- |
| k=21 |     |     |     | k=21 |     |     |     |
| k=31 |     |     | 0.6 | k=31 |     |     |     |
0.7
| −S k=41                                                |     |     | −S                           | k=41                          |     |     |     |
| ------------------------------------------------------ | --- | --- | ---------------------------- | ----------------------------- | --- | --- | --- |
| k=51                                                   |     |     |                              | k=51                          |     |     |     |
| (cid:12)(cid:12)(cid:12) (cid:12)(cid:12)(cid:12) k=61 |     |     | 0.5 (cid:12)(cid:12)(cid:12) | (cid:12)(cid:12)(cid:12) k=61 |     |     |     |
| 0.6 k=71                                               |     |     |                              | k=71                          |     |     |     |
| k=81                                                   |     |     | 0.4                          | k=81                          |     |     |     |
| k=91                                                   |     |     |                              | k=91                          |     |     |     |
0.5
| 0.2 | 0.4 0.6 0.8 | 1.0 1.2 | 1.4 | 0.2 0.4 | 0.6 0.8 | 1.0 1.2 | 1.4 |
| --- | ----------- | ------- | --- | ------- | ------- | ------- | --- |
|     | CFL         |         |     |         | CFL     |         |     |
Implicit, pressure-centered with d = 1 Implicit, fully-centered (d = d = 0)
|      |     | p   |     |      |     | p u |     |
| ---- | --- | --- | --- | ---- | --- | --- | --- |
| 9    |     |     | 8   |      |     |     |     |
| k=1  |     |     |     | k=1  |     |     |     |
| k=11 |     |     |     | k=11 |     |     |     |
8
| k=21                                                     |     |     | 7                          | k=21                          |     |     |     |
| -------------------------------------------------------- | --- | --- | -------------------------- | ----------------------------- | --- | --- | --- |
| 2 (cid:12)(cid:12)(cid:12) (cid:12)(cid:12)(cid:12) k=31 |     |     | 2 (cid:12)(cid:12)(cid:12) | (cid:12)(cid:12)(cid:12) k=31 |     |     |     |
7
| k=41        |     |     | 6           | k=41 |     |     |     |
| ----------- | --- | --- | ----------- | ---- | --- | --- | --- |
| 2/1S)k(A2/1 |     |     | 2/1S)k(A2/1 |      |     |     |     |
| 6 k=51      |     |     |             | k=51 |     |     |     |
| k=61        |     |     | 5           | k=61 |     |     |     |
| k=71        |     |     |             | k=71 |     |     |     |
5
| k=81 |     |     |     | k=81 |     |     |     |
| ---- | --- | --- | --- | ---- | --- | --- | --- |
4
| 4 k=91                                            |     |     |                          | k=91                     |     |     |     |
| ------------------------------------------------- | --- | --- | ------------------------ | ------------------------ | --- | --- | --- |
| −S                                                |     |     | −S                       |                          |     |     |     |
| 3                                                 |     |     | 3                        |                          |     |     |     |
| (cid:12)(cid:12)(cid:12) (cid:12)(cid:12)(cid:12) |     |     | (cid:12)(cid:12)(cid:12) | (cid:12)(cid:12)(cid:12) |     |     |     |
| 2                                                 |     |     | 2                        |                          |     |     |     |
1
1
| 0.2 | 0.4 0.6 0.8 | 1.0 1.2 | 1.4 | 0.2 0.4 | 0.6 0.8 | 1.0 1.2 | 1.4 |
| --- | ----------- | ------- | --- | ------- | ------- | ------- | --- |
|     | CFL         |         |     |         | CFL     |         |     |
Implicit-explicit, pressure-centered with d p = 1/2 Implicit-explicit, fully-centered
S 1/2A(k)S1/2
Figure 3.2: Norm of the matrices − with explicit, implicit and implicit-
|||·|||2
| explicit | time integrations. |     |     |     |     |     |     |
| -------- | ------------------ | --- | --- | --- | --- | --- | --- |

CHAPTER 3. STUDY OF THE STABILITY OF THE STAGGERED SCHEMES FOR
| THE ONE | DIMENSIONAL | WAVE SYSTEM |     |     | 74  |
| ------- | ----------- | ----------- | --- | --- | --- |
1600
1400
))t(
|     |     | Explicit,fully-upwindwithd | =d =1/2 |     |     |
| --- | --- | -------------------------- | ------- | --- | --- |
|     | h   |                            | p u     |     |     |
U(E 1200
|     |     | Explicit,pressure-centeredwithd | =1  |     |     |
| --- | --- | ------------------------------- | --- | --- | --- |
p
|     |     | Implicit,fully-upwindwithd      | p =d u =1/2 |     |     |
| --- | --- | ------------------------------- | ----------- | --- | --- |
|     |     | Implicit,pressure-centeredwithd | =1          |     |     |
p
1000
Implicit,fully-centered
|     |     | ImEx,pressure-centeredwithd | =1/2 |     |     |
| --- | --- | --------------------------- | ---- | --- | --- |
p
800 ImEx,fully-centered
|     |     | 0.00 0.02 | 0.04 0.06 | 0.08 0.10 |     |
| --- | --- | --------- | --------- | --------- | --- |
t
Figure 3.3: Test case on energy dissipation: energy of the numerical solution with respect to
time.

CHAPTER 3. STUDY OF THE STABILITY OF THE STAGGERED SCHEMES FOR
| THE ONE | DIMENSIONAL | WAVE SYSTEM |     |     |     |     | 75  |
| ------- | ----------- | ----------- | --- | --- | --- | --- | --- |
| 1.0     |             |             |     | 1.5 |     |     |     |
1.0
| 0.5 | Explicit,fully-upwindwithdp=du=1/2 |     |     |     |     |     |     |
| --- | ---------------------------------- | --- | --- | --- | --- | --- | --- |
|     | Explicit,pressure-centeredwithdp=1 |     |     | 0.5 |     |     |     |
ImEx,pressure-centeredwithdp=1/2
Implicit,fully-upwindwithdp=du=1/2
| hp 0.0 |     |     |     | hp 0.0 | ImEx,fully-centered |     |     |
| ------ | --- | --- | --- | ------ | ------------------- | --- | --- |
Implicit,pressure-centeredwithdp=1
exactsolution
|     | Implicit,fully-centered |     |     | 0.5 |     |     |     |
| --- | ----------------------- | --- | --- | --- | --- | --- | --- |
−
| 0.5 | exactsolution |     |     |     |     |     |     |
| --- | ------------- | --- | --- | --- | --- | --- | --- |
| −   |               |     |     | 1.0 |     |     |     |
−
| − 1.0 |     |     |     | 1.5 |     |     |     |
| ----- | --- | --- | --- | --- | --- | --- | --- |
−
| 0.0 | 0.2 0.4 | 0.6 0.8 | 1.0 | 0.0 0.2 | 0.4 | 0.6 | 0.8 1.0 |
| --- | ------- | ------- | --- | ------- | --- | --- | ------- |
|     | x       |         |     |         | x   |     |         |
1.0
1.5
1.0
| 0.5 | Explicit,fully-upwindwithdp=du=1/2 |     |     |     |     |     |     |
| --- | ---------------------------------- | --- | --- | --- | --- | --- | --- |
|     | Explicit,pressure-centeredwithdp=1 |     |     | 0.5 |     |     |     |
ImEx,pressure-centeredwithdp=1/2
Implicit,fully-upwindwithdp=du=1/2
| hu 0.0 |     |     |     | hu 0.0 | ImEx,fully-centered |     |     |
| ------ | --- | --- | --- | ------ | ------------------- | --- | --- |
Implicit,pressure-centeredwithdp=1
exactsolution
|     | Implicit,fully-centered |     |     | − 0.5 |     |     |     |
| --- | ----------------------- | --- | --- | ----- | --- | --- | --- |
0.5
| −   | exactsolution |     |     |     |     |     |     |
| --- | ------------- | --- | --- | --- | --- | --- | --- |
1.0
−
| 1.0 |         |         |     | 1.5     |                                  |     |         |
| --- | ------- | ------- | --- | ------- | -------------------------------- | --- | ------- |
| −   |         |         |     | −       |                                  |     |         |
| 0.0 | 0.2 0.4 | 0.6 0.8 | 1.0 | 0.0 0.2 | 0.4                              | 0.6 | 0.8 1.0 |
|     | x       |         |     |         | x                                |     |         |
|     |         |         |     | 0.6     | ImEx,pressure-centeredwithdp=1/2 |     |         |
0.1
ImEx,fully-centered
0.4
exactsolution
0.2
0.0
Explicit,fully-upwindwithdp=du=1/2
| −hC                                    |     |     |     | −hC   |     |     |     |
| -------------------------------------- | --- | --- | --- | ----- | --- | --- | --- |
| b Explicit,pressure-centeredwithdp=1   |     |     |     | 0.0 b |     |     |     |
| 0.1 Implicit,fully-upwindwithdp=du=1/2 |     |     |     |       |     |     |     |
| −                                      |     |     |     | 0.2   |     |     |     |
| Implicit,pressure-centeredwithdp=1     |     |     |     | −     |     |     |     |
| Implicit,fully-centered                |     |     |     | 0.4   |     |     |     |
| 0.2                                    |     |     |     | −     |     |     |     |
| − exactsolution                        |     |     |     |       |     |     |     |
0.6
−
| 0.0 | 0.2 0.4 | 0.6 0.8 | 1.0 | 0.0 0.2 | 0.4 | 0.6 | 0.8 1.0 |
| --- | ------- | ------- | --- | ------- | --- | --- | ------- |
|     | x       |         |     |         | x   |     |         |
3
2
2
Explicit,fully-upwindwithdp=du=1/2
1
|     | Explicit,pressure-centeredwithdp=1 |     |     | 1   |     |     |     |
| --- | ---------------------------------- | --- | --- | --- | --- | --- | --- |
ImEx,pressure-centeredwithdp=1/2
Implicit,fully-upwindwithdp=du=1/2
| +      |     |     |     | + 0  |                     |     |     |
| ------ | --- | --- | --- | ---- | ------------------- | --- | --- |
| hC 0 b |     |     |     | hC b | ImEx,fully-centered |     |     |
Implicit,pressure-centeredwithdp=1
exactsolution
|     | Implicit,fully-centered |     |     | − 1 |     |     |     |
| --- | ----------------------- | --- | --- | --- | --- | --- | --- |
| 1   | exactsolution           |     |     |     |     |     |     |
| −   |                         |     |     | 2   |     |     |     |
−
3
| 2   |     |     |     | −   |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- |
−
| 0.0 | 0.2 0.4 | 0.6 0.8 | 1.0 | 0.0 0.2 | 0.4 | 0.6 | 0.8 1.0 |
| --- | ------- | ------- | --- | ------- | --- | --- | ------- |
|     | x       |         |     |         | x   |     |         |
Figure 3.4: Riemann problem : pressure p , velocity u and characteristic variables and +
|     |     |     | h   | h   |     | Ch− | Ch  |
| --- | --- | --- | --- | --- | --- | --- | --- |
(defined in Definition 3.4.2) obtained with implicit and explicit staggered schemes (left) and
with implicit-explicit staggered schemes (right) at time t = 0.1. (cid:98) (cid:98)

CHAPTER 3. STUDY OF THE STABILITY OF THE STAGGERED SCHEMES FOR
| THE ONE | DIMENSIONAL | WAVE SYSTEM |     |     | 76  |
| ------- | ----------- | ----------- | --- | --- | --- |
2.00 Explicit,fully-upwindwithdp=du=1/2 2.00 Explicit,fully-upwindwithdp=du=1/2
1.75 Explicit,pressure-centeredwithdp=1 1.75 Explicit,pressure-centeredwithdp=1
Implicit,fully-upwindwithdp=du=1/2 Implicit,fully-upwindwithdp=du=1/2
| 1.50 |     |     | 1.50 |     |     |
| ---- | --- | --- | ---- | --- | --- |
Implicit,pressure-centeredwithdp=1 Implicit,pressure-centeredwithdp=1
| 1.25 |     |     | 1.25 |     |     |
| ---- | --- | --- | ---- | --- | --- |
∞k−hCk Implicit,fully-centered ∞k−hCk Implicit,fully-centered
| 1.00 |                                  |     | 1.00 |     |     |
| ---- | -------------------------------- | --- | ---- | --- | --- |
| b    | ImEx,pressure-centeredwithdp=1/2 |     | b    |     |     |
ImEx,fully-centered
| 0.75      |           |           | 0.75      |                                    |           |
| --------- | --------- | --------- | --------- | ---------------------------------- | --------- |
| 0.50      |           |           | 0.50      |                                    |           |
| 0.25      |           |           | 0.25      |                                    |           |
| 0.00      |           |           | 0.00      |                                    |           |
| 0.00 0.02 | 0.04 0.06 | 0.08 0.10 | 0.00 0.02 | 0.04 0.06                          | 0.08 0.10 |
|           | t         |           |           | t                                  |           |
| 3.4       |           |           |           | Explicit,fully-upwindwithdp=du=1/2 |           |
2.30
Explicit,pressure-centeredwithdp=1
3.2 Explicit,fully-upwindwithdp=du=1/2
|     |     |     | 2.25 | Implicit,fully-upwindwithdp=du=1/2 |     |
| --- | --- | --- | ---- | ---------------------------------- | --- |
3.0 Explicit,pressure-centeredwithdp=1
Implicit,pressure-centeredwithdp=1
Implicit,fully-upwindwithdp=du=1/2 2.20
| ∞khCk 2.8                            |     |     | ∞khCk  | Implicit,fully-centered |     |
| ------------------------------------ | --- | --- | ------ | ----------------------- | --- |
| + Implicit,pressure-centeredwithdp=1 |     |     | +      |                         |     |
| b                                    |     |     | 2.15 b |                         |     |
2.6 Implicit,fully-centered
| 2.4 ImEx,pressure-centeredwithdp=1/2 |     |     | 2.10 |     |     |
| ------------------------------------ | --- | --- | ---- | --- | --- |
ImEx,fully-centered
| 2.2       |           |           | 2.05      |                                    |           |
| --------- | --------- | --------- | --------- | ---------------------------------- | --------- |
| 2.0       |           |           | 2.00      |                                    |           |
| 0.00 0.02 | 0.04 0.06 | 0.08 0.10 | 0.00 0.02 | 0.04 0.06                          | 0.08 0.10 |
|           | t         |           |           | t                                  |           |
| 30        |           |           | 6         | Explicit,fully-upwindwithdp=du=1/2 |           |
Explicit,pressure-centeredwithdp=1
| 25  | Explicit,fully-upwindwithdp=du=1/2 |     | 5   |     |     |
| --- | ---------------------------------- | --- | --- | --- | --- |
Implicit,fully-upwindwithdp=du=1/2
Explicit,pressure-centeredwithdp=1
| 20     |                                    |     | 4    | Implicit,pressure-centeredwithdp=1 |     |
| ------ | ---------------------------------- | --- | ---- | ---------------------------------- | --- |
| )−hC   | Implicit,fully-upwindwithdp=du=1/2 |     | )−hC |                                    |     |
| b      |                                    |     | b    | Implicit,fully-centered            |     |
| (VT 15 | Implicit,pressure-centeredwithdp=1 |     | (VT  |                                    |     |
3
Implicit,fully-centered
10
|     | ImEx,pressure-centeredwithdp=1/2 |     | 2   |     |     |
| --- | -------------------------------- | --- | --- | --- | --- |
ImEx,fully-centered
5
1
0
| 0.00 0.02 | 0.04 0.06 | 0.08 0.10 | 0.00 0.02 | 0.04 0.06 | 0.08 0.10 |
| --------- | --------- | --------- | --------- | --------- | --------- |
|           | t         |           |           | t         |           |
| 40        |           |           | 9.75      |           |           |
Explicit,fully-upwindwithdp=du=1/2
| 35  |     |     | 9.50 | Explicit,pressure-centeredwithdp=1 |     |
| --- | --- | --- | ---- | ---------------------------------- | --- |
Explicit,fully-upwindwithdp=du=1/2
Implicit,fully-upwindwithdp=du=1/2
| 30  | Explicit,pressure-centeredwithdp=1 |     | 9.25 |     |     |
| --- | ---------------------------------- | --- | ---- | --- | --- |
Implicit,pressure-centeredwithdp=1
| )+      | Implicit,fully-upwindwithdp=du=1/2 |     | )+ 9.00 |                         |     |
| ------- | ---------------------------------- | --- | ------- | ----------------------- | --- |
| hC 25 b |                                    |     | hC b    | Implicit,fully-centered |     |
| (VT     | Implicit,pressure-centeredwithdp=1 |     | (VT     |                         |     |
8.75
| 20  | Implicit,fully-centered          |     |      |     |     |
| --- | -------------------------------- | --- | ---- | --- | --- |
|     | ImEx,pressure-centeredwithdp=1/2 |     | 8.50 |     |     |
15
ImEx,fully-centered
8.25
10
8.00
| 0.00 0.02 | 0.04 0.06 | 0.08 0.10 | 0.00 0.02 | 0.04 0.06 | 0.08 0.10 |
| --------- | --------- | --------- | --------- | --------- | --------- |
|           | t         |           |           | t         |           |
Figure 3.5: Riemann problem : l ∞ -norm and total variation of the characteristic variables Ch−
+
and with respect to time obtained with staggered schemes. On the left, all the schemes are
Ch
plotted, while on the right only the explicit and implicit schemes are plotted. (cid:98)
(cid:98)

CHAPTER 3. STUDY OF THE STABILITY OF THE STAGGERED SCHEMES FOR
| THE ONE | DIMENSIONAL | WAVE SYSTEM |     | 77  |
| ------- | ----------- | ----------- | --- | --- |
1.00
0.98
0.96
))t(
h
U(E 0.94
Explicit,fully-upwindwithdp=du=1/2
|     | 0.92 Explicit,pressure-centeredwithdp=1 |     |     |     |
| --- | --------------------------------------- | --- | --- | --- |
Implicit,fully-upwindwithdp=du=1/2
Implicit,pressure-centeredwithdp=1
|     | 0.90 Implicit,fully-centered |     |     |     |
| --- | ---------------------------- | --- | --- | --- |
ImEx,pressure-centeredwithdp=1/2
|     | 0.88 ImEx,fully-centered |                |           |     |
| --- | ------------------------ | -------------- | --------- | --- |
|     | 0.00                     | 0.02 0.04 0.06 | 0.08 0.10 |     |
t
Figure 3.6: Riemann problem : energy of the numerical solution with respect to time.

Chapter 4
Hodge-Helmholtz decomposition
and de Rham complexes: continuous
and discrete aspects
Contents
4.1 Introduction . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 80
4.2 Continuous de Rham complexes and Hodge-Helmholtz decompos-
ition . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 80
4.2.1 A discussion on the case of the 2D de Rham complex. . . . . . . . . . 80
4.2.2 The de Rham formalism and harmonic forms . . . . . . . . . . . . . . 83
4.3 Discrete de Rham complexes . . . . . . . . . . . . . . . . . . . . . . 85
4.4 The N´ed´elec-Raviart-Thomas de Rham staggered approximation
space . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 87
4.4.1 Mesh and notation . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 87
4.4.2 The Raviart-Thomas space . . . . . . . . . . . . . . . . . . . . . . . . 87
4.4.3 Properties of the discrete complex . . . . . . . . . . . . . . . . . . . . 91
4.5 Conclusion . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 102
79

| CHAPTER    | 4.  | HODGE-HELMHOLTZ |     |     | DECOMPOSITION |     |         |     | AND DE RHAM |     |
| ---------- | --- | --------------- | --- | --- | ------------- | --- | ------- | --- | ----------- | --- |
| COMPLEXES: |     | CONTINUOUS      |     | AND | DISCRETE      |     | ASPECTS |     |             | 80  |
4.1 Introduction
In chapter 2 we have shown that the low Mach number limit of the Euler equations is formally
equivalent to the long time limit of a particular wave system. This limit is identified with a
particular Hodge-Helmholtz decomposition, on which boundary conditions are imposed. We
postulate now that these kinds of decompositions are natural byproducts of the de Rham
complexes formalism; in this chapter we aim at introducing gently this concept and ultimately
motivate its use for the discretization of problems which require structure preservation. To do
| so, this | chapter | is separated |     | as follow: |     |     |     |     |     |     |
| -------- | ------- | ------------ | --- | ---------- | --- | --- | --- | --- | --- | --- |
1) In section 4.2, a particular and simplified de Rham complex is presented in order to show
its core properties. Then these fundamental properties are shown to naturally result in
a Hodge-Helmholtz decomposition. Once this particular example is understood, a more
| generic | definition |     | can | be considered. |     |     |     |     |     |     |
| ------- | ---------- | --- | --- | -------------- | --- | --- | --- | --- | --- | --- |
2) Then, in section 4.3, the principles of discrete de Rham complexes are introduced, and
in particular we briefly touch on a necessary condition that makes a discrete de Rham
| complex |     | relevant. |     |     |     |     |     |     |     |     |
| ------- | --- | --------- | --- | --- | --- | --- | --- | --- | --- | --- |
3) section 4.4 will unravel thoroughly our choice of the discrete de Rham complex and detail
| its | useful | byproducts. |     |     |     |     |     |     |     |     |
| --- | ------ | ----------- | --- | --- | --- | --- | --- | --- | --- | --- |
4) Finally we will summarize the main ideas of the chapter in section 4.5.
| 4.2 | Continuous |     | de  | Rham | complexes |     |     | and | Hodge-Helmholtz | de- |
| --- | ---------- | --- | --- | ---- | --------- | --- | --- | --- | --------------- | --- |
composition
| 4.2.1 | A discussion |     | on  | the case | of  | the | 2D de | Rham | complex |     |
| ----- | ------------ | --- | --- | -------- | --- | --- | ----- | ---- | ------- | --- |
In this section we introduce in a semi-formal discussion the concept of de Rham complexes in
simple terms in 2 space dimensions. The aim of the following is not to provide rigorous proofs
but simply to show how the definition of complexes naturally lead to byproducts such as the
| Hodge-Helmholtz |     | decomposition. |     |     |     |     |     |     |     |     |
| --------------- | --- | -------------- | --- | --- | --- | --- | --- | --- | --- | --- |
Now, let Ω be an open set of R2, H1(Ω), H(div;Ω), the Sobolev spaces for which the
gradient (in the sense of distributions), respectively, the divergence, is L2 integrable and
−
∇ = ( ∂ ,∂ )t. One of the canonical de Rham complexes in 2 dimensions reads [83]:
| ⊥   | y   | x   |     |     |     |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
−
∇⊥
|     |     |     | H1(Ω)(Ω) |     |     | H(div;Ω) |     | div | L2(Ω). | (4.1) |
| --- | --- | --- | -------- | --- | --- | -------- | --- | --- | ------ | ----- |
−
|     |     |     |     |     | −−→ |     |     | −−−→ |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | ---- | --- | --- |
The arrow characterizes the fundamental property of a complex: an inclusion of the range of
| the first | operator | in the | kernel | of the | second | operator |     |     |     |     |
| --------- | -------- | ------ | ------ | ------ | ------ | -------- | --- | --- | --- | --- |
H1(Ω)
|     |     |     |     | ∇ ⊥     |     | ker(    | div | H(div;Ω) | ),  | (4.2) |
| --- | --- | --- | --- | ------- | --- | ------- | --- | -------- | --- | ----- |
|     |     |     |     |         |     | ⊂       | −   | |        |     |       |
|     |     |     |     | (cid:2) |     | (cid:3) |     |          |     |       |

| CHAPTER    |     | 4. HODGE-HELMHOLTZ |     |     |     | DECOMPOSITION |     |         | AND DE | RHAM |     |     |
| ---------- | --- | ------------------ | --- | --- | --- | ------------- | --- | ------- | ------ | ---- | --- | --- |
| COMPLEXES: |     | CONTINUOUS         |     |     | AND | DISCRETE      |     | ASPECTS |        |      |     | 81  |
which is no more than a formalization of the well-known calculus formula
|     |     |     |     |     |     | div(∇ |     | 0.  |     |     |     |     |
| --- | --- | --- | --- | --- | --- | ----- | --- | --- | --- | --- | --- | --- |
⊥ ) =
| This inclusion |     | property | enables | to  | define | the | following | space |     |     |     |     |
| -------------- | --- | -------- | ------- | --- | ------ | --- | --------- | ----- | --- | --- | --- | --- |
H1(Ω)
|     |     |     |     | ker( | div |            | )/∇ |         | ,       |     |     |     |
| --- | --- | --- | --- | ---- | --- | ---------- | --- | ------- | ------- | --- | --- | --- |
|     |     |     |     |      | −   | | H(div;Ω) |     | ⊥       |         |     |     |     |
|     |     |     |     |      |     |            |     | (cid:2) | (cid:3) |     |     |     |
which, we will call the space of harmonic forms. It measures the gap between divergence-free
fields and fields that are rotated gradient of functions. In fact, we will say that the complex
is exact [84] if this gap is reduced to 0 . In this case, any divergence-free field will actually
|        |      |             |          |         |               | { } |     |         |         |          |     |       |
| ------ | ---- | ----------- | -------- | ------- | ------------- | --- | --- | ------- | ------- | -------- | --- | ----- |
| derive | from | the rotated | gradient |         | of a function |     | ;   |         |         |          |     |       |
|        |      | exactness   | of the   | complex |               |     | ∇   | H1(Ω)   | = ker(  | div      | ).  | (4.3) |
|        |      |             |          |         |               |     | ⊥   |         |         | H(div;Ω) |     |       |
|        |      |             |          |         |               | ←→  |     |         | −       | |        |     |       |
|        |      |             |          |         |               |     |     | (cid:2) | (cid:3) |          |     |       |
The exactness of the complex is strongly linked to the topological properties of the underlying
domain on which it is defined, we will elaborate on this point later. In parallel, it should be
noted that the differential operator ∇ can be considered as a linear operator defined on a dense
subspace of L2(Ω): H1(Ω). This makes it a densely defined unbounded operator on L2(Ω) with
H1(Ω):
domain ’unbounded’ because it is defined on a subspace and ’densely defined’ because
thissubspaceisdenseinL2(Ω).
|     |     |     |     | Inthesameway, |     |     |     | div canalsobeconsideredasasuchdensely |     |     |     |     |
| --- | --- | --- | --- | ------------- | --- | --- | --- | ------------------------------------- | --- | --- | --- | --- |
−
defined unbounded operator onL2(Ω;R2)withdomainH(div;Ω). Hence, followingthisremark,
since div is a linear operator defined on a subspace of L2(Ω;R2) and ker( div ) is
H(div;Ω)
| −   |     |     |     |     |     |     |     |     |     |     | − | |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
closed in L2(Ω;R2), we have by the classical orthogonality theorem in Hilbert spaces, that
L2
L2(Ω;R2)
|     |     |     |     | = ker( |     | div H(div;Ω) | )⊥  | ker( | div H(div;Ω) | ) ⊥ , |     | (4.4) |
| --- | --- | --- | --- | ------ | --- | ------------ | --- | ---- | ------------ | ----- | --- | ----- |
|     |     |     |     |        | −   | |            |     | ⊕    | − |          |       |     |       |
L2
| where | the orthogonality |     | is  | meant | in the |     | scalar | product | sense |     |     |     |
| ----- | ----------------- | --- | --- | ----- | ------ | --- | ------ | ------- | ----- | --- | --- | --- |
−
|     |     |     |     |     | v,w       | L2(Ω)2    | :=  | v wdx. |     |     |     |     |
| --- | --- | --- | --- | --- | --------- | --------- | --- | ------ | --- | --- | --- | --- |
|     |     |     |     |     | (cid:104) | (cid:105) |     | ·      |     |     |     |     |
(cid:90)Ω
In [84, Chapter 3] it is shown that the adjoint of div is ∇ using the theory of
|     |     |     |     |     |     |     | −   | | H(div;Ω) |     | | H 1(Ω) |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | ---------- | --- | -------- | --- | --- |
0
densely defined unbounded closed operators (here, closed means with closed graph). For the
sake of clarity, we do not go through this step, however we can quickly convince ourselves of
this property by remarking that if v H(div;Ω) and ϕ H1(Ω), using integration by parts,
0
|                  |     |      |             |     | ∈    |       |        |     | ∈    |     |     |     |
| ---------------- | --- | ---- | ----------- | --- | ---- | ----- | ------ | --- | ---- | --- | --- | --- |
| [84, Proposition |     | 3.11 | and Theorem |     | 3.12 | p 28] | yields |     |      |     |     |     |
|                  |     |      | div(v)ϕdx   |     | =    | v     | ∇ϕdx   | v   | n, ϕ | ,   |     |     |
(cid:105)H1/2(∂Ω)
|     |     |     | −         |     |     |           | ·   | −(cid:104) | ·   |     |     |     |
| --- | --- | --- | --------- | --- | --- | --------- | --- | ---------- | --- | --- | --- | --- |
|     |     |     | (cid:90)Ω |     |     | (cid:90)Ω |     |            |     |     |     |     |
=0
| or equivalently |     |            |     |     |       |     |            |     | (cid:124)(cid:123)(cid:122)(cid:125) |           |     |       |
| --------------- | --- | ---------- | --- | --- | ----- | --- | ---------- | --- | ------------------------------------ | --------- | --- | ----- |
|                 |     | v H(div;Ω) |     | ϕ   | H1(Ω) |     | div(v),    | ϕ   | = v,                                 | ∇ϕ        | ,   | (4.5) |
|                 |     |            |     |     |       | 0   |            |     | L2(Ω)                                | L2(Ω)2    |     |       |
|                 |     | ∀ ∈        |     | ∀   | ∈     |     | (cid:104)− |     | (cid:105) (cid:104)                  | (cid:105) |     |       |

CHAPTER 4. HODGE-HELMHOLTZ DECOMPOSITION AND DE RHAM
COMPLEXES: CONTINUOUS AND DISCRETE ASPECTS 82
which is the definition of an adjoint in Hilbert spaces. Moreover if v ker( div ) then
H(div;Ω)
∈ − |
(4.5) shows that
ker( div ) = ∇ H1(Ω) . (4.6)
− | H(div;Ω) ⊥ 0
Inserting (4.6) in (4.4) yields (cid:2) (cid:3)
L2
L2(Ω;R2) = ∇ H1(Ω) ⊥ ker( div ). (4.7)
0 ⊕ − | H(div;Ω)
(cid:2) (cid:3)
We observe, by the inclusion property (4.2), that ∇ [H1(Ω)] is closed linear subspace of
⊥
ker( div ), so invoking again the orthogonality theorem in Hilbert spaces yields
H(div;Ω)
− |
L2
ker( div ) = ∇ [H1(Ω)]⊥ ∇ [H1(Ω)] ⊥. (4.8)
H(div;Ω) ⊥ ⊥
− | ⊕
(cid:16) (cid:17)
It is possible to show thanks to a Poincar´e inequality ([85, discussion p 22]) that the range of
the rotated gradient is closed; yielding from (4.8) :
L2
ker( div ) = ∇ H1(Ω) ⊥ ∇ H1(Ω) ⊥. (4.9)
H(div;Ω) ⊥ ⊥
− | ⊕
(cid:16) (cid:17)
(cid:2) (cid:3) (cid:2) (cid:3)
Using the symbol to denote the relation between spaces that are isomorphic, we have
≡
∇ ⊥ H1(Ω) ⊥ ker( div H(div;Ω) )/∇ ⊥ H1(Ω) . (4.10)
≡ − |
(cid:2) (cid:3) (cid:2) (cid:3)
Plugging (4.10) in (4.9) yields
L2
ker( div ) = ∇ H1(Ω) ⊥ ker( div )/∇ H1(Ω) . (4.11)
H(div;Ω) ⊥ H(div;Ω) ⊥
− | ⊕ − |
(cid:20) (cid:21)
(cid:2) (cid:3) (cid:2) (cid:3)
Gathering (4.11) in (4.7) gives the general Hodge-Helmholtz decomposition :
Theorem 4.2.1 (Hodge-Helmholtz decomposition with complexes). Let Ω R2 a bounded
⊂
open set. We have the following Hodge-Helmholtz decomposition
L2
L2(Ω,R2) = ∇ H1(Ω) ⊥ L ,
0 ⊕ Ψ
where (cid:2) (cid:3)
L2
L := ∇ H1(Ω) ⊥ ker( div )/∇ H1(Ω) ,
Ψ ⊥ H(div;Ω) ⊥
⊕ − |
(cid:20) (cid:21)
(cid:2) (cid:3) (cid:2) (cid:3)
and
div[L ] = 0 .
Ψ
{ }
Now thatwe have explainedtheformalism ofcomplexes andits practical consequencessuch
as the HHD, we must address the space of harmonic forms. Indeed, throughout chapter 2, we
have mainly insisted on decompositions of the form
L2
V = ∇[H]⊥ L div[L ] = 0 ,
Ψ Ψ
⊕ − { }

| CHAPTER    |     | 4. HODGE-HELMHOLTZ |     |     |     | DECOMPOSITION |         |     | AND | DE RHAM |     |
| ---------- | --- | ------------------ | --- | --- | --- | ------------- | ------- | --- | --- | ------- | --- |
| COMPLEXES: |     | CONTINUOUS         |     |     | AND | DISCRETE      | ASPECTS |     |     |         | 83  |
asin(4.7)whichmightleadtothinkthattheevocationofharmonic forms isuseless, sincethey
are included in the divergence kernel (see Theorem 4.2.1). However it is these very harmonic
forms that ensures the relevance of a de Rham complex at the discrete scale.
| 4.2.2 | The | de  | Rham | formalism |     | and | harmonic |     | forms |     |     |
| ----- | --- | --- | ---- | --------- | --- | --- | -------- | --- | ----- | --- | --- |
So now, with this example in mind, we introduce more rigorously the complex and its funda-
| mental | concept, | the | harmonic | forms. |     |     |     |     |     |     |     |
| ------ | -------- | --- | -------- | ------ | --- | --- | --- | --- | --- | --- | --- |
Definition 4.2.1 ((Hilbert) Complex, see [84, Chapter 2, p 16] ). A (Hilbert) complex is
(Vk,dk)
a sequence of couples (closed densely defined) linear operators/(Hilbert) spaces k M ,
| M N | such | that | dk : Vk | Vk+1 | and |     |            |     |     |     | ∈   |
| --- | ---- | ---- | ------- | ---- | --- | --- | ---------- | --- | --- | --- | --- |
| ⊂   |      |      |         | −→   |     |     |            |     |     |     |     |
|     |      |      |         |      | dk  | Vk  | ker(dk+1). |     |     |     |     |
(4.12)
⊂
(cid:104) (cid:105)
| It is said | to  | be exact | iff |     |     |     |             |     |     |     |     |
| ---------- | --- | -------- | --- | --- | --- | --- | ----------- | --- | --- | --- | --- |
|            |     |          |     |     | dk  | Vk  | = ker(dk+1) | .   |     |     |     |
schem(cid:2)atic(cid:3)form
| It can | be put | under | the following |     |      |      |     |      |      |     |     |
| ------ | ------ | ----- | ------------- | --- | ---- | ---- | --- | ---- | ---- | --- | --- |
|        |        |       | dk−2          |     |      | dk−1 |     | dk   | dk+1 |     |     |
|        |        |       | ...           |     | Vk 1 |      | Vk  | Vk+1 |      | ... |     |
−
|     |     |     | −−−→ |     |     | −−−→ |     | −→  | −−−→ |     |     |
| --- | --- | --- | ---- | --- | --- | ---- | --- | --- | ---- | --- | --- |
The interest of introducing this definition appears when defining systematically the spaces
k(Ω)
of k harmonic forms which measures, as previously, the default of exactness of the
| −   |     |     | H   |     |     |     |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
complex
Definition 4.2.2 (Harmonic forms, see [84, Chapter 4 p 35] ). Given a complex (Vk,dk) ,
k M
∈
| the space | of  | k harmonic |     | forms | are given | by  |     |     |     |     |     |
| --------- | --- | ---------- | --- | ----- | --------- | --- | --- | --- | --- | --- | --- |
−
|     |     |     |     |     | k(Ω) | ker(dk+1)/dk |     | Vk  | .   |     |     |
| --- | --- | --- | --- | --- | ---- | ------------ | --- | --- | --- | --- | --- |
:=
H
|     |     |     |     |     |     |     |     | (cid:104) | (cid:105) |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --------- | --------- | --- | --- |
By convention the first space of the sequence (Vk) is indexed by the superscript 0 and its
k M
∈
| associated | space | of  | harmonic | forms | is  | defined | as         |     |     |     |     |
| ---------- | ----- | --- | -------- | ----- | --- | ------- | ---------- | --- | --- | --- | --- |
|            |       |     |          |       |     | 0(Ω)    | := ker(d0) |     |     |     |     |
H
| Then | in  | the case | (4.1), |     |     |     |     |     |     |     |     |
| ---- | --- | -------- | ------ | --- | --- | --- | --- | --- | --- | --- | --- |
∇⊥
|     |     |     |     | H1(Ω) |     | H(div;Ω) |     | div | L2(Ω), |     |     |
| --- | --- | --- | --- | ----- | --- | -------- | --- | --- | ------ | --- | --- |
−
|            |     |          |       |      | −−→ |     |     | −−−→ |     |     |     |
| ---------- | --- | -------- | ----- | ---- | --- | --- | --- | ---- | --- | --- | --- |
| the spaces | of  | harmonic | forms | are: |     |     |     |      |     |     |     |
• 0(Ω)
|     |      | := ker(∇ | ⊥ H1(Ω) | ),       |     |         |         |     |     |     |     |
| --- | ---- | -------- | ------- | -------- | --- | ------- | ------- | --- | --- | --- | --- |
| H   |      |          | |       |          |     |         |         |     |     |     |     |
| •   | 1(Ω) | := ker(  | div     |          | )/∇ | H1(Ω)   | ,       |     |     |     |     |
|     |      |          |         | H(div;Ω) |     | ⊥       |         |     |     |     |     |
| H   |      |          | −       | |        |     |         |         |     |     |     |     |
|     |      |          |         |          |     | (cid:2) | (cid:3) |     |     |     |     |

| CHAPTER    | 4.      | HODGE-HELMHOLTZ |                 |     |          | DECOMPOSITION |         |     | AND | DE RHAM |     |
| ---------- | ------- | --------------- | --------------- | --- | -------- | ------------- | ------- | --- | --- | ------- | --- |
| COMPLEXES: |         | CONTINUOUS      |                 | AND | DISCRETE |               | ASPECTS |     |     |         | 84  |
| •          | 2(Ω) := | L2(Ω)/(         | div[H(div;Ω)]). |     |          |               |         |     |     |         |     |
| H          |         |                 | −               |     |          |               |         |     |     |         |     |
The depth and fundamentality of this notion is finally understood when we note that the space
of k harmonic forms are related to topological properties of the domain: in the case of the
−
k(Ω))
complex (4.1), the spaces ( k are finite dimensional and the dimension of each of them
H
|     |     |     |     |     |     |     |     | Betti’s | numbers | (b  |     |
| --- | --- | --- | --- | --- | --- | --- | --- | ------- | ------- | --- | --- |
is linked to topological features of the domain through i ) 0 i 2 [86]. In the
≤≤
case of the sequence (4.1), a simple example shows that some of these links are very natural: if
a H1 function is gradient free, then it is constant on each connected component of its domain
| of definition, | so  | that | straightforwardly |     |     |     |     |     |     |     |     |
| -------------- | --- | ---- | ----------------- | --- | --- | --- | --- | --- | --- | --- | --- |
0(Ω)
dim := b where b is the number of connected components of the domain.
|     | H   |     | 0   | 0   |     |     |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
Nonetheless, it is not as evident to get a feeling of the link between the dimension of 1(Ω)
H
| and the | topology | of   | the domain | since,  | in   | fact       |     |          |         |             |     |
| ------- | -------- | ---- | ---------- | ------- | ---- | ---------- | --- | -------- | ------- | ----------- | --- |
|         | dim      | 1(Ω) | := b       | , where | b is | the number |     | of holes | through | the domain, |     |
|         |          | H    | 1          |         | 1    |            |     |          |         |             |     |
and finally
dim 2(Ω) := b , where b is 0 for any bounded domain in R2 and 1 for the torus.
|     |     |     | 2   | 2   |     |     |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
H
Remark 4.2.2. dim 2(Ω) = 0 for bounded domains leads in particular to the surjectivity of
H
| the operator |     | div | :   |     |               |     |     |        |     |     |     |
| ------------ | --- | --- | --- | --- | ------------- | --- | --- | ------ | --- | --- | --- |
|              | −   |     |     |     |               |     |     | L2(Ω), |     |     |     |
|              |     |     |     |     | div[H(div;Ω)] |     | =   |        |     |     |     |
−
| which is | equivalent |     | to the continuous |     | inf-sup | condition |     | [87]. |     |     |     |
| -------- | ---------- | --- | ----------------- | --- | ------- | --------- | --- | ----- | --- | --- | --- |
Remark 4.2.3. There exists another de Rham complex in two space dimensions, similar to
| (4.1) which | reads | [83] |     |       |     |           |     |      |       |     |        |
| ----------- | ----- | ---- | --- | ----- | --- | --------- | --- | ---- | ----- | --- | ------ |
|             |       |      |     | H1(Ω) | ∇   | H(curl;Ω) |     | curl | L2(Ω) |     | (4.13) |
|             |       |      |     |       | −→  |           |     | −−→  |       |     |        |
where curlu := ∂ uy ∂ ux and H(curl;Ω) the Sobolev space for which the distributional curl
|     |     | x   | y   |     |     |     |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
−
is L2(Ω). Equivalently, if Ω is an open set of R3, the complex (4.1) becomes
|     |     | H1(Ω) | ∇   |          |     | rot |          |     | div  | L2(Ω), |     |
| --- | --- | ----- | --- | -------- | --- | --- | -------- | --- | ---- | ------ | --- |
|     |     |       |     | H(rot;Ω) |     |     | H(div;Ω) |     | −    |        |     |
|     |     |       | −→  |          |     | −−→ |          |     | −−−→ |        |     |
where
|     |     |     |     |     |     |     | uz  | uy  |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
|     |     |     |     |     |     | ∂   | y   | ∂ z |     |     |     |
−
|     |     |     |     | rotu |     | ∂   | ux  | ∂ uz |     |     |     |
| --- | --- | --- | --- | ---- | --- | --- | --- | ---- | --- | --- | --- |
|     |     |     |     |      | :=  |     | z   | x    |     |     |     |
|     |     |     |     |      |     |    | −   |      |    |     |     |
|     |     |     |     |      |     | ∂   | uy  | ∂ ux |     |     |     |
|     |     |     |     |      |     |     | x   | y    |     |     |     |
−
|     |     |     |     |     |     |    |     |     |    |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
L2.
and H(rot;Ω) the Sobolev space for which the distributional rot is Its harmonic forms are
| also related | to  | the topology | of  | the | underlying | domain |     | of definition | Ω:  |     |     |
| ------------ | --- | ------------ | --- | --- | ---------- | ------ | --- | ------------- | --- | --- | --- |
dim 0(Ω) := b is the number of connected components of the domain,
0
H

| CHAPTER    | 4.  | HODGE-HELMHOLTZ |     |     | DECOMPOSITION |     |         |     | AND | DE RHAM |     |
| ---------- | --- | --------------- | --- | --- | ------------- | --- | ------- | --- | --- | ------- | --- |
| COMPLEXES: |     | CONTINUOUS      |     | AND | DISCRETE      |     | ASPECTS |     |     |         | 85  |
then
|     | dim | 1(Ω) | :=  | b is | the number | of  | of tunnels | through |     | the domain, |     |
| --- | --- | ---- | --- | ---- | ---------- | --- | ---------- | ------- | --- | ----------- | --- |
1
H
and
|     |     | dim 2(Ω) | :=  | b is | the number | of  | voids | enclosed | by  | the domain. |     |
| --- | --- | -------- | --- | ---- | ---------- | --- | ----- | -------- | --- | ----------- | --- |
2
H
Finally
|     |     | 3(Ω) | b    | is  | for bounded | domains |     | in R3 | and 1 | on the torus. |     |
| --- | --- | ---- | ---- | --- | ----------- | ------- | --- | ----- | ----- | ------------- | --- |
|     | dim |      | := 3 | 0   |             |         |     |       |       |               |     |
H
Remark 4.2.4. Note that the complex (4.1) is usually given under the following enriched form
|     |     |     | id    |     | ∇⊥       |     |      | div |       | 0      |     |
| --- | --- | --- | ----- | --- | -------- | --- | ---- | --- | ----- | ------ | --- |
|     |     | R   | H1(Ω) |     | H(div;Ω) |     | −    |     | L2(Ω) | 0 ,    |     |
|     |     |     | −→    |     | −−→      |     | −−−→ |     |       | →− { } |     |
whichyieldsan exactcomplexinthecaseofaconnecteddomain, asopposedto (4.1). Indeedthis
)/R,
complex yields, by definition, a space of 0 harmonic forms equal to ker(∇ which
⊥ H1(Ω)
|     |     |     |     |     | −   |     |     |     |     | |   |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
simplify to 0 in the case of a connected domain. This choice of presentation is however only
{ }
| made for     | the sake | of clarity |      | and does | not affect | the | following. |     |     |     |     |
| ------------ | -------- | ---------- | ---- | -------- | ---------- | --- | ---------- | --- | --- | --- | --- |
| 4.3 Discrete |          | de         | Rham |          | complexes  |     |            |     |     |     |     |
The formalism of complexes encapsulates a lot of fundamental and natural mathematical
properties that arise at the continuous level and replicating in simulations some of these
byproducts is essential. Indeed, in a lot of applications, it is necessary to preserve involution
constraints [88] linked to the inclusion property (4.12), differential constraints (as for Maxwell
or incompressible Navier-Stokes equations) or even replicate a discrete inf-sup conditions [87].
These observations are not new, and in fact, embracing this philosophy has been very fertile
in the development of approximation spaces and new numerical methods that are compatible
with these continuous differential constraints: one can cite mimetic finite difference methods
[89,90,91,92],compatiblediscretizations[93,94,95],theDiscretedeRhamframework[96,97],
finiteelementexteriorcalculus[85,84]oralsointhediscontinuousGalerkinframework[83,98].
The main difficulty, and not the least, when designing a discrete de Rham complex of type
(4.1), resides in showing that the discrete complex ensures an harmonic forms isomorphism and
that these discrete spaces are satisfying approximations of the continuous harmonic forms: the
discrete harmonic forms k (Ω) should be isomorphic to the continuous ones k(Ω),
M
|     |     |     | H   |     |       |     |       |     |     | H   |     |
| --- | --- | --- | --- | --- | ----- | --- | ----- | --- | --- | --- | --- |
|     |     |     |     |     | k     |     | k(Ω). |     |     |     |     |
|     |     |     |     |     | M (Ω) |     |       |     |     |     |     |
|     |     |     |     |     | H     | ≡   | H     |     |     |     |     |

| CHAPTER    | 4.  | HODGE-HELMHOLTZ |     |     |          | DECOMPOSITION |         | AND | DE RHAM |     |     |
| ---------- | --- | --------------- | --- | --- | -------- | ------------- | ------- | --- | ------- | --- | --- |
| COMPLEXES: |     | CONTINUOUS      |     | AND | DISCRETE |               | ASPECTS |     |         |     | 86  |
Since these spaces depend in particular on topological features of the domain, typically, in two
and three dimensions, the number of holes and tunnels in the domain, it is not obvious how to
deal with such isomorphism in a systematic manner for a general domain.
We avoid this problem by basing our staggered scheme on a preexisting finite element complex,
theN´ed´elec-Raviart-Thomas complexwhichreads[84,99]onquadrangularandtriangular
| meshes, | Figure | 4.1: |     |        |     |        |     |        |     |     |     |
| ------- | ------ | ---- | --- | ------ | --- | ------ | --- | ------ | --- | --- | --- |
|         |        |      |     |        | ∇⊥  | RT1(Ω) | div |        |     |     |     |
|         |        |      |     | cG1(Ω) |     |        | −   | dG0(Ω) |     |     |     |
|         |        |      |     |        | −→  |        | −→  |        |     |     |     |
Figure 4.1: The H(div;Ω) N´ed´elec-Raviart-Thomas complex in 2 space dimensions
−
cG1(Ω)
In particular this discrete complex is included in the complex (4.1) in the sense that
⊂
H1(Ω) where cG1(Ω) is, on triangular meshes, the space of continuous polynomials of order
1 on each variable x,y and, on quadrangular meshes, the space of continuous polynomials of
|              |     |      |         |         |      | RT1(Ω) |          |     | dG0(Ω) | L2(Ω) |       |
| ------------ | --- | ---- | ------- | ------- | ---- | ------ | -------- | --- | ------ | ----- | ----- |
| total degree | at  | most | 1, then | we also | have |        | H(div;Ω) | and |        |       | where |
|              |     |      |         |         |      |        | ⊂        |     |        | ⊂     |       |
dG0(Ω) is the space of cellwise constant funtions. These spaces verify an important property of
commutationbetweenthedifferentialoperatorsandboundedinterpolationoperators(I ) :
|        |               |     |         |          |     |     |     |     |     |     | k 0 k 2 |
| ------ | ------------- | --- | ------- | -------- | --- | --- | --- | --- | --- | --- | ------- |
|        |               |     |         |          |     |     |     |     |     |     | ≤ ≤     |
| indeed | the following |     | diagram | commutes |     |     |     |     |     |     |         |
∇⊥
|     |     |     |     | H1(Ω) |     | H(div;Ω) | div | L2(Ω) |     |     |     |
| --- | --- | --- | --- | ----- | --- | -------- | --- | ----- | --- | --- | --- |
−
|     |     |     |     |                | −→  |                 | −→  |                |     |     |     |
| --- | --- | --- | --- | -------------- | --- | --------------- | --- | -------------- | --- | --- | --- |
|     |     |     |     | I              |     | I               |     | I              | .   |     |     |
|     |     |     |     | 0              |     |                 | 1   | 2              |     |     |     |
|     |     |     |     |               |     |                |     |               |     |     |     |
|     |     |     |     |               |     |                |     |               |     |     |     |
|     |     |     |     |               |     |                |     |               |     |     |     |
|     |     |     |     | (cid:121) 1(Ω) | ∇⊥  | RT(cid:121)1(Ω) | div | (cid:121) 0(Ω) |     |     |     |
|     |     |     |     | cG             |     |                 | −   | dG             |     |     |     |
|     |     |     |     |                | −→  |                 | −→  |                |     |     |     |
The interpolation operators (I k ) 0 k 2 can be constructed by regularizing the canonical inter-
| polation | operators | [100, | Section | 5.4], ≤ | ≤ [101]. |     |     |     |     |     |     |
| -------- | --------- | ----- | ------- | ------- | -------- | --- | --- | --- | --- | --- | --- |

| CHAPTER    |             | 4.      | HODGE-HELMHOLTZ          |          |             | DECOMPOSITION |       |         | AND         | DE        | RHAM |         |     |
| ---------- | ----------- | ------- | ------------------------ | -------- | ----------- | ------------- | ----- | ------- | ----------- | --------- | ---- | ------- | --- |
| COMPLEXES: |             |         | CONTINUOUS               |          | AND         | DISCRETE      |       | ASPECTS |             |           |      |         | 87  |
| 4.4        |             | The     | N´ed´elec-Raviart-Thomas |          |             |               |       | de      | Rham        | staggered |      | approx- |     |
|            |             | imation | space                    |          |             |               |       |         |             |           |      |         |     |
| Our        | proposition |         | follows                  | the line | of research |               | [66], | in the  | sense that: |           |      |         |     |
•
|     | we  | base | our approximation |     | on  | a finite | element | ground, |     |     |     |     |     |
| --- | --- | ---- | ----------------- | --- | --- | -------- | ------- | ------- | --- | --- | --- | --- | --- |
• the vectorial unknown is located, in some way we will define, on the faces of the mesh.
Bycontrast, weintroducethefollowingoriginality: ourstaggereddiscretizationisbasedonthe
de Rham complex of N´ed´elec-Raviart-Thomas given in Figure 4.1, where the spaces cG1(Ω)
which is H1 conforming, dG0(Ω) which is L2 conforming are, respectively, the continuous
|     |     | −   |     |     |     |     | −   |     |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
polynomials of order 1 and the space of cellwise constant functions. Finally, the main space
RT1(Ω),
of interest here is the lowest order Raviart-Thomas finite element space, which is the
smallest normal-conforming finite element space for which the divergence is cellwise constant
[102]. We detail in the following its definition but first we introduce the notation we will
| repeatedly |     | use. |              |     |     |     |     |     |     |     |     |     |     |
| ---------- | --- | ---- | ------------ | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 4.4.1      |     | Mesh | and notation |     |     |     |     |     |     |     |     |     |     |
R2,
Given a domain Ω is the mesh on this surface, will denote the set of cells of ,
|     |     |     | ⊂   | M   |     |     |     |     | C   |     |     |     | M   |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
K a cell of this mesh and ν(K) the number of neighbouring cells of K. These cells will be
∈ C
either triangular or quadrangular shaped (so ν(K) is either equal to 3 or 4). ν will denote
max
max[ν(K)]. σ ∂K will denote a face of a cell K and , int, b are respectively the set
| K   |     |     | ∈   |     |     |     |     |     | F F | F   |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
of∈aCll
faces of the mesh which are oriented, the set of interior faces and the set of boundary
faces. For a fixed cell K and a face σ ∂K, n is the outward unit normal to σ. For
K,σ
|     |     |     |     | ∈ C |     |     | ⊂   |     |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
a given interior face σ, n is a unit normal associated to this face, we denote L and K the
|     |     |     |     | σ   |     |     |     |     |     |     |     | σ   | σ   |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
cells such that L K = σ and L is the exterior cell with respect to the normal n . This set
|     |     |     | σ σ |     | σ   |     |     |     |     |     |     | σ   |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
∩
of neighbouring cells will be denoted by (σ) = L ,K . In the case where σ b then we
|     |     |     |     |     |     |     |     | σ   | σ   |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
|     |     |     |     |     |     | C   |     | {   | }   |     |     | ∈ F |     |
have (σ) = K h will denote min K and σ = max σ . Finally, we define the following
|          | C   | {   | σ } |     |       | | |  | | | max |      | | |   |      |       |     |        |
| -------- | --- | --- | --- | --- | ----- | ---- | ------- | ---- | ----- | ---- | ----- | --- | ------ |
|          |     |     |     |     | K     |      |         |      | σ     |      |       |     |        |
| notation |     |     |     |     | ∈C    |      |         |      | ∈F    |      |       |     |        |
|          |     |     | σ   | ,   | K     | (σ), | ε (σ)   | := n | n     |      | 1,1 . |     | (4.14) |
|          |     |     |     |     |       |      | K       |      | K,σ σ |      |       |     |        |
|          |     |     | ∀ ∈ | F   | ∀ ∈ C |      |         |      | ·     | ∈ {− | }     |     |        |
For the boundary faces, the associated normal n is oriented in such way that it matches the
σ
| outgoing |       | normal         |                |     |       |      |           |     |         |       |      |     |     |
| -------- | ----- | -------------- | -------------- | --- | ----- | ---- | --------- | --- | ------- | ----- | ---- | --- | --- |
|          |       |                | σ              | b,  | K     | (σ), | n         | =   | n (ε    | (σ) = | 1) . |     |     |
|          |       |                |                |     |       |      | K,σ       |     | σ K     |       |      |     |     |
|          |       |                | ∀              | ∈ F | ∀     | ∈ C  |           |     |         |       |      |     |     |
| 4.4.2    |       | The            | Raviart-Thomas |     | space |      |           |     |         |       |      |     |     |
| The      | local | Raviart-Thomas |                |     | space | on a | reference |     | element |       |      |     |     |
Following Ciarlet [103] a finite element is defined by a triplet (Kˆ,Vˆ,(θ ) ) where Kˆ is a
|     |     |     |     |     |     |     |     |     |     |     | i 1 i | v   |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | ----- | --- | --- |
|     |     |     |     |     |     | Rd  |     |     |     |     | ≤     | ≤   | Vˆ  |
compact, connected, Lipschitz subset of with non-empty interior, typica ll y a polygon, is
Rp
a vector space of functions, here of dimension v, generally polynomials with values in and

CHAPTER 4. HODGE-HELMHOLTZ DECOMPOSITION AND DE RHAM
COMPLEXES: CONTINUOUS AND DISCRETE ASPECTS 88
(cid:126)n (cid:126)n = n
L,σ σ Kσ
K L
σ σ
σ
Figure 4.2: Representation of (σ) the set of cells neighbouring a given face σ such that
C
ε (σ) = 1 and ε (σ) = 1.
Kσ Lσ
−
lastly, a family of v linear forms (θ ) acting on Vˆ . In particular, it is often required that
i 1 i v
− ≤≤
the family of linear forms is unisolvent:
p Vˆ, 1 i v, θ (p) = 0 = p = 0.
i
∀ ∈ ∀ ≤ ≤ ⇒
Thus, the Raviart-Thomas finite element is defined first on a reference element. The linear
forms are defined on the faces σˆ ∂Kˆ as
∈
1
vˆ Vˆ θ (vˆ) := vˆ nˆ dΓ (4.15)
∀ ∈
σˆ
σˆ ·
σˆ
| | (cid:90)
σˆ
where for a given face σˆ , σˆ denotes its (d 1)-Lebesgue measure. Then a basis is given on a
| | −
reference element, for example
• On the reference triangle (see Figure 4.4)
K := (xˆ,yˆ),0 xˆ 1, 0 yˆ 1 xˆ ,
{ ≤ ≤ ≤ ≤ − }
the basis functions reads(cid:98)
xˆ xˆ 1 xˆ
, − , . (4.16)
yˆ 1 yˆ yˆ
(cid:18) − (cid:19) (cid:18) (cid:19) (cid:18) (cid:19)
• On the reference quadrangle (see Figure 4.5)
K := [0,1]2,
the basis functions are (cid:98)
xˆ 0 xˆ 1 0
, , − , . (4.17)
0 yˆ 0 yˆ 1
(cid:18) (cid:19) (cid:18) (cid:19) (cid:18) (cid:19) (cid:18) − (cid:19)

| CHAPTER    | 4. HODGE-HELMHOLTZ |     |     | DECOMPOSITION |     |         | AND | DE RHAM |     |
| ---------- | ------------------ | --- | --- | ------------- | --- | ------- | --- | ------- | --- |
| COMPLEXES: | CONTINUOUS         |     | AND | DISCRETE      |     | ASPECTS |     |         | 89  |
They are built in such way that the linear forms (4.15) verify for any faces σˆ,fˆ ∂Kˆ and basis
∈
| function | Ψˆ (associated | to  | the face | σˆ) on the | reference | element: |     |     |     |
| -------- | -------------- | --- | -------- | ---------- | --------- | -------- | --- | --- | --- |
σˆ
|     |     |     |     | θ (Ψˆ | ) = | δ ,   |     |     | (4.18) |
| --- | --- | --- | --- | ----- | --- | ----- | --- | --- | ------ |
|     |     |     |     | fˆ    | σˆ  | fˆ,σˆ |     |     |        |
(Ψˆ
where δ is the Kronecker’symbol. With σˆ ) defined in (4.16), (4.17), the local space
| fˆ,σˆ |     |     |     |     | σˆ  | ∂Kˆ |     |     |     |
| ----- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
∈
| on the reference | element | is  | defined | as: |         |      |     |     |     |
| ---------------- | ------- | --- | ------- | --- | ------- | ---- | --- | --- | --- |
|                  |         |     | RT1(Kˆ) | :=  | span Ψˆ | , σˆ | ∂Kˆ |     |     |
σˆ
∈
|     |     |     |     |     | (cid:26) |     | (cid:27) |     |     |
| --- | --- | --- | --- | --- | -------- | --- | -------- | --- | --- |
ThePiolatransform: akeytooltodefinetheRaviart-Thomasonaphysicalelement
We then construct a local space on an arbitrary element K. The classical way in finite element
is to define the physical space on a cell K through a simple change of variable: if we denote
T : Kˆ K the transformation from the reference element to the physical element of
K
−→
interest (illustrated in Figure 4.3), then the local physical space V (K) is defined usually from
h
(Kˆ)
| the reference | space | V h | as  |     |     |     |     |     |     |
| ------------- | ----- | --- | --- | --- | --- | --- | --- | --- | --- |
(Kˆ)
|     |     | V   | (K) := | x        | pˆ(T 1(x)), | pˆ  | V   | .        |     |
| --- | --- | --- | ------ | -------- | ----------- | --- | --- | -------- | --- |
|     |     |     | h      |          | K−          |     | h   |          |     |
|     |     |     |        | →        |             | ∈   |     |          |     |
|     |     |     |        | (cid:26) |             |     |     | (cid:27) |     |
However, this procedure is not adapted to define the Raviart-Thomas finite element in K: in
particular it does not ensure that the basis functions on the physical element K defined by this
change of variable will verify (4.18) on K. Consequently, a particular transformation is needed
in order to preserve the normal component of the function on the boundary of the physical
element: the contravariant Piola transform. From a function vˆ RT1(Kˆ), it defines a function
∈
v as [104]:
1
|     |     | v(x) | = ( vˆ)(x) | :=  |           | B(x)vˆ(T |     | 1(x)), | (4.19) |
| --- | --- | ---- | ---------- | --- | --------- | -------- | --- | ------ | ------ |
|     |     |      |            |     | det(B(x)) |          | K−  |        |        |
P
|     |     |     |     |     | |   | |   |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
1(x)).
where B(x) := ∇ xˆ T K (T K− We can regroup the fundamental properties of this transform-
| ation in | the following | lemma: |     |     |     |     |     |     |     |
| -------- | ------------- | ------ | --- | --- | --- | --- | --- | --- | --- |
Lemma 4.4.1 (see [105, discussion p 2430 ] ). Let v = vˆ, vˆ RT1(Kˆ) a function defined
|     |     |     |     |     |     |     | P   | ∈   |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
by the contravariant Piola transformation (4.19). Then the contravariant Piola transform
∂Kˆ
| i) preserves | the | normal | component: | for | σ ∂K,σ | = T | K (σˆ) with | σˆ  |     |
| ------------ | --- | ------ | ---------- | --- | ------ | --- | ----------- | --- | --- |
|              |     |        |            |     | ∈      |     |             | ∈   |     |
dΓˆ.
|               |     |             |     | v           | n σ dΓ = | vˆ nˆ    | σˆ          |     |     |
| ------------- | --- | ----------- | --- | ----------- | -------- | -------- | ----------- | --- | --- |
|               |     |             |     | ·           |          | ·        |             |     |     |
|               |     |             |     | σ           |          | σˆ       |             |     |     |
|               |     |             |     | (cid:90)    |          | (cid:90) |             |     |     |
| ii) preserves | the | divergence: |     |             |          |          |             |     |     |
|               |     |             |     | div x (v)dx | =        | div      | xˆ (vˆ)dxˆ. |     |     |
|               |     |             |     | K           |          | Kˆ       |             |     |     |
|               |     |             |     | (cid:90)    |          | (cid:90) |             |     |     |

CHAPTER 4. HODGE-HELMHOLTZ DECOMPOSITION AND DE RHAM
COMPLEXES: CONTINUOUS AND DISCRETE ASPECTS 90
Then we define the Raviart-Thomas space on K as
RT1(K) := span Ψ σˆ ∂Kˆ ,
σˆ
P ∈
(cid:26) (cid:27)
(cid:98)
with (Ψ ) the basis of RT1(Kˆ).
σˆ σˆ ∂Kˆ
∈
(cid:98)
T
K
Kˆ K
Figure 4.3: Representation of the transformation T from a reference element Kˆ to a physical
K
element K.
The Global Raviart-Thomas space
Finally,theglobalspaceissimplydefinedbygluingthelocalspaceswithacontinuityconstraint
on the linear forms. In other words, we impose a normal-continuity type constraint on each
face:
RT1(Ω) := v L2(Ω)d, K ,v RT1(K), σ int, [[v n ]]dΓ = 0 ,
K σ
(cid:26) ∈ ∀ ∈ C | ∈ ∀ ∈ F (cid:90) σ · (cid:27)
where
[[v n ]] (x) := lim (v(x+tn ) n v(x tn ) n ).
σ σ σ σ σ σ
· t 0+ · − − ·
→
A basis of the global Raviart-Thomas space is given by the basis functions (Ψ ) that are
σ σ
∈F
such that
1
f,σ θ (Ψ ) := Ψ n dΓ = δ . (4.20)
f σ σ f σ,f
∀ ∈ F f ·
f
| | (cid:90)
Note that [[v n ]] = [[v ( n )]] so that it is independent of the orientation of the normal
σ σ σ σ
· · −
n . Finally, the Raviart-Thomas elements verify the following fundamental property:
σ
u RT1(Ω) div(u)dx = σ ε (σ)u . (4.21)
K σ
∀ ∈ | |
K
(cid:90) σ ∂K
(cid:88)∈
where we will denote u := u(x ) n , x is the center of mass of the face σ
σ σ σ σ
·

| CHAPTER    |     | 4. HODGE-HELMHOLTZ |     |     |          | DECOMPOSITION |     |         | AND DE RHAM |     |
| ---------- | --- | ------------------ | --- | --- | -------- | ------------- | --- | ------- | ----------- | --- |
| COMPLEXES: |     | CONTINUOUS         |     | AND | DISCRETE |               |     | ASPECTS |             | 91  |
div
−
RT1(Ω)
dG0(Ω)
p
|     |     |     |     | u   |     |     |     |     | h   |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
h
Figure4.4: Triangles: Representationofthefiniteelementapproximationspacesforthevelocity
| and the  | pressure   |       |        |          |     |         |     |     |     |     |
| -------- | ---------- | ----- | ------ | -------- | --- | ------- | --- | --- | --- | --- |
| 4.4.3    | Properties |       | of the | discrete |     | complex |     |     |     |     |
| Harmonic |            | forms |        |          |     |         |     |     |     |     |
k
Since it is a complex, we can define the spaces of discrete harmonic forms (Definition 4.2.2):
−
•
|     | 0 (Ω) | := ker(∇   |       | ),      |          |          |     |     |     |     |
| --- | ----- | ---------- | ----- | ------- | -------- | -------- | --- | --- | --- | --- |
|     | M     |            | ⊥ cG1 |         |          |          |     |     |     |     |
|     | H     |            | |     |         |          |          |     |     |     |     |
| •   | 1     |            |       |         | cG1      |          |     |     |     |     |
|     | (Ω)   | := ker(    | div   | RT1 )/∇ | ⊥        | ,        |     |     |     |     |
|     | H M   |            | −     | |       |          |          |     |     |     |     |
|     |       |            |       |         | (cid:20) | (cid:21) |     |     |     |     |
| •   |       |            |       |         | RT1      |          |     |     |     |     |
|     | 2 (Ω) | := dG0(Ω)/ |       | div     |          | .        |     |     |     |     |
M
|     | H   |     |     | −   |     |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
(cid:20) (cid:21)
These harmonic forms consist in fact of polynomial spaces that are good approximations of the
continuous harmonic forms. Moreover, they are isomorphic to the continuous harmonic spaces:
Theorem 4.4.1 (Harmonicformsisomorphism[84,Section7.6]). The discrete harmonic forms
k of the N´ed´elec-Raviart-Thomas complex are isomorphic to the harmonic forms of
| ( M (Ω)) | k   |     |     |     |     |     |     |     |     |     |
| -------- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
H
the continuous complex (4.1) ( k(Ω)) on both triangular and quadrangular meshes:
k
H
|     |     |     |     | k   | 0,1,2 |     | k   |     | k(Ω). |     |
| --- | --- | --- | --- | --- | ----- | --- | --- | --- | ----- | --- |
M (Ω)
|     |     |     |     | ∀   | ∈ { | }   | H   | ≡ H |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
Remark 4.4.2. In 3 space dimensions the N´ed´elec-Raviart-Thomas complex [99] reads
|     |     |     | cG1(Ω) |     | ∇ N1(Ω) |     | rot RT1(Ω) |     | div dG0(Ω), |     |
| --- | --- | --- | ------ | --- | ------- | --- | ---------- | --- | ----------- | --- |
−
|     |     |     |     | −→  |     | −→  |     | −→  |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
N1(Ω)
with the first order N´ed´elec finite element space [106]. Similar properties on harmonic
| forms | isomorphism |     | stands | in a domain |     | in three | space | dimensions. |     |     |
| ----- | ----------- | --- | ------ | ----------- | --- | -------- | ----- | ----------- | --- | --- |

| CHAPTER    | 4.  | HODGE-HELMHOLTZ |     |     | DECOMPOSITION |     |         | AND DE RHAM |     |
| ---------- | --- | --------------- | --- | --- | ------------- | --- | ------- | ----------- | --- |
| COMPLEXES: |     | CONTINUOUS      |     | AND | DISCRETE      |     | ASPECTS |             | 92  |
div
|     |     |     | RT1(Ω) |     |     | −   |     | dG0(Ω) |     |
| --- | --- | --- | ------ | --- | --- | --- | --- | ------ | --- |
p
h
u h
Figure 4.5: Quadrangles: Representation of the finite element approximation spaces for the
| velocity | and the      | pressure |           |     |     |     |     |     |     |
| -------- | ------------ | -------- | --------- | --- | --- | --- | --- | --- | --- |
| Discrete | differential |          | operators |     |     |     |     |     |     |
The discrete complex enables to define naturally discrete differential operators: we recall that
| the complex | of  | interest | in two | space  | dimensions |        | is        |     |     |
| ----------- | --- | -------- | ------ | ------ | ---------- | ------ | --------- | --- | --- |
|             |     |          |        |        | ∇⊥         |        | div       |     |     |
|             |     |          |        | cG1(Ω) |            | RT1(Ω) | − dG0(Ω). |     |     |
|             |     |          |        |        | −→         |        | −→        |     |     |
For reasons that will become clear in the next chapters, in the following cG1(Ω), dG0(Ω) will
| be equipped | with | the | usual L2 | scalar | product |     |       |     |     |
| ----------- | ---- | --- | -------- | ------ | ------- | --- | ----- | --- | --- |
|             |      |     |          |        | q, ϕ    | :=  | qϕdx, |     |     |
L2(Ω)
|     |     |     |     | (cid:104) | (cid:105) |     | (cid:90)Ω |     |     |
| --- | --- | --- | --- | --------- | --------- | --- | --------- | --- | --- |
whereas RT1(Ω) will be equipped with a discrete scalar product that approximates the natural
L2(Ω)2 scalarproduct. Wewilldenotethis, to-be-defined(inchapter5),discretescalarproduct
by
|     |     |     |     |     |     | ., .                | .   |     |     |
| --- | --- | --- | --- | --- | --- | ------------------- | --- | --- | --- |
|     |     |     |     |     |     | (cid:104) (cid:105) | h   |     |     |
Discrete adjoints: By definition of a complex, the diagram shows in particular that
|     |     |     |     |     | RT1(Ω) |     | dG0(Ω), |     |     |
| --- | --- | --- | --- | --- | ------ | --- | ------- | --- | --- |
div
|              |       |         |     | −       |         |     | ⊂           |     |     |
| ------------ | ----- | ------- | --- | ------- | ------- | --- | ----------- | --- | --- |
|              |       | RT1(Ω), |     | dG0(Ω)  | (cid:2) |     | (cid:3)     |     |     |
| as a result, | for v |         |     | q       | , the   | sca | lar product |     |     |
|              |       | ∈       |     | ∈       |         |     |             |     |     |
|              |       |         |     | div(v), | q       | =   | div(v)qdx   |     |     |
L2(Ω)
|     |     |     |     | (cid:104)− | (cid:105) |     | −   |     |     |
| --- | --- | --- | --- | ---------- | --------- | --- | --- | --- | --- |
(cid:90)Ω
| is defined, | which | enables | the | following |     |     |     |     |     |
| ----------- | ----- | ------- | --- | --------- | --- | --- | --- | --- | --- |

| CHAPTER    |     | 4.  | HODGE-HELMHOLTZ |     |     | DECOMPOSITION |     |         |     | AND DE | RHAM |     |     |
| ---------- | --- | --- | --------------- | --- | --- | ------------- | --- | ------- | --- | ------ | ---- | --- | --- |
| COMPLEXES: |     |     | CONTINUOUS      |     | AND | DISCRETE      |     | ASPECTS |     |        |      |     | 93  |
Definition 4.4.1 (Gradient on dG0(Ω)). The discrete gradient ( div) is defined on dG0(Ω)
∗
−
by duality:
|     |       | g           | dG0(Ω), |          | v RT1(Ω)   |        | (         | div)    | g,v       | := g, div(v) |           | .   |     |
| --- | ----- | ----------- | ------- | -------- | ---------- | ------ | --------- | ------- | --------- | ------------ | --------- | --- | --- |
|     |       |             |         |          |            |        |           |         | ∗         | h            | L2(Ω)     |     |     |
|     |       | ∀           | ∈       |          | ∀ ∈        |        | (cid:104) | −       | (cid:105) | (cid:104) −  | (cid:105) |     |     |
|     | Also, | the diagram |         | tells by | definition | that   |           |         |           |              |           |     |     |
|     |       |             |         |          |            | cG1(Ω) |           | RT1(Ω), |           |              |           |     |     |
∇ ⊥
⊂
|     |              |      | RT1(Ω), |                | cG1(Ω | (cid:2) |           | (cid:3)     |     |     |     |     |     |
| --- | ------------ | ---- | ------- | -------------- | ----- | ------- | --------- | ----------- | --- | --- | --- | --- | --- |
| as  | a result,for | v    |         |                | ϕ     | ),      | the sca   | lar product |     |     |     |     |     |
|     |              |      | ∈       |                | ∈     |         |           |             |     |     |     |     |     |
|     |              |      |         |                |       |         | ∇         | ϕ, v        |     |     |     |     |     |
|     |              |      |         |                |       |         | ⊥         | h           |     |     |     |     |     |
|     |              |      |         |                |       |         | (cid:104) | (cid:105)   |     |     |     |     |     |
| is  | defined,     | this | enables | the following: |       |         |           |             |     |     |     |     |     |
Definition 4.4.2 (Curl on Raviart-Thomas space). The discrete curl, (∇ ⊥ ) is defined on
∗
| RT1(Ω) |     | by duality: |         |     |     |        |     |           |         |                     |     |           |     |
| ------ | --- | ----------- | ------- | --- | --- | ------ | --- | --------- | ------- | ------------------- | --- | --------- | --- |
|        |     |             | RT1(Ω), |     |     | cG1(Ω) |     |           |         |                     |     |           |     |
|        |     |             | u       |     | ϕ   |        |     | (∇ ⊥      | ) ∗ u,ϕ | L2(Ω) := u,∇        | ⊥ ϕ | h .       |     |
|        |     |             | ∀ ∈     |     | ∀ ∈ |        |     | (cid:104) |         | (cid:105) (cid:104) |     | (cid:105) |     |
k
Remark 4.4.3. Classically a differential operator acting on a function that is will lose
C
|     |     |     |     |     |     | k 1. |     |     |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | ---- | --- | --- | --- | --- | --- | --- | --- |
an order of regularity and thus be By contrast, it is remarkable that, for a cellwise
C −
constant function p, p, its discrete gradient, actually takes values in RT1(Ω) (and is
|     |     |     | (   | div) | ∗   |     |     |     |     |     |     |     |     |
| --- | --- | --- | --- | ---- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
−
consequently normal-conform), hence yielding a gain in regularity. A similar remark can be
| made | on  | the discrete | curl | on  | the Raviart-Thomas |     |     | space. |     |     |     |     |     |
| ---- | --- | ------------ | ---- | --- | ------------------ | --- | --- | ------ | --- | --- | --- | --- | --- |
Second order discrete differential operators: With this formalism, we can go even fur-
ther and define second order differential operators that will mimic again continuous differential
operators. Indeed, by duality we have two adjoint complexes; the original one
|     |     |     |     |     | cG1(Ω) | ∇⊥  | RT1(Ω) |     | div dG0(Ω), |     |     |     |     |
| --- | --- | --- | --- | --- | ------ | --- | ------ | --- | ----------- | --- | --- | --- | --- |
−
|     |     |     |     |     |     | −→  |     | −→  |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
and the one defined by Definition 4.4.2 and Definition 4.4.1 that ’goes’ in the other direction:
|     |     |     |     |     | cG1(Ω) | (∇⊥)∗ | RT1(Ω) | (   | div)∗ | dG0(Ω). |     |     |     |
| --- | --- | --- | --- | --- | ------ | ----- | ------ | --- | ----- | ------- | --- | --- | --- |
−
|     |              |     |        |        |               | ←−  |       |         | ←−  |     |     |     |     |
| --- | ------------ | --- | ------ | ------ | ------------- | --- | ----- | ------- | --- | --- | --- | --- | --- |
| It  | is a complex |     | in the | sense  | of Definition |     | 4.2.1 | because |     |     |     |     |     |
|     | •            |     | dG0(Ω) | RT1(Ω) |               |     |       |         |     |     |     |     |     |
( div) ∗ and similarly for a Lagrange finite element q we have for
|     | −   |     |     | ⊂   |     |     |     |     |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
p dG0(Ω),
|     | any |           |     | using   | twice | the definition |     | of the | adjoints: |     |     |     |     |
| --- | --- | --------- | --- | ------- | ----- | -------------- | --- | ------ | --------- | --- | --- | --- | --- |
|     |     | ∈ (cid:2) |     | (cid:3) |       |                |     |        |           |     |     |     |     |
(∇ ⊥ ) ∗ ( div) ∗ p,q L2(Ω) := ( div) ∗ p,∇ ⊥ q h := p, div∇ ⊥ q L2(Ω) = 0,
|     |     | (cid:104) | −       |        | (cid:105) |          | (cid:104) − |        | (cid:105) | (cid:104) − |     | (cid:105) |     |
| --- | --- | --------- | ------- | ------ | --------- | -------- | ----------- | ------ | --------- | ----------- | --- | --------- | --- |
|     |     |           |         | dG0(Ω) |           | (∇       |             |        |           |             |     |           |     |
|     | so  | that (    | div) ∗  |        | ker       |          | ⊥ ) ∗       | RT1(Ω) |           |             |     |           |     |
|     |     | −         |         |        | ⊂         |          | |           |        |           |             |     |           |     |
|     |     |           |         |        |           | (cid:16) |             |        | (cid:17)  |             |     |           |     |
|     |     |           | (cid:2) |        | (cid:3)   |          |             |        |           |             |     |           |     |

| CHAPTER    | 4.         | HODGE-HELMHOLTZ |             |     |      | DECOMPOSITION |     |         | AND | DE RHAM |     |     |
| ---------- | ---------- | --------------- | ----------- | --- | ---- | ------------- | --- | ------- | --- | ------- | --- | --- |
| COMPLEXES: |            | CONTINUOUS      |             |     | AND  | DISCRETE      |     | ASPECTS |     |         |     | 94  |
| • By       | definition | of              | the adjoint |     | (∇ ) | (RT1(Ω)       |     | cG1(Ω). |     |         |     |     |
⊥ ∗
⊂
| Then |     |     |     |     |     | (cid:2) |     | (cid:3) |     |     |     |     |
| ---- | --- | --- | --- | --- | --- | ------- | --- | ------- | --- | --- | --- | --- |
RT1(Ω)
|                |         |        |           | u,Ψ  |        |           | (           | div) ∗ div(u),Ψ |       | h         |     | (4.22) |
| -------------- | ------- | ------ | --------- | ---- | ------ | --------- | ----------- | --------------- | ----- | --------- | --- | ------ |
|                |         |        | ∀         |      | ∈      |           | (cid:104) − |                 |       | (cid:105) |     |        |
| is proprely    | defined | and    |           |      |        |           |             |                 |       |           |     |        |
|                |         |        |           | u,Ψ  | RT1(Ω) |           | ∇           | (∇ )            | (u),Ψ | ,         |     | (4.23) |
|                |         |        |           |      |        |           |             | ⊥ ⊥             | ∗     | h         |     |        |
|                |         |        |           | ∀    | ∈      |           | (cid:104)   |                 |       | (cid:105) |     |        |
| is too. Adding |         | (4.23) | to (4.22) | will | lead   | to the    | discrete    | equivalent      |       | of        |     |        |
|                |         |        |           | ∆u   |        | ∇div(u)+∇ |             | curl(u),        |       |           |     |        |
|                |         |        |           |      | =      |           |             | ⊥               |       |           |     |        |
with∆uthevectorialLaplacian, whichiscalledtheHodge-LaplacianoftheRaviart-Thomas
space:
RT1(Ω)
Definition 4.4.3 (Hodge-Laplacian). The Hodge-Laplacian ∆u of u is defined
∈
R2
| for a domain | Ω   |     | as  |     |     |     |     |     |     |     |     |     |
| ------------ | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
⊂
RT1(Ω)
|     | Ψ   |     |     | ∆u,Ψ      | h         | :=  | ( div) | ∗ div(u)+∇ |     | ⊥ (∇ ⊥ ) ∗ (u),Ψ | .   |     |
| --- | --- | --- | --- | --------- | --------- | --- | ------ | ---------- | --- | ---------------- | --- | --- |
|     | ∀   | ∈   |     | (cid:104) | (cid:105) |     | −      |            |     |                  |     |     |
(cid:29)h
(cid:28)
Its exact expression will depend on the definition of the discrete scalar product .. .
h
(cid:104) (cid:105)
Discrete Hodge-Helmholtz decomposition with boundary conditions
As discussed in section 4.2, a de Rham complex yields naturally a Hodge-Helmholtz decompos-
ition (HHD) of the type Theorem 4.2.1. In fact from there it is possible to derive other types
of HHD, of which, ones with boundary conditions Theorem 2.3.1. In this case we want to show
| that at | the discrete | scale | a   | decomposition |     | of  | the following |     | type |     |     |        |
| ------- | ------------ | ----- | --- | ------------- | --- | --- | ------------- | --- | ---- | --- | --- | ------ |
|         |              |       |     |               |     | u = | ∇ϕ+u          | ,   |      |     |     | (4.24) |
Ψ
stands, with
|     |     |     |     | div(u | Ψ ) = | 0 and | u Ψ | n ∂Ω = | u b n | ∂Ω , |     |     |
| --- | --- | --- | --- | ----- | ----- | ----- | --- | ------ | ----- | ---- | --- | --- |
|     |     |     |     |       |       |       | ·   |        | ·     |      |     |     |
|     |     |     |     |       |       |       |     | |      |       | |    |     |     |
and ϕ is a field and u,u are vector fields. Formally, at the continuous level, a way to show
Ψ
the existence of such decomposition is to first get a variational formulation defining ϕ. In this
| aim, we | apply minus |     | the divergence |     | operator |        | to (4.24), | yielding |     |     |     |        |
| ------- | ----------- | --- | -------------- | --- | -------- | ------ | ---------- | -------- | --- | --- | --- | ------ |
|         |             |     |                |     |          | div(u) | =          | ∆ϕ.      |     |     |     | (4.25) |
|         |             |     |                |     |          | −      |            | −        |     |     |     |        |
Then multiplying by a test function g (4.25) and integrating by parts we get:
∂ϕ
|     |     | ∇ϕ        | ∇gdx |     |          | gdΓ | =         | u ∇gdx |     | u ngdΓ.  |     | (4.26) |
| --- | --- | --------- | ---- | --- | -------- | --- | --------- | ------ | --- | -------- | --- | ------ |
|     |     |           | ·    | −   | ∂n       |     |           | ·      | −   | ·        |     |        |
|     |     | (cid:90)Ω |      |     | ∂Ω       |     | (cid:90)Ω |        |     | ∂Ω       |     |        |
|     |     |           |      |     | (cid:90) |     |           |        |     | (cid:90) |     |        |

CHAPTER 4. HODGE-HELMHOLTZ DECOMPOSITION AND DE RHAM
COMPLEXES: CONTINUOUS AND DISCRETE ASPECTS 95
But u n = u n , so
Ψ ∂Ω b ∂Ω
· | · |
∇ϕ n = u n u n .
∂Ω ∂Ω b ∂Ω
· | · | − · |
Plugging this in (4.26) we obtain the following variational formula
Find ϕ H1 g H1 ∇ϕ ∇gdx = u ∇gdx gu ndΓ. (4.27)
b
∈ ∀ ∈ · · − ·
(cid:90)Ω (cid:90)Ω (cid:90) ∂Ω
Alternatively, if u is H(div), then (4.27) can be rephrased as
Find ϕ H1 g H1 ∇ϕ ∇gdx = (divu)gdx+ g(u n u n) dΓ.
b
∈ ∀ ∈ · − · − ·
(cid:90)Ω (cid:90)Ω (cid:90) ∂Ω
(4.28)
In this context, the existence of such ϕ can be proven with the Lax-Milgram lemma modulo
the compatibility condition stands on the boundary velocity u [43]
b
u ndΓ = 0. (4.29)
b
·
(cid:90)
∂Ω
Atthecontinuouslevel, theformalismofcomplexesenablestoproperlydefinetheoperators
contributing to (4.27) or (4.28). At the discrete level, however, adaptations should be made in
order to obtain the relevant equivalent of (4.27) or (4.28).
Resolution of the discrete system in the original complex: Staying in the original
complex;
cG1(Ω) ∇⊥ RT1(Ω)
−
div dG0(Ω),
−→ −→
in order to obtain the decomposition while trying to get divu = 0 is in fact a dead end.
Ψ
Indeed, for any g dG0, we have
∈
gdivu = gdiv(u u )
Ψ ϕ
−
(cid:90)Ω (cid:90)Ω
= gdivu gdiv(u )
ϕ
−
(cid:90)Ω (cid:90)Ω
= gdivu u div g
ϕ ∗
− ·
(cid:90)Ω (cid:90)Ω
= gdivu+ ( div) ϕ ( div) g.
∗ ∗
− · −
(cid:90)Ω (cid:90)Ω
Therefore, the divergence of u is zero provided we have:
Ψ
( div) ϕ ( div) g = gdivu.
∗ ∗
− · − −
(cid:90)Ω (cid:90)Ω
This last equality matches with a Laplace variational formulation on ϕ, but it clearly does not
include any boundary condition. It matches neither (4.27), nor (4.28).

| CHAPTER    | 4.  | HODGE-HELMHOLTZ |     |     | DECOMPOSITION |     |         |     | AND DE | RHAM |     |
| ---------- | --- | --------------- | --- | --- | ------------- | --- | ------- | --- | ------ | ---- | --- |
| COMPLEXES: |     | CONTINUOUS      |     | AND | DISCRETE      |     | ASPECTS |     |        |      | 96  |
Increasing the number of degrees of freedom: In order to take into account correctly
the von Neumann boundary conditions, we propose to add to dG0 some degrees of freedom on
the boundary faces. This procedure is inspired by [107, 83], in which it is proposed to consider
discontinuous versions of the Raviart-Thomas and N´ed´elec finite elements. By their very
definition, discontinuous Raviart-Thomas finite elements lose the normal conformity, which
means that the divergence in the sense of the distributions will include both the classical
divergence in the cells and a normal component jump across each face, multiplied by a Dirac
mass. In our case, we only need to assume ’non-conformity’ on boundary faces so that the
divergence includes a normal jump across the boundary faces; which can be seen as the jump
between the value of the trace, and a boundary value that we would like to impose.
The construction of the approximation space for the extraction of the potential (4.28)
is strongly inspired by what was done in [45, Section 3.2.2]. In fact, we define the space of
| potentials | as  |     |     |     |           |     |          |     |     |     |     |
| ---------- | --- | --- | --- | --- | --------- | --- | -------- | --- | --- | --- | --- |
|            |     |     |     | Φ   | := dG0(Ω) |     | dG0(∂Ω), |     |     |     |     |
×
which means that we have one degree of freedom on each cell, and one degree of freedom on
each boundary face, representing the boundary trace. For an element g Φ, we will denote by
∈
g for K the value associated with the cell and by g for f b the value associated with
| K        |        |     |         |             |           |      |         | f           |     |     |     |
| -------- | ------ | --- | ------- | ----------- | --------- | ---- | ------- | ----------- | --- | --- | --- |
|          | ∈ C    |     |         |             |           |      |         |             | ∈ F |     |     |
| boundary | faces. | The | space Φ | is equipped |           | with | a graph | L2 product: |     |     |     |
|          |        |     | p,q Φ   |             | p,q       | :=   | K p     | q +         | f p | q . |     |
|          |        |     |         |             | Φ         |      | K       | K           | f   | f   |     |
|          |        | ∀   | ∈       | (cid:104)   | (cid:105) | K    | | |     |             | | | |     |     |
f b
|     |     |     |     |     |     | (cid:80) ∈C |     |     | (cid:80)∈F |     |     |
| --- | --- | --- | --- | --- | --- | ----------- | --- | --- | ---------- | --- | --- |
We now wish to define the gradient of each of the function basis. For this, we start by defining
an extension of the (minus) divergence operator of the original complex to the space Φ, which
| we denote | by  | div: |     |     |     |     |     |     |     |     |     |
| --------- | --- | ---- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
−
|         |        |              |          |        |     |          | divu        | =   | divu |     |     |
| ------- | ------ | ------------ | -------- | ------ | --- | -------- | ----------- | --- | ---- | --- | --- |
|         |        |              |          | RT1(Ω) |     |          |             |     | K    |     |     |
|         |        |              | u        |        |     |          | −           | K − |      |     |     |
|         |        |              |          |        |     |          | divu        | =(u | n) . |     |     |
|         |        |              | ∀ ∈      |        |     | (cid:40) | (cid:1)f    |     | f    |     |     |
|         |        |              |          |        |     | (cid:0)− |             |     | ·    |     |     |
|         |        |              |          |        |     | (cid:0)  | (cid:1)     |     |      |     |     |
| We then | define | the discrete | gradient |        | div | as       | the adjoint | of  | div: |     |     |
∗
|            |       |           |     |                   | −   |               |           |     | −            |                   |     |
| ---------- | ----- | --------- | --- | ----------------- | --- | ------------- | --------- | --- | ------------ | ----------------- | --- |
|            |       |           |     | Gradient(cid:0)on |     | dG(cid:1)0(Ω) | dG0(∂Ω)). |     |              |                   |     |
| Definition | 4.4.4 | (Discrete |     |                   |     |               |           |     | The extended | discrete gradient |     |
×
|     | is defined | by  | duality | by: |     |     |     |     |     |     |     |
| --- | ---------- | --- | ------- | --- | --- | --- | --- | --- | --- | --- | --- |
div ∗
−
| (cid:0) (cid:1) |          |          | RT1(Ω)   |     |                |           |         |           |            |           |     |
| --------------- | -------- | -------- | -------- | --- | -------------- | --------- | ------- | --------- | ---------- | --------- | --- |
|                 |          | u        |          | ϕ   | Φ              | u,        | div     | ∗ϕ h      | = divu,ϕ   | Φ .       |     |
|                 |          | ∀ ∈      |          | ∀   | ∈              | (cid:104) | −       | (cid:105) | (cid:104)− | (cid:105) |     |
|                 |          |          |          |     |                |           | (cid:0) | (cid:1)   |            |           |     |
| Then            | this new | operator | verifies |     | the following: |           |         |           |            |           |     |
Lemma 4.4.2. Let div ∗ the operator given by Definition 4.4.4, div ∗ has a kernel
|     |     |     | −   |     |     |     |     |     |     | −   |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
reduced to the functions constant on each connected component of Ω, and:
|     |     |     | (cid:0) (cid:1) |     |     |      |        |     |     | (cid:0) (cid:1) |     |
| --- | --- | --- | --------------- | --- | --- | ---- | ------ | --- | --- | --------------- | --- |
|     |     |     | dim             |     | div | ∗(Φ) | = # +# | b   | b , |                 |     |
0
|     |     |     |     | −        |     |          | C   | F   | −   |     |     |
| --- | --- | --- | --- | -------- | --- | -------- | --- | --- | --- | --- | --- |
|     |     |     |     | (cid:16) |     | (cid:17) |     |     |     |     |     |
(cid:0) (cid:1)
where b is the zeroth Betti number which is equal to the number of connected components of
0
Ω.

CHAPTER 4. HODGE-HELMHOLTZ DECOMPOSITION AND DE RHAM
COMPLEXES: CONTINUOUS AND DISCRETE ASPECTS 97
Proof. Let g Φ such that div ∗g = 0 then it implies that
∈ −
(cid:0) (cid:1)
Ψ RT1(Ω) div ∗g, Ψ
h
= 0. (4.30)
∀ ∈ (cid:104) − (cid:105)
(cid:0) (cid:1)
SincebydefinitionofΦ,g = α 1 + α 1 forsomecoefficients(α ) ,(α ) R
K K f f K K f f b
K f b ∈C ∈F ⊂
(wherewedenote1 ,1 thei(cid:80)n ∈ d C icatorfun(cid:80)c ∈ t F ionof,respectively,anycellK/facef),(4.30)yields
K f
Ψ RT1(Ω) α K div ∗ 1 K ,Ψ h + α f div ∗ 1 f ,Ψ h = 0.
∀ ∈ (cid:104) − (cid:105) (cid:104) − (cid:105)
K f b
(cid:80) ∈C (cid:0) (cid:1) (cid:80)∈F (cid:0) (cid:1)
Using the definition of div ∗ Definition 4.4.4, we get
−
(cid:0) (cid:1)
α divΨdx+ α Ψ n dΓ = 0. (4.31)
K f f
− ·
K (cid:90) K f b (cid:90) f
(cid:80) ∈C (cid:80)∈F
Suppose that Ψ is Ψ , the Raviart-Thomas basis function associated with an interior face σ.
σ
Then, denoting Ω the i th connected component of Ω in which σ is contained, then
i
−
α divΨ dx = 0.
K σ
K (σ) Ωi (cid:90) K
∈C ∩
(cid:80)
Since
divΨ = Ψ n = σ ε (σ), (4.32)
σ σ K
· | |
K ∂K
(cid:90) (cid:90)
we have
σ α = 0.
K
| |
K (σ) Ωi
∈C ∩
(cid:80)
As a consequence, on a given connected component of Ω, Ω
i
K,L Ω α = α . (4.33)
i K L
∀ ∈ C ∩
We now restart from (4.31), and test with a Raviart-Thomas basis function associated with a
boundary face σ ∂Ω , where Ω still the i th connected component of Ω:
i i
∈ −
α divΨ dx+α Ψ n dΓ = 0,
−
Kσ σ σ σ
·
σ
(cid:90)
Kσ
(cid:90)
σ
yielding by (4.32) α = α , so that all the coefficients (α ) are equal in each connected
σ Kσ K K
∈C
component Ω i . Finally, applying the rank theorem to div ∗ gives:
−
(cid:0) (cid:1)
dimΦ = # +# b = range div ∗+ker div ∗ = range div ∗+b 0 ,
C F − − −
(cid:0) (cid:1) (cid:0) (cid:1) (cid:0) (cid:1)
which ends the proof.

| CHAPTER    |     | 4.  | HODGE-HELMHOLTZ |     |     |              | DECOMPOSITION |     |         | AND | DE  | RHAM |     |
| ---------- | --- | --- | --------------- | --- | --- | ------------ | ------------- | --- | ------- | --- | --- | ---- | --- |
| COMPLEXES: |     |     | CONTINUOUS      |     |     | AND DISCRETE |               |     | ASPECTS |     |     |      | 98  |
We are now interested in the resolution of (4.27), but in the approximation space Φ. Given
u RT1(Ω),
∈
| Find | ϕ   | Φ   | g Φ |           | div | ∗ϕ  | div | ∗g        | = u,      | div | ∗g          | gu ndΓ. | (4.34) |
| ---- | --- | --- | --- | --------- | --- | --- | --- | --------- | --------- | --- | ----------- | ------- | ------ |
|      |     |     |     |           |     |     |     | h         |           |     | h           | b       |        |
|      | ∈   | ∀   | ∈   | (cid:104) | −   | · − |     | (cid:105) | (cid:104) | −   | (cid:105) − | ·       |        |
(cid:90) ∂Ω
|     |      |      |     | (cid:0) | (cid:1) | (cid:0) | (cid:1) |     |     | (cid:0) | (cid:1) |     |     |
| --- | ---- | ---- | --- | ------- | ------- | ------- | ------- | --- | --- | ------- | ------- | --- | --- |
| We  | then | have |     |         |         |         |         |     |     |         |         |     |     |
Proposition 4.4.1 (Existence of ϕ). Suppose that (4.29) holds and Ω is connected. Then the
| variational |     | formulation |     | (4.34) | has | a unique |     | solution | up  | to a constant. |     |     |     |
| ----------- | --- | ----------- | --- | ------ | --- | -------- | --- | -------- | --- | -------------- | --- | --- | --- |
Proof.
For proving the proposition, we want to use the Lax-Milgram theorem. However, as
we have Neumann boundary conditions, the bilinear form is not coercive. In parallel, the
kernel of div was studied in Lemma 4.4.2, and is equal to span 1 since Ω is connected
|     |     |     | ∗   |     |     |     |     |     |     |     | Ω   |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
|     |     | −   |     |     |     |     |     |     |     |     | {   | }   |     |
. Therefore, instead of considering the formulation (4.34), we consider the same problem, but
|     |     | (cid:0) | (cid:1) |     |     |     |     |     |     |     |     |     |     |
| --- | --- | ------- | ------- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
1
| posed | on  | the | quotient | space | Φ/span |     | :   |     |     |     |     |     |     |
| ----- | --- | --- | -------- | ----- | ------ | --- | --- | --- | --- | --- | --- | --- | --- |
Ω
|     |     |     |     |     |        | {      | }   |     |     |        |     |     |     |
| --- | --- | --- | --- | --- | ------ | ------ | --- | --- | --- | ------ | --- | --- | --- |
|     |     |     |     |     | Find ϕ | Φ/span |     | 1   | g   | Φ/span | 1   |     |     |
|     |     |     |     |     |        |        |     | Ω   |     |        | Ω   |     |     |
|     |     |     |     |     |        | ∈      | {   | }   | ∀ ∈ |        | { } |     |     |
(4.35)
|     |     |     |             | div ∗ϕ, |     | div ∗g | =           | u,          | div ∗g |               | gu  | ndΓ. |     |
| --- | --- | --- | ----------- | ------- | --- | ------ | ----------- | ----------- | ------ | ------------- | --- | ---- | --- |
|     |     |     | (cid:104) − |         | −   |        | (cid:105) h | (cid:104) − |        | (cid:105) h − | b   | ·    |     |
∂Ω
(cid:90)
|     |     |     | (cid:0) | (cid:1) | (cid:0) | (cid:1) |     | (cid:0) | (cid:1) |     |     |     |     |
| --- | --- | --- | ------- | ------- | ------- | ------- | --- | ------- | ------- | --- | --- | --- | --- |
This is possible, because the right hand side is invariant under a constant shift: if g Φ, then
∈
|     |     |     | 1   |     |     |     |     |     |     |     | 1   |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
for any λ span Ω ; by Lemma 4.4.2, ker( div ∗) = span Ω so that
|     |           | ∈       | {       | }         |          |        |     | −       |           | {       | }             |          |     |
| --- | --------- | ------- | ------- | --------- | -------- | ------ | --- | ------- | --------- | ------- | ------------- | -------- | --- |
|     |           |         |         |           |          |        |     | (cid:0) | (cid:1)   |         |               |          |     |
|     | u,        | div     | ∗(g+λ)  |           |          | (g+λ)u |     | ndΓ.=   | u,        | div     | ∗g            | gu ndΓ   |     |
|     | (cid:104) | −       |         | (cid:105) | h −      |        | b   | ·       | (cid:104) | −       | (cid:105) h − | b ·      |     |
|     |           |         |         |           | ∂Ω       |        |     |         |           |         |               | ∂Ω       |     |
|     |           |         |         |           | (cid:90) |        |     |         |           |         |               | (cid:90) |     |
|     |           | (cid:0) | (cid:1) |           |          |        |     |         |           | (cid:0) | (cid:1)       |          |     |
|     |           |         |         |           |          |        |     |         |           | λ       | u b ndΓ       |          |     |
|     |           |         |         |           |          |        |     |         | −         |         | ·             |          |     |
∂Ω
(cid:90)
|     |     |     |     |     |     |     |     |     | u,        |     | ∗g          | gu ndΓ |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --------- | --- | ----------- | ------ | --- |
|     |     |     |     |     |     |     |     |     | =         | div | h           | b      |     |
|     |     |     |     |     |     |     |     |     | (cid:104) | −   | (cid:105) − | ·      |     |
(cid:90) ∂Ω
|         |       |        |           |       |            |       |           |     |     | (cid:0)     | (cid:1) |     |     |
| ------- | ----- | ------ | --------- | ----- | ---------- | ----- | --------- | --- | --- | ----------- | ------- | --- | --- |
| because |       | (4.29) | holds.    | Then, | in (4.35): |       |           |     |     |             |         |     |     |
|         | • the | left   | hand side | is    | a bilinear | form, | coercive, |     | and | continuous, |         |     |     |
•
|     | the | right | hand | side is | a linear | form | that | is continuous. |     |     |     |     |     |
| --- | --- | ----- | ---- | ------- | -------- | ---- | ---- | -------------- | --- | --- | --- | --- | --- |
Following Lax-Milgram theorem, (4.35) has a unique solution, and so (4.34) has a unique
| solution |     | up to | a constant, |     | which | ends the | proof. |     |     |     |     |     |     |
| -------- | --- | ----- | ----------- | --- | ----- | -------- | ------ | --- | --- | --- | --- | --- | --- |
We now obtain a discrete counterpart of the Hodge-Helmholtz decomposition of The-
RT1(Ω):
| orem | 2.3.1, | but | on  |     |     |     |     |     |     |     |     |     |     |
| ---- | ------ | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |

| CHAPTER    |     | 4.  | HODGE-HELMHOLTZ |     |     |              | DECOMPOSITION |         |     | AND | DE RHAM |     |     |
| ---------- | --- | --- | --------------- | --- | --- | ------------ | ------------- | ------- | --- | --- | ------- | --- | --- |
| COMPLEXES: |     |     | CONTINUOUS      |     |     | AND DISCRETE |               | ASPECTS |     |     |         |     | 99  |
Theorem 4.4.4 (Discrete Hodge-Helmholtz Decomposition with boundary conditions). Let
Ω R2, an open connected Lipschitz set and u such that (4.29) is verified. Then, for all
b
⊂
| u   | RT1(Ω) | it  | exists | a unique | decomposition: |     |     |     |     |     |     |     |     |
| --- | ------ | --- | ------ | -------- | -------------- | --- | --- | --- | --- | --- | --- | --- | --- |
∈
|     |     |     |     |     |     | u = | div | ∗ϕ+u | ,   |     |     |     | (4.36) |
| --- | --- | --- | --- | --- | --- | --- | --- | ---- | --- | --- | --- | --- | ------ |
Ψ
−
|       |     |     |        |          |       |               | (cid:0) | (cid:1) |     |     |     |     |     |
| ----- | --- | --- | ------ | -------- | ----- | ------------- | ------- | ------- | --- | --- | --- | --- | --- |
| where |     | div | is the | operator | given | by Definition |         | 4.4.4   | and |     |     |     |     |
∗
−
|     | (cid:0) | (cid:1) |     |         |     |       |        |     |     |     |       |       |        |
| --- | ------- | ------- | --- | ------- | --- | ----- | ------ | --- | --- | --- | ----- | ----- | ------ |
|     | u       | U       | v   | RT1(Ω), |     | K     | div(v) |     | 0,  | σ   | b v n | u n , |        |
|     | Ψ       | Ψ       | :=  |         |     |       |        | K = |     |     | σ =   | b σ   | (4.37) |
|     |         | ∈       | {   | ∈       |     | ∀ ∈ C |        |     | ∀   | ∈ F | ·     | · }   |        |
| and | ϕ       | Φ/span  | 1   | .       |     |       |        |     |     |     |       |       |        |
Ω
|     | ∈   |     | {   | }   |     |     |     |     |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
Proof. For a given u RT1(Ω), if the condition (4.29) holds, then because Ω is connected by
∈
Proposition 4.4.1, (4.34) has a unique solution ϕ up to a constant. Then we define
|     |     |     |     |     |     | u := | u   | div | ∗ϕ. |     |     |     |     |
| --- | --- | --- | --- | --- | --- | ---- | --- | --- | --- | --- | --- | --- | --- |
Ψ
|     |     |     |     |     |     |     | −   | −       |         |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | ------- | ------- | --- | --- | --- | --- |
|     |     |     |     |     |     |     |     | (cid:0) | (cid:1) |     |     |     |     |
Following this definition of u , and the fact that ϕ is solution of (4.35), we have
Ψ
|     |     |     | g   | Φ/span | 1   |     | u ,       |     | ∗g        |     | gu ndΓ. |     |        |
| --- | --- | --- | --- | ------ | --- | --- | --------- | --- | --------- | --- | ------- | --- | ------ |
|     |     |     |     |        |     | Ω   | Ψ         | div | h         | =   | b       |     | (4.38) |
|     |     |     | ∀   | ∈      | {   | }   | (cid:104) | −   | (cid:105) |     | ·       |     |        |
(cid:90) ∂Ω
(cid:0) (cid:1)
|     |     |     |     |     |     | 1   |     |     | b,  |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
We test the above equality with g = where σ and find for the left hand side of (4.38),
|     |     |     |     |           |         | σ         |           | ∈ F    |           |     |       |     |     |
| --- | --- | --- | --- | --------- | ------- | --------- | --------- | ------ | --------- | --- | ----- | --- | --- |
|     |     |     |     |           |         | 1         | 1         |        |           |     |       |     |     |
|     |     |     |     | u ,       | div     | ∗         | =         | , divu |           | = σ | u n , |     |     |
|     |     |     |     | Ψ         |         | σ h       | σ         |        | Ψ Φ       |     | Ψ σ   |     |     |
|     |     |     |     | (cid:104) | −       | (cid:105) | (cid:104) | −      | (cid:105) | |   | | ·   |     |     |
|     |     |     |     |           | (cid:0) | (cid:1)   |           |        |           |     |       |     |     |
by definition of div and div ∗. The right hand side of (4.38) tested with 1 is
σ
|     |     |     | −   |         | −   |         |     |     |     |     |     |     |     |
| --- | --- | --- | --- | ------- | --- | ------- | --- | --- | --- | --- | --- | --- | --- |
|     |     |     |     | (cid:0) |     | (cid:1) |     |     |     |     |     |     |     |
|     |     |     |     |         |     | 1       | u n | = σ | u n | .   |     |     |     |
|     |     |     |     |         |     | σ       | b   |     | b   | σ   |     |     |     |
|     |     |     |     |         |     |         | ·   | | | | ·   |     |     |     |     |
(cid:90) ∂Ω
|     |     |     |     |     | 1   |     |     | b   |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
Therefore, (4.38) tested with σ for any σ gives (u ) = u n σ . Now testing (4.38) by
|     |     |     |     |     |     |     | ∈ F |     | Ψ   | σ   | b · |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
1
K yields
|     |     |     |     |     |     | u ,       | div | 1   | =         | 0   |     |     |     |
| --- | --- | --- | --- | --- | --- | --------- | --- | --- | --------- | --- | --- | --- | --- |
|     |     |     |     |     |     | Ψ         |     | ∗ K | h         |     |     |     |     |
|     |     |     |     |     |     | (cid:104) | −   |     | (cid:105) |     |     |     |     |
which, by definition of div yields div(u(cid:0) Ψ ) K =(cid:1)0. This ends the proof.
−
Remark 4.4.5. While Theorem 4.4.4 is stated for a domain Ω connected, the proof can be ad-
apted to a multiply-connected domain by noticing that Proposition 4.4.1 stands in such domains.
Indeed posing the problem (4.35) in Φ quotiented by the functions constant in each connected
component of Ω will lead to this extended result. In pratice, for our purpose, examining the
| case | of multiply-connected |     |     |     | domains | is not | of much | interest. |     |     |     |     |     |
| ---- | --------------------- | --- | --- | --- | ------- | ------ | ------- | --------- | --- | --- | --- | --- | --- |
Mimickingargumentsgivenin[45],wecangofurtherandgivethefollowingcharacterization
of U :
Ψ

| CHAPTER    |     | 4.  | HODGE-HELMHOLTZ |     |     |     | DECOMPOSITION |     |         | AND | DE RHAM |     |     |
| ---------- | --- | --- | --------------- | --- | --- | --- | ------------- | --- | ------- | --- | ------- | --- | --- |
| COMPLEXES: |     |     | CONTINUOUS      |     |     | AND | DISCRETE      |     | ASPECTS |     |         |     | 100 |
Proposition 4.4.2 (Characterization of the space of U ). Suppose that the hypothesis of
Ψ
Theorem 4.4.4 stand: we denote by 0 and 0 the ϕ and Ψ components of the decomposition
|     |          |        |       |       |            | ϕ   |        | Ψ   |             |     |     |     |     |
| --- | -------- | ------ | ----- | ----- | ---------- | --- | ------ | --- | ----------- | --- | --- | --- | --- |
| of  | the null | vector | field | given | by Theorem |     | 4.4.4. | We  | also denote | by: |     |     |     |
•
ψ k the continuous finite element that is equal to one on an interior point k, and 0 in the
|     | other | points, |     |     |     |     |     |     |     |     |     |     |     |
| --- | ----- | ------- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
• ψ the continuous finite element that is equal to one on all the points of the connected
∂
k
|     | component |     | ∂   | of the | boundary | of  | Ω, and | 0 in | the other | points. |     |     |     |
| --- | --------- | --- | --- | ------ | -------- | --- | ------ | ---- | --------- | ------- | --- | --- | --- |
k
RT1(Ω),
We consider u and u ϕ , u Ψ its decomposition provided by Theorem 4.4.4, then
∈
|     |     |     |     |     | u   | 0   | ∇   | span     | ψ ,ψ | .        |     |     | (4.39) |
| --- | --- | --- | --- | --- | --- | --- | --- | -------- | ---- | -------- | --- | --- | ------ |
|     |     |     |     |     |     | Ψ Ψ |     | ⊥        | k    | ∂        |     |     |        |
|     |     |     |     |     |     | −   | ∈   |          | {    | k}       |     |     |        |
|     |     |     |     |     |     |     |     | (cid:20) |      | (cid:21) |     |     |        |
Proof. Fora given u , we defineU = div ∗ϕ, u RT1(Ω) ,given by the decomposition
|     |     |     |     | b   |     | ϕ   |     |     |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
|     |     |     |     |     |     |     | −   |     | ∈   |     |     |     |     |
Theorem 4.4.4. This space is, followin(cid:110)g Lemma 4.4.2 in the ca(cid:111)se of a connected domain, of
|     |     |     |     |     |     |     | (cid:0) | (cid:1) |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | ------- | ------- | --- | --- | --- | --- | --- |
dimension # +# b 1. Then we remark that if ψ span ψ ,ψ , then div ∇ ψ = 0
|     |     |     |     |     |     |     |     |     |     | k   | ∂   |     | ⊥   |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
|     |     | C   |     | F − |     |     |     |     | ∈   | {   | k}  |     |     |
and ∇ ψ n = 0 on any boundary face f. As a consequence, the discrete Hodge-Helmholtz
|     | ⊥   |     | f   |     |     |     |     |     |     |     |     | (cid:0) | (cid:1) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | ------- | ------- |
·
decomposition of 0 +∇ ψ has a zero ϕ component. This means that:
|     |     |     |     | Ψ ⊥ |     |     |          |      |          |     |     |     |        |
| --- | --- | --- | --- | --- | --- | --- | -------- | ---- | -------- | --- | --- | --- | ------ |
|     |     |     |     |     | 0   | +∇  | span     | ψ ,ψ |          | U . |     |     | (4.40) |
|     |     |     |     |     | Ψ   |     | ⊥        | k    | ∂        | Ψ   |     |     |        |
|     |     |     |     |     |     |     |          | {    | k} ⊂     |     |     |     |        |
|     |     |     |     |     |     |     | (cid:20) |      | (cid:21) |     |     |     |        |
If we denote by r the number of holes of the domain Ω, then r+1 is the number of connected
components of ∂Ω. Let # the number of points on the mesh , then, remarking that the
|     |     |     |     |     | P   |     |     |     |     | M   |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
number of boundary points is equal to the number of boundary faces, the number of interior
|     |     |     |     | b.  |     |     |     |     |     |     |     | b)  |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
points is # # Since there are as many ψ as interior points (# # and as many
|     |     | P         | − F        |     |     |          |     | k       |      |     | P − F |     |     |
| --- | --- | --------- | ---------- | --- | --- | -------- | --- | ------- | ---- | --- | ----- | --- | --- |
| ψ ∂ | as  | connected | components |     | of  | ∂Ω (r+1) |     | We then | have |     |       |     |     |
k
b+r+1,
|     |          |               |           |      | dimspan | ψ        | ,ψ   | = #      | #   |          |       |     |      |
| --- | -------- | ------------- | --------- | ---- | ------- | -------- | ---- | -------- | --- | -------- | ----- | --- | ---- |
|     |          |               |           |      |         | { k      | ∂ k} | P        | − F |          |       |     |      |
| and | as       | a consequence |           |      |         |          |      |          |     |          |       |     |      |
|     |          |               |           | dim∇ |         | span     | ψ ,ψ |          | = # | # b+r.   |       |     |      |
|     |          |               |           |      | ⊥       |          | k    | ∂ k}     |     |          |       |     |      |
|     |          |               |           |      |         |          | {    |          | P − | F        |       |     |      |
|     |          |               |           |      |         | (cid:20) |      | (cid:21) |     |          |       |     |      |
| In  | parallel | we            | also      | have |         |          |      |          |     |          |       |     |      |
|     | dimU     |               | dimRT1(Ω) |      | dimU    |          |      |          |     | b        |       |     | b+1. |
|     |          | Ψ =           |           |      |         | ϕ        | = #  | #        | +#  | 1        | = # # | #   |      |
|     |          |               |           |      | −       |          | F    | −        | C F | −        | F − C | − F |      |
|     |          |               |           |      |         |          |      | (cid:16) |     | (cid:17) |       |     |      |
and the Euler relation [108, Theorem 2.27 p.140] states in 2 spaces dimensions that
|     |     |     |     |     |     | #   | #   | +#  | = 1 | r,  |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
|     |     |     |     |     |     | C   | − F | P   | −   |     |     |     |     |

CHAPTER 4. HODGE-HELMHOLTZ DECOMPOSITION AND DE RHAM
COMPLEXES: CONTINUOUS AND DISCRETE ASPECTS 101
which leads to
dimU = # +r 1 # b+1 = # # b+r = dim∇ span ψ ,ψ .
Ψ P − − F P − F ⊥ { k ∂ k}
(cid:20) (cid:21)
This last equality, combined with (4.40) gives (4.39).
Remark 4.4.6. We remark that the link between the continuous Hodge-Helmholtz decompos-
ition Theorem 2.3.1 and its discrete counterpart Theorem 4.4.4 is purely geometrical. The
properties that are preserved can indeed be stated in terms of decomposition into a divergence
free and a curl free component, with similar boundary conditions nonetheless we do not have
any property of approximation on this discrete decomposition, namely if uh RT1(Ω) is an
∈
approximation of a H1 vector field u, we do not have any proof of componentwise convergence
(e.g. uh u ).
ϕ ϕ
→
Remark 4.4.7. As noted in section 1.3, a line of research on staggered schemes has based the
approximationofthevectorialunknownontheCrouzeix-Raviartspaceoneachofitscomponents
[55, 56, 57, 65, 66, 67, 68]. However, we are not aware of the existence of a Hodge-Helmholtz
decomposition for these finite element spaces but it is now clear that the HHD is a key tool
in our analysis, so that the choice of this space is not, in our case, relevant as long as such
decomposition is not exhibited. In fact, Crouzeix-Raviart elements are also widely used for
incompressible Navier-Stokes equations, it is known that they suffer from a lack of orthogonality
[109, 110]. To make it clear, at the continuous level, two fundamental properties arise on the
incompressible Stokes equations:
ν∆u+∇p = f
−
, (4.41)

 div(u) = 0
with no-slip boundary conditions u ∂Ω = 0. f being a source term, for example the gravity, and
|
ν the viscosity. These two properties are
(i) the well-known inf-sup condition [87]
q div(u)dx
β R+ inf sup (cid:90)Ω β > 0, (4.42)
∃ ∈ ∗ q ∈ L2 0 (Ω)/ { 0 } v ∈ H 0 1(Ω)d/ { 0 } (cid:107) q (cid:107) L2 0 (Ω)(cid:107) ∇u (cid:107) L2(Ω)d ≥
(ii) an invariance property, called pressure robustness[109] linked to the Hodge-Helmholtz
decomposition; if the source term f is modified by a purely potential ∇ϕ then only the
pressure is affected:
p p+ϕ
f f+∇ϕ then the solution .
−→ u −→ u
(cid:18) (cid:19) (cid:18) (cid:19)
A discrete inf-sup condition with a constant independent of the mesh is known to be essential for
adiscretizationoftheincompressible(Navier-)Stokesequations, incontrastthesecondcondition

| CHAPTER    | 4. HODGE-HELMHOLTZ |     | DECOMPOSITION | AND DE RHAM |     |
| ---------- | ------------------ | --- | ------------- | ----------- | --- |
| COMPLEXES: | CONTINUOUS         | AND | DISCRETE      | ASPECTS     | 102 |
is generally not verified for classical mixed finite element methods [111] and in particular when
| the discrete | velocity space | is chosen | as the Crouzeix-Raviart | elements. |     |
| ------------ | -------------- | --------- | ----------------------- | --------- | --- |
4.5 Conclusion
In this chapter we have presented the formalism of de Rham complexes and shown that it
is the natural set up behind Hodge-Helmholtz decompositions. Replicating at the discrete
level its principles is possible modulo an harmonic forms isomorphism stands and the discrete
harmonic forms are relevant approximation of the continuous ones. Showing the commutation
of differential operators and canonical interpolation operators is not obvious, besides, because
ofitsintricaterelationwithtopologicalfeaturesofthecomputationaldomain,thisisomorphism
| is not evident | to prove | for a discrete | complex. |     |     |
| -------------- | -------- | -------------- | -------- | --- | --- |
As a consequence, we base our staggered discretization on a preexisting finite element complex
for which the harmonic forms isomorphism is already established: the N´ed´elec-Raviart-Thomas
complex. Notably, staggering is done with the Raviart-Thomas space and the pressure and
densityarepiecewiseconstant. Finally, wehaveexploredthenaturalbyproductsofthediscrete
complex,ofwhich: thedefinitionofdiscretedifferentialadjoints,discretesecondorderoperators
and a discrete HHD with boundary conditions. These fundamental features will be exploited
in the following chapter in order to study the long time consistency of a class of staggered
schemes.

Chapter 5
Development of a class of long time
consistent staggered schemes on the
wave system
Contents
5.1 Introduction . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 105
5.2 The Raviart-Thomas staggered scheme for the multi-dimensional
wave system. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 107
5.2.1 Deriving the ’centred’ finite volume scheme . . . . . . . . . . . . . . . 107
5.2.2 Adding the appropriate, curl-preserving, diffusion . . . . . . . . . . . . 112
5.3 L2 Stability analysis . . . . . . . . . . . . . . . . . . . . . . . . . . 115
−
5.3.1 Tools for stability. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 115
5.3.2 Stability results . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 119
5.3.3 Proof of Proposition 5.3.1 . . . . . . . . . . . . . . . . . . . . . . . . . 121
5.3.4 Proof of Proposition 5.3.3 . . . . . . . . . . . . . . . . . . . . . . . . . 126
5.4 Discrete long time behaviour . . . . . . . . . . . . . . . . . . . . . . 127
5.5 Discussion on some preexisting staggered schemes through low
Mach number asymptotics . . . . . . . . . . . . . . . . . . . . . . . . 136
5.6 Numerical Results . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 138
5.6.1 Numerical long time . . . . . . . . . . . . . . . . . . . . . . . . . . . . 139
5.6.2 On the necessity of a stationary preserving diffusion . . . . . . . . . . 142
5.7 Conclusion . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 147
104

| CHAPTER   |     | 5.  | DEVELOPMENT |     |        | OF   | A CLASS | OF     | LONG | TIME | CONSISTENT |     |     |
| --------- | --- | --- | ----------- | --- | ------ | ---- | ------- | ------ | ---- | ---- | ---------- | --- | --- |
| STAGGERED |     |     | SCHEMES     |     | ON THE | WAVE |         | SYSTEM |      |      |            |     | 105 |
5.1 Introduction
The main focus of this chapter is developing a class of staggered schemes that is long time
| consistent |     | on the | wave | system: |     |     |     |     |     |     |     |     |     |
| ---------- | --- | ------ | ---- | ------- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
1
|     |     |     |     |     | p                 |     | 0   |     | div | p        |          |     |       |
| --- | --- | --- | --- | --- | ----------------- | --- | --- | --- | --- | -------- | -------- | --- | ----- |
|     |     |     |     | ∂   |                   | +   |     | ρ   |     |          | = 0,     |     | (5.1) |
|     |     |     |     | τ   | u                 |    |     | 0   |    | u        |          |     |       |
|     |     |     |     |     |                   |     | κ ∇ |     | 0   |          |          |     |       |
|     |     |     |     |     | (cid:18) (cid:19) |     | 0   |     |     | (cid:18) | (cid:19) |     |       |
|     |     |     |     |     |                   |    |     |     |    |          |          |     |       |
where p is the pressure and u is the velocity, κ ,ρ are parameters such that the wave velocity
|     | κ   |     |     |     |     |     |     | 0 0 |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
0
isc2 = . Thislinearhyperbolicsystemisenrichedwithweaklyenforcedboundaryconditions
0 ρ
0
| of    | inlet/outlet |       | type:     |      |              |     |           |            |     |     |     |       |       |
| ----- | ------------ | ----- | --------- | ---- | ------------ | --- | --------- | ---------- | --- | --- | --- | ----- | ----- |
|       |              |       |           |      |              |     |           | u          | n+u | n   | c   |       |       |
|       |              |       | 1         |      |              |     |           | 1          |     | b   | 0   |       |       |
|       |              |       |           |      |              |     |           | ·          |     | · + | (p  | p )   |       |
|       |              |       |           | u n  |              |     | ρ         |            |     |     |     | b     |       |
|       |              |       | ρ         | ·    |              | =   |           | 0          | 2   |     | 2 − | ,     | (5.2) |
|       |              |       | 0         |      |              |     |          | p+p        |     | c   |     |      |       |
|       |              |       |          |      |             |     |           |            | b   | 0   |     |       |       |
|       |              |       | κ         | 0 pn |              |     | κ 0       |            | n+  | (u  | n u | b n)n |       |
|       |              |       |           |      | inlet/outlet |     |           | 2          |     | 2 · | −   | ·     |       |
|       |              |       |           |      |              |     |          |            |     |     |     |      |       |
|       |              |       |          |      |             |     |          |            |     |     |     |      |       |
| where | p            | and u | are given |      | pressure     | and | velocity, | satisfying |     |     |     |       |       |
|       | b            |       | b         |      |              |     |           |            |     |     |     |       |       |
|       |              |       |           |      |              |     | u         | ndΓ        | = 0 |     |     |       | (5.3) |
b
·
|     |         |       |     |     |     |     | (cid:90) ∂Ω |     |     |     |     |     |     |
| --- | ------- | ----- | --- | --- | --- | --- | ----------- | --- | --- | --- | --- | --- | --- |
| or  | of wall | type: |     |     |     |     |             |     |     |     |     |     |     |
1
|     |     |     |     |     | u    | n   |     |            | 0    |      |          |     |       |
| --- | --- | --- | --- | --- | ---- | --- | --- | ---------- | ---- | ---- | -------- | --- | ----- |
|     |     |     |     |     | ρ 0  | ·   | =   |            |      |      | .        |     | (5.4) |
|     |     |     |     |     |     |    |     | κ          | pn+c | u nn |          |     |       |
|     |     |     |     |     | κ pn |     |     | (cid:20) 0 |      | 0    | (cid:21) |     |       |
|     |     |     |     |     | 0    |     |     |            |      | ·    |          |     |       |
wall
|     |     |     |     |     |    |    |     |     |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
The discussion in chapter 2 put in light several criteria that our staggered discretization
| should | verify | for | a consistent |     | long | time | behaviour: |     |     |     |     |     |     |
| ------ | ------ | --- | ------------ | --- | ---- | ---- | ---------- | --- | --- | --- | --- | --- | --- |
1 Our staggered discretization should enable the identification of the limit with Hodge-
Helmholtz Decomposition, as a discrete equivalent of Theorem 2.3.1: this is done in
|     | Theorem |     | 4.4.4 |     |     |     |     |     |     |     |     |     |     |
| --- | ------- | --- | ----- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
2 Structure preservation, in particular keeping the invariance, at the discrete scale, of the
divergence-free part associated with the discrete HHD of the velocity as in Lemma 2.3.1.
3 Last but not least, to even get a chance to mimick Lemma 2.3.2 at the discrete level we
will need to enforce the relative energy dissipation: at the discrete scale it is not obvious,
|     | especially |     | in explicit |     | time integration. |     |     |     |     |     |     |     |     |
| --- | ---------- | --- | ----------- | --- | ----------------- | --- | --- | --- | --- | --- | --- | --- | --- |
Stability was studied in chapter 3 in one space dimension and the analysis hints towards
three possible time integrations: an explicit time integration with diffusion on both equations,

| CHAPTER   |     | 5.  | DEVELOPMENT |     | OF     | A CLASS | OF LONG | TIME | CONSISTENT |     |
| --------- | --- | --- | ----------- | --- | ------ | ------- | ------- | ---- | ---------- | --- |
| STAGGERED |     |     | SCHEMES     |     | ON THE | WAVE    | SYSTEM  |      |            | 106 |
a (only) pressure centred ImEx time stepping and a fully-centred implicit one. As far as the
theoretical interest goes, the explicit time integration presents the strongest one: finding the
equilibrium between point 3 and point 2 will force us to introduce a new numerical diffusion
operator.
As for point 1 , chapter 4 shows that a discrete Hodge-Helmholtz decomposition is seen as
a byproduct of a discrete de Rham complex so we base our staggered discretization on the
| N´ed´elec-Raviart-Thomas |      |         |              | finite | element | complex.       |     |     |     |     |
| ------------------------ | ---- | ------- | ------------ | ------ | ------- | -------------- | --- | --- | --- | --- |
|                          | This | chapter | is separated |        | in six  | main sections: |     |     |     |     |
1) First,insection5.2,thefiniteelementdeRhamcomplexwhichconstitutesthefoundation
of our staggered scheme is defined and detailed. Then we describe which ingredients are
used to derive our scheme, and in fine how to derive it in order to go from finite element
spaces to a finite volume-like approximation. Then we describe our proposition in order
|     | to  | obtain | simultaneously |     | stability | 3   | and structure | preservation | 2 . |     |
| --- | --- | ------ | -------------- | --- | --------- | --- | ------------- | ------------ | --- | --- |
L2-stability
2) In section 5.3 the of the different time stepping strategies are analyzed in
order to obtain energy dissipation criterion for the study of the long time limit.
3) Next,insection5.4thediscreterelativeenergydissipationisexaminedinordertoidentify
|     | the | long time | limit. |     |     |     |     |     |     |     |
| --- | --- | --------- | ------ | --- | --- | --- | --- | --- | --- | --- |
4) Then in section 5.5, this work is compared with pre-existing staggered schemes.
5) Lastly, in section 5.6 numerical simulations that validate the properties of the scheme are
presented.
|     | 6) The | findings | of  | this chapter | are | summed | up in section | 5.7. |     |     |
| --- | ------ | -------- | --- | ------------ | --- | ------ | ------------- | ---- | --- | --- |

| CHAPTER   | 5.          | DEVELOPMENT    |      |        | OF     | A CLASS   |        | OF LONG |        | TIME | CONSISTENT |        |     |
| --------- | ----------- | -------------- | ---- | ------ | ------ | --------- | ------ | ------- | ------ | ---- | ---------- | ------ | --- |
| STAGGERED | SCHEMES     |                |      | ON THE | WAVE   |           | SYSTEM |         |        |      |            |        | 107 |
| 5.2       | The         | Raviart-Thomas |      |        |        | staggered |        |         | scheme |      | for the    | multi- |     |
|           | dimensional |                | wave |        | system |           |        |         |        |      |            |        |     |
Details on the discretization spaces can be found in chapter 4, we simply recall that we base
our approximation of the velocity and the pressure on the N´ed´elec-Raviart-Thomas complex,
| which reads | ([84, | 99]) | on quadrangular |        |     | and       | triangular | meshes: |         |     |     |     |     |
| ----------- | ----- | ---- | --------------- | ------ | --- | --------- | ---------- | ------- | ------- | --- | --- | --- | --- |
|             |       |      |                 | cG1(Ω) |     | ∇⊥ RT1(Ω) |            | div     | dG0(Ω), |     |     |     |     |
−
|     |     |     |     |     |     | −→  |     | −→  |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
Thus,thediscretevelocityisintheRaviart-Thomasspace(RT1(Ω))andthepressureiscellwise
(dG0(Ω)).
constant
| 5.2.1 | Deriving | the | ’centred’ |     | finite |     | volume | scheme |     |     |     |     |     |
| ----- | -------- | --- | --------- | --- | ------ | --- | ------ | ------ | --- | --- | --- | --- | --- |
Given finite element spaces we can derive the discrete wave system approximation; in our case
we impose the boundary conditions weakly. To do so, the derivation will be Discontinuous
Galerkin inspired [112], in order to impose fluxes with numerical traces for boundary faces.
We recall that the normal associated with boundary faces are oriented in such way that they
are outer normals and denote b , b the boundary faces on which are imposed,
|     |     |     |     | Fwall |     | Finlet/outlet |     |     |     |     |     |     |     |
| --- | --- | --- | --- | ----- | --- | ------------- | --- | --- | --- | --- | --- | --- | --- |
respectively, wall boundary conditions (5.4), and inlet/outlet boundary conditions (5.2).
| Pressure | equation |     |     |     |     |     |     |     |     |     |     |     |     |
| -------- | -------- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
In the derivation, we first suppose that the velocity and the pressure are regular, u (Ω)d
∞
∈ C
and p (Ω). We then multiply the pressure equation by a test function ϕ dG0(Ω) and
∞
| ∈         | C      |         |     |     |     |     |     |     |     |     | ∈   |     |     |
| --------- | ------ | ------- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| integrate | on the | domain: |     |     |     |     |     |     |     |     |     |     |     |
1
|             |           |           |           | ∂         | pϕdx+ |             | div(u)ϕdx   |           | =   | 0,         |         |     | (5.5) |
| ----------- | --------- | --------- | --------- | --------- | ----- | ----------- | ----------- | --------- | --- | ---------- | ------- | --- | ----- |
|             |           |           |           |           | τ     | ρ           |             |           |     |            |         |     |       |
|             |           |           |           | (cid:90)Ω |       |             | 0 (cid:90)Ω |           |     |            |         |     |       |
| integrating | by parts, |           | and using | that      | ϕ     | is cellwise |             | constant: |     |            |         |     |       |
|             | 1         |           |           |           | 1     |             |             |           | 1   |            |         |     |       |
|             |           | div(u)ϕdx |           | =         |       |             | u ∇ϕdx+     |           |     |            | u nϕdΓ, |     | (5.6) |
|             | ρ         |           |           |           | −ρ    |             |             |           | ρ   |            |         |     |       |
|             | 0         | (cid:90)Ω |           |           | 0     |             | K ·         |           | 0   |            | ∂K ·    |     |       |
|             |           |           |           |           |       | K (cid:90)  |             |           |     | K (cid:90) |         |     |       |
|             |           |           |           |           |       | (cid:88)∈C  |             | =0        |     | (cid:88)∈C |         |     |       |
(cid:124)(cid:123)(cid:122)(cid:125)
| which yields | by  | using | (5.6) | in (5.5) |     |     |     |     |     |     |     |     |     |
| ------------ | --- | ----- | ----- | -------- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
1
|     |     |     |     | ∂ τ pϕdx+ |     |     |            | u   | nϕdΓ | = 0. |     |     | (5.7) |
| --- | --- | --- | --- | --------- | --- | --- | ---------- | --- | ---- | ---- | --- | --- | ----- |
|     |     |     |     |           |     | ρ   |            | ·   |      |      |     |     |       |
|     |     |     |     | (cid:90)Ω |     | 0   |            | ∂K  |      |      |     |     |       |
|     |     |     |     |           |     |     | K (cid:90) |     |      |      |     |     |       |
(cid:88)∈C
Now when the unknowns are not regular anymore but u RT1(Ω) and p dG0(Ω), the
|     |     |     |     |     |     |     |     |     |     | ∈   | ∈   |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
normal trace (u n ) is, because of the properties of the Raviart-Thomas finite element,
|     |     | K,σ | σ   |     |     |     |     |     |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
|     |     | ·   | |   |     |     |     |     |     |     |     |     |     |     |
replaced by a constant by face normal trace u n . Then, taking ϕ as the indicator function
· K,σ
of a fixed cell K in (5.7) gives, denoting u σ := u(x σ ) n σ (x σ center of gravity of the face
|     |     | ∈ C |     |     |     |     |     |     | ·   |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |

| CHAPTER   |     | 5. DEVELOPMENT |     |     | OF       | A CLASS | OF     | LONG | TIME | CONSISTENT |     |     |     |
| --------- | --- | -------------- | --- | --- | -------- | ------- | ------ | ---- | ---- | ---------- | --- | --- | --- |
| STAGGERED |     | SCHEMES        |     | ON  | THE WAVE |         | SYSTEM |      |      |            |     |     | 108 |
σ),
1
|     |     |     |     |     | K ∂ p | +   | σ   | ε (σ)u | =   | 0,  |     |     | (5.8) |
| --- | --- | --- | --- | --- | ----- | --- | --- | ------ | --- | --- | --- | --- | ----- |
|     |     |     |     |     | τ K   |     |     | K      | σ   |     |     |     |       |
|     |     |     |     |     | | |   | ρ 0 | | | |        |     |     |     |     |       |
σ ∂K
(cid:88)∈
| because |     | for ϕ = 1 | and | u a | Raviart-Thomas |     | function: |     |     |     |     |     |     |
| ------- | --- | --------- | --- | --- | -------------- | --- | --------- | --- | --- | --- | --- | --- | --- |
K
|     |            | u           | n1  | dΓ = | u ndΓ       | =   | ε         | (σ)        | u n | dΓ = |           | σ ε (σ)u | .   |
| --- | ---------- | ----------- | --- | ---- | ----------- | --- | --------- | ---------- | --- | ---- | --------- | -------- | --- |
|     |            |             | K   |      |             |     | K         |            | σ   |      |           | K        | σ   |
|     |            | ·           |     |      | ·           |     |           |            | ·   |      |           | | |      |     |
|     | L          | (cid:90) ∂L |     |      | (cid:90) ∂K | σ   | ∂K        | (cid:90) σ |     |      | σ ∂K      |          |     |
|     | (cid:88)∈C |             |     |      |             |     | (cid:88)∈ |            |     |      | (cid:88)∈ |          |     |
Since we will aslo have to treat boundary conditions,the normal trace u will be replaced by a
σ
numerical normal trace u , yielding the following updated formulation of (5.8)
σ
1
|     |     |     |     | (cid:99)  | ∂ pϕdx+ |     | σ    | ε (σ)u | =   | 0.  |     |     | (5.9) |
| --- | --- | --- | --- | --------- | ------- | --- | ---- | ------ | --- | --- | --- | --- | ----- |
|     |     |     |     |           | τ       |     |      | K      | σ   |     |     |     |       |
|     |     |     |     |           |         | ρ   | |    | |      |     |     |     |     |       |
|     |     |     |     | (cid:90)Ω |         | 0   | σ ∂K |        |     |     |     |     |       |
(cid:88)∈
(cid:99)
| where | we  | define the | numerical |     | trace u | as: |     |     |     |     |     |     |     |
| ----- | --- | ---------- | --------- | --- | ------- | --- | --- | --- | --- | --- | --- | --- | --- |
σ
Definition 5.2.1 (Velocity flux). The normal velocity flux u is defined
|     |     |     |     |      | (cid:99) |       |      |     | σ             |     |     |     |     |
| --- | --- | --- | --- | ---- | -------- | ----- | ---- | --- | ------------- | --- | --- | --- | --- |
|     |     |     |     |      | u        |       | if σ |     | int           |     |     |     |     |
|     |     |     |     |      |          | σ     |      |     | (cid:99)      |     |     |     |     |
|     |     |     |     |      |          |       | ∈    |     | F             |     |     |     |     |
|     |     |     |     |      | u σ +(u  | b ) σ |      |     |               |     |     |     |     |
|     |     |     |     | u := |         |       | if σ |     | b             |     |     |     |     |
|     |     |     |     | σ    |          |       |      |     | Finlet/outlet |     |     |     |     |
|     |     |     |     |      |  2      |       | ∈    |     |               |     |     |     |     |
|     |     |     |     |      |        |       | if σ |     | b             |     |     |     |     |
0
|     |     |     |     |     |     |     | ∈   |     | Fwall |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | ----- | --- | --- | --- | --- |
(cid:99)

 
| Velocity |     | equation |     |     |     |     |     |     |     |     |     |     |     |
| -------- | --- | -------- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
For the velocity equation, we suppose as previously in a first place that the unknown are in
|     |     |     |     | u   | (Ω)d | p   |     |     |     |     |     |     |     |
| --- | --- | --- | --- | --- | ---- | --- | --- | --- | --- | --- | --- | --- | --- |
fact infinitely regular: ∞ and ∞ (Ω). We then proceed similarly: we multiply by
|     |     |     |     | ∈ C |     | ∈ C |     |     |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
a test function in Ψ RT1(Ω) the velocity equation and integrate over the domain:
∈
|             |     |           |           |       | ∂ u       | Ψdx+κ             | ∇p          | Ψdx | 0.       |        |          |     |        |
| ----------- | --- | --------- | --------- | ----- | --------- | ----------------- | ----------- | --- | -------- | ------ | -------- | --- | ------ |
|             |     |           |           |       | τ         |                   | 0           |     | =        |        |          |     | (5.10) |
|             |     |           |           |       | ·         |                   |             | ·   |          |        |          |     |        |
|             |     |           |           |       | (cid:90)Ω |                   | (cid:90)Ω   |     |          |        |          |     |        |
| Integrating |     | again     | by parts, | we    | get:      |                   |             |     |          |        |          |     |        |
|             |     | κ         | ∇p        | Ψdx=κ |           |                   | p div(Ψ)dx+ |     |          | pΨ ndΓ |          |     |        |
|             |     | 0         |           | ·     | 0         | −                 |             |     |          | ·      |          |     |        |
|             |     | (cid:90)Ω |           |       |           | K                 |             |     | ∂K       |        |          |     |        |
|             |     |           |           |       | K         | (cid:20) (cid:90) |             |     | (cid:90) |        | (cid:21) |     |        |
(cid:88)∈C
|     |     |     |     |     | = κ        | p        | div(Ψ)dx+κ |     |            | [[pΨ        | n   | ]]dΓ |        |
| --- | --- | --- | --- | --- | ---------- | -------- | ---------- | --- | ---------- | ----------- | --- | ---- | ------ |
|     |     |     |     |     | − 0        |          |            |     | 0          |             | ·   | f    | (5.11) |
|     |     |     |     |     |            | K        |            |     |            | f           |     |      |        |
|     |     |     |     |     | K          | (cid:90) |            |     | f          | int(cid:90) |     |      |        |
|     |     |     |     |     | (cid:88)∈C |          |            |     | ∈(cid:88)F |             |     |      |        |
|     |     |     |     |     | +κ         |          | pΨ n       | dΓ, |            |             |     |      |        |
|     |     |     |     |     | 0          |          |            | f   |            |             |     |      |        |
f ·
|     |     |     |     |     |     | f b(cid:90) |     |     |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | ----------- | --- | --- | --- | --- | --- | --- | --- |
(cid:88)∈F

| CHAPTER   | 5.         | DEVELOPMENT |          |                   | OF A | CLASS  | OF    | LONG | TIME | CONSISTENT |     |        |
| --------- | ---------- | ----------- | -------- | ----------------- | ---- | ------ | ----- | ---- | ---- | ---------- | --- | ------ |
| STAGGERED |            | SCHEMES     |          | ON THE            | WAVE | SYSTEM |       |      |      |            |     | 109    |
| Using the | regularity |             | of p and | normal-regularity |      |        | of Ψ, | we   | get  |            |     |        |
|           |            |             |          | κ                 |      | [[pΨ   | n     | ]]dΓ | = 0. |            |     | (5.12) |
|           |            |             |          | 0                 |      |        |       | f    |      |            |     |        |
·
int(cid:90) f
f ∈(cid:88)F
| By combining |     | (5.12)    | and (5.11) | in         | (5.10)     | we           | have |     |            |      |         |        |
| ------------ | --- | --------- | ---------- | ---------- | ---------- | ------------ | ---- | --- | ---------- | ---- | ------- | ------ |
|              |     | ∂ u       | Ψdx        | κ          |            | p div(Ψ)dx+κ |      |     |            | pΨ n | dΓ = 0. | (5.13) |
|              |     | τ         |            | 0          |            |              |      | 0   |            | f    |         |        |
|              |     |           | ·          | −          |            |              |      |     |            | ·    |         |        |
|              |     | (cid:90)Ω |            | K          | (cid:90) K |              |      | f   | b(cid:90)  | f    |         |        |
|              |     |           |            | (cid:88)∈C |            |              |      |     | (cid:88)∈F |      |         |        |
RT1(Ω)
In the case, where u and p are not regular but respectively and dG0(Ω), we modify
| (5.13) as | follow: |           |     |            |          |              |     |     |             |        |         |        |
| --------- | ------- | --------- | --- | ---------- | -------- | ------------ | --- | --- | ----------- | ------ | ------- | ------ |
|           |         | ∂ τ u     | Ψdx | κ 0        |          | p div(Ψ)dx+κ |     | 0   |             | pΨ n f | dΓ = 0, | (5.14) |
|           |         |           | ·   | −          |          |              |     |     |             | ·      |         |        |
|           |         | (cid:90)Ω |     |            | K        |              |     |     |             | f      |         |        |
|           |         |           |     | K          | (cid:90) |              |     |     | f b(cid:90) |        |         |        |
|           |         |           |     | (cid:88)∈C |          |              |     |     | (cid:88)∈F  |        |         |        |
(cid:98)
where we replaced the trace of the regular pressure, p by a constant by face numerical trace
σ
|
| p which | we define | as: |     |     |     |     |     |     |     |     |     |     |
| ------- | --------- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
σ
|     | p Kσ | +p b | c 0 |     |     |     |     |     |     |     |     |     |
| --- | ---- | ---- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
1. p := + ((u ) u ) for inlet/outlet boundary conditions (5.2),
| (cid:98) σ |      |      |     | b σ           | σ   |            |     |       |     |     |     |     |
| ---------- | ---- | ---- | --- | ------------- | --- | ---------- | --- | ----- | --- | --- | --- | --- |
|            |      | 2    | 2   | −             |     |            |     |       |     |     |     |     |
| 2. p       | := p | +c u | for | wall boundary |     | conditions |     | (5.4) | .   |     |     |     |
| (cid:98)σ  | Kσ   | 0    | σ   |               |     |            |     |       |     |     |     |     |
(cid:98)
dG0(Ω),
From now on p so taking as test function Ψ = Ψ σ the basis function associated
∈
σ
| with a fixed | face | in  | (5.14),    | we have  |     |      |     |            |          |             |     |        |
| ------------ | ---- | --- | ---------- | -------- | --- | ---- | --- | ---------- | -------- | ----------- | --- | ------ |
|              |      |     | κ          | pdiv(Ψ   | σ   | )dx= | κ   |            | p K      | div(Ψ σ )dx |     |        |
|              |      | −   | 0          |          |     |      | − 0 |            |          |             |     |        |
|              |      |     |            | K        |     |      |     |            | K        |             |     |        |
|              |      |     | K          | (cid:90) |     |      | K   | (σ)        | (cid:90) |             |     |        |
|              |      |     | (cid:88)∈C |          |     |      |     | (cid:88)∈C |          |             |     |        |
|              |      |     |            |          |     | =    | κ   |            | p        | Ψ n         | dΓ  |        |
|              |      |     |            |          |     |      | 0   |            | K        | σ f         |     | (5.15) |
|              |      |     |            |          |     |      | −   |            | ∂K       | ·           |     |        |
|              |      |     |            |          |     |      | K   | (σ)        | (cid:90) |             |     |        |
(cid:88)∈C
|     |     |     |     |     |     | =   | κ σ   |     | ε K | (σ)p K . |     |     |
| --- | --- | --- | --- | --- | --- | --- | ----- | --- | --- | -------- | --- | --- |
|     |     |     |     |     |     |     | − 0 | | |   |     |          |     |     |
K (σ)
(cid:88)∈C
| Plugging | the | expression | (5.15) | in (5.14) |              | yields |      |     |              |      |       |        |
| -------- | --- | ---------- | ------ | --------- | ------------ | ------ | ---- | --- | ------------ | ---- | ----- | ------ |
|          |     | ∂ u        | dx     | κ σ       |              | ε      | (σ)p | +κ  |              | pΨ n | dΓ    |        |
|          |     | τ          | Ψ σ    | 0         |              | K      | K    | 0   |              | σ    | f = 0 | (5.16) |
|          |     |            | ·      | − |       | |            |        |      |     |              | ·    |       |        |
|          |     | (cid:90)Ω  |        |           |              |        |      |     | b(cid:90)    | f    |       |        |
|          |     |            |        |           | K (cid:88)∈C | (σ)    |      |     | f (cid:88)∈F |      |       |        |
(cid:98)

| CHAPTER      | 5.  | DEVELOPMENT |     |       | OF     | A CLASS |        | OF LONG | TIME | CONSISTENT |     |
| ------------ | --- | ----------- | --- | ----- | ------ | ------- | ------ | ------- | ---- | ---------- | --- |
| STAGGERED    |     | SCHEMES     |     | ON    | THE    | WAVE    | SYSTEM |         |      |            | 110 |
| Mass-lumping |     | Finally,    |     | since | u(x,τ) | =       | u (τ)Ψ | (x),    |      |            |     |
|              |     |             |     |       |        |         | f      | f       |      |            |     |
f (cid:88)∈F
|     |     |     | σ   |     | ∂         | u Ψ | dx = | ∂ u | Ψ         | Ψ dx | (5.17) |
| --- | --- | --- | --- | --- | --------- | --- | ---- | --- | --------- | ---- | ------ |
|     |     |     |     |     |           | τ   | σ    | τ   | f         | f σ  |        |
|     |     |     | ∀   | ∈ F |           | ·   |      |     |           | ·    |        |
|     |     |     |     |     | (cid:90)Ω |     |      | f   | (cid:90)Ω |      |        |
(cid:88)∈F
M
| Then denoting |     | the | square | matrix |     | of size | # such | that |     |     |     |
| ------------- | --- | --- | ------ | ------ | --- | ------- | ------ | ---- | --- | --- | --- |
F
|     |     |     |     |     | σ,f |     | M = | Ψ   | Ψ dx |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | ---- | --- | --- |
|     |     |     |     |     |     |     | σ,f | f   | σ    |     |     |
|     |     |     |     | ∀   | ∈   | F   |     |     | ·    |     |     |
(cid:90)Ω
R#
| and U := | (u ) |     | , (5.17) |     | is equivalent |     | to      |     |     |     |     |
| -------- | ---- | --- | -------- | --- | ------------- | --- | ------- | --- | --- | --- | --- |
|          | f    | f ∈ | F        |     |               |     |         |     |     |     |     |
|          |      |     |          |     |               |     | ∂ [MU], |     |     |     |     |
τ
withManon-diagonalmatrix. Thisisnotreconciliablewithouraimtodevelopafinitevolume
scheme. Wecancircumventthisdifficultywithamass-lumping; thegeneralideaistotransform
the matrix M into a diagonal matrix: in [52] this procedure is applied for nodal finite elements
M˜
by replacing the mass matrix M by a diagonal matrix for which the diagonal entries are
the sum of all the terms on each rows. The authors of [52] are able to show that the obtained
measure can be associated to a geometric control volume around a node. In [53, 113], the mass-
lumping is defined by approximating the (L2)d scalar product using the linear forms (4.18) as
integrationnodes. Theweightsofintegrationassociatedaredeterminedbysearchingadiagonal
|     | M˜  |     |     |     |     |     |     |     |     | (L2)d |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | ----- | --- |
matrix such that the bilinear form associated is consistent with the original scalar
product for constant fields. Also, in [54], similar results are recovered with a Petrov-Galerkin
approach.
In our case the mass-lumping procedure will mimic the approach from [52]:
Definition 5.2.2 (Mass-lumped scalar product). Denote M the mass matrix of the Raviart-
Thomas space. We define the mass-lumping by replacing the initial mass matrix M by a lumped
M˜
| matrix | :   |     |     |     |     |     |     |     |     |       |     |
| ------ | --- | --- | --- | --- | --- | --- | --- | --- | --- | ----- | --- |
|        |     |     |     |     |     |     |     | M   | if  | σ = f |     |
σ,j
M˜
|     |     |     | σ,f |     |     | σ,f := | j            |     |      | .   |     |
| --- | --- | --- | --- | --- | --- | ------ | ------------ | --- | ---- | --- | --- |
|     |     |     | ∀   | ∈ F |     |        |  (cid:88)∈F |     |      |     |     |
|     |     |     |     |     |     |        |              | 0   | else |     |     |

| So that | the lumped | scalar |     | product | .,. | verifies: |     |     |     |     |     |
| ------- | ---------- | ------ | --- | ------- | --- | --------- | --- | --- | --- | --- | --- |
|         |            |        |     |         |     | h         |    |     |     |     |     |
(cid:104) (cid:105)
RT1(Ω),
|     |     |     | u   |     |     |     | σ   | u,Ψ σ               | h := | D σ u σ , |     |
| --- | --- | --- | --- | --- | --- | --- | --- | ------------------- | ---- | --------- | --- |
|     |     |     | ∀   | ∈   |     | ∀   | ∈ F | (cid:104) (cid:105) | |    | |         |     |
with
|     |     |     |     |     | D   | :=  | Ψ   | Ψ dx. |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | ----- | --- | --- | --- |
|     |     |     |     |     | |   | σ | |     | f · σ |     |     |     |
(cid:90)Ω
f
(cid:88)∈F

| CHAPTER   | 5.  | DEVELOPMENT |     |     | OF  | A CLASS |        | OF LONG | TIME | CONSISTENT |     |
| --------- | --- | ----------- | --- | --- | --- | ------- | ------ | ------- | ---- | ---------- | --- |
| STAGGERED |     | SCHEMES     |     | ON  | THE | WAVE    | SYSTEM |         |      |            | 111 |
Gathering the definition of the mass-lumping and the definition of the pressure trace p (1),
σ
(2) in (5.16) yields the following discrete ’centred’ staggered velocity equation
(cid:98)
c
0
|          |           |             |     | D ∂ | u +κ | σ [[p]] | +   | σ [[u | n]] | = 0, | (5.18) |
| -------- | --------- | ----------- | --- | --- | ---- | ------- | --- | ----- | --- | ---- | ------ |
|          |           |             |     | σ   | τ σ  | 0       | σ   |       |     | σ    |        |
|          |           |             |     | | | |      | | |     |     | 2 | | | ·   |      |        |
| with the | following | definitions |     |     |      |         |     |       |     |      |        |
Definition 5.2.3 (Discrete pressure jump ). We define the discrete pressure gradient for p
| dG0(Ω) | [[p]] | as follow: |     |     |     |     |     |     |     |     |     |
| ------ | ----- | ---------- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
σ
∈
|     |     |       |      |            |       | ε (σ)p |     | if σ |     | int            |     |
| --- | --- | ----- | ---- | ---------- | ----- | ------ | --- | ---- | --- | -------------- | --- |
|     |     |       |      |            | K (σ) | K      | K   |      |     |                |     |
|     |     |       |      |            | ∈pC   | p      |     | ∈    |     | F              |     |
|     |     |       |      |            | b     | Kσ     |     |      |     | b              |     |
|     |     | [[p]] | σ := |  (cid:80) |       | −      |     | if σ |     |                |     |
|     |     |       |      |            |       | 2      |     | ∈    |     | Fi nlet/outlet |     |

|     |     |     |     |    |     | 0   |     | if σ |     | b   |     |
| --- | --- | --- | --- | --- | --- | --- | --- | ---- | --- | --- | --- |
Fwall
∈

Definition 5.2.4 (Discret  e normal velocity jump ). We define the discrete normal velocity
| jump [[u | n]] | for u | RT1(Ω) |     | as follow: |     |     |     |     |     |     |
| -------- | --- | ----- | ------ | --- | ---------- | --- | --- | --- | --- | --- | --- |
σ
· ∈
|     |     |     |         |     |     | 0   | if  | σ   |     | int            |     |
| --- | --- | --- | ------- | --- | --- | --- | --- | --- | --- | -------------- | --- |
|     |     |     |         |     |     |     |     | ∈   |     | F              |     |
|     |     |     | [[u n]] | :=  | (u  | ) u | if  | σ   |     | b              |     |
|     |     |     |         | σ   |    | b σ | σ   |     |     | Fi nlet/outlet |     |
|     |     |     | ·       |     |     | −   |     | ∈   |     |                |     |
|     |     |     |         |     |    | 2u  | if  | σ   |     | b              |     |
|     |     |     |         |     |    | σ   |     |     |     | Fw all         |     |
∈
 
| Preservation |     | of a | discrete | curl | on  | the | centred | scheme |     |     |     |
| ------------ | --- | ---- | -------- | ---- | --- | --- | ------- | ------ | --- | --- | --- |
OnthetorusΩ = T2 thecontinuouswavesystem(5.1)inducesthepreservationofthefollowing
| differential | operator |     |     |         |     |     |       |        |     |     |        |
| ------------ | -------- | --- | --- | ------- | --- | --- | ----- | ------ | --- | --- | ------ |
|              |          |     |     |         |     |     |       |        | uy  | ux, |        |
|              |          |     | ∂   | τ curlu | h = | 0   | curlu | h := ∂ | x   | ∂ y | (5.19) |
−
which is fondamentally linked to the preservation of a divergence-free part of the velocity
Lemma 2.3.1. In fact, the centred velocity equation (5.18) induces the conservation at the
| discrete | scale of | (5.19) | thanks | to  | the | formalism | of  | complexes: |     |     |     |
| -------- | -------- | ------ | ------ | --- | --- | --------- | --- | ---------- | --- | --- | --- |
Lemma 5.2.1 (Curl preservation ). Suppose that Ω = T2 and let (∇ ) the discrete curl
⊥ ∗
definedinDefinition4.4.2withthescalarproduct ., . givenbythemass-lumpedscalarproduct
h
|     |     |     |     |     |     |     | (cid:104) | (cid:105) |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --------- | --------- | --- | --- | --- |
Definition 5.2.2. Then, the discrete velocity equation (5.16) yields the preservation of this
| discrete | curl |     |     |     |     |        |     |           |     |     |     |
| -------- | ---- | --- | --- | --- | --- | ------ | --- | --------- | --- | --- | --- |
|          |      |     |     |     | ∂   | (∇ ) u | = 0 | in cG1(Ω) |     |     |     |
|          |      |     |     |     | τ   | ⊥ ∗    |     |           |     |     |     |
Proof. The derivation of the scheme shows in particular that (5.18) is equivalent on the bidi-
| mensional | torus | to  |     |           |      |             |            |           |     |     |        |
| --------- | ----- | --- | --- | --------- | ---- | ----------- | ---------- | --------- | --- | --- | ------ |
|           |       |     |     | ∂         | u, Ψ |             |            | pdiv(Ψ)dx | =   | 0,  | (5.20) |
|           |       |     |     |           | τ    | h           |            |           |     |     |        |
|           |       |     |     | (cid:104) |      | (cid:105) − |            |           |     |     |        |
|           |       |     |     |           |      | K           | (cid:90) K |           |     |     |        |
(cid:88)∈C

| CHAPTER   | 5.  | DEVELOPMENT |     |        | OF        | A CLASS | OF     | LONG | TIME | CONSISTENT |     |     |
| --------- | --- | ----------- | --- | ------ | --------- | ------- | ------ | ---- | ---- | ---------- | --- | --- |
| STAGGERED |     | SCHEMES     |     | ON THE | WAVE      |         | SYSTEM |      |      |            |     | 112 |
| taking Ψ  | = ∇ | ϕ for       | ϕ   | cG1(Ω) | in (5.20) | leads   | to     |      |      |            |     |     |
⊥
∈
|     |     |     |           | ∂ u, ∇ | ϕ         |     | pdiv(∇ |     | ϕ)dx | 0.  |     |        |
| --- | --- | --- | --------- | ------ | --------- | --- | ------ | --- | ---- | --- | --- | ------ |
|     |     |     |           | τ      | ⊥         | h   |        |     | ⊥    | =   |     | (5.21) |
|     |     |     | (cid:104) |        | (cid:105) | −   |        |     |      |     |     |        |
K
|     |     |     |     |     |     | K (cid:88)∈C | (cid:90)  | =0                 |           |     |     |     |
| --- | --- | --- | --- | --- | --- | ------------ | --------- | ------------------ | --------- | --- | --- | --- |
|     |     |     |     |     |     |              | (cid:124) | (cid:123)(cid:122) | (cid:125) |     |     |     |
Now, we note that (5.21) is equivalent, since u(x,τ) = u (τ)Ψ (x), to
|     |     |     |     |     |     |     |     |     | σ   | σ   |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
σ
(cid:88)∈F
|     |     |     |     |     |     | ∂ u Ψ     | , ∇ | ϕ =       | 0,  |     |     | (5.22) |
| --- | --- | --- | --- | --- | --- | --------- | --- | --------- | --- | --- | --- | ------ |
|     |     |     |     |     |     | τ σ       | σ ⊥ | h         |     |     |     |        |
|     |     |     |     |     |     | (cid:104) |     | (cid:105) |     |     |     |        |
σ
(cid:88)∈F
since only the d.o.fs are time dependent. Using Definition 4.4.2 this yields
|     |     |     | ϕ   | cG1(Ω) |     | ∂   | u (∇      | ) Ψ | ,ϕ        | = 0, |     |     |
| --- | --- | --- | --- | ------ | --- | --- | --------- | --- | --------- | ---- | --- | --- |
|     |     |     |     |        |     |     | τ σ       | ⊥ ∗ | σ L2(Ω)   |      |     |     |
|     |     |     | ∀   | ∈      |     |     | (cid:104) |     | (cid:105) |      |     |     |
σ
(cid:88)∈F
| so finally | by linearity |     |     |     |     |     |     |     |     |     |     |     |
| ---------- | ------------ | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
cG1(Ω)
|                 |     |     |       | ϕ   |     | ∂           | (∇ ⊥ )    | ∗ u ,ϕ    |                 | = 0, |     |     |
| --------------- | --- | --- | ----- | --- | --- | ----------- | --------- | --------- | --------------- | ---- | --- | --- |
|                 |     |     | ∀     | ∈   |     | (cid:104) τ |           |           | (cid:105) L2(Ω) |      |     |     |
|                 |     |     |       |     |     |             | (cid:104) | (cid:105) |                 |      |     |     |
| which concludes |     | the | proof |     |     |             |           |           |                 |      |     |     |
Remark 5.2.1. This is known that the existence of non-trivial stationary states is equivalent
for some systems to the preservation of differential operators [114, Theorem 2.3], which in our
case is the curl of the velocity. Preserving non-trivial stationary states (e.g for non-uniform
divergence-free fields) will prove to be of paramount importance for the long time convergence.
| 5.2.2 | Adding | the | appropriate, |     |     | curl-preserving, |     |     | diffusion |     |     |     |
| ----- | ------ | --- | ------------ | --- | --- | ---------------- | --- | --- | --------- | --- | --- | --- |
The conclusions established in the chapter 3 show that a stabilization operator is needed on the
velocity equation if we hope for energy dissipation with explicit time stepping: in the following
we develop numerical diffusion operators that are appropriate to maintaining the preservation
| of the curl | Lemma    | 5.2.1. |     |     |     |     |     |     |     |     |     |     |
| ----------- | -------- | ------ | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Pressure    | equation |        |     |     |     |     |     |     |     |     |     |     |
The stabilization of pressure equation does not cause any trouble: in order to both preserve
pressures with null gradient and dissipate energy, it is sufficient to select a second order dis-
sipative operator that preserves constant pressure. This yields the choice of the classical finite
volume Laplacian and therefore, the upwind version of (5.9) would be:
1
|     |     | K   | ∂ p | +   |           | σ ε (σ)u | =   | d c | σ         | ε (σ)[[p]] | ,   |     |
| --- | --- | --- | --- | --- | --------- | -------- | --- | --- | --------- | ---------- | --- | --- |
|     |     |     | τ K |     |           | K        | σ   | 1 0 |           | K          | σ   |     |
|     |     | |   | |   | ρ   |           | | |      |     |     | |         | |          |     |     |
|     |     |     |     | 0   | σ ∂K      |          |     | σ   | ∂K        |            |     |     |
|     |     |     |     |     | (cid:88)∈ |          |     |     | (cid:88)∈ |            |     |     |
(cid:99)

| CHAPTER   | 5.  | DEVELOPMENT |     |        | OF   | A CLASS | OF  | LONG | TIME | CONSISTENT |     |
| --------- | --- | ----------- | --- | ------ | ---- | ------- | --- | ---- | ---- | ---------- | --- |
| STAGGERED |     | SCHEMES     |     | ON THE | WAVE | SYSTEM  |     |      |      |            | 113 |
with d R+ and u and [[p]] given by, respectively, Definition 5.2.1 and Definition 5.2.3.
| 1   |     | σ   |     | σ   |     |     |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
∈
∗
(cid:99)
| Velocity | equation |     |     |     |     |     |     |     |     |     |     |
| -------- | -------- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
Inthecaseofthevelocity, thechoiceisalittlemoretricky. Lemma2.3.1motivatestheneedfor
our scheme to be able to preserve divergence free velocities. As a matter of fact, a stablization
obtained with a discrete Laplacian is very natural and would lead to the formal equation
|     |     |     |     |     | ∂ u+κ | ∇p  | =   | c h∆u |     |     |     |
| --- | --- | --- | --- | --- | ----- | --- | --- | ----- | --- | --- | --- |
|     |     |     |     |     | τ     | 0   |     | 0     |     |     |     |
wherec2 := κ0 andhisalengthrelatedtothemesh. Howeverin2dimensions,thisisequivalent,
| 0             | ρ0        |     |        |           |     |               |     |     |         |          |        |
| ------------- | --------- | --- | ------ | --------- | --- | ------------- | --- | --- | ------- | -------- | ------ |
| for a regular | solution, |     | to     |           |     |               |     |     |         |          |        |
|               |           |     | ∂      | u+κ       | ∇p  | c h ∇div(u)+∇ |     |     | curl(u) | ,        |        |
|               |           |     | τ      | 0         | =   | 0             |     |     | ⊥       |          | (5.23) |
|               |           |     |        |           |     | (cid:18)      |     |     |         | (cid:19) |        |
| applying      | the curl  | to  | (5.23) | we obtain |     |               |     |     |         |          |        |
curl(∂ u)+κ curl(∇p) = c h curl(∇div(u))+curl ∇ curl(u) , (5.24)
|     |     | τ   | 0   |     | 0   |     |     |     |     | ⊥   |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
(cid:18) (cid:19)
(cid:16) (cid:17)
commuting curl(∂ u) = ∂ (curlu) and curl(∇) = 0 implies together that (5.24) becomes
|     |     | τ   | τ   |           |     |       |          |           |     |          |        |
| --- | --- | --- | --- | --------- | --- | ----- | -------- | --------- | --- | -------- | ------ |
|     |     |     |     | ∂ (curlu) |     | = c h | curl     | ∇ curl(u) |     | ,        | (5.25) |
|     |     |     |     | τ         |     | 0     |          | ⊥         |     |          |        |
|     |     |     |     |           |     |       | (cid:16) |           |     | (cid:17) |        |
which differs from (5.19). In order to recover (5.19), we can suppress ∇ curl(u) from ∆u,
⊥
yielding a diffusion in ∇div(u). The discrete equivalent of ∇div(u) is in fact already defined
in chapter 4 as the composition of the divergence div and the discrete gradient ( div) of
− ∗
Definition 4.4.1 obtained with the lumped-scalar product ., . h of Definition 5.2.2:
|     |     |     |     |     |             |        |       |             | (cid:104) | (cid:105) |     |
| --- | --- | --- | --- | --- | ----------- | ------ | ----- | ----------- | --------- | --------- | --- |
|     |     |     |     |     | (           | div) ∗ | divu, | Ψ σ         |           |           |     |
|     |     |     |     |     | (cid:104) − |        |       | (cid:105) h |           |           |     |
Then defining a diffusion analogous to h∇div(u) can be done by first introducing the operator
|     |     | div | : RT1(Ω) |     |     |     |     | dG0(Ω) |     |     |     |
| --- | --- | --- | -------- | --- | --- | --- | --- | ------ | --- | --- | --- |
−→
(5.26)
|     |     | (cid:102) |     |     |     | (cid:93) |     | 1   |     |            |     |
| --- | --- | --------- | --- | --- | --- | -------- | --- | --- | --- | ---------- | --- |
|     |     |           |     | u   |     | (divu)   |     | :=  |     | σ ε (σ)u . |     |
|     |     |           |     |     |     |          | K   | ∂K  |     | K σ        |     |
|     |     |           |     |     | −→  |          |     |     |     | | |        |     |
|     |     |           |     |     |     |          |     | |   | | σ | ∂K         |     |
(cid:88)∈
The image of a Raviart-Thomas function by (5.26) is cellwise constant, enabling the definition
| of the following |     | composition |     |     |     |     |     |     |     |     |     |
| ---------------- | --- | ----------- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
(cid:93)
|     |     |     |     | Ψ   | RT1(Ω) | c   | ( div)      | d   | ivu, | Ψ ,       |     |
| --- | --- | --- | --- | --- | ------ | --- | ----------- | --- | ---- | --------- | --- |
|     |     |     |     |     |        | 0   |             | ∗   |      | h         |     |
|     |     |     |     | ∀ ∈ |        |     | (cid:104) − |     |      | (cid:105) |     |

| CHAPTER   | 5.  | DEVELOPMENT |     |     | OF  | A CLASS |        | OF  | LONG | TIME | CONSISTENT |     |     |
| --------- | --- | ----------- | --- | --- | --- | ------- | ------ | --- | ---- | ---- | ---------- | --- | --- |
| STAGGERED |     | SCHEMES     |     | ON  | THE | WAVE    | SYSTEM |     |      |      |            |     | 114 |
K
which is formally equivalent to (5.23) with h = h = | | . Then, we define
K
∂K
|     |     |     |     |     |     |     |     |     | |   | |   |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
Definition 5.2.5 (Grad-Div diffusion operator). Let ( div) the discrete gradient defined
∗
−
in Definition 4.4.1 with the scalar product ., . given by the mass-lumped scalar product
|            |        |     |          |     |           |           | (cid:104) | (cid:105) h |         |     |     |     |     |
| ---------- | ------ | --- | -------- | --- | --------- | --------- | --------- | ----------- | ------- | --- | --- | --- | --- |
| Definition | 5.2.2. | The | grad-div |     | numerical | diffusion |           | is          | defined | as: |     |     |     |
(cid:93)
|              |     |     |     | Ψ     | RT1(Ω)   |            | (           | div) | d ivu, | Ψ .       |     |     | (5.27) |
| ------------ | --- | --- | --- | ----- | -------- | ---------- | ----------- | ---- | ------ | --------- | --- | --- | ------ |
|              |     |     |     |       |          |            |             |      | ∗      | h         |     |     |        |
|              |     |     |     | ∀     | ∈        |            | (cid:104) − |      |        | (cid:105) |     |     |        |
| It is equal, | for | Ψ = | Ψ a | basis | function | associated |             | to   | a face | σ, to     |     |     |        |
σ
|     |     |             | (cid:93) |     |             |              | (cid:94) |     |      |       | (cid:94)  |     |     |
| --- | --- | ----------- | -------- | --- | ----------- | ------------ | -------- | --- | ---- | ----- | --------- | --- | --- |
|     |     | ( div)      | ∗ divu,  | Ψ   | σ =         | σ [[div(u)]] |          | σ = | σ    | ε K   | (σ)div(u) | .   |     |
|     |     | (cid:104) − |          |     | (cid:105) h | | |          |          |     | −| | |       |           | K   |     |
|     |     |             |          |     |             |              |          |     |      | K (σ) |           |     |     |
(cid:88)∈C
Adding (5.27) to (5.18) we obtain the discrete velocity equation with diffusion:
c
|     |     |     |     |     |         |     | 0     |     |     | (cid:93)  |       |     |        |
| --- | --- | --- | --- | --- | ------- | --- | ----- | --- | --- | --------- | ----- | --- | ------ |
|     |     |     | D ∂ | u + | σ [[p]] | +   | σ [[u | n]] | =   | d c σ [[d | ivu]] | ,   | (5.28) |
|     |     |     | σ   | τ σ |         | σ   | 2     |     | σ   | 2 0       | σ     |     |        |
|     |     |     | | | |     | | |     |     | | |   | ·   |     | | |       |       |     |        |
R+
| d 2 | a diffusion |     | coefficient. |     |     |     |     |     |     |     |     |     |     |
| --- | ----------- | --- | ------------ | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
∈
∗
| Preservation |     | of a | discrete | curl | on  | the | ’upwind’ |     | scheme |     |     |     |     |
| ------------ | --- | ---- | -------- | ---- | --- | --- | -------- | --- | ------ | --- | --- | --- | --- |
With this new diffusion we can replicate Lemma 5.2.1 on the ’upwind’ velocity equation (5.28).
T2
Lemma 5.2.2 (Curl preservation diffusion operator). Suppose that Ω = and let (∇ )
⊥ ∗
the discrete curl defined in Definition 4.4.2 with the scalar product ., . given by the mass-
|     |     |     |     |     |     |     |     |     |     |     | (cid:104) (cid:105) h |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --------------------- | --- | --- |
lumped scalar product Definition 5.2.2. Then, the discrete velocity equation (5.28) yields the
| preservation | of  | this | discrete | curl |     |      |     |             |     |     |     |     |     |
| ------------ | --- | ---- | -------- | ---- | --- | ---- | --- | ----------- | --- | --- | --- | --- | --- |
|              |     |      |          |      | ∂   | (∇ ) | u = | 0 in cG1(Ω) |     |     |     |     |     |
|              |     |      |          |      | τ   | ⊥    | ∗   |             |     |     |     |     |     |
RT1(Ω)
| Proof. | (5.28) is | equivalent |     | for Ψ |     |     | on the | torus | to  |     |     |     |     |
| ------ | --------- | ---------- | --- | ----- | --- | --- | ------ | ----- | --- | --- | --- | --- | --- |
∈
(cid:93)
|     |     |     | ∂ u,      | Ψ           |     | pdiv(Ψ)dx |     | =   | ( div)      | divu, | Ψ         | ,   | (5.29) |
| --- | --- | --- | --------- | ----------- | --- | --------- | --- | --- | ----------- | ----- | --------- | --- | ------ |
|     |     |     | τ         | h           |     |           |     |     |             | ∗     | h         |     |        |
|     |     |     | (cid:104) | (cid:105) − |     | K         |     |     | (cid:104) − |       | (cid:105) |     |        |
K (cid:90)
(cid:88)∈C
| By Definition |     | 4.4.1 | (5.29) | is equivalent |     | to  |     |     |     |     |     |     |     |
| ------------- | --- | ----- | ------ | ------------- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
(cid:93)
|     |     | ∂   | u, Ψ |     |     | pdiv(Ψ)dx |     | =   |     | divu div(Ψ)dx. |     |     | (5.30) |
| --- | --- | --- | ---- | --- | --- | --------- | --- | --- | --- | -------------- | --- | --- | ------ |
τ h
|             |     | (cid:104) |          | (cid:105) − |            |      |        | −       |            |            |     |     |     |
| ----------- | --- | --------- | -------- | ----------- | ---------- | ---- | ------ | ------- | ---------- | ---------- | --- | --- | --- |
|             |     |           |          | K           | (cid:90) K |      |        |         | K (cid:90) | K          |     |     |     |
|             |     |           |          | (cid:88)∈C  |            |      |        |         | (cid:88)∈C |            |     |     |     |
| Then taking | as  | test      | function | Ψ           | = ∇        | ϕ in | (5.30) | implies | in         | particular |     |     |     |
⊥
|     |     |     |     |     | cG1(Ω) |     |     | (cid:93) |       |      |     |     |     |
| --- | --- | --- | --- | --- | ------ | --- | --- | -------- | ----- | ---- | --- | --- | --- |
|     |     |     |     | ϕ   |        |     |     | d ivu    | div(∇ | ϕ) = | 0.  |     |     |
|     |     |     |     | ∀ ∈ |        |     |     |          |       | ⊥    |     |     |     |
K
|     |     |     |     |     |     | K   | (cid:90) |     |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | -------- | --- | --- | --- | --- | --- | --- |
(cid:88)∈C

| CHAPTER   |      | 5.      | DEVELOPMENT |           |     | OF    | A CLASS |        | OF LONG |     | TIME CONSISTENT |     |     |
| --------- | ---- | ------- | ----------- | --------- | --- | ----- | ------- | ------ | ------- | --- | --------------- | --- | --- |
| STAGGERED |      |         | SCHEMES     |           | ON  | THE   | WAVE    | SYSTEM |         |     |                 |     | 115 |
| The       | rest | follows | by          | the proof | of  | Lemma | 5.2.1.  |        |         |     |                 |     |     |
Remark 5.2.2. Notice that apart from the mass-lumping, the derivation of the scheme does
not use directly mesh shape properties so that the scheme is written in the same manner in
triangles and quadrilaterals. Similarly in 3 space dimensions, the only term that depends on
| the | geometry | is  | the | mass lumping. |     |     |     |     |     |     |     |     |     |
| --- | -------- | --- | --- | ------------- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
Lemma 5.2.1 and Lemma 5.2.2 would however state the preservation of a different operator in
| 3D; | the | 3D rotational. |     |     |     |     |     |     |     |     |     |     |     |
| --- | --- | -------------- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
−
Remark 5.2.3. We see here that to conserve the divergence-free part of the velocity in the
scheme we have to insert diffusion terms with stationarity preserving features. In the case of a
bi-dimensional torus, we ended up preserving the curl(2d) of the velocity which is the involution
| preserved |     | in [114]. |           |     |          |     |     |     |     |     |     |     |     |
| --------- | --- | --------- | --------- | --- | -------- | --- | --- | --- | --- | --- | --- | --- | --- |
| 5.3       | L2  |           | Stability |     | analysis |     |     |     |     |     |     |     |     |
−
In this section we wish to obtain an energy dissipation criterion for the scheme derived previ-
| ously, | which | we  | can | write | in the | following | manner |     |     |     |     |     |     |
| ------ | ----- | --- | --- | ----- | ------ | --------- | ------ | --- | --- | --- | --- | --- | --- |
1
|     |     |     |     | K ∂ p |     |                | σ ε | (σ)u     | d       | c              | σ ε (σ)[[p]]    | ,   |        |
| --- | --- | --- | --- | ----- | --- | -------------- | --- | -------- | ------- | -------------- | --------------- | --- | ------ |
|     |     |     |     | τ     | K + |                | K   | σ        | = 1     | 0              | K               | σ   |        |
|     |     |     | |   | |     | ρ   |                | | | |          |         |                | | |             |     |        |
|     |     |     |    |       |     | 0              |     |          |         |                |                 |     |        |
|     |     |     |     |       |     | σ (cid:88)∈ ∂K |     |          |         | σ (cid:88)∈ ∂K |                 |     |        |
|     |     |     |    |       |     |                |     |          |         |                |                 |     | (5.31) |
|     |     |     |    |       |     |                |     | (cid:99) |         |                |                 |     |        |
|     |     |     |   |       |     |                |     | c        |         |                |                 |     |        |
|     |     |     |    |       |     |                |     | 0        |         |                | (cid:94)        |     |        |
|     |     |     | D   | ∂ u   | +κ  | σ [[p]]        | =   | σ        | [[u n]] | +d             | c σ [[d iv(u)]] | ,   |        |
|     |     |     | |   | σ | τ | σ 0 | | |            | σ − | 2 |      | | ·     | σ              | 2 0 | |         | σ   |        |

 

|     |     |     |    | (cid:94) |     |     |     |     |     |     |     |     |     |
| --- | --- | --- | --- | -------- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
with u σ , σ [[p]] σ , σ [[div(u)]] σ and σ [[u n]] σ are given by, respectively, Definition 5.2.1,
|            |     | | |    | |          | |   |        |            | | | | ·      |      |     |     |     |     |
| ---------- | --- | ------ | ---------- | --- | ------ | ---------- | --- | ------ | ---- | --- | --- | --- | --- |
| Definition |     | 5.2.3, | Definition |     | 5.2.5, | Definition |     | 5.2.4. | Also |     |     |     |     |
• (cid:99)
|     | d   | = d | = 0 in | Euler | implicit | integration, |              |     |     |     |     |     |     |
| --- | --- | --- | ------ | ----- | -------- | ------------ | ------------ | --- | --- | --- | --- | --- | --- |
|     | 1   | 2   |        |       |          |              |              |     |     |     |     |     |     |
|     | • d | d   | d      |       |          |              |              |     |     |     |     |     |     |
|     | 1   | = 2 | = in   | Euler | explicit | time         | integration, |     |     |     |     |     |     |
• and finally d = d,d = 0 and the pressure gradient is implicited in semi-implicit time
|     |     |     | 1   |     | 2   |     |     |     |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
integration.
WefocusfirstonEulerexplicitintegrationandimplicitintegration,whereaswecouldexpectit,
fortheformerweobtainthestabilityunderaCFLconditionandforthelatteranunconditional
stability.
| 5.3.1 |      | Tools   | for       | stability   |         |     |            |     |     |        |       |        |     |
| ----- | ---- | ------- | --------- | ----------- | ------- | --- | ---------- | --- | --- | ------ | ----- | ------ | --- |
| We    | will | use the | following |             | lemmas: |     |            |     |     |        |       |        |     |
|       |      |         |           |             |         |     |            |     |     | RT1(Ω) |       | dG0(Ω) |     |
| Lemma |      | 5.3.1   | (Discrete | integration |         |     | by parts). | Let | u   |        | and p |        |     |
|       |      |         |           |             |         |     |            |     | ∈   |        | ∈     |        |     |

CHAPTER 5. DEVELOPMENT OF A CLASS OF LONG TIME CONSISTENT
STAGGERED SCHEMES ON THE WAVE SYSTEM 116
i)
σ ε (σ)u p + σ [[p]] u = σ p u ,
| |
K σ K
| |
σ σ
| |
Kσ σ
K σ ∂K σ int σ b
(cid:88)∈C (cid:88)∈ ∈(cid:88)F (cid:88)∈F
which is the discrete equivalent of
div(u)pdx+ u ∇pdx = u npdΓ.
· ·
(cid:90)Ω (cid:90)Ω (cid:90) ∂Ω
ii)
1
σ ε (σ)[[p]] p = σ [[p]]2 + σ [[p2]] ,
K σ K σ σ
| | − | | 2 | |
K σ ∂K σ σ b
(cid:88)∈C (cid:88)∈ (cid:88)∈F (cid:88)∈F
which is, up to a metric, the discrete equivalent of
∂p
∆ppdx = ∇p 2dx+ pdΓ.
− | | ∂n
(cid:90)Ω (cid:90)Ω (cid:90) ∂Ω
iii)
σ [[d (cid:93) ivu]] u = ∂K (d (cid:93) ivu)2 + σ (d (cid:93) ivu) u ,
| |
σ σ
− | |
K
| |
Kσ σ
σ int K σ b
∈(cid:88)F (cid:88)∈C (cid:88)∈F
which is the discrete equivalent, again up to a metric, of
∇div(u) udx = (div(u))2dx+ u n div(u)dΓ.
· − ·
(cid:90)Ω (cid:90)Ω (cid:90) ∂Ω
Proof. point i) is obtained by switching the sum on the cells to a sum on the faces:
σ ε (σ)u p = σ ε (σ)u p
K σ K K σ K
| | | |
K σ ∂K σ K (σ)
(cid:88)∈C (cid:88)∈ (cid:88)∈F (cid:88)∈C
= σ u ε (σ)p
| | σ K K (5.32)
σ K (σ)
(cid:88)∈F (cid:88)∈C
= σ u ε (σ)p + σ u ε (σ)p ,
σ K K σ K K
| | | |
σ int K (σ) σ b K (σ)
∈(cid:88)F (cid:88)∈C (cid:88)∈F (cid:88)∈C
(5.32) yields i) by definition of the discrete gradient Definition 5.2.3 and of (σ) for a boundary
C
face σ
σ ε (σ)u p = σ [[p]] u + σ p u , (5.33)
| |
K σ K
− | |
σ σ
| |
Kσ σ
K σ ∂K σ int σ b
(cid:88)∈C (cid:88)∈ ∈(cid:88)F (cid:88)∈F
which gives i).
As for point ii),we first recall that for a R, (b ) Rm
i 1 i m
∈ ≤≤ ∈
m m
2 (b a)a = (b2 a2 (b a)2). (5.34)
i − i − − i −
i=1 i=1
(cid:80) (cid:80)

CHAPTER 5. DEVELOPMENT OF A CLASS OF LONG TIME CONSISTENT
STAGGERED SCHEMES ON THE WAVE SYSTEM 117
In parallel, for a given cell K, by definition of the pressure gradient Definition 5.2.3 we have
1
σ ε (σ)[[p]] p = σ (p p )p + σ ((p ) p )p .
| |
K σ K
| |
Lσ
−
K K
2 | |
b σ
−
K K (5.35)
σ ∂K σ ∂K int σ ∂K b
(cid:88)∈ ∈ (cid:88)∩F ∈(cid:88)∩F
Using (5.34) in (5.35) we obtain
1
σ ε (σ)[[p]] p = σ p2 p2 (p p )2
| |
K σ K
2 | |
Lσ
−
K
−
Lσ
−
K
σ ∂K σ ∂K int (cid:18) (cid:19)
(cid:88)∈ ∈ (cid:88)∩F
(5.36)
1
+ σ (p )2 p2 ((p ) p )2 ,
b σ K b σ K
4 | | − − −
σ ∂K b (cid:18) (cid:19)
∈(cid:88)∩F
using again the definition of the pressure gradient (Definition 5.2.3) (5.36) can be written
1
σ ε K (σ)[[p]] σ p K = σ ε K (σ)[[p2]] σ [[p]]2 σ . (5.37)
| | 2 | | −
(cid:32) (cid:33)
σ ∂K σ ∂K
(cid:88)∈ (cid:88)∈
Summing on the cells K we have
1
σ ε K (σ)[[p]] σ p K = σ ε K (σ)[[p2]] σ [[p]]2 σ . (5.38)
| | 2 | | −
(cid:32) (cid:33)
K σ ∂K K σ ∂K
(cid:88)∈C (cid:88)∈ (cid:88)∈C (cid:88)∈
But
1 1 1
σ ε (σ)[[p2]] = σ [[p2]] ε (σ)+ σ ε (σ)[[p2]] . (5.39)
K σ σ K K σ
2 | | 2 | | 2 | |
K σ ∂K σ int K (σ) σ b
(cid:88)∈C (cid:88)∈ ∈(cid:88)F (cid:88)∈C (cid:88)∈F
=0
Using (5.39) in (5.38) yields ii). (cid:124) (cid:123)(cid:122) (cid:125)
Point iii) is similar to point i), we detail how the metric factor emerges. By Definition 4.4.1
of the grad-div diffusion we have
(cid:93) (cid:93)
σ [[divu]] u = σ ε (σ)(divu) u
σ σ K K σ
| | − | |
σ int σ intK (σ)
∈(cid:88)F ∈(cid:88)F (cid:88)∈C
(cid:93) (cid:93)
= σ ε (σ)(divu) u + σ (divu) u
− | |
K K σ
| |
Kσ σ
K σ ∂K σ b
(cid:88)∈C (cid:88)∈ (cid:88)∈F
(cid:93) (cid:93)
= (divu) σ ε (σ)u + σ (divu) u .
−
K
| |
K σ
| |
Kσ σ
K σ ∂K σ b
(cid:88)∈C (cid:88)∈ (cid:88)∈F
(cid:93)
iii) follows by using the definition of (divu) .
K

| CHAPTER   |     | 5.      | DEVELOPMENT |           |        | OF             | A CLASS | OF     | LONG   | TIME     | CONSISTENT |     |     |
| --------- | --- | ------- | ----------- | --------- | ------ | -------------- | ------- | ------ | ------ | -------- | ---------- | --- | --- |
| STAGGERED |     | SCHEMES |             |           | ON THE | WAVE           |         | SYSTEM |        |          |            |     | 118 |
| Lemma     |     | 5.3.2   | (Inverse    | Poincar´e |        | inequalities). |         |        | i) for | p dG0(Ω) |            |     |     |
∈
|     |       |     |        |      |              |                | σ [[p]] | )2  | σ     |              | σ [[p]]2, |     |     |
| --- | ----- | --- | ------ | ---- | ------------ | -------------- | ------- | --- | ----- | ------------ | --------- | --- | --- |
|     |       |     |        |      |              | (              |         | σ   | 2 max |              |           | σ   |     |
|     |       |     |        |      |              |                | | |     | ≤   | | |   |              | | |       |     |     |
|     |       |     |        |      | K (cid:88)∈C | σ (cid:88)∈ ∂K |         |     |       | σ (cid:88)∈F |           |     |     |
|     | where | we  | recall | that | σ max        | := max[        | σ       | ]   |       |              |           |     |     |
|     |       |     |        |      | | |          | σ              | | |     |     |       |              |           |     |     |
∈F
|     |         | RT1(Ω), |     |     | (cid:93) |         |               |     |       |      |     |     |     |
| --- | ------- | ------- | --- | --- | -------- | ------- | ------------- | --- | ----- | ---- | --- | --- | --- |
|     | ii) For | u       |     | let | d ivu    | defined | as Definition |     | 5.2.5 | then |     |     |     |
∈
|     |       |     |            |     |            | (cid:93) |        |     |            |     | (cid:93)  |     |     |
| --- | ----- | --- | ---------- | --- | ---------- | -------- | ------ | --- | ---------- | --- | --------- | --- | --- |
|     |       |     |            |     |            | σ 2[[d   | ivu]]2 | 2ν  |            | ∂K  | 2(d ivu)2 | ,   |     |
|     |       |     |            |     |            |          |        | σ   | max        |     |           | K   |     |
|     |       |     |            |     |            | | |      |        | ≤   |            | |   | |         |     |     |
|     |       |     |            |     | σ int      |          |        |     | K          |     |           |     |     |
|     |       |     |            |     | ∈(cid:88)F |          |        |     | (cid:88)∈C |     |           |     |     |
|     | where | ν   | := maxν(K) |     |            |          |        |     |            |     |           |     |     |
max
K
∈C
| Proof. | For | i) we | have |              |                |           |     |       |              |                |         |     |        |
| ------ | --- | ----- | ---- | ------------ | -------------- | --------- | --- | ----- | ------------ | -------------- | ------- | --- | ------ |
|        |     |       |      |              |                |           | )2  |       |              |                | [[p]]2, |     |        |
|        |     |       |      |              |                | ( σ [[p]] | σ   | σ max |              |                | σ       |     | (5.40) |
|        |     |       |      |              |                | | |       |     | ≤ | | |              |                | | |     | σ   |        |
|        |     |       |      | K (cid:88)∈C | σ (cid:88)∈ ∂K |           |     |       | K (cid:88)∈C | σ (cid:88)∈ ∂K |         |     |        |
but
|     |     |     |       |                      |     | [[p]]2 |         |          |            | [[p]]2 |            | [[p]]2   |        |
| --- | --- | --- | ----- | -------------------- | --- | ------ | ------- | -------- | ---------- | ------ | ---------- | -------- | ------ |
|     |     | σ   |       |                      | σ   |        | = σ     | 2        |            | σ      | +          | σ        |        |
|     |     | |   | | max |                      | |   | | σ    | | | max |          |            | | |    | σ          | | | σ    |        |
|     |     |     |       |                      |     |        |         | (cid:34) |            |        |            | (cid:35) |        |
|     |     |     | K     | σ                    | ∂K  |        |         | σ        | int        |        | σ          | b        |        |
|     |     |     |       | (cid:88)∈C (cid:88)∈ |     |        |         |          | ∈(cid:88)F |        | (cid:88)∈F |          | (5.41) |
|     |     |     |       |                      |     |        | 2 σ     |          | σ [[p]]2.  |        |            |          |        |
|     |     |     |       |                      |     |        |         | max      |            | σ      |            |          |        |
|     |     |     |       |                      |     |        | ≤ | |   |          | | |        |        |            |          |        |
σ
(cid:88)∈F
| So  | using | (5.41)  | in (5.40)      | we  | obtain | i). |     |     |     |     |     |     |     |
| --- | ----- | ------- | -------------- | --- | ------ | --- | --- | --- | --- | --- | --- | --- | --- |
|     | Then  | for ii) | by definition, |     |        |     |     |     |     |     |     |     |     |
2
|     |      |         |                  | 2[[d       | (cid:93) ivu]]2 |                      |                  | 2          |                  |            | (cid:93)       |          |        |
| --- | ---- | ------- | ---------------- | ---------- | --------------- | -------------------- | ---------------- | ---------- | ---------------- | ---------- | -------------- | -------- | ------ |
|     |      |         |                  | σ          |                 | :=                   |                  | σ          |                  | ε K        | (σ)(d ivu)     | K ,      | (5.42) |
|     |      |         |                  | | |        |                 | σ                    |                  | | |        |                  | −          |                |          |        |
|     |      |         |                  |            |                 |                      |                  | (cid:32)   |                  |            |                | (cid:33) |        |
|     |      |         | σ ∈(cid:88)F int |            |                 |                      | σ ∈(cid:88)F int |            | K (cid:88)∈C (σ) |            |                |          |        |
| we  | know | that (a | b)2              | 2(a2+b2)   |                 | so                   | (5.42)           | yields     |                  |            |                |          |        |
|     |      |         | −                | ≤          |                 |                      |                  |            |                  |            |                |          |        |
|     |      |         |                  |            |                 | 2[[d (cid:93) ivu]]2 |                  |            | 2                |            | (cid:93) ivu)2 |          |        |
|     |      |         |                  |            | σ               |                      |                  | 2          | σ                |            | (d             | .        |        |
|     |      |         |                  |            | | |             |                      | σ ≤              |            | | |              |            |                | K        | (5.43) |
|     |      |         |                  | σ          | int             |                      |                  | σ int      |                  | K (σ)      |                |          |        |
|     |      |         |                  | ∈(cid:88)F |                 |                      |                  | ∈(cid:88)F |                  | (cid:88)∈C |                |          |        |
Now,
|     |       |     |            |     | 2          | (cid:93) ivu)2 |     |            | (cid:93) ivu)2 |     |            | 2,  |        |
| --- | ----- | --- | ---------- | --- | ---------- | -------------- | --- | ---------- | -------------- | --- | ---------- | --- | ------ |
|     |       |     |            |     | σ          | (d             |     | =          | (d             |     |            | σ   | (5.44) |
|     |       |     |            | |   | |          |                | K   |            |                | K   |            | | | |        |
|     |       |     | σ          | int | K          | (σ)            |     | K          |                | σ   | ∂K int     |     |        |
|     |       |     | ∈(cid:88)F |     | (cid:88)∈C |                |     | (cid:88)∈C |                | ∈   | (cid:88)∩F |     |        |
| but | since | K   |            |     |            |                |     |            |                |     |            |     |        |
|     |       | ∀ ∈ | C          |     |            |                |     |            |                |     |            |     |        |
|     |       |     |            |     |            | σ 2            |     | ∂K         | 2 ν            | ∂K  | 2,         |     |        |
max
|     |     |     |     |     |           | | | ≤ |           | |   | | ≤ | |   | |   |     |     |
| --- | --- | --- | --- | --- | --------- | ----- | --------- | --- | --- | --- | --- | --- | --- |
|     |     |     |     |     | σ ∂K      |       | σ ∂K      |     |     |     |     |     |     |
|     |     |     |     |     | (cid:88)∈ |       | (cid:88)∈ |     |     |     |     |     |     |

| CHAPTER   |     | 5.      | DEVELOPMENT |     | OF  | A    | CLASS  | OF LONG | TIME | CONSISTENT |     |     |
| --------- | --- | ------- | ----------- | --- | --- | ---- | ------ | ------- | ---- | ---------- | --- | --- |
| STAGGERED |     | SCHEMES |             | ON  | THE | WAVE | SYSTEM |         |      |            |     | 119 |
(5.44) yields
|     |     |     |     |            |            | (cid:93) |     |     |            | (cid:93)  |     |        |
| --- | --- | --- | --- | ---------- | ---------- | -------- | --- | --- | ---------- | --------- | --- | ------ |
|     |     |     |     |            | σ 2        | (d ivu)2 |     | ν   | ∂K         | 2(d ivu)2 | .   | (5.45) |
|     |     |     |     |            |            |          | K   | max |            |           | K   |        |
|     |     |     |     |            | | |        |          | ≤   |     | | |        |           |     |        |
|     |     |     |     | σ int      | K (σ)      |          |     | K   |            |           |     |        |
|     |     |     |     | ∈(cid:88)F | (cid:88)∈C |          |     |     | (cid:88)∈C |           |     |        |
ii).
| Plugging   |     | (5.45)        | in (5.43) | gives      |               |            |        |     |            |     |     |        |
| ---------- | --- | ------------- | --------- | ---------- | ------------- | ---------- | ------ | --- | ---------- | --- | --- | ------ |
| 5.3.2      |     | Stability     | results   |            |               |            |        |     |            |     |     |        |
| Definition |     | 5.3.1         | (Energy   |            | dissipation). | We         | define | the | energy     |     |     |        |
|            |     |               |           |            |               |            |        | p2  |            | u2  |     |        |
|            |     |               |           |            |               |            |        | K   |            | σ.  |     |        |
|            |     |               |           | E(U        | ) :=          |            | K ρ κ  | +   | D σ        |     |     | (5.46) |
|            |     |               |           |            | h             |            | | | 0  | 0 2 | |          | | 2 |     |        |
|            |     |               |           |            |               | K          |        |     | σ          |     |     |        |
|            |     |               |           |            |               | (cid:88)∈C |        |     | (cid:88)∈F |     |     |        |
| We         | say | that a scheme |           | dissipates | the           | energy     | if     |     |            |     |     |        |
|            |     |               |           |            | E(Un+1)       | E(Un)      |        |     |            |     |     |        |
|            |     |               |           |            | h             | −          | h      | +   | σ Φ        | 0,  |     |        |
σ
|     |     |     |     |     |     | δτ  |     |     | | | ≤ |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | ----- | --- | --- | --- |
|     |     |     |     |     |     |     |     | σ   | b     |     |     |     |
(cid:88)∈F
| with | Φ   | a flux that | depends |     | on the solution |     | at time | step | n or n+1. |     |     |     |
| ---- | --- | ----------- | ------- | --- | --------------- | --- | ------- | ---- | --------- | --- | --- | --- |
σ
In order to study the energy dissipation of the numerical schemes, we make two simplifica-
tionsonthetreatmentoftheboundaryterms: firstly, thevelocitydegreesoffreedomlocatedon
boundary faces will now be strongly imposed. Then for all time step if σ is a boundary face u σ
|     |     |     |     |     |     |     |     |     |     | u   | n   |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
the normal component of the velocity at this face will be equal to b σ , the boundary velocity
·
we impose. Secondly, we will consider that only inlet/outlet boundary conditions are applied
on the whole boundary ∂Ω. The results in the general case where a wall flux is imposed on a
part of the boundary is obtained by simply changing in the proof (u ) = u and (p ) = p
|     |     |     |     |     |     |     |     |     |     | b   | σ σ | b σ Kσ |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | ------ |
−
| for | σ   | b . Moreover, |     | we  | denote |     |     |     |     |     |     |     |
| --- | --- | ------------- | --- | --- | ------ | --- | --- | --- | --- | --- | --- | --- |
Fwall
∈
|     |     |     |     |     | h   | = min(min |     | D ,min | K ) |     |     | (5.47) |
| --- | --- | --- | --- | --- | --- | --------- | --- | ------ | --- | --- | --- | ------ |
σ
|     |     |     |     |     |     |     | σ | | | K | | | |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
|     |     |     |     |     |     |     | ∈F  | ∈C  |     |     |     |     |
Proposition 5.3.1 (Energy dissipation of the Euler Explicit scheme). We consider the follow-
ing scheme:
|     |     | pn+1 | pn  |     |     |     |     |     |     |     |     |     |
| --- | --- | ---- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
1
|     |     | K K    | − K | +   | σ ε       | (σ)un | = dc |           | σ ε (σ)[[pn]] |     | K ,   | (5.48a) |
| --- | --- | ------ | --- | --- | --------- | ----- | ---- | --------- | ------------- | --- | ----- | ------- |
|     |     |        |     |     |           | K     | σ    | 0         | K             |     | σ     |         |
|     |     | | | δτ |     | ρ   | | |       |       |      |           | | |           |     | ∀ ∈ C |         |
|     |    |        |     | 0   | σ ∂K      |       |      | σ ∂K      |               |     |       |         |
|     |     |        |     |     | (cid:88)∈ |       |      | (cid:88)∈ |               |     |       |         |
 
|     |    |     |     | un  | +1 un |     |          |      | (cid:94)     |     |        |         |
| --- | --- | --- | --- | --- | ----- | --- | -------- | ---- | ------------ | --- | ------ | ------- |
|     |   |     |     | D σ | σ     | +κ  | σ [[pn]] | = dc | σ [[d ivun]] |     | σ int, | (5.48b) |
|     |    |     |     | σ   | −     | 0   |          | σ    | 0            | σ   |        |         |
|     |    |     | |   | |   | δτ    | |   | |        |      | | |          |     | ∀ ∈ F  |         |
|     |     |     |     |     |       |     |          |      | un           |     | b.     |         |
|     |    |     |     |     |       |     |          |      | =            | u b | n σ σ  | (5.48c) |
|     |    |     |     |     |       |     |          |      | σ            | ·   | ∀ ∈ F  |         |
 

 

CHAPTER 5. DEVELOPMENT OF A CLASS OF LONG TIME CONSISTENT
STAGGERED SCHEMES ON THE WAVE SYSTEM 120
The scheme (5.48) dissipates the energy in the sense of Definition 5.3.1 with boundary flux:
σ b Φn := κ pn un dc (d (cid:94) ivun) un dc ρ κ [[(pn)2]] ,
∀ ∈ F
σ 0 Kσ σ
−
0 Kσ σ
−
0 0 0 σ
under the CFL condition:
c δτmax( ∂K ) 1
0 | | ,
h ≤ 1 (5.49)
2dν +
max
d
where h is defined by (5.47).
Proposition 5.3.2 (Energy dissipation for the Pressure centred ImEx scheme). We consider
the following scheme
pn+1 pn 1
K K − K + σ ε (σ)un = dc σ ε (σ)[[pn]] K , (5.50a)
| | δτ ρ | | K σ 0 | | K σ ∀ ∈ C
 0
σ ∂K σ ∂K
 (cid:88)∈ (cid:88)∈
  un+1 un
    | D σ | σ δτ − σ +κ 0 | σ | [[pn+1]] σ = 0 ∀ σ ∈ F int, (5.50b)
un = u n σ b. (5.50c)
   σ b · σ ∀ ∈ F



the scheme dissipates the energy in the sense of Definition 5.3.1 with boundary flux
dc ρ κ
σ b Φn := κ pn un δτc2(divun) un 0 0 0 [[(pn)2]] +δτdc κ (∆pn) un,
∀ ∈ F σ 0 Kσ σ − 0 Kσ σ − 2 σ 0 0 Kσ σ
under the CFL condition: (cid:103)
δτc σ 1
0 | | max ,
h ≤ 3 (5.51)
ν (1+d)+
max
2d
where h is defined by (5.47).
Proposition 5.3.3 (Energy dissipation for Euler implicit time integration). We consider now
the following scheme
pn+1 pn 1
K K − K + σ ε (σ)un+1 = 0 K , (5.52a)
| | δτ ρ | | K σ ∀ ∈ C
 0
σ ∂K
 (cid:88)∈
  un+1 un
    | D σ | σ δτ − σ +κ 0 | σ | [[pn+1]] σ = 0 ∀ σ ∈ F int, (5.52b)
un = u n σ b. (5.52c)
   σ b · σ ∀ ∈ F



i) For general inlet/outlet conditions, the scheme (5.52) dissipates the energy uncondition-
ally with the boundary flux:
σ b Φn+1 := κ pn+1un+1.
∀ ∈ F σ 0 Kσ σ

CHAPTER 5. DEVELOPMENT OF A CLASS OF LONG TIME CONSISTENT
STAGGERED SCHEMES ON THE WAVE SYSTEM 121
ii) With periodic or Neumann boundary conditions, given Un, (5.52) ensures a unique solu-
h
tion Un+1.
h
Remark 5.3.1. Notice that the semi-implicit scheme does not require the inversion of a matrix
and can be rewritten in a fully explicit form. Indeed using
δτ 1
pn+1 = pn + K div(un) +dc σ ε (σ)[[pn]] ,
K K K − ρ | | K 0 | | K σ
(cid:32) 0 (cid:33)
| | σ ∂K
(cid:88)∈
1
we have seen that the system can be rewritten, by denoting ∆pn := σ ε (σ)[[pn]] ,
K K | | K σ
| | σ ∂K
in the following form:
(cid:88)∈
(cid:103)
pn+1 pn 1
K K − K + σ ε (σ)un = dc σ ε (σ)[[pn]] K ,
| | δτ ρ | | K σ 0 | | K σ ∀ ∈ M
 0
σ ∂K σ ∂K
    D un σ +1 − un σ +κ σ [[p (cid:88)∈ n]] = δτc2 σ [[div(un) (cid:88) ] ∈ ] δτdc κ σ [[∆pn]] σ int, (5.53)
   | σ | δτ 0 | | σ 0| | σ − 0 0 | | σ ∀ ∈ F

  un σ = u b n σ σ b. (cid:103)
 · ∀ ∈ F





5.3.3 Proof of Proposition 5.3.1
Step 1: we make appear an equation on the global energy
We multiply (5.48a) by pn ρ κ
K 0 0
pn+1 pn 1
K K − Kρ κ pn + σ ε (σ)unρ κ pn = dc σ ε (σ)[[pn]] ρ κ pn , (5.54)
| | δτ 0 0 K ρ | | K σ 0 0 K 0 | | K σ 0 0 K
0
σ ∂K σ ∂K
(cid:88)∈ (cid:88)∈
and (5.48b) by un
σ
D
un
σ
+1
−
un
σun+κ σ [[pn]] un = dc σ [[d (cid:94) ivun]] un. (5.55)
| σ | δτ σ 0 | | σ σ 0 | | σ σ
Using 2(a b)b = a2 b2 (a b)2 we have
− − − −
pn+1 pn K (pn+1)2 (pn )2 K
K K − Kρ κ pn = | |ρ κ K K | |ρ κ (pn+1 pn )2,
| | δτ 0 0 K δτ 0 0 2 − 2 − 2δτ 0 0 K − K
 (cid:18) (cid:19) (5.56)
un+1 un D (un+1)2 (un)2 D
   | D σ | σ δτ − σun σ = | δτ σ | σ 2 − 2 σ − | 2δ σ τ |(un σ +1 − un σ )2.
(cid:18) (cid:19)




CHAPTER 5. DEVELOPMENT OF A CLASS OF LONG TIME CONSISTENT
STAGGERED SCHEMES ON THE WAVE SYSTEM 122
Using (5.56) in (5.54) yields:
K (pn+1)2 (pn )2 K 1
| |ρ κ K K | |ρ κ (pn+1 pn )2+ σ ε (σ)unρ κ pn
δτ 0 0 2 − 2 − 2δτ 0 0 K − K ρ | | K σ 0 0 K
(cid:18) (cid:19) 0 σ ∂K (5.57)
= dc σ ε (σ)[[pn]] ρ κ pn ,
(cid:88)∈
0 K σ 0 0 K
| |
σ ∂K
(cid:88)∈
and (5.56) in (5.55) yields:
| D σ | (un σ +1)2 (un σ )2 | D σ |(un+1 un)2+κ σ [[pn]] un = dc σ [[d (cid:94) ivun]] un . (5.58)
δτ 2 − 2 − 2δτ σ − σ 0 | | σ σ 0 | | σ σ
(cid:18) (cid:19)
Summing (5.57) on the cells K gives
∈ C
K (pn+1)2 (pn )2
| |ρ κ K K +κ σ ε (σ)unpn
δτ 0 0 2 − 2 0 | | K σ K
K (cid:18) (cid:19) K σ ∂K
(cid:88)∈C (cid:88)∈C (cid:88)∈
(5.59)
K
= dc σ ε (σ)[[pn]] pn ρ κ + | |ρ κ (pn+1 pn )2 ,
0 | | K σ K 0 0 2δτ 0 0 K − K
(cid:32) (cid:33)
K σ ∂K
(cid:88)∈C (cid:88)∈
and similarly, summing (5.58) on the face σ int
∈ F
D (un+1)2 (un)2
| σ | σ σ +κ σ [[pn]] un
δτ 2 − 2 0 | | σ σ
σ int (cid:18) (cid:19) σ int
∈(cid:88)F ∈(cid:88)F
(5.60)
= dc σ [[d (cid:94) ivun]] un+ | D σ |(un+1 un)2 .
0 | | σ σ 2δτ σ − σ
(cid:32) (cid:33)
σ int
∈(cid:88)F
By summing (5.59) and (5.60) and using point i) of Lemma 5.3.1:
E(Un+1) E(Un)
h δτ − h +κ 0 | σ | pn Kσ un σ = L, (5.61)
σ b
(cid:88)∈F
where
K
L = dc σ ε (σ)[[pn]] pn ρ κ + | |ρ κ (pn+1 pn )2
0 | | K σ K 0 0 2δτ 0 0 K − K
(cid:32) (cid:33)
K σ ∂K
(cid:88)∈C (cid:88)∈
(5.62)
+ dc [[d (cid:94) ivun]] un+ | D σ |(un+1 un)2 .
0 σ σ 2δτ σ − σ
(cid:32) (cid:33)
σ int
∈(cid:88)F
Step 2: we make appear the diffusion thanks to discrete integrations by parts

CHAPTER 5. DEVELOPMENT OF A CLASS OF LONG TIME CONSISTENT
STAGGERED SCHEMES ON THE WAVE SYSTEM 123
Now we use ii) of Lemma 5.3.1 so that
dc ρ κ
dc σ ε (σ)[[pn]] pn ρ κ = dc ρ κ σ [[pn]]2 + 0 0 0 σ [[(pn)2]] ,
0 K σ K 0 0 0 0 0 σ σ
| | − | | 2 | |
K σ ∂K σ σ b
(cid:88)∈C (cid:88)∈ (cid:88)∈F (cid:88)∈F
(5.63)
point iii) of Lemma 5.3.1 yields
dc [[d (cid:94) ivun]] un = dc ∂K (d (cid:94) ivun)2 +dc σ (d (cid:94) ivun) un. (5.64)
0 σ σ
−
0
| |
K 0
| |
Kσ σ
σ int K σ b
∈(cid:88)F (cid:88)∈C (cid:88)∈F
So that plugging (5.63) and (5.64) in (5.62) yields
dc ρ κ K
L = dc ρ κ σ [[pn]]2 + 0 0 0 σ [[(pn)2]] + | |ρ κ (pn+1 pn )2
− 0 0 0 | | σ 2 | | σ 2δτ 0 0 K − K
σ σ b K
(cid:88)∈F (cid:88)∈F (cid:88)∈C
(5.65)
dc ∂K (d (cid:94) ivun)2 +dc ∂K (d (cid:94) ivun) u + | D σ |(un+1 un)2.
− 0 | | K 0 | | Kσ σ 2δτ σ − σ
K σ b σ int
(cid:88)∈C (cid:88)∈F ∈(cid:88)F
Using (5.65) in (5.61) we obtain then
E(Un+1) E(Un)
h δτ − h + | σ | Φn σ = R, (5.66)
σ b
(cid:88)∈F
with the boundary flux:
Φn := κ pn un dc (d (cid:94) ivun) un dc 0 ρ 0 κ 0 [[(pn)2]] ,
σ 0 Kσ σ
−
0 Kσ σ
− 2
σ
with R equal to
R = dc ∂K (d (cid:94) ivun)2 dc ρ κ σ [[pn]]2
0 K 0 0 0 σ
− | | − | |
K σ
(cid:88)∈C (cid:88)∈F
diffusiveterms
(5.67)
(cid:124) (cid:123)(cid:122) (cid:125)
K D
+ | |ρ κ (pn+1 pn )2+ | σ |(un+1 un)2.
2δτ 0 0 K − K 2δτ σ − σ
K σ int
(cid:88)∈C ∈(cid:88)F
non-negativeterms
Step 3: we bound (cid:124)the non-negative term(cid:123)s(cid:122) in the rest R as facto(cid:125)rs of the diffusive
terms

CHAPTER 5. DEVELOPMENT OF A CLASS OF LONG TIME CONSISTENT
STAGGERED SCHEMES ON THE WAVE SYSTEM 124
Obviously, by squaring the equations (5.48), we have
2
| K |ρ κ (pn+1 pn )2 = δτ ρ κ dc σ ε (σ)[[pn]] | ∂K |(d (cid:94) ivun)
 2δτ 0 0 K − K 2 K 0 0 (cid:32) 0 | | K σ − ρ 0 K (cid:33)
| | σ ∂K
  (cid:88)∈



  2
 | D σ |(un+1 un)2 = δτ dc σ [[d (cid:94) ivun]] κ σ [[pn]] ,
   2δτ σ − σ 2 | D σ | (cid:32) 0 | | σ − 0 | | σ (cid:33)



so us  ing (a b)2 2(a2+b2) for a,b R
− ≤ ∈
| K |ρ κ (pn+1 pn )2 δτ ρ κ dc σ ε (σ)[[pn]] 2 + | ∂K |(d (cid:94) ivun) 2 , (5.68)
2δτ 0 0 K − K ≤ K 0 0 0 | | K σ ρ K
| | (cid:34) (cid:18) σ ∂K (cid:19) (cid:18) 0 (cid:19) (cid:35)
(cid:88)∈
and
| D σ |(un+1 un)2 δτ dc σ [[d (cid:94) ivun]] 2 + κ σ [[pn]] ) 2 . (5.69)
2δτ σ − σ ≤ D 0 | | σ 0 | | σ
σ (cid:34) (cid:35)
| | (cid:18) (cid:19) (cid:18) (cid:19)
But using Jensen’s inequality on the square function we have
2
dc σ ε (σ)[[pn]] d2c2ν(K) σ 2[[pn]]2 d2c2ν σ 2[[pn]]2. (5.70)
0 | | K σ ≤ 0 | | σ ≤ 0 max | | σ
(cid:18) σ ∂K (cid:19) σ ∂K σ ∂K
(cid:88)∈ (cid:88)∈ (cid:88)∈
Plugging (5.70) in (5.68) we have
| K |ρ κ (pn+1 pn )2 δτ ρ κ d2c2ν σ 2[[pn]]2 + | ∂K |(d (cid:94) ivun) 2 . (5.71)
2δτ 0 0 K − K ≤ K 0 0 0 max | | σ ρ K
| | (cid:34) σ ∂K (cid:18) 0 (cid:19) (cid:35)
(cid:88)∈
Summing (5.71) on the cells K and (5.69) on the faces σ int and summing the result we
∈ F
obtain the following bound:
R δτ ρ κ d2c2ν ( σ [[pn]] )2+ | ∂K | 2 d (cid:94) ivun 2
≤ K 0 0 0 max | | σ ρ2 K
K (cid:88)∈C | | (cid:32) σ (cid:88)∈ ∂K 0 (cid:16) (cid:17) (cid:33)
+ δτ d2c2 σ 2[[d (cid:94) ivun]]2 +κ2 σ 2[[pn]]2 (5.72)
D 0| | σ 0| | σ
σ (cid:32) (cid:33)
σ int | |
∈(cid:88)F
dc ρ κ σ [[pn]]2 dc ∂K d (cid:94) ivun 2 .
0 0 0 σ 0 K
− | | − | |
σ (cid:88)∈F K (cid:88)∈C (cid:16) (cid:17)
Now using the discrete inverse Poincar´e inequalities from Lemma 5.3.2, yields ,with h =

CHAPTER 5. DEVELOPMENT OF A CLASS OF LONG TIME CONSISTENT
STAGGERED SCHEMES ON THE WAVE SYSTEM 125
min(min D ,min K ), σ := max σ
σ σ K max
| | | | | | σ | |
∈F
δτ δτ σ c 2dν
ρ κ d2c2ν ( σ [[pn]] )2 | | max 0 maxρ κ dc σ [[pn]]2, (5.73)
K 0 0 0 max | | σ ≤ h 0 0 0 | | σ
K | | σ ∂K σ
(cid:88)∈C (cid:88)∈ (cid:88)∈F
and
δτ d2c2 σ 2[[d (cid:94) ivun]]2 dc ∂K 2(d (cid:94) ivun)2 dc δτ 2ν . (5.74)
D 0| | σ ≤ 0 | | K 0 h max
σ
σ int | | K
∈(cid:88)F (cid:88)∈C
Gathering (5.73) and (5.74) in (5.72) we obtain
c δτ σ 1 c δτ σ
R dc ρ κ σ [[pn]]2 0 | | max 1 +dc ρ κ σ [[pn]]2 0 | | max 2dν
≤ 0 0 0 | | σ h d − 0 0 0 | | σ h max
(cid:32) (cid:33)
σ σ
(cid:88)∈F (cid:88)∈F
(5.75)
+dc d (cid:94) ivun 2 ∂K 2 c 0 δτ 1 ∂K +dc ∂K 2(d (cid:94) ivun)2 dc δτ 2ν .
0 K(cid:32) | | h d −| | (cid:33) 0 | | K 0 h max
K (cid:88)∈C (cid:16) (cid:17) K (cid:88)∈C
Reorganizing this gives
R ρ κ σ [[pn]]2 + ∂K (d (cid:94) ivun)2 c 0 δτmax( | ∂K | ) 2dν + 1 1 dc .
≤ 0 0 | | σ | | K h max d − 0
(cid:32) (cid:33)
(cid:20) σ K (cid:21) (cid:18) (cid:19)
(cid:88)∈F (cid:88)∈C
(5.76)
Using (5.76) in (5.66) yields
E(Un+1) E(Un)
h − h + σ Φn
δτ | | σ ≤
σ ∂Ω
(cid:88)∈
ρ κ σ [[pn]]2 + ∂K (d (cid:94) ivun)2 c 0 δτmax( | ∂K | ) 2dν + 1 1 dc .
0 0 | | σ | | K h max d − 0
(cid:32) (cid:33)
(cid:20) σ K (cid:21) (cid:18) (cid:19)
(cid:88)∈F (cid:88)∈C
So that energy dissipation occurs when
c δτmax( ∂K ) 1
0 | | 2dν + 1 0,
max
h d − ≤
(cid:18) (cid:19)
yielding the CFL (5.49).
Remark 5.3.2. We can easily optimize the best coefficient for d by searching the global min-
1 1 1
imum of f(x) = + 2xν , f (x) = + 2ν this means f (x) = 0 iff x = .
x max (cid:48) −x2 max (cid:48) √2ν
max
1 1
f (x) 0 when 0 < x and f (x) 0 when x so the maximum is reached
(cid:48) (cid:48)
≥ ≤ √2ν ≤ ≥ √2ν
max max

CHAPTER 5. DEVELOPMENT OF A CLASS OF LONG TIME CONSISTENT
STAGGERED SCHEMES ON THE WAVE SYSTEM 126
1 1
in d = with value .
√2ν 2√2ν
max max
The proof of Proposition 5.3.2 shares similar arguments with the previous section. Details
are given in Appendix A.
5.3.4 Proof of Proposition 5.3.3
Proof of point i): we multiply (5.52a) by pn+1ρ κ
K 0 0
pn+1 pn 1
K K − Kρ κ pn+1+ σ ε (σ)un+1ρ κ pn+1 = 0, (5.77)
| | δτ 0 0 K ρ | | K σ 0 0 K
0
σ ∂K
(cid:88)∈
and (5.52b) by un+1
σ
un+1 un
D σ − σun+1+κ σ [[pn+1]] un+1 = 0, (5.78)
| σ | δτ σ 0 | | σ σ
Using 2(a b)b = a2 b2+(a b)2 we have
− − −
pn+1 pn K (pn+1)2 (pn )2 K
K K − Kρ κ pn = | |ρ κ K K + | |ρ κ (pn+1 pn )2,
| | δτ 0 0 K δτ 0 0 2 − 2 2δτ 0 0 K − K
 (cid:18) (cid:19) (5.79)
un+1 un D (un+1)2 (un)2 D
   | D σ | σ δτ − σun σ = | δτ σ | σ 2 − 2 σ + | 2δ σ τ |(un σ +1 − un σ )2.
(cid:18) (cid:19)


Using (5.79) in (5.77) yields:
K (pn+1)2 (pn )2 K 1
| |ρ κ K K | |ρ κ (pn+1 pn )2+ σ ε (σ)un+1ρ κ pn+1 = 0,
δτ 0 0 2 − 2 − 2δτ 0 0 K − K ρ | | K σ 0 0 K
(cid:18) (cid:19) 0 σ ∂K
(cid:88)∈
(5.80)
and (5.79) in (5.78) yields:
D (un+1)2 (un)2 D
| σ | σ σ + | σ |(un+1 un)2+κ σ [[pn+1]] un+1 = 0. (5.81)
δτ 2 − 2 2δτ σ − σ 0 | | σ σ
(cid:18) (cid:19)
Summing (5.80) on the cells K gives
∈ C
K (pn+1)2 (pn )2
| |ρ κ K K +κ σ ε (σ)un+1pn+1
δτ 0 0 2 − 2 0 | | K σ K
K (cid:18) (cid:19) K σ ∂K
(cid:88)∈C (cid:88)∈C (cid:88)∈
(5.82)
K
= | |ρ κ (pn+1 pn )2,
− 2δτ 0 0 K − K
K
(cid:88)∈C

| CHAPTER        | 5.  | DEVELOPMENT |        |     | OF       | A CLASS | OF     | LONG | TIME | CONSISTENT |     |
| -------------- | --- | ----------- | ------ | --- | -------- | ------- | ------ | ---- | ---- | ---------- | --- |
| STAGGERED      |     | SCHEMES     |        | ON  | THE WAVE |         | SYSTEM |      |      |            | 127 |
| and similarly, |     | summing     | (5.81) | on  | the face | σ       | int    |      |      |            |     |
∈ F
|     |     |            |     | D        | (un+1)2 | (un)2 |          |            |          |      |     |
| --- | --- | ---------- | --- | -------- | ------- | ----- | -------- | ---------- | -------- | ---- | --- |
|     |     |            |     | σ        | σ       |       | σ        |            | [[pn+1]] | un+1 |     |
|     |     |            | |   | |        |         |       | +κ       |            | σ        |      |     |
|     |     |            |     | δτ       | 2       | −     | 2        | 0          | | |      | σ σ  |     |
|     |     | σ          | int | (cid:18) |         |       | (cid:19) | σ          | int      |      |     |
|     |     | ∈(cid:88)F |     |          |         |       |          | ∈(cid:88)F |          |      |     |
(5.83)
D
|     |     |     |     |     |     |     | σ |(un+1 | un)2. |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | -------- | ----- | --- | --- | --- |
|     |     |     |     |     | =   | |   |          |       |     |     |     |
|     |     |     |     |     |     |     | 2δτ σ    |       | σ   |     |     |
|     |     |     |     |     | −   |     |          | −     |     |     |     |
|     |     |     |     |     | σ   | int |          |       |     |     |     |
∈(cid:88)F
By summing (5.82) and (5.83) and using again point i) of Lemma 5.3.1:
|     |     |     | E(Un+1) |     | E(Un) |     |     |     |        |     |        |
| --- | --- | --- | ------- | --- | ----- | --- | --- | --- | ------ | --- | ------ |
|     |     |     |         |     |       | h   |     | pn  | 1un    |     |        |
|     |     |     |         | h   | −     | +κ  |     | σ   | + +1 = | R,  | (5.84) |
|     |     |     |         |     | δτ    |     | 0   | K   | σ      |     |        |
|     |     |     |         |     |       |     |     | | | | σ      |     |        |
|     |     |     |         |     |       |     | σ b |     |        |     |        |
(cid:88)∈F
| with the | boundary |     | flux: |     |     |     |     |     |     |     |     |
| -------- | -------- | --- | ----- | --- | --- | --- | --- | --- | --- | --- | --- |
pn+1un+1,
|     |     |     |     |     | Φn+1 | :=  | κ    |     |     |     |     |
| --- | --- | --- | --- | --- | ---- | --- | ---- | --- | --- | --- | --- |
|     |     |     |     |     | σ    |     | 0 Kσ | σ   |     |     |     |
and
|         |        |           |            | K     |      |     |       |            | D      |           |        |
| ------- | ------ | --------- | ---------- | ----- | ---- | --- | ----- | ---------- | ------ | --------- | ------ |
|         |        |           |            |       | (pn  | +1  | pn )2 |            | σ |(un | +1 un )2. |        |
|         |        | R         | =          | |     | |ρ κ |     |       |            | |      |           | (5.85) |
|         |        |           | −          | 2 δτ  | 0 0  | K − | K −   |            | 2δ τ σ | − σ       |        |
|         |        |           | K          |       |      |     |       | σ int      |        |           |        |
|         |        |           | (cid:88)∈C |       |      |     |       | ∈(cid:88)F |        |           |        |
| Since R | 0 this | concludes |            | point | i).  |     |       |            |        |           |        |
|         | ≤      |           |            |       |      |     |       |            |        |           | Φn+1   |
Now for point ii), in the periodic case, or with Neumann boundary conditions
σ
σ b
(cid:88)∈F
| vanishes | so that | (5.84) | becomes |     |         |     |       |      |     |     |        |
| -------- | ------- | ------ | ------- | --- | ------- | --- | ----- | ---- | --- | --- | ------ |
|          |         |        |         |     | E(Un+1) |     | E(Un) |      |     |     |        |
|          |         |        |         |     |         | h   | h     |      |     |     | (5.86) |
|          |         |        |         |     |         |     | −     | = R, |     |     |        |
δτ
Un+1,
in the implicit case we have to solve a linear system to obtain the solution we thus show
h
that given Un it exists a unique Un+1 associated to it through the scheme. This is simply done
h
h
by showing that the kernel of the matrix A associated to the scheme (5.52) by the relation
AUn+1 = Un is reduced to 0 : let Un = 0 then it implies, by (5.85) that R = E(Un+1)
|          |           |     |     |     |     | h   |     |     |     |     | h   |
| -------- | --------- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
|          |           |     |     | {   | }   |     |     |     |     | −   |     |
| yielding | in (5.86) |     |     |     |     |     |     |     |     |     |     |
2E(Un+1)
|     |     |     |     |     |     |     | =   | 0.  |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
h
| By definition |     | of E(.), | Un+1 | = 0. | This | ends | the proof | of existence. |     |     |     |
| ------------- | --- | -------- | ---- | ---- | ---- | ---- | --------- | ------------- | --- | --- | --- |
h
| 5.4 | Discrete |     | long | time | behaviour |     |     |     |     |     |     |
| --- | -------- | --- | ---- | ---- | --------- | --- | --- | --- | --- | --- | --- |
In this section, we show a discrete version of the relative energy dissipation Lemma 2.3.2 and
more importantly, we prove, that every solution built with our numerical scheme converges in
| long time | to  | this limit. | This | is  | detailed | in the | following | theorem: |     |     |     |
| --------- | --- | ----------- | ---- | --- | -------- | ------ | --------- | -------- | --- | --- | --- |

| CHAPTER   |     | 5.  | DEVELOPMENT |     |     | OF       | A CLASS | OF LONG |     | TIME | CONSISTENT |     |     |
| --------- | --- | --- | ----------- | --- | --- | -------- | ------- | ------- | --- | ---- | ---------- | --- | --- |
| STAGGERED |     |     | SCHEMES     |     | ON  | THE WAVE |         | SYSTEM  |     |      |            |     | 128 |
Theorem 5.4.1 (Semi-discrete long time limit). : Suppose that Ω R2 is an open connected
⊂
p
set. Let U p dG0(Ω), u RT1(Ω) the solution of the semi-discrete wave system
=
|     |     |     | u   | ∈   |     |     | ∈   |     |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
(cid:18) (cid:19)
(5.31) with:
|     | • initial | conditions |     | p(τ | = 0,x) | :=  | p (x) | and u(τ = | 0,x) | := u (x) |     |     |     |
| --- | --------- | ---------- | --- | --- | ------ | --- | ----- | --------- | ---- | -------- | --- | --- | --- |
|     |           |            |     |     |        |     | 0     |           |      | 0        |     |     |     |
•
|     | weakly | imposed |     | wall | conditions | on  | b   |     |     |     |     |     |     |
| --- | ------ | ------- | --- | ---- | ---------- | --- | --- | --- | --- | --- | --- | --- | --- |
Fwall
p
|     | •   |     |     | b   |     |     |     |     |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
a state U = weakly imposed on b which is such that, p ,u are time
|     |     |     | b   | u   |     |     |     | Finlet/outlet |     |     |     | b   | b   |
| --- | --- | --- | --- | --- | --- | --- | --- | ------------- | --- | --- | --- | --- | --- |
b
(cid:18) (cid:19)
|     | -independent |     | and |            |     |     |            |                   |     |           |     |       |     |
| --- | ------------ | --- | --- | ---------- | --- | --- | ---------- | ----------------- | --- | --------- | --- | ----- | --- |
|     |              |     | p   | is uniform |     | and | u verifies | the compatibility |     | condition |     |       |     |
|     |              |     |     | b          |     |     | b          |                   |     |           |     | (5.3) |     |
Suppose moreover that the semi-discrete scheme (5.31) induces energy dissipation.
|     | Then, | the solution |     | U(τ) | converges |     | in long | time and | its long | time | limit |     |     |
| --- | ----- | ------------ | --- | ---- | --------- | --- | ------- | -------- | -------- | ---- | ----- | --- | --- |
p
|     |     |     |     |     |     |     | U = | ∞   | ,   |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
∞
u
|     |     |     |     |     |     |     |     | (cid:18) ∞ (cid:19) |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | ------------------- | --- | --- | --- | --- | --- |
is such that
|     | • p | is uniform |     | equal | to p | in dG0(Ω) |     |     |     |     |     |     |     |
| --- | --- | ---------- | --- | ----- | ---- | --------- | --- | --- | --- | --- | --- | --- | --- |
|     |     | ∞          |     |       | b    |           |     |     |     |     |     |     |     |
|     | •   |            |     |       |      | RT1(Ω)    |     |     |     |     |     |     |     |
u is equal to (u ) , in the divergence-free part of the initial velocity u with
|     |          | ∞   |            | 0   | Ψ             |     |      |             |         |        |     |     | 0   |
| --- | -------- | --- | ---------- | --- | ------------- | --- | ---- | ----------- | ------- | ------ | --- | --- | --- |
|     | boundary |     | conditions |     | u , extracted |     | with | (4.36) with | Theorem | 4.4.4. |     |     |     |
b
In order to conclude in the semi-discrete setting, we will need the following intermediary
results
(p,u)t
Lemma 5.4.1 (Invariance of u ). Let solution of (5.31), then u from the decomposi-
|        |        |              |     |          |       | Ψ      |      |            |     |     |     | Ψ   |     |
| ------ | ------ | ------------ | --- | -------- | ----- | ------ | ---- | ---------- | --- | --- | --- | --- | --- |
| tion   | (4.36) | is invariant |     | in       | time: |        |      |            |     |     |     |     |     |
|        |        |              |     |          |       |        | ∂ (u | 0.         |     |     |     |     |     |
|        |        |              |     |          |       |        | τ    | Ψ ) =      |     |     |     |     |     |
| Proof. |        | The velocity |     | equation | from  | (5.31) | can  | be written |     |     |     |     |     |
|        |        |              |     |          |       |        |      | p+p        |     | c   |     |     |     |
b 0
|     | ∂ τ u,Ψ   | h           | κ 0       | p div(Ψ)dx+κ |     | 0   | Ψ   | n κ 0    |     | (u b | n   | u n) dΓ  |     |
| --- | --------- | ----------- | --------- | ------------ | --- | --- | --- | -------- | --- | ---- | --- | -------- | --- |
|     | (cid:104) | (cid:105) − |           |              |     |     |     | · 2      | −   | 2    | · − | ·        |     |
|     |           |             | (cid:90)Ω |              |     |     | ∂Ω  | (cid:18) |     |      |     | (cid:19) |     |
(cid:90)
(5.87)
|     |     |     |     |     |     | d c | divu | div(Ψ)dx. |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | ---- | --------- | --- | --- | --- | --- | --- |
|     |     |     |     |     | =   | 2 0 |      |           |     |     |     |     |     |
−
(cid:90)Ω
(cid:102)
So that taking in (5.87) Ψ div a test function in the kernel of the divergence operator such that
n
| Ψ   | div | ∂Ω = 0 | we obtain, |     | because | of  | Definition | 4.4.1 |     |     |     |     |     |
| --- | --- | ------ | ---------- | --- | ------- | --- | ---------- | ----- | --- | --- | --- | --- | --- |
·
|
|     |     |     |     |     |     |     | ∂ τ u,Ψ   | h = 0.        |     |     |     |     | (5.88) |
| --- | --- | --- | --- | --- | --- | --- | --------- | ------------- | --- | --- | --- | --- | ------ |
|     |     |     |     |     |     |     | (cid:104) | div (cid:105) |     |     |     |     |        |

| CHAPTER   | 5.      | DEVELOPMENT |     |        | OF A | CLASS  | OF LONG |     | TIME CONSISTENT |     |
| --------- | ------- | ----------- | --- | ------ | ---- | ------ | ------- | --- | --------------- | --- |
| STAGGERED | SCHEMES |             |     | ON THE | WAVE | SYSTEM |         |     |                 | 129 |
Nowweintroduce0 RT1(Ω)thefieldequalto0everywhere. Let0 from(4.36)itsdivergence
Ψ
∈
free part such that 0 n = u n on the boundary. This field is time-independent so that
|                    |            | Ψ      |             | b    |             |           |           |      |             |        |
| ------------------ | ---------- | ------ | ----------- | ---- | ----------- | --------- | --------- | ---- | ----------- | ------ |
|                    |            |        | ·           | ·    |             |           |           |      |             |        |
| (5.88) is          | equivalent | to     |             |      |             |           |           |      |             |        |
|                    |            |        |             |      | ∂ τ (u      | 0 Ψ ),Ψ   | div h     | = 0, |             | (5.89) |
|                    |            |        |             |      | (cid:104)   | −         | (cid:105) |      |             |        |
| but by Proposition |            | 4.4.2, |             |      |             |           |           |      |             |        |
|                    |            |        | ∂           | τ (u | 0 Ψ ),Ψ div | h =       | ∂ τ (u Ψ  | 0 Ψ  | ),Ψ div h , | (5.90) |
|                    |            |        | (cid:104)   | −    |             | (cid:105) | (cid:104) | −    | (cid:105)   |        |
| which leads,       | because    | of     | Proposition |      | 4.4.2,      | to        |           |      |             |        |
|                    |            |        |             |      | ∂           | (u        | 0 0.      |      |             |        |
|                    |            |        |             |      | τ           | Ψ         | Ψ ) =     |      |             | (5.91) |
−
| By time-independency |           |           | of 0 | Ψ we obtain   | the       | result     |              |       |       |     |
| -------------------- | --------- | --------- | ---- | ------------- | --------- | ---------- | ------------ | ----- | ----- | --- |
| In the               | following | we        | use  | the following |           | denotation |              |       |       |     |
| Definition           | 5.4.1     | (Relative |      | state).       | We define |            | the relative | state | as:   |     |
|                      |           |           |      | prel          |           |            | urel         |       |       |     |
|                      |           |           |      |               | := p p    | b and      | :=           | u (u  | 0 ) Ψ |     |
|                      |           |           |      |               | −         |            |              | −     |       |     |
and
p p
|     |     |     |     |     | Urel |          | b    |          |     |     |
| --- | --- | --- | --- | --- | ---- | -------- | ---- | -------- | --- | --- |
|     |     |     |     |     | :=   |          | −    |          |     |     |
|     |     |     |     |     |      | u        | (u ) |          |     |     |
|     |     |     |     |     |      |          | 0    | Ψ        |     |     |
|     |     |     |     |     |      | (cid:18) | −    | (cid:19) |     |     |
Lemma 5.4.2 (Flux annihilation). With u σ given by Definition 5.2.1 and σ [[p]] σ by Defini-
| |
| tion 5.2.3. | Then, | we have |     |     |     |     |     |     |     |     |
| ----------- | ----- | ------- | --- | --- | --- | --- | --- | --- | --- | --- |
(cid:99)
|     |     |     |            | prel      |       |      |            |         | urel |     |
| --- | --- | --- | ---------- | --------- | ----- | ---- | ---------- | ------- | ---- | --- |
|     |     |     |            |           | σ ε   | (σ)u | +          | σ [[p]] | = 0. |     |
|     |     |     |            | K         | | | K |      | σ          | | |     | σ σ  |     |
|     |     |     | K          | σ         | ∂K    |      | σ          |         |      |     |
|     |     |     | (cid:88)∈C | (cid:88)∈ |       |      | (cid:88)∈F |         |      |     |
(cid:99)

CHAPTER 5. DEVELOPMENT OF A CLASS OF LONG TIME CONSISTENT
STAGGERED SCHEMES ON THE WAVE SYSTEM 130
Proof. First using the definition of u we expand the following expression
σ
prel σ ε (σ (cid:99) )u = prel σ ε (σ)u
K K σ K K σ
| | | |
(cid:34)
K σ ∂K K σ ∂K int
(cid:88)∈C (cid:88)∈ (cid:88)∈C ∈ (cid:88)∩F
(cid:99) + σ prel u σ +(u b ) σ
| |
Kσ
2
(cid:35)
σ ∂K b
∈(cid:88)∩F
= preldiv(u)dx σ prelu
− | |
Kσ σ
K
K (cid:90) σ b
(cid:88)∈C (cid:88)∈F
(5.92)
u +(u )
+ σ prel σ b σ
| |
Kσ
2
σ b
(cid:88)∈F
= preldiv(u)dx
K
K (cid:90)
(cid:88)∈C
u +(u )
+ σ prel − σ b σ .
| |
Kσ
2
σ b
(cid:88)∈F
Using K ,div((u ) ) = 0 we have K ,div(urel) = div(u) , so this yields in (5.92)
0 Ψ K K K
∀ ∈ C ∀ ∈ C
u +(u )
prel σ ε (σ)u = preldiv(urel)dx+ σ prel − σ b σ . (5.93)
K
| |
K σ
| |
Kσ
2
K
K σ ∂K K (cid:90) σ b
(cid:88)∈C (cid:88)∈ (cid:88)∈C (cid:88)∈F
(cid:99)
On the other hand, by Definition 5.2.3 of [[p]] , we have
σ
p p
σ [[p]] urel= σ [[p]] urel + σ urel b − Kσ
σ σ σ σ σ
| | | | | | 2
σ σ int σ b
(cid:88)∈F ∈(cid:88)F (cid:88)∈F
= σ [[prel]] urel σ urel(p p )
| |
σ σ
− | |
σ Kσ
−
b
σ int σ b
∈(cid:88)F (cid:88)∈F
(5.94)
p p
+ σ urel Kσ − b
σ
| | 2
σ b
(cid:88)∈F
p p
= preldiv(urel)dx+ σ urel Kσ − b .
σ
− | | 2
K
K (cid:90) σ b
(cid:88)∈C (cid:88)∈F

CHAPTER 5. DEVELOPMENT OF A CLASS OF LONG TIME CONSISTENT
STAGGERED SCHEMES ON THE WAVE SYSTEM 131
Summing the resulting expression of (5.94) and (5.93) we have
u +(u )
prel σ ε (σ)u + σ [[p]] urel= σ prel − σ b σ
K
| |
K σ
| |
σ σ
| |
Kσ
2
K σ ∂K σ σ b
(cid:88)∈C (cid:88)∈ (cid:88)∈F (cid:88)∈F
(cid:99)
p p
+ σ urel Kσ − b
σ
| | 2
σ b
(cid:88)∈F
urel
= σ prel
− | |
Kσ
2
σ b
(cid:88)∈F
prel
+ σ urel Kσ
σ
| | 2
σ b
(cid:88)∈F
=0.
(cid:94)
Lemma 5.4.3 (Semi-discreteRelativeenergy). Letu , σ [[p]] , σ [[div(urel)]] and σ [[u n]]
σ σ σ σ
| | | | | | ·
given by, respectively, Definition 5.2.1, Definition 5.2.3, Definition 5.2.5, Definition 5.2.4.
Then the relative energy verifies the following equat(cid:99)ion:
∂ E(Urel) = d c 0 ρ 0 κ 0 σ [[prel]]2 d c ∂K di (cid:94) v(urel) 2
τ − 1 2 | | σ − 2 0 | | K
σ int K
∈(cid:88)F (cid:88)∈C
(5.95)
d c ρ κ c
1 0 0 0 σ (p p )2 0 σ (urel)2 c σ (urel)2.
− 2 | |
Kσ
−
b
− 2 | |
σ
−
0
| |
σ
σ ∈Fi b (cid:88)nlet/outlet σ ∈Fi b (cid:88)nlet/outlet σ ∈ (cid:88) Fw b all
Proof. We multiply the pressure equation by ρ κ prel = ρ κ (p p ) of (5.31)
0 0 K 0 0 K − b
1
K ∂ (p )ρ κ prel + σ ε (σ)u ρ κ prel = d c σ ε (σ)[[p]] ρ κ prel. (5.96)
| | τ K 0 0 K ρ | | K σ 0 0 K 1 0 | | K σ 0 0 K
0
σ ∂K σ ∂K
(cid:88)∈ (cid:88)∈
(cid:99)
By stationarity of p , ∂ (p ) = ∂ (p p ), so that summing on K (5.96) becomes:
b τ K τ K b
− ∈ C
∂ ( prel 2 )
K ρ κ τ K +κ σ ε (σ)u prel = d c ρ κ σ ε (σ)[[p]] prel,
0 0 0 K σ K 1 0 0 0 K σ K
| | 2 | | | |
K (cid:0) (cid:1) K σ ∂K K σ ∂K
(cid:88)∈C (cid:88)∈C (cid:88)∈ (cid:88)∈C (cid:88)∈
(cid:99)
(5.97)

CHAPTER 5. DEVELOPMENT OF A CLASS OF LONG TIME CONSISTENT
STAGGERED SCHEMES ON THE WAVE SYSTEM 132
p is constant so we have σ [[p]] = σ [[prel]] which yields that (5.97) becomes:
b σ σ
| | | |
∂ ( prel 2 )
K ρ κ τ K +κ σ ε (σ)u prel = d c ρ κ σ ε (σ)[[prel]] prel.
0 0 0 K σ K 1 0 0 0 K σ K
| | 2 | | | |
K (cid:0) (cid:1) K σ ∂K K σ ∂K
(cid:88)∈C (cid:88)∈C (cid:88)∈ (cid:88)∈C (cid:88)∈
(cid:99)
(5.98)
Similarly, we multiply the velocity equation of (5.31) by urel
σ
D ∂ u urel +κ σ [[p]] urel = c 0 σ [[u n]] urel +d c σ [[d (cid:94) iv(u)]] urel. (5.99)
σ τ σ σ 0 σ σ σ σ 2 0 σ σ
| | | | − 2 | | · | |
By stationarity of (u ) , ∂ (u ) = ∂ (u ((u ) )) ), so that summing (5.99) on σ ,
0 Ψ τ σ τ σ 0 Ψ σ
− ∈ F
D ∂ τ ( ur σ el 2 ) +κ σ [[p]] urel = c 0 σ [[u n]] urel +d c σ [[d (cid:94) iv(u)]] urel.
σ 0 σ σ σ σ 2 0 σ σ
| | 2 | | − 2 | | · | |
σ (cid:0) (cid:1) σ σ b σ
(cid:88)∈F (cid:88)∈F (cid:88)∈F (cid:88)∈F
(5.100)
(cid:94) (cid:94)
By definition div((u ) ) = 0 so div(u) = div(urel) which gives thanks to (5.100)
0 Ψ
D ∂ τ ( ur σ el 2 ) +κ σ [[p]] urel = c 0 σ [[u n]] urel +d c σ [[di (cid:94) v(urel)]] urel.
σ 0 σ σ σ σ 2 0 σ σ
| | 2 | | − 2 | | · | |
σ (cid:0) (cid:1) σ σ b σ
(cid:88)∈F (cid:88)∈F (cid:88)∈F (cid:88)∈F
(5.101)
Finally summing (5.98) and (5.101) using Lemma 5.4.2 we have:
c
∂ E(Urel) = d c ρ κ σ ε (σ)[[prel]] prel 0 σ [[u n]] urel
τ 1 0 0 0 K σ K σ σ
| | − 2 | | ·
K σ ∂K σ b
(cid:88)∈C (cid:88)∈ (cid:88)∈F
(5.102)
(cid:94) (cid:94)
+d c σ [[div(urel)]] urel d c σ div(u) urel,
2 0 | | σ σ − 2 0 | | Kσ σ
σ int σ b
∈(cid:88)F (cid:88)∈F
since urel := u (u ) for σ an inlet/oulet boundary face and urel := u for a wall boundary
σ σ b σ σ σ
−
face we have that, using the definition of the velocity jump on the boundary ;
[[u n]] = u (u ) for inlet/outlet boundary faces,
σ σ b σ
· −
and
[[u n]] = 2u for σ b ,
· σ σ ∈ Fwall

CHAPTER 5. DEVELOPMENT OF A CLASS OF LONG TIME CONSISTENT
STAGGERED SCHEMES ON THE WAVE SYSTEM 133
so that
c c
0 σ [[u n]] urel = 0 σ (urel)2 c σ (urel)2. (5.103)
σ σ σ 0 σ
− 2 | | · − 2 | | − | |
σ (cid:88)∈F b σ ∈Fi b n(cid:88)let/outlet σ ∈ (cid:88) Fw b all
Finally, using the integration by parts Lemma 5.3.1 iii) on the velocity divergence integral, we
get
(cid:94) (cid:94) (cid:94) 2
d c σ [[div(urel)]] urel d c σ div(u) = d c div(urel) , (5.104)
2 0 | | σ σ − 2 0 | | Kσ − 2 0 K
σ int σ b K
∈(cid:88)F (cid:88)∈F (cid:88)∈C
and mimicking ii) from Lemma 5.3.1 on the pressure Laplacian
d c ρ κ σ ε (σ)[[prel]] prel=d c ρ κ σ ε (σ)[[prel]] prel
1 0 0 0 K σ K 1 0 0 0 K σ K
| | | |
K σ ∂K K σ ∂K, int
(cid:88)∈C (cid:88)∈ (cid:88)∈C ∈ (cid:88)∩F
d c ρ κ
1 0 0 0 σ (p p )2
− 2 | |
b
−
Kσ
σ ∈Fi b n(cid:88)let/outlet
(5.105)
c ρ κ
= d 0 0 0 σ [[prel]]2
1 σ
− 2 | |
σ int
∈(cid:88)F
d c ρ κ
1 0 0 0 σ (p p )2.
− 2 | |
b
−
Kσ
σ ∈Fi b n(cid:88)let/outlet
Gathering (5.105), (5.104) and (5.103) in (5.102) we obtain the result.
Lastly, we have
Lemma 5.4.4. Suppose that v RT1(Ω) is such that
∈
v = div ∗ϕ for some ϕ Φ := dG0(Ω) dG0(∂Ω)
− ∈ ×
and (cid:0) (cid:1)
K div(v) = 0, v n = 0 (5.106)
K ∂Ω
∀ ∈ C · |
Then
v = 0
Proof. If v = div ∗ϕ then by Definition 4.4.4,
−
(cid:0) (cid:1)
ω RT1(Ω) v,ω
h
= div ∗ϕ,ω
h
= div(ω)ϕdx+ ϕω ndΓ. (5.107)
∀ ∈ (cid:104) (cid:105) (cid:104) − (cid:105) − ·
(cid:90)Ω (cid:90) ∂Ω
(cid:0) (cid:1)
Then if div(ω) = 0 such that ω n = 0, (5.107) yields
∂Ω
· |
v ω RT1(Ω), div(ω) = 0, ω n = 0 . (5.108)
∂Ω
⊥ { ∈ · }

CHAPTER 5. DEVELOPMENT OF A CLASS OF LONG TIME CONSISTENT
STAGGERED SCHEMES ON THE WAVE SYSTEM 134
But (5.108) combined with hypothesis (5.106) imply together that
v ω RT1(Ω), div(ω) = 0, ω n = 0 ω RT1(Ω), div(ω) = 0, ω n = 0 := 0 ,
∂Ω ∂Ω ⊥
∈ { ∈ · | }∩{ ∈ · | } { }
so that v = 0.
We can now combine these intermediary results in order to prove Theorem 5.4.1:
proof of Theorem 5.4.1. The energy dissipation implies
E(Urel(τ)) c R , (5.109)
τ−→+ ∈
∗+
→ ∞
which implies that the sequence is bounded. Also, E(.) is a norm , therefore, since the spaces
are finite dimensional, we can extract a subsquence that converges to a limit. Combining this
with Lemma 5.4.3, it implies that in the limit of this subsequence τ + we have:
→ ∞
d c 0 ρ 0 κ 0 σ [[prel]]2 +d c ∂K di (cid:94) v(urel) 2
1 2 | | σ 2 0 | | K
σ int K
∈(cid:88)F (cid:88)∈C
d c ρ κ c
+ 1 0 0 0 σ (p p )2+ 0 σ (urel)2+c σ (urel)2 = 0.
2 | |
Kσ
−
b
2 | |
σ 0
| |
σ
σ ∈Fi b n(cid:88)let/outlet σ ∈Fi b n(cid:88)let/outlet σ ∈ (cid:88) Fw b all
(5.110)
Since each term is non-negative, they are in fact each equal to 0. Leading to urel n = 0
∂Ω
and prel = 0. It is important to note that because of the stationarity preserving p · rop | erty of
∂Ω
the sc|heme, we have by Lemma 5.4.1
σ ∂ (u ) = 0, (5.111)
τ Ψ σ
∀ ∈ F
which leads to urel = u (u ) = u +u (u ) = u +(u ) (u ) = u . Then,
0 Ψ ϕ Ψ 0 Ψ ϕ 0 Ψ 0 Ψ ϕ
− − −
(cid:94) 2
Case d > 0, d > 0: (5.110) leads to ∂K div(urel) = 0 yielding to div(urel) = 0 .
1 2 | | K
K
Then, because urel n
∂Ω
= 0 and urel = u (cid:88)
ϕ
∈C = div ∗ϕ then Lemma 5.4.4 yields urel = 0.
(5.110) yields prel = · 0, | giving the conclusion. −
(cid:0) (cid:1)

CHAPTER 5. DEVELOPMENT OF A CLASS OF LONG TIME CONSISTENT
STAGGERED SCHEMES ON THE WAVE SYSTEM 135
Case d > 0, d = 0: In this case we write the system verified by Urel:
1 2
K 1
K ∂ prel + | | div(urel) = σ ε (σ)urel +d c σ [[prel]] ,
| | τ K ρ K 2ρ | | K σ 1 0 | | σ
 0 0 σ K ∂Ω σ ∂K

∈(cid:88)∩ (cid:88)∈



D ∂ urel +κ σ [[prel]] = 0 σ int.
σ τ σ 0 σ
| | | | ∀ ∈ F



 (cid:94)
Nowreplacingthegradientjump[[p]] byurel inpointi)ofLemma5.3.2andreplacingdiv(u)
σ σ K
by prel in point ii) of Lemma 5.3.2 we can pass to the limit (of the extracted subsequence) in
K
the system thanks to the relative energy dissipation. Then, the derivative in time vanishes (by
time independency of the limit ). We also know from (5.110) that urel n = 0 and prel = 0,
so that all terms depending on the boundary disappear, combined with · d | ∂Ω = 0 it yields| ∂Ω
2
K
| | div(urel) = d c σ [[prel]] K ,
K 1 0 σ
ρ | | ∀ ∈ C
 0 σ ∂K

(cid:88)∈



κ σ [[prel]] = 0 σ int.
0 σ
| | ∀ ∈ F



The second equationimplies that prel is constant in the limit and is thus equal to 0. As a
consequence, the pressure Laplacian on the first equation vanishes leading to a divergence free
urel in the limit. Recalling that urel = u , urel n = 0 and K div(urel) = 0 we
ϕ ∂Ω K
obtain by Lemma 5.4.4 urel = 0. · | ∀ ∈ C
Then for both of these cases we have in the limit of the extracted subsequence E(Urel) = 0,
then, because (5.109) stands, we can show by a classical compacity argument that the whole
sequence converges the this limit.
Remark 5.4.2. The fully discrete setting follows from the same arguments up to the hypo-
thesis that the time integration yields energy dissipation. Indeed, the proofs of Lemma 5.4.1,
Lemma 5.4.2 are almost identical in the fully discrete setting, no matter the choice of time
stepping. In parallel, if the time integration is Explicit or ImEx the scheme yields energy dis-
sipation under CFL condition, as shown in section 5.3 and unconditionnally in Implicit time
integration. In this fully discrete setting, the relative energy balance is very similar to (5.95),
up to terms of the form Un+1 Un 2, which vanishes in the limit n + . Since only
(cid:107) rel − rel(cid:107) −→ ∞
the limit energy balance matters, we see that under energy dissipation conditions the proof of
Theorem 5.4.1 stands in fully discrete.
Remark 5.4.3. It is not detailed here because it is not useful in our case but we can actually
prove the case d = 0,d > 0 easily since in the limit it will yield div(urel) = 0 so by the sta-
1 2
tionarity preservation qualities of the grad-div diffusion, we will obtain also a constant relative
pressure equal to 0.
Remark 5.4.4. Theorem 5.4.1 shows that we have in fact commutation between the long time
limit and the limit when the mesh size tends to 0.

| CHAPTER   |            | 5. DEVELOPMENT |         | OF A        | CLASS OF | LONG      | TIME | CONSISTENT |         |     |
| --------- | ---------- | -------------- | ------- | ----------- | -------- | --------- | ---- | ---------- | ------- | --- |
| STAGGERED |            | SCHEMES        | ON      | THE WAVE    | SYSTEM   |           |      |            |         | 136 |
| 5.5       | Discussion |                | on some | preexisting |          | staggered |      | schemes    | through |     |
|           | low        | Mach           | number  | asymptotics |          |           |      |            |         |     |
In the following section we will compare the schemes developped in this chapter with staggered
schemes that can be found in the literature. Since most of them are presented on Euler’ system
(or Navier-Stokes), we will focus on their asymptotic expansion in Mach number in order to
recover our formalism. It is fundamental to recall that most of the staggered finite volume
schemes found are based on vectorial Crouzeix-Raviart [44] (simplicial meshes) or vectorial
Rannacher-Turek [63] (quad/hexa meshes) as the discretization basis of the staggered velocity.
These finite element spaces consist in fields that are affine on each component and have the
global constraint of being continuous in the middle point of each face of the mesh. Normal con-
tinuity is thus not ensured in every case; it is known that they are neither H(div) conforming
−
H1
nor conforming. In this context, we are not aware of the existence of a discrete Hodge-
−
Helmholtz Decomposition for these finite elements spaces and we are much less able to charac-
terize the corresponding discrete (non-conforming) de Rham complex resembling the N´ed´elec-
| Raviart-Thomas |     | complex | of our | discretization. |     |     |     |     |     |     |
| -------------- | --- | ------- | ------ | --------------- | --- | --- | --- | --- | --- | --- |
Besides,thediscreteconvergenceinlongtimeisaconsequenceof,inparticular,thepreserva-
tion of Lemma 2.3.1 at the discrete level; in this aim stationarity preserving numerical diffusion
isrequiredaswellasHodge-orthogonalityofthevelocityunknownspace. Asalreadydiscussed,
wearenotabletoreplicatethesecondpointforvelocitiesinCrouzeix-Raviart/Rannacher-Turek
finite element space. Naturally, we are therefore not able to characterize, in general, the long
| time | consistency | with | this choice | of velocity | approximation |     | space. |     |     |     |
| ---- | ----------- | ---- | ----------- | ----------- | ------------- | --- | ------ | --- | --- | --- |
However, we emphasized two other fundamental points for accurate discrete long time be-
haviour to occur: energy dissipation and stationarity preserving diffusion; we can investigate
thesepropertiesonthelowMachnumberasymptoticsoftheaforementionedstaggeredschemes.
As a short reminder we detailed in Appendix C the low Mach number asymptotic analysis for
| two | of the schemes, |     | the others | will follow | equivalently. |     |     |     |     |     |
| --- | --------------- | --- | ---------- | ----------- | ------------- | --- | --- | --- | --- | --- |
The multiple dimensions case is far more interesting than the one dimensional case because; on
theonehandwestillhavedifferenttimeintegrationstrategieswhereasontheotherhandmulti-
dimensionality will allow for much richer spatial diffusive operators. With this in mind, we will
arrange staggered schemes found in three categories depending on the time integration but
without necessarily any details on the upwinding. So that we have the schemes asymptotically
| consistent | with: |     |     |     |     |     |     |     |     |     |
| ---------- | ----- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
•
Explicit discretizations: The low Mach number asymptotic analysis on the explicit
scheme of [69] proposed on Shallow Water equations gives (see Proposition C.1.1 in Ap-

CHAPTER 5. DEVELOPMENT OF A CLASS OF LONG TIME CONSISTENT
STAGGERED SCHEMES ON THE WAVE SYSTEM 137
pendix C)
(h˜(1) )n+1 (h˜(1) )n h˜(0) h˜(0) σ˜
K − K + σ˜ (u(0))n n = γ∆τ σ˜ | | (h˜(1) )n (h˜(1) )n ,
 

∆τ | K˜ | σ (cid:88)∈ ∂K | | σ · K,σ | K˜ | σ (cid:88)∈ ∂K | | | D˜ σ | (cid:16) L − K (cid:17)
  (cid:101)

       (u( σ 0) )n+ ∆ 1 − τ (u( σ 0) )n + D | ˜ σ˜ | (h˜( L 1) )n − (h˜( K 1) )n n K,σ =
σ
| | (cid:104) (cid:105)
(cid:101) (cid:101)∂L ∂K
        2α∆τh˜(0) (cid:34) (cid:18) | |(cid:102) L˜ | | (cid:16) (u( σ 0) )n − (u( K 0) )n (cid:17) · n L,σ (cid:19) − (cid:18) | |(cid:103) K˜ | | (cid:16) (u( σ 0) )n − (u( K 0) )n (cid:17) · n K,σ (cid:19) (cid:35) n K,σ .

  (cid:101) (cid:101) (cid:101) (cid:101)
 So this scheme has diffusion on both equations that depends on pressure AND velocity
terms, it nonetheless differs from our stabilization since it is not a grad-divergence oper-
ator. Indeed, from (C.2), (C.3), we are not able to conclude that if
1
div(u) := σ u(x ) n = 0 (5.112)
K σ K,σ
K | | ·
| | σ ∂K
(cid:88)∈
with x barycenter of the face σ, for all cells K implies that the diffusion vanishes. In
σ
fact, if we restrict to a 2d Cartesian grid on a periodic domain we see that for a fixed face
σ
∂L ∂K
| | (u(0))n (u(0) )n n | | (u(0))n (u(0) )n n n
(cid:34)
L˜ σ − K · L,σ − K˜ σ − K · K,σ
(cid:35)
K,σ
(cid:18) |(cid:102)| (cid:16) (cid:17) (cid:19) (cid:18) |(cid:103)| (cid:16) (cid:17) (cid:19)
(cid:101) (cid:101) (cid:101) (cid:101)
4
= (u(0))n (u(0) )n n (u(0))n (u(0) )n n n
∆x σ − K · L,σ − σ − K · K,σ K,σ
(cid:34) (cid:35)
(cid:18) (cid:19) (cid:18) (cid:19)
(cid:16) (cid:17) (cid:16) (cid:17)
(cid:101) (cid:101) (cid:101) (cid:101)
4
= (u(0))n (u(0) )n ( n ) (u(0))n (u(0) )n n n
∆x σ − K · − K,σ − σ − K · K,σ K,σ
(cid:34) (cid:35)
(cid:18) (cid:19) (cid:18) (cid:19)
(cid:16) (cid:17) (cid:16) (cid:17)
(cid:101) (cid:101) (cid:101) (cid:101)
4
= (u(0) )n n 2(u(0))n n +(u(0) )n n n .
∆x L · K,σ − σ · K,σ K · K,σ K,σ
(cid:34) (cid:35)
While the exact expres(cid:101)sion of this opera(cid:101)tor depends on t(cid:101)he reconstruction in the cell u˜ ,
K
we observe that it is formally very close to the normal component of a discrete Laplacian;
is not evident that this diffusion will preserve discrete divergence free velocities (5.112)
• Implicit (fully-centred) discretizations: we mention the implicit schemes of [49, 66,
81], which are, by their implicit nature, energy dissipative and since no spatial diffusive
terms are needed to obtain this property, they do not suffer from the loss of stationarity
preservation.
• ImEx discretizations: [65] results on the wave system in a fully-centred ImEx scheme
whichisveryreminescentoftheImExintegration(5.53)butdiffersfromitbytheamount

| CHAPTER   |     | 5. DEVELOPMENT |     |     | OF A     | CLASS  | OF LONG | TIME CONSISTENT |     |
| --------- | --- | -------------- | --- | --- | -------- | ------ | ------- | --------------- | --- |
| STAGGERED |     | SCHEMES        |     | ON  | THE WAVE | SYSTEM |         |                 | 138 |
of diffusive terms: there is no discrete diffusion on the pression equation. The one dimen-
sionalstudyinchapter3showsthatsuchschemeswithonecentredtermareinsufficiently
diffusive and authorize the energy to increase for some particular modes. Since this class
of schemes is not energy-dissipative, which is a necessary condition in order complete our
study on the long time behaviour, we are, as a consequence, not able to hint on its low
|     | Mach      | number | accuracy. |         |     |     |     |     |     |
| --- | --------- | ------ | --------- | ------- | --- | --- | --- | --- | --- |
| 5.6 | Numerical |        |           | Results |     |     |     |     |     |
Inthissectionweexplorethenumericallongtimebehaviourofthestaggeredschemesintroduced
SolverLab.
in the chapter. The numerical scheme is implemented in the C++ code This
code is originally designed for cell based finite volume solvers, some node based finite element
methods can also be found and most importantly it offers a huge flexibility with regards to the
managment of meshes; quadrangular, triangular or even generic meshes are easily generated on
non-trivialgeometries. Thisworkmarkstheintroductionoffacebasedmethodsintheplatform;
in this aim a new class for the staggered discretization of the wave system was implemented
with the addition of discrete differential operators such as the discrete staggered gradient,
divergence but also the grad-div operator for Raviart-Thomas elements. Also, Raviart-Thomas
basis functions were coded in order to compute the mass-lumping and an interpolation from
faces to cells formula for post-processing as been coded. It is given by integrating on the cell
| the Raviart-Thomas |     |     | unknown |     | u RT1(Ω): |     |     |     |     |
| ------------------ | --- | --- | ------- | --- | --------- | --- | --- | --- | --- |
∈
|     |     |     |        | 1    |     | 1   |      |                     |         |
| --- | --- | --- | ------ | ---- | --- | --- | ---- | ------------------- | ------- |
|     |     |     | (u h ) | K := | udx | =   | σ (x | σ x K )ε K (σ)u σ , | (5.113) |
|     |     |     |        | K    |     | K   | | |  | −                   |         |
K
|     |     |     |     | |   | | (cid:90) | | | σ | ∂K  |     |     |
| --- | --- | --- | --- | --- | ---------- | ----- | --- | --- | --- |
(cid:88)∈
(cid:102)
where x σ is the barycenter of the face σ and x K the center of mass of K. Finally, adapted
flux formulations of boundary terms were implemented. As far as the time integration goes,
the staggered scheme for the wave system is available in Euler explicit method and centred
implicit method.
The numerical simulations are ran with d = 1/2, CFL=0.6 and ρ = κ = 1, and the
0 0
Un+1 Un 10
numerical results correspond to steady states obtain with a precision of 10
| − | 2 ≈ −
.
| with | 2 is | the 2 | Euclidian | norm. |     |     |     |     |     |
| ---- | ---- | ----- | --------- | ----- | --- | --- | --- | --- | --- |
|      | | |  | −     |           |       |     |     |     |     |     |
In order to test the accuracy of the solver developed in this chapter, we tested (5.48)
on the following test case: Ω the domain in which the equations stand is an annulus, in polar
R2.
coordinates, Ω = [r ,r ] [0,2π] We simulate the scattering of a wave in Ω with
|     |     |     | 0   | 1   |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
|     |     |     |     | ×   | ⊂   |     |     |     |     |
r = 0.8 and r = 6 where wall boundary conditions (5.4) on the inner circle are imposed
| 0   |     | 1   |     |     |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
weakly on r = r 0 and weakly inlet/outlet boundary conditions (5.2) are imposed on the outer
|     | r   | r   |     | p   | u   | (1,0)t |     |     |     |
| --- | --- | --- | --- | --- | --- | ------ | --- | --- | --- |
circle on = 1 , where b = 0 and b = are given pressure and velocity. Conveniently,

| CHAPTER   | 5.   | DEVELOPMENT |           | OF A          | CLASS  | OF LONG | TIME | CONSISTENT |     |
| --------- | ---- | ----------- | --------- | ------------- | ------ | ------- | ---- | ---------- | --- |
| STAGGERED |      | SCHEMES     | ON        | THE WAVE      | SYSTEM |         |      |            | 139 |
| the long  | time | limit for   | this test | case is known | and    | is:     |      |            |     |
|           |      |             |           | p             | (r,θ)  |         |      |            |     |
|           |      |             |           |               | exact  | = 0     |      |            |     |
|           |      |            |           |               |        |         |      |           |     |
r2
0 cos(2θ)
|     |     |    |         |     |     | 1     |     |  . | (5.114) |
| --- | --- | --- | ------- | --- | --- | ----- | --- | --- | ------- |
|     |     |     |         |     | r2  | − r 2 |     |     |         |
|     |     |    |         |     |     |       |     |    |         |
|     |     |    | u (r,θ) | =   | 1  |       |     |   |         |
|     |     |    | exact   | r2  | r2  |       |     |    |         |
2
|     |     |    |     | 1   | − 0  | r         |     |   |     |
| --- | --- | --- | --- | --- | ----- | --------- | --- | --- | --- |
|     |     |    |     |     |       | 0 sin(2θ) |     |    |     |
|     |     |    |     |     |      |           |     |   |     |
|     |     |     |     |     |      | − r 2     |     |    |     |
|     |     |    |     |     |      |           |     |   |     |
|     |     |    |     |     |       |           |     |    |     |
This test case illustrates Theorem 5.4.1 in the sense that, under energy dissipation conditions,
the choice of velocity staggering and relevant grad-div diffusion operator gathered together
yield the convergence to the accurate long time limit as detailed in subsection 5.6.1. These
conditionsarestatedtobesufficientinTheorem5.4.1, howeverinsubsection5.6.2weconstruct
an example of a diffusive operator that does not preserve divergence-free velocity and challenge
itsbehaviouronthesametest: theresultsindicatesthatthestationaritypreservationcondition
| might be | in fact   | necessary. |      |     |     |     |     |     |     |
| -------- | --------- | ---------- | ---- | --- | --- | --- | --- | --- | --- |
| 5.6.1    | Numerical | long       | time |     |     |     |     |     |     |
Qualitative results on the long time limit: the velocity isocontours
In Figure 5.2 and Figure 5.1, we compare the iso-values of the 2-Euclidean norms of: the face-
to-cell interpolation of long time numerical limit obtained on a mesh composed of unstructured
20 64 squares (Figure 5.2), and the face-to-cell interpolation of the known analytical long
×
| time limit(5.114) |     | on the | same mesh | (Figure | 5.1). |     |     |     |     |
| ----------------- | --- | ------ | --------- | ------- | ----- | --- | --- | --- | --- |
Figure 5.2 shows that the numerical solution captured with the numerical scheme proposed
in this chapter is qualitatively consistent with the exact long time limit Figure 5.1, which is,
as we will see in subsection 5.6.2, not evident (see also [39, Chapter 5 p 103] for examples
of numerical schemes that do not capture this tendency on the iso-contours). Also, results in
explicit and implicit time integration are, for the eye assessment, identical, confirming that,
up to the appropriate choice of stationarity preserving diffusion the problem of capturing the
| accurate | long | time stationary | state | is essentially |     | spatial. |     |     |     |
| -------- | ---- | --------------- | ----- | -------------- | --- | -------- | --- | --- | --- |

| CHAPTER   | 5. DEVELOPMENT | OF A                  | CLASS OF LONG | TIME CONSISTENT |     |
| --------- | -------------- | --------------------- | ------------- | --------------- | --- |
| STAGGERED | SCHEMES        | ON THE WAVE           | SYSTEM        |                 | 140 |
|           | Figure         | 5.1: Norm isocontours | of the exact  | long time limit |     |

| CHAPTER   | 5. DEVELOPMENT |      |        | OF A        | CLASS  | OF       | LONG TIME      |     | CONSISTENT |        |     |
| --------- | -------------- | ---- | ------ | ----------- | ------ | -------- | -------------- | --- | ---------- | ------ | --- |
| STAGGERED | SCHEMES        |      | ON     | THE WAVE    | SYSTEM |          |                |     |            |        | 141 |
|           | (a) Explicit   |      | scheme |             |        |          |                | (b) | Implicit   | scheme |     |
|           | Figure         | 5.2: | Norm   | isocontours | of     | the long | time numerical |     | solutions  |        |     |
Quantitative results on the long time limit: the mesh convergence
Using the exact long time limit (5.114), we are able to perform a convergence analysis on the
stationary state obtained with the explicit scheme: the analysis is performed on 4 meshes, of
size5 8, 5 16, 10 32, 20 64. Figure5.3showsthatthediscretizationproposedconvergesto
| ×   | ×   | ×   | ×   |     |     |     |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
the exact long time limit with an order slighly inferior to 1, showing that we have commutation
| between | the limit | in mesh | size | and the | long time | limit. |     |     |     |     |     |
| ------- | --------- | ------- | ---- | ------- | --------- | ------ | --- | --- | --- | --- | --- |
Remark 5.6.1. On the cylinder test, the quadrangular mesh has the particularity that it is
actually not affine, meaning that for any cell K the transformation between the reference
|     | Kˆ  |     |     |     |     | ∈ C |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
element T : : K has a non-constant gradient. At first sight, an unfortunate consequence
K −→
of this is that the divergence of any Raviart-Thomas is not constant by cell, thus disrupting an
important property intrinsic to the definition of the discrete N´ed´elec-Raviart-Thomas complex:
div(RT1(Ω) dG0(Ω). Luckily, it is possible to modify the basis functions in order to recover
⊂
this fundamental property [115]: for the basis functions on the reference element, knowning the
| associated | physical | element           | K   | the modified         | basis | yields |          |     |     |     |     |
| ---------- | -------- | ----------------- | --- | -------------------- | ----- | ------ | -------- | --- | --- | --- | --- |
|            |          |                   |     |                      |       |        | 1 βxˆ(xˆ |     | 1)  |     |     |
|            |          | Ψˆmodified(xˆ,yˆ) |     | = Ψˆoriginal(xˆ,yˆ)+ |       |        |          |     |     | ,   |     |
−
|     |     |     |     |     |     |     | K γyˆ(yˆ     |     | 1)       |     |     |
| --- | --- | --- | --- | --- | --- | --- | ------------ | --- | -------- | --- | --- |
|     |     |     |     |     |     |     | | | (cid:18) | −   | (cid:19) |     |     |
where
|     |     |     |     | det(∇T |       | α+βxˆ+γyˆ. |     |     |     |     |     |
| --- | --- | --- | --- | ------ | ----- | ---------- | --- | --- | --- | --- | --- |
|     |     |     |     |        | K ) = |            |     |     |     |     |     |

| CHAPTER   |     | 5. DEVELOPMENT |     |     | OF  | A    | CLASS  | OF  | LONG | TIME         | CONSISTENT |     |     |
| --------- | --- | -------------- | --- | --- | --- | ---- | ------ | --- | ---- | ------------ | ---------- | --- | --- |
| STAGGERED |     | SCHEMES        |     | ON  | THE | WAVE | SYSTEM |     |      |              |            |     | 142 |
|           |     | Implicit,RT1   |     |     |     |      |        |     |      | Implicit,RT1 |            |     |     |
10−1
|     |     | Explicit,RT1 |     |     |     |     |     |      |     | Explicit,RT1 |     |     |     |
| --- | --- | ------------ | --- | --- | --- | --- | --- | ---- | --- | ------------ | --- | --- | --- |
| 2   |     |              |     |     |     |     |     | 2    |     |              |     |     |     |
| k   |     |              |     |     |     |     |     | k    |     |              |     |     |     |
| xe  |     |              |     |     |     |     |     | xe   |     |              |     |     |     |
| ) x |     |              |     |     |     |     |     | )    |     |              |     |     |     |
| u(  |     |              |     |     |     |     |     | u( y |     |              |     |     |     |
10−1
ssllooppee==00..99
| −    |     |     |     |     |     |     |     | −    | slope=0.8 |     |     |     |     |
| ---- | --- | --- | --- | --- | --- | --- | --- | ---- | --------- | --- | --- | --- | --- |
| h    |     |     |     |     |     |     |     | h    |           |     |     |     |     |
| )    |     |     |     |     |     |     |     | )    |           |     |     |     |     |
| u( x |     |     |     |     |     |     |     | u( y |           |     |     |     |     |
| k    |     |     |     |     |     |     |     | k    |           |     |     |     |     |
10−2
|     | 10−3 |        |        |             | 10−2 |     |      |        | 10−3  |       |                    | 10−2 |     |
| --- | ---- | ------ | ------ | ----------- | ---- | --- | ---- | ------ | ----- | ----- | ------------------ | ---- | --- |
|     |      |        |        | h           |      |     |      |        |       |       | h                  |      |     |
|     | (a)  | Error  | on the | x component |      |     |      |        | (b)   | Error | on the y component |      |     |
|     |      |        |        | −           |      |     |      |        |       |       | −                  |      |     |
|     |      | Figure | 5.3:   | Convergence |      | in  | mesh | to the | exact | long  | time solution      |      |     |
However, from a pure computational point of view, the difference in the results between the
ones obtained with the enriched basis and the others obtained with the usual basis is basically
| non-existent. |     | In fact, | if  | the basis | functions |     | (Ψ  | )   | verify |     |     |     |     |
| ------------- | --- | -------- | --- | --------- | --------- | --- | --- | --- | ------ | --- | --- | --- | --- |
|               |     |          |     |           |           |     | σ   | σ   |        |     |     |     |     |
∈F
1
|     |     |     |     | f,σ |     | θ(Ψ ) | :=  |            | Ψ n | dΓ = | δ , |     |     |
| --- | --- | --- | --- | --- | --- | ----- | --- | ---------- | --- | ---- | --- | --- | --- |
|     |     |     |     |     |     | σ     | f   | f          | σ   | f    | f,σ |     |     |
|     |     |     |     | ∀   | ∈ F |       |     | f          | ·   |      |     |     |     |
|     |     |     |     |     |     |       | |   | | (cid:90) |     |      |     |     |     |
then their actual expressions will only affect the mass-lumping, which, for the long time study,
has no effects.
Remark 5.6.2. Interestingly, the variational formulation established in order to prove the
Hodge-Decomposition is actually independent of the number of neighbours of each cell, showing
that the problem can be solvable for any type of mesh. It is thus possible to construct a vector
of u¯ R# such that (Bu¯) = ε¯ (σ) σ u¯ = 0 (B R# # the ’divergence’ matrix)
|     | F   |     |     | K   | σ   | ∂K K |     | σ   |     | C×  | F   |     |     |
| --- | --- | --- | --- | --- | --- | ---- | --- | --- | --- | --- | --- | --- | --- |
| ∈   |     |     |     |     |     |      | |   | |   |     | ∈   |     |     |     |
with u¯ = (u n) for all σ ∈b. In parallel, the numerical results confirm this remark.
|     | σ   | b   | σ   |     |            |     |     |     |     |     |     |     |     |
| --- | --- | --- | --- | --- | ---------- | --- | --- | --- | --- | --- | --- | --- | --- |
|     |     | ·   |     |     | ∈(cid:80)F |     |     |     |     |     |     |     |     |
Generally, the proof of existence of the HHD suggests that it might exists spaces that verify
these hypothesis for more general meshes. In this regard, it might be able to adapt ideas from
| the Virtual |     | Element | Method    | [116, | 117].           |     |     |            |     |           |     |     |     |
| ----------- | --- | ------- | --------- | ----- | --------------- | --- | --- | ---------- | --- | --------- | --- | --- | --- |
| 5.6.2       | On  | the     | necessity |       | of a stationary |     |     | preserving |     | diffusion |     |     |     |
In this section we explore numerically the structure preservation condition (Lemma 2.3.1)
on the numerical diffusion. On the one hand we have extensively motivated the need of
stationarity preserving diffusion operator and on the other hand we have shown that it brings
the expected numerical results. However we want to drive the point home by showing how a
different and more naive operator can degrade the properties of the scheme.

| CHAPTER   |     | 5.  | DEVELOPMENT |     |     | OF A | CLASS  | OF LONG | TIME | CONSISTENT |     |     |
| --------- | --- | --- | ----------- | --- | --- | ---- | ------ | ------- | ---- | ---------- | --- | --- |
| STAGGERED |     |     | SCHEMES     | ON  | THE | WAVE | SYSTEM |         |      |            |     | 143 |
So, a very natural diffusive operator to use in the first place would have been the vec-
torial Laplacian. At first glance, defining a Laplacian on staggered grid is not trivial, moreover
it is also not evident how to extract a relevant diffusion operator inspired by such Laplacian.
For the former, we have shown in chapter 4 that the complex yields the discrete equivalent of
| vectorial |     | Laplacian |     |     |     |             |     |     |        |     |     |     |
| --------- | --- | --------- | --- | --- | --- | ----------- | --- | --- | ------ | --- | --- | --- |
|           |     |           |     |     | ∆u  | = ∇div(u)+∇ |     |     | curlu, |     |     |     |
⊥
whichiscalledtheHodge-LaplacianoftheRaviart-ThomasspaceDefinition4.4.3. Theexact
|     |     |     | ∇ (∇ | (u),Ψ |     |     |     |     |     |     |     |     |
| --- | --- | --- | ---- | ----- | --- | --- | --- | --- | --- | --- | --- | --- |
expression of ⊥ ⊥ ) ∗ h depends on the mass-lumping procedure chosen: we derive
|     |     |     | (cid:104) |     | (cid:105) |     |     |     |     |     |     |     |
| --- | --- | --- | --------- | --- | --------- | --- | --- | --- | --- | --- | --- | --- |
rapidly a formula for the discrete masss-lumped scalar product Definition 5.2.2. We denote in
| the | following | for | any | node n | and face | σ      |         |     |        |     |     |     |
| --- | --------- | --- | --- | ------ | -------- | ------ | ------- | --- | ------ | --- | --- | --- |
|     |           |     |     |        | ε        | (n) := | sign((x | x   | ) n ), |     |     |     |
|     |           |     |     |        | σ        |        |         | n σ | ⊥σ     |     |     |     |
|     |           |     |     |        |          |        |         | −   | ·      |     |     |     |
|     | x         |     |     |        |          |        | n       | x   |        |     | σ   |     |
where n is the point associated to the node and σ the barycenter of the face and
D
σ
|     |     |     |     |     |     |     | (cid:63)σ := | | |. |     |     |     | (5.115) |
| --- | --- | --- | --- | --- | --- | --- | ------------ | ---- | --- | --- | --- | ------- |
|     |     |     |     |     |     | |   | |            | σ    |     |     |     |         |
| |
Lemma 5.6.1. Let (∇ ) the discrete curl defined in Definition 4.4.2 with the scalar product
⊥ ∗
., . given by the mass-lumped scalar product Definition 5.2.2. Let ϕ the basis function of
|     | h   |     |     |     |     |     |     |     |     | n   |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
(cid:104) (cid:105)
cG1(Ω) at the node n. Then we have in both triangles and quadrangular meshes:
|     |     |     | u   | RT1(Ω) |     | ϕ ,(∇     | ) u       | =     |     | (cid:63)f ε (n)u , |     |     |
| --- | --- | --- | --- | ------ | --- | --------- | --------- | ----- | --- | ------------------ | --- | --- |
|     |     |     |     |        |     | n         | ⊥ ∗       | L2(Ω) |     | f f                |     |     |
|     |     |     | ∀   | ∈      |     | (cid:104) | (cid:105) |       | |   | |                  |     |     |
f (n)
∈(cid:88)F
| where |     | (n) denotes |     | the set of | faces | touching | the | node n. |     |     |     |     |
| ----- | --- | ----------- | --- | ---------- | ----- | -------- | --- | ------- | --- | --- | --- | --- |
F
n
Figure 5.4: Illustration of the rotated gradient ∇ around a node and its associated volume
⊥
(cid:63)n

| CHAPTER   | 5.            | DEVELOPMENT |         |             | OF     | A CLASS   | OF     | LONG      | TIME |           | CONSISTENT |     |
| --------- | ------------- | ----------- | ------- | ----------- | ------ | --------- | ------ | --------- | ---- | --------- | ---------- | --- |
| STAGGERED |               | SCHEMES     |         | ON THE      | WAVE   |           | SYSTEM |           |      |           |            | 144 |
| Proof.    | By Definition |             | 4.4.2,  | for ϕ       | cG1(Ω) | and       | u      | RT1(Ω)    |      |           |            |     |
|           |               |             |         |             | ∈      |           |        | ∈         |      |           |            |     |
|           |               | ϕ           | n ,(∇ ⊥ | ) u L2(Ω)   | :=     | ∇ ⊥       | ϕ n ,u | h =       | u f  | ∇ ⊥ ϕ     | n ,Ψ f h , |     |
|           |               | (cid:104)   |         | ∗ (cid:105) |        | (cid:104) |        | (cid:105) |      | (cid:104) | (cid:105)  |     |
f
(cid:88)∈F
whether the mesh is composed of triangles of quads we have that the support of ϕ is the union
n
of the cells which have the node in their boundary so by Definition 5.2.2
|     |     |            | u ∇       | ϕ   | ,Ψ        | =          | u   | ∇         | ϕ ,Ψ |           |          |     |
| --- | --- | ---------- | --------- | --- | --------- | ---------- | --- | --------- | ---- | --------- | -------- | --- |
|     |     |            | f         | ⊥ n | f h       |            | f   | ⊥         | n f  | h         |          |     |
|     |     |            | (cid:104) |     | (cid:105) |            |     | (cid:104) |      | (cid:105) |          |     |
|     |     | f          |           |     |           | f          | (n) |           |      |           |          |     |
|     |     | (cid:88)∈F |           |     |           | ∈(cid:88)F |     |           |      |           |          |     |
|     |     |            |           |     |           | =          | D   | u         | (∇ ϕ | ) n       | (x ) .   |     |
|     |     |            |           |     |           |            |     | f f       | ⊥    | n f       | f        |     |
|     |     |            |           |     |           |            | |   | |         |      | ·         |          |     |
|     |     |            |           |     |           |            |     | (cid:20)  |      |           | (cid:21) |     |
f ∈(cid:88)F (n)
Computations given in Appendix B show that in both triangular and quadrangular meshes:
1
|     |     |     | ∇   | ϕ (x |     | n   | sign(n | (x  | x   |      | .   |     |
| --- | --- | --- | --- | ---- | --- | --- | ------ | --- | --- | ---- | --- | --- |
|     |     |     |     | ⊥ n  | f ) | f = |        | ⊥σ  | n   | f )) |     |     |
|     |     |     |     |      | ·   |     |        | ·   | −   | f    |     |     |
| |
| yielding | the result. |     |     |     |     |     |     |     |     |     |     |     |
| -------- | ----------- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
Formally
|     |     |     |     | ϕ ,(∇     |     | u         |     | (∇                 | udx,  |     |     |     |
| --- | --- | --- | --- | --------- | --- | --------- | --- | ------------------ | ----- | --- | --- | --- |
|     |     |     |     | n         | ⊥   | ) ∗ L2(Ω) |     |                    | ⊥ ) ∗ |     |     |     |
|     |     |     |     | (cid:104) |     | (cid:105) | ≈   |                    |       |     |     |     |
|     |     |     |     |           |     |           |     | (cid:90) (cid:63)n |       |     |     |     |
where (cid:63)n is a dual volume associated with node n (see Figure 5.4). Here it could be defined as
ϕ
the measure of the support of n . By analogy with the integration of the divergence we can
(cid:94)
| extrapolate | the | metric | we need | to  | put | to define | (∇  | ) u: | indeed |     |     |     |
| ----------- | --- | ------ | ------- | --- | --- | --------- | --- | ---- | ------ | --- | --- | --- |
⊥ ∗
1
|     |     |     |     | (divu) |     | =   |     | div(u)dx, |     |     |     |     |
| --- | --- | --- | --- | ------ | --- | --- | --- | --------- | --- | --- | --- | --- |
K
|     |     |     |     |     |     |     | K K          |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | ------------ | --- | --- | --- | --- | --- |
|     |     |     |     |     |     |     | | | (cid:90) |     |     |     |     |     |
and
|     |     |     |     | (cid:94) |     |     | 1   |           |     |     |     |     |
| --- | --- | --- | --- | -------- | --- | --- | --- | --------- | --- | --- | --- | --- |
|     |     |     |     | (divu)   |     | =   |     | div(u)dx, |     |     |     |     |
K
∂K
|     |     |     |     |     |     |     | | | (cid:90) | K   |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | ------------ | --- | --- | --- | --- | --- |
(cid:94)
so we see that to define the analogue for (∇ ) u we need to normalize by the measure of the
|     |     |     |     |     |     |     | ⊥ ∗ |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
boundary of the element on which it’s integrated, in this case, the boundary of (cid:63)n.

| CHAPTER   | 5. DEVELOPMENT |     |     |     | OF A | CLASS  | OF  | LONG | TIME CONSISTENT |     |     |
| --------- | -------------- | --- | --- | --- | ---- | ------ | --- | ---- | --------------- | --- | --- |
| STAGGERED | SCHEMES        |     | ON  | THE | WAVE | SYSTEM |     |      |                 |     | 145 |
(cid:94)
Definition 5.6.1. We define (∇ ) u thanks to the previous lemma and remark as:
|     |     |     |          | ⊥   | ∗     |              |           |       |           |     |     |
| --- | --- | --- | -------- | --- | ----- | ------------ | --------- | ----- | --------- | --- | --- |
|     |     |     | (cid:94) |     |       | 1            |           |       |           |     |     |
|     |     |     | ((∇      | )   | u) := |              | ϕ         | ,(∇ ) | u         |     |     |
|     |     |     |          | ⊥ ∗ | n     | ∂((cid:63)n) | n         | ⊥     | ∗ L2(Ω)   |     |     |
|     |     |     |          |     |       |              | (cid:104) |       | (cid:105) |     |     |
|     |     |     |          |     |       | |            | |         |       |           |     |     |
1
|     |     |     |     |                |            |     | ε (n) | (cid:63)f | u .      |     |     |
| --- | --- | --- | --- | -------------- | ---------- | --- | ----- | --------- | -------- | --- | --- |
|     |     |     |     |                |            |     | f     |           | f        |     |     |
|     |     |     |     | ≈ ∂((cid:63)n) |            |     |       | | |       |          |     |     |
|     |     |     |     | |              | |(cid:18)f | (n) |       |           | (cid:19) |     |     |
∈(cid:88)F
This enables to define a diffusion operator inspired by the Hodge-Laplacian:
Definition 5.6.2 (Hodge-Laplacian numerical diffusion). Let ., . , the mass-lumped scalar
h
|     |     |     |     |     |     |     |     |     | (cid:104) (cid:105) |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | ------------------- | --- | --- |
product Definition 5.2.2. We define a Hodge-Laplacian numerical diffusion on the
|                |     |       |              |     |     | R2  |         |     | RT1(Ω) |     |     |
| -------------- | --- | ----- | ------------ | --- | --- | --- | ------- | --- | ------ | --- | --- |
| Raviart-Thomas |     | space | for a domain |     | Ω   | as, | for any | Ψ,u |        |     |     |
|                |     |       |              |     | ⊂   |     |         |     | ∈      |     |     |
(cid:94)
|     |     | ∆u,Ψ      |           | :=  | ( div)   | div(u)+∇  |     | (∇  | ) (u),Ψ | ,         |     |
| --- | --- | --------- | --------- | --- | -------- | --------- | --- | --- | ------- | --------- | --- |
|     |     |           |           | h   |          | ∗         |     | ⊥   | ⊥ ∗     |           |     |
|     |     | (cid:104) | (cid:105) |     | −        |           |     |     |         |           |     |
|     |     |           |           |     | (cid:28) |           |     |     |         | (cid:29)h |     |
|     |     | (cid:103) |           |     |          | (cid:103) |     |     |         |           |     |
(cid:94)
with div(u) given by Definition 5.2.5 and (∇ ⊥ ) (u) given by Definition 5.6.1. Let, for q defined
∗
at the nodes
| (cid:103) |     |     |     |     | [[q]] |             | ε (n)q | ,   |     |     |         |
| --------- | --- | --- | --- | --- | ----- | ----------- | ------ | --- | --- | --- | ------- |
|           |     |     |     |     | ⊥σ    | :=          | σ      | n   |     |     | (5.116) |
|           |     |     |     |     |       | n (cid:88)⊂ | ∂σ     |     |     |     |         |
RT1(Ω)
then the Hodge Laplacian diffusion reads in the fully discrete setting, for u and Ψ
σ
∈
| a Raviart-Thomas |     | basis | function |     |     |     |     |     |     |     |     |
| ---------------- | --- | ----- | -------- | --- | --- | --- | --- | --- | --- | --- | --- |
(cid:94)
(cid:94)
|     |     |     | ∆u,Ψ      | σ h       | = σ [[div(u)]]+ |     |     | (cid:63)σ [[(∇ | ⊥ ) (u)]] | ⊥ . |     |
| --- | --- | --- | --------- | --------- | --------------- | --- | --- | -------------- | --------- | --- | --- |
|     |     |     | (cid:104) | (cid:105) | | |             |     | |   | |              | ∗         |     |     |
Remark 5.6.3. This s(cid:103)ection is actually the only part of the numerical scheme that is
dimension-dependent, since the N´ed´elec-Raviart-Thomas complex is different in three space di-
mensions: first the operator preceding the divergence in the sequence is not the rotated gradient
anymore, secondly the space preceding the Raviart-Thomas space is not the space of polynomials
of degree one but is now the N´ed´elec finite element space for H(rot;Ω) fields. Indeed the three
| dimensions | N´ed´elec-Raviart-Thomas |       |     |     | reads | as follows: |        |     |            |     |     |
| ---------- | ------------------------ | ----- | --- | --- | ----- | ----------- | ------ | --- | ---------- | --- | --- |
|            |                          | G1(Ω) |     | ∇   | N0(Ω) | rot         | RT1(Ω) |     | div dG0(Ω) |     |     |
−
|     |     |     |     | −→  |     | −→  |     |     | −→  |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
SothatinthreespacedimensionstheHodge-LaplacianoftheRaviart-Thomasspaceisdefined
| as a discrete | equivalent |     | of  |     |     |       |          |     |     |     |     |
| ------------- | ---------- | --- | --- | --- | --- | ----- | -------- | --- | --- | --- | --- |
|               |            |     |     | ∆u  | =   | ∇divu | rotrotu, |     |     |     |     |
−
instead of
|     |     |     |     | ∆u  | = ∇divu+∇ |     |     | curlu, |     |     |     |
| --- | --- | --- | --- | --- | --------- | --- | --- | ------ | --- | --- | --- |
⊥
in two dimensions. Conveniently, since we are interested by acoustic stabilization, only the
grad-div diffusion is relevant here and is thereby ’indifferent’ to the dimension of the domain.

| CHAPTER   | 5. DEVELOPMENT | OF A        | CLASS OF LONG | TIME CONSISTENT |     |
| --------- | -------------- | ----------- | ------------- | --------------- | --- |
| STAGGERED | SCHEMES        | ON THE WAVE | SYSTEM        |                 | 146 |
In the following numerical example, the numerical diffusion is now defined as in Defini-
tion 5.6.2.
Qualitative degradation of the long time limit: the velocity isocontours
The iso-contours of the norm of numerical long time velocity obtained with the explicit scheme
enrichedwiththeHodge-LaplaciandiffusionDefinition5.6.2areshowninFigure5.5: thisfigure
now illustrates blatantly the necessity of a stationarity perserving diffusion. By contrast with
Figure 5.2, the isocontours are in this case very far from the ones of the exact long time limit
(5.114).
(a) Exact long time limit (b) Hodge-Laplacian diffusion solution
Figure 5.5: Comparison of the norm isocontours of the exact limit and the numerical one
| obtained | with Hodge-Laplacian | diffusion |     |     |     |
| -------- | -------------------- | --------- | --- | --- | --- |
Quantitative degradation of the long time limit: the mesh convergence
Usingagainmeshesofsize5 8, 5 16, 10 32, 20 64, Figure5.6demonstratesthebreakdown
|     |     | × × | × × |     |     |
| --- | --- | --- | --- | --- | --- |
in mesh convergence of the explicit scheme stabilized with the Hodge-Laplacian numerical dif-
fusion, in comparison with the explicit scheme with the suited grad-div numerical diffusion,
indicating that the full Hodge-Laplacian diffusion is inadapted for problems requiring the pre-
| servation | of the wave kernel. |     |     |     |     |
| --------- | ------------------- | --- | --- | --- | --- |

| CHAPTER   |     | 5. DEVELOPMENT |     |     | OF  | A CLASS     | OF  | LONG TIME | CONSISTENT |     |     |
| --------- | --- | -------------- | --- | --- | --- | ----------- | --- | --------- | ---------- | --- | --- |
| STAGGERED |     | SCHEMES        |     | ON  | THE | WAVE SYSTEM |     |           |            |     | 147 |
100
| 2   | 100 |     |     |     |     |     | 2   |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
|     | k   |     |     |     |     |     | k   |     |     |     |     |
| xe) |     |     |     |     |     |     | xe) |     |     |     |     |
x
| u(  |     |     |     |     |     |     | u( y |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | ---- | --- | --- | --- | --- |
|     | −   |     |     |     |     |     | −    |     |     |     |     |
10−1
| h    |      | slope=0.9 |                          |     |     |     | h    |           |                          |     |     |
| ---- | ---- | --------- | ------------------------ | --- | --- | --- | ---- | --------- | ------------------------ | --- | --- |
| )    |      |           |                          |     |     |     | )    | slope=0.8 |                          |     |     |
| u( x | 10−1 |           |                          |     |     |     | u( y |           |                          |     |     |
|      | k    |           |                          |     |     |     | k    |           |                          |     |     |
|      |      |           |                          |     |     | RT1 |      |           |                          |     | RT1 |
|      |      |           | ExplicitHodge-Laplacian, |     |     |     |      |           | ExplicitHodge-Laplacian, |     |     |
|      |      |           | Explicit,RT1             |     |     |     |      |           | Explicit,RT1             |     |     |
10−2
|     |     | 10−3      |        |             | 10−2 |         |        | 10−3       |                    | 10−2 |     |
| --- | --- | --------- | ------ | ----------- | ---- | ------- | ------ | ---------- | ------------------ | ---- | --- |
|     |     |           |        | h           |      |         |        |            | h                  |      |     |
|     |     | (a) Error | on the | x component |      |         |        | (b) Error  | on the y component |      |     |
|     |     |           |        | −           |      |         |        |            | −                  |      |     |
|     |     | Figure    | 5.6:   | Convergence |      | in mesh | to the | exact long | time solution      |      |     |
5.7 Conclusion
In this chapter we constructed a class of staggered numerical schemes that are long time con-
sistent on the wave system for both quadrilateral and triangular meshes 2D. Our study is based
on approximation spaces that can be fitted in a de Rham complex. Using a mass-lumping on
the Raviart-Thomas mass matrix we were able to define a staggered finite volume scheme for
| the | first | order wave | system | in  | semi-discrete | setting, | then: |     |     |     |     |
| --- | ----- | ---------- | ------ | --- | ------------- | -------- | ----- | --- | --- | --- | --- |
• In explicit time integration we had to introduce a special operator of the form ∇(divu)
in order to get both stability results and preservation of the long time behaviour. This
operator exhibits convenient properties such as: natural discrete integration by parts and
an inverse Poincar´e inequality that enables to get a natural bound of the divergence in
the stability analysis. A convenient quality of this operator is that it does not actu-
ally require a dual mesh, as usual classical staggered schemes do, while still being truly
multidimensional [114, 118](as opposed to classical collocated solvers).
• ImEx time integration was also introduced and its stability was proven using arguments
|     | from | the explicit |     | time integration |     | proof. |     |     |     |     |     |
| --- | ---- | ------------ | --- | ---------------- | --- | ------ | --- | --- | --- | --- | --- |
•
In implicit time integration our methodology enables us to get another point of view
in the understanding of the efficiency near the incompressible regime of more classical
|     | staggered | schemes. |     |     |     |     |     |     |     |     |     |
| --- | --------- | -------- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
Since the choice of the velocity staggering enabled us to prove the existence of a long time
limit thanks to a peculiar HHD, then, for all of the schemes we were able to show theoretically
the convergence in long time limit to the relevant stationary state and we could reenact this
behaviour numerically with the scattering of wave in a 2D cylinder. This test case is very
convenient because it enables us to experiment with non-structured meshes on a non-trivial

| CHAPTER   | 5. DEVELOPMENT | OF A        | CLASS OF LONG | TIME CONSISTENT |     |
| --------- | -------------- | ----------- | ------------- | --------------- | --- |
| STAGGERED | SCHEMES        | ON THE WAVE | SYSTEM        |                 | 148 |
topology. Moreover the exact long time limit is known in this particular case, so we recover
qualitative results by looking at iso-contours of the velocity norm but also quantitative results
| with convergence | estimates. |     |     |     |     |
| ---------------- | ---------- | --- | --- | --- | --- |
Even though the equivalence between the low Mach limit and the long time behaviour
on the wave system is formal, we conjecture through this analysis that there is not much
room as for the choice of stabilization operator in order to get a great low Mach number
behaviour: indeed, using the de Rham formalism we defined the full Hodge-Laplacian on the
Raviart-Thomas space; this operator is composed of the grad-div operator introduced for
explicit stability combined with another part that does not necessarily preserve divergence
free fields. Through the same numerical test case we were able to show that adding this other
part in order to complete the full Hodge-Laplacian will in fact breakdown the convergence
to the relevant long time limit. Showing, as a consequence, the importance of a stationarity
preserving diffusion operator. We will pursue this work with an extension of the scheme on
| Euler’s barotropic | system. |     |     |     |     |
| ------------------ | ------- | --- | --- | --- | --- |

| Chapter    | 6      |     |                 |     |     |     |
| ---------- | ------ | --- | --------------- | --- | --- | --- |
| Extension  | of     | the | Raviart-Thomas  |     |     |     |
| staggered  | scheme |     | to compressible |     |     |     |
| barotropic | flows  |     |                 |     |     |     |
6.1 Introduction
Contents
6.1 Introduction . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 150
6.2 The Raviart-Thomas staggered scheme for the two-dimensional
Euler barotropic equations . . . . . . . . . . . . . . . . . . . . . . . . 152
6.2.1 Deriving the ’centred’ scheme . . . . . . . . . . . . . . . . . . . . . . . 152
6.2.2 Adding the appropriate diffusion . . . . . . . . . . . . . . . . . . . . . 155
6.2.3 Numerical treatment of the boundary conditions . . . . . . . . . . . . 158
6.3 Conservation . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 160
6.4 Low Mach number analysis . . . . . . . . . . . . . . . . . . . . . . . 162
6.5 Discussion on some preexisting staggered schemes . . . . . . . . . 165
6.5.1 Low Mach number behaviour . . . . . . . . . . . . . . . . . . . . . . . 165
6.5.2 Discrete entropy dissipation . . . . . . . . . . . . . . . . . . . . . . . 166
6.5.3 Other computational questions . . . . . . . . . . . . . . . . . . . . . . 167
6.6 Numerical results . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 168
6.6.1 1D Riemann problems . . . . . . . . . . . . . . . . . . . . . . . . . . . 168
6.6.2 Cylinder scattering . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 169
6.6.3 PropagationofalowMachnumberacousticwavethroughastationary
|     | low Mach number | vortex | . . . . . . . | . . . . . . . | . . . . . . . | . . . . 170 |
| --- | --------------- | ------ | ------------- | ------------- | ------------- | ----------- |
6.7 Conclusion . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 172
150

CHAPTER 6. EXTENSION OF THE RAVIART-THOMAS STAGGERED SCHEME TO
| COMPRESSIBLE |     |     | BAROTROPIC |     |     | FLOWS |     |     |     |     |     | 151 |
| ------------ | --- | --- | ---------- | --- | --- | ----- | --- | --- | --- | --- | --- | --- |
In this chapter, we extend to the Euler barotropic equations the de Rham staggered dis-
cretization introduced previously. The equations are defined, in space on a set Ω of Rd and in
| time | on  | the interval |     | [0,T] with | finite | time | T:        |     |     |     |     |     |
| ---- | --- | ------------ | --- | ---------- | ------ | ---- | --------- | --- | --- | --- | --- | --- |
|      |     |              |     |            |        | ∂    | ρ+div(ρu) |     | = 0 |     |     |     |
t

|     |     |     |     |     |    | ∂ (ρu)+div(ρu |     | u)+∇p |     | = 0 |     | (6.1) |
| --- | --- | --- | --- | --- | --- | ------------- | --- | ----- | --- | --- | --- | ----- |
|     |     |     |     |     |   | t             |     |       |     |     |     |       |
|     |     |     |     |     |    |               |     | ⊗     |     |     |     |       |
 
ργ;γ
|     |     |     |     |     |     |     | p(ρ) = |     | 1   |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | ------ | --- | --- | --- | --- | --- |
|     |     |     |     |     |    |     |        | ≥   |     |     |     |     |

 

with ρ is the density, the equation of state p = f(ρ) links the pressure p to the density and u
is the the velocity. Classically, the Euler barotropic equations can be put under the form of a
| non-linear |     | system | of  | conservation |     | laws : |       |     |      |     |       |       |
| ---------- | --- | ------ | --- | ------------ | --- | ------ | ----- | --- | ---- | --- | ----- | ----- |
|            |     |        | ∂   | U+divF(U)    |     | = 0,   | U : Ω | Rd  | [0,+ | [   | Rd+1, | (6.2) |
t
|     |     |     |     |     |     |     |     | ⊂   | ×   | ∞   | −→  |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
with
ρut
|                            |     |     |     |     |          |                                                      |           |     | R(d+1)   |     | d.  |       |
| -------------------------- | --- | --- | --- | --- | -------- | ---------------------------------------------------- | --------- | --- | -------- | --- | --- | ----- |
|                            |     |     |     |     | F(U)     | =                                                    |           |     |          |     |     | (6.3) |
|                            |     |     |     |     |          | ρu                                                   | u+pI      |     |          | ×   |     |       |
|                            |     |     |     |     |          |                                                      |           | d   | ∈        |     |     |       |
|                            |     |     |     |     |          | (cid:18)                                             | ⊗         |     | (cid:19) |     |     |       |
| Undertheassumptionthatp(ρ) |     |     |     |     |          | >                                                    |           |     |          |     |     |       |
|                            |     |     |     |     | (cid:48) | 0thesystemishyperbolicwitheigenvaluesinthenormalized |           |     |          |     |     |       |
| direction                  |     | n   |     |     |          |                                                      |           |     |          |     |     |       |
|                            |     |     |     | λ   | := u     | n, λ                                                 | := u n+c, |     | λ :=     | u n | c,  | (6.4) |
|                            |     |     |     | 1   |          | 2                                                    |           |     | 3        |     |     |       |
|                            |     |     |     |     |          | ·                                                    | ·         |     |          | ·   | −   |       |
where c := p(ρ). With this in mind, we recall the challenges this chapter is tackling :
(cid:48)
1 The s(cid:112)taggered discretization implies inherently that a vectorial unknown (typically the
velocity field) is approximated at the faces while the discrete density and pressure are
located at the center of the cells. Defining properly product of variables located at the
|     | cells | with | variables | located |     | at the | faces is | then | necessary. |     |     |     |
| --- | ----- | ---- | --------- | ------- | --- | ------ | -------- | ---- | ---------- | --- | --- | --- |
2 A consequence of 1 is that it becomes also non-trivial to understand conservation for
|     | terms | such | as  |     |     |     |        |     |     |     |     |     |
| --- | ----- | ---- | --- | --- | --- | --- | ------ | --- | --- | --- | --- | --- |
|     |       |      |     |     |     |     | div(ρu |     | u). |     |     |     |
⊗
3 Finally, the staggered discretization was first motivated by the approximation of the low
Mach number regime, so that using, in part, the underlying Hodge-Helmholtz decompos-
ition of the de Rham scheme, we want to conclude on the low Mach number behaviour of
|     | this      | scheme | in    | some sense |      | we will | define.      |     |     |        |     |     |
| --- | --------- | ------ | ----- | ---------- | ---- | ------- | ------------ | --- | --- | ------ | --- | --- |
|     | To answer |        | these | questions, | this | chapter | is conducted |     | as  | follow | :   |     |
1) In section 6.2 is presented our proposition for the treatment of 1 and from there how
is derived the numerical scheme on Euler equations. The choice of sensible numerical
diffusion operators and the treatment of boundary conditions is also developed.

CHAPTER 6. EXTENSION OF THE RAVIART-THOMAS STAGGERED SCHEME TO
| COMPRESSIBLE |     |     | BAROTROPIC |     |     | FLOWS |     |     |     |     | 152 |
| ------------ | --- | --- | ---------- | --- | --- | ----- | --- | --- | --- | --- | --- |
2) Then, in section 6.3, the problem of defining conservation in the framework of the de
|     | Rham | staggered |     | discretization |     | 2 is treated | .   |     |     |     |     |
| --- | ---- | --------- | --- | -------------- | --- | ------------ | --- | --- | --- | --- | --- |
3) In section 6.4, the low Mach number of the Raviart-Thomas staggered scheme 3 is
|     | addressed |     | using elements |     | developed | throughout |     | this thesis. |     |     |     |
| --- | --------- | --- | -------------- | --- | --------- | ---------- | --- | ------------ | --- | --- | --- |
4) Next, as in chapter 3 and chapter 5, the Raviart-Thomas staggered scheme is briefly put
in perspective with respect to other staggered discretizations in section 6.5.
5) In section 6.6, we put to the test the conservation and low Mach number properties of
|     | the | de Rham | scheme | with | numerical | simulations. |     |     |     |     |     |
| --- | --- | ------- | ------ | ---- | --------- | ------------ | --- | --- | --- | --- | --- |
6) Finally, the findings and conclusions of this chapter are summed up in section 6.7.
| 6.2 | The    | Raviart-Thomas |       |            |     | staggered |     | scheme | for the | two-dimen- |     |
| --- | ------ | -------------- | ----- | ---------- | --- | --------- | --- | ------ | ------- | ---------- | --- |
|     | sional |                | Euler | barotropic |     | equations |     |        |         |            |     |
We recall briefly here the approximation spaces used and the different definitions needed in
the following. Details on the discretization spaces can be found in chapter 4. We base our
approximation of the discrete unknowns on the N´ed´elec-Raviart-Thomas complex, which reads
| ([84, | 99]) | on quadrangular |     | and | triangular | meshes |     |         |     |     |     |
| ----- | ---- | --------------- | --- | --- | ---------- | ------ | --- | ------- | --- | --- | --- |
|       |      |                 |     |     |            | ∇⊥     |     | div     |     |     |     |
|       |      |                 |     |     | cG1(Ω)     | RT1(Ω) | −   | dG0(Ω). |     |     |     |
|       |      |                 |     |     |            | −→     | −→  |         |     |     |     |
Thus, the discrete momentum q := ρu is in the Raviart-Thomas space (RT1(Ω)) and the
pressure and density are cellwise constant (dG0(Ω)). In the following, we will denote for any
face f
|     | ∈   | F   |     |     |      |             |     |     |     |     |     |
| --- | --- | --- | --- | --- | ---- | ----------- | --- | --- | --- | --- | --- |
|     |     |     |     |     |      | v(x+tn)+v(x |     | tn) |     |     |     |
|     |     |     |     |     | v := | lim         |     | −   | ,   |     |     |
|     |     |     |     | {{  | }}f  | 0+          |     | 2   |     |     |     |
t
→
v.
| for   | any field | or  | vector | field     |     |        |     |     |     |     |     |
| ----- | --------- | --- | ------ | --------- | --- | ------ | --- | --- | --- | --- | --- |
| 6.2.1 | Deriving  |     | the    | ’centred’ |     | scheme |     |     |     |     |     |
We present here the derivation of the ’centred’ scheme in a discontinuous Galerkin inspired
manner, with the momentum approximated in the Raviart-Thomas space and the density
being constant by cell. Contrary to the derivation on the wave system the presentation of the
boundary conditions is postponed to subsection 6.2.3: the domain is considered periodic until
| subsection |     | 6.2.3. |     |     |     |     |     |     |     |     |     |
| ---------- | --- | ------ | --- | --- | --- | --- | --- | --- | --- | --- | --- |

CHAPTER 6. EXTENSION OF THE RAVIART-THOMAS STAGGERED SCHEME TO
| COMPRESSIBLE |     | BAROTROPIC |     | FLOWS |     |     |     |     |     |     | 153 |
| ------------ | --- | ---------- | --- | ----- | --- | --- | --- | --- | --- | --- | --- |
Mass equation
As for the approximation of the wave system, we multiply the mass equation by a test function
ϕ dG0(Ω) where we first suppose that q (Ω)d and ρ (Ω) and integrate on the
|     |     |     |     |     |     | ∞   |     |     | ∞   |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| ∈   |     |     |     |     | ∈   | C   |     |     | ∈ C |     |     |
domain Ω:
|     |     |     | (∂        | t ρ)ϕdx+ |     | div(q)ϕdx |     | =   | 0.  |     | (6.5) |
| --- | --- | --- | --------- | -------- | --- | --------- | --- | --- | --- | --- | ----- |
|     |     |     | (cid:90)Ω |          |     | K         |     |     |     |     |       |
|     |     |     |           |          | K   | (cid:90)  |     |     |     |     |       |
(cid:88)∈C
So, integrating by parts, as for the derivation of the scheme on the wave system;
|     |            |          | div(q)ϕdx | =          |          | q        | ∇ϕdx+ |          | q nϕdΓ | .        | (6.6) |
| --- | ---------- | -------- | --------- | ---------- | -------- | -------- | ----- | -------- | ------ | -------- | ----- |
|     |            | K        |           |            | −        | K        | ·     |          | ∂K ·   |          |       |
|     | K          | (cid:90) |           | K          | (cid:20) | (cid:90) |       | (cid:90) |        | (cid:21) |       |
|     | (cid:88)∈C |          |           | (cid:88)∈C |          |          | =0    |          |        |          |       |
(cid:124)(cid:123)(cid:122)(cid:125)
| Hence, plugging | (6.6) | in  | (6.5) we | obtain |     |     |        |     |     |     |     |
| --------------- | ----- | --- | -------- | ------ | --- | --- | ------ | --- | --- | --- | --- |
|                 |       |     | (∂       | ρ)ϕdx+ |     |     | q nϕdΓ |     | 0.  |     |     |
|                 |       |     |          | t      |     |     |        | =   |     |     |     |
·
|     |     |     | (cid:90)Ω |     |              | ∂K       |     |     |     |     |     |
| --- | --- | --- | --------- | --- | ------------ | -------- | --- | --- | --- | --- | --- |
|     |     |     |           |     | K (cid:88)∈C | (cid:90) |     |     |     |     |     |
RT1(Ω)
In the case where q and ρ are not regular but respectively and dQ0(Ω), the normal
trace q n is still defined and in fact, because of the properties of the Raviart-Thomas, it is
·
|     |     |     |     |     | 1   |     |     |     | RT1(Ω) |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | ------ | --- | --- |
constant by face. So then, by taking ϕ = K for K a cell, for q we obtain
∈
|     |     |     |     | K ∂ t ρ | +    | σ q | n     | = 0, |     |     |     |
| --- | --- | --- | --- | ------- | ---- | --- | ----- | ---- | --- | --- | --- |
|     |     |     |     | | | K   |      | | | | · K,σ |      |     |     |     |
|     |     |     |     |         | σ ∂K |     |       |      |     |     |     |
(cid:88)∈
with q n the normal component of the Raviart-Thomas q at the face σ.
| · K,σ    |          |     |     |     |     |     |     |     |     |     |     |
| -------- | -------- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Momentum | equation |     |     |     |     |     |     |     |     |     |     |
Analogously, we multiply the momentum equation, where q,p,ρ are supposed for now
RT1(Ω)
(Ω) regular, by a test function Ψ and integrate on the domain
| C ∞ − |     |           |          | ∈   |            |          |     |          |          |     |       |
| ----- | --- | --------- | -------- | --- | ---------- | -------- | --- | -------- | -------- | --- | ----- |
|       |     |           |          |     |            | q        | q   |          |          |     |       |
|       |     |           | ∂ q Ψdx+ |     | div        | ⊗        | +I  | p        | Ψdx = 0. |     | (6.7) |
|       |     |           | t        |     |            |          |     | d        |          |     |       |
|       |     |           | ·        |     |            |          | ρ   | ·        |          |     |       |
|       |     | (cid:90)Ω |          | K   | (cid:90) K | (cid:20) |     | (cid:21) |          |     |       |
(cid:88)∈C

CHAPTER 6. EXTENSION OF THE RAVIART-THOMAS STAGGERED SCHEME TO
COMPRESSIBLE BAROTROPIC FLOWS 154
q q
Let G := ⊗ +I p, integrating by parts yields
d
ρ
divG Ψdx= G : ∇Ψdx+ (G Ψ) n dΓ
· − · ·
K (cid:34) K ∂K (cid:35)
K (cid:90) K (cid:90) (cid:90)
(cid:88)∈C (cid:88)∈C
= G : ∇Ψdx+ [[(G n) Ψ]]dΓ (6.8)
− · ·
K f
K (cid:90) f int(cid:90)
(cid:88)∈C ∈(cid:88)F
+ (G n) ΨdΓ,
· ·
f
f b(cid:90)
(cid:88)∈F
where we used the symmetry of G for (G Ψ) n = (G n) Ψ. Now by regularity assumptions
· · · ·
on the unknowns
[[(G n) Ψ]]dΓ = (G n) [[Ψ]]dΓ. (6.9)
· · · ·
f f
f int(cid:90) f int(cid:90)
∈(cid:88)F ∈(cid:88)F
So that gathering, (6.9) and (6.8) in (6.7) yields
∂ q Ψdx G : ∇Ψdx+ (G n) [[Ψ]]dΓ+ (G n) ΨdΓ = 0.
t
· − · · · ·
(cid:90)Ω K (cid:90) K f int(cid:90) f f b(cid:90) f
(cid:88)∈C ∈(cid:88)F (cid:88)∈F
(6.10)
In the case where the unknowns are such that q RT1(Ω), ρ,p dG0(Ω), the trace is replaced
∈ (cid:92) ∈
by a constant by face numerical trace (G n) (G n) defined as
f f
· | −→ · |
(cid:92) q
f int (G n) := q n +n p . (6.11)
∀ ∈ F · | f · (cid:26)(cid:26) ρ (cid:27)(cid:27)f {{ }}f
The boundary faces will be treated in the following section. Then by taking Ψ = Ψ the
σ
Raviart-Thomas basis function associated with a face σ we have, using the mass-lumping
∈ F
of Definition 5.2.2 that (6.10) becomes :
(cid:91) (cid:91)
D ∂ q G : ∇Ψ dx+ G n [[Ψ ]]dΓ+ G n Ψ dΓ = 0,
σ t σ σ σ σ
| | − · · · ·
K f f
K (σ)(cid:90) f int(cid:90) f b(cid:90)
(cid:88)∈C ∈(cid:88)F (cid:88)∈F
which is in fact equivalent to the more concise formulation:
D ∂ q + q ⊗ q +I p : ∇Ψ dx+ G (cid:91) n Ψ dΓ = 0. (6.12)
| σ | t σ − ρ d σ · · σ
K (cid:20) (cid:90) K(cid:18) (cid:19) (cid:90) ∂K (cid:21)
(cid:88)∈C

CHAPTER 6. EXTENSION OF THE RAVIART-THOMAS STAGGERED SCHEME TO
| COMPRESSIBLE |     |        | BAROTROPIC |             |     | FLOWS     |     |     |     |     |     |     | 155 |
| ------------ | --- | ------ | ---------- | ----------- | --- | --------- | --- | --- | --- | --- | --- | --- | --- |
| 6.2.2        |     | Adding | the        | appropriate |     | diffusion |     |     |     |     |     |     |     |
Mass equation
We inoculate diffusive terms through a classical upwinding of the numerical flux
|     |     |     |     |     | K   | ∂ ρ | +   | σ q(cid:92)n |     | = 0, |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | ------------ | --- | ---- | --- | --- | --- |
|     |     |     |     |     |     | t K |     |              | K,σ |      |     |     |     |
|     |     |     |     |     | |   | |   |     | | |          | ·   |      |     |     |     |
σ ∂K
(cid:88)∈
| where | now | the    | numerical | trace      | is  | modified | to   |       |       |        |      |     |     |
| ----- | --- | ------ | --------- | ---------- | --- | -------- | ---- | ----- | ----- | ------ | ---- | --- | --- |
|       |     |        |           |            |     |          |      | u     | n K,σ | +c max |      |     |     |
|       |     |        | σ         | q(cid:92)n |     | = q      | n    | + | · | |     |        | (ρ ρ | ) , |     |
|       |     |        |           |            | K,σ |          | K,σ  |       |       |        | Lσ   | K   |     |
|       |     |        | ∀ ∈       | F          | ·   |          | ·    |       | 2     |        | −    |     |     |
|       |     |        |           |            |     |          | q    | n K,σ |       |        |      |     |     |
| with  | c   | := max | p(ρ       | ) and      | u   | n        | := · |       | where | ρ =    | ρ .  |     |     |
|       | max |        | (cid:48)  | K          |     | K,σ      |      |       |       | σ      | }}σ  |     |     |
|       |     | K      |           |            |     | ·        |      | ρ σ   |       |        | {{   |     |     |
∈C
(cid:112)
| Remark |     | 6.2.1. | This | is equivalent |     | to the | finite | volume | formulation: |     |     |     |     |
| ------ | --- | ------ | ---- | ------------- | --- | ------ | ------ | ------ | ------------ | --- | --- | --- | --- |
u +c
|     |     |     |       |       |           |            |     |           | σ   | max |                |     |        |
| --- | --- | --- | ----- | ----- | --------- | ---------- | --- | --------- | --- | --- | -------------- | --- | ------ |
|     |     |     | K ∂ t | ρ K + |           | σ ε K (σ)q | σ = |           | | | |     | σ ε K (σ)[[ρ]] | σ , | (6.13) |
|     |     |     | | |   |       | |         | |          |     |           |     | 2   | | |            |     |        |
|     |     |     |       |       | σ ∂K      |            |     | σ ∂K      |     |     |                |     |        |
|     |     |     |       |       | (cid:88)∈ |            |     | (cid:88)∈ |     |     |                |     |        |
q
σ
with q σ the Raviart-Thomas d.o.f at the face σ and u σ = . This formulation will turn
ρ
σ
out to be useful for: the low Mach number analysis section 6.4 and the following proposition
| Proposition |     | 6.2.1. |     |     |     |     |     |     |     |     |     |     |     |
| ----------- | --- | ------ | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
Integrating with Euler explicit method (6.13), we note instantly the following proposition
| on  | the positivity |     | of the | density: |     |     |     |     |     |     |     |     |     |
| --- | -------------- | --- | ------ | -------- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
Proposition 6.2.1 (Positivity of the discrete density). Let (ρn) the sequence of constant
n N
⊂
| by  | cell density |     | built with | the | following | relation: |     | for | all cell | K   |     |     |     |
| --- | ------------ | --- | ---------- | --- | --------- | --------- | --- | --- | -------- | --- | --- | --- | --- |
|     |              |     |            |     |           |           |     |     |          | ∈   | C   |     |     |
K
|     |     | | |(ρn+1 |     | ρn )+ |           | σ ε (σ)un |     | ρn  | =         | dn  | σ ε (σ)[[ρn]] | ,   | (6.14) |
| --- | --- | -------- | --- | ----- | --------- | --------- | --- | --- | --------- | --- | ------------- | --- | ------ |
|     |     |          | K   | K     |           | K         | σ   | }}σ |           | σ   | K             | σ   |        |
|     |     | δt       | −   |       |           | | |       | {{  |     |           |     | | |           |     |        |
|     |     |          |     |       | σ ∂K      |           |     |     | σ         | ∂K  |               |     |        |
|     |     |          |     |       | (cid:88)∈ |           |     |     | (cid:88)∈ |     |               |     |        |
then if
un
|     |     |     |     |     |     | σ   |     | dn  | | σ |, |     |     |     | (6.15) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | ------ | --- | --- | --- | ------ |
σ
|     |        |      |               |     |     | ∀        | ∈ F | ≥   | 2   |     |     |     |     |
| --- | ------ | ---- | ------------- | --- | --- | -------- | --- | --- | --- | --- | --- | --- | --- |
| and | δt the | time | step verifies |     |     |          |     |     |     |     |     |     |     |
|     |        |      |               |     |     | δt maxdn | max | ∂K  |     |     |     |     |     |
σ
|     |     |     |     |     |     | σ   | K     | |   | |   | 1   |     |     |        |
| --- | --- | --- | --- | --- | --- | --- | ----- | --- | --- | --- | --- | --- | ------ |
|     |     |     |     |     |     | ∈F  |       | ∈C  | <   | ,   |     |     | (6.16) |
|     |     |     |     |     |     |     | min K |     |     | 2   |     |     |        |
|     |     |     |     |     |     |     | K |   | |   |     |     |     |     |        |
∈C
| the | relation | (6.14) | preserves |     | positivity | :   |      |     |      |      |     |     |     |
| --- | -------- | ------ | --------- | --- | ---------- | --- | ---- | --- | ---- | ---- | --- | --- | --- |
|     |          |        |           |     | K          |     | ρn > |     | ρn+1 | > 0. |     |     |     |
|     |          |        |           |     |            |     | K    | 0 = |      |      |     |     |     |
|     |          |        |           |     | ∀          | ∈ C |      | ⇒   | K    |      |     |     |     |

CHAPTER 6. EXTENSION OF THE RAVIART-THOMAS STAGGERED SCHEME TO
COMPRESSIBLE BAROTROPIC FLOWS 156
Proof. Equation (6.14) can be rewritten as
δt σ ε (σ)un δt σ ε (σ)un
ρn+1 = ρn 1 | | K σ +dn + | | dn K σ ρn , (6.17)
K K − K 2 σ K σ − 2 Lσ
(cid:32) (cid:33) (cid:32) (cid:33)
σ ∂K | | (cid:18) (cid:19) σ ∂K | |
(cid:88)∈ (cid:88)∈
un
Suppose now ρn > 0, then on the one hand, because of the condition (6.15) : dn | σ |
σ
≥ 2 ≥
ε (σ)un ε (σ)un
K σ, we have dn K σ 0 and as a consequence
σ
2 − 2 ≥
δt σ ε (σ)un
| | dn K σ ρn 0. (6.18)
K σ − 2 Lσ ≥
(cid:32) (cid:33)
σ ∂K | |
(cid:88)∈
On the other hand, again with (6.15)
ε (σ)un un
0 K σ +dn | σ | +dn 2dn, (6.19)
σ σ σ
≤ 2 ≤ 2 ≤
so (6.19) with the CFL (6.16) yields
δt σ ε (σ)un
1 | | K σ +dn > 0. (6.20)
− K 2 σ
σ ∂K | | (cid:18) (cid:19)
(cid:88)∈
Combining (6.18) and (6.20) yields
δt σ ε (σ)un
ρn+1 ρn 1 | | K σ +dn > 0.
K ≥ K − K 2 σ
(cid:32) (cid:33)
σ ∂K | | (cid:18) (cid:19)
(cid:88)∈
Remark 6.2.2. In our case, defining un σ := ρ q σ n n σ with ρ σ := {{ ρn }}σ and dn σ = | un σ| + 2 cn max >
un
| σ| ,(6.13) is equivalent to (6.14), so that the discrete mass equation preserves the positivity of
2
the density. Moreover, since we have a strict inequality
un
dn > | σ |
σ
2
the positivity related CFL is relaxed to
δt maxdn max ∂K
σ σ K | | 1
∈F ∈C ,
min K ≤ 2
K | |
∈C
or in this case
δt (max un +cn ) max ∂K
σ | σ | max K | |
∈F ∈C 1
min K ≤
K | |
∈C

CHAPTER 6. EXTENSION OF THE RAVIART-THOMAS STAGGERED SCHEME TO
COMPRESSIBLE BAROTROPIC FLOWS 157
Momentum equation
As we are now getting used to, the vectorial equation is much more difficult to stabilize,
especially because vectorial unknowns offer richer differential operators, thereby making their
approximation and preservation more subtle. As a consequence, we will split the numerical
diffusion in an acoustic part and a transport part: the former, which is of the order of the
sound propagation speed, will respect the acoustic kernel, and will be based on the grad-div
acoustic stationarity preserving diffusion Definition 5.2.5: for a fixed Raviart-Thomas basis
function Ψ
σ
q RT1(Ω) ( div )d (cid:93) ivq,Ψ = σ [[d (cid:93) ivq]] ,
∗ σ L2 σ
∀ ∈ (cid:104) − (cid:105) | |
where the operator div is defined as in Definition 4.4.1.
∗
On the transport scale, no special structure preservation is needed, any naive diffusive
operator will suffice; classical upwinding are often resulting on a discrete operator consistent
with a metric (related to the mesh) multiplied by a Laplacian. Thanks to the Hodge-Laplacian
we are able to define in the staggered framework a similar numerical diffusion (see Defini-
tion 5.6.2): since the transport is of the order of the fluid velocity we integrate in the discrete
momentum equation a numerical diffusion following the Hodge-Laplacian diffusion on this
scale Definition 5.6.2: for a fixed Raviart-Thomas basis function Ψ ;
σ
u u (cid:94)
q RT1(Ω) max ∆q,Ψ = max ( div) div(q)+∇ (∇ ) (q),Ψ ,
∀ ∈ 2 (cid:104) (cid:105) h 2 − ∗ ⊥ ⊥ ∗
(cid:28) (cid:29)h
(cid:103) (cid:103)
with
1 q
u := max dx ,
max
K ∈C(cid:12)| K | (cid:90) K ρ (cid:12)2
(cid:12) (cid:12)
where . is the 2 Euclidian norm. (cid:12) (cid:12)
| | 2 − (cid:12) (cid:12)
By adding these diffusion operators in (6.12), we then obtain the following discrete mo-
mentum equation
D ∂ q + q ⊗ q +I p : ∇Ψ dx+ G (cid:91) n Ψ dΓ
σ t σ 2 σ σ
| | − ρ · ·
K (cid:20) (cid:90) K(cid:18) (cid:19) (cid:90) ∂K (cid:21)
(cid:88)∈C (6.21)
= c max ( div )d (cid:93) ivq,Ψ + u max ∆q,Ψ ,
∗ σ h σ h
2 (cid:104) − (cid:105) 2 (cid:104) (cid:105)
with the numerical trace (6.11). (cid:103)
Remark 6.2.3. Denoting, the discrete convection term
q q q
C(q) = ⊗ : ∇Ψ dx+ q n Ψ dΓ ,
σ σ σ
− ρ · ρ ·
K (cid:20) (cid:90) K (cid:90) ∂K (cid:26)(cid:26) (cid:27)(cid:27) (cid:21)
(cid:88)∈C

CHAPTER 6. EXTENSION OF THE RAVIART-THOMAS STAGGERED SCHEME TO
COMPRESSIBLE BAROTROPIC FLOWS 158
(6.21) yields the equivalent formulation
D ∂ q +C(q) + σ [[p]] = u max +c max σ [[d (cid:94) iv(q)]] + u max (cid:63)σ [[( (cid:94) ∇ ) q]] , (6.22)
| σ | t σ σ | | σ 2 | | σ 2 | | ⊥ ∗ ⊥σ
with the notations used in chapter 5:
D (cid:94) 1
σ
| (cid:63)σ | := | σ |, [[h]] ⊥σ := ε σ (n)h n , (∇ ⊥ ) ∗ q := ∂((cid:63)n) ε f (n) | (cid:63)f | q f
| | n ∂σ (cid:18) (cid:19)n | | f (n)
(cid:88)⊂ ∈(cid:88)F
Formulation (6.22) will turn out to be convenient for low Mach number analysis.
6.2.3 Numerical treatment of the boundary conditions
Since the numerical scheme we proposed is Discontinuous Galerkin (DG) inspired, it is sensible
to try and replicate the treatment of the boundary conditions of these methods in our case
[119, 112, 120]. The difference in the treatment of the boundary conditions with the wave
system lies in the fact that, for the latter, boundary fluxes are directly derived using the
explicit formula of a solution of a Riemann problem. For a non-linear system such as the Euler
system, the determination of a boundary flux is much less evident. In particular, if we aim at
imposing a state U or imposing a wall boundary condition, it is not obvious how we should
b
weight the information coming from the inside of the domain Ω with respect to information at
the boundary ∂Ω. We will use two types of boundary conditions; wall boundary conditions and
inlet/oulet boundary conditions. For the former, we will use a general approximation which
relies on the following definition:
Definition 6.2.1 (Generic wall flux). Let F(U ,U ,n) a numerical flux. Let
L R
1 0t
P (n) :=
wall 0 I 2n n
d
(cid:18) − ⊗ (cid:19)
with 0 R2 with each component equal to 0 and 0t its transposed vector. We define the wall
∈
flux associated with the given numerical flux as
F (U,n) := F(U,P (n)U,n).
wall wall
We can use for example the wall Roe flux
FRoe(U,n) := FRoe(U,P (n)U,n),
wall wall
where FRoe is the classical Roe flux [12], which will turn out to be relevant for the low Mach
number behaviour. As for the inlet/outlet conditions, a reasonable approximation of the phys-
ical flux at the boundary is the Steger-Warming flux, defined as;
Definition 6.2.2 (Steger-Warming flux). The Steger-Warming flux F (U,U ,n) is
SW b
F (U,U ,n) := F(U ) n+A+(U ,n)(U U ),
SW b b b b
· −

CHAPTER 6. EXTENSION OF THE RAVIART-THOMAS STAGGERED SCHEME TO
COMPRESSIBLE BAROTROPIC FLOWS 159
with F(U ) the flux (6.3) A+(U ,n) the positive part of the Jacobian of the system,
b b
0 nt
A(U,n) := .
 
c2n u nu u n+u nI
d
− · ⊗ ·
 
TheexactformulationofthepositivepartoftheJacobianisthengivenbysimpledefinitions
and computations: for a diagonalizable matrix A the positive part A+ is defined as
A+ := Pdiag(λ+)P 1,
−
where diag(λ+) is the diagonal matrix, with the positive part of the eigenvalues of the original
matrix as diagonal entries. Then for Euler barotropic equations (6.1), this yields the following
formula
(u n c)+ 0 0
· −
 
A+(U,n) := P 0 (u n)+ 0 P 1,
−
 · 
 
 
 0 0 (u n+c)+ 
 · 
 
with, denoting t is the vector n rotated by π/2 in two dimensions,
u n+c nt
·
2c −2c
 
1 0 1
P := P − 1 :=  u t tt .
   − · 
u cn t u+cn  
−  
   u n c nt 
 · − 
 − 2c 2c 
 
 
Finally, wegatherthedefinitionsoftheboundaryfluxeswiththedefinitionofthenumerical
scheme on a periodic domain to obtain the discrete-in-space system in a bounded domain:
K ∂ ρ + σ q(cid:92)n dΓ = 0, K (6.23a)
t K K,σ
| | | | · ∀ ∈ C
 σ ∂K
(cid:88)∈

          | D σ | ∂ t q σ + K (cid:88)∈C (σ) (cid:34) (cid:90) K(cid:18) q ⊗ ρ q +I 2 p (cid:19) : ∇Ψ σ dx+ (cid:90) ∂K G (cid:91) · n · Ψ σ dΓ (cid:35) = (6.23b)
  c max ( div )d (cid:93) ivq,Ψ + u max ∆q,Ψ , σ
 ∗ σ h σ h
  2 (cid:104) − (cid:105) 2 (cid:104) (cid:105) ∀ ∈ F




where   (cid:103)

CHAPTER 6. EXTENSION OF THE RAVIART-THOMAS STAGGERED SCHEME TO
| COMPRESSIBLE |       | BAROTROPIC |        |     | FLOWS    |      |            |     |     |     |       | 160 |
| ------------ | ----- | ---------- | ------ | --- | -------- | ---- | ---------- | --- | --- | --- | ----- | --- |
| Definition   | 6.2.3 | (Mass      | flux). |     | The mass | flux | is defined |     | as  |     |       |     |
|              |       |            |        |     | u        | n    | +c         |     |     |     |       |     |
|              |       |            | q      | n   |          | K,σ  | max        | (ρ  | ρ   | if  | σ int |     |
|              |       |            |        | K,σ | + | ·    | |    |            | Lσ  | K ) |     |       |     |
|              |       |            |        | ·   |          | 2    |            |     | −   |     | ∈ F   |     |

|     | q(cid:92) n | :=  |   |     |     |      |      |     |     |      | b                |     |
| --- | ----------- | --- | --- | --- | --- | ---- | ---- | --- | --- | ---- | ---------------- | --- |
|     | K,σ         |     |    |     | [F  | (U,U | ,n)] |     |     | if σ |                  |     |
|     | ·           |     |   |     | SW  |      | b    | [ρ] |     |      | ∈ Fi nlet/outlet |     |


b
|     |     |     |     |     | [F   | (U,n)] | (=  | 0)  |     | if  | σ       |     |
| --- | --- | --- | --- | --- | ---- | ------ | --- | --- | --- | --- | ------- | --- |
|     |     |     |   |     | wall |        | [ρ] |     |     |     | ∈ Fwall |     |

 
,n)] 
where [F (U,U given by the first component of Definition 6.2.2 and
|         | SW   | b     | [ρ] |     |                 |     |     |            |        |     |     |     |
| ------- | ---- | ----- | --- | --- | --------------- | --- | --- | ---------- | ------ | --- | --- | --- |
| [F (U,U | ,n)] | given | by  | the | first component |     | of  | Definition | 6.2.1. |     |     |     |
| wall    | b    | [ρ]   |     |     |                 |     |     |            |        |     |     |     |
and
| Definition | 6.2.4 | (Momentum |     |     | flux). The | momentum |     | flux | is defined |     |     |     |
| ---------- | ----- | --------- | --- | --- | ---------- | -------- | --- | ---- | ---------- | --- | --- | --- |
q
|     |     |     |     | q   | n   |     | +n p |     | if  | σ   | int |     |
| --- | --- | --- | --- | --- | --- | --- | ---- | --- | --- | --- | --- | --- |
}}f
|     |     |     |     |     | · ρ              |                   | {{  |     |     | ∈ F |     |     |
| --- | --- | --- | --- | --- | ---------------- | ----------------- | --- | --- | --- | --- | --- | --- |
|     |     |     |     |    | (cid:26)(cid:26) | (cid:27)(cid:27)f |     |     |     |     |     |     |

(cid:91) 
|     |     | G   | n := |   |     |      |       |     |      | b                |     |     |
| --- | --- | --- | ---- | --- | --- | ---- | ----- | --- | ---- | ---------------- | --- | --- |
|     |     |     | ·    |    | [F  | (U,U | ,n)]  |     | if σ |                  |     |     |
|     |     |     |      |   | SW  |      | b [q] |     |      | ∈ Fi nlet/outlet |     |     |

b
|     |     |     |     |    | [F   | (U,n)] |     |     | if  | σ    |     |     |
| --- | --- | --- | --- | --- | ---- | ------ | --- | --- | --- | ---- | --- | --- |
|     |     |     |     |   | wall |        | [q] |     |     | ∈ Fw | all |     |

 

where [F (U,U ,n)] isgiven once again by the second and third components of Defini-
|     | SW  | b   | [q] |     |     |     |     |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
tion 6.2.2. [F (U,U ,n)] is given by the second and third components of Definition 6.2.1.
|     | wall |     | b   | [q] |     |     |     |     |     |     |     |     |
| --- | ---- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
6.3 Conservation
Let us address now the question of conservation. The proposed staggered scheme is based on
a Galerkin framework which states that the equations of interest are projected on a discrete
space. The discrete density equation is given by a classical flux formulation, yielding, by the
very definition of these fluxes, a local conservation of the discrete density. The case of the
| discrete | momentum |     | equation | is  | more demanding: |     |     |     |     |     |     |     |
| -------- | -------- | --- | -------- | --- | --------------- | --- | --- | --- | --- | --- | --- | --- |
Proposition 6.3.1 (Momentum global conservation). Suppose that the domain is periodic,
| then the | momentum         |     | equation    | (6.23b)   | yields:           |      |           |          |                   |      |     |     |
| -------- | ---------------- | --- | ----------- | --------- | ----------------- | ---- | --------- | -------- | ----------------- | ---- | --- | --- |
|          |                  |     |             |           | 1                 |      |           |          | 0                 |      |     |     |
|          |                  |     | ∂           | t q,      |                   | = 0, |           | ∂ t q,   |                   | = 0, |     |     |
|          |                  |     |             |           | 0                 |      |           |          | 1                 |      |     |     |
|          |                  |     | (cid:28)    | (cid:18)  | (cid:19)(cid:29)h |      | (cid:28)  | (cid:18) | (cid:19)(cid:29)h |      |     |     |
| which    | are the discrete |     | equivalents |           | of                |      |           |          |                   |      |     |     |
|          |                  |     |             |           | qxdx              |      |           | qydx     |                   |      |     |     |
|          |                  |     |             |           | ∂                 | = 0, |           | ∂        | =                 | 0.   |     |     |
|          |                  |     |             |           | t                 |      |           | t        |                   |      |     |     |
|          |                  |     |             | (cid:90)Ω |                   |      | (cid:90)Ω |          |                   |      |     |     |

CHAPTER 6. EXTENSION OF THE RAVIART-THOMAS STAGGERED SCHEME TO
COMPRESSIBLE BAROTROPIC FLOWS 161
Proof. (6.23b) is equivalent to
Ψ RT1(Ω) ∂ q, Ψ + q ⊗ q +I p : ∇Ψdx+ G (cid:91) n ΨdΓ
t h 2
∀ ∈ (cid:104) (cid:105) ρ · ·
K (cid:88)∈C (cid:34) (cid:90) K(cid:18) (cid:19) (cid:90) ∂K (cid:35) (6.24)
= c max ( div )d (cid:93) ivq,Ψ + u max ∆q,Ψ .
∗ h h
2 (cid:104) − (cid:105) 2 (cid:104) (cid:105)
But notice that (cid:103)
q q 1
⊗ +I p : ∇ dx = 0, (6.25)
2
ρ 0
K (cid:90) K(cid:18) (cid:19) (cid:20)(cid:18) (cid:19)(cid:21)
(cid:88)∈C
and by Definition 4.4.1
(cid:93) 1 (cid:93) 1
( div )divq, = divq, div = 0 (6.26)
∗ h
(cid:104) − 0 (cid:105) − 0
(cid:18) (cid:19) (cid:28) (cid:18) (cid:19)(cid:29)L2
and also by Definition 4.4.2
(cid:94) 1 (cid:94) 1
∇ (∇ ) (q), = (∇ ) (q),(∇ ) = 0, (6.27)
⊥ ⊥ ∗ 0 ⊥ ∗ ⊥ ∗ 0
(cid:28) (cid:18) (cid:19)(cid:29)h (cid:28) (cid:18) (cid:19)(cid:29)L2
sothat,gathering(6.26)(6.27)andbythedefinitionoftheHodgeLaplaciannumericaldiffusion
(Definition 5.6.2) we get
c max ( div )d (cid:93) ivq, 1 = 0, u max ∆q, 1 = 0. (6.28)
∗ h h
2 (cid:104) − 0 (cid:105) 2 (cid:104) 0 (cid:105)
(cid:18) (cid:19) (cid:18) (cid:19)
(cid:103)
1
Then, taking Ψ = as test function in (6.24) and using (6.25) and (6.28) we obtain
0
(cid:18) (cid:19)
1 (cid:91) 1
∂ q, + G n dΓ = 0, (6.29)
t
0 · · 0
(cid:28) (cid:18) (cid:19)(cid:29)h K (cid:90) ∂K (cid:18) (cid:19)
(cid:88)∈C
But on a periodic domain
(cid:91) 1 (cid:92) 1
G n dΓ = G n [[ ]]dΓ = 0
f
· · 0 · · 0
K (cid:90) ∂K (cid:18) (cid:19) f (cid:90) f (cid:18) (cid:19)
(cid:88)∈C (cid:88)∈F
yielding
1
∂ q, = 0.
t
0
(cid:28) (cid:18) (cid:19)(cid:29)h
0
The result with Ψ = is obtained similarly. This ends the proof.
1
(cid:18) (cid:19)

CHAPTER 6. EXTENSION OF THE RAVIART-THOMAS STAGGERED SCHEME TO
COMPRESSIBLE BAROTROPIC FLOWS 162
6.4 Low Mach number analysis
We are now interested in concluding on the low Mach number behaviour of the proposed
Raviart-Thomas staggered scheme. It is very clear by now that the de Rham discretization
is mainly motivated by the existence of an underlying Hodge-Helmholtz decomposition. This
property was used in chapter 5 to show long time consistency on the first order linear wave
system. In particular, in this section, we finally show how this long time consistency of the
Raviart-Thomas staggered scheme implies in some sense (to be defined) low Mach number
accuracy. Let us first clarify what we mean by a low Mach number accurate numerical scheme
:
u
Definition 6.4.1 (Low Mach number accurate scheme). Let M := 0 be the Mach number
c
0
and let the expansion
N
ϕ˜(x˜,t˜,M) = Mnϕ˜(n)(x˜,t˜,τ)+ (MN+1), (6.30)
O
n=0
(cid:88)
Suppose that the initial and boundary conditions are well-prepared in the sense of Defini-
tion 1.2.3 and assume the following hypothesis
Hyp 1 a low Mach number asymptotic analysis using the expansion (6.30) on the discrete un-
knowns shows that the discrete first order pressure is coupled with the zeroth-order mo-
mentum through a discretization consistent with a wave system,
Hyp 2 the resulting discretization of this wave system is able to capture the accurate long time
limit.
Then if Hyp 1 is true, we say that a numerical scheme is low Mach number accurate if Hyp 2
is verified.
For the following, we use the characteristic lengths ρ , x and u , c2 = p(ρ ) and the
0 0 0 0 (cid:48) 0
dimensionless variables defined in chapter 1, subsection 1.2.2. In addition, we denote
q K σ (cid:63)σ c u
q := K := | | σ˜ := | | (cid:63)σ˜ := | | c˜ := max u˜ := max.
ρ u | | x2 | | x | | x max c max u
0 0 0 0 0 0 0
These(cid:101)notations are (cid:101) used in the following basic results:
Lemma 6.4.1 (Dimensionless density equation). The dimensionless mass equation of the
staggered scheme (6.23a) on a periodic domain reads
u˜ c˜
K˜ ∂ ρ˜ + σ˜ ε (σ)q = | σ | + max σ˜ ε (σ)[[ρ˜]] (6.31)
| |
t˜ K
| |
K σ
2 2M | |
K σ
σ ∂K σ ∂K(cid:18) (cid:19)
(cid:88)∈ (cid:88)∈
(cid:101)
Proof. Usingthedefinitionofthecharacteristicquantitiesρ ,x andu werewrite(6.13)under
0 0 0

CHAPTER 6. EXTENSION OF THE RAVIART-THOMAS STAGGERED SCHEME TO
COMPRESSIBLE BAROTROPIC FLOWS 163
the following form
ρ ρ u x ρ x u˜ c˜
K 0∂ ρ˜ + 0 0 0 σ˜ ε (σ)q = 0 0 u | σ | +c max σ˜ ε (σ)[[ρ˜]] . (6.32)
| |t t˜ K x2 | | K σ x2 0 2 0 2 | | K σ
0 0 σ ∂K 0 σ ∂K(cid:18) (cid:19)
(cid:88)∈ (cid:88)∈
(cid:101) (cid:101)
t
0
Multiplying (6.32) by we obtain the result.
ρ x
0 0
Now onto the momentum equation,
Lemma 6.4.2 (Dimensionless momentum equation). The dimensionless staggered momentum
equation (6.23b) is given in a periodic domain by:
1
D˜ ∂ q +C(q) + σ˜ [[p˜]]
| σ | t˜ σ σ γM2| | σ
(6.33)
= u˜ max + c˜ max (cid:101) σ˜ [[d (cid:94) iv (cid:101) (q)]] + u˜ max ˜(cid:63)σ [[( (cid:94) ∇ ) q]] ,
2 2M | | σ 2 | | ⊥ ∗ ⊥σ
(cid:18) (cid:19)
c2ρ (cid:101) (cid:101)
with γ := 0 0,
p
0
q q q˜
C(q˜) := ⊗ : ∇ Ψ dx˜+ q˜ n Ψ dΓ˜
x˜ σ σ
K (σ)(cid:20)
−
(cid:90)
K˜ ρ˜
(cid:90)
∂˜K ·
(cid:26)(cid:26)
ρ˜
(cid:27)(cid:27)
·
(cid:21)
(cid:88)∈C (cid:101) (cid:101)
Proof. Using characteristic quantities ρ ,u ,p and x we have by (6.22)
0 0 0 0
ρ u x2
0 0 0 D˜ ∂ q +ρ u2x C(q) +x p σ˜ [[p˜]]
t | σ | t˜ σ 0 0 0 σ 0 0 | | σ
0
(6.34)
= x ρ u u u˜ max +c c˜ max (cid:101) σ˜ [[d (cid:94) iv(q)]] + (cid:101) u x ρ u u˜ max (cid:63)σ˜ [[( (cid:94) ∇ ) q]] ,
0 0 0 0 2 0 2 | | σ 0 0 0 0 2 | | ⊥ ∗ ⊥σ
(cid:18) (cid:19)
ρ u x2 ρ u x2(cid:101) (cid:101)
Dividing (6.34) by 0 0 0 and using 0 0 0 = ρ u2x we obtain the result.
t t 0 0 0
0 0
Combining these results we get the following proposition:
Proposition 6.4.1 (Discrete wave system coupling the first order density and the zeroth order
velocity). Suppose that the initial conditions are well-prepared in the sense of Definition 1.2.3,
using the asymptotic expansion (6.30) on each variables, the staggered scheme (6.23a) (6.23b)
yields that
ρ˜(0) is independent of x˜,t˜,τ˜, (6.35)
and it is asymptotically consistent, on a periodic domain, with the following discrete wave

CHAPTER 6. EXTENSION OF THE RAVIART-THOMAS STAGGERED SCHEME TO
COMPRESSIBLE BAROTROPIC FLOWS 164
system coupling the first order pressure and the zeroth order momentum:
K˜ ∂ p˜ (1) +γ c˜(0) 2 σ˜ ε (σ)ρ˜(0)u˜(0) =
c˜(0)
σ˜ ε (σ)[[p˜(1)]] ,
| | τ K | | K σ 2 | | K σ
 σ ∂K σ ∂K
 (cid:0) (cid:1) (cid:88)∈ (cid:88)∈



  1 c˜(0) (cid:94)
D˜ ∂ ρ˜(0)u˜ (0) + σ˜ [[ p˜(1)]] = σ˜ [[div(ρ˜(0)u˜(0))]] .
σ τ σ σ σ
| | γ| | 2 | |





Proof. Weinject in the dimensionless mass equation the expansion on Mach number (6.30) in
1
each variables. From that, (6.33) yields at order
O M2
(cid:18) (cid:19)
σ σ˜ [[p˜(0)]] = 0, (6.36)
σ
∀ ∈ F | |
so that p(0) is uniform in space . As a consequence ρ˜(0) is also uniform. Using this information
1
with (6.31) at order we obtain
O M
(cid:18) (cid:19)
c˜(0)
K˜ ∂ ρ˜ (0) = σ˜ [[ρ˜(0)]] = 0, (6.37)
| | τ K 2 | | σ
σ ∂K
(cid:88)∈
so ρ˜(0)(x˜,t˜,τ˜) = ρ˜(0)(t˜). Because of the well-preparedness Definition 1.2.3, the latter is time-
independent, yielding (6.35) and, invoking again (6.31), but here at order (1) we obtain
O
c˜(0) d
K˜ ∂ ρ˜ (1) + σ˜ ε (σ)ρ˜(0)u˜(0) = σ˜ ε (σ)[[ρ˜(1)]] K˜ ρ˜ (0). (6.38)
| | τ K | | K σ 2 | | K σ −| |dt˜ K
σ ∂K σ ∂K
(cid:88)∈ (cid:88)∈
=0
1 (cid:124) (cid:123)(cid:122) (cid:125)
Finally, using (6.33) at order we get
O M
(cid:18) (cid:19)
1 c˜(0) (cid:94)
D˜ ∂ ρ˜(0)u˜ (0) + σ˜ [[p(1)]] = σ˜ [[div(ρ˜(0)u˜ (0) )]] . (6.39)
σ τ σ σ σ σ
| | γ| | 2 | |
Since p˜ = p˜(ρ˜ ) we can adapt (2.2.1) to remark that
K K
2
p˜ (1) = p˜(ρ˜ (0) )ρ(1) = p˜(ρ˜(0))ρ(1), γ c˜(0) = p(ρ˜(0)) (6.40)
K (cid:48) K K (cid:48) K (cid:48)
(cid:16) (cid:17)
where the last equality stands because ρ˜(0) is constant. Then gathering (6.39) and (6.38)
multiplied by p˜(ρ˜(0)) yields, because of (6.40), the result.
(cid:48)
Gathering this proposition with the definition of the boundary fluxes we obtain the
following results on the numerical scheme:

CHAPTER 6. EXTENSION OF THE RAVIART-THOMAS STAGGERED SCHEME TO
COMPRESSIBLE BAROTROPIC FLOWS 165
Theorem 6.4.1 (Low Mach number accuracy of the de Rham scheme). The Raviart-Thomas
staggered scheme is low Mach number accurate in the sense of Definition 6.4.1.
Proof. It is noted in [43, Section 2.1, 2.2] that the Steger-Warming Definition 6.2.2 and (Roe-
)wall fluxes Definition 6.2.1 are consistent in low Mach number asymptotic with the fluxes
γ c˜(0) 2 (ρ˜u˜)(0) n γ c˜(0) 2
(ρ˜u˜)(0)+(ρ˜
b
u˜
b
)(0)
n+
c˜(0)
(p˜(1) 0)
1 · =  2 · 2 − , (6.41)
 (cid:0) γ (cid:1)p(1)n  (cid:0) 1p˜(1 (cid:1) )+0 n+ c˜(0) ((ρ˜u˜)(0) (ρ˜ b u˜ b )(0)) nn
inlet/outlet  γ 2 2 − · 
   
 
and
γ c˜(0) 2 (ρ˜u˜)(0) n 0
1 · = 1 , (6.42)
 (cid:0) (cid:1)p(1)n   p(1)n+c˜(0)(ρ˜u˜)(0) nn 
γ γ ·
wall
   
which are the boundary fluxes we use for the wave system discretization developed in chapter 5
with 1 := γ c˜(0) 2 and κ := 1 such that c2 := κ 0 = c˜(0) 2 . Remarking this implies that in
ρ 0 γ 0 ρ
0 0
order to check(cid:0) the(cid:1)low Mach number consistency it is s(cid:0)uffic(cid:1)ient to take a look at the interior.
In this regards, Proposition 6.4.1 shows that the de Rham discretization is indeed consistent
with a wave sytem on the first order pressure and the zeroth order momentum. This answers to
point Hyp 1 , then chapter 5 shows that the fully-discrete scheme, combined with the wall and
inlet/oulet boundary fluxes (6.41), (6.42) converges to the accurate long time limit for Euler
Explicit, ImEx (under CFL conditions) and Implicit method when
u˜ (0) ndΓ˜ = 0,
b ·
(cid:90)
∂Ω
which stands under the condition of well-preparedness; showing that Hyp 2 is verified: this
concludes the proof.
6.5 Discussion on some preexisting staggered schemes
In this section, we conduct a similar discussion to the ones made in chapter 3 and chapter 5 on
the proposed Raviart-Thomas staggered scheme and how it compares to preexisting staggered
schemes.
6.5.1 Low Mach number behaviour
The low Mach number behaviour of multiple preexisting staggered schemes was adressed in
the previous chapters (chapter 3, section 3.5 and chapter 5, section 5.5), we recall here the
fundamental takeaways:

CHAPTER 6. EXTENSION OF THE RAVIART-THOMAS STAGGERED SCHEME TO
| COMPRESSIBLE |     | BAROTROPIC |     |     | FLOWS |     |     |     |     |     | 166 |
| ------------ | --- | ---------- | --- | --- | ----- | --- | --- | --- | --- | --- | --- |
1) First, we recall that most of the discretizations found [82, 69, 49, 66, 81, 65] do not
benefit from a Hodge-Helmholtz decomposition (HHD) (or at least we are not aware of
| such | property). |     |     |     |     |     |     |     |     |     |     |
| ---- | ---------- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
2) Inparallel,lowMachnumberasymptoticanalysisshowsthatingeneralclassicalstaggered
discretizations [66, 65] lead to ImEx fully centred or Implicit schemes; the former is not
energystableonthewavesystemwhilethelatterisnaturally. Bothsufferfromoscillations
which appear to be damped by adding diffusion. Note that, as largely discussed in this
thesis, in multiple space dimensions these diffusion operators should be chosen carefully,
| especially |     | if they | are | carried | by the | acoustic | scale. |     |     |     |     |
| ---------- | --- | ------- | --- | ------- | ------ | -------- | ------ | --- | --- | --- | --- |
3) Lastly,discretizationsresultingonexplicitschemes[69]onthewavesystemwereidentified
in the litterature. However, the diffusion chosen, see section 5.5, chapter 5, does not seem
to preserve the stationary state as needed for us to conclude on a potential accurate low
| Mach | number |     | behaviour. |     |     |     |     |     |     |     |     |
| ---- | ------ | --- | ---------- | --- | --- | --- | --- | --- | --- | --- | --- |
By contrast, the discretization we propose results through low Mach number asymptotics on a
long time consistent (as shown in Theorem 6.4.1) wave system since it checks all the necessary
features: the ability to identify the limit is yielded by a natural HHD on the velocity space
RT1(Ω) and, in parallel, a well-designed stationarity preserving diffusion enables to obtain
energy dissipation properties while preserving the invariance of the divergence-free part of the
HHD.
| 6.5.2 | Discrete | entropy |     | dissipation |     |     |     |     |     |     |     |
| ----- | -------- | ------- | --- | ----------- | --- | --- | --- | --- | --- | --- | --- |
Entropy types inequalites or stability properties in the staggered framework have been ob-
tained for multiple set ups through the Crouzeix-Raviart/Rannacher-Turek velocity staggering
([55, 56, 57, 66, 67, 68, 69]). It is founded on the postulate that the convection term has a
particular form which links the mass flux to the convection flux. Since the framework proposed
here is quite different, we are not able to replicate these arguments for the convection term
| (6.22) | nor we | are, more | generally, |     | able to | prove | entropy dissipation. |     |     |     |     |
| ------ | ------ | --------- | ---------- | --- | ------- | ----- | -------------------- | --- | --- | --- | --- |
Nonetheless we note that the discrete grad-div introduced in this thesis can also readily
be defined in the Crouzeix-Raviart/Rannacher-Turek context using their natural discrete
divergence;
1
|     | div(u) |     | :=  |        | σ u(x ) | n   | with x barycenter |     | of the face | σ   |     |
| --- | ------ | --- | --- | ------ | ------- | --- | ----------------- | --- | ----------- | --- | --- |
|     |        | K   | K   |        | σ       | K,σ | σ                 |     |             |     |     |
|     |        |     |     |        | | |     | ·   |                   |     |             |     |     |
|     |        |     | |   | | σ ∂K |         |     |                   |     |             |     |     |
(cid:88)⊂
| so that | the grad-div |     | diffusion | can        | be adapted | as       |          |     |     |     |     |
| ------- | ------------ | --- | --------- | ---------- | ---------- | -------- | -------- | --- | --- | --- | --- |
|         |              |     |           | (cid:94)   |            | (cid:94) | (cid:94) |     |     |     |     |
|         |              |     | σ         | [[div(u)]] | n :=       | (div(u)  | div(u)   | )n  |     |     |     |
|         |              |     |           |            | σ K,σ      |          | Lσ       | K   | K,σ |     |     |
|         |              |     | | |       |            |            |          | −        |     |     |     |     |
With this operator, we proposed in Appendix D an Euler fully explicit time integration using
the classical staggered set up (Crouzeix-Raviart/Rannacher-Turek staggered velocity) that is
globally entropic under CFL condition. While it is not the core of this thesis, this proposition

CHAPTER 6. EXTENSION OF THE RAVIART-THOMAS STAGGERED SCHEME TO
| COMPRESSIBLE |     | BAROTROPIC |     | FLOWS |     |     |     | 167 |
| ------------ | --- | ---------- | --- | ----- | --- | --- | --- | --- |
shows that other staggered discretizations can benefit from the grad-div stabilization we intro-
duced. Again, we recall that this choice of discretization, to our knowledge, does not yield any
Hodge-Helmholtz decomposition and as a consequence we are not able to identify the long time
| limit on | the wave | system.       |     |           |     |     |     |     |
| -------- | -------- | ------------- | --- | --------- | --- | --- | --- | --- |
| 6.5.3    | Other    | computational |     | questions |     |     |     |     |
One of the most remarkable advantages of the Raviart-Thomas staggered scheme, apart from
its strong low Mach number theoretical basis, resides in the reduction of degrees of freedom
by comparison with Crouzeix-Raviart type discretizations. Indeed these finite element spaces
intrinsically use the full velocity vector at the face, yielding, in the bidimensional case, a
number of velocity d.o.fs that is twice the number of faces and in three space dimensions, three
times the number of faces. In contrast, the Raviart-Thomas based discretization only needs
the normal component whether it be for one, two or three dimensions simulations. Naturally,
it implies that the number of velocity/momentum discrete unknowns is, on the one hand,
independent of the space dimension (even though a three dimensional domain will generally
demand a finer mesh) and on the other hand, equal to the number of faces. As a consequence,
it is three times smaller than in the Crouzeix-Raviart/Rannacher-Turek framework in 3
dimensions. Again, the underlying space of discretization does not enable the identification of
| the long | time | limit on | the wave | system. |     |     |     |     |
| -------- | ---- | -------- | -------- | ------- | --- | --- | --- | --- |
As mentioned in chapter 5, adapting the scheme for the three space dimensions case,
does not require any additional effort as far as the acoustic part is concerned. In general,
for the proposed Raviart-Thomas staggered scheme to be able to tackle the case of three
space dimensions with analogous operators, the following adaptations are needed: as noted in
chapter 5, the N´ed´elec-Raviart-Thomas in 3D differs from the 2D case, and reads
|     |     |     |        | ∇ N1(Ω) |     | rot RT1(Ω) | div       |     |
| --- | --- | --- | ------ | ------- | --- | ---------- | --------- | --- |
|     |     |     | cG1(Ω) |         |     |            | − dG0(Ω). |     |
|     |     |     |        | −→      | −→  |            | −→        |     |
So that in three space dimensions the Hodge-Laplacian of the Raviart-Thomas space is
| defined | as a | discrete equivalent |     | of  |         |          |     |     |
| ------- | ---- | ------------------- | --- | --- | ------- | -------- | --- | --- |
|         |      |                     |     | ∆u  | = ∇divu | rotrotu, |     |     |
−
instead of
|     |     |     |     | ∆u  | = ∇divu+∇ | ⊥   | curlu, |     |
| --- | --- | --- | --- | --- | --------- | --- | ------ | --- |
in 2D. So the diffusion term should be adapted in 3 space dimensions. Since the grad-div
operator is identically defined in 3 or 2 dimensions, we have to compute in particular for any
| RT1(Ω) |     | N1(Ω), |     |              |           |           |            |     |
| ------ | --- | ------ | --- | ------------ | --------- | --------- | ---------- | --- |
| Ψ      |     | and Φ  |     | the operator |           |           |            |     |
| ∈      |     | ∈      |     |              |           |           |            |     |
|        |     |        |     | (rot) ∗ Ψ,Φ  | L2(Ω)3    | :=        | Ψ,rotΦ h . |     |
|        |     |        |     | (cid:104)    | (cid:105) | (cid:104) | (cid:105)  |     |
Remark 6.5.1. Note that since we have defined the Hodge-Laplacian on the Raviart-Thomas
| space, | we should | be able | to treat | viscosity | terms | of the | form |     |
| ------ | --------- | ------- | -------- | --------- | ----- | ------ | ---- | --- |

CHAPTER 6. EXTENSION OF THE RAVIART-THOMAS STAGGERED SCHEME TO
COMPRESSIBLE BAROTROPIC FLOWS 168
∇div(u) and ∆u,
which are readily defined in the context developed here. In practice we have not examined this
possibility.
6.6 Numerical results
We are now interested in testing the properties of Raviart-Thomas staggered scheme, with a
particular emphasis on the conservation and the low Mach number behaviour. The implement-
ation is done in the same simulation C++ code SolverLab used to develop the staggered
discretization on the wave system (see details on this implementation in the introduction of
section 5.6, chapter 5). The implementation of the scheme for Euler equations inherits, for
practical reasons, from the class of the wave system introduced in chapter 5. Because of this
inheritance, most of the discrete differential operators are reusable on the Euler barotropic
system. Nonetheless, in order to treat it fully, additional operators were coded, of which, the
convection term, a class managing the pressure law and the gradients of the Raviart-Thomas
basis functions on quadrangular, but also for future applications, for triangular meshes. For
post-processing of the momentum we used the interpolation formula (5.113) implemented for
chapter 5 that takes a face-based field and gives a cell-based interpolation. The numerical in-
tegration of the convection term is done with a trapezoid formula on the reference element (see
details on this choice in Appendix F) while two time-integrations are available in SolverLab:
either the Euler explicit time stepping or the Semi-Implicit time stepping (where only acoustic
terms are implicited). In both cases CFL conditions are required, we constraint the explicit
time stepping with a CFL of the type δt h whereas in semi-implicit it is constrained by
≈ u+c
a CFL of the form δt h , where h is length | rel | ated to the size of the mesh. Unless specified,
≈ u
the numerical results sho|w|n are obtained with the ImEx Euler time stepping. Finally, the law
used for the pressure is p(ρ) = ρ2.
6.6.1 1D Riemann problems
In[121]itisshownthatnon-conservativeschemesconvergetoawrongsolutionandthatthegap
with the accurate solution is supported on the space around the discontinuities. Subsequently,
in order to challenge the conservativity of the scheme, we first investigate the behaviour of
the scheme in configurations where such discontinuities arise: classical tests which encapsulate
these conditions are Riemann problems. On the Euler barotropic system an exact solution can
be computed numerically by solving a Newton method (see Appendix E) which enables us to
compare the scheme with the exact solution from the Newton solution. All the simulations are
run with final time T := 0.03 and 400 cells. In Figure 6.1, Figure 6.2 the numerical and exact
solutions of three Riemann problems are shown, solved in explicit time stepping Figure 6.1 and
in ImEx time stepping Figure 6.2.
It is known that in the case of the Euler barotropic system the solution of a Riemann problem

CHAPTER 6. EXTENSION OF THE RAVIART-THOMAS STAGGERED SCHEME TO
| COMPRESSIBLE |     | BAROTROPIC |     |     | FLOWS |     |     |     |     |     |     | 169 |
| ------------ | --- | ---------- | --- | --- | ----- | --- | --- | --- | --- | --- | --- | --- |
consists of the combinaison of shocks and rarefaction waves by contrast with the full Euler
system where contact discontinuities constitutes a supplementary type of elementary solutions.
Hence, the solution of a Riemann problem, if it exists, can be expressed as the composition of
rarefaction waves and shocks. We put a special emphasis on cases where the solution contains
a shock: a 1 shock wave and a 2 rarefaction (Figure 6.1a in explicit, Figure 6.2a in semi-
|     | −   |     |     |     | −   |     |     |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
implicit) or, a 1 rarefaction wave and a 2 shock (Figure 6.1b in explicit, Figure 6.2b in semi-
|     |     | −   |     |     |     | −   |     |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
implicit), a 1 shock and a 2 shock wave (Figure 6.1c in explicit, Figure 6.2c in semi-implicit).
|     | −   |     |     | −   |     |     |     |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
Ineachoftheseteststheshocksarewell-resolvedandnogapappearsaroundthediscontinuities
between the exact solution and the one obtained with the numerical scheme. Unsurprisingly,
| the semi-implicit |          | version |            | is a litlle | bit | more diffusive. |     |     |     |     |     |     |
| ----------------- | -------- | ------- | ---------- | ----------- | --- | --------------- | --- | --- | --- | --- | --- | --- |
| 6.6.2             | Cylinder |         | scattering |             |     |                 |     |     |     |     |     |     |
In this section, we explore the low Mach number properties of the numerical scheme with the
R2,
cylinder test case: the domain is, in polar coordinates, Ω = [0,2π] [r ,r ] we simulate
|     |     |     |     |     |     |     |     |     | ×   | 0 1 | ⊂   |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
|     |     |     |     |     | r   | 0.8 | r   |     |     |     |     |     |
the scattering of a flow in Ω with 0 = and 1 = 6 where wall boundary conditions are
imposed on the inner circle, namely on r = r , whereas inlet/outlet boundary conditions are
0
imposedontheoutercircler = r . Concerninginlet/outletconditions,wehaveρ = 2andwith
|     |     |     |     |     | 1   |     |     |     |     |     | b   |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
u = M p(ρ ), u = (u ,0)t are the given density and velocity and M the Mach number, is
| 0 b         | (cid:48) | 0 b      |       | 0   |         |          |     |              |     | b   |     |     |
| ----------- | -------- | -------- | ----- | --- | ------- | -------- | --- | ------------ | --- | --- | --- | --- |
| a parameter | of       | the test | case. | The | initial | velocity | is  | simply given | by  |     |     |     |
(cid:112)
|               |     |       |           |          |      | ρ            | (x) =     | ρ ,       |       |          |     |        |
| ------------- | --- | ----- | --------- | -------- | ---- | ------------ | --------- | --------- | ----- | -------- | --- | ------ |
|               |     |       |           |          |      | 0            |           | b         |       |          |     | (6.43) |
|               |     |       |           |          |      | u            | (x) =     | 0.        |       |          |     |        |
|               |     |       |           |          |      | (cid:26) 0   |           |           |       |          |     |        |
| Conveniently, |     | a     | reference | solution |      | for this     | test case | is known  | and   | is:      |     |        |
|               |     |       |           |          |      | r2           | 2         | r2        | r4    | 2        |     |        |
|               |     | ρ     | (r,θ)     | =        | ρ +ρ |              | 1         | 0 cos(2θ) | 0     | M2,      |     |        |
|               |     | exact |           |          | b    | b            |           |           |       |          | b   |        |
|               |     |       |           |          |      | r2           | r2        | r2        | − 2r4 |          |     |        |
|               |     |      |           |          |      | (cid:18) 1 − | 0(cid:19) | (cid:18)  |       | (cid:19) |     |        |
 

|     |     |   |     |     |     |     |     | r 2       |     |     |     |        |
| --- | --- | --- | --- | --- | --- | --- | --- | --------- | --- | --- | --- | ------ |
|     |     |    |     |     |     |     |     | 0 cos(2θ) |     |     | .   | (6.44) |
|     |     |   |     |     |     |     |     | 1         |     |     |     |        |
|     |     |   |     |     |     | r2  |     | − r2      |     |     |     |        |

|     |     |     | u   | (r,θ) | =   | u    | 1  |     |     |    |     |     |
| --- | --- | --- | --- | ----- | --- | ---- | --- | --- | --- | --- | --- | --- |
|     |     |     |     | exact |     | b r2 | r2  |     |     |     |     |     |
2
|     |     |    |     |     |     | 1 − | 0  | r         |     |    |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --------- | --- | --- | --- | --- |
|     |     |   |     |     |     |     |     | 0 sin(2θ) |     |     |     |     |
|     |     |    |     |     |     |     |    |           |     |    |     |     |
|     |     |   |     |     |     |     |    | − r 2     |     |    |     |     |
|     |     |    |     |     |     |     |    |           |     |    |     |     |
 

An initial qualitative result is obtained by looking, as for the wave system, at the isocontours
|     |     |     |     |     |     |     |     | T 100s; |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | ------- | --- | --- | --- | --- |
of the velocity norm on a stationary state obtained at = in Figure 6.3 are shown the
isocontours of the reference solution in comparison with the discrete velocity obtained with the
scheme on a mesh with n = 20 subdivisions in the radial direction and n = 32 subdivisions
|     |     |     |     | r   |     |     |     |     |     |     | θ   |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
in the angular direction. The left one on Figure 6.3, obtained with the scheme (6.23a)(6.23b),
demonstrates the relevant qualitative behaviour, while the isocontours on the right are ob-
tained with (6.23a)(6.23b) enriched with a full Hodge-Laplacian diffusion Definition 5.6.2 on

CHAPTER 6. EXTENSION OF THE RAVIART-THOMAS STAGGERED SCHEME TO
| COMPRESSIBLE |     | BAROTROPIC |     | FLOWS |     |     |     |     |     | 170 |
| ------------ | --- | ---------- | --- | ----- | --- | --- | --- | --- | --- | --- |
theacousticscale. Asseeninsubsection5.6.2, chapter5, theseresultsillustratetheimportance
| of putting | the | appropriate | diffusion | operator |     | on the | acoustic | scale. |     |     |
| ---------- | --- | ----------- | --------- | -------- | --- | ------ | -------- | ------ | --- | --- |
This is not sufficient to conclude on the accuracy of the scheme, especially since a reference
solution (6.44) is available on this test case. In the aim to better characterize this appropriate
L2
behaviour, Figure 6.4 shows, for a fixed reference Mach number M b = 1e 4, the error
|     |     |     |     |     |     |     |     |     | −   | −   |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
L2
(with respect to this reference solution) of the discrete density and the error of each com-
−
ponent of the discrete momentum with respect to mesh sizes given by n = 6 n = 12,
|     |     |     |     |     |     |     |     |     | r   | θ   |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
×
n = 10 n = 20, n = 20 n = 40, n = 40 n = 80. The rate of convergence to the
| r         | θ        |        | r       | θ        | r    |     | θ   |     |     |     |
| --------- | -------- | ------ | ------- | -------- | ---- | --- | --- | --- | --- | --- |
|           | ×        |        |         | ×        |      |     | ×   |     |     |     |
| reference | solution | (6.44) | is very | close to | one. |     |     |     |     |     |
Finally,inFigure6.5,theevolutionoftheL2 normofthediscretedivergenceofthemomentum
−
| and the | L2 norm | of  | the discrete | gradient | of  | the density |     |     |     |     |
| ------- | ------- | --- | ------------ | -------- | --- | ----------- | --- | --- | --- | --- |
−
|           | 2           |     | [[ρ˜]]2 |             |       |        |           |           | 2 (div(q˜))2 |     |
| --------- | ----------- | --- | ------- | ----------- | ----- | ------ | --------- | --------- | ------------ | --- |
| ∇         | h ρ˜ :=     | σ   |         | (= ( div)   | ∗ ρ˜, | ( div) | ∗ ρ˜ h    | ), divq˜  | :=           | dx  |
| (cid:107) | (cid:107)L2 | |   | | σ     | (cid:104) − |       | −      | (cid:105) | (cid:107) | (cid:107)L2  |     |
(cid:90)Ω
σ
(cid:88)∈F
areplottedwithrespecttotheMachnumberonafixedmeshofsize10 20. Figure6.5confirms
×
| that, | as it should; |     |     |           |     |         |     |      |     |     |
| ----- | ------------- | --- | --- | --------- | --- | ------- | --- | ---- | --- | --- |
|       |               |     |     | ∇ρ˜= (M2) |     | div(q˜) | =   | (M). |     |     |
|       |               |     |     | O         |     |         |     | O    |     |     |
6.6.3 Propagation of a low Mach number acoustic wave through a stationary
|     | low Mach |     | number | vortex |     |     |     |     |     |     |
| --- | -------- | --- | ------ | ------ | --- | --- | --- | --- | --- | --- |
Finally, we examine the ability of the de Rham scheme to take on incompressible problems
while propagating waves at low Mach number. It is motivated by [75], where it is shown
that some classical collocated solvers widly used for low Mach number flows approximation,
such as the Roe-Turkel scheme, are unable to propagate low Mach number acoustic waves.
The test case used in [75] to put in exergue this behaviour is the problem of a steady vortex
crossed by a low Mach number wave: the domain is a rectangle Ω := [ 0.1;1.1] [0,1]
|     |     |     |     |     |     |     |     |     | −   | ×   |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
meshed with 380 300 rectangular cells, with 380 subdivisions in the x direction and 300
|     |     | ×   |     |     |     |     |     |     | −   |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
y
in the direction. Periodic conditions are applied to the top and bottom boundaries and
−
inlet/outlet conditions are applied to the left and right boundaries. The initial condition is the
sum of two compactly supported functions carrying each a particular behaviour. On the one
hand the incompressible phenomena will be carried by a vortex centred in (x ,y ) := (0.5,0.5)
c c
| characterized | by  | a reference |     | Mach number | M   | :   |     |     |     |     |
| ------------- | --- | ----------- | --- | ----------- | --- | --- | --- | --- | --- | --- |
ref

CHAPTER 6. EXTENSION OF THE RAVIART-THOMAS STAGGERED SCHEME TO
| COMPRESSIBLE |     |         | BAROTROPIC |      |     | FLOWS |     |     |     |     | 171 |
| ------------ | --- | ------- | ---------- | ---- | --- | ----- | --- | --- | --- | --- | --- |
|              | √(x | xc)2+(y |            | yc)2 |     |       |     |     |     |     |     |
For r¯= − − 1 this function is equal to 0, else it is defined as:
0.45
≥
2
|     |     |     |     |        |     |     |     | M ref |     | 2α2   |     |
| --- | --- | --- | --- | ------ | --- | --- | --- | ----- | --- | ----- | --- |
|     |     |     |     | ρ      | (x) | =   | ρ 1 |       | e   | ,     |     |
|     |     |     |     | vortex |     |     | 0   |       |     | −1−r¯ |     |
− λ
|     |     |     |    |     |     |     | (cid:18) | (cid:18) max(cid:19) |     | (cid:19) |        |
| --- | --- | --- | --- | --- | --- | --- | -------- | -------------------- | --- | -------- | ------ |
|     |     |     |    |     |     |     |          |                      |     |          | (6.45) |

 
|     |     |     |    |        |     |     | y                 |     | 2α  | α2       |     |
| --- | --- | --- | --- | ------ | --- | --- | ----------------- | --- | --- | -------- | --- |
|     |     |     |    | u      | (x) | u   |                   |     |     | e −1−r¯, |     |
|     |     |     |     | vortex |     | =   | 0                 |     |     |          |     |
|     |     |     |     |        |     |     | x                 | λ r | (1  | r¯2)     |     |
|     |     |     |     |        |     |     | (cid:18) (cid:19) | max | 0   |          |     |
|     |     |     |   |        |     |     | −                 |     | −   |          |     |

 
|     | ρ   | u   | M  |     | p(ρ α |     | λ   |     |     |     |     |
| --- | --- | --- | --- | --- | ----- | --- | --- | --- | --- | --- | --- |
where 0 = 2, 0 := ref (cid:48) 0 ), = 2. max is chosen so that the maximal velocity norm is
u :
0
|     |     |     |     | (cid:112) |     |     | α2  |     |     |     |     |
| --- | --- | --- | --- | --------- | --- | --- | --- | --- | --- | --- | --- |
2αr¯
|     |     |     | λ   | =   | m ax  | e−1−r¯m | 2 r¯   | :=  | α2+       | 1+α4. |     |
| --- | --- | --- | --- | --- | ----- | ------- | ------ | --- | --------- | ----- | --- |
|     |     |     | max |     |       |         | ax max |     |           |       |     |
|     |     |     |     | 1   | r¯ 2  |         |        |     | −         |       |     |
|     |     |     |     |     | − max |         |        |     | (cid:113) |       |     |
(cid:112)
On the other hand, the acoustic phenomena is supported by a low Mach number acoustic wave
computed thanks to the Riemann invariants of the system: for x¯ = x/0.05 1 it is equal to
|     |            |         |     |     |     |     |     |     |     | | | | | ≥ |     |
| --- | ---------- | ------- | --- | --- | --- | --- | --- | --- | --- | --------- | --- |
| 0,  | else it is | defined | as: |     |     |     |     |     |     |           |     |
1
|     |     |     |     | ρ    | (x) = |     | ρ 1+M    |     | e 1 −1−x¯2 | ,        |     |
| --- | --- | --- | --- | ---- | ----- | --- | -------- | --- | ---------- | -------- | --- |
|     |     |     |     | wave |       |     | 0        | ref |            |          |     |
|     |     |     |    |      |       |     | (cid:18) |     |            | (cid:19) |     |
(6.46)

 
|     |     |     |    |      |       |     | 2        |      |     |            |     |
| --- | --- | --- | --- | ---- | ----- | --- | -------- | ---- | --- | ---------- | --- |
|     |     |     |    | u    | (x) = |     | p(ρ      | (x)  |     | p(ρ ) .    |     |
|     |     |     |     | wave |       |     | (cid:48) | wave |     | (cid:48) 0 |     |
|     |     |     |     |      |       | γ   | 1        |      | −   |            |     |
−
|     |     |     |   |     |     |     | (cid:16)(cid:112) |     |     | (cid:112) (cid:17) |     |
| --- | --- | --- | --- | --- | --- | --- | ----------------- | --- | --- | ------------------ | --- |

| Then | the | initial | condition  | is  | defined | as  |         |        |      |     |     |
| ---- | --- | ------- | ------------ | --- | ------- | --- | ------- | ------ | ---- | --- | --- |
|      |     |         |              |     | ρ (x)   | :=  | ρ (x)+ρ |        | (x), |     |     |
|      |     |         |              |     | 0       |     | wave    | vortex |      |     |     |

|     |     |     |     |     | u   | (x) | u (x)+u |        | (x). |     |     |
| --- | --- | --- | --- | --- | --- | --- | ------- | ------ | ---- | --- | --- |
|     |     |     |     |     |  0 | =   | wave    | vortex |      |     |     |

At least three things are expected to be seen in order to confirm the accuracy of the
numerical scheme in these conditions: firstly, it is expected that the vortex remains stationary.
Secondly, the low Mach number wave should be purely transported, and finally the low Mach
| number | wave | must | not | interact | with | the | vortex when | crossing |     | it. |     |
| ------ | ---- | ---- | --- | -------- | ---- | --- | ----------- | -------- | --- | --- | --- |
In Figure 6.6, the Mach number obtained with the explicit scheme is shown at different times;
at the initial time the Mach number clearly illustrates the two compactly supported functions
constituting the low Mach number wave (6.46) and steady vortex (6.45), at time 0.25s the
acoustic wave is shown crossing the vortex and at time 0.5s it is clear that the wave has been
transportedthroughthesteadylowMachnumbervortexandhasleftitunaffected,asexpected.
| Finally, | at  | time | 3.5s | we show | that the | vortex | is preserved. |     |     |     |     |
| -------- | --- | ---- | ---- | ------- | -------- | ------ | ------------- | --- | --- | --- | --- |

CHAPTER 6. EXTENSION OF THE RAVIART-THOMAS STAGGERED SCHEME TO
| COMPRESSIBLE | BAROTROPIC |     | FLOWS |     | 172 |
| ------------ | ---------- | --- | ----- | --- | --- |
6.7 Conclusion
In this chapter we have extended the Raviart-Thomas staggered scheme to Euler barotropic
equations. Combining ideas from Discontinuous Galerkin methods with mass-lumping tech-
niques and diffusion operators defined thanks to the underlying finite element de Rham com-
plex we have derived a staggered scheme on the density and the momentum. Then with this
| discretization | we were able | to address | three core | points: |     |
| -------------- | ------------ | ---------- | ---------- | ------- | --- |
• Defining the staggered variable as the momentum, gathered with the DG inspired frame-
work, enabledustonaturallytreatvariablesandoperatorsscatteredatdifferentlocations
| of the | mesh. |     |     |     |     |
| ------ | ----- | --- | --- | --- | --- |
•
Conservation is achieved by using the Galerkin formalism on which the discretization
rely on. The choice of numerical diffusion and, for the convection term, the centered
flux, is compatible with the preservation of uniform fields, yielding naturally a global
| conservation | of the | momentum. |     |     |     |
| ------------ | ------ | --------- | --- | --- | --- |
•
By using the formal link between the low Mach number limit and the long time limit of
a wave system we show that the Raviart-Thomas staggered scheme is low Mach number
asymptotically consistent with the discrete wave system we have developed in chapter 5.
This discretization, provided that it verifies natural stability constraints, converges to the
accurate long time limit, yielding, for the Euler equations, a low Mach number accurate
| staggered | scheme. |     |     |     |     |
| --------- | ------- | --- | --- | --- | --- |
We challenged these properties with various numerical tests: Riemann problems for the ro-
bustness with respect to conservation and the cylinder scattering test for low Mach number
accuracy. Also, we have numerically illustrated that the Raviart-Thomas staggered scheme is
able to propagate a low Mach number acoustic wave through a steady incompressible vortex
| without damaging | it. |     |     |     |     |
| ---------------- | --- | --- | --- | --- | --- |

CHAPTER 6. EXTENSION OF THE RAVIART-THOMAS STAGGERED SCHEME TO
| COMPRESSIBLE | BAROTROPIC |        | FLOWS          |     |                  |     | 173 |
| ------------ | ---------- | ------ | -------------- | --- | ---------------- | --- | --- |
|              |            | (a) (ρ | L ,u L )=(1,2) | and | (ρ R ,u R )=(10, | 3)  |     |
−
|     |     | (b) (ρ | ,u )=(9,1.5) | and | (ρ ,u )=(1,  | 3)  |     |
| --- | --- | ------ | ------------ | --- | ------------ | --- | --- |
|     |     |        | L L          |     | R R          | −   |     |
|     |     | (c) (ρ | ,u )=(13,0)  | and | (ρ ,u )=(13, | 10) |     |
|     |     | L      | L            |     | R R          | −   |     |
fully explicit
Figure 6.1: Resolution of various 1d Riemann problems with time-integration

CHAPTER 6. EXTENSION OF THE RAVIART-THOMAS STAGGERED SCHEME TO
| COMPRESSIBLE | BAROTROPIC |        | FLOWS          |     |                  |     | 174 |
| ------------ | ---------- | ------ | -------------- | --- | ---------------- | --- | --- |
|              |            | (a) (ρ | L ,u L )=(1,2) | and | (ρ R ,u R )=(10, | 3)  |     |
−
|     |     | (b) (ρ | ,u )=(9,1.5) | and | (ρ ,u )=(1,  | 3)  |     |
| --- | --- | ------ | ------------ | --- | ------------ | --- | --- |
|     |     |        | L L          |     | R R          | −   |     |
|     |     | (c) (ρ | ,u )=(13,0)  | and | (ρ ,u )=(13, | 10) |     |
|     |     | L      | L            |     | R R          | −   |     |
semi-implicit
Figure 6.2: Resolution of various 1d Riemann problems with time-integration

CHAPTER 6. EXTENSION OF THE RAVIART-THOMAS STAGGERED SCHEME TO
| COMPRESSIBLE | BAROTROPIC | FLOWS              |        |     | 175 |
| ------------ | ---------- | ------------------ | ------ | --- | --- |
|              |            | Reference solution | (6.44) |     |     |
Grad-div diffusion on the acoustic scale, Hodge-Laplacian diffusion on the acoustic scale
|     | Figure | 6.3: Isocontours | of the velocity | norm |     |
| --- | ------ | ---------------- | --------------- | ---- | --- |

CHAPTER 6. EXTENSION OF THE RAVIART-THOMAS STAGGERED SCHEME TO
|     | COMPRESSIBLE |     |     | BAROTROPIC |     | FLOWS |     |     |     |     |     | 176 |
| --- | ------------ | --- | --- | ---------- | --- | ----- | --- | --- | --- | --- | --- | --- |
RT1
Semi-Implicit,
|     |     |     |     |     |        | Explicit, RT1 |     |     |     |     |     |     |
| --- | --- | --- | --- | --- | ------ | ------------- | --- | --- | --- | --- | --- | --- |
|     |     |     |     |     | 6 10−9 |               |     |     |     |     |     |     |
×
2L
k
|     |     |     |     |     | 4 10−9 |     |     |     |     |     |     |     |
| --- | --- | --- | --- | --- | ------ | --- | --- | --- | --- | --- | --- | --- |
|     |     |     |     | xe  | ×      |     |     |     |     |     |     |     |
ρ
−
|     |     |     |     |     | 3 10−9 | slope=0.9 |     |     |     |     |     |     |
| --- | --- | --- | --- | --- | ------ | --------- | --- | --- | --- | --- | --- | --- |
|     |     |     |     | h   | ×      |           |     |     |     |     |     |     |
ρ
k
|     |     |     |     |     | 2 10−9 |     |     |     |     |     |     |     |
| --- | --- | --- | --- | --- | ------ | --- | --- | --- | --- | --- | --- | --- |
×
|     |     |     |     |     | 10−1 | 2 10−1 | 3 1 0−1 4 10−1 | 6   | 10−1 | 100 |     |     |
| --- | --- | --- | --- | --- | ---- | ------ | -------------- | --- | ---- | --- | --- | --- |
|     |     |     |     |     |      | ×      | × ×            | ×   |      |     |     |     |
h
L2
|     |     |     |     |     |     | Error | on ρ |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | ----- | ---- | --- | --- | --- | --- | --- |
6 10− 5
×
|     |       | Semi-Implicit, |     | RT1 |     |     |          |     | Semi-Implicit, |     | RT1 |     |
| --- | ----- | -------------- | --- | --- | --- | --- | -------- | --- | -------------- | --- | --- | --- |
|     |       |                | RT1 |     |     |     |          |     |                | RT1 |     |     |
| 2L  |       | Explicit,      |     |     |     |     | 2L 6 10− | 5   | Explicit,      |     |     |     |
| 4   | 10− 5 |                |     |     |     |     | ×        |     |                |     |     |     |
| k × |       |                |     |     |     |     | k        |     |                |     |     |     |
| xe  |       |                |     |     |     |     | xe       |     |                |     |     |     |
)
| x 3   | 10− 5 |     |     |     |     |     | )     |     |     |     |     |     |
| ----- | ----- | --- | --- | --- | --- | --- | ----- | --- | --- | --- | --- | --- |
| uρ( × |       |     |     |     |     |     | uρ( y | 5   |     |     |     |     |
4 10−
×
| −   |     | slope=0.9 |     |     |     |     | −   |     | slope=0.8 |     |     |     |
| --- | --- | --------- | --- | --- | --- | --- | --- | --- | --------- | --- | --- | --- |
3 10− 5
| 2     | 10− 5 |     |     |     |     |     | ×     |     |     |     |     |     |
| ----- | ----- | --- | --- | --- | --- | --- | ----- | --- | --- | --- | --- | --- |
| h ×   |       |     |     |     |     |     | h     |     |     |     |     |     |
| )     |       |     |     |     |     |     | )     |     |     |     |     |     |
| uρ( x |       |     |     |     |     |     | uρ( y |     |     |     |     |     |
2 10− 5
| k   |     |     |     |     |     |     | k × |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
10− 5
|     |     |     | 1   | 14  | 1   | 1   |     |     |     | 1   | 14  | 1 1 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
10− 1 2 10− 3 1 0− 10− 6 10− 100 10− 1 2 10− 3 1 0− 10− 6 10− 100
|     |     |       | ×   | × h   | ×   | ×   |     |     |     | ×        | × h   | × × |
| --- | --- | ----- | --- | ----- | --- | --- | --- | --- | --- | -------- | ----- | --- |
|     |     | Error | L2  | on qx |     |     |     |     |     | Error L2 | on qy |     |
Figure 6.4: Convergence in mesh with a fixed Mach number M = 1e 4
b
−

CHAPTER 6. EXTENSION OF THE RAVIART-THOMAS STAGGERED SCHEME TO
| COMPRESSIBLE |     | BAROTROPIC        |     |     | FLOWS |     |     |                   |     |     | 177 |
| ------------ | --- | ----------------- | --- | --- | ----- | --- | --- | ----------------- | --- | --- | --- |
|              |     | Semi-Implicit,RT1 |     |     |       |     |     | Semi-Implicit,RT1 |     |     |     |
10−3
10−2
|     |      | Explicit,RT1 |     |     |     |     |      | Explicit,RT1 |     |     |     |
| --- | ---- | ------------ | --- | --- | --- | --- | ---- | ------------ | --- | --- | --- |
|     | 10−5 |              |     |     |     |     | 10−3 |              |     |     |     |
)vid(H
1H| 10−7
10−4
|
| h   |       |                    |     |     |     | h      |      |           |     |     |     |
| --- | ----- | ------------------ | --- | --- | --- | ------ | ---- | --------- | --- | --- | --- |
| ˜ρ  | 10−9  | ssllooppee==22..00 |     |     |     | )u˜˜ρ( |      | slope=1.0 |     |     |     |
| |   |       |                    |     |     |     |        | 10−5 |           |     |     |     |
|     | 10−11 |                    |     |     |     |        | |    |           |     |     |     |
10−6
10−13
10−7
10−7 10−6 10−5 10−4 10−3 10−2 10−1 10−7 10−6 10−5 10−4 10−3 10−2 10−1
|     |           |           | Mach | number      |        |                   |        |                     | Mach number    |        |       |
| --- | --------- | --------- | ---- | ----------- | ------ | ----------------- | ------ | ------------------- | -------------- | ------ | ----- |
|     | (a)       | ∇ρ˜ w.r.t | the  | Mach        | number |                   |        | divρ˜u˜             |                |        |       |
|     |           | L2        |      |             |        |                   | (b)    | L2                  | w.r.t the Mach | number |       |
|     | (cid:107) | (cid:107) |      |             |        |                   |        | (cid:107) (cid:107) |                |        |       |
|     |           | Figure    | 6.5: | Convergence |        | in Mach           | number | on a fixed          | mesh           |        |       |
|     |           |           |      |             |        | Initial condition |        |                     |                |        |       |
| t   | = 0.25 s  |           |      |             |        | t = 0.5           | s      |                     |                | t =    | 3.5 s |
Figure 6.6: Propagation of a low Mach acoustic wave through a stationary vortex: snapshots
| of the | Mach | number | at different |     | times. |     |     |     |     |     |     |
| ------ | ---- | ------ | ------------ | --- | ------ | --- | --- | --- | --- | --- | --- |

Chapter 7
Conclusion
This work comes within the scope of the research of numerical tool that is able to simulate low
Mach number flows in a conservative manner. In this aim, a proposition of staggered scheme
for Euler barotropic equations based on a de Rham complex is introduced. The methodology
unraveled in this work, while applied to the construction of a staggered discretization, is not
specific to the approximation spaces chosen here. Indeed, as far as the low Mach number
behaviour is concerned, this thesis, combined with the pre-existing literature, establishes a
clear understanding of the Euler barotropic system. Let us sum up the main takeaways:
In chapter 2, we recall that a formal link exists between the low Mach number limit on
Euler barotropic equations and a linear wave system. Indeed, a double time scale asymptotic
analysis in Mach number enables to demonstrate the formal equivalence between the low Mach
number limit and the long time behaviour of a first-order linear wave system. In this simplified
set up, the ability to identify the limit is entirely determined by the existence of a particular
decompositionofthevelocityfield: aHodge-Helmholtzdecompositionwhichincludesboundary
conditions on the divergence free part.
Then, it should be understood that the design of a class of low Mach number accurate
numerical scheme can be made by constructing a numerical scheme that is long time consistent
onthewavesystem. Inthisregard, itisshownthatsufficientconditionsforanumericalscheme
to be long time consistent are: existence of an underlying Hodge-Helmholtz decomposition,
preservation of continuous structures such as the stationary states of the system and finally
relative energy dissipation with a peculiar stationary state.
Energy dissipation is easily investigated in the particular case of the one dimension, so
that potential candidates of time steppings and numerical diffusion operators can be readily
identified in these conditions. In parallel, in one space dimension, all staggered schemes are
identical with regard to the space of approximation. Hence, in chapter 3 we develop a thorough
studyontheonedimensionalwavesystemofenergydissipationandoscillationspropertiesofthe
staggereddiscretizationwhencombinedwithdifferenttimeintegrationstrategiesandnumerical
diffusions. It is concluded that classical staggered schemes yield either a fully centered Implicit
scheme or an ImEx pressure centered scheme; the former is energy dissipative while the latter
179

CHAPTER 7. CONCLUSION 180
is not. Both schemes suffer from oscillations that can be damped by adding diffusion terms.
We introduce also fully explicit time integrations for which it is proven that energy dissipation
can only be achieved by adding diffusion on both the pressure and velocity equations.
In an abstract setting, the Hodge-Helmholtz decomposition can be seen as the byproduct
of an underlying de Rham complex; chapter 4 is dedicated to a gentle introduction of this link.
A canonical de Rham complex in two space dimensions and the natural consequences of this
formalism are presented. The staggered scheme proposed in this thesis is based on a discrete
equivalentofthisdeRhamcomplex: theN´ed´elec-Raviart-Thomasfiniteelementcomplex. From
there, weareabletoshowtheexistenceofadiscreteHodge-Helmholtzdecompositionincluding
boundary conditions by tweaking the original discrete complex in such way that it takes into
account the boundary.
More generally, this chapter shows that a discrete de Rham complex leads to a natural
definition of differential operators and most importantly, adjoints operators. Consequently, a
discrete Laplacian, called Hodge-Laplacian can be defined by composing and combining these
discrete differential operators. In fact, given the differential operators and scalar products
contained in the definition of a discrete de Rham complex, one can define by duality these
adjoint operators. The very existence of the complex and its link with a continuous one ensures
therelevanceofthesedefinitions. Also,wenoticethat,nomatterthediscretedeRhamcomplex
involved, one can readily extend these operators in order to solve the fitting Poisson problem
leadingtotheexistenceofadiscreteHodge-Helmholtzdecompositionwithboundaryconditions.
Using approximation spaces based on a discrete de Rham complex, one can derive a numer-
icalschemeforthemulti-dimensionalwavesystem. Inchapter5thisderivationisproposedina
Galerkinfashion,withastaggeredRaviart-Thomasvelocityandacellwiseconstantpressure. It
is complemented with the proposition of a new numerical diffusion operator obtained from the
differential operators given by the discrete de Rham complex. Indeed, introducing numerical
diffusion is, in general, desirable for stability reasons, that might be related either to the time
stepping or, in a non-linear setting, to the occurence of shocks. In this context, a numerical
diffusion in grad-divergence is defined thanks to the complex and its natural adjoint operators.
It is demonstrated that this diffusion operator preserves a discrete curl, or equivalently station-
ary states. Moreover, stability properties are established on the fully discrete scheme, notably
in explicit time stepping using this new grad-div diffusion operator.
Thus, on the one hand, the discretization space yields a Hodge-Helmholtz decomposition
with boundary conditions and on the other hand the numerical scheme induces energy dis-
sipation and preservation of stationary states thanks to the appropriate diffusion operator.
Combining these elements, we are able to show long time consistency on the wave system of
the proposed staggered discretization.
Besides, we have constructed, thanks to the aforementioned Hodge-Laplacian, a numerical
diffusion that does not preserve stationary state of the wave system. The loss of this property
implied, as shown in numerical simulations, the loss of convergence to the accurate long time
limit. This conclusion hints that, contrary to popular belief, using staggered schemes do not

CHAPTER 7. CONCLUSION 181
necessarily imply low Mach number precision.
Again, stepping back, we notice that given a discrete de Rham complex, one can define a
stationarity preserving and energy dissipating numerical diffusion with the composition of the
discrete wave operator and its adjoint.
Following a (partially discontinuous) Galerkin approach detailed in chapter 6, one can
extend the de Rham complex based staggered discretization from the wave system to the Euler
barotropic equations. Therefore, the density is approximated by a cellwise constant function
whereasthediscretemomentumisintheRaviart-Thomasfiniteelementspace. InthisGalerkin
framework, the resulting upwind discrete mass equation preserves the positivity of the density
while it is sufficient for the discrete momentum space to contain uniform vector fields to ensure
global momentum conservation. In parallel, the low Mach number precision is obtained by
making sure that the diffusion appearing on the acoustic scale respects the acoustic kernel; the
new grad-div staggered diffusion is the only one depending on the sound speed, so that the
doubletimescaleasymptoticanalysisinMachnumbershowsconsistencywiththediscretization
introduced on the first order wave system.
A few improvements for our proposition of staggered scheme should be aimed, starting
with a discrete entropy dissipation property. Indeed, discrete entropy type inequalities have
beenprovenforCrouzeix-Raviart/Rannacher-Turek(CR/RT)staggereddiscretizationsbutthe
arguments used are not available in the case of the proposed staggered scheme. However, the
grad-div numerical diffusion introduced in this thesis is readily defined for CR/RT staggered
discretizations, and this operator can be used to obtain global entropy inequalities in Euler
explicit time integration as shown in Appendix D.
Also while the discrete density is locally conserved, the obtained momentum conservation is
only global; it is very desirable to obtain a local conservation property for the latter in order to
simulateequationsarisingfromconservationprinciples. Athree-dimensionaladaptationispos-
sible by basing the discretization on the three-dimensional N´ed´elec-Raviart-Thomas complex;
discrete differential operators and a Hodge-Helmholtz decomposition can be redefined fittingly
following the methodology introduced in chapter 4, chapter 5 and chapter 6.
Finally, the organic sequel of this work is to extend this scheme to the full Euler system,
whichintroducesnewchallenges. Inparticular, becauseofthestaggering, classicalconvexcom-
bination arguments [122, 13, 123] are not available and thus, defining a scheme that conserves
total energy while preserving positivity of the internal energy is a subtle problem.
As a final note, one should have in mind that the natural formalism for the de Rham com-
plexes is the language of differential forms introduced by [124], and in particular the choice of
the Raviart-Thomas elements for the velocity/momentum staggering can be seen as approx-
imating 1 differential forms in 2D or 2 differential forms in 3D. With this in mind, we infer
− −
that the approach introduced in this thesis, notably adapting results obtained in chapter 4 and
chapter 5, might benefit other systems that require preservation of differential constraints, such
as Maxwell equations or magnetohydrodynamics models.

Appendix A
Proof of Energy Dissipation in ImEx
We use formulation (5.53) and thus we denote in what follows
1
∆pn := σ ε (σ)[[pn]] .
K K | | K σ
| | σ ∂K
(cid:88)∈
(cid:103)
Step 1: we make appear an equation on the global energy
Following Step1 of the proof for the Explicit scheme (5.48), we multiply the pressure
equation of (5.53) ρ κ pn
0 0 K
pn+1 pn 1
K K − Kρ κ pn + σ ε (σ)unρ κ pn = dc σ ε (σ)[[pn]] ρ κ pn , (A.1)
| | δτ 0 0 K ρ | | K σ 0 0 K 0 | | K σ 0 0 K
0
σ ∂K σ ∂K
(cid:88)∈ (cid:88)∈
and the equation on the velocity of (5.53) by un
σ
un+1 un
D σ − σun+κ σ [[pn]] un = δτc2 σ [[div(un)]] un δτdc κ σ [[∆pn]] un, (A.2)
| σ | δτ σ 0 | | σ σ 0| | σ σ − 0 0 | | σ σ
Using (5.56) in (A.1) yields: (cid:103)
K (pn+1)2 (pn )2 K 1
| |ρ κ K K | |ρ κ (pn+1 pn )2+ σ ε (σ)unρ κ pn =
δτ 0 0 2 − 2 − 2δτ 0 0 K − K ρ | | K σ 0 0 K
(cid:18) (cid:19) 0 σ ∂K
dc σ ε (σ)[[pn]] ρ κ pn
(cid:88)∈
0 K σ 0 0 K
| |
σ ∂K
(cid:88)∈
(A.3)
183

APPENDIX A. PROOF OF ENERGY DISSIPATION IN IMEX 184
, and (5.56) in (A.2) yields:
D (un+1)2 (un)2 D
σ σ σ | σ |(un+1 un)2+κ σ [[pn]] un
δτ 2 − 2 − 2δτ σ − σ 0 | | σ σ
(cid:18) (cid:19) (A.4)
= δτc2 σ [[div(un)]] un δτdc κ σ [[∆pn]] un.
0| | σ σ − 0 0 | | σ σ
Summing (A.3) on the cells K gives (cid:103)
∈ C
K (pn+1)2 (pn )2
| |ρ κ K K +κ σ ε (σ)unpn
δτ 0 0 2 − 2 0 | | K σ K
K (cid:18) (cid:19) K σ ∂K
(cid:88)∈C (cid:88)∈C (cid:88)∈
(A.5)
K
= dc σ [[pn]] pn ρ κ + | |ρ κ (pn+1 pn )2 ,
0 | | σ K 0 0 2δτ 0 0 K − K
(cid:32) (cid:33)
K σ ∂K
(cid:88)∈C (cid:88)∈
and similarly, summing (A.4) on the face σ int
∈ F
D (un+1)2 (un)2
| σ | σ σ +κ σ [[pn]] un
δτ 2 − 2 0 | | σ σ
σ int (cid:18) (cid:19) σ int
∈(cid:88)F ∈(cid:88)F
(A.6)
D
= δτc2 σ [[div(un)]] un δτdc κ σ [[∆pn]] un+ | σ |(un+1 un)2 .
0| | σ σ − 0 0 | | σ σ 2δτ σ − σ
(cid:32) (cid:33)
σ int
∈(cid:88)F
(cid:103)
By summing (A.5) and (A.6) and using point i) of Lemma 5.3.1:
E(Un+1) E(Un)
h δτ − h +κ 0 | σ | pn Kσ un σ = L, (A.7)
σ b
(cid:88)∈F
. where
K
L = dc σ ε (σ)[[pn]] pn ρ κ + | |ρ κ (pn+1 pn )2
0 | | K σ K 0 0 2δτ 0 0 K − K
(cid:32) (cid:33)
K σ ∂K
(cid:88)∈C (cid:88)∈
(A.8)
D
+ δτc2 σ [[div(un)]] un δτdc κ σ [[∆pn]] un+ | σ |(un+1 un)2 .
0| | σ σ − 0 0 | | σ σ 2δτ σ − σ
(cid:32) (cid:33)
σ int
∈(cid:88)F
(cid:103)
Step 2: we make appear appear the diffusive terms by integration by parts Now

APPENDIX A. PROOF OF ENERGY DISSIPATION IN IMEX 185
by integration by parts with pointii) of Lemma 5.3.1 we have again
dc ρ κ
dc σ ε (σ)[[pn]] pn ρ κ = dc ρ κ σ [[pn]]2 + 0 0 0 σ [[(pn)2]] ,
0 K σ K 0 0 0 0 0 σ σ
| | − | | 2 | |
K σ ∂K σ σ b
(cid:88)∈C (cid:88)∈ (cid:88)∈F (cid:88)∈F
(A.9)
and iii) of Lemma 5.3.1
δτc2 σ [[div(un)]] un = δτc2 K div(un)2+δτc2 σ (divun) un, (A.10)
0 | | σ σ − 0 | | 0 | | Kσ σ
σ int K σ b
∈(cid:88)F (cid:88)∈C (cid:88)∈F
by applying iii) of Lemma 5.3.1 with p replaced by ∆pn
K K
δτdc κ σ [[∆pn]] un = δτdc κ K (∆pn) d(cid:103)iv(un) δτdc κ σ (∆pn) un.
−
0 0
| |
σ σ 0 0
| |
K K
−
0 0
| |
Kσ σ
σ int K σ b
∈(cid:88)F (cid:88)∈C (cid:88)∈F
(cid:103) (cid:103) (cid:103)
(A.11)
Gathering (A.9), (A.10) and (A.11) in (A.8) we obtain
dc ρ κ K
L= dc ρ κ σ [[pn]]2 + 0 0 0 σ [[(pn)2]] + | |ρ κ (pn+1 pn )2
− 0 0 0 | | σ 2 | | σ 2δτ 0 0 K − K
σ σ b K
(cid:88)∈F (cid:88)∈F (cid:88)∈C
δτc2 K div(un)2+δτc2 σ (divun) un
− 0 | | 0 | | Kσ σ
K σ b
(cid:88)∈C (cid:88)∈F
(A.12)
+δτdc κ K (∆pn) div(un) δτdc κ σ (∆pn) un
0 0
| |
K K
−
0 0
| |
Kσ σ
K σ b
(cid:88)∈C (cid:88)∈F
(cid:103) (cid:103)
D
+ | σ |(un+1 un)2
2δτ σ − σ
σ int
∈(cid:88)F
Using (A.12) in (A.7) yields
E(Un+1) E(Un)
h − h + σ Φn = R, (A.13)
δτ | | σ
σ ∂Ω
(cid:88)∈
with the boundary flux:
dc ρ κ
Φn := κ pn un δτc2(divun) un 0 0 0 [[(pn)2]] +δτdc κ (∆pn) un,
σ 0 Kσ σ − 0 Kσ σ − 2 σ 0 0 Kσ σ
(cid:103)

APPENDIX A. PROOF OF ENERGY DISSIPATION IN IMEX 186
with now R equal to
R = c δτc 0 | ∂K | ∂K d (cid:94) iv(un) 2 dc ρ κ σ [[pn]]2
− K | | − 0 0 0 | | σ
K | | σ
(cid:88)∈C (cid:88)∈F
diffusiveterms
(cid:124) (cid:123)(cid:122) (cid:125)
K D
+ | |ρ κ (pn+1 pn )2+ | σ |(un+1 un)2+δτdc ρ κ K (∆pn) div(un) .
2δτ 0 0 K − K 2δτ σ − σ 0 0 0 | | K K
K σ int K
(cid:88)∈C ∈(cid:88)F (cid:88)∈C
(cid:103)
non-negativeterms unsignedterms
(cid:124) (cid:123)(cid:122) (cid:125) (cid:124) (cid:123)(cid:122) (cid:125)
(A.14)
Step 3: we bound the non-negative and unsigned terms in the rest R as a factors
of the diffusive terms
First we notice from the pressure equation of (5.53) that:
K K δτ 1 2
| |ρ κ (pn+1 pn )2= | |ρ κ K div(un) +dc K ∆pn
2δτ 0 0 K − K 2δτ 0 0 K − ρ | | K 0 | | K
K K (cid:20)| |(cid:18) 0 (cid:19)(cid:21)
(cid:88)∈C (cid:88)∈C
(cid:103)
K 1 2
= | |ρ κ δτ2 div(un) +dc ∆pn
2δτ 0 0 − ρ K 0 K
K (cid:20) 0 (cid:21)
(cid:88)∈C
(cid:103)
= | K |ρ κ δτ2 1 div(un)2 2dc 0 div(un) ∆pn +d2c2∆pn 2 .
2δτ 0 0 ρ2 K − ρ K K 0 K
K (cid:20) 0 0 (cid:21)
(cid:88)∈C
(cid:103) (cid:103)
(A.15)

APPENDIX A. PROOF OF ENERGY DISSIPATION IN IMEX 187
This gives
K
| |ρ κ (pn+1 pn )2
2δτ 0 0 K − K
K
= c2 0dδτ K div(un)2 δτd (cid:88) c ∈C κ K div(un) ∆pn + | K |ρ κ δτd2c2∆pn 2 .
2 | | K − 0 0 | | K K 2 0 0 0 K
K K K
(cid:88)∈C (cid:88)∈C (cid:88)∈C
(cid:103) (cid:103)
(A.16)
Using the expression (A.16) we have
K
| |ρ κ (pn+1 pn )2+δτdc κ K ∆pn div(un)
2δτ 0 0 K − K 0 0 | | K K
K K
(cid:88)∈C (cid:88)∈C
(cid:103) (A.17)
= c ρ κ d δτc 0 ( σ [[pn]] )2+ c 0 δτc 0 | ∂K | ∂K d (cid:94) ivun 2 ,
0 0 0 2 K | | σ 2 K | | K
K | | σ K | |
(cid:88)∈C (cid:88)∈F (cid:88)∈C
meaning with an Inverse Poincar´e inequality i) from Lemma 5.3.2 on the pressure Laplacian in
(A.17)
K
| |ρ κ (pn+1 pn )2+δτdc ρ κ K (∆pn) div(un)
2δτ 0 0 K − K 0 0 0 | | K K
K K
(cid:88)∈C (cid:88)∈C
(cid:103) (A.18)
c ρ κ δτc 0 | σ | max dν σ [[pn]]2 + c 0 δτc 0 | ∂K | ∂K d (cid:94) ivun 2 .
≤ 0 0 0 h max | | σ 2 K | | K
σ K | |
(cid:88)∈F (cid:88)∈C
In parallel the velocity equation from (5.53) gives
2
D δτ
| σ |(un+1 un)2 = κ σ [[pn]] +δτc2 σ [[div(un)]] δτdc κ σ [[∆pn]] .
2δτ σ − σ 2 D − 0 | | σ 0| | σ − 0 0 | | σ
σ (cid:32) (cid:33)
σ int σ int | |
∈(cid:88)F ∈(cid:88)F
(cid:103)
(A.19)
Now a discrete Jensen inequality yields for three real numbers a,b,c,
(a+b+c)2 3(a2+b2+c2). (A.20)
≤ 0

APPENDIX A. PROOF OF ENERGY DISSIPATION IN IMEX 188
Applying (A.20) in (A.19) implies
D 3δτ
| σ |(un+1 un)2 κ2( σ [[pn]] )2+δτ2c4 σ 2[[div(un)]]2 +δτ2d2c2κ2 σ 2[[∆pn]]2 ,
2δτ σ − σ ≤ 2 D 0 | | σ | | σ 0 0| | σ
σ (cid:32) (cid:33)
σ int σ int | |
∈(cid:88)F ∈(cid:88)F
(cid:103)
(A.21)
then, first with h = min(min D ,min K )
σ σ K
| | | |
3δτ 3δτc σ
κ2( σ [[pn]] )2 ρ cκ 0 | | max σ [[pn]]2. (A.22)
2 D 0 | | σ ≤ 0 0 2h | | σ
σ
σ int | | σ int
∈(cid:88)F ∈(cid:88)F
By Inverse Poincar´e inequalities ii) of Lemma 5.3.2 on the grad-div velocity
3δτ 3δτ
δτ2c4 σ 2[[div(un)]]2 δτ2c4 ν ∂K 2div(un)2 , (A.23)
2 D | | σ ≤ h max | | K
σ
σ int | | K
∈(cid:88)F (cid:88)
∂K 2 (cid:94) 2
using div(un)2 = | | div(un) (A.23) becomes
K K 2 K
| |
δτ2c4 3δτ σ 2[[div(un)]]2 δτc 0 | ∂K |3ν c δτc 0 | ∂K | 2 ∂K d (cid:94) iv(un) 2 . (A.24)
2 D | | σ ≤ h max K | | K
σ
σ int | | K (cid:18) | | (cid:19)
∈(cid:88)F (cid:88)
(cid:94)
Finally, using the inverse Poincar´e inequality ii) with ∆pn instead of div(un) we obtain
K K
3δτ δτ2d2c2κ2 σ 2[[∆pn]]2 3δτ δτ2d2c2(cid:103) κ2(2ν ) ∂K 2∆pn 2 , (A.25)
2 D 0 0| | σ ≤ 2h 0 0 max | | K
σ
σ int | | K
∈(cid:88)F (cid:88)
(cid:103) (cid:103)
but
2
∂K 2∆pn 2 := ∂K 2 1 σ ε (σ)[[pn]]
| | K | | K 2 | | K σ
K K (cid:20)| | σ ∂K (cid:21)
(cid:88) (cid:88) (cid:88)∈
(cid:103)
∂K 2 2
| |max σ ε K (σ)[[pn]] σ (A.26)
≤ h2 | |
K (cid:20)σ ∂K (cid:21)
(cid:88) (cid:88)∈
∂K 2
| |maxν σ 2[[pn]]2,
≤ h2 max | | σ
K σ ∂K
(cid:88) (cid:88)∈
2
where we used again discrete Jensen inequality σ ε (σ)[[pn]]
K σ
| | ≤
(cid:20)σ ∂K (cid:21)
(cid:88)∈
ν σ 2[[pn]]2 for the last inequality. Then, using in (A.26) the first inverse
max σ
| |
K σ ∂K
(cid:88) (cid:88)∈

APPENDIX A. PROOF OF ENERGY DISSIPATION IN IMEX 189
Poincar´e inequality i) of Lemma 5.3.2 we obtain
∂K 2∆pn 2 (| ∂K | 2 max2ν σ ) σ [[pn]]2. (A.27)
| | K ≤ h2 max | | max | | σ
K σ
(cid:88) (cid:88)∈F
(cid:103)
Plugging (A.27) in (A.25) we have
3δτ 3δτ ∂K 2
δτ2d2c2κ2 σ 2[[∆pn]]2 δτ2d2c2κ2(2ν )(| |max2ν σ ) σ [[pn]]2.
2 D 0 0| | σ ≤ 2h 0 0 max h2 max | | max | | σ
σ
σ int | | σ
∈(cid:88)F (cid:88)∈F
(cid:103)
(A.28)
Gathering (A.28) (A.24) (A.22) in (A.21) we have
D 3δτc σ
| σ |(un+1 un)2 ρ cκ 0 | | max σ [[pn]]2
2δτ σ − σ ≤ 0 0 2h | | σ
σ int σ int
∈(cid:88)F ∈(cid:88)F
+ δτc 0 | ∂ h K | max 3ν max c δτc 0 K | ∂K | 2 | ∂K | d (cid:94) iv(un) 2 K (A.29)
K (cid:18) | | (cid:19)
(cid:88)∈C
δτ σ δτc ∂K
+ | | max 6( | | max dκ ν )2 σ [[pn]]2.
h h 0 max | | σ
σ int
∈(cid:88)F
Gathering (A.29) and (A.18) in (A.14), we thus obtain that the rest R is bounded by:
δτc σ 13 δτc ∂K
R dc ρ κ 0 | | max ν + +6d( | | max ν )2 1 σ [[pn]]2
≤ 0 0 0 h max d2 h max − | | σ
(cid:32) (cid:33)
(cid:18) (cid:19) σ
(cid:88)∈F
+c δτc 0 | ∂K | 1 +3ν δτc 0 | ∂K | δτc 0 | ∂K | max 1 ∂K d (cid:94) iv(un) 2 ,
K 2 max K h − | | K
(cid:32) (cid:33)
K | | (cid:18) | | (cid:19)
(cid:88)∈C
we find that the velocity term is non-positive under the condition
δτc ∂K 1
0 | | max ,
h ≤ √6ν
max
while under this CFL:
δτc σ 3 δτc ∂K δτc σ 3
0 | | max ν + +6d( | | max ν )2 1 0 | | max (ν (1+d)+ ) 1,
max max max
h 2d h − ≤ h 2d −
(cid:18) (cid:19)

| APPENDIX     | A.   | PROOF           | OF  | ENERGY | DISSIPATION |      | IN  | IMEX | 190 |
| ------------ | ---- | --------------- | --- | ------ | ----------- | ---- | --- | ---- | --- |
| the pressure | term | is non-positive |     | under  | the         | CFL: |     |      |     |
|              |      |                 |     | δτc    | σ           |      | 1   |      |     |
0 max
|     |     |     |     |     | | | |          |     | .   | (A.30) |
| --- | --- | --- | --- | --- | --- | -------- | --- | --- | ------ |
|     |     |     |     |     | h   | ≤        |     | 3   |        |
|     |     |     |     |     |     | ν (1+d)+ |     |     |        |
max
2d
| which concludes |     | the proof.        |     |     |     |        |                      |     |     |
| --------------- | --- | ----------------- | --- | --- | --- | ------ | -------------------- | --- | --- |
|                 |     |                   |     |     |     |        | 3                    |     | 3   |
| RemarkA.0.1.    |     | Noticenowthatg(x) |     |     | ν   | (1+x)+ | hasaglobalminimuminx |     |     |
|                 |     |                   |     |     | =   | max    |                      | =   |     |
|                 |     |                   |     |     |     |        | 2x                   |     | 2ν  |
(cid:114) max
1
| so that | the maximum |     | of the cfl |     | is  |     | .   | Obviously |     |
| ------- | ----------- | --- | ---------- | --- | --- | --- | --- | --------- | --- |
(A.30)
|     |     |     |     |     |     | ν +√6ν |     |     |     |
| --- | --- | --- | --- | --- | --- | ------ | --- | --- | --- |
|     |     |     |     |     |     | max    | max |     |     |
|     |     |     |     |     | 1   |        | 1   |     |     |
,
|          |      |            |        | ν     | +√6ν   |       | √6ν |     |     |
| -------- | ---- | ---------- | ------ | ----- | ------ | ----- | --- | --- | --- |
|          |      |            |        | max   |        | max ≤ | max |     |     |
| thus the | more | penalizing | cfl is |       |        |       |     |     |     |
|          |      |            |        | δτc 0 | ∂K max |       | 1   |     |     |
|          |      |            |        |       | | |    |       |     | .   |     |
h
|     |     |     |     |     |     | ≤ ν max | +√6ν | max |     |
| --- | --- | --- | --- | --- | --- | ------- | ---- | --- | --- |

| Appendix    |     | B   |        |        |         |
| ----------- | --- | --- | ------ | ------ | ------- |
| Computation |     |     | of the | lumped | rotated |
gradient
| We detail | here briefly | the formula |        |        |     |
| --------- | ------------ | ----------- | ------ | ------ | --- |
|           |              |             | (∇ ϕ ) | n (x ) |     |
|           |              |             | ⊥ n    | σ σ    |     |
·
| for ϕ a basis | function | of cG1(Ω) associated | to a | node n : |     |
| ------------- | -------- | -------------------- | ---- | -------- | --- |
n
| ∇ ϕ on | triangles |     |     |     |     |
| ------ | --------- | --- | --- | --- | --- |
⊥ n
v
σ
n
w
Given a cell K and a node n ∂K, we have, in triangles, that the basis function associated
∈
with the node n is equal to ϕ n (x) = λ n (x) where λ n is the barycentric coordinate associated
K
|               |              | |      | det(x   | v,x w)    |     |
| ------------- | ------------ | ------ | ------- | --------- | --- |
| with the node | n. Recalling | that λ | (x) = − | − we have |     |
|               |              | n      | 2       | K         |     |
| |
1
|     |     |     | ∂ ϕ = | (v w)x1, |     |
| --- | --- | --- | ----- | -------- | --- |
x2 n
|     |     | −   | 2 K | −   |     |
| --- | --- | --- | --- | --- | --- |
|     |     |     | |   | |   |     |
1
w)x2,
|     |     |     | ∂ x1 ϕ n = | (v  |     |
| --- | --- | --- | ---------- | --- | --- |
|     |     |     | 2 K        | −   |     |
|     |     |     | | |        |     |     |
192

| APPENDIX | B.  | COMPUTATION |     |     | OF  | THE | LUMPED |     | ROTATED |     | GRADIENT |     | 193 |
| -------- | --- | ----------- | --- | --- | --- | --- | ------ | --- | ------- | --- | -------- | --- | --- |
v w
| leading | to ∇ | ϕ = | −   | . Meaning | that |     |     |     |     |     |     |     |     |
| ------- | ---- | --- | --- | --------- | ---- | --- | --- | --- | --- | --- | --- | --- | --- |
|         | ⊥    | n   | K   |           |      |     |     |     |     |     |     |     |     |
| |
|     |     |     |       |      | v   | w   |      |     |          |     | 2   |     |     |
| --- | --- | --- | ----- | ---- | --- | --- | ---- | --- | -------- | --- | --- | --- | --- |
|     |     | (∇  | ϕ ) n | (x ) | = − | n   | = (v | w)  | n        |     |     | ,   |     |
|     |     | ⊥   | n     | σ σ  |     |     | σ    |     | σ        |     |     |     |     |
|     |     |     | ·     |      | 2 K | ·   |      | −   | · 2det(n |     | v,w | v)  |     |
|     |     |     |       |      | |   | |   |      |     |          | −   |     | −   |     |
notice now that det(n v,w v) = det(v w,n v) = σ det(v w,n )(n (n v)). But
|     |     |     | −   | −   |     | −   |     | −   | | | | −   | ⊥σ  | ⊥σ · − |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | ------ | --- |
det(u,v ⊥ ) = u v meaning det(n v,w v) = σ (v w) n σ sign(n ⊥σ (n v)) so
|     |     | ·     |      |     | −   | −   | |     | | − | ·    |      | ·   | −      |     |
| --- | --- | ----- | ---- | --- | --- | --- | ----- | --- | ---- | ---- | --- | ------ | --- |
|     |     |       |      | v   | w   |     |       |     |      |      | 1   |        |     |
|     | (∇  | ϕ     | n (x |     |     | n   | (v w) | n   |      |      |     |        |     |
|     |     | ⊥ n ) | σ σ  | ) = | −   | σ = |       |     | σ    |      |     |        |     |
|     |     | ·     |      |     | K · |     | −     | ·   | σ (v | w) n | (n  | (n v)) |     |
σ ⊥σ
|     |     |     |          |     | | | |     |          |     | | | − | ·    |     | · − |     |
| --- | --- | --- | -------- | --- | --- | --- | -------- | --- | ----- | ---- | --- | --- | --- |
|     |     |     |          |     |     |     | 1        |     |       |      | 1   |     |     |
|     |     |     | = sign(n |     | (n  | v)) | = sign(n |     | (x    | x )) | .   |     |     |
|     |     |     |          | ⊥σ  |     |     |          |     | ⊥σ n  | σ    |     |     |     |
|     |     |     |          |     | · − |     | σ        |     | ·     | −    | σ   |     |     |
|     |     |     |          |     |     | |   | |        |     |       |      | | | |     |     |
∇ ⊥ ϕ n on quadrangles In quadrangles the basis function associated with the node n in a
cell K for which n ∂K is expressed thanks to the transformation from the reference element
⊂
| to the | element | K.  |     |     |     |     |     |     |     |     |     |     |     |
| ------ | ------- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
w,(1,1)
M
2
v,(0,1)
M
1
σ
M
3
u,(1,0)
M
0
n,(0,0)
| Following | the | figure | we define |       |       |           |     |               |     |     |        |     |     |
| --------- | --- | ------ | --------- | ----- | ----- | --------- | --- | ------------- | --- | --- | ------ | --- | --- |
|           |     | T      | (xˆ,yˆ)   | := (1 | xˆ)(1 | yˆ)n+xˆ(1 |     | yˆ)u+xˆyˆw+(1 |     |     | xˆ)yˆv |     |     |
K
|     |     |     |     |     | −   | −   |     | −   |     |     | −   |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
nisthusassociatedwiththebasisfunctionof(0,0)inthereferenceelement,meaningϕˆ (xˆ,yˆ)
|     |     |     |     |     |     |     |     |     |     |     |     | 0   | =   |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
(1 xˆ)(1 yˆ). Now ∇ ϕ (x) := ∇ ϕˆ (T 1(x)) = (∇ T 1(x)) t∇ ϕˆ (T 1(x)). Since T is
|     |     |     | x   | n   | x   | 0 K− |     |     | x − | −   | xˆ 0 | −   |     |
| --- | --- | --- | --- | --- | --- | ---- | --- | --- | --- | --- | ---- | --- | --- |
| −   | −   |     |     |     |     |      |     |     |     |     |      |     |     |
invertible, ∇ T 1(x) = (∇ T(T 1(x))) 1. From the definition of the transformation
|     | x   | −   |     | xˆ  | −   | −   |      |      |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | ---- | ---- | --- | --- | --- | --- | --- |
|     |     |     |     |     | ∇   | T = | µ(y) | ν(x) | ,   |     |     |     |     |
xˆ
|     |     |     |     |     |     |     | (cid:0) |     | (cid:1) |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | ------- | --- | ------- | --- | --- | --- | --- |

APPENDIX B. COMPUTATION OF THE LUMPED ROTATED GRADIENT 194
with µ(y) = (1 yˆ)M +yˆM and ν(x) = (1 xˆ)M +xˆM . This gives
0 2 1 2
− −
1
(∇ T) t = ν(x) ,µ(y) ,
xˆ − det(µ(y),ν(x)) ⊥ ⊥
(cid:0) (cid:1)
so
∇ ϕ (x) = (∇ T 1(x)) t∇ ϕˆ (T 1(x))
x n x − − xˆ 0 −
1 1 yˆ 1
= ν(x) ,µ(y) − = (1 yˆ)ν(x) +(1 xˆ)µ(y) .
det(µ(y),ν(x)) ⊥ ⊥ 1 xˆ det(µ(y),ν(x)) − ⊥ − ⊥
(cid:18) − (cid:19) (cid:18) (cid:19)
(cid:0) (cid:1)
Leading to
1
∇ ϕ (x) = (yˆ 1)ν(x)+(xˆ 1)µ(y) , (B.1)
⊥x n
det(µ(y),ν(x)) − −
(cid:18) (cid:19)
since for a vector of R2 (q ) = q. Finally
⊥ ⊥
−
1 M M +M
∇ ϕ (x ) n = 3 0 2 n ,
⊥ n σ · σ
det(
M 0 +M 2,M
)(cid:18)
− 2 − 2
(cid:19)
· σ (B.2)
3
2
changing the order of the vectors in the determinant we multiply by minus one. Now we notice
that :
M +M M +M
det( 0 2,M ) = det(M + 0 2,M )
3 3 3
2 2
M +M M +M
= det(M + 0 2,n ) σ sign(M n ) = (M + 0 2 ) n σ sign(M n ),
3 ⊥σ 3 ⊥σ 3 σ 3 ⊥σ
2 | | · 2 · | | ·
(B.3)
Plugging (B.3) in (B.1) yields as in the triangular case
1
∇ ϕ (x ) n = sign(n (x x ))
⊥ n σ
·
σ ⊥σ
·
n
−
σ
σ
| |

|     | Appendix |      | C   |        |     |           |            |     |         |     |     |     |
| --- | -------- | ---- | --- | ------ | --- | --------- | ---------- | --- | ------- | --- | --- | --- |
|     | Low      | Mach |     | number |     |           | asymptotic |     |         |     |     |     |
|     | analysis |      | of  | some   |     | staggered |            |     | schemes |     |     |     |
We show two examples of low Mach number asymptotic analysis on two staggered schemes.
For both schemes, h K /ρ K is the height of the water/density at the cell K and plays a similar
role to the density in the Euler barotropic equations, D σ h σ := D K,σ h K + D L,σ h L , where
|     |     |     |     |     |     |     |     | | | | |   | |   | | | |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
K,L are the cells touching the face σ as in Figure C.1. u is the full velocity vector at the face
σ
σ, D is the dual cell associated to this face. Then e ∂D is a face on the boundary of D .
|     | σ   |     |     |     |     |     |     | σ   |     |     |     | σ   |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
⊂
For a fixed cell K and a fixed face σ ∂K, L is the only cell such that L K = σ so that n
K,σ
|     |     |     |     |     | ⊂   |     |     |     |     | ∩   |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
is oriented towards L. Similarly, for a fixed face σ and a fixed dual face e ∂D σ the only
|     |     |     |     |     |     |     |     |     |     |     | σ e(cid:48) |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | ----------- | --- |
⊂
face such that D D σ(cid:48) = e. For any cell of the mesh C dual or primal and n ∂C any dual
|     |     | σ   | ∩ e |     |     |     |     |     |     |     | ⊂   |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
of primal face n C,n denotes the unit outward normal at the face n and cell C.
C.1 First example: the explicit scheme of Duran .A, Vila. J-P
|     |     | and Baraille |     | .   | R   |     |     |     |     |     |     |     |
| --- | --- | ------------ | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
The explicit scheme of [69] is given on the shallow waters equations by
|     | hn+1 | hn    |      |     |     |     |     |     |     |     |     |     |
| --- | ---- | ----- | ---- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
|     | K K  | − K + | σ    | n n | = 0 |     |     |     |     |     |     |     |
|     |      |       |      | σ   | K,σ |     |     |     |     |     |     |     |
|     | | |  | ∆t    | | |F | ·   |     |     |     |     |     |     |     |     |
|    |      | σ     | ∂K   |     |     |     |     |     |     |     |     |     |
(cid:88)⊂
 

 
|    | hn+1un+1 | hnun |     |                       |        |         |     |                 |          |     |          |     |
| --- | -------- | ---- | --- | --------------------- | ------ | ------- | --- | --------------- | -------- | --- | -------- | --- |
|    | D σ      | σ σ  | σ + |                       | un ( n | n )++un |     | ( n n           | )− +∆t   | D   | hn ∇Φn , | = 0 |
|     | σ        | −    |     |                       | σ e    | Dσ,e    | σ   | (cid:48) e Dσ,e |          | σ   | σ σ ∗    |     |
|     | | |      | ∆t   |     |                       | F ·    |         |     | e F ·           |          | |   | |        |     |
|   |          |      | e   | ⊂(cid:88) ∂Dσ(cid:16) |        |         |     |                 | (cid:17) |     |          |     |

 


(C.1)
196

| APPENDIX  | C. LOW MACH | NUMBER | ASYMPTOTIC | ANALYSIS | OF SOME |     |
| --------- | ----------- | ------ | ---------- | -------- | ------- | --- |
| STAGGERED | SCHEMES     |        |            |          |         | 197 |
K
D
K,σ
D
|     |     |     | σ L,σ |     |     |     |
| --- | --- | --- | ----- | --- | --- | --- |
n L
K,σ
e
D
σ(cid:48) e
FigureC.1: RepresentationofadualcellassociatedtoafaceσforCrouzeix-Raviart/Rannacher-
| Turek velocity | staggering |     |     |     |     |     |
| -------------- | ---------- | --- | --- | --- | --- | --- |

APPENDIX C. LOW MACH NUMBER ASYMPTOTIC ANALYSIS OF SOME
STAGGERED SCHEMES 198
Then,
σ
n n
K,σ
= hn
K
(un
σ
n
K,σ
)++hn
L
(un
σ
n
K,σ
)− Πn
K,σ
n
K,σ
F · · · − ·
with (cid:0) (cid:1)
σ
Πn := γ∆t | | hn(Φn Φn )n
K,σ D σ L − K K,σ
σ
| |
where γ is a coefficient and
σ
∇Φn σ , ∗ := D | | Φ n L , ∗ − Φ n K , ∗ n K,σ
σ
| |
(cid:0) (cid:1)
where
∂K
Φ n K , ∗ := Φn K − 2αg∆t| K | ((hu)n σ − (hu)n K ) · n K,σ
| |
so that, we in fact have the much clearer definition
| D σ | ∇Φ n σ , ∗ = | σ | (Φn L− Φn K )n K,σ
∂L ∂K
σ 2αg∆t | | ((hu)n (hu)n) n | | ((hu)n (hu)n ) n n
−| | L σ − L · L,σ − K σ − K · K,σ K,σ
(cid:34) (cid:35)
| | | |
and it is sufficient, in order to obtain entropy dissipation, to define for any approximation of
the discharge in a cell K (hu)
K
(hu) := λ (hu) (C.2)
K K K
with
2
0 if σ (hu) n = 0
| | K · K,σ
 σ (cid:88)⊂ ∂K (cid:104) (cid:105)


  1
λ K :=   
   |
σ
|
(hu)
σ ·
n
K,σ
2 2 (C.3)
 σ (cid:88)⊂ ∂K (cid:104) (cid:105) 
2

     
 
  

σ (cid:88)⊂ ∂K |
σ
| (cid:104)
(hu)
K ·
n
K,σ (cid:105)   

The flux
e
is defined such th at the following discrete conservation law stands
F
D
| σ | hn+1 hn + n n = 0
∆t σ − σ F e · Dσ,e
(cid:0) (cid:1)
e
⊂(cid:88)
∂Dσ
Thanks to this dual mass balance equation, manipulating the mass and discharge equations,
we recover the discrete evolution of the velocity
| D ∆ σ t | un σ +1 − un σ + un σ h e (cid:48) n − +1 un σ ( F e n · n Dσ,e )−+∆t | D σ |h h n+ n σ 1 ∇Φn σ , ∗ = 0
(cid:0) (cid:1)
e
⊂(cid:88)
∂Dσ (cid:18) σ (cid:19) σ
From that we have :

| APPENDIX  |     | C. LOW  | MACH |     | NUMBER |     | ASYMPTOTIC |     | ANALYSIS |     | OF  | SOME |     |
| --------- | --- | ------- | ---- | --- | ------ | --- | ---------- | --- | -------- | --- | --- | ---- | --- |
| STAGGERED |     | SCHEMES |      |     |        |     |            |     |          |     |     |      | 199 |
LemmaC.1.1(Dimensionlessheight( density)equation). Thedimensionlessheightequation
≈
| of (C.1) | is given | by    |     |     |           |     |          |             |          |     |           |          |       |
| -------- | -------- | ----- | --- | --- | --------- | --- | -------- | ----------- | -------- | --- | --------- | -------- | ----- |
|          |          | h˜n+1 |     | h˜n |           |     |          |             |          |     |           |          |       |
|          |          | 1     |     |     | 1         |     | h˜n      |             | )++h˜n   |     |           |          |       |
|          |          |       | K − | K + |           | σ˜  |          | (un         | n        |     | (un n     | )−       |       |
|          |          | M     | τ   |     | K˜        |     | K        | σ           | K,σ      | L   | σ         | K,σ      |       |
|          |          |       | ∆   |     |           | |   | |        | ·           |          |     | ·         |          |       |
|          |          |       |     |     | | | σ     | ∂K  |          |             |          |     |           |          |       |
|          |          |       |     |     | (cid:88)⊂ |     | (cid:16) |             |          |     |           | (cid:17) |       |
|          |          |       |     |     | γ∆τ       | 1   |          | (cid:101)σ˜ |          |     |           |          |       |
|          |          |       |     |     |           |     |          |             | h˜n h˜n  | h˜n | (cid:101) |          |       |
|          |          |       |     |     |           |     | σ˜       | | |         |          |     | = 0.      |          | (C.4) |
|          |          |       |     |     | M         | K˜  |          | D˜          | σ L      | K   |           |          |       |
|          |          |       |     | −   |           |     | |        | | σ         |          | −   |           |          |       |
|          |          |       |     |     | |         | | σ | ∂K       | | |         | (cid:16) |     | (cid:17)  |          |       |
(cid:88)⊂
Proof. Using h , (cid:96) and u as scaling variables and assuming that Φ = gh, we get
|         |            | 0 0    |     | 0           |               |                  |          |             |            |        |           |          |     |
| ------- | ---------- | ------ | --- | ----------- | ------------- | ---------------- | -------- | ----------- | ---------- | ------ | --------- | -------- | --- |
|         |            | h˜n+1  | h˜n |             |               |                  |          |             |            |        |           |          |     |
|         | h          |        |     | h u         | (cid:96) 1    |                  |          |             |            |        |           |          |     |
|         | 0          | K      | − K | + 0         | 0 0           |                  | σ˜       | h˜n (un     | n          | )++h˜n | (un       | n )−     |     |
|         |            |        |     |             |               |                  |          | K           | σ K,σ      |        | L σ       | K,σ      |     |
|         | t          | ∆      | t˜  | (cid:96)2   | K˜            |                  | | |      |             | ·          |        |           | ·        |     |
|         | 0          |        |     |             | 0             | σ ∂K             |          |             |            |        |           |          |     |
|         |            |        |     |             | |             | | (cid:88)⊂      | (cid:16) |             |            |        |           | (cid:17) |     |
|         |            |        |     | gh2         | (cid:96)2 t 1 |                  |          | (cid:101)σ˜ |            |        | (cid:101) |          |     |
|         |            |        |     | 0           | 0 0           |                  | σ˜ γ∆t˜  |             | h˜n        | h˜n    | h˜n       |          |     |
|         |            |        |     |             |               |                  |          | |           | | σ        | L      | K =       | 0        |     |
|         |            |        |     | − (cid:96)4 | K˜            |                  | | |      | D˜          |            | −      |           |          |     |
|         |            |        |     |             | 0             |                  |          |             | σ          |        |           |          |     |
|         |            |        |     |             | |             | | σ (cid:88)⊂ ∂K |          | |           | | (cid:16) |        | (cid:17)  |          |     |
| where t | = (cid:96) | /u and | c   | = √gh       | so that       |                  |          |             |            |        |           |          |     |
|         | 0 0        | 0      | 0   |             | 0             |                  |          |             |            |        |           |          |     |
|         |            | h˜n+1  | h˜n |             |               |                  |          |             |            |        |           |          |     |
1
|     |     | K   | −    | K + |     | σ˜  | h˜n (un | n   | )++h˜n |     | (un n | )−  |     |
| --- | --- | --- | ---- | --- | --- | --- | ------- | --- | ------ | --- | ----- | --- | --- |
|     |     |     |      |     |     |     | K       | σ   | K,σ    | L   | σ K,σ |     |     |
|     |     |     | ∆ t˜ | K˜  |     | | | |         | ·   |        |     | ·     |     |     |
σ ∂K
|     |     |     |     | |   | | (cid:88)⊂  |               | (cid:16) |           |            |     |           | (cid:17) |     |
| --- | --- | --- | --- | --- | ------------ | ------------- | -------- | --------- | ---------- | --- | --------- | -------- | --- |
|     |     |     |     | gh  | t2           | 1             |          | (cid:101) | σ˜         |     | (cid:101) |          |     |
|     |     |     |     |     | 0 0          |               | σ˜       | γ ∆t˜     | h˜n        | h˜n | h˜n       |          |     |
|     |     |     |     |     |              |               |          | |         | | σ        | L   | K =       | 0        |     |
|     |     |     |     | −   | (cid:96)2 K˜ |               | | |      | D˜        |            | −   |           |          |     |
|     |     |     |     |     | 0            |               |          |           | σ          |     |           |          |     |
|     |     |     |     |     | |            | | σ (cid:88)⊂ | ∂K       | |         | | (cid:16) |     | (cid:17)  |          |     |
with
|          |     |             |     |          |      | gh t2     | c2      |     |     |     |     |     |     |
| -------- | --- | ----------- | --- | -------- | ---- | --------- | ------- | --- | --- | --- | --- | --- | --- |
|          |     |             |     |          |      | 0         |         |     | 1   |     |     |     |     |
|          |     |             |     |          |      |           | 0 =     | 0 = |     |     |     |     |     |
|          |     |             |     |          |      | (cid:96)2 | u2      | M2  |     |     |     |     |     |
|          |     |             |     |          |      | 0         |         | 0   |     |     |     |     |     |
| Denoting | τ = | (cid:96) /c | the | acoustic | time | step,     | we have |     |     |     |     |     |     |
0 0 0
|          |     |         |      | ∆t  |          | ∆t  | u        | ∆t          | ∆t     |       |     |          |     |
| -------- | --- | ------- | ---- | --- | -------- | --- | -------- | ----------- | ------ | ----- | --- | -------- | --- |
|          |     |         | ∆t˜= |     | =        |     | = 0      |             | = M    | = M∆τ |     |          |     |
|          |     |         |      | t   | (cid:96) | /u  | c        | (cid:96) /c | τ      |       |     |          |     |
|          |     |         |      |     | 0        | 0 0 | 0        | 0 0         |        | 0     |     |          |     |
| Then, we | get |         |      |     |          |     |          |             |        |       |     |          |     |
|          |     | 1 h˜n+1 |      | h˜n | 1        |     |          |             |        |       |     |          |     |
|          |     |         | K    | K   |          |     | h˜n      | (un         | )++h˜n |       | (un |          |     |
|          |     |         | −    | +   |          | σ˜  |          |             | n      |       | n   | )−       |     |
|          |     | M       | ∆ τ  |     | K˜       |     | K        | σ           | K,σ    | L     | σ   | K,σ      |     |
|          |     |         |      |     |          | |   | |        | ·           |        |       | ·   |          |     |
|          |     |         |      |     | | | σ    | ∂K  | (cid:16) |             |        |       |     | (cid:17) |     |
(cid:88)⊂
|            |     |     |      |     | γ∆t˜   |           |     | (cid:101)σ˜ |          |     |           |     |     |
| ---------- | --- | --- | ---- | --- | ------ | --------- | --- | ----------- | -------- | --- | --------- | --- | --- |
|            |     |     |      |     |        | 1         |     |             |          |     | (cid:101) |     |     |
|            |     |     |      |     |        |           | σ˜  | | | h˜n     | h˜n      | h˜n | = 0       |     |     |
|            |     |     |      |     | M2 K˜  |           |     | D˜          | σ L      | K   |           |     |     |
|            |     |     |      | −   |        |           | | | |             |          | −   |           |     |     |
|            |     |     |      |     | |      | | σ       | ∂K  | | σ |       |          |     |           |     |     |
|            |     |     |      |     |        | (cid:88)⊂ |     |             | (cid:16) |     | (cid:17)  |     |     |
| using ∆t˜= | M∆τ | we  | have | the | result |           |     |             |          |     |           |     |     |
Similarly,
LemmaC.1.2(Dimensionlessvelocityequation). Thedimensionlessvelocityequationof (C.1)

APPENDIX C. LOW MACH NUMBER ASYMPTOTIC ANALYSIS OF SOME
STAGGERED SCHEMES 200
is given by
1 un σ +1 − un σ + 1 e˜ u˜n σ e (cid:48) − u˜n σ h˜u˜ γ ∆τ | σ˜ | h˜n(h˜n h˜n ) −
M ∆τ D˜ | | h˜n+1 − M D˜ σ L− K
(cid:101) (cid:101) | σ | e ⊂ (cid:80) ∂Dσ 1 (cid:18) h˜n σ σ˜ (cid:19)(cid:18) (cid:110)(cid:110) (cid:111)(cid:111) | σ | (cid:19)
+ σ | | h˜n h˜n n =
M2h˜n+1 D˜ L− K K,σ
σ σ
h˜n ∆τ ∂L | | (cid:16) ∂K (cid:17)
2 σ α | | (h˜un (h˜u˜)n n | | (h˜un (h˜u˜)n n n
h˜n
σ
+1 M
(cid:34) (cid:18) |(cid:102)
L˜
| (cid:16)
σ − L
(cid:17)
· L,σ
(cid:19)
−
(cid:18) |(cid:103)
K˜
| (cid:16)
σ − K
(cid:17)
· K,σ
(cid:19) (cid:35)
K,σ
(cid:101) (cid:101)
Proof. Using h , (cid:96) and u as scaling variables and assuming that Φ = gh, we get for the
0 0 0
discharge equation, which is given by
un σ + ∆ 1 − t un σ + | D 1 σ | e
⊂(cid:88)
∂Dσ (cid:18) un σ h e (cid:48) n σ − +1 un σ (cid:19) ( F e n · n Dσ,e )−+∆t h h n σ + n σ 1 ∇Φn σ , ∗ = 0
that, using the definition of ∇ Φn,
D ∗
u 0 un σ +1 − un σ + u 0 1 e˜ u˜n σ e (cid:48) − u˜n σ h u h˜u˜ γ t 0 gh2 0∆t˜ | σ˜ | h˜n(h˜n h˜n ) −
t
0 (cid:101)h˜n
∆
gh
t˜
(cid:101) σ˜
h
0
(cid:96)
0 |
D˜
σ | e ⊂ (cid:80) ∂Dσ
| |
(cid:18)
h˜n
σ
+1
(cid:19)(cid:18)
0 0
(cid:110)(cid:110) (cid:111)(cid:111)
− (cid:96)
0 |
D˜
σ |
σ L− K
(cid:19)
+ σ 0 | | h˜n h˜n n =
h˜n
σ
+1 (cid:96)
0
D˜
σ
L− K K,σ
h˜n gh | σ˜ | (cid:16) t u (cid:17) ∂L ∂K
2 σ α 0 | | ∆t˜ 0 0 | | (h˜un (h˜u˜)n n | | (h˜un (h˜u˜)n n n
h˜n
σ
+1 (cid:96)
0 |
D˜
σ |
(cid:96)
0 (cid:34) (cid:18) |(cid:102)
L˜
| (cid:16)
σ − L
(cid:17)
· L,σ
(cid:19)
−
(cid:18) |(cid:103)
K˜
| (cid:16)
σ − K
(cid:17)
· K,σ
(cid:19) (cid:35)
K,σ
(cid:101) (cid:101)
Since c2 = gh , it gives
0 0
un σ +1 − un σ + 1 e˜ u˜n σ e (cid:48) − u˜n σ h˜u˜ γ ∆t˜ | σ˜ | h˜n(h˜n h˜n ) −
∆t˜ D˜ | | h˜n+1 − M2 D˜ σ L− K
(cid:101) 1 h˜(cid:101)n σ˜ | σ | e ⊂ (cid:80) ∂Dσ (cid:18) σ (cid:19)(cid:18) (cid:110)(cid:110) (cid:111)(cid:111) | σ | (cid:19)
+ σ | | h˜n h˜n n =
M2h˜n+1 D˜ L− K K,σ
σ σ
h˜n ∆t˜ | | ∂ (cid:16) L (cid:17) ∂K
2 σ α | | (h˜un (h˜u˜)n n | | (h˜un (h˜u˜)n n n
h˜n
σ
+1 M2
(cid:34) (cid:18) |(cid:102)
L˜
| (cid:16)
σ − L
(cid:17)
· L,σ
(cid:19)
−
(cid:18) |(cid:103)
K˜
| (cid:16)
σ − K
(cid:17)
· K,σ
(cid:19) (cid:35)
K,σ
(cid:101) (cid:101)
With ∆t˜= M∆τ, we get the result
Combining these lemmas we have
Proposition C.1.1 (Discrete wave system coupling the first order height and the zeroth order
velocity). The scheme (C.1) is asymptotically consistent with the following discrete wave system

APPENDIX C. LOW MACH NUMBER ASYMPTOTIC ANALYSIS OF SOME
STAGGERED SCHEMES 201
coupling the first order height and the zeroth order velocity :
(h˜(1) )n+1 (h˜(1) )n h˜(0) h˜(0) σ˜
K − K + σ˜ (u(0))n n = γ∆τ σ˜ | | (h˜(1) )n (h˜(1) )n
 

∆τ | K˜ | σ (cid:88)⊂ ∂K | | σ · K,σ | K˜ | σ (cid:88)⊂ ∂K | | | D˜ σ | (cid:16) L − K (cid:17)
  (cid:101)

       (u( σ 0) )n+ ∆ 1 − τ (u( σ 0) )n + D | ˜ σ˜ | (h˜( L 1) )n − (h˜( K 1) )n n K,σ =
σ
| | (cid:104) (cid:105)
(cid:101) (cid:101)∂L ∂K
        2α∆τh˜(0) (cid:34) (cid:18) | |(cid:102) L˜ | | (cid:16) (u( σ 0) )n − (u( K 0) )n (cid:17) · n L,σ (cid:19) − (cid:18) | |(cid:103) K˜ | | (cid:16) (u( σ 0) )n − (u( K 0) )n (cid:17) · n K,σ (cid:19) (cid:35) n K,σ

  (cid:101) (cid:101) (cid:101) (cid:101)
Proof. Injecting, the development as power of M of h˜ and u˜, Lemma C.1.2 at order 1/M2 gives
σ˜
| | (h˜(0) )n (h˜(0) )n = 0
D˜ L − K
σ
| | (cid:16) (cid:17)
so that (h˜(0))n is uniform. Then, Lemma C.1.2 at order 1/M gives
(u(0) )n+1 (u(0) )n σ˜
σ − σ + | | (h˜(1) )n (h˜(1) )n n =
∆τ D˜ L − K K,σ
σ
| | (cid:104) (cid:105)
(cid:101) ∂L(cid:101) ∂K
2α∆τ | | (h˜un (h˜u˜)n n | | (h˜un (h˜u˜)n n n
(cid:34)
L˜ σ − L · L,σ − K˜ σ − K · K,σ
(cid:35)
K,σ
(cid:18) |(cid:102)| (cid:16) (cid:17) (cid:19) (cid:18) |(cid:103)| (cid:16) (cid:17) (cid:19)
(cid:101) (cid:101)
and Lemma C.1.1 at order 1/M gives
(h˜(1) )n+1 (h˜(1) )n h˜(0) h˜(0) σ˜
K − K + σ˜ (u(0))n n = γ∆τ σ˜ | | (h˜(1) )n (h˜(1) )n
∆τ K˜ | | σ · K,σ K˜ | | D˜ L − K
| | σ (cid:88)⊂ ∂K | | σ (cid:88)⊂ ∂K | σ | (cid:16) (cid:17)
(cid:101)
(C.5)
leading to the result.
So this scheme has diffusion on both equations that depends on pressure AND velocity
terms, it nonetheless differs from our stabilization since it is not a grad-divergence operator.
Moreover, from (C.2), (C.3), we are not able to conclude that if
1
div(u) := σ u(x ) n = 0
K σ K,σ
K | | ·
| | σ ∂K
(cid:88)⊂
with x barycenter of the face σ, for all cells K implies that the diffusion vanishes. In fact, if
σ
we restrict to a 2d cartesian grid such that ∆x = ∆y on a periodic domain we see that for a
fixed face σ
∂L ∂K
| | (u(0))n (u(0) )n n | | (u(0))n (u(0) )n n n =
(cid:34)
L˜ σ − K · L,σ − K˜ σ − K · K,σ
(cid:35)
K,σ
(cid:18) |(cid:102)| (cid:16) (cid:17) (cid:19) (cid:18) |(cid:103)| (cid:16) (cid:17) (cid:19)
(cid:101) (cid:101) (cid:101) (cid:101)

APPENDIX C. LOW MACH NUMBER ASYMPTOTIC ANALYSIS OF SOME
STAGGERED SCHEMES 202
4
(u(0))n (u(0) )n n (u(0))n (u(0) )n n n =
∆x σ − K · L,σ − σ − K · K,σ K,σ
(cid:34) (cid:35)
(cid:18) (cid:19) (cid:18) (cid:19)
(cid:16) (cid:17) (cid:16) (cid:17)
(cid:101) (cid:101) (cid:101) (cid:101)
4
(u(0))n (u(0) )n ( n ) (u(0))n (u(0) )n n n =
∆x σ − K · − K,σ − σ − K · K,σ K,σ
(cid:34) (cid:35)
(cid:18) (cid:19) (cid:18) (cid:19)
(cid:16) (cid:17) (cid:16) (cid:17)
(cid:101) (cid:101) (cid:101) (cid:101)
4
(u(0) )n n 2(u(0))n n +(u(0) )n n n
∆x L · K,σ − σ · K,σ K · K,σ K,σ
(cid:34) (cid:35)
(cid:101) (cid:101) (cid:101)
C.2 Second example: the implicit scheme of Herbin .R, Kheriji.
W and Latch´e .J-C
For example, the implicit scheme of [66] is given by
K
| | ρn+1 ρn + σ n+1 n = 0
∆t K − K | |F σ · K,σ
 σ ∂K
    | D ∆ σ t | (cid:0) ρn σ +1un σ +1 (cid:1) − ρn σ (cid:88)⊂ un σ + F e n+1 · n Dσ,e un e +1+ | D σ | ( ∇ p)n σ +1 = | D σ | (∆u) n σ +1

(cid:0) (cid:1)
e
⊂(cid:88)
∂Dσ



(C.6)
Then,
n+1 n := ρn+1(un+1 n )++ρn+1(un+1 n )
F σ · K,σ K σ · K,σ L σ · K,σ −
F e n+1 · n Dσ,e un e +1 := un σ +1 F e n+1 · n Dσ,e + +un σ e (cid:48) +1 F e n+1 · n Dσ,e −
To define the mass flux through a d(cid:0)ual face e ou(cid:1)tward D (cid:0)the objective(cid:1)is that the flux
e σ
F
satisfies a finite volume discretization over diamond cell
D
| σ | ρn+1 ρn + n+1 n = 0
∆t σ − σ F e · Dσ,e
(cid:0) (cid:1)
e
⊂(cid:88)
∂Dσ
The pressure gradient is given by
σ
( p)n+1 = | | pn+1 pn+1 n
∇ σ D L − K K,σ
σ
| |
(cid:0) (cid:1)
and the stabilization term is given by
D (∆u) n+1 = νhd 2 un+1 un+1
| σ | σ e− σ e (cid:48) − σ
e ⊂(cid:88) ∂Dσ0 (cid:16) (cid:17)

APPENDIX C. LOW MACH NUMBER ASYMPTOTIC ANALYSIS OF SOME
STAGGERED SCHEMES 203
where
1
νhd 2 = n+1
e− e
2|F |
Lemma C.2.1 (Dimensionlessdensityequation). The dimensionless density equation of (C.6)
is given by
1 ρ˜n+1 ρ˜n 1
M K ∆ − τ K + K˜ | σ˜ | ρ˜n K un σ +1 · n K,σ + +ρ˜n L un σ +1 · n K,σ − = 0
| | σ (cid:88)⊂ ∂K (cid:16) (cid:0) (cid:1) (cid:0) (cid:1) (cid:17)
(cid:101) (cid:101)
Proof. Using ρ , (cid:96) and u as scaling variables
0 0 0
ρ ρ˜n+1 ρ˜n ρ u (cid:96) 1
t 0 K ∆ − t˜ K + 0 (cid:96)2 0 0 K˜ | σ˜ | ρ˜n K un σ +1 · n K,σ + +ρ˜n L un σ +1 · n K,σ − = 0
0 0 | | σ (cid:88)⊂ ∂K (cid:16) (cid:0) (cid:1) (cid:0) (cid:1) (cid:17)
(cid:101) (cid:101)
where t = (cid:96) /u
0 0 0
ρ˜n+1 ρ˜n 1
K ∆ − t˜ K + K˜ | σ˜ | h˜n K un σ +1 · n K,σ + +ρ˜n L un σ +1 · n K,σ − = 0
| | σ (cid:88)⊂ ∂K (cid:16) (cid:0) (cid:1) (cid:0) (cid:1) (cid:17)
(cid:101) (cid:101)
Denoting τ = (cid:96) /c the acoustic time step, we have
0 0 0
∆t ∆t u ∆t ∆t
∆t˜= = = 0 = M = M∆τ
t (cid:96) /u c (cid:96) /c τ
0 0 0 0 0 0 0
LemmaC.2.2(Dimensionlessvelocityequation). Thedimensionlessvelocityequationof (C.6)
is given by
1 ρ˜n+1 ρ˜n 1
M K ∆ − τ K + K˜ | σ˜ | ρ˜n K un σ +1 · n K,σ + +ρ˜n L un σ +1 · n K,σ − = 0
| | σ (cid:88)⊂ ∂K (cid:16) (cid:0) (cid:1) (cid:0) (cid:1) (cid:17)
(cid:101) (cid:101)
Proof. Using ρ , (cid:96) and u as scaling variables. The velocity equation is given by
0 0 0
un σ +
∆
1 −
t
un σ +
D
1
σ (cid:32)
un σ e (cid:48) +1
ρn
σ
−
+1
un σ +1
(cid:33)
(
F e
n
·
n
Dσ,e
)−+
ρn
σ
1
+1
(
∇
p)n
σ
+1 =
ρn
σ
1
+1
(∆u) n
σ
+1
| | e
⊂(cid:88)
∂Dσ
so that, using the definitions
u 0 un σ +1 − un σ + u 0 1 e˜ u˜n σ e (cid:48) +1 − u˜n σ +1 ρ u ρ˜u˜ n+1 −
t
0
∆t˜ ρ
0
(cid:96)
0
|
D˜
σ
|
e
⊂
∂Dσ
| |
(cid:32)
ρ˜n
σ
+1
(cid:33)
(cid:16)
0 0 {{ }}
(cid:17)
(cid:101) (cid:101) (cid:80)
p σ˜ 1 ρ u2
+ 0 | | p˜n+1 p˜n+1 n = 0 0 e˜ ρ˜u˜ n+1 u˜n+1 u˜n+1
(cid:96) 0 | D˜ σ | ρ 0 ρ˜n σ +1 (cid:0) L − K (cid:1) K,σ | D˜ σ | (cid:96) 0 |2ρ 0 ρ˜n σ +1 e ⊂(cid:88) ∂Dσ | |(cid:12) (cid:12) (cid:12) {{ }} (cid:12) (cid:12) (cid:12) (cid:16) σ e (cid:48) − σ (cid:17)
(cid:12) (cid:12)
(cid:12) (cid:12)

APPENDIX C. LOW MACH NUMBER ASYMPTOTIC ANALYSIS OF SOME
STAGGERED SCHEMES 204
ρ p u
Let γ = c2 0 and recall that c2 = 0 and M = 0 it yields
0p 0 ρ c
0 0 0
un σ +1 − un σ + 1 e˜ u˜n σ e (cid:48) +1 − u˜n σ +1 ρ˜u˜ n+1 −
∆t˜
|
D˜
σ
|
e
⊂
∂Dσ
| |
(cid:32)
ρ˜n
σ
+1
(cid:33)
(cid:16)
{{ }}
(cid:17)
(cid:101) (cid:101) (cid:80)
1 σ˜ 1
+ | | p˜n+1 p˜n+1 n = e˜ ρ˜u˜ n+1 u˜n+1 u˜n+1
γM2 | D˜ σ | ρ˜n σ +1 (cid:0) L − K (cid:1) K,σ | D˜ σ | 2ρ˜n σ +1| e ⊂(cid:88) ∂Dσ | |(cid:12) (cid:12) (cid:12) {{ }} (cid:12) (cid:12) (cid:12) (cid:16) σ e (cid:48) − σ (cid:17)
(cid:12) (cid:12)
With ∆t˜= M∆τ, we get (cid:12) (cid:12)
M
1 un σ +
∆
1 −
τ
un σ +
D˜
1
|
e˜
|
u˜n σ
ρ˜
e (cid:48)
n
−
+1
u˜n σ (
{{
ρ˜u˜
}}
)−
| σ | e ⊂ ∂Dσ (cid:18) σ (cid:19)
(cid:101) (cid:101) (cid:80) (C.7)
1 σ˜ 1
+ | | (p˜n p˜n )n = e˜ ρ˜u˜ n+1 u˜n+1 u˜n+1
γM2 | D˜ σ | ρ˜n σ +1 L− K K,σ | D˜ σ | 2ρ˜n σ +1| e ⊂(cid:88) ∂Dσ | |(cid:12) (cid:12)
(cid:12)
{{ }} (cid:12) (cid:12)
(cid:12)
(cid:16) σ e (cid:48) − σ (cid:17)
(cid:12) (cid:12)
(cid:12) (cid:12)
Gathering these lemmas we get
Proposition C.2.1 (Discrete wave system coupling the first order height and the zeroth order
velocity). The scheme (C.6) is asymptotically consistent with the following discrete wave system
coupling the first order height and the zeroth order velocity :
(ρ˜ (1) )n+1 (ρ˜ (1) )n ρ˜(0)
K − K + σ˜ (u(0))n+1 n = 0
 ∆τ K˜ | | σ · K,σ


| | σ (cid:88)⊂ ∂K
  (cid:101)

  (u(0) )n+1 (u(0) )n σ˜
(ρ˜ ( σ 0) )n+1 σ ∆ − τ σ + γ | D˜ | (ρ˜ ( L 1) )n+1 − (ρ˜ ( K 1) )n+1 n K,σ = 0
   | σ | (cid:104) (cid:105)
  (cid:101) (cid:101)

Proof. If we suppose for clarity that p(ρ) = ρ then, injecting, the expansion as power of M of
ρ˜and u˜, Lemma C.2.2 at order 1/M2 gives
σ˜
| | (ρ˜ (0) )n+1 (ρ˜ (0) )n+1 = 0
D˜ L − K
σ
| | (cid:16) (cid:17)
so that (ρ˜(0))n+1 is uniform. Then, Lemma C.2.2 at order 1/M gives
(u(0) )n+1 (u(0) )n σ˜
σ − σ + | | (ρ˜ (1) )n+1 (ρ˜ (1) )n+1 n = 0
∆τ γ D˜ L − K K,σ
σ
| | (cid:104) (cid:105)
(cid:101) (cid:101)

| APPENDIX    | C. LOW   | MACH     | NUMBER | ASYMPTOTIC |              | ANALYSIS | OF SOME |       |
| ----------- | -------- | -------- | ------ | ---------- | ------------ | -------- | ------- | ----- |
| STAGGERED   | SCHEMES  |          |        |            |              |          |         | 205   |
| and (C.2.1) | at order | 1/M      | gives  |            |              |          |         |       |
|             |          | (1)      | (1)    |            |              |          |         |       |
|             |          | (ρ˜ )n+1 | (ρ˜ )n | ρ˜(0)      |              |          |         |       |
|             |          | K        | K      |            | σ˜ (u(0))n+1 | n        |         |       |
|             |          |          | −      | +          |              | σ        | K,σ = 0 | (C.8) |
|             |          |          | ∆τ     | K˜         | | |          | ·        |         |       |
| | σ (cid:88)⊂ ∂K
(cid:101)
| leading to | the result. |     |     |     |     |     |     |     |
| ---------- | ----------- | --- | --- | --- | --- | --- | --- | --- |

| Appendix         |     | D            |       |           |     |     |
| ---------------- | --- | ------------ | ----- | --------- | --- | --- |
| Entropy          |     | for Explicit |       |           |     |     |
| Crouzeix-Raviart |     |              |       | staggered |     |     |
| discretizations  |     |              | using | ∇         | div |     |
stabilization
In this section, we show that adapting the grad-div diffusion introduced in chapter 5 to the
Crouzeix-Raviart/Rannacher-Turekstaggeredschemeenablestheobtentionofadiscreteglobal
entropy inequality for an Euler explicit time integration. Two main arguments are important
here: first, the approximation space and the choice of a particular convection term enables to
write a discrete kinetic energy balance [65, 66, 67, 68]. Then, the novelty in this set up is to
introduce the appropriate acoustic dependent diffusion; a grad-div type operator is added in
the momentum equation while a classical density Laplacian type operator is added in the mass
equation.
In this context, ρ is the density at the cell K , D h := D h + D h , where K,L
|     | K   |     |     | σ σ | K,σ K | L,σ L |
| --- | --- | --- | --- | --- | ----- | ----- |
|     |     |     |     | | | | | |   | | |   |
are the cells touching the face σ as in Figure D.1. u is the full velocity vector at the face σ,
σ
D σ is the dual cell associated to this face. Then e ∂D σ is a face on the boundary of D σ . For
⊂
a fixed cell K and a fixed face σ ∂K, L is the only cell such that L K = σ so that n K,σ
|     |     | ⊂   |     |     |     | ∩   |
| --- | --- | --- | --- | --- | --- | --- |
is oriented towards L. Similarly, for a fixed face σ and a fixed dual face e ∂D σ the only
σ e(cid:48)
⊂
face such that D D = e. For any cell of the mesh C dual or primal and n ∂C any dual
|     | σ   | σ(cid:48) |     |     |     |     |
| --- | --- | --------- | --- | --- | --- | --- |
|     | ∩   | e         |     |     |     | ⊂   |
of primal face n denotes the unit outward normal at the face n and cell C. ˜ will be the
C,n
F
| set of dual | faces e ∂D | σ . |     |     |     |     |
| ----------- | ---------- | --- | --- | --- | --- | --- |
σ
|     | ∈   | ∈ F |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- |
207

| APPENDIX        | D. ENTROPY | FOR EXPLICIT       | CROUZEIX-RAVIART | STAGGERED |     |
| --------------- | ---------- | ------------------ | ---------------- | --------- | --- |
| DISCRETIZATIONS | USING      | ∇div STABILIZATION |                  |           | 208 |
K
D
K,σ
D
σ L,σ
|     |     |     | n L |     |     |
| --- | --- | --- | --- | --- | --- |
K,σ
e
D
σ(cid:48) e
Figure D.1: Representation of a dual cell associated to a face σ for Crouzeix-
| Raviart/Rannacher-Turek |     | velocity staggering |     |     |     |
| ----------------------- | --- | ------------------- | --- | --- | --- |

| APPENDIX        |           | D. ENTROPY |       | FOR      | EXPLICIT |               | CROUZEIX-RAVIART |     |     |     |     | STAGGERED |     |
| --------------- | --------- | ---------- | ----- | -------- | -------- | ------------- | ---------------- | --- | --- | --- | --- | --------- | --- |
| DISCRETIZATIONS |           |            | USING |          | ∇div     | STABILIZATION |                  |     |     |     |     |           | 209 |
| D.1             | The       | numerical  |       |          | scheme   |               |                  |     |     |     |     |           |     |
| The             | numerical | scheme     |       | we study | is       |               |                  |     |     |     |     |           |     |
K
| |(ρn+1 |     | ρn )+ |     | n n | =   | 0,  |     |     |     |     |     |     |     |
| ------ | --- | ----- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| |      | K   | K     |     | σ   | K,σ |     |     |     |     |     |     |     |     |
| δτ     |     | −     |     | F · |     |     |     |     |     |     |     |     |     |
σ ∂K
|    |     |     | (cid:88)∈ |     |     |     |     |     |     |     |     |     |     |
| --- | --- | --- | --------- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
 
|      |     |     |     |     |     |     |     |     |     |     |     | ρn cn |          |
| ----- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | ----- | -------- |
|   D |     |     |     |     |     |     |     |     |     |     |     |       | (cid:94) |
 | σ |((ρ u )n+1 (ρ u )n)+ n n un + σ [[pn]] n = (cid:107) (cid:107)∞ max σ [[d ivun]] n ,
|     | σ   | σ   | σ   | σ   |     | e   | Dσ,e | e   |     | σ K,σ |     |       | σ K,σ |
| --- | --- | --- | --- | --- | --- | --- | ---- | --- | --- | ----- | --- | ----- | ----- |
| δτ  |     | −   |     |     |     | F · |      | |   | |   |       |     | 2 | | |       |
e ∂Dσ
|   |     |     |     |     | ⊂(cid:88) |     |     |     |     |     |     |     |     |
| --- | --- | --- | --- | --- | --------- | --- | --- | --- | --- | --- | --- | --- | --- |

 

(D.1)
where
|     |     |     |     |        | (cid:94) |     | 1     |     |     |     |     |     |       |
| --- | --- | --- | --- | ------ | -------- | --- | ----- | --- | --- | --- | --- | --- | ----- |
|     |     |     |     | div(u) |          | :=  |       | σ   | u   | n , |     |     | (D.2) |
|     |     |     |     |        |          | K   | ∂K    |     | σ   | K,σ |     |     |       |
|     |     |     |     |        |          |     |       | |   | | · |     |     |     |       |
|     |     |     |     |        |          |     | | | σ | ∂K  |     |     |     |     |       |
(cid:88)∈
|     | cn  |               |     | p(ρ), |          | ρn                   | maxρn. |     |      |           |     |     |     |
| --- | --- | ------------- | --- | ----- | -------- | -------------------- | ------ | --- | ---- | --------- | --- | --- | --- |
| and | :=  |               | max |       |          |                      | :=     |     | Also | we denote |     |     |     |
|     | max | [minρn,maxρn] |     |       | (cid:48) | (cid:107) (cid:107)∞ |        |     |      |           |     |     |     |
ρ
∈
(cid:112)
|     |     |       |     |      |     |        |       | un  | n   | +cn |      |       |     |
| --- | --- | ----- | --- | ---- | --- | ------ | ----- | --- | --- | --- | ---- | ----- | --- |
|     |     | n     |     | un   |     | ρn     |       | σ   | K,σ | max |      |       |     |
|     |     |       | n   | := σ | n   |        |       | | · | |   |     | σ (ρ | ρ ),  |     |
|     |     | F σ · | K,σ | | |  | σ · | K,σ {{ | }}σ − |     | 2   |     | | |  | L − K |     |
and
|     |      |     |     | n     | un   |          | n       | + un    | n       |         | un           |     |     |
| --- | ---- | --- | --- | ----- | ---- | -------- | ------- | ------- | ------- | ------- | ------------ | --- | --- |
|     |      |     |     | n     | Dσ,e | :=       |         |         | +       | −       | ,            |     |     |
|     |      |     |     | F e · |      | e        | G D σ,e | σ       | G D     | σ,e     | σ e (cid:48) |     |     |
| and | also |     |     |       |      | (cid:0)  |         | (cid:1) | (cid:0) | (cid:1) |              |     |     |
|     |      |     |     |       |      | ρ        | if      | un      | n       | 0       |              |     |     |
|     |      |     |     | ρupw  |      | K        |         | σ       | K,σ     |         | .            |     |     |
|     |      |     |     | σ     | :=   |          |         |         | ·       | ≥       |              |     |     |
|     |      |     |     |       |      |          | ρ       |         | else    |         |              |     |     |
|     |      |     |     |       |      | (cid:26) | Lσ      |         |         |         |              |     |     |
The dual density and the fluxes are chosen such that the following dual mass conservation
stands
|     |     |     |     |     | D σ |      |      |     |      |      |     |     |       |
| --- | --- | --- | --- | --- | --- | ---- | ---- | --- | ---- | ---- | --- | --- | ----- |
|     |     |     |     | |   | |   | ρn+1 | ρn + |     | n    | = 0. |     |     | (D.3) |
|     |     |     |     |     |     | σ    | σ    |     | Dσ,e |      |     |     |       |
|     |     |     |     |     | δt  | −    |      |     | G    |      |     |     |       |
e ∂Dσ
|     |     |     |     |     | (cid:0) |     | (cid:1) | ∈(cid:88) |     |     |     |     |     |
| --- | --- | --- | --- | --- | ------- | --- | ------- | --------- | --- | --- | --- | --- | --- |
In order to obtain this property we introduce the following flux (see [57] for details)
|     |     |     |     |     |     | n       | ωn  | n    | dΓ, |     |     |     |     |
| --- | --- | --- | --- | --- | --- | ------- | --- | ---- | --- | --- | --- | --- | --- |
|     |     |     |     |     |     | Dσ,e := |     | Dσ,e |     |     |     |     |     |
|     |     |     |     |     | G   |         |     | ·    |     |     |     |     |     |
e
(cid:90)
| where | ω is | a H(div;Ω) |     | field such | that |     |     |     |     |     |     |     |     |
| ----- | ---- | ---------- | --- | ---------- | ---- | --- | --- | --- | --- | --- | --- | --- | --- |
n
|     |     |     |     |     | div(ω)dx |     | =   |     | n K,σ | .   |     |     |     |
| --- | --- | --- | --- | --- | -------- | --- | --- | --- | ----- | --- | --- | --- | --- |
|     |     |     |     |     |          |     |     | F   | σ ·   |     |     |     |     |
K
|     |     |     |     |     | (cid:90) |     |     | σ ∂K |     |     |     |     |     |
| --- | --- | --- | --- | --- | -------- | --- | --- | ---- | --- | --- | --- | --- | --- |
(cid:88)∈
We do not discuss the existence of such field. On simplicial meshes and quadrangular/hexa-
hedral meshes, it actually can be constructed by using the Raviart-Thomas elements. In the

APPENDIX D. ENTROPY FOR EXPLICIT CROUZEIX-RAVIART STAGGERED
DISCRETIZATIONS USING ∇div STABILIZATION 210
following, the domain is assumed to be the d dimensional torus in order to simplify the com-
−
putations. Denote h = min(σ D ,K K ), ρ˜= min(ρ). The main result we want to
σ
∈ F| | ∈ C| |
address is:
Theorem D.1.1 (Global entropy dissipation). Let ψ in C2(R,R) function such that ρψ (ρ)
(cid:48)
−
ψ(ρ) = p(ρ). Suppose ψ strictly convex, then the scheme given by (D.1) yields
K D un+1 2 un 2
| | ψ(ρn+1) ψ(ρn ) + | σ | ρn+1| σ | ρn | σ | 0
δτ K − K δτ σ 2 − σ 2 ≤
K (cid:20) (cid:21) σ (cid:20) (cid:21)
(cid:88)∈C (cid:88)∈F
under the following CFL condition:
cn δτ ∂K 3 2 δt3ν
min max | | max 2 ρn + , max max ,
(cid:34) h (cid:107) (cid:107)∞(cid:32) ρ˜n+1 ρ˜n (cid:33) h ρ˜n+1 σ |G σ,e |
∈F
δt ∂K 3 (cn )2 δt ∂K
| | max cn (cn )2 , max | | max 2 un 2 1
h max
(cid:32)
2ρ˜n+1 min (ψ
(cid:48)(cid:48)
) max
(cid:33)
ρ˜n min (ψ
(cid:48)(cid:48)
)min un
σ 2
h (cid:107) (cid:107)∞(cid:35) ≤
[minρn,maxρn] [minρn,maxρn] σ | |
∈F
(D.4)
For that we introduce the following lemmas
1
Lemma D.1.1. Let div(ρu) := σ u(x ) n ρ , Then denoting u :=
K K | | σ · K,σ { } σ (cid:107) (cid:107)∞
| | σ ∂K
max u we have (cid:88)∈ (cid:8) (cid:9)
σ 2
σ | |
∈F
div(ρu) 2 | ∂K | max 2 ρn 2 ∂K di (cid:94) v(u) 2 + | ∂K | max u 2 σ [[ρ]]2.
(cid:107) (cid:107)L2(Ω) ≤ h (cid:107) (cid:107)∞ | | K h (cid:107) (cid:107)∞ | | σ
K σ
(cid:88)∈C (cid:88)∈F
ρ +ρ ρ ρ
Proof. Using
K Lσ
= ρ
K
−
Lσ
and the triangular inequality twice we have
K
2 − 2
ρ ρ
σ u n ρ = ρ σ u n σ u n
K
−
Lσ
| | σ · K,σ { } σ K | | σ · K,σ − | | σ · K,σ 2
(cid:12)σ ∂K (cid:12) (cid:12) σ ∂K σ ∂K (cid:12)
(cid:12) (cid:88)∈ (cid:8) (cid:9) (cid:12) (cid:12) (cid:88)∈ (cid:88)∈ (cid:12)
(cid:12) (cid:12) (cid:12) (cid:12) (D.5)
(cid:12) (cid:12) (cid:12) (cid:12)
[[ρ]]
σ
ρ σ u n + u σ | |,
σ K,σ
≤ (cid:107) (cid:107)∞ | | · (cid:107) (cid:107)∞ | | 2
(cid:12)σ ∂K (cid:12) σ ∂K
(cid:12) (cid:88)∈ (cid:12) (cid:88)∈
(cid:12) (cid:12)
(cid:12) (cid:12)

APPENDIX D. ENTROPY FOR EXPLICIT CROUZEIX-RAVIART STAGGERED
DISCRETIZATIONS USING ∇div STABILIZATION 211
so applying the square function to each side of the inequality (D.5) yields
2 [[ρ]] 2
σ
| σ | u σ · n K,σ { ρ } σ ≤ (cid:107) ρ (cid:107)∞ | σ | u σ · n K,σ + (cid:107) u (cid:107)∞ | σ | | 2 | . (D.6)
(cid:18)σ ∂K (cid:19) (cid:20) (cid:12)σ ∂K (cid:12) σ ∂K (cid:21)
(cid:88)∈ (cid:8) (cid:9) (cid:12) (cid:88)∈ (cid:12) (cid:88)∈
(cid:12) (cid:12)
(cid:12) 1 (cid:12)
Then (a+b)2 2a2+2b2 gives in (D.6) multiplied by
≤ K 2
| |
2 2
1 1
σ u n ρ 2 ρ 2 σ u n
K 2 | | σ · K,σ { } σ ≤ K 2 (cid:107) (cid:107)∞ | | σ · K,σ
| | (cid:18)σ ∂K (cid:19) | | (cid:18)σ ∂K (cid:19)
(cid:88)∈ (cid:8) (cid:9) (cid:88)∈
(D.7)
u 2 2
+(cid:107) (cid:107)∞ σ [[ρ]]
σ
.
2 | || |
(cid:20)σ ∂K (cid:21)
(cid:88)∈
2
Finally, using a discrete Jensen’s inequality σ [[ρ]] ∂K σ 2[[ρ]]2 in
σ max σ
| || | ≤ | | | |
(cid:20)σ ∂K (cid:21) σ ∂K
(D.7)
(cid:88)∈ (cid:88)∈
2 2
1 1
σ u n ρ 2 ρ 2 σ u n
K 2 | | σ · K,σ { } σ ≤ K 2 (cid:107) (cid:107)∞ | | σ · K,σ
| | (cid:18)σ ∂K (cid:19) | | (cid:18)σ ∂K (cid:19)
(cid:88)∈ (cid:8) (cid:9) (cid:88)∈
(D.8)
1 u 2
+
K 2
(cid:107)
2
(cid:107)∞
|
∂K
| max |
σ
|
2[[ρ]]2
σ
.
| | σ ∂K
(cid:88)∈
Multiplying by K (D.8) and summing on the cells K , we obtain the result.
| | ∈ C
We will also use these two lemmas proved in chapter 5 :
Lemma D.1.2 (Discrete integration by parts). Let u CRd(Ω) and p dG0(Ω), then
∈ ∈
i)
σ [[p]] n u = K pdivu + σ p u n .
| |
σ K,σ
·
σ
− | |
K
| |
Kσ σ
·
K,σ
σ int K σ ∂Ω
∈(cid:88)F (cid:88) (cid:88)∈
ii)
σ [[d (cid:93) ivu]] n u = ∂K (d (cid:93) ivu)2 + σ (d (cid:93) ivu) u n .
| |
σ K,σ
·
σ
− | |
K
| |
Kσ σ
·
K,σ
σ int K σ ∂Ω
∈(cid:88)F (cid:88) (cid:88)∈
Lemma D.1.3 (Inverse Poincar´e inequalities). Let u CRd(Ω) and p dG0(Ω), then
∈ ∈
i) for
2
σ [[p]] 2 ∂K σ [[p]]2,
σ max σ
| | ≤ | | | |
K (cid:18)σ ∂K (cid:19) σ
(cid:88)∈C (cid:88)∈ (cid:88)∈F

APPENDIX D. ENTROPY FOR EXPLICIT CROUZEIX-RAVIART STAGGERED
DISCRETIZATIONS USING ∇div STABILIZATION 212
(cid:93)
ii) let divu defined as (D.2) then
σ 2[[d (cid:93) ivu]]2 2 ∂K ∂K (d (cid:93) ivu)2 ,
σ max K
| | ≤ | | | |
σ int K
∈(cid:88)F (cid:88)∈C
D.2 Properties of the momentum equation
Thanks to the dual mass conservation, a discrete kinetic equation can be written
Lemma D.2.1 (Discrete kinetic transport equation, see [65]). If ρ the dual density interpol-
σ
ation and are chosen such that (D.3) stands, then a discrete kinetic equation is verified:
G
Dσ,e
D
| σ |((ρ u )n+1 (ρ u )n)+ un un
δτ σ σ − σ σ G Dσ,e {{ }}e · σ
(cid:34) (cid:35)
e
∈(cid:88)
∂Dσ
= | D σ |(ρn+1| un σ +1 | 2 ρn| un σ | 2 )+ un σ e (cid:48) · un σ | D σ |ρn+1 un+1 un 2
δτ σ 2 − σ 2 G Dσ,e 2 − 2δτ σ | σ − σ |
e
∈(cid:88)
∂Dσ
This yields the following kinetic energy global balance:
Proposition D.2.1 (Global kinetic energy equation). We have the following equality
D un+1 2 un 2
| σ |(ρn+1| σ | ρn | σ | )+ σ [[pn]] n un = R ,
δτ σ 2 − σ 2 | | σ K,σ · σ u
σ σ
(cid:88)∈F (cid:88)∈F
where
R
u
=
−
(cid:107) ρn (cid:107)∞
2
cn max
|
∂K
|
(d (cid:94) ivun)2
K −
|G D
2
σ,e |
|
un
σ e (cid:48) −
un
σ |
2+ |
2
D
δ
σ
τ
|ρn
σ
+1
|
un
σ
+1
−
un
σ |
2.
K e σ
(cid:88)∈C (cid:88)∈F (cid:88)∈F
Proof. Multiply the momentum equation of (D.1) by un gives
σ
D
| σ |((ρ u )n+1 (ρ u )n)+ un un+ σ [[pn]] n un =
δτ σ σ − σ σ G Dσ,e { } e · σ | | σ K,σ · σ
(cid:34) (cid:35)
e
∈(cid:88)
∂Dσ
(cid:8) (cid:9) (D.9)
(cid:107)
ρn
(cid:107)∞
2
cn
max
|
σ
|
[[d (cid:94) ivun]]
σ
n
K,σ ·
un
σ
+ |G D
2
σ,e |(un
σ e (cid:48) −
un
σ
)
·
un
σ
.
e
∈(cid:88)
∂Dσ

APPENDIX D. ENTROPY FOR EXPLICIT CROUZEIX-RAVIART STAGGERED
DISCRETIZATIONS USING ∇div STABILIZATION 213
Now with 2(un un) un = un 2 un 2 un un 2
σ e (cid:48) − σ · σ | σ e (cid:48) | −| σ | −| σ e (cid:48) − σ |
|G Dσ,e |(un un) un. = |G Dσ,e |( un 2 un 2) |G Dσ,e | un un 2
2 σ e (cid:48) − σ · σ 4 | σ e (cid:48) | −| σ | − 4 | σ e (cid:48) − σ |
e
∈(cid:88)
∂Dσ e
∈(cid:88)
∂Dσ e
∈(cid:88)
∂Dσ
(D.10)
Combining (D.10) and Lemma D.2.1 (D.9) becomes :
D un+1 2 un 2
| σ |(ρn+1| σ | ρn| σ | )+ G Dσ,e un un |G Dσ,e |( un 2 un 2) + σ [[pn]] n un
δτ σ 2 − σ 2 2 σ e (cid:48) · σ − 4 | σ e (cid:48) | −| σ | | | σ K,σ · σ
e
∈(cid:88)
∂Dσ (cid:18) (cid:19)
= (cid:107) ρn (cid:107)∞
2
cn max
|
σ
|
[[d (cid:94) ivun]]
σ
n
K,σ ·
un
σ −
|G D
4
σ,e |
|
un
σ e (cid:48) −
un
σ |
2+ |
2
D
δ
σ
τ
|ρn
σ
+1
|
un
σ
+1
−
un
σ |
2.
e
∈(cid:88)
∂Dσ
(D.11)
Since the domain is considered periodic,
G Dσ,e un un |G Dσ,e |( un 2 un 2) = 0
2 σ e (cid:48) · σ − 4 | σ e (cid:48) | −| σ |
σ
(cid:88)∈F
e
∈(cid:88)
∂Dσ (cid:18) (cid:19)
(D.12)
|G Dσ,e | un un 2 = |G Dσ,e | un un 2
4 | σ e (cid:48) − σ | 2 | σ e (cid:48) − σ |
σ (cid:88)∈F e ∈(cid:88) ∂Dσ e(cid:88)
∈F
˜
Summing (D.11) on σ and using (D.12) we get
∈ F
D un+1 2 un 2
| σ |(ρn+1| σ | ρn | σ | ) σ [[pn]] n un
δτ σ 2 − σ 2 | | σ K,σ · σ
σ σ
(cid:88)∈F (cid:88)∈F
= (cid:107) ρn (cid:107)∞
2
cn max
|
σ
|
[[d (cid:94) ivun]]
σ
n
K,σ ·
un
σ −
|G D
2
σ,e |
|
un
σ e (cid:48) −
un
σ |
2+ |
2
D
δ
σ
τ
|ρn
σ
+1
|
un
σ
+1
−
un
σ |
2.
σ (cid:88)∈F e(cid:88)
∈F
˜ σ (cid:88)∈F
(D.13)
The final result is obtained by integrating by parts the Grad-Div term in (D.13) with
Lemma 5.3.1 ii).

APPENDIX D. ENTROPY FOR EXPLICIT CROUZEIX-RAVIART STAGGERED
DISCRETIZATIONS USING ∇div STABILIZATION 214
We now focus on the non-negative kinetic energy term arising from the Euler forward
discretisation.
Lemma D.2.2. Let
D
L := | σ |ρn+1 un+1 un 2 (D.14)
u 2δτ σ | σ − σ |
σ
(cid:88)∈F
then
L u ≤ δt | ∂K h | maxcn max
(cid:32)
3 4 (cid:107) ρ ρ ˜n n + (cid:107) 1 ∞
(cid:33)
| ∂K | div (cid:94) (un) K 2 (cid:107) ρn (cid:107)∞ cn max
K
(cid:88)∈C
δt 3ν
+ max ( ) 2 [[un]] 2
h ρ˜n+1 | G σ,e − | | e |
(cid:32) (cid:33)
e(cid:88) ˜
∈F
δt ∂K 3
+ | | max cn (cn )3 σ [[ρn]]2
h max 2ρ˜n+1 max | | σ
(cid:32) (cid:33)
σ
(cid:88)∈F
Proof. Using the dual density balance equation (D.3) the convection term can be rewritten as
| D σ |((ρ u )n+1 (ρ u )n)+ un = | D σ |ρn+1(un+1 un)+
un
σ e (cid:48) −
un
σ .
δτ σ σ − σ σ G Dσ,e {{ }}e δτ σ σ − σ G Dσ,e 2
e
∈(cid:88)
∂Dσ e
∈(cid:88)
∂Dσ
As a consequence the momentum equation of (D.1) can be rewritten as
D
| δτ σ |ρn σ +1(un σ +1 − un σ ) = ( G σ,e ) − [[un]] e
e
∈(cid:88)
∂Dσ
(D.15)
ρn cn
(cid:94)
σ [[pn]]
σ
n
K,σ
+ (cid:107) (cid:107)∞ max σ [[div(un)]]
σ
n
K,σ
−| | 2 | |
Using (D.15) in the definition of R we have
u
2
δt ρn cn (cid:94)
L u = 2 D σ ρn σ +1 (cid:12) ( G σ,e ) − [[un]] e −| σ | [[pn]] σ n K,σ + (cid:107) (cid:107)∞ 2 max | σ | [[div(un)]] σ n K,σ (cid:12)
σ (cid:88)∈F | | (cid:12)
(cid:12)
e ∈(cid:88) ∂Dσ (cid:12)
(cid:12)
(cid:12) (cid:12)
Then, since, by the sta(cid:12)ndard Jensen inequality on the square function (cid:12)
2
ρn cn
(cid:94)
( σ,e ) − [[un]] e σ [[pn]] σ n K,σ + (cid:107) (cid:107)∞ max σ [[div(un)]] σ n K,σ
(cid:12) G −| | 2 | | (cid:12)
(cid:12)
(cid:12)
e ∈(cid:88) ∂Dσ (cid:12)
(cid:12)
(cid:12) (cid:12)
(cid:12) 2 2 ρn cn
(cid:94)
(cid:12) 2
3 ( σ,e ) − [[un]] e + σ [[pn]] σ n K,σ + (cid:107) (cid:107)∞ max σ [[div(un)]] σ n K,σ ,
≤ G | | 2 | |
(cid:34)(cid:12)
(cid:12)
e
∈(cid:88)
∂Dσ (cid:12)
(cid:12)
(cid:12)
(cid:12)
(cid:12)
(cid:12)
(cid:12)
(cid:12)
(cid:12)
(cid:12)
(cid:35)
(cid:12) (cid:12) (cid:12) (cid:12) (cid:12) (cid:12)
(cid:12) (cid:12) (cid:12) (cid:12) (cid:12) (cid:12)

| APPENDIX        | D.  | ENTROPY |       | FOR  | EXPLICIT      |     | CROUZEIX-RAVIART |     |     |     | STAGGERED |     |
| --------------- | --- | ------- | ----- | ---- | ------------- | --- | ---------------- | --- | --- | --- | --------- | --- |
| DISCRETIZATIONS |     |         | USING | ∇div | STABILIZATION |     |                  |     |     |     |           | 215 |
we have
|     |     |     |     |     | 3      | δt  |       |     |       |     |     |        |
| --- | --- | --- | --- | --- | ------ | --- | ----- | --- | ----- | --- | --- | ------ |
|     |     |     |     | L   |        | (L  | ) +(L | )   | +(L ) | ,   |     | (D.16) |
|     |     |     |     | u   | 2ρ˜n+1 | h   | u 1   | u 2 | u 3   |     |     |        |
≤
|     |     |     |     |     |     | (cid:20) |     |     |     | (cid:21) |     |     |
| --- | --- | --- | --- | --- | --- | -------- | --- | --- | --- | -------- | --- | --- |
with
2
|     |     |     |     | (L  |       |     |     | [[un]]  | ,   |     |     |     |
| --- | --- | --- | --- | --- | ----- | --- | --- | ------- | --- | --- | --- | --- |
|     |     |     |     | u   | ) 1 = |     | (   | σ,e ) − | e   |     |     |     |
G
|     |     |     |     |     |     | (cid:12)                |               |           | (cid:12) |     |     |     |
| --- | --- | --- | --- | --- | --- | ----------------------- | ------------- | --------- | -------- | --- | --- | --- |
|     |     |     |     |     |     | σ (cid:88)∈F (cid:12) e | ∈(cid:88) ∂Dσ |           | (cid:12) |     |     |     |
|     |     |     |     |     |     | (cid:12)                |               |           | (cid:12) |     |     |     |
| and |     |     |     |     |     | (cid:12)                |               |           | (cid:12) |     |     |     |
|     |     |     |     |     |     | (cid:12)                |               |           | (cid:12) |     |     |     |
|     |     |     |     |     | (L  |                         | σ             | 2[[pn]]2, |          |     |     |     |
|     |     |     |     |     | u   | ) 2 =                   |               | σ         |          |     |     |     |
| |
σ (cid:88)∈F
and
|     |     |     |     |     | ρn  | cn  | 2   |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
(cid:94)
|     |     |     | (L  | ) = | (cid:107) | (cid:107)∞ max |          | σ 2[[div(un)]]2 |     | .   |     |     |
| --- | --- | --- | --- | --- | --------- | -------------- | -------- | --------------- | --- | --- | --- | --- |
|     |     |     |     | u 3 |           |                |          |                 |     | σ   |     |     |
|     |     |     |     |     |           | 2              |          | | |             |     |     |     |     |
|     |     |     |     |     | (cid:18)  |                | (cid:19) | σ               |     |     |     |     |
(cid:88)∈F
| Then, by | Jensen | inequality |     |      |     |     |     |       |           |     |     |        |
| -------- | ------ | ---------- | --- | ---- | --- | --- | --- | ----- | --------- | --- | --- | ------ |
|          |        |            |     | (L ) | 2ν  |     | (   | ) 2   | [[un]] 2, |     |     | (D.17) |
|          |        |            |     | u    | 1   | max |     | σ,e − | e         |     |     |        |
|          |        |            |     |      | ≤   |     | | G | | |   | |         |     |     |        |
˜
e(cid:88)
∈F
| In parallel, | since   | by mean | value             | theorem; |           |          |       |     |                 |     |          |     |
| ------------ | ------- | ------- | ----------------- | -------- | --------- | -------- | ----- | --- | --------------- | --- | -------- | --- |
|              | [[pn]]  |         |                   |          | p(ρ)[[ρ]] |          |       | cn  |                 |     | p(ρ)     |     |
|              |         | σ       |                   | max      |           | (cid:48) | σ and |     | :=              | max | (cid:48) |     |
|              |         |         | ≤ ρ [minρn,maxρn] |          |           |          |       | max | ρ [minρn,maxρn] |     |          |     |
|              |         |         | ∈                 |          |           |          |       |     | ∈               |     |          |     |
| we have      | on (L u | )       |                   |          |           |          |       |     |                 |     |          |     |
2
|     |     |     |     | (L ) | ∂K  |     | (cn | )4  | σ [[ρn]]2 |     |     | (D.18) |
| --- | --- | --- | --- | ---- | --- | --- | --- | --- | --------- | --- | --- | ------ |
|     |     |     |     | u    | 2   | max | max |     |           | σ   |     |        |
|     |     |     |     |      | ≤ | |     | |   | |   | |         |     |     |        |
σ
(cid:88)∈F
| and finally | with | the Inverse |        | Poincar´e | inequality |           | Lemma          | 5.3.2 | :        |     |     |        |
| ----------- | ---- | ----------- | ------ | --------- | ---------- | --------- | -------------- | ----- | -------- | --- | --- | ------ |
|             |      |             |        | ∂K        |            |           |                |       | (cid:94) | 2   |     |        |
|             |      |             |        |           | max        | ρn        | cn             | )2    | iv(un)   |     |     |        |
|             |      |             | (L u ) | 3 |       | |          | (         |                |       | ∂K d     |     |     | (D.19) |
|             |      |             |        | ≤         | 2          | (cid:107) | (cid:107)∞ max |       | | |      | K   |     |        |
K (cid:88)∈C
Combining the bounds (D.17), (D.18) (D.19) in (D.16) we obtain the result.
| D.3    | Properties |           | of      | the | mass   | equation |             |     |        |     |     |     |
| ------ | ---------- | --------- | ------- | --- | ------ | -------- | ----------- | --- | ------ | --- | --- | --- |
| Now we | derive     | a similar | balance |     | on the | global   | ’potential’ |     | energy |     |     |     |
C2(R,R)
Proposition D.3.1 (Global ’potential’ energy estimate). Let ψ a function such that
| ρψ (ρ) | ψ(ρ) = | p(ρ), | then |     |     |     |     |     |     |     |     |     |
| ------ | ------ | ----- | ---- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
(cid:48) −
K
|     |     |     | |(ψ(ρn+1) |     | ψ(ρn |     |     | div(un) |     | p(ρn |     |     |
| --- | --- | --- | --------- | --- | ---- | --- | --- | ------- | --- | ---- | --- | --- |
|     |     |     | |         |     |      | ))+ |     | K       |     | )    | = R |     |
|     |     |     | δτ        | K   | −    | K   |     | | |     | K   | K    | ρ   |     |
|
|     |     | K          |     |     |     |     | K          |     |     |     |     |     |
| --- | --- | ---------- | --- | --- | --- | --- | ---------- | --- | --- | --- | --- | --- |
|     |     | (cid:88)∈C |     |     |     |     | (cid:88)∈C |     |     |     |     |     |

APPENDIX D. ENTROPY FOR EXPLICIT CROUZEIX-RAVIART STAGGERED
DISCRETIZATIONS USING ∇div STABILIZATION 216
with
cn 1
R = max σ ψ (ρn )[[ρn]]2 σ (un n ) ψ (ρupw )(ρupw ρn )2
ρ − 2 | | (cid:48)(cid:48) K,Lσ σ − 2 | | σ · K,σ − (cid:48)(cid:48) σ,K σ − K
K σ ∂K K σ ∂K
(cid:88)∈C (cid:88)∈ (cid:88)∈C (cid:88)∈
K
+ 2 | δτ |ψ (cid:48)(cid:48) (ρn K ,n+1 )(ρn K +1 − ρn K )2
K
(cid:88)∈C
where for K a fixed cell and σ ∂K,
∈
• ρn [min(ρn ,ρn ),max(ρn ,ρn )],
K,Lσ ∈ K Lσ K Lσ
• ρn,n+1 [max(ρn ,ρn+1),max(ρn ,ρn+1)],
K ∈ K K K K
• and ρupw [min(ρn ,ρupw ),max(ρn ,ρupw )]
σ,K ∈ K σ K σ
Proof. We multiply the density equation of (D.1) by ψ (ρn ) so that it becomes
(cid:48) K
K cn
| δτ |(ρn K +1 − ρn K )ψ (cid:48) (ρn K )+ | σ | un σ · n K,σ ρu σ pwψ (cid:48) (ρn K ) = m 2 ax | σ | ε K (σ)[[ρn]] σ ψ (cid:48) (ρn K )
σ ∂K σ ∂K
(cid:88)∈ (cid:88)∈
We separate it in three terms T +T = T
1 2 3
K
T 1 = | δτ |(ρn K +1 − ρn K )ψ (cid:48) (ρn K )
T = σ un n ρupwψ (ρn )
2 σ K,σ σ (cid:48) K
| | ·
σ ∂K
(cid:88)∈
cn
T = max σ ε (σ)[[ρn]] ψ (ρn )
3 K σ (cid:48) K
2 | |
σ ∂K
(cid:88)∈
By Taylor expansion
K K
T 1 = | δτ |(ψ(ρn K +1) − ψ(ρn K )) − 2 | δτ |ψ (cid:48)(cid:48) (ρn K ,n+1 )(ρn K +1 − ρn K )2. (D.20)
Then we rewrite T under the following form
2
T = σ un n ρupwψ (ρn )= σ un n ψ(ρupw)
2 σ K,σ σ (cid:48) K σ K,σ σ
| | · | | ·
σ ∂K σ ∂K
(cid:88)∈ (cid:88)∈
(D.21)
+ σ un n (ρupwψ (ρn ) ψ(ρupw)).
σ K,σ σ (cid:48) K σ
| | · −
σ ∂K
(cid:88)∈

| APPENDIX        |     | D. ENTROPY |       | FOR  | EXPLICIT |               | CROUZEIX-RAVIART |     |     |     | STAGGERED |     |     |
| --------------- | --- | ---------- | ----- | ---- | -------- | ------------- | ---------------- | --- | --- | --- | --------- | --- | --- |
| DISCRETIZATIONS |     |            | USING | ∇div |          | STABILIZATION |                  |     |     |     |           |     | 217 |
Now since p(ρn ) = ρn ψ (ρn ) ψ(ρn ), 0 = p(ρn ) (ρn ψ (ρn ) ψ(ρn )), so
|     |     | K   | K (cid:48) | K   | K      |         | K        | K   | (cid:48) | K   | K   |     |     |
| --- | --- | --- | ---------- | --- | ------ | ------- | -------- | --- | -------- | --- | --- | --- | --- |
|     |     |     |            | −   |        |         |          | −   |          | −   |     |     |     |
|     |     |     | ρupwψ      | (ρn | ψ(ρupw | )=ρupwψ |          | (ρn | ψ(ρupw   |     |     |     |     |
|     |     |     | σ          | )   |        | σ       | σ        | )   |          | σ ) |     |     |     |
|     |     |     | (cid:48)   | K − |        |         | (cid:48) | K   | −        |     |     |     |     |
(D.22)
|            |       |     |       |          |            |     | +p(ρn | ) (ρn | ψ (ρn      | )   | ψ(ρn | )), |     |
| ---------- | ----- | --- | ----- | -------- | ---------- | --- | ----- | ----- | ---------- | --- | ---- | --- | --- |
|            |       |     |       |          |            |     | K     |       | K (cid:48) | K   | K    |     |     |
|            |       |     |       |          |            |     |       | −     |            | −   |      |     |     |
| and since, | using | yet | again | a Taylor | expansion: |     |       |       |            |     |      |     |     |
1
ψ(ρn ) ψ(ρupw)+ψ (ρn )(ρupw ρn ) = ψ (ρupw )(ρupw ρn )2. (D.23)
|           |        | K   |            | σ         | (cid:48) | K      | σ K      |           | (cid:48)(cid:48) | σ,K  | σ     | K   |     |
| --------- | ------ | --- | ---------- | --------- | -------- | ------ | -------- | --------- | ---------------- | ---- | ----- | --- | --- |
|           |        |     | −          |           |          |        | −        | −2        |                  |      |       | −   |     |
| Combining | (D.23) |     | and (D.22) | in        | (D.21)   | yields |          |           |                  |      |       |     |     |
|           |        |     | T =        |           | σ un     | n      | ψ(ρupw)+ |           | σ                | un n | pn    |     |     |
|           |        |     | 2          |           | σ        | K,σ    | σ        |           |                  | σ    | K,σ K |     |     |
|           |        |     |            |           | | |      | ·      |          |           | |                | | ·  |       |     |     |
|           |        |     |            | σ ∂K      |          |        |          | σ         | ∂K               |      |       |     |     |
|           |        |     |            | (cid:88)∈ |          |        |          | (cid:88)∈ |                  |      |       |     |     |
(D.24)
|     |     |     |     | 1   |     |      | (ρupw                |        |     |        |     |     |     |
| --- | --- | --- | --- | --- | --- | ---- | -------------------- | ------ | --- | ------ | --- | --- | --- |
|     |     |     |     |     | σ   | un n | ψ                    | )(ρupw |     | ρn )2. |     |     |     |
|     |     |     |     |     |     | σ    | K,σ (cid:48)(cid:48) | σ,K    | σ   | K      |     |     |     |
|     |     |     |     | −2  | |   | | ·  |                      |        |     | −      |     |     |     |
|     |     |     |     | σ   | ∂K  |      |                      |        |     |        |     |     |     |
(cid:88)∈
| Finally, | with | a Taylor | expansion |     | on  | T we | obtain |     |     |     |     |     |     |
| -------- | ---- | -------- | --------- | --- | --- | ---- | ------ | --- | --- | --- | --- | --- | --- |
3
cn
|     |     | T   | = max |     | σ (ρn | ρn  | )ψ (ρn       | )   |     |     |     |     |     |
| --- | --- | --- | ----- | --- | ----- | --- | ------------ | --- | --- | --- | --- | --- | --- |
|     |     | 3   |       |     | Lσ    |     | K (cid:48) K |     |     |     |     |     |     |
|     |     |     | 2     |     | | |   | −   |              |     |     |     |     |     |     |
σ ∂K
(cid:88)∈
(D.25)
|     |     |     | cn  |     |          |              |     | ψ   | (ρn                   | )   |     |          |     |
| --- | --- | --- | --- | --- | -------- | ------------ | --- | --- | --------------------- | --- | --- | -------- | --- |
|     |     |     | max |     |          | (σ)[[ψ(ρn)]] |     |     | (cid:48)(cid:48) K,Lσ | (ρn |     | )2       |     |
|     |     |     | =   |     | σ ε      |              |     |     |                       |     | ρ   | ,        |     |
|     |     |     | 2   |     | | |      | K            | σ   | −   | 2                     | Lσ  | −   | K        |     |
|     |     |     |     |     | (cid:32) |              |     |     |                       |     |     | (cid:33) |     |
σ ∂K
(cid:88)∈

APPENDIX D. ENTROPY FOR EXPLICIT CROUZEIX-RAVIART STAGGERED
DISCRETIZATIONS USING ∇div STABILIZATION 218
So, gathering (D.25), (D.24) and (D.20) gives
K cn
| |(ψ(ρn+1) ψ(ρn ))+ σ un n ψ(ρupw)+ K div(un) pn + max σ ε (σ)[[ψ(ρn)]]
δτ K − K | | σ · K,σ σ | | K K 2 | | K K,σ
σ ∂K σ ∂K
(cid:88)∈ (cid:88)∈
cn K
= − m 2 ax | σ | ψ (cid:48)(cid:48) (ρn K,Lσ )(ρn Lσ − ρ K )2+ 2 | δτ |ψ (cid:48)(cid:48) (ρn K ,n+1 )(ρn K +1 − ρn K )2
σ ∂K
(cid:88)∈
1
+ σ un n ψ (ρupw )(ρupw ρn )2.
2 | | σ · K,σ (cid:48)(cid:48) σK σ − K
σ ∂K
(cid:88)∈
(D.26)
Note that because of conservativity and the fact that the domain has no boundary
cn
max σ ε (σ)[[ψ(ρn)]] = 0 σ un n ψ(ρupw) = 0, (D.27)
K K,σ σ K,σ σ
2 | | | | ·
K σ ∂K K σ ∂K
(cid:88)∈C (cid:88)∈ (cid:88)∈C (cid:88)∈
also because the density
ρupw
is upwinded:
σ
1 1
σ un n ψ (ρupw )(ρupw ρn )2 = σ (un n ) ψ (ρupw )(ρupw ρn )2
2 | | σ · K,σ (cid:48)(cid:48) σK σ − K −2 | | σ · K,σ − (cid:48)(cid:48) σK σ − K
σ ∂K σ ∂K
(cid:88)∈ (cid:88)∈
(D.28)
Summing on the cells (D.26) and using (D.27) and (D.28), we obtain the results
We take a look at the non-negative reminder in the ”potential” estimate,
Lemma D.3.1. Let
K
L ρ := 2 | δτ |ψ (cid:48)(cid:48) (ρn K ,n+1 )(ρn K +1 − ρn K )2, (D.29)
K
(cid:88)∈C

APPENDIX D. ENTROPY FOR EXPLICIT CROUZEIX-RAVIART STAGGERED
DISCRETIZATIONS USING ∇div STABILIZATION 219
then the following bound stands:
L (cn max )2δt | ∂K | max 2 ρn 2 ∂K div (cid:94) (un) 2
ρ ≤ ρ˜n h (cid:107) (cid:107)∞ | | K
K
(cid:88)∈C
(cn )2δt ∂K
+ max | | max 2 un 2 σ [[ρn]]2
ρ˜n h (cid:107) (cid:107)∞ | | σ
σ
(cid:88)∈F
(cn )2δt ∂K
+ max | | max (cn )2 σ [[ρn]]2.
ρ˜n h max | | σ
σ
(cid:88)∈F
Proof. for the first point we have
2
δτ un +cn
(ρn+1 ρn )2 = ( )2 σ un n ρn + | σ | max σ [[ρn]] , (D.30)
K − K K − | | σ · K,σ { } σ 2 | | σ
(cid:32) (cid:33)
| | σ ∂K σ ∂K
(cid:88)∈ (cid:8) (cid:9) (cid:88)∈
but by Jensen inequality for the square function:
2
un +cn
σ un n ρn + | σ | max σ [[ρn]]
− | | σ · K,σ { } σ 2 | | σ
(cid:32) (cid:33)
σ ∂K σ ∂K
(cid:88)∈ (cid:8) (cid:9) (cid:88)∈ (D.31)
( un +cn )2 2
2 K 2div(ρu)2 + | σ | max σ [[ρn]] .
≤ | | K 2 | | σ
(cid:18)σ ∂K (cid:19)
(cid:88)∈
Combining (D.31) and (D.30) in the definition of R yields, because ψ is convex
ρ
2
K δτ ( un +cn )2
L ρ ≤ 2 | δτ |ψ (cid:48)(cid:48) (ρn K ,n+1 )( K )2 2 | K | 2div(ρu)2 K + (cid:107) (cid:107)∞ 2 max (cid:107) σ | [[ρn]] σ
(cid:32) (cid:33)
K | | (cid:18) σ ∂K (cid:19)
(cid:88)∈C (cid:88)∈
(D.32)
2
( un +cn )2 1
≤(cid:107) ψ (cid:48)(cid:48) (cid:107) L∞(ρ˜n,ρ¯n) δt (cid:107) div(ρnun) (cid:107) 2 L2 + (cid:107) (cid:107)∞ 4 max K (cid:107) σ | [[ρn]] σ .
(cid:32) (cid:33)
(cid:20) K | | σ ∂K (cid:21)
(cid:88)∈C (cid:88)∈
Then Lemma D.1.1 on div(ρnun) 2 and by Inverse Poincar´e inequality in
(cid:107) (cid:107)L2

| APPENDIX        | D.  | ENTROPY |       | FOR  | EXPLICIT      |     | CROUZEIX-RAVIART |     |     |     | STAGGERED |     |     |
| --------------- | --- | ------- | ----- | ---- | ------------- | --- | ---------------- | --- | --- | --- | --------- | --- | --- |
| DISCRETIZATIONS |     |         | USING | ∇div | STABILIZATION |     |                  |     |     |     |           |     | 220 |
2
1
|     |     | σ [[ρn]] |     | , (D.32) | becomes |     |     |     |     |     |     |     |     |
| --- | --- | -------- | --- | -------- | ------- | --- | --- | --- | --- | --- | --- | --- | --- |
σ
K (cid:107) |
| K | |      | (cid:32) σ ∂K |            | (cid:33)                     |     |       |     |                       |     |          |      |     |     |     |
| ---------- | ------------- | ---------- | ---------------------------- | --- | ----- | --- | --------------------- | --- | -------- | ---- | --- | --- | --- |
| (cid:88)∈C | (cid:88)∈     |            |                              |     |       |     |                       |     |          |      |     |     |     |
|            |               |            |                              |     | δt ∂K |     |                       |     |          |      |     |     |     |
|            |               |            |                              |     |       | max |                       |     | (cid:94) |      | 2   |     |     |
|            |               | L ψ        |                              |     | | |   | 2   | ρn 2                  | ∂K  | div      | (un) |     |     |     |
|            |               | ρ          | (cid:48)(cid:48) L∞(ρ˜n,ρ¯n) |     | h     |     |                       |     |          | K    |     |     |     |
|            |               | ≤(cid:107) | (cid:107)                    |     |       |     | (cid:107) (cid:107) ∞ | |   | |        |      |     |     |     |
K
(cid:88)∈C
|     |     |           |                              |     | δt ∂K |           |             |     |         |     |     |     |     |
| --- | --- | --------- | ---------------------------- | --- | ----- | --------- | ----------- | --- | ------- | --- | --- | --- | --- |
|     |     |           |                              |     |       | max       |             |     | [[ρn]]2 |     |     |     |     |
|     |     | + ψ       |                              |     | | |   |           | un 2        | σ   |         |     |     |     |     |
|     |     |           | (cid:48)(cid:48) L∞(ρ˜n,ρ¯n) |     | h     |           |             |     | σ       |     |     |     |     |
|     |     | (cid:107) | (cid:107)                    |     |       | (cid:107) | (cid:107) ∞ | |   | |       |     |     |     |     |
σ
(cid:88)∈F
|     |     |     |     |     | δt ∂K |     | un  | +cn | )2  |     |     |     |     |
| --- | --- | --- | --- | --- | ----- | --- | --- | --- | --- | --- | --- | --- | --- |
max (
|     |     | + ψ       |                              |     | | | |     | (cid:107) (cid:107)∞ | max |     | σ [[ρn]]2 | .   |     |     |
| --- | --- | --------- | ---------------------------- | --- | --- | --- | -------------------- | --- | --- | --------- | --- | --- | --- |
|     |     |           | (cid:48)(cid:48) L∞(ρ˜n,ρ¯n) |     | h   |     |                      |     |     |           | σ   |     |     |
|     |     | (cid:107) | (cid:107)                    |     |     |     |                      | 2   |     | | |       |     |     |     |
σ
(cid:88)∈F
|     |     |     |     |     |     |     |     |     |     |     | (cn | )2  |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
max
by derivating ρψ ψ = p(ρ) we obtain ρψ = p(ρ) so, ψ L∞(ρ˜n,ρ¯n) . Combining
|     |     | (cid:48) − |     |     |     | (cid:48)(cid:48) | (cid:48) |     | (cid:107) (cid:48)(cid:48) (cid:107) |     | ≤   | ρ˜n |     |
| --- | --- | ---------- | --- | --- | --- | ---------------- | -------- | --- | ------------------------------------ | --- | --- | --- | --- |
|     | un  | +cn        | )2  |     |     |                  |          |     |                                      |     |     |     |     |
(
this with, (cid:107) (cid:107)∞ max ( un )2+(cn )2 we obtain the bound.
max
|     |     | 2   |     | ≤ (cid:107) | (cid:107)∞ |     |     |     |     |     |     |     |     |
| --- | --- | --- | --- | ----------- | ---------- | --- | --- | --- | --- | --- | --- | --- | --- |
Combining these estimates we are able to prove the global entropy dissipation theorem:
proof of Theorem D.1.1. Summing the ’potential’ energy balance given by Proposition D.3.1
| and the | global     | kinetic   | balance | given | by   | Proposition |     | D.2.1   | gives |     |       |     |     |
| ------- | ---------- | --------- | ------- | ----- | ---- | ----------- | --- | ------- | ----- | --- | ----- | --- | --- |
|         |            | K         |         |       |      |             | D   |         | un+1  | 2   | un 2  |     |     |
|         |            | |(ψ(ρn+1) |         |       | ψ(ρn |             | σ   | |(ρn+1| | σ     | ρn  | σ     |     |     |
|         |            | |         |         |       | ))+  |             | |   |         |       | |   | | | ) | = R |     |
|         |            | δτ        | K       | −     | K    |             | δτ  | σ       | 2     | −   | σ 2   |     |     |
|         | K          |           |         |       |      | σ           |     |         |       |     |       |     |     |
|         | (cid:88)∈C |           |         |       |      | (cid:88)∈F  |     |         |       |     |       |     |     |
where
cn
|            |           |           |                       |            | 1     |            |           |            |       | (ρupw              |        |         |     |
| ---------- | --------- | --------- | --------------------- | ---------- | ----- | ---------- | --------- | ---------- | ----- | ------------------ | ------ | ------- | --- |
| R = max    |           | σ ψ       | (ρn                   | )[[ρn]]2   |       |            |           | σ (un      | n     | ) ψ                | )(ρupw | ρn )2+L |     |
|            |           |           | (cid:48)(cid:48) K,Lσ |            | σ     |            |           | σ          | K,σ   | − (cid:48)(cid:48) | σ,K    | σ K     | ρ   |
| − 2        |           | | |       |                       |            | − 2   |            | |         | |          | ·     |                    |        | −       |     |
| K          | σ         | ∂K        |                       |            |       | K          | σ ∂K      |            |       |                    |        |         |     |
| (cid:88)∈C | (cid:88)∈ |           |                       |            |       | (cid:88)∈C | (cid:88)∈ |            |       |                    |        |         |     |
|            |           | ρn        | cn                    |            |       |            |           |            |       |                    |        |         |     |
|            |           |           |                       |            |       | (cid:94)   |           |            | D σ,e |                    |        |         |     |
|            |           | (cid:107) | (cid:107)∞ max        |            | ∂K (d | ivun)2     |           | |G         | | un  |                    | un 2+L |         |     |
|            |           |           |                       |            |       |            | K         |            |       | σ (cid:48)         | σ      | u       |     |
|            |           | −         | 2                     |            | | |   |            | −         |            | 2 |   | e −                | |      |         |     |
|            |           |           |                       | K          |       |            |           | e          |       |                    |        |         |     |
|            |           |           |                       | (cid:88)∈C |       |            |           | (cid:88)∈F |       |                    |        |         |     |
(D.33)

APPENDIX D. ENTROPY FOR EXPLICIT CROUZEIX-RAVIART STAGGERED
DISCRETIZATIONS USING ∇div STABILIZATION 221
with L given by (D.29), L given by (D.14). Then Lemma D.3.1 gives
ρ u
L (cn max )2δt | ∂K | max 2 ρn 2 ∂K div (cid:94) (un) 2
ρ ≤ ρ˜n h (cid:107) (cid:107)∞ | | K
K
(cid:88)∈C
(cn )2δt ∂K
+ m ρ˜ a n x | h | max 2 (cid:107) un (cid:107) 2 ∞ | σ | [[ρn]]2 σ (D.34)
σ
(cid:88)∈F
(cn )2δt ∂K
+ max | | max (cn )2 σ [[ρn]]2.
ρ˜n h max | | σ
σ
(cid:88)∈F
and by Lemma D.2.2
L u ≤ δt | ∂K h | maxcn max
(cid:32)
3 4 (cid:107) ρ ρ ˜n n + (cid:107) 1 ∞
(cid:33)
| ∂K | div (cid:94) (un) K 2 (cid:107) ρn (cid:107)∞ cn max
K
(cid:88)∈C
δt 3ν
+ h ρ˜n m + a 1 x | ( G σ,e ) − | 2 | [[un]] e | 2 (D.35)
(cid:32) (cid:33)
e(cid:88) ˜
∈F
δt ∂K 3
+ | | max cn (cn )3 σ [[ρn]]2
h max 2ρ˜n+1 max | | σ
(cid:32) (cid:33)
σ
(cid:88)∈F
Also
cn
max σ ψ (ρn )[[ρn]]2 min (ψ )cn σ [[ρn]]2. (D.36)
− 2 | | (cid:48)(cid:48) K,Lσ σ ≤ −[minρn,maxρn] (cid:48)(cid:48) max | | σ
K σ ∂K σ
(cid:88)∈C (cid:88)∈ (cid:88)∈F
Gathering (D.36), (D.35) and (D.34) in (D.33) yields
R Aacoustic,u+Aacoustic,ρ+Ctransport,u+Ctransport,ρ, (D.37)
≤
with
Aacoustic,u = δt | ∂K | maxcn ρn 3 + 2 1 ∂K d (cid:94) iv(un) 2 cn ρn ,
(cid:34)
h max(cid:107) (cid:107)∞(cid:32) 4ρ˜n+1 ρ˜n
(cid:33)
− 2
(cid:35)
| | K max(cid:107) (cid:107)∞
K
(cid:88)∈C
δt3ν 1
Ctransport,u = max [[un]]2,
h ρ˜n+1 |G σ,e |− 2 |G σ,e | e
(cid:32) (cid:33)
e(cid:88) ˜
∈F
and with
δt ∂K 3
Aacoustic,ρ = cn max | σ | [[ρn]]2 σ | h | max cn max (cid:32) 2ρ˜n+1 (cn max )2 (cid:33) −[minρ m n, i m n axρn] (ψ (cid:48)(cid:48) ) ,
σ (cid:20) (cid:21)
(cid:88)∈F

| APPENDIX        | D. ENTROPY |         | FOR  | EXPLICIT      | CROUZEIX-RAVIART |     | STAGGERED |     |     |
| --------------- | ---------- | ------- | ---- | ------------- | ---------------- | --- | --------- | --- | --- |
| DISCRETIZATIONS |            | USING   | ∇div | STABILIZATION |                  |     |           |     | 222 |
|                 | (cn        | )2δt ∂K |      |               |                  | 1   |           |     |     |
Ctransport,ρ m a x max un 2 σ [[ρn]]2 σ (un n ψ (ρu p w )(ρu pw ρn )2.
|     | =    | |   | | 2 |                     |     | σ   | σ K,σ ) − | (cid:48)(cid:48) | σ K |
| --- | ---- | --- | --- | ------------------- | --- | --- | --------- | ---------------- | --- |
|     | ρ˜ n |     | h   | (cid:107) (cid:107) | | | | −2  | | | ·     | σ, K             | −   |
∞
|     |     |     |     |     | σ (cid:88)∈F | K (cid:88)∈C σ (cid:88)∈ | ∂K  |     |     |
| --- | --- | --- | --- | --- | ------------ | ------------------------ | --- | --- | --- |
ThedefinitionsofAacoustic,u,Aacoustic,ρ,Ctransport,u,Ctransport,ρ yieldsR 0undertheCFL
≤
| (D.4): this | ends the | prood. |     |     |     |     |     |     |     |
| ----------- | -------- | ------ | --- | --- | --- | --- | --- | --- | --- |
Remark D.3.1. CFL (D.4) should be improved by bettering the bounds given by Lemma D.3.1
| and Lemma | D.2.2. |     |     |     |     |     |     |     |     |
| --------- | ------ | --- | --- | --- | --- | --- | --- | --- | --- |

| Appendix |     |     |            | E   |     |     |      |         |     |     |         |     |     |
| -------- | --- | --- | ---------- | --- | --- | --- | ---- | ------- | --- | --- | ------- | --- | --- |
| Solution |     |     |            | of  | the |     | 1d   | Riemann |     |     | Problem |     |     |
| in       |     | the | barotropic |     |     |     | case |         |     |     |         |     |     |
In this section we detail the resolution of the Riemann problem for Euler barotropic equations.
| In  | this | aim, we | briefly | adapt | the | work from | Toro      | [14] | to the | following |     |     |     |
| --- | ---- | ------- | ------- | ----- | --- | --------- | --------- | ---- | ------ | --------- | --- | --- | --- |
|     |      |         |         |       |     | ∂         | ρ+div(ρu) |      | = 0    |           |     |     |     |
t

|     |     |     |     |     |    | ∂ t (ρu)+div(ρu |     | u)+∇p |     | = 0 |     |     | (E.1) |
| --- | --- | --- | --- | --- | --- | --------------- | --- | ----- | --- | --- | --- | --- | ----- |
|     |     |     |     |     |   |                 |     | ⊗     |     |     |     |     |       |

 
|     |     |     |     |     |     |     | p(ρ) = | ργ;γ | 1   |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | ------ | ---- | --- | --- | --- | --- | --- |
|     |     |     |     |     |    |     |        | ≥    |     |     |     |     |     |
 

|     | ρ   |     |     |     | equation  |     | p   | f(ρ) |     |     | p   |     | u   |
| --- | --- | --- | --- | --- | ----------- | --- | --- | ---- | --- | --- | --- | --- | --- |
with is the density, the of state = links the pressure to the density and
is the the velocity. Classically, the Euler barotropic equations can be put under the form of a
| non-linear |     | system | of  | conservation |     | laws : |       |     |      |     |       |     |       |
| ---------- | --- | ------ | --- | ------------ | --- | ------ | ----- | --- | ---- | --- | ----- | --- | ----- |
|            |     |        |     |              |     |        |       | Rd  |      |     | Rd+1, |     |       |
|            |     |        |     | ∂ U+divF(U)  |     | = 0,   | U : Ω |     | [0,+ | [   |       |     | (E.2) |
t
|     |     |     |     |     |     |     |     | ⊂   | ×   | ∞   | −→  |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
with
ρut
|     |     |     |     |     |      |     |      |     | R(d+1) |     | d.  |     |       |
| --- | --- | --- | --- | --- | ---- | --- | ---- | --- | ------ | --- | --- | --- | ----- |
|     |     |     |     |     | F(U) | :=  |      |     |        | ×   |     |     | (E.3) |
|     |     |     |     |     |      | ρu  | u+pI |     | ∈      |     |     |     |       |
d
|     |          |     |        |         |        | (cid:18) | ⊗   |     | (cid:19) |     |     |     |       |
| --- | -------- | --- | ------ | ------- | ------ | -------- | --- | --- | -------- | --- | --- | --- | ----- |
| The | jacobian |     | of the | flux is | given  | by       |     |     |          |     |     |     |       |
|     |          |     |        |         |        |          |     | 0   |          |     | nt  |     |       |
|     |          |     | A(U,n) | :=      | ∇(F(U) | n)       | =   |     |          |     |     | .   | (E.4) |
|     |          |     |        |         |        | ·        |    |     |          |     | nI  |    |       |
|     |          |     |        |         |        |          | c2n | u   | nu       | u   | n+u |     |       |
d
|     |     |     |     |     |     |     |     | −   | ·   | ⊗   | ·   |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
|     |     |     |     |     |     |     |    |     |     |     |     |    |     |
Undertheassumptionthatp(ρ) > 0thesystemishyperbolic (Definition1.2.1)witheigenvalues
(cid:48)
| in  | the normalized |     | direction |     | n    |      |           |     |      |     |     |     |     |
| --- | -------------- | --- | --------- | --- | ---- | ---- | --------- | --- | ---- | --- | --- | --- | --- |
|     |                |     |           | λ   | := u | n, λ | := u n+c, |     | λ := | u n | c,  |     |     |
|     |                |     |           |     | 1    | 2    |           |     | 3    |     |     |     |     |
|     |                |     |           |     |      | ·    | ·         |     |      | ·   | −   |     |     |
224

APPENDIX E. SOLUTION OF THE 1D RIEMANN PROBLEM IN THE BAROTROPIC
CASE 225
where c := p(ρ). The eigenvalues are generally supposed to be either linearly degenerate or
(cid:48)
genuinely non-linear
(cid:112)
Definition E.0.1 (Genuinely non-linear fields). The eigenvalue λ (U) is genuinely non-linear
j
if it has multiplicity one and if
U, r ker[A(U) λ (U)I] ∇λ (U) r(U) = 0.
j j
∀ ∀ ∈ − · (cid:54)
Definition E.0.2 (Linearly degenerate fields). The eigenvalue λ (U) is linearly degenerate if
j
U, r ker[A(U) λ (U)I] ∇λ (U) r(U) = 0.
j j
∀ ∀ ∈ − ·
Let us recall the Rankine-Hugoniot conditions
Lemma E.0.1 (Rankine-Hugoniotconditions,[122,Lemma1p3]). Let C be a 1 curve in R2
C
defined by x = ξ(t),ξ 1, that cuts the open set Ω in two open sets Ω := x Ω,x < ξ(t)
∈ C − { ∈ }
and Ω := x Ω,x > ξ(t) . Consider a function U defined on Ω that is of class 1 in Ω
+
{ ∈ } C −
and Ω . Then, U solves (E.1) in the sense of distributions in Ω if and only if U is a classical
+
solution in Ω and Ω , and the Rankine–Hugoniot jump relation
+
−
F(U ) F(U ) = σ(U U ) on C Ω σ = ξ , (E.5)
+ + (cid:48)
− − − − ∩
stands.
As a final reminder, an entropy-flux pair Definition 1.2.2 for (E.1) is given by
u 2 p(s) u 2 p(s)
η := ρ| | +ρ ds, G := ρ| | +ρ ds+p(ρ) u. (E.6)
2 s2 2 s2
(cid:32) (cid:33)
(cid:90) (cid:90)
E.1 The Riemann problem in one dimension
In one space dimension, the Riemann problem is defined by solving (E.1) with an initial con-
ρ ρ
dition composed of two constant states U = L and U = R separated by a
L R
u u
L R
(cid:18) (cid:19) (cid:18) (cid:19)
discontinuity. In this 1D context, the eigenvalues are
λ := u c, λ := u+c.
1 2
−
In the general case, the solution of the Riemann problem is composed of three states
(ρ ,u ),(ρ ,u ),(ρ ,u ) linked by two elementary waves. Those waves are
L L ∗ ∗ R R
• for genuinely non-linear fields Definition E.0.1,
– either rarefaction waves
Definition E.1.1 (Rarefaction wave). We say that two constant states U ,U can
(cid:63) (cid:5)
be joined by a j rarefaction wave if there exists some 1 path U(t) for t [t ,t ],
1 2
− C ∈

APPENDIX E. SOLUTION OF THE 1D RIEMANN PROBLEM IN THE BAROTROPIC
| CASE |     |           |     |     |     |     |     |     |     |     |     |     | 226 |
| ---- | --- | --------- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
|      |     | such that |     |     |     |     |     |     |     |     |     |     |     |
dU(t)
|     |     |             |     |     | =   | r (U(t)), |     | U(t | ) = U ,    | U(t | ) = U . |     | (E.7) |
| --- | --- | ----------- | --- | --- | --- | --------- | --- | --- | ---------- | --- | ------- | --- | ----- |
|     |     |             |     | ds  |     | j         |     |     | 1 (cid:63) | 2   | (cid:5) |     |       |
|     |     | – or shocks |     |     |     |           |     |     |            |     |         |     |       |
DefinitionE.1.2(Shockwave). TheshockwavelinkingU toU isadiscontinuous
(cid:63) (cid:5)
|     |     | auto-similar | field | of  | the form | :   |     |     |          |     |     |     |     |
| --- | --- | ------------ | ----- | --- | -------- | --- | --- | --- | -------- | --- | --- | --- | --- |
|     |     |              |       |     |          |     |     | U   | , x < σt |     |     |     |     |
(cid:63)
|     |     |     |     |     |     | U(x,t) | =   |     |        |     |     |     |     |
| --- | --- | --- | --- | --- | --- | ------ | --- | --- | ------ | --- | --- | --- | --- |
|     |     |     |     |     |     |        |     | U   | x > σt |     |     |     |     |
(cid:5)
(cid:26)
|     |     | with σ satisfying |     | the | Rankine-Hugoniot |     |     | conditions | Lemma | E.0.1. |     |     |     |
| --- | --- | ----------------- | --- | --- | ---------------- | --- | --- | ---------- | ----- | ------ | --- | --- | --- |
•
for linearly degenerate fields Definition E.0.2: contact discontinuities which can be de-
|     |         |       | Rankine-Hugoniot |     |     | conditions |     |     |     |     |     |     |     |
| --- | ------- | ----- | ---------------- | --- | --- | ---------- | --- | --- | --- | --- | --- | --- | --- |
|     | scribed | using |                  |     |     |            |     |     |     |     |     |     |     |
Definition E.1.3 (Contact discontinuity). A j contact discontinuity linking U to U
|     |     |     |     |     |     |     |     |     |     |     |     |     | (cid:63) (cid:5) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | ---------------- |
−
|     | is  | a discontinuous | auto-similar |     |     | field of | the      | form     | :      |     |     |     |     |
| --- | --- | --------------- | ------------ | --- | --- | -------- | -------- | -------- | ------ | --- | --- | --- | --- |
|     |     |                 |              |     |     |          |          | U ,      | x < σt |     |     |     |     |
|     |     |                 |              |     |     | U(x,t)   | =        | (cid:63) |        |     |     |     |     |
|     |     |                 |              |     |     |          |          | U        | x > σt |     |     |     |     |
|     |     |                 |              |     |     |          | (cid:26) | (cid:5)  |        |     |     |     |     |
where
– the wave speed σ, given in (E.5) is equal to the characteristic speed λ (U) and
j
|     |     | – the associated |     | field is  | linearly | degenerate |     | Definition | E.0.2.  |     |         |     |      |
| --- | --- | ---------------- | --- | --------- | -------- | ---------- | --- | ---------- | ------- | --- | ------- | --- | ---- |
|     | For | barotropic       |     | equations |          | however,   |     |            | one can |     | readily | see | that |
contact discontinuities do not exist (all fields are genuinely non-linear) so the element-
ary waves can either be rarefaction waves or shocks, depending on the initial values. We now
| detail | how | their nature | is  | uncovered. |     |     |     |     |     |     |     |     |     |
| ------ | --- | ------------ | --- | ---------- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| E.1.1  |     | Shocks       |     |            |     |     |     |     |     |     |     |     |     |
In the case where the elementary wave linking U to U is a shock, we can derive explicit jump
|                  |     |         |            |         |     |     |     | (cid:63) | (cid:5) |     |     |     |     |
| ---------------- | --- | ------- | ---------- | ------- | --- | --- | --- | -------- | ------- | --- | --- | --- | --- |
| conditions       |     | between | both       | states: |     |     |     |          |         |     |     |     |     |
| Rankine-Hugoniot |     |         | conditions |         |     |     |     |          |         |     |     |     |     |
Following Rankine-Hugoniot conditions Lemma E.0.1, we have for discontinuous solution of
(E.1):
|     |     |     |     |                   | ρ         | u ρ      | u =             | s(ρ      | ρ )                         |         |     |     |       |
| --- | --- | --- | --- | ----------------- | --------- | -------- | --------------- | -------- | --------------------------- | ------- | --- | --- | ----- |
|     |     |     |     |                   | (cid:63)  | (cid:63) | (cid:5) (cid:5) | (cid:63) | (cid:5)                     | ,       |     |     | (E.8) |
|     |     |     |     |                   |           | −        |                 | −        |                             |         |     |     |       |
|     |     |     |     | ρ u2              | ρ (u      | )2+p     | p               | = s(ρ    | u ρ u                       | )       |     |     |       |
|     |     |     |     | (cid:63) (cid:63) | − (cid:5) | (cid:5)  | (cid:63) −      | (cid:5)  | (cid:63) (cid:63) − (cid:5) | (cid:5) |     |     |       |
with s the speed of propagation of the discontinuity. Dividing the first equation of (E.8) by
(ρ (cid:63) ρ ) we obtain an explicit formula for s. Hence plugging this expression in the second
− (cid:5)

APPENDIX E. SOLUTION OF THE 1D RIEMANN PROBLEM IN THE BAROTROPIC
| CASE     |          |           |          |         |     |          |            |           |     | 227   |
| -------- | -------- | --------- | -------- | ------- | --- | -------- | ---------- | --------- | --- | ----- |
| equation | of (E.8) | we obtain |          |         |     |          |            |           |     |       |
|          |          |           |          |         |     |          | 1          | 1         |     |       |
|          |          |           | u        | u =     |     | (p       | p )(       | ).        |     | (E.9) |
|          |          |           | (cid:63) | (cid:5) |     | (cid:63) | (cid:5) ρ  | ρ         |     |       |
|          |          |           |          | −       | ±   | −        | − (cid:63) | − (cid:5) |     |       |
(cid:114)
| Entropy | jump conditions |     |     |     |     |     |     |     |     |     |
| ------- | --------------- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
Given an entropy-flux pair (η,G), the Rankine-Hugoniot conditions Lemma E.0.1 can be up-
dated for an entropy solution to the following weak Rankine-Hugoniot conditions
|     |     | G(U | G(U |     | s(η(U |     | η(U | )), | s ξ ,      |        |
| --- | --- | --- | --- | --- | ----- | --- | --- | --- | ---------- | ------ |
|     |     |     | + ) |     | )     | +   | )   |     | = (cid:48) | (E.10) |
|     |     |     | −   | −   | ≤     |     | −   | −   |            |        |
In order to verify the entropy condition (E.10) with the entropy-flux pair (E.6), only one of the
| signs in | (E.9) is correct: |     |          |           |     |              |           |     |     |        |
| -------- | ----------------- | --- | -------- | --------- | --- | ------------ | --------- | --- | --- | ------ |
|          |                   |     |          |           |     |              | 1         | 1   |     |        |
|          |                   |     | u        | u         | =   | (p           | p )(      | ).  |     | (E.11) |
|          |                   |     | (cid:63) | − (cid:5) |     | − (cid:63) − | (cid:5) ρ | − ρ |     |        |
(cid:63) (cid:5)
(cid:114)
Note also that the square root is supposed real (otherwise the velocity would be imaginary), so
| by monotony | of p(ρ)     | (E.11) | makes               | sense | only | p          | > p .    |             |     |     |
| ----------- | ----------- | ------ | ------------------- | ----- | ---- | ---------- | -------- | ----------- | --- | --- |
|             |             |        |                     |       |      | (cid:5)    | (cid:63) |             |     |     |
| E.1.2       | Rarefaction |        | waves               |       |      |            |          |             |     |     |
| Rarefaction | waves       | can be | fully characterized |       |      | by Riemann |          | invariants. |     |     |
Definition E.1.4 (Riemann invariants). A j Riemann invariant is a function ω(U) that
−
verifies
|     |     |     | r ker[A(U) |     |     | λ (U)I] | ∇ω  | r   | = 0. |     |
| --- | --- | --- | ---------- | --- | --- | ------- | --- | --- | ---- | --- |
j
|           |        |        | ∀ ∈ |         | −          |     |           | ·   |     |     |
| --------- | ------ | ------ | --- | ------- | ---------- | --- | --------- | --- | --- | --- |
| For Euler | system | (E.1), | the | Riemann | invariants |     | are given | by  |     |     |
c(s)
|     |     |     |     | ω(U) | :=  | u   | ds. |     |     |     |
| --- | --- | --- | --- | ---- | --- | --- | --- | --- | --- | --- |
|     |     |     |     |      |     | ±   | s   |     |     |     |
(cid:90)
Through the j rarefaction (E.7), the j Riemann invariant is constant, yielding for a
|               | −           |      |          |           |      | −      |          |      |     |        |
| ------------- | ----------- | ---- | -------- | --------- | ---- | ------ | -------- | ---- | --- | ------ |
| 1 rarefaction | wave        |      |          |           |      |        |          |      |     |        |
| −             |             |      |          | ρ(cid:63) |      |        | ρ(cid:5) |      |     |        |
|               |             |      |          |           | c(s) |        |          | c(s) |     |        |
|               |             |      | u +      |           |      | ds = u | +        | ds,  |     | (E.12) |
|               |             |      | (cid:63) |           | s    |        | (cid:5)  | s    |     |        |
|               |             |      |          | (cid:90)  |      |        | (cid:90) |      |     |        |
| and a 2       | rarefaction | wave |          |           |      |        |          |      |     |        |
−
|       |            |     |          | ρ(cid:63) |      |        | ρ(cid:5) |      |     |        |
| ----- | ---------- | --- | -------- | --------- | ---- | ------ | -------- | ---- | --- | ------ |
|       |            |     |          |           | c(s) |        |          | c(s) |     |        |
|       |            |     | u        |           |      | ds = u |          | ds.  |     | (E.13) |
|       |            |     | (cid:63) |           |      |        | (cid:5)  |      |     |        |
|       |            |     | −        |           | s    |        | −        | s    |     |        |
|       |            |     |          | (cid:90)  |      |        | (cid:90) |      |     |        |
| E.1.3 | The Newton |     | solver   |           |      |        |          |      |     |        |
Applying (E.11) and (E.12) to U = U , U = U we get the following relation
|     |     |     |     | (cid:63) | L   | (cid:5) | ∗     |     |     |        |
| --- | --- | --- | --- | -------- | --- | ------- | ----- | --- | --- | ------ |
|     |     |     |     | u        | u   | = f(p   | ,W ), |     |     | (E.14) |
|     |     |     |     |          | L   | ∗       | ∗ L   |     |     |        |
−

APPENDIX E. SOLUTION OF THE 1D RIEMANN PROBLEM IN THE BAROTROPIC
CASE 228
with
1 1
(p p)( ) if p > p
∗
◦− ρ − ρ ◦
f(p,W ) :=  (cid:114) ◦ . (E.15)

◦    ρ c(s) ρ◦ c(s)

ds ds else
s − s
 (cid:90) (cid:90)

Similarly, applying (E.11) and (E .13) to U = U , U = U we get
 (cid:63) ∗ (cid:5) R
u u = f(p ,W ), (E.16)
∗ R ∗ R
−
with f defined as in (E.15). Then summing (E.14) and (E.16) yields
u u = f(p ,W )+f(p ,W ). . (E.17)
R L ∗ R ∗ L
−
Also, substracting (E.14) to (E.16) gives
u +u f(p ,W ) f(p ,W )
u ∗ = L R + ∗ R − ∗ L .
2 2
It all boils down to finding the root p of equation (E.17) which can be solved thanks to a
∗
Newton scheme. As a consequence we need the explicit formula of ∂ f(p,W). For example,
p
p 1
when p(ρ) = Aργ = we have ρ = ( )γ, c(ρ) = Aγργ 1:
⇒ A −
(cid:112)
1 1
∂ (p p)( ) if p > p
p − L − p L 1 − p 1 L
 (cid:34)(cid:118) ( )γ ( )γ (cid:35)
∂ p f(p,W L ) :=    (cid:117) (cid:117) A A
   (cid:116)2 ∂ p γA( p ) γ− γ 1 c(ρ L ) else
γ 1 A −
 −
(cid:34)(cid:32)(cid:114) (cid:33)(cid:35)



now  
1 1
∂ (p p)( )
p − L − p L 1 − p 1
(cid:34) ( )γ ( )γ) (cid:35)
1 1 A A
∂ (p p)( ) =
p − L − p L 1 − p 1 2f(p,W)
(cid:34)(cid:118) ( )γ ( )γ (cid:35)
(cid:117) A A
(cid:117)
(cid:116)
1 1 1 1
(p p)∂ ( )+( )
− L − p p L 1 − p 1 p L 1 − p 1 ,
( )γ ( )γ) ( )γ ( )γ)
A A A A
=
2f(p,W)
− 1 1 (p L − p)p− ( γ 1+1) +( p L 1 1 − p 1 1 )
γAγ ( )γ ( )γ)
A A
2f(p,W)

APPENDIX E. SOLUTION OF THE 1D RIEMANN PROBLEM IN THE BAROTROPIC
CASE 229
and
1
|       |       | 2                         | p γ−1   |                    |     | 2 √γA |              | γ−1      | A2γ γ+1 |
| ----- | ----- | ------------------------- | ------- | ------------------ | --- | ----- | ------------ | -------- | ------- |
|       |       | ∂                         | γA(     | c(ρ                |     |       | ∂ p          |          | p− ,    |
|       |       | p                         | ) γ     | L )                | =   |       | p            | 2γ =     | 2γ      |
|       | γ     | 1                         | A       | −                  | γ   | 1     | γ − 1        |          | 1       |
|       |       | (cid:34)(cid:32)(cid:114) |         | (cid:33)(cid:35)   |     | A     | 2 γ (cid:34) | (cid:35) | γ2      |
|       |       | −                         |         |                    |     | −     |              |          |         |
| which | gives | all the necessary         | formula | for the resolution |     | of    | (E.17).      |          |         |

| Appendix    |     |                | F   |     |     |     |     |            |           |     |     |      |
| ----------- | --- | -------------- | --- | --- | --- | --- | --- | ---------- | --------- | --- | --- | ---- |
| Computation |     |                |     |     | of  | the |     | convection |           |     |     | term |
| for         | the | Raviart-Thomas |     |     |     |     |     |            | staggered |     |     |      |
scheme
Using finite element spaces has been very convenient in order to prove long time consistency
of the wave system. Nonetheless, once we treat the Euler equations, this formalism implies the
| computation |     | of the term |     |     |     |     |     |     |     |     |     |     |
| ----------- | --- | ----------- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
|             |     |             | q   | q   |     |     |     |     |     | q   |     |     |
+I
|     |              |                    | ⊗   |     | d p      | : ∇Ψ σ | dx+      | q          | n                | +np                      | Ψ   | σ dΓ     |
| --- | ------------ | ------------------ | --- | --- | -------- | ------ | -------- | ---------- | ---------------- | ------------------------ | --- | -------- |
|     |              | −                  |     | ρ   |          |        |          |            | ·                | ρ                        | ·   |          |
|     |              | (cid:34) K(cid:18) |     |     | (cid:19) |        |          | ∂K(cid:18) | (cid:26)(cid:26) | (cid:27)(cid:27)(cid:19) |     | (cid:35) |
|     | K (cid:88)∈C | (σ) (cid:90)       |     |     |          |        | (cid:90) |            |                  |                          |     |          |
On affine meshes, even if the coefficients are not analytically known a priori as in [57], it causes
little to no problem since the gradient of the basis functions is constant. By opposition, on
meshes where cells lead to a non-constant transformation’s gradient, the Piola transform will
lead to new difficulties. Indeed using the contravariant Piola transform, we define the basis
K
| functions | on  | the physical | element |     | with | the | following | formula: |     |     |     |     |
| --------- | --- | ------------ | ------- | --- | ---- | --- | --------- | -------- | --- | --- | --- | --- |
σ
|     |     |     |     | (x) |     |             | B(x)Ψ |     | (T 1(x)) |     |     |     |
| --- | --- | --- | --- | --- | --- | ----------- | ----- | --- | -------- | --- | --- | --- |
|     |     |     |     | Ψ σ | =   | | |         |       | σˆ  | K−       |     |     |     |
|     |     |     |     |     |     | det( B (x)) |       |     |          |     |     |     |
|     |     |     |     |     | |   |             | |     |     |          |     |     |     |
(cid:98)
where
|         |           |       |           |     | B(x)         | := [∇ | T ](T  | 1(x)) |     |     |     |     |
| ------- | --------- | ----- | --------- | --- | ------------ | ----- | ------ | ----- | --- | --- | --- | --- |
|         |           |       |           |     |              |       | xˆ K   | K−    |     |     |     |     |
| So, for | a general | basis | function, |     | its gradient | is    | of the | form  |     |     |     |     |
B(x)
|     |     |     |     |       |      |            |           | Ψˆ  | 1(xˆ)    |          |     |     |
| --- | --- | --- | --- | ----- | ---- | ---------- | --------- | --- | -------- | -------- | --- | --- |
|     |     |     | ∇   | x Ψ σ | = ∇x | σ          |           |     | σˆ (T K− |          |     |     |
|     |     |     |     |       |      | || |       | det(B(x)) |     |          |          |     |     |
|     |     |     |     |       |      | (cid:18) | |           | |   |          | (cid:19) |     |     |
For a generic quadrangular cell, the gradient of the basis functions are especially tedious to
integrate, and in particular, a classical mid-point formula on the cell is insufficient to preserve
q
simple stationary states such as constant fields. To tell the truth, if is a constant vector,
231

APPENDIX F. COMPUTATION OF THE CONVECTION TERM FOR THE
RAVIART-THOMAS STAGGERED SCHEME 232
ρ = 1, the challenge is to integrate
G : ∇Ψ dx
σ
K
(cid:90)
with G a constant matrix. But
1 1 1
∇ ( BΨˆ (T 1(x))) = ∂ BΨˆ (T 1(x)) + ∂ BΨˆ (T 1(x))
x
detB
σ K− xj
detB
σ K−
detB
xj σ K−
(cid:20) (cid:21)i,j (cid:20) (cid:21)(cid:18) (cid:19)i (cid:20)(cid:18) (cid:19)i(cid:21)
2
1 1
= ∂ (detB) BΨˆ (T 1(x)) + ∂ B Ψˆ (T 1(x))
− detB
xj σ K−
detB
xj i,k σ K− k
(cid:18) (cid:19) (cid:18) (cid:19)i (cid:20) (cid:21)
2
1 1 1
= ∂ (detB) BΨˆ (T 1(x)) + ∂ B Ψˆ (T 1(x)) + B∇ Ψˆ (T 1(x))
− detB
xj σ K−
detB
xj i,k σ K− k
detB
x σ K−
(cid:18) (cid:19) (cid:18) (cid:19)i (cid:18) (cid:19)i,j
whereweuseEinsteinconventionontherepeatedindices. Theexactexpressionsofthetermsis
not releveant for this discussion; we can illustrate the problem intrinsic to the Piola transform
with the last term. The computation of the integral is problematic mainly because of terms in
B 1,
−
1 1
GG : B(x)∇ Ψˆ (T 1(x))dx = G : B(x)∇ Ψˆ (T 1(x))B 1(x)dx,
det(B(x))
x σ K−
det(B(x))
xˆ σ K− −
K K
(cid:90) | | (cid:90) | |
now by definition det(B(x))B 1(x) = Com[B(x)] where Com[B(x)] is the comatrix of B(x),
−
so that applying the change of variable x = T (xˆ) will yield:
K
1 1
G : B(x)∇ Ψˆ (T 1(x))dx = G : B(x)∇ Ψˆ (T 1(x))Com[B(x)]dxˆ.
K
det(B(x))
x σ K−
Kˆ
det(B(x))
xˆ σ K−
(cid:90) | | (cid:90)
Onaffinemeshesdet(B(x)) = det(B)isnecessarilyconstantinthecell,onagenericquadrangle
however it is not, so the computation boils down to integrating a rational function instead of
classically integrating a polynomial. Integrating on the reference element with the trapezoid
formula enabled us to recover the preservation of constant state for the meshes we used in
application.

Bibliography
[1] Neil E Todreas and Mujid S Kazimi. Nuclear systems volume I: Thermal hydraulic fun-
| damentals. | CRC press, | 2021. |     |     |     |     |
| ---------- | ---------- | ----- | --- | --- | --- | --- |
[2] Edwige Godlewski and Pierre-Arnaud Raviart. Numerical approximation of hyperbolic
systems of conservation laws, volume 118. Springer Science & Business Media, 2013.
[3] Bruno Despr´es. Numerical methods for Eulerian and Lagrangian conservation laws.
| Birkh¨auser, | 2017. |     |     |     |     |     |
| ------------ | ----- | --- | --- | --- | --- | --- |
[4] Sergiu Klainerman and Andrew Majda. Singular limits of quasilinear hyperbolic systems
with large parameters and the incompressible limit of compressible fluids. Communica-
| tions on | pure and | applied Mathematics, | 34(4):481–524, |     | 1981. |     |
| -------- | -------- | -------------------- | -------------- | --- | ----- | --- |
[5] Steven Schochet. Fast singular limits of hyperbolic PDEs. Journal of differential equa-
tions,
| 114(2):476–512, |     | 1994. |     |     |     |     |
| --------------- | --- | ----- | --- | --- | --- | --- |
[6] Emmanuel Grenier. Oscillatory perturbations of the Navier-Stokes equations. Journal de
| Math´ematiques | Pures | et Appliqu´ees, | 76(6):477–498, |     | 1997. |     |
| -------------- | ----- | --------------- | -------------- | --- | ----- | --- |
[7] P-L Lions and Nader Masmoudi. Incompressible limit for a viscous compressible fluid.
| Journal | de math´ematiques | pures | et appliqu´ees, | 77(6):585–627, |     | 1998. |
| ------- | ----------------- | ----- | --------------- | -------------- | --- | ----- |
[8] Jonathan Jung and Vincent Perrier. Behavior of the discontinuous Galerkin method for
compressible flows at low Mach number on triangles and tetrahedrons. SIAM Journal on
| Scientific | Computing, | 46(1):A452–A482, | 2024. |     |     |     |
| ---------- | ---------- | ---------------- | ----- | --- | --- | --- |
[9] Yousef Saad. Iterative methods for sparse linear systems. SIAM, 2003.
[10] CliftonWall,CharlesDPierce,andParvizMoin. Asemi-implicitmethodforresolutionof
acousticwavesinlowMachnumberflows. Journal of Computational Physics,181(2):545–
563, 2002.
[11] Giacomo Dimarco, Rapha¨el Loub`ere, and Marie-H´el`ene Vignal. Study of a new asymp-
SIAM journal
toticpreservingschemefortheEulersysteminthelowMachnumberlimit.
| on Scientific | Computing, | 39(5):A2099–A2128, |     | 2017. |     |     |
| ------------- | ---------- | ------------------ | --- | ----- | --- | --- |
[12] Philip L Roe. Approximate Riemann solvers, parameter vectors, and difference schemes.
| Journal | of computational | physics, | 43(2):357–372, | 1981. |     |     |
| ------- | ---------------- | -------- | -------------- | ----- | --- | --- |
234

BIBLIOGRAPHY 235
[13] Amiram Harten, Peter D Lax, and Bram van Leer. On upstream differencing and
Godunov-typeschemesforhyperbolicconservationlaws. SIAM review,25(1):35–61,1983.
[14] EleuterioFToro. Riemann solvers and numerical methods for fluid dynamics: a practical
introduction.
|     | Springer |     | Science & | Business | Media, | 2013. |
| --- | -------- | --- | --------- | -------- | ------ | ----- |
[15] Francis H Harlow. MAC numerical calculation of time-dependent viscous incompressible
| flow of fluid | with | free | surface. Phys. | Fluid, | 8:12, | 1965. |
| ------------- | ---- | ---- | -------------- | ------ | ----- | ----- |
[16] J Eddie Welch, Francis Harvey Harlow, John P Shannon, and Bart J Daly. The MAC
method-a computing technique for solving viscous, incompressible, transient fluid-flow
problems involving free surfaces. Technical report, Los Alamos National Lab.(LANL),
| Los Alamos, | NM  | (United | States), | 1965. |     |     |
| ----------- | --- | ------- | -------- | ----- | --- | --- |
[17] Randall J LeVeque. Numerical methods for conservation laws, volume 132. Springer,
1992.
Finitevolumemethodsforhyperbolicproblems,volume31.
[18] RandallJLeVeque. Cambridge
| university | press, | 2002. |     |     |     |     |
| ---------- | ------ | ----- | --- | --- | --- | --- |
[19] G Volpe. Performance of compressible flow codes at low Mach numbers. AIAA journal,
| 31(1):49–56, | 1993. |     |     |     |     |     |
| ------------ | ----- | --- | --- | --- | --- | --- |
[20] Herv´e Guillard and C´ecile Viozat. On the behaviour of upwind schemes in the low Mach
| number limit. | Computers |     | & fluids, | 28(1):63–86, |     | 1999. |
| ------------- | --------- | --- | --------- | ------------ | --- | ----- |
[21] Herv´eGuillardandAngeloMurrone. OnthebehaviorofupwindschemesinthelowMach
number limit: II. Godunov type schemes. Computers & fluids, 33(4):655–675, 2004.
[22] Herv´eGuillard. OnthebehaviorofupwindschemesinthelowMachnumberlimit.IV:P0
approximationontriangularandtetrahedralcells. Computers & fluids,38(10):1969–1972,
2009.
[23] Herv´e Guillard and Boniface Nkonga. On the behaviour of upwind schemes in the low
Mach number limit: A review. Handbook of Numerical Analysis, 18:203–231, 2017.
[24] KeithWilliamMortonandPhilipLRoe. Vorticity-preservingLaxWendroff-typeschemes
for the system wave equation. SIAM Journal on Scientific Computing, 23(1):170–192,
2001.
[25] St´ephane Dellacherie. Analysis of Godunov type schemes applied to the compressible
Euler system at low Mach number. Journal of Computational Physics, 229(4):978–1016,
2010.
Annual Review
[26] Eli Turkel. Preconditioning techniques in computational fluid dynamics.
| of Fluid | Mechanics, | 31(1):385–416, |     | 1999. |     |     |
| -------- | ---------- | -------------- | --- | ----- | --- | --- |
[27] E Turkel, A Fiterman, and Bram Van Leer. Preconditioning and the limit to the incom-
| pressible | flow equations. |     | Technical | report, | 1993. |     |
| --------- | --------------- | --- | --------- | ------- | ----- | --- |

BIBLIOGRAPHY 236
[28] Xue-songLiandChun-weiGu. Anall-speedRoe-typeschemeanditsasymptoticanalysis
of low Mach number behaviour. Journal of Computational Physics, 227(10):5144–5159,
2008.
[29] Xue-song Li, Chun-wei Gu, and Jian-zhong Xu. Development of Roe-type scheme for
all-speed flows based on preconditioning method. Computers & Fluids, 38(4):810–817,
2009.
[30] Felix Rieper. A low-Mach number fix for Roe’s approximate Riemann solver. Journal of
Computational Physics, 230(13):5263–5287, 2011.
[31] Kai Oßwald, Alexander Siegmund, Philipp Birken, Volker Hannemann, and Andreas
Meister. L2 Roe: a low dissipation version of Roe’s approximate Riemann solver for
low Mach numbers. International Journal for Numerical Methods in Fluids, 81(2):71–86,
2016.
[32] Christophe Chalons, Mathieu Girardin, and Samuel Kokh. An all-regime Lagrange-
Projection like scheme for the gas dynamics equations on unstructured meshes. Commu-
nications in Computational Physics, 20(1):188–233, 2016.
[33] Floraine Cordier, Pierre Degond, and Anela Kumbaro. An asymptotic-preserving all-
speed scheme for the Euler and Navier–Stokes equations. Journal of Computational
Physics, 231(17):5685–5704, 2012.
[34] Pierre Degond, Shi Jin, and J Yuming. Mach-number uniform asymptotic-preserving
gaugeschemesforcompressibleflows. Bulletin-Institute of Mathematics Academia Sinica,
2(4):851, 2007.
[35] Pierre Degond and Min Tang. All speed scheme for the low Mach number limit of the
isentropic Euler equations. Communications in Computational Physics, 10(1):1–31, 2011.
[36] Jeffrey Haack, Shi Jin, and Jian-Guo Liu. An all-speed asymptotic-preserving method
for the isentropic Euler and Navier-Stokes equations. Communications in Computational
Physics, 12(4):955–980, 2012.
[37] Klaus Kaiser, Jochen Schu¨tz, Ruth Sch¨obel, and Sebastian Noelle. A new stable splitting
for the isentropic Euler equations. Journal of scientific computing, 70:1390–1407, 2017.
[38] St´ephane Dellacherie. Checkerboard modes and wave equation. In Proceedings of AL-
GORITMY, volume 2009, pages 71–80, 2009.
[39] Ibtissem Lannabi. Analysis of spurious oscillations problem of Finite Volume Methods
for low Mach number flows in fluid mechanics. PhD thesis, Universit´e de Pau et des Pays
de l’Adour, 2024.
[40] Jonathan Jung, Ibtissem Lannabi, and Vincent Perrier. On the convergence of the
Godunov scheme with a centered discretization of the pressure gradient. In Interna-
tional Conference on Finite Volumes for Complex Applications, pages 201–208. Springer,
2023.

BIBLIOGRAPHY 237
[41] Felix Rieper. Influence of cell geometry on the behaviour of the first-order Roe scheme
in the low Mach number regime. Finite Volumes for Complex Applications V, pages
| 625–632, | 2008. |     |     |     |     |     |     |     |
| -------- | ----- | --- | --- | --- | --- | --- | --- | --- |
[42] St´ephane Dellacherie, Pascal Omnes, and Felix Rieper. The influence of cell geometry
Journal of Computational
| on the   | Godunov            | scheme | applied |       | to the linear | wave | equation. |     |
| -------- | ------------------ | ------ | ------- | ----- | ------------- | ---- | --------- | --- |
| Physics, | 229(14):5315–5338, |        |         | 2010. |               |      |           |     |
[43] Jonathan Jung and Vincent Perrier. Steady low Mach number flows: identification of
the spurious mode and filtering method. Journal of Computational Physics, 468:111462,
2022.
[44] MichelCrouzeixandP-ARaviart. Conformingandnonconformingfiniteelementmethods
forsolvingthestationaryStokesequationsI. Revue fran¸caise d’automatique informatique
| recherche | op´erationnelle. |     | Math´ematique, |     | 7(R3):33–75, |     | 1973. |     |
| --------- | ---------------- | --- | -------------- | --- | ------------ | --- | ----- | --- |
[45] Jonathan Jung and Vincent Perrier. A curl preserving finite volume scheme by space
Journal of
velocity enrichment. application to the low Mach number accuracy problem.
| Computational |     | Physics, | 515:113252, |     | 2024. |     |     |     |
| ------------- | --- | -------- | ----------- | --- | ----- | --- | --- | --- |
[46] Roy A Nicolaides. Analysis and convergence of the mac scheme. I. the linear problem.
| SIAM Journal |     | on Numerical |     | Analysis, | 29(6):1579–1591, |     | 1992. |     |
| ------------ | --- | ------------ | --- | --------- | ---------------- | --- | ----- | --- |
[47] Roy A. Nicolaides and X. Wu. Analysis and convergence of the MAC scheme. II. Navier-
| Stokes | equations. | Mathematics |     | of  | Computation, |     | 65(213):29–44, | 1996. |
| ------ | ---------- | ----------- | --- | --- | ------------ | --- | -------------- | ----- |
[48] Vincenzo Casulli and D Greenspan. Pressure method for the numerical solution of tran-
sient, compressible fluid flows. International Journal for Numerical Methods in Fluids,
| 4(11):1001–1012, |     | 1984. |     |     |     |     |     |     |
| ---------------- | --- | ----- | --- | --- | --- | --- | --- | --- |
[49] RaadIIssa, ADGosman, andAPWatkins. Thecomputationofcompressibleandincom-
pressible recirculating flows by a non-iterative implicit scheme. Journal of Computational
| Physics, | 62(1):66–82, |     | 1986. |     |     |     |     |     |
| -------- | ------------ | --- | ----- | --- | --- | --- | --- | --- |
[50] Robert Eymard, Thierry Gallou¨et, Raphaele Herbin, and Jean-Claude Latch´e. Con-
SIAM Journal on
| vergence  | of the    | MAC              | scheme | for | the compressible |     | Stokes equations. |     |
| --------- | --------- | ---------------- | ------ | --- | ---------------- | --- | ----------------- | --- |
| Numerical | Analysis, | 48(6):2218–2246, |        |     | 2010.            |     |                   |     |
[51] ThierryGallou¨et,RaphaeleHerbin,andJean-ClaudeLatch´e. Aconvergentfiniteelement-
finite volume scheme for the compressible Stokes problem. part I: The isothermal case.
| Mathematics |     | of Computation, |     | 78(267):1333–1352, |     |     | 2009. |     |
| ----------- | --- | --------------- | --- | ------------------ | --- | --- | ----- | --- |
[52] V Selmin and Luca Formaggia. Unified construction of finite element and finite volume
discretizations for compressible flows. International Journal for Numerical Methods in
| Engineering, |     | 39(1):1–32, | 1996. |     |     |     |     |     |
| ------------ | --- | ----------- | ----- | --- | --- | --- | --- | --- |

| BIBLIOGRAPHY |     |     |     |     |     |     |     | 238 |
| ------------ | --- | --- | --- | --- | --- | --- | --- | --- |
[53] Jacques Baranger, Jean-Fran¸cois Maitre, and Fabienne Oudin. Connection between fi-
nite volume and mixed finite element methods. ESAIM: Mathematical Modelling and
| Numerical | Analysis, | 30(4):445–465, |     | 1996. |     |     |     |     |
| --------- | --------- | -------------- | --- | ----- | --- | --- | --- | --- |
[54] Franc¸ois Dubois, Isabelle Greff, and Charles Pierre. Raviart–thomas finite elements
|                    |     |       | ESAIM: |     | Mathematical | Modelling | and Numerical | Analysis, |
| ------------------ | --- | ----- | ------ | --- | ------------ | --------- | ------------- | --------- |
| of petrov–Galerkin |     | type. |        |     |              |           |               |           |
| 53(5):1553–1576,   |     | 2019. |        |     |              |           |               |           |
[55] Thierry Gallou¨et, Laura Gastaldo, Raphaele Herbin, and Jean-Claude Latch´e. An un-
conditionally stable pressure correction schemefor the compressible barotropic Navier-
Stokes equations. ESAIM: Mathematical Modelling and Numerical Analysis, 42(2):303–
331, 2008.
[56] Laura Gastaldo, Rapha`ele Herbin, and Jean-Claude Latch´e. A discretization of the phase
mass balance in fractional step algorithms for the drift-flux model. IMA Journal of
| Numerical | Analysis, | 31(1):116–146, |     | 2011. |     |     |     |     |
| --------- | --------- | -------------- | --- | ----- | --- | --- | --- | --- |
[57] G Ansanay-Alex, F Babik, JC Latch´e, and D Vola. An L2-stable approximation of the
Navier–Stokes convection operator for low-order non-conforming finite elements. Inter-
national Journal for Numerical Methods in Fluids, 66(5):555–580, 2011.
[58] Katsushi Ohmori and Teruo Ushijima. A technique of upstream type applied to a linear
nonconforming finite element approximation of convective diffusion equations. RAIRO.
| Analyse | num´erique, |                |     |       |     |     |     |     |
| ------- | ----------- | -------------- | --- | ----- | --- | --- | --- | --- |
|         |             | 18(3):309–332, |     | 1984. |     |     |     |     |
[59] PhilippeAngot, V´ıtDolejˇs´ı, MiloslavFeistauer, andJiˇr´ıFelcman. Analysisofacombined
barycentricfinitevolume—nonconformingfiniteelementmethodfornonlinearconvection-
| diffusion | problems. | Applications |     | of Mathematics, |     | 43(4):263–310, | 1998. |     |
| --------- | --------- | ------------ | --- | --------------- | --- | -------------- | ----- | --- |
[60] Robert Eymard, Danielle Hilhorst, and Martin Vohral´ık. A combined finite volume–
nonconforming/mixed-hybrid finite element scheme for degenerate parabolic problems.
| Numerische | Mathematik, |     | 105(1):73–131, |     | 2006. |     |     |     |
| ---------- | ----------- | --- | -------------- | --- | ----- | --- | --- | --- |
[61] F Schieweck and L Tobiska. A nonconforming finite element method of upstream type
applied to the stationary Navier-Stokes equation. ESAIM: Mathematical Modelling and
| Numerical | Analysis, |                |     |       |     |     |     |     |
| --------- | --------- | -------------- | --- | ----- | --- | --- | --- | --- |
|           |           | 23(4):627–647, |     | 1989. |     |     |     |     |
[62] F Schieweck and L Tobiska. An optimal order error estimate for an upwind discretization
of the Navier-Stokes equations. Numerical Methods for Partial Differential Equations:
| An International |     | Journal, | 12(4):407–421, |     | 1996. |     |     |     |
| ---------------- | --- | -------- | -------------- | --- | ----- | --- | --- | --- |
[63] Rolf Rannacher and Stefan Turek. Simple nonconforming quadrilateral Stokes element.
Numerical Methods for Partial Differential Equations, 8(2):97–111, 1992.
[64] Stefan Turek. Efficient solvers for incompressible flow problems: An algorithmic and
computational approache, volume 6. Springer Science & Business Media, 1999.

| BIBLIOGRAPHY |     |     |     |     |     | 239 |
| ------------ | --- | --- | --- | --- | --- | --- |
[65] Raphaele Herbin, Jean-Claude Latch´e, and Trung Tan Nguyen. Consistent explicit
staggered schemes for compressible flows part I: the barotropic Euler equations. 2013.
[66] Raphaele Herbin, Walid Kheriji, and J-C Latch´e. On some implicit and semi-implicit
staggered schemes for the shallow water and Euler equations. ESAIM: Mathematical
| Modelling | and Numerical | Analysis, |     |     |     |     |
| --------- | ------------- | --------- | --- | --- | --- | --- |
48(6):1807–1857, 2014.
[67] Dionysis Grapsas, Rapha`ele Herbin, Walid Kheriji, and Jean-Claude Latch´e. An uncon-
ditionally stable staggered pressure correction scheme for the compressible Navier-Stokes
equations. The SMAI journal of computational mathematics, 2:51–97, 2016.
[68] Aubin Brunel, Rapha`ele Herbin, and Jean-Claude Latch´e. A staggered scheme for the
compressible Euler equations on general 3d meshes. arXiv preprint arXiv:2209.06474,
2022.
[69] Arnaud Duran, Jean-Paul Vila, and R´emy Baraille. Energy-stable staggered schemes for
the shallow water equations. Journal of computational Physics, 401:109051, 2020.
[70] Thierry Goudon, Julie Llobell, and Sebastian Minjeaud. A staggered scheme for the
Euler equations. In Finite Volumes for Complex Applications VIII-Hyperbolic, Elliptic
and Parabolic Problems: FVCA 8, Lille, France, June 2017 8, pages 91–99. Springer,
2017.
|     |     |     |     |     | SIAM Journal | on  |
| --- | --- | --- | --- | --- | ------------ | --- |
[71] Roy A Nicolaides. Direct discretization of planar div-curl problems.
| Numerical | Analysis, | 29(1):32–56, | 1992. |     |     |     |
| --------- | --------- | ------------ | ----- | --- | --- | --- |
[72] BernhardMu¨ller. Low-Mach-numberasymptoticsoftheNavier-Stokesequations. Journal
| of Engineering | Mathematics, | 34(1):97–109, |     | 1998. |     |     |
| -------------- | ------------ | ------------- | --- | ----- | --- | --- |
[73] Thomas Gali´e, Jonathan Jung, Ibtissem Lannabi, and Vincent Perrier. Extension of an
all-Mach Roe scheme able to deal with low Mach acoustics to full Euler system. ESAIM:
| Proceedings | and Surveys, | 76:35–51, | 2024. |     |     |     |
| ----------- | ------------ | --------- | ----- | --- | --- | --- |
[74] St´ephane Dellacherie, Jonathan Jung, Pascal Omnes, and P-A Raviart. Construction of
modifiedGodunov-typeschemesaccurateatanyMachnumberforthecompressibleEuler
|         | Mathematical | Models | and Methods | in Applied Sciences, |                   |       |
| ------- | ------------ | ------ | ----------- | -------------------- | ----------------- | ----- |
| system. |              |        |             |                      | 26(13):2525–2615, | 2016. |
[75] PascalBruel, SimonDelmas, JonathanJung, andVincentPerrier. AlowMachcorrection
able to deal with low Mach acoustics. Journal of Computational Physics, 378:723–759,
2019.
[76] Vyacheslav Ivanovich Lebedev. Difference analogues of orthogonal decompositions, basic
differential operators and some boundary problems of mathematical physics. I. USSR
Computational Mathematics and Mathematical Physics, 4(3):69–92, 1964.
[77] RobertD.RichtmyerandE.H.Dill. Differencemethodsforinitial-valueproblems. Physics
| Today, | 12(4):50–50, | 1959. |     |     |     |     |
| ------ | ------------ | ----- | --- | --- | --- | --- |

BIBLIOGRAPHY 240
[78] AmyLBauer,Rapha¨elLoubere,andBurtonWendroff. Onstabilityofstaggeredschemes.
| SIAM | Journal | on Numerical |     | Analysis, | 46(2):996–1011, |     | 2008. |
| ---- | ------- | ------------ | --- | --------- | --------------- | --- | ----- |
[79] Robert D. Richtmyer and K. W. Morton. Difference Methods for initial-value problems,
volume4ofInterscienceTractsinPureandAppliedMathematics. Intersciencepublishers,
| John | Wiley | & Sons, | 1967. |     |     |     |     |
| ---- | ----- | ------- | ----- | --- | --- | --- | --- |
[80] Peter D. Lax and Robert D. Richtmyer. Survey of the stability of linear finite difference
equations. In Selected Papers Volume I, pages 125–151. Springer, 2005.
[81] Rapha`ele Herbin, J-C Latch´e, and Khaled Saleh. Low Mach number limit of some
Mathematics of Computation,
| staggered          | schemes |     | for compressible |     | barotropic | flows. |     |
| ------------------ | ------- | --- | ---------------- | --- | ---------- | ------ | --- |
| 90(329):1039–1087, |         |     | 2021.            |     |            |        |     |
[82] Martin Parisot and Jean-Paul Vila. Centered-potential regularization for the advection
upstream splitting method. SIAM J. Numer. Anal., 54(5):3083–3104, 2016.
[83] Vincent Perrier. discrete de Rham complex involving a discontinuous finite element space
forvelocities: thecaseofperiodicstraighttriangularandCartesianmeshes.arXivpreprint
| arXiv:2404.19545, |     | 2024.   |        |         |                    |     |             |
| ----------------- | --- | ------- | ------ | ------- | ------------------ | --- | ----------- |
| [84] Douglas      | N   | Arnold. | Finite | element | exterior calculus. |     | SIAM, 2018. |
[85] Douglas Arnold, Richard Falk, and Ragnar Winther. Finite element exterior calculus:
from Hodge theory to numerical stability. Bulletin of the American mathematical society,
| 47(2):281–354, |     | 2010. |     |     |     |     |     |
| -------------- | --- | ----- | --- | --- | --- | --- | --- |
[86] Martin Costabel and Alan McIntosh. On Bogovski˘ı and regularized Poincar´e integral
operators for de Rham complexes on Lipschitz domains. Mathematische Zeitschrift,
| 265(2):297–320, |     | 2010. |     |     |     |     |     |
| --------------- | --- | ----- | --- | --- | --- | --- | --- |
[87] Vivette Girault and Pierre-Arnaud Raviart. Finite element approximation of the Navier-
| Stokes | equations. | Springer, |     | 1979. |     |     |     |
| ------ | ---------- | --------- | --- | ----- | --- | --- | --- |
Archive for
[88] Constantine M Dafermos. Quasilinear hyperbolic systems with involutions.
| Rational | Mechanics |     | and Analysis, |     | 94(4):373–389, | 1986. |     |
| -------- | --------- | --- | ------------- | --- | -------------- | ----- | --- |
[89] Konstantin Lipnikov, Gianmarco Manzini, and Mikhail Shashkov. Mimetic finite differ-
ence method. Journal of Computational Physics, 257:1163–1227, 2014.
[90] Len G Margolin, Mikhail Shashkov, and Piotr K Smolarkiewicz. A discrete operator
calculus for finite difference approximations. Computer methods in applied mechanics
| and engineering, |     | 187(3-4):365–383, |     |     | 2000. |     |     |
| ---------------- | --- | ----------------- | --- | --- | ----- | --- | --- |
[91] JamesMHymanandMikhailShashkov. Naturaldiscretizationsforthedivergence, gradi-
ent, andcurlonlogicallyrectangulargrids. Computers & Mathematics with Applications,
| 33(4):81–104, |     | 1997. |     |     |     |     |     |
| ------------- | --- | ----- | --- | --- | --- | --- | --- |

BIBLIOGRAPHY 241
[92] Franco Brezzi, Annalisa Buffa, and Konstantin Lipnikov. Mimetic finite differences for
elliptic problems. ESAIM: Mathematical Modelling and Numerical Analysis, 43(2):277–
295, 2009.
[93] Ren´e Beltman, MJH Anthonissen, and Barry Koren. Conservative polytopal mimetic
Journal of Computational
| discretization | of the       | incompressible |              | Navier–Stokes |       | equations. |     |
| -------------- | ------------ | -------------- | ------------ | ------------- | ----- | ---------- | --- |
| and Applied    | Mathematics, |                | 340:443–473, |               | 2018. |            |     |
[94] J´erˆome Bonelle. Compatible discrete operator schemes on polyhedral meshes for elliptic
| and Stokes | equations. | Ph. | D. Thesis, | 2014. |     |     |     |
| ---------- | ---------- | --- | ---------- | ----- | --- | --- | --- |
[95] J´erˆome Bonelle, Alexandre Ern, and Riccardo Milani. Compatible discrete operator
schemes for the steady incompressible Stokes and Navier–Stokes equations. In Interna-
tional Conference on Finite Volumes for Complex Applications, pages 93–101. Springer,
2020.
[96] Daniele A Di Pietro and J´erˆome Droniou. An arbitrary-order discrete de Rham complex
on polyhedral meshes: Exactness, poincar´e inequalities, and consistency. Foundations of
| Computational | Mathematics, |     | 23(1):85–164, |     | 2023. |     |     |
| ------------- | ------------ | --- | ------------- | --- | ----- | --- | --- |
[97] Daniele A Di Pietro, J´erˆome Droniou, and Silvano Pitassi. Cohomology of the discrete
Calcolo,
de Rham complex on domains of general topology. 60(2):1–25, 2023.
[98] Vincent Perrier. Development of discontinuous Galerkin methods for hyperbolic systems
that preserve a curl or a divergence constraint: the case of linear systems. arXiv preprint
| arXiv:2405.04347, |     | 2024. |     |     |     |     |     |
| ----------------- | --- | ----- | --- | --- | --- | --- | --- |
[99] Alexandre Ern and Jean-Luc Guermond. Theory and practice of finite elements, volume
| 159. Springer, | 2004. |     |     |     |     |     |     |
| -------------- | ----- | --- | --- | --- | --- | --- | --- |
[100] DouglasNArnold,RichardSFalk,andRagnarWinther. Finiteelementexteriorcalculus,
homological techniques, and applications. Acta numerica, 15:1–155, 2006.
[101] Snorre H Christiansen. Stability of hodge decompositions in finite element spaces of
differential forms in arbitrary dimension. Numerische Mathematik, 107(1):87–106, 2007.
[102] P-ARaviartandJean-MarieThomas. Primalhybridfiniteelementmethodsfor2ndorder
|          |            | Mathematics |     | of computation, |     |                  |       |
| -------- | ---------- | ----------- | --- | --------------- | --- | ---------------- | ----- |
| elliptic | equations. |             |     |                 |     | 31(138):391–413, | 1977. |
[103] Philippe G Ciarlet. Basic error estimates for elliptic problems. Handbook of Numerical
| Analysis, | Vol II, Finite | Element |     | Methods | (Part | 1), 1991. |     |
| --------- | -------------- | ------- | --- | ------- | ----- | --------- | --- |
[104] Marie E Rognes, Robert C Kirby, and Anders Logg. Efficient assembly of H(div) and
H(curl) conforming finite elements. SIAM Journal on Scientific Computing, 31(6):4130–
4151, 2010.
[105] Douglas N Arnold, Daniele Boffi, and Richard S Falk. Quadrilateral H(div) finite ele-
| ments. | SIAM Journal | on  | Numerical | Analysis, | 42(6):2429–2451, |     | 2005. |
| ------ | ------------ | --- | --------- | --------- | ---------------- | --- | ----- |

BIBLIOGRAPHY 242
[106] Jean-Claude N´ed´elec. Mixed finite elements in R. Numerische Mathematik, 35(3):315–
341, 1980.
[107] Martin Werner Licht. Complexes of discrete distributional differential forms and their
homology theory. Foundations of Computational Mathematics, 17(4):1085–1122, 2017.
[108] Allen Hatcher. Algebraic Topology. Cambridge University Press, 2001.
[109] Alexander Linke. On the role of the Helmholtz decomposition in mixed methods for
incompressibleflowsandanewvariationalcrime. Computermethodsinappliedmechanics
| and engineering, | 268:782–800, | 2014. |     |     |     |
| ---------------- | ------------ | ----- | --- | --- | --- |
[110] Volker John, Alexander Linke, Christian Merdon, Michael Neilan, and Leo G Rebholz.
On the divergence constraint in mixed finite element methods for incompressible flows.
| SIAM review, | 59(3):492–544, | 2017. |     |     |     |
| ------------ | -------------- | ----- | --- | --- | --- |
[111] Keith J Galvin, Alexander Linke, Leo G Rebholz, and Nicholas E Wilson. Stabilizing
poor mass conservation in incompressible flow problems with large irrotational forcing
and application to thermal convection. Computer Methods in Applied Mechanics and
| Engineering, | 237:166–176, | 2012. |     |     |     |
| ------------ | ------------ | ----- | --- | --- | --- |
[112] Bernardo Cockburn. Discontinuous Galerkin methods. ZAMM-Journal of Applied Math-
ematics and Mechanics/Zeitschrift fu¨r Angewandte Mathematik und Mechanik: Applied
| Mathematics | and Mechanics, | 83(11):731–754, |     | 2003. |     |
| ----------- | -------------- | --------------- | --- | ----- | --- |
[113] Xiaozhe Hu, Carmen Rodrigo, Francisco J Gaspar, and Ludmil T Zikatanov. A non-
conforming finite element method for the Biot’s consolidation model in poroelasticity.
Journal of Computational and Applied Mathematics, 310:143–154, 2017.
[114] Wasilij Barsukow. Stationarity preserving schemes for multi-dimensional linear systems.
| Mathematics | of Computation, | 88(318):1621–1645, |     | 2019. |     |
| ----------- | --------------- | ------------------ | --- | ----- | --- |
[115] Daniele Boffi and Lucia Gastaldi. Some remarks on quadrilateral mixed finite elements.
| Computers | & structures, | 87(11-12):751–757, |     | 2009. |     |
| --------- | ------------- | ------------------ | --- | ----- | --- |
[116] LourencoBeir˜aodaVeiga,FrancoBrezzi,AndreaCangiani,GianmarcoManzini,LDona-
Math-
tella Marini, and Alessandro Russo. Basic principles of virtual element methods.
ematical Models and Methods in Applied Sciences, 23(01):199–214, 2013.
[117] L Beir˜ao da Veiga, Franco Brezzi, Luisa Donatella Marini, and Alessandro Russo. The
hitchhiker’s guide to the virtual element method. Mathematical models and methods in
| applied | sciences, 24(08):1541–1573, |     | 2014. |     |     |
| ------- | --------------------------- | --- | ----- | --- | --- |
[118] Wasilij Barsukow. Truly multi-dimensional all-speed schemes for the Euler equations on
| cartesian | grids. Journal | of Computational |     | Physics, 435:110216, | 2021. |
| --------- | -------------- | ---------------- | --- | -------------------- | ----- |
[119] Bernardo Cockburn, George E Karniadakis, and Chi-Wang Shu. The development of dis-
continuous Galerkin methods. In Discontinuous Galerkin methods: theory, computation
| and applications, | pages | 3–50. Springer, | 2000. |     |     |
| ----------------- | ----- | --------------- | ----- | --- | --- |

BIBLIOGRAPHY 243
[120] Bernardo Cockburn, George E Karniadakis, and Chi-Wang Shu. Discontinuous Galerkin
methods: theory, computation and applications, volume 11. Springer Science & Business
Media, 2012.
[121] Thomas Y Hou and Philippe G LeFloch. Why nonconservative schemes converge to
wrong solutions: error analysis. Mathematics of computation, 62(206):497–530, 1994.
[122] Franc¸ois Bouchut. Nonlinear stability of finite Volume Methods for hyperbolic conserva-
tion laws: And Well-Balanced schemes for sources. Springer Science & Business Media,
2004.
[123] ChristopheBerthonandChristopheChalons. Afullywell-balanced,positiveandentropy-
satisfying Godunov-type method for the shallow-water equations. Mathematics of Com-
putation, 85(299):1281–1307, 2016.
[124] E´lie Cartan. Sur certaines expressions diff´erentielles et le probl`eme de pfaff. In Annales
scientifiques de l’E´cole normale sup´erieure, volume 16, pages 239–332, 1899.