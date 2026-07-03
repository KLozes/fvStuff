| A finite     | volume   | scheme   | accurate      | at low   | Mach number | on  |
| ------------ | -------- | -------- | ------------- | -------- | ----------- | --- |
| quadrangular |          | mesh     | by space      | velocity | enrichment  |     |
|              |          | Jonathan | Jung, Vincent | Perrier  |             |     |
| To cite this | version: |          |               |          |             |     |
Jonathan Jung, Vincent Perrier. A finite volume scheme accurate at low Mach number on quadrangular
meshbyspacevelocityenrichment. ECCOMAS2024-9thEuropeanCongressonComputationalMethodsin
AppliedSciencesandEngineering,Jun2024,Lisbonne,Portugal. ⟨hal-04836529⟩
|     |     | HAL | Id: hal-04836529 |     |     |     |
| --- | --- | --- | ---------------- | --- | --- | --- |
https://hal.science/hal-04836529v1
Submittedon13Dec2024
HAL is a multi-disciplinary open access archive L’archiveouvertepluridisciplinaireHAL,estdes-
for the deposit and dissemination of scientific re- tinée au dépôt et à la diffusion de documents scien-
searchdocuments,whethertheyarepublishedornot. tifiquesdeniveaurecherche,publiésounon,émanant
Thedocumentsmaycomefromteachingandresearch des établissements d’enseignement et de recherche
institutionsinFranceorabroad,orfrompublicorpri- français ou étrangers, des laboratoires publics ou
| vateresearchcenters. |     |     | privés. |     |     |     |
| -------------------- | --- | --- | ------- | --- | --- | --- |
DistributedunderaCreativeCommonsCCBY4.0-Attribution-InternationalLicense

The9thEuropeanCongressonComputationalMethodsinAppliedSciencesandEngineering
ECCOMASCongress2024
3-–7June2024,Lisbon,Portugal
| A finite     | volume | scheme | accurate | at          | low Mach   | number | on  |
| ------------ | ------ | ------ | -------- | ----------- | ---------- | ------ | --- |
| quadrangular |        | mesh   | by space | velocity    | enrichment |        |     |
|              |        | J.     | Jung1,   | V. Perrier2 |            |        |     |
1 LMA, E2S-UPPA, and Cagire team, Inria Bordeaux Sud-Ouest, Avenue de
|     |     | l’Universit´e, | 64013 | Pau Cedex, | France. |     |     |
| --- | --- | -------------- | ----- | ---------- | ------- | --- | --- |
2
Cagire team, Inria Bordeaux Sud-Ouest and LMA, E2S-UPPA, Avenue de
|     |     | l’Universit´e, | 64013 | Pau Cedex, | France. |     |     |
| --- | --- | -------------- | ----- | ---------- | ------- | --- | --- |
Keywords: Low Mach number flows, Finite volume scheme, Compressible solver, Finite
| element method, | fluid | mechanics |     |     |     |     |     |
| --------------- | ----- | --------- | --- | --- | --- | --- | --- |
Classical finite volume schemes for compressible Euler system are not accurate on quad-
rangular mesh at low Mach number in the sense that they do not converge to the in-
compressible limit when the Mach number tends to zero [1]. The spurious mode that
jeopardizes the convergence can be identified and corresponds to the long time limit of a
first order wave system whose properties and discretization depend on the scheme used
for the compressible system [2]. Then, the low Mach number accuracy can be analysed
by studying the long time solution of the associated wave system.
In this presentation, we will propose to enrich the velocity space approximation by adding
adivergencefreeelementineachcell. UsingthisnewapproximationspaceandaGodunov’
numerical flux, the long-time limit of the wave system discretization corresponds to the
divergence free component of the initial condition (as expected). This divergence free
component will be defined via a discrete Hodge-Helmoltz decomposition.
As a consequence, for Euler equations, using this new approximation space and a numer-
ical flux that degenerates in the low Mach number limit to a Godunov’ flux for the wave
system, the numerical scheme provides pressure fluctuations of order Mach square and a
divergence of the velocity field of order Mach in the low Mach number limit.
This work can be seen as an extension to quadrangles of the properties obtained on
triangles [3].
REFERENCES
[1] H. Guillard and C. Viozat, On the behaviour of upwind schemes in the low Mach
| number | limit. Computers | & fluids, | 28(1), | pp. 63–86, | 1999. |     |     |
| ------ | ---------------- | --------- | ------ | ---------- | ----- | --- | --- |
[2] J. Jung and V. Perrier, Steady low Mach number flows: identification of the spurious
mode and filtering method. Journal of Computational Physics, 468, 111462, 2022.
[3] H. Guillard, On the behavior of upwind schemes in the low Mach number limit. IV:
P0 approximation on triangular and tetrahedral cells. Computers & fluids, 38(10),
| pp. 1969-1972, | 2009. |     |     |     |     |     |     |
| -------------- | ----- | --- | --- | --- | --- | --- | --- |