Journal of Computational Physics 555 (2026) 114764
Contents lists available at ScienceDirect
Journal of Computational Physics
journal homepage: www.elsevier.com/locate/jcp
A unified framework for non-linear reconstruction schemes in a
compact stencil. Part 2: Learning operators from neural networks
Minsheng Huang a, Xi Deng b,∗, Omar K. Matar b, Wenjun Ying a,c
aSchool of Mathematical Sciences, Shanghai Jiao Tong University, Shanghai, PR China
bDepartment of Chemical Engineering, Imperial College London, SW7 2AZ, United Kingdom
cInstitute of Natural Sciences, MOE-LSC, Shanghai Jiao Tong University, Shanghai, PR China
a r t i c l e i n f o a b s t r a c t
Keywords: In the preceding part, a family of ROUND (Reconstruction Operators on Unified Normalised-
ROUND variable Diagram) schemes was introduced within a unified framework that includes existing non-
WENO3-NN linear reconstruction schemes developed in a compact stencil. While the initial ROUND schemes
Neural network
demonstrated promising performance, there remains potential to enhance their spectral properties
Low dissipation
by designing alternative functions within this framework. Recently, Bezgin et al.[1] proposed the
Shock capturing
WENO3-NN scheme, which employs neural networks as weighting functions in Weighted Essen-
tially Non-Oscillatory (WENO) schemes. Although WENO3-NN exhibits improved accuracy and
reduced numerical dissipation, it incurs higher computational costs compared to conventional
schemes. In this study, we first unify the WENO3-NN scheme within the normalised-variable di-
agram framework. This analysis reveals that WENO3-NN deviates from the second-order upwind
scheme and tends toward the second-order central scheme as the normalised variable decreases
below zero. For less smooth functions, the cell-face reconstruction loss dominates over the ideal
weight loss, driving WENO3-NN toward central behaviour. Motivated by this insight, we propose
a new ROUND formulation that mirrors the behaviour of WENO3-NN while including a greater
contribution from the central scheme to further enhance the dissipation property. The resulting
ROUND scheme preserves the key features of WENO3-NN on the normalised-variable diagram
while maintaining algorithmic simplicity. Building on this formulation, we further develop a
low-dissipative ROUND (LD-ROUND) scheme by introducing controlled anti-dissipation errors.
Comparative study of numerical error versus CPU cost demonstrates that both ROUND and LD-
ROUND achieve substantial efficiency gains over WENO3-NN. For example, using the same grid,
LD-ROUND produces numerical errors and CPU costs that are both less than half of those gener-
ated by WENO3-NN. Finally, we validate the proposed schemes on a range of compressible single-
and multi-phase flow problems, demonstrating their superior performance. In several benchmark
tests, the LD-ROUND scheme achieves higher resolution than classical fifth-order WENO schemes.
1. Introduction
Solving hyperbolic partial differential equation systems–such as time-dependent shallow water equations and Euler equations
for compressible flows–requires numerical schemes that can simultaneously capture discontinuities and resolve fine-scale flow fea-
tures. Balancing these competing demands, namely ensuring sufficient numerical dissipation near discontinuities while maintaining
∗ Corresponding author.
E-mail address: x.deng@imperial.ac.uk (X. Deng).
https://doi.org/10.1016/j.jcp.2026.114764
Received 7 July 2025; Received in revised form 20 October 2025; Accepted 6 February 2026
Available online 11 February 2026
0021-9991/© 2026 The Author(s). Published by Elsevier Inc. This is an open access article under the CC BY license
(h ttp://creativecommons.org/licenses/by/4.0/ ).

M. Huang, X. Deng, O.K. Matar et al. Journal of Computational Physics 555 (2026) 114764
low dissipation in smooth regions, remains a persistent challenge in the design of high-fidelity numerical methods. Over the past
decades, high-order numerical methods such as the Finite Difference Method (FDM) [2–4] and the Flux Reconstruction (FR) method
[5–7], which includes the Discontinuous Galerkin (DG) approach [8], have been developed for solving systems of hyperbolic partial
differential equations. Nevertheless, the Finite Volume Method (FVM), particularly developed in a compact stencil, remains widely
adopted in computational fluid dynamics (CFD) for engineering applications such as unstructured-grid-based methods [9,10]. More-
over, FVM in a compact stencil continues to serve as a robust candidate for reconstruction in some high-order schemes such as the
Multi-dimensional Optimal Order Detection (MOOD) method [11,12] and hybrid DG/FV approaches [13–15].
However, designing an accurate and non-oscillatory reconstruction scheme within a compact stencil is non-trivial. One key chal-
lenge arises from Godunov’s theorem, which states that linear schemes cannot simultaneously achieve second-order (or higher)
accuracy and maintain monotonicity. Additionally, as noted by Liu and Shen[16], a three-cell stencil cannot effectively distinguish
between smooth regions and discontinuities. Over recent decades, numerous methodologies have been proposed to develop non-
oscillatory, second-order (or higher) accurate schemes within compact stencils. Representative approaches include Total Variation
Diminishing (TVD) schemes [17,18], Normalised Variable Diagram (NVD) schemes [19–21], and third-order Weighted Essentially
Non-Oscillatory (WENO) schemes [22,23], among others. A comprehensive review of these compact-stencil-based schemes is provided
in the preceding part [24]. The limitations of these schemes have also been highlighted in the review by Deng[24]. For instance,
TVD schemes that satisfy the symmetry condition often suffer from accuracy degradation and distortion of the transported scalar
profile. WENO schemes that rely on smoothness indicators face challenges, as three-cell-based indicators are insufficient for reliably
distinguishing between smooth regions and discontinuities.
Machine learning (ML) methods have recently emerged as powerful tools for enhancing numerical schemes in CFD. Relevant
studies are generally categorised into three major types: data-driven surrogate models, physics-constrained surrogate models, and
ML-assisted numerical methods [25]. Although both data-driven and physics-constrained surrogate models are innovative and ex-
hibit significant potential, they often face challenges in end-to-end direct simulations, particularly during long-term simulations where
error accumulation becomes critical and the scenarios are not explicitly included during the training process. In contrast, ML-assisted
methods augment traditional numerical schemes rather than entirely replacing them, frequently demonstrating comparable or su-
perior resolution, accuracy, efficiency, and generality [26]. Specifically, ML-assisted numerical schemes have garnered substantial
attention due to their capability to significantly enhance numerical performance [27–30], which is also the primary focus of this
work. Pioneering research by Stevens and Colonius [31] proposed an enhanced fifth-order WENO5 scheme, improving the resolution
of the classic WENO5-JS scheme [23] by perturbing its non-linear weights. Building on their earlier work [32], Bezgin et al. [1]
developed the WENO3-NN scheme within a compact stencil, utilizing a neural network to predict optimal non-linear weights, thus
refining the conventional three-stencil WENO approach. Subsequently, Shahane et al.[33] adapted rational functions with enhanced
representational properties to create the rational-WENO3-NN scheme, achieving higher efficiency and reduced errors in the low-
resolution regime. Inspired by the WENO3-NN framework, Nogueira et al. [34] introduced the WENO5-NN scheme, where a neural
network regressor directly determines non-linear weights. Expanding on the same technique as the WENO3-NN scheme, Zhang et al.
[35] extended the third-order Weighted Compact Non-linear Scheme (WCNS) as WCNS3-NN scheme, achieving the scale-invariant
property and ENO-property. Based on previous studies and the Multi-Resolution (MR) strategy [36,37], Fan et al.[26] designed the
WCNS3-MR-NN scheme, which incorporates multi-resolution capabilities and achieves optimal accuracy across various grid resolu-
tions. In this study, we focus on the WENO3-NN scheme, which is developed within a compact stencil and demonstrates enhanced
accuracy and dissipation properties compared to conventional methods. However, these improvements come at the cost of increased
computational expense relative to traditional schemes.
In the preceding part [24], a new family of high-resolution, non-oscillatory schemes termed ROUND (Reconstruction Operator
on Unified Normalised-variable Diagram) was proposed. The ROUND scheme is formulated within the Unified Normalised-variable
Diagram (UND), which provides a unified framework for existing second-order polynomial-based and non-polynomial-based shock-
capturing schemes developed in a compact stencil. The ROUND scheme has been further extended to unstructured grids using the
virtual upwind point method [38], and to non-uniform grids as demonstrated in Deng et al.[39]. The efficacy of ROUND schemes
has been demonstrated through numerical simulations of high-speed compressible flows [40] and combustion processes [41]. Despite
the promising performance of the initial ROUND schemes, there remains potential to construct alternative functions within the UND
framework to further enhance the spectral properties of ROUND schemes. Therefore, this work unifies the WENO3-NN scheme into
the normalised-variable diagram and designs a new ROUND formulation with improved accuracy. The behaviour of WENO3-NN on
the normalised-variable diagram reveals the mechanism underlying its improved dissipation properties. As the normalised variable
decreases below zero, WENO3-NN diverges from the second-order upwind scheme and approaches the second-order central scheme.
Negative values of the normalised variable indicate regions containing the critical points, where the reconstruction loss becomes more
dominant than the ideal weight loss in the WENO3-NN scheme. Building on insights from WENO3-NN, we introduce a novel ROUND
formulation that emulates its behaviour while incorporating a stronger central scheme component to further improve dissipation
properties. This new scheme retains the essential features of WENO3-NN on the normalised-variable diagram and remains algorith-
mically simple. Extending this approach, we develop a low-dissipative variant, LD-ROUND, by introducing controlled anti-dissipation
mechanisms. A comparative analysis of numerical error versus computational cost reveals that both ROUND and LD-ROUND offer
significant efficiency improvements over WENO3-NN. The proposed schemes are validated across a range of compressible single- and
multi-phase flow problems, demonstrating superior performance. In several benchmark cases, the resolution achieved by LD-ROUND
even surpasses classical fifth-order WENO schemes.
The remainder of this paper is organised as follows. Section 2 provides a brief introduction to the finite volume ROUND
scheme. In Section 3, we review the WENO3-NN scheme and reveal its mechanism for reducing reconstruction errors using the
2

M. Huang, X. Deng, O.K. Matar et al. Journal of Computational Physics 555 (2026) 114764
normalised-variable diagram. Section 4 presents a novel ROUND formulation inspired by the behaviour of WENO3-NN, along with
a low-dissipation variant incorporating controlled anti-dissipation errors. The theoretical properties of the proposed schemes are
discussed in Section 5. Section 6 reports a wide range of numerical experiments to evaluate the performance of the new methods.
Finally, concluding remarks are provided in Section 7.
2.  ROUND schemes in a compact stencil
2.1.  Finite volume method
To illustrate the finite volume method, we consider the one-dimensional scalar conservation law, expressed as
| 𝜕𝜙 𝜕𝑓(𝜙) |     |     |     |     |     |     | (2.1) |
| -------- | --- | --- | --- | --- | --- | --- | ----- |
| +        | =0, |     |     |     |     |     |       |
𝜕𝑡 𝜕𝑥
where 𝜙(𝑥,𝑡) is the solution function and 𝑓(𝜙) the flux function. Following the standard finite volume method, the computational
domain is decomposed into 𝑁 non-overlapping finite volume cells,  ∶𝑥∈[𝑥 ,𝑥 𝑖+1],𝑖=1,2,…𝑁, of uniform grid spacing ℎ=
𝑖 𝑖−1
|                                            |     |              |     | 2   | 2   |     |     |
| ------------------------------------------ | --- | ------------ | --- | --- | --- | --- | --- |
| 𝑥 −𝑥 . We define the cell average value 𝜙̄ |     |  over cell  |  as |     |     |     |     |
| 𝑖+1 𝑖−1                                    |     | 𝑖            | 𝑖   |     |     |     |     |
2 2
| 1      | 𝑥 𝑖+1       |     |     |     |     |     |       |
| ------ | ----------- | --- | --- | --- | --- | --- | ----- |
| 𝜙̄(𝑡)= | 2 𝜙(𝑥,𝑡)𝑑𝑥. |     |     |     |     |     | (2.2) |
𝑖 ℎ∫
|     | 𝑥 𝑖−1 |     |     |     |     |     |     |
| --- | ----- | --- | --- | --- | --- | --- | --- |
2
Then the cell average of each cell, 𝜙̄, can be updated by the following semi-discrete form
𝑖
𝑑 𝜙̄
| 𝑖 =− 1 | (𝑓̃ −𝑓̃ 𝑖−1), |     |     |     |     |     | (2.3) |
| ------ | ------------- | --- | --- | --- | --- | --- | ----- |
𝑖+1
| 𝑑 𝑡 ℎ | 2 2 |     |     |     |     |     |     |
| ----- | --- | --- | --- | --- | --- | --- | --- |
where 𝑓̃ is the numerical flux across cell boundaries and is computed by using a Riemann solver
| 𝑓̃ =𝑓 𝑅 | 𝑖 𝑒𝑚𝑎𝑛𝑛(𝜙𝐿 ,𝜙𝑅 ). |     |     |     |     |     | (2.4) |
| ------- | ----------------- | --- | --- | --- | --- | --- | ----- |
| 𝑖+1     | 1 𝑖+1 𝑖+1         |     |     |     |     |     |       |
2 𝑖 +
|     | 2 2 2 |     |     |     |     |     |     |
| --- | ----- | --- | --- | --- | --- | --- | --- |
The left-side cell face value 𝜙𝐿  and right-side cell face value 𝜙𝑅  are computed from the reconstructions over left- and right-biased
|     | 𝑖+1 |     | 𝑖+1 |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- |
stencils, respectively. The Riem2ann flux usually can be formulated2  in a canonical form as follows:
|              | (                 |       |          | )   |     |     |       |
| ------------ | ----------------- | ----- | -------- | --- | --- | --- | ----- |
| 𝑅 𝑖 𝑒𝑚𝑎𝑛𝑛(𝜙𝐿 | ,𝜙𝑅 1 𝑓(𝜙𝐿 )+𝑓(𝜙𝑅 |       | |(𝜙𝐿 −𝜙𝑅 |     |     |     | (2.5) |
| 𝑓            | )=                | )−|𝛼̃ | 𝑖+1      | ) , |     |     |       |
| 𝑖 + 1        | 𝑖+1 𝑖+1 2 𝑖+1     | 𝑖+1   | 2 𝑖+1    | 𝑖+1 |     |     |       |
| 2            | 2 2 2             | 2     | 2        | 2   |     |     |       |
where 𝛼̃  stands for a characteristic speed in a hyperbolic equation. The semi-discrete formulation (2.3) can be advanced using
𝑖+1
explicit tim2e integration schemes, with this study opting for the 3rd-order Runge-Kutta scheme. The main task is approximating the
left-side and right-side cell face values via a reconstruction scheme. The following sections will give details of how to reconstruct the
left-side value 𝜙𝐿  through the reconstruction operator 𝜙𝐿 =𝐿[𝜙̄ ,𝜙̄,𝜙̄ ] in a compact stencil. For clarity, the superscript
𝑖−1 𝑖 𝑖+1
|     | 𝑖+1 |     | 𝑖+1∕2 | 𝑖   |     |     |     |
| --- | --- | --- | ----- | --- | --- | --- | --- |
𝐿 will be omitted i2n the following section. Reconstructing the right-side value can then be achieved symmetrically. Unless otherwise
specified, the HLLC Riemann solver [42] is employed in this work to approximate the numerical fluxes following the reconstruction
step.
2.2.  ROUND methodology
According to Godunov’s theorem, a non-linear reconstruction scheme is required to obtain non-oscillatory solutions with accuracy
higher than the first order. Instead of directly devising a non-linear reconstruction function, the ROUND methodology constructs the
non-linear reconstruction function in the normalised variable space. We define the normalised variable 𝜙̂ and normalised reconstruc-
𝑖
tion value 𝜙̂  as follows
𝑖+1∕2
−𝜙̄
| 𝜙̄ −𝜙̄ | 𝜙 𝑖+1∕2          | 𝑖−1 |     |     |     |     |       |
| ------ | ---------------- | --- | --- | --- | --- | --- | ----- |
| 𝜙̂ = 𝑖 | 𝑖−1 , 𝜙̂ =       | .   |     |     |     |     | (2.6) |
| 𝑖 𝜙̄   | −𝜙̄ 𝑖+1∕2 𝜙̄ −𝜙̄ |     |     |     |     |     |       |
| 𝑖+1    | 𝑖−1 𝑖+1          | 𝑖−1 |     |     |     |     |       |
According to Proposition 1 in [24], the location of the first-order critical point within the compact stencil can affect the value of 𝜙̂. Let
𝑖
𝝓∈ℝ𝑛𝑠 collect the cell averages on the stencil and define the unified/normalised scalar by 𝜙∶=[𝝓]. The normalised reconstruction
operator ̂
[⋅] of cell 𝑖 is defined by
𝑖
| ̂ ∶ℝ→ℝ, | 𝜙̂ =̂ [𝜙̂ ]. |     |     |     |     |     | (2.7) |
| -------- | ------------- | --- | --- | --- | --- | --- | ----- |
| 𝑖        | 𝑖+1 𝑖 𝑖       |     |     |     |     |     |       |
2
This differs from the conventional reconstruction operator [⋅]:
𝑖
|  ∶ℝ𝑛𝑠 | = [𝝓̄       |     |     |     |     |     | (2.8) |
| ------ | ------------ | --- | --- | --- | --- | --- | ----- |
| 𝑖 →ℝ,  | 𝜙 𝑖+1 𝑖 𝑖 ], |     |     |     |     |     |       |
2
where 𝑛  denotes the stencil size. In this work, we consider the three-cell-stencil scheme (𝑛 =3) and 𝝓̄ =(𝜙̄ ,𝜙̄,𝜙̄ ).  For
| 𝑠   |     |     |     |     | 𝑠   | 𝑖 𝑖−1 𝑖 | 𝑖+1 |
| --- | --- | --- | --- | --- | --- | ------- | --- |
clarity, we denote the normalised operator by ̂  and the conventional operator by , respectively. Once 𝜙̂  is obtained via
|     |     | 𝑖   |     |     | 𝑖   | 𝑖+1∕2 |     |
| --- | --- | --- | --- | --- | --- | ----- | --- |
3

M. Huang, X. Deng, O.K. Matar et al. Journal of Computational Physics 555 (2026) 114764
the operator ̂ , the physical cell boundary value can be recovered through a denormalisation process. The normalisation approach
𝑖
offers several advantages. First, it provides a unified framework for existing schemes developed within compact stencils. Second,
as demonstrated in [24], designing the reconstruction function in the normalised variable space enables the scheme to achieve
scale-invariant properties. Some common linear reconstruction schemes, such as the first-order upwind scheme 𝑈𝑊 =𝜙̄ , the
𝑖 𝑖−1
|     |     |     |     |     |     |     |     | 3𝜙̄ | −1𝜙̄ |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | ---- | --- |
first-order downwind scheme 𝐷𝑊 =𝜙̄ , the second-order upwind scheme 𝑈𝑊2= , the second-order central scheme
|     |     |     | 𝑖   | 𝑖+1 |     |     |     | 𝑖 2 𝑖 | 2 𝑖−1 |     |
| --- | --- | --- | --- | --- | --- | --- | --- | ----- | ----- | --- |
𝐶𝐷2= 1𝜙̄ +1𝜙̄ , and third-order linear upwind scheme 𝑃3=−1𝜙̄ +5𝜙̄ +1𝜙̄  can be projected into their corresponding
| 𝑖                                      | 𝑖    | 𝑖+1  |     |     |     | 𝑖   | 𝑖−1 | 𝑖 𝑖+1 |     |       |
| -------------------------------------- | ---- | ---- | --- | --- | --- | --- | --- | ----- | --- | ----- |
| normalised reconstruction operators as | 2    | 2    |     |     |     |     | 6   | 6 3   |     |       |
|                                        | ̂𝑈𝑊 | =𝜙̂, |     |     |     |     |     |       |     | (2.9) |
𝑖
𝑖
|     | ̂𝐷𝑊 | =1, |     |     |     |     |     |     |     | (2.10) |
| --- | ---- | --- | --- | --- | --- | --- | --- | --- | --- | ------ |
𝑖
3
|     | ̂𝑈𝑊2= | 𝜙̂,    |     |     |     |     |     |     |     | (2.11) |
| --- | ------ | ------ | --- | --- | --- | --- | --- | --- | --- | ------ |
|     | 𝑖      | 2 𝑖    |     |     |     |     |     |     |     |        |
|     | ̂𝐶𝐷2= | 1 𝜙̂   | 1   |     |     |     |     |     |     |        |
|     |        | +      | ,   |     |     |     |     |     |     | (2.12) |
|     | 𝑖      | 2 𝑖    | 2   |     |     |     |     |     |     |        |
|     | ̂𝑃3=  | 5 𝜙̂ 1 |     |     |     |     |     |     |     | (2.13) |
|     |        | 𝑖 +    | .   |     |     |     |     |     |     |        |
|     | 𝑖      | 6 3    |     |     |     |     |     |     |     |        |
According to the truncation error analysis in [24], the dissipation and anti-dissipation regions are defined in the normalised-variable
diagram. To design non-linear reconstruction operators that achieve both non-oscillatory behaviour and high-resolution accuracy,
Deng[24] further projects existing non-linear schemes–such as TVD, classical WENO, and non-polynomial-based methods like THINC
(Tangent of Hyperbola for INterface Capturing) [43]–into the normalised space. Based on the behaviour of these schemes in this
space, a Unified Normalised-variable Diagram (UND) is proposed, which encompasses both the TVD and ENO regions. This framework
enables the direct construction of a family of high-resolution, non-oscillatory schemes by designing appropriate functions within the
UND. These schemes are termed as ROUND schemes. Despite the promising performance of the initial ROUND schemes demonstrated
by Deng[24,38], there remains potential to construct alternative functions within the UND framework to further enhance the spectral
properties of ROUND schemes.
3.  Review of the WENO3 scheme trained by deep neural network
3.1.  WENO methodology
It is known that the third-order linear scheme can be expressed as a linear combination of the second-order upwind scheme and
the second-order central scheme:
|     | 𝑃3=𝑑 | 𝑈𝑊2+𝑑 | 𝐶𝐷2, |     |     |     |     |     |     | (3.1) |
| --- | ----- | ------ | ----- | --- | --- | --- | --- | --- | --- | ----- |
|     | 𝑖     | 0 𝑖    | 1 𝑖   |     |     |     |     |     |     |       |
where {𝑑 ,𝑑 }={1,2}, which are ideal weights to achieve uniform third-order accuracy.
|                                                                                                                                 | 0   | 1   |     |     |     |     |     |     |     |     |
| ------------------------------------------------------------------------------------------------------------------------------- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| To suppress numerical oscillations introduced by high-order linear interpolation, WENO schemes replace the linear coefficients  |     | 3 3 |     |     |     |     |     |     |     |     |
} by non-linear weighting functions {𝜔 }, resulting in the WENO reconstruction operator 𝑊𝐸𝑁𝑂 defined as
| {𝑑 0 ,𝑑 | 1   |     |     |     | 0 ,𝜔 1 |     |     |     |     |     |
| ------- | --- | --- | --- | --- | ------ | --- | --- | --- | --- | --- |
𝑖
|     | 𝑊𝐸𝑁𝑂=𝜔 | 𝑈𝑊2+𝜔 | 𝐶𝐷2, |     |     |     |     |     |     | (3.2) |
| --- | ------- | ------ | ----- | --- | --- | --- | --- | --- | --- | ----- |
|     | 𝑖       | 0      | 𝑖 1   | 𝑖   |     |     |     |     |     |       |
For the third-order WENO schemes, the nonlinear weighting functions {𝜔 ,𝜔 } are designed to meet the following conditions:
0 1
| Consistency Condition: 𝜔 |     |     |            | =1;    |     |     |     |     |     |     |
| ------------------------ | --- | --- | ---------- | ------ | --- | --- | --- | --- | --- | --- |
| •                        |     |     | 0 +𝜔       | 1      |     |     |     |     |     |     |
| • Accuracy Condition: 𝜔  |     |     | −𝑑 =𝑂(ℎ2), | 𝑘=0,1; |     |     |     |     |     |     |
|                          |     |     | 𝑘 𝑘        |        |     |     |     |     |     |     |
• Essentially Non-oscillatory Condition: 𝜔 →0 if and only if the corresponding stencil contains a discontinuity.
𝑘
To satisfy the above conditions, the classical WENO weighting functions are formulated based on the so-called smoothness indicators,
defined as
|     |      |           | [       | ]2    |     |         | [    | ]2  |     |       |
| --- | ---- | --------- | ------- | ----- | --- | ------- | ---- | --- | --- | ----- |
|     |      | ∑ 2 𝑥 𝑖+1 | 𝜕𝑙𝑈    | 𝑊2    | ∑ 2 | 𝑥 𝑖+1   | 𝜕𝑙𝐶 | 𝐷2  |     |       |
|     | 𝐼𝑆 = |           | 2 ℎ2𝑙−1 | 𝑖 𝑑𝑥, | 𝛽 = | 2 ℎ2𝑙−1 | 𝑖    | 𝑑𝑥. |     | (3.3) |
|     | 0    | ∫         | 𝜕𝑥𝑙     |       | 1   | ∫       | 𝜕𝑥𝑙  |     |     |       |
|     |      | 𝑙=1 𝑥     |         |       | 𝑙=1 | 𝑥       |      |     |     |       |
|     |      | 𝑖−1       |         |       |     | 𝑖−1     |      |     |     |       |
|     |      | 2         |         |       |     | 2       |      |     |     |       |
The explicit forms of smoothness indicators for second-order linear schemes are given as follows:
|     | =(𝜙̄ | −𝜙̄)2,𝐼𝑆 | =(𝜙̄ | −𝜙̄ )2. |     |     |     |     |     | (3.4) |
| --- | ---- | -------- | ---- | ------- | --- | --- | --- | --- | --- | ----- |
|     | 𝐼𝑆 0 | 𝑖−1      | 𝑖 1  | 𝑖 𝑖+1   |     |     |     |     |     |       |
The weighting functions in WENO-JS proposed in [23] are given as
|     |       | 𝜁 𝑘   | 𝑑     | 𝑘        |     |     |     |     |     | (3.5) |
| --- | ----- | ----- | ----- | -------- | --- | --- | --- | --- | --- | ----- |
|     | 𝜔 𝑘 = | ,     | 𝜁 𝑘 = | , 𝑘=0,1, |     |     |     |     |     |       |
|     |       | ∑1 𝜁  | (𝜖+𝐼𝑆 | )2       |     |     |     |     |     |       |
|     |       | 𝑗=0 𝑗 |       | 𝑘        |     |     |     |     |     |       |
where the parameter 𝜖 is a small positive constant introduced to prevent division by zero.
To reduce the excessive numerical dissipation of the WENO-JS scheme, Borges et al. [44] proposed a novel WENO variant, known
as WENO-Z, by introducing a higher-order reference smoothness indicator. This new indicator is constructed through a specific
4

M. Huang, X. Deng, O.K. Matar et al. Journal of Computational Physics 555 (2026) 114764
| Fig. 1. Schematic of the WENO3-NN architecture. The stencil 𝑆={𝜙̄ |     |     | ,𝜙̄,𝜙̄ |     |
| ----------------------------------------------------------------- | --- | --- | ------ | --- |
} is initially processed through the Delta layer, which produces
| Galilean invariant features Δ for 𝑗∈{1,2,3,4}. The resulting hidden representation is then fed into three fully connected layers, each consisting  |     |     | 𝑖−1 𝑖 𝑖+1 |     |
| -------------------------------------------------------------------------------------------------------------------------------------------------- | --- | --- | --------- | --- |
𝑗
of 16 nodes. During training, the neural network outputs the predicted weights 𝜔NN. In practical applications, the final weights 𝜔̃NN are obtained
| through an ENO layer. |     |     | 𝑘   | 𝑘   |
| --------------------- | --- | --- | --- | --- |
combination of the classical smoothness indicators used in WENO-JS. Afterwards, the third-order WENO-Z3 scheme [45] is presented
by considering the reference smoothness indicator of the form
( 𝜏 )
| 𝜁 =𝑑 1+ | 𝑍 , 𝑗=0,1, |     |     | (3.6) |
| ------- | ---------- | --- | --- | ----- |
| 𝑗 𝑗     | 𝐼𝑆 +𝜖      |     |     |       |
𝑘
(3.7)
| 𝜏 𝑍 =|𝐼𝑆 | −𝐼𝑆 1|. |     |     |     |
| -------- | ------- | --- | --- | --- |
0
The WENO-Z3 smoothness indicator assigns a larger weight to the stencil containing discontinuities than the WENO-JS scheme. Other
contributions aimed at enhancing the performance of WENO schemes developed in a three-cell-based compact stencil can be found in
the literature, such as [46–50]. Nevertheless, the development of high-resolution, non-oscillatory schemes within a compact stencil
remains a challenging task.
3.2.  Basic formulation and WENO3-NN architecture
Unlike traditional approaches that rely on smoothness indicators to design WENO weighting functions, WENO3-NN [1] and its
variant [33] construct the weighting functions by training a neural network to minimise a chosen loss function. Following Bezgin
et al.[1], we briefly introduce the basic formulation and architecture of the WENO3-NN in this subsection. The main idea behind
the WENO3-NN is to train a neural-network-based surrogate model to predict the weights 𝜔  and 𝜔 , which is given as the following
0 1
formulation:
| 𝜔𝑁𝑁 =𝑁𝑁(𝜙̄ | ,𝜙̄,𝜙̄ ;Θ). |     |     | (3.8) |
| ---------- | ----------- | --- | --- | ----- |
| 𝑘          | 𝑖−1 𝑖 𝑖+1   |     |     |       |
Here, 𝑁𝑁 denotes a neural network and Θ presents the parameters of the neural network. Different from the traditional WENO
schemes, the weight of the target cell is determined by the trained neural network.
The architecture of the network is illustrated in Fig. 1.  The network architecture comprises an input layer, three hidden layers with
16 nodes each, and a softmax output layer. This configuration results in a total of 658 trainable parameters. The input stencil values
𝑆 ={𝜙̄ ,𝜙̄,𝜙̄ } are first processed through a Delta layer, which extracts physically meaningful features. This transformation is
| 𝑖 𝑖−1 𝑖 𝑖+1 |     |     |     |     |
| ----------- | --- | --- | --- | --- |
designed to enforce Galilean invariance, thereby incorporating prior knowledge into the model and enhancing its generalisation
capability. In [1], the following four features are used:
Δ̃ =|𝜙̄ −𝜙̄ 𝑖−1|,Δ̃ =|𝜙̄ −𝜙̄ 𝑖|,Δ̃ =|𝜙̄ −𝜙̄ 𝑖−1|,Δ̃ =|𝜙̄ −2𝜙̄ +𝜙̄ 𝑖−1|. (3.9)
| 1 𝑖 | 2 𝑖+1 | 3 𝑖+1 | 4 𝑖+1 𝑖 |     |
| --- | ----- | ----- | ------- | --- |
All features are then normalised into amplitudes of the first and second-order derivatives of the local input data:
| Δ =Δ̃ ∕max(Δ̃ | ,Δ̃ ,𝜖), for 𝑗={1,2,3,4}. |     |     | (3.10) |
| ------------- | ------------------------- | --- | --- | ------ |
| 𝑖 𝑗           | 1 2                       |     |     |        |
Here, 𝜖 is a small positive constant. After the Delta layer, the neural network consists of three hidden layers with 16 nodes in each
and a softmax output layer. To prevent the spurious oscillations, an ENO layer is constructed to restore the ENO property, which is
essentially a cut-off function defined as
{
| 𝜔𝑁     | 𝑁𝜓        | 0, if 𝜔𝑁 𝑁 <𝑐eno, |     |        |
| ------ | --------- | ----------------- | --- | ------ |
| 𝜔̃𝑁𝑁 = | 𝑘 𝑘 , 𝜓 = | 𝑘                 |     | (3.11) |
| 𝑘 ∑1   | 𝑗         | otherwise.        |     |        |
|        | 𝜔𝑁𝑁𝜓      | 1,                |     |        |
𝑗=0 𝑗 𝑗
5

M. Huang, X. Deng, O.K. Matar et al. Journal of Computational Physics 555 (2026) 114764
Here 𝑐eno  is the threshold and is set as 2×10−4, as recommanded in [1]. The WENO3-NN scheme is designed to predict suitable
weights in smooth and discontinuous stencils. In regions of perfect smoothness, the predicted weights 𝜔̃  should closely approximate
𝑘
the ideal weights 𝑑 , whereas in stencils containing discontinuities, 𝜔̃  should be driven toward zero. Therefore, the following loss
|     | 𝑘   |     | 𝑘   |     |     |
| --- | --- | --- | --- | --- | --- |
function is proposed in [1] to train the parameters Θ:
|  =   |             |     |     |     |     |
| ------ | ------------ | --- | --- | --- | --- |
| +𝛽     | +𝛽 𝑊||Θ|| 2, |     |     |     |     |
| 𝑙𝑜𝑠𝑠 𝑟 | 𝑑 𝑑          |     |     |     |     |
𝑁𝑏
1 ∑
|  =   | (𝛾𝑠)𝛼(𝜙𝑁𝑁,𝑠−𝜙𝑠 | )2, |     |     |        |
| ----- | -------------- | --- | --- | --- | ------ |
| 𝑟 𝑁   | 𝑖 𝑖+1          | 𝑖+1 |     |     | (3.12) |
| 𝑏 𝑠=1 | 2              | 2   |     |     |        |
𝑁𝑏 1
| 1 ∑      | ∑                 |     |     |     |     |
| -------- | ----------------- | --- | --- | --- | --- |
|  =      | [1−(𝛾𝑠)𝛼](𝜔𝑁𝑁,𝑠−𝑑 | )2. |     |     |     |
| 𝑑 𝑁      | 𝑖                 | 𝑘 𝑘 |     |     |     |
| 𝑏 𝑠=1𝑘=0 |                   |     |     |     |     |
The total loss, denoted as  , is comprised of three components: the cell face reconstruction loss , the ideal weight loss  , and
|     | 𝑙𝑜𝑠𝑠 |     |     | 𝑟   | 𝑑   |
| --- | ---- | --- | --- | --- | --- |
the 𝐿  regularization of the parameters Θ to prevent overfitting. The expressions for  and   calculate the average loss across
| 2   |     |     | 𝑟   | 𝑑   |     |
| --- | --- | --- | --- | --- | --- |
samples, with 𝑠 representing the input sample 𝑆 and 𝑁  denoting the number of mini-batch samples. 𝛾𝑠 is an indicator to measure
|     |     | 𝑖 𝑏 |     | 𝑖   |     |
| --- | --- | --- | --- | --- | --- |
the well-resolvedness of the function in the target stencil 𝑆. For stencils with perfectly smooth functions (e.g., linear functions), 𝛾𝑠
|     |     | 𝑖   |     |     | 𝑖   |
| --- | --- | --- | --- | --- | --- |
tends to zero, which allows the neural network to approximate the non-linear weights closely to their ideal counterparts. Conversely,
for a stencil 𝑆 containing discontinuities, 𝛾𝑠 tends to 1, emphasizing  and ensuring the network output 𝜙𝑁𝑁,𝑠 closely approximates
| 𝑖   |     | 𝑖   | 𝑟   | 𝑖+1 |     |
| --- | --- | --- | --- | --- | --- |
the interpolated value 𝜙𝑠 . For intermediate cases, 0<𝛾𝑠<1 enables the neural network to automatically2  balance reconstruction
𝑖+1 𝑖
loss and ideal weight loss d2uring training. Therefore, 𝛾𝑠 proposed in [1] is defined as
𝑖
| |𝜙̄       | −2𝜙̄ +𝜙̄ 𝑖+1|     |             |     |     |        |
| --------- | ----------------- | ----------- | --- | --- | ------ |
| 𝛾𝑠=       | 𝑖−1 𝑖             | , 𝜖 =10−15. |     |     | (3.13) |
| 𝑖 |𝜙̄ −𝜙̄ | 𝑖−1|+|𝜙̄ −𝜙̄ 𝑖|+𝜖 | 𝛾           |     |     |        |
| 𝑖         | 𝑖+1               | 𝛾           |     |     |        |
Furthermore, the 𝛼 in Eq. (3.12) is a parameter to control the scale separation mechanism. When 𝛼=0, the cell face reconstruction
loss becomes dominant. As 𝛼 approaches one, the weights obtained from the neural network will become closer to the ideal weights.
In this work, we focus on the WENO3-NN scheme with parameters 𝛼=0.1 and 𝛽 =0.1, as these values have been shown to
𝑑
improve the dissipative properties of the scheme in [1].  It’s noted that the choice of these parameters does not guarantee optimality.
To pursue optimal performance, more advanced techniques such as Bayesian optimisation [51] or reinforcement learning frameworks
[52] could be employed. However, due to the inherent complexity of high-speed compressible flows and the varying definitions of
optimisation objectives, the resulting schemes inevitably represent a compromise among low dissipation, low dispersion, and effective
discontinuity capturing.
It’s also noted that the architecture of the neural network used in this study is not optimal and may be overparameterised.
Determining the optimal number of parameters in deep neural networks involves balancing model capacity, generalisation ability,
and computational efficiency. It’s noted that recent work by [33] explores the use of a lightweight neural network to train WENO
schemes. Their MLP consists of three hidden layers with four neurons each, resulting in only 105 trainable parameters–significantly
fewer than the 658 parameters in the current network. However, our primary focus is on the numerical properties of the trained
WENO schemes rather than on minimising the number of parameters of neural networks. Specifically, the dissipation characteristics
and discontinuity-capturing capabilities of the neural-network-trained WENO schemes are of greater importance for the development
of our new ROUND formulation. While overparameterisation may exist, we prioritise achieving desirable numerical behaviour over
architectural minimalism.
3.3.  Training procedure of WENO3-NN
The training of WENO3-NN closely follows the methodology outlined in [1]. The training dataset comprises canonical functions
designed to emulate local solution features of hyperbolic conservation laws. A summary of the dataset is provided in Table 1. It
includes polynomials of degree up to three, jump discontinuities, sawtooth functions, and trigonometric functions. Polynomials and
hyperbolic tangent functions are evaluated over the domain [−1,1], while all other functions are evaluated over [0,1]. The dataset
is discretized with a spatial resolution of Δ𝑥=0.01. For jump and sawtooth functions, only stencils containing a discontinuity are
included. During training, the dataset is split into training and validation subsets using a validation split of 0.1. The neural network
is trained using the Adam optimizer [53] for 100 epochs, with a fixed learning rate of 10−3 and a mini-batch size of 100.
3.4.  Unifying the WENO3-NN on the normalised-variable diagram
Due to the complexity of the neural network architecture, it is not feasible to directly normalise the reconstruction operator of the
WENO3-NN scheme. To investigate the behaviour of WENO3-NN on the normalised variable diagram, we assume that the WENO3-NN
scheme exhibits a certain degree of scale independence. Using this assumption, we can generate the normalised variable diagram
for the WENO3-NN schemes by setting 𝜙̄ =0 and 𝜙̄ =1, while varying 𝜙̄ in increments of Δ𝜙̄ =0.02 as input to the trained
|     |     | 𝑖−1 𝑖+1 | 𝑖   | 𝑖   |     |
| --- | --- | ------- | --- | --- | --- |
neural network. The normalised variable diagram for WENO3-NN is presented in Fig. 2. It can be seen that when 𝜙̂ is close to 0.5,
𝑖
the WENO3-NN scheme retrieves the third-order linear scheme. This result is expected: for the smooth function without the critical
6

M. Huang, X. Deng, O.K. Matar et al. Journal of Computational Physics 555 (2026) 114764
Table 1
Summary of training functions and sampling parameters. The notation
 denotes a uniform distribution, while  denotes a Bernoulli distri-
bution.
|     |     |  Function 𝑓(𝑥)             |  Parameters           |  Number of samples   |     |
| --- | --- | -------------------------- | --------------------- | -------------------- | --- |
|     |     | ∑ 𝑛 𝑎 𝑥 𝑘                  | 𝑎 ∈  ( −1 ,1 ) ∀ 𝑘   |  4 0 0 0  for each 𝑘 |     |
|     |     | 𝑘 =0 𝑘                     | 𝑘                     |                      |     |
|     |     | 𝑢 𝑙( 𝑥 < 0 .5 )+𝑢 𝑟(𝑥>0.5) | 𝑢 ,𝑢 ∈  (− 1 0 ,1 0) |  8 0 0 0             |     |
|     |     |                            | 𝑙 𝑟                   |  4000                |     |
|     |     | (−1)𝑎𝑥+𝛿(𝑥>0.5)            | 𝑎∈(0.5),𝛿∈(0.5,2)   |                      |     |
|     |     |                            | 𝑘∈(2,20)             |  4000                |     |
sin(𝑘𝜋𝑥)
|     |     | tanh(𝑘𝑥) | 𝑘∈(5,30) |  4000 |     |
| --- | --- | -------- | --------- | ----- | --- |
Fig. 2. The 𝜙̂ −𝜙̂  relationship of WENO3-NN schemes on the normalised-variable diagram. (a) mainly shows the range of 𝜙̂ >0, and (b)
| 𝑖   | 𝑖+1∕2 |     |     | 𝑖   |     |
| --- | ----- | --- | --- | --- | --- |
mainly shows 𝜙̂ <0.
𝑖
| point, the Taylor expansion at 𝑥 |     |  yields |     |     |     |
| -------------------------------- | --- | ------- | --- | --- | --- |
𝑐
| 1      | 𝜙′′(𝑥 )ℎ  |     |     |     |        |
| ------ | --------- | --- | --- | --- | ------ |
| 𝜙̂ = − | 𝑐 +𝑂(ℎ2). |     |     |     | (3.14) |
| 𝑖 2    | 4𝜙′(𝑥 )   |     |     |     |        |
𝑐
Moreover, for smooth functions, the ideal weight loss   dominates in Eq. (3.12), which drives WENO3-NN toward the linear third-
𝑑
order scheme. From Fig. 2(a), we can also observe that when 𝜙̂ is close to 1.0, WENO3-NN retrieves the second-order central scheme,
𝑖
whereas it retrieves the second-order upwind scheme when 𝜙̂ is close to 0.0. This behaviour is the same as that of traditional WENO
𝑖
schemes on the normalised-variable diagram, as previously shown in [24]. When 0<𝜙̂ <1, the plot of WENO3-NN in the normalised
𝑖
variable diagram predominantly lies below the line 𝜙̂ =𝑘𝜙̂, where 𝑘=2.5. According to the analysis in [54], this suggests that
𝑖+1∕2 𝑖
the CFL number should not be larger than 0.4 for WENO3-NN to strictly satisfy the convection boundedness criterion across the
discontinuities.
An important feature distinguishing WENO3-NN from other non-linear schemes developed in a compact stencil is shown in
Fig. 2(b). Fig. 2(b) shows that as 𝜙̂ decreases further below zero, the behaviour of WENO3-NN gradually deviates from that of
𝑖
the second-order upwind scheme and tends toward the second-order central scheme. According to the analysis in [24], a negative
value of 𝜙̂ indicates that cell 𝑖 contains a critical point. As 𝜙̂ decreases below zero further, the local function around cell 𝑖 becomes
| 𝑖   |     |     | 𝑖   |     |     |
| --- | --- | --- | --- | --- | --- |
less smooth. Based on Eq. (3.12), the cell-face reconstruction loss  becomes more dominant than the ideal weight loss  , which
|     |     |     | 𝑟   |     | 𝑑   |
| --- | --- | --- | --- | --- | --- |
drives the WENO3-NN scheme toward the second-order central scheme in order to minimize reconstruction errors. Consequently, the
spectral property analysis reveals that WENO3-NN has better dissipation properties compared to other non-linear schemes developed
in a compact stencil.
4.  Learning ROUND operators from neural networks
Although WENO3-NN demonstrates improved accuracy and dissipation characteristics compared to existing non-linear schemes
developed using compact stencils, it incurs a higher computational cost than conventional numerical methods. To address this, we
7

M. Huang, X. Deng, O.K. Matar et al. Journal of Computational Physics 555 (2026) 114764
Fig. 3. The normalised operator of the ROUND scheme in the normalised-variable diagram.
adopt the ROUND methodology in this section and design the reconstruction operator directly on the unified normalised-variable
diagram. The proposed ROUND schemes aim to retain the key features of WENO3-NN on the normalised-variable diagram while
remaining algorithmically simple and computationally efficient.
4.1.  Formulation of improved ROUND schemes
Based on the behaviour of the WENO3-NN scheme observed in the normalised-variable diagram, we design piecewise linear
functions using the normalised variable to preserve the key characteristics of WENO3-NN while maintaining algorithmic simplicity.
The resulting ROUND scheme, which retains the essential features of WENO3-NN, is formulated as
|     | ⎧max{min{1 | 𝜙̂  | +1 ,−3 | 𝜙̂ −9 },3 | 𝜙̂ }, | if 𝜙̂ ≤0, |     |     |     |
| --- | ---------- | --- | ------ | --------- | ----- | --------- | --- | --- | --- |
|     |            | 2   | 𝑖 2 2  | 𝑖 5       | 2 𝑖   | 𝑖         |     |     |     |
̂𝑅 𝑂𝑈𝑁𝐷[𝜙̂ ]= ⎪ min{min{5 𝜙̂ ,1 +5 𝜙̂ }, 3 𝜙̂ +4 7 }, if 0<𝜙̂ ≤1, (4.1)
| 𝑖   | 𝑖 ⎨  |     | 𝑖   | 𝑖 𝑖 |     | 𝑖          |     |     |     |
| --- | ---- | --- | --- | --- | --- | ---------- | --- | --- | --- |
|     |      | 2   | 3 6 | 5 0 | 5 0 |            |     |     |     |
|     | ⎪1𝜙̂ | +1, |     |     |     | otherwise. |     |     |     |
⎩2 𝑖 2
The behaviour of the above ROUND formulation on the normalised variable diagram is shown in Fig. 3. Similar to the WENO3-NN
scheme, the proposed ROUND formulation recovers the third-order linear scheme when 𝜙̂ is close to 0.5, tends toward the second-
𝑖
order central scheme as 𝜙̂ approaches 1.0, and approaches the second-order upwind scheme when 𝜙̂ is near 0.0. The plot of the
𝑖 𝑖
ROUND scheme is also constrained to lie below the line 𝜙̂ =2.5𝜙̂ for 0<𝜙̂ <1. According to the analysis in [54], this constraint
|     |     |     |     |     | 𝑖+1∕2 | 𝑖   | 𝑖   |     |     |
| --- | --- | --- | --- | --- | ----- | --- | --- | --- | --- |
ensures that the convection boundedness criterion is strictly satisfied for a maximum CFL number of 0.4. More importantly, the
ROUND scheme deviates from the second-order upwind scheme and tends toward the second-order central scheme as 𝜙̂ decreases
𝑖
further below zero. This behaviour mirrors a key feature of the WENO3-NN scheme, which helps reduce reconstruction errors for
less smooth regions. It is noteworthy that the ROUND scheme employs a greater proportion of the second-order central scheme
than WENO3-NN when 𝜙̂ <0. This design choice is motivated by numerical experiments, indicating that incorporating more central
𝑖
schemes in less smooth regions can further enhance accuracy and improve dissipation properties.
4.2.  Formulation of low-dissipative ROUND schemes
The numerical error analysis presented in [24], along with the behaviour of the THINC scheme on the normalised-variable diagram,
demonstrates that the dissipation in ROUND schemes can be decreased by introducing controlled anti-dissipation errors. Therefore,
we propose the following formulation to further mitigate numerical dissipation while maintaining accuracy.
|                 | ⎧          | m a x { m i n | { 1 𝜙 ̂ + 1 | , − 3 𝜙̂ −   | 9 }, 3 𝜙 ̂ | } ,         | i f   | 𝜙 ̂ ≤ 0 , |       |
| --------------- | ---------- | ------------- | ----------- | ------------ | ---------- | ----------- | ----- | --------- | ----- |
|                 |            |               | 2 𝑖 2       | 2 𝑖          | 5 2        | 𝑖           |       | 𝑖         |       |
|                 | ⎪          | 5 ̂           | 1           | 5 𝜙̂         |            | 5 𝜙̂        |       | 𝜙̂ 1      |       |
|                 | ⎪          | m in { 𝜙      | , 𝜔 ( +     | ) +          | (1 − 𝜔 )   | ( )}        | i f   | 0 < ≤ ,   |       |
| ̂𝐿𝐷−𝑅𝑂𝑈𝑁𝐷[𝜙̂]= |            | 2 𝑖           | 𝑖 ,1 3      | 6 𝑖          | 𝑖, 1       | 2 𝑖         |       | 𝑖 2       | (4.2) |
| 𝑖               | 𝑖 ⎨min{3𝜙̂ |               | +47,𝜔       | (1+5𝜙̂)+(1−𝜔 |            | )(3𝜙̂ +47)} | if 1  | <𝜙̂ ≤1,   |       |
|                 | ⎪          | 50            | 𝑖 50 𝑖,2    | 3 6          | 𝑖          | 𝑖,2 50 𝑖 50 |       | 2 𝑖       |       |
⎪1𝜙̂ +1
otherwise.
|     | ⎩2  | 𝑖 2 |     |     |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
8

M. Huang, X. Deng, O.K. Matar et al. Journal of Computational Physics 555 (2026) 114764
Fig. 4. The normalised operator of the low-dissipative ROUND scheme in the normalised-variable diagram.
Table 2
Comparisons of floating point operations between two
ROUND schemes and four other WENO-type schemes.
Scheme Floating Point Operations
WENO3-JS 19
WENO5-JS 55
WENO3-NN [1] 2139
WENO3-NN-Rational [33] 508
ROUND 11
LD ROUND 20
Here, 𝜔 =𝜔(𝜙̂) is a general multi-quadric function with respect to 𝜙̂, which is defined as
𝑖 𝑖 𝑖 𝑖
𝜔 = 1 , 𝑘=0,1, (4.3)
𝑖,𝑘 ( )2
1+𝜃 (𝜙̂ −1)4
𝑘 𝑖 2
where 𝜃 is the shape parameter of the general multi-quadric function and controls the magnitude of the introduced anti-dissipation
𝑘
errors. Based on the analysis in the previous work [24], a larger value of 𝜃 corresponds to a greater magnitude of anti-dissipation
errors. In this work, we choose 𝜃 =180 and 𝜃 =600 in order to minimise dissipation while avoiding excessive anti-dissipation
1 2
errors. It should be noted that using values of 𝜃 beyond this is not recommended, as it may introduce excessive anti-dissipation
errors and potentially lead to numerical instability. We term the ROUND formulation of Eq. (4.2) as LD-ROUND (Low-Dissipative
ROUND). The behaviour of LD-ROUND on the normalised-variable diagram is presented in Fig. 4. Unlike the ROUND formulation in
Eq. (4.1), a portion of the LD-ROUND plot lies above the line representing the third-order linear scheme for 0<𝜙̂ <1. As shown by
𝑖
the numerical error analysis in [24] and the behaviour of the THINC scheme on the normalised-variable diagram, the region above
the line of the third-order linear scheme corresponds to anti-dissipation. Thus, LD-ROUND can produce low-dissipative simulation
results by introducing the controlled anti-dissipation errors for 0<𝜙̂ <1 while retaining key features of the WENO3-NN scheme for
𝑖
𝜙̂ <0. The low-dissipation property of LD-ROUND will be demonstrated in the following numerical tests.
𝑖
5. Properties of the proposed schemes
5.1. Floating point operations
As illustrated in Table 2, we compare the floating point operations (FLOPs) of the two ROUND schemes with four existing WENO-
type schemes. Following the previous work in [1], the FLOPs for the WENO3-JS, WENO5-JS, and WENO3-NN schemes are provided by
JAX. The ROUND scheme requires approximately half the FLOPs of WENO3-JS and only 5% of the cost of WENO3-NN. In contrast,
the LD-ROUND scheme demands twice the computational cost of the basic ROUND scheme, making it comparable to WENO3-JS
while still being nearly half the cost of WENO5-JS and less than 10% of the cost of WENO3-NN. To provide a more comprehensive
9

M. Huang, X. Deng, O.K. Matar et al. Journal of Computational Physics 555 (2026) 114764
Fig. 5. The spectral property of the ROUND and LD-ROUND scheme: (a) and (b) show the real and imaginary parts of the modified wavenumber,
respectively, plotted against the reduced wavenumber, following [55].
comparison, we have included the WENO scheme trained using a lightweight neural network, as explored in the recent work by [33].
Their network architecture consists of three hidden layers with four neurons each, resulting in only 105 trainable parameters and
508 floating-point operations (FLOPs). However, it is important to note that the FLOPs of the lightweight network still exceed those
of classical WENO schemes and the newly formulated ROUND schemes. It’s also noted that this comparison is limited to currently
available architectures. There remains potential for further reducing FLOPs by optimising the neural network structure, which could
lead to even more efficient models.
5.2. Spectral property
The approximate dispersion relation (ADR) method described in [55] is applied to study the spectral property of the proposed
scheme. As illustrated in Fig. 5, the spectral behaviour of the ROUND and LD-ROUND schemes are compared against representa-
tive conventional WENO5-JS and ML-assisted WENO3-NN schemes. The real and imaginary components of the complex modified
wavenumber (left and right panels, respectively) characterise dispersion and dissipation properties.
Regarding dispersion characteristics, all schemes exhibit comparable performance at low wavenumbers. However, both ROUND
schemes demonstrate reduced dispersive errors at high wavenumbers relative to WENO5-JS and WENO3-NN. In dissipation analysis,
the ROUND scheme exhibits significantly lower numerical dissipation than both reference schemes in the high-wavenumber regime.
The LD-ROUND scheme manifests dissipation behaviour similar to WENO5-JS at low wavenumbers while demonstrating superior
dissipation characteristics at elevated wavenumbers.
As established in [1], the WENO3-NN scheme achieves near-zero dissipation near the cutoff wavenumber due to its data-driven
formulation. Capitalising on this fundamental insight, the proposed two ROUND schemes attain enhanced dissipation properties,
similarly approaching vanishing dissipation near the cutoff wavenumber while maintaining favourable dispersion properties.
6. Numerical experiments
In this section, we perform a wide range of benchmark tests to evaluate and compare the numerical performance of the ROUND,
LD-ROUND, and other classical high-order WENO schemes. Unless otherwise specified, the CFL number is set to 0.4, and HLLC is
used as the approximate Riemann solver for all benchmark tests.
6.1. Accuracy tests
The accuracy of the proposed schemes is initially assessed through the solution of the linear advection equation with smooth
initial profiles. Numerical errors and convergence rates are quantified based on these solutions, using representative smooth initial
profiles:
𝜙 (𝑥)=sin(𝜋𝑥), 𝑥∈[−1,1]. (6.1)
0
10

M. Huang, X. Deng, O.K. Matar et al. Journal of Computational Physics 555 (2026) 114764
Table 3
Numerical errors, convergence rates and CPU cost.
|     |  Schemes |  Meshes 𝐿 errors | 𝐿 order |  errors  order  |  CPU cost (𝑠) |
| --- | -------- | ---------------- | ------- | --------------- | ------------- |
|     |          |                  |         | 𝐿 𝐿             |               |
|     |          | 1                | 1       | ∞ ∞             |               |
|     |          |  320 2.38×10−4   |  -      | 2.95×10−3  -    |  2.6          |
|     | WENO3-NN |  640 4.40×10−5   |  2.44   | 9.48×10−4  1.64 |  9.6          |
|     |          |  1280 8.37×10−6  |  2.39   | 3.12×10−4  1.60 |  37.0         |
|     |          |  320             |  -      |  -              |  0.6          |
|     |          | 1.72×10−4        |         | 2.27×10−3       |               |
|     | ROUND    |  640 3.25×10−5   |  2.40   | 7.35×10−4  1.63 |  2.2          |
|     |          |  1280 6.05×10−6  |  2.43   | 2.38×10−4  1.63 |  8.5          |
|     |          |  320             |  -      |  -              |  0.7          |
|     |          | 1.15×10−4        |         | 9.39×10−4       |               |
|     | LD-ROUND |  640 1.92×10−5   |  2.58   | 2.82×10−4  1.74 |  2.8          |
|     |          |  1280 3.21×10−6  |  2.58   | 8.58×10−5  1.72 |  11.3         |
Fig. 6. Variations of 𝐿  and 𝐿  errors with CPU cost. The number is calculated in a log scale.
1 ∞
| We define the 𝐿  and 𝐿 |  errors as |     |     |     |     |
| ---------------------- | ---------- | --- | --- | --- | --- |
| 1                      | ∞          |     |     |     |     |
𝑁
| 1 ∑        | |𝜙𝑠 −𝜙𝑒 | |𝜙𝑠 −𝜙𝑒       |     |     |     |
| ---------- | ------- | ------------- | --- | --- | --- |
| 𝐿 1,error= | 𝑖|, 𝐿   | ∞,error= m ax | 𝑖|. |     |     |
| 𝑁          | 𝑖       | 𝑖= 1,… ,𝑁 𝑖   |     |     |     |
𝑖=1
where 𝜙𝑠 and 𝜙𝑒 denote the numerical and exact solutions of cell 𝑖, respectively. The simulation is carried out up to time 𝑡=20. The
𝑖 𝑖
corresponding numerical errors, convergence rates and CPU costs are reported in Table 3. Compared with the WENO3-NN scheme,
the ROUND and LD-ROUND schemes exhibit smaller 𝐿  and 𝐿  errors. Specifically, the ROUND scheme achieves a comparable
1 ∞
convergence order, yielding better accuracy than the WENO3-NN scheme. The CPU cost of the ROUND scheme is approximately
25% of that incurred by the WENO3-NN scheme. Moreover, the LD-ROUND scheme demonstrates superior performance in both
convergence order and accuracy relative to the ROUND and WENO3-NN schemes. It shows that using the same grid, LD-ROUND
produces numerical errors and CPU costs that are both less than half of those generated by WENO3-NN.
To further evaluate computational efficiency, we analyse how numerical error changes with CPU time. The results are shown in
Fig. 6. These figures illustrate that both proposed schemes substantially improve computational efficiency compared to the WENO3-
NN scheme. Although the ROUND scheme has the lowest computational cost, it tends to produce higher numerical errors compared
to the LD-ROUND scheme.
6.2.  Advection of complex profiles
To evaluate the capability of the proposed schemes in handling varying degrees of smoothness, we consider the propagation of
a complex wave [23]. This benchmark is commonly used to assess both the non-oscillatory behaviour and numerical dissipation of
11

M. Huang, X. Deng, O.K. Matar et al. Journal of Computational Physics 555 (2026) 114764
Fig. 7. Numerical solutions for advection of complex wave at 𝑡=2.0 with 200 uniform cells.
shock-capturing schemes. The initial condition for the advected scalar field is given as in [23]:
| 1                                   |     |     | ≤0.1,      |     |     |
| ----------------------------------- | --- | --- | ---------- | --- | --- |
| ⎧ (𝐺(𝑥,𝛽,𝑧−𝛿)+𝐺(𝑥,𝛽,𝑧+𝛿)+4𝐺(𝑥,𝛽,𝑧)) |     |     | if |𝑥+0.7| |     |     |
6
| ⎪1  |     |     | if |𝑥+0.3| ≤0.1, |     |     |
| --- | --- | --- | ---------------- | --- | --- |
⎪
| 𝜙 (𝑥)= 1−|10𝑥−1| |     |     | if |𝑥−0.1| ≤0.1, | ,   | (6.2) |
| ---------------- | --- | --- | ---------------- | --- | ----- |
0 ⎨
| ⎪1(𝐹(𝑥,𝜇,𝑎−𝛿)+𝐹(𝑥,𝜇,𝑎+𝛿)+4𝐹(𝑥,𝜇,𝑎)) |     |     | if |𝑥−0.5| ≤0.1, |     |     |
| ----------------------------------- | --- | --- | ---------------- | --- | --- |
⎪6
| ⎩0  |     |     | otherwise |     |     |
| --- | --- | --- | --------- | --- | --- |
where functions 𝐹 and 𝐺 are defined as
√
| 𝐺(𝑥,𝛽,𝑧)=𝑒−𝛽(𝑥−𝑧)2,𝐹(𝑥,𝜇,𝑎)= |     | ( 1−𝜇2(𝑥−𝑎)2,0 | )   |     | (6.3) |
| ---------------------------- | --- | -------------- | --- | --- | ----- |
|                              |     | max            | ,   |     |       |
and the coefficients are
|               |                  | (     | 36𝛿2) |     | (6.4) |
| ------------- | ---------------- | ----- | ----- | --- | ----- |
| 𝑎=0.5, 𝑧=0.7, | 𝛿=0.005, 𝜇=10.0, | 𝛽=log | .     |     |       |
2
The computational domain is defined as [−1,1] with periodic boundary conditions. Numerical results at time 𝑡=2.0 are presented
using a 200-cell mesh, as shown in Fig. 7.
In Fig. 7(a), we compare the performance of the proposed ROUND scheme with that of the WENO3-NN scheme. Both schemes
effectively resolve discontinuities while suppressing spurious oscillations. The ROUND scheme delivers results comparable to or even
superior to those of the WENO3-NN scheme, particularly for initial profiles with local maxima, such as the Gaussian, triangular,
and semi-elliptic functions. In Fig. 7(b), the LD-ROUND scheme exhibits significantly enhanced resolution of local maxima due to its
low-dissipation design. Furthermore, it shows slightly better performance in capturing the square wave discontinuity compared to
both the WENO3-NN and ROUND schemes.
6.3.  Lax problems
To check the ability of the proposed numerical scheme to capture relatively strong shock, we solve the Lax problem, originally
introduced in [56], with the following initial condition:
{
|             | (0.445,0.698,3.528) | 0⩽𝑥⩽0.5 |     |     |       |
| ----------- | ------------------- | ------- | --- | --- | ----- |
| (𝜌 ,𝑢 ,𝑝 )= |                     |         |     |     | (6.5) |
| 0 0 0       | (0.5,0,0.571)       | 0.5<𝑥⩽1 |     |     |       |
The computation is performed using 200 uniformly spaced cells. The numerical solutions obtained from the proposed schemes at
𝑡=0.15 are presented in Fig. 8, compared with the WENO3-NN scheme. The proposed schemes effectively capture the contact dis-
continuity and accurately resolve the shock wave without producing numerical oscillations. As illustrated in Fig. 8(a), the ROUND
scheme performs comparably to the WENO3-NN method. In comparison, the LD-ROUND scheme shown in Fig. 8(b) resolves the
contact discontinuity with slightly less numerical diffusion.
12

M. Huang, X. Deng, O.K. Matar et al. Journal of Computational Physics 555 (2026) 114764
Fig. 8. Numerical solution for the Lax problem at 𝑡=0.16 with 200 uniform cells.
6.4. Shock-entropy wave interactions
This test case examines the interaction between a Mach 3 shock wave and a localised density perturbation, producing a complex
flow field that features both smooth regions and sharp discontinuities. Commonly known as the shock-density wave interaction
problem [57], it serves as a standard benchmark for evaluating the performance of numerical schemes in capturing mixed flow
features. The initial conditions are specified as follows:
{
(3.857148,2.629369,10.333333) 0⩽𝑥⩽0.1
(𝜌 ,𝑢 ,𝑝 )= (6.6)
0 0 0 (1+0.2sin(50𝑥−25),0,1) 0.1<𝑥⩽1
The computational domain is defined as [0,1]. The computations are performed using 400 uniform cells and are integrated up to
a time of 𝑡=0.18. As shown in Fig. 9, the exact solution, indicated by the black solid line, is computed using the classical 5th-order
WENO5-JS scheme with 2000 mesh cells. It is evident that the flow structures produced by the WENO3-NN scheme are significantly
smeared out, whereas the proposed schemes are proficient in capturing small-scale waves with higher accuracy. We further compare
the performance of the ROUND and LD-ROUND schemes in Fig. 9(a) and (b). These comparisons reveal that the LD-ROUND scheme
resolves density perturbations more effectively and captures the peaks of the waves with greater accuracy.
6.5. Blast waves
The interacting blast waves problem, introduced in [58], presents a demanding test scenario involving multiple interactions
between strong shock and rarefaction waves. In this study, it is used to assess the effectiveness of the proposed schemes in resolving
intense shocks and intricate flow features. The initial condition is specified as follows:
⎧(1, 0, 1000), if 0≤𝑥<0.1,
(𝜌 , 𝑢 , 𝑝 )= ⎪ (1, 0, 0.01), if 0.1≤𝑥<0.9, (6.7)
0 0 0 ⎨
⎪(1, 0, 100), if 0.9≤𝑥<1.
⎩
Reflective boundary conditions are applied at both ends of the computational domain. The simulation uses 400 uniformly spaced
cells, with time integration carried out to 𝑡=0.038. The reference solution, depicted as a solid line, is obtained using the WENO5
scheme on a highly refined mesh. According to the literature [1], the WENO3-NN scheme provides more accurate results compared
to the traditional WENO3-JS scheme. As shown in Fig. 10, the proposed schemes achieve a performance that is comparable to, and
in some aspects better than, the WENO3-NN scheme. Notably, the ROUND scheme resolves a density profile that is similar to the
WENO3-NN scheme, as illustrated in Fig. 10(a), while requiring less computational cost. Although the LD-ROUND scheme requires
slightly higher costs than the ROUND scheme, it delivers low-dissipative results, particularly noticeable near the valley and at the
right density peak in Fig. 10(b).
13

M. Huang, X. Deng, O.K. Matar et al. Journal of Computational Physics 555 (2026) 114764
Fig. 9. Numerical solutions for the shock density wave interaction at 𝑡=0.18 with 400 uniform cells.
Fig. 10. Numerical solutions for blast-wave interaction problem with 400 uniform cells.
6.6.  2D Riemann problem
This work simulates test case 3 of the 2D Riemann problems described in [59]. The computational domain is set as [0,1]×[0,1].
The initial condition is specified as follows:
⎧(1.5,
|     | 0.0, | 0.0, 1.5) | [0.8,1]×[0.8,1], |     |
| --- | ---- | --------- | ---------------- | --- |
⎪
| (         | ) ⎪(0.5323,  | 1.206, 0.0, 0.3)     | [0.0,0.8]×[0.8,1],   |       |
| --------- | ------------ | -------------------- | -------------------- | ----- |
| 𝜌 , 𝑢 , 𝑣 | , 𝑝 =        |                      |                      | (6.8) |
| 0 0       | 0 0 ⎨(0.138, | 1.206, 1.206, 0.029) | [0.0,0.8]×[0.0,0.8], |       |
⎪
|     | ⎪(0.5323, | 0.0, 1.206, 0.3) | [0.8,1.0]×[0.0,0.8], |     |
| --- | --------- | ---------------- | -------------------- | --- |
⎩
The simulation is set to run up until 𝑡=0.8 with 600×600 uniform meshes.
14

M. Huang, X. Deng, O.K. Matar et al. Journal of Computational Physics 555 (2026) 114764
Fig. 11. Density contours of 2D Riemann problem from ROUND(upper left), LD ROUND(upper right), WENO3-NN(lower left) and WENO5-JS(lower-
right) at the time 𝑡=0.8 with resolution of 600×600. The figure is plotted with 43 density contours between 0.135 to 1.75.
As illustrated in Fig. 11, the results of the two proposed schemes are shown in the upper panels, while the WENO5-JS and WENO3-
NN schemes are presented in the lower panels. The ability to resolve small-scale structures is highly sensitive to the numerical
dissipation introduced by each scheme. It can be observed that the ROUND scheme captures slightly fewer but still comparable
vortical structures relative to the WENO3-NN scheme. Due to its higher-order accuracy, the WENO5-JS scheme demonstrates lower
numerical dissipation, allowing it to resolve more intricate flow features. Notably, the proposed LD-ROUND scheme reproduces even
finer small-scale flow structures than the fifth-order WENO5-JS scheme. Furthermore, it is important to emphasise that both proposed
schemes maintain flow symmetry throughout the simulation.
6.7. Double Mach reflection
The double Mach reflection problem described in [58] is simulated here. This benchmark is always adopted to evaluate the ability
of capturing the shock waves as well as resolving the small-scale vortical structure along the slip line in the re-circulation zone, which
is generated by the Kelvin-Helmholtz instabilities. The initial condition is given as
( 𝜌 , 𝑢 , 𝑣 , 𝑝 ) = { (8.0, 7.145, −4.125, 116.8333) 𝑦≥1.732(𝑥−0.1667) (6.9)
0 0 0 0 (1.4, 0.0, 0.0, 1.0) 𝑦<1.732(𝑥−0.1667)
15

M. Huang, X. Deng, O.K. Matar et al. Journal of Computational Physics 555 (2026) 114764
Fig. 12. Numerical solution for the double Mach 10 refection problem. The figure shows the density contour with 30 contours between 1.4 and
22.4. The zoomed perspective of the re-circulation zone of the original WENO5, WENO3-NN, ROUND and LD-ROUND schemes are presented with
the grid resolution of Δ𝑥=Δ𝑦=1∕400.
The computational domain is [0,4]×[0,1] and the final simulation time is 𝑡=0.2. The solutions are computed by the uniform grid size,
Δ𝑥=Δ𝑦=1∕400. Initially, a right-moving Mach 10 shock is imposed at 𝑥=0.1667 with 60◦ angle relative to 𝑥-axis. At the bottom,
a reflecting wall boundary is prescribed from 𝑥=0.1667 to 𝑥=4, while a post-shock condition is enforced from 𝑥=0 to 𝑥=0.1667.
In this simulation, the HLL Riemann solver [42] is employed to mitigate the carbuncle phenomenon, as the HLLC solver is known to
exhibit numerical instabilities in the presence of strong shock waves [60].
As shown in Fig. 12, the density contour of the zoomed perspective of the re-circulation zone calculated by the original WENO5-JS,
WENO3-NN and two proposed ROUND schemes are presented. From the numerical results, the proposed schemes are able to capture
the shocks without spurious oscillation. The WENO3-NN, WENO5-JS, and ROUND schemes produce comparable results, whereas the
LD-ROUND scheme resolves the more complex vortical structures along the primary slip line than the other three schemes.
6.8. Rayleigh-Taylor instability
The Rayleigh-Taylor instability (RTI), originally investigated by Xu and Shu [61], serves as a fundamental benchmark for evaluat-
ing numerical methods for the compressible Euler equations [4,62]. Numerical schemes with lower dissipation are expected to yield
16

M. Huang, X. Deng, O.K. Matar et al. Journal of Computational Physics 555 (2026) 114764
Fig. 13. Rayleigh-Taylor instability at t = 1.95. 43 equally spaced density contours between 0.9 and 2.2 are presented. Resolution 128×512.
richer and more detailed flow structures in the simulation of inviscid Euler equations. The initial condition is specified as follows:
{ (2,0,𝑣(𝑥),2𝑦+1), 𝑦≥0.5,
(𝜌 ,𝑢 ,𝑣 ,𝑝 )= (6.10)
0 0 0 0 (1,0,𝑣(𝑥),𝑦+1.5), otherwise.
Here, 𝑣(𝑥) represents the initial velocity in the 𝑦-direction, perturbed in the 𝑥-direction to initiate mixing between the heavy and light
fluids. Following the definition in [61], 𝑣(𝑥) is given by:
𝑣(𝑥)=−0.025𝑐cos(8𝜋𝑥),
where the sound speed is defined as 𝑐= √ 𝛾𝑝∕𝜌, with the ratio of specific heats set to 𝛾=5∕3. The computational domain is [0,0.25]×
[0,1]. Reflective boundary conditions are imposed on the left and right boundaries. At the bottom and top boundaries, fixed primitive
states are prescribed as (𝜌,𝑢,𝑣,𝑝)=(1,0,0,2.5) and (𝜌,𝑢,𝑣,𝑝)=(2,0,0,1), respectively. According to the definition in Eq. (6.10), the
initial interface is located at 𝑦=0.5, separating the heavy and light fluids and forming a contact discontinuity. The simulation is
performed on a uniform grid with resolution 128×512, corresponding to a uniform grid spacing of Δ𝑥=Δ𝑦=1∕512.
Fig. 13 shows the solution at 𝑡=1.95 with 43 equally spaced density contours ranging from 0.9 to 2.2. Numerical results for all
three-stencil WENO variants–including WENO3-JS, WENO3-Z, and WENO3-NN1/NN2–have been previously reported in [1]. In this
example, we further compare the proposed schemes with five-stencil WENO methods, specifically WENO5-JS and WENO5-Z [44].
The ROUND scheme demonstrates improved performance compared to the WENO5-JS method, successfully capturing significantly
finer flow structures, although it is still slightly less detailed than the WENO5-Z scheme. In contrast, the LD-ROUND scheme benefits
from a low-dissipation design, allowing it to resolve mushroom-shaped structures with greater complexity and capture more intricate
vortical patterns. These results indicate a lower level of dissipation than that of both the three-stencil schemes reported in the study
by Bezgin et al.[1] and the five-stencil WENO5-JS/Z schemes. Notably, both the ROUND and LD-ROUND schemes reproduce the
upturned tail of the light fluid, which is typically observed only in high-order schemes, such as those discussed by [4] and [63]. This
highlights their effectiveness in resolving subtle flow features.
6.9. Richtmyer-Meshkov instability
The Richtmyer-Meshkov instability (RMI), first analytically investigated by Richtmyer[64] and experimentally studied by
Meshkov[65], is employed to evaluate the capability of the proposed schemes in resolving small-scale flow structures near dis-
continuities, which are particularly sensitive to numerical dissipation. The RMI originates from an impulsive acceleration acting on a
contact discontinuity, typically induced by a shock wave interaction. The computational domain is defined as [−0.5,0.5]×[0,4], and
17

M. Huang, X. Deng, O.K. Matar et al. Journal of Computational Physics 555 (2026) 114764
Fig. 14. Numerical solution for the Richtmyer-Meshkov instability at time 𝑡=2.
the initial conditions are specified as follows:
|     | ⎧(1.0,0.0,0.0,1.0) | if 𝑦>=2−0.5cos(2𝜋𝑥), |
| --- | ------------------ | -------------------- |
⎪
| (𝜌 ,𝑢 ,𝑣 ,𝑝 | )= (0.25,0.0,0.0,1.0) | if 0.5≤𝑦≤2−0.5cos(2𝜋𝑥), |
| ----------- | --------------------- | ----------------------- |
| 0 0 0       | 0 ⎨                   |                         |
|             | ⎪(8∕3,0.0,0.0,4.5)    | otherwise,              |
⎩
The simulation is carried out to time 𝑡=2.0 with a uniform grid size of Δ𝑥=Δ𝑦=1∕200. Wall boundary conditions are imposed on
the left and right boundaries, while Dirichlet boundary conditions are applied at the top and bottom of the computational domain.
The comparison of density contours obtained using the proposed two schemes is illustrated in Fig. 14. These contours are compared
with those from the three-stencil WENO3-NN scheme and the five-stencil WENO5-JS scheme. The ROUND scheme shows comparable
capabilities for capturing shocks and vortices when compared to the WENO3-NN scheme, with only minor differences noted in the
vortices. In contrast, the LD-ROUND scheme not only effectively resolves contact discontinuities but also captures finer and more
intricate vortical structures compared to the WENO5-JS scheme. These results indicate that the ROUND scheme outperforms the
WENO3-NN scheme in terms of overall resolution quality, achieving similar results at a lower computational cost. Meanwhile, the
LD-ROUND scheme demonstrates low dissipation and superior capability in resolving small-scale flow features.
6.10.  Shock vortex interaction
This test case, first proposed in [23], investigates the interaction between a stationary shock wave and a vortex. A Mach 1.1 shock
wave is initially positioned at 𝑥=0.5 and is oriented normally to the 𝑥-axis. Based on the pre-shock conditions to the left of the
shock, the post-shock conditions to the right are determined using the Rankine-Hugoniot relations. Additionally, a vortex, centred at
(𝑥 ,𝑦 )=(0.25,0.5), is introduced as a perturbation in the initial flow field, which is given as
𝑐 𝑐
| ⎧Δ𝑢=𝜖𝜏𝑒𝛼(  | 1−𝜏2) sin𝜃 |     |
| ---------- | ---------- | --- |
| ⎪Δ𝑣=−𝜖𝜏𝑒𝛼( | 1−𝜏2)      |     |
cos𝜃
| ⎪           | ( 1−𝜏2) |     |
| ----------- | ------- | --- |
| ⎨ =−(𝛾−1)𝜖2 | 2 𝛼     |     |
| ⎪Δ𝑇         | 𝑒       |     |
4 𝛼 𝛾
⎪
⎩Δ𝑆=0
Here, 𝑇 and 𝑆 denote the temperature and entropy, respectively. The parameter 𝜏 is calculated by
| 𝑟       | √            |              |
| ------- | ------------ | ------------ |
| 𝜏= , 𝑟= | (𝑥−𝑥 )2+(𝑦−𝑦 | )2, 𝑟 =0.05. |
|         | 𝑐            | 𝑐 𝑐          |
𝑟 𝑐
𝜖 indicates the strength of the vortex and 𝛼 represents the decay rate of the vortex. As in the same definition as [23], we take the
value of 𝜖=0.3 and 𝛼=0.204. The whole computational region is [0,2]×[0,1] with a uniform grid size of ℎ=1∕100.
As illustrated in Fig. 15, we compare the pressure profiles computed by the proposed schemes with those obtained using two
WENO-type schemes: WENO3-NN and WENO5-JS. The reference solution is generated with a much finer mesh using the WENO5-
JS scheme. It can be observed that the ROUND scheme yields more accurate results at the shock front compared to the other two
18

M. Huang, X. Deng, O.K. Matar et al. Journal of Computational Physics 555 (2026) 114764
Fig. 15. Pressure distribution along the 𝑦=0.5. The comparison of the ROUND (red circle), LD-ROUND(blue diamond) with WENO5-JS(green
triangle), WENO5-JS(brown square) and reference solution (solid line). (For interpretation of the references to colour in this figure legend, the
reader is referred to the web version of this article.)
Fig. 16. Schlieren images of the flow structure of the supersonic jet. Each time instant compares the WENO5-JS scheme with the proposed two
schemes. PS: Primary shock, VS: Vortex-induced shock, MV: Main vortex, ML: Main layer, EW: Expansion wave, KHI: Kelvin-Helmholtz instability,
TP: Triple point, R: Rearward facing shock, SS: Slipstream.
WENO-type schemes. Furthermore, the LD-ROUND scheme provides the most accurate solution near the shock front at 𝑥≈0.55 and
captures the lowest pressure value around 𝑥≈0.95, demonstrating its superior accuracy.
6.11. Supersonic planar jet
In this benchmark test, a simulation of a supersonic planar jet is performed to evaluate the ability of the proposed scheme to capture
complex flow characteristics, including shock waves, shear layers, and vortex structures [40,66]. The initial setup follows Case 3 in
[66]. The calculation is carried out on a uniform grid with a mesh spacing of Δ𝑥=Δ𝑦=1.5625×10−4. At the inlet, a fixed flow
condition is imposed on the central place of the left boundary, specified as (𝜌,𝑢,𝑣,𝑝)=(1.625kg∕m3, 486.14m∕s, 0m∕s, 141855Pa),
corresponding to a Mach 1.4 under-expanded jet flow.
19

M. Huang, X. Deng, O.K. Matar et al. Journal of Computational Physics 555 (2026) 114764
Fig. 17. Schlieren images of the density fields for two-dimensional Mach 6 air-helium interaction at time 𝑡=0.006,0.012,0.018,0.024. In each time
instant, the upper half represents the ROUND or LD-ROUND scheme while the lower half denotes the result of WENO5-JS.
Fig. 16 illustrates the evolution of the inlet supersonic jet. At each time instant, the numerical solutions obtained by the two
proposed schemes are shown in the upper half of each panel, while the results from the WENO5-JS scheme are presented in the
lower half. These results are also comparable to Fig. 7(a) in [66], which was computed using a large-eddy simulation (LES) with a
high-order TCD-WENO scheme [67] on a uniform grid of size 8.33×10−5.
Initially, a strong primary shock (PS) forms near the inlet boundary. Following the shock, a reflected expansion wave (EW) is
generated and propagates downstream. Subsequently, the effects of viscosity and expansion cause the mixing layer (ML) to roll up,
leading to the formation of the main vortex (MV). Due to entrainment within the cores of the primary vortices, peak velocities arise
in the recirculation region of the main vortex, where the local velocity exceeds the speed of sound, resulting in the generation of
vortex-induced shocks (VS1 and VS2). As time progresses, the mixing layer becomes unstable due to the large velocity difference
across the shear layer, which is indicative of Kelvin-Helmholtz instability (KHI) [68]. In the later stages, the rearward-facing shock
(R), vortex-induced shock (VS2), and the Mach stem merge at triple points (TP).
There are two primary interactions of particular interest. First, the interaction between the mixing layer and the expansion waves
accelerates the instability, giving rise to a series of counter-rotating vortical structures. Second, the development of vortices along the
shock wave interaction between VS1 and VS2 initially generates small vortices behind the shock waves, which subsequently evolve
into larger vortical structures. These structures are highly sensitive to the numerical dissipation of the schemes used. As shown in
Fig. 16, the ROUND scheme produces results comparable to those of WENO5-JS. Furthermore, the LD-ROUND scheme captures even
more complex and intricate vortical structures along the mixing layer and within the region between VS1 and VS2. It outperforms
the results obtained by the ROUND_A and ROUND_A+ schemes [40], as well as the TCD-WENO scheme [66], which requires a finer
grid resolution.
6.12. Mach 6 air-helium shock bubble interaction
A compressible multi-phase flow example, initially proposed by [69] describes the interaction between a right-moving Mach 6
shock wave in air and a cylindrical helium bubble. The flow dynamics are governed by the inviscid compressible Euler equations.
The interface between the two phases is captured using a sharp interface method based on the level set function [70,71]. To further
demonstrate the robustness and accuracy of the proposed schemes, both methods are also applied within the framework of a multi-
resolution approach for compressible multi-phase flow [72]. Numerical results for the same problem have been reported previously
in [72,73]. The initial conditions are specified as follows:
⎧(1,0,0,1,1.4) pre-shock air,
(𝜌,𝑢,𝑣,𝑝,𝛾)= ⎪ (5.268,5.752,0,41.83,1.4) post-shock air, (6.11)
⎨
⎪(0.138,0,0,1,1.667) helium bubble,
⎩
20

M. Huang, X. Deng, O.K. Matar et al. Journal of Computational Physics 555 (2026) 114764
and the level set function 𝜑(𝑥,𝑦,𝑡) at the beginning is defined as
√
𝜑(𝑥,𝑦,0)=−0.025+ (𝑥−0.15)2+(𝑦−0.0445)2. (6.12)
Here, the region with 𝜑>0 corresponds to air, while 𝜑<0 represents the helium bubble centered at (0.15,0.0445) with a radius of
0.025. The initial shock wave is positioned at 𝑥=0.1. The dimensional reference quantities are defined as 1 atmosphere, 1kg∕m3,
and 1m. The computational domain is [0,0.356]×[0,0.089], corresponding to an aspect ratio of 4∶1. Solid wall boundary conditions
are imposed on the upper and lower boundaries, while inflow and outflow conditions are applied at the left and right boundaries,
respectively. Different from [72], we set the maximum level of resolution to 𝐿 =3, with the finest resolution corresponding to a
max
uniform grid of 1024×256 cells.
The numerical results are depicted as density Schlieren images captured at various time intervals, as illustrated in Fig. 17. The
findings from the proposed schemes are compared with those obtained using the high-order WENO5-JS scheme. For each time instant,
the upper panel displays results from the ROUND or LD-ROUND schemes, while the lower panel shows the results from the WENO5-
JS scheme. As seen in Fig. 17, the results align well with previous studies [69,72,73]. Despite the relatively coarse resolution, the
proposed schemes produce results that are comparable to those of the WENO5-JS scheme. The vortex structures forming around the
helium bubble are sensitive to numerical dissipation. At times 𝑡=0.018 and 𝑡=0.024, the ROUND scheme demonstrates slightly more
dissipation, while the LD-ROUND scheme captures more complex vortical structures than the WENO5-JS scheme.
7. Concluding remarks
Following the preceding part, this study projects the WENO3-NN scheme into the normalised-variable framework to reveal the
mechanism behind its enhanced dissipation properties and to guide the development of a more accurate ROUND formulation. Anal-
ysis of WENO3-NN within this framework shows that as the normalised variable drops below zero, indicating regions near critical
points, the scheme transitions from a second-order upwind scheme to a second-order central scheme. With these insights, we design
a new ROUND scheme that captures the advantageous traits of WENO3-NN while amplifying the central scheme contribution to
further reduce dissipation. This formulation maintains simplicity and is extended to a low-dissipation variant, LD-ROUND, through
the introduction of controlled anti-dissipation effects. Comparative evaluations of numerical accuracy and computational efficiency
demonstrate that both ROUND and LD-ROUND significantly outperform WENO3-NN. Extensive validation on compressible single-
and multi-phase flow problems confirms the robustness of the proposed methods, with LD-ROUND achieving resolution that exceeds
those of classical fifth-order WENO schemes in several benchmark scenarios.
The proposed methodology is also applicable to variants of the WENO3-NN scheme that utilise rational neural networks [33],
as well as to other non-linear schemes employing neural networks trained within a compact stencil, such as WCNS3-MR-NN [26].
By adopting this approach, these schemes that use neural networks can benefit from simplified algorithmic structures, reduced
computational cost, and potentially improved accuracy. It is also worth noting that the proposed ROUND formulation can be extended
to unstructured grids by employing the virtual upwind point method, following the approach outlined in [38].
CRediT authorship contribution statement
Minsheng Huang: Writing – original draft, Visualization, Validation, Software, Methodology, Investigation, Formal analysis;
Xi Deng: Writing – original draft, Visualization, Validation, Supervision, Software, Methodology, Investigation, Formal analysis,
Conceptualization; Omar K. Matar: Writing – review & editing, Supervision, Resources, Project administration, Methodology, Funding
acquisition, Conceptualization; Wenjun Ying: Writing – review & editing, Supervision, Resources, Funding acquisition.
Data availability
Data will be made available on request.
Declaration of competing interest
The authors declare that they have no known competing financial interests or personal relationships that could have appeared to
influence the work reported in this paper.
Acknowledgements
XD and OKM acknowledge the funding provided by the Engineering and Physical Sciences Research Council and First Light Fusion
through the ICL-FLF Prosperity Partnership (grant number EP/X025373/1). WY acknowledges the funding provided by the National
Key R&D Program of China (Project No. 2020YFA0712000), the Shanghai Science and Technology Innovation Action Plan in Basic
Research Area (Project No. 22JC1401700) and the Fundamental Research Funds for the Central Universities. Dr. Abd Essamade Saufi
from FLF and Dr. Francesco Miniati from Mach42 are acknowledged for the fruitful discussion. MH would give special thanks to Dr.
Wenhua Ma for inspirational discussions.
21

M. Huang, X. Deng, O.K. Matar et al. Journal of Computational Physics 555 (2026) 114764
References
[1] D.A. Bezgin, S.J. Schmidt, N.A. Adams, WENO3-NN: a maximum-order three-point data-driven weighted essentially non-oscillatory scheme, J. Comput. Phys.
452 (2022) 110920.
[2] M.M. Rai, P. Moin, Direct numerical simulation of transition and turbulence in a spatially evolving boundary layer, J. Comput. Phys. 109 (2) (1993) 169–192.
[3] C.-W. Shu, High order weighted essentially nonoscillatory schemes for convection dominated problems, SIAM Rev. 51 (1) (2009) 82–126.
[4] L. Fu, X.Y. Hu, N.A. Adams, A family of high-order targeted ENO schemes for compressible-fluid simulations, J. Comput. Phys. 305 (2016) 333–359.
[5] H.T. Huynh, A flux reconstruction approach to high-order schemes including discontinuous Galerkin methods, in: 18th AIAA Computational Fluid Dynamics
Conference, 2007, p. 4079.
[6] P.E. Vincent, P. Castonguay, A. Jameson, A new class of high-order energy stable flux reconstruction schemes, J. Sci. Comput. 47 (1) (2011) 50–72.
[7] F.D. Witherden, A.M. Farrington, P.E. Vincent, PyFR: an open source framework for solving advection–diffusion type problems on streaming architectures using
the flux reconstruction approach, Comput. Phys. Commun. 185 (11) (2014) 3028–3040.
[8] P. Lesaint, P.-A. Raviart, On a finite element method for solving the neutron transport equation, Publ. Math. Inf. Rennes (S4) (1974) 1–40.
[9] H. Jasak, A. Jemcov, Z. Tukovic, et al., OpenFOAM: A C++ library for complex physics simulations, in: International Workshop on Coupled Methods in
Numerical Dynamics, vol. 1000, IUC Dubrovnik Croatia, 2007, pp. 1–20.
[10] F. Palacios, J. Alonso, K. Duraisamy, M. Colonno, J. Hicken, A. Aranake, A. Campos, S. Copeland, T. Economon, A. Lonkar, T. Lukaczyk, T. Taylor, Stanford
university unstructured (SU^2): an open-source integrated computational environment for multi-physics simulation and design, in: 51st AIAA Aerospace Sciences
Meeting Including the New Horizons Forum and Aerospace Exposition, American Institute of Aeronautics and Astronautics, 2013. https://doi.org/10.2514/6.
2013-287. https://doi.org/10.2514/6.2013-287
[11] S. Clain, S. Diot, R. Loubère, A high-order finite volume method for systems of conservation laws–multi-dimensional optimal order detection (MOOD), J.
Comput. Phys. 230 (10) (2011) 4028–4050.
[12] S. Diot, S. Clain, R. Loubère, Improved detection criteria for the multi-dimensional optimal order detection (MOOD) on unstructured meshes with very high-order
polynomials, Comput. Fluids 64 (2012) 43–63.
[13] M. Sonntag, C.-D. Munz, Efficient parallelization of a shock capturing for discontinuous Galerkin methods using finite volume sub-cells, J. Sci. Comput. 70 (3)
(2017) 1262–1289.
[14] M. Dumbser, O. Zanotti, R. Loubère, S. Diot, A posteriori subcell limiting of the discontinuous Galerkin finite element method for hyperbolic conservation laws,
J. Comput. Phys. 278 (2014) 47–75.
[15] Z.-H. Jiang, X. Deng, F. Xiao, C. Yan, J. Yu, S. Lou, Hybrid discontinuous Galerkin/finite volume method with subcell resolution for shocked flows, AIAA J.
(2021) 1–18.
[16] S. Liu, Y. Shen, Discontinuity-detecting method for a four-point stencil and its application to develop a third-order hybrid-WENO scheme, J. Sci. Comput. 81
(3) (2019) 1732–1766.
[17] A. Harten, High resolution schemes for hyperbolic conservation laws, J. Comput. Phys. 49 (1983) 357–393.
[18] A. Harten, On a class of high resolution total-variation-stable finite-difference schemes, SIAM J. Numer. Anal. 21 (1) (1984) 1–23.
[19] B.P. Leonard, Universal limiter for transient interpolation modeling of the advective transport equations: the ULTIMATE conservative difference scheme, NASA
Tech. Memorandum 100916 (1988) 115.
[20] B.P. Leonard, The ULTIMATE conservative difference scheme applied to unsteady one-dimensional advection, Comput. Methods Appl. Mech. Eng. 88 (1) (1991)
17–74.
[21] B.P. Leonard, Simple high-accuracy resolution program for convective modelling of discontinuities, Int. J. Numer. Methods Fluids 8 (10) (1988) 1291–1318.
[22] X.-D. Liu, S. Osher, T. Chan, Weighted essentially non-oscillatory schemes, J. Comput. Phys. 115 (1) (1994) 200–212.
[23] G.-S. Jiang, C.-W. Shu, Efficient implementation of weighted ENO schemes, J. Comput. Phys. 126 (1) (1996) 202–228.
[24] X. Deng, A unified framework for non-linear reconstruction schemes in a compact stencil. Part 1: beyond second order, J. Comput. Phys. 481 (2023) 112052.
[25] H. Wang, Y. Cao, Z. Huang, Y. Liu, P. Hu, X. Luo, Z. Song, W. Zhao, J. Liu, J. Sun, et al., Recent advances on machine learning for computational fluid dynamics:
a survey, arXiv:2408.12171 (2024).
[26] S. Fan, J. Qin, Y. Dong, Y. Jiang, X. Deng, WCNS3-MR-NN: a machine learning-based shock-capturing scheme with accuracy-preserving and high-resolution
properties, J. Comput. Phys. 532 (2025) 113973.
[27] Y. Feng, T. Liu, A characteristic-featured shock wave indicator on unstructured grids based on training an artificial neuron, J. Comput. Phys. 443 (2021)
110446.
[28] M. Huang, L. Cheng, W. Ying, X. Deng, F. Xiao, A low-dissipation reconstruction scheme for compressible single- and multi-phase flows based on artificial neural
networks, J. Comput. Phys. 530 (2025) 113894.
[29] X. Wen, W.S. Don, Z. Gao, J.S. Hesthaven, An edge detector based on artificial neural network with application to hybrid compact-WENO finite difference
scheme, J. Sci. Comput. 83 (3) (2020) 49.
[30] T. Kossaczká, A.D. Jagtap, M. Ehrhardt, Deep smoothness weighted essentially non-oscillatory method for two-dimensional hyperbolic conservation laws: a
deep learning approach for learning smoothness indicators, Phys. Fluids 36 (3) (2024) 036603.
[31] B. Stevens, T. Colonius, Enhancement of shock-capturing methods via machine learning, Theor. Comput. Fluid Dyn. 34 (4) (2020) 483–496.
[32] D.A. Bezgin, S.J. Schmidt, N.A. Adams, A data-driven physics-informed finite-volume scheme for nonclassical undercompressive shocks, J. Comput. Phys. 437
(2021) 110324.
[33] S. Shahane, S. Chammas, D.A. Bezgin, A.B. Buhendwa, S.J. Schmidt, N.A. Adams, S.H. Bryngelson, Y.-F. Chen, Q. Wang, F. Sha, L. Zepeda-Núñez, Rational-WENO:
a lightweight, physically-consistent three-point weighted essentially non-oscillatory scheme, arXiv:2409.09217 (2024).
[34] X. Nogueira, J. Fernández-Fidalgo, L. Ramos, I. Couceiro, L. Ramírez, Machine learning-based WENO5 scheme, Comput. Math. Appl. 168 (2024) 84–99.
[35] Z. Zhang, Y. Dong, Y. Zou, H. Zhang, X. Deng, A data-driven scale-invariant weighted compact nonlinear scheme for hyperbolic conservation laws, Commun.
Comput. Phys. 35 (4) (2024) 1120–1154.
[36] J. Zhu, C.-W. Shu, A new type of multi-resolution WENO schemes with increasingly higher order of accuracy, J. Comput. Phys. 375 (2018) 659–683.
[37] J. Qin, Y. Chen, Y. Lin, X. Deng, On construction of shock-capturing boundary closures for high-order finite difference method, Comput. Fluids 255 (2023).
[38] X. Deng, A new open-source library based on novel high-resolution structure-preserving convection schemes, J. Comput. Sci. 74 (2023) 102150.
[39] X. Deng, Z.-h. Jiang, C. Yan, Efficient ROUND schemes on non-uniform grids applied to discontinuous Galerkin schemes with Godunov-type finite volume
sub-cell limiting, J. Comput. Phys. 522 (2025) 113575.
[40] L. Cheng, X. Deng, B. Xie, An accurate and practical numerical solver for simulations of shock, vortices and turbulence interaction problems, Acta Astronaut.
210 (2023) 1–13.
[41] X. Deng, J.C. Massey, N. Swaminathan, Large-eddy simulation of bluff-body stabilized premixed flames with low-dissipative, structure-preserving convection
schemes, AIP Adv. 13 (5) (2023) 055014.
[42] E.F. Toro, Riemann Solvers and Numerical Methods for Fluid Dynamics: A Practical Introduction, Springer Science & Business Media, 2013.
[43] F. Xiao, Y. Honma, T. Kono, A simple algebraic interface capturing scheme using hyperbolic tangent function, Int. J. Numer. Methods Fluids 48 (2005)
1023–1040. https://doi.org/10.1002/fld.975
[44] R. Borges, M. Carmona, B. Costa, W.S. Don, An improved weighted essentially non-oscillatory scheme for hyperbolic conservation laws, J. Comput. Phys. 227
(6) (2008) 3191–3211.
[45] W.-S. Don, R. Borges, Accuracy of the weighted essentially non-oscillatory conservative finite difference schemes, J. Comput. Phys. 250 (2013) 347–372.
[46] S. Liu, Y. Shen, B. Chen, F. Zeng, Novel local smoothness indicators for improving the third-order WENO scheme, Int. J. Numer. Methods Fluids 87 (2) (2018)
51–69.
22

M. Huang, X. Deng, O.K. Matar et al. Journal of Computational Physics 555 (2026) 114764
[47] Y. Wang, Y. Du, K. Zhao, L. Yuan, A low-dissipation third-order weighted essentially nonoscillatory scheme with a new reference smoothness indicator, Int. J.
Numer. Methods Fluids 92 (9) (2020) 1212–1234.
[48] Y. Ha, C.H. Kim, H. Yang, J. Yoon, Construction of an improved third-order WENO scheme with a new smoothness indicator, J. Sci. Comput. 82 (3) (2020)
1–23.
[49] W. Xiaoshuai, Z. Yuxin, A high-resolution hybrid scheme for hyperbolic conservation laws, Int. J. Numer. Methods Fluids 78 (3) (2015) 162–187.
[50] W. Xu, W. Wu, An improved third-order WENO-Z scheme, J. Sci. Comput. 75 (3) (2018) 1808–1841.
[51] Y. Feng, F.S. Schranner, J. Winter, N.A. Adams, A multi-objective Bayesian optimization environment for systematic design of numerical schemes for compressible
flow, J. Comput. Phys. 468 (2022) 111477.
[52] Y. Feng, F.S. Schranner, J. Winter, N.A. Adams, A deep reinforcement learning framework for dynamic optimization of numerical schemes for compressible
flow simulations, J. Comput. Phys. 493 (2023) 112436.
[53] D. Kinga, J.B. Adam, et al., A method for stochastic optimization, in: International Conference on Learning Representations (ICLR), 5, California;, 2015.
[54] X. Deng, Z.-h. Jiang, O.K. Matar, C. Yan, On the convection boundedness of numerical schemes across discontinuities, Comput. Fluids (2025) 106645.
[55] S. Pirozzoli, On the spectral properties of shock-capturing schemes, J. Comput. Phys. 219 (2) (2006) 489–497.
[56] P.D. Lax, Weak solutions of nonlinear hyperbolic equations and their numerical computation, Commun. Pure Appl. Math. 7 (1) (1954) 159–193.
[57] C.-W. Shu, S. Osher, Efficient implementation of essentially non-oscillatory shock-capturing schemes, II, in: Upwind and High-Resolution Schemes, Springer,
1989, pp. 328–374.
[58] P. Woodward, P. Colella, The numerical simulation of two-dimensional fluid flow with strong shocks, J. Comput. Phys. 54 (1) (1984) 115–173.
[59] P.D. Lax, X.-D. Liu, Solution of two-dimensional riemann problems of gas dynamics by positive schemes, SIAM J. Sci. Comput. 19 (2) (1998) 319–340.
[60] J.J. Quirk, A Contribution to the Great Riemann Solver Debate, Springer, 1997.
[61] Z. Xu, C.-W. Shu, Anti-diffusive flux corrections for high order finite difference WENO schemes, J. Comput. Phys. 205 (2) (2005) 458–485.
[62] H. Wakimura, S. Takagi, F. Xiao, Symmetry-preserving enforcement of low-dissipation method based on boundary variation diminishing principle, Comput.
Fluids 233 (2022).
[63] S. Takagi, L. Fu, H. Wakimura, F. Xiao, A novel high-order low-dissipation TENO-THINC scheme for hyperbolic conservation laws, J. Comput. Phys. 452 (2022)
110899.
[64] R.D. Richtmyer, Taylor instability in shock acceleration of compressible fluids, Commun. Pure Appl. Math. 13 (2) (1960) 297–319.
[65] E.E. Meshkov, Instability of the interface of two gases accelerated by a shock wave, Fluid Dyn. 4 (5) (1969) 101–104.
[66] H. Zhang, Z. Chen, X. Jiang, Z. Huang, The starting flow structures and evolution of a supersonic planar jet, Comput. Fluids 114 (2015) 98–109.
[67] D.J. Hill, D.I. Pullin, Hybrid tuned center-difference-WENO method for large eddy simulations in the presence of strong shocks, J. Comput. Phys. 194 (2) (2004)
435–450.
[68] A. Rikanati, O. Sadot, G. Ben-Dor, D. Shvarts, T. Kuribayashi, K. Takayama, Shock-wave Mach-reflection slip-stream instability: a secondary small-scale turbulent
mixing phenomenon, Phys. Rev. Lett. 96 (2006) 174503.
[69] A. Bagabir, D. Drikakis, Mach number effects on shock-bubble interaction, Shock Waves 11 (3) (2001) 209–218.
[70] S. Osher, J.A. Sethian, Fronts propagating with curvature-dependent speed: algorithms based on Hamilton-Jacobi formulations, J. Comput. Phys. 79 (1) (1988)
12–49.
[71] S. Osher, R. Fedkiw, Level Set Methods and Dynamic Implicit surfaces, Vol. 153, New York, NY: Springer Nature, New York, NY, 1st ed., New York, NY, 2006.
[72] L. Han, X.Y. Hu, N.A. Adams, Adaptive multi-resolution method for compressible multi-phase flows with sharp interface model and pyramid data structure, J.
Comput. Phys. 262 (2014) 131–152.
[73] X.Y. Hu, B.C. Khoo, N.A. Adams, F.L. Huang, A conservative interface method for compressible flows, J. Comput. Phys. 219 (2) (2006) 553–578.
23