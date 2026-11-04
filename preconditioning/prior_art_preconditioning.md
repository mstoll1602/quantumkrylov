# Prior-art check: preconditioning QLSAs for KKT systems of PDE-constrained optimization

Checked on 2026-11-04 against the draft `paper/quantum_preconditioning.tex`
(Stoll, Pearson; working draft, November 2026), contributions C1–C5.

Verification markers used below:
- **[read]**: I read the relevant sections of the full text (arXiv HTML or PDF).
- **[partial]**: I read part of the full text. Which parts is stated.
- **[abs]**: abstract, metadata or search snippets only. Treat any detail as unconfirmed.

---

## 1. Verdict

**Go, but with substantial repositioning.** Nothing I found treats saddle-point/KKT
systems of PDE-constrained optimization with a QLSA and a parameter-robust
preconditioner. The *saddle-point* parts of C2 and the *quantum time-periodic* part
of C4 still look new. But three building blocks that the draft presents as
contributions are already in the literature:

1. **The commuting lemma (C2) is not new for Q1.**
   - Deiml–Peterseim (Math. Comp. 2026, arXiv:2403.19512), Thm. 6.2 and its
     commutative diagram: the gradient is applied on the coarse level and then
     transferred by an inclusion of discontinuous Q1 spaces with *orthonormal
     columns* (an isometry). Normalization grows linearly in L, i.e.
     α_G ~ √L in our notation. **[read]**
   - Jin–Liu–Ma–Yu (arXiv:2505.06866): the same "commuting identity" for finite
     differences. **[partial]**
   - **Yang–Yu–Zhang (arXiv:2610.06987, 4 Oct 2026)**: Lemma 3.3 gives
     E_L P_ℓ = T⁽⁰⁾E_ℓ and C_L P_ℓ = T C_ℓ with *both* transfers isometries, for
     *values and gradients*, for Q1 reaction–diffusion K + cM, i.e. our A_β.
     Block-encoding normalization 12dL. **[read]** This is our Lemma 3.x and
     Proposition 3.x for Q1, in an L²-orthonormal (Legendre) basis instead of
     Gauss-point values.

   For Q1 the two bases differ by an element-wise orthogonal 2×2 map, so the
   statements are equivalent. What remains ours: Q_k/P_k generality, the
   quadrature-point formulation, and the joint use with the KKT factorization.
2. **C3 (β-robust multilevel weights) is partly pre-empted.**
   - Yang–Yu–Zhang Thm. 2.1 proves, for weights ω_ℓ = (h_ℓ⁻² + c)^(-1/2) with
     *exact level mass scaling* M_ℓ^(-1/2), the bound 1/L ≤ λ ≤ 12dL for all
     c ≥ 0.
   - So κ ≤ 12dL², uniformly in the reaction coefficient c. With c = β^(-1/2) this
     is a *proven* β-uniform (up to L²) bound for A_β = M + √β K. **[read]**
   - Our s = 1 weights with Jacobi scaling are different. They target O(1)
     instead of O(L²), but have a 3^d floor and no proof. They must be compared
     with YYZ's weights, and the TODO "prove Assumption 3.x" now has a
     published partial answer.
3. **The √3 bound of C4 is classical for λ = iω.**
   - Krendl–Simoncini–Zulehner, Numer. Math. 124 (2013) 183–213, Thm. 5: for
     exactly our block [[M, √ν(K − iωM)], [√ν(K + iωM), −M]] and
     P = blkdiag(P, P) with P = M + √ν(K + ωM), the spectrum of P⁻¹𝓜 lies in
     [−1, −1/√3] ∪ [1/√3, 1]. **[read]**
   - Our Theorem 4.1 is this result extended to complex λ with Re λ ≥ 0, which
     covers the backward-Euler symbols. That extension is a modest
     generalization and must be attributed as such.
   - The same preconditioner D = M + √λ(K + kωσM) is already in Kolmbauer–Langer
     (SISC 2012), with a cruder bound κ ≤ 3 (as reported by Kollmann–Kolmbauer).
     **[read, via KK12 report]**
   - Note: KK12, which the draft cites for this, is **not** the λ-robust
     reference. Its preconditioner diag(K_ν + M_kωσ, …) has no mass term and
     degrades as the cost parameter → 0, as its own conclusions state. **[read]**

**The "Gap" paragraph is factually outdated.** "Quantum approaches to
PDE-constrained optimization have so far been variational heuristics" is no
longer true:
- Sato et al. (arXiv:2511.14420): fault-tolerant, explicit block encodings for
  PDE-constrained optimization.
- Huang et al. (arXiv:2610.26843, Oct 2026): a QLSA applied to the state–costate
  optimality system of finite-horizon LQR, with a horizon-independent
  condition-number bound.

Neither treats PDE KKT systems with mesh- and β-robust preconditioning, so the gap
survives in a narrower form.

**Contributions at risk:**
- **High:** C3 (needs comparison with YYZ and a proof).
- **High:** the "commuting lemma" and the √L subnormalization as a stand-alone
  claim in C2 (must be credited to DP24 / JLMY / YYZ).
- **Medium–high:** the √3 theorem in C4 (credit KSZ13; claim only the Re λ ≥ 0
  extension).
- **Medium:** C1 (every design principle appears somewhere already; present it
  as a synthesis).

**Still new as far as I can tell:**
- the two-sided preconditioned KKT system R̃ᵀ𝒜R̃ for QLSAs;
- the J-factorization 𝒢*𝒥𝒢 with a signed-permutation reflection, so that an
  *indefinite* saddle-point operator gets a well-scaled block encoding;
- the effective-condition-number theorem for the KKT system;
- the QFT-in-time treatment of time-periodic optimal control on a quantum
  computer, with the scaled-reflection coupling J_k = s_k U_k.

---

## 2. Works by topic

### 2.1 Quantum solvers for saddle-point / KKT / indefinite block systems

**Systematic search.** I queried the arXiv API (all categories) for
"saddle point"/"saddle-point" together with quantum/QLSA/preconditioner/block
encoding, for "KKT" in quant-ph, and for "block encoding" with "indefinite". I
also ran web searches for Stokes/Darcy/mixed FEM on quantum computers.

**Result.** I found **no** paper that block-encodes and preconditions a PDE
saddle-point (velocity–pressure, KKT, mixed FEM) system for a QLSA. Quantum
Navier–Stokes/Stokes works avoid the monolithic saddle-point system:
- Jin et al. (arXiv:2607.00281) use artificial compressibility with
  Schrödingerisation, i.e. no saddle point **[abs]**;
- arXiv:2406.00280 and arXiv:2603.18222 use classical pressure projection with
  a QLSA for the pressure Poisson equation **[abs]**;
- arXiv:2509.08807 (QNSS, spectral): the abstract does not say which linear
  system is solved **[abs]**.

**Overlap with C2:** none. Positioning: "to our knowledge, the first
preconditioned QLSA treatment of a PDE saddle-point system".

- **Sünderhauf, Generalized QSVT**, arXiv:2312.00723 (2023). **[abs]**
  - QSVT for arbitrary matrices with complex polynomials, and an eigenvalue
    variant for Hermitian matrices that allows polynomials of *indefinite
    parity*.
  - Overlap: low. It concerns the draft's statement in Section 2 that mixed
    parity "doubles α". That is not true in general for Hermitian inputs (see
    also QETU-type constructions). Soften the sentence.
- **Orsucci–Dunjko**, arXiv:2101.11868 (Quantum 2021), positive-definite QLSs
  with quadratically improved κ-dependence under extra structure. **[abs]**
  - Overlap: low.
  - Meaning for us: SPD structure can be exploited, but saddle-point systems are
    indefinite. This motivates why the J-factorization keeps the *indefinite*
    form instead of squaring it (normal equations would square κ).

### 2.2 Quantum algorithms for PDE-constrained optimization and optimal control

- **Sato, Kato, Yano, Ito, Yamamoto**, "Explicit block-encoding for partial
  differential equation-constrained optimization", arXiv:2511.14420 (v1 Nov
  2025, v2 May 2026). **[abs]**
  - A fully coherent pipeline: a quantum PDE solver (LCHS, per search snippets)
    prepares the state, and a quantum optimizer queries a block-encoded
    objective, so no tomography is needed.
  - Examples: Black–Scholes calibration, wave-equation material design.
  - Overlap with C2/C5: low–medium. Same application class (design/control/
    inference), but it is an outer-loop optimizer over a parameterized PDE, not
    the KKT system, and has no preconditioning.
  - *Meaning:* it contradicts the "only variational heuristics" sentence. Cite
    it as a complementary, black-box-optimizer approach.
- **Yano, Sato**, "Quantum framework for parameterizing PDEs via diagonal
  block-encoding", arXiv:2603.01358 (2026). **[abs, title only]** Same group,
  related. Overlap: low.
- **Huang, Álvarez Yuste, Guseynov, Liu, González-Conde**, "Quantum algorithms for
  optimal control: large-scale finite-horizon linear-quadratic regulation",
  arXiv:2610.26843 (v1 19 Oct 2026, v2 1 Nov 2026). **[partial: Sections 1–5.5]**
  - The state–costate system (Hamiltonian matrix) is written as a *multiple-
    shooting* linear system and solved by a QLSA.
  - κ(L_ex) is horizon-independent under stabilizability and detectability, but
    κ = O(h⁻¹) under time refinement.
  - No preconditioning. ODE (mass–spring chain), not PDE; no mesh or β analysis.
  - Not presented as symmetric/KKT. No Fourier in time.
  - Overlap: medium for the *message* ("QLSA on an optimality system"), low
    technically. **Must cite** as the closest quantum optimal-control-via-QLSA
    work. Position: we handle PDE constraints, mesh- and β-robustness, and the
    symmetric indefinite KKT form.
- **Clayton, Leng, Yang, Qiao, Lin, Wu**, "Differentiable quantum computing for
  large-scale linear control", NeurIPS 2024, arXiv:2411.01391. **[abs]**
  Infinite-horizon LQR via a quantum Lyapunov solver and policy gradients.
  Overlap: low.
- **Krovi**, arXiv:2404.07303 (2024), quadratic classical Hamiltonians and LQR via
  quantum ODE solvers. **[abs]** Overlap: low.
- **Variational PDE-constrained optimization.** **[abs]** Overlap: low. These are
  the "variational heuristics" the draft refers to; keep the citation.
  - Surana, Gnanasekaran, "Variational quantum framework for PDE constrained
    optimization", arXiv:2405.16651: bi-level, VQLS plus a black-box optimizer.
  - Gnanasekaran, Surana, Zhu, arXiv:2410.13688: the nonlinear version via
    Carleman linearization.
- **Quantum inverse problems.** **[abs]** Overlap: none.
  - Shao–Xiang, "Quantum regularized least squares solver with parameter
    estimate" (arXiv:1812.09934): Tikhonov via HHL on an extended matrix.
  - Chakraborty–Morolia–Peduri, "Quantum regularized least squares", Quantum
    2023.
- **Stein et al.**, "Exponential quantum speedup for simulation-based optimization
  applications", arXiv:2305.08482. **[abs]** QSVT combined with QAOA for
  "LinQuSO". No PDE KKT. Overlap: none.

### 2.3 Preconditioning for QLSAs (general)

- **Clader, Jacobs, Sprouse**, PRL 110, 250504 (2013). SPAI preconditioning
  inside HHL. Overlap with C1: low (historical). **[metadata verified]**
- **Shao, Xiang**, "Quantum circulant preconditioner for linear system of
  equations", PRA 98, 062321 (2018), arXiv:1807.04563. **[abs]**
  - A circulant preconditioner diagonalized by the QFT.
  - Overlap with C4: low–medium. It is the earliest QFT-diagonalized
    preconditioner in a QLSA. Cite it in Section 4.
- **Tong, An, Wiebe, Lin**, PRA 104, 032422 (2021). Fast inversion of a dominant
  part, giving preconditioned QLSAs. Overlap: low (cited). **[metadata verified]**
- **Bagherimehrab, Nakaji, Wiebe, Brennen, Sanders, Aspuru-Guzik**, "Fast quantum
  algorithm for differential equations", arXiv:2306.11802 (v3 Sep 2025); Phys.
  Rev. Research, accepted 19 Oct 2026, DOI 10.1103/2mlq-dbfg. **[abs + APS
  accepted page]**
  - A wavelet basis with a diagonal preconditioner gives N-independent κ.
  - Overlap: low–medium for C1. A multilevel basis is the same idea as BPX in a
    different guise.
- **Lapworth, Sünderhauf**, "Preconditioned block encodings for quantum linear
  systems", Quantum Sci. Technol. 10(4) (2025) 045064, arXiv:2502.20908 (v3 Oct
  2025). **[abs]**
  - Encoding A and P separately and multiplying them loses through
    subnormalization; encoding the classical product works better.
  - Introduces "preamplified quantum multiplication".
  - Overlap with C1/D2: **medium-high**. D2 is their main empirical finding.
    Cite it as the source of D2, and mention preamplified multiplication as an
    alternative to our factorization approach.
- **Deiml, Peterseim**, "Quantum realization of the finite element method", Math.
  Comp. 95 (2026) no. 362, 2637–2668 (online 14 Aug 2025), arXiv:2403.19512 (v4
  Jul 2025). **[read: Sections 3–6, conclusion, hardware]**
  - BPX via the generating system F, P = FFᵀ (Sec. 5.2).
  - Factorization S = C_Lᵀ(D_A ⊗ I)C_L (6.1), with C_L the gradient into an
    L²-orthonormal basis of discontinuous Q1.
  - FᵀSF = C_Fᵀ(D_A ⊗ I)C_F (6.2).
  - Thm. 6.2 (commutative diagram): C_L F|_{V_ℓ} = T̃_{ℓ,L} 2^{−ℓ(2−d)/2} C_ℓ,
    with T_{ℓ,L} the inclusion Q_ℓ → Q_L. It "has orthonormal columns",
    normalization and subnormalization 1.
  - Thm. 6.3 (Cartesian): subnormalization < (β/α) γ̃(U_A) d(L + π²/4),
    normalization 4γ(U_A)dL.
  - Def. 4.1: κ_eff ≥ γ̃(U_S) κ(S), i.e. our κ_α.
  - Kernel of FᵀSF is ignored by the solver (our D4).
  - Claimed advantage from d = 2. Remark 6.7: with h ≈ tol^{1/s}, advantage when
    d > s.
  - No mass/reaction term, no saddle point, no optimal control, no time
    dependence. Diffusion coefficient A may vary.
  - Overlap: **high** with the multilevel block encoding in C2 (Lemma 3.x,
    Prop. 3.x, √L) and with D1/D2/D4 in C1. **Must cite prominently.** Present
    our Lemma as "the analogue of [DP, Thm 6.2] including the mass term, in
    quadrature form, for Q_k/P_k".
  - Their bibliography also lists "Deiml, Peterseim, *Quantum multilevel
    preconditioning*, Oberwolfach Rep. (2025)", a short report (not checked).
- **Jin, Liu, Ma, Yu**, "Quantum preconditioning method for finite difference
  discretizations of the Poisson equation via Schrödingerization",
  arXiv:2505.06866 (v1 May 2025, v2 Jun 2026). **[partial: Sections 1–4]**
  - BPX with a "commuting identity" (Lemma 4.3): D_L P_ℓ = (1/r)E_{ℓ→L} D_ℓ,
    with E = I ⊗ 1_r, so T is a scaled isometry.
  - Exact block encoding of the symmetrically preconditioned FD Poisson matrix,
    normalization O(d²(L+1)). Solved via Schrödingerization, ε⁻¹ query
    complexity.
  - No mass term, no saddle point.
  - Overlap: **medium-high** with the C2 commuting lemma (FD version). Cite.
- **Yang, Yu, Zhang**, "A structure-aware construction of quantum BPX
  preconditioners for Q1 finite element discretizations of reaction–diffusion
  problems", arXiv:2610.06987 (v1, 4 Oct 2026). **[read: Sections 1–4,
  conclusions, references]**
  - −Δu + cu, constant c ≥ 0, Q1 on (0,1)^d.
  - Exact Gram factorizations M_ℓ = E_ℓᵀE_ℓ, A_ℓ = C_ℓᵀC_ℓ in a cellwise
    L²-orthonormal basis (Lemma 3.2).
  - **Lemma 3.3: E_L P_ℓ = T⁽⁰⁾E_ℓ and C_L P_ℓ = T C_ℓ, both isometries.**
  - Weights ω_ℓ = (h_ℓ⁻² + c)^(-1/2), S_ℓ = ω_ℓ P_ℓ M_ℓ^(-1/2).
  - **Thm. 2.1: 1/L ≤ λ(SᵀKS) ≤ 12dL for every c ≥ 0** (proof by L²-projection
    stable splitting).
  - Thm. 3.1: α = 12dL. They reject the product encoding, since
    α^prod ≥ κ(K)/L (our (3.x)).
  - Thm. 4.1/4.2: state preparation and functional estimation with explicit
    dL³ / dL² costs.
  - Cites DP and JLMY as precursors. No saddle point, no optimal control, no
    time dependence. Outlook: shifted systems for fractional operators.
  - Overlap: **high** with C2 (Lemma/Prop. for Q1, including the mass term) and
    **medium-high** with C3 (a proven c-uniform bound with log loss).
  - **Must cite and compare.** (i) Our A_β is their K + cM with c = β^(-1/2),
    scaled. (ii) Their weights equal our s = 0 weights with Jacobi scaling
    replaced by M_ℓ^(-1/2); their L² loss matches our observation that plain
    BPX degrades as β → 0. (iii) Our θ-weights try to remove the L-loss, so
    that has to be shown or proven. (iv) Their exact-mass scaling may remove our
    3^d floor, which is exactly the "better mass factor" in our outlook.
- **Fressart, Nowak, Spillane**, "Quantum domain decomposition for preconditioning
  the finite element method", arXiv:2605.26090 (May 2026). **[partial]**
  - Two-level additive Schwarz with BPX local solves, built on DP's factorized
    form.
  - Identity C_L R^(i)ᵀ = ℛ^(i)ᵀ C_L^(i) (Eq. 14). Normalization O(NdL).
  - Poisson only. Overlap: low–medium (same factorized-encoding philosophy).
    Cite as related.
- **Suresh**, "A preconditioned quantum linear solver for periodic two-phase heat
  conduction with mesh-independent condition number", arXiv:2610.41690 (30 Oct
  2026). **[abs]**
  - Folds the closed-form inverse square root of the homogeneous operator into
    the FE factors. The preconditioned operator becomes
    (norm-preserving gradient-type factor) × diag × transpose, with
    subnormalization 1, implemented by a QFT in *space*.
  - Overlap: medium with C1/D2 and with the "factor so every factor is well
    scaled" idea. SPD only. Cite.
- **Nie, Lai, An**, "Pauli-structured preconditioning for QLS solvers",
  arXiv:2606.01733 (Jun 2026). **[abs]**
  - Normalization overhead of composed encodings cancels preconditioning;
    regrouping Pauli expansions fixes it in some cases.
  - Overlap: medium with D2. Cite.
- **Pechan, Golden, O'Malley**, "Block encoding the 3D heterogeneous Poisson
  equation with application to fracture flow", arXiv:2508.07125 (v2 Mar 2026).
  **[abs]** Separate encoding of A and P does not reduce the effective κ, which
  they call "a key barrier". Overlap: medium with D2. Also Golden, O'Malley,
  Viswanathan, arXiv:2205.06379 (preconditioners for hydrological QLSs) **[abs]**.
- **Low, Su**, arXiv:2410.18178, optimal queries to state preparation with
  "block preconditioning". **[abs]** Overlap: low (different meaning of
  "preconditioning").
- **Dalzell, Li, Su**, "Faster quantum linear system solver beyond the condition
  number", arXiv:2607.07691 (Jul/Oct 2026). **[abs]**
  - Solvers whose cost depends on an instance-dependent κ_eff (norms like
    ‖(A†A)^(-t/2)|x⟩‖) rather than κ.
  - Overlap: concerns C1/D1. "Only κ_α counts; clustering does not help" holds
    for worst-case QSVT inversion, but this line of work shows that
    instance-dependent quantities can matter. Qualify D1 ("for QSVT-based
    inversion with worst-case guarantees").
- **Jiang, Park, Xu**, "A polynomial dimension-dependence analysis of BPX
  preconditioners", arXiv:2512.06166 (v3 Oct 2026). **[abs]** Explicit polynomial
  d-dependence of BPX for P1, motivated by quantum FEM. Useful for our C_d
  discussion.
- **X. Li**, "A Quantum Path to Partial Differential Equations", arXiv:2607.09639,
  lecture notes. **[abs]** No saddle point, control or preconditioning.

### 2.4 Quantum interior point methods (Newton/KKT systems)

- **Kerenidis, Prakash**, "A quantum interior point method for LPs and SDPs", ACM
  TQC 1 (2020) 1–32. **[verified via citation in Augustino et al.]**
- **Augustino, Nannicini, Terlaky, Zuluaga**, "Quantum interior point methods for
  semidefinite optimization", Quantum 7, 1110 (2023), arXiv:2112.06025. **[abs]**
- **Dalzell, Clader, Salton, Berta, Lin, Bader, Stamatopoulos, Schuetz, Brandão,
  Katzgraber, Zeng**, "End-to-end resource analysis for quantum interior-point
  methods and portfolio optimization", PRX Quantum 4, 040325 (2023),
  arXiv:2211.12489. **[partial]**
  - Solves the full Newton system (II-QIPM) or a null-space reduced system.
  - "A simple preconditioning method compatible with the QIPM": in numerical
    experiments it lowers κ_F by more than an order of magnitude, with no
    asymptotic bound.
- **Wu, Mohammadisiahroudi, Augustino, Yang, Terlaky**, "An inexact feasible
  quantum interior point method for linearly constrained quadratic
  optimization", Entropy 25 (2023) 330, arXiv:2301.05357. **[abs]**
  - The orthogonal subspace system (OSS) gives a square system on which inexact
    solves keep feasibility.
  - The QP KKT system is handled by a null-space reformulation, not by
    saddle-point preconditioning.
- **Mohammadisiahroudi, Wu, Augustino, Carr, Terlaky**, "Improvements to quantum
  interior point method for linear optimization", arXiv:2310.07574. **[abs]**
  Iterative refinement with a preconditioned normal equation system gives a
  uniform κ bound.
- **Wu, Yang, Terlaky**, "A preconditioned inexact infeasible quantum interior
  point method for linear optimization", arXiv:2412.11307 (Dec 2024). **[abs]**
  Optimal-partition-based (Chai–Toh) preconditioning: κ goes from O(1/μ²) to
  O(1/μ).
- **Mohammadisiahroudi, Wu, Sampourmahani, Harkness, Terlaky**, "Quantum Interior
  Point Methods: A Review of Developments and An Optimally Scaling Framework",
  arXiv:2512.06224 (v1, 6 Dec 2025; symposium paper). **[abs]**
- **Overlap with us:** low. These works precondition *finite-dimensional
  LP/QP/SDP* Newton systems (normal equations/OSS, μ-dependence). There are no
  PDE structure, mesh-robustness or saddle-point block preconditioners of
  Schöberl–Zulehner type.
- *Meaning:* our outlook item 3 (control constraints, semismooth Newton/IPM)
  is a genuine bridge. Cite Wu–Yang–Terlaky and the review there.

### 2.5 Parabolic / time-periodic problems, Fourier in time, circulant structure

- **No quantum work found** that applies a QFT *in time* to a time-periodic
  (circulant) all-at-once system, or that uses ParaDiag/α-circulant structure in
  a QLSA.
  - Searches: arXiv API "quantum + circulant + linear system", "quantum +
    time-periodic"; web searches for quantum ParaDiag, Floquet/harmonic balance
    with QLSA, quantum all-at-once heat.
  - Closest:
    - Shao–Xiang 2018 (circulant preconditioner via QFT, Sec. 2.3).
    - Wan et al., "Block-encoding based quantum algorithm for linear systems with
      displacement structures", PRA 104, 062414 (2021), arXiv:1912.12212
      (circulant/Toeplitz block encodings). **[abs]**
    - Suresh 2026 (QFT in *space*, periodic cell).
    - Huang et al., "Classical combinations of quantum states for solving banded
      circulant linear systems", arXiv:2309.11451. **[abs]**
    - History-state/Dyson-series all-at-once ODE solvers (Berry–Costa, Quantum
      2024; Fang–Lin–Tong time marching) do not exploit periodicity. **[abs]**
    - Linden–Montanaro–Shao (heat equation; cited).
  - **Overlap with C4 (quantum side): low.** The QFT-in-time decoupling into
    frequency-wise KKT blocks, and the scaled reflection J_k, look new.
- **Classical multiharmonic / time-periodic preconditioning**, which answers the
  specific question:
  - **Krendl, Simoncini, Zulehner**, "Stability estimates and structural spectral
    properties of saddle point problems", Numer. Math. 124(1) (2013) 183–213,
    arXiv:1202.3330. **[read, Section 5.1]**
    - The reduced system [[M, √ν(K − iωM)], [√ν(K + iωM), −M]] (their (39)).
    - P = blkdiag(P, P) with **P = M + √ν(K + ωM)**, i.e. exactly our
      A_k = (1 + √β|λ|)M + √βK for λ = iω.
    - **Theorem 5: inf-sup constant γ ≥ 1/√3 and ‖B‖ ≤ 1, so
      spec(P⁻¹𝓜) ⊂ [−1, −1/√3] ∪ [1/√3, 1], symmetric about 0.**
    - Proof via the auxiliary operator H = [[I, (1−i)I], [(1+i)I, −I]].
    - Also the time-periodic Stokes control (Thm. 6, interval
      [−1, ½(1−√3)] ∪ [µ₃, ½(1+√5)] for a different system).
    - **This is the √3 constant of C4.** Our only extension is complex λ with
      Re λ ≥ 0 (backward Euler; KSZ use the ω real argument
      ‖y‖²_K + ω‖y‖²_M ≤ … with iω).
  - **Kolmbauer, Langer**, "A robust preconditioned MinRes solver for distributed
    time-periodic eddy current optimal control problems", SIAM J. Sci. Comput.
    34(6) (2012) B785–B809.
    - **[not read; content as reported in KK12's NuMa report]**: block-diagonal
      P = diag(D, D, λ⁻¹D, λ⁻¹D) with D = √λK_ν + √λM_kωσ + M, robust in
      h, k, ω, σ, ν, λ, with "a rigorous condition number bound κ ≤ 3".
  - **Kollmann, Kolmbauer**, "A preconditioned MinRes solver for time-periodic
    parabolic optimal control problems", Numer. Linear Algebra Appl. 20(5)
    (2013) 761–784 (online 2012). **[read, NuMa-Report 2011-06 version]**
    - Uses C = diag(K_ν + M_kωσ, …), *without* the mass term.
    - Robust in h, kω, ν, σ but **not** in the cost parameter λ (their
      conclusion: "for very small λ the convergence rate deteriorates").
    - Constants are not explicit (c̄/c via Zulehner's Thm. 2.6).
    - The extension to control constraints and partial observation/control
      domains is its main point.
    - It is the wrong primary reference for a λ-robust √3 bound.
  - **Kollmann, Kolmbauer, Langer, Wolfmayr, Zulehner**, "A robust finite element
    solver for a multiharmonic parabolic optimal control problem", Comput. Math.
    Appl. 65(3) (2013) 469–486. **[metadata only; full text not accessible:
    publisher 403, JKU preprint server has an expired certificate]**
    - By its abstract, the main contribution is a robust block preconditioner for
      the multiharmonic optimality system.
    - **I could not confirm which constant it proves.** Most likely the same D
      with a bound of the KL12/KSZ13 type.
    - Check the paper before attributing a specific constant.
  - **Langer, Wolfmayr**, "Multiharmonic finite element analysis of a
    time-periodic parabolic optimal control problem", J. Numer. Math. 21 (2013),
    doi 10.1515/jnum-2013-0011. **[abs]** Robust AMLI-type preconditioners. Also
    Wolfmayr's PhD thesis (JKU 2014), not checked.
  - Zulehner, "Nonstandard norms and robust estimates for saddle point
    problems", SIAM J. Matrix Anal. Appl. 32 (2011) 536–560: the ω = 0 case and
    the framework.
- **Classical ParaDiag for parabolic optimal control.** **[abs or partial]**
  - Wu–Zhou (ESAIM COCV 2020); Bouillon–Samaey–Meerbergen arXiv:2302.06406
    (α-circulant, generalized ParaDiag for optimal control); Lin–Wu–Zhang
    arXiv:2109.12524; Hon–Fung–Lin arXiv:2410.22686.
  - Hon–Fung–Lin, **[partial]**: the ε-circulant RBD preconditioner uses the
    complex λ_k per frequency, not |λ_k|, and has GMRES bounds without a 1/√3.
  - None of the works I inspected uses blkdiag(A_k, A_k) with
    A_k = (1 + √β|λ_k|)M + √βK for complex λ_k with Re λ_k ≥ 0, i.e. our
    "|D_t|" preconditioner. The Re λ ≥ 0 extension is not exhaustively checked:
    PRESB-type analyses by Axelsson, Neytcheva, Liang and others for complex /
    time-harmonic control may contain it.
  - Also cite McDonald–Pestana–Wathen 2018 (already cited) and the
    absolute-value circulant preconditioner idea for the |D_t| remark.

### 2.6 J-factorizations / indefinite factorizations for block encodings; Hermitian dilations

- I found **no** quantum paper that writes an indefinite saddle-point matrix as
  𝒢*𝒥𝒢 with a signature/reflection 𝒥 to control subnormalization. Searches
  combined "block encoding" with indefinite, signature matrix, J-factorization,
  hyperbolic, and G^T J G.
- Closest structural relatives:
  - DP24 (6.1): S = C_Lᵀ(D_A ⊗ I)C_L, with SPD middle factor D_A.
  - Fressart et al. (same form).
  - Suresh 2026: Vᵀ diag V, SPD.
  - Gram/LCU factorizations in quantum chemistry ("double factorization"),
    which are definite.
- Hermitian dilation [[0, A], [A*, 0]] is the standard HHL device for
  non-Hermitian A. It is unnecessary with QSVT, and for KKT matrices it would
  hide rather than exploit the saddle-point structure.
- **Overlap with C2 (J-factorization): low.** This looks new for QLSAs. The
  classical idea (J-orthogonal/hyperbolic factorizations; LDLᵀ with D = ±1) is
  old and is cited as such in the draft.

---

## 3. Closest competitors (ranked) and how to position

1. **Yang–Yu–Zhang, arXiv:2610.06987 (Oct 2026).**
   - Same matrix A_β (reaction–diffusion), the same isometric commuting transfers
     for values and gradients, the same √L / L normalization, and a proven
     c-uniform bound.
   - *Must cite.* Position: we use their kind of encoding as a *building block*
     inside a two-sided KKT preconditioner, plus the J-factorization for an
     indefinite system. Our quadrature form extends to Q_k/P_k.
   - Our s = 1 weights must be compared with their (h_ℓ⁻² + c)^(-1/2) M_ℓ^(-1/2)
     weights, in both theory and numerics. Adopting their exact level-mass
     scaling may remove the 3^d floor.
   - If the s = 1 weights do not beat L² in theory, drop C3 as a stand-alone
     contribution or restate it as "numerically O(1)".
2. **Deiml–Peterseim, Math. Comp. 2026.**
   - Origin of the factorized BPX block encoding (commuting diagram, isometric
     inclusion, normalization ∝ dL) and of κ_eff = subnormalization × κ.
   - *Must cite prominently* in the C1 and C2 text. Reconcile our ε⁻¹-with-d = 3
     advantage claim with their d = 2 claim. They take h ≈ tol (H²; Remark 6.7);
     we take h ~ ε^{1/2} (L²/QoI superconvergence). State the cost-model
     difference explicitly.
3. **Krendl–Simoncini–Zulehner, Numer. Math. 2013** (with Kolmbauer–Langer SISC
   2012 and KKLWZ CAMWA 2013).
   - The √3 bound for exactly our frequency-wise blocks when λ = iω.
   - Restate Theorem 4.1 as "a slight extension of [KSZ13, Thm 5] to
     Re λ_k ≥ 0", keep our short proof, and remove "we prove a √3 bound" from
     the abstract as a headline result.
   - The novelty in C4 is the quantum part: QFT decoupling in superposition,
     the scaled reflection J_k, and that N_t and τ do not enter κ_α.
4. **Jin–Liu–Ma–Yu, arXiv:2505.06866.** The commuting identity in FD form. Cite in
   the same sentence as DP24 and YYZ.
5. **Huang et al., arXiv:2610.26843** and **Sato et al., arXiv:2511.14420.**
   Non-variational quantum optimal control / PDE-constrained optimization. Cite
   them in a rewritten "Gap" paragraph: "existing approaches either use
   variational heuristics [SG24, GSZ24], an outer quantum optimizer over a
   block-encoded objective [Sato+25], or QLSAs for ODE-constrained LQR without
   preconditioning [Huang+26]. None addresses the mesh- and β-dependence of the
   KKT system of a PDE-constrained problem."
6. **Lapworth–Sünderhauf (QST 2025), Pechan–Golden–O'Malley, Nie–Lai–An,
   Suresh.** Prior statements of D2 ("separately encoded factors lose"). Cite in
   Section 2 as the sources of D2. Suresh is a further SPD instance of "fold the
   preconditioner into the factors".
7. **QIPM line** (Wu–Yang–Terlaky 2024; review 2512.06224; Dalzell et al. 2023):
   only for the outlook (control constraints).

---

## 4. What is still new, contribution by contribution

- **C1 (design principles).** Each principle is already stated somewhere:
  - D1 κ_eff: DP24 Def. 4.1, LS25, YYZ;
  - D2: LS25, Pechan+, Nie–Lai–An, DP24, YYZ;
  - D3 additive vs multiplicative: implicit in DP24's choice of BPX;
  - D4 kernels: DP24, YYZ;
  - D5 QoIs: MP16, DP24, YYZ Thm 4.2.

  New is only the *synthesis for indefinite saddle-point systems* (e.g. "every
  factor of a product must be well scaled" applied to a J-factorization). Present
  C1 as a collection with citations, not as a contribution. Qualify D1 in light
  of Dalzell–Li–Su 2026 and the parity remark in light of generalized QSVT.
- **C2 (Poisson control).**
  - **New:** (i) applying QLSA preconditioning to the KKT system at all; (ii) the
    two-sided R̃ᵀ𝒜R̃ with Theorem 3.x (eigenvalue intervals via Ostrowski,
    consistency, x = R̃𝓑⁺R̃ᵀf); (iii) the J-factorization 𝓑 = 𝒢ᵀ𝒥𝒢 with 𝒥 a
    two-qubit signed permutation, giving an indefinite operator a well-scaled
    encoding, and the demonstration that the product encoding fails
    (α ≥ (c₀/2)κ(A_β)); (iv) right-hand side and outputs through 𝒢 (D5 for the
    KKT setting).
  - **Not new:** the commuting lemma for Q1, the isometric transfers, and
    α_G = O(√L). New only in the Q_k/P_k and quadrature-point generality, which
    is a modest extension.
- **C3 (level weights).** At risk. YYZ already prove β-uniform (up to L²) bounds
  with natural weights. What can stay:
  - the observation that plain Jacobi-scaled BPX degrades as β → 0 while
    θ-weights keep κ(CA_β) = O(1) numerically;
  - the "truncation does not work" observation.

  To be a contribution it needs a proof of an L-independent, β-robust bound, or
  at least a numerical comparison showing that it beats YYZ's L-growth.
  Classical robust multilevel theory for I − ε²Δ should also be checked before
  writing (not covered here).
- **C4 (time-periodic heat control).**
  - **New:** the quantum treatment. QFT in time on the KKT system, all
    frequencies in superposition, the frequency-dependent scaled reflection 𝒥_k
    with √r_k value weighting (‖J_k‖ ≤ 1), and κ_α independent of N_t and τ.
    Also "the level weights are no longer optional" (β_k → small at high
    frequencies).
  - **Not new:** the √3 bound for λ = iω (KSZ13).
  - **Small extension:** Re λ ≥ 0 (backward Euler, possibly other A-stable
    symbols), and the time-domain form with |D_t|.
- **C5 (numerics and complexity).** Fine. The complexity discussion must
  reconcile with DP24 (d = 2 claim) and YYZ (dL³ state preparation, dL² + dL³/ε
  functionals). The "rough in time" argument for heat control appears new.

---

## 5. Bibliography corrections

| Key | Draft | Correct (verified via Crossref/arXiv/publisher) |
|---|---|---|
| KK12 | NLA (2012), doi 10.1002/nla.1842 | M. Kollmann, M. Kolmbauer, *Numer. Linear Algebra Appl.* **20**(5) (2013) 761–784 (online 27 Apr 2012), doi 10.1002/nla.1842. Rename key KK13. **Not** the λ-robust reference (see 2.5). |
| KKLWZ13 | CAMWA 65 (2013) 469–486 | Correct: *Comput. Math. Appl.* **65**(3) (2013) 469–486, doi 10.1016/j.camwa.2012.06.012. |
| QIPMreview | "M. Mohammadisiahroudi et al." | M. Mohammadisiahroudi, Z. Wu, P. Sampourmahani, A. Harkness, T. Terlaky, "Quantum Interior Point Methods: A Review of Developments and An Optimally Scaling Framework", arXiv:2512.06224 (2025). Title in draft is fine. |
| BNWA23 | 4 authors, arXiv only | **Six authors**: M. Bagherimehrab, K. Nakaji, N. Wiebe, **G. K. Brennen, B. C. Sanders**, A. Aspuru-Guzik, "Fast quantum algorithm for differential equations", *Phys. Rev. Research*, accepted 19 Oct 2026, doi 10.1103/2mlq-dbfg; arXiv:2306.11802. Volume and article number not yet assigned. Rename key (e.g. BNWBSA26). |
| LS25 | arXiv:2502.20908 | L. Lapworth, C. Sünderhauf, "Preconditioned block encodings for quantum linear systems", *Quantum Sci. Technol.* **10**(4) (2025) 045064, doi 10.1088/2058-9565/ae0f4b. |
| MP16 | PRA 93 (2016) 032324 | Correct. |
| TAWL21 | PRA 104 (2021) 032422 | Correct. The published title includes "fast Green's-function computation"; the arXiv title omits it. |
| CJS13 | PRL 110 (2013) 250504 | Correct. |
| DP24 | arXiv:2403.19512 | M. Deiml, D. Peterseim, "Quantum realization of the finite element method", *Math. Comp.* **95** (2026), no. 362, 2637–2668 (online 14 Aug 2025), doi 10.1090/mcom/4137. Rename key DP26 or keep DP24 with the journal data. |
| SZ07 | SIMAX 29 (2007) 752–773 | Correct (issue 3). |

**Add:**
- Krendl–Simoncini–Zulehner, Numer. Math. 124(1) (2013) 183–213, doi
  10.1007/s00211-012-0507-3.
- Kolmbauer–Langer, SIAM J. Sci. Comput. 34(6) (2012) B785–B809, doi
  10.1137/110842533.
- Yang–Yu–Zhang arXiv:2610.06987.
- Jin–Liu–Ma–Yu arXiv:2505.06866.
- Huang et al. arXiv:2610.26843.
- Sato et al. arXiv:2511.14420.
- Shao–Xiang PRA 98, 062321 (2018).
- Suresh arXiv:2610.41690.
- Fressart–Nowak–Spillane arXiv:2605.26090.
- Pechan–Golden–O'Malley arXiv:2508.07125.
- Nie–Lai–An arXiv:2606.01733.
- Wu–Yang–Terlaky arXiv:2412.11307.
- Dalzell–Li–Su arXiv:2607.07691.
- Zulehner SIMAX 32 (2011) 536–560.
- Optionally Montanaro–Pallister is already present. Add Jiang–Park–Xu
  arXiv:2512.06166 for the d-dependence of BPX.

---

## 6. Must-read before revising

1. Yang–Yu–Zhang, arXiv:2610.06987: Thm. 2.1 (proof), Lemma 3.3, Thm. 3.1, and
   the weights (2.9). Directly decides the fate of C3.
2. Deiml–Peterseim, arXiv:2403.19512 v4 / Math. Comp.: Sections 5–6 (Thm. 6.2,
   6.3, 6.5, Remarks 6.6–6.7).
3. Krendl–Simoncini–Zulehner 2013, Section 5.1 (Theorem 5 and Remark 1).
4. KKLWZ13 (CAMWA) and Kolmbauer–Langer 2012 (SISC): which constants they prove.
   I could not access KKLWZ13.
5. Jin–Liu–Ma–Yu, arXiv:2505.06866 v2 (Lemma 4.3 and Theorem 4.1).
6. Huang et al., arXiv:2610.26843, and Sato et al., arXiv:2511.14420 (for the
   rewritten "Gap").
7. Lapworth–Sünderhauf (QST 2025): preamplified multiplication, as a possible
   alternative to our factorization.
8. Classical robust multilevel methods for I − ε²Δ: not covered here; needed for
   C3.

---

## 7. Sources

- DP24: https://arxiv.org/abs/2403.19512 ; https://arxiv.org/html/2403.19512v4 ; doi https://doi.org/10.1090/mcom/4137
- Jin–Liu–Ma–Yu: https://arxiv.org/abs/2505.06866 ; https://arxiv.org/html/2505.06866v2
- Yang–Yu–Zhang: https://arxiv.org/abs/2610.06987 ; https://arxiv.org/html/2610.06987v1
- Fressart–Nowak–Spillane: https://arxiv.org/abs/2605.26090
- Suresh: https://arxiv.org/abs/2610.41690
- Pechan–Golden–O'Malley: https://arxiv.org/abs/2508.07125
- Nie–Lai–An: https://arxiv.org/abs/2606.01733
- Lapworth–Sünderhauf: https://arxiv.org/abs/2502.20908 ; doi https://doi.org/10.1088/2058-9565/ae0f4b
- Bagherimehrab et al.: https://arxiv.org/abs/2306.11802 ; https://journals.aps.org/prresearch/accepted/10.1103/2mlq-dbfg
- Tong–An–Wiebe–Lin: https://arxiv.org/abs/2008.13295 ; https://journals.aps.org/pra/abstract/10.1103/PhysRevA.104.032422
- Clader–Jacobs–Sprouse: https://arxiv.org/abs/1301.2340
- Montanaro–Pallister: https://arxiv.org/abs/1512.05903
- Shao–Xiang circulant: https://arxiv.org/abs/1807.04563
- Wan et al. displacement structures: https://arxiv.org/abs/1912.12212
- Dalzell–Li–Su: https://arxiv.org/abs/2607.07691
- Low–Su: https://arxiv.org/abs/2410.18178
- Sünderhauf generalized QSVT: https://arxiv.org/abs/2312.00723
- Orsucci–Dunjko: https://arxiv.org/abs/2101.11868
- Jiang–Park–Xu: https://arxiv.org/abs/2512.06166
- X. Li lecture notes: https://arxiv.org/abs/2607.09639
- Sato et al.: https://arxiv.org/abs/2511.14420
- Huang et al. LQR: https://arxiv.org/abs/2610.26843 ; https://arxiv.org/html/2610.26843
- Clayton et al.: https://arxiv.org/abs/2411.01391
- Krovi: https://arxiv.org/abs/2404.07303
- Surana–Gnanasekaran: https://arxiv.org/abs/2405.16651 ; Gnanasekaran et al.: https://arxiv.org/abs/2410.13688
- Stein et al. QuSO: https://arxiv.org/abs/2305.08482
- Shao–Xiang regularized LS: https://arxiv.org/abs/1812.09934
- QIPM review: https://arxiv.org/abs/2512.06224
- Wu–Yang–Terlaky: https://arxiv.org/abs/2412.11307
- Wu et al. IF-QIPM QP: https://arxiv.org/abs/2301.05357
- Mohammadisiahroudi et al. improvements: https://arxiv.org/abs/2310.07574
- Dalzell et al. QIPM: https://arxiv.org/abs/2211.12489
- Augustino et al. SDP QIPM: https://arxiv.org/abs/2112.06025
- Krendl–Simoncini–Zulehner: https://arxiv.org/abs/1202.3330 ; doi https://doi.org/10.1007/s00211-012-0507-3
- Kollmann–Kolmbauer NuMa-Report 2011-06: https://numa.jku.at/media/filer_public/0b/b2/0bb2a440-b8d9-4191-a0ac-12197f5be5db/2011-06.pdf ; doi https://doi.org/10.1002/nla.1842
- KKLWZ13: https://research.jku.at/en/publications/a-robust-finite-element-solver-for-a-multiharmonic-parabolic-opti/ ; doi https://doi.org/10.1016/j.camwa.2012.06.012
- Kolmbauer–Langer: doi https://doi.org/10.1137/110842533
- Langer–Wolfmayr: https://www.degruyterbrill.com/document/doi/10.1515/jnum-2013-0011/pdf
- Bouillon–Samaey–Meerbergen: https://arxiv.org/abs/2302.06406
- Hon–Fung–Lin: https://arxiv.org/abs/2410.22686 ; https://arxiv.org/html/2410.22686v2
- Lin–Wu–Zhang: https://arxiv.org/abs/2109.12524
- Quantum Stokes via Schrödingerisation: https://arxiv.org/abs/2607.00281
- QNSS: https://arxiv.org/abs/2509.08807
- Crossref metadata API (bibliographic checks): https://api.crossref.org/works/<doi>
