# Review of `proof_bpx_robust.tex` (robustness of the weighted BPX factor for M + sqrt(beta) K)

Reviewer: Claude (adversarial check), 2026-11-04.
Scripts and raw outputs: `/tmp/claude-1000/-home-martin-Work-research-codes-QuantumKrylov/ed397104-e7df-47b8-bcd7-dbdae8e25de3/scratchpad/proofcheck/`
(`common.py`, `check1_lions_nodal_lambda.py`, `check2_jackson_scs.py`, `check3_N_kappa.py`,
`check4_heat_trunc_limit.py`, `check5_CH_bounds.py`, `check6_floor.py`, `out*.txt`).
They import `fe` and `prolong` from `code/qp_bpx_sparse.py`, with bytecode writing disabled. No repository file was changed.

Notation as in the draft: nu = sqrt(beta), t_l = nu 4^l, mu_l = 1 + t_l, lambda_l = mu_l / theta~_l.

---

## 1. Verdict

**Correct, with minor fixes.** The two cited lemmas are classical and true for Q1 on uniform Cartesian grids with homogeneous Dirichlet conditions. Given them, every step of the argument is right, and so are the constants c0, c1 and c1/c0 <= 3d 3^d K0 K1. No step fails numerically.

The findings fall into four groups:

- **The paper, not the proof.** The paper's own definitions (eq. (theta), and theta_{l,k} in the heat section) use **unnormalized** weights. The proof and the heat corollary need the **normalized** weights theta~_l = theta_l / theta_L.
  - With the paper's weights, spec(C_k A_k) is contained in theta_{L,k} [c0, c1].
  - c1/c0 stays uniform, but c0 does not. It degrades whenever sqrt(beta_k) h_L^{-2} < 1.
  - Theorem heatB uses c0 itself (kappa_alpha <= sqrt3 alpha_G^2 / c0), and so does the Poisson cost bound. The paper should switch to normalized weights. This does not affect Prop. alphaG, because theta~_l <= 1.
- **Wording and constants.**
  - The remark that the factor 3d 3^d "is the mass-matrix floor seen in the numerics" is misleading.
  - C_I = 4d 3^d is valid but loose; the sharp value is 12d.
  - A few statements are imprecise.
- **Scope.** The proof works verbatim for every s > 0, including s > 1 and the truncated limit s = infinity. This contradicts the paper's sentence "Truncation does not work". Truncation is robust in theory, only with worse constants.
- **Novelty.** The result is essentially known.
  - Griebel–Oswald 1995 prove robust (truncated) nodal BPX/MDS for -eps Delta + I with Q1 on (0,1)^d.
  - Bramble–Pasciak–Vassilevski 2000 and Bacuta–Hayes–Jacavage 2022 give an additive multilevel preconditioner whose level weights gamma_j - gamma_{j+1} are, in the resolved regime, within level-independent factors of the s = 1 weights.
  - What is new is the elementary self-contained proof, the family s in (0,1], and explicit (pessimistic) constants.

---

## 2. Step-by-step findings

| # | Step | Severity | Finding |
|---|------|----------|---------|
| 1 | Lemma lions (statement and proof) | ok | Correct. The minimizer y_k = D_k^{-1} R_k^T C^{-1} x and the Cauchy–Schwarz step are fine. Checked numerically (rel. err <= 6e-15, d = 1, 2, all s). Attribution: see section 3.3. |
| 2 | Application R_l = P_l, D_l = (delta_l/theta~_l) I | ok | P_l from `prolong` is the natural embedding: P_l^T M_L P_l = M_l and P_l^T K_L P_l = K_l to 4e-16, and P_L = I (d = 1, 2, 3). The map y_l -> u_l is bijective, sum P_l y_l = x iff sum u_l = u. P_L = I gives spanning. |
| 3 | (eq:twosided) | ok | Order reversal under inversion is correct: c1^{-1} A <= C^{-1} <= c0^{-1} A. |
| 4 | Lemma nodal | ok (minor) | Mass eigenvalues lie in ((h/3)^d, h^d); the diagonal formula for delta_l is exact; 1 <= (1 + 3d t)/(1 + t) <= 3d. The strict inequalities need y != 0. |
| 5 | (eq:NC) and c0, c1, c1/c0 | ok | Per-term factors (2/3)^d and 3d 2^d are right. The sandwich is nearly attained: x^T C^{-1} x / N(u) lies in [0.68, 6.00] for d = 1 (bound [0.667, 6]) and in [0.45, 23.75] for d = 2, L = 5 (bound [0.444, 24]). c0 = 1/(3d 2^d K0), c1 = (3/2)^d K1 and c1/c0 = 3d 3^d K0 K1 are correct. |
| 6 | "the factor 3d 3^d ... is the mass-matrix floor seen in the numerics" | minor (misleading) | The 3d part is an artifact of bounding delta_l by mu_l (diagonal of K vs h^{-2}). The sandwich is sharp, but its extremes do not align with those of N vs a, so kappa is much smaller. The observed unresolved floor is not 3^d but **3^d/(1 - 4^{-s})**: kappa(C M) -> 4.000 (d = 1, s = 1), 5.99 (d = 1, s = 1/2), 11.93 (d = 2, s = 1, L = 6), 17.3 (d = 2, s = 1/2). This also corrects the paper's statement "approach the Jacobi constant 3^d" (Sec. 3 bullets and Sec. 6). In 3D the floor for s = 1 should approach 36, not 27, as L grows. |
| 7 | Regime split, lambda_l formulas, t_{l*} < 4 | ok | All formulas are right. t_{l*} = 4 t_{l*-1} < 4 when l* > 1. When l* = 1 there are no l < l* levels and t_1 may be large; this is handled correctly downstream. At t_L = 1 exactly, l* = L and lambda_L = 2. The two regimes join continuously at t_L = 1. |
| 8 | Lemma lambda (a) | ok | theta~_l <= 1 in both regimes, so lambda_l >= 1 > t_l when t_l < 1. Holds for all s >= 0. Verified on a grid (L <= 20, nu in [1e-14, 100]). |
| 9 | Lemma lambda (b) | ok (not sharp) | The bound 4/3 + 1/(1 - 4^{-s}) is valid. The numerical maximum of S is exactly 1/(1 - 4^{-s}): 1.3333 for s = 1 (bound 2.667) and 2.000 for s = 1/2 (bound 3.333). Valid for **every s > 0**, also s > 1. For s -> infinity (truncation), S <= 4/3 + 1. |
| 10 | s = 0 statement in Lemma lambda | minor | "S ~ L" and "S ~ l* ~ log_4(1/nu)" are right in order. More precisely, S is between (number of levels with t_l < 1)/2 and that number + 4/3. Numerically S(s = 0) = L in the unresolved limit, and 16.49 for L = 20 with 17 mass-dominated levels. |
| 11 | Lemma jackson (cited) | ok | True for Q1, Dirichlet, uniform dyadic grids; sources in section 3.1. Numerically C_Q(L) = 0.333, 0.406, ..., 0.4441, 0.4442 (d = 1, L = 1..11; converged) and 0.167, 0.274, 0.330, 0.384, 0.406, 0.423 (d = 2, L = 1..6; increments shrink by about 0.8 per level, so the limit is about 0.5). Dirichlet is essential: with Q_0 = 0, the l = 1 term needs Poincaré, and the lemma is false on Neumann spaces (constants). |
| 12 | Lemma scs (cited) and C_I = 4d 3^d | ok / minor | Strengthened Cauchy–Schwarz holds with rate 2^{-(l-k)/2}, and this rate is sharp in 1D. Numerically sup c(k,l) 2^{(l-k)/2} = 2.63 (d = 1, L = 10) and 3.46 (d = 2, L = 6) over k < l. The harmonicity argument is right: every d-linear function is harmonic. C_I = 4d 3^d is correct but loose. The sharp value is **C_I = 12d**, since h^2 lambda_max(K_l, M_l) < 12d; numerically 11.9999, 23.96, 34.98 for d = 1, 2, 3. It equals 4d 3^d only for d = 1. |
| 13 | Lemma H1 | ok | The diagonal and off-diagonal estimates and the Schur-test constant sum_{j>=1} 2^{-j/2} = 1/(sqrt2 - 1) are right. C_H = C_I + 2 C_S C_I^{1/2}/(sqrt2 - 1) is correct. Actual sup of nu|u|_1^2 / sum t_l ||u_l||^2: 22.7 (d = 1), 29.8 (d = 2, L = 5). The formula gives 56.0 (d = 1) and 213.6 / 105.8 (d = 2 with C_I = 72 / 24). |
| 14 | Lemma upper | ok | Mass part ||sum u_l||^2 <= S sum lambda_l ||u_l||^2; stiffness via H1 and (a). K1 = S_bound + C_H. Actual max a/N: 22.6 (d = 1, L = 10) and 29.5 (d = 2, L = 5), below the theoretical 58.7 and 216. |
| 15 | Lemma lower | ok | Q_{l*}u lies in V_{l*} and is a valid component; the components telescope to Q_L u = u. The l* = 1 case is fine, and l* = L reduces to u_L = u. In the unresolved case, lambda_L = 1 + t_L < 2. Valid for all s >= 0, because levels below l* get zero components. Actual max N/a (= K0): <= 1.28 for d = 1, 2, versus the theoretical K0 = max(5, 2 C_Q) = 5. |
| 16 | Theorem constants and kappa bound | ok | With computed constants the bound is about 2.6e3 (d = 1) and 5.8e4 (d = 2; 2.9e4 with C_I = 12d). The observed maxima are 6.14 (d = 1, L = 10) and 11.5 (d = 2, L = 5). The bound is valid but pessimistic by 10^2.6 to 10^3.7, as the draft says. |
| 17 | "same condition number bound for the unnormalized weights" | ok, but **gap in application** | True for kappa: C changes by the scalar theta_L. The paper, however, also uses c0 itself: eq. (naivealpha), Theorem B and Theorem heatB (kappa_alpha <= sqrt3 alpha_G^2 / c0), with alpha_G bounded independently of theta. With the paper's weights, c0 becomes theta_L c0 = min(1, t_L)^s c0, so the bounds on kappa_alpha degrade by t_L^{-s} in the unresolved regime. The paper's eq. (theta) should use normalized weights. |
| 18 | Remark s = 0 | ok (add ref.) | The growth is real, not just a weakness of the bound. Kolev–Xu–Zhu 2016, Lemma 4.4 and Remark 4.5 show lambda_max ~ |log h| cannot be improved for unweighted levels. Numerically K1(s = 0) = L in the unresolved limit; kappa(C A) for s = 0, d = 1, L = 10 grows 4.92, 7.33, 9.80, 12.16 as beta goes 1e-6 -> 1e-12. The YYZ comparison is right: their weights are s = 0 (see section 4). |
| 19 | Remark sharpness | ok | 3d 3^d = 243 for d = 3 is right. |
| 20 | Corollary heat (algebra) | ok | delta_{l,k} = r_k delta_l(nu_k) with nu_k = nu/r_k = sqrt(beta_k); theta_{l,k} = theta_l(nu_k); C_k = r_k^{-1} C(nu_k); C_k A_k = C(nu_k) A_{nu_k}. Checked numerically: spectra identical. |
| 21 | Corollary heat vs. the paper's definitions | **error as applied to the paper** | The corollary is correct for normalized weights, but the paper defines theta_{l,k} = min(1, sqrt(beta_k) h_l^{-2})^s, which is **unnormalized**. Then spec(C_k A_k) is contained in theta_{L,k} [c0, c1], not uniform in k once t_{L,k} = sqrt(beta_k) h_L^{-2} < 1. Since sqrt(beta_k) is about tau/2 for high frequencies, this happens when tau < about 2 h_L^2. Numerically (d = 2, L = 4, beta = 1e-2): r = 1e3 gives [9.2e-3, 6.6e-2] with paper weights versus [0.36, 2.57] normalized; r = 1e5 gives [6.7e-5, 7.4e-4] versus [0.26, 2.88]. The hypothesis of Theorem heatB ("c0, c1 independent of k") then fails. Fix: use theta~_{l,k}. |
| 22 | Paper: "Truncation does not work" | minor (paper) | The proof covers s -> infinity, i.e. dropping levels with t_l < 1, which is exactly Griebel–Oswald's robust truncated BPX. Truncation is therefore robust in theory; it only has worse constants at intermediate beta. Numerically kappa(trunc) = 19.23 vs 6.14 (s = 1) at d = 1, L = 10, beta = 1e-4, and 10.68 vs 2.90 at d = 2, L = 5, beta = 1e-4. It is better than s = 1 in the deeply unresolved 2D cases (8.67 vs 11.48). |

No hidden assumptions beyond those stated: nested uniform dyadic Q1 grids, coarse elements being unions of fine ones, homogeneous Dirichlet conditions (needed in the Jackson lemma at l = 1), and C_Q, C_S, C_I independent of L.

---

## 3. Status of the cited lemmas

V = verified by reading the text (preprint, scan or OCR). I = inferred from an abstract, a secondary source, or an argument.

### 3.1 Lemma jackson: sum_l h_l^{-2} ||(Q_l - Q_{l-1}) u||^2 <= C_Q |u|_1^2

**True as stated for Q1, uniform grids, Dirichlet.**

- **Griebel & Oswald, "Tensor product type subspace splittings and multilevel iterative methods for anisotropic problems", Adv. Comput. Math. 4 (1995) 171–206** [V, preprint text]. This is the best match for Q1.
  - Setting: bilinear elements on (0,1)^d, H^1_0, form int (grad u)^T C grad v + c0 u v.
  - Proof of Theorem 1, eq. (23)/(24): stability of {H^1_0(0,1)} = sum_j {W_{j,0}, a^{2j}(.,.)_{L2}} with W_j = V_j minus V_{j-1} (L2-orthogonal), via spline approximation and Besov characterizations. The tensor-product extension uses Props. 1–3.
  - With c1 = c2 = 1 and c0 = 0 this gives both directions of the norm equivalence. The text states that the result extends to d > 2 and c0 = 0.
- **Bornemann & Yserentant, "A basic norm equivalence for the theory of multilevel methods", Numer. Math. 64 (1993) 455–476** [V, scan]. The Jackson direction is Thm. 7.5, eq. (7.14), and Thm. 7.6, eq. (7.15); the upper direction is Thm. 4.1. The setting is simplicial P1. The proof uses only ||u - Q_k u|| <= c K(4^{-k}, u), a Fourier argument and zero extension, so it carries over verbatim to Q1 [I].
- **Xu, SIAM Rev. 34 (1992) 581–613** [V, OCR]. Lemma 6.7 (proof in Appendix, Prop. 8.6) gives sum ||(Q_k - Q_{k-1})v||_{H^1}^2 ~ ||v||_{H^1}^2 for P1, using H^{1+alpha} regularity. The h_k^{-2} form follows from the inverse estimate on V_k plus orthogonality. Xu credits Oswald and Zhang for the general case.
- **Dahmen & Kunoth, Numer. Math. 63 (1992) 315–344** [V, OCR]. Thm. 3.2 shows that Jackson plus Bernstein give O(1) conditioning. Section 6 (Props. 6.1/6.2) covers shift-invariant spaces on R^d including tensor-product B-splines, but not the Dirichlet cube.
- Oswald's book (1994) was not read.
- Numerical confirmation: C_Q about 0.444 (d = 1) and below 0.5 (d = 2), bounded in L.

Suggested citation in the draft: Griebel–Oswald 1995 (Adv. Comput. Math.), Thm. 1 and eq. (24); Bornemann–Yserentant 1993, Thm. 7.5; Oswald 1994.

### 3.2 Lemma scs: |(grad u, grad v)| <= C_S 2^{-(l-k)/2} |u|_1 h_l^{-1} ||v||, plus the inverse inequality

**True as stated.**

- **Xu, SIAM Rev. 34 (1992), Lemma 6.1** [V, page image]: a(u,v) <~ gamma^{j-i} h_j^{-1} ||u||_A ||v|| for u in M_i, v in M_j, i <= j, with gamma = 2^{-1/2} for dyadic refinement. This is exactly the draft's form.
  - The proof is Green's formula on coarse elements with Delta u = 0, plus a boundary-strip estimate. Xu credits Bank–Dupont and Yserentant, Numer. Math. 49 (1986).
  - For Q1, the step "grad u constant on K" is replaced by a polynomial trace-inverse estimate. Harmonicity still holds because every multilinear function is harmonic [I, checked by hand].
- Restated as **Kolev–Xu–Zhu, J. Sci. Comput. 67 (2016), Lemma 4.3** [V].
- BPX 1990 and Yserentant 1986 were not read: ams.org blocked the downloads.
- Inverse inequality: the sharp constant for Q1 on uniform grids is **C_I = 12d** (1D generalized eigenvalues (12/h^2) sin^2(phi/2)/(2 + cos phi), Kronecker sum).
- Numerics: C_S about 2.63 (d = 1) and 3.46 (d = 2); C_I = 11.9999, 23.96, 34.98 for d = 1, 2, 3.

### 3.3 Lemma lions (fictitious-space / additive Schwarz identity)

- **Griebel & Oswald, "On the abstract theory of additive and multiplicative Schwarz algorithms", Numer. Math. 70 (1995) 163–180**, eq. (18) and Thm. 1(b) [V, preprint]. They attribute it to Matsokin–Nepomnyaschikh 1985, Bjørstad–Mandel 1991, Xu 1992 and Zhang 1992.
- **Xu & Zikatanov, J. AMS 15 (2002) 573–597, Lemma 2.4** [I, via L. Chen's notes].
- P.-L. Lions 1988 could not be confirmed to contain this identity. Prefer "Griebel–Oswald 1995; Xu–Zikatanov 2002" over "Lions".

---

## 4. Novelty relative to the literature

The result, i.e. robust additive multilevel preconditioning of -eps Delta + I uniformly in eps and h, is **known**.

- **Griebel & Oswald, Adv. Comput. Math. 4 (1995), Theorems 1 and 2** [V]:
  - Thm. 1, eq. (19), with c1 = c2 = eps and c0 = 1 gives exactly ||u||^2 + eps|u|_1^2 ~ sum_l (1 + eps 4^l) ||(Q_l - Q_{l-1})u||^2, uniformly in eps.
  - Thm. 2, eqs. (26) and (27): nodal-basis BPX and multilevel diagonal scaling (MDS, i.e. Jacobi with respect to the full form) for Q1 on (0,1)^d, with condition numbers independent of c_k and J. Levels coarser than j0 = [log_4(1/eps)], i.e. h_l^2 > eps (mass-dominated), are **dropped**.
  - This is the truncated version (s = infinity) of the draft's theorem. GO also note that the full splitting without this choice is not robust.
- **Bramble, Pasciak & Vassilevski, Math. Comp. 69 (2000) 463–480** [I: abstract and citation in Bacuta et al.]: B = sum_j (eps h_j^{-2} + c)^{-1} (Q_j - Q_{j-1}) is a uniform preconditioner for -eps Delta + c.
- **Bacuta, Hayes & Jacavage, arXiv:2104.11337 (Comput. Math. Appl. 2022)** [V, Sec. 2 and eq. (8.4)]:
  - By summation by parts, B = gamma_J I + sum_{j<J} (gamma_j - gamma_{j+1}) Q_j with gamma_j = (eps h_j^{-2} + c)^{-1}.
  - Replacing Q_j by lumped-mass "Q~_j" gives B = gamma_J D_J^{-1} + sum_j (gamma_j - gamma_{j+1}) E_j D_j^{-1} E_j^T, which is uniform by their Lemmas 2.1–2.2.
  - In our scaling, gamma_j - gamma_{j+1} = 3 t_j / ((1 + t_j)(1 + 4 t_j)), which is about 3 t_j for t_j << 1 and about 3/(4 t_j) for t_j >> 1.
  - In the resolved regime these level weights are within level-independent factors of the draft's s = 1 weights theta_l/delta_l: t_l for mass-dominated levels, about 1/t_l for the others. The two preconditioners are then spectrally equivalent there. In the unresolved regime, BPV/BHJ damp the coarse levels by an additional factor t_L; both choices are robust.
  - **So s = 1 is, up to Jacobi versus lumped mass and normalization, the known BPV/BHJ preconditioner.**
- **Kolev, Xu & Zhu, J. Sci. Comput. 67 (2016), Lemma 4.4 and Remark 4.5** [V]: with naive weights (omega h_k^{-2} + rho)^{-1} on V_k (the s = 0 case), lambda_max(BA) <~ |log h|, and this cannot be improved. This confirms the draft's s = 0 remark.
- **Multiplicative and AMLI methods** [I, abstracts only]: Olshanskii–Reusken, Computing 65 (2000) 193–202 (robust V/W-cycle); Kraus–Wolfmayr (robust AMLI for reaction–diffusion).
- **Yang–Yu–Zhang, arXiv:2610.06987, Theorem 2.1** [V, full text]:
  - Operator -Delta + c. B = sum omega_l^2 P_l M_l^{-1} P_l^T with omega_l^2 = (h_l^{-2} + c)^{-1}.
  - In our notation this is **s = 0 with exact M_l^{-1}** instead of Jacobi. Mass-dominated levels are not damped.
  - The bound is 1/L <= lambda <= 12dL, so kappa <= 12dL^2 uniformly in c.
  - The proof sums the per-level estimate ||v - Q_l v|| <= (h_l/pi)||grad v|| (no Jackson sum) and uses Cauchy–Schwarz over the L levels (no strengthened Cauchy–Schwarz). Both ends lose a factor L, even at c = 0 where BPX is O(1).
  - Their operator satisfies B K_h = C_N A for s = 0, so its exact condition number is K0 K1(s = 0). Numerically it is <= 8.9 (d = 1, L = 10) and <= 9.8 (d = 2, L = 5), compared with their bounds 1200 and 600.
  - With Lemmas jackson and scs, their bound improves to O(min(L, log(1/nu)) + C_H), and with damped weights to O(1).

**What remains new in the draft:**
- the one-parameter family theta^s with s in (0,1] (actually any s > 0), in the Jacobi-scaled nodal form needed for the block encoding;
- a short self-contained proof with explicit constants;
- the normalization that makes [c0, c1] beta-independent, which the heat problem needs.

The paper should cite GO95 (Adv. Comput. Math.), BPV 2000 and Bacuta et al., and present the theorem as an adaptation, not a new robustness result.

---

## 5. Numerical checks

Setup: `.venv/bin/python`, dense linear algebra, d = 1 (L <= 11) and d = 2 (L <= 6).

1. **Lions identity**: x^T C^{-1} x versus the minimum-norm decomposition via lstsq. Relative error <= 6e-15 for d = 1 (L = 6) and d = 2 (L = 4), nu in {1, 1e-3, 1e-6}, s in {0, 1/2, 1}.
2. **Embedding**: P_l^T M_L P_l = M_l and P_l^T K_L P_l = K_l to 4e-16; P_L = I (d = 1, 2, 3).
3. **Lemma nodal**: eig(M_l)/h^d lies in [0.3349, 0.9984] (d = 1, l = 5), [0.1122, 0.9968] (d = 2) and [0.0462, 0.9258] (d = 3, l = 3), inside (3^{-d}, 1). The formula for delta_l is exact, and delta/((2h/3)^d mu) lies in [1, 3d].
4. **Lemma lambda** (L <= 20, nu in [1e-14, 100], 161 values): t_l <= lambda_l always; t_{l*} < 4 always.

   | s | max S | bound |
   |---|-------|-------|
   | 1 | 1.3333 | 2.667 |
   | 1/2 | 2.0000 | 3.333 |
   | 0 | 19.99 (L = 20) | none (grows like L) |

5. **Jackson constant C_Q(L)** (generalized eigenvalue of sum 4^l (G_l - G_{l-1}) against K):

   | d | L = 1 | 2 | 3 | 4 | 5 | 6 | 8 | 10 | 11 |
   |---|---|---|---|---|---|---|---|---|---|
   | 1 | .3333 | .4064 | .4296 | .4378 | .4410 | .4425 | .4437 | .4441 | .4442 |
   | 2 | .1667 | .2740 | .3297 | .3840 | .4058 | .4234 | | | |

   Bounded and saturating. The draft's K0 = max(5, 2 C_Q) is 5.
6. **SCS and inverse constants**: C_S (k < l) = 2.63 (d = 1, L = 10) and 3.46 (d = 2, L = 6). In 1D, c(k,l) 2^{(l-k)/2} is flat at about 2.62 for l - k = 2..6, so the rate 1/2 is sharp. C_I = 11.9999, 23.96, 34.98 (sharp 12d; the draft's 4d 3^d gives 12, 72, 324).
7. **K0 = max N/a and K1 = max a/N** (N computed exactly via C_N = sum P_l (lambda_l M_l)^{-1} P_l^T), over beta in {1, ..., 1e-12}:

   | case | max K0 | max K1 | max K0 K1 | theory K0 | theory K1 |
   |------|--------|--------|-----------|-----------|-----------|
   | d = 1, L = 10, s = 1 | 1.28 | 22.6 | 19.6 | 5 | 58.7 |
   | d = 2, L = 5, s = 1 | 1.26 | 29.5 | 20.1 | 5 | 216 (108 with C_I = 12d) |
   | d = 2, L = 5, s = 0 | 0.998 | 29.5 | 9.8 | 5 | grows like L; = L in the unresolved limit |

8. **kappa(C A_beta)**: d = 1, L = 10 and d = 2, L = 4 reproduce Table 2 of the paper to all printed digits (e.g. 6.14 for d = 1, s = 1, beta = 1e-4; 10.51 for d = 2, L = 4, s = 1, beta = 1e-10). The theorem's bound kappa <= 3d 3^d K0 K1 evaluates to 2641 (d = 1, s = 1) and 58,390 (d = 2, s = 1; 29,272 with C_I = 12d), versus observed maxima 6.14 and 11.48. The bound is valid but pessimistic by a factor of 400 to 5000.
9. **Unresolved floor**: kappa(C M) with nu -> 0 converges to 3^d/(1 - 4^{-s}), not 3^d:

   | d | s | L | kappa(C M) | 3^d/(1 - 4^{-s}) |
   |---|---|---|------------|------------------|
   | 1 | 1 | 10 | 4.000 | 4 |
   | 1 | 1/2 | 10 | 5.990 | 6 |
   | 2 | 1 | 6 | 11.93 | 12 |
   | 2 | 1/2 | 6 | 17.31 | 18 |

   For comparison, kappa(diag(M)^{-1} M) = 3.000 (d = 1) and 8.97 (d = 2).
10. **Heat corollary** (r in {1, 10, 1e3, 1e5}, beta in {1e-2, 1e-6}): C_k A_k with normalized weights has exactly the spectrum of C(nu_k) A_{nu_k}, and it stays inside a k-independent interval (d = 2, L = 4: within [0.26, 2.89]). With the paper's unnormalized theta_{l,k}, the interval collapses by theta_{L,k}, e.g. to [6.6e-7, 7.4e-6] for d = 2, beta = 1e-6, r = 1e5.
11. **Truncation versus s = 1** (kappa):

    | case | beta | s = 0 | s = 1 | s = 2 | truncation |
    |------|------|-------|-------|-------|------------|
    | d = 1, L = 10 | 1e-4 | 5.02 | 6.14 | 9.72 | 19.23 |
    | d = 1, L = 10 | 1e-10 | 9.80 | 3.61 | 6.13 | 13.14 |
    | d = 2, L = 5 | 1e-4 | 3.96 | 2.90 | 4.86 | 10.68 |
    | d = 2, L = 5 | 1e-12 | 36.57 | 11.48 | 9.25 | 8.67 |

    Truncation is bounded, as the theory says, but worse at intermediate beta.

---

## 6. Suggested corrections

**(a) Lemma lambda, extend the range of s.** Replace "For $s\in(0,1]$" in Lemma lambda, Lemma upper and the Theorem by "For $s>0$". Add after the proof of Lemma lambda:

> The same proof applies verbatim for every $s>0$. Letting $s\to\infty$ (i.e. $\tilde\theta_\ell=0$ for $t_\ell<1$ when $t_L\ge1$, and only level $L$ when $t_L<1$) gives $S\le 4/3+1$. Hence the truncated multilevel method of Griebel and Oswald is covered too.

Numerically, the maximum of S is 1/(1 - 4^{-s}); the stated bound remains valid.

**(b) Lemma scs note, sharper inverse constant.** Replace "the inverse inequality holds with $C_I=4d\cdot3^d$, by the eigenvalue bounds ... ; an elementwise argument gives a smaller constant" with:

> For $Q_1$ on a uniform grid, $|v|_{H^1}^2\le 12d\,h_\ell^{-2}\|v\|^2$. In 1D the generalized eigenvalues of $(K_1,M_1)$ are $\frac{12}{h^2}\frac{\sin^2(\phi/2)}{2+\cos\phi}<\frac{12}{h^2}$, and in $d$ dimensions they are sums of $d$ such values. So $C_I=12d$, which is sharp as $h\to0$.

Then C_H = 12d + 2 C_S sqrt(12d)/(sqrt2 - 1).

**(c) Sentence after (eq:NC).** Replace "The factor $3d\cdot3^d$ comes only from Lemma~\ref{lem:nodal}, i.e. from Jacobi scaling and the nodal basis. It is the ``mass-matrix floor'' seen in the numerics." with:

> The factor $3d\cdot3^d$ comes only from Lemma~\ref{lem:nodal}, i.e. from Jacobi scaling and the nodal basis. The comparison~\eqref{eq:NC} itself is nearly sharp, but its extremes are not attained by the same functions as those of $N(u)/a(u,u)$, so the actual $\kappa(CA_\beta)$ is much smaller. In the unresolved limit $\nu\to0$ one observes $\kappa(CM)\to3^d/(1-4^{-s})$.

**(d) Lemma lambda, s = 0 sentence.** Replace with:

> For $s=0$, (a) still holds, and $\tfrac12m\le S\le m+\tfrac43$, where $m=\#\{\ell:t_\ell<1\}$, i.e. $m=L$ in the unresolved regime and $m=\ell^*-1=\lceil\log_4(1/\nu)\rceil-1$ in the resolved one.

**(e) Remark s = 0.** Add:

> The growth is genuine: for unweighted levels, $\lambda_{\max}(CA)\gtrsim m$ (test $u\in V_1$ split equally over the $m$ mass-dominated levels; cf. Kolev--Xu--Zhu, J. Sci. Comput. 67 (2016), Remark 4.5). Yang--Yu--Zhang (arXiv:2610.06987, Thm. 2.1) use $\omega_\ell^2=(h_\ell^{-2}+c)^{-1}$ with $M_\ell^{-1}$, i.e. our $s=0$ with exact mass inverse. Their bound $12dL^2$ loses one factor $L$ at each end through per-level estimates. With Lemmas~\ref{lem:jackson} and~\ref{lem:scs} it improves to $O(\min(L,\log(1/\nu))+C_H)$.

**(f) Theorem, last sentence.** Replace "The same condition number bound holds for the unnormalized weights $\theta_\ell$ of the paper." with:

> For the unnormalized weights $\theta_\ell$, $\operatorname{spec}(CA_\beta)\subset\theta_L[c_0,c_1]$ with $\theta_L=\min(1,t_L)^s$. The condition number bound is the same, but $c_0$ degrades by $t_L^{s}$ in the unresolved regime. Since the cost bounds of the paper involve $c_0$ itself (with $\alpha_G$ independent of $\theta$), the paper should use the normalized weights $\tilde\theta_\ell\le1$. They leave Proposition~\ref{prop:alphaG} unchanged.

**(g) Corollary heat.** Add:

> Note that the paper's $\theta_{\ell,k}=\min(1,\sqrt{\beta_k}h_\ell^{-2})^s$ must be replaced by $\tilde\theta_{\ell,k}=\big(\min(1,\sqrt{\beta_k}h_\ell^{-2})/\min(1,\sqrt{\beta_k}h_L^{-2})\big)^s$. Otherwise $\operatorname{spec}(C_kA_k)\subset\theta_{L,k}[c_0,c_1]$, which is not uniform in $k$ when $\sqrt{\beta_k}h_L^{-2}<1$ (high frequencies with $\tau\lesssim2h_L^2$).

**(h) Corresponding edits in `quantum_preconditioning.tex`** (not made; for the authors):
- **eq. (theta):** use theta~_l = (min(1, sqrt(beta) h_l^{-2}) / min(1, sqrt(beta) h_L^{-2}))^s. Also, s in {0, 1} is inconsistent with Table 2, which uses s = 1/2.
- **Heat section:** make the same change to theta_{l,k}.
- **Bullets after Assumption 3.4 and Sec. 6:** replace "approaches the constant $3^d$ of the Jacobi-scaled $Q_1$ mass matrix" with "approaches about $3^d/(1-4^{-s})$ (4 in 1D and 12 in 2D for $s=1$)". The worst 3D value 20.5 will grow towards about 36 for larger L at tiny beta.
- **"Truncation does not work":** replace with "Truncation (dropping the mass-dominated levels) is also robust in theory [GO95, Thm. 2], but its constants are worse at intermediate $\beta$ (e.g. 19.2 vs 6.1 for $d=1$, $L=10$, $\beta=10^{-4}$)."
- **Prior work:** cite Griebel–Oswald, Adv. Comput. Math. 4 (1995) 171–206, Thms. 1–2; Bramble–Pasciak–Vassilevski, Math. Comp. 69 (2000) 463–480; Bacuta–Hayes–Jacavage (arXiv:2104.11337), eq. (2.6)/(8.4). Present Theorem robust as an adaptation of these to the Jacobi-scaled factored form, with s-family and explicit constants.

**(i) Sources for the cited lemmas** (replace "Sources to check"):
- Lemma jackson: Griebel–Oswald, Adv. Comput. Math. 4 (1995), proof of Thm. 1, eqs. (23)–(24) (Q1, Dirichlet, (0,1)^d); Bornemann–Yserentant, Numer. Math. 64 (1993), Thm. 7.5 (P1; the proof transfers); Oswald 1994; Xu 1992, Lemma 6.7.
- Lemma scs: Xu, SIAM Rev. 34 (1992), Lemma 6.1 (after Yserentant 1986); the proof transfers to Q1 because multilinear functions are harmonic.
- Lemma lions: Griebel–Oswald, Numer. Math. 70 (1995), eq. (18) and Thm. 1(b); Xu–Zikatanov, JAMS 15 (2002), Lemma 2.4.
