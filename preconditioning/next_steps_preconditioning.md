# Next steps: "Preconditioning for quantum linear solvers in PDE-constrained optimization"

For Martin and John. Status as of 2026-11-04 (draft at commit `23e9124` plus uncommitted edits).

## Where we stand

**Proven and checked to machine precision**
- Schöberl–Zulehner bound for Poisson control. This is known; we use it.
- Two-sided theorem (Thm 3.5): the nonzero spectrum of 𝓑 = R̃ᵀ𝒜R̃ lies in ±[c₀/√2, c₁].
- J-factorization 𝓑 = 𝒢ᵀ𝒥𝒢, with 𝒥 a signed permutation (one two-qubit gate).
- Commuting lemma: F·P_ℓ = T_ℓ·F_ℓ with an element-local isometry. Checked for Q1–Q3 and for P1/P2 on nested meshes.
- Time-periodic heat equation: QFT in time, √3 bound per frequency, and the coupling matrix as a scaled reflection with ‖𝒥_k‖ ≤ 1.
- BPX robustness in β (Theorem 3.5, normalized weights, all s > 0; Q1, uniform grids). Proof in `paper/proof_bpx_robust.tex`, agent-checked; essentially known (Griebel–Oswald 1995), new is the factored form, s-family, normalization, constants.

**Only numerical so far**
- 3D grids up to n ≈ 30k.

**Only argued so far**
- The subnormalization α_G ≤ C_d√L (Prop. 3.8, proof sketch with a pessimistic constant C₃ ≈ 45).
- Hence the main Theorem 3.9 is still a *target*.
- No end-to-end cost statement yet.

**Known problem**
- In the heat case, the weighted BPX makes the condition number over *all* frequencies worse than the worst single frequency (d = 1, L = 6: 7.5 vs 4.6).
- Normalizing by the finest-level weight does **not** fix it (tested).
- Cause: the spectral bounds c₀,ₖ and c₁,ₖ drift with the frequency k, and QSVT sees max_k c₁,ₖ / min_k c₀,ₖ.

## Target

A paper for SISC or SIMAX (NLA framing) whose core claims are *theorems with explicit constants*, backed by measured subnormalizations and an honest end-to-end cost. A venue like Quantum or PRX Quantum only makes sense if Phase 2 includes resource counts.

Realistic effort: about 2–3 months part time. Suggested leads are proposals, not assignments.

---

## Phase 0: housekeeping (1–2 days)

- [x] **Rescue the Q_k scripts.** Now in `code/` as `qp_elements_qk.py`, `qp_elements_qk_extra.py` and `qp_elements_q2_kron.py` (2026-11-04).
- [ ] **Rescue the triangle scripts** (P1/P2 on red-refined meshes). Still only in the temporary scratchpad (`tri.py`, `tri_kron.py`).
- [ ] **Verify the bibliography.** Entries written from memory: KK12 (year), KKLWZ13 (pages), the QIPM review (authors), BNWA23 (journal version, PRR 2025?). Also check every arXiv number.
- [x] **Serious prior-art check** done (2026-11-04): `prior_art_preconditioning.md`. Verdict: go, but reposition. The commuting lemma for Q1 (Deiml–Peterseim, Yang–Yu–Zhang), the β-uniform BPX bound (Yang–Yu–Zhang, up to L²) and the √3 bound for λ = iω (Krendl–Simoncini–Zulehner 2013) are known. Still new: the two-sided KKT preconditioning, the 𝒢ᵀ𝒥𝒢 J-factorization, and the quantum time-periodic treatment. **Rewrite the contributions and the "Gap" paragraph accordingly.** Topics covered:
  - quantum KKT and saddle-point solvers;
  - quantum optimal control and PDE-constrained optimization;
  - preconditioned block encodings since 2025;
  - follow-ups to Deiml–Peterseim.
  - *Go/no-go:* if someone has done saddle-point preconditioning for QLSAs, reposition before investing more.
- [x] `qk_rational_fair.py` moved, unfinished, to `../archive/rational_krylov/` together with all rational-Krylov material (2026-11-04). The rebuilt `quantum_krylov_note.pdf` was committed.
- [ ] Give John access to the GitHub repo.

## Phase 1: make the core claims hard (3–5 weeks)

### 1.1 Measure α_G (most important; suggested lead: Martin, with Claude for scripts)
- [ ] Build the block encoding of G as an explicit unitary matrix for d = 1 (L ≤ 6) and d = 2 (L ≤ 4), using exactly the factors a circuit uses:
  - LCUs of shifts;
  - Kronecker products;
  - the weighted stack over values and gradients;
  - the isometries T_ℓ, completed to unitaries;
  - the concatenation over levels.
- [ ] Check unitarity and that the leading block equals G/α_G. Tabulate α_G against ‖G‖ for varying L, d, β.
- [ ] Decide whether the √L factor is real. If so, try a weighted concatenation or a strengthened Cauchy–Schwarz argument to remove it.
- [ ] **Done when:** Prop. 3.8 becomes a theorem with an explicit constant, and Table 1 reports κ_α = α_G²/σ_min instead of the ideal ‖G‖²/σ_min.

### 1.2 β-robustness of the multilevel factor (suggested lead: John and Martin)
- [x] Literature: known. Griebel–Oswald, Adv. Comput. Math. 4 (1995), Thms 1–2 (truncated, s → ∞); Bramble–Pasciak–Vassilevski 2000; Bacuta–Hayes–Jacavage 2021 (≈ our s = 1). Plain BPX log growth is genuine (Kolev–Xu–Zhu 2016).
- [x] Proof written for all s > 0 with normalized weights (`paper/proof_bpx_robust.tex`), checked by an agent (`paper/proof_bpx_robust_review.md`). Paper now states it as Theorem 3.5 and uses normalized weights θ_ℓ/θ_L.
- [ ] Read the sources ourselves: GO95 eqs. (23)–(24) for the Jackson lemma; Xu 1992 Lemma 6.1 and the Q1 transfer of strengthened Cauchy–Schwarz.
- [ ] Decide whether the proof goes into the paper (appendix) or stays a supplement.
- [ ] **Mass-matrix floor** 3^d/(1−4^{-s}) (4, 12, 36 for s = 1), from Jacobi scaling of the Q1 mass matrix; it matters for the heat problem when τ ≲ 2h². Remedy found, see `paper/nla_for_quantum.tex`, section "Removing the mass-matrix floor", and `code/qp_mass_floor.py`:
  - level scaling D_ℓ = (M_ℓ + √β diag(K_ℓ) I)⁻¹ instead of δ_ℓ⁻¹ (floor → 4/3; resolved cases unchanged);
  - encoding via the isometry U_val = F_val M^{-1/2} and sine-basis rotations, as in Yang–Yu–Zhang Secs. 3.4–3.5 (prior art for the encoding); measured κ_α: 3D unresolved 24.5 → 1.3, resolved within ~15% of Jacobi. Tight value encoding via U_val also improves the current Jacobi route by up to 2.4× in 3D (relevant for item 1.1).
  - [ ] Adapt the robustness proof to D_ℓ.
  - [ ] Gate count per call to G (sine transforms, angle arithmetic, 3-call amplification) vs the gain in κ_α.
  - [ ] Decide whether the paper switches to D_ℓ or keeps Jacobi and mentions D_ℓ as a remedy.
- [ ] Rerun 3D at larger L / smaller β to see the floor approach 36.

### 1.3 Fix the cross-frequency scaling for the heat equation (suggested lead: Martin)
- [ ] Add per-frequency scalars that align the spectra of the frequency blocks, e.g. scale 𝓑_k by an a-priori estimate of 1/c₁,ₖ (Gershgorin or Fourier-symbol bound). The global κ then ≈ √3 · max_k(c₁,ₖ/c₀,ₖ).
- [ ] Check that the scalars are cheap on the device: diagonal in k, a smooth function of sin(πk/N_t).
- [ ] Rerun `qp_heat_periodic_check.py freq` (d = 1, L = 6 and 8; d = 2, L = 4) and fill Table 5.
- [ ] Restate Thm 4.3 with the per-frequency scaling.

### 1.4 Attribution of the heat bound (suggested lead: John)
- [ ] Read Kollmann–Kolmbauer (NLA 2012) and Kollmann–Kolmbauer–Langer–Wolfmayr–Zulehner (CMA 2013). Is our √3 bound for periodic backward Euler new, or a special case? Attribute accordingly.

### 1.5 Small additions
- [ ] Q2 numerics in the paper: Kronecker check, commuting lemma, κ(CA_β). Currently only in a remark.
- [ ] Rerun the killed jobs: 3D L = 5 with β = 10⁻⁸; 2D L = 8, 9 (`qp_bpx_sparse.py`).

## Phase 2: end-to-end cost and positioning (3–5 weeks)

### 2.1 End-to-end theorem (suggested lead: Martin)
- [ ] Preparing g (weighted ŷ at Gauss points): which classes of ŷ allow polylog(n) circuits? Cite or prove.
- [ ] Success probability of QSVT postselection, amplitude amplification, and estimation of ‖z‖.
- [ ] QoI readout (amplitude estimation, O(1/ε)).
- [ ] Error budget across polynomial approximation, discretization and readout. A single theorem: total cost to compute a QoI to accuracy ε.

### 2.2 Analytic resource counts (suggested lead: Martin with Claude)
- [ ] Toffoli and rotation counts per building block as formulas in L, d, ε: shifts and comparators about ℓ Toffolis; ramp state; level control; QFT for heat.
- [ ] Multiply out with the measured α_G from 1.1. Give a table for d = 2, 3 and realistic L, plus wall-clock estimates against classical multigrid–MINRES.
- [ ] Expected outcome: no practical advantage at feasible sizes. Say so, and identify the dominant constant.
- [ ] Compiled counts (Qualtran / Azure Resource Estimator) only if we target a quantum venue.

### 2.3 Positioning and baselines (suggested lead: John)
- [ ] **Reduced versus full KKT system.** Explain why the all-at-once saddle-point system is preferable on a quantum device: the reduced Hessian βM + M K⁻¹ M K⁻¹ M needs nested inverses, whose costs multiply. This is a strong argument for our approach and missing from the draft.
- [ ] One-sided preconditioning, or applying Deiml–Peterseim block by block, as baselines. Why do they lose?
- [ ] Reconcile our "advantage only for d = 3" count with Deiml–Peterseim's "from d = 2".

**Decision point after Phase 2:** choose the venue (SISC/SIMAX or Quantum), depending on how strong 2.1–2.2 turn out.

## Phase 3: writing (2–3 weeks)

- [ ] **Underexplained for NLA readers:**
  - sparse block encodings ("entries computable by a circuit");
  - LCU in general (one example now in Sec. 2);
  - the role of the success probability and amplitude amplification;
  - why concatenating levels costs √L;
  - |D_t| and the √r_k weight in the heat section.
- [ ] Figures:
  - a diagram of the factorization 𝒢ᵀ𝒥𝒢 (what each factor does);
  - κ against L and β (instead of only tables);
  - a register/circuit sketch.
- [ ] Possibly an appendix: "quantum linear algebra for numerical analysts" (reuse from the quantum Krylov note).
- [ ] Reproducibility: README for the `qp_*.py` scripts, fixed versions, a script-to-table map.
- [ ] Shorten the Introduction and Contributions once the results are final. Check notation consistency.

## Phase 4: not for this paper (follow-ups)

- **Control constraints and interior point methods** (John's area): the (2,2) block changes every Newton step. This connects to quantum IPMs. A natural second paper.
- **Initial-value heat problem.** The α-circulant preconditioner trades preconditioner quality against encoding cost (κ/α). A toy model shows a bounded optimum at α ≈ 0.1–0.25, where the classical choice is α ≈ 10⁻⁸. This is a genuinely new NLA question, and a strong standalone paper.
- Non-Cartesian and unstructured meshes (sparse-access oracles), variable coefficients, boundary control.

## Risks

| risk | consequence | mitigation |
|---|---|---|
| Prior art found in Phase 0 | the novelty claim weakens | reposition early; the heat and J-factorization parts are the most likely to stay new |
| α_G much larger than ‖G‖ | the main theorem has bad constants | change the construction (Kronecker stack instead of sparse access); report honestly |
| End-to-end cost shows no advantage | weaker "quantum" story | frame the paper as "the NLA that makes quantum solvers well posed" (already the framing) |

## Files

All paths are relative to this folder (`preconditioning/`).

- Paper: `paper/quantum_preconditioning.tex` / `.pdf`; supplement `paper/nla_for_quantum.tex` (how block encodings are realized)
- Checks (in `code/`; keep them together, they import each other):
  - `qp_poisson_check.py` (Tables 1, 4)
  - `qp_bpx_weights.py` (Table 2)
  - `qp_bpx_sparse.py` (Table 3)
  - `qp_heat_periodic_check.py` (Table 5; `global` / `freq`)
  - `qp_elements_qk.py`, `qp_elements_qk_extra.py` (Q_k generality: factorization, commuting lemma, isometry, under-integration, mapped geometry; Remark 3.1)
  - `qp_elements_q2_kron.py` (direct 2D Q2 assembly against the Kronecker form)
  - `qp_block_encoding_rules.py` (explicit unitary realizations of all block-encoding rules; supplement `paper/nla_for_quantum.tex`)
- Background: `../krylov_note/` (note, `dead_ends.md`), `../archive/rational_krylov/` (rational Krylov, parked; `prior_art.md`)
