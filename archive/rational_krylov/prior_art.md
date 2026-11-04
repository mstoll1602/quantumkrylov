# Prior-art check: rational Krylov / Loewner methods from quantum time-series data

Checked on 2026-11-03 against the idea in `manuscript_outline.md`: use
resolvent filters $(z_k - H)^{-1}v_0$, build the Loewner pencil from
$r(z) = v_0^*(z-H)^{-1}v_0$, and estimate $r$ as a Laplace transform of
$g(t) = v_0^*e^{-itH}v_0$ by random-time Monte Carlo.

**Verdict.** The core mechanism is **not new classically**: a basis of filtered
states at chosen (possibly complex) energies, with matrices computed from the
autocorrelation signal alone, is classical *filter diagonalization* (FD). FD has
also been brought to quantum computers, most recently with prolate windows and
error bounds. A paper is still possible, but only with a narrower claim
(Section 4 below). Reviewers will read our method as "filter diagonalization
with Lorentzian windows", and we must position against that from the start.

---

## 1. Quantum filter diagonalization (Parrish–McMahon, arXiv:1909.08925, 2019)

- The basis is **plain time-propagated states** $e^{-i2\pi k H/\kappa}|\Phi\rangle$,
  $k=-k_{\max},\dots,k_{\max}$ (their Eq. 5), i.e. real-time Krylov, possibly
  with several reference states. "Lorentzian" appears only in the line
  broadening of their plotted spectra.
- They state that classical FD uses "a basis of either time-propagated
  reference states or a Chebyshev expansion of $(H-\lambda)^{-\alpha}$ acting on
  reference states". They implement the first variant and relate the second to
  quantum inverse iteration.
- **Overlap with us:** low. QFD is real-time Krylov; it does not use resolvent
  filters.

## 2. Classical filter diagonalization (Neuhauser 1990/1994; Wall–Neuhauser 1995; Mandelshtam–Taylor 1995/1997; review Mandelshtam 2001)

From the review (Prog. NMR Spectrosc. 38 (2001) 159–196, Eqs. 66–74):

- **Fourier basis at reference frequencies $w_j$**:
  $|w_j) = \sum_{n=0}^{M-1} e^{in\tau(w_j-\Omega)}|0)
        = \frac{1-e^{iM\tau(w_j-\Omega)}}{1-e^{i\tau(w_j-\Omega)}}\,|0)$.
  This is a **truncated resolvent of the unitary** $U=e^{-i\tau\Omega}$. For
  complex $w_j$ (damping) and $M\to\infty$ it becomes
  $(I - z_jU)^{-1}|0)$, i.e. rational Krylov **in $U$**.
- **Matrices from the signal alone.** The matrix elements of $U^p$ in this
  basis are linear functionals of $c_n$ (Eq. 72), with denominators
  $1-e^{i\tau(w_{j'}-w_j)}$. These are divided differences, i.e. the
  Loewner-type structure we derived.
- Using $U^{(0)}$ and $U^{(1)}$ gives a unitary pencil, so **no $H$ is
  needed**.
- Regularization exists: the **Regularized Resolvent Transform (RRT)**.
- FD belongs to the harmonic-inversion family together with Padé (Belkić's
  "fast Padé transform"), linear prediction and Prony. Their relative accuracy
  depends on the setting (Main, Dando, Belkić, Taylor 2000; arXiv:nlin/0104064).
- Neuhauser's original FD filtered states at regularly spaced reference
  energies; resolvent/Green's-function filters $(H-\lambda)^{-1}$ appear as a
  variant.
- **Overlap with us:** high for the mechanism. FD is discrete-time (z-transform
  on a grid, finite window); we use a continuous-time Laplace transform, i.e.
  the resolvent of $H$ rather than of $U$.

## 3. Filter diagonalization on quantum computers

- **Quantum prolate diagonalization (QPD)**, Stroschein, Castaldo, Reiher,
  arXiv:2507.15148, *Phys. Rev. Research* 8 (2026) 033177. An eigenvalue
  problem built from the autocorrelation function with **prolate spheroidal
  windows**, rigorous **error bounds** (a sharp accuracy transition set by the
  observation time and spectral density), Heisenberg-limited accuracy for
  ground and excited states, and robustness to imperfect state preparation.
  **Closest competitor. Read it fully before writing.**
  Method details, from the HTML version:
  - Trial vectors $v_l=\int f_l(\tau)\,e^{-i\tau H}v_0\,d\tau$, with prolate
    spheroidal windows $f_l$ on $[-T,T]$.
  - $(V^*V)_{sl}=\iint f_s(\tau)C(\tau-t)f_l(t)\,dt\,d\tau$, and
    $V^*HV$ is the same integral with $-i\,\partial_\tau C$. Both come from
    the signal only.
  - The signal is sampled on a Nyquist grid $t_k=\pi k/W_s$ over
    $[-2T,2T]$.
  - No statistical shot-noise model. No resolvents, Laplace transforms,
    complex energies, Loewner pencils or rational approximation.
  - Error bounds come from prolate band concentration.

  **In these terms, our method is FD with one-sided exponential windows**
  $f_k(t)=e^{iz_kt}$ on $[0,\infty)$: $\int_0^\infty e^{iz_kt}e^{-iHt}v_0\,dt=i(z_k-H)^{-1}v_0$.
  - **Risk:** prolates are *optimally* concentrated for a given $T$, while
    Lorentzian windows have heavy $1/(E-a)^2$ tails. QPD may well beat us in
    accuracy per unit of evolution time.
  - **Our possible edge:** a closed-form Loewner pencil, rational
    approximation theory, no grid, a shared-shot random-time estimator with
    an explicit shot-noise analysis, and decoherence as a pole shift.
  - **Must-have experiment:** a head-to-head comparison with QPD at equal
    total evolution time and equal shots.
- **Cortes–Gray** (arXiv:2109.06868, Phys. Rev. A 2022) relate real-time
  quantum Krylov to FD with a non-orthogonal basis of "filter energies" and
  count the quantum calls needed.
- **Quantum inverse iteration**, Kyriienko (arXiv:1901.09988, npj Quantum Inf.
  2020). It writes $(H-\lambda)^{-1}$ as a Fourier-type sum of time
  evolutions and uses it in a **power iteration** (one shift, no subspace or
  Loewner pencil).

## 4. Fault-tolerant contour and rational methods (different setting)

- Usuda (QIC 2020); Takahira, Ohashi, Sogabe, Usuda (arXiv:2106.08076);
  Jiang–An (arXiv:2601.11959): $f(A)|\psi\rangle$ via a Cauchy integral, with
  resolvents from quantum linear solvers or QSVT, combined by LCU.
- Quantum rational transformation via LCHS (arXiv:2408.07742): resolvents
  from integrals of Hamiltonian simulations, Zolotarev sign function.
- Futamura–Ye–Sakurai (arXiv:2112.05395): eigenvalue density via a contour
  integral with HHL.
- **Overlap with us:** these apply rational functions to a *state* on
  error-corrected hardware. They do not build subspace or Loewner pencils
  from near-term time-series data.

---

## 5. What could still be new (our claim, to be tested)

1. **Continuous-time Laplace formulation.** We use resolvents of $H$, not of
   $U$, with arbitrary complex poles. The Loewner pencil is exactly rational
   Krylov for $H$, and there is no time grid (so no Nyquist or aliasing
   constraint on $\Delta t$).
2. **Random-time Monte Carlo estimator** with shots shared across all poles,
   plus a variance and shot-count analysis. FD assumes a given signal on a
   grid, and QPD's sampling model needs checking.
3. **NLA theory.** Convergence via rational approximation (pole selection,
   Zolotarev-type rates) and perturbation theory for the noisy, thresholded
   Loewner pencil, analogous to Epperly–Lin–Nakatsukasa. The FD literature is
   largely empirical and signal-processing oriented; QPD has bounds for
   prolate windows, which must be compared.
4. **Decoherence as a pole shift** without noise amplification. Classical FD
   handles decaying signals natively; the quantum-noise framing and its
   comparison with rescaling may be new. QPD's robustness claims must be
   checked.

## 6. Must-read before writing

- Stroschein, Castaldo, Reiher, arXiv:2507.15148 (QPD); method and error
  bounds.
- Mandelshtam, review (2001), and Mandelshtam–Taylor, J. Chem. Phys. 107,
  6756 (1997); Wall–Neuhauser, J. Chem. Phys. 102, 8011 (1995).
- Cortes–Gray, arXiv:2109.06868 (FD connection to quantum Krylov).
- Kyriienko, arXiv:1901.09988 (resolvent via time evolutions).
- Loewner framework (Mayo–Antoulas 2007) and rational Krylov theory (Ruhe;
  Güttel's survey) for the NLA side.

## Sources

- QFD: https://arxiv.org/abs/1909.08925
- FDM review (Mandelshtam 2001): https://www.chem.uci.edu/~mandelsh/publ/review.pdf
- Wall–Neuhauser 1995: https://pubs.aip.org/aip/jcp/article/102/20/8011/293261/Extraction-through-filter-diagonalization-of
- Harmonic inversion comparison: https://arxiv.org/html/nlin/0104064
- QPD: https://arxiv.org/abs/2507.15148
- Cortes–Gray: https://arxiv.org/pdf/2109.06868
- Kyriienko: https://arxiv.org/abs/1901.09988v2
- Takahira et al.: https://arxiv.org/pdf/2106.08076 ; Jiang–An: https://arxiv.org/abs/2601.11959
- Rational transformation via LCHS: https://arxiv.org/html/2408.07742
- Futamura–Ye–Sakurai: https://arxiv.org/abs/2112.05395
