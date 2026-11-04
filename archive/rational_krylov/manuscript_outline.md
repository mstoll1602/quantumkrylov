# Manuscript outline: rational Krylov methods from quantum time-series data

*Working titles:* "Rational Krylov subspace methods for quantum eigenvalue
estimation" / "A Loewner framework for quantum Krylov diagonalization"

## 1. Core idea

The rational Krylov basis is $w_k = (z_k - H)^{-1} v_0$, with poles
$z_k = a_k + i\eta_k$ and $\eta_k > 0$. The projected pencil depends only on the
scalar function $r(z) = v_0^*(z-H)^{-1}v_0$, and forms a **Loewner pencil**:
$$
G_{jk} = \frac{\overline{r(z_j)} - r(z_k)}{z_k - \bar z_j}, \qquad
M_{jk} = z_k G_{jk} - \overline{r(z_j)} .
$$
Both matrices are Hermitian by construction, also for noisy $r$. Verified
numerically to $10^{-15}$.

$r$ is a **linear functional of the signal** $g(t) = v_0^* e^{-itH} v_0$ that
real-time Krylov already measures:
$$
r(z) = -i\int_0^\infty e^{izt} g(t)\,dt
     = \frac{-i}{\eta}\,\mathbb{E}_{t\sim \mathrm{Exp}(\eta)}\big[e^{iat} g(t)\big],
\qquad z = a + i\eta .
$$
One random time per shot and one Hadamard test give an unbiased estimate. The
same shots serve **all** poles with the same $\eta$.

## 2. Claimed advantages (to be established)

1. **No new circuits**: only Hadamard tests of $e^{-itH}$, the same as real-time
   Krylov.
2. **No measurements of $H$**: $M$ follows from $r$, so the $L$ Pauli-term
   measurements of $h_k$ are not needed.
3. **Approximation power**: rational filters with poles near the target beat
   polynomial and trigonometric filters (Zolotarev-type rates).
4. **Damping**: if $g(t) \to e^{-\gamma t} g(t)$, the data are exactly
   $r(z + i\gamma)$. Correcting for this is a pole shift, with no
   amplification of noise. Real-time Krylov must instead rescale $g_k$ by
   $e^{\gamma k\Delta t}$, which amplifies noise.
5. **Classical toolbox**: Loewner and model-reduction methods (Mayo–Antoulas,
   SVD truncation, AAA), rational Krylov theory (Ruhe, Güttel), and pole
   selection.

## 3. Pilot results (`qk_rational_pilot.py`, `fig_rational_pilot.pdf`)

Setup: Ising chain, $N = 10$, $h = 1$, start vector $|{+}\dots{+}\rangle$
with $|\gamma_0|^2 = 0.47$. Errors are relative to $\lambda_{\max} - \lambda_0$.

- **E1, noiseless.** Loewner with poles in the lower part of the spectrum,
  $\eta = 0.5$, reaches $8\cdot10^{-12}$ at $m = 8$. Real-time Krylov reaches
  $5\cdot10^{-4}$ and Lanczos $3\cdot10^{-4}$. Caveat: exact $r$ is infinite-time
  data.
- **E2, equal shot budget** (median of 10 runs):

  | method | $10^4$ shots | $10^5$ | $10^6$ |
  |---|---|---|---|
  | real time, best of $m = 8, 12, 16$ | $1.5\cdot10^{-2}$ | $7.5\cdot10^{-3}$ | $3.3\cdot10^{-3}$ |
  | Loewner, naive poles ($a_k \in [-\alpha, v_0^*Hv_0]$) | fails | fails | $2.2\cdot10^{-2}$ |
  | Loewner, poles near a rough estimate, $m = 8$ | $4.0\cdot10^{-3}$ | $4.5\cdot10^{-3}$ | $7.0\cdot10^{-4}$ |

  **Pole placement decides everything.** With poles near a rough estimate
  ($\lambda_0 + 0.3$), Loewner is about 2–5× more accurate. Caveats:
  - The shots needed to obtain the estimate are not counted.
  - The mean evolution time per shot is 2.0 for Loewner versus 0.76 for real
    time, so at equal *total evolution time* the advantage shrinks to
    roughly 1.5–3×.
  - The exponential time distribution has an unbounded tail; truncate it.
- **E3, damping, noiseless.** Both methods recover the undamped result exactly
  when $\gamma$ is known (pole shift versus rescaling). Without correction,
  Loewner fails ($2.5\cdot10^{-1}$) and real time degrades ($2\cdot10^{-3}$
  to $10^{-2}$). **The decisive experiment is still missing**: damping plus shot
  noise, where rescaling amplifies noise and the pole shift broadens the
  filter.
- **Unexplained.** With $\eta = 1$ and poles near the estimate, the error
  stalls at about $5.5\cdot10^{-3}$ from $10^5$ to $10^6$ shots. It may be a
  bias, or the threshold rule.

## 4. Planned contributions

**Theory.**
1. Loewner identity and the Monte Carlo estimator, with variance at most
   $\sim 1/(\eta^2 N_s)$ per pole and correlations across poles.
2. Convergence via rational approximation on the spectral interval; dependence
   on pole distance and gap.
3. Perturbation theory for the noisy Loewner pencil under thresholding,
   analogous to Epperly–Lin–Nakatsukasa. Is the energy error linear in the
   noise?
4. Pole selection: a two-stage scheme (cheap estimate, then poles), and
   optimal $\eta$ as a trade-off between filter sharpness and variance.
5. Damping as a pole shift, and the comparison of noise amplification.

**Algorithms.** A two-stage method, shared-shot estimation, a truncated time
distribution, SVD-based Loewner regularization, and adaptive (greedy, AAA-like)
pole choice.

**Numerics.**
- Ising across phases ($h = 0.5, 1, 2$) and $N$ up to about 14–16 (statevector).
- A Heisenberg or Hubbard model and a small molecule.
- Comparisons at equal shots **and** equal total evolution time.
- Damping plus shot noise.

## 5. Risks and prior art to check first

- **Quantum filter diagonalization** (Parrish–McMahon, arXiv:1909.08925, and
  the PRX Quantum follow-up): filtered states from time evolution, mention of
  Lorentzians. How close is it to resolvent filters?
- **Classical filter diagonalization and harmonic inversion** (Neuhauser;
  Mandelshtam–Taylor): their known equivalence with Padé approximation.
  Position the Loewner view against it.
- Green's functions and resolvents estimated from quantum time series.
- **Physics limit.** Resolving a gap $\delta$ still needs evolution times of
  about $1/\delta$. Gains are in constants and robustness, not in scaling.
- **Conditioning stays bad**: $\kappa(G) \approx 5\cdot10^{17}$ at $m = 8$
  with good poles. Is the noise behaviour better than for real time
  nonetheless?

## 6. Proposed structure

1. Introduction
2. Quantum Krylov background (short; reuse the note)
3. Rational Krylov and the Loewner pencil from $r(z)$
4. Estimating $r(z)$ from time-series data (Monte Carlo, variance,
   truncation, damping)
5. Convergence (rational approximation, pole selection)
6. Stability under noise (thresholded Loewner pencil)
7. Numerical experiments
8. Discussion

*Target venue:* SIAM J. Sci. Comput. or SIMAX (NLA framing), alternatively
Quantum.

## 7. Next steps

1. Read QFD and the filter-diagonalization/Padé literature, and decide how
   novel the work is.
2. Run the damping-plus-noise experiment; resolve the $\eta = 1$ plateau.
3. Do cost accounting at equal total evolution time, including the stage that
   produces the estimate.
4. Repeat on a second model with a smaller gap.
