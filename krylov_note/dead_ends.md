# Ideas that did not work (or only partly)

Notes on approaches we examined for quantum Krylov methods that fail, or work
only in restricted settings. The aim is to avoid revisiting them and to keep the
reasons on record.

**Setting.** The device prepares Krylov vectors $v_k$, e.g. $v_k = U^k v_0$ with
$U = e^{-i\Delta t H}$. It returns noisy estimates of inner products only. Each
number is accurate to an additive $\varepsilon \sim N_s^{-1/2}$, where $N_s$ is
the number of shots; in practice $\varepsilon \approx 10^{-3}$. The Gram matrix
$G = V^*V$ has $\kappa(G) = \kappa(V)^2$, which grows exponentially with $m$.
The core difficulty is the noisy, ill-conditioned pencil $(M, G)$
(note, Sections 3 and 5).

---

## 1. Orthogonalizing the Krylov vectors on the device

**Idea.** Run Gram–Schmidt or the Lanczos recurrence on the device, so that
$G = I$.

**Why it fails.**
- The coefficients $q_j^* v_k$ exist only as noisy numbers obtained by
  measurement. They are never available inside the device.
- The map $v_k \mapsto q_k$ depends on $v_k$ through the coefficients and the
  normalization. It is not multiplication by a fixed unitary matrix, which is
  all the device can do.
- The device holds one vector, and vectors cannot be copied. So
  $q_0,\dots,q_{k-1}$ cannot be kept to subtract them.

**What remains.** Orthogonalization in coefficient space: $G = R^*R$ and
$Q = VR^{-1}$. This is CholeskyQR. It squares the condition number and fails
once $\kappa(V)^2 \gtrsim 1/\varepsilon$, which is far earlier than the
classical $1/\mathbf{u}$.

**Status.** Dead in this form. Two variants prepare orthonormal Lanczos vectors
directly: block-encoded Lanczos polynomials, and the operator-level recursion
of arXiv:2607.09476 (not yet read). Both still need recurrence coefficients
computed from noisy data. Their loss of orthogonality, analysed in Paige's
style with $\varepsilon$ in place of $\mathbf{u}$, is an open question.

## 2. Oblivious sketching of $\mathbb{C}^n$ to avoid squaring $\kappa$

**Idea.** Sketched Rayleigh–Ritz (Nakatsukasa–Tropp): work with $SV$ instead of
$G$, where $S \in \mathbb{C}^{s\times n}$ is random and $s \approx 2m$.
Householder QR of $SV$ has $\kappa(SV) \approx \kappa(V)$, so nothing is
squared. On the device, an entry $s_i^* v_k$ is an inner product that a
Hadamard test can estimate.

**Why it fails.**
- **Exponentially small entries.** For a random unit vector,
  $|s_i^* v_k|^2 \approx 1/n = 2^{-N}$, but the noise is additive, of size
  $\varepsilon$. The rescaling by $\sqrt{n/s}$ that a subspace embedding needs
  amplifies the noise by about $2^{N/2}$.
- **No preparation.** Random vectors in $\mathbb{C}^n$ have no short circuit.
- **Not fixable by a better sketch.** Any sketch chosen without knowledge of
  $V$ must work for every subspace. By concentration of measure, typical
  subspaces have overlaps of size $n^{-1/2}$ with any fixed probe vectors.

**Status.** Dead for oblivious sketches.

## 3. Hadamard probe vectors and the SRHT

**Idea.** Use rows of the Walsh–Hadamard matrix
$F = h\otimes\cdots\otimes h$ as sketch vectors, as in probing for diagonal
estimation.

**Result.**
- It removes the preparation problem: $F e_b$ is a product state that needs
  $N$ gates.
- It does not remove the magnitude problem. The rows form an orthonormal basis,
  so $\sum_b |(Fe_b)^* v|^2 = 1$ and the typical overlap is again $1/n$.
- A Hadamard sketch is just a change of basis. It helps only if the vectors are
  concentrated on few rows of $F$. In that case, measuring $Fv$ samples those
  rows directly, and no sketch is needed.
- The subsampled randomized Hadamard transform (SRHT) is designed to *spread*
  vectors evenly over all rows, which is exactly the wrong property here.

**Evidence.** Ising chain, $N = 10$, $n = 1024$. The table gives the number of
basis vectors that carry 99% of the squared norm (script `code/qk_concentration.py`):

| field $h$ | vector | standard basis | Hadamard basis |
|---|---|---|---|
| 0.5 | ground state | 83 | 467 |
| 1.0 (critical) | ground state | 793 | 206 |
| 1.0 | $v_{10}$ | 664 | 348 |
| 2.0 | ground state | 971 | 16 |
| 2.0 | $v_{10}$ | 976 | 45 |

**Contrast.** Trace and diagonal estimation *do* work with random or Hadamard
probes on a device (e.g. the density of states in Summer et al.,
arXiv:2303.13476). There the target $\operatorname{tr}(A)/n$ is $O(1)$. A
Krylov sketch needs $O(1)$ overlaps with one particular subspace.

**Status.** Dead as a fix for conditioning.

## 4. Sketches that depend on $V$ (cross-Gram matrices)

**Idea.** Take sketch vectors correlated with $v_0$, e.g.
$s_i = W_i v_0$ with cheap unitaries $W_i$. Their overlaps are $O(1)$.

**Result.**
- $SV = W^*V$ is a cross-Gram matrix, and $(W^*HV, W^*V)$ is a Petrov–Galerkin
  (oblique) projection.
- In the ideal case, $W$ is orthonormal with
  $\operatorname{range}(V) \subseteq \operatorname{range}(W)$. Then
  $\kappa(W^*V) = \kappa(V)$, with no squaring.
- But preparing such a $W$ is the orthogonalization problem again (item 1).

**Status.** Not dead, but no construction that helps is known. Open.

## 5. Monte Carlo over the Pauli terms of $h_k$

**Idea.** Estimate $h_k = \sum_\ell \alpha_\ell\, v_0^* P_\ell U^k v_0$ by
drawing $\ell$ with probability $|\alpha_\ell|/\alpha$.

**Result.** The estimator is unbiased and needs
$O(\alpha^2/\varepsilon^2)$ shots, with $\alpha = \sum_\ell |\alpha_\ell|$.
This matches an optimal fixed allocation of shots over the $L$ terms, so there
is **no gain in shot count**. Sampling over terms pays off only where the terms
cannot be enumerated or must be compiled:
- qDRIFT for the time evolution itself;
- continuous time integrals (Lin–Tong).

**Status.** Correct but not an improvement for $h_k$.

## 6. Sample-based subspaces (row sampling by measurement)

**Idea.** Measurement samples index $b$ with probability $|v_b|^2$, i.e.
squared-row-norm or leverage sampling. Collect the observed indices $B$ and
solve for $H(B,B)$ classically (sample-based diagonalization and its Krylov
variant, arXiv:2501.09702).

**Result.**
- No conditioning problem, since $G = I$.
- It is variational: $\lambda_{\min}(H(B,B)) \ge \lambda_0$ by Cauchy
  interlacing.
- It needs the ground state to be concentrated in *some* basis, and which basis
  depends on the problem (table in item 3). Near criticality neither the
  standard nor the Hadamard basis works well. In general $|B|$ grows
  exponentially with $N$.

**Status.** Works in a restricted regime. Not a general fix.

## 7. Monte Carlo for the eigenvalue itself

**Idea.** Sample over time indices or terms to estimate the Ritz value directly.

**Why it fails.** Rayleigh–Ritz is nonlinear in the data. Monte Carlo over
terms applies only to linear functionals, so here it can be used only entry by
entry. Methods built on linear functionals of the time signal, e.g. Lin–Tong's
cumulative spectral function, can sample one random time per shot, but they are
a different algorithm.

**Status.** Dead for Rayleigh–Ritz.

---

## Lessons for the next direction

1. Anything we measure must be an **$O(1)$ quantity**, read off with additive
   noise. Exponentially small overlaps are hopeless.
2. **Structure in coefficient space** is the main lever: Toeplitz $G$,
   Szegő/Levinson recursions, and regularization that preserves structure.
3. **Linear functionals of the time signal** $g(t) = v_0^* e^{-itH} v_0$ are
   cheap. They allow Monte Carlo over $t$, and their noise is controlled by the
   $\ell_1$ norm of the kernel.
4. Point 3 is what makes **resolvent-based (rational) methods** attractive:
   $v_0^*(z-H)^{-1}v_0 = -i\int_0^\infty e^{izt} g(t)\,dt$ for $\operatorname{Im} z > 0$
   is exactly such a linear functional.
