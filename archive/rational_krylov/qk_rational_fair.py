"""Fair comparison: rational Krylov (Loewner) versus real-time Krylov from Hadamard-test data.

Extends qk_rational_pilot.py (same conventions: start vector, relative error
|theta - lambda_0| / (lambda_max - lambda_0), thresholded Rayleigh-Ritz after
diagonal scaling, exact binomial +-1 shot simulation restricted to the
eigenpairs that carry weight in v0).

Methods (both use only Hadamard tests of exp(-itH), optionally followed by a
Pauli string):
  RT   real-time Krylov, v_k = U^k v0, U = exp(-i dt H), dt = pi/(lam_max-lam_min);
       Toeplitz G, M from g_k, h_k; h_k measured term by term (shots ~ |alpha_l|).
  LW   Loewner pencil of r(z) = v0^*(z-H)^{-1}v0 at poles z = a + i*eta, estimated
       by Monte Carlo with one random time per shot from Exp(eta) truncated at
       T_max = KAPPA/eta (density renormalised, so the estimator is unbiased for
       the truncated integral r_T). The same shots serve all poles.
  2S   two-stage: stage 1 = RT with small m1 on a fraction f of the budget gives
       theta_1; stage 2 = LW with poles theta_1 + linspace(-1, 2, m). Both stages
       are charged to the budget. lambda_0 is never used.

Cost: (a) total shots, (b) total evolution time = sum over shots of the time
of that shot's circuit. Decoherence: every Hadamard-test expectation at time t
is damped by exp(-gamma t); RT rescales by exp(+gamma t_k), LW evaluates the
Loewner pencil at z + i*gamma. Thresholds: rule tau = C_RULE*sigma (sigma = a
priori bound on the entrywise noise std of the scaled G) and an oracle
(best multiplier over a grid, chosen per configuration on the median).

Parts: P0 noiseless/truncation checks, P1 the eta=1 plateau diagnosis,
P2 main sweep (models x gamma x budgets, multiprocessing), figures and tables.
Results are printed and pickled to rational_fair_results.pkl (scratch).
"""
import os
import sys
import time
import pickle
import numpy as np

T0 = time.time()
SEED = 20261103
C_GRID = np.array([0.1, 0.3, 1, 3, 10, 30, 100, 300, 1000])
C_RULE = 10.0
I_RULE = int(np.where(C_GRID == C_RULE)[0][0])
KAPPA = 8.0                           # T_max = KAPPA / eta
OFFS_M = (4, 6, 8)                    # Loewner pole counts
OFFS = np.concatenate([np.linspace(-1.0, 2.0, m) for m in OFFS_M])   # pole real parts - centre
M_RT = (4, 6, 8, 12, 16, 20)          # real-time dimensions
M1, F_STAGE1 = 4, (0.1, 0.2, 0.3)     # stage-1 dimension and budget fractions


# ---------------------------------------------------------------- models
class Model:
    """Hamiltonian as Pauli strings; spectral data restricted to the support of v0.

    Each Pauli string P acts as P e_b = phase[b] e_{b xor mask}."""

    def __init__(self, name, N, terms, v0):
        self.name, self.N = name, N
        n = 2**N
        idx = np.arange(n)
        H = np.zeros((n, n))
        for c, mask, ph in terms:
            H[idx ^ mask, idx] += c * ph
        lam, V = np.linalg.eigh(H)
        del H
        cw = V.T @ v0
        sup = cw**2 > 1e-14
        self.lam_all0, self.lam_all1 = lam[0], lam[-1]
        self.lam, self.c = lam[sup], cw[sup]          # v0 = sum_j c_j x_j (+ orthogonal part)
        self.w = self.c**2
        Vs = V[:, sup]
        del V
        self.coef = np.array([t[0] for t in terms])
        self.alpha = np.abs(self.coef).sum()
        # B[l, j] = v0^T P_l x_j, so v0^* P_l exp(-itH) v0 = sum_j B[l,j] c_j exp(-i t lam_j)
        PV = np.empty((len(terms), n))
        for l, (c, mask, ph) in enumerate(terms):
            PV[l, idx ^ mask] = ph * v0
        self.B = PV @ Vs
        assert np.allclose(self.coef @ self.B, self.lam * self.c, atol=1e-10)
        self.lam0 = self.lam_all0
        self.spread = self.lam_all1 - self.lam_all0
        self.dt = np.pi / self.spread
        self.h0 = np.sum(self.lam * self.w)
        self.gap = self.lam[1] - self.lam[0]
        self.S = sup.sum()
        self._grid = None

    def err(self, th):
        return np.abs(np.asarray(th) - self.lam0) / self.spread

    def g(self, t):
        return np.exp(-1j * np.outer(t, self.lam)) @ self.w

    def terms_t(self, t):
        """v0^* P_l exp(-itH) v0 for all l (len(t) x L)."""
        return (np.exp(-1j * np.outer(t, self.lam)) * self.c) @ self.B.T

    def r_exact(self, z, T=np.inf, gamma=0.0):
        """r_T(z) = -i int_0^T exp(izt) exp(-gamma t) g(t) dt (closed form)."""
        z = np.asarray(z)[:, None] + 1j * gamma
        d = z - self.lam[None, :]
        f = 1.0 if np.isinf(T) else 1 - np.exp(1j * d * T)
        return np.sum(self.w * f / d, axis=1)

    # cubic Hermite table of g_c(t) = exp(i*cen*t) g(t) (frequencies in [-Lam, Lam])
    def g_fast(self, t, Tmax):
        if self._grid is None or self._grid[0] < Tmax:
            cen = 0.5 * (self.lam[0] + self.lam[-1])
            Lam = 0.5 * (self.lam[-1] - self.lam[0])
            hgrid = 0.1 / Lam
            tg = np.arange(0, Tmax + 2 * hgrid, hgrid)
            E = np.exp(-1j * np.outer(tg, self.lam - cen))
            self._grid = (Tmax, cen, hgrid, E @ self.w, E @ (-1j * (self.lam - cen) * self.w))
        _, cen, hg, G0, G1 = self._grid
        u = t / hg
        i = np.floor(u).astype(int)
        s = u - i
        s2, s3 = s * s, s * s * s
        gc = ((2 * s3 - 3 * s2 + 1) * G0[i] + (s3 - 2 * s2 + s) * hg * G1[i]
              + (-2 * s3 + 3 * s2) * G0[i + 1] + (s3 - s2) * hg * G1[i + 1])
        return gc * np.exp(-1j * cen * t)


def build_model(kind, N, h, hz=0.0):
    n = 2**N
    idx = np.arange(n)
    zsign = lambda s: 1.0 - 2.0 * ((idx >> (N - 1 - s)) & 1)
    terms = []
    if kind == "ising":                  # -J sum Z Z - h sum X (- hz sum Z), J = 1, start |+...+>
        for i in range(N - 1):
            terms.append((-1.0, 0, zsign(i) * zsign(i + 1)))
        for i in range(N):
            terms.append((-h, 1 << (N - 1 - i), np.ones(n)))
        if hz:
            for i in range(N):
                terms.append((-hz, 0, zsign(i)))
        v0 = np.ones(n) / np.sqrt(n)
        name = f"Ising N={N} h={h:g}" + (f" hz={hz:g}" if hz else "")
    else:                                # Heisenberg XXZ: sum XX + YY + h ZZ, start Neel |0101...>
        for i in range(N - 1):
            mask = (1 << (N - 1 - i)) | (1 << (N - 2 - i))
            zz = zsign(i) * zsign(i + 1)
            terms += [(1.0, mask, np.ones(n)), (1.0, mask, -zz), (h, 0, zz)]
        v0 = np.zeros(n); v0[int("01" * (N // 2), 2)] = 1.0
        name = f"XXZ N={N} Delta={h:g} (Neel)"
    return Model(name, N, terms, v0)


# ---------------------------------------------------------------- linear algebra
def ritz_multi(M, G, taus):
    """Smallest Ritz value of the thresholded, diagonally scaled pencil for each tau."""
    out = np.full(len(taus), np.nan)
    d = np.real(np.diag(G))
    if np.any(d <= 0) or not np.all(np.isfinite(G)):
        return out
    D = 1 / np.sqrt(d)
    G, M = D[:, None] * G * D[None, :], D[:, None] * M * D[None, :]
    s, W = np.linalg.eigh((G + G.conj().T) / 2)
    Mh = (M + M.conj().T) / 2
    last = None
    for i, tau in enumerate(taus):
        keep = s > tau
        if not keep.any():
            continue
        if last is not None and keep.sum() == last[0]:
            out[i] = last[1]; continue
        Wk = W[:, keep] / np.sqrt(s[keep])
        th = np.linalg.eigvalsh(Wk.conj().T @ Mh @ Wk)[0]
        out[i] = th; last = (keep.sum(), th)
    return out


def toeplitz_herm(sig, m):
    j, k = np.meshgrid(np.arange(m), np.arange(m), indexing="ij")
    T = sig[np.abs(k - j)]
    return np.where(k >= j, T, np.conj(T))


def loewner(z, r):
    G = (np.conj(r)[:, None] - r[None, :]) / (z[None, :] - np.conj(z)[:, None])
    M = z[None, :] * G - np.conj(r)[:, None]
    return M, G


def loewner_sigma(z_eff, G, s_r):
    """A priori bound on the entrywise noise std of the diagonally scaled Loewner G."""
    d = np.abs(np.real(np.diag(G)))
    den = np.abs(z_eff[None, :] - np.conj(z_eff)[:, None])
    return np.max(2 * s_r / den / np.sqrt(d[:, None] * d[None, :]))


def errs_loewner(model, z, r, gamma, s_r, correct=True, taus=None):
    z_eff = z + 1j * gamma if correct else z
    M, G = loewner(z_eff, r)
    if taus is None:
        d = np.real(np.diag(G))
        if np.any(d <= 0):
            return np.full(len(C_GRID), np.inf)
        taus = C_GRID * loewner_sigma(z_eff, G, s_r)
    e = model.err(ritz_multi(M, G, taus))
    return np.where(np.isnan(e), np.inf, e)


# ---------------------------------------------------------------- shot simulators
def binom_mean(rng, x, shots):
    """Mean of `shots` +-1 outcomes with expectation x (exact binomial sampling)."""
    return 2 * rng.binomial(shots, (1 + np.clip(x, -1, 1)) / 2) / shots - 1


def rt_alloc(model, m, budget, kind):
    """Shot allocation of RT: Nq shots per g-number, Nl ~ Nq |alpha_l|/alpha per Pauli term.

    Returns Nq, Nl, total shots, total evolution time, max time."""
    tk = model.dt * np.arange(m)
    frac = np.abs(model.coef) / model.alpha
    if kind == "shots":
        Nq = max(int(budget // (4 * m - 3)), 1)
    else:                                # evolution-time budget
        Nq = max(int(budget // (4 * tk[1:].sum())), 1)
    Nl = np.maximum(np.round(Nq * frac).astype(int), 1)
    shots = 2 * (m - 1) * Nq + (2 * m - 1) * Nl.sum()
    ttot = tk[1:].sum() * (2 * Nq + 2 * Nl.sum())
    return Nq, Nl, shots, ttot, tk[-1]


def rt_run(model, m, budget, kind, gamma, rng):
    """Real-time Krylov; returns errors over C_GRID (rescaled, unrescaled), cost, raw Ritz (rule)."""
    Nq, Nl, shots, ttot, tmax = rt_alloc(model, m, budget, kind)
    tk = model.dt * np.arange(m)
    damp = np.exp(-gamma * tk)
    gk = model.g(tk) * damp
    trm = model.terms_t(tk) * damp[:, None]                     # m x L
    g_est = np.empty(m, complex)
    g_est[0] = 1.0
    g_est[1:] = binom_mean(rng, gk[1:].real, Nq) + 1j * binom_mean(rng, gk[1:].imag, Nq)
    re = binom_mean(rng, trm.real, Nl[None, :]) @ model.coef
    im = binom_mean(rng, trm.imag, Nl[None, :]) @ model.coef
    im[0] = 0.0
    h_est = re + 1j * im
    out = []
    for corr in (True, False):
        sc = 1 / damp if corr else np.ones(m)
        G = toeplitz_herm(g_est * sc, m)
        M = toeplitz_herm(h_est * sc, m)
        sigma = np.sqrt(2 / Nq) * sc[-1]
        th = ritz_multi(M, G, C_GRID * sigma)
        e = model.err(th)
        out.append(np.where(np.isnan(e), np.inf, e))
        if corr:
            th_rule = th[I_RULE]
    return out[0], out[1], shots, ttot, tmax, th_rule


def lw_cost(eta, Tmax):
    """Mean evolution time per shot of the truncated exponential distribution."""
    q = np.exp(-eta * Tmax)
    return 1 / eta - Tmax * q / (1 - q)


def lw_mc(model, eta, gamma, groups, rng, chunk=100_000):
    """Monte Carlo estimate of r_T at poles centre + OFFS for several groups.

    groups: list of (centre, sorted list of shot counts). One shot stream (truncated
    Exp(eta) times, one Hadamard test of Re or Im of exp(-gamma t) g(t) per shot) is
    shared; group g uses its first N shots. Returns list of arrays (len(Ns), len(OFFS))."""
    Tmax = KAPPA / eta
    norm = (1 - np.exp(-eta * Tmax)) / eta
    Nmax = max(max(Ns) for _, Ns in groups)
    bounds = sorted(set(range(0, Nmax, chunk)) | {N for _, Ns in groups for N in Ns} | {Nmax})
    acc = [np.zeros(len(OFFS), complex) for _ in groups]
    res = [[] for _ in groups]
    q = 1 - np.exp(-eta * Tmax)
    for s0, s1 in zip(bounds[:-1], bounds[1:]):
        c = s1 - s0
        t = -np.log1p(-q * rng.random(c)) / eta
        gt = model.g_fast(t, Tmax) * np.exp(-gamma * t)
        use_re = rng.random(c) < 0.5
        x = np.where(use_re, gt.real, gt.imag)
        Y = np.where(rng.random(c) < (1 + x) / 2, 2.0, -2.0) * np.where(use_re, 1.0, 1j)
        Eoff = np.exp(1j * np.outer(OFFS, t))
        for gi, (cen, Ns) in enumerate(groups):
            if s0 >= Ns[-1]:
                continue
            acc[gi] += Eoff @ (Y * np.exp(1j * cen * t))
            if s1 in Ns:
                res[gi].append(-1j * norm * acc[gi] / s1)
    return [np.array(r) for r in res], 2 * norm / np.sqrt(1.0)   # s_r = 2*norm/sqrt(N)


def split_offs(r):
    """Split an array over OFFS into the pole sets m in OFFS_M."""
    out, p = {}, 0
    for m in OFFS_M:
        out[m] = r[..., p:p + m]; p += m
    return out
