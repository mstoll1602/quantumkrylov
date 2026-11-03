"""Pilot: rational Krylov (Loewner pencil) from quantum time-series data.

Transverse-field Ising chain with N spins (n = 2^N), start vector |+...+>.
Rational Krylov basis w_k = (z_k - H)^{-1} v0, z_k = a_k + i*eta. Its projected
pencil is a Loewner pencil of the scalar r(z) = v0^* (z - H)^{-1} v0:
    G_jk = (conj(r_j) - r_k) / (z_k - conj(z_j)),   M_jk = z_k G_jk - conj(r_j).
r(z) = -i int_0^inf e^{izt} g(t) dt with g(t) = v0^* e^{-itH} v0 is estimated by
Monte Carlo: per shot draw t ~ Exp(eta), run one Hadamard test (Re or Im of g(t))
and reuse every shot for all poles. Compared with real-time Krylov (G, M Toeplitz
from g_k, h_k with h_k measured term by term over the Pauli strings).

Experiments:
  E1  noiseless convergence in m (Lanczos, real-time, Loewner),
  E2  shot noise, equal total number of shots; Loewner with naive poles
      (real parts between -alpha and v0^*Hv0) and with poles around a rough
      energy estimate (the cost of obtaining that estimate is not counted),
  E3  signal damping g(t) -> exp(-gamma t) g(t) (decoherence model), noiseless;
      both methods can undo known damping (shifted poles vs rescaled g_k, h_k).
"""
import numpy as np
import matplotlib.pyplot as plt

rng = np.random.default_rng(0)

N, J, hx = 10, 1.0, 1.0
n = 2**N
I2 = np.eye(2)
X = np.array([[0, 1], [1, 0]], dtype=float)
Z = np.diag([1.0, -1.0])


def op(single, site):
    out = np.array([[1.0]])
    for s in range(N):
        out = np.kron(out, single if s == site else I2)
    return out


paulis, coeffs = [], []
for i in range(N - 1):
    paulis.append(op(Z, i) @ op(Z, i + 1)); coeffs.append(-J)
for i in range(N):
    paulis.append(op(X, i)); coeffs.append(-hx)
coeffs = np.array(coeffs)
H = sum(c * P for c, P in zip(coeffs, paulis))
alpha = np.sum(np.abs(coeffs))

lam, V = np.linalg.eigh(H)
spread = lam[-1] - lam[0]
v0 = np.ones(n) / np.sqrt(n)
w = (V.T @ v0)**2
sup = w > 1e-14                       # eigenpairs seen by v0
lam_s, w_s, V_s = lam[sup], w[sup], V[:, sup]
h0 = np.sum(lam_s * w_s)              # v0^* H v0, cheap to measure
print(f"n={n}  lambda_0={lam[0]:.4f}  |gamma_0|^2={w[0]:.3f}  support={sup.sum()}  "
      f"alpha={alpha}  v0^*Hv0={h0:.2f}")


def g_of_t(t):
    return np.exp(-1j * np.outer(t, lam_s)) @ w_s


def ritz_threshold(M, G, tau):
    """Smallest Ritz value after diagonal scaling and discarding eigenvalues of G below tau."""
    d = np.real(np.diag(G))
    if np.any(d <= 0):
        return np.nan
    D = 1 / np.sqrt(d)
    G, M = D[:, None] * G * D[None, :], D[:, None] * M * D[None, :]
    s, W = np.linalg.eigh((G + G.conj().T) / 2)
    keep = s > tau
    if not keep.any():
        return np.nan
    W = W[:, keep] / np.sqrt(s[keep])
    return np.linalg.eigvalsh(W.conj().T @ ((M + M.conj().T) / 2) @ W)[0]


def toeplitz_herm(sig, m):
    T = np.empty((m, m), dtype=complex)
    for j in range(m):
        for k in range(m):
            T[j, k] = sig[k - j] if k >= j else np.conj(sig[j - k])
    return T


def loewner(z, r):
    G = (np.conj(r)[:, None] - r[None, :]) / (z[None, :] - np.conj(z)[:, None])
    M = z[None, :] * G - np.conj(r)[:, None]
    return M, G


def poles(m, eta):
    """Real parts between the cheap lower bound -alpha and v0^*Hv0; no knowledge of lambda_0."""
    return np.linspace(-alpha, h0, m) + 1j * eta


def poles_near(m, eta, centre):
    """Real parts around a rough estimate of lambda_0."""
    return centre + np.linspace(-1.0, 2.0, m) + 1j * eta


def lanczos_ritz(m):
    Q = np.zeros((n, m)); q = v0.copy()
    for j in range(m):
        Q[:, j] = q
        x = H @ q
        x -= Q[:, :j + 1] @ (Q[:, :j + 1].T @ x)
        x -= Q[:, :j + 1] @ (Q[:, :j + 1].T @ x)
        q = x / np.linalg.norm(x)
    return np.linalg.eigvalsh(Q.T @ H @ Q)[0]


dt = np.pi / spread
err = lambda th: abs(th - lam[0]) / spread

# ---------------- E1: noiseless convergence ----------------
ms = np.arange(2, 25)
E1 = {"Lanczos": [], "real time": [], "Loewner eta=0.5": [], "Loewner eta=1": []}
for m in ms:
    E1["Lanczos"].append(err(lanczos_ritz(m)))
    tk = dt * np.arange(m)
    gk = g_of_t(tk)
    hk = np.exp(-1j * np.outer(tk, lam_s)) @ (lam_s * w_s)
    E1["real time"].append(err(ritz_threshold(toeplitz_herm(hk, m), toeplitz_herm(gk, m), 1e-13)))
    for eta, key in [(0.5, "Loewner eta=0.5"), (1.0, "Loewner eta=1")]:
        z = poles(m, eta)
        r = np.array([np.sum(w_s / (zk - lam_s)) for zk in z])
        E1[key].append(err(ritz_threshold(*loewner(z, r), 1e-13)))
print("E1 (noiseless, m=8/16/24):", {k: [f"{v[i]:.1e}" for i in (6, 14, 22)] for k, v in E1.items()})


# ---------------- shot-noise simulators ----------------
def hadamard_estimate(x, shots):
    """Mean of `shots` +-1 outcomes with expectation x (exact binomial sampling)."""
    return 2 * rng.binomial(shots, (1 + x) / 2) / shots - 1


def realtime_noisy(m, Ntot, gamma=0.0):
    """Real-time Krylov with Ntot shots spread evenly over the 4(m-1)+1 real numbers."""
    Q = 4 * (m - 1) + 1
    Nq = max(Ntot // Q, 1)
    tk = dt * np.arange(m)
    damp = np.exp(-gamma * tk)
    Uv = (V_s * np.exp(-1j * np.outer(tk, lam_s))[:, None, :]) @ (V_s.T @ v0)   # m x n
    gk = (Uv @ v0) * damp
    terms = np.array([[(P @ v0) @ u for P in paulis] for u in Uv]) * damp[:, None]  # m x L
    g_est = np.empty(m, complex); h_est = np.empty(m, complex)
    g_est[0] = 1.0
    for k in range(1, m):
        g_est[k] = hadamard_estimate(gk[k].real, Nq) + 1j * hadamard_estimate(gk[k].imag, Nq)
    Nl = np.maximum((Nq * np.abs(coeffs) / alpha).astype(int), 1)
    for k in range(m):
        re = sum(c * hadamard_estimate(x.real, s) for c, x, s in zip(coeffs, terms[k], Nl))
        im = 0.0 if k == 0 else sum(c * hadamard_estimate(x.imag, s) for c, x, s in zip(coeffs, terms[k], Nl))
        h_est[k] = re + 1j * im
    tau = 10 / np.sqrt(Nq)
    th = ritz_threshold(toeplitz_herm(h_est, m), toeplitz_herm(g_est, m), tau)
    mean_time = np.mean(tk[1:])
    return th, mean_time


def loewner_noisy(z_all, eta, Ntot, gamma=0.0, correct=True, chunk=20000):
    """Monte Carlo estimate of r at every pole set in the dict z_all from one set of Ntot shots.

    All poles must have imaginary part eta (the rate of the exponential time distribution)."""
    a_all = np.concatenate([z.real for z in z_all.values()])
    acc = np.zeros(a_all.size, complex)
    for start in range(0, Ntot, chunk):
        c = min(chunk, Ntot - start)
        t = rng.exponential(1 / eta, c)
        gt = g_of_t(t) * np.exp(-gamma * t)
        use_re = rng.random(c) < 0.5
        x = np.where(use_re, gt.real, gt.imag)
        s = np.where(rng.random(c) < (1 + x) / 2, 1.0, -1.0)
        Y = 2 * s * np.where(use_re, 1.0, 1j)
        acc += np.exp(1j * np.outer(a_all, t)) @ Y
    r_all = -1j * acc / Ntot / eta
    out, pos = {}, 0
    for key, z in z_all.items():
        m = z.size
        r = r_all[pos:pos + m]; pos += m
        z_eff = z + 1j * gamma if correct else z
        M, G = loewner(z_eff, r)
        dmin = np.min(np.real(np.diag(G)))
        tau = 10 * np.sqrt(2) / (eta**2 * np.sqrt(Ntot) * max(dmin, 1e-300))
        out[key] = ritz_threshold(M, G, tau)
    return out, 1 / eta


# ---------------- E2: equal shot budget ----------------
Ns = [10**4, 10**5, 10**6]
m_rt, eta2 = [8, 12, 16], 0.5
centre = lam[0] + 0.3                 # rough estimate, about the accuracy of a cheap real-time run
z_sets = {("Loewner naive poles", 8): poles(8, eta2),
          ("Loewner poles near estimate", 4): poles_near(4, eta2, centre),
          ("Loewner poles near estimate", 8): poles_near(8, eta2, centre)}
nrep = 10
E2 = {}
for Ntot in Ns:
    for m in m_rt:
        vals = [err(realtime_noisy(m, Ntot)[0]) for _ in range(nrep)]
        E2[("real time", m, Ntot)] = np.nanmedian(vals)
    reps = [loewner_noisy(z_sets, eta2, Ntot)[0] for _ in range(nrep)]
    for key in z_sets:
        E2[key + (Ntot,)] = np.nanmedian([err(r[key]) for r in reps])
print("E2 median relative error (equal total shots):")
for (meth, m, Ntot), v in sorted(E2.items(), key=lambda x: (x[0][0], x[0][1], x[0][2])):
    print(f"  {meth:28s} m={m:2d} Ntot={Ntot:.0e}: {v:.2e}")
print(f"  mean evolution time per shot: real time (m=12) {np.mean(dt*np.arange(1, 12)):.2f}, "
      f"Loewner eta={eta2} {1/eta2:.2f}")

# ---------------- E3: signal damping, noiseless ----------------
gammas = np.array([0.0, 0.02, 0.05, 0.1, 0.2, 0.4])
E3 = {"real time m=12": [], "real time m=12, rescaled": [],
      "Loewner m=8, uncorrected": [], "Loewner m=8, corrected": []}
for gam in gammas:
    m = 12
    tk = dt * np.arange(m)
    damp = np.exp(-gam * tk)
    gk = g_of_t(tk) * damp
    hk = (np.exp(-1j * np.outer(tk, lam_s)) @ (lam_s * w_s)) * damp
    E3["real time m=12"].append(err(ritz_threshold(toeplitz_herm(hk, m), toeplitz_herm(gk, m), 1e-13)))
    E3["real time m=12, rescaled"].append(
        err(ritz_threshold(toeplitz_herm(hk / damp, m), toeplitz_herm(gk / damp, m), 1e-13)))
    z = poles(8, 0.5)
    r_damped = np.array([np.sum(w_s / (zk + 1j * gam - lam_s)) for zk in z])   # = exact r(z + i gamma)
    E3["Loewner m=8, uncorrected"].append(err(ritz_threshold(*loewner(z, r_damped), 1e-13)))
    E3["Loewner m=8, corrected"].append(err(ritz_threshold(*loewner(z + 1j * gam, r_damped), 1e-13)))
print("E3 (damping gamma =", gammas, "):")
for k, v in E3.items():
    print(f"  {k:26s}", " ".join(f"{x:.1e}" for x in v))

# ---------------- figure ----------------
plt.rcParams.update({"font.size": 9, "font.family": "serif"})
fig, ax = plt.subplots(1, 3, figsize=(9.6, 2.7))
for (k, v), sty in zip(E1.items(), ["k-", "o-", "s-", "^-"]):
    ax[0].semilogy(ms, v, sty, ms=3, label=k)
ax[0].set_xlabel("subspace dimension $m$"); ax[0].set_ylabel("relative error")
ax[0].set_title("E1: noiseless"); ax[0].legend(fontsize=6.5); ax[0].set_ylim(1e-15, 1)
for m, c in zip(m_rt, ["C0", "C1", "C2"]):
    ax[1].loglog(Ns, [E2[("real time", m, N_)] for N_ in Ns], "o-", color=c, ms=3, label=f"real time $m={m}$")
for key, c, sty in zip(z_sets, ["C3", "C4", "C5"], ["x:", "s--", "s--"]):
    lab = ("Loewner naive" if "naive" in key[0] else "Loewner near est.") + f" $m={key[1]}$"
    ax[1].loglog(Ns, [E2[key + (N_,)] for N_ in Ns], sty, color=c, ms=3, label=lab)
ax[1].set_xlabel("total shots"); ax[1].set_title("E2: equal shot budget"); ax[1].legend(fontsize=6)
for (k, v), sty in zip(E3.items(), ["o-", "d--", "x:", "s-"]):
    ax[2].semilogy(gammas, np.maximum(v, 1e-16), sty, ms=3, label=k)
ax[2].set_xlabel(r"damping rate $\gamma$"); ax[2].set_title("E3: signal damping")
ax[2].legend(fontsize=6.5)
fig.tight_layout()
fig.savefig("fig_rational_pilot.pdf")
