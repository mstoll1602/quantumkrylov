"""Small illustration for the note on quantum Krylov subspace methods.

Transverse-field Ising chain with N spins (n = 2^N), start vector |+...+>.
Compares
  * Lanczos (full reorthogonalisation) on K_m(H, v0),
  * real-time Krylov:  v_k = exp(-i k dt H) v0,  Gram/projected matrices are
    Toeplitz and built from the 2m-1 numbers g_k = v0^* U^k v0, h_k = v0^* H U^k v0,
  * the same with additive noise on g_k, h_k (mimicking shot noise), with and
    without thresholding of the Gram matrix.
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


H = np.zeros((n, n))
for i in range(N - 1):
    H -= J * op(Z, i) @ op(Z, i + 1)
for i in range(N):
    H -= hx * op(X, i)

lam, V = np.linalg.eigh(H)
spread = lam[-1] - lam[0]
v0 = np.ones(n) / np.sqrt(n)          # |+...+>
c = V.T @ v0                          # coefficients in the eigenbasis
print("n =", n, " lambda_0 =", lam[0], " |u0^* v0|^2 =", abs(c[0])**2)

mmax = 30
dt = np.pi / spread

# signals g_k = v0^* U^k v0 and h_k = v0^* H U^k v0, k = 0..mmax-1
ks = np.arange(mmax)
phase = np.exp(-1j * dt * np.outer(ks, lam))
g = phase @ (abs(c)**2)
h = phase @ (lam * abs(c)**2)


def toeplitz_herm(sig, m):
    """T[j,k] = sig[k-j] for k>=j and conj(sig[j-k]) otherwise."""
    T = np.empty((m, m), dtype=complex)
    for j in range(m):
        for k in range(m):
            T[j, k] = sig[k - j] if k >= j else np.conj(sig[j - k])
    return T


def ritz_threshold(M, G, tau):
    """Smallest eigenvalue of the pencil (M, G) after discarding eigenvalues of G below tau."""
    s, W = np.linalg.eigh(G)
    keep = s > tau
    if not keep.any():
        return np.nan
    W = W[:, keep] / np.sqrt(s[keep])
    return np.linalg.eigvalsh(W.conj().T @ M @ W)[0]


def ritz_plain(M, G):
    """Smallest real part among the eigenvalues of G^{-1} M (no regularisation)."""
    ev = np.linalg.eigvals(np.linalg.solve(G, M))
    return np.min(ev.real)


def lanczos_ritz(m):
    Q = np.zeros((n, m))
    q = v0.copy()
    for j in range(m):
        Q[:, j] = q
        w = H @ q
        w -= Q[:, :j + 1] @ (Q[:, :j + 1].T @ w)
        w -= Q[:, :j + 1] @ (Q[:, :j + 1].T @ w)
        q = w / np.linalg.norm(w)
    return np.linalg.eigvalsh(Q.T @ H @ Q)[0]


ms = np.arange(1, mmax + 1)
err = {k: [] for k in ["lanczos", "rt", "rt6", "rt3", "rt3plain"]}
condG = []
nrep = 50
for m in ms:
    G, M = toeplitz_herm(g, m), toeplitz_herm(h, m)
    condG.append(np.linalg.cond(G))
    err["lanczos"].append(lanczos_ritz(m) - lam[0])
    err["rt"].append(ritz_threshold(M, G, 1e-14) - lam[0])
    for key, eps, tau in [("rt6", 1e-6, 1e-5), ("rt3", 1e-3, 1e-2), ("rt3plain", 1e-3, None)]:
        vals = []
        for _ in range(nrep):
            dg = eps * (rng.standard_normal(mmax) + 1j * rng.standard_normal(mmax)) / np.sqrt(2)
            dh = eps * spread * (rng.standard_normal(mmax) + 1j * rng.standard_normal(mmax)) / np.sqrt(2)
            dg[0] = dg[0].real
            dh[0] = dh[0].real
            Gn, Mn = toeplitz_herm(g + dg, m), toeplitz_herm(h + dh, m)
            th = ritz_plain(Mn, Gn) if tau is None else ritz_threshold(Mn, Gn, tau)
            vals.append(abs(th - lam[0]))
        err[key].append(np.median(vals))

for key in err:
    err[key] = np.abs(np.array(err[key])) / spread

plt.rcParams.update({"font.size": 9, "font.family": "serif"})
fig, ax = plt.subplots(1, 2, figsize=(6.4, 2.5))
ax[0].semilogy(ms, err["lanczos"], "k-", label="Lanczos")
ax[0].semilogy(ms, err["rt"], "o-", ms=3, color="C0", label=r"real time, exact")
ax[0].semilogy(ms, err["rt6"], "s-", ms=3, color="C1", label=r"real time, $\varepsilon=10^{-6}$")
ax[0].semilogy(ms, err["rt3"], "^-", ms=3, color="C2", label=r"real time, $\varepsilon=10^{-3}$")
ax[0].semilogy(ms, err["rt3plain"], "x:", ms=3, color="C3", label=r"$\varepsilon=10^{-3}$, no threshold")
ax[0].set_xlabel("subspace dimension $m$")
ax[0].set_ylabel(r"$|\theta_{\min}-\lambda_0|/(\lambda_{\max}-\lambda_0)$")
ax[0].set_ylim(1e-14, 1e2)
ax[0].legend(fontsize=6.5, loc="lower left")
ax[1].semilogy(ms, condG, "o-", ms=3, color="C0")
ax[1].set_xlabel("subspace dimension $m$")
ax[1].set_ylabel(r"$\kappa_2(G)$, real-time basis")
fig.tight_layout()
fig.savefig("fig_qk.pdf")
print("final errors:", {k: v[-1] for k, v in err.items()})
print("cond(G) at m=10,20,30:", condG[9], condG[19], condG[29])
