"""Numerical checks for the time-periodic heat control section.

Control of y_t - Laplace(y) = u on (0,1)^d x (0,T), y = 0 on the boundary,
y(0) = y(T). Q1 in space (h = 2^-L), periodic backward Euler in time
(N_t = 2^Lt steps, tau = T/N_t). Scaled reduced KKT matrix
    A = [[I(x)M, sqrt(beta) Kt^T], [sqrt(beta) Kt, -I(x)M]],  Kt = D_t (x) M + I (x) K,
with D_t = (I - S)/tau circulant. The DFT in time diagonalises D_t with eigenvalues
lam_k = (1 - exp(-2 pi i k/N_t))/tau, and A splits into Hermitian blocks
    A_k = [[M, sqrt(beta)(K + conj(lam_k) M)], [sqrt(beta)(K + lam_k M), -M]].
We check
  1. kappa(A) and the bound |eig(P^{-1} A)| in [1/sqrt3, 1] for the global
     preconditioner P = blkdiag(Pt, Pt), Pt = I(x)(M + sqrt(beta)K) + sqrt(beta)|D_t|(x)M,
     i.e. per frequency A_k := r_k M + sqrt(beta) K, r_k = 1 + sqrt(beta)|lam_k|,
  2. kappa(C_k A_k) for weighted BPX per frequency (plain: s = 0, weighted: s = 1),
  3. B_k = Rt_k^T A_k Rt_k = Gt_k^T J_k Gt_k with ||J_k|| <= 1, and kappa of the
     nonzero spectrum of B = blkdiag_k(B_k).
"""
import numpy as np
from qp_bpx_sparse import fe, prolong
from qp_poisson_check import matrices


def time_ops(Nt, T=1.0):
    tau = T / Nt
    S = np.roll(np.eye(Nt), 1, axis=0)          # cyclic shift: (S y)_k = y_{k-1}
    Dt = (np.eye(Nt) - S) / tau
    lam = (1 - np.exp(-2j * np.pi * np.arange(Nt) / Nt)) / tau
    w, V = np.linalg.eigh(Dt.T @ Dt)
    absDt = V @ np.diag(np.sqrt(np.maximum(w, 0))) @ V.T
    return Dt, absDt, lam


def pencil_abs_eigs(A, P):
    Li = np.linalg.inv(np.linalg.cholesky(P))
    return np.abs(np.linalg.eigvalsh(Li @ A @ Li.T))


def bpx_k(L, d, sb, r, s, normalise=False):
    """Weighted BPX factor for A_k = r M + sqrt(beta) K. With normalise=True the
    level weights are divided by the finest one (theta_L = 1), so that all
    frequencies share the scale of a Jacobi-scaled finest level."""
    theta = lambda l: min(1.0, sb / r * 4.0**l) ** s
    blocks = []
    for l in range(1, L + 1):
        Ml, Kl = fe(l, d)
        delta = r * Ml[0, 0] + sb * Kl[0, 0]
        t = theta(l) / theta(L) if normalise else theta(l)
        blocks.append(prolong(l, L, d).toarray() * np.sqrt(t / delta))
    return np.hstack(blocks)


def global_check(L, d, Lt, beta):
    sb = np.sqrt(beta)
    M, K = (X.toarray() for X in fe(L, d))
    Nt = 2**Lt
    Dt, absDt, _ = time_ops(Nt)
    It = np.eye(Nt)
    MM = np.kron(It, M)
    Kt = np.kron(Dt, M) + np.kron(It, K)
    A = np.block([[MM, sb * Kt.T], [sb * Kt, -MM]])
    Pt = np.kron(It, M + sb * K) + sb * np.kron(absDt, M)
    ev = np.abs(np.linalg.eigvalsh(A))
    e = pencil_abs_eigs(A, np.kron(np.eye(2), Pt))
    return ev.max() / ev.min(), e.min(), e.max()


def frequency_check(L, d, Lt, beta, s, normalise=False):
    sb = np.sqrt(beta)
    M, K = (X.toarray() for X in fe(L, d))
    F0, Fg = matrices(L, d)
    _, _, lam = time_ops(2**Lt)
    kapC, eig_all, err_fac, normJ = [], [], 0.0, 0.0
    # lam_{Nt-k} = conj(lam_k): those blocks are complex conjugates with the same spectrum
    for lk in lam[:2**(Lt - 1) + 1]:
        r = 1 + sb * abs(lk)
        Ak = r * M + sb * K
        R = bpx_k(L, d, sb, r, s, normalise)
        e = np.linalg.eigvalsh(R.T @ Ak @ R)
        e = e[e > 1e-10 * e.max()]
        kapC.append(e.max() / e.min())
        # B_k and its factorisation
        B0 = R.T @ M @ R
        B1 = sb * R.T @ (K + lk * M) @ R
        Bk = np.block([[B0, B1.conj().T], [B1, -B0]])
        Gv = np.sqrt(r) * F0 @ R
        Gg = np.vstack([beta**0.25 * F @ R for F in Fg])
        Jv = np.array([[1, sb * np.conj(lk)], [sb * lk, -1]]) / r
        normJ = max(normJ, np.linalg.norm(Jv, 2))
        GJG = np.block([[Jv[0, 0] * Gv.T @ Gv, Jv[0, 1] * Gv.T @ Gv + Gg.T @ Gg],
                        [Jv[1, 0] * Gv.T @ Gv + Gg.T @ Gg, Jv[1, 1] * Gv.T @ Gv]])
        err_fac = max(err_fac, np.abs(GJG - Bk).max() / np.abs(Bk).max())
        eb = np.abs(np.linalg.eigvalsh(Bk))
        eig_all.append(eb[eb > 1e-10 * eb.max()])
    eig_all = np.concatenate(eig_all)
    return max(kapC), eig_all.max() / eig_all.min(), err_fac, normJ


if __name__ == "__main__":
    import sys
    parts = sys.argv[1:] or ["global", "freq"]
    if "global" in parts:
        print("1. global Schoeberl-Zulehner-type preconditioner (kappa bound sqrt3 = 1.732)")
        for d, L, Lt in [(1, 5, 5), (2, 3, 4)]:
            for beta in [1.0, 1e-2, 1e-4, 1e-6]:
                kA, emin, emax = global_check(L, d, Lt, beta)
                print(f"d={d} L={L} Nt={2**Lt:3d} beta={beta:7.0e}  kappa(A)={kA:9.2e}  "
                      f"|eig(P^-1 A)| in [{emin:.3f}, {emax:.3f}]  kappa={emax / emin:.3f}", flush=True)
    if "freq" in parts:
        print("\n2./3. per frequency: max_k kappa(C_k A_k) and kappa(B) over all k for plain,\n"
              "weighted and normalised-weighted BPX; factorisation error, max ||J_k||")
        for d, L, Lt in [(1, 6, 6), (1, 8, 8), (2, 4, 5)]:
            for beta in [1.0, 1e-2, 1e-4, 1e-6]:
                out = [frequency_check(L, d, Lt, beta, s, nm)
                       for s, nm in [(0, False), (1, False), (1, True)]]
                print(f"d={d} L={L} Nt={2**Lt:4d} beta={beta:7.0e}  "
                      f"plain: kC={out[0][0]:5.2f} kB={out[0][1]:5.2f}  "
                      f"weighted: kC={out[1][0]:5.2f} kB={out[1][1]:5.2f}  "
                      f"normalised: kC={out[2][0]:5.2f} kB={out[2][1]:5.2f}  "
                      f"fac err={max(o[2] for o in out):.1e}  max||J||={out[1][3]:.3f}", flush=True)
