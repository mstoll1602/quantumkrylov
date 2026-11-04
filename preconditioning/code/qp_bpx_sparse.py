"""kappa(C A_beta) for the weighted BPX factor on larger grids, d = 2, 3 (sparse).

A_beta = M + sqrt(beta) K (Q1, uniform grid, h = 2^-L), C = R R^T with
R = [omega_l P_l]_l and omega_l^2 = theta_l / delta_l, theta_l = min(1, sqrt(beta)/h_l^2)^s
(s = 0 is plain diagonally scaled BPX). The extreme eigenvalues of C A_beta are
those of the symmetric-definite pencil (C, A_beta^{-1}); we compute them by
Lanczos, applying A_beta^{-1} by CG preconditioned with C.
"""
import sys
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import LinearOperator, cg, eigsh


def fe1d(l):
    n, h = 2**l - 1, 2.0**-l
    e = np.ones(n)
    M = sp.diags([e[:-1], 4 * e, e[:-1]], [-1, 0, 1]) * (h / 6)
    K = sp.diags([-e[:-1], 2 * e, -e[:-1]], [-1, 0, 1]) / h
    return M.tocsr(), K.tocsr()


def kron_all(mats):
    out = mats[0]
    for A in mats[1:]:
        out = sp.kron(out, A, format="csr")
    return out


def fe(l, d):
    M1, K1 = fe1d(l)
    M = kron_all([M1] * d)
    K = sum(kron_all([K1 if j == i else M1 for j in range(d)]) for i in range(d))
    return M, K


def prolong1d(l):
    nc, nf = 2**l - 1, 2**(l + 1) - 1
    j = np.arange(nc)
    rows = np.r_[2 * j + 1, 2 * j, 2 * j + 2]
    vals = np.r_[np.ones(nc), 0.5 * np.ones(2 * nc)]
    return sp.csr_matrix((vals, (rows, np.r_[j, j, j])), shape=(nf, nc))


def prolong(l, L, d):
    P = sp.identity(2**l - 1, format="csr")
    for k in range(l, L):
        P = prolong1d(k) @ P
    return kron_all([P] * d)


def kappa(L, d, beta, s):
    sb = np.sqrt(beta)
    M, K = fe(L, d)
    A = (M + sb * K).tocsr()
    blocks = []
    for l in range(1, L + 1):
        Ml, Kl = fe(l, d)
        delta = Ml[0, 0] + sb * Kl[0, 0]
        theta = min(1.0, sb * 4.0**l) ** s
        blocks.append(prolong(l, L, d) * np.sqrt(theta / delta))
    R = sp.hstack(blocks, format="csr")
    RT = R.T.tocsr()
    n = A.shape[0]
    C = LinearOperator((n, n), matvec=lambda x: R @ (RT @ x), dtype=float)

    def ainv(x):
        y, info = cg(A, x, M=C, rtol=1e-12, maxiter=2000)
        assert info == 0
        return y

    Ainv = LinearOperator((n, n), matvec=ainv, dtype=float)
    Aop = LinearOperator((n, n), matvec=lambda x: A @ x, dtype=float)
    lmax = eigsh(C, k=1, M=Ainv, Minv=Aop, which="LA", tol=1e-6, return_eigenvectors=False)[0]
    lmin = eigsh(C, k=1, M=Ainv, Minv=Aop, which="SA", tol=1e-6, return_eigenvectors=False)[0]
    return lmax / lmin


if __name__ == "__main__":
    d = int(sys.argv[1])
    Ls = [int(a) for a in sys.argv[2:]]
    for L in Ls:
        for beta in [1.0, 1e-2, 1e-4, 1e-6, 1e-8]:
            ks = [kappa(L, d, beta, s) for s in (0, 1)]
            print(f"d={d} L={L} n={(2**L - 1)**d:7d} beta={beta:7.0e}  "
                  f"plain: {ks[0]:6.2f}  weighted(s=1): {ks[1]:6.2f}", flush=True)
