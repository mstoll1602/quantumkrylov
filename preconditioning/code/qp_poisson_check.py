"""Numerical checks for the draft section on quantum preconditioning (Poisson control).

Distributed control of -Laplace(y) = u on (0,1)^d, homogeneous Dirichlet data,
Q1 elements on a uniform grid with h = 2^-L. Scaled reduced KKT matrix
    As = [[M, sqrt(beta) K], [sqrt(beta) K, -M]],  A_beta = M + sqrt(beta) K.
We check
  1. kappa(As) without preconditioning,
  2. the Schoeberl-Zulehner bound: eig(P^{-1} As) in [-1,-1/sqrt2] u [1/sqrt2,1]
     for P = blkdiag(A_beta, A_beta),
  3. kappa(C A_beta) for diagonally scaled BPX, C = R R^T, R = [P_l D_l^{-1/2}]_l,
  4. kappa_eff of the two-sided system B = Rt^T As Rt (nonzero spectrum) and the
     factorisation B = Gt^T J Gt (J is a reflection, J^2 = I, by construction),
  5. (1D) the norm loss of the naive product F_grad P_l versus the commuted
     form Pi_l F_grad,l, which decides the block-encoding subnormalisation.
"""
import numpy as np

g = 1 / (2 * np.sqrt(3))
XI = np.array([0.5 - g, 0.5 + g])          # 2-point Gauss on [0,1]


def factors1d(L):
    """F0, F1 with F0^T F0 = M1 and F1^T F1 = K1 (interior nodes, 2-pt Gauss)."""
    h, ne = 2.0**-L, 2**L
    n = ne - 1
    F0 = np.zeros((2 * ne, n))
    F1 = np.zeros((2 * ne, n))
    w = np.sqrt(h / 2)
    for e in range(ne):
        for q in range(2):
            r = 2 * e + q
            for node, phi, dphi in [(e - 1, 1 - XI[q], -1 / h), (e, XI[q], 1 / h)]:
                if 0 <= node < n:
                    F0[r, node] = w * phi
                    F1[r, node] = w * dphi
    return F0, F1


def matrices(L, d):
    F0, F1 = factors1d(L)
    if d == 1:
        return F0, [F1]
    return np.kron(F0, F0), [np.kron(F1, F0), np.kron(F0, F1)]


def prolong1d(l):
    """Linear interpolation from level l to level l+1 (interior nodes)."""
    nc, nf = 2**l - 1, 2**(l + 1) - 1
    P = np.zeros((nf, nc))
    for j in range(nc):
        P[2 * j + 1, j] = 1.0
        P[2 * j, j] = P[2 * j + 2, j] = 0.5
    return P


def prolong(l, L, d):
    P = np.eye(2**l - 1)
    for k in range(l, L):
        P = prolong1d(k) @ P
    return P if d == 1 else np.kron(P, P)


def kmat(L, d, sb):
    F0, Fg = matrices(L, d)
    M = F0.T @ F0
    K = sum(F.T @ F for F in Fg)
    return M, K, F0, Fg


def bpx(L, d, sb):
    blocks = []
    for l in range(1, L + 1):
        Ml, Kl, _, _ = kmat(l, d, sb)
        Dl = np.diag(Ml + sb * Kl)
        blocks.append(prolong(l, L, d) / np.sqrt(Dl)[None, :])
    return np.hstack(blocks)


def nonzero(ev, rel=1e-10):
    a = np.abs(ev)
    return a[a > rel * a.max()]


def run(L, d, beta):
    sb = np.sqrt(beta)
    M, K, F0, Fg = kmat(L, d, sb)
    n = M.shape[0]
    Ab = M + sb * K
    As = np.block([[M, sb * K], [sb * K, -M]])
    ev = np.linalg.eigvalsh(As)
    kap_A = np.abs(ev).max() / np.abs(ev).min()
    # Schoeberl-Zulehner: eig of Lc^{-1} As Lc^{-T}, Ab = Lc Lc^T
    Li = np.linalg.inv(np.linalg.cholesky(Ab))
    Pi = np.kron(np.eye(2), Li)
    sz = np.abs(np.linalg.eigvalsh(Pi @ As @ Pi.T))
    # BPX
    R = bpx(L, d, sb)
    e_c = nonzero(np.linalg.eigvalsh(R.T @ Ab @ R))
    kap_C = e_c.max() / e_c.min()
    # two-sided system B = Rt^T As Rt and its factorisation Gt^T J Gt, where
    # G = [G0; G1] and J = [[E0, E1], [E1, -E0]] with E0, E1 the coordinate
    # projectors onto the G0 and G1 rows: G^T E0 G = G0^T G0, G^T E1 G = G1^T G1
    G0 = F0 @ R
    G1 = np.vstack([beta**0.25 * F @ R for F in Fg])
    G = np.vstack([G0, G1])
    B0, B1 = R.T @ M @ R, sb * (R.T @ K @ R)
    B = np.block([[B0, B1], [B1, -B0]])
    err_fac = max(np.abs(G0.T @ G0 - B0).max(), np.abs(G1.T @ G1 - B1).max()) / np.abs(B).max()
    e_b = nonzero(np.linalg.eigvalsh(B))
    return dict(n=n, kap_A=kap_A, sz_min=sz.min(), sz_max=sz.max(), kap_C=kap_C,
                kap_B=e_b.max() / e_b.min(), normG=np.linalg.norm(G, 2),
                err_fac=err_fac)


def norm_loss_1d(L):
    """||F1|| ||P_l|| / ||F1 P_l|| for each level (naive product vs. exact norm)."""
    _, F1 = factors1d(L)
    out = []
    for l in range(1, L):
        P = prolong(l, L, 1)
        out.append(np.linalg.norm(F1, 2) * np.linalg.norm(P, 2) / np.linalg.norm(F1 @ P, 2))
    return out


if __name__ == "__main__":
    betas = [1.0, 1e-2, 1e-4, 1e-6, 1e-8]
    for d, Ls in [(1, [4, 6, 8]), (2, [3, 4, 5])]:
        print(f"\n=== d = {d} ===")
        print(f"{'L':>2} {'n':>5} {'beta':>7} {'kap(As)':>10} {'SZ eig |.|':>15} "
              f"{'kap(CA_b)':>9} {'kap_eff(B)':>10} {'||G||':>6} {'fac err':>8}", flush=True)
        for L in Ls:
            for beta in betas:
                r = run(L, d, beta)
                print(f"{L:2d} {r['n']:5d} {beta:7.0e} {r['kap_A']:10.2e} "
                      f"[{r['sz_min']:.3f},{r['sz_max']:.3f}] {r['kap_C']:9.2f} "
                      f"{r['kap_B']:10.2f} {r['normG']:6.2f} {r['err_fac']:8.1e}", flush=True)
    print("\nnorm loss ||F1|| ||P_l|| / ||F1 P_l||, d = 1, L = 8, l = 1..7:")
    print(np.round(norm_loss_1d(8), 1))
