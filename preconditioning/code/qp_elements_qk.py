"""Qk elements in 1D (and 2D by Kronecker), interior DOFs (homogeneous Dirichlet).

Checks: F0^T F0 = M, F1^T F1 = K (exact Galerkin matrices), F P_l = T_l F_l for
values and derivatives, T_l^T T_l = I, composition of one-level isometries.
Optional geometry map x = Phi(xi) = xi + c xi(1-xi) on [0,1] (quadratic,
'isoparametric' for k>=2 with geometry fixed on the coarsest level).
"""
import numpy as np
from numpy.polynomial.legendre import leggauss

np.set_printoptions(precision=2)


def gauss01(n):
    x, w = leggauss(n)
    return (x + 1) / 2, w / 2


def lagrange(k, t):
    """values and derivatives of the k+1 Lagrange basis functions (equispaced
    nodes m/k on [0,1]) at points t. Returns (len(t), k+1) arrays."""
    nodes = np.linspace(0, 1, k + 1)
    V = np.vander(nodes, k + 1, increasing=True)
    C = np.linalg.inv(V)                       # columns: monomial coefficients
    Tm = np.vander(t, k + 1, increasing=True)
    dT = np.zeros_like(Tm)
    for p in range(1, k + 1):
        dT[:, p] = p * t ** (p - 1)
    return Tm @ C, dT @ C


class Geo:
    def __init__(self, c=0.0):
        self.c = c

    def J(self, xi):                            # dPhi/dxi
        return 1 + self.c * (1 - 2 * xi)


def factors(k, l, nq, geo):
    """F0 (values), F1 (x-derivatives), weighted by sqrt(w |J|), level l."""
    ne, h = 2**l, 2.0**-l
    n = k * ne - 1
    t, w = gauss01(nq)
    phi, dphi = lagrange(k, t)
    F0 = np.zeros((ne * nq, n)); F1 = np.zeros((ne * nq, n))
    for e in range(ne):
        xi = (e + t) * h
        J = geo.J(xi)
        for q in range(nq):
            r = e * nq + q
            for m in range(k + 1):
                j = k * e + m - 1
                if 0 <= j < n:
                    F0[r, j] = np.sqrt(w[q] * h * J[q]) * phi[q, m]
                    F1[r, j] = np.sqrt(w[q] * h / J[q]) * dphi[q, m] / h
    return F0, F1


def exact_MK(k, l, geo, nq=25):
    """Galerkin matrices with (practically) exact integration."""
    A0, A1 = factors(k, l, nq, geo)
    return A0.T @ A0, A1.T @ A1


def prolong(k, l, L):
    """interpolation of coarse Qk functions at fine nodes (interior DOFs)."""
    nc, nf = k * 2**l - 1, k * 2**L - 1
    xf = np.arange(1, nf + 1) / (k * 2**L)
    P = np.zeros((nf, nc))
    hc = 2.0**-l
    e = np.minimum((xf / hc).astype(int), 2**l - 1)
    t = xf / hc - e
    phi, _ = lagrange(k, t)
    for i in range(nf):
        for m in range(k + 1):
            j = k * e[i] + m - 1
            if 0 <= j < nc:
                P[i, j] += phi[i, m]
    return P


def transfer(k, l, L, nq, geo, kind):
    """T_l elementwise: weighted coarse Gauss values -> weighted fine values.
    kind 'val': weights sqrt(w h J), 'grad': sqrt(w h / J)."""
    t, w = gauss01(nq)
    r = 2**(L - l)
    hc, hf = 2.0**-l, 2.0**-L
    Vc = np.vander(t, k + 1, increasing=True)            # P_k at coarse pts
    tf = np.concatenate([(s + t) / r for s in range(r)])  # fine pts, coarse local coord
    Vf = np.vander(tf, k + 1, increasing=True)
    blocks = []
    for E in range(2**l):
        xc = (E + t) * hc
        xf = (E + tf) * hc
        if kind == 'val':
            wc = np.sqrt(w * hc * geo.J(xc)); wf = np.sqrt(np.tile(w, r) * hf * geo.J(xf))
        else:
            wc = np.sqrt(w * hc / geo.J(xc)); wf = np.sqrt(np.tile(w, r) * hf / geo.J(xf))
        blocks.append(wf[:, None] * (Vf @ np.linalg.pinv(Vc)) / wc[None, :])
    from scipy.linalg import block_diag
    return block_diag(*blocks)


def mx(A):
    return np.abs(A).max()


def check(k, L, nq, geo, label):
    F0, F1 = factors(k, L, nq, geo)
    M, K = exact_MK(k, L, geo)
    print(f"{label}: k={k} nq={nq} L={L}: |F0'F0-M|/|M|={mx(F0.T@F0-M)/mx(M):.1e} "
          f"|F1'F1-K|/|K|={mx(F1.T@F1-K)/mx(K):.1e}")
    worst = dict(v=0, g=0, iso=0, isog=0, rng=0, rngg=0, ls=0)
    for l in range(1, L):
        P = prolong(k, l, L)
        F0l, F1l = factors(k, l, nq, geo)
        Tv = transfer(k, l, L, nq, geo, 'val')
        Tg = transfer(k, l, L, nq, geo, 'grad')
        worst['v'] = max(worst['v'], mx(F0 @ P - Tv @ F0l) / mx(F0 @ P))
        worst['g'] = max(worst['g'], mx(F1 @ P - Tg @ F1l) / mx(F1 @ P))
        worst['iso'] = max(worst['iso'], mx(Tv.T @ Tv - np.eye(Tv.shape[1])))
        worst['isog'] = max(worst['isog'], mx(Tg.T @ Tg - np.eye(Tg.shape[1])))
        # isometry restricted to the relevant range (range of F_l)
        for F, T, key in [(F0l, Tv, 'rng'), (F1l, Tg, 'rngg')]:
            Q, _ = np.linalg.qr(F)
            worst[key] = max(worst[key], mx(Q.T @ T.T @ T @ Q - np.eye(Q.shape[1])))
        # same T for values and derivatives?
        worst['ls'] = max(worst['ls'], mx(Tv - Tg))
    print(f"   max_l |F0 P - T F0l|={worst['v']:.1e}  |F1 P - T F1l|={worst['g']:.1e}  "
          f"|T'T-I|: val {worst['iso']:.1e} grad {worst['isog']:.1e}  "
          f"on range(F_l): val {worst['rng']:.1e} grad {worst['rngg']:.1e}  |Tval-Tgrad|={worst['ls']:.1e}")


def check_lsq(k, L, nq):
    """No T can exist if the quadrature is not unisolvent: least-squares residual."""
    geo = Geo(0)
    F0, F1 = factors(k, L, nq, geo)
    M, K = exact_MK(k, L, geo)
    res = 0
    for l in range(1, L):
        P = prolong(k, l, L); F0l, _ = factors(k, l, nq, geo)
        Tls = (F0 @ P) @ np.linalg.pinv(F0l)
        res = max(res, mx(F0 @ P - Tls @ F0l) / mx(F0 @ P))
    print(f"under-integration k={k} nq={nq}: |F0'F0-M|/|M|={mx(F0.T@F0-M)/mx(M):.1e}, "
          f"rank F0={np.linalg.matrix_rank(F0)} of {F0.shape[1]}, best LSQ |F0 P - T F0l|={res:.1e}")


def check2d(k, L):
    geo = Geo(0)
    nq = k + 1
    F0, F1 = factors(k, L, nq, geo)
    M, K = exact_MK(k, L, geo)
    Fv, Fx, Fy = np.kron(F0, F0), np.kron(F1, F0), np.kron(F0, F1)
    M2, K2 = np.kron(M, M), np.kron(K, M) + np.kron(M, K)
    e_fac = max(mx(Fv.T @ Fv - M2) / mx(M2), mx(Fx.T @ Fx + Fy.T @ Fy - K2) / mx(K2))
    e_c = e_i = 0
    for l in range(1, L):
        P1 = prolong(k, l, L); P = np.kron(P1, P1)
        A0, A1 = factors(k, l, nq, geo)
        T1 = transfer(k, l, L, nq, geo, 'val'); T = np.kron(T1, T1)
        for F, Fl in [(Fv, np.kron(A0, A0)), (Fx, np.kron(A1, A0)), (Fy, np.kron(A0, A1))]:
            e_c = max(e_c, mx(F @ P - T @ Fl) / mx(F @ P))
        e_i = max(e_i, mx(T.T @ T - np.eye(T.shape[1])))
    print(f"2D Q{k}, L={L}: factorization {e_fac:.1e}, commuting (val,dx,dy) {e_c:.1e}, |T'T-I| {e_i:.1e}")


def check_composition(k, L):
    """T_{l->L} = product of one-level maps I (x) That_1."""
    geo = Geo(0); nq = k + 1
    for l in [1]:
        T = transfer(k, l, L, nq, geo, 'val')
        C = np.eye(2**l * nq)
        for m in range(l, L):
            C = transfer(k, m, m + 1, nq, geo, 'val') @ C
        print(f"composition k={k} l={l} L={L}: |T - prod one-level| = {mx(T - C):.1e}")


if __name__ == "__main__":
    L = 6
    for k in [1, 2, 3]:
        check(k, L, k + 1, Geo(0), "Cartesian")
    check(2, L, 4, Geo(0), "Cartesian overintegrated")
    for k in [2, 3]:
        check_lsq(k, L, k)
    for k in [1, 2, 3]:
        check2d(k, 4)
    for k in [2, 3]:
        check_composition(k, 5)
    print("\n-- geometry map x = xi + 0.5 xi(1-xi) (nonaffine, fixed on coarsest level) --")
    for k in [1, 2, 3]:
        check(k, L, k + 1, Geo(0.5), "mapped")
