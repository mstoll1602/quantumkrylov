"""Level weights theta_l = min(1, sqrt(beta)/h_l^2)^s for the BPX factor (Table 2 of
quantum_preconditioning.tex). s = 0 is plain diagonally scaled BPX."""
import numpy as np

from qp_poisson_check import kmat, prolong, nonzero

def kappa(L, d, beta, s):
    """BPX with Jacobi weights times theta_l = min(1, sqrt(beta)/h_l^2)^s."""
    sb = np.sqrt(beta)
    M, K, _, _ = kmat(L, d, sb)
    Ab = M + sb * K
    bl = []
    for l in range(1, L + 1):
        Ml, Kl, _, _ = kmat(l, d, sb)
        delta = (Ml + sb * Kl)[0, 0]
        theta = min(1.0, sb * 4.0**l) ** s
        bl.append(prolong(l, L, d) * np.sqrt(theta / delta))
    R = np.hstack(bl)
    e = nonzero(np.linalg.eigvalsh(R.T @ Ab @ R))
    return e.max() / e.min()

if __name__ == "__main__":
    for d, Ls in [(1, [6, 8, 10]), (2, [4, 5])]:
        for L in Ls:
            for beta in [1, 1e-2, 1e-4, 1e-6, 1e-8, 1e-10]:
                ks = [kappa(L, d, beta, s) for s in (0, 0.5, 1)]
                print(f"d={d} L={L:2d} beta={beta:7.0e}  s=0: {ks[0]:6.2f}  s=1/2: {ks[1]:6.2f}  s=1: {ks[2]:6.2f}", flush=True)
