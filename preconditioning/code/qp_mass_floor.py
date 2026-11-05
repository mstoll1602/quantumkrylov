"""Mass-matrix floor: kappa(CA), subnormalization and kappa_alpha for three level scalings
(Jacobi = paper; YYZ = Yang-Yu-Zhang arXiv:2610.06987; hybrid D_l = (M_l + sqrt(beta) diag(K_l))^{-1}).
See paper/nla_for_quantum.tex, section "Removing the mass-matrix floor".

Raw encodings as in the paper, plus "tight" encodings for the Jacobi and hybrid routes that
use the YYZ isometry U1 = E1 M1^{-1/2} (alpha 1) everywhere:
  values:      F_val D^1/2 = U_val * [M^{1/2} D^{1/2}]                     (diag in sine basis)
  gradient q:  F_grad^(q) D^1/2 = [U1 x .. x C1 x .. x U1] * [prod_{p!=q} mu_p^{1/2} D^{1/2}]
               alpha_q^2 = (4/h) * max_modes( prod_{p!=q} mu_p * D )
Columns: jacobi(raw), jacobi(tight), yyz, hybrid(raw), hybrid(tight)."""
import sys, numpy as np
sys.dont_write_bytecode = True
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from qp_bpx_sparse import fe, prolong

def eig1d(l):
    h = 2.0**-l; phi = np.pi*np.arange(1, 2**l)/2**l
    return h/3*(2+np.cos(phi)), 2/h*(1-np.cos(phi))

def modes(l, d):
    """per-mode arrays: mu (M eig), prod_{p!=q} mu_p for each q."""
    m1, _ = eig1d(l)
    grids = np.meshgrid(*([m1]*d), indexing='ij')
    mus = [g.ravel() for g in grids]
    mu = np.prod(mus, axis=0)
    others = [mu/mq for mq in mus]
    return mu, others

def msqrt(X, p):
    w, V = np.linalg.eigh(X); return (V*np.maximum(w, 0)**p) @ V.T

def run(d, L, beta, prec):
    nu = np.sqrt(beta); tL = nu*4.0**L
    M, K = fe(L, d); M = M.toarray(); K = K.toarray(); A = M + nu*K
    n = len(M); C = np.zeros((n, n)); a2 = {}
    for l in range(1, L+1):
        h = 2.0**-l; t = nu*4.0**l
        P = prolong(l, L, d).toarray()
        Ml, Kl = fe(l, d); Ml = Ml.toarray(); kd = Kl[0, 0]
        mu, others = modes(l, d)
        th = (min(1, t)/min(1, tL))
        raw_v2 = (4*h/3)**d; raw_g2 = d*(4/h)*(4*h/3)**(d-1)
        if prec == 'jacobi':
            Dd = 1/(Ml[0, 0] + nu*kd) * np.ones_like(mu)
            C += th*Dd[0] * P @ P.T
            add = {'raw': th*Dd[0]*(raw_v2 + nu*raw_g2),
                   'tight': th*((mu*Dd).max() + nu*sum((4/h)*(o*Dd).max() for o in others))}
        elif prec == 'yyz':
            C += 1/(1+t) * P @ np.linalg.inv(Ml) @ P.T
            add = {'raw': 1/(1+t)*(1 + nu*12*d/h**2)}
        else:
            c = nu*kd; Dd = 1/(mu + c)
            C += th * P @ np.linalg.inv(Ml + c*np.eye(len(Ml))) @ P.T
            g2 = mu*Dd
            add = {'raw': th*(g2.max() + nu*min(12*d/h**2*g2.max(), raw_g2/(mu.min()+c))),
                   'tight': th*(g2.max() + nu*sum((4/h)*(o*Dd).max() for o in others))}
        for k, v in add.items(): a2[k] = a2.get(k, 0) + v
    Ch = msqrt(C, 0.5)
    lam = np.linalg.eigvalsh(Ch @ A @ Ch)
    kC, G2 = lam.max()/lam.min(), lam.max()
    Z = np.zeros_like(A); Cb = np.block([[Ch, Z], [Z, Ch]])
    ev = np.abs(np.linalg.eigvalsh(Cb @ np.block([[M, nu*K], [nu*K, -M]]) @ Cb))
    smin = ev[ev > 1e-10*ev.max()].min()
    return kC, {k: (v/G2, v/smin) for k, v in a2.items()}

cols = [('jacobi', 'raw'), ('jacobi', 'tight'), ('yyz', 'raw'), ('hybrid', 'raw'), ('hybrid', 'tight')]
print("columns: jacobi-raw jacobi-tight yyz hybrid-raw hybrid-tight")
print(f"{'d':>2}{'L':>3}{'beta':>8} | kappa(CA) jac/yyz/hyb | alpha_G^2/||G||^2 (5 cols) | kappa_alpha (5 cols)")
for d, L in [(1, 8), (2, 5), (3, 3)]:
    for beta in [1, 1e-2, 1e-4, 1e-6, 1e-8, 1e-12]:
        res = {p: run(d, L, beta, p) for p in ('jacobi', 'yyz', 'hybrid')}
        row = f"{d:>2}{L:>3}{beta:>8.0e} |" + "".join(f"{res[p][0]:6.1f}" for p in ('jacobi', 'yyz', 'hybrid')) + " |"
        row += "".join(f"{res[p][1][e][0]:6.1f}" for p, e in cols) + " |"
        row += "".join(f"{res[p][1][e][1]:6.1f}" for p, e in cols)
        print(row, flush=True)
