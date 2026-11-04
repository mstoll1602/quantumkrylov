"""Gram consistency and local transfer maps for Qk with under-integration and
mapped geometry (supplements qp_elements_qk.py)."""
import numpy as np
from qp_elements_qk import *

def gram(k, L, nq, geo):
    F0, F1 = factors(k, L, nq, geo)
    out = []
    for l in range(1, L):
        P = prolong(k, l, L); A0, A1 = factors(k, l, nq, geo)
        e0 = mx(P.T @ F0.T @ F0 @ P - A0.T @ A0) / mx(A0.T @ A0)
        e1 = mx(P.T @ F1.T @ F1 @ P - A1.T @ A1) / mx(A1.T @ A1)
        # local T on one coarse element: does a (k+1)x? map exist? residual of local LSQ
        out.append((l, e0, e1))
    return out

L = 6
for k, nq, c in [(2, 2, 0), (3, 3, 0), (2, 3, 0.5), (3, 4, 0.5)]:
    print(f"k={k} nq={nq} c={c}: Gram consistency |P'F'FP - F_l'F_l|/|.| per level (val, grad):")
    for l, e0, e1 in gram(k, L, nq, Geo(c)):
        print(f"   l={l}: {e0:.1e} {e1:.1e}")

# local T for under-integration: element-local LSQ residual on one coarse element (l=L-1)
for k, nq in [(2, 2), (3, 3)]:
    t, w = gauss01(nq)
    Vc = np.vander(t, k + 1, increasing=True)
    tf = np.concatenate([(s + t) / 2 for s in range(2)])
    Vf = np.vander(tf, k + 1, increasing=True)
    Tls = Vf @ np.linalg.pinv(Vc)
    print(f"under-int k={k} nq={nq}: local residual |Vf - T Vc| = {mx(Vf - Tls @ Vc):.1e} (no local T)")
