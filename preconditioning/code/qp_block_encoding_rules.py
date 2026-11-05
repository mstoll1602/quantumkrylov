"""Explicit unitary realizations of the block-encoding combination rules.

Companion to paper/nla_for_quantum.tex. Every construction is built as a dense
unitary matrix and checked: unitarity and the leading block. Ordering convention:
ancilla/index registers are the OUTER Kronecker factors, the system is innermost,
so the "leading block" (all ancillas in state 0) is the top-left n x n block.
"""
import numpy as np

rng = np.random.default_rng(0)


def crand(m, n):
    return rng.standard_normal((m, n)) + 1j * rng.standard_normal((m, n))


def psd_sqrt(X):
    w, V = np.linalg.eigh(X)
    return V @ np.diag(np.sqrt(np.maximum(w, 0))) @ V.conj().T


def dilation(A, alpha):
    """Halmos dilation: unitary with leading block A/alpha (one ancilla qubit)."""
    n = A.shape[0]
    X = A / alpha
    I = np.eye(n)
    return np.block([[X, psd_sqrt(I - X @ X.conj().T)],
                     [psd_sqrt(I - X.conj().T @ X), -X.conj().T]])


def unitary_with_first_column(p):
    """Householder-type unitary V with V e_0 = p (p a unit vector)."""
    k = len(p)
    Q, _ = np.linalg.qr(np.column_stack([p, crand(k, k - 1)]))
    return Q * (p[0] / Q[0, 0]) if abs(Q[0, 0]) > 0 else Q   # fix phase of first column


def err_unitary(U):
    return np.abs(U.conj().T @ U - np.eye(U.shape[0])).max()


def report(name, U, block, target):
    print(f"{name:34s} unitary err {err_unitary(U):.1e}   block err {np.abs(block - target).max():.1e}")


n, m = 4, 3
A, B = crand(n, n), crand(n, n)
aA, aB = 1.3 * np.linalg.norm(A, 2), 1.1 * np.linalg.norm(B, 2)
UA, UB = dilation(A, aA), dilation(B, aB)
report("single matrix (Halmos dilation)", UA, UA[:n, :n], A / aA)

# products: separate ancillas, ordering (ancA, ancB, system)
I2 = np.eye(2)
UA_t = np.einsum('aibj,cd->acibdj', UA.reshape(2, n, 2, n), I2).reshape(4 * n, 4 * n)
UB_t = np.einsum('cidj,ab->acibdj', UB.reshape(2, n, 2, n), I2).reshape(4 * n, 4 * n)
W = UA_t @ UB_t
report("product AB (fresh ancillas)", W, W[:n, :n], A @ B / (aA * aB))
W_bad = UA @ UB
print(f"{'product AB, SAME ancilla (wrong)':34s} block err {np.abs(W_bad[:n, :n] - A @ B / (aA * aB)).max():.1e}")

# Kronecker product: U_A (x) U_B with ordering (ancA, sysA, ancB, sysB) -> (ancA, ancB, sysA, sysB)
C = crand(m, m); aC = 1.2 * np.linalg.norm(C, 2); UC = dilation(C, aC)
K = np.einsum('aibj,ckdl->acikbdjl', UA.reshape(2, n, 2, n), UC.reshape(2, m, 2, m)).reshape(4 * n * m, 4 * n * m)
report("Kronecker A (x) C", K, K[:n * m, :n * m], np.kron(A, C) / (aA * aC))

# LCU: sum_j c_j A_j with k = 4 terms (phases of c_j absorbed into U_j)
k = 4
As = [crand(n, n) for _ in range(k)]
als = [1.2 * np.linalg.norm(X, 2) for X in As]
cs = crand(k, 1).ravel()
Us = [np.exp(1j * np.angle(c)) * dilation(X, a) for c, X, a in zip(cs, As, als)]
s = sum(abs(c) * a for c, a in zip(cs, als))
p = np.sqrt(np.array([abs(c) * a for c, a in zip(cs, als)]) / s)
V = unitary_with_first_column(p)
SEL = np.zeros((k * 2 * n, k * 2 * n), complex)
for j, U in enumerate(Us):
    SEL[j * 2 * n:(j + 1) * 2 * n, j * 2 * n:(j + 1) * 2 * n] = U
VI = np.kron(V, np.eye(2 * n))
W = VI.conj().T @ SEL @ VI
report("LCU sum c_j A_j", W, W[:n, :n], sum(c * X for c, X in zip(cs, As)) / s)

# stacking [A_1; ...; A_k]: W = SEL (V (x) I), p_j = alpha_j / N
N = np.sqrt(sum(a**2 for a in als))
Us0 = [dilation(X, a) for X, a in zip(As, als)]
SEL0 = np.zeros_like(SEL)
for j, U in enumerate(Us0):
    SEL0[j * 2 * n:(j + 1) * 2 * n, j * 2 * n:(j + 1) * 2 * n] = U
Vs = unitary_with_first_column(np.array(als) / N)
W = SEL0 @ np.kron(Vs, np.eye(2 * n))
W6 = W.reshape(k, 2, n, k, 2, n)
stack_block = W6[:, 0, :, 0, 0, :].reshape(k * n, n)        # rows: (j, anc=0), cols: (j=0, anc=0)
report("stacking [A_1; ...; A_k]", W, stack_block, np.vstack(As) / N)

# concatenation [A_1, ..., A_k]: W = (V^H (x) I) SEL, p_j = alpha_j / N
W = np.kron(Vs, np.eye(2 * n)).conj().T @ SEL0
W6 = W.reshape(k, 2, n, k, 2, n)
conc_block = W6[0, 0, :, :, 0, :].reshape(n, k * n)        # rows: (j=0, anc=0), cols: (j, anc=0)
report("concatenation [A_1, ..., A_k]", W, conc_block, np.hstack(As) / N)

# sparse matrices: A / (amax sqrt(s_r s_c)) = T_L^H T_R with isometries T_L, T_R
nn = 8
S = np.zeros((nn, nn), complex)
for i in range(nn):
    for j in [i - 1, i, i + 1, (i + 3) % nn]:
        if 0 <= j < nn:
            S[i, j] = complex(*rng.standard_normal(2))
sr = int((S != 0).sum(1).max()); sc = int((S != 0).sum(0).max())
amax = np.abs(S).max()
TL = np.zeros((nn * nn * 2, nn), complex)            # columns psi_i in C^n (x) C^n (x) C^2
TR = np.zeros((nn * nn * 2, nn), complex)            # columns phi_j
idx = lambda i, j, f: (i * nn + j) * 2 + f
for i in range(nn):
    cols = list(np.nonzero(S[i])[0])
    cols += [j for j in range(nn) if j not in cols][:sr - len(cols)]   # pad to s_r entries
    for j in cols:
        a = S[i, j] / amax
        TL[idx(i, j, 0), i] = np.conj(a) / np.sqrt(sr)
        TL[idx(i, j, 1), i] = np.sqrt(1 - abs(a)**2) / np.sqrt(sr)
for j in range(nn):
    rows = list(np.nonzero(S[:, j])[0])
    rows += [i for i in range(nn) if i not in rows][:sc - len(rows)]   # pad to s_c entries
    for i in rows:
        TR[idx(i, j, 0), j] = 1 / np.sqrt(sc)
def complete(T):
    """Unitary whose first columns are the orthonormal columns of T. QR of [T, X]
    reproduces T up to unimodular column factors, so we put T back in place."""
    Q = np.linalg.qr(np.column_stack([T, crand(T.shape[0], T.shape[0] - T.shape[1])]))[0]
    Q[:, :T.shape[1]] = T
    return Q


err_iso = lambda T: np.abs(T.conj().T @ T - np.eye(T.shape[1])).max()
UL, UR = complete(TL), complete(TR)
W = UL.conj().T @ UR
print(f"{'sparse: isometries T_L, T_R':34s} |T_L^H T_L - I| {err_iso(TL):.1e}  |T_R^H T_R - I| {err_iso(TR):.1e}")
report(f"sparse (s_r={sr}, s_c={sc})", W, W[:nn, :nn], S / (amax * np.sqrt(sr * sc)))
print(f"   alpha = amax*sqrt(s_r s_c) = {amax * np.sqrt(sr * sc):.2f}  vs ||A|| = {np.linalg.norm(S, 2):.2f}")

# amplitude amplification: rotation by 2*theta in span{g, h}, sin(theta) = a
nA = 8
A = crand(nA, nA) + 6 * np.eye(nA)                  # moderately conditioned
U = dilation(A, 3.0 * np.linalg.norm(A, 2))           # loose alpha -> small success amplitude
b = crand(nA, 1).ravel(); b /= np.linalg.norm(b)
u = np.concatenate([b, np.zeros(nA)])                 # e_0 (x) b, ancilla outermost
psi = U @ u
Pi = np.diag(np.concatenate([np.ones(nA), np.zeros(nA)]))
a = np.linalg.norm(Pi @ psi)
th = np.arcsin(a)
g = Pi @ psi / a
hh = (psi - Pi @ psi) / np.linalg.norm(psi - Pi @ psi)
S_Pi = np.eye(2 * nA) - 2 * Pi
S_psi = U @ (np.eye(2 * nA) - 2 * np.outer(u, u.conj())) @ U.conj().T
Q = -S_psi @ S_Pi
kopt = int(np.floor(np.pi / (4 * th)))
err, v = 0.0, psi.copy()
for k in range(1, kopt + 1):
    v = Q @ v
    err = max(err, np.abs(v - (np.sin((2 * k + 1) * th) * g + np.cos((2 * k + 1) * th) * hh)).max())
print(f"{'amplitude amplification':34s} a = {a:.3f}, rounds = {kopt}, rotation err {err:.1e}, "
      f"success prob {a**2:.3f} -> {np.linalg.norm(Pi @ v)**2:.3f}")
