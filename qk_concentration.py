"""Basis dependence of sampling, for the sketching section of the note.

Transverse-field Ising chain with N spins (n = 2^N), start vector |+...+>.
For the ground state u0 and the real-time Krylov vector v_10 = U^10 v0 we count
how many basis vectors carry 99% of the squared norm, once in the standard
basis and once in the Walsh-Hadamard basis F = h (x) ... (x) h. Sampling-based
methods only work if this number is small.
"""
import numpy as np

N, J = 10, 1.0
n = 2**N
I2 = np.eye(2)
X = np.array([[0, 1], [1, 0]], dtype=float)
Z = np.diag([1.0, -1.0])
h = np.array([[1, 1], [1, -1]]) / np.sqrt(2)


def op(single, site):
    out = np.array([[1.0]])
    for s in range(N):
        out = np.kron(out, single if s == site else I2)
    return out


def count99(w):
    """Number of largest entries of w that carry 99% of ||w||^2."""
    p = np.sort(np.abs(w)**2)[::-1]
    return np.searchsorted(np.cumsum(p), 0.99) + 1


F = np.array([[1.0]])
for _ in range(N):
    F = np.kron(F, h)

v0 = np.ones(n) / np.sqrt(n)          # |+...+>
ZZ = sum(op(Z, i) @ op(Z, i + 1) for i in range(N - 1))
XX = sum(op(X, i) for i in range(N))
for hx in [0.5, 1.0, 2.0]:
    H = -J * ZZ - hx * XX
    lam, V = np.linalg.eigh(H)
    dt = np.pi / (lam[-1] - lam[0])
    v10 = V @ (np.exp(-1j * 10 * dt * lam) * (V.T @ v0))
    for name, w in [("u0", V[:, 0]), ("v_10", v10)]:
        print(f"h={hx}  {name:5s} standard: {count99(w):4d}   Hadamard: {count99(F @ w):4d}   (of {n})")
