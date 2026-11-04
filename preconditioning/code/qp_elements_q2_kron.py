"""Direct 2D Q2 assembly (9-node elements, 3x3 Gauss) vs. Kronecker form."""
import numpy as np
L = 3; ne = 2**L; h = 1.0 / ne
nn = 2 * ne + 1                      # 1D nodes incl. boundary
g, w = np.polynomial.legendre.leggauss(3); g = (g + 1) / 2; w = w / 2
def phi(x):  return np.array([2*(x-.5)*(x-1), -4*x*(x-1), 2*x*(x-.5)])
def dphi(x): return np.array([4*x-3, -8*x+4, 4*x-1])
# 1D assembly
M1 = np.zeros((nn, nn)); K1 = np.zeros((nn, nn))
for e in range(ne):
    idx = [2*e, 2*e+1, 2*e+2]
    for xq, wq in zip(g, w):
        p, dp = phi(xq), dphi(xq) / h
        M1[np.ix_(idx, idx)] += h * wq * np.outer(p, p)
        K1[np.ix_(idx, idx)] += h * wq * np.outer(dp, dp)
# direct 2D assembly with 9-node elements; global index (i, j) -> j*nn + i (x inner)
N = nn * nn; M = np.zeros((N, N)); K = np.zeros((N, N))
for ey in range(ne):
    for ex in range(ne):
        ix = [2*ex, 2*ex+1, 2*ex+2]; iy = [2*ey, 2*ey+1, 2*ey+2]
        idx = [jy * nn + jx for jy in iy for jx in ix]
        for xq, wx in zip(g, w):
            for yq, wy in zip(g, w):
                px, py = phi(xq), phi(yq); dx, dy = dphi(xq) / h, dphi(yq) / h
                P  = np.kron(py, px)                       # 9 basis values
                Gx = np.kron(py, dx); Gy = np.kron(dy, px)  # gradients
                wt = h * h * wx * wy
                M[np.ix_(idx, idx)] += wt * np.outer(P, P)
                K[np.ix_(idx, idx)] += wt * (np.outer(Gx, Gx) + np.outer(Gy, Gy))
inner = np.arange(1, nn - 1)
I2 = np.ix_([jy*nn+jx for jy in inner for jx in inner], [jy*nn+jx for jy in inner for jx in inner])
m1, k1 = M1[np.ix_(inner, inner)], K1[np.ix_(inner, inner)]
print("n per direction:", len(inner), "= 2^(L+1)-1 =", 2**(L+1)-1)
print("|M - M1(x)M1|        :", np.abs(M[I2] - np.kron(m1, m1)).max())
print("|K - (K1(x)M1+M1(x)K1)|:", np.abs(K[I2] - (np.kron(k1, m1) + np.kron(m1, k1))).max())
print("element mass * 30/h:\n", np.round(30/h*M1[:3,:3] , 6))
