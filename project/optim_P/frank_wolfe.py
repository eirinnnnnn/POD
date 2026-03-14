import numpy as np
from scipy.optimize import linear_sum_assignment, minimize_scalar

# ============================================================
# 1) Input: (16,7) eBCH generator matrix
# ============================================================
G_str = [
    "1100010111000000",
    "0110011100100000",
    "0011001110010000",
    "1101110000001000",
    "1010111000000100",
    "1001011100000010",
    "1000101110000001",
]

G = np.array([[int(c) for c in row] for row in G_str], dtype=float)
K, N = G.shape
assert (K, N) == (7, 16)

# ============================================================
# 2) Standard polar transform with bit-reversal:
#    G_N = B_N F^{⊗n}
# ============================================================
def kron_power(F, n):
    M = np.array([[1.0]])
    for _ in range(n):
        M = np.kron(M, F)
    return M

def bit_reverse_perm(N):
    m = int(np.log2(N))
    perm = []
    for i in range(N):
        b = format(i, f"0{m}b")
        perm.append(int(b[::-1], 2))
    return perm

F = np.array([[1.0, 0.0],
              [1.0, 1.0]])

n = int(np.log2(N))
assert 2**n == N

B_N = np.eye(N)[bit_reverse_perm(N)]
G_N = B_N @ kron_power(F, n)

# ============================================================
# 3) Bhattacharyya parameters p_i
#    Here I use BEC(eps=0.5), then reorder with bit-reversal
# ============================================================
def bec_bhattacharyya(N, eps=0.5):
    Z = [eps]
    while len(Z) < N:
        newZ = []
        for z in Z:
            z_bad = 2*z - z*z
            z_good = z*z
            newZ.extend([z_bad, z_good])
        Z = newZ
    return np.array(Z, dtype=float)

Z_plain = bec_bhattacharyya(N, eps=0.5)
br = bit_reverse_perm(N)
p = Z_plain[br]  # bit-reversal-consistent order

# ============================================================
# 4) Objective:
#    f(P) = sum_{i=1}^{N-1} (p_i - p_{i+1}) ||Y_i(P)||_*
# ============================================================
def objective(P, G, G_N, p):
    N = G_N.shape[0]
    total = 0.0
    for i in range(1, N):
        E_i = G_N[:, :i]
        Y_i = G @ P.T @ E_i
        s = np.linalg.svd(Y_i, compute_uv=False)
        total += (p[i-1] - p[i]) * np.sum(s)
    return total

# ============================================================
# 5) Gradient / subgradient matrix C
#    C = sum_i (p_i - p_{i+1}) G^T U_i V_i^T E_i^T
# ============================================================
def build_C(P, G, G_N, p):
    N = G_N.shape[0]
    C = np.zeros((N, N), dtype=float)

    for i in range(1, N):
        E_i = G_N[:, :i]
        Y_i = G @ P.T @ E_i
        U, s, Vt = np.linalg.svd(Y_i, full_matrices=False)
        UVt = U @ Vt
        C += (p[i-1] - p[i]) * (G.T @ UVt @ E_i.T)

    return C

# ============================================================
# 6) Hungarian for minimization
#    Solve min_S Tr(C^T S)
# ============================================================
def hungarian_min(C):
    row_ind, col_ind = linear_sum_assignment(C)   # minimize directly
    S = np.zeros_like(C)
    S[row_ind, col_ind] = 1.0
    return S, col_ind

# ============================================================
# 7) Build permutation matrix from a 0-based permutation
# ============================================================
def perm_to_matrix(perm):
    N = len(perm)
    P = np.zeros((N, N), dtype=float)
    for r, c in enumerate(perm):
        P[r, c] = 1.0
    return P

# ============================================================
# 8) Frank-Wolfe minimization
# ============================================================
def frank_wolfe_min(
    G, G_N, p,
    max_iter=50,
    tol=1e-8,
    init="identity",
    verbose=True
):
    N = G.shape[1]

    if init == "identity":
        P = np.eye(N)
    elif init == "random":
        perm0 = np.random.permutation(N)
        P = np.eye(N)[perm0]
    else:
        raise ValueError("init must be 'identity' or 'random'")

    history = []

    for k in range(max_iter):
        f_val = objective(P, G, G_N, p)
        C = build_C(P, G, G_N, p)

        # minimization oracle
        S, perm_S = hungarian_min(C)
        D = S - P

        # FW directional derivative <grad, D>
        fw_dir = np.sum(C * D)
        history.append((f_val, fw_dir))

        if verbose:
            print(f"[iter {k:02d}] f(P) = {f_val:.10f}, <grad,D> = {fw_dir:.10e}")

        if abs(fw_dir) < tol:
            break

        # line search: minimize objective((1-gamma)P + gamma S)
        def phi(gamma):
            P_new = (1.0 - gamma) * P + gamma * S
            return objective(P_new, G, G_N, p)

        res = minimize_scalar(
            phi,
            bounds=(0.0, 1.0),
            method="bounded",
            options={"xatol": 1e-4}
        )
        gamma = float(res.x)

        # fallback step size if needed
        if gamma < 1e-8:
            gamma = 2.0 / (k + 2.0)

        P = (1.0 - gamma) * P + gamma * S

    # final rounding:
    # nearest permutation in Frobenius norm = maximize Tr(P^T S)
    row_ind, col_ind = linear_sum_assignment(-P)
    P_round = np.zeros_like(P)
    P_round[row_ind, col_ind] = 1.0

    return P, P_round, col_ind, history

# ============================================================
# 9) GF(2) utilities to extract pivot set I(P)
# ============================================================
def gf2_matmul(A, B):
    return ((A.astype(np.uint8) @ B.astype(np.uint8)) % 2).astype(np.uint8)

def gf2_rref_pivots(A):
    A = A.copy().astype(np.uint8)
    m, n = A.shape
    row = 0
    pivots = []

    for col in range(n):
        pivot = None
        for r in range(row, m):
            if A[r, col]:
                pivot = r
                break
        if pivot is None:
            continue

        if pivot != row:
            A[[row, pivot]] = A[[pivot, row]]

        for r in range(m):
            if r != row and A[r, col]:
                A[r, :] ^= A[row, :]

        pivots.append(col)
        row += 1
        if row == m:
            break

    return pivots

def pivot_set_from_perm(perm, G, G_N):
    P = perm_to_matrix(perm).astype(np.uint8)
    M = gf2_matmul(gf2_matmul(G.astype(np.uint8), P.T), G_N.astype(np.uint8))
    return gf2_rref_pivots(M)

# ============================================================
# 10) Run solver
# ============================================================
P_soft, P_final, perm_final, history = frank_wolfe_min(
    G, G_N, p,
    max_iter=50,
    tol=1e-8,
    init="identity",
    verbose=True
)

perm_final = perm_final.tolist()
I_final = pivot_set_from_perm(perm_final, G, G_N)
sum_final = float(sum(p[i] for i in I_final))

print("\n=== New minimization-based FW result ===")
print("perm (0-based) =", perm_final)
print("pivot set I =", I_final)
print("sum_{i in I} p_i =", sum_final)

# ============================================================
# 11) Compare with the previous permutation
# ============================================================
prev_perm = [12, 9, 11, 3, 4, 13, 7, 15, 14, 8, 1, 2, 0, 10, 5, 6]
I_prev = pivot_set_from_perm(prev_perm, G, G_N)
sum_prev = float(sum(p[i] for i in I_prev))

print("\n=== Previous permutation ===")
print("perm (0-based) =", prev_perm)
print("pivot set I =", I_prev)
print("sum_{i in I} p_i =", sum_prev)