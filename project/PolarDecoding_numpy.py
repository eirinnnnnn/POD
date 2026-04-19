import os
import numpy as np
import numba as nb

# def load_binary_matrix(path: str):
#     assert os.path.isfile(path)
#     with open(path, 'r') as din:
#         c = din.read().splitlines()
#     row, col = map(int, c[0].split())
#     matrix = [int(s, 2) for s in c[1:row+1]]
#     return row, col, matrix
def load_binary_matrix(path: str):
    assert os.path.isfile(path)
    with open(path, 'r') as din:
        c = din.read().splitlines()
    row, col = map(int, c[0].split())
    matrix = []
    for line in c[1:row+1]:
        if ' ' in line:
            matrix.append(int(''.join(line.split()), 2))
        else:
            matrix.append(int(line, 2))
    return row, col, matrix


def binary_matrix_transprose(matrix: list[int], col: int):
    m = [0] * col
    for r in matrix:
        for cdx in range(col):
            bit = (r >> (col - 1 - cdx)) & 1
            m[cdx] = (m[cdx] << 1) | bit
    return m

def inv_eirin(perm: list[int]):
    p = np.array(perm)
    p = np.argsort(perm).tolist()
    print(p)
    return p


# ── packed-bit helpers (Python, used in __debug__ only) ───────────────────────

def _int_to_packed(val: int, words: int) -> np.ndarray:
    """Convert a Python int to a uint64[words] packed row (LSB-first)."""
    row = np.zeros(words, dtype=np.uint64)
    for w in range(words):
        row[w] = np.uint64((val >> (w * 64)) & 0xFFFFFFFFFFFFFFFF)
    return row


def _packed_to_int(row: np.ndarray) -> int:
    """Convert a uint64[words] packed row back to a Python int."""
    val = 0
    for w, word in enumerate(row):
        val |= int(word) << (w * 64)
    return val


def _packed_mat_to_int_list(mat: np.ndarray) -> list:
    return [_packed_to_int(mat[i]) for i in range(mat.shape[0])]


# ── numba JIT functions ────────────────────────────────────────────────────────
# Each N-bit row is stored as uint64[words] (LSB-first packing).
# words = ceil(N / 64)
# bit b of row r:  (mat[r, b//64] >> (b%64)) & 1
# set bit b:        mat[r, b//64] |= 1 << (b%64)

@nb.njit(cache=True, inline='always')
def _ctz64(x):
    """Count trailing zeros of a non-zero uint64."""
    c = nb.uint64(0)
    while (x & nb.uint64(1)) == nb.uint64(0):
        c += nb.uint64(1)
        x >>= nb.uint64(1)
    return c


@nb.njit(cache=True)
def _compute_kernel_bha(N, words, bha_init, ops, na_arr, nb_arr):
    """Kernel construction + bha update.

    kernel[i] starts as identity with bit (N-1-i) set (anti-diagonal,
    matching kernel_col_idx = [[N-1],[N-2],...,[0]]).
    XOR step replaces kernel_col_idx[nb].extend(kernel_col_idx[na]).
    """
    kernel = np.zeros((N, words), dtype=np.uint64)
    for i in range(N):
        b = N - 1 - i
        kernel[i, b // 64] = nb.uint64(1) << nb.uint64(b % 64)

    bha = bha_init.copy()

    for i in range(len(ops)):
        if ops[i]:
            for w in range(words):
                kernel[nb_arr[i], w] ^= kernel[na_arr[i], w]
            prod = bha[na_arr[i]] * bha[nb_arr[i]]
            bha[na_arr[i]] = bha[na_arr[i]] + bha[nb_arr[i]] - prod
            bha[nb_arr[i]] = prod

    return kernel, bha


@nb.njit(cache=True)
def _build_GP_T(kernel, permutation, N, words):
    """Build GP_T from kernel and permutation.

    For each kernel row k_idx, iterates over set bits x (original columns).
    dst = permutation[x]; sets bit (N-1-k_idx) in GP_T[dst].
    """
    GP_T = np.zeros((N, words), dtype=np.uint64)
    for k_idx in range(N):
        rb = N - 1 - k_idx                    # bit position of row_bit
        rb_w = rb // 64
        rb_mask = nb.uint64(1) << nb.uint64(rb % 64)

        for w in range(words):
            temp = kernel[k_idx, w]
            while temp != nb.uint64(0):
                lsb = temp & (~temp + nb.uint64(1))
                x = w * 64 + int(_ctz64(lsb))  # original column index
                dst = permutation[x]
                GP_T[dst, rb_w] |= rb_mask
                temp ^= lsb

    return GP_T


@nb.njit(cache=True)
def _build_GPH_T(GP_T, H_T_cols, H_T_nnz, H_rows, words):
    """Build GPH_T using precomputed H_T column indices."""
    GPH_T = np.zeros((H_rows, words), dtype=np.uint64)
    for h_idx in range(H_rows):
        for j in range(H_T_nnz[h_idx]):
            col = H_T_cols[h_idx, j]
            for w in range(words):
                GPH_T[h_idx, w] ^= GP_T[col, w]
    return GPH_T


@nb.njit(cache=True)
def _gaussian_elimination_v1(matrix, N, words):
    """Forward Gaussian Elimination (rows below pivot only).

    Scans bit positions 0 (LSB) to N-1 (MSB), matching the original
    shift = 1 → 1<<N loop.

    Returns (matrix_in_place, frozen) where frozen[words] is the packed
    bitmask of pivot columns (== lead_bit in original).
    """
    H_rows = matrix.shape[0]
    frozen = np.zeros(words, dtype=np.uint64)

    l = 0
    for b in range(N):
        if l >= H_rows:
            break
        w       = b // 64
        bit_pos = nb.uint64(b % 64)
        bit_mask = nb.uint64(1) << bit_pos

        t_index = -1
        for i in range(l, H_rows):
            if matrix[i, w] & bit_mask:
                t_index = i
                break

        if t_index == -1:
            continue

        if l != t_index:
            for k in range(words):
                matrix[l, k], matrix[t_index, k] = matrix[t_index, k], matrix[l, k]

        frozen[w] |= bit_mask

        pivot_row = l
        for mdx in range(t_index + 1, H_rows):
            if matrix[mdx, w] & bit_mask:
                for k in range(words):
                    matrix[mdx, k] ^= matrix[pivot_row, k]

        l += 1

    return matrix, frozen


# ── v2-specific JIT helpers ───────────────────────────────────────────────────

@nb.njit(cache=True)
def _compute_PH_bha_v2(H_np, permutation, N, words, bha_init, ops, na_arr, nb_arr):
    """Build PH = H[permutation] then XOR rows and update bha."""
    PH = np.empty((N, words), dtype=np.uint64)
    for i in range(N):
        for w in range(words):
            PH[i, w] = H_np[permutation[i], w]

    bha = bha_init.copy()

    for i in range(len(ops)):
        if ops[i]:
            for w in range(words):
                PH[nb_arr[i], w] ^= PH[na_arr[i], w]
            prod = bha[na_arr[i]] * bha[nb_arr[i]]
            bha[na_arr[i]] = bha[na_arr[i]] + bha[nb_arr[i]] - prod
            bha[nb_arr[i]] = prod

    return PH, bha


@nb.njit(cache=True)
def _packed_transpose_nxn(arr, N, words):
    """Packed N×N bit matrix transpose (matches binary_matrix_transprose).

    Input:  arr[N, words]  — row r, bit c stored at arr[r, c//64] bit c%64
    Output: res[N, words]  — res[c] bit r = arr[r] bit c
    """
    res = np.zeros((N, words), dtype=np.uint64)
    for r in range(N):
        for c in range(N):
            bit = (arr[r, c // 64] >> nb.uint64(c % 64)) & nb.uint64(1)
            if bit:
                res[c, r // 64] |= nb.uint64(1) << nb.uint64(r % 64)
    return res


# ── PolarDecodeSystem ─────────────────────────────────────────────────────────

class PolarDecodeSystem():
    """
    PolarDecodeSystem could initial by H matrix file
    matrix_T is restore matrix in transprose
    """
    def __init__():
        raise RuntimeError("Please initial by PolarDecodeSystem.from_*")

    @classmethod
    def from_G_matrix_file(cls,
                           level: int,
                           matrix_path: str = None,
                           bha_init: float = 0.198997):
        raise RuntimeError("TODO")

    @classmethod
    def from_H_matrix_file(cls,
                           level: int,
                           matrix_path: str,
                           bha_init: float = 0.198997):
        obj = object.__new__(cls)
        row, col, obj.H_ori = load_binary_matrix(matrix_path)
        obj.N_ori, obj.H_ori_T = row, binary_matrix_transprose(obj.H_ori, col)
        obj.init_system(level=level, bha_init=bha_init)
        return obj

    def init_system(self, level: int, bha_init: float):
        # init once
        self.N, self.level = 2**level, level
        self.shorten = self.N - self.N_ori
        assert self.shorten >= 0, \
            f'level({self.level}), N({self.N}) is not large enough for N_ori({self.N_ori})'

        self.words = (self.N + 63) // 64   # uint64 words per packed row

        # ── H_T and H (keep Python lists for __debug__ prints) ───────────────
        self.H_T = [m << self.shorten for m in self.H_ori_T]
        self.H   = [m << self.shorten for m in self.H_ori]
        for i in range(self.shorten - 1, -1, -1):
            self.H_T.append(1 << i)
            self.H.append(1 << i)

        if __debug__:
            print(f'{"H^T":=^{self.N}}')
            for i in self.H_T:
                print(f'{i:0{self.N}b}')

        if __debug__:
            print(f'{"H":=^{len(self.H_T)}}')
            for i in self.H:
                print(f'{i:0{len(self.H_T)}b}')

        self.bha_init = bha_init
        self.bha_init_array = [self.bha_init] * self.N_ori + [self.bha_init] * self.shorten  # TODO for shorten

        # ── default state ─────────────────────────────────────────────────────
        self.operation_array = "0" * (self.level * self.N // 2)
        self.permutation = list(range(self.N - 1, -1, -1))

        # ── numpy arrays for JIT ──────────────────────────────────────────────

        # index_pair → two int32 arrays
        na_list, nb_list = [], []
        for stage in range(self.level - 1, -1, -1):
            for index in range(self.N):
                if index & (1 << stage):
                    na_list.append(index ^ (1 << stage))
                    nb_list.append(index)
        self.na_arr = np.array(na_list, dtype=np.int32)
        self.nb_arr = np.array(nb_list, dtype=np.int32)

        # bha_init as float64 (matches bha_init_array)
        self.bha_init_np = np.array(self.bha_init_array, dtype=np.float64)

        # bha_init uniform for v2 (matches [self.bha_init] * self.N)
        self.bha_init_uniform_np = np.full(self.N, self.bha_init, dtype=np.float64)

        # H_T column indices precomputed (GP_T index = N - lsb.bit_length())
        H_rows  = len(self.H_T)
        max_nnz = max(bin(h).count('1') for h in self.H_T)
        self.H_T_cols = np.zeros((H_rows, max_nnz), dtype=np.int32)
        self.H_T_nnz  = np.zeros(H_rows, dtype=np.int32)
        for i, h in enumerate(self.H_T):
            j = 0
            temp = h
            while temp:
                lsb = temp & -temp
                self.H_T_cols[i, j] = self.N - lsb.bit_length()
                j += 1
                temp ^= lsb
            self.H_T_nnz[i] = j

        # H matrix packed (for set_system_parameter_v2)
        self.H_np = np.zeros((self.N, self.words), dtype=np.uint64)
        for i, h in enumerate(self.H):
            for w in range(self.words):
                self.H_np[i, w] = np.uint64((h >> (w * 64)) & 0xFFFFFFFFFFFFFFFF)

        self.kernel_col_idx = None
        self.bha = None
        self.bha_info = None

        # warm-up JIT compilation
        self.set_system_parameter_v1(self.permutation, self.operation_array)

    # ── hot path v1 ───────────────────────────────────────────────────────────

    def set_system_parameter_v1(self, permutation: list[int], operation_array: str):
        # compute ( G * P * H )^T

        # record
        self.permutation = permutation
        self.operation_array = operation_array

        if __debug__:
            print("permutation:", permutation)

        perm = np.asarray(permutation, dtype=np.int32)
        ops  = (np.frombuffer(operation_array.encode(), dtype=np.uint8) - ord('0')).astype(np.uint8)

        # ── kernel + bha ──────────────────────────────────────────────────────
        # kernel row i has initial bit at column N-1-i, which _build_GP_T maps
        # to channel permutation[N-1-i].  Permute bha_init accordingly so each
        # kernel row starts with the BHA value of its actual output channel.
        bha_init_permuted = self.bha_init_np[perm[::-1]]
        kernel, bha = _compute_kernel_bha(
            self.N, self.words, bha_init_permuted, ops, self.na_arr, self.nb_arr)

        if __debug__:
            N = self.N
            kernel_ints = _packed_mat_to_int_list(kernel)
            print(f'{"kernel":=^{N}}  {"bha"}')
            for k, b in zip(kernel_ints, bha.tolist()):
                print(f'{k:0{N}b}  {b:8f}')

        # ── GP_T ──────────────────────────────────────────────────────────────
        GP_T = _build_GP_T(kernel, perm, self.N, self.words)

        if __debug__:
            N = self.N
            GP_T_ints = _packed_mat_to_int_list(GP_T)
            print(f'{"GP":=^{N}}')
            for i in binary_matrix_transprose(GP_T_ints, N):
                print(f'{i:0{N}b}')

        # ── GPH_T ─────────────────────────────────────────────────────────────
        H_rows = len(self.H_T)
        GPH_T = _build_GPH_T(GP_T, self.H_T_cols, self.H_T_nnz, H_rows, self.words)

        if __debug__:
            N = self.N
            print(f'{"H_T":=^{N}}')
            for i in self.H_T:
                print(f'{i:0{N}b}')

            print(f'{"GPH_T":=^{N}}')
            for i in _packed_mat_to_int_list(GPH_T):
                print(f'{i:0{N}b}')

        # ── Gaussian Elimination ──────────────────────────────────────────────
        GPH_T, frozen_np = _gaussian_elimination_v1(GPH_T, self.N, self.words)

        frozen = _packed_to_int(frozen_np)
        N = self.N
        self.bha_info = [
            (idx, float(bha[idx]))
            for idx in range(N)
            if not ((frozen >> (N - 1 - idx)) & 1)
        ]

        if __debug__:
            print(f'{"GE GPH_T":=^{N}}')
            for i in _packed_mat_to_int_list(GPH_T):
                print(f'{i:0{N}b}')

            print(f'{"info":=^{N}}')
            print(f'{frozen:0{N}b}')
            for i, v in self.bha_info:
                print(f'{i:3d} : {v:8f}')
            print(f'sum bha_info: {sum([i[1] for i in self.bha_info])}')

    # ── hot path v2 ───────────────────────────────────────────────────────────

    def set_system_parameter_v2(self, permutation: list[int], operation_array: str):
        # record
        self.permutation = permutation
        self.operation_array = operation_array

        # compute ( G * P * H )^T
        N    = self.N
        ops  = (np.frombuffer(operation_array.encode(), dtype=np.uint8) - ord('0')).astype(np.uint8)
        perm = np.asarray(permutation, dtype=np.int32)

        # ── PH + bha ──────────────────────────────────────────────────────────
        PH, bha = _compute_PH_bha_v2(
            self.H_np, perm, N, self.words,
            self.bha_init_uniform_np, ops, self.na_arr, self.nb_arr)

        if __debug__:
            PH_ints = _packed_mat_to_int_list(PH)
            print(f'{"PH":=^{len(self.H_T)}}')
            for i in PH_ints:
                print(f'{i:0{len(self.H_T)}b}')

        # ── transpose PH → GPH_T ──────────────────────────────────────────────
        GPH_T = _packed_transpose_nxn(PH, N, self.words)

        if __debug__:
            print(f'{"GPH_T":=^{N}}')
            for i in _packed_mat_to_int_list(GPH_T):
                print(f'{i:0{N}b}')

        # ── Gaussian Elimination ──────────────────────────────────────────────
        GPH_T, frozen_np = _gaussian_elimination_v1(GPH_T, N, self.words)

        frozen = _packed_to_int(frozen_np)
        self.bha_info = [
            (idx, float(bha[idx]))
            for idx in range(N)
            if not ((frozen >> (N - 1 - idx)) & 1)
        ]

        if __debug__:
            print(f'{"GE GPH_T":=^{N}}')
            for i in _packed_mat_to_int_list(GPH_T):
                print(f'{i:0{N}b}')

            print(f'{"info":=^{N}}')
            print(f'{frozen:0{N}b}')
            for i, v in self.bha_info:
                print(f'{i:3d} : {v:8f}')

    # ── kept as Python method for API compatibility ───────────────────────────

    def Gaussian_Elimination(self, matrix):
        H_rows = len(matrix)
        matrix_np = np.zeros((H_rows, self.words), dtype=np.uint64)
        for i, val in enumerate(matrix):
            for w in range(self.words):
                matrix_np[i, w] = np.uint64((val >> (w * 64)) & 0xFFFFFFFFFFFFFFFF)

        result_np, frozen_np = _gaussian_elimination_v1(matrix_np, self.N, self.words)

        result = _packed_mat_to_int_list(result_np)
        lead_bit = _packed_to_int(frozen_np)
        return result, lead_bit


if __name__ == '__main__':
    matrix_path = "matrix/eRS_64_36_HT.matrix"
    S = PolarDecodeSystem.from_H_matrix_file(level=6, matrix_path=matrix_path)
    permutation = [34, 53, 25, 13, 31, 4, 11, 6, 63, 29, 35, 50, 61, 38, 3, 42, 24, 5, 36, 17, 57, 26, 8, 56, 33, 19, 23, 0, 39, 22, 43, 28, 21, 52, 58, 2, 15, 49, 45, 59, 62, 46, 10, 47, 32, 60, 1, 40, 14, 54, 30, 51, 9, 55, 44, 37, 12, 48, 41, 27, 16, 18, 7, 20]
    P = "011111001100111011010110010111110101111011011111111111010110001011111111111111111000110011111111111101110000000010001000100101001111111100000000000000000000000100000000100001000000000000000100"
    def o2n(P):
        return "".join([P[i:i+S.N//2] for i in range(0, len(P), S.N//2)][::-1])


    P = o2n(P)
    permutation = permutation[::-1]


    # ori
    global_best_op = '111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111110000000000000000000000000000000000000000000000000000000000000000'
    global_best_perm = [51, 50, 49, 48, 47, 46, 45, 44, 55, 54, 53, 52, 27, 26, 25, 24, 31, 30, 29, 28, 39, 38, 37, 36, 59, 58, 57, 56, 15, 14, 13, 12, 43, 42, 41, 40, 23, 22, 21, 20, 35, 34, 33, 32, 11, 10, 9, 8, 19, 18, 17, 16, 7, 6, 5, 4, 3, 2, 1, 0, 63, 62, 61, 60]

    # best -2
    global_best_op = "111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111101111111111111111111111111111111011110100110011101101111111000011011111101011000111110101011101"
    global_best_perm = [51, 46, 49, 28, 7, 50, 45, 36, 55, 54, 57, 32, 63, 58, 9, 40, 19, 30, 17, 12, 23, 38, 5, 56, 3, 14, 53, 20, 11, 10, 13, 8, 43, 22, 41, 16, 39, 42, 21, 4, 35, 34, 1, 52, 15, 2, 25, 48, 31, 18, 29, 60, 47, 6, 37, 0, 59, 62, 33, 44, 27, 26, 61, 24]

    # best -3
    # global_best_op = "111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111001011111111111111111111111111110111001101110111111101101111011000100101011101100010011011111111"
    # global_best_perm = [51, 42, 29, 28, 39, 22, 45, 4, 35, 14, 53, 32, 11, 54, 9, 60, 31, 18, 49, 40, 47, 38, 37, 44, 3, 26, 57, 56, 63, 2, 61, 24, 43, 50, 17, 16, 7, 46, 21, 36, 55, 62, 33, 52, 27, 34, 25, 12, 19, 30, 41, 48, 23, 6, 5, 20, 59, 10, 1, 0, 15, 58, 13, 8]


    permutation = permutation[::-1]
    def target_block(loop):
        for i in range(loop):
            # S.set_system_parameter_v1(permutation, P)
            S.set_system_parameter_v1(global_best_perm, global_best_op)
    
    if __debug__:
        target_block(1)
        x = sorted([i[1] for i in S.bha_info])
        x = [sum(x[:len(x)-i]) for i in range(len(x))]
        print(x)
        exit(0)

    if 0:
        from line_profiler import LineProfiler
        lp = LineProfiler()
        lp.add_function(PolarDecodeSystem.set_system_parameter_v1)
        lp.run('target_block(1000)')
        lp.print_stats()
    if 1:
        import cProfile
        cProfile.run("target_block(1000)")
    # target_block(1000)
