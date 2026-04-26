#!/usr/bin/env python3
import argparse
import math
import numpy as np


# 3GPP TS 38.212 Table 5.3.1.2-1.
# Q_0^{Nmax-1}, Nmax = 1024.
# Order: ascending reliability, i.e. least reliable -> most reliable.
NR_5G_POLAR_SEQUENCE_1024 = [
    0, 1, 2, 4, 8, 16, 32, 3, 5, 64, 9, 6, 17, 10, 18, 128,
    12, 33, 65, 20, 256, 34, 24, 36, 7, 129, 66, 512, 11, 40, 68, 130,
    19, 13, 48, 14, 72, 257, 21, 132, 35, 258, 26, 513, 80, 37, 25, 22,
    136, 260, 264, 38, 514, 96, 67, 41, 144, 28, 69, 42, 516, 49, 74, 272,
    160, 520, 288, 528, 192, 544, 70, 44, 131, 81, 50, 73, 15, 320, 133, 52,
    23, 134, 384, 76, 137, 82, 56, 27, 97, 39, 259, 84, 138, 145, 261, 29,
    43, 98, 515, 88, 140, 30, 146, 71, 262, 265, 161, 576, 45, 100, 640, 51,
    148, 46, 75, 266, 273, 517, 104, 162, 53, 193, 152, 77, 164, 768, 268,
    274, 518, 54, 83, 57, 521, 112, 135, 78, 289, 194, 85, 276, 522, 58,
    168, 139, 99, 86, 60, 280, 89, 290, 529, 524, 196, 141, 101, 147, 176,
    142, 530, 321, 31, 200, 90, 545, 292, 322, 532, 263, 149, 102, 105, 304,
    296, 163, 92, 47, 267, 385, 546, 324, 208, 386, 150, 153, 165, 106, 55,
    328, 536, 577, 548, 113, 154, 79, 269, 108, 578, 224, 166, 519, 552, 195,
    270, 641, 523, 275, 580, 291, 59, 169, 560, 114, 277, 156, 87, 197, 116,
    170, 61, 531, 525, 642, 281, 278, 526, 177, 293, 388, 91, 584, 769, 198,
    172, 120, 201, 336, 62, 282, 143, 103, 178, 294, 93, 644, 202, 592, 323,
    392, 297, 770, 107, 180, 151, 209, 284, 648, 94, 204, 298, 400, 608, 352,
    325, 533, 155, 210, 305, 547, 300, 109, 184, 534, 537, 115, 167, 225, 326,
    306, 772, 157, 656, 329, 110, 117, 212, 171, 776, 330, 226, 549, 538, 387,
    308, 216, 416, 271, 279, 158, 337, 550, 672, 118, 332, 579, 540, 389, 173,
    121, 553, 199, 784, 179, 228, 338, 312, 704, 390, 174, 554, 581, 393, 283,
    122, 448, 353, 561, 203, 63, 340, 394, 527, 582, 556, 181, 295, 285, 232,
    124, 205, 182, 643, 562, 286, 585, 299, 354, 211, 401, 185, 396, 344, 586,
    645, 593, 535, 240, 206, 95, 327, 564, 800, 402, 356, 307, 301, 417, 213,
    568, 832, 588, 186, 646, 404, 227, 896, 594, 418, 302, 649, 771, 360, 539,
    111, 331, 214, 309, 188, 449, 217, 408, 609, 596, 551, 650, 229, 159, 420,
    310, 541, 773, 610, 657, 333, 119, 600, 339, 218, 368, 652, 230, 391, 313,
    450, 542, 334, 233, 555, 774, 175, 123, 658, 612, 341, 777, 220, 314, 424,
    395, 673, 583, 355, 287, 183, 234, 125, 557, 660, 616, 342, 316, 241, 778,
    563, 345, 452, 397, 403, 207, 674, 558, 785, 432, 357, 187, 236, 664, 624,
    587, 780, 705, 126, 242, 565, 398, 346, 456, 358, 405, 303, 569, 244, 595,
    189, 566, 676, 361, 706, 589, 215, 786, 647, 348, 419, 406, 464, 680, 801,
    362, 590, 409, 570, 788, 597, 572, 219, 311, 708, 598, 601, 651, 421, 792,
    802, 611, 602, 410, 231, 688, 653, 248, 369, 190, 364, 654, 659, 335, 480,
    315, 221, 370, 613, 422, 425, 451, 614, 543, 235, 412, 343, 372, 775, 317,
    222, 426, 453, 237, 559, 833, 804, 712, 834, 661, 808, 779, 617, 604, 433,
    720, 816, 836, 347, 897, 243, 662, 454, 318, 675, 618, 898, 781, 376, 428,
    665, 736, 567, 840, 625, 238, 359, 457, 399, 787, 591, 678, 434, 677, 349,
    245, 458, 666, 620, 363, 127, 191, 782, 407, 436, 626, 571, 465, 681, 246,
    707, 350, 599, 668, 790, 460, 249, 682, 573, 411, 803, 789, 709, 365, 440,
    628, 689, 374, 423, 466, 793, 250, 371, 481, 574, 413, 603, 366, 468, 655,
    900, 805, 615, 684, 710, 429, 794, 252, 373, 605, 848, 690, 713, 632, 482,
    806, 427, 904, 414, 223, 663, 692, 835, 619, 472, 455, 796, 809, 714, 721,
    837, 716, 864, 810, 606, 912, 722, 696, 377, 435, 817, 319, 621, 812, 484,
    430, 838, 667, 488, 239, 378, 459, 622, 627, 437, 380, 818, 461, 496, 669,
    679, 724, 841, 629, 351, 467, 438, 737, 251, 462, 442, 441, 469, 247, 683,
    842, 738, 899, 670, 783, 849, 820, 728, 928, 791, 367, 901, 630, 685, 844,
    633, 711, 253, 691, 824, 902, 686, 740, 850, 375, 444, 470, 483, 415, 485,
    905, 795, 473, 634, 744, 852, 960, 865, 693, 797, 906, 715, 807, 474, 636,
    694, 254, 717, 575, 913, 798, 811, 379, 697, 431, 607, 489, 866, 723, 486,
    908, 718, 813, 476, 856, 839, 725, 698, 914, 752, 868, 819, 814, 439, 929,
    490, 623, 671, 739, 916, 463, 843, 381, 497, 930, 821, 726, 961, 872, 492,
    631, 729, 700, 443, 741, 845, 920, 382, 822, 851, 730, 498, 880, 742, 445,
    471, 635, 932, 687, 903, 825, 500, 846, 745, 826, 732, 446, 962, 936, 475,
    853, 867, 637, 907, 487, 695, 746, 828, 753, 854, 857, 504, 799, 255, 964,
    909, 719, 477, 915, 638, 748, 944, 869, 491, 699, 754, 858, 478, 968, 383,
    910, 815, 976, 870, 917, 727, 493, 873, 701, 931, 756, 860, 499, 731, 823,
    922, 874, 918, 502, 933, 743, 760, 881, 494, 702, 921, 501, 876, 847, 992,
    447, 733, 827, 934, 882, 937, 963, 747, 505, 855, 924, 734, 829, 965, 938,
    884, 506, 749, 945, 966, 755, 859, 940, 830, 911, 871, 639, 888, 479, 946,
    750, 969, 508, 861, 757, 970, 919, 875, 862, 758, 948, 977, 923, 972, 761,
    877, 952, 495, 703, 935, 978, 883, 762, 503, 925, 878, 735, 993, 885, 939,
    994, 980, 926, 764, 941, 967, 886, 831, 947, 507, 889, 984, 751, 942, 996,
    971, 890, 509, 949, 973, 1000, 892, 950, 863, 759, 1008, 510, 979, 953,
    763, 974, 954, 879, 981, 982, 927, 995, 765, 956, 887, 985, 997, 986, 943,
    891, 998, 766, 511, 988, 1001, 951, 1002, 893, 975, 894, 1009, 955, 1004,
    1010, 957, 983, 958, 987, 1012, 999, 1016, 767, 989, 1003, 990, 1005, 959,
    1011, 1013, 895, 1006, 1014, 1017, 1018, 991, 1020, 1007, 1015, 1019,
    1021, 1022, 1023
]


def gf2_matmul(A, B):
    return (A @ B) % 2


def kronecker_power_F(m):
    F = np.array([[1, 0], [1, 1]], dtype=np.uint8)
    G = np.array([[1]], dtype=np.uint8)
    for _ in range(m):
        G = np.kron(G, F) % 2
    return G.astype(np.uint8)


def bit_reverse(i, m):
    """
    Reverse the m-bit binary representation of i.

    Example:
        m=4, i=1=0001 -> 1000=8.
    """
    rev = 0
    for _ in range(m):
        rev = (rev << 1) | (i & 1)
        i >>= 1
    return rev


def bit_reversal_matrix(m):
    """
    Return B_m with the convention

        (B_m @ A)[i, :] = A[bit_reverse(i), :].

    Hence Gp = B_m F^{otimes m} when include_bitreversal=True.
    """
    n = 1 << m
    B = np.zeros((n, n), dtype=np.uint8)
    for i in range(n):
        B[i, bit_reverse(i, m)] = 1
    return B


def build_polar_transform(n, include_bitreversal=False, check_involutory=True):
    """
    Build the polar transform Gp used in the experiment.

    include_bitreversal=False:
        Gp = F^{otimes m}

    include_bitreversal=True:
        Gp = B_m F^{otimes m}

    For both conventions here, Gp is involutory over GF(2), so Gp^{-1}=Gp.
    The theoretical p,q routine below relies on this involutory property.
    """
    if n <= 0 or (n & (n - 1)) != 0:
        raise ValueError("n must be a power of two.")

    m = int(math.log2(n))
    Fm = kronecker_power_F(m)

    if include_bitreversal:
        Bm = bit_reversal_matrix(m)
        Gp = gf2_matmul(Bm, Fm)
    else:
        Gp = Fm

    if check_involutory:
        identity = np.eye(n, dtype=np.uint8)
        if not np.array_equal(gf2_matmul(Gp, Gp), identity):
            raise RuntimeError("Expected Gp^{-1}=Gp, but Gp @ Gp != I.")

    return Gp


def selector_matrix(I, n):
    S = np.zeros((len(I), n), dtype=np.uint8)
    for row, idx in enumerate(I):
        S[row, idx] = 1
    return S


def transposition_matrix(n, a, b):
    """
    Matrix P_(a,b) for right action on row vectors: c -> cP.
    """
    P = np.eye(n, dtype=np.uint8)
    P[:, [a, b]] = P[:, [b, a]]
    return P


def gf2_rref(A):
    R = A.copy().astype(np.uint8)
    m, n = R.shape
    pivots = []
    row = 0

    for col in range(n):
        pivot = None
        for r in range(row, m):
            if R[r, col]:
                pivot = r
                break

        if pivot is None:
            continue

        if pivot != row:
            R[[row, pivot]] = R[[pivot, row]]

        for r in range(m):
            if r != row and R[r, col]:
                R[r, :] ^= R[row, :]

        pivots.append(col)
        row += 1

        if row == m:
            break

    return R, pivots


def support(vec):
    return [i for i, x in enumerate(vec.tolist()) if x % 2 == 1]


def pw_reliability_order(n, beta=2 ** 0.25):
    """
    Non-3GPP fallback for n up to 2048.
    Order returned is least reliable -> most reliable under polarization weight.
    Larger PW score means more reliable.
    """
    m = int(math.log2(n))
    order = []
    for i in range(n):
        score = 0.0
        for j in range(m):
            if (i >> j) & 1:
                score += beta ** j
        order.append((score, i))
    order.sort(key=lambda x: (x[0], x[1]))
    return [i for _, i in order]


def reliability_order(n, ranking):
    """
    Return reliability order from least reliable to most reliable.
    """
    if ranking == "5g":
        if n > 1024:
            raise ValueError(
                "3GPP 5G NR polar reliability sequence is defined only up to Nmax=1024. "
                "Use --ranking pw for a non-standard n>1024 fallback."
            )
        order = [i for i in NR_5G_POLAR_SEQUENCE_1024 if i < n]
        if len(order) != n:
            raise RuntimeError("5G reliability order length mismatch.")
        return order

    if ranking == "pw":
        if n > 2048:
            raise ValueError("PW fallback in this script is limited to n<=2048.")
        return pw_reliability_order(n)

    raise ValueError(f"Unknown ranking: {ranking}")


def default_info_set(n, k, ranking="5g"):
    """
    Choose the best k channels.

    reliability_order returns least -> most reliable,
    so the information set is the last k entries, sorted increasingly
    for the row-selector S_I.
    """
    order = reliability_order(n, ranking)
    if len(order) != n:
        raise RuntimeError("Reliability order length mismatch.")
    return sorted(order[-k:])


def normalize_info_set(n, k, I=None, ranking="5g"):
    if I is None:
        return default_info_set(n, k, ranking=ranking)

    I = sorted(I)
    if len(I) != k:
        raise ValueError(f"Given I has size {len(I)}, but k={k}.")
    if len(set(I)) != len(I):
        raise ValueError("I contains duplicated indices.")
    if any(i < 0 or i >= n for i in I):
        raise ValueError("I must be a subset of [0,n-1].")
    return I


def compute_theoretical_pq(n, I, Gp, a, b):
    """
    Compute theoretical T, q, beta, residual, p for transposition (a,b).

    This implementation intentionally uses the actual current Gp, not UPO
    support shortcuts. Therefore it works for both:
        Gp = F^{otimes m}
        Gp = B_m F^{otimes m}

    Since the chosen Gp is involutory:
        Gp^{-1} = Gp.

    Hence:
        v   = S_I Gp(e_a + e_b), determined by c_a + c_b.
        y^T = (e_a + e_b)^T Gp^{-1}, determined by row_a(Gp)+row_b(Gp).
    """
    I = sorted(I)
    I_set = set(I)

    c_a = Gp[:, a]
    c_b = Gp[:, b]
    col_diff = (c_a ^ c_b).astype(np.uint8)

    y = (Gp[a, :] ^ Gp[b, :]).astype(np.uint8)

    # T is in coordinate labels, not row numbers.
    T = [i for i in I if col_diff[i] == 1]

    if len(T) == 0:
        return {
            "T": [],
            "q": None,
            "beta": None,
            "residual": None,
            "p": None,
            "predicted_pivots": I.copy(),
            "case": "unchanged because T is empty",
            "col_diff_support": support(col_diff),
            "y_support": support(y),
        }

    q = max(T)

    beta = 1
    for i in T:
        beta ^= int(y[i])

    residual = np.zeros(n, dtype=np.uint8)

    # y_{I^c}
    for j in range(n):
        if j not in I_set:
            residual[j] = y[j]

    # + beta e_q
    if beta == 1:
        residual[q] ^= 1

    residual_support = support(residual)
    if not residual_support:
        raise RuntimeError(
            "Theoretical residual is zero. This should not happen because W has rank k."
        )

    p = min(residual_support)
    predicted = sorted((I_set - {q}) | {p})

    return {
        "T": T,
        "q": q,
        "beta": beta,
        "residual": residual,
        "p": p,
        "predicted_pivots": predicted,
        "case": "one-swap formula",
        "col_diff_support": support(col_diff),
        "y_support": support(y),
    }


def direct_rref_pivots_for_transposition(n, I, Gp, a, b):
    """
    Directly compute W_(a,b)=S_I Gp P_(a,b) Gp^{-1} and its RREF pivots.

    The current Gp is involutory, so Gp^{-1}=Gp.
    """
    S = selector_matrix(I, n)
    P = transposition_matrix(n, a, b)

    W = gf2_matmul(gf2_matmul(gf2_matmul(S, Gp), P), Gp)
    R, pivots = gf2_rref(W)

    return W, R, pivots


def verify(n, k, a, b, I=None, ranking="5g", include_bitreversal=False, verbose=True):
    if n <= 0 or (n & (n - 1)) != 0:
        raise ValueError("n must be a power of two.")
    if n > 2048:
        raise ValueError("This script supports n up to 2048.")
    if not (0 < k <= n):
        raise ValueError("k must satisfy 0 < k <= n.")
    if a == b:
        raise ValueError("a and b must be distinct.")
    if not (0 <= a < n and 0 <= b < n):
        raise ValueError("a,b must be in [0,n-1].")

    I = normalize_info_set(n, k, I=I, ranking=ranking)
    Gp = build_polar_transform(n, include_bitreversal=include_bitreversal)

    W, R, actual_pivots = direct_rref_pivots_for_transposition(n, I, Gp, a, b)

    theory = compute_theoretical_pq(n, I, Gp, a, b)
    predicted_pivots = theory["predicted_pivots"]
    ok = actual_pivots == predicted_pivots

    if verbose:
        print("=" * 80)
        print(f"n = {n}, k = {k}, ranking = {ranking}")
        print(f"include_bitreversal = {include_bitreversal}")
        print(f"I = {I}")
        print(f"transposition (a,b) = ({a},{b})")
        print("-" * 80)
        print(f"supp(c_a+c_b) = {theory['col_diff_support']}")
        print(f"supp(y) = supp(row_a+row_b) = {theory['y_support']}")
        print(f"T = {theory['T']}")
        print(f"q = {theory['q']}")
        print(f"beta = {theory['beta']}")
        if theory["residual"] is not None:
            print(f"supp(residual r) = {support(theory['residual'])}")
        else:
            print("residual r = None")
        print(f"p = {theory['p']}")
        print("-" * 80)
        print(f"Predicted pivots = {predicted_pivots}")
        print(f"Actual RREF pivots = {actual_pivots}")
        print(f"Match? {ok}")
        print("-" * 80)
        print("RREF(W_pi) =")
        print(R.astype(int))
        print("=" * 80)

    return {
        "ok": ok,
        "I": I,
        "W": W,
        "RREF": R,
        "actual_pivots": actual_pivots,
        "theory": theory,
        "include_bitreversal": include_bitreversal,
    }


def first_frozen_leakage_p(n, I, Gp, a, b):
    """
    Candidate-level p from y_{I^c} only.

    This ignores T, q, beta.
    It only checks:
        p0 = min supp(y_{I^c})
    where y^T = row_a(Gp) + row_b(Gp).

    This uses the actual current Gp, either F^{otimes m}
    or B_m F^{otimes m}; no UPO support shortcut is assumed.
    """
    I_set = set(I)

    y = (Gp[a, :] ^ Gp[b, :]).astype(np.uint8)
    frozen_support = [j for j in range(n) if j not in I_set and y[j] == 1]

    if len(frozen_support) == 0:
        return None

    return min(frozen_support)


def verify_candidate(n, I, Gp, a, b, target_p=None, target_q=None):
    """
    Exact verifier for one candidate (a,b).

    It checks:
      1. theoretical q,p from compute_theoretical_pq
      2. direct RREF pivots
      3. whether the candidate matches desired target_p / target_q
    """
    theory = compute_theoretical_pq(n, I, Gp, a, b)
    _, R, actual_pivots = direct_rref_pivots_for_transposition(n, I, Gp, a, b)

    predicted_pivots = theory["predicted_pivots"]
    formula_matches_rref = (actual_pivots == predicted_pivots)

    p_ok = (target_p is None) or (theory["p"] == target_p)
    q_ok = (target_q is None) or (theory["q"] == target_q)

    return {
        "a": a,
        "b": b,
        "T": theory["T"],
        "q": theory["q"],
        "beta": theory["beta"],
        "residual_support": support(theory["residual"]) if theory["residual"] is not None else None,
        "p": theory["p"],
        "predicted_pivots": predicted_pivots,
        "actual_pivots": actual_pivots,
        "formula_matches_rref": formula_matches_rref,
        "target_p_ok": p_ok,
        "target_q_ok": q_ok,
        "ok": formula_matches_rref and p_ok and q_ok,
        "RREF": R,
        "col_diff_support": theory["col_diff_support"],
        "y_support": theory["y_support"],
    }


def generate_candidates_for_p(
    n,
    k,
    target_p,
    target_q=None,
    I=None,
    ranking="5g",
    include_bitreversal=False,
    verify_rref=True,
    max_print=None,
):
    """
    Enumerate all transpositions (a,b) and find candidates for desired incoming pivot p.

    There are two levels:

    weak_candidates:
        These satisfy only
            target_p = min supp(y_{I^c})
        i.e. p is the first frozen leakage of y.
        This does NOT guarantee that final RREF incoming pivot is target_p.

    exact_candidates:
        These satisfy the full theoretical formula:
            theory["p"] == target_p
        and optionally
            theory["q"] == target_q.

    verified_candidates:
        These are exact_candidates also checked by direct RREF computation.

    The generator uses actual c_a+c_b and row_a+row_b from the chosen Gp.
    Therefore it works with both include_bitreversal=False and True.
    """
    if n <= 0 or (n & (n - 1)) != 0:
        raise ValueError("n must be a power of two.")
    if n > 2048:
        raise ValueError("This script supports n up to 2048.")
    if not (0 < k <= n):
        raise ValueError("k must satisfy 0 < k <= n.")
    if not (0 <= target_p < n):
        raise ValueError("target_p must be in [0,n-1].")

    I = normalize_info_set(n, k, I=I, ranking=ranking)
    I_set = set(I)

    if target_p in I_set:
        raise ValueError(
            "target_p is already in I. For an incoming pivot, target_p should be in I^c."
        )

    Gp = build_polar_transform(n, include_bitreversal=include_bitreversal)

    weak_candidates = []
    exact_candidates = []
    verified_candidates = []

    for a in range(n):
        for b in range(a + 1, n):
            # Weak filter: only check the first frozen leakage of actual y.
            p0 = first_frozen_leakage_p(n, I, Gp, a, b)

            if p0 == target_p:
                weak_candidates.append((a, b))

            # Full theoretical p/q computation.
            theory = compute_theoretical_pq(n, I, Gp, a, b)

            if theory["p"] != target_p:
                continue

            if target_q is not None and theory["q"] != target_q:
                continue

            exact_candidates.append((a, b, theory))

            if verify_rref:
                checked = verify_candidate(
                    n=n,
                    I=I,
                    Gp=Gp,
                    a=a,
                    b=b,
                    target_p=target_p,
                    target_q=target_q,
                )

                if checked["ok"]:
                    verified_candidates.append(checked)
                else:
                    print("WARNING: formula/RREF mismatch or target mismatch")
                    print(checked)
            else:
                verified_candidates.append({
                    "a": a,
                    "b": b,
                    "T": theory["T"],
                    "q": theory["q"],
                    "beta": theory["beta"],
                    "residual_support": support(theory["residual"]) if theory["residual"] is not None else None,
                    "p": theory["p"],
                    "predicted_pivots": theory["predicted_pivots"],
                    "col_diff_support": theory["col_diff_support"],
                    "y_support": theory["y_support"],
                })

    print("=" * 80)
    print(f"n={n}, k={k}, ranking={ranking}")
    print(f"include_bitreversal={include_bitreversal}")
    print(f"I={I}")
    print(f"target_p={target_p}, target_q={target_q}")
    print("-" * 80)
    print(f"weak candidates from first frozen leakage only: {len(weak_candidates)}")
    print(f"exact theoretical candidates: {len(exact_candidates)}")
    print(f"RREF-verified candidates: {len(verified_candidates)}")
    print("=" * 80)

    if max_print is not None:
        for item in verified_candidates[:max_print]:
            print(
                f"(a,b)=({item['a']},{item['b']}), "
                f"T={item['T']}, q={item['q']}, beta={item['beta']}, "
                f"residual_support={item['residual_support']}, p={item['p']}, "
                f"pivots={item.get('actual_pivots', item.get('predicted_pivots'))}"
            )

    return {
        "I": I,
        "target_p": target_p,
        "target_q": target_q,
        "include_bitreversal": include_bitreversal,
        "weak_candidates": weak_candidates,
        "exact_candidates": exact_candidates,
        "verified_candidates": verified_candidates,
    }


def parse_info_set(s):
    if s is None or not s.strip():
        return None
    return sorted(int(x.strip()) for x in s.split(","))


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Verify theoretical p,q and exact RREF(W_pi), or inverse-search "
            "transpositions producing a desired incoming pivot p."
        )
    )
    parser.add_argument(
        "--mode",
        choices=["verify", "inverse"],
        default="verify",
        help="verify: check one transposition. inverse: search candidates for target p."
    )
    parser.add_argument("--n", type=int, required=True, help="Code length n=2^m.")
    parser.add_argument("--k", type=int, required=True, help="Dimension k.")

    parser.add_argument("--a", type=int, default=None, help="First transposition index.")
    parser.add_argument("--b", type=int, default=None, help="Second transposition index.")

    parser.add_argument(
        "--target-p",
        type=int,
        default=None,
        help="Desired incoming pivot p for inverse mode."
    )
    parser.add_argument(
        "--target-q",
        type=int,
        default=None,
        help="Optional desired outgoing pivot q for inverse mode."
    )

    parser.add_argument(
        "--I",
        type=str,
        default=None,
        help="Optional comma-separated information set. If omitted, choose best k channels."
    )
    parser.add_argument(
        "--ranking",
        choices=["5g", "pw"],
        default="5g",
        help=(
            "Reliability ranking. '5g' uses 3GPP NR sequence for n<=1024. "
            "'pw' is a non-standard polarization-weight fallback for n<=2048."
        )
    )
    parser.add_argument(
        "--bit-reversal",
        action="store_true",
        help="Use Gp = B_m F^{otimes m}. If omitted, use Gp = F^{otimes m}."
    )
    parser.add_argument(
        "--no-rref-verify",
        action="store_true",
        help="In inverse mode, skip direct RREF verification and only use the formula."
    )
    parser.add_argument(
        "--max-print",
        type=int,
        default=20,
        help="Maximum number of inverse candidates to print."
    )

    args = parser.parse_args()
    I = parse_info_set(args.I)

    if args.mode == "verify":
        if args.a is None or args.b is None:
            raise ValueError("--a and --b are required in verify mode.")

        verify(
            args.n,
            args.k,
            args.a,
            args.b,
            I=I,
            ranking=args.ranking,
            include_bitreversal=args.bit_reversal,
            verbose=True,
        )

    elif args.mode == "inverse":
        if args.target_p is None:
            raise ValueError("--target-p is required in inverse mode.")

        generate_candidates_for_p(
            n=args.n,
            k=args.k,
            target_p=args.target_p,
            target_q=args.target_q,
            I=I,
            ranking=args.ranking,
            include_bitreversal=args.bit_reversal,
            verify_rref=(not args.no_rref_verify),
            max_print=args.max_print,
        )


if __name__ == "__main__":
    main()
