#!/usr/bin/env python3
"""
Generate a static polar-code generator matrix using the 5G NR reliability sequence.

Usage:
    python make_static_polar_matrix.py --n 128 --k 64

Output:
    static_polar_128_64.matrix

Assumed .matrix format:
    k n
    row_0 entries separated by spaces
    row_1 entries separated by spaces
    ...
"""

import argparse
from pathlib import Path


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


def bit_reverse(i: int, m: int) -> int:
    out = 0
    for _ in range(m):
        out = (out << 1) | (i & 1)
        i >>= 1
    return out


def is_power_of_two(n: int) -> bool:
    return n > 0 and (n & (n - 1)) == 0


def arikan_kernel_power(n: int) -> list[list[int]]:
    """
    Return F^{⊗m}, where

        F = [[1, 0],
             [1, 1]]

    over GF(2).
    """
    if not is_power_of_two(n):
        raise ValueError("--n must be a power of two.")

    G = [[1]]
    F = [[1, 0], [1, 1]]

    size = 1
    while size < n:
        new_G = [[0 for _ in range(2 * size)] for _ in range(2 * size)]

        for a in range(2):
            for b in range(2):
                if F[a][b] == 0:
                    continue
                for i in range(size):
                    for j in range(size):
                        new_G[a * size + i][b * size + j] = G[i][j]

        G = new_G
        size *= 2

    return G


def bit_reversal_permuted_arikan_generator(n: int) -> list[list[int]]:
    """
    Return

        G_N = B_N F^{⊗m}

    meaning row-bit-reversal followed by the Arikan kernel.
    """
    m = n.bit_length() - 1
    F_power = arikan_kernel_power(n)

    G = []
    for i in range(n):
        reversed_i = bit_reverse(i, m)
        G.append(F_power[reversed_i][:])

    return G


def bit_reversal_permuted_arikan_inverse(n: int) -> list[list[int]]:
    """
    Return

        G_N^{-1} = (B_N F^{⊗m})^{-1}
                 = F^{⊗m} B_N.

    Since F^{-1} = F over GF(2), F^{⊗m} is also self-inverse.
    Right multiplication by B_N means column-bit-reversal.
    """
    m = n.bit_length() - 1
    F_power = arikan_kernel_power(n)

    G_inv = [[0 for _ in range(n)] for _ in range(n)]

    for col in range(n):
        reversed_col = bit_reverse(col, m)
        for row in range(n):
            G_inv[row][col] = F_power[row][reversed_col]

    return G_inv


def get_5g_information_set(n: int, k: int) -> list[int]:
    """
    Select the k most reliable bit-channel indices for length n.

    The 5G NR reliability sequence is ordered from less reliable to more reliable.
    For length n, keep only indices < n, then take the last k entries.
    """
    if n > 1024:
        raise ValueError("The 5G NR sequence only supports n <= 1024.")
    if not is_power_of_two(n):
        raise ValueError("--n must be a power of two.")
    if not (0 < k <= n):
        raise ValueError("--k must satisfy 0 < k <= n.")

    q_n = [i for i in NR_5G_POLAR_SEQUENCE_1024 if i < n]

    if len(q_n) != n:
        raise ValueError(
            f"Reliability sequence filtering failed: expected {n} entries, got {len(q_n)}."
        )

    info_set = q_n[-k:]
    return sorted(info_set)


def get_frozen_set(n: int, info_set: list[int]) -> list[int]:
    info = set(info_set)
    return [i for i in range(n) if i not in info]


def make_static_polar_generator(n: int, info_set: list[int]) -> list[list[int]]:
    """
    Return

        G_static = G_N[I, :].
    """
    G_n = bit_reversal_permuted_arikan_generator(n)
    return [G_n[i] for i in info_set]


def make_static_polar_parity_check(n: int, frozen_set: list[int]) -> list[list[int]]:
    """
    Return parity-check matrix H for the static polar code.

    For frozen index f,

        (c G_N^{-1})_f = 0.

    Hence the corresponding parity-check row is

        h_f = (G_N^{-1})[:, f]^T.
    """
    G_inv = bit_reversal_permuted_arikan_inverse(n)

    H = []
    for f in frozen_set:
        row = [G_inv[i][f] for i in range(n)]
        H.append(row)

    H_T = [list(row) for row in zip(*H)]
    

    return H_T


def write_matrix_file(path: Path, matrix: list[list[int]]) -> None:
    """
    Assumed .matrix format:

        rows cols
        row_0
        row_1
        ...

    Each row is written as space-separated 0/1 entries.
    """
    rows = len(matrix)
    cols = len(matrix[0]) if rows > 0 else 0

    with path.open("w", encoding="utf-8") as f:
        f.write(f"{rows} {cols}\n")
        for row in matrix:
            f.write(" ".join(str(x & 1) for x in row))
            f.write("\n")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate static polar G and H matrices using the 5G NR reliability sequence."
    )
    parser.add_argument("--n", type=int, required=True, help="Polar block length, power of two, <= 1024.")
    parser.add_argument("--k", type=int, required=True, help="Dimension / number of information bits.")
    parser.add_argument(
        "--output-dir",
        type=str,
        default=".",
        help="Directory for output .matrix files.",
    )

    args = parser.parse_args()

    n = args.n
    k = args.k

    info_set = get_5g_information_set(n, k)
    frozen_set = get_frozen_set(n, info_set)

    G_static = make_static_polar_generator(n, info_set)
    H_static = make_static_polar_parity_check(n, frozen_set)

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    G_path = output_dir / f"static_polar_{n}_{k}.matrix"
    H_path = output_dir / f"static_polar_{n}_{k}_H.matrix"

    write_matrix_file(G_path, G_static)
    write_matrix_file(H_path, H_static)

    print(f"Wrote {G_path}")
    print(f"Wrote {H_path}")
    print(f"n = {n}, k = {k}")
    print(f"Information set I = {info_set}")
    print(f"Frozen set F = {frozen_set}")


if __name__ == "__main__":
    main()