# Generate the permutation matrix from the provided permutation
import numpy as np

perm = [
31,26,14,3,51,54,61,7,17,9,53,29,63,22,55,15,
35,27,30,19,13,21,62,23,25,57,58,59,28,10,1,6,
5,49,60,33,11,37,38,39,41,52,42,43,44,45,46,47,
48,12,50,18,20,40,8,32,16,56,24,0,36,4,34,2
]

n = len(perm)
P = np.zeros((n, n), dtype=int)
for i, j in enumerate(perm):
    P[i, j] = 1

# Save to txt
path = "eBCH_m6_t11_perm_fw.txt"
with open(path, "w") as f:
    for row in P:
        f.write(" ".join(map(str, row)) + "\n")
