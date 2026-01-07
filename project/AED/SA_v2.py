import random
import os,sys,math

STAGE = 7 
N, K = 2**STAGE, 99 
FILE_NAME = os.path.basename(sys.argv[1])[:-4]
NEIGHBOR_LIST = [1]*N + [2]*int(N*(N-1)/2) + [3]*int(N*(N-1)*(N-2)/6) + [4]*int(N*(N-1)*(N-2)*(N-3)/24)
NEIGHBOR_SIZE = len(NEIGHBOR_LIST)
NEIGHBOR_SIZE = 10000000
P_IDX = -8
T_INIT, T_ALPHA, T_COUNT = 0.8, 0.99, 2000
pSize = N
opSize = STAGE*N>>1

# best_gene = [[ 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1],
# [65, 15, 11, 162, 230, 188, 34, 138, 161, 85, 99, 224, 225, 221, 206, 143, 31, 81, 24, 106, 218, 227, 177, 255, 211, 4, 104, 210, 142, 168, 209, 214, 83, 75, 155, 158, 139, 251, 90, 237, 163, 38, 176, 97, 79, 197, 160, 171, 67, 91, 152, 131, 100, 107, 244, 222, 182, 50, 114, 102, 223, 133, 105, 249, 37, 5, 62, 178, 89, 250, 242, 42, 66, 217, 121, 215, 68, 126, 58, 45, 184, 19, 74, 167, 226, 39, 187, 29, 232, 112, 213, 86, 109, 196, 169, 72, 55, 32, 118, 192, 94, 115, 36, 98, 122, 84, 110, 146, 51, 43, 245, 35, 229, 2, 20, 233, 216, 48, 157, 200, 93, 59, 80, 27, 159, 247, 207, 236, 8, 203, 14, 186, 175, 220, 70, 241, 165, 49, 127, 60, 44, 123, 145, 69, 180, 253, 235, 30, 87, 179, 136, 147, 252, 174, 150, 195, 108, 199, 46, 191, 12, 47, 78, 181, 76, 243, 137, 113, 212, 52, 254, 132, 190, 151, 77, 173, 156, 154, 117, 234, 148, 13, 194, 22, 204, 166, 82, 116, 231, 238, 3, 135, 41, 88, 248, 228, 119, 202, 198, 153, 141, 71, 189, 25, 183, 111, 172, 6, 16, 201, 149, 56, 208, 40, 10, 125, 128, 7, 33, 101, 240, 96, 23, 124, 129, 1, 239, 205, 140, 170, 64, 144, 246, 134, 73, 219, 193, 54, 21, 92, 120, 26, 164, 17, 9, 103, 95, 63, 185, 53, 28, 130, 61, 0, 57, 18]]

with open("%s.ini"%FILE_NAME, 'r') as din:
    gold_ini = din.read().split('\n')[:-1]

def getOneStepNeighbor(gene):
    _gene = [[i for i in gene[0]], [i for i in gene[1]]]
    if random.randint(0, 1):
        # random op
        idx_x = random.randint(0, opSize-1)
        _gene[0][idx_x] ^= 1
    else:
        idx_a, idx_b = 0, 0
        while idx_a == idx_b:
            idx_a = random.randint(0, pSize-1)
            idx_b = random.randint(0, pSize-1)
        _gene[1][idx_a], _gene[1][idx_b] = _gene[1][idx_b], _gene[1][idx_a]
    return _gene

def genData(gene):
    # matrix
    with open("%s.matrix"%FILE_NAME, 'w') as dout:
        dout.write("%d %d\n"%(pSize, pSize))
        for i in range(pSize):
            dout.write("%s\n"%("".join(["0" if gene[1][i] != j else "1" for j in range(pSize)])))

    # ini
    with open("%s.ini"%FILE_NAME, 'w') as dout:
        for line in gold_ini:
            if line.startswith("operationArray"):
                dout.write("operationArray                 = " + "".join(["%d"%i for i in gene[0]]) + '\n')
            else:
                dout.write(line+'\n')

def geneValue(gene):
    genData(gene)
    os.system("./AdjustPolarDecoder_bha -ini %s.ini >log_%s.log"%(FILE_NAME,FILE_NAME))
    with open("log_%s.log"%FILE_NAME ,'r') as din:
        line = din.read().split('\n')

    bha_ = [float(l[6:]) for l in line[N:N+K]]
    bha_ = sorted(bha_)
    bha_sum = [sum(bha_[:l+1]) for l in range(K)]
    return bha_sum

exit_flag = 0
try:
    print(geneValue(best_gene))
    print("gen by best gene")
    exit_flag = 1
except:
    pass
if exit_flag:
    exit(0)

# while True:
for x in range(int(sys.argv[2],10)):

    # gene = (op string, permutation list)
    geneInit = [[random.randint(0, 1) for i in range(STAGE*N>>1)], [i for i in range(pSize)]]
    random.shuffle(geneInit[1])
    curentPair = (geneValue(geneInit)[P_IDX], geneInit)
    unchange_cnt = 0
    t, t_cnt = T_INIT, 0

    while unchange_cnt < (NEIGHBOR_SIZE):
        tmpGene = curentPair[1]
        step = random.choice(NEIGHBOR_LIST)
        for s_idx in range(step):
            tmpGene = getOneStepNeighbor(tmpGene)
        tmpGeneV = geneValue(tmpGene)[P_IDX]

        delta = tmpGeneV - curentPair[0]
        if delta > 0:
            if random.random() < math.exp(-delta/t):
                nextPair = (tmpGeneV, tmpGene)
            else:
                nextPair = curentPair
        else:
            nextPair = (tmpGeneV, tmpGene)

        if nextPair[0] == curentPair[0]:
            unchange_cnt = unchange_cnt+1
        else:
            curentPair = nextPair
            unchange_cnt = 0

        if t_cnt == T_COUNT:
            t_cnt = 0
            t = max(t*T_ALPHA, 0.00000000000000001)        
            print(t, nextPair[0])
            with open("generat_log_%s.log"%FILE_NAME ,'w') as dout:
                gene = curentPair[1]
                bha_sum = geneValue(gene)
                dout.write("%s\n"%(", ".join(["%f"%x for x in bha_sum])))
                dout.write("op %s\n"%(", ".join(["%d"%i for i in gene[0]])))
                dout.write("[%s]\n"%(", ".join(["%d"%i for i in gene[1]])))
        else:
            t_cnt = t_cnt + 1
        
    with open("final_log_%s.log"%FILE_NAME ,'a+t') as dout:
        gene = curentPair[1]
        bha_sum = geneValue(gene)
        dout.write("%s\n"%(", ".join(["%f"%x for x in bha_sum])))
        dout.write("op %s\n"%(", ".join(["%d"%i for i in gene[0]])))
        dout.write("[%s]\n"%(", ".join(["%d"%i for i in gene[1]])))
