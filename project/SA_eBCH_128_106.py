import os
import PolarDecoding_numpy as PD
import random, math, time, heapq


def target_block(run_dir: str = ".", seed: int = None, P_IDX: int = -2):

    os.makedirs(run_dir, exist_ok=True)
    watchdog_log = os.path.join(run_dir, "watchdog.log")
    result_log   = os.path.join(run_dir, "result.log")

    if seed is not None:
        random.seed(seed)

    matrix_path = "./POD/eBCH_m7_t3_H.matrix"
    # matrix_path = "../data/extend_BCH_128_106_H.matrix"
    S = PD.PolarDecodeSystem.from_H_matrix_file(level=7, matrix_path=matrix_path)

    NEIGHBOR_WEIGHTS = [S.N, int(S.N*(S.N-1)/2), int(S.N*(S.N-1)*(S.N-2)/6), int(S.N*(S.N-1)*(S.N-2)*(S.N-3)/24)]
    UNCHANGE_LIMIT = sum(NEIGHBOR_WEIGHTS)

    T_INIT            = 1e-17
    T_ALPHA           = 0.99
    T_TEMP_INTERVAL   = 2000   # temperature step independent of watchdog
    WATCHDOG_INTERVAL = 50000  # status report interval (steps)

    def bha2value(bha_info):
        v = [i[1] for i in bha_info]
        n = len(v) + P_IDX  # P_IDX is negative, so take the smallest (len+P_IDX) elements
        return sum(heapq.nsmallest(n, v))

    def getOneStepNeighbor(param):
        _param = [param[0], [i for i in param[1]]]
        # if random.random() < 0.5:
        #     # flip one bit in operation_array
        #     idx_x = random.randint(0, (S.level*S.N//2)-1)
        #     flip = '0' if param[0][idx_x] == '1' else '1'
        #     _param[0] = param[0][:idx_x]+flip+param[0][idx_x+1:]
        # else:
        # swap two elements in permutation
        idx_a, idx_b = random.sample(range(S.N), 2)
        _param[1][idx_a], _param[1][idx_b] = _param[1][idx_b], _param[1][idx_a]
        return _param

    start_time  = time.time()

    for restart in range(5):
        param_init = [f'{random.randint(0, 2**(S.level*S.N//2)):0{S.level*S.N//2}b}', list(range(S.N))]
        random.shuffle(param_init[1])
        permm = [
            0, 122, 121, 123, 50, 51, 84, 124, 25, 85, 102, 52, 41, 26, 106, 125, 61, 103, 81, 86, 73, 42, 76, 53, 94, 107, 37, 27, 111, 62, 89, 126, 13, 82, 119, 104, 100, 74, 79, 87, 35, 77, 117, 43, 115, 95, 45, 54, 19, 38, 47, 108, 70, 112, 97, 28, 67, 90, 21, 63, 31, 14, 56, 127, 7, 120, 49, 83, 24, 101, 40, 105, 60, 80, 72, 75, 93, 36, 110, 88, 12, 118, 99, 78, 34, 116, 114, 44, 18, 46, 69, 96, 66, 20, 30, 55, 6, 48, 23, 39, 59, 71, 92, 109, 11, 98, 33, 113, 17, 68, 65, 29, 5, 22, 58, 91, 10, 32, 16, 64, 4, 57, 9, 15, 3, 8, 2, 1]
        # 0, 127, 126, 124, 120, 112, 96, 64,
        # 125, 122, 116, 104, 80, 32,
        # 61, 123, 118, 108, 88, 48,
        # 93, 58, 113, 98, 68, 8,
        # 13, 27, 55, 111, 94, 60,
        # 117, 106, 84, 40,
        # 77, 26, 49, 99, 70, 12,
        # 21, 43, 87, 46, 89, 50,
        # 97, 66, 4,
        # 5, 11, 23, 47, 95, 62,
        # 121, 114, 100, 72, 16,
        # 29, 59, 119, 110, 92, 56,
        # 109, 90, 52, 101, 74, 20,
        # 37, 75, 22, 41, 83, 38,
        # 73, 18, 33, 67, 6,
        # 9, 19, 39, 79, 30, 57,
        # 115, 102, 76, 24,
        # 45, 91, 54, 105, 82, 36,
        # 69, 10, 17, 35, 71, 14,
        # 25, 51, 103, 78, 28,
        # 53, 107, 86, 44,
        # 85, 42, 81, 34, 65, 2, 1, 3, 7, 15, 31, 63
        # ]

        param_init = [
        '1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111',
        # PD.inv_eirin(permm)
        permm
        ]
        
        param_init[1] = param_init[1][::-1]
        global_best = (param_init, float('inf'))
        S.set_system_parameter_v1(param_init[1], param_init[0])

        curentPair   = (param_init, bha2value(S.bha_info))
        run_best     = curentPair
        unchange_cnt = 0
        total_steps  = 0
        t            = T_INIT
        t_cnt        = 0
        watchdog_cnt = 0

        while unchange_cnt < UNCHANGE_LIMIT:
            tmpParam = curentPair[0]
            step = random.choices([i+1 for i in range(len(NEIGHBOR_WEIGHTS))], weights=NEIGHBOR_WEIGHTS)[0]
            for s_idx in range(step):
                tmpParam = getOneStepNeighbor(tmpParam)
            S.set_system_parameter_v1(tmpParam[1], tmpParam[0])
            tmpParamV = bha2value(S.bha_info)

            delta = tmpParamV - curentPair[1]
            if delta > 0:
                if random.random() < math.exp(-delta/t):
                    nextPair = (tmpParam, tmpParamV)
                else:
                    nextPair = curentPair
            else:
                nextPair = (tmpParam, tmpParamV)

            total_steps += 1
            if nextPair is curentPair:
                unchange_cnt += 1
            else:
                curentPair = nextPair
                unchange_cnt = 0

            # track run best and global best
            if nextPair[1] <= run_best[1]:
                run_best = nextPair
                if run_best[1] < global_best[1]:
                    global_best = run_best

            # temperature schedule (independent counter)
            t_cnt += 1
            if t_cnt >= T_TEMP_INTERVAL:
                t_cnt = 0
                t = max(t * T_ALPHA, 1e-17)

            # watchdog: report elapsed time, global best, run best
            watchdog_cnt += 1
            if watchdog_cnt >= WATCHDOG_INTERVAL:
                watchdog_cnt = 0
                elapsed = time.time() - start_time
                summary = (
                    f"[restart {restart+1:d}/5]  "
                    f"elapsed={elapsed:8.1f}s  "
                    f"steps={total_steps:10d}  "
                    f"unchange={unchange_cnt:7d}/{UNCHANGE_LIMIT:d}  "
                    f"t={t:10.3e}  "
                    f"current={curentPair[1]:12.6f}  "
                    f"run_best={run_best[1]:12.6f}  "
                    f"global_best={global_best[1]:12.6f}\n"
                )
                print(summary, end='', flush=True)
                with open(watchdog_log, 'w') as f:
                    f.write(summary)
                    f.write(f"run_best_op={run_best[0][0]}\n")
                    f.write(f"run_best_perm={run_best[0][1]}\n")
                    f.write(f"global_best_op={global_best[0][0]}\n")
                    f.write(f"global_best_perm={global_best[0][1]}\n")
                    S.set_system_parameter_v1(global_best[0][1], global_best[0][0])
                    bha_list = sorted([i[1] for i in S.bha_info])
                    bha_sum = [sum(bha_list[:l+1]) for l in range(len(bha_list))]
                    f.write("%s\n"%(", ".join(["%f"%x for x in bha_sum[::-1]])))

        # save final result of this restart
        elapsed = time.time() - start_time
        with open(result_log, 'a') as f:
            f.write(f"=== restart {restart+1} | elapsed={elapsed:.1f}s ===\n")
            f.write(f"value={run_best[1]:.6f}\n")
            f.write(f"op={run_best[0][0]}\n")
            f.write(f"perm={run_best[0][1]}\n")
            f.write("\n")
        print(f"[restart {restart+1}/5] done  best={run_best[1]:.6f}  global_best={global_best[1]:.6f}", flush=True)

    # append global best summary
    with open(result_log, 'a') as f:
        f.write(f"=== GLOBAL BEST ===\n")
        f.write(f"value={global_best[1]:.6f}\n")
        f.write(f"op={global_best[0][0]}\n")
        f.write(f"perm={global_best[0][1]}\n")
        f.write("\n")


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--run-dir', default='.')
    parser.add_argument('--P_IDX', type=int, default=-2)
    parser.add_argument('--seed', type=int, default=None)
    args = parser.parse_args()
    target_block(run_dir=args.run_dir, seed=args.seed, P_IDX=args.P_IDX)
