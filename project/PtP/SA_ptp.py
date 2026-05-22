import os
import PolarDecoding_numpy as PD
import random, math, time, heapq


def target_block(run_dir: str = ".", seed: int = None, P_IDX: int = -2):

    os.makedirs(run_dir, exist_ok=True)
    watchdog_log = os.path.join(run_dir, "watchdog.log")
    result_log   = os.path.join(run_dir, "result.log")

    if seed is not None:
        random.seed(seed)

    matrix_path = "static_polar_32_16_H.matrix"
    S = PD.PolarDecodeSystem.from_H_matrix_file(level=5, matrix_path=matrix_path)

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
        if random.random() < 0.5:
            # flip one bit in operation_array
            idx_x = random.randint(0, (S.level*S.N//2)-1)
            flip = '0' if param[0][idx_x] == '1' else '1'
            _param[0] = param[0][:idx_x]+flip+param[0][idx_x+1:]
        else:
            # swap two elements in permutation
            idx_a, idx_b = random.sample(range(S.N), 2)
            _param[1][idx_a], _param[1][idx_b] = _param[1][idx_b], _param[1][idx_a]
        return _param

    start_time  = time.time()

    for restart in range(5):
        param_init = [f'{random.randint(0, 2**(S.level*S.N//2)):0{S.level*S.N//2}b}', list(range(S.N))]
        random.shuffle(param_init[1])

        param_init = [
        '11111111111111111111111111111111111111111111111111111111111111111111111111111111',
        [23, 7, 18, 2, 30, 11, 5, 27, 14, 0, 21, 9, 31, 16, 4, 25,
 12, 28, 1, 19, 6, 24, 10, 3, 29, 15, 22, 8, 26, 13, 20, 17]
        ]
        
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
            if nextPair[1] < run_best[1]:
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
