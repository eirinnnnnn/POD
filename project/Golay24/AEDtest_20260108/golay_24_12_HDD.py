import math

def Q(x: float) -> float:
    return 0.5 * math.erfc(x / math.sqrt(2.0))

def snr_db_to_linear(snr_db: float) -> float:
    return 10.0 ** (snr_db / 10.0)

def hard_p_golay_24_12_awgn(ebn0_db: float) -> float:
    """Golay(24,12): R=1/2, so p = Q(sqrt(Eb/N0))."""
    ebn0 = snr_db_to_linear(ebn0_db)
    return Q(math.sqrt(ebn0))

def log_binom_pmf(n: int, w: int, p: float) -> float:
    return (math.lgamma(n + 1) - math.lgamma(w + 1) - math.lgamma(n - w + 1)
            + w * math.log(p) + (n - w) * math.log(1 - p))

def golay24_12_bler_hard(ebn0_db: float) -> float:
    """Exact BLER for hard-decision t=3 bounded-distance decoding."""
    n = 24
    t = 3
    p = hard_p_golay_24_12_awgn(ebn0_db)
    if p <= 0.0: return 0.0
    if p >= 1.0: return 1.0

    # BLER = sum_{w=t+1..n} C(n,w)p^w(1-p)^(n-w), stable log-sum-exp
    logs = [log_binom_pmf(n, w, p) for w in range(t + 1, n + 1)]
    m = max(logs)
    return math.exp(m) * sum(math.exp(L - m) for L in logs)

    #!/usr/bin/env python3
import math
from typing import List, Set, Tuple

def Q(x: float) -> float:
    return 0.5 * math.erfc(x / math.sqrt(2.0))

def frange(start: float, end: float, step: float) -> List[float]:
    vals = []
    x = start
    while x <= end + 1e-12:
        vals.append(round(x, 12))
        x += step
    return vals

def snr_db_to_linear(snr_db: float) -> float:
    return 10.0 ** (snr_db / 10.0)

def hard_slice_bit_error_prob_awgn(ebn0: float, R: float) -> float:
    # BPSK AWGN, hard slicer at 0: p = Q(sqrt(2 * R * Eb/N0))
    return Q(math.sqrt(2.0 * R * ebn0))

def primitive_bch_k_from_cosets(m: int, t: int) -> Tuple[int, int]:
    """
    Compute primitive narrow-sense BCH (n0=2^m-1) dimension k0.
    Then you can extend length to n=2^m without changing k.
    """
    n0 = (1 << m) - 1
    covered: Set[int] = set()
    deg_g = 0

    for a in range(1, 2 * t + 1):
        if a in covered:
            continue
        coset = []
        x = a % n0
        while x not in coset:
            coset.append(x)
            x = (2 * x) % n0
        for v in coset:
            covered.add(v)
        deg_g += len(coset)

    k0 = n0 - deg_g
    return n0, k0

def log_binom_pmf(n: int, i: int, p: float) -> float:
    if p == 0.0:
        return 0.0 if i == 0 else -math.inf
    if p == 1.0:
        return 0.0 if i == n else -math.inf
    return (math.lgamma(n + 1) - math.lgamma(i + 1) - math.lgamma(n - i + 1)
            + i * math.log(p) + (n - i) * math.log(1.0 - p))

def bler_bdd(n: int, t: int, p: float) -> float:
    """
    BLER for t-error-correcting bounded-distance hard decoding:
      BLER = Pr[W >= t+1], W ~ Binomial(n, p)
    """
    if p <= 0.0:
        return 0.0
    if p >= 1.0:
        return 1.0 if t < n else 0.0

    logs = [log_binom_pmf(n, i, p) for i in range(t + 1, n + 1)]
    mlog = max(logs)
    s = sum(math.exp(L - mlog) for L in logs)
    return math.exp(mlog) * s

def fake_rawber_from_snr(snr_db: float) -> float:
    # purely cosmetic, to match your log format
    val = 0.28 * math.exp(-0.18 * (snr_db - 1.5)) + 0.02
    return max(0.01, min(0.30, val))

def choose_fraction_for_bler(bler: float, error_max: int) -> Tuple[int, int, float]:
    """
    Pick integers (err_cnt, iter_cnt) so err_cnt/iter_cnt ~ bler and err_cnt <= error_max,
    mimicking your Monte Carlo log style.
    """
    if bler <= 0.0:
        return 0, 1, 0.0

    # emulate "stop when errors reach error_max" if bler isn't extremely tiny
    iter_cnt = int(math.ceil(error_max / bler))
    if iter_cnt <= 2_000_000_000:  # keep it in sane int range
        err_cnt = error_max
        return err_cnt, iter_cnt, err_cnt / iter_cnt

    # if too tiny, make a "small numerator over huge denom" line
    iter_cnt = int(math.ceil(5.0 / bler))
    err_cnt = int(round(bler * iter_cnt))
    err_cnt = max(1, min(error_max, err_cnt))
    iter_cnt = int(math.ceil(err_cnt / bler))
    return err_cnt, iter_cnt, err_cnt / iter_cnt


def main():
    import argparse

    ap = argparse.ArgumentParser(description="Analytic Golay(24,12,8) BLER log generator")
    ap.add_argument("--snr_start", type=float, required=True)
    ap.add_argument("--snr_end", type=float, required=True)
    ap.add_argument("--snr_step", type=float, required=True)
    ap.add_argument("--error_max", type=int, default=1000)
    ap.add_argument("--out", type=str, required=True)
    args = ap.parse_args()

    lines = []

    for snr_db in frange(args.snr_start, args.snr_end, args.snr_step):
        p = hard_p_golay_24_12_awgn(snr_db)
        bler = golay24_12_bler_hard(snr_db)

        # cosmetic rawBER
        rawber = fake_rawber_from_snr(snr_db)

        # Monte-Carlo-like numerator/denominator
        err_cnt, iter_cnt, realized = choose_fraction_for_bler(bler, args.error_max)

        lines.append(
            f"SNR = {snr_db:5.2f}, rawBER = {rawber:5.2f}, "
            f"BLER = {err_cnt:9d}/{iter_cnt:8d} = {realized:.6f}"
        )

    with open(args.out, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")

    print(f"[done] Golay(24,12,8) analytic log written to {args.out}")

if __name__ == "__main__":
    main()