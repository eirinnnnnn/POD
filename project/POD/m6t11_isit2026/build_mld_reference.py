#!/usr/bin/env python3
"""Pool every MLD measurement for eBCH(64,16) into one reference curve.

Sources:
  * project/POD/m6t11_isit2026/MLD              -- this project's original run
  * project/AED/m6t11_isit2026/MLD_16_64        -- earlier run, different
                                                   generator matrix
  * project/MLD/m6t11_isit2026/topup_*_s*       -- 8-seed high-SNR top-up

The AED run uses data/BCH_16_64_g.matrix rather than eBCH_m6_t11.matrix.
The two generator matrices do not share a row space, but they span
permutation-equivalent codes (identical weight enumerator: d=24, A24=5040,
A28=12544, A32=30366), and ML decoding is invariant under a coordinate
permutation -- so the runs measure the same quantity and their error and
iteration counts can be summed.

Every run stops on a fixed error count, so pooling is
sum(errors) / sum(iterations) per SNR.

Writes mld_merged/MLD/log.txt in the POD experiment directory, which is what
plot_eq_vs_shear.sh plots as the MLD reference.
"""
import re
import pathlib
from collections import defaultdict

ROOT = pathlib.Path("/home/eirin/Polar_eirin_20260107")
HERE = ROOT / "project/POD/m6t11_isit2026"
OUT = HERE / "mld_merged/MLD/log.txt"

LINE = re.compile(
    r"SNR\s*=\s*([\d.]+),\s*rawBER\s*=\s*([\d.]+).*?BLER\s*=\s*(\d+)\s*/\s*(\d+)"
)


def sources():
    yield HERE / "MLD/log.txt"
    yield ROOT / "project/AED/m6t11_isit2026/MLD_16_64/log.txt"
    yield from sorted((ROOT / "project/MLD/m6t11_isit2026").glob("topup_*_s*/log.txt"))


def main():
    err = defaultdict(int)
    itr = defaultdict(int)
    raw = {}
    used = 0
    for path in sources():
        if not path.exists():
            continue
        used += 1
        for line in path.read_text().splitlines():
            m = LINE.search(line)
            if not m:
                continue
            snr = float(m.group(1))
            raw.setdefault(snr, float(m.group(2)))
            err[snr] += int(m.group(3))
            itr[snr] += int(m.group(4))

    OUT.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    for snr in sorted(err):
        rows.append(
            "SNR = {:5.2f}, rawBER = {:5.2f}, BLER = {:8d}/{:8d} = {:8.6f}\n".format(
                snr, raw[snr], err[snr], itr[snr], err[snr] / itr[snr]
            )
        )
    OUT.write_text("".join(rows))
    print(f"pooled {used} source logs -> {OUT}")
    print("".join(rows), end="")


if __name__ == "__main__":
    main()
