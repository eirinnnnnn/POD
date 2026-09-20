#!/usr/bin/env python3
"""
Plot BLER curves from PolarAED run logs.

Each run directory <name>/ holds a log.txt written by main.cpp, one line per
SNR point. The directory name is the legend entry.

Usage:
    python3 plot_bler.py --dir ../codes/n128_k100_highsnr \
        --runs SC SCL32 AE32SC_LTA AE32SC_UTL AE4SCL8_UTL \
        --title "(128,100) polar, high-SNR DE/GA design" --out fig3.png
"""

import argparse
import os
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

LINE_RE = re.compile(
    r"SNR\s*=\s*([+-]?\d+(?:\.\d+)?)\s*,.*?BLER\s*=\s*(\d+)\s*/\s*(\d+)\s*=\s*"
    r"([+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)")

# Style: colour = decoder family, dash = automorphism family.
STYLE = {
    "SC":          dict(color="0.35", ls="-",  marker="o", label="SC"),
    "SCL32":       dict(color="tab:orange", ls="-", marker="s", label="SCL 32"),
    "SCL8":        dict(color="tab:orange", ls="--", marker="s", label="SCL 8"),
    "OSD1":        dict(color="black", ls=(0, (1, 1)), marker="*",
                        label="OSD order 1"),
    "AE32SC_LTA":  dict(color="tab:red", ls=":", marker="v",
                        label="AE32-SC  LTA"),
    "AE32SC_UTL":  dict(color="tab:blue", ls="-", marker="^",
                        label="AE32-SC  UTL"),
    "AE32SC_PU":   dict(color="tab:green", ls="-", marker="D",
                        label="AE32-SC  BLTA (P$\\cdot$U)"),
    "AE32SC_RND":  dict(color="tab:purple", ls="--", marker="x",
                        label="AE32-SC  random A"),
    "AE4SCL8_UTL": dict(color="tab:cyan", ls="-.", marker="P",
                        label="AE4-SCL8  UTL"),
}

# The constant-effective-list-size family: M branches x list size L, M*L = 32.
# SCL32 is its M = 1 endpoint and AE32-SC its L = 1 endpoint, so one colour
# ramp over M tells the "trade list depth for branch diversity" story.
EFF32 = [("SCL32", 1, 32), ("AE2SCL16_UTL", 2, 16), ("AE4SCL8_UTL", 4, 8),
         ("AE8SCL4_UTL", 8, 4), ("AE16SCL2_UTL", 16, 2), ("AE32SC_UTL", 32, 1)]
EFF32_MARKERS = ["s", "o", "D", "^", "v", "P"]


def eff32_style():
    import matplotlib.cm as cm
    out = {}
    for k, (name, M, L) in enumerate(EFF32):
        out[name] = dict(color=cm.viridis(k / (len(EFF32) - 1.0)), ls="-",
                         marker=EFF32_MARKERS[k],
                         label="$M{=}%d$, SCL$_{%d}$" % (M, L)
                               if L > 1 else "$M{=}32$, SC")
    out["SCL32"]["label"] = "$M{=}1$, SCL$_{32}$  (plain SCL32)"
    return out


def read_log(path):
    pts = []
    with open(path) as f:
        for line in f:
            m = LINE_RE.search(line)
            if m:
                snr, err, it, bler = m.groups()
                pts.append((float(snr), int(err), int(it), float(bler)))
    # keep the last measurement per SNR (logs are appended to)
    by_snr = {}
    for snr, err, it, bler in pts:
        by_snr[snr] = (err, it, bler)
    return sorted((s, v[0], v[1], v[2]) for s, v in by_snr.items())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True)
    ap.add_argument("--runs", nargs="+", required=True)
    ap.add_argument("--title", default="")
    ap.add_argument("--out", default="bler.png")
    ap.add_argument("--ymin", type=float, default=1e-5)
    ap.add_argument("--eff32", action="store_true",
                    help="colour the constant-effective-list-size family by M")
    args = ap.parse_args()

    style = dict(STYLE)
    if args.eff32:
        style.update(eff32_style())

    fig, ax = plt.subplots(figsize=(6.8, 5.2))
    for name in args.runs:
        path = os.path.join(args.dir, name, "log.txt")
        if not os.path.exists(path):
            print("  [skip] %s (no log)" % name)
            continue
        pts = read_log(path)
        if not pts:
            print("  [skip] %s (empty log)" % name)
            continue
        st = dict(style.get(name, {}))
        st.setdefault("label", name)
        st.setdefault("marker", "o")
        x = [p[0] for p in pts]
        y = [max(p[3], 1e-12) for p in pts]
        ax.semilogy(x, y, ms=5, lw=1.6, **st)
        tail = ", ".join("%.2f:%.2e(%d/%d)" % (p[0], p[3], p[1], p[2])
                         for p in pts)
        print("  %-14s %s" % (name, tail))

    ax.set_xlabel("$E_b/N_0$ [dB]")
    ax.set_ylabel("BLER")
    ax.set_ylim(args.ymin, 1.0)
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=9)
    if args.title:
        ax.set_title(args.title, fontsize=11)
    fig.tight_layout()
    out = os.path.join(args.dir, args.out)
    fig.savefig(out, dpi=150)
    print("-> %s" % out)


if __name__ == "__main__":
    main()
