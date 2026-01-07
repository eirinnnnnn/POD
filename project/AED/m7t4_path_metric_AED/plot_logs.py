#!/usr/bin/env python3
# plot_logs.py
import argparse
import re
from pathlib import Path
from typing import List, Tuple, Optional

import matplotlib.pyplot as plt


LINE_RE = re.compile(
    r"SNR\s*=\s*([+-]?\d+(?:\.\d+)?)\s*,.*?BLER\s*=\s*\d+\s*/\s*\d+\s*=\s*([+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)"
)

def parse_log(log_path: Path) -> Tuple[List[float], List[float]]:
    """
    Parse SNR and BLER pairs from a log.txt.
    Returns (snr_list, bler_list).
    """
    snrs: List[float] = []
    blers: List[float] = []

    text = log_path.read_text(encoding="utf-8", errors="replace").splitlines()
    for line in text:
        m = LINE_RE.search(line)
        if not m:
            continue
        snr = float(m.group(1))
        bler = float(m.group(2))
        snrs.append(snr)
        blers.append(bler)

    # sort by SNR (in case the log isn't ordered)
    pairs = sorted(zip(snrs, blers), key=lambda x: x[0])
    if not pairs:
        return [], []
    snrs_sorted, blers_sorted = zip(*pairs)
    return list(snrs_sorted), list(blers_sorted)


def folder_label(folder: Path) -> str:
    # Use last path component as the curve label
    return folder.name if folder.name else str(folder)


def main():
    ap = argparse.ArgumentParser(
        description="Plot BLER vs SNR from multiple folders, each containing log.txt."
    )
    ap.add_argument(
        "folders",
        nargs="+",
        help="Folder paths that each contain a log.txt",
    )
    ap.add_argument(
        "--log-name",
        default="log.txt",
        help="Log filename inside each folder (default: log.txt)",
    )
    ap.add_argument(
        "--title",
        default="BLER vs SNR",
        help="Plot title",
    )
    ap.add_argument(
        "--out",
        default="bler_vs_snr.png",
        help="Output image path (e.g., .png or .pdf)",
    )
    ap.add_argument(
        "--show",
        action="store_true",
        help="Show the plot window (in addition to saving).",
    )
    ap.add_argument(
        "--no-grid",
        action="store_true",
        help="Disable grid lines.",
    )

    args = ap.parse_args()

    plt.figure()
    any_curve = False

    for folder_str in args.folders:
        folder = Path(folder_str).expanduser().resolve()
        log_path = folder / args.log_name

        if not log_path.exists():
            print(f"[warn] missing: {log_path}")
            continue

        snrs, blers = parse_log(log_path)
        if not snrs:
            print(f"[warn] no parsed points in: {log_path}")
            continue

        # Avoid zeros on log-scale; clip to a tiny positive value if needed
        blers_plot = [max(b, 1e-15) for b in blers]

        plt.semilogy(snrs, blers_plot, marker="o", linewidth=1.5, label=folder_label(folder))
        any_curve = True
        print(f"[info] {folder_label(folder)}: {len(snrs)} points from {log_path}")

    if not any_curve:
        raise SystemExit("[error] No curves were plotted. Check folder paths and log format.")

    plt.xlabel("SNR (dB)")
    plt.ylabel("BLER")
    plt.title(args.title)
    if not args.no_grid:
        plt.grid(True, which="both", linestyle="--", linewidth=0.5, alpha=0.6)
    plt.legend()
    plt.tight_layout()

    out_path = Path(args.out).expanduser().resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path, dpi=300)
    print(f"[ok] saved: {out_path}")

    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
