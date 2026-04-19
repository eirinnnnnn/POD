#!/usr/bin/env python3
# plot_logs.py
import argparse
import re
from pathlib import Path
from typing import List, Tuple, Dict

import matplotlib.pyplot as plt

LINE_RE = re.compile(
    r"SNR\s*=\s*([+-]?\d+(?:\.\d+)?)\s*,.*?BLER\s*=\s*\d+\s*/\s*\d+\s*=\s*([+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)"
)

# ===== POD legend mapping =====
AED_SC_RE  = re.compile(r"^AED(\d+)SC$", re.IGNORECASE)
AED_SCL_RE = re.compile(r"^AED(\d+)SCL(\d+)$", re.IGNORECASE)
SCL_RE     = re.compile(r"^SCL(\d+)$", re.IGNORECASE)
OSD_RE     = re.compile(r"^OSD(\d+)$", re.IGNORECASE)

# =========================
# Publication-grade palette
# =========================
# Deep, high-contrast, print-friendly colors (NO gray/cyan/purple/pink)
C_BLUE   = "#1f4ed8"   # royal blue
C_ORANGE = "#d97706"   # amber
C_GREEN  = "#15803d"   # forest green
C_RED    = "#b91c1c"   # deep red
C_BROWN  = "#7c2d12"   # chestnut
C_OLIVE  = "#4d7c0f"   # olive
C_TEAL   = "#0f766e"   # dark teal
C_GOLD   = "#ca8a04"   # gold
C_MAROON = "#7f1d1d"   # wine (still "red-family" but distinct)
C_NAVY   = "#1e3a8a"   # navy

COLOR_POOL = [
    C_BLUE, C_ORANGE, C_GREEN, C_RED,
    C_BROWN, C_OLIVE, C_TEAL, C_GOLD, C_MAROON, C_NAVY
]

# Enumerate BOTH (color, linestyle) within a family
POD_SC_LINESTYLES  = ["--", "-.", ":"]
POD_SCL_LINESTYLES = ["-.", ":", "--"]

def _combos(colors: List[str], linestyles: List[str]) -> List[Tuple[str, str]]:
    """Return list of (color, linestyle) pairs in stable order."""
    return [(c, ls) for ls in linestyles for c in colors]

# POD families: marker fixed per family, but enumerate (color, linestyle)
POD_SC_COMBOS = _combos(
    colors=[C_BROWN, C_BLUE, C_GREEN, C_RED, C_OLIVE, C_NAVY],
    linestyles=POD_SC_LINESTYLES,
)
POD_SCL_COMBOS = _combos(
    colors=[C_RED, C_GREEN, C_BROWN, C_BLUE, C_OLIVE, C_TEAL],
    linestyles=POD_SCL_LINESTYLES,
)

# (ii) Hard-fix these curves across ALL figures
# Keys are folder names (case-insensitive handling is done in get_decoder_style)
POD_FIXED: Dict[str, Tuple[str, str, str]] = {
    # folder      (linestyle, marker, color)
    "AED4SC":     ("--", "v", C_BROWN),
    "AED8SC":     ("-.", "v", C_BLUE),
    "AED8SCL2":   ("-",  "D", C_RED),
    "AED8SCL8":   ("-",  "D", "m"),
}

def get_decoder_style(folder_name: str) -> Tuple[str, str, str]:
    """
    Return (linestyle, marker, color) based on decoder type.

    Requirements implemented:
      (i) Within same POD family: enumerate both color and linestyle; marker fixed per family.
      (ii) POD_FIXED entries are identical across figures.
      (iii) Avoid light colors (cyan/purple/pink) and also avoid gray by construction.
    """
    name = folder_name.strip()
    up = name.upper()

    # Fixed POD styles first
    if up in POD_FIXED:
        return POD_FIXED[up]

    # Baselines
    if up == "HD":
        return "-", "o", C_NAVY
    if up == "SC":
        return "--", "s", C_ORANGE
    if up == "MLD":
        return "-.", None, "black"  # keep black for reference curve

    # SCL family: keep marker consistent, color varies with L
    # m = SCL_RE.match(name)
    # if m:
    #     L = int(m.group(1))
    #     scl_color_map = {
    #         2:  C_GREEN,
    #         4:  C_MAROON,
    #         8:  C_RED,
    #         16: C_BROWN,
    #         32: C_GREEN,
    #         64: C_NAVY,
    #     }
    #     color = scl_color_map.get(L, COLOR_POOL[L % len(COLOR_POOL)])
    #     return "-.", "P", color
    # Canonical SCL baseline ladder (explicit, cross-figure fixed)
    m = SCL_RE.match(name)
    if m:
        L = int(m.group(1))

        SCL_CANONICAL = {
            4:  ("--", "P", "tab:orange"),
            8:  ("-.", "P", "tab:pink"),
            32: (":",  "P", "tab:green"),
            64: ("-",  "P", "tab:gray"),
        }

        if L in SCL_CANONICAL:
            return SCL_CANONICAL[L]

        # fallback for any other SCL
        return "-.", "P", C_OLIVE


    # OSD family: keep marker consistent; avoid gray
    m = OSD_RE.match(name)
    if m:
        order = int(m.group(1))
        osd_colors = ["b", C_OLIVE, C_BROWN, C_BLUE, C_GREEN, C_RED, C_GOLD, C_NAVY]
        color = osd_colors[(order - 1) % len(osd_colors)]
        return "-.", "^", color

    # POD_M-SC (AED{M}SC): marker fixed = 'v', enumerate (color, linestyle)
    m = AED_SC_RE.match(name)
    if m:
        M = int(m.group(1))
        color, ls = POD_SC_COMBOS[M % len(POD_SC_COMBOS)]
        return ls, "v", color

    # POD_M-SCL_L (AED{M}SCL{L}): marker fixed = 'D', enumerate (color, linestyle)
    m = AED_SCL_RE.match(name)
    if m:
        M, L = int(m.group(1)), int(m.group(2))
        # deterministic mixing of (M, L) into an index
        idx = (1000 * M + L) % len(POD_SCL_COMBOS)
        color, ls = POD_SCL_COMBOS[idx]
        return ls, "D", color

    # Fallback: still avoid gray
    return "-", "X", C_GOLD


def parse_log(log_path: Path) -> Tuple[List[float], List[float]]:
    snrs, blers = [], []
    for line in log_path.read_text(encoding="utf-8", errors="replace").splitlines():
        m = LINE_RE.search(line)
        if m:
            snrs.append(float(m.group(1)))
            blers.append(float(m.group(2)))
    pairs = sorted(zip(snrs, blers), key=lambda x: x[0])
    return (list(zip(*pairs))[0], list(zip(*pairs))[1]) if pairs else ([], [])


def folder_label(folder: Path) -> str:
    name = folder.name

    m = AED_SC_RE.match(name)
    if m:
        M = m.group(1)
        return rf"$\mathrm{{POD}}_{{{M}}}\!-\!\mathrm{{SC}}$"

    m = AED_SCL_RE.match(name)
    if m:
        M, L = m.group(1), m.group(2)
        return rf"$\mathrm{{POD}}_{{{M}}}\!-\!\mathrm{{SCL}}_{{{L}}}$"

    m = SCL_RE.match(name)
    if m:
        L = m.group(1)
        return rf"$\mathrm{{SCL}}_{{{L}}}$"

    m = OSD_RE.match(name)
    if m:
        L = m.group(1)
        return rf"$\mathrm{{OSD}}_{{{L}}}$"

    return name


def main():
    ap = argparse.ArgumentParser(description="Plot BLER vs SNR from folders.")
    ap.add_argument("folders", nargs="+")
    ap.add_argument("--log-name", default="log.txt")
    ap.add_argument("--title", default="BLER vs SNR")
    ap.add_argument("--out", default="bler_vs_snr.png")
    ap.add_argument("--show", action="store_true")
    ap.add_argument("--no-grid", action="store_true")
    args = ap.parse_args()

    plt.figure()
    any_curve = False

    for f in args.folders:
        folder = Path(f).expanduser().resolve()
        log_path = folder / args.log_name
        if not log_path.exists():
            print(f"[warn] missing {log_path}")
            continue

        snrs, blers = parse_log(log_path)
        if not snrs:
            print(f"[warn] no data in {log_path}")
            continue

        blers = [max(b, 1e-15) for b in blers]
        label = folder_label(folder)

        if folder.name.upper() == "MLD":
            plt.semilogy(
                snrs, blers,
                linestyle="--", linewidth=1.2,
                color="black", marker=None,
                label=label, zorder=10
            )
        else:
            # linestyle, marker, color = get_decoder_style(folder.name)
            # plt.semilogy(
            #     snrs, blers,
            #     linestyle=linestyle,
            #     marker=marker,
            #     color=color,
            #     linewidth=1.0,
            #     markersize=5,
            #     markerfacecolor="none",
            #     label=label
            # )
            linestyle, marker, color = get_decoder_style(folder.name)
            is_landmark = folder.name.upper() in ("AED8SCL2", "AED8SCL8", "AED16SC")

            plt.semilogy(
                snrs, blers,
                linestyle=linestyle,
                marker=marker,
                color=color,
                linewidth=1.6 if is_landmark else 1.0,
                markersize=10 if is_landmark else 5,
                markerfacecolor="none",
                markeredgecolor= color if is_landmark else None,
                markeredgewidth=1.2 if is_landmark else None,
                zorder=5 if is_landmark else 3,
                label=label
            )


        print(f"[info] {label}: {len(snrs)} points")
        any_curve = True

    if not any_curve:
        raise SystemExit("No curves plotted.")

    plt.xlabel("$\mathrm{E_b/N_0}$ (dB)")
    plt.ylabel("BLER")
    plt.title(args.title)
    if not args.no_grid:
        plt.grid(True, which="both", linestyle="--", linewidth=0.5, alpha=0.6)

    plt.legend(
        loc="lower left",
        ncol=2,
        frameon=True,
        fontsize=10,
        columnspacing=1.2,
        handletextpad=0.6,
        borderaxespad=0.5
    )

    plt.tight_layout()
    out_path = Path(args.out).expanduser().resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path, dpi=300)
    print(f"[ok] saved to {out_path}")
    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
