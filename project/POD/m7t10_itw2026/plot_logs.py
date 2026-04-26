#!/usr/bin/env python3
import argparse
import re
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
import numpy as np

LINE_RE = re.compile(
    r"SNR\s*=\s*([+-]?\d+(?:\.\d+)?)\s*,.*?BLER\s*=\s*\d+\s*/\s*\d+\s*=\s*([+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)"
)

# Supported names:
#   HD
#   SC
#   SCL4, SCL8, ...
#   POD4SC, POD8SCL4, ...
#   AED4SC, AED8SCL4, ...   (legacy alias -> treated as POD)
#   MLD
#   OSD1, OSD2, ...
HD_RE = re.compile(r"^HD$", re.IGNORECASE)
SC_RE = re.compile(r"^SC$", re.IGNORECASE)
SCL_RE = re.compile(r"^SCL(\d+)$", re.IGNORECASE)
POD_SC_RE = re.compile(r"^(?:POD|AED)(\d+)SC$", re.IGNORECASE)
POD_SCL_RE = re.compile(r"^(?:POD|AED)(\d+)SCL(\d+)$", re.IGNORECASE)
MLD_RE = re.compile(r"^MLD$", re.IGNORECASE)
OSD_RE = re.compile(r"^OSD(\d+)$", re.IGNORECASE)

# Effective list sizes to be represented consistently across all plots.
EFF_SIZES = [1, 4, 8, 16, 32, 64, 128, 256]

# Color = decoder family / type.
# For POD/PED, this is only the fallback base color.
FAMILY_COLOR: Dict[str, str] = {
    "HD": "tab:blue",
    "SC": "tab:orange",
    "SCL": "tab:orange",
    "POD_SC": "tab:green",
    "POD_SCL": "tab:green",
    "MLD": "black",
    "OSD": "tab:blue",
    "UNKNOWN": "tab:blue",
}

# Dark green palette for PED/POD ensemble size M.
# The first one starts from the normal matplotlib tab:green color.
PED_M_COLOR: Dict[int, str] = {
    4: "#2ca02c",   # normal tab:green
    # 32: "#238823",
    # 16: "#1f6f1f",
    32: "#104510",
    64: "#104510",
    # 64: "#185c18",
    # 64: "#104510",
}

# Marker = effective list size.
EFF_MARKER: Dict[int, str] = {
    1: "o",
    4: "s",
    8: "^",
    16: "D",
    32: "v",
    64: "P",
    128: "X",
    256: "*",
}

# Fallback marker when a decoder does not have one of the canonical effective sizes.
FAMILY_FALLBACK_MARKER: Dict[str, str] = {
    "HD": "o",
    "SC": "o",
    "SCL": "o",
    "POD_SC": "o",
    "POD_SCL": "o",
    "MLD": "X",
    "OSD": "^",
    "UNKNOWN": "o",
}

# Linestyle = broad family.
LINESTYLE_MAP: Dict[str, object] = {
    "HD": "-",
    "SC": "--",
    "SCL": "--",
    "POD_SC": "-",
    "POD_SCL": "-",
    "MLD": "-",
    "OSD": (0, (3, 1, 1, 1)),
    "UNKNOWN": "-",
}

EMPH_COLOR_BY_FAMILY: Dict[str, str] = {
    "POD_SCL": "tab:red",
    "UNKNOWN": "tab:blue"
}

# Emphasis should not destroy the semantic color encoding.
# Therefore emphasized curves keep their original color, but become visually stronger.
def apply_emphasis(style: Dict[str, object]) -> Dict[str, object]:
    emph = dict(style)
    emph["linewidth"] = max(float(style["linewidth"]), 2.0)
    emph["markersize"] = max(float(style["markersize"]), 10.0)
    emph["markeredgecolor"] = "tab:red"
    emph["color"] = "tab:red"
    emph["markeredgewidth"] = 1.8
    emph["zorder"] = 6
    return emph


def parse_log(log_path: Path) -> Tuple[List[float], List[float]]:
    snrs, blers = [], []
    text = log_path.read_text(encoding="utf-8", errors="replace")
    for line in text.splitlines():
        m = LINE_RE.search(line)
        if m:
            snrs.append(float(m.group(1)))
            blers.append(float(m.group(2)))

    pairs = sorted(zip(snrs, blers), key=lambda x: x[0])
    if not pairs:
        return [], []
    xs, ys = zip(*pairs)
    return list(xs), list(ys)


def decode_info(name: str) -> Dict[str, object]:
    s = name.strip()

    if HD_RE.match(s):
        return {
            "family": "HD",
            "eff_size": 1,
            "M": None,
            "L": None,
            "label": r"$\mathrm{HD}$",
        }

    if SC_RE.match(s):
        return {
            "family": "SC",
            "eff_size": 1,
            "M": None,
            "L": 1,
            "label": r"$\mathrm{SC}$",
        }

    m = SCL_RE.match(s)
    if m:
        L = int(m.group(1))
        return {
            "family": "SCL",
            "eff_size": L,
            "M": None,
            "L": L,
            "label": rf"$\mathrm{{SCL}}_{{{L}}}$",
        }

    m = POD_SC_RE.match(s)
    if m:
        M = int(m.group(1))
        return {
            "family": "POD_SC",
            "eff_size": M,
            "M": M,
            "L": 1,
            "label": rf"$\mathrm{{PED}}_{{{M}}}\!-\!\mathrm{{SC}}$",
        }

    m = POD_SCL_RE.match(s)
    if m:
        M = int(m.group(1))
        L = int(m.group(2))
        return {
            "family": "POD_SCL",
            "eff_size": M * L,
            "M": M,
            "L": L,
            "label": rf"$\mathrm{{PED}}_{{{M}}}\!-\!\mathrm{{SCL}}_{{{L}}}$",
        }

    if MLD_RE.match(s):
        return {
            "family": "MLD",
            "eff_size": None,
            "M": None,
            "L": None,
            "label": r"$\mathrm{MLD}$",
        }

    m = OSD_RE.match(s)
    if m:
        order = int(m.group(1))
        return {
            "family": "OSD",
            "eff_size": order,
            "M": None,
            "L": None,
            "label": rf"$\mathrm{{OSD}}_{{{order}}}$",
        }

    return {
        "family": "UNKNOWN",
        "eff_size": None,
        "M": None,
        "L": None,
        "label": s,
    }


def color_of(info: Dict[str, object]) -> str:
    family = str(info["family"])

    # For PED/POD curves:
    #   hue = PED family, fixed as green
    #   darkness = ensemble size M
    if family in {"POD_SC", "POD_SCL"}:
        M = info.get("M", None)
        if isinstance(M, int):
            return PED_M_COLOR.get(M, FAMILY_COLOR.get(family, "tab:green"))

    return FAMILY_COLOR.get(family, "tab:gray")


def style_of(name: str) -> Dict[str, object]:
    info = decode_info(name)
    family = str(info["family"])
    eff_size = info["eff_size"]

    color = color_of(info)
    marker = EFF_MARKER.get(eff_size, FAMILY_FALLBACK_MARKER.get(family, "o"))
    linestyle = LINESTYLE_MAP.get(family, "-")

    linewidth = 1.8 if family in {"HD", "MLD"} else 1.3
    markersize = 6 if family not in {"MLD"} else 7

    return {
        "color": color,
        "marker": marker,
        "linestyle": linestyle,
        "linewidth": linewidth,
        "markersize": markersize,
        "label": info["label"],
        "family": family,
        "eff_size": eff_size,
        "M": info.get("M", None),
        "L": info.get("L", None),
    }


def sort_key(folder_name: str) -> Tuple[int, int, int, int, str]:
    info = decode_info(folder_name)
    family_order = {
        "HD": 0,
        "SC": 1,
        "SCL": 2,
        "POD_SC": 3,
        "POD_SCL": 4,
        "MLD": 5,
        "OSD": 6,
        "UNKNOWN": 7,
    }

    eff = info["eff_size"]
    M = info.get("M", None)
    L = info.get("L", None)

    eff_rank = EFF_SIZES.index(eff) if eff in EFF_SIZES else 999
    M_rank = int(M) if isinstance(M, int) else 999
    L_rank = int(L) if isinstance(L, int) else 999

    return (
        family_order.get(str(info["family"]), 999),
        eff_rank,
        M_rank,
        L_rank,
        folder_name.upper(),
    )


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Plot BLER vs SNR from decoder log folders.")
    p.add_argument("folders", nargs="+", help="Decoder result folders.")
    p.add_argument("--title", default="BLER vs SNR", help="Plot title.")
    p.add_argument("--savepath", default="./bler_vs_snr.png", help="Output figure path.")
    p.add_argument(
        "--snr",
        nargs=3,
        type=float,
        metavar=("START", "STOP", "STEP"),
        default=(2.0, 4.6, 0.25),
        help="SNR axis setup: START STOP STEP",
    )
    p.add_argument(
        "--emph",
        nargs="*",
        default=[],
        metavar="FOLDER_NAME",
        help="Emphasize the specified folder names, matched against folder basename.",
    )
    return p


def main() -> None:
    args = build_parser().parse_args()

    snr_start, snr_stop, snr_step = args.snr
    xticks = np.arange(snr_start, snr_stop + 0.5 * snr_step, snr_step)

    plt.figure(figsize=(8.0, 5.6))
    plotted = False
    emph_names = {name.strip().upper() for name in args.emph}

    for folder_str in sorted(args.folders, key=lambda s: sort_key(Path(s).name)):
        folder = Path(folder_str).expanduser().resolve()
        log_path = folder / "log.txt"
        if not log_path.exists():
            print(f"[warn] missing {log_path}")
            continue

        snrs, blers = parse_log(log_path)
        if not snrs:
            print(f"[warn] no valid BLER lines in {log_path}")
            continue

        st = style_of(folder.name)
        is_emph = folder.name.upper() in emph_names
        if is_emph:
            st = apply_emphasis(st)

        blers = [max(b, 1e-15) for b in blers]

        plt.semilogy(
            snrs,
            blers,
            label=st["label"],
            color=st["color"],
            linestyle=st["linestyle"],
            marker=st["marker"],
            linewidth=st["linewidth"],
            markersize=st["markersize"],
            markerfacecolor="none",
            markeredgecolor=st.get("markeredgecolor", st["color"]),
            markeredgewidth=st.get("markeredgewidth", 1.0),
            zorder=st.get("zorder", 3),
        )

        eff = st["eff_size"]
        M = st["M"]
        L = st["L"]

        eff_text = f", eff={eff}" if eff is not None else ""
        ped_text = f", M={M}, L={L}" if M is not None else ""
        emph_text = ", emphasized" if is_emph else ""

        print(f"[info] {folder.name}: {len(snrs)} points{eff_text}{ped_text}{emph_text}")
        plotted = True

    if not plotted:
        raise SystemExit("No curves were plotted.")

    plt.xlabel(r"$\mathrm{E_b/N_0}$ (dB)")
    plt.ylabel("BLER")
    plt.title(args.title)
    plt.xlim(snr_start, snr_stop)
    plt.xticks(xticks)
    plt.grid(True, which="both", linestyle="--", linewidth=0.5, alpha=0.6)

    plt.legend(
        loc="lower left",
        ncol=2,
        frameon=True,
        fontsize=10,
        columnspacing=1.0,
        handletextpad=0.5,
        borderaxespad=0.4,
    )

    out_path = Path(args.savepath).expanduser().resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    print(f"[ok] saved to {out_path}")


if __name__ == "__main__":
    main()