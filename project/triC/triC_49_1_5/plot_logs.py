#!/usr/bin/env python3
import argparse
import re
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
import numpy as np

LINE_RE = re.compile(
    r"\b(SNR|p)\b\s*=\s*([+-]?\d+(?:\.\d+)?)\s*,.*?BLER\s*=\s*\d+\s*/\s*\d+\s*=\s*([+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)"
)

XLABEL_BY_CHANNEL: Dict[str, str] = {
    "SNR": r"$\mathrm{E_b/N_0}$ (dB)",
    "p": r"crossover probability $p$",
}

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


def parse_log(log_path: Path) -> Tuple[List[float], List[float], str]:
    xs, ys, channel = [], [], ""
    text = log_path.read_text(encoding="utf-8", errors="replace")
    for line in text.splitlines():
        m = LINE_RE.search(line)
        if m:
            channel = channel or m.group(1)
            xs.append(float(m.group(2)))
            ys.append(float(m.group(3)))

    pairs = sorted(zip(xs, ys), key=lambda x: x[0])
    if not pairs:
        return [], [], channel
    xs, ys = zip(*pairs)
    return list(xs), list(ys), channel


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
    p = argparse.ArgumentParser(description="Plot BLER vs SNR (AWGN) or vs p (BSC) from decoder log folders.")
    p.add_argument("folders", nargs="+", help="Decoder result folders.")
    p.add_argument("--title", default=None, help="Plot title.")
    p.add_argument("--savepath", default=None, help="Output figure path.")
    p.add_argument(
        "--snr",
        nargs=3,
        type=float,
        metavar=("START", "STOP", "STEP"),
        default=(2.0, 4.6, 0.25),
        help="SNR axis setup: START STOP STEP (used when logs report SNR).",
    )
    p.add_argument(
        "--p",
        nargs=3,
        type=float,
        metavar=("START", "STOP", "STEP"),
        default=None,
        help="Crossover-probability axis setup: START STOP STEP "
             "(used when logs report p; defaults to the data range if omitted).",
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

    plt.figure(figsize=(8.0, 5.6))
    plotted = False
    emph_names = {name.strip().upper() for name in args.emph}
    channel_kind = ""
    all_xs: List[float] = []

    # for folder_str in sorted(args.folders, key=lambda s: sort_key(Path(s).name)):
    for folder_str in args.folders:
        folder = Path(folder_str).expanduser().resolve()
        log_path = folder / "log.txt"
        if not log_path.exists():
            print(f"[warn] missing {log_path}")
            continue

        xs, blers, folder_channel = parse_log(log_path)
        if not xs:
            print(f"[warn] no valid BLER lines in {log_path}")
            continue

        if folder_channel:
            if channel_kind and channel_kind != folder_channel:
                print(f"[warn] {folder.name}: channel '{folder_channel}' differs from "
                      f"previously seen '{channel_kind}'; using '{channel_kind}' for axis setup")
            else:
                channel_kind = folder_channel
        all_xs.extend(xs)

        st = style_of(folder.name)
        is_emph = folder.name.upper() in emph_names
        if is_emph:
            st = apply_emphasis(st)

        blers = [max(b, 1e-15) for b in blers]

        plt.semilogy(
            xs,
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

        print(f"[info] {folder.name}: {len(xs)} points{eff_text}{ped_text}{emph_text}")
        plotted = True

    if not plotted:
        raise SystemExit("No curves were plotted.")

    channel_kind = channel_kind or "SNR"
    if channel_kind == "p":
        if args.p is not None:
            x_start, x_stop, x_step = args.p
            xticks = np.arange(x_start, x_stop + 0.5 * x_step, x_step)
        else:
            x_start, x_stop = min(all_xs), max(all_xs)
            xticks = sorted(set(all_xs))
        default_title = "BLER vs p"
        default_savepath = "./bler_vs_p.png"
    else:
        x_start, x_stop, x_step = args.snr
        xticks = np.arange(x_start, x_stop + 0.5 * x_step, x_step)
        default_title = "BLER vs SNR"
        default_savepath = "./bler_vs_snr.png"

    plt.xlabel(XLABEL_BY_CHANNEL.get(channel_kind, channel_kind))
    plt.ylabel("BLER")
    plt.title(args.title or default_title)
    plt.xlim(x_start, x_stop)
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

    out_path = Path(args.savepath or default_savepath).expanduser().resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    print(f"[ok] saved to {out_path}")


if __name__ == "__main__":
    main()