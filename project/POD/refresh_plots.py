#!/usr/bin/env python3
"""Re-render every comparison figure from whatever has landed in the logs.

    python3 refresh_plots.py              # all figures
    python3 refresh_plots.py m8t24 m8t12  # only these projects

Each project keeps its own plot_logs.py (they differ in colour schemes), so
decoder families and effective list sizes are read through that project's
decode_info().  A curve is drawn only if its log has at least one SNR line;
HD / OSD1 / MLD are added as shared reference curves and are never paired.
The m6t11 set is rendered by its own wrapper, plot_eq_vs_shear.sh.
"""
import importlib.util
import pathlib
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

ROOT = pathlib.Path(__file__).resolve().parent

VARIANT_RE = re.compile(r"^(?P<dec>.+)_(?P<var>P|GL|GLlo|GLhi)$")


def load_module(proj):
    path = ROOT / proj / "plot_logs.py"
    spec = importlib.util.spec_from_file_location(f"plot_logs_{proj}", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def n_points(d):
    log = d / "log.txt"
    if not log.exists():
        return 0
    return sum(1 for line in log.open(errors="replace") if "SNR" in line)


def select(proj, mod, pred, refs):
    """Folder names of decoder curves whose (family, eff) satisfy pred."""
    d = ROOT / proj
    names = []
    for f in sorted(d.iterdir()):
        m = VARIANT_RE.match(f.name)
        if not (m and f.is_dir()) or n_points(f) == 0:
            continue
        info = mod.decode_info(f.name)
        if info["family"] in ("SC", "SCL", "POD_SC", "POD_SCL", "CPC_SC", "CPC_SCL") and pred(info):
            names.append(f.name)
    return names + [r for r in refs if n_points(d / r) > 0]


ALL = lambda i: True
GE = lambda k: (lambda i: (i["eff_size"] or 0) >= k)
IN = lambda *ks: (lambda i: i["eff_size"] in ks)

# (project, output name, predicate, references, title)
P_GL = r"$P$ (solid) vs $P_{GL}$ optimum (dashed)"
T12 = r"$P$ (solid), $P_{GL}$ 2-3dB (dashed), $P_{GL}$ 3.5-6dB (dotted)"
FIGS = [
    # m7t10 figures carry the _han suffix
    ("m7t10_itw2026", "eBCH_m7_t10_P_vs_GLopt_han.png",
     lambda i: (i["eff_size"] or 0) < 127, ["HD", "OSD1"],
     r"eBCH(128,64): " + P_GL + r", HD, OSD$_1$"),
    ("m7t10_itw2026", "eBCH_m7_t10_PED127_P_GL_HD_OSD1_han.png",
     lambda i: (i.get("M") == 127 or i["eff_size"] == 128) and i.get("M") != 32,
     ["HD", "OSD1"],
     r"eBCH(128,64) t=10: PED$_{127}$ / CPC$_{127}$ and eff-128 splits, $P$ (solid) vs $P_{GL}$ (dashed), HD, OSD$_1$"),
    ("m7t10_itw2026", "eBCH_m7_t10_eff128_P_GL_HD_OSD1_han.png",
     lambda i: i["eff_size"] in (127, 128) and not i["family"].startswith("CPC"), ["HD", "OSD1"],
     r"eBCH(128,64) t=10, eff 128 splits (PED$_{127}$-SC included): $P$ (solid) vs $P_{GL}$ (dashed), HD, OSD$_1$"),
    ("m8t24_itw2026", "eBCH_m8_t24_P_vs_GLopt.png", ALL, ["HD", "OSD1"],
     r"eBCH(256,91) t=24: " + P_GL),
    ("m8t24_itw2026", "eBCH_m8_t24_eff8plus.png", GE(8), ["HD", "OSD1"],
     r"eBCH(256,91) t=24, eff$\geq$8: " + P_GL),
    ("m8t24_itw2026", "eBCH_m8_t24_eff8_eff32_HD.png", IN(8, 32), ["HD"],
     r"eBCH(256,91) t=24, eff 8 & 32: $P$ (solid) vs $P_{GL}$ (dashed), HD"),
    ("m8t24_itw2026", "eBCH_m8_t24_eff32_eff128_OSD1.png", IN(32, 128), ["OSD1"],
     r"eBCH(256,91) t=24, eff 32 & 128: $P$ (solid) vs $P_{GL}$ (dashed), OSD$_1$"),
    ("m8t19_itw2026", "eBCH_m8_t19_P_vs_GLopt.png", ALL, ["HD", "OSD1"],
     r"eBCH(256,123) t=19: " + P_GL),
    ("m8t19_itw2026", "eBCH_m8_t19_eff16plus.png", GE(16), ["HD", "OSD1"],
     r"eBCH(256,123) t=19, eff$\geq$16: " + P_GL),
    ("m8t19_itw2026", "eBCH_m8_t19_eff8_eff32_HD.png", IN(8, 32), ["HD"],
     r"eBCH(256,123) t=19, eff 8 & 32: $P$ (solid) vs $P_{GL}$ (dashed), HD"),
    ("m8t19_itw2026", "eBCH_m8_t19_eff32_eff128_OSD1.png", IN(32, 128), ["OSD1"],
     r"eBCH(256,123) t=19, eff 32 & 128: $P$ (solid) vs $P_{GL}$ (dashed), OSD$_1$"),
    ("m8t12_itw2026", "eBCH_m8_t12_P_vs_GLopt.png", ALL, ["HD", "OSD1"],
     r"eBCH(256,163) t=12: " + T12),
    ("m8t12_itw2026", "eBCH_m8_t12_eff32plus.png", GE(32), ["HD", "OSD1"],
     r"eBCH(256,163) t=12, eff$\geq$32: " + T12),
    ("m8t12_itw2026", "eBCH_m8_t12_eff32_HD.png", IN(32), ["HD"],
     r"eBCH(256,163) t=12, eff 32: " + T12 + ", HD"),
    ("m8t12_itw2026", "eBCH_m8_t12_eff128_OSD1.png", IN(128), ["OSD1"],
     r"eBCH(256,163) t=12, eff 128: " + T12 + r", OSD$_1$"),
    ("m8t12_itw2026", "eBCH_m8_t12_eff512_OSD1.png", IN(512), ["OSD1"],
     r"eBCH(256,163) t=12, eff 512: " + T12 + r", OSD$_1$"),
]


def render(fig, mods):
    proj, out, pred, refs, title = fig
    mod = mods[proj]
    names = select(proj, mod, pred, refs)
    if not names:
        return f"[skip] {proj}/{out}: nothing to plot"
    r = subprocess.run(
        [sys.executable, "plot_logs.py", *names, "--title", title, "--savepath", f"./{out}"],
        cwd=ROOT / proj, capture_output=True, text=True)
    ok = any(l.startswith("[ok]") for l in r.stdout.splitlines())
    return (f"[ok]   {proj}/{out}  ({len(names)} curves)" if ok
            else f"[FAIL] {proj}/{out}: {r.stderr.strip()[-200:]}")


def main():
    want = [a for a in sys.argv[1:] if not a.startswith("-")]
    figs = [f for f in FIGS if not want or any(f[0].startswith(w) for w in want)]
    mods = {p: load_module(p) for p in {f[0] for f in figs}}
    if not want or any("m6t11".startswith(w) or w.startswith("m6") for w in want):
        r = subprocess.run(["./plot_eq_vs_shear.sh"], cwd=ROOT / "m6t11_isit2026",
                           capture_output=True, text=True)
        n = sum(1 for l in r.stdout.splitlines() if l.startswith("[ok]"))
        print(f"[ok]   m6t11_isit2026: {n} figures via plot_eq_vs_shear.sh")
    with ThreadPoolExecutor(max_workers=4) as ex:
        for line in ex.map(lambda f: render(f, mods), figs):
            print(line, flush=True)


if __name__ == "__main__":
    main()
