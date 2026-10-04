#!/usr/bin/env python3
"""
Plots of the finite-difference (FD) approximations of q and p on the boundary, compared with the
analytical values, from the files written by export_fd_q_p:

    qFD_k=<k>_a=<a>_N=<N>_h=<h>.txt      pFD_k=..._h=....txt
    header : # k=.. a=.. N=.. h=.. scheme=..
    rows   : theta  Re(fd) Im(fd)  Re(exact) Im(exact)  abs_error

For every (k, a, N) and for q and p, one figure  fd_compare_<q|p>_k=<k>_a=<a>_N=<N>.png  with
    top-left     : Re  -> analytical curve (black) + FD values for every h (markers)
    top-right    : Im  -> same
    bottom-left  : |Re(fd) - Re(exact)| versus theta (log scale), one curve per h
    bottom-right : |Im(fd) - Im(exact)| versus theta (log scale), one curve per h

Usage:  python plot_fd.py [folder] [--skip-existing | --no-skip-existing]
"""
import argparse
import glob
import os
import re
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np

# ---- True = skip figures whose PNG already exists (override with --no-skip-existing)
SKIP_EXISTING = True

parser = argparse.ArgumentParser(description="Plot FD approximations of q and p vs analytical values.")
parser.add_argument("folder", nargs="?", default=".", help="folder containing the qFD/pFD .txt files")
parser.add_argument("--skip-existing", action=argparse.BooleanOptionalAction, default=SKIP_EXISTING)
args = parser.parse_args()
folder, skip_existing = args.folder, args.skip_existing

FD_RE = re.compile(
    r"^(?P<name>q|p)FD_k[=_](?P<k>.+?)_a[=_](?P<a>.+?)_N[=_](?P<N>\d+)_h[=_](?P<h>.+?)\.txt$"
)


def num(s):
    return float(s.replace("_", "."))


def fmt(x):
    return f"{x:g}"


def read_fd(path):
    """-> h, scheme, array with columns theta, Re fd, Im fd, Re ex, Im ex (sorted by theta)"""
    h, scheme, rows = None, "centered", []
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                m = re.search(r"\bh=(\S+)", line)
                if m:
                    h = float(m.group(1))
                m = re.search(r"scheme=(\S+)", line)
                if m:
                    scheme = m.group(1)
                continue
            t = line.split()
            if len(t) == 6:
                try:
                    rows.append([float(x) for x in t[:5]])
                except ValueError:
                    pass
    a = np.array(rows)
    if a.size:
        a = a[np.argsort(a[:, 0])]
    return h, scheme, a


def plot_group(name, key, items):
    k, a, N = key
    tag = f"k={fmt(k)}_a={fmt(a)}_N={N}"
    out = os.path.join(folder, f"fd_compare_{name}_{tag}.png")
    if skip_existing and os.path.exists(out):
        print("skipped (already exists):", out)
        return False

    items = sorted(items, key=lambda it: -it[0])           # largest h first
    items = [it for it in items if it[2].size and np.isfinite(it[2][:, 1:]).any()]
    if not items:
        print(f"[{name}, {tag}] no finite data, skipped")
        return False
    cmap = plt.get_cmap("viridis")
    colors = [cmap(i / max(len(items) - 1, 1)) for i in range(len(items))]

    fig, axes = plt.subplots(2, 2, figsize=(14, 9), sharex=True)
    th, ex = items[-1][2][:, 0], items[-1][2][:, 3:5]      # analytical values (same in every file)
    for j, part in enumerate(("Real part", "Imaginary part")):
        axes[0, j].plot(th, ex[:, j], "k-", lw=2, label=f"analytical {name}")
    for (h, scheme, d), c in zip(items, colors):
        lab = f"FD, h={h:.0e}"
        for j in range(2):
            axes[0, j].plot(d[:, 0], d[:, 1 + j], "o", ms=5, color=c, mfc="none", label=lab if j == 0 else None)
            err = np.abs(d[:, 1 + j] - d[:, 3 + j])
            err = np.where(err > 0, err, np.nan)            # log scale: drop exact zeros
            axes[1, j].semilogy(d[:, 0], err, "o-", ms=4, lw=1, color=c, label=lab if j == 0 else None)
    schemes = ", ".join(sorted({it[1] for it in items}))
    info = f"k={fmt(k)}, a={fmt(a)}, N={N}, scheme: {schemes}"
    for j, part in enumerate(("Real part", "Imaginary part")):
        axes[0, j].set_title(f"{name}: {part} - FD vs analytical\n{info}")
        axes[0, j].set_ylabel(f"{name}")
        axes[1, j].set_title(f"|{part[:2]}(FD) - {part[:2]}(analytical)|")
        axes[1, j].set_ylabel("absolute error")
        axes[1, j].set_xlabel(r"$\theta$ (rad)")
        for i in (0, 1):
            axes[i, j].grid(alpha=0.3, which="both")
    axes[0, 0].legend(fontsize=8, ncol=2)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    print("saved", out)
    return True


if __name__ == "__main__":
    groups = defaultdict(list)   # (name, (k,a,N)) -> [(h, scheme, data)]
    for path in glob.glob(os.path.join(folder, "*FD_k*.txt")):
        m = FD_RE.match(os.path.basename(path))
        if not m:
            continue
        h, scheme, d = read_fd(path)
        if h is None:
            h = num(m["h"])
        groups[(m["name"], (num(m["k"]), num(m["a"]), int(m["N"])))].append((h, scheme, d))
    if not groups:
        print("No qFD/pFD files found in", os.path.abspath(folder))
    made = 0
    for (name, key), items in sorted(groups.items()):
        made += bool(plot_group(name, key, items))
    if made:
        plt.show()
    else:
        print("Nothing new to plot.")