#!/usr/bin/env python3
#!/usr/bin/env python3
"""
Plots for etudie_convergence_airy (airy_check.cpp): checks the Airy-type truncation
estimate  N_series ~= ka + C*(ka)^(1/3)  for the series defining u+ at the illuminated
point theta = 0, across several target tolerances.

Reads:
  outputs/airyResu_<prefix>_k=<k>_a=<a>.txt
      header: # k=.. a=.. ka=.. Neff=..
      rows  : n  diff_n           (diff_n = |S_n - S_{n-1}|)
  outputs/airySummary_<prefix>.txt
      header: # k a ka tol n_star
      rows  : k a ka tol n_star    (long format: one row per (k, tol) pair; n_star = -1 if
                                     the series broke down before reaching that tolerance)

Produces, for each prefix found:
  airy_diff_<prefix>.png    : |S_n - S_{n-1}| vs n (log y), one curve per k, with a dotted
                              horizontal line at each tested tolerance
  airy_nstar_<prefix>.png   : n_star vs ka (log-log), one color per tolerance: measured points
                              (markers) plus the best-fit Airy curve ka + C*(ka)^(1/3) for that
                              tolerance (line), with the fitted C value given in the legend

Usage:  python plot_airy.py [folder] [--skip-existing | --no-skip-existing]
"""
import argparse
import glob
import os
import re
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np

SKIP_EXISTING = True

parser = argparse.ArgumentParser(description="Plot the Airy convergence test results.")
parser.add_argument("folder", nargs="?", default=".", help="folder containing the airy*.txt files")
parser.add_argument("--skip-existing", action=argparse.BooleanOptionalAction, default=SKIP_EXISTING)
args = parser.parse_args()
folder, skip_existing = args.folder, args.skip_existing

RESU_RE = re.compile(r"^airyResu_(?P<prefix>.+)_k[=_](?P<k>.+?)_a[=_](?P<a>.+?)\.txt$")


def num(s):
    return float(s.replace("_", "."))


def fmt(x):
    return f"{x:g}"


def read_diff(path):
    n, diff, ka = [], [], None
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                m = re.search(r"\bka=(\S+)", line)
                if m:
                    ka = float(m.group(1))
                continue
            t = line.split()
            if len(t) == 2:
                try:
                    n.append(float(t[0])); diff.append(float(t[1]))
                except ValueError:
                    pass
    return np.array(n), np.array(diff), ka


def read_summary(path):
    """-> structured array with columns k, a, ka, tol, n_star"""
    rows = []
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            t = line.split()
            if len(t) == 5:
                try:
                    rows.append([float(x) for x in t])
                except ValueError:
                    pass
    return np.array(rows)


def plot_diff(prefix, items, tol_values):
    out = os.path.join(folder, f"airy_diff_{prefix}.png")
    if skip_existing and os.path.exists(out):
        print("skipped (already exists):", out)
        return False
    items = sorted(items, key=lambda it: it[0])  # by k
    fig, ax = plt.subplots(figsize=(9, 6))
    cmap = plt.get_cmap("viridis")
    colors = [cmap(i / max(len(items) - 1, 1)) for i in range(len(items))]
    for (k, path), c in zip(items, colors):
        n, diff, ka = read_diff(path)
        diff = np.where(diff > 0, diff, np.nan)
        ax.semilogy(n, diff, "-", color=c, lw=1.3, label=f"k={fmt(k)} (ka={fmt(ka)})")
    for tol in sorted(tol_values):
        ax.axhline(tol, color="red", ls=":", lw=1)
        ax.text(ax.get_xlim()[1], tol, f" tol={tol:g}", color="red", fontsize=8,
                 va="center", ha="left", clip_on=False)
    ax.set_xlabel("n (number of terms)")
    ax.set_ylabel(r"$|S_n - S_{n-1}|$")
    ax.set_title(f"Convergence of u+ at the illuminated point ($\\theta=0$)\nprefix: {prefix}")
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=8, ncol=2)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    print("saved", out)
    return True


def plot_nstar(prefix, rows):
    out = os.path.join(folder, f"airy_nstar_{prefix}.png")
    if skip_existing and os.path.exists(out):
        print("skipped (already exists):", out)
        return False

    tol_values = sorted(set(rows[:, 3]), reverse=True)  # loosest tol first
    fig, ax = plt.subplots(figsize=(9, 6.5))
    cmap = plt.get_cmap("plasma")
    colors = [cmap(i / max(len(tol_values) - 1, 1)) for i in range(len(tol_values))]

    for tol, color in zip(tol_values, colors):
        sel = rows[rows[:, 3] == tol]
        sel = sel[np.argsort(sel[:, 2])]
        ka, n_star = sel[:, 2], sel[:, 4]
        valid = n_star > 0
        if (~valid).any():
            print(f"[{prefix}, tol={tol:g}] n_star not reached for ka = {ka[~valid].tolist()}")
        ax.loglog(ka[valid], n_star[valid], "o", color=color, ms=6,
                   label=f"measured, tol={tol:g}")
        if valid.sum() >= 2:
            C_fit = np.median((n_star[valid] - ka[valid]) / ka[valid] ** (1. / 3.))
            ka_fine = np.linspace(ka[valid].min(), ka[valid].max(), 200)
            ax.loglog(ka_fine, ka_fine + C_fit * ka_fine ** (1. / 3.), "-", color=color, lw=1.8,
                       label=rf"fit tol={tol:g}: $N=ka+{C_fit:.2f}\,(ka)^{{1/3}}$")
            print(f"[{prefix}, tol={tol:g}] best-fit C (median over k) = {C_fit:.3f}")

    ax.set_xlabel("ka")
    ax.set_ylabel(r"$N$ needed for convergence")
    ax.set_title(f"Truncation order vs ka: measured vs best-fit Airy estimate\nprefix: {prefix}")
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    print("saved", out)
    return True


if __name__ == "__main__":
    groups = defaultdict(list)  # prefix -> [(k, path), ...]
    for path in glob.glob(os.path.join(folder, "airyResu_*.txt")):
        m = RESU_RE.match(os.path.basename(path))
        if m:
            groups[m["prefix"]].append((num(m["k"]), path))

    made = 0
    summaries = {}
    for path in glob.glob(os.path.join(folder, "airySummary_*.txt")):
        prefix = os.path.basename(path)[len("airySummary_"):-len(".txt")]
        rows = read_summary(path)
        if rows.size:
            summaries[prefix] = rows

    for prefix, items in sorted(groups.items()):
        tol_values = set(summaries[prefix][:, 3]) if prefix in summaries else set()
        made += bool(plot_diff(prefix, items, tol_values))

    for prefix, rows in sorted(summaries.items()):
        made += bool(plot_nstar(prefix, rows))

    if not groups and not summaries:
        print("No airyResu_*/airySummary_*.txt files found in", os.path.abspath(folder))
    if made:
        plt.show()
    else:
        print("Nothing new to plot.")