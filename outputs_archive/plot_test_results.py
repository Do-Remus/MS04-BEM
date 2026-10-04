#!/usr/bin/env python3
"""
Plots for the partial sums of q and p produced by main.cpp.

Expected file names (as written by main.cpp):
  final approximations : {q|p}Resu_k=<k>_a=<a>_N=<N>.txt
  convergence (1/theta): {q|p|solExt}Resu_k=<k>_a=<a>_theta=<theta>_N=<N>.txt
  (solExt = partial sums of u_N^+ at the exterior points; convergence plots only)
(the '=' may also be '_' and the decimal point '_', as in some upload systems)

For every (k, a, N) combination found, and for both q and p, it produces:
  convergence_<name>_k=<k>_a=<a>_N=<N>.png : x = n, y = partial sum (Re | Im), all thetas superposed
  error_<name>_k=<k>_a=<a>_N=<N>.png       : x = n, y = error (log scale), for the names in ERROR_NAMES
  final_<name>_k=<k>_a=<a>_N=<N>.png       : x = theta, y = final value (Re and Im)
  error_<name>_k=<k>_a=<a>_N=<N>.png       : (solExt) error vs n on a log scale
k, a and N also appear in the figure titles.

Usage:  python plot_results.py [folder] [--skip-existing | --no-skip-existing]

  --skip-existing      only plot figures whose PNG does not exist yet (data is not even read)
  --no-skip-existing   always (re)plot everything and overwrite the PNGs
Default is set by SKIP_EXISTING below.
"""
import argparse
import glob
import os
import re
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np

# ---- Option: True = skip plots whose PNG already exists, False = always replot.
#      (can be overridden on the command line with --skip-existing / --no-skip-existing)
SKIP_EXISTING = True

parser = argparse.ArgumentParser(description="Plot q and p results.")
parser.add_argument("folder", nargs="?", default=".", help="folder containing the .txt files")
parser.add_argument("--skip-existing", action=argparse.BooleanOptionalAction,
                    default=SKIP_EXISTING,
                    help="skip figures whose PNG already exists")
args = parser.parse_args()
folder = args.folder
skip_existing = args.skip_existing

NAME_RE = re.compile(
    r"^(?P<name>q|p|solExt|solExtFar)Resu_k[=_](?P<k>.+?)_a[=_](?P<a>.+?)"
    r"(?:_theta[=_](?P<theta>.+?))?_N[=_](?P<N>\d+)\.txt$"
)


# Names for which an error-vs-n figure is also produced (add "q", "p" to extend it)
ERROR_NAMES = {"solExt"}

# Names whose points lie on the boundary r = a, where the exact value is known:
# u+(a, theta) = -exp(-i k a cos(theta)).  For these, the error is measured against it.
# (Only valid if the series is evaluated AT the boundary point, not at an exterior point.)
EXACT_BOUNDARY = {"solExt"}

# name used in file names -> label shown on the figures
DISPLAY = {"q": "q", "p": "p", "solExt": "u_N^+"}


def num(s):
    """'0_1' or '0.100000' -> float"""
    return float(s.replace("_", "."))


def fmt(x):
    """Compact number for titles/file names: 5.000000 -> 5, 0.100000 -> 0.1"""
    return f"{x:g}"


def parse_name(path):
    m = NAME_RE.match(os.path.basename(path))
    if not m:
        return None
    return dict(
        name=m["name"],
        k=num(m["k"]),
        a=num(m["a"]),
        N=int(m["N"]),
        theta=num(m["theta"]) if m["theta"] else None,
    )


def try_floats(tokens):
    try:
        return [float(t) for t in tokens]  # float() accepts 'nan' / 'inf'
    except ValueError:
        return None


# ---------------------------------------------------------------- readers
def read_convergence(path):
    """Lines 'n real imag'; any other line (header, ...) is ignored."""
    rows = []
    with open(path) as f:
        for line in f:
            t = line.split()
            v = try_floats(t) if len(t) == 3 else None
            if v:
                rows.append(v)
    a = np.array(rows)
    return a[:, 0], a[:, 1], a[:, 2]


def read_final(path):
    """Header 'k =.. a=..', then records: '(x;y)' line + 'theta real imag' line."""
    rows = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("("):
                continue
            t = line.split()
            v = try_floats(t) if len(t) == 3 else None
            if v:
                rows.append(v)
    a = np.array(rows)
    a = a[np.argsort(a[:, 0])]
    return a[:, 0], a[:, 1], a[:, 2]


# ---------------------------------------------------------------- plots
def plot_convergence(name, key, files):
    k, a, N = key
    tag = f"k={fmt(k)}_a={fmt(a)}_N={N}"
    info = rf"k={fmt(k)}, a={fmt(a)}, N={N}"
    out = os.path.join(folder, f"convergence_{name}_{tag}.png")
    if skip_existing and os.path.exists(out):
        print("skipped (already exists):", out)
        return False
    files = sorted(files, key=lambda fi: fi[1]["theta"])
    fig, axes = plt.subplots(1, 2, figsize=(14, 5), sharex=True)
    for path, meta in files:
        n, re_, im_ = read_convergence(path)
        lab = rf"$\theta$={meta['theta']:.3f}"
        axes[0].plot(n, re_, lw=0.8, label=lab)
        axes[1].plot(n, im_, lw=0.8, label=lab)
        bad = np.where(~np.isfinite(re_) | ~np.isfinite(im_))[0]
        if bad.size:
            print(f"[{name}, {tag}] {os.path.basename(path)}: first NaN at n = {int(n[bad[0]])}")
    for ax, part in zip(axes, ("Real part", "Imaginary part")):
        ax.set_xlabel("n (number of terms)")
        ax.set_ylabel(f"{DISPLAY.get(name, name)}: partial sum ({part})")
        ax.set_title(f"Convergence of {DISPLAY.get(name, name)} - {part}\n{info}")
        ax.grid(alpha=0.3)
    axes[1].legend(fontsize=8, ncol=2)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    print("saved", out)
    return True


def plot_error(name, key, files):
    """Error of the partial sums vs n (log scale), all thetas superposed.
    Left : |S_n - u_exact| with u_exact = -exp(-i k a cos(theta)) for names in EXACT_BOUNDARY,
           otherwise |S_n - S_ref| with S_ref = last finite partial sum.
    Right: |S_n - S_(n-1)|, difference between two successive partial sums.
    Values that are exactly 0 (machine precision reached) are not drawn on the log axis."""
    k, a, N = key
    tag = f"k={fmt(k)}_a={fmt(a)}_N={N}"
    info = rf"k={fmt(k)}, a={fmt(a)}, N={N}"
    label = DISPLAY.get(name, name)
    out = os.path.join(folder, f"error_{name}_{tag}.png")
    if skip_existing and os.path.exists(out):
        print("skipped (already exists):", out)
        return False
    files = sorted(files, key=lambda fi: fi[1]["theta"])
    fig, axes = plt.subplots(1, 2, figsize=(14, 5), sharex=True)
    for path, meta in files:
        n, re_, im_ = read_convergence(path)
        S = re_ + 1j * im_
        good = np.where(np.isfinite(S))[0]
        if good.size == 0:
            print(f"[{name}, {tag}] {os.path.basename(path)}: no finite value, skipped")
            continue
        if name in EXACT_BOUNDARY:
            # exact value on the boundary r = a:  u+(a, theta) = -exp(-i k a cos(theta))
            ref = -np.exp(-1j * k * a * np.cos(meta["theta"]))
            err = np.abs(S - ref)
        else:
            ref_idx = good[-1]                 # no exact solution: last finite partial sum
            if ref_idx != len(S) - 1:
                print(f"[{name}, {tag}] {os.path.basename(path)}: reference taken at n = "
                      f"{int(n[ref_idx])} (NaN afterwards)")
            err = np.abs(S - S[ref_idx])
            err[ref_idx:] = np.nan
        err[err == 0] = np.nan
        diff = np.abs(np.diff(S, prepend=np.nan))
        diff[diff == 0] = np.nan
        lab = rf"$\theta$={meta['theta']:.3f}"
        axes[0].semilogy(n, err, lw=0.9, label=lab)
        axes[1].semilogy(n, diff, lw=0.9, label=lab)
    if name in EXACT_BOUNDARY:
        axes[0].set_title(f"Error of {label} vs n: " + r"$|S_n + e^{-ika\cos\theta}|$" + f"\n{info}")
    else:
        axes[0].set_title(f"Error of {label} vs n: " + r"$|S_n - S_{ref}|$" + f"\n{info}")
    axes[1].set_title(f"Successive differences of {label}: " + r"$|S_n - S_{n-1}|$" + f"\n{info}")
    axes[0].set_ylabel("error (log scale)")
    axes[1].set_ylabel("difference (log scale)")
    for ax in axes:
        ax.set_xlabel("n (number of terms)")
        ax.grid(alpha=0.3, which="both")
    axes[1].legend(fontsize=8, ncol=2)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    print("saved", out)
    return True


def plot_final(name, key, path):
    k, a, N = key
    tag = f"k={fmt(k)}_a={fmt(a)}_N={N}"
    info = rf"k={fmt(k)}, a={fmt(a)}, N={N}"
    out = os.path.join(folder, f"final_{name}_{tag}.png")
    if skip_existing and os.path.exists(out):
        print("skipped (already exists):", out)
        return False
    theta, re_, im_ = read_final(path)
    if not np.isfinite(re_).any():
        print(f"[{name}, {tag}] all final values are NaN - plot will be empty")
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(theta, re_, "o-", label="Real part")
    ax.plot(theta, im_, "s--", label="Imaginary part")
    ax.set_xlabel(r"$\theta$ (rad)")
    ax.set_ylabel(f"{DISPLAY.get(name, name)}(θ)")
    ax.set_title(f"Final approximation of {DISPLAY.get(name, name)}\n{info}")
    ax.grid(alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    print("saved", out)
    return True


if __name__ == "__main__":
    conv = defaultdict(list)   # (name, (k,a,N)) -> [(path, meta)]
    final = {}                 # (name, (k,a,N)) -> path
    for path in glob.glob(os.path.join(folder, "*Resu_*.txt")):
        meta = parse_name(path)
        if meta is None:
            continue
        key = (meta["name"], (meta["k"], meta["a"], meta["N"]))
        if meta["theta"] is None:
            final[key] = path
        else:
            conv[key].append((path, meta))

    if not conv and not final:
        print("No matching {q|p}Resu_k=..._a=..._N=... files found in", os.path.abspath(folder))
    made = 0
    for (name, key), files in sorted(conv.items()):
        made += bool(plot_convergence(name, key, files))
        if name in ERROR_NAMES:
            made += bool(plot_error(name, key, files))
    for (name, key), path in sorted(final.items()):
        made += bool(plot_final(name, key, path))
    if made:
        plt.show()
    else:
        print("Nothing new to plot.")