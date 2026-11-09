import glob, sys
import numpy as np
import matplotlib.pyplot as plt

def lire(motif="svd_blocs_rank*.txt"):
    blocs = []
    for fichier in glob.glob(motif):
        with open(fichier) as f:
            for ligne in f:
                t = ligne.split()
                if not t:
                    continue
                typ, m, n, nl, nc, diag, nb = t[0], *map(int, t[1:7])
                sig = np.array(t[7:7 + nb], dtype=float)
                blocs.append(dict(type=typ, m=m, n=n, niv=(nl, nc), diag=diag, sigma=sig))
    return blocs

def apparier(blocs, nb_paires=20, tol_rel=0.15):
    """Pour chaque bloc admissible (m, n), cherche un bloc dense de taille proche."""
    adm = [b for b in blocs if b["type"] == "ADMISSIBLE"]
    den = [b for b in blocs if b["type"] == "DENSE"]
    paires, utilises = [], set()
    # on prend d'abord les plus gros admissibles
    for a in sorted(adm, key=lambda b: -b["m"] * b["n"]):
        cands = [(abs(d["m"] - a["m"]) / a["m"] + abs(d["n"] - a["n"]) / a["n"], i)
                 for i, d in enumerate(den) if i not in utilises]
        if not cands:
            break
        err, i = min(cands)
        if err < 2 * tol_rel:
            paires.append((a, den[i]))
            utilises.add(i)
        if len(paires) == nb_paires:
            break
    return paires

if __name__ == "__main__":
    blocs = lire(sys.argv[1] if len(sys.argv) > 1 else "svd_blocs_rank*.txt")
    paires = apparier(blocs)
    if not paires:
        sys.exit("Aucune paire admissible/dense de taille comparable trouvée.")

    ncols = min(4, len(paires))                 # max 4 plots per row
    nrows = int(np.ceil(len(paires) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 3.6 * nrows),
                             squeeze=False)
    axes = axes.ravel()

    for ax, (a, d) in zip(axes, paires):
        for b, c, lab in ((a, "tab:blue", "admissible"), (d, "tab:red", "dense")):
            s = b["sigma"] / b["sigma"][0]
            ax.semilogy(np.arange(len(s)), s, "o-", ms=3, color=c,
                        label=f'{lab} {b["m"]}×{b["n"]}')
        ax.set_xlabel("k")
        ax.set_ylabel(r"$\sigma_k/\sigma_0$")
        ax.grid(True, which="both", alpha=.3)
        ax.legend(fontsize=8)

    # hide unused axes in the last row
    for ax in axes[len(paires):]:
        ax.set_visible(False)

    plt.tight_layout()
    plt.savefig("svd_blocs.png", dpi=150)
    plt.show()