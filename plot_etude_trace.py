#!/usr/bin/env python3
"""
Trace les resultats de l'etude TP2_etude (temps de construction de A et
erreur relative sur p en fonction du nombre de segments N et du nombre de
points de quadrature nq).

Entree  : outputs/etude_temps_erreur_A.txt   (genere par le programme C++)
Sortie  : outputs/figures/*.png

Usage :
    python3 trace_etude_A.py
    python3 trace_etude_A.py -i outputs/etude_temps_erreur_A.txt -o outputs/figures --show
"""

import argparse
import os

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm, Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.ticker import NullFormatter, ScalarFormatter

# Indices des colonnes du fichier
C_NDEM, C_N, C_PAS, C_NQ, C_TC, C_TE, C_EC, C_EE = range(8)


# ----------------------------------------------------------------------------
# Lecture
# ----------------------------------------------------------------------------
def charge(fichier):
    data = np.loadtxt(fichier, comments="#", ndmin=2)  # "nan" est lu correctement
    if data.shape[1] != 8:
        raise ValueError(f"8 colonnes attendues, {data.shape[1]} trouvees")
    return data


def parametres_entete(fichier):
    """Recupere la ligne '# k=... rayon=...' pour l'afficher dans les titres."""
    with open(fichier, encoding="utf-8") as f:
        for ligne in f:
            if ligne.startswith("# k="):
                return ligne[2:].strip()
    return ""


def grille(data, colonne):
    """Matrice [iN, iNq] de la colonne demandee, + vecteurs N et nq tries."""
    Ns = np.unique(data[:, C_N])
    nqs = np.unique(data[:, C_NQ])
    Z = np.full((len(Ns), len(nqs)), np.nan)
    for ligne in data:
        i = np.searchsorted(Ns, ligne[C_N])
        j = np.searchsorted(nqs, ligne[C_NQ])
        Z[i, j] = ligne[colonne]
    return Ns, nqs, Z


# ----------------------------------------------------------------------------
# Outils de trace
# ----------------------------------------------------------------------------
def courbes_colorees(ax, x, Y, labels_couleur, cmap, titre_cbar, log_cbar=False):
    """Une courbe par ligne de Y, coloree selon labels_couleur (avec colorbar)."""
    vals = np.asarray(labels_couleur, dtype=float)
    norm = (LogNorm if log_cbar else Normalize)(vals.min(), vals.max())
    cm = plt.get_cmap(cmap)
    for y, v in zip(Y, vals):
        ax.plot(x, y, "o-", ms=3, lw=1.2, color=cm(norm(v)))
    sm = ScalarMappable(norm=norm, cmap=cm)
    cb = plt.colorbar(sm, ax=ax)
    cb.set_label(titre_cbar)


def courbes_legende(ax, x, Y, labels, cmap, prefixe):
    """Une courbe par ligne de Y, avec legende."""
    cm = plt.get_cmap(cmap)
    n = len(labels)
    for i, (y, lab) in enumerate(zip(Y, labels)):
        ax.plot(x, y, "o-", ms=3, lw=1.2, color=cm(i / max(n - 1, 1)),
                label=f"{prefixe} = {int(lab)}")
    ax.legend(fontsize=8, loc="center left", bbox_to_anchor=(1.01, 0.5))


def pente_ref(ax, x, y_ref, x0, pente, label):
    """Droite de reference y = y_ref * (x/x0)^pente."""
    ax.plot(x, y_ref * (x / x0) ** pente, "k--", lw=1, alpha=0.6, label=label)


def mise_en_forme(ax, xlabel, ylabel, titre, xlog=True, ylog=True):
    if xlog:
        ax.set_xscale("log")
    if ylog:
        ax.set_yscale("log")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(titre, fontsize=10)
    ax.grid(True, which="both", alpha=0.3)


def ticks_nq(ax, nqs):
    """Graduations entieres lisibles pour un axe nq (meme en echelle log)."""
    pas = max(1, int(round(len(nqs) / 8)))
    ax.set_xticks(nqs[::pas])
    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.xaxis.set_minor_formatter(NullFormatter())


def sauve(fig, dossier, nom, show):
    chemin = os.path.join(dossier, nom)
    fig.tight_layout()
    fig.savefig(chemin, dpi=150)
    print("  ->", chemin)
    if not show:
        plt.close(fig)


# ----------------------------------------------------------------------------
# Figures
# ----------------------------------------------------------------------------
def fig_temps(data, dossier, sous_titre, show):
    Ns, nqs, T = grille(data, C_TC)

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    # (a) temps vs N, une courbe par nq
    ax = axes[0]
    courbes_colorees(ax, Ns, T.T, nqs, "viridis", "nq")
    j = len(nqs) // 2
    ok = np.isfinite(T[:, j])
    if ok.any():
        i0 = np.where(ok)[0][-1]
        pente_ref(ax, Ns, T[i0, j], Ns[i0], 2, r"$\propto N^2$")
        ax.legend(fontsize=8, loc="upper left")
    mise_en_forme(ax, "Nombre de segments N", "Temps de construction de A (s)",
                  "Construction de A (cached) : temps vs N")

    # (b) temps vs nq, une courbe par N
    ax = axes[1]
    courbes_legende(ax, nqs, T, Ns, "plasma", "N")
    i = len(Ns) - 1
    ok = np.isfinite(T[i, :])
    if ok.any():
        j0 = np.where(ok)[0][-1]
        pente_ref(ax, nqs, T[i, j0], nqs[j0], 2, r"$\propto nq^2$")
    mise_en_forme(ax, "Points de quadrature nq", "Temps de construction de A (s)",
                  "Construction de A (cached) : temps vs nq")
    ticks_nq(ax, nqs)

    fig.suptitle(sous_titre, fontsize=8)
    sauve(fig, dossier, "temps_construction_A.png", show)


def fig_erreur(data, dossier, sous_titre, show, col=C_EC, version="cached",
               nom="erreur_p.png"):
    """Erreur relative sur p vs N et vs nq (version 'cached' ou 'exacte')."""
    Ns, nqs, E = grille(data, col)

    # On ne garde que les maillages pour lesquels la colonne est renseignee
    # (la version exacte n'est calculee que pour N <= nbSegmentsMaxExact).
    lignes = np.where(np.isfinite(E).any(axis=1))[0]
    if len(lignes) == 0:
        print(f"  (pas de donnees pour la version {version} : {nom} ignoree)")
        return
    Ns, E = Ns[lignes], E[lignes]

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    # (a) erreur vs N, une courbe par nq
    ax = axes[0]
    courbes_colorees(ax, Ns, E.T, nqs, "viridis", "nq")
    j = len(nqs) // 2
    ok = np.isfinite(E[:, j])
    if ok.any() and len(Ns) > 1:
        i0 = np.where(ok)[0][0]
        pente_ref(ax, Ns, E[i0, j], Ns[i0], -1, r"$\propto N^{-1}$ (ordre 1)")
        pente_ref(ax, Ns, E[i0, j], Ns[i0], -2, r"$\propto N^{-2}$ (ordre 2)")
        ax.legend(fontsize=8, loc="lower left")
    mise_en_forme(ax, "Nombre de segments N",
                  r"Erreur relative $\|p-p_{ana}\|/\|p_{ana}\|$",
                  f"Erreur sur p ({version}) vs N")

    # (b) erreur vs nq, une courbe par N
    ax = axes[1]
    courbes_legende(ax, nqs, E, Ns, "plasma", "N")
    mise_en_forme(ax, "Points de quadrature nq",
                  r"Erreur relative $\|p-p_{ana}\|/\|p_{ana}\|$",
                  f"Erreur sur p ({version}) vs nq", xlog=False)
    ticks_nq(ax, nqs)

    fig.suptitle(sous_titre, fontsize=8)
    sauve(fig, dossier, nom, show)


def fig_cartes(data, dossier, sous_titre, show):
    """Cartes de chaleur dans le plan (N, nq)."""
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    for ax, col, titre, cmap, cbl in (
        (axes[0], C_TC, "Temps de construction de A (s)", "magma", "temps (s)"),
        (axes[1], C_EC, "Erreur relative sur p", "viridis_r", "erreur"),
    ):
        Ns, nqs, Z = grille(data, col)
        Zp = np.where(Z > 0, Z, np.nan)
        if not np.isfinite(Zp).any():
            continue
        im = ax.imshow(Zp, origin="lower", aspect="auto", cmap=cmap,
                       norm=LogNorm(np.nanmin(Zp), np.nanmax(Zp)))
        ax.set_xticks(range(len(nqs)))
        ax.set_xticklabels([int(q) for q in nqs], fontsize=7)
        ax.set_yticks(range(len(Ns)))
        ax.set_yticklabels([int(n) for n in Ns], fontsize=8)
        ax.set_xlabel("Points de quadrature nq")
        ax.set_ylabel("Nombre de segments N")
        ax.set_title(titre, fontsize=10)
        plt.colorbar(im, ax=ax).set_label(cbl)

    fig.suptitle(sous_titre, fontsize=8)
    sauve(fig, dossier, "cartes_temps_erreur.png", show)


def fig_exact_vs_cached(data, dossier, sous_titre, show):
    """Comparaison exacte / cached (seulement la ou la version exacte existe)."""
    Ns, nqs, TE = grille(data, C_TE)
    _, _, TC = grille(data, C_TC)
    _, _, EE = grille(data, C_EE)
    _, _, EC = grille(data, C_EC)

    lignes = np.where(np.isfinite(TE).any(axis=1))[0]
    if len(lignes) == 0:
        print("  (pas de donnees pour la version exacte : figure ignoree)")
        return

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    # (a) speedup = t_exact / t_cached
    ax = axes[0]
    speedup = TE / TC
    courbes_colorees(ax, Ns[lignes], speedup[lignes].T, nqs, "viridis", "nq")
    ax.axhline(1, color="k", lw=0.8, ls=":")
    mise_en_forme(ax, "Nombre de segments N", r"Speedup $t_{exact}/t_{cached}$",
                  "Gain de la version cached", ylog=False)

    # (b) ecart relatif entre les erreurs exacte et cached
    ax = axes[1]
    ecart = np.abs(EC - EE) / EE
    ecart = np.where(ecart > 0, ecart, np.nan)
    courbes_colorees(ax, Ns[lignes], ecart[lignes].T, nqs, "viridis", "nq")
    mise_en_forme(ax, "Nombre de segments N",
                  r"$|err_{cached}-err_{exact}|/err_{exact}$",
                  "Impact du cache sur l'erreur sur p")

    fig.suptitle(sous_titre, fontsize=8)
    sauve(fig, dossier, "exact_vs_cached.png", show)


# ----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-i", "--input", default="outputs/etude_temps_erreur_A.txt")
    ap.add_argument("-o", "--outdir", default="outputs/figures")
    ap.add_argument("--show", action="store_true", help="afficher les figures a l'ecran")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    data = charge(args.input)
    sous_titre = parametres_entete(args.input)

    print(f"{len(data)} lignes lues dans {args.input}")
    print("Figures :")
    fig_temps(data, args.outdir, sous_titre, args.show)
    fig_erreur(data, args.outdir, sous_titre, args.show,
               col=C_EC, version="cached", nom="erreur_p.png")
    fig_erreur(data, args.outdir, sous_titre, args.show,
               col=C_EE, version="exacte", nom="erreur_p_exacte.png")
    fig_cartes(data, args.outdir, sous_titre, args.show)
    fig_exact_vs_cached(data, args.outdir, sous_titre, args.show)

    if args.show:
        plt.show()


if __name__ == "__main__":
    main()