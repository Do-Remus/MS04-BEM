#!/usr/bin/env python3
"""
Trace les resultats de l'etude TP2_etude_CG : temps du gradient conjugue en
fonction du nombre de segments N et du nombre de points de quadrature nq.

Entree  : outputs/etude_temps_CG.txt   (genere par le programme C++)
Sortie  : outputs/figures/cg_*.png

Usage :
    python3 trace_etude_CG.py
    python3 trace_etude_CG.py -i outputs/etude_temps_CG.txt -o outputs/figures --show
"""

import argparse
import os

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LogNorm, Normalize
from matplotlib.ticker import NullFormatter, ScalarFormatter

# Colonnes : nbSegments_demande nbSegments_reel pasMaillage nq
#            temps_assemblage temps_CG erreur_p
C_NDEM, C_N, C_PAS, C_NQ, C_TASS, C_TCG, C_ERR = range(7)


# ----------------------------------------------------------------------------
# Lecture
# ----------------------------------------------------------------------------
def charge(fichier):
    data = np.loadtxt(fichier, comments="#", ndmin=2)
    if data.shape[1] != 7:
        raise ValueError(f"7 colonnes attendues, {data.shape[1]} trouvees")
    return data


def parametres_entete(fichier):
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
        Z[np.searchsorted(Ns, ligne[C_N]), np.searchsorted(nqs, ligne[C_NQ])] = ligne[colonne]
    return Ns, nqs, Z


# ----------------------------------------------------------------------------
# Outils de trace
# ----------------------------------------------------------------------------
def courbes_colorees(ax, x, Y, valeurs, cmap, titre_cbar):
    """Une courbe par ligne de Y, coloree selon `valeurs` (avec colorbar)."""
    vals = np.asarray(valeurs, dtype=float)
    norm = Normalize(vals.min(), vals.max())
    cm = plt.get_cmap(cmap)
    for y, v in zip(Y, vals):
        ax.plot(x, y, "o-", ms=3, lw=1.2, color=cm(norm(v)))
    plt.colorbar(ScalarMappable(norm=norm, cmap=cm), ax=ax).set_label(titre_cbar)


def courbes_legende(ax, x, Y, labels, cmap, prefixe):
    """Une courbe par ligne de Y, avec legende a l'exterieur du graphique."""
    cm = plt.get_cmap(cmap)
    n = len(labels)
    for i, (y, lab) in enumerate(zip(Y, labels)):
        ax.plot(x, y, "o-", ms=3, lw=1.2, color=cm(i / max(n - 1, 1)),
                label=f"{prefixe} = {int(lab)}")
    ax.legend(fontsize=8, loc="center left", bbox_to_anchor=(1.01, 0.5))


def pente_ref(ax, x, y_ref, x0, pente, label):
    ax.plot(x, y_ref * (x / x0) ** pente, "k--", lw=1, alpha=0.6, label=label)


def derniere_valeur_valide(v):
    ok = np.where(np.isfinite(v) & (v > 0))[0]
    return ok[-1] if len(ok) else None


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
def fig_cg_temps(data, dossier, sous_titre, show):
    """Temps du CG vs N (une courbe par nq) et vs nq (une courbe par N)."""
    Ns, nqs, T = grille(data, C_TCG)
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    # (a) temps vs N
    ax = axes[0]
    courbes_colorees(ax, Ns, T.T, nqs, "viridis", "nq")
    j = len(nqs) // 2
    i0 = derniere_valeur_valide(T[:, j])
    if i0 is not None:
        pente_ref(ax, Ns, T[i0, j], Ns[i0], 2, r"$\propto N^2$ (1 produit matrice-vecteur)")
        pente_ref(ax, Ns, T[i0, j], Ns[i0], 3, r"$\propto N^3$")
        ax.legend(fontsize=8, loc="upper left")
    mise_en_forme(ax, "Nombre de segments N", "Temps du gradient conjugué (s)",
                  "Gradient conjugué : temps vs N")

    # (b) temps vs nq
    ax = axes[1]
    courbes_legende(ax, nqs, T, Ns, "plasma", "N")
    mise_en_forme(ax, "Points de quadrature nq", "Temps du gradient conjugué (s)",
                  "Gradient conjugué : temps vs nq", xlog=False)
    ticks_nq(ax, nqs)

    fig.suptitle(sous_titre, fontsize=8)
    sauve(fig, dossier, "cg_temps.png", show)


def fig_cg_cartes(data, dossier, sous_titre, show):
    """Cartes de chaleur : temps du CG, et part du CG dans (assemblage + CG)."""
    Ns, nqs, TCG = grille(data, C_TCG)
    _, _, TASS = grille(data, C_TASS)
    part = TCG / (TCG + TASS)

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    for ax, Z, titre, cmap, cbl, log in (
        (axes[0], TCG, "Temps du gradient conjugué (s)", "magma", "temps (s)", True),
        (axes[1], part, "Part du CG dans (assemblage + CG)", "viridis", "t_CG / (t_ass + t_CG)", False),
    ):
        Zp = np.where(Z > 0, Z, np.nan)
        if not np.isfinite(Zp).any():
            continue
        norm = LogNorm(np.nanmin(Zp), np.nanmax(Zp)) if log else Normalize(0, 1)
        im = ax.imshow(Zp, origin="lower", aspect="auto", cmap=cmap, norm=norm)
        ax.set_xticks(range(len(nqs)))
        ax.set_xticklabels([int(q) for q in nqs], fontsize=7)
        ax.set_yticks(range(len(Ns)))
        ax.set_yticklabels([int(n) for n in Ns], fontsize=8)
        ax.set_xlabel("Points de quadrature nq")
        ax.set_ylabel("Nombre de segments N")
        ax.set_title(titre, fontsize=10)
        plt.colorbar(im, ax=ax).set_label(cbl)

    fig.suptitle(sous_titre, fontsize=8)
    sauve(fig, dossier, "cg_cartes.png", show)


def fig_cg_vs_assemblage(data, dossier, sous_titre, show):
    """Compare le temps du CG au temps d'assemblage de A, pour quelques nq."""
    Ns, nqs, TCG = grille(data, C_TCG)
    _, _, TASS = grille(data, C_TASS)

    # trois valeurs de nq : min, milieu, max
    choix = sorted({0, len(nqs) // 2, len(nqs) - 1})
    cm = plt.get_cmap("viridis")

    fig, ax = plt.subplots(figsize=(7.5, 5))
    for c, j in enumerate(choix):
        couleur = cm(c / max(len(choix) - 1, 1))
        ax.plot(Ns, TASS[:, j], "o-", ms=3, color=couleur, label=f"assemblage, nq = {int(nqs[j])}")
        ax.plot(Ns, TCG[:, j], "s--", ms=3, color=couleur, label=f"CG, nq = {int(nqs[j])}")
    mise_en_forme(ax, "Nombre de segments N", "Temps (s)",
                  "Assemblage de A vs gradient conjugué")
    ax.legend(fontsize=8)

    fig.suptitle(sous_titre, fontsize=8)
    sauve(fig, dossier, "cg_vs_assemblage.png", show)


def fig_cg_erreur(data, dossier, sous_titre, show):
    """Erreur sur p : verifie que le CG a bien converge pour chaque (N, nq)."""
    Ns, nqs, E = grille(data, C_ERR)
    fig, ax = plt.subplots(figsize=(8.5, 5))
    courbes_legende(ax, nqs, E, Ns, "plasma", "N")
    mise_en_forme(ax, "Points de quadrature nq",
                  r"Erreur relative $\|p-p_{ana}\|/\|p_{ana}\|$",
                  "Erreur sur p (controle de la convergence du CG)", xlog=False)
    ticks_nq(ax, nqs)

    fig.suptitle(sous_titre, fontsize=8)
    sauve(fig, dossier, "cg_erreur.png", show)


# ----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-i", "--input", default="outputs/etude_temps_CG.txt")
    ap.add_argument("-o", "--outdir", default="outputs/figures")
    ap.add_argument("--show", action="store_true", help="afficher les figures a l'ecran")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    data = charge(args.input)
    sous_titre = parametres_entete(args.input)

    print(f"{len(data)} lignes lues dans {args.input}")
    print("Figures :")
    fig_cg_temps(data, args.outdir, sous_titre, args.show)
    fig_cg_cartes(data, args.outdir, sous_titre, args.show)
    fig_cg_vs_assemblage(data, args.outdir, sous_titre, args.show)
    fig_cg_erreur(data, args.outdir, sous_titre, args.show)

    if args.show:
        plt.show()


if __name__ == "__main__":
    main()