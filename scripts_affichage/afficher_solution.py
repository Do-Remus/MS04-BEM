
#!/usr/bin/env python3
"""
Affiche et anime les solutions complexes exportées par ecrire_solution().

Format attendu (7 colonnes) :
x y uExact_re uExact_im uCached_re uCached_im erreur

Le script recherche par défaut le fichier outputs/u_*.txt le plus récent.
Les figures et le GIF sont enregistrés dans un sous-dossier de outputs/.

Dépendances :
    python -m pip install numpy matplotlib pillow

Exécution depuis la racine du projet :
    python scripts_affichage/afficher_solution.py
    python scripts_affichage/afficher_solution.py --input outputs/u_k10.000000.txt
"""

from __future__ import annotations

import argparse
import glob
import re
from pathlib import Path
from typing import Optional

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.tri import Triangulation
from matplotlib.path import Path as MplPath
from matplotlib.patches import PathPatch
from matplotlib.collections import LineCollection
import numpy as np



def trouver_fichier(
    racine: Path, fichier_argument: Optional[str]
) -> Path:
    """Trouve le fichier demandé ou le fichier u_*.txt le plus récent."""
    if fichier_argument:
        chemin = Path(fichier_argument)
        if not chemin.is_absolute():
            chemin = racine / chemin

        if not chemin.is_file():
            raise FileNotFoundError(f"Fichier introuvable : {chemin}")

        return chemin

    motif = str(racine / "outputs" / "u_*.txt")
    candidats = [Path(p) for p in glob.glob(motif)]

    if not candidats:
        raise FileNotFoundError(
            f"Aucun fichier outputs/u_*.txt trouvé depuis {racine}."
        )

    return max(candidats, key=lambda p: p.stat().st_mtime)


def lire_donnees(fichier: Path):
    """Lit les sept colonnes produites par ecrire_solution()."""
    donnees = np.loadtxt(fichier, comments="#", ndmin=2)

    if donnees.shape[1] < 7:
        raise ValueError(
            f"{fichier} contient {donnees.shape[1]} colonnes ; "
            "7 sont attendues."
        )

    x = donnees[:, 0]
    y = donnees[:, 1]

    u_exact = donnees[:, 2] + 1j * donnees[:, 3]
    u_cached = donnees[:, 4] + 1j * donnees[:, 5]
    erreur = donnees[:, 6]

    return x, y, u_exact, u_cached, erreur


def extraire_k(nom_fichier: str, k_argument: Optional[float]) -> float:
    """Extrait k d'un nom comme u_k10.000000_R1.000000_....txt."""
    if k_argument is not None:
        return k_argument

    match = re.search(
        r"_k([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)",
        nom_fichier,
    )

    if match:
        return float(match.group(1))

    raise ValueError(
        "Impossible de déduire k du nom du fichier. "
        "Fournissez --k VALEUR."
    )

def lire_maillage(fichier: Path):
    """
    Lit outputs/maillage.txt.
    Chaque ligne contient : x1 y1 x2 y2.
    """
    segments = np.loadtxt(fichier, comments="#", ndmin=2)

    if segments.shape[1] != 4:
        raise ValueError(
            f"{fichier} doit contenir exactement 4 colonnes."
        )

    p1 = segments[:, 0:2]
    p2 = segments[:, 2:4]

    return segments, p1, p2

def construire_contour(p1, p2):
    """
    Reconstruit le contour en supposant que les segments sont
    déjà rangés dans l'ordre et raccordés les uns aux autres.
    """
    points = np.vstack((p1[:1], p2))

    if not np.allclose(points[-1], points[0]):
        raise ValueError(
            "Le contour n'est pas fermé ou les segments ne sont "
            "pas dans l'ordre."
        )

    return points


def ajouter_zone_noire(ax, contour):
    """Remplit l'intérieur du contour en noir."""
    ax.fill(
        contour[:, 0],
        contour[:, 1],
        color="black",
        linewidth=0,
        zorder=10,
    )

def enregistrer_figure(
    x, y, valeurs, titre, sortie, triangulation,
    contour, cmap="viridis", symetrique=False
):
    """Enregistre un champ dans une image indépendante."""
    fig, ax = plt.subplots(figsize=(8, 7), constrained_layout=True)

    valeurs = np.asarray(valeurs)

    if symetrique:
        vmax = float(np.nanmax(np.abs(valeurs)))
        if not np.isfinite(vmax) or vmax == 0:
            vmax = 1.0
        vmin = -vmax
    else:
        vmin = None
        vmax = None

    image = ax.tricontourf(
        triangulation,
        valeurs,
        levels=60,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
    )

    ajouter_zone_noire(ax, contour)

    fig.colorbar(image, ax=ax, label=titre)
    ax.set_title(titre)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_aspect("equal", adjustable="box")

    fig.savefig(sortie, dpi=180)
    plt.close(fig)


def enregistrer_animation(
    x, y, u_cached, k, sortie, nb_frames, fps, contour
):
    """
    Anime les parties réelle et imaginaire de :

        u(x,t) = u_cached(x) * exp(-i*k*t)

    On suppose une vitesse d'onde c = 1, donc omega = k.
    Les couleurs restent fixes pendant toute l'animation.
    """
    triangulation = Triangulation(x, y)

    periode = 2.0 * np.pi / abs(k) if k != 0 else 2.0 * np.pi
    temps = np.linspace(
        0.0, periode, nb_frames, endpoint=False
    )

    vmax = float(np.max(np.abs(u_cached)))
    if not np.isfinite(vmax) or vmax == 0:
        vmax = 1.0

    champ_initial = u_cached * np.exp(-1j * k * temps[0])

    fig, axes = plt.subplots(
        1, 2, figsize=(13, 5.5), constrained_layout=True
    )

    # tripcolor avec shading="gouraud" interpole les valeurs
    # entre les sommets des triangles pour un rendu continu.
    re_plot = axes[0].tripcolor(
        triangulation,
        champ_initial.real,
        shading="gouraud",
        cmap="RdBu_r",
        vmin=-vmax,
        vmax=vmax,
    )

    im_plot = axes[1].tripcolor(
        triangulation,
        champ_initial.imag,
        shading="gouraud",
        cmap="RdBu_r",
        vmin=-vmax,
        vmax=vmax,
    )

    ajouter_zone_noire(axes[0], contour)
    ajouter_zone_noire(axes[1], contour)

    fig.colorbar(re_plot, ax=axes[0], label="Partie réelle")
    fig.colorbar(im_plot, ax=axes[1], label="Partie imaginaire")

    axes[0].set_title(r"$\Re(u_{\mathrm{cached}}e^{-ikt})$")
    axes[1].set_title(r"$\Im(u_{\mathrm{cached}}e^{-ikt})$")

    for ax in axes:
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_aspect("equal", adjustable="box")

    titre = fig.suptitle("")

    def actualiser(indice):
        phase = np.exp(-1j * k * temps[indice])
        champ = u_cached * phase

        re_plot.set_array(champ.real)
        im_plot.set_array(champ.imag)

        titre.set_text(
            f"Champ cached dépendant du temps — "
            f"t = {temps[indice]:.5g}"
        )

        return re_plot, im_plot, titre

    animation = FuncAnimation(
        fig,
        actualiser,
        frames=nb_frames,
        interval=1000 / fps,
        blit=False,
    )

    animation.save(sortie, writer=PillowWriter(fps=fps))
    plt.close(fig)



def main():
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument(
        "--input",
        help="Chemin du fichier de solution .txt",
    )
    parser.add_argument(
        "--k",
        type=float,
        help="Nombre d'onde k (déduit du nom sinon)",
    )
    parser.add_argument(
        "--frames",
        type=int,
        default=48,
        help="Nombre d'images GIF (défaut : 48)",
    )
    parser.add_argument(
        "--fps",
        type=int,
        default=12,
        help="Images par seconde (défaut : 12)",
    )

    args = parser.parse_args()

    if args.frames < 2:
        parser.error("--frames doit être supérieur ou égal à 2.")
    if args.fps < 1:
        parser.error("--fps doit être supérieur ou égal à 1.")

    def etape(numero, total, message):
        print(f"\n[{numero}/{total}] {message}", flush=True)

    total_etapes = 5

    # Étape 1 : préparation des chemins
    etape(1, total_etapes, "Préparation des dossiers")

    racine = Path(__file__).resolve().parent.parent
    dossier_outputs = racine / "outputs"
    dossier_outputs.mkdir(parents=True, exist_ok=True)

    # Étape 2 : lecture des données
    etape(2, total_etapes, "Lecture du fichier de solution")

    fichier = trouver_fichier(racine, args.input)
    print(f"    Fichier : {fichier}", flush=True)

    x, y, u_exact, u_cached, erreur = lire_donnees(fichier)
    k = extraire_k(fichier.name, args.k)

    print(f"    Nombre de points : {len(x)}", flush=True)
    print(f"    Nombre d'onde k : {k:g}", flush=True)

    # Lecture du maillage
    print("\n    Lecture du maillage...", flush=True)

    fichier_maillage = dossier_outputs / "maillage.txt"
    if not fichier_maillage.is_file():
        raise FileNotFoundError(
            f"Maillage introuvable : {fichier_maillage}"
        )

    segments, p1, p2 = lire_maillage(fichier_maillage)
    contour = construire_contour(p1, p2)
    print(f"    Nombre de segments : {len(segments)}", flush=True)

    # Préparation des sorties
    dossier_sortie = (
        dossier_outputs / f"{fichier.stem}_visualisations"
    )
    dossier_sortie.mkdir(parents=True, exist_ok=True)

    print(f"    Dossier de sortie : {dossier_sortie}", flush=True)

    triangulation = Triangulation(x, y)

    champs = [
        (
            u_exact.real,
            r"Partie réelle de $u_{\mathrm{exact}}$",
            "uExact_Re.png",
            "RdBu_r",
            True,
        ),
        (
            u_exact.imag,
            r"Partie imaginaire de $u_{\mathrm{exact}}$",
            "uExact_Im.png",
            "RdBu_r",
            True,
        ),
        (
            np.abs(u_exact),
            r"Norme de $u_{\mathrm{exact}}$",
            "uExact_Norm.png",
            "viridis",
            False,
        ),
        (
            u_cached.real,
            r"Partie réelle de $u_{\mathrm{cached}}$",
            "uCached_Re.png",
            "RdBu_r",
            True,
        ),
        (
            u_cached.imag,
            r"Partie imaginaire de $u_{\mathrm{cached}}$",
            "uCached_Im.png",
            "RdBu_r",
            True,
        ),
        (
            np.abs(u_cached),
            r"Norme de $u_{\mathrm{cached}}$",
            "uCached_Norm.png",
            "viridis",
            False,
        ),
        (
            erreur,
            r"Erreur absolue",
            "Erreur.png",
            "magma",
            False,
        ),
    ]

    # Étape 3 : génération des images
    etape(3, total_etapes, "Génération des figures PNG")

    nombre_figures = len(champs)

    for indice, (valeurs, titre, nom, cmap, symetrique) in enumerate(
        champs, start=1
    ):
        print(
            f"    Figure {indice}/{nombre_figures} : {nom}",
            flush=True,
        )

        sortie = dossier_sortie / nom

        enregistrer_figure(
            x,
            y,
            valeurs,
            titre,
            sortie,
            triangulation,
            contour,
            cmap=cmap,
            symetrique=symetrique,
        )

        print(f"        Enregistrée : {sortie}", flush=True)

    # Étape 4 : génération du GIF
    etape(4, total_etapes, "Génération de l'animation GIF")

    sortie_gif = (
        dossier_sortie / f"{fichier.stem}_animation.gif"
    )

    print(f"    Nombre d'images : {args.frames}", flush=True)
    print(f"    Fréquence : {args.fps} images/s", flush=True)
    print("    Calcul des images et écriture du GIF...", flush=True)

    enregistrer_animation(
        x,
        y,
        u_cached,
        k,
        sortie_gif,
        args.frames,
        args.fps,
        contour
    )

    print(f"    GIF enregistré : {sortie_gif}", flush=True)

    # Étape 5 : bilan
    etape(5, total_etapes, "Traitement terminé")

    print(f"    Fichier de solution : {fichier}", flush=True)
    print(f"    Dossier des résultats : {dossier_sortie}", flush=True)
    print(f"    Nombre de points : {len(x)}", flush=True)
    print(f"    Nombre de figures : {nombre_figures}", flush=True)
    print(f"    Nombre d'images GIF : {args.frames}", flush=True)
    print(f"    Nombre d'onde : k = {k:g}", flush=True)
    print("    Convention temporelle : exp(-i*k*t), avec c = 1.",
          flush=True)


if __name__ == "__main__":
    main()

