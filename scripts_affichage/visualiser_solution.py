from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.tri import Triangulation
import numpy as np


EXPORTS_SOLUTION = {
    "uExact_Re.png": ("u_exact_re", "Partie réelle de la solution exacte", "RdBu_r", True),
    "uExact_Im.png": ("u_exact_im", "Partie imaginaire de la solution exacte", "RdBu_r", True),
    "uExact_Norm.png": ("u_exact_norm", "Norme de la solution exacte", "viridis", False),
    "uCached_Re.png": ("u_cached_re", "Partie réelle de la solution calculée", "RdBu_r", True),
    "uCached_Im.png": ("u_cached_im", "Partie imaginaire de la solution calculée", "RdBu_r", True),
    "uCached_Norm.png": ("u_cached_norm", "Norme de la solution calculée", "viridis", False),
    "Erreur.png": ("erreur", "Erreur absolue", "magma", False),
}


def lire_donnees(fichier):
    donnees = np.loadtxt(fichier, comments="#", ndmin=2)
    if donnees.shape[1] < 7:
        raise ValueError(
            f"{fichier} contient {donnees.shape[1]} colonnes ; 7 sont attendues."
        )
    x, y = donnees[:, 0], donnees[:, 1]
    u_exact = donnees[:, 2] + 1j * donnees[:, 3]
    u_cached = donnees[:, 4] + 1j * donnees[:, 5]
    erreur = donnees[:, 6]
    return x, y, u_exact, u_cached, erreur


def lire_contour(fichier_maillage):
    segments = np.loadtxt(fichier_maillage, comments="#", ndmin=2)
    if segments.shape[1] != 4:
        raise ValueError(f"{fichier_maillage} doit contenir 4 colonnes.")
    p1 = segments[:, 0:2]
    p2 = segments[:, 2:4]
    contour = np.vstack((p1[:1], p2))
    if not np.allclose(contour[-1], contour[0]):
        raise ValueError(
            "Le contour du maillage n'est pas fermé ou les segments ne sont pas ordonnés."
        )
    return contour


def ajouter_zone_obstacle(ax, contour):
    ax.fill(contour[:, 0], contour[:, 1], color="black", linewidth=0, zorder=10)


def enregistrer_figure(
    x, y, valeurs, titre, sortie, triangulation, contour,
    cmap="viridis", symetrique=False,
):
    fig, ax = plt.subplots(figsize=(8, 7), constrained_layout=True)
    valeurs = np.asarray(valeurs)

    if symetrique:
        vmax = float(np.nanmax(np.abs(valeurs)))
        if not np.isfinite(vmax) or vmax == 0:
            vmax = 1.0
        vmin = -vmax
    else:
        vmin = vmax = None

    image = ax.tricontourf(
        triangulation, valeurs, levels=60, cmap=cmap, vmin=vmin, vmax=vmax
    )
    ajouter_zone_obstacle(ax, contour)
    fig.colorbar(image, ax=ax, label=titre)
    ax.set_title(titre)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_aspect("equal", adjustable="box")
    fig.savefig(sortie, dpi=180)
    plt.close(fig)


def enregistrer_animation(x, y, u_cached, k, sortie, nb_frames, fps, contour):
    if len(x) < 3:
        raise ValueError("Au moins trois points sont nécessaires pour trianguler le champ.")

    triangulation = Triangulation(x, y)
    periode = 2.0 * np.pi / abs(k) if k != 0 else 2.0 * np.pi
    temps = np.linspace(0.0, periode, nb_frames, endpoint=False)

    vmax = float(np.max(np.abs(u_cached)))
    if not np.isfinite(vmax) or vmax == 0:
        vmax = 1.0

    champ_initial = u_cached * np.exp(1j * k * temps[0])
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5), constrained_layout=True)

    re_plot = axes[0].tripcolor(
        triangulation, champ_initial.real, shading="gouraud",
        cmap="RdBu_r", vmin=-vmax, vmax=vmax,
    )
    im_plot = axes[1].tripcolor(
        triangulation, champ_initial.imag, shading="gouraud",
        cmap="RdBu_r", vmin=-vmax, vmax=vmax,
    )
    ajouter_zone_obstacle(axes[0], contour)
    ajouter_zone_obstacle(axes[1], contour)
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
        champ = u_cached * np.exp(1j * k * temps[indice])
        re_plot.set_array(champ.real)
        im_plot.set_array(champ.imag)
        titre.set_text(f"Champ calculé dépendant du temps — t = {temps[indice]:.5g}")
        return re_plot, im_plot, titre

    animation = FuncAnimation(
        fig, actualiser, frames=nb_frames,
        interval=1000 / fps, blit=False,
    )
    animation.save(sortie, writer=PillowWriter(fps=fps))
    plt.close(fig)


def visualiser_solution(
    fichier_solution,
    fichier_maillage,
    dossier_sortie,
    exports=None,
    k=None,
    frames=48,
    fps=12,
):
    fichier_solution = Path(fichier_solution)
    fichier_maillage = Path(fichier_maillage)
    dossier_sortie = Path(dossier_sortie)
    dossier_sortie.mkdir(parents=True, exist_ok=True)

    if frames < 2:
        raise ValueError("frames doit être supérieur ou égal à 2.")
    if fps < 1:
        raise ValueError("fps doit être supérieur ou égal à 1.")

    if exports is None:
        exports = list(EXPORTS_SOLUTION.keys()) + ["animation"]
    exports = set(exports)

    x, y, u_exact, u_cached, erreur = lire_donnees(fichier_solution)
    contour = lire_contour(fichier_maillage)
    triangulation = Triangulation(x, y)

    valeurs_disponibles = {
        "u_exact_re": u_exact.real,
        "u_exact_im": u_exact.imag,
        "u_exact_norm": np.abs(u_exact),
        "u_cached_re": u_cached.real,
        "u_cached_im": u_cached.imag,
        "u_cached_norm": np.abs(u_cached),
        "erreur": erreur,
    }

    sorties = []
    for nom, (cle, titre, cmap, symetrique) in EXPORTS_SOLUTION.items():
        if nom not in exports:
            continue
        sortie = dossier_sortie / nom
        enregistrer_figure(
            x, y, valeurs_disponibles[cle], titre, sortie,
            triangulation, contour, cmap=cmap, symetrique=symetrique,
        )
        print(f"Figure créée : {sortie}")
        sorties.append(sortie)

    if "animation" in exports:
        if k is None:
            import re

            # Chercher k dans le nom du fichier, puis dans les dossiers parents.
            chemins_a_tester = [fichier_solution.stem]
            chemins_a_tester.extend(
                parent.name for parent in fichier_solution.parents
            )

            for nom in chemins_a_tester:
                match = re.search(
                    r"(?:^|_)k([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)",
                    nom,
                )
                if match:
                    k = float(match.group(1))
                    break
        if k is None:
            raise ValueError(
                "Impossible de déduire k du nom du fichier. "
                "Passez k=... à visualiser_solution()."
            )

        sortie = dossier_sortie / f"{fichier_solution.stem}_animation.gif"
        enregistrer_animation(
            x, y, u_cached, float(k), sortie, frames, fps, contour
        )
        print(f"Animation créée : {sortie}")
        sorties.append(sortie)

    return sorties
