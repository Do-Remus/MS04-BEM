from pathlib import Path
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.tri import Triangulation
import numpy as np


EXPORTS_SOLUTION = {
    "uExact_Re.png": (
        "u_exact_re",
        "Partie réelle de la solution exacte",
        "RdBu_r",
        True,
    ),
    "uExact_Im.png": (
        "u_exact_im",
        "Partie imaginaire de la solution exacte",
        "RdBu_r",
        True,
    ),
    "uExact_Norm.png": (
        "u_exact_norm",
        "Norme de la solution exacte",
        "viridis",
        False,
    ),
    "uCached_Re.png": (
        "u_cached_re",
        "Partie réelle de la solution calculée",
        "RdBu_r",
        True,
    ),
    "uCached_Im.png": (
        "u_cached_im",
        "Partie imaginaire de la solution calculée",
        "RdBu_r",
        True,
    ),
    "uCached_Norm.png": (
        "u_cached_norm",
        "Norme de la solution calculée",
        "viridis",
        False,
    ),
    "uExactTot_Re.png": (
        "u_exact_tot_re",
        "Partie réelle du champ total exact",
        "RdBu_r",
        True,
    ),
    "uExactTot_Im.png": (
        "u_exact_tot_im",
        "Partie imaginaire du champ total exact",
        "RdBu_r",
        True,
    ),
    "uExactTot_Norm.png": (
        "u_exact_tot_norm",
        "Norme du champ total exact",
        "viridis",
        False,
    ),
    "uCachedTot_Re.png": (
        "u_cached_tot_re",
        "Partie réelle du champ total calculé",
        "RdBu_r",
        True,
    ),
    "uCachedTot_Im.png": (
        "u_cached_tot_im",
        "Partie imaginaire du champ total calculé",
        "RdBu_r",
        True,
    ),
    "uCachedTot_Norm.png": (
        "u_cached_tot_norm",
        "Norme du champ total calculé",
        "viridis",
        False,
    ),
    "Erreur.png": (
        "erreur",
        "Erreur absolue",
        "magma",
        False,
    ),
}


EXPORTS_TOTAUX = {
    "uExactTot_Re.png",
    "uExactTot_Im.png",
    "uExactTot_Norm.png",
    "uCachedTot_Re.png",
    "uCachedTot_Im.png",
    "uCachedTot_Norm.png",
}


def lire_donnees(fichier):
    donnees = np.loadtxt(fichier, comments="#", ndmin=2)

    if donnees.shape[1] < 7:
        raise ValueError(
            f"{fichier} contient {donnees.shape[1]} colonnes ; "
            "7 sont attendues."
        )

    x, y = donnees[:, 0], donnees[:, 1]
    u_exact = donnees[:, 2] + 1j * donnees[:, 3]
    u_cached = donnees[:, 4] + 1j * donnees[:, 5]
    erreur = donnees[:, 6]

    return x, y, u_exact, u_cached, erreur


def trouver_nombre_onde(fichier_solution):
    """Recherche k dans le nom du fichier et dans ses dossiers parents."""
    fichier_solution = Path(fichier_solution)

    noms = [fichier_solution.stem]
    noms.extend(parent.name for parent in fichier_solution.parents)

    motif = re.compile(
        r"(?:^|_)k([-+]?(?:\d+(?:\.\d*)?|\.\d+)"
        r"(?:[eE][-+]?\d+)?)(?=_|$)"
    )

    for nom in noms:
        match = motif.search(nom)

        if match:
            k = float(match.group(1))

            if k <= 0:
                raise ValueError(
                    f"Le nombre d'onde détecté doit être positif : {k}."
                )

            print(f"Nombre d'onde détecté : k = {k} (dans « {nom} »)")
            return k

    return None


def lire_contour(fichier_maillage):
    segments = np.loadtxt(fichier_maillage, comments="#", ndmin=2)

    if segments.shape[1] != 4:
        raise ValueError(
            f"{fichier_maillage} doit contenir 4 colonnes."
        )

    p1 = segments[:, 0:2]
    p2 = segments[:, 2:4]
    contour = np.vstack((p1[:1], p2))

    if not np.allclose(contour[-1], contour[0]):
        raise ValueError(
            "Le contour du maillage n'est pas fermé ou les segments "
            "ne sont pas ordonnés."
        )

    return contour


def ajouter_zone_obstacle(ax, contour):
    ax.fill(
        contour[:, 0],
        contour[:, 1],
        color="black",
        linewidth=0,
        zorder=10,
    )


def enregistrer_figure(
    x,
    y,
    valeurs,
    titre,
    sortie,
    triangulation,
    contour,
    cmap="viridis",
    symetrique=False,
):
    fig, ax = plt.subplots(
        figsize=(8, 7),
        constrained_layout=True,
    )

    valeurs = np.asarray(valeurs)

    if symetrique:
        vmax = float(np.nanmax(np.abs(valeurs)))

        if not np.isfinite(vmax) or vmax == 0:
            vmax = 1.0

        vmin = -vmax
    else:
        vmin = vmax = None

    image = ax.tricontourf(
        triangulation,
        valeurs,
        levels=60,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
    )

    ajouter_zone_obstacle(ax, contour)

    fig.colorbar(image, ax=ax, label=titre)
    ax.set_title(titre)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_aspect("equal", adjustable="box")

    fig.savefig(sortie, dpi=180)
    plt.close(fig)


def enregistrer_animation(
    x,
    y,
    u_cached_tot,
    omega,
    sortie,
    nb_frames,
    fps,
    contour,
):
    if len(x) < 3:
        raise ValueError(
            "Au moins trois points sont nécessaires pour trianguler le champ."
        )

    triangulation = Triangulation(x, y)

    periode = 2.0 * np.pi / abs(omega) if omega != 0 else 2.0 * np.pi
    temps = np.linspace(0.0, periode, nb_frames, endpoint=False)

    vmax = float(np.max(np.abs(u_cached_tot)))

    if not np.isfinite(vmax) or vmax == 0:
        vmax = 1.0

    # Convention temporelle : exp(+i omega t).
    champ_initial = u_cached_tot * np.exp(1j * omega * temps[0])

    fig, axes = plt.subplots(
        1,
        2,
        figsize=(13, 5.5),
        constrained_layout=True,
    )

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

    ajouter_zone_obstacle(axes[0], contour)
    ajouter_zone_obstacle(axes[1], contour)

    fig.colorbar(re_plot, ax=axes[0], label="Partie réelle")
    fig.colorbar(im_plot, ax=axes[1], label="Partie imaginaire")

    axes[0].set_title(r"$\Re(u_{\mathrm{tot}} e^{i\omega t})$")
    axes[1].set_title(r"$\Im(u_{\mathrm{tot}} e^{i\omega t})$")

    for ax in axes:
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_aspect("equal", adjustable="box")

    titre = fig.suptitle("")

    def actualiser(indice):
        champ = u_cached_tot * np.exp(1j * omega * temps[indice])

        re_plot.set_array(champ.real)
        im_plot.set_array(champ.imag)

        titre.set_text(
            f"Champ total calculé — "
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

    animation.save(
        sortie,
        writer=PillowWriter(fps=fps),
    )

    plt.close(fig)


def visualiser_solution(
    fichier_solution,
    fichier_maillage,
    dossier_sortie,
    exports=None,
    k=None,
    omega=None,
    frames=48,
    fps=12,
):
    fichier_solution = Path(fichier_solution)
    fichier_maillage = Path(fichier_maillage)
    dossier_sortie = Path(dossier_sortie)

    dossier_sortie.mkdir(
        parents=True,
        exist_ok=True,
    )

    if frames < 2:
        raise ValueError("frames doit être supérieur ou égal à 2.")

    if fps < 1:
        raise ValueError("fps doit être supérieur ou égal à 1.")

    if exports is None:
        exports = list(EXPORTS_SOLUTION.keys()) + ["animation"]

    exports = set(exports)

    besoin_champ_total = bool(exports & EXPORTS_TOTAUX)
    besoin_animation = "animation" in exports
    besoin_de_k = besoin_champ_total or besoin_animation

    # Trouver k dans le nom du fichier ou de ses dossiers parents.
    if k is None and besoin_de_k:
        k = trouver_nombre_onde(fichier_solution)

    if besoin_de_k and k is None:
        raise ValueError(
            "Impossible de trouver k dans le nom du fichier de solution "
            "ou de ses dossiers parents. Passez k=... à "
            "visualiser_solution()."
        )

    if k is not None:
        k = float(k)

        if k <= 0:
            raise ValueError("Le nombre d'onde k doit être positif.")

    # omega est indépendant de k en général.
    # Si aucune valeur n'est fournie, on suppose une vitesse de propagation
    # normalisée à 1, donc omega = k.
    if omega is None and besoin_animation:
        omega = k

    if omega is not None:
        omega = float(omega)

        if omega <= 0:
            raise ValueError("La pulsation omega doit être positive.")

    x, y, u_exact, u_cached, erreur = lire_donnees(
        fichier_solution
    )

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

    # Champ incident : u_inc(x, y) = exp(-i k x).
    # On suppose ici que les champs lus sont les champs diffusés.
    u_inc = None
    u_cached_tot = None

    if besoin_champ_total or besoin_animation:
        u_inc = np.exp(-1j * k * x)
        u_cached_tot = u_cached + u_inc

    if besoin_champ_total:
        u_exact_tot = u_exact + u_inc

        valeurs_disponibles.update({
            "u_exact_tot_re": u_exact_tot.real,
            "u_exact_tot_im": u_exact_tot.imag,
            "u_exact_tot_norm": np.abs(u_exact_tot),
            "u_cached_tot_re": u_cached_tot.real,
            "u_cached_tot_im": u_cached_tot.imag,
            "u_cached_tot_norm": np.abs(u_cached_tot),
        })

    sorties = []

    for nom, (cle, titre, cmap, symetrique) in EXPORTS_SOLUTION.items():
        if nom not in exports:
            continue

        sortie = dossier_sortie / nom

        enregistrer_figure(
            x,
            y,
            valeurs_disponibles[cle],
            titre,
            sortie,
            triangulation,
            contour,
            cmap=cmap,
            symetrique=symetrique,
        )

        print(f"Figure créée : {sortie}")
        sorties.append(sortie)

    if besoin_animation:
        sortie = (
            dossier_sortie
            / f"{fichier_solution.stem}_animation.gif"
        )

        enregistrer_animation(
            x,
            y,
            u_cached_tot,
            omega,
            sortie,
            frames,
            fps,
            contour,
        )

        print(f"Animation créée : {sortie}")
        sorties.append(sortie)

    return sorties