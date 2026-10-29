from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib import colormaps
from matplotlib.patches import Rectangle

from io_utils import lire_csv, lire_milieux_segments, lire_points


def visualiser_clusters(
    fichier_maillage,
    fichier_clusters,
    type_maillage="segments",
    dossier_sortie=None,
    fps=1,
):
    fichier_maillage = Path(fichier_maillage)
    fichier_clusters = Path(fichier_clusters)
    dossier_sortie = Path(dossier_sortie or fichier_clusters.parent / "visualisations")
    dossier_sortie.mkdir(parents=True, exist_ok=True)

    print(f"Lecture des clusters : {fichier_clusters}", flush=True)

    if type_maillage == "segments":
        points = lire_milieux_segments(fichier_maillage)
    elif type_maillage == "points":
        points = lire_points(fichier_maillage)
    else:
        raise ValueError("type_maillage doit être 'segments' ou 'points'.")

    arbre = lire_csv(
        fichier_clusters,
        {"niveau", "taille", "x_min", "x_max", "y_min", "y_max"},
    )

    idx_racine = arbre["niveau"].idxmin()
    taille_racine = int(arbre.loc[idx_racine, "taille"])
    if len(points) != taille_racine:
        raise ValueError(
            f"Le maillage contient {len(points)} points, mais la racine "
            f"contient {taille_racine} indices. Vérifie que les fichiers "
            "proviennent du même calcul."
        )

    x = points["x"].to_numpy()
    y = points["y"].to_numpy()
    arbre["niveau"] = arbre["niveau"].astype(int)
    niveaux = sorted(arbre["niveau"].unique())
    niveau_max = max(niveaux)

    xmin = min(float(points["x"].min()), float(arbre["x_min"].min()))
    xmax = max(float(points["x"].max()), float(arbre["x_max"].max()))
    ymin = min(float(points["y"].min()), float(arbre["y_min"].min()))
    ymax = max(float(points["y"].max()), float(arbre["y_max"].max()))
    marge = 0.05 * max(xmax - xmin, ymax - ymin, 1.0)

    fig, ax = plt.subplots(figsize=(9, 9))

    def mise_a_jour(niveau_courant):
        ax.clear()
        ax.scatter(
            x, y, s=12, color="0.78",
            label="Milieux des segments" if type_maillage == "segments" else "Points",
            zorder=1,
        )

        visibles = arbre[arbre["niveau"] <= niveau_courant]
        for _, noeud in visibles.iterrows():
            niveau = int(noeud["niveau"])
            couleur = colormaps["viridis"](
                0.12 + 0.72 * niveau / max(niveau_max, 1)
            )
            ax.add_patch(Rectangle(
                (noeud["x_min"], noeud["y_min"]),
                noeud["x_max"] - noeud["x_min"],
                noeud["y_max"] - noeud["y_min"],
                fill=False,
                edgecolor=couleur,
                linewidth=1.7,
                alpha=0.95,
                zorder=3,
            ))

        nb = int((arbre["niveau"] == niveau_courant).sum())
        ax.set_title(
            f"Arbre de clusters — niveau {niveau_courant}\n"
            f"{nb} cluster(s) à ce niveau"
        )
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_xlim(xmin - marge, xmax + marge)
        ax.set_ylim(ymin - marge, ymax + marge)
        ax.set_aspect("equal", adjustable="box")
        ax.grid(True, alpha=0.2)
        ax.legend(loc="upper right")

    nom = fichier_clusters.stem
    chemin_gif = dossier_sortie / f"{nom}.gif"
    chemin_png = dossier_sortie / f"{nom}.png"

    animation = FuncAnimation(
        fig, mise_a_jour, frames=niveaux,
        interval=1000 // max(fps, 1), repeat=True, blit=False,
    )
    animation.save(chemin_gif, writer=PillowWriter(fps=max(fps, 1)))
    mise_a_jour(niveau_max)
    fig.savefig(chemin_png, dpi=200, bbox_inches="tight")
    plt.close(fig)

    print(f"Animation créée : {chemin_gif}")
    print(f"Image créée : {chemin_png}")
    return [chemin_gif, chemin_png]
