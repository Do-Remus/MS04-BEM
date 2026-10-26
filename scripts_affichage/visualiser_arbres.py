
import os
import pandas as pd
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib import colormaps


# ============================================================
# Configuration
# ============================================================

FICHIER_MAILLAGE = "outputs/cercle_maillage.txt"
FICHIER_CLUSTERS = "outputs/arbre_clusters.csv"
FICHIER_BLOCS = "outputs/arbre_blocs.csv"

DOSSIER_SORTIE = "outputs"
FPS = 1

os.makedirs(DOSSIER_SORTIE, exist_ok=True)


# ============================================================
# Lecture du maillage
# ============================================================

def lire_milieux_segments(fichier):
    """
    Lit les segments x1 y1 x2 y2 et calcule leurs milieux.
    Chaque ligne du fichier correspond à un segment du maillage.
    """
    if not os.path.isfile(fichier):
        raise FileNotFoundError(
            f"Fichier de maillage introuvable : {fichier}"
        )

    segments = pd.read_csv(
        fichier,
        sep=r"\s+",
        header=None,
        names=["x1", "y1", "x2", "y2"],
    )

    if segments.empty:
        raise ValueError(f"Le maillage {fichier} est vide.")

    if segments.isna().any().any():
        raise ValueError(
            f"Le fichier {fichier} contient des lignes invalides."
        )

    return pd.DataFrame({
        "x": 0.5 * (segments["x1"] + segments["x2"]),
        "y": 0.5 * (segments["y1"] + segments["y2"]),
    })


def lire_csv(fichier, colonnes_requises):
    """Lit un CSV et vérifie les colonnes nécessaires."""
    if not os.path.isfile(fichier):
        raise FileNotFoundError(
            f"Fichier introuvable : {fichier}"
        )

    df = pd.read_csv(fichier)

    if df.empty:
        raise ValueError(f"Le fichier {fichier} est vide.")

    manquantes = set(colonnes_requises) - set(df.columns)

    if manquantes:
        raise ValueError(
            f"Colonnes manquantes dans {fichier} : "
            f"{sorted(manquantes)}\n"
            f"Colonnes présentes : {list(df.columns)}"
        )

    return df


# ============================================================
# 1. Visualisation spatiale de l'arbre de clusters
# ============================================================

def visualiser_clusters():
    print("Lecture du maillage et de l'arbre de clusters...",
          flush=True)

    points = lire_milieux_segments(FICHIER_MAILLAGE)

    arbre = lire_csv(
        FICHIER_CLUSTERS,
        {
            "niveau", "taille",
            "x_min", "x_max", "y_min", "y_max",
        },
    )

    # Vérification de la cohérence entre le maillage et la racine.
    idx_racine = arbre["niveau"].idxmin()
    taille_racine = int(arbre.loc[idx_racine, "taille"])

    if len(points) != taille_racine:
        raise ValueError(
            f"Le maillage contient {len(points)} segments, "
            f"mais la racine contient {taille_racine} indices.\n"
            "Vérifie que le maillage exporté est celui utilisé "
            "pour construire l'arbre."
        )

    x = points["x"].to_numpy()
    y = points["y"].to_numpy()

    niveau_max = int(arbre["niveau"].max())
    niveaux = sorted(arbre["niveau"].astype(int).unique())

    xmin = min(points["x"].min(), arbre["x_min"].min())
    xmax = max(points["x"].max(), arbre["x_max"].max())
    ymin = min(points["y"].min(), arbre["y_min"].min())
    ymax = max(points["y"].max(), arbre["y_max"].max())

    marge = 0.05 * max(xmax - xmin, ymax - ymin, 1.0)

    fig, ax = plt.subplots(figsize=(9, 9))

    def mise_a_jour(niveau_courant):
        ax.clear()

        ax.scatter(
            x, y,
            s=12,
            color="black",
            label="Milieux des segments",
            zorder=5,
        )

        visibles = arbre[arbre["niveau"] <= niveau_courant]

        for _, noeud in visibles.iterrows():
            niveau = int(noeud["niveau"])

            couleur = colormaps["viridis"](
                niveau / max(niveau_max, 1)
            )

            largeur = noeud["x_max"] - noeud["x_min"]
            hauteur = noeud["y_max"] - noeud["y_min"]

            ax.add_patch(Rectangle(
                (noeud["x_min"], noeud["y_min"]),
                largeur,
                hauteur,
                fill=False,
                edgecolor=couleur,
                linewidth=1.4,
            ))

        nb = len(arbre[arbre["niveau"] == niveau_courant])

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

    animation = FuncAnimation(
        fig,
        mise_a_jour,
        frames=niveaux,
        interval=1000 // FPS,
        repeat=True,
        blit=False,
    )

    chemin_gif = os.path.join(
        DOSSIER_SORTIE, "arbre_clusters_spatial.gif"
    )
    chemin_png = os.path.join(
        DOSSIER_SORTIE, "arbre_clusters_spatial.png"
    )

    animation.save(
        chemin_gif,
        writer=PillowWriter(fps=FPS),
    )

    mise_a_jour(niveau_max)
    fig.savefig(chemin_png, dpi=200, bbox_inches="tight")
    plt.close(fig)

    print(f"Animation créée : {chemin_gif}")
    print(f"Image créée : {chemin_png}")


# ============================================================
# 2. Visualisation de l'arbre des blocs
# ============================================================

def visualiser_blocs():
    print("Lecture de l'arbre des blocs...", flush=True)

    arbre = lire_csv(
        FICHIER_BLOCS,
        {
            "niveau", "type",
            "i_debut", "i_fin",
            "j_debut", "j_fin",
        },
    )

    # Normalisation des types.
    arbre["type"] = (
        arbre["type"].astype(str).str.strip().str.lower()
    )

    correspondances = {
        "interne": "interne",
        "bloc_interne": "interne",
        "bloc interne": "interne",
        "admissible": "admissible",
        "bloc_admissible": "admissible",
        "bloc admissible": "admissible",
        "dense": "dense",
        "bloc_dense": "dense",
        "bloc dense": "dense",
    }

    arbre["type"] = arbre["type"].map(
        lambda t: correspondances.get(t, t)
    )

    couleurs = {
        "interne": "gray",
        "admissible": "green",
        "dense": "red",
        "vide": "blue",
    }

    types_inconnus = set(arbre["type"]) - set(couleurs)

    if types_inconnus:
        raise ValueError(
            f"Types de blocs inconnus : {sorted(types_inconnus)}.\n"
            f"Types acceptés : {sorted(couleurs)}"
        )

    # Vérification des intervalles.
    colonnes_indices = [
        "i_debut", "i_fin", "j_debut", "j_fin"
    ]

    for colonne in colonnes_indices:
        arbre[colonne] = pd.to_numeric(
            arbre[colonne], errors="raise"
        )

    if (arbre["i_fin"] < arbre["i_debut"]).any():
        raise ValueError("Un bloc possède un intervalle de lignes vide.")

    if (arbre["j_fin"] < arbre["j_debut"]).any():
        raise ValueError("Un bloc possède un intervalle de colonnes vide.")

    if (arbre[colonnes_indices] < 0).any().any():
        raise ValueError("Les indices des blocs doivent être positifs ou nuls.")

    niveau_max = int(arbre["niveau"].max())
    niveaux = sorted(arbre["niveau"].astype(int).unique())

    N_lignes = int(arbre["i_fin"].max())
    N_colonnes = int(arbre["j_fin"].max())

    fig, ax = plt.subplots(figsize=(10, 9))

    def mise_a_jour(niveau_courant):
        ax.clear()

        visibles = arbre[arbre["niveau"] <= niveau_courant]

        for _, bloc in visibles.iterrows():
            i0 = bloc["i_debut"]
            i1 = bloc["i_fin"]
            j0 = bloc["j_debut"]
            j1 = bloc["j_fin"]

            type_bloc = bloc["type"]

            ax.add_patch(Rectangle(
                (j0, i0),
                j1 - j0,
                i1 - i0,
                facecolor=couleurs[type_bloc],
                edgecolor="black",
                linewidth=0.6,
                alpha=0.45,
            ))

        ax.set_xlim(0, N_colonnes)
        ax.set_ylim(N_lignes, 0)
        ax.set_aspect("auto")

        ax.set_xlabel("Indices des colonnes")
        ax.set_ylabel("Indices des lignes")
        ax.set_title(
            f"Arbre des blocs — niveau {niveau_courant}"
        )
        ax.grid(False)

        legendes = [
            Rectangle(
                (0, 0), 1, 1,
                facecolor=couleurs[type_bloc],
                edgecolor="black",
                alpha=0.45,
                label=nom,
            )
            for nom, type_bloc in [
                ("Interne", "interne"),
                ("Admissible", "admissible"),
                ("Dense", "dense"),
            ]
        ]

        ax.legend(handles=legendes, loc="upper right")

    animation = FuncAnimation(
        fig,
        mise_a_jour,
        frames=niveaux,
        interval=1000 // FPS,
        repeat=True,
        blit=False,
    )

    chemin_gif = os.path.join(
        DOSSIER_SORTIE, "arbre_blocs.gif"
    )
    chemin_png = os.path.join(
        DOSSIER_SORTIE, "arbre_blocs.png"
    )

    animation.save(
        chemin_gif,
        writer=PillowWriter(fps=FPS),
    )

    mise_a_jour(niveau_max)
    fig.savefig(chemin_png, dpi=200, bbox_inches="tight")
    plt.close(fig)

    print(f"Animation créée : {chemin_gif}")
    print(f"Image créée : {chemin_png}")


# ============================================================
# Programme principal
# ============================================================

if __name__ == "__main__":
    print("=== Visualisation spatiale des clusters ===", flush=True)
    visualiser_clusters()

    print("\n=== Visualisation de l'arbre des blocs ===", flush=True)
    visualiser_blocs()

    print("\nToutes les visualisations sont terminées.", flush=True)
