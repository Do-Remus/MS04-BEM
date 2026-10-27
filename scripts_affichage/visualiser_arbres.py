
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

FICHIER_MAILLAGE = "outputs/maillage.txt"
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

    # --------------------------------------------------------
    # 1. Normalisation des types
    # --------------------------------------------------------

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
    }

    types_inconnus = set(arbre["type"]) - set(couleurs)

    if types_inconnus:
        raise ValueError(
            f"Types de blocs inconnus : {sorted(types_inconnus)}"
        )

    # --------------------------------------------------------
    # 2. Vérification des indices
    # --------------------------------------------------------

    colonnes_indices = [
        "i_debut", "i_fin", "j_debut", "j_fin"
    ]

    for colonne in colonnes_indices:
        arbre[colonne] = pd.to_numeric(
            arbre[colonne], errors="raise"
        )

        if not (arbre[colonne] % 1 == 0).all():
            raise ValueError(
                f"Les indices de {colonne} doivent être entiers."
            )

        arbre[colonne] = arbre[colonne].astype(int)

    if (arbre[colonnes_indices] < 0).any().any():
        raise ValueError("Indices négatifs détectés.")

    if (arbre["i_fin"] <= arbre["i_debut"]).any():
        raise ValueError("Intervalle de lignes vide détecté.")

    if (arbre["j_fin"] <= arbre["j_debut"]).any():
        raise ValueError("Intervalle de colonnes vide détecté.")

    arbre["niveau"] = pd.to_numeric(
        arbre["niveau"], errors="raise"
    ).astype(int)

    # --------------------------------------------------------
    # 3. Lecture de est_diagonal
    # --------------------------------------------------------

    def convertir_booleen(valeur):
        if isinstance(valeur, bool):
            return valeur

        if pd.isna(valeur):
            raise ValueError("Valeur est_diagonal vide.")

        valeur = str(valeur).strip().lower()

        if valeur in {"true", "1", "oui", "yes"}:
            return True

        if valeur in {"false", "0", "non", "no"}:
            return False

        raise ValueError(
            f"Valeur est_diagonal non reconnue : {valeur!r}"
        )

    if "est_diagonal" in arbre.columns:
        arbre["est_diagonal"] = (
            arbre["est_diagonal"].map(convertir_booleen)
        )
    else:
        arbre["est_diagonal"] = (
            (arbre["i_debut"] == arbre["j_debut"])
            & (arbre["i_fin"] == arbre["j_fin"])
        )
        print(
            "Colonne est_diagonal absente : déduction "
            "depuis les intervalles.",
            flush=True,
        )

    # Vérification de cohérence entre le booléen et les intervalles.
    diagonal_intervalles = (
        (arbre["i_debut"] == arbre["j_debut"])
        & (arbre["i_fin"] == arbre["j_fin"])
    )

    incoherents = arbre[
        arbre["est_diagonal"] != diagonal_intervalles
    ]

    if not incoherents.empty:
        print(
            f"Attention : {len(incoherents)} nœud(s) ont un "
            "est_diagonal incompatible avec leurs intervalles.",
            flush=True,
        )

    # --------------------------------------------------------
    # 4. Préparation de l'affichage
    # --------------------------------------------------------

    niveaux = sorted(arbre["niveau"].unique())
    niveau_max = max(niveaux)

    N_lignes = int(arbre["i_fin"].max())
    N_colonnes = int(arbre["j_fin"].max())

    fig, ax = plt.subplots(figsize=(10, 9))

    # Hypothèse : les nœuds de type "interne" sont subdivisés,
    # tandis que les blocs admissibles et denses sont terminaux.
    feuilles = arbre[
        arbre["type"].isin(["admissible", "dense"])
    ]

    internes = arbre[arbre["type"] == "interne"]

    def dessiner_contour(i0, i1, j0, j1, couleur,
                         largeur=0.7, style="-"):
        ax.add_patch(Rectangle(
            (j0, i0),
            j1 - j0,
            i1 - i0,
            fill=False,
            edgecolor=couleur,
            linewidth=largeur,
            linestyle=style,
        ))

    def dessiner_feuille(bloc, miroir=False):
        i0 = int(bloc["i_debut"])
        i1 = int(bloc["i_fin"])
        j0 = int(bloc["j_debut"])
        j1 = int(bloc["j_fin"])

        if miroir:
            # Transposition des indices du bloc.
            i0, i1, j0, j1 = j0, j1, i0, i1

        couleur = couleurs[bloc["type"]]

        ax.add_patch(Rectangle(
            (j0, i0),
            j1 - j0,
            i1 - i0,
            facecolor=couleur,
            edgecolor=couleur if miroir else "black",
            linewidth=0.6,
            alpha=0.13 if miroir else 0.75,
            hatch="//" if miroir else None,
        ))

    # --------------------------------------------------------
    # 5. Animation
    # --------------------------------------------------------

    def mise_a_jour(niveau_courant):
        ax.clear()

        visibles = arbre[arbre["niveau"] <= niveau_courant]

        feuilles_visibles = feuilles[
            feuilles["niveau"] <= niveau_courant
        ]

        internes_visibles = internes[
            internes["niveau"] <= niveau_courant
        ]

        # 1) Blocs terminaux stockés.
        for _, bloc in feuilles_visibles.iterrows():
            dessiner_feuille(bloc)

        # 2) Miroirs uniquement des blocs terminaux non diagonaux.
        for _, bloc in feuilles_visibles.iterrows():
            if bool(bloc["est_diagonal"]):
                continue

            dessiner_feuille(bloc, miroir=True)

        # 3) Contours des nœuds internes, dessinés en dernier.
        for _, bloc in internes_visibles.iterrows():
            dessiner_contour(
                int(bloc["i_debut"]),
                int(bloc["i_fin"]),
                int(bloc["j_debut"]),
                int(bloc["j_fin"]),
                couleur="black",
                largeur=0.55,
            )
        
        # 4) Contours bleus de tous les nœuds diagonaux visibles,
        # y compris les nœuds internes.
        diagonaux_bleus = arbre[
            arbre["est_diagonal"]
            & (
                (arbre["niveau"] == niveau_courant)
                | (
                    (arbre["niveau"] < niveau_courant)
                    & arbre["type"].isin(["admissible", "dense"])
                )
            )
        ]

        for _, bloc in diagonaux_bleus.iterrows():
            dessiner_contour(
                int(bloc["i_debut"]),
                int(bloc["i_fin"]),
                int(bloc["j_debut"]),
                int(bloc["j_fin"]),
                couleur="blue",
                largeur=1.6,
            )

        nb_feuilles = len(feuilles_visibles)
        nb_admissibles = int(
            (feuilles_visibles["type"] == "admissible").sum()
        )
        nb_denses = int(
            (feuilles_visibles["type"] == "dense").sum()
        )
        nb_diagonaux = int(
            visibles["est_diagonal"].sum()
        )

        ax.set_xlim(0, N_colonnes)
        ax.set_ylim(N_lignes, 0)
        ax.set_aspect("equal", adjustable="box")

        ax.set_xlabel("Indices des colonnes")
        ax.set_ylabel("Indices des lignes")

        ax.set_title(
            f"Arbre des blocs symétrique — niveau {niveau_courant}\n"
            f"{nb_feuilles} bloc(s) terminal(aux) visible(s) "
            f"({nb_admissibles} admissibles, {nb_denses} denses)\n"
            f"{nb_diagonaux} nœud(s) diagonal(aux) visible(s)"
        )

        legendes = [
            Rectangle(
                (0, 0), 1, 1,
                facecolor=couleurs["admissible"],
                edgecolor="black",
                alpha=0.75,
                label="Bloc admissible stocké",
            ),
            Rectangle(
                (0, 0), 1, 1,
                facecolor=couleurs["dense"],
                edgecolor="black",
                alpha=0.75,
                label="Bloc dense stocké",
            ),
            Rectangle(
                (0, 0), 1, 1,
                fill=False,
                edgecolor="black",
                label="Nœud interne",
            ),
            Rectangle(
                (0, 0), 1, 1,
                fill=False,
                edgecolor="blue",
                linewidth=1.6,
                label="Bloc diagonal",
            ),
            Rectangle(
                (0, 0), 1, 1,
                facecolor="gray",
                edgecolor="gray",
                hatch="//",
                alpha=0.13,
                label="Miroir déduit par symétrie",
            ),
        ]

        ax.legend(handles=legendes, loc="upper right")
        ax.grid(False)

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
