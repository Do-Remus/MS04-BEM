from pathlib import Path

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.patches import Rectangle

from io_utils import lire_csv


def _convertir_booleen(valeur):
    if isinstance(valeur, bool):
        return valeur
    if pd.isna(valeur):
        raise ValueError("Valeur booléenne vide dans le CSV.")
    valeur = str(valeur).strip().lower()
    if valeur in {"true", "1", "oui", "yes"}:
        return True
    if valeur in {"false", "0", "non", "no"}:
        return False
    raise ValueError(f"Valeur booléenne non reconnue : {valeur!r}")


def visualiser_blocs(fichier_blocs, dossier_sortie=None, fps=1):
    fichier_blocs = Path(fichier_blocs)
    dossier_sortie = Path(dossier_sortie or fichier_blocs.parent / "visualisations")
    dossier_sortie.mkdir(parents=True, exist_ok=True)

    arbre = lire_csv(
        fichier_blocs,
        {"niveau", "type", "i_debut", "i_fin", "j_debut", "j_fin", "est_symetrique"},
    ).copy()

    arbre["type"] = arbre["type"].astype(str).str.strip().str.lower()
    correspondances = {
        "interne": "interne", "bloc_interne": "interne", "bloc interne": "interne",
        "admissible": "admissible", "bloc_admissible": "admissible",
        "bloc admissible": "admissible",
        "dense": "dense", "bloc_dense": "dense", "bloc dense": "dense",
    }
    arbre["type"] = arbre["type"].map(lambda t: correspondances.get(t, t))
    couleurs = {"interne": "gray", "admissible": "green", "dense": "red"}

    inconnus = set(arbre["type"]) - set(couleurs)
    if inconnus:
        raise ValueError(f"Types de blocs inconnus : {sorted(inconnus)}")

    colonnes_indices = ["i_debut", "i_fin", "j_debut", "j_fin"]
    for col in colonnes_indices:
        arbre[col] = pd.to_numeric(arbre[col], errors="raise")
        if not (arbre[col] % 1 == 0).all():
            raise ValueError(f"Les indices de {col} doivent être entiers.")
        arbre[col] = arbre[col].astype(int)

    if (arbre[colonnes_indices] < 0).any().any():
        raise ValueError("Indices négatifs détectés.")
    if (arbre["i_fin"] <= arbre["i_debut"]).any():
        raise ValueError("Intervalle de lignes vide détecté.")
    if (arbre["j_fin"] <= arbre["j_debut"]).any():
        raise ValueError("Intervalle de colonnes vide détecté.")

    arbre["niveau"] = pd.to_numeric(arbre["niveau"], errors="raise").astype(int)
    arbre["est_symetrique"] = arbre["est_symetrique"].map(_convertir_booleen)
    valeurs_symetrie = arbre["est_symetrique"].unique()
    if len(valeurs_symetrie) != 1:
        raise ValueError("est_symetrique doit avoir la même valeur pour tous les blocs.")
    est_symetrique = bool(valeurs_symetrie[0])

    if "est_diagonal" in arbre.columns:
        arbre["est_diagonal"] = arbre["est_diagonal"].map(_convertir_booleen)
    else:
        arbre["est_diagonal"] = (
            (arbre["i_debut"] == arbre["j_debut"])
            & (arbre["i_fin"] == arbre["j_fin"])
        )

    niveaux = sorted(arbre["niveau"].unique())
    niveau_max = max(niveaux)
    n_lignes = int(arbre["i_fin"].max())
    n_colonnes = int(arbre["j_fin"].max())
    feuilles = arbre[arbre["type"].isin(["admissible", "dense"])]
    internes = arbre[arbre["type"] == "interne"]

    fig, ax = plt.subplots(figsize=(10, 9))

    def dessiner_contour(i0, i1, j0, j1, couleur, largeur=0.7, style="-"):
        ax.add_patch(Rectangle(
            (j0, i0), j1 - j0, i1 - i0,
            fill=False, edgecolor=couleur,
            linewidth=largeur, linestyle=style,
        ))

    def dessiner_feuille(bloc, miroir=False):
        i0, i1 = int(bloc["i_debut"]), int(bloc["i_fin"])
        j0, j1 = int(bloc["j_debut"]), int(bloc["j_fin"])
        if miroir:
            i0, i1, j0, j1 = j0, j1, i0, i1
        couleur = couleurs[bloc["type"]]
        ax.add_patch(Rectangle(
            (j0, i0), j1 - j0, i1 - i0,
            facecolor=couleur,
            edgecolor=couleur if miroir else "black",
            linewidth=0.6,
            alpha=0.13 if miroir else 0.75,
            hatch="//" if miroir else None,
        ))

    def mise_a_jour(niveau_courant):
        ax.clear()
        feuilles_visibles = feuilles[feuilles["niveau"] <= niveau_courant]
        internes_visibles = internes[internes["niveau"] <= niveau_courant]

        for _, bloc in feuilles_visibles.iterrows():
            dessiner_feuille(bloc)

        if est_symetrique:
            for _, bloc in feuilles_visibles.iterrows():
                if not bool(bloc["est_diagonal"]):
                    dessiner_feuille(bloc, miroir=True)

        for _, bloc in internes_visibles.iterrows():
            dessiner_contour(
                int(bloc["i_debut"]), int(bloc["i_fin"]),
                int(bloc["j_debut"]), int(bloc["j_fin"]),
                couleur="black", largeur=0.65,
            )

        diagonaux = arbre[
            arbre["est_diagonal"]
            & (
                (arbre["niveau"] == niveau_courant)
                | (
                    (arbre["niveau"] < niveau_courant)
                    & arbre["type"].isin(["admissible", "dense"])
                )
            )
        ]
        if est_symetrique:
            for _, bloc in diagonaux.iterrows():
                dessiner_contour(
                    int(bloc["i_debut"]), int(bloc["i_fin"]),
                    int(bloc["j_debut"]), int(bloc["j_fin"]),
                    couleur="blue", largeur=1.6,
                )

        nb_admissibles = int((feuilles_visibles["type"] == "admissible").sum())
        nb_denses = int((feuilles_visibles["type"] == "dense").sum())
        nature = "symétrique" if est_symetrique else "non symétrique / rectangulaire"

        ax.set_xlim(0, n_colonnes)
        ax.set_ylim(n_lignes, 0)
        ax.set_aspect("auto")
        ax.set_xlabel("Indices des colonnes")
        ax.set_ylabel("Indices des lignes")
        ax.set_title(
            f"Arbre des blocs {nature} — niveau {niveau_courant}\n"
            f"{len(feuilles_visibles)} bloc(s) terminal(aux) visible(s) "
            f"({nb_admissibles} admissibles, {nb_denses} denses)"
        )

        legendes = [
            Rectangle((0, 0), 1, 1, facecolor="green", edgecolor="black",
                      alpha=0.75, label="Bloc admissible stocké"),
            Rectangle((0, 0), 1, 1, facecolor="red", edgecolor="black",
                      alpha=0.75, label="Bloc dense stocké"),
            Rectangle((0, 0), 1, 1, fill=False, edgecolor="black",
                      label="Nœud interne"),
        ]
        if est_symetrique:
            legendes.extend([
                Rectangle((0, 0), 1, 1, fill=False, edgecolor="blue",
                          linewidth=1.6, label="Bloc diagonal"),
                Rectangle((0, 0), 1, 1, facecolor="gray", edgecolor="gray",
                          hatch="//", alpha=0.13,
                          label="Miroir déduit par symétrie"),
            ])
        ax.legend(handles=legendes, loc="upper right")
        ax.grid(False)

    nom = fichier_blocs.stem
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
