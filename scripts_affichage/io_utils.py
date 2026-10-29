from pathlib import Path
import pandas as pd


# Alias centralisés : on peut adapter la nomenclature ici sans toucher
# aux fonctions de visualisation.
ALIAS_FICHIERS = {
    "solution": ["solution.txt", "u_*.txt"],
    "maillage": ["maillage.txt"],
    "maillage_solution": ["maillage_solution.txt"],
    "clusters_bem": ["arbre_clusters_BEM.csv", "arbre_clusters_bem.csv"],
    "blocs_bem": ["arbre_blocs_BEM.csv", "arbre_blocs_bem.csv"],
    "clusters_solution": ["arbre_clusters_solution.csv"],
    "blocs_solution": ["arbre_blocs_solution.csv"],
}


def trouver_fichier(dossier, cle, obligatoire=False):
    """Trouve le premier fichier correspondant aux alias de la clé."""
    dossier = Path(dossier)
    motifs = ALIAS_FICHIERS.get(cle)
    if motifs is None:
        raise KeyError(f"Clé de fichier inconnue : {cle}")

    for motif in motifs:
        candidats = sorted(dossier.glob(motif))
        fichiers = [p for p in candidats if p.is_file()]
        if fichiers:
            return fichiers[0]

    if obligatoire:
        attendus = ", ".join(motifs)
        raise FileNotFoundError(
            f"Aucun fichier correspondant à [{attendus}] dans {dossier}"
        )
    return None


def lire_csv(fichier, colonnes_requises):
    """Lit un CSV et vérifie les colonnes nécessaires."""
    fichier = Path(fichier)
    if not fichier.is_file():
        raise FileNotFoundError(f"Fichier introuvable : {fichier}")

    df = pd.read_csv(fichier)
    if df.empty:
        raise ValueError(f"Le fichier {fichier} est vide.")

    manquantes = set(colonnes_requises) - set(df.columns)
    if manquantes:
        raise ValueError(
            f"Colonnes manquantes dans {fichier} : {sorted(manquantes)}. "
            f"Colonnes présentes : {list(df.columns)}"
        )
    return df


def lire_milieux_segments(fichier):
    """Lit x1 y1 x2 y2 et calcule les milieux des segments."""
    fichier = Path(fichier)
    if not fichier.is_file():
        raise FileNotFoundError(f"Fichier de maillage introuvable : {fichier}")

    segments = pd.read_csv(
        fichier, sep=r"\s+", header=None,
        names=["x1", "y1", "x2", "y2"], comment="#"
    )
    if segments.empty or segments.isna().any().any():
        raise ValueError(f"Le fichier de maillage {fichier} est vide ou invalide.")

    return pd.DataFrame({
        "x": 0.5 * (segments["x1"] + segments["x2"]),
        "y": 0.5 * (segments["y1"] + segments["y2"]),
    })


def lire_points(fichier):
    """Lit un fichier x y, avec ou sans compteur de points en première ligne."""
    fichier = Path(fichier)
    if not fichier.is_file():
        raise FileNotFoundError(f"Fichier de points introuvable : {fichier}")

    with fichier.open("r", encoding="utf-8") as f:
        premiere_ligne = f.readline().strip()

    try:
        nombre_points = int(premiere_ligne)
        a_un_compteur = True
    except ValueError:
        a_un_compteur = False

    if a_un_compteur:
        points = pd.read_csv(
            fichier, sep=r"\s+", header=None, skiprows=1,
            names=["x", "y"], comment="#"
        )
        if len(points) != nombre_points:
            raise ValueError(
                f"{fichier} annonce {nombre_points} points, "
                f"mais en contient {len(points)}."
            )
    else:
        points = pd.read_csv(
            fichier, sep=r"\s+", header=None,
            names=["x", "y"], comment="#"
        )

    if points.empty or points.isna().any().any():
        raise ValueError(f"Le fichier de points {fichier} est vide ou invalide.")
    return points
