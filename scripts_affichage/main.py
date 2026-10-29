from pathlib import Path

from io_utils import trouver_fichier
from visualiser_solution import visualiser_solution, EXPORTS_SOLUTION
from visualiser_clusters import visualiser_clusters
from visualiser_blocs import visualiser_blocs


RACINE = Path(__file__).resolve().parent.parent
DOSSIER_OUTPUTS = RACINE / "outputs"


EXPORTS = {
    "1": ("Solution : images et animation", "solution"),
    "2": ("Clusters BEM", "clusters_bem"),
    "3": ("Blocs BEM", "blocs_bem"),
    "4": ("Clusters de la solution", "clusters_solution"),
    "5": ("Blocs de reconstruction", "blocs_solution"),
}


def trouver_solution(dossier):
    return trouver_fichier(dossier, "solution", obligatoire=False)


def lister_calculs():
    if not DOSSIER_OUTPUTS.is_dir():
        raise FileNotFoundError(f"Dossier outputs introuvable : {DOSSIER_OUTPUTS}")

    calculs = []
    for dossier in sorted(DOSSIER_OUTPUTS.iterdir()):
        if dossier.is_dir() and trouver_solution(dossier) is not None:
            calculs.append((dossier, trouver_solution(dossier)))

    # Compatibilité avec l'ancien format où les fichiers étaient directement
    # dans outputs/ plutôt que dans un sous-dossier de calcul.
    if not calculs and trouver_fichier(DOSSIER_OUTPUTS, "solution", obligatoire=False):
        calculs.append((
            DOSSIER_OUTPUTS,
            trouver_fichier(DOSSIER_OUTPUTS, "solution", obligatoire=True),
        ))
    return calculs


def choisir_calcul(calculs):
    print("\n=== Calculs disponibles ===")
    for i, (dossier, solution) in enumerate(calculs, start=1):
        print(f"{i}. {dossier.name}  [solution : {solution.name}]")

    while True:
        saisie = input("\nNuméro du calcul (q pour quitter) : ").strip().lower()
        if saisie == "q":
            return None
        if saisie.isdigit() and 1 <= int(saisie) <= len(calculs):
            return calculs[int(saisie) - 1]
        print("Choix invalide.")


def choisir_exports():
    print("\n=== Que souhaitez-vous visualiser ? ===")
    for numero, (libelle, _) in EXPORTS.items():
        print(f"{numero}. {libelle}")
    print("A. Tout générer")

    while True:
        saisie = input("Choix séparés par des virgules (ex. 1,3,5) : ").strip().lower()
        if saisie == "a":
            return {cle for _, cle in EXPORTS.values()}
        choix = {x.strip() for x in saisie.split(",") if x.strip()}
        invalides = choix - set(EXPORTS)
        if not invalides and choix:
            return {EXPORTS[x][1] for x in choix}
        print("Choix invalide : indiquez un ou plusieurs numéros proposés.")


def choisir_exports_solution():
    print("\n=== Exports de la solution ===")
    options = list(EXPORTS_SOLUTION.keys()) + ["animation"]
    descriptions = {
        "uExact_Re.png": "Partie réelle de la solution exacte",
        "uExact_Im.png": "Partie imaginaire de la solution exacte",
        "uExact_Norm.png": "Norme de la solution exacte",
        "uCached_Re.png": "Partie réelle de la solution calculée",
        "uCached_Im.png": "Partie imaginaire de la solution calculée",
        "uCached_Norm.png": "Norme de la solution calculée",
        "uExactTot_Re.png": "Partie réelle du champ total exact",
        "uExactTot_Im.png": "Partie imaginaire du champ total exact",
        "uExactTot_Norm.png": "Norme du champ total exact",
        "uCachedTot_Re.png": "Partie réelle du champ total calculé",
        "uCachedTot_Im.png": "Partie imaginaire du champ total calculé",
        "uCachedTot_Norm.png": "Norme du champ total calculé",
        "Erreur.png": "Erreur absolue",
        "animation": "Animation temporelle GIF",
    }
    for i, cle in enumerate(options, start=1):
        print(f"{i}. {descriptions[cle]}")
    print("A. Tous les exports")

    while True:
        saisie = input("Choix séparés par des virgules (ex. 1,4,8) : ").strip().lower()
        if saisie == "a":
            return set(options)
        choix = {x.strip() for x in saisie.split(",") if x.strip()}
        if choix and all(x.isdigit() and 1 <= int(x) <= len(options) for x in choix):
            return {options[int(x) - 1] for x in choix}
        print("Choix invalide.")


def executer_visualisations(dossier, fichier_solution, exports):
    dossier_sortie = dossier / "visualisations"
    dossier_sortie.mkdir(parents=True, exist_ok=True)
    print(f"\nDossier du calcul : {dossier}")
    print(f"Les visualisations seront enregistrées dans : {dossier_sortie}")

    if "solution" in exports:
        maillage = trouver_fichier(dossier, "maillage", obligatoire=False)
        if maillage is None:
            print("Maillage absent : visualisation de la solution ignorée.")
        else:
            exports_solution = choisir_exports_solution()
            visualiser_solution(
                fichier_solution=fichier_solution,
                fichier_maillage=maillage,
                dossier_sortie=dossier_sortie,
                exports=exports_solution,
                frames=48,
                fps=12,
            )

    if "clusters_bem" in exports:
        maillage = trouver_fichier(dossier, "maillage")
        arbre = trouver_fichier(dossier, "clusters_bem")
        if maillage and arbre:
            visualiser_clusters(maillage, arbre, "segments", dossier_sortie)
        else:
            print("Clusters BEM ignorés : maillage ou CSV manquant.")

    if "blocs_bem" in exports:
        arbre = trouver_fichier(dossier, "blocs_bem")
        if arbre:
            visualiser_blocs(arbre, dossier_sortie)
        else:
            print("Blocs BEM ignorés : CSV manquant.")

    if "clusters_solution" in exports:
        maillage = trouver_fichier(dossier, "maillage_solution")
        arbre = trouver_fichier(dossier, "clusters_solution")
        if maillage and arbre:
            visualiser_clusters(maillage, arbre, "points", dossier_sortie)
        else:
            print("Clusters solution ignorés : maillage de solution ou CSV manquant.")

    if "blocs_solution" in exports:
        arbre = trouver_fichier(dossier, "blocs_solution")
        if arbre:
            visualiser_blocs(arbre, dossier_sortie)
        else:
            print("Blocs solution ignorés : CSV manquant.")


def main():
    try:
        calculs = lister_calculs()
    except Exception as exc:
        print(f"Erreur : {exc}")
        return

    if not calculs:
        print(
            f"Aucun dossier de calcul trouvé dans {DOSSIER_OUTPUTS}.\n"
            "Chaque dossier de calcul doit contenir solution.txt ou u_*.txt."
        )
        return

    selection = choisir_calcul(calculs)
    if selection is None:
        print("Arrêt demandé.")
        return

    dossier, solution = selection
    exports = choisir_exports()
    executer_visualisations(dossier, solution, exports)
    print("\nToutes les visualisations sélectionnées sont terminées.")


if __name__ == "__main__":
    main()
