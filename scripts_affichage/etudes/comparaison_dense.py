
#!/usr/bin/env python3

import re
import argparse
from pathlib import Path
from collections import defaultdict

import numpy as np
import matplotlib.pyplot as plt


# ============================================================
# CONFIGURATION
# ============================================================

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = Path(__file__).resolve().parents[2]

DEFAULT_LOG_DIR = PROJECT_DIR / "log"
OUTPUT_DIR = SCRIPT_DIR / "figures"

DEFAULT_NB_TESTS = 20


# ============================================================
# EXTRACTION DES DONNEES
# ============================================================

def extraire_valeur(texte, motif, nom, fichier):

    resultat = re.search(motif, texte, re.MULTILINE)

    if resultat is None:
        raise ValueError(
            f"Impossible de trouver '{nom}' dans {fichier.name}"
        )

    return float(resultat.group(1))


def extraire_resultat_gmres(texte, nom, fichier):

    motif = (
        rf"^\s*{re.escape(nom)}\s+"
        r"(\d+)\s+"
        r"([0-9.eE+-]+)\s+"
        r"(OUI|NON)\s*$"
    )

    resultat = re.search(motif, texte, re.MULTILINE)

    if resultat is None:
        raise ValueError(
            f"Impossible de trouver les résultats de '{nom}' "
            f"dans {fichier.name}"
        )

    return {
        "iterations": int(resultat.group(1)),
        "residu": float(resultat.group(2)),
        "convergence": resultat.group(3) == "OUI"
    }


def lire_log(fichier):

    texte = fichier.read_text(
        encoding="utf-8",
        errors="replace"
    )

    # --------------------------------------------------------
    # Pas de maillage
    # --------------------------------------------------------

    match = re.search(
        r"_hM([0-9]+(?:\.[0-9]+)?)",
        fichier.name
    )

    if match is None:
        raise ValueError(
            f"Impossible d'extraire hM du nom : {fichier.name}"
        )

    hM = float(match.group(1))

    # --------------------------------------------------------
    # Tolérance ACA
    # --------------------------------------------------------

    tolerance_aca = extraire_valeur(
        texte,
        r"Tolerance ACA\s*:\s*([0-9.eE+-]+)",
        "Tolerance ACA",
        fichier
    )

    # --------------------------------------------------------
    # Mémoire
    # --------------------------------------------------------

    memoire_dense = extraire_valeur(
        texte,
        r"Matrice dense\s*:\s*([0-9.eE+-]+)\s*MiB",
        "Matrice dense",
        fichier
    )

    memoire_H = extraire_valeur(
        texte,
        r"H total estimee?\s*:\s*([0-9.eE+-]+)\s*MiB",
        "H total estimee",
        fichier
    )

    # --------------------------------------------------------
    # Temps moyens
    # --------------------------------------------------------

    temps_dense = extraire_valeur(
        texte,
        r"Temps moyen dense\s*:\s*([0-9.eE+-]+)\s*s",
        "Temps moyen dense",
        fichier
    )

    temps_H = extraire_valeur(
        texte,
        r"Temps moyen H\s*:\s*([0-9.eE+-]+)\s*s",
        "Temps moyen H",
        fichier
    )

    # --------------------------------------------------------
    # Erreurs relatives
    # --------------------------------------------------------

    erreur_moyenne = extraire_valeur(
        texte,
        r"Erreur moyenne\s*:\s*([0-9.eE+-]+)",
        "Erreur moyenne",
        fichier
    )

    ecart_type = extraire_valeur(
        texte,
        r"Ecart-type\s*:\s*([0-9.eE+-]+)",
        "Ecart-type",
        fichier
    )

    match_nb = re.search(
        r"Nombre de tests\s*:\s*(\d+)",
        texte
    )

    nb_tests = (
        int(match_nb.group(1))
        if match_nb
        else DEFAULT_NB_TESTS
    )

    # --------------------------------------------------------
    # Résultats GMRES
    # --------------------------------------------------------

    gmres_dense = extraire_resultat_gmres(
        texte,
        "GMRES dense",
        fichier
    )

    gmres_H = extraire_resultat_gmres(
        texte,
        "GMRES H sans Schwarz",
        fichier
    )

    gmres_H_schwarz = extraire_resultat_gmres(
        texte,
        "GMRES H avec Schwarz",
        fichier
    )

    return {
        "hM": hM,
        "tolerance_aca": tolerance_aca,

        "memoire_dense": memoire_dense,
        "memoire_H": memoire_H,

        "temps_dense": temps_dense,
        "temps_H": temps_H,

        "erreur_moyenne": erreur_moyenne,
        "ecart_type": ecart_type,
        "nb_tests": nb_tests,

        "iterations_dense": gmres_dense["iterations"],
        "iterations_H": gmres_H["iterations"],
        "iterations_H_schwarz": gmres_H_schwarz["iterations"],

        "residu_dense": gmres_dense["residu"],
        "residu_H": gmres_H["residu"],
        "residu_H_schwarz": gmres_H_schwarz["residu"],

        "convergence_dense": gmres_dense["convergence"],
        "convergence_H": gmres_H["convergence"],
        "convergence_H_schwarz": gmres_H_schwarz["convergence"],

        "fichier": fichier.name
    }


# ============================================================
# REGROUPEMENT DES RESULTATS
# ============================================================

def agreger_groupe(groupe):

    nb_fichiers = len(groupe)

    # --------------------------------------------------------
    # Mémoire
    # --------------------------------------------------------

    memoire_dense = np.mean([
        r["memoire_dense"] for r in groupe
    ])

    memoire_H = np.mean([
        r["memoire_H"] for r in groupe
    ])

    # --------------------------------------------------------
    # Temps : moyenne et min/max entre les runs
    # --------------------------------------------------------

    temps_dense = np.array([
        r["temps_dense"] for r in groupe
    ])

    temps_H = np.array([
        r["temps_H"] for r in groupe
    ])

    # --------------------------------------------------------
    # Erreurs : combinaison des statistiques
    # --------------------------------------------------------

    effectif_total = sum(
        r["nb_tests"] for r in groupe
    )

    erreur_moyenne = sum(
        r["nb_tests"] * r["erreur_moyenne"]
        for r in groupe
    ) / effectif_total

    variance_numerateur = 0.0

    for r in groupe:

        n = r["nb_tests"]
        mu = r["erreur_moyenne"]
        sigma = r["ecart_type"]

        variance_numerateur += (
            (n - 1) * sigma**2
            + n * (mu - erreur_moyenne)**2
        )

    erreur_ecart_type = (
        np.sqrt(
            variance_numerateur / (effectif_total - 1)
        )
        if effectif_total > 1
        else 0.0
    )

    # --------------------------------------------------------
    # Nombre d'itérations GMRES
    # --------------------------------------------------------

    iterations_dense = np.array([
        r["iterations_dense"] for r in groupe
    ])

    iterations_H = np.array([
        r["iterations_H"] for r in groupe
    ])

    iterations_H_schwarz = np.array([
        r["iterations_H_schwarz"] for r in groupe
    ])

    # --------------------------------------------------------
    # Résidus relatifs
    # --------------------------------------------------------

    residus_dense = np.array([
        r["residu_dense"] for r in groupe
    ])

    residus_H = np.array([
        r["residu_H"] for r in groupe
    ])

    residus_H_schwarz = np.array([
        r["residu_H_schwarz"] for r in groupe
    ])

    return {
        "nb_fichiers": nb_fichiers,

        "memoire_dense": memoire_dense,
        "memoire_H": memoire_H,

        "temps_dense": np.mean(temps_dense),
        "temps_dense_min": np.min(temps_dense),
        "temps_dense_max": np.max(temps_dense),

        "temps_H": np.mean(temps_H),
        "temps_H_min": np.min(temps_H),
        "temps_H_max": np.max(temps_H),

        "erreur_moyenne": erreur_moyenne,
        "erreur_ecart_type": erreur_ecart_type,
        "nb_tests": effectif_total,

        "iterations_dense": np.mean(iterations_dense),
        "iterations_dense_min": np.min(iterations_dense),
        "iterations_dense_max": np.max(iterations_dense),

        "iterations_H": np.mean(iterations_H),
        "iterations_H_min": np.min(iterations_H),
        "iterations_H_max": np.max(iterations_H),

        "iterations_H_schwarz": np.mean(iterations_H_schwarz),
        "iterations_H_schwarz_min": np.min(iterations_H_schwarz),
        "iterations_H_schwarz_max": np.max(iterations_H_schwarz),

        "residu_dense": np.mean(residus_dense),
        "residu_H": np.mean(residus_H),
        "residu_H_schwarz": np.mean(residus_H_schwarz),

        "convergence_dense": all(
            r["convergence_dense"] for r in groupe
        ),
        "convergence_H": all(
            r["convergence_H"] for r in groupe
        ),
        "convergence_H_schwarz": all(
            r["convergence_H_schwarz"] for r in groupe
        )
    }


def regrouper_resultats(resultats):

    groupes = defaultdict(list)

    for resultat in resultats:

        cle = (
            resultat["tolerance_aca"],
            resultat["hM"]
        )

        groupes[cle].append(resultat)

    donnees = {}

    for cle, groupe in sorted(groupes.items()):

        donnees[cle] = agreger_groupe(groupe)

    return donnees


# ============================================================
# FONCTIONS GENERIQUES POUR LES GRAPHIQUES
# ============================================================

def obtenir_abscisses(donnees, variable_x):

    return np.array([
        d[variable_x] for d in donnees
    ])


def configurer_axe_x(ax, variable_x):

    if variable_x == "hM":

        ax.set_xlabel(r"Pas de maillage $h_M$")

    else:

        ax.set_xlabel(r"Tolérance ACA")

    ax.set_xscale("log")


# ============================================================
# GRAPHIQUE 1 : MEMOIRE
# ============================================================

def plot_memoire(donnees, dossier, variable_x):

    x = obtenir_abscisses(donnees, variable_x)

    dense = [d["memoire_dense"] for d in donnees]
    hierarchique = [d["memoire_H"] for d in donnees]

    fig, ax = plt.subplots(figsize=(9, 6))

    ax.plot(
        x, dense,
        marker="o",
        linestyle="-",
        label="Matrice dense"
    )

    ax.plot(
        x, hierarchique,
        marker="s",
        linestyle="-",
        label="H-matrice"
    )

    configurer_axe_x(ax, variable_x)

    ax.set_ylabel("Mémoire (MiB)")
    ax.set_title("Comparaison mémoire : dense / H-matrice")

    ax.grid(True, which="both", linestyle="--", alpha=0.5)
    ax.legend()

    fig.tight_layout()

    fig.savefig(
        dossier / "01_memoire.png",
        dpi=300,
        bbox_inches="tight"
    )

    plt.close(fig)


# ============================================================
# GRAPHIQUE 2 : PERFORMANCE
# ============================================================

def plot_performance(donnees, dossier, variable_x):

    x = obtenir_abscisses(donnees, variable_x)

    temps_dense = np.array([
        d["temps_dense"] for d in donnees
    ])

    temps_H = np.array([
        d["temps_H"] for d in donnees
    ])

    erreur_dense = np.array([
        [
            d["temps_dense"] - d["temps_dense_min"]
            for d in donnees
        ],
        [
            d["temps_dense_max"] - d["temps_dense"]
            for d in donnees
        ]
    ])

    erreur_H = np.array([
        [
            d["temps_H"] - d["temps_H_min"]
            for d in donnees
        ],
        [
            d["temps_H_max"] - d["temps_H"]
            for d in donnees
        ]
    ])

    fig, ax = plt.subplots(figsize=(9, 6))

    ax.errorbar(
        x, temps_dense,
        yerr=erreur_dense,
        marker="o",
        linestyle="-",
        capsize=5,
        label="Matrice dense"
    )

    ax.errorbar(
        x, temps_H,
        yerr=erreur_H,
        marker="s",
        linestyle="-",
        capsize=5,
        label="H-matrice"
    )

    configurer_axe_x(ax, variable_x)

    ax.set_ylabel("Temps moyen par produit Ax (s)")
    ax.set_title("Performance des produits matrice-vecteur")

    ax.grid(True, which="both", linestyle="--", alpha=0.5)
    ax.legend()

    fig.tight_layout()

    fig.savefig(
        dossier / "02_performance.png",
        dpi=300,
        bbox_inches="tight"
    )

    plt.close(fig)


# ============================================================
# GRAPHIQUE 3 : ERREUR RELATIVE
# ============================================================

def plot_erreur(donnees, dossier, variable_x):

    x = obtenir_abscisses(donnees, variable_x)

    erreur = np.array([
        d["erreur_moyenne"] for d in donnees
    ])

    ecart_type = np.array([
        d["erreur_ecart_type"] for d in donnees
    ])

    fig, ax = plt.subplots(figsize=(9, 6))

    ax.errorbar(
        x,
        erreur,
        yerr=ecart_type,
        marker="o",
        linestyle="-",
        capsize=5,
        label="Erreur relative moyenne ± écart-type"
    )

    configurer_axe_x(ax, variable_x)

    ax.set_ylabel("Erreur relative")
    ax.set_title("Précision de l'approximation H-matrice")

    ax.set_yscale("log")

    ax.grid(True, which="both", linestyle="--", alpha=0.5)
    ax.legend()

    fig.tight_layout()

    fig.savefig(
        dossier / "03_erreur.png",
        dpi=300,
        bbox_inches="tight"
    )

    plt.close(fig)


# ============================================================
# GRAPHIQUE 4 : NOMBRE D'ITERATIONS GMRES
# ============================================================

def plot_iterations(donnees, dossier, variable_x):

    x = obtenir_abscisses(donnees, variable_x)

    # --------------------------------------------------------
    # Itérations moyennes
    # --------------------------------------------------------

    iterations_dense = np.array([
        d["iterations_dense"] for d in donnees
    ])

    iterations_H = np.array([
        d["iterations_H"] for d in donnees
    ])

    iterations_H_schwarz = np.array([
        d["iterations_H_schwarz"] for d in donnees
    ])

    # --------------------------------------------------------
    # Barres d'erreur min/max entre les runs
    # --------------------------------------------------------

    erreur_dense = np.array([
        [
            d["iterations_dense"] - d["iterations_dense_min"]
            for d in donnees
        ],
        [
            d["iterations_dense_max"] - d["iterations_dense"]
            for d in donnees
        ]
    ])

    erreur_H = np.array([
        [
            d["iterations_H"] - d["iterations_H_min"]
            for d in donnees
        ],
        [
            d["iterations_H_max"] - d["iterations_H"]
            for d in donnees
        ]
    ])

    erreur_H_schwarz = np.array([
        [
            d["iterations_H_schwarz"] - d["iterations_H_schwarz_min"]
            for d in donnees
        ],
        [
            d["iterations_H_schwarz_max"] - d["iterations_H_schwarz"]
            for d in donnees
        ]
    ])

    # --------------------------------------------------------
    # Tracé
    # --------------------------------------------------------

    fig, ax = plt.subplots(figsize=(9, 6))

    ax.errorbar(
        x,
        iterations_dense,
        yerr=erreur_dense,
        marker="o",
        linestyle="-",
        capsize=5,
        label="GMRES dense"
    )

    ax.errorbar(
        x,
        iterations_H,
        yerr=erreur_H,
        marker="s",
        linestyle="-",
        capsize=5,
        label="GMRES H sans Schwarz"
    )

    ax.errorbar(
        x,
        iterations_H_schwarz,
        yerr=erreur_H_schwarz,
        marker="^",
        linestyle="-",
        capsize=5,
        label="GMRES H avec Schwarz"
    )

    configurer_axe_x(ax, variable_x)

    ax.set_ylabel("Nombre d'itérations GMRES")
    ax.set_title("Comparaison du nombre d'itérations GMRES")

    ax.grid(True, which="both", linestyle="--", alpha=0.5)
    ax.legend()

    fig.tight_layout()

    fig.savefig(
        dossier / "04_iterations.png",
        dpi=300,
        bbox_inches="tight"
    )

    plt.close(fig)


# ============================================================
# GENERATION DES QUATRE GRAPHIQUES
# ============================================================

def generer_graphiques(donnees, dossier, variable_x):

    dossier.mkdir(
        parents=True,
        exist_ok=True
    )

    plot_memoire(donnees, dossier, variable_x)
    plot_performance(donnees, dossier, variable_x)
    plot_erreur(donnees, dossier, variable_x)
    plot_iterations(donnees, dossier, variable_x)


# ============================================================
# AFFICHAGE DES RESULTATS
# ============================================================

def afficher_tableau(donnees, titre, variable_x):

    print(f"\n================ {titre} ================\n")

    nom_x = "hM" if variable_x == "hM" else "Tol ACA"

    print(
        f"{nom_x:>12} | "
        f"{'Dense MiB':>10} | "
        f"{'H MiB':>10} | "
        f"{'Dense (s)':>12} | "
        f"{'H (s)':>12} | "
        f"{'Dense it.':>10} | "
        f"{'H it.':>10} | "
        f"{'H+Schwarz':>10}"
    )

    print("-" * 115)

    for d in donnees:

        valeur_x = d[variable_x]

        print(
            f"{valeur_x:12.5g} | "
            f"{d['memoire_dense']:10.4f} | "
            f"{d['memoire_H']:10.4f} | "
            f"{d['temps_dense']:12.6e} | "
            f"{d['temps_H']:12.6e} | "
            f"{d['iterations_dense']:10.1f} | "
            f"{d['iterations_H']:10.1f} | "
            f"{d['iterations_H_schwarz']:10.1f}"
        )

    print("\n--- Résidus relatifs moyens ---")

    for d in donnees:

        print(
            f"{nom_x} = {d[variable_x]:.5g} | "
            f"Dense = {d['residu_dense']:.6e} | "
            f"H = {d['residu_H']:.6e} | "
            f"H + Schwarz = {d['residu_H_schwarz']:.6e}"
        )

    print("\n--- Convergence ---")

    for d in donnees:

        print(
            f"{nom_x} = {d[variable_x]:.5g} | "
            f"Dense : {'OUI' if d['convergence_dense'] else 'NON'} | "
            f"H : {'OUI' if d['convergence_H'] else 'NON'} | "
            f"H + Schwarz : "
            f"{'OUI' if d['convergence_H_schwarz'] else 'NON'}"
        )


# ============================================================
# ETUDE 1 : VARIATION DU MAILLAGE
# ============================================================

def etude_maillage(donnees_par_couple):

    groupes = defaultdict(list)

    # Regroupement par tolérance ACA
    for (tolerance_aca, hM), donnees in donnees_par_couple.items():

        groupes[tolerance_aca].append({
            "hM": hM,
            **donnees
        })

    print("\n\n############################################")
    print("ETUDE 1 : VARIATION DU PAS DE MAILLAGE")
    print("############################################")

    for tolerance_aca, donnees in sorted(groupes.items()):

        donnees.sort(key=lambda d: d["hM"])

        nom_dossier = f"tol_ACA_{tolerance_aca:.0e}"

        dossier_sortie = (
            OUTPUT_DIR
            / "maillage"
            / nom_dossier
        )

        print("\n============================================")
        print(f"Tolérance ACA = {tolerance_aca:.1e}")
        print(f"Dossier : {dossier_sortie}")
        print("============================================")

        afficher_tableau(
            donnees,
            f"VARIATION DU MAILLAGE (tol ACA = {tolerance_aca:.1e})",
            "hM"
        )

        generer_graphiques(
            donnees,
            dossier_sortie,
            "hM"
        )

        print("\nGraphiques enregistrés dans :")
        print(dossier_sortie.resolve())


# ============================================================
# ETUDE 2 : VARIATION DE LA TOLERANCE ACA
# ============================================================

def etude_tolerance_aca(donnees_par_couple):

    groupes = defaultdict(list)

    # Regroupement par pas de maillage
    for (tolerance_aca, hM), donnees in donnees_par_couple.items():

        groupes[hM].append({
            "tolerance_aca": tolerance_aca,
            **donnees
        })

    print("\n\n############################################")
    print("ETUDE 2 : VARIATION DE LA TOLERANCE ACA")
    print("############################################")

    for hM, donnees in sorted(groupes.items()):

        donnees.sort(key=lambda d: d["tolerance_aca"])

        nom_dossier = f"h_{hM:.0e}"

        dossier_sortie = (
            OUTPUT_DIR
            / "toleranceACA"
            / nom_dossier
        )

        print("\n============================================")
        print(f"Pas de maillage hM = {hM:.1e}")
        print(f"Dossier : {dossier_sortie}")
        print("============================================")

        afficher_tableau(
            donnees,
            f"VARIATION DE LA TOLERANCE ACA (hM = {hM:.1e})",
            "tolerance_aca"
        )

        generer_graphiques(
            donnees,
            dossier_sortie,
            "tolerance_aca"
        )

        print("\nGraphiques enregistrés dans :")
        print(dossier_sortie.resolve())


# ============================================================
# MAIN
# ============================================================

def main():

    parser = argparse.ArgumentParser(
        description=(
            "Etude des performances, de la précision et "
            "des itérations des H-matrices."
        )
    )

    parser.add_argument(
        "--log-dir",
        type=Path,
        default=DEFAULT_LOG_DIR,
        help="Dossier contenant les fichiers .log"
    )

    args = parser.parse_args()

    log_dir = args.log_dir

    if not log_dir.exists():

        raise FileNotFoundError(
            f"Dossier de logs introuvable : {log_dir}"
        )

    fichiers = sorted(log_dir.glob("*.log"))

    if not fichiers:

        raise RuntimeError(
            f"Aucun fichier .log trouvé dans {log_dir}"
        )

    print(f"\nNombre de fichiers trouvés : {len(fichiers)}")

    resultats = []

    for fichier in fichiers:

        try:

            resultat = lire_log(fichier)
            resultats.append(resultat)

            print(
                f"Lu : {fichier.name} "
                f"(hM = {resultat['hM']}, "
                f"tol ACA = {resultat['tolerance_aca']:.1e})"
            )

        except (ValueError, OSError) as erreur:

            print(f"[ERREUR] {erreur}")

    if not resultats:

        raise RuntimeError("Aucun fichier valide.")

    # Agrégation commune aux deux études
    donnees_par_couple = regrouper_resultats(resultats)

    # --------------------------------------------------------
    # Etude 1 : variation du maillage
    # --------------------------------------------------------

    etude_maillage(donnees_par_couple)

    # --------------------------------------------------------
    # Etude 2 : variation de la tolérance ACA
    # --------------------------------------------------------

    etude_tolerance_aca(donnees_par_couple)

    print("\n\n============================================")
    print("Toutes les études sont terminées.")
    print(f"Figures enregistrées dans : {OUTPUT_DIR.resolve()}")
    print("============================================")


if __name__ == "__main__":
    main()
