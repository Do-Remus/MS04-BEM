#!/bin/bash

# ============================================================
# Lance le programme Helmholtz avec MPI et sauvegarde le log
# ============================================================

CONFIG="config.txt"
EXEC="bin/helmholtz"
LOG_DIR="log"

# ------------------------------------------------------------
# Vérifications
# ------------------------------------------------------------

if [ ! -f "$CONFIG" ]; then
    echo "Erreur : $CONFIG introuvable."
    exit 1
fi

if [ ! -f "$EXEC" ]; then
    echo "Erreur : $EXEC introuvable."
    exit 1
fi

# ------------------------------------------------------------
# Nombre de coeurs MPI
# ------------------------------------------------------------

read -p "Nombre de coeurs MPI : " NCORES

if ! [[ "$NCORES" =~ ^[1-9][0-9]*$ ]]; then
    echo "Erreur : le nombre de coeurs doit être un entier positif."
    exit 1
fi

# ------------------------------------------------------------
# Fonction pour récupérer un paramètre dans config.txt
#
# Accepte par exemple :
#   k=10
#   k = 10
# ------------------------------------------------------------

get_param() {
    local param="$1"

    grep -E "^[[:space:]]*$param[[:space:]]*=" "$CONFIG" \
        | head -n 1 \
        | sed -E 's/^[^=]*=[[:space:]]*//' \
        | tr -d '[:space:]'
}

# ------------------------------------------------------------
# Lecture des paramètres
# ------------------------------------------------------------

k=$(get_param "k")
rayon=$(get_param "rayon")
L=$(get_param "L")
pasMaillage=$(get_param "pasMaillage")
pasSolution=$(get_param "pasSolution")
idxTroncature=$(get_param "idxTroncature")
ordre=$(get_param "ordre")
maxIterGradConj=$(get_param "maxIterGradConj")
tolGradConj=$(get_param "tolGradConj")

# ------------------------------------------------------------
# Conversion de la tolérance en notation scientifique
# uniquement pour rendre le nom plus compact.
# ------------------------------------------------------------

tol_log=$(awk -v t="$tolGradConj" 'BEGIN {
    printf "%.0e", t
}')

# Retire le + dans e+00 éventuel
tol_log=${tol_log/e+/e}

# ------------------------------------------------------------
# Nom de base du run
# ------------------------------------------------------------

BASE="k${k}_R${rayon}_L${L}_hM${pasMaillage}_hS${pasSolution}_N${idxTroncature}_q${ordre}_tol${tol_log}_it${maxIterGradConj}_core${NCORES}"

# ------------------------------------------------------------
# Création du dossier log
# ------------------------------------------------------------

mkdir -p "$LOG_DIR"

# ------------------------------------------------------------
# Recherche du prochain numéro de run
# ------------------------------------------------------------

RUN=1

while true; do
    RUN_NAME=$(printf "%s_run%03d.log" "$BASE" "$RUN")
    LOG_FILE="$LOG_DIR/$RUN_NAME"

    if [ ! -e "$LOG_FILE" ]; then
        break
    fi

    ((RUN++))
done

# ------------------------------------------------------------
# Affichage des informations
# ------------------------------------------------------------

echo
echo "=============================================="
echo "          Lancement du calcul MPI"
echo "=============================================="
echo "Coeurs MPI       : $NCORES"
echo "Configuration    : $CONFIG"
echo "Exécutable       : $EXEC"
echo "Log              : $LOG_FILE"
echo "=============================================="
echo

# ------------------------------------------------------------
# Lancement
#
# tee permet :
#   - d'afficher la sortie dans le terminal
#   - de la copier simultanément dans le fichier log
# ------------------------------------------------------------

mpirun -n "$NCORES" "$EXEC" 2>&1 | tee "$LOG_FILE"

# ------------------------------------------------------------
# Code de retour de mpirun
# ------------------------------------------------------------

STATUS=${PIPESTATUS[0]}

echo
echo "=============================================="

if [ "$STATUS" -eq 0 ]; then
    echo "Calcul terminé avec succès."
else
    echo "Calcul terminé avec une erreur (code $STATUS)."
fi

echo "Log : $LOG_FILE"
echo "=============================================="

exit "$STATUS"