
#!/usr/bin/env bash

# ============================================================
# Lance Helmholtz avec MPI et sauvegarde des logs
#
# Usage :
#   ./run.sh
#   ./run.sh -debug
#
# Mode debug :
#   - mesure périodique du PSS et du RSS des processus Helmholtz
#   - sauvegarde des mesures dans un CSV
#   - résumé des pics observés
# ============================================================

set -o pipefail

CONFIG="config.txt"
EXEC="bin/helmholtz"
LOG_DIR="log"
PROCESS_NAME="$(basename "$EXEC")"

DEBUG=false
STATUS=1

# ------------------------------------------------------------
# Gestion des options
# ------------------------------------------------------------

if [[ "$#" -gt 1 ]]; then
    printf 'Usage : %s [-debug]\n' "$0"
    exit 1
fi

if [[ "$#" -eq 1 ]]; then
    if [[ "$1" == "-debug" ]]; then
        DEBUG=true
    else
        printf 'Option inconnue : %s\n' "$1"
        printf 'Usage : %s [-debug]\n' "$0"
        exit 1
    fi
fi

# ------------------------------------------------------------
# Vérifications
# ------------------------------------------------------------

if [[ ! -f "$CONFIG" ]]; then
    printf 'Erreur : %s introuvable.\n' "$CONFIG"
    exit 1
fi

if [[ ! -x "$EXEC" ]]; then
    printf 'Erreur : %s introuvable ou non exécutable.\n' "$EXEC"
    exit 1
fi

if [[ "$DEBUG" == true ]]; then
    for cmd in mpirun pgrep awk tee; do
        if ! command -v "$cmd" >/dev/null 2>&1; then
            printf 'Erreur : commande %s introuvable.\n' "$cmd"
            exit 1
        fi
    done

    if [[ ! -d /proc/self ]]; then
        printf 'Erreur : le mode debug nécessite Linux (/proc).\n'
        exit 1
    fi
fi

# ------------------------------------------------------------
# Nombre de rangs MPI
# ------------------------------------------------------------

read -r -p "Nombre de coeurs MPI : " NCORES

if [[ ! "$NCORES" =~ ^[1-9][0-9]*$ ]]; then
    printf 'Erreur : le nombre de coeurs doit être un entier positif.\n'
    exit 1
fi

# ------------------------------------------------------------
# Lecture des paramètres de config.txt
# ------------------------------------------------------------

get_param() {
    local param="$1"

    awk -F '=' -v key="$param" '
        {
            name = $1
            gsub(/^[[:space:]]+|[[:space:]]+$/, "", name)

            if (name == key) {
                value = substr($0, index($0, "=") + 1)
                gsub(/^[[:space:]]+|[[:space:]]+$/, "", value)
                print value
                exit
            }
        }
    ' "$CONFIG"
}

k=$(get_param "k")
rayon=$(get_param "rayon")
L=$(get_param "L")
pasMaillage=$(get_param "pasMaillage")
pasSolution=$(get_param "pasSolution")
idxTroncature=$(get_param "idxTroncature")
ordre=$(get_param "ordre")
maxIterGradConj=$(get_param "maxIterGradConj")
tolGradConj=$(get_param "tolGradConj")
maillage=$(get_param "maillage")

# ------------------------------------------------------------
# Conversion de la tolérance
# ------------------------------------------------------------

tol_log=$(awk -v t="$tolGradConj" 'BEGIN {
    if (t ~ /^[0-9.eE+-]+$/ && t + 0 > 0)
        printf "%.0e", t
    else
        print "inconnue"
}')

tol_log="${tol_log/e+/e}"

# ------------------------------------------------------------
# Nom du run
# ------------------------------------------------------------

BASE="${maillage}_k${k}_R${rayon}_L${L}_hM${pasMaillage}_hS${pasSolution}_N${idxTroncature}_q${ordre}_tol${tol_log}_it${maxIterGradConj}_core${NCORES}"

mkdir -p "$LOG_DIR" || {
    printf 'Erreur : impossible de créer %s.\n' "$LOG_DIR"
    exit 1
}

RUN=1

while true; do
    RUN_NAME=$(printf '%s_run%03d' "$BASE" "$RUN")

    LOG_FILE="$LOG_DIR/${RUN_NAME}.log"
    MEM_FILE="$LOG_DIR/${RUN_NAME}_memory.csv"
    MEM_SUMMARY="$LOG_DIR/${RUN_NAME}_memory_summary.txt"

    if [[ ! -e "$LOG_FILE" &&
          ! -e "$MEM_FILE" &&
          ! -e "$MEM_SUMMARY" ]]; then
        break
    fi

    ((RUN++))
done

# ------------------------------------------------------------
# Affichage des informations
# ------------------------------------------------------------

printf '\n==============================================\n'
printf '          Lancement du calcul MPI\n'
printf '==============================================\n'
printf 'Rangs MPI       : %s\n' "$NCORES"
printf 'Configuration   : %s\n' "$CONFIG"
printf 'Exécutable      : %s\n' "$EXEC"
printf 'Log             : %s\n' "$LOG_FILE"

if [[ "$DEBUG" == true ]]; then
    printf 'Mode debug      : activé\n'
    printf 'Mesures mémoire : %s\n' "$MEM_FILE"
    printf 'Résumé mémoire  : %s\n' "$MEM_SUMMARY"
else
    printf 'Mode debug      : désactivé\n'
fi

printf '==============================================\n\n'

# ------------------------------------------------------------
# Exécution normale
# ------------------------------------------------------------

if [[ "$DEBUG" == false ]]; then

    mpirun -n "$NCORES" "$EXEC" 2>&1 | tee "$LOG_FILE"
    STATUS=${PIPESTATUS[0]}

else

    # --------------------------------------------------------
    # Initialisation des mesures
    # --------------------------------------------------------

    printf 'elapsed_s,n_processes,sum_pss_kib,sum_rss_kib,max_rss_kib\n' \
        > "$MEM_FILE"

    START_TIME=$SECONDS

    PEAK_PSS_KIB=0
    PEAK_RSS_KIB=0
    PEAK_SINGLE_RSS_KIB=0
    PEAK_PROCESSES=0
    SAMPLES=0

    # --------------------------------------------------------
    # Lancement MPI en arrière-plan
    # --------------------------------------------------------

    mpirun -n "$NCORES" "$EXEC" \
        > >(tee "$LOG_FILE") 2>&1 &

    MPI_PID=$!

    # --------------------------------------------------------
    # Surveillance mémoire
    # --------------------------------------------------------

    while kill -0 "$MPI_PID" 2>/dev/null; do

        NPROC=0
        SUM_PSS=0
        SUM_RSS=0
        MAX_RSS=0

        for PID in $(pgrep -x "$PROCESS_NAME" 2>/dev/null); do

            FILE="/proc/$PID/smaps_rollup"

            if [[ -r "$FILE" ]]; then

                RSS=""
                PSS=""

                read -r RSS PSS < <(
                    awk '
                        /^Rss:/ { rss = $2 }
                        /^Pss:/ { pss = $2 }

                        END {
                            if (rss != "" && pss != "")
                                print rss, pss
                        }
                    ' "$FILE" 2>/dev/null
                )

                if [[ "$RSS" =~ ^[0-9]+$ &&
                      "$PSS" =~ ^[0-9]+$ ]]; then

                    ((NPROC++))
                    ((SUM_RSS += RSS))
                    ((SUM_PSS += PSS))

                    if (( RSS > MAX_RSS )); then
                        MAX_RSS=$RSS
                    fi
                fi
            fi
        done

        ELAPSED=$((SECONDS - START_TIME))

        printf '%s,%s,%s,%s,%s\n' \
            "$ELAPSED" "$NPROC" "$SUM_PSS" "$SUM_RSS" "$MAX_RSS" \
            >> "$MEM_FILE"

        ((SAMPLES++))

        if (( SUM_PSS > PEAK_PSS_KIB )); then
            PEAK_PSS_KIB=$SUM_PSS
        fi

        if (( SUM_RSS > PEAK_RSS_KIB )); then
            PEAK_RSS_KIB=$SUM_RSS
        fi

        if (( MAX_RSS > PEAK_SINGLE_RSS_KIB )); then
            PEAK_SINGLE_RSS_KIB=$MAX_RSS
        fi

        if (( NPROC > PEAK_PROCESSES )); then
            PEAK_PROCESSES=$NPROC
        fi

        sleep 0.5
    done

    # Récupération du code de retour de mpirun
    wait "$MPI_PID"
    STATUS=$?

    # --------------------------------------------------------
    # Résumé mémoire
    # --------------------------------------------------------

    {
        printf 'Résumé de la mesure mémoire\n'
        printf 'Exécutable : %s\n' "$EXEC"
        printf 'Rangs MPI demandés : %s\n' "$NCORES"
        printf 'Nombre maximal de processus détectés : %s\n' \
            "$PEAK_PROCESSES"
        printf 'Nombre d échantillons : %s\n\n' "$SAMPLES"

        awk \
            -v p="$PEAK_PSS_KIB" \
            -v r="$PEAK_RSS_KIB" \
            -v m="$PEAK_SINGLE_RSS_KIB" '
            BEGIN {
                printf "Pic PSS cumulé observé : %.2f MiB\n", p / 1024
                printf "Pic RSS cumulé observé : %.2f MiB\n", r / 1024
                printf "Pic RSS d un rang : %.2f MiB\n", m / 1024
            }
        '

        printf '\n'
        printf 'PSS : pages partagées réparties proportionnellement.\n'
        printf 'RSS cumulé : peut compter plusieurs fois les pages partagées.\n'
        printf 'Mesure échantillonnée toutes les 0,5 s.\n'
        printf 'Les processus auxiliaires MPI peuvent ne pas être inclus.\n'

    } > "$MEM_SUMMARY"

    cat "$MEM_SUMMARY"
fi

# ------------------------------------------------------------
# Résultat final
# ------------------------------------------------------------

printf '\n==============================================\n'

if [[ "$STATUS" -eq 0 ]]; then
    printf 'Calcul terminé avec succès.\n'
else
    printf 'Calcul terminé avec une erreur (code %s).\n' "$STATUS"
fi

printf 'Log : %s\n' "$LOG_FILE"

if [[ "$DEBUG" == true ]]; then
    printf 'Mesures mémoire : %s\n' "$MEM_FILE"
    printf 'Résumé mémoire : %s\n' "$MEM_SUMMARY"
fi

printf '==============================================\n'

exit "$STATUS"
