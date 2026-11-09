#include "../information.hpp"

void affichage_parametres()
{
    affichage_ligne("Nombre d'onde k", k);
}

void affichage_approximations()
{
    affichage_ligne("Indice de troncature N", idxTroncature);
    affichage_ligne("Ordre de quadrature", ordre);
    affichage_ligne("Rang maximal ACA", rang_max);
    affichage_ligne("Tolerance ACA", toleranceAca);
    affichage_ligne("Facteur du pas green", facteur_step_green);
    affichage_ligne("Pas du cache de green", green_cache_step);
    affichage_ligne("Taille du cache", green_cache_size);
    affichage_ligne("Rang FFM", ordreFMM);

    const Real memoire = green_cache_size * sizeof(Complex) / (1024.0 * 1024.0);
    affichage_ligne("Memoire cache Green", std::to_string(memoire) + " MiB");
}

void affichage_erreur_relative_matrice(Real erreur)
{
    affichage_ligne("Erreur relative matrice H", erreur);
}

void affichage_maillage(const Point &centre, unsigned int nbSegments, const std::string &nomFichier)
{
    affichage_ligne("Centre", centre);
    affichage_ligne("Rayon", rayon);
    affichage_ligne("Pas", pasMaillage);
    affichage_ligne("Nombre de segments", nbSegments);
    affichage_ligne("Maillage exporte", nomFichier);
}

void affichage_erreur_relative_p(Real erreur)
{
    affichage_ligne("Erreur relative L2 p", erreur);
}

void affichage_maillage_solution(unsigned int nbPointsSolution)
{
    affichage_ligne("Domaine", std::string("[") +
                                   std::to_string(-L / 2.0) + ", " +
                                   std::to_string(L / 2.0) + "] x [" +
                                   std::to_string(-L / 2.0) + ", " +
                                   std::to_string(L / 2.0) + "]");
    affichage_ligne("Taille du domaine", std::to_string(L) + " x " + std::to_string(L));
    affichage_ligne("Pas", pasSolution);
    affichage_ligne("Distance a la frontiere", delta);
    affichage_ligne("Rayon minimal autorise", rayon + delta);
    affichage_ligne("Nombre de points", nbPointsSolution);
}

void affichage_erreurs_solution(Real erreurL2, Real erreurLinf)
{
    affichage_ligne("Erreur L2", erreurL2);
    affichage_ligne("Erreur Linf", erreurLinf);
}

void affichage_statistiques_matrice_H(const ArbreBlocs &arbre)
{
    if (!arbre)
    {
        affichage_ligne("Matrice H", "Arbre de blocs vide");
        return;
    }

    struct Statistiques
    {
        std::size_t nb_blocs = 0;
        std::size_t nb_blocs_internes = 0;
        std::size_t nb_blocs_admissibles = 0;
        std::size_t nb_blocs_denses = 0;

        std::size_t coefficients_denses = 0;
        std::size_t coefficients_UV = 0;

        std::size_t rang_maximal = 0;
        std::size_t somme_rangs = 0;
    };

    Statistiques stats;

    // ============================================================
    // 1. Parcours de l'arbre de blocs
    // ============================================================

    std::function<void(const Bloc &)> parcourir =
        [&](const Bloc &bloc)
    {
        ++stats.nb_blocs;

        switch (bloc.type)
        {
        case TypeBloc::Interne:
        {
            ++stats.nb_blocs_internes;

            for (const auto &ligne : bloc.fils)
            {
                for (const auto &fils : ligne)
                {
                    if (fils)
                        parcourir(*fils);
                }
            }

            break;
        }

        case TypeBloc::Dense:
        {
            ++stats.nb_blocs_denses;

            stats.coefficients_denses +=
                static_cast<std::size_t>(bloc.dense.size());

            break;
        }

        case TypeBloc::Admissible:
        {
            ++stats.nb_blocs_admissibles;

            const std::size_t rang =
                static_cast<std::size_t>(bloc.U.cols());

            stats.rang_maximal =
                std::max(stats.rang_maximal, rang);

            stats.somme_rangs += rang;

            stats.coefficients_UV +=
                static_cast<std::size_t>(bloc.U.size()) +
                static_cast<std::size_t>(bloc.V.size());

            break;
        }
        }
    };

    parcourir(*arbre);

    // ============================================================
    // 2. Dimensions et nombre de coefficients
    // ============================================================

    const std::size_t M =
        arbre->cluster_lignes
            ? arbre->cluster_lignes->indices.size()
            : 0;

    const std::size_t N =
        arbre->cluster_colonnes
            ? arbre->cluster_colonnes->indices.size()
            : 0;

    const std::size_t coefficients_H =
        stats.coefficients_denses + stats.coefficients_UV;

    const long double coefficients_denses_complets =
        static_cast<long double>(M) *
        static_cast<long double>(N);

    const long double octets_par_coefficient =
        sizeof(Matrice::Scalar);

    // ============================================================
    // 3. Estimation de la mémoire
    // ============================================================

    // Mémoire des objets Bloc.
    // sizeof(Bloc) inclut les matrices Eigen (objets uniquement),
    // les pointeurs vers les clusters et les unique_ptr des fils.
    const long double memoire_structure =
        static_cast<long double>(stats.nb_blocs) *
        sizeof(Bloc);

    // Mémoire des coefficients réellement alloués par Eigen.
    const long double memoire_coefficients =
        static_cast<long double>(coefficients_H) *
        octets_par_coefficient;

    // Mémoire totale estimée.
    const long double memoire_H =
        memoire_structure + memoire_coefficients;

    // Mémoire d'une matrice dense complète.
    const long double memoire_dense =
        coefficients_denses_complets *
        octets_par_coefficient;

    constexpr long double octets_par_MiB =
        1024.0L * 1024.0L;

    const double memoire_structure_MiB =
        static_cast<double>(memoire_structure / octets_par_MiB);

    const double memoire_coefficients_MiB =
        static_cast<double>(memoire_coefficients / octets_par_MiB);

    const double memoire_H_MiB =
        static_cast<double>(memoire_H / octets_par_MiB);

    const double memoire_dense_MiB =
        static_cast<double>(memoire_dense / octets_par_MiB);

    // ============================================================
    // 4. Statistiques complémentaires
    // ============================================================

    const double pourcentage_stockage =
        coefficients_denses_complets > 0.0L
            ? 100.0 * static_cast<double>(coefficients_H) /
                  static_cast<double>(coefficients_denses_complets)
            : 0.0;

    const double pourcentage_memoire =
        memoire_dense > 0.0L
            ? 100.0 * static_cast<double>(memoire_H / memoire_dense)
            : 0.0;

    const double rang_moyen =
        stats.nb_blocs_admissibles > 0
            ? static_cast<double>(stats.somme_rangs) /
                  static_cast<double>(stats.nb_blocs_admissibles)
            : 0.0;

    // ============================================================
    // 5. Affichage
    // ============================================================

    affichage_ligne("Taille de la matrice",
                    std::to_string(M) + " x " + std::to_string(N));

    affichage_ligne("Nombre total de blocs", stats.nb_blocs);
    affichage_ligne("Blocs internes", stats.nb_blocs_internes);
    affichage_ligne("Blocs admissibles", stats.nb_blocs_admissibles);
    affichage_ligne("Blocs denses", stats.nb_blocs_denses);

    affichage_ligne("Coefs denses", stats.coefficients_denses);
    affichage_ligne("Coefs U/V", stats.coefficients_UV);
    affichage_ligne("Coefs stockes", coefficients_H);

    affichage_ligne("Rang maximal ACA", stats.rang_maximal);
    affichage_ligne("Rang moyen ACA", std::to_string(rang_moyen));

    affichage_ligne("Taux de stockage coefficients",
                    std::to_string(pourcentage_stockage) + " %");

    // Mémoire
    affichage_ligne("Memoire structure arbre",
                    std::to_string(memoire_structure_MiB) + " MiB");

    affichage_ligne("Memoire coefficients",
                    std::to_string(memoire_coefficients_MiB) + " MiB");

    affichage_ligne("Memoire totale estimee H",
                    std::to_string(memoire_H_MiB) + " MiB");

    affichage_ligne("Memoire matrice dense",
                    std::to_string(memoire_dense_MiB) + " MiB");

    affichage_ligne("Taux memoire total",
                    std::to_string(pourcentage_memoire) + " %");
}

void afficher_stats_FMM(const StructureFMM &fmm)
{
    Index nombre_feuilles = 0;
    Index nombre_internes = 0;
    Index N = 0;

    constexpr double OCTETS_PAR_MIB = 1024.0 * 1024.0;

    // ============================================================
    // 1. Mémoire des noeuds
    // ============================================================

    const std::size_t memoire_noeuds =
        fmm.noeuds.size() * sizeof(NoeudFMM);

    std::size_t memoire_basse_indices = 0;
    std::size_t memoire_haute_indices = 0;

    std::size_t memoire_coefficients = 0;

    for (const NoeudFMM &noeud : fmm.noeuds)
    {
        if (noeud.est_feuille())
        {
            ++nombre_feuilles;
            N += noeud.indices.size();
        }
        else
        {
            ++nombre_internes;
        }

        memoire_basse_indices +=
            noeud.indices.size() * sizeof(Index);

        memoire_haute_indices +=
            noeud.indices.capacity() * sizeof(Index);

        memoire_coefficients +=
            static_cast<std::size_t>(
                noeud.multipole.size() + noeud.local.size()) *
            sizeof(Complex);
    }

    // ============================================================
    // 2. Blocs denses
    // ============================================================

    Index nombre_coefficients_denses = 0;

    std::size_t memoire_basse_dense = 0;
    std::size_t memoire_haute_dense = 0;

    for (const BlocDenseFMM &bloc : fmm.blocs_denses)
    {
        nombre_coefficients_denses +=
            static_cast<Index>(bloc.coefficients.size());

        // Eigen alloue généralement exactement la taille demandée.
        memoire_basse_dense +=
            static_cast<std::size_t>(bloc.coefficients.size()) * sizeof(Complex);

        memoire_haute_dense +=
            static_cast<std::size_t>(bloc.coefficients.size()) * sizeof(Complex);
    }

    // ============================================================
    // 3. Opérateurs de translation compacts
    // ============================================================

    auto calculer_memoire_translations =
        [](const std::vector<TranslationFMM> &translations)
    {
        std::size_t memoire_basse = 0;
        std::size_t memoire_haute = 0;

        Index coefficients = 0;

        for (const TranslationFMM &translation : translations)
        {
            coefficients +=
                static_cast<Index>(translation.coefficients.size());

            memoire_basse +=
                static_cast<std::size_t>(
                    translation.coefficients.size()) *
                sizeof(Complex);

            memoire_haute +=
                static_cast<std::size_t>(
                    translation.coefficients.size()) *
                sizeof(Complex);
        }

        return std::make_tuple(
            memoire_basse,
            memoire_haute,
            coefficients);
    };

    const auto [memoire_basse_M2M, memoire_haute_M2M, coefficients_M2M] =
        calculer_memoire_translations(fmm.translations_M2M);

    const auto [memoire_basse_M2L, memoire_haute_M2L, coefficients_M2L] =
        calculer_memoire_translations(fmm.translations_M2L);

    const auto [memoire_basse_L2L, memoire_haute_L2L, coefficients_L2L] =
        calculer_memoire_translations(fmm.translations_L2L);

    const std::size_t memoire_basse_translations_coefficients =
        memoire_basse_M2M + memoire_basse_M2L + memoire_basse_L2L;

    const std::size_t memoire_haute_translations_coefficients =
        memoire_haute_M2M + memoire_haute_M2L + memoire_haute_L2L;

    const Index nombre_coefficients_translations =
        coefficients_M2M + coefficients_M2L + coefficients_L2L;

    const std::size_t memoire_basse_translations_objets =
        (fmm.translations_M2M.size() + fmm.translations_M2L.size() + fmm.translations_L2L.size()) * sizeof(TranslationFMM);

    const std::size_t memoire_haute_translations_objets =
        (fmm.translations_M2M.capacity() + fmm.translations_M2L.capacity() + fmm.translations_L2L.capacity()) * sizeof(TranslationFMM);

    const std::size_t memoire_basse_translations =
        memoire_basse_translations_objets + memoire_basse_translations_coefficients;

    const std::size_t memoire_haute_translations =
        memoire_haute_translations_objets + memoire_haute_translations_coefficients;

    // ============================================================
    // 4. Fonctions de base et quadratures
    // ============================================================

    const std::size_t memoire_bases =
        static_cast<std::size_t>(fmm.bases_points.size()) * sizeof(Complex);

    const Index nombre_coefficients_bases =
        static_cast<Index>(fmm.bases_points.size());

    const std::size_t memoire_offsets_basse =
        fmm.offsets_quadratures.size() * sizeof(Index);

    const std::size_t memoire_offsets_haute =
        fmm.offsets_quadratures.capacity() * sizeof(Index);

    // ============================================================
    // 5. Interactions
    // ============================================================

    const std::size_t memoire_basse_proches =
        fmm.interactions_proches.size() * sizeof(InteractionFMM);

    const std::size_t memoire_haute_proches =
        fmm.interactions_proches.capacity() * sizeof(InteractionFMM);

    const std::size_t memoire_basse_interactions_M2L =
        fmm.interactions_M2L.size() * sizeof(InteractionFMM);

    const std::size_t memoire_haute_interactions_M2L =
        fmm.interactions_M2L.capacity() * sizeof(InteractionFMM);

    // ============================================================
    // 6. Objets des blocs denses
    // ============================================================

    const std::size_t memoire_basse_blocs =
        fmm.blocs_denses.size() * sizeof(BlocDenseFMM);

    const std::size_t memoire_haute_blocs =
        fmm.blocs_denses.capacity() * sizeof(BlocDenseFMM);

    // ============================================================
    // 7. Mémoire totale estimée
    // ============================================================

    const std::size_t memoire_basse =
        memoire_noeuds + memoire_basse_indices + memoire_coefficients + memoire_basse_dense + memoire_basse_blocs + memoire_basse_translations + memoire_bases + memoire_offsets_basse + memoire_basse_proches + memoire_basse_interactions_M2L;

    const std::size_t memoire_haute =
        memoire_noeuds + memoire_haute_indices + memoire_coefficients + memoire_haute_dense + memoire_haute_blocs + memoire_haute_translations + memoire_bases + memoire_offsets_haute + memoire_haute_proches + memoire_haute_interactions_M2L;

    // ============================================================
    // 8. Comparaison avec une matrice dense
    // ============================================================

    const double memoire_dense =
        static_cast<double>(N) * N * sizeof(Complex);

    const double basse_Mo =
        static_cast<double>(memoire_basse) / OCTETS_PAR_MIB;

    const double haute_Mo =
        static_cast<double>(memoire_haute) / OCTETS_PAR_MIB;

    const double dense_Mo =
        memoire_dense / OCTETS_PAR_MIB;

    const double blocs_Mo =
        static_cast<double>(memoire_basse_dense) / OCTETS_PAR_MIB;

    const double translations_M2M_Mo =
        static_cast<double>(memoire_basse_M2M) / OCTETS_PAR_MIB;

    const double translations_M2L_Mo =
        static_cast<double>(memoire_basse_M2L) / OCTETS_PAR_MIB;

    const double translations_L2L_Mo =
        static_cast<double>(memoire_basse_L2L) / OCTETS_PAR_MIB;

    const double bases_Mo =
        static_cast<double>(memoire_bases) / OCTETS_PAR_MIB;

    const double translations_Mo =
        static_cast<double>(memoire_basse_translations) / OCTETS_PAR_MIB;

    const double taux_basse =
        memoire_dense > 0.0
            ? 100.0 * static_cast<double>(memoire_basse) / memoire_dense
            : 0.0;

    const double taux_haute =
        memoire_dense > 0.0
            ? 100.0 * static_cast<double>(memoire_haute) / memoire_dense
            : 0.0;

    const double gain_basse =
        memoire_basse > 0
            ? memoire_dense / static_cast<double>(memoire_basse)
            : 0.0;

    const double gain_haute =
        memoire_haute > 0
            ? memoire_dense / static_cast<double>(memoire_haute)
            : 0.0;

    // ============================================================
    // 9. Affichage des statistiques
    // ============================================================

    affichage_ligne("Nombre de noeuds", fmm.noeuds.size());
    affichage_ligne("Nombre de feuilles", nombre_feuilles);
    affichage_ligne("Nombre de noeuds internes", nombre_internes);
    affichage_ligne("Nombre de niveaux", fmm.nombre_niveaux);

    affichage_ligne("Nombre de degres de liberte", N);

    affichage_ligne("Interactions proches", fmm.interactions_proches.size());
    affichage_ligne("Interactions M2L", fmm.interactions_M2L.size());

    affichage_ligne("Nombre de blocs denses", fmm.blocs_denses.size());
    affichage_ligne("Coefficients denses stockes", nombre_coefficients_denses);

    affichage_ligne("Translations M2M", fmm.translations_M2M.size());
    affichage_ligne("Translations M2L", fmm.translations_M2L.size());
    affichage_ligne("Translations L2L", fmm.translations_L2L.size());

    affichage_ligne("Coefficients translations", nombre_coefficients_translations);

    affichage_ligne("Coefficients bases", nombre_coefficients_bases);
    affichage_ligne("Points quadrature indexes", fmm.offsets_quadratures.size());

    affichage_ligne("Nombre d'onde k", fmm.k);
    affichage_ligne("Ordre FMM", fmm.ordre);
    affichage_ligne("Parametre eta", fmm.eta);

    // ============================================================
    // 10. Mémoire détaillée
    // ============================================================

    affichage_ligne(
        "Memoire noeuds",
        std::to_string(memoire_noeuds / OCTETS_PAR_MIB) + " MiB");

    affichage_ligne(
        "Memoire indices",
        std::to_string(memoire_basse_indices / OCTETS_PAR_MIB) + " / " + std::to_string(memoire_haute_indices / OCTETS_PAR_MIB) + " MiB");

    affichage_ligne(
        "Memoire multipoles et locaux",
        std::to_string(memoire_coefficients / OCTETS_PAR_MIB) + " MiB");

    affichage_ligne(
        "Memoire blocs denses",
        std::to_string(blocs_Mo) + " MiB");

    affichage_ligne(
        "Memoire translations M2M",
        std::to_string(translations_M2M_Mo) + " MiB");

    affichage_ligne(
        "Memoire translations M2L",
        std::to_string(translations_M2L_Mo) + " MiB");

    affichage_ligne(
        "Memoire translations L2L",
        std::to_string(translations_L2L_Mo) + " MiB");

    affichage_ligne(
        "Memoire totale translations",
        std::to_string(translations_Mo) + " MiB");

    affichage_ligne(
        "Memoire bases quadrature",
        std::to_string(bases_Mo) + " MiB");

    affichage_ligne(
        "Memoire offsets quadrature",
        std::to_string(memoire_offsets_basse / OCTETS_PAR_MIB) + " / " + std::to_string(memoire_offsets_haute / OCTETS_PAR_MIB) + " MiB");

    affichage_ligne(
        "Memoire interactions proches",
        std::to_string(memoire_basse_proches / OCTETS_PAR_MIB) + " / " + std::to_string(memoire_haute_proches / OCTETS_PAR_MIB) + " MiB");

    affichage_ligne(
        "Memoire interactions M2L",
        std::to_string(memoire_basse_interactions_M2L / OCTETS_PAR_MIB) + " / " + std::to_string(memoire_haute_interactions_M2L / OCTETS_PAR_MIB) + " MiB");

    // ============================================================
    // 11. Mémoire totale et comparaison
    // ============================================================

    affichage_ligne("Memoire FMM basse", std::to_string(basse_Mo) + " MiB");
    affichage_ligne("Memoire FMM haute", std::to_string(haute_Mo) + " MiB");
    affichage_ligne("Memoire matrice dense", std::to_string(dense_Mo) + " MiB");

    affichage_ligne("Taux stockage bas", std::to_string(taux_basse) + " %");
    affichage_ligne("Taux stockage haut", std::to_string(taux_haute) + " %");

    affichage_ligne("Gain memoire bas", gain_basse);
    affichage_ligne("Gain memoire haut", gain_haute);
}
