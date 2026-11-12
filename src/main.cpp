#include "wrapper/commun.hpp"

int main(int argc, char **argv)
{
    // -- Initialisation MPI --

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    // -- Initialisations globales --

    // - Random -
    srand(time(NULL));

    // - Parametères -
    const string config = "config.txt";
    get_config(config);

    std::string end_output_folder = "_k" + formater_nombre(k) +
                                    "_R" + formater_nombre(rayon) +
                                    "_L" + formater_nombre(L) +
                                    "_hM" + formater_nombre(pasMaillage) +
                                    "_hS" + formater_nombre(pasSolution) +
                                    "_d" + formater_nombre(delta) +
                                    "_N" + std::to_string(idxTroncature) +
                                    "_q" + std::to_string(ordre);
    std::string output_folder;

#define FINAL_CODE // FINAL_CODE or TESTS_CODE

#ifdef FINAL_CODE
    /* ------------------------- */
    /* -- Area for Final code -- */
    /* ------------------------- */

    /* === Début du programme === */

    Timers timers;
    timers.start("programme");

    if (mpi_rank == 0)
    {
        std::cout << "\n=== Début du programme ===" << std::endl;
        affichage_parametres();
    }

    // /* === Initialisation cache de green === */

    run("Initialisation de Green", timers, "green_cache", initialiser_green_cache, affichage_approximations);

    /* === Création Maillage === */

    Maillage maillage = run("Préparation du maillage de l'obstacle", timers, "maillage", [&]() -> Maillage
                            {
                            if (nom_maillage == "narval")
                            {
                                return maillage_narval(
                                    1.0,         // echelle
                                    pasMaillage, // pas
                                    true         // centrer
                                );
                            }
                            else if (nom_maillage == "mur_acoustique")
                            {
                                return maillage_mur_acoustique(
                                    0.5,         // largeur
                                    4.0,         // hauteur
                                    8,           // nbPointes
                                    0.6,         // profondeur
                                    pasMaillage, // pas
                                    true         // centrer
                                );
                            }
                            else if (nom_maillage == "balle_golf")
                            {
                                return maillage_balle_de_golf(
                                    rayon,
                                    12,         // nbAlveoles
                                    0.12,       // profondeur
                                    0.8,        // remplissage
                                    pasMaillage // pas
                                );
                            }
                            else if (nom_maillage == "rectangle")
                            {
                                return maillage_rectangle(
                                    2.0,        // largeur
                                    1.0,        // hauteur
                                    pasMaillage // pas
                                );
                            }
                            else if (nom_maillage == "ellipse")
                            {
                                return maillage_ellipse(
                                    1.0,        // demi-largeur
                                    0.5,        // demi-hauteur
                                    pasMaillage // pas
                                );   
                            }
                            else if (nom_maillage == "rectangle")
                            {
                                return maillage_rectangle(
                                    2.0,        // largeur
                                    1.0,        // hauteur
                                    pasMaillage // pas
                                );
                            }
                            else if (nom_maillage == "ellipse")
                            {
                                return maillage_ellipse(
                                    1.0,        // demi-largeur
                                    0.5,        // demi-hauteur
                                    pasMaillage // pas
                                );   
                            }
                            return maillage_cercle(); }, [&](const Maillage &maillage)
                            {
                            const std::filesystem::path dossier_sortie = std::filesystem::path("outputs") / (maillage.nom + end_output_folder);
                            std::filesystem::create_directories(dossier_sortie);
                            output_folder = dossier_sortie.string();

                            maillage.export_maillage((dossier_sortie / "maillage.txt").string());
                            affichage_maillage(Point(0, 0), maillage.size(), (dossier_sortie / "maillage.txt").string()); });

    /* === Calcul de la quadrature de Legendre === */

    const std::vector<QuadratureSegment> quadrature_maillage = run("Calcul de la quadrature", timers, "quadrature",
                                                                   [&]()
                                                                   {
                                                                       return get_quadrature_maillage(maillage, get_legendre_data(ordre));
                                                                   });

    /* === Réorganisation de la quadrature pour le SIMD === */

    // const QuadratureSIMD quadrature_simd = run("Préparation quadrature SIMD", timers, "quadrature_simd",
    //                                            [&]()
    //                                            {
    //                                                return quadrature_SIMD(maillage, quadrature_maillage);
    //                                            });

    /* === Construction de l'arbre de clusters du maillage === */

    Arbre arbre_clusters_maillage = run("Construction de l'arbre de clusters", timers, "arbre_clusters", [&]()
                                        { return construire_arbre_maillage(maillage, N_leaves); }, [&](const Arbre &arbre)
                                        { ecrire_arbre_clusters(arbre, output_folder + "/arbre_clusters_BEM.csv"); });

    /* === Construction de l'arbre de blocs H-matrice === */

    ArbreBlocs arbre_blocs_BEM = run("Construction de l'arbre de blocs", timers, "arbre_blocs_H", [&]()
                                     {
                                     ArbreBlocs arbre = construire_arbre_blocs(arbre_clusters_maillage, parametre_admissibilite, facteur_admissibilite, true);
                                     remplir_blocs_bem(*arbre, maillage, quadrature_maillage, rang_max, toleranceAca);
                                     return arbre; }, [&](const ArbreBlocs &arbre)
                                     { 
                                    IntervallesClusters intervalles = construire_intervalles_clusters(arbre_clusters_maillage);
                                    affichage_statistiques_matrice_H(arbre); 
                                    ecrire_arbre_blocs(arbre, output_folder + "/arbre_blocs_BEM.csv", intervalles, intervalles, true); });

    /* === Calcul b de la FV de l'équation intégrale === */

    Vecteur vect_b = run("Calcul du second membre de la FV", timers, "b",
                         [&]()
                         {
                             return calculer_second_membre(quadrature_maillage);
                         });

    /* === Construction du préconditionneur par blocs diagonaux === */

    AccesMatrice acces_A = [&](Eigen::Index i, Eigen::Index j) -> Complex
    {
        return matrice_bem(i, j, maillage, quadrature_maillage);
    };

    PreconditionneurBlocs preconditionneur_blocs = run("Construction preconditionneur blocs", timers, "preconditionnement", [&]()
                                                       { return construire_inverses_blocs_diagonaux(arbre_blocs_BEM, maillage, quadrature_maillage); });

    /* === Inversion du système Ap = b (GMRES préconditionné par blocs) === */

    Vecteur vect_p_GMRES_BLOC = run("Resolution systeme GMRES blocs", timers, "resolution GMRES blocs", [&]()
                                    {
                                    ResultatGMRES resultat = resoudre_GMRES_BLOC(
                                        arbre_blocs_BEM,
                                        preconditionneur_blocs,
                                        vect_b,
                                        tolGradConj,
                                        maxIterGradConj,
                                        40,  // restart GMRES
                                        true // matrice symétrique
                                    );

                                    if (mpi_rank == 0)
                                    {
                                        affichage_ligne("Itérations", resultat.iterations);
                                        affichage_ligne("Erreur", resultat.erreur);
                                        affichage_ligne("Info", resultat.info);
                                    }

                                    return resultat.solution; });

    /* === Inversion du système Ap = b (méthode itérative) === */

    // Vecteur vect_p_COCG = run("Resolution systeme", timers, "resolution COCG",
    //                           [&]()
    //                           {
    //                               return COCG(maillage, vect_b, quadrature_maillage, arbre_blocs_BEM, tolGradConj, maxIterGradConj, true);
    //                           });

    /* === Erreur avec le vecteur de p sur les milieux des bords === */

    Vecteur vect_p_cached = std::move(vect_p_GMRES_BLOC);

    if (affichageErreur)
        run("Validation de p", timers, "erreur_p", [&]()
            { return calculer_erreur_relative_p(maillage, vect_p_cached); }, [](const Real erreurP)
            { affichage_erreur_relative_p(erreurP); });

    /* === Maillage pour la solution === */

    const std::vector<Point> pointsSolution = run("Maillage pour la solution", timers, "maillage_solution", [&]()
                                                  { return construire_maillage_solution(maillage); }, [&](const std::vector<Point> &pointsSolution)
                                                  {
            if (mpi_rank == 0)
                export_maillage_solution(pointsSolution, output_folder + "/maillage_solution.txt");
            affichage_maillage_solution(pointsSolution.size()); });

    /* === Construction de l'arbre de clusters des points de solution === */

    Arbre arbre_clusters_solution = run("Construction de l'arbre de clusters de la solution", timers, "arbre_clusters_solution", [&]()
                                        { return construire_arbre(pointsSolution, N_leaves); }, [&](const Arbre &arbre)
                                        { ecrire_arbre_clusters(arbre, output_folder + "/arbre_clusters_solution.csv"); });

    /* === Construction de l'arbre de blocs H-matrice de reconstruction === */

    ArbreBlocs arbre_blocs_solution = run("Construction de l'arbre de blocs de reconstruction", timers, "arbre_blocs_solution", [&]()
                                          {
            ArbreBlocs arbre = construire_arbre_blocs(arbre_clusters_solution, arbre_clusters_maillage, parametre_admissibilite, facteur_admissibilite, false);
            remplir_blocs_solution(*arbre, pointsSolution, quadrature_maillage, rang_max, toleranceAca);
            return arbre; }, [&](const ArbreBlocs &arbre)
                                          {
            const IntervallesClusters intervalles_lignes = construire_intervalles_clusters(arbre_clusters_solution);
            const IntervallesClusters intervalles_colonnes = construire_intervalles_clusters(arbre_clusters_maillage);
            affichage_statistiques_matrice_H(arbre);
            ecrire_arbre_blocs(arbre, output_folder + "/arbre_blocs_solution.csv", intervalles_lignes, intervalles_colonnes, false); });

    /* === Reconstruction finale === */

    const Vecteur u_cached = run("Reconstruction de la solution", timers, "reconstruction", [&]()
                                 { return appliquer_H(arbre_blocs_solution, vect_p_cached, false); });

    /* === Validation de la solution  === */

    run("Validation de la solution", timers, "erreurs_sol", [&]()
        { return calculer_erreurs_solution(pointsSolution, u_cached); }, [&](const ErreursSolution &erreurs)
        {
            ecrire_solution(pointsSolution, u_cached, erreurs, output_folder + "/solution.txt");
            if (affichageErreur)
                affichage_erreurs_solution(erreurs.erreurL2, erreurs.erreurMax); });

    // -- End time --

    timers.end("programme");

    if (mpi_rank == 0)
        timers.print_all();

    afficher_memoire();

    // -- liberation mémoire --
    MPI_Win_free(&green_cache_win);
#endif

#ifdef TESTS_CODE

    /* ====================================================== */
    /* INITIALISATION POUR LES TESTS                          */
    /* ====================================================== */

    if (mpi_size != 1)
    {
        if (mpi_rank == 0)
            std::cerr << "Les tests doivent être exécutés avec un seul processus MPI."
                      << std::endl;
    }
    else
    {
        Timers timers;
        timers.start("programme_tests");

        std::cout << "\n=== INITIALISATION DES TESTS ===" << std::endl;
        affichage_parametres();

        /* === Initialisation du cache de Green === */

        run("Initialisation de Green", timers, "green_cache", initialiser_green_cache, affichage_approximations);

        /* === Création du maillage === */

        Maillage maillage = run(
            "Préparation du maillage",
            timers,
            "maillage",
            [&]() -> Maillage
            {
                if (nom_maillage == "narval")
                    return maillage_narval(1.0, pasMaillage, true);

                else if (nom_maillage == "mur_acoustique")
                    return maillage_mur_acoustique(0.5, 4.0, 8, 0.6, pasMaillage, true);

                else if (nom_maillage == "balle_golf")
                    return maillage_balle_de_golf(rayon, 12, 0.12, 0.8, pasMaillage);

                else if (nom_maillage == "rectangle")
                    return maillage_rectangle(2.0, 1.0, pasMaillage);

                else if (nom_maillage == "ellipse")
                    return maillage_ellipse(1.0, 0.5, pasMaillage);

                return maillage_cercle();
            },
            [&](const Maillage &maillage)
            {
                const std::filesystem::path dossier_sortie = std::filesystem::path("outputs") / (maillage.nom + end_output_folder);

                std::filesystem::create_directories(dossier_sortie);
                output_folder = dossier_sortie.string();

                maillage.export_maillage((dossier_sortie / "maillage.txt").string());
            });

        const Index N = static_cast<Index>(maillage.size());

        std::cout << "Nombre de segments : " << N << std::endl;

        /* === Calcul de la quadrature === */

        const std::vector<QuadratureSegment> quadrature_maillage =
            run("Calcul de la quadrature",
                timers,
                "quadrature",
                [&]()
                {
                    return get_quadrature_maillage(maillage, get_legendre_data(ordre));
                });

        /* === Construction de l'arbre de clusters === */

        Arbre arbre_clusters_maillage = run("Construction de l'arbre de clusters",
                                            timers,
                                            "arbre_clusters",
                                            [&]()
                                            {
                                                return construire_arbre_maillage(maillage, N_leaves);
                                            });

        /* === Construction de l'arbre de blocs H-matrice === */

        ArbreBlocs arbre_blocs_BEM = run("Construction de l'arbre de blocs", timers, "arbre_blocs_H", [&]()
                                         { ArbreBlocs arbre = construire_arbre_blocs(arbre_clusters_maillage, parametre_admissibilite, facteur_admissibilite, true);
                                           remplir_blocs_bem(*arbre, maillage, quadrature_maillage, rang_max, toleranceAca);
                                           return arbre; }, [&](const ArbreBlocs &arbre)
                                         { affichage_statistiques_matrice_H(arbre); });

        /* === Construction du préconditionneur de Schwarz === */

        PreconditionneurBlocs preconditionneur_blocs =
            run("Construction du préconditionneur de Schwarz",
                timers,
                "preconditionneur_schwarz",
                [&]()
                {
                    return construire_inverses_blocs_diagonaux(
                        arbre_blocs_BEM,
                        maillage,
                        quadrature_maillage);
                });

        /* === Construction de la matrice dense === */

        timers.start("construction_matrice_dense");

        Matrice A_dense(N, N);

        for (Index i = 0; i < N; ++i)
        {
            for (Index j = i; j < N; ++j)
            {
                const Complex aij = matrice_bem(i, j, maillage, quadrature_maillage);

                A_dense(i, j) = aij;
                A_dense(j, i) = aij;
            }
        }

        timers.end("construction_matrice_dense");

        /* === Construction du second membre === */

        Vecteur vect_b(N);

        timers.start("construction_second_membre");

        for (Index i = 0; i < N; ++i)
        {
            vect_b(i) = integ_simple(
                [](const Point &P)
                {
                    return -u_inc(P);
                },
                quadrature_maillage[i]);
        }

        timers.end("construction_second_membre");

        /* ====================================================== */
        /* TEST 1 : ERREUR RELATIVE MATRICE H / MATRICE DENSE    */
        /* ====================================================== */

        {
            std::cout << "\n=== TEST 1 : Erreur relative H / dense ===" << std::endl;

            const int nb_tests = 20;

            std::cout << "Nombre de vecteurs : " << nb_tests << std::endl;

            /* === Génération et comparaison === */

            std::mt19937_64 generateur(123456789);
            std::normal_distribution<double> distribution(0.0, 1.0);

            std::vector<double> erreurs;
            erreurs.reserve(nb_tests);

            timers.start("comparaison_produits");

            for (int test = 0; test < nb_tests; ++test)
            {
                Vecteur x(N);

                for (Index i = 0; i < N; ++i)
                {
                    x(i) = Complex(distribution(generateur), distribution(generateur));
                }

                const Vecteur y_dense = A_dense * x;
                const Vecteur y_H = appliquer_H(arbre_blocs_BEM, x, true);

                const double denominateur = y_dense.norm();

                if (denominateur == 0.0)
                    continue;

                const double erreur = (y_H - y_dense).norm() / denominateur;

                erreurs.push_back(erreur);
            }

            timers.end("comparaison_produits");

            /* === Affichage des résultats === */

            if (!erreurs.empty())
            {
                const double erreur_min = *std::min_element(erreurs.begin(), erreurs.end());
                const double erreur_max = *std::max_element(erreurs.begin(), erreurs.end());
                const double erreur_moyenne = std::accumulate(erreurs.begin(), erreurs.end(), 0.0) / erreurs.size();

                double variance = 0.0;

                for (double erreur : erreurs)
                    variance += std::pow(erreur - erreur_moyenne, 2);

                variance /= erreurs.size();

                std::cout << "\n--- Résultats ---" << std::endl;
                std::cout << std::scientific << std::setprecision(6);

                std::cout << "Nombre de tests       : "
                          << erreurs.size() << std::endl;

                std::cout << "Erreur minimale      : "
                          << erreur_min << std::endl;

                std::cout << "Erreur moyenne       : "
                          << erreur_moyenne << std::endl;

                std::cout << "Erreur maximale      : "
                          << erreur_max << std::endl;

                std::cout << "Ecart-type           : "
                          << std::sqrt(variance) << std::endl;
            }
        }

        /* ====================================================== */
        /* TEST 2 : PERFORMANCE DES PRODUITS MATRICE-VECTEUR      */
        /* ====================================================== */

        {
            std::cout << "\n=== TEST 2 : Performance des produits Ax ==="
                      << std::endl;

            const int nb_iterations = 100;
            const int nb_vecteurs = 10;

            std::cout << "Nombre de vecteurs : " << nb_vecteurs << std::endl;
            std::cout << "Nombre d'itération : " << nb_iterations << std::endl;

            // Génération des vecteurs avant les mesures
            std::mt19937_64 generateur_perf(987654321);
            std::normal_distribution<double> distribution_perf(0.0, 1.0);

            std::vector<Vecteur> vecteurs(nb_vecteurs);

            for (int k = 0; k < nb_vecteurs; ++k)
            {
                vecteurs[k].resize(N);

                for (Index i = 0; i < N; ++i)
                {
                    vecteurs[k](i) = Complex(distribution_perf(generateur_perf), distribution_perf(generateur_perf));
                }
            }

            // Phase de chauffe
            Vecteur y_dense = A_dense * vecteurs[0];
            Vecteur y_H = appliquer_H(arbre_blocs_BEM, vecteurs[0], true);

            // Vérification de cohérence
            const double erreur_chauffe = (y_H - y_dense).norm() / y_dense.norm();

            // Benchmark matrice dense
            double temps_dense = 0.0;

            {
                const auto debut = std::chrono::steady_clock::now();

                for (int iter = 0; iter < nb_iterations; ++iter)
                {
                    const Vecteur &x = vecteurs[iter % nb_vecteurs];
                    y_dense = A_dense * x;
                }

                const auto fin = std::chrono::steady_clock::now();

                temps_dense = std::chrono::duration<double>(fin - debut).count();
            }

            // Benchmark H-matrice
            double temps_H = 0.0;

            {
                const auto debut = std::chrono::steady_clock::now();

                for (int iter = 0; iter < nb_iterations; ++iter)
                {
                    const Vecteur &x = vecteurs[iter % nb_vecteurs];
                    y_H = appliquer_H(arbre_blocs_BEM, x, true);
                }

                const auto fin = std::chrono::steady_clock::now();

                temps_H = std::chrono::duration<double>(fin - debut).count();
            }

            // Résultats
            std::cout << std::scientific << std::setprecision(6);

            std::cout << "\n--- Résultats ---" << std::endl;

            std::cout << "Nombre de produits : "
                      << nb_iterations << std::endl;

            std::cout << "Temps total dense  : "
                      << temps_dense << " s" << std::endl;

            std::cout << "Temps total H      : "
                      << temps_H << " s" << std::endl;

            std::cout << "Temps moyen dense  : "
                      << temps_dense / nb_iterations << " s" << std::endl;

            std::cout << "Temps moyen H      : "
                      << temps_H / nb_iterations << " s" << std::endl;

            std::cout << "Speedup dense / H  : "
                      << temps_dense / temps_H << std::endl;

            std::cout << "Erreur de contrôle : "
                      << erreur_chauffe << std::endl;
        }

        /* ====================================================== */
        /* TEST 3 : COMPARAISON MEMOIRE DENSE / H-MATRICE         */
        /* ====================================================== */

        {
            std::cout << "\n=== TEST 3 : Comparaison mémoire ===" << std::endl;

            // Mémoire de la matrice dense
            const std::size_t memoire_dense = static_cast<std::size_t>(N) * static_cast<std::size_t>(N) * sizeof(Complex);

            // Comptage mémoire H-matrice
            std::size_t memoire_coefficients_H = 0;
            std::size_t memoire_structure_H = 0;

            std::size_t nb_blocs_denses = 0;
            std::size_t nb_blocs_aca = 0;
            std::size_t nb_noeuds = 0;

            auto compter_memoire_H = [&](auto &&self, const Bloc &bloc) -> void
            {
                ++nb_noeuds;

                // Structure du bloc
                memoire_structure_H += sizeof(Bloc);

                if (bloc.type == TypeBloc::Dense)
                {
                    ++nb_blocs_denses;

                    memoire_coefficients_H += static_cast<std::size_t>(bloc.dense.size()) * sizeof(Complex);
                }
                else if (bloc.type == TypeBloc::Admissible)
                {
                    ++nb_blocs_aca;

                    memoire_coefficients_H += static_cast<std::size_t>(bloc.U.size() + bloc.V.size()) * sizeof(Complex);
                }

                // Parcours des enfants
                for (int i = 0; i < 2; ++i)
                {
                    for (int j = 0; j < 2; ++j)
                    {
                        if (bloc.fils[i][j])
                            self(self, *bloc.fils[i][j]);
                    }
                }
            };

            compter_memoire_H(compter_memoire_H, *arbre_blocs_BEM);

            // Mémoire totale estimée de la H-matrice
            const std::size_t memoire_H = memoire_coefficients_H + memoire_structure_H;

            // Conversion en Mio
            const double dense_MiB = static_cast<double>(memoire_dense) / (1024.0 * 1024.0);
            const double H_coefficients_MiB = static_cast<double>(memoire_coefficients_H) / (1024.0 * 1024.0);
            const double H_structure_MiB = static_cast<double>(memoire_structure_H) / (1024.0 * 1024.0);
            const double H_total_MiB = static_cast<double>(memoire_H) / (1024.0 * 1024.0);

            std::cout << std::fixed << std::setprecision(4);

            std::cout << "\n--- Résultats ---" << std::endl;

            std::cout << "Matrice dense       : "
                      << dense_MiB << " MiB" << std::endl;

            std::cout << "H coefficients      : "
                      << H_coefficients_MiB << " MiB" << std::endl;

            std::cout << "H structure         : "
                      << H_structure_MiB << " MiB" << std::endl;

            std::cout << "H total estime      : "
                      << H_total_MiB << " MiB" << std::endl;

            std::cout << "\n--- Statistiques H ---" << std::endl;

            std::cout << "Nombre de noeuds    : "
                      << nb_noeuds << std::endl;

            std::cout << "Blocs denses        : "
                      << nb_blocs_denses << std::endl;

            std::cout << "Blocs ACA           : "
                      << nb_blocs_aca << std::endl;

            std::cout << "\nRatio memoire dense / H : "
                      << dense_MiB / H_total_MiB << std::endl;
        }

        /* ====================================================== */
        /* TEST 4 : COMPARAISON DES RESOLUTIONS                    */
        /* ====================================================== */

        {
            std::cout << "\n=== TEST 4 : Comparaison des résolutions ==="
                      << std::endl;

            // Paramètres communs aux trois résolutions
            const double tolerance_gmres = tolGradConj;
            const Eigen::Index max_iterations = maxIterGradConj;
            const Eigen::Index restart = 40;

            std::cout << "Taille du système    : " << N << std::endl;
            std::cout << "Tolerance            : " << tolerance_gmres << std::endl;
            std::cout << "Iterations max       : " << max_iterations << std::endl;
            std::cout << "Restart GMRES        : " << restart << std::endl;

            /* ================================================== */
            /* 1. GMRES SUR LA MATRICE DENSE                       */
            /* ================================================== */

            std::cout << "\n--- GMRES dense ---" << std::endl;

            ResultatGMRES resultat_dense;

            timers.start("resolution_GMRES_dense");

            {
                Eigen::GMRES<Matrice, Eigen::IdentityPreconditioner> gmres;

                gmres.set_restart(restart);
                gmres.setMaxIterations(max_iterations);
                gmres.setTolerance(tolerance_gmres);

                gmres.compute(A_dense);

                resultat_dense.solution = gmres.solve(vect_b);
                resultat_dense.iterations = gmres.iterations();
                resultat_dense.erreur = gmres.error();
                resultat_dense.info = gmres.info();
            }

            timers.end("resolution_GMRES_dense");

            // Résidu calculé dans le système original
            const double residu_dense =
                (A_dense * resultat_dense.solution - vect_b).norm() / std::max(vect_b.norm(), 1e-30);

            /* ================================================== */
            /* 2. GMRES SUR LA H-MATRICE SANS PRECONDITIONNEMENT  */
            /* ================================================== */

            std::cout << "\n--- GMRES H sans preconditionnement ---"
                      << std::endl;

            ResultatGMRES resultat_H;

            timers.start("resolution_GMRES_H");

            resultat_H = resoudre_GMRES_Id(
                arbre_blocs_BEM,
                vect_b,
                tolerance_gmres,
                max_iterations,
                restart,
                true);

            timers.end("resolution_GMRES_H");

            // Résidu calculé dans le système original
            const Vecteur residu_H =
                appliquer_H(
                    arbre_blocs_BEM,
                    resultat_H.solution,
                    true) -
                vect_b;

            const double norme_residu_H =
                residu_H.norm() / std::max(vect_b.norm(), 1e-30);

            /* ================================================== */
            /* 3. GMRES SUR LA H-MATRICE AVEC SCHWARZ             */
            /* ================================================== */

            std::cout << "\n--- GMRES H avec Schwarz ---"
                      << std::endl;

            ResultatGMRES resultat_H_schwarz;

            timers.start("resolution_GMRES_H_Schwarz");

            resultat_H_schwarz = resoudre_GMRES_BLOC(
                arbre_blocs_BEM,
                preconditionneur_blocs,
                vect_b,
                tolerance_gmres,
                max_iterations,
                restart,
                true);

            timers.end("resolution_GMRES_H_Schwarz");

            // Attention : on calcule le résidu dans le système original
            // et non dans le système préconditionné.
            const Vecteur residu_H_schwarz =
                appliquer_H(
                    arbre_blocs_BEM,
                    resultat_H_schwarz.solution,
                    true) -
                vect_b;

            const double norme_residu_H_schwarz =
                residu_H_schwarz.norm() / std::max(vect_b.norm(), 1e-30);

            /* ================================================== */
            /* 4. AFFICHAGE DES RESULTATS                         */
            /* ================================================== */

            std::cout << std::scientific << std::setprecision(6);

            std::cout << "\n=== RESULTATS COMPARATIFS ===" << std::endl;

            std::cout << std::left
                      << std::setw(32) << "Methode"
                      << std::right
                      << std::setw(14) << "Iterations"
                      << std::setw(18) << "Residu relatif"
                      << std::setw(12) << "Convergence"
                      << std::endl;

            std::cout << std::string(76, '-') << std::endl;

            auto afficher_resultat =
                [&](const std::string &nom,
                    const ResultatGMRES &resultat,
                    const double residu)
            {
                std::cout << std::left
                          << std::setw(32) << nom
                          << std::right
                          << std::setw(14) << resultat.iterations
                          << std::setw(18) << residu
                          << std::setw(12)
                          << (resultat.info == Eigen::Success ? "OUI" : "NON")
                          << std::endl;
            };

            afficher_resultat(
                "GMRES dense",
                resultat_dense,
                residu_dense);

            afficher_resultat(
                "GMRES H sans Schwarz",
                resultat_H,
                norme_residu_H);

            afficher_resultat(
                "GMRES H avec Schwarz",
                resultat_H_schwarz,
                norme_residu_H_schwarz);

            /* ================================================== */
            /* 5. COMPARAISON DES REDUCTIONS D'ITERATIONS          */
            /* ================================================== */

            std::cout << "\n=== GAINS EN ITERATIONS ===" << std::endl;

            if (resultat_H.iterations > 0)
            {
                std::cout << "Dense / H sans Schwarz : "
                          << static_cast<double>(resultat_dense.iterations) / resultat_H.iterations
                          << std::endl;
            }

            if (resultat_H_schwarz.iterations > 0)
            {
                std::cout << "Dense / H avec Schwarz : "
                          << static_cast<double>(resultat_dense.iterations) / resultat_H_schwarz.iterations
                          << std::endl;

                if (resultat_H.iterations > 0)
                {
                    std::cout << "H sans / avec Schwarz : "
                              << static_cast<double>(resultat_H.iterations) / resultat_H_schwarz.iterations
                              << std::endl;
                }
            }
        }

        timers.end("programme_tests");

        timers.print_all();

        afficher_memoire();

        MPI_Win_free(&green_cache_win);
    }

#endif

    // -- Fin --
    MPI_Finalize();

    return 0;
}
