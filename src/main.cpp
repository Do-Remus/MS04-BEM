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

#define TESTS_CODE // FINAL_CODE or TESTS_CODE

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

    /* === Initialisation cache de green === */

    run("Initialisation de Green", timers, "green_cache", initialiser_green_cache, affichage_approximations);

    /* === Création Maillage === */

    // Maillage du Cercle

    Maillage maillage = run("Préparation du maillage de l'obstacle", timers, "maillage", maillage_cercle,
                            [&](const Maillage &maillage)
                            {
                                maillage.export_maillage("outputs/cercle_maillage.txt");
                                affichage_maillage(Point(0, 0), maillage.size(), "outputs/cercle_maillage.txt");
                            });

    // Maillage du Narval

    // Maillage maillage = run("Préparation du maillage de l'obstacle", timers, "maillage", [&]()
    //                         { return maillage_narval(
    //                               1.0,         // echelle
    //                               pasMaillage, // pas
    //                               true         // centrer
    //                           ); }, [&](const Maillage &maillage)
    //                         {
    //                             maillage.export_maillage("outputs/narval_maillage.txt");

    //                             affichage_maillage(Point(0, 0), maillage.size(), "outputs/narval_maillage.txt"
    //                         ); });

    // Maillage du mur acoustique

    // Maillage maillage = run("Préparation du maillage de l'obstacle", timers, "maillage", [&]()
    //                         { return maillage_mur_acoustique(
    //                               0.5,         // largeur
    //                               4.0,         // hauteur
    //                               8,           // nbPointes
    //                               0.6,         // profondeur
    //                               pasMaillage, // pas
    //                               true         // centrer
    //                           ); }, [&](const Maillage &maillage)
    //                         {
    //                             maillage.export_maillage("outputs/mur_acoustique_maillage.txt");

    //                             affichage_maillage(Point(0, 0), maillage.size(), "outputs/mur_acoustique_maillage.txt"); });

    // Maillage de la balle de golf

    // Maillage maillage = run("Préparation du maillage de l'obstacle", timers, "maillage", [&]()
    //                         { return maillage_balle_de_golf(
    //                               rayon,
    //                               12,         // nbAlveoles
    //                               0.12,       // profondeur
    //                               0.8,        // remplissage
    //                               pasMaillage // pas
    //                           ); }, [&](const Maillage &maillage)
    //                         {
    //                             maillage.export_maillage("outputs/balle_de_golf_maillage.txt");

    //                             affichage_maillage(Point(0, 0), maillage.size(), "outputs/balle_de_golf_maillage.txt"); });

    /* === Calcul de la quadrature de Legendre === */

    const std::vector<QuadratureSegment> quadrature_maillage = run("Calcul de la quadrature", timers, "quadrature",
                                                                   [&]()
                                                                   {
                                                                       return get_quadrature_maillage(maillage, get_legendre_data(ordre));
                                                                   });

    /* === Réorganisation de la quadrature pour le SIMD === */

    const QuadratureSIMD quadrature_simd = run("Préparation quadrature SIMD", timers, "quadrature_simd",
                                               [&]()
                                               {
                                                   return construire_quadrature_SIMD(maillage, quadrature_maillage);
                                               });

    /* === Calcul b de la FV de l'équation intégrale === */

    Vecteur vect_b = run("Calcul du second membre de la FV", timers, "b",
                         [&]()
                         {
                             return calculer_second_membre(quadrature_maillage);
                         });

    /* === Inversion du système Ap = b (méthode itérative) === */

    Vecteur vect_p_cached = run("Resolution systeme", timers, "resolution",
                                [&]()
                                {
                                    return COCG(maillage, vect_b, quadrature_maillage, quadrature_simd, tolGradConj, maxIterGradConj);
                                });

    /* === Erreur avec le vecteur de p sur les milieux des bords === */

    run("Validation de p", timers, "erreur_p", [&]()
        { return calculer_erreur_relative_p(maillage, vect_p_cached); }, [](const Real erreurP)
        { affichage_erreur_relative_p(erreurP); });

    /* === Maillage pour la solution === */

    const std::vector<Point> pointsSolution = run("Maillage pour la solution", timers, "maillage_solution", [&]()
                                                  { return construire_maillage_solution(maillage); }, [&](const std::vector<Point> &pointsSolution)
                                                  { affichage_maillage_solution(pointsSolution.size()); });

    /* === Reconstruction finale === */

    const Vecteur u_cached = run("Reconstruction de la solution", timers, "reconstruction", [&]()
                                 { return reconstruire_solution(pointsSolution, maillage, quadrature_maillage, vect_p_cached); });

    /* === Validation de la solution  === */

    run("Validation de la solution", timers, "erreurs_sol", [&]()
        { return calculer_erreurs_solution(pointsSolution, u_cached); }, [&](const ErreursSolution &erreurs)
        { ecrire_solution(pointsSolution, u_cached, erreurs);
          affichage_erreurs_solution(erreurs.erreurL2, erreurs.erreurMax); });

    // -- End time --

    timers.end("programme");

    if (mpi_rank == 0)
        timers.print_all();

    // -- liberation mémoire --
    MPI_Win_free(&green_cache_win);
#endif

#ifdef TESTS_CODE
    /* -------------------- */
    /* -- Area for tests -- */
    /* -------------------- */

    Timers timers;

    // Matrice A = matrice_aleatoire(maxRank, N);
    Maillage maillage = run("Préparation du maillage de l'obstacle", timers, "maillage", maillage_cercle,
                            [&](const Maillage &maillage)
                            {
                                maillage.export_maillage("outputs/cercle_maillage.txt");
                                affichage_maillage(Point(0, 0), maillage.size(), "outputs/cercle_maillage.txt");
                            });

    // const std::vector<QuadratureSegment> quadrature_maillage = run("Calcul de la quadrature", timers, "quadrature",
    //                                                                [&]()
    //                                                                {
    //                                                                    return get_quadrature_maillage(maillage, get_legendre_data(ordre));
    //                                                                });
    // const Eigen::Index N = maillage.size();
    // const Eigen::Index maxRank = N;
    // // const Eigen::Index rankTroncature = maxRank;
    // for (Eigen::Index rankTroncature = 10; rankTroncature < N + 1; rankTroncature += 30)
    // {
    //     std::cout << " N =" << N << std::endl;
    //     Matrice A = get_matrice_bem(maillage, quadrature_maillage);
    //     // Vérification des dimensions
    //     std::cout << "Dimensions de A : "
    //               << A.rows() << " x " << A.cols() << std::endl;

    //     std::cout << "SVD..." << std::endl;

    //     // Pour une matrice complexe, utiliser JacobiSVD
    //     timers.start("svd");
    //     Eigen::JacobiSVD<Matrice> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
    //     std::cout << "Troncature..." << std::endl;

    //     const auto &sigma = svd.singularValues();

    //     const Eigen::Index r = std::min(rankTroncature, static_cast<Eigen::Index>(sigma.size()));

    //     if (r < 0)
    //         throw std::invalid_argument("Rang de troncature invalide");

    //     const Real sigmaMax = sigma.size() > 0 ? sigma(0) : 0.0;

    //     const Real tol = std::numeric_limits<Real>::epsilon() * std::max(A.rows(), A.cols()) * sigmaMax;

    //     Eigen::Index rangReel = 0;

    //     for (Eigen::Index i = 0; i < sigma.size(); ++i)
    //         if (sigma(i) > tol)
    //             ++rangReel;

    //     timers.end("svd");

    //     // SVD tronquée
    //     Matrice U_r = svd.matrixU().leftCols(r);
    //     Matrice V_r = svd.matrixV().leftCols(r);
    //     auto sigma_r = sigma.head(r);

    //     Matrice A_r = U_r * sigma_r.asDiagonal() * V_r.adjoint();

    //     const Real normeA = A.norm();
    //     const Real erreurAbsolue = (A - A_r).norm();
    //     const Real erreurRelative =
    //         normeA > 0.0 ? erreurAbsolue / normeA : 0.0;

    //     std::cout << "ACA..." << std::endl;

    //     // -- Test de ACA pivotage complet
    //     timers.start("ACA_pivot_complet");
    //     ACA_decomposition aca = full_pivoted_cross_approximation(A, maxRank, 1e-8);
    //     timers.end("ACA_pivot_complet");
    //     std::cout << "ACA Partiel..." << std::endl;

    //     // -- Test de ACA pivotage partiel
    //     timers.start("ACA_pivot_partiel");
    //     ACA_decomposition aca_partial = partial_pivoted_cross_approximation(A, maxRank, 1e-8);
    //     timers.end("ACA_pivot_partiel");
    //     std::cout << "Préaffichage..." << std::endl;

    //     Matrice A_aca = aca.U * aca.V.transpose();
    //     Matrice A_aca_partial = aca_partial.U * aca_partial.V.transpose();

    //     // -- Erreurs
    //     const Real erreurAbsolueACA = (A - A_aca).norm();
    //     const Real erreurRelativeACA =
    //         normeA > 0.0 ? erreurAbsolueACA / normeA : 0.0;

    //     const Real erreurAbsolueACAPartial = (A - A_aca_partial).norm();
    //     const Real erreurRelativeACAPartial =
    //         normeA > 0.0 ? erreurAbsolueACAPartial / normeA : 0.0;

    //     // -- Comparaison de mémoire SVD tronqué
    //     const std::size_t originalBytes = static_cast<std::size_t>(A.size()) * sizeof(Complex);

    //     const std::size_t truncatedBytes = static_cast<std::size_t>(U_r.size()) * sizeof(Complex) + static_cast<std::size_t>(sigma_r.size()) * sizeof(Real) + static_cast<std::size_t>(V_r.size()) * sizeof(Complex);

    //     const double originalMB = originalBytes / (1024.0 * 1024.0);
    //     const double truncatedMB = truncatedBytes / (1024.0 * 1024.0);

    //     const double savedBytes = static_cast<double>(originalBytes) - truncatedBytes;

    //     const double savedMB = savedBytes / (1024.0 * 1024.0);

    //     const double savedPercent = originalBytes > 0 ? 100.0 * savedBytes / originalBytes : 0.0;

    //     std::cout << std::scientific << std::setprecision(6);

    //     // -- Comparaison de mémoire ACA
    //     auto memoireACA = [](const ACA_decomposition &d) -> std::size_t
    //     {
    //         return static_cast<std::size_t>(
    //                    d.U.size() + d.V.size()) *
    //                sizeof(Complex);
    //     };

    //     const std::size_t memoireACAComplet = memoireACA(aca);
    //     const std::size_t memoireACAPartiel = memoireACA(aca_partial);

    //     auto MiB = [](std::size_t bytes) -> double
    //     {
    //         return static_cast<double>(bytes) / (1024.0 * 1024.0);
    //     };

    //     auto tauxCompression = [&](std::size_t bytes) -> double
    //     {
    //         return bytes > 0
    //                    ? 100.0 * (1.0 -
    //                               static_cast<double>(bytes) /
    //                                   static_cast<double>(originalBytes))
    //                    : 0.0;
    //     };

    //     std::cout << "\n========== ANALYSE SVD ==========\n";
    //     std::cout << "Dimensions             : "
    //               << A.rows() << " x " << A.cols() << '\n';
    //     std::cout << "Rang demande            : " << maxRank << '\n';
    //     std::cout << "Rang numerique estime   : " << rangReel << '\n';
    //     std::cout << "Rang de troncature      : " << r << '\n';
    //     std::cout << "Tolerance de rang       : " << tol << '\n';
    //     std::cout << "Temps de calcul         : ";
    //     timers.print("svd");
    //     std::cout << '\n';

    //     std::cout << "\n========== ERREUR ==========\n";
    //     std::cout << "Norme de A              : " << normeA << '\n';
    //     std::cout << "Norme de A_r            : " << A_r.norm() << '\n';
    //     std::cout << "Erreur absolue          : " << erreurAbsolue << '\n';
    //     std::cout << "Erreur relative         : " << erreurRelative << '\n';
    //     std::cout << "================================\n";

    //     std::cout << "\n========== MEMOIRE SVD ==========\n";
    //     std::cout << "Original matrix       : "
    //               << originalMB << " MiB\n";
    //     std::cout << "Truncated SVD factors : "
    //               << truncatedMB << " MiB\n";
    //     std::cout << "Memory saved          : "
    //               << savedMB << " MiB\n";
    //     std::cout << "Space reduction       : "
    //               << savedPercent << " %\n";
    //     std::cout << "====================================\n";

    //     std::cout << "\n========== ANALYSE ACA ==========\n";

    //     std::cout << "\n--- ACA : pivotage complet ---\n";
    //     std::cout << "Rang obtenu             : " << aca.rang << '\n';
    //     std::cout << "Erreur absolue          : " << erreurAbsolueACA << '\n';
    //     std::cout << "Erreur relative         : " << erreurRelativeACA << '\n';
    //     std::cout << "Memoire des facteurs    : "
    //               << MiB(memoireACAComplet) << " MiB\n";
    //     std::cout << "Reduction memoire       : "
    //               << tauxCompression(memoireACAComplet) << " %\n";
    //     std::cout << "Temps de calcul         : ";
    //     timers.print("ACA_pivot_complet");
    //     std::cout << '\n';

    //     std::cout << "\n--- ACA : pivotage partiel ---\n";
    //     std::cout << "Rang obtenu             : " << aca_partial.rang << '\n';
    //     std::cout << "Erreur absolue          : "
    //               << erreurAbsolueACAPartial << '\n';
    //     std::cout << "Erreur relative         : "
    //               << erreurRelativeACAPartial << '\n';
    //     std::cout << "Memoire des facteurs    : "
    //               << MiB(memoireACAPartiel) << " MiB\n";
    //     std::cout << "Reduction memoire       : "
    //               << tauxCompression(memoireACAPartiel) << " %\n";

    //     std::cout << "Temps de calcul         : ";
    //     timers.print("ACA_pivot_partiel");
    //     std::cout << '\n';

    //     std::cout << "=================================\n";
    // }

    const Index N_leafs = 100;
    const double eta = 3.0;

    Arbre arbre_clusters = run("Construction de l'arbre de clusters", timers, "arbre_clusters", [&]()
                               { return construire_arbre_maillage(maillage, N_leafs); }, [&](const Arbre &arbre)
                               { ecrire_arbre_clusters(arbre, "outputs/arbre_clusters.csv"); });

    std::vector<Index> permutation = run("Construction de la permutation", timers, "permutation",
                                         [&]()
                                         {
                                             return construire_permutation(arbre_clusters);
                                         });

    std::vector<Index> permutation_inverse = run("Construction de la permutation inverse", timers, "permutation_inverse",
                                                 [&]()
                                                 {
                                                     return construire_permutation_inverse(permutation);
                                                 });

    IntervallesClusters intervalles = run("Construction des intervalles des clusters", timers, "intervalles_clusters",
                                          [&]()
                                          {
                                              return construire_intervalles_clusters(arbre_clusters);
                                          });

    ArbreBlocs arbre_blocs = run("Construction de l'arbre des blocs", timers, "arbre_blocs", [&]()
                                 { return construire_arbre_blocs(arbre_clusters, eta); }, [&](const ArbreBlocs &arbre)
                                 { ecrire_arbre_blocs(arbre, "outputs/arbre_blocs.csv", intervalles); });

#endif

    // -- Fin --
    MPI_Finalize();

    return 0;
}