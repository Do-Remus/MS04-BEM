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
                                    return COCG(maillage, vect_b, quadrature_maillage, tolGradConj, maxIterGradConj);
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

    const std::vector<Complex> u_cached = run("Reconstruction de la solution", timers, "reconstruction", [&]()
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
#endif

#ifdef TESTS_CODE
    /* -------------------- */
    /* -- Area for tests -- */
    /* -------------------- */
#endif

    // -- Fin --
    MPI_Win_free(&green_cache_win);
    MPI_Finalize();

    return 0;
}