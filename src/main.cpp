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

    // Maillage maillage = run("Préparation du maillage de l'obstacle", timers, "maillage",
    //                         [&]()
    //                         {
    //                             return maillage_narval(
    //                                 1.0,    // echelle
    //                                 0.008,  // pas
    //                                 true    // centrer
    //                             );
    //                         },
    //                         [&](const Maillage &maillage)
    //                         {
    //                             maillage.export_maillage("outputs/narval_maillage.txt");

    //                             affichage_maillage(Point(0, 0), maillage.size(), "outputs/narval_maillage.txt"
    //                         );
    //                         });

    // Maillage du mur acoustique

    // Maillage maillage = run("Préparation du maillage de l'obstacle", timers, "maillage",
    //                         [&]()
    //                         {
    //                             return maillage_mur_acoustique(
    //                                 0.5,   // largeur
    //                                 4.0,   // hauteur
    //                                 8,     // nbPointes
    //                                 0.6,   // profondeur
    //                                 0.01,  // pas
    //                                 true   // centrer
    //                             );
    //                         },
    //                         [&](const Maillage &maillage)
    //                         {
    //                             maillage.export_maillage("outputs/mur_acoustique_maillage.txt");

    //                             affichage_maillage(Point(0, 0), maillage.size(), "outputs/mur_acoustique_maillage.txt");
    //                         });

    // Maillage de la balle de golf

    // Maillage maillage = run("Préparation du maillage de l'obstacle", timers, "maillage",
    //                         [&]()
    //                         {
    //                             return maillage_balle_de_golf(
    //                                 rayon,
    //                                 12,    // nbAlveoles
    //                                 0.12,  // profondeur
    //                                 0.8,   // remplissage
    //                                 0.02   // pas
    //                             );
    //                         },
    //                         [&](const Maillage &maillage)
    //                         {
    //                             maillage.export_maillage("outputs/balle_de_golf_maillage.txt");

    //                             affichage_maillage(Point(0, 0), maillage.size(), "outputs/balle_de_golf_maillage.txt");
    //                         });

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
                             return calculer_second_membre(quadrature_maillage, maillage.size());
                         });

    /* === Inversion du système Ap = b (méthode itérative) === */

    Vecteur vect_p_cached = run("Resolution systeme", timers, "resolution",
                                [&]()
                                {
                                    return gradConjMatrixFree(maillage, vect_b, quadrature_maillage, tolGradConj, maxIterGradConj);
                                });

    /* === Erreur avec le vecteur de p sur les milieux des bords === */

    if (mpi_rank == 0)
    {
        Vecteur vect_p_ana(maillage.size());
        for (unsigned int i = 0; i < maillage.size(); i++)
        {
            vect_p_ana[i] = p_analytique(maillage[i].milieu, idxTroncature);
        }

        Real erreurPCached = (vect_p_cached - vect_p_ana).norm() / vect_p_ana.norm();

        std::cout << "    Erreur p cached / p analytique L2 : "
                  << erreurPCached
                  << std::endl;
    }

    /* === Maillage pour la solution === */

    if (mpi_rank == 0)
        std::cout << "\n=== Maillage pour la solution ===" << std::endl;

    timers.start("reconstruction");

    if (mpi_rank == 0)
    {
        // -- Maillage pour la solution --

        vector<Point> pointsSolution;

        const Real xmin = -L / 2.0;
        const Real xmax = L / 2.0;
        const Real ymin = -L / 2.0;
        const Real ymax = L / 2.0;

        // Cercle de reference (sert uniquement a la solution analytique du cercle)
        const Real distanceMin = rayon + delta;

        // -- Boite englobante de la frontiere (pour accelerer les tests) --

        Real bxmin = std::numeric_limits<Real>::max(), bxmax = -bxmin;
        Real bymin = bxmin, bymax = -bxmin;
        for (unsigned int i = 0; i < maillage.size(); i++)
        {
            bxmin = std::min<Real>(bxmin, maillage[i].P1.x);
            bxmax = std::max<Real>(bxmax, maillage[i].P1.x);
            bymin = std::min<Real>(bymin, maillage[i].P1.y);
            bymax = std::max<Real>(bymax, maillage[i].P1.y);
        }

        // -- Point exterieur a l'obstacle (lancer de rayon) ET a au moins delta de la frontiere --

        auto exterieurValide = [&](const Point &M) -> bool
        {
            // Loin de la boite englobante : exterieur et eloigne de la frontiere
            if (M.x < bxmin - delta || M.x > bxmax + delta ||
                M.y < bymin - delta || M.y > bymax + delta)
                return true;

            bool dedans = false;
            Real d2min = std::numeric_limits<Real>::max();

            for (unsigned int i = 0; i < maillage.size(); i++)
            {
                const Point &A = maillage[i].P1;
                const Point &B = maillage[i].P2;

                if ((A.y > M.y) != (B.y > M.y) &&
                    M.x < (B.x - A.x) * (M.y - A.y) / (B.y - A.y) + A.x)
                    dedans = !dedans;

                const Real ex = B.x - A.x;
                const Real ey = B.y - A.y;
                const Real len2 = ex * ex + ey * ey;
                Real t = (len2 > 0.0) ? ((M.x - A.x) * ex + (M.y - A.y) * ey) / len2 : 0.0;
                t = std::max<Real>(0.0, std::min<Real>(1.0, t));
                const Real dx = M.x - (A.x + t * ex);
                const Real dy = M.y - (A.y + t * ey);
                d2min = std::min(d2min, dx * dx + dy * dy);
            }

            return !dedans && d2min >= delta * delta;
        };

        // -- 1) Grille reguliere, filtree par la vraie forme de l'obstacle --

        for (Real x = xmin; x <= xmax; x += pasSolution)
        {
            for (Real y = ymin; y <= ymax; y += pasSolution)
            {
                const Point M(x, y);
                if (exterieurValide(M))
                    pointsSolution.emplace_back(M);
            }
        }

        // -- 2) Couches de points le long de la frontiere (pres de l'obstacle) --
        // Tous les pasProche (abscisse curviligne), on place des points sur la
        // normale sortante a des distances multiples * pasMaillage (>= delta).

        const Real pasProche = 0.02;
        const std::vector<Real> multiplesProches = {2.0, 4.0, 8.0, 16.0, 32.0};
        const unsigned int nbPointsGrille = pointsSolution.size();

        Real arc = pasProche; // force un point sur le premier segment
        for (unsigned int i = 0; i < maillage.size(); i++)
        {
            const Real ex = maillage[i].P2.x - maillage[i].P1.x;
            const Real ey = maillage[i].P2.y - maillage[i].P1.y;
            const Real len = std::sqrt(ex * ex + ey * ey);
            arc += len;

            if (arc < pasProche || len <= 0.0)
                continue;
            arc = 0.0;

            // Normale sortante (contour anti-horaire)
            const Real nx = ey / len;
            const Real ny = -ex / len;
            const Point &C = maillage[i].milieu;

            for (const Real m : multiplesProches)
            {
                const Real d = std::max<Real>(delta, m * pasMaillage);
                const Point M(C.x + d * nx, C.y + d * ny);

                if (M.x < xmin || M.x > xmax || M.y < ymin || M.y > ymax)
                    continue;

                if (exterieurValide(M))
                    pointsSolution.emplace_back(M);
            }
        }

        const unsigned int nbPointsSolution = pointsSolution.size();

        // -- Affichage --

        affichage_maillage_solution(nbPointsGrille, nbPointsSolution);

        /* === Reconstruction fnale === */

        std::cout << "\n=== Reconstruction de la solution ===" << std::endl;

        // -- Creation du fichier de resultats --

        const string filename = string("outputs/cercle_u") + "_k" + std::to_string(k) + "_R" + std::to_string(rayon) + "_L" + std::to_string(L) + "_hM" + std::to_string(pasMaillage) + "_hS" + std::to_string(pasSolution) + "_d" + std::to_string(delta) + "_N" + std::to_string(idxTroncature) + "_q" + std::to_string(ordre) + ".txt";
        ofstream file(filename);

        if (!file.is_open())
        {
            std::cout << "ERROR: Le fichier " << filename << " n'a pas pu être ouvert" << endl;
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }

        // -- Construction solution approchée et des erreurs --

        Real erreurTotaleL2 = 0.0;
        Real erreurTotaleMax = 0.0;
        unsigned int nbPointsErreur = 0; // points ou la solution analytique du cercle est definie

        for (const Point &Pj : pointsSolution)
        {
            Complex uCached = 0.0;

            // -- Solution analytique --

            // (solution analytique du cercle : evitee la ou elle n'a pas de sens)
            const bool horsCercle = std::hypot(Pj.x, Pj.y) >= distanceMin;
            const Complex uExact = horsCercle ? u_N_plus_analytique(Pj, rayon, idxTroncature) : Complex(0.0);

            // -- Interpolation cte --

            for (unsigned int i = 0; i < maillage.size(); i++)
            {
                uCached += vect_p_cached[i] * integ_simple([&Pj](const Point &Q)
                                                           { return green_cached(Pj, Q); },
                                                           quadrature_maillage[i]);
            }

            // Erreurs locales
            const Real erreurTotale = horsCercle ? erreur_relative(uCached, uExact) : 0.0;

            // Accumulation
            if (horsCercle)
            {
                erreurTotaleL2 += erreurTotale * erreurTotale;
                erreurTotaleMax = std::max(erreurTotaleMax, erreurTotale);
                nbPointsErreur++;
            }

            // -- Ecriture Fichier --

            file
                << Pj.x << " "
                << Pj.y << " "

                << uExact.real() << " "
                << uExact.imag() << " "

                << uCached.real() << " "
                << uCached.imag() << " "

                << erreurTotale

                << "\n";
        }

        file.close();

        // -- Erreurs L2 --

        erreurTotaleL2 = std::sqrt(erreurTotaleL2 / static_cast<Real>(std::max(1u, nbPointsErreur)));

        // -- Affichage --

        affichage_erreurs_solution(erreurTotaleL2, erreurTotaleMax);
    }

    timers.end("reconstruction");

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