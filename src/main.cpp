#include "global/global.hpp"
#include "BEM/free_matrix/commun.hpp"
#include "geometrie/commun.hpp"
#include "functions/green.hpp"
#include "functions/solutions_analytiques.hpp"

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    int rank;
    int size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

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

    /* === Début Programme === */

    if (rank == 0)
    {
        std::cout << "\n=== Parametres du probleme ===" << std::endl;
        std::cout << "  Nombre d'onde k          : " << k << std::endl;
    }

    MPI_Barrier(MPI_COMM_WORLD);
    const Real start = MPI_Wtime();

    /* === Initialisation Cache Green === */

    MPI_Barrier(MPI_COMM_WORLD);
    const Real startGreenCache = MPI_Wtime();

    if (rank == 0)
    {
        std::cout << "\n=== Initialisation de Green ===" << std::endl;
    }

    initialiser_green_cache();

    MPI_Barrier(MPI_COMM_WORLD);
    const Real endGreenCache = MPI_Wtime();

    const Real tempsGreenCached = endGreenCache - startGreenCache;

    if (rank == 0)
    {
        std::cout << "  Indice de troncature N   : " << idxTroncature << std::endl;
        std::cout << "  Ordre de quadrature      : " << ordre << std::endl;
        std::cout << "  Pas large du cache       : " << green_cache_step_large << std::endl;
        std::cout << "  Pas fin du cache         : " << green_cache_step_fine << std::endl;
        std::cout << "  Cutoff du cache          : " << green_cache_cutoff << std::endl;
        std::cout << "  Taille du cache          : " << green_cache_size << std::endl;
        std::cout << "  Mémoire cache Green      : " << green_cache_size * sizeof(Complex) / (1024.0 * 1024.0) << " MiB" << std::endl;
        std::cout << "  Temps de construction    : " << tempsGreenCached << " s" << endl;
    }

    /* === Création Maillage === */

    MPI_Barrier(MPI_COMM_WORLD);
    const Real startMaillage = MPI_Wtime();

    if (rank == 0)
    {
        std::cout << "\n=== Préparation du maillage de l'obstacle ===" << std::endl;
    }

    Point O(0, 0);
    Cercle cercle(rayon, O);

    Maillage maillage;
    maillage.ajoute_cercle(pasMaillage, cercle);

    // Autres maillages:

    /* //balle de golf
    vector<Point> pts = pointsBalleGolf(1.0,    // rayon R
                                    10,     // nombre d'alvéoles
                                    0.12,   // profondeur = 12 % de R (amplifié)
                                    0.8,    // remplissage angulaire d'une alvéole
                                    pasMaillage); // pas du maillage
    */

    /* // Narval
    const vector<Point> pts = pointsNarval(1.0, pasMaillage);
    Maillage maillage(pts);
    */

    const unsigned int nbSegmentsMaillage = maillage.size();

    MPI_Barrier(MPI_COMM_WORLD);
    const Real endMaillage = MPI_Wtime();

    const Real tempsMaillage = endMaillage - startMaillage;

    if (rank == 0)
    {
        // -- Export du maillage (un seul processus ecrit) --
        maillage.export_maillage("outputs/cercle_maillage.txt");

        std::cout << "  Centre                 : " << O << std::endl;
        std::cout << "  Rayon                  : " << rayon << std::endl;
        std::cout << "  Pas                    : " << pasMaillage << std::endl;
        std::cout << "  Nombre de segments     : " << nbSegmentsMaillage << std::endl;
        std::cout << "  Maillage exporte       : outputs/cercle_maillage.txt" << std::endl;
        std::cout << "  Temps de fabrication   : " << tempsMaillage << " s" << std::endl;
    }

    /* === Récupération Coefficients Legendre === */

    MPI_Barrier(MPI_COMM_WORLD);
    const Real startQuad = MPI_Wtime();

    if (rank == 0)
    {
        std::cout << "\n=== Calcul de la quadrature ===" << std::endl;
    }

    const LegendreData &legendreData = get_legendre_data(ordre);
    const std::vector<QuadratureSegment> quadrature_maillage = get_quadrature_maillage(maillage, legendreData);

    MPI_Barrier(MPI_COMM_WORLD);
    const Real endQuad = MPI_Wtime();

    const Real tempsQuad = endQuad - startQuad;

    if (rank == 0)
    {
        std::cout << "  Temps de fabrication   : " << tempsQuad << " s" << std::endl;
    }

    /* === Calcul b de la FV de l'équation intégrale === */

    MPI_Barrier(MPI_COMM_WORLD);
    const Real startCalcB = MPI_Wtime();

    if (rank == 0)
    {
        std::cout << "\n=== Calcul du second membre de la FV ===" << std::endl;
    }

    Vecteur b_local(nbSegmentsMaillage, 0.0);
    for (unsigned int i = rank; i < nbSegmentsMaillage; i += size)
    {
        b_local[i] = integ_simple([](const Point &Q)
                                  { return -u_inc(Q); },
                                  quadrature_maillage[i]);
    }

    Vecteur vect_b(nbSegmentsMaillage, 0.0);

    MPI_Allreduce(b_local.data(), vect_b.data(), static_cast<int>(nbSegmentsMaillage), mpi_complex_type(), MPI_SUM, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);
    const Real endCalcB = MPI_Wtime();

    const Real tempsCalcB = endCalcB - startCalcB;

    if (rank == 0)
    {
        std::cout << "  Temps de fabrication   : " << tempsCalcB << " s" << std::endl;
    }

    /* === Inversion de la matrice (méthode itérative) === */

    MPI_Barrier(MPI_COMM_WORLD);
    const Real startResolutionCached = MPI_Wtime();

    if (rank == 0)
    {
        std::cout << "\n=== Resolution systeme ===" << std::endl;
    }

    Vecteur vect_p_cached = gradConjMatrixFree(maillage, vect_b, quadrature_maillage, tolGradConj, maxIterGradConj);

    MPI_Barrier(MPI_COMM_WORLD);
    const Real endResolutionCached = MPI_Wtime();

    const Real tempsResolutionCached = endResolutionCached - startResolutionCached;

    if (rank == 0)
    {
        std::cout << "  Temps de résolution   : " << tempsResolutionCached << std::endl;
    }

    /* === Erreur avec le vecteur de p sur les milieux des bords === */

    if (rank == 0)
    {
        Vecteur vect_p_ana(nbSegmentsMaillage);
        for (unsigned int i = 0; i < nbSegmentsMaillage; i++)
        {
            vect_p_ana[i] = p_analytique(maillage[i].milieu, idxTroncature);
        }

        Real erreurPCached = (vect_p_cached - vect_p_ana).norm() / vect_p_ana.norm();

        std::cout << "\n=== Analyse de p ===" << std::endl;

        std::cout << "  Erreur p cached / p analytique L2 : "
                  << erreurPCached
                  << std::endl;
    }

    // -- Reconstruction de la solution --

    MPI_Barrier(MPI_COMM_WORLD);
    const Real startReconstruction = MPI_Wtime();

    if (rank == 0)
    {

        std::cout << "\n=== Reconstruction de la solution ===" << std::endl;

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
        for (unsigned int i = 0; i < nbSegmentsMaillage; i++)
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

            for (unsigned int i = 0; i < nbSegmentsMaillage; i++)
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
        for (unsigned int i = 0; i < nbSegmentsMaillage; i++)
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

        std::cout << "  Points de grille    : " << nbPointsGrille << std::endl;
        std::cout << "  Points pres du bord : " << pointsSolution.size() - nbPointsGrille << std::endl;

        const unsigned int nbPointsSolution = pointsSolution.size();

        // -- Affichage --

        std::cout << "\n=== Maillage de solution ===" << std::endl;
        std::cout << "  Domaine             : [" << xmin << ", " << xmax << "] x ["
                  << ymin << ", " << ymax << "]" << std::endl;
        std::cout << "  Taille du domaine   : " << L << " x " << L << std::endl;
        std::cout << "  Pas                 : " << pasSolution << std::endl;
        std::cout << "  Distance à la frontière : " << delta << std::endl;
        std::cout << "  Rayon minimal autorisé   : " << distanceMin << std::endl;
        std::cout << "  Nombre de points    : " << nbPointsSolution << std::endl;

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

            for (unsigned int i = 0; i < nbSegmentsMaillage; i++)
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

        // -- Prints --

        std::cout << "\n=== Erreurs ===" << std::endl;

        std::cout << "    P0 erreur L² : Totale = "
                  << erreurTotaleL2
                  << endl;

        std::cout << "    P0 erreur L∞ : Totale = "
                  << erreurTotaleMax
                  << endl;
    }

    MPI_Barrier(MPI_COMM_WORLD);
    const Real endReconstruction = MPI_Wtime();

    const Real tempsReconstruction = endReconstruction - startReconstruction;

    // -- End time --

    MPI_Barrier(MPI_COMM_WORLD);
    const Real end = MPI_Wtime();

    const Real tempsTotal = end - start;

    if (rank == 0)
    {
        std::cout << "\n=== Temps ===" << std::endl;

        std::cout << "    execution totale      : "
                  << tempsTotal << " s" << std::endl;

        std::cout << "    creating maillage     : "
                  << tempsMaillage << " s" << std::endl;

        std::cout << "    creating quadrature   : "
                  << tempsQuad << " s" << std::endl;

        std::cout << "    creation de b         : "
                  << tempsCalcB << " s" << std::endl;

        std::cout << "    creation cache        : "
                  << tempsGreenCached << " s" << std::endl;

        std::cout << "    resolution GC         : "
                  << tempsResolutionCached << " s" << std::endl;

        std::cout << "    reconstruction u      : "
                  << tempsReconstruction << " s" << std::endl;
    }
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