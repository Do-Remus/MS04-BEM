#include <iostream>
#include <time.h>
#include "src/config/config.hpp"
#include "src/config/constantes.hpp"
#include "src/config/external.hpp"
#include "src/headers/utils.hpp"

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    int rank;
    int size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // -- Initialisation des paramètres --

    const string config = "config.txt";
    get_config(config);

    // -- Affichage Parametres --

    if (rank == 0)
    {
        std::cout << "\n=== Parametres du probleme ===" << std::endl;
        std::cout << "  Nombre d'onde k          : " << k << std::endl;
    }

    // -- Start time --

    MPI_Barrier(MPI_COMM_WORLD);
    const Real start = MPI_Wtime();

    // -- Initialisation cache de Green --

    MPI_Barrier(MPI_COMM_WORLD);
    const Real startGreenCache = MPI_Wtime();

    initialiser_green_cache();

    MPI_Barrier(MPI_COMM_WORLD);
    const Real endGreenCache = MPI_Wtime();

    const Real tempsGreenCached = endGreenCache - startGreenCache;

    if (rank == 0)
    {
        std::cout << "\n=== Approximations ===" << std::endl;
        std::cout << "  Indice de troncature N   : " << idxTroncature << std::endl;
        std::cout << "  Ordre de quadrature      : " << ordre << std::endl;
        std::cout << "  Pas large du cache       : " << green_cache_step_large << std::endl;
        std::cout << "  Pas fin du cache         : " << green_cache_step_fine << std::endl;
        std::cout << "  Cutoff du cache          : " << green_cache_cutoff << std::endl;
        std::cout << "  Taille du cache          : " << green_cache_size << std::endl;
        std::cout << "  Mémoire cache Green      : " << green_cache_size * sizeof(Complex) / (1024.0 * 1024.0) << " MiB" << std::endl;
    }

    // -- Création maillage --

    Point O(0, 0);
    Cercle cercle(rayon, O);
    Maillage maillage;
    maillage.ajoute_cercle(pasMaillage, cercle);
    const unsigned int nbSegmentsMaillage = maillage.size();

    if (rank == 0)
    {
        std::cout << "\n=== Maillage du cercle ===" << std::endl;
        std::cout << "  Centre              : " << O << std::endl;
        std::cout << "  Rayon               : " << rayon << std::endl;
        std::cout << "  Pas                 : " << pasMaillage << std::endl;
        std::cout << "  Nombre de segments  : " << nbSegmentsMaillage << std::endl;
    }

    // -- Récupération des coefficients de Legendre --

    MPI_Barrier(MPI_COMM_WORLD);
    const Real startQuad = MPI_Wtime();

    if (rank == 0)
    {
        std::cout << "\n=== Calcul des points de quadrature ===" << std::endl;
    }

    const LegendreData &legendreData = get_legendre_data(ordre);
    const std::vector<QuadratureSegment> quadrature_maillage = get_quadrature_maillage(maillage, legendreData);

    MPI_Barrier(MPI_COMM_WORLD);
    const Real endQuad = MPI_Wtime();

    const Real tempsQuad = endQuad - startQuad;

    // -- Calcul du vecteur de b de la FV de l'équation intégrale --

    MPI_Barrier(MPI_COMM_WORLD);
    const Real startCalcB = MPI_Wtime();

    if (rank == 0)
    {
        std::cout << "\n=== Calcul du second membre ===" << std::endl;
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

    // -- Inversion de la matrice (méthode itérative) --

    MPI_Barrier(MPI_COMM_WORLD);

    if (rank == 0)
    {
        std::cout << "\n=== Resolution systeme ===" << std::endl;
    }

    const Real startResolutionCached = MPI_Wtime();

    Vecteur vect_p_cached = gradConjMatrixFree(maillage, vect_b, quadrature_maillage, tolGradConj, maxIterGradConj);

    MPI_Barrier(MPI_COMM_WORLD);
    const Real endResolutionCached = MPI_Wtime();

    const Real tempsResolutionCached = endResolutionCached - startResolutionCached;

    // -- Erreur avec le vecteur de p sur les milieux des bords --

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

    Real tempsCached = 0.0;

    if (rank == 0)
    {

        std::cout << "\n=== Reconstruction de la solution ===" << std::endl;

        // -- Maillage pour la solution --

        vector<Point> pointsSolution;

        const Real xmin = -L / 2.0;
        const Real xmax = L / 2.0;
        const Real ymin = -L / 2.0;
        const Real ymax = L / 2.0;

        const Real distanceMin = rayon + delta;

        for (Real x = xmin; x <= xmax; x += pasSolution)
        {
            for (Real y = ymin; y <= ymax; y += pasSolution)
            {
                const Real distance =
                    std::sqrt(x * x + y * y);

                // On conserve uniquement les points
                // suffisamment éloignés du cercle.
                if (distance >= distanceMin)
                {
                    pointsSolution.emplace_back(x, y);
                }
            }
        }

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

        const string filename = string("outputs/u") + "_k" + std::to_string(k) + "_R" + std::to_string(rayon) + "_L" + std::to_string(L) + "_hM" + std::to_string(pasMaillage) + "_hS" + std::to_string(pasSolution) + "_d" + std::to_string(delta) + "_N" + std::to_string(idxTroncature) + "_q" + std::to_string(ordre) + ".txt";
        ofstream file(filename);

        if (!file.is_open())
        {
            std::cout << "ERROR: Le fichier " << filename << " n'a pas pu être ouvert" << endl;
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }

        // -- Construction solution approchée et des erreurs --

        Real erreurTotaleL2 = 0.0;
        Real erreurTotaleMax = 0.0;

        for (const Point &Pj : pointsSolution)
        {
            Complex uCached = 0.0;

            // -- Solution analytique --

            const Complex uExact = u_N_plus_analytique(Pj, rayon, idxTroncature);

            // -- Interpolation cte --

            // Green cached
            const Real startCached = MPI_Wtime();
            ;

            for (unsigned int i = 0; i < nbSegmentsMaillage; i++)
            {
                uCached += vect_p_cached[i] * integ_simple([&Pj](const Point &Q)
                                                           { return green_cached_vec(Pj, Q); },
                                                           quadrature_maillage[i]);
            }

            const Real endCached = MPI_Wtime();

            tempsCached += endCached - startCached;

            // Erreurs locales
            const Real erreurTotale = erreur_relative(uCached, uExact);

            // Accumulation
            erreurTotaleL2 += erreurTotale * erreurTotale;
            erreurTotaleMax = std::max(erreurTotaleMax, erreurTotale);

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

        erreurTotaleL2 = std::sqrt(erreurTotaleL2 / static_cast<Real>(nbPointsSolution));

        // -- Prints --

        std::cout << "\n=== Erreurs ===" << std::endl;

        std::cout << "    P0 erreur L² : Totale = "
                  << erreurTotaleL2
                  << endl;

        std::cout << "    P0 erreur L∞ : Totale = "
                  << erreurTotaleMax
                  << endl;
    }

    // -- End time --

    MPI_Barrier(MPI_COMM_WORLD);
    const Real end = MPI_Wtime();

    Real seconds = end - start;

    if (rank == 0)
    {
        std::cout << "\n=== Temps ===" << std::endl;

        std::cout << "    execution        : "
                  << seconds << " s" << std::endl;

        std::cout << "    creating quadrature   : "
                  << tempsQuad << " s" << std::endl;

        std::cout << "    creation de b         : "
                  << tempsCalcB << " s" << std::endl;

        std::cout << "    creation cache        : "
                  << tempsGreenCached << " s" << std::endl;

        std::cout << "    resolution GC         : "
                  << tempsResolutionCached << " s" << std::endl;

        std::cout << "    reconstruction u      : "
                  << tempsCached << " s" << std::endl;
    }

    // -- Fin --
    MPI_Win_free(&green_cache_win);
    MPI_Finalize();

    return 0;
}