#include <iostream>
#include <time.h>
#include "src/config/config.hpp"
#include "src/config/constantes.hpp"
#include "src/config/external.hpp"
#include "src/headers/utils.hpp"

Vecteur produit_A_cached_blocked(const Maillage &maillage, const Vecteur &x, const std::vector<QuadratureSegment> &quadrature_maillage)
{
    const std::size_t N = maillage.size();

    int rank;
    int size;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    Vecteur y_local(N, 0.0);
    constexpr std::size_t BLOCK = 8; // Modifiable pour optimiser

    // Buffers réutilisés
    std::vector<std::complex<double>> yi(BLOCK);
    std::vector<std::complex<double>> yj(BLOCK);

    for (std::size_t ib = rank * BLOCK; ib < N; ib += size * BLOCK)
    {
        const std::size_t i_end = std::min(ib + BLOCK, N);
        const std::size_t ni = i_end - ib;

        // Remise à zéro du bloc i
        std::fill(yi.begin(), yi.begin() + ni, std::complex<double>(0.0));

        for (std::size_t jb = ib; jb < N; jb += BLOCK)
        {
            const std::size_t j_end = std::min(jb + BLOCK, N);
            const std::size_t nj = j_end - jb;

            // Remise à zéro du bloc j
            std::fill(yj.begin(), yj.begin() + nj, std::complex<double>(0.0));

            for (std::size_t i = ib; i < i_end; ++i)
            {
                const QuadratureSegment &qS = quadrature_maillage[i];

                const std::size_t j_start = std::max(jb, i + 1);

                for (std::size_t j = j_start; j < j_end; ++j)
                {
                    const std::complex<double> aij = integ_double([](const Point &Q,
                                                                     const Point &P)
                                                                  { return green_cached_vec(Q, P); },
                                                                  qS,
                                                                  quadrature_maillage[j]);

                    yi[i - ib] += aij * x[j];
                    yj[j - jb] += aij * x[i];
                }
            }

            // Une seule écriture contiguë du bloc j
            for (std::size_t j = jb; j < j_end; ++j)
            {
                y_local[j] += yj[j - jb];
            }
        }

        // Écriture contiguë du bloc i
        for (std::size_t i = ib; i < i_end; ++i)
        {
            y_local[i] += yi[i - ib];
        }
    }

    Vecteur y(N, 0.0);

    MPI_Allreduce(y_local.data(), y.data(), static_cast<int>(N), MPI_C_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);

    return y;
}

Vecteur produit_A_cached(const Maillage &maillage, const Vecteur &x, const std::vector<QuadratureSegment> &quadrature_maillage)
{
    const unsigned int N = maillage.size();

    int rank;
    int size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    /*
     * Chaque processus possède son propre vecteur local.
     *
     * Il contient les contributions calculées par ce processus
     * pour toutes les composantes.
     */
    Vecteur y_local(N, 0.0);

    /*
     * Répartition cyclique des lignes.
     *
     * rank 0 : 0, size, 2*size, ...
     * rank 1 : 1, size+1, 2*size+1, ...
     */
    for (unsigned int i = rank; i < N; i += size)
    {
        const Segment &S = maillage[i];
        const QuadratureSegment &qS = quadrature_maillage[i];

        std::complex<double> aii =
            integ_simple([&S](const Point &Q)
                         { return integ_simple_log(Q, S); },
                         qS);

        aii += integ_double([](const Point &Q, const Point &P)
                            { return green_reguliere_cached_vec(Q, P); },
                            qS,
                            qS);

        y_local[i] += aii * x[i];

        for (unsigned int j = i + 1; j < N; ++j)
        {
            const std::complex<double> aij = integ_double([](const Point &Q, const Point &P)
                                                          { return green_cached_vec(Q, P); },
                                                          qS,
                                                          quadrature_maillage[j]);

            y_local[i] += aij * x[j];
            y_local[j] += aij * x[i];
        }
    }

    /*
     * Chaque rang possède maintenant une partie des contributions.
     *
     * On somme les contributions de tous les rangs.
     */
    Vecteur y(N, 0.0);

    MPI_Allreduce(y_local.data(), y.data(), static_cast<int>(N), MPI_C_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);

    return y;
}

Vecteur gradConjMatrixFree(const Maillage &maillage, const Vecteur &b, const std::vector<QuadratureSegment> &quadrature_maillage, double tol, unsigned int maxIter)
{
    const std::size_t n = b.size();
    Vecteur x(n, 0.0);

    Vecteur r = b; // car x0 = 0
    Vecteur d = r;

    const double bnorm = b.norm();

    if (bnorm == 0.0)
        return x;

    double relativeResidual = 1.0;

    std::complex<double> rho = r.produitBilineaire(r);

    for (unsigned int iter = 0; iter < maxIter; ++iter)
    {
        Vecteur Ad = produit_A_cached_blocked(maillage, d, quadrature_maillage);
        const std::complex<double> denom = d.produitBilineaire(Ad);
        const double scale = d.norm() * Ad.norm();

        if (scale == 0.0 || std::abs(denom) < 1e-20 * scale)
        {
            int rank;
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);

            if (rank == 0)
            {
                std::cerr
                    << "COCG : breakdown, "
                    << "denominateur nul "
                    << "a l'iteration "
                    << iter
                    << ", résidu relatif = "
                    << relativeResidual
                    << std::endl;
            }

            return x;
        }

        const std::complex<double> alpha = rho / denom;

        for (std::size_t i = 0; i < n; ++i)
        {
            x[i] += alpha * d[i];
            r[i] -= alpha * Ad[i];
        }

        relativeResidual = r.norm() / bnorm;

        if (relativeResidual < tol)
        {
            int rank;
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);

            if (rank == 0)
            {
                std::cout
                    << "COCG converge en "
                    << iter + 1
                    << " iterations, "
                    << "résidu relatif = "
                    << relativeResidual
                    << std::endl;
            }

            return x;
        }

        const std::complex<double> rhoNew = r.produitBilineaire(r);

        if (std::abs(rho) < 1e-30)
        {
            int rank;
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);

            if (rank == 0)
            {
                std::cerr
                    << "COCG : breakdown de rho "
                    << "a l'iteration "
                    << iter
                    << std::endl;
            }

            return x;
        }

        const std::complex<double> beta = rhoNew / rho;

        for (std::size_t i = 0; i < n; ++i)
        {
            d[i] = r[i] + beta * d[i];
        }

        rho = rhoNew;
    }

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0)
    {
        std::cout
            << "COCG : nombre maximal "
            << "d'iterations atteint."
            << std::endl;
    }

    return x;
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    int rank;
    int size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // initialisation des paramètres
    const string config = "config.txt";
    get_config(config);

    // -- Affichage Parametres --
    if (rank == 0)
    {
        std::cout << "\n=== Parametres du probleme ===" << std::endl;
        std::cout << "  Nombre d'onde k          : " << k << std::endl;

        std::cout << "\n=== Approximations ===" << std::endl;
        std::cout << "  Indice de troncature N   : " << idxTroncature << std::endl;
        std::cout << "  Ordre de quadrature      : " << ordre << std::endl;
        std::cout << "  Pas du cache de Green    : " << green_cache_step << std::endl;
        std::cout << "  Nombre max d'indices     : " << max_index << std::endl;
    }

    // -- Start time --

    MPI_Barrier(MPI_COMM_WORLD);
    const double start = MPI_Wtime();

    // -- Initialisation cache de Green --

    MPI_Barrier(MPI_COMM_WORLD);
    const double startGreenCache = MPI_Wtime();

    initialiser_green_cache();

    MPI_Barrier(MPI_COMM_WORLD);
    const double endGreenCache = MPI_Wtime();

    const double tempsGreenCached = endGreenCache - startGreenCache;

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
    const double startQuad = MPI_Wtime();

    if (rank == 0)
    {
        std::cout << "\n=== Calcul des points de quadrature ===" << std::endl;
    }

    const LegendreData &legendreData = get_legendre_data(ordre);
    const std::vector<QuadratureSegment> quadrature_maillage = get_quadrature_maillage(maillage, legendreData);

    MPI_Barrier(MPI_COMM_WORLD);
    const double endQuad = MPI_Wtime();

    const double tempsQuad = endQuad - startQuad;

    // -- Calcul du vecteur de b de la FV de l'équation intégrale --

    MPI_Barrier(MPI_COMM_WORLD);
    const double startCalcB = MPI_Wtime();

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

    MPI_Allreduce(b_local.data(), vect_b.data(), static_cast<int>(nbSegmentsMaillage), MPI_C_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);
    const double endCalcB = MPI_Wtime();

    const double tempsCalcB = endCalcB - startCalcB;

    // -- Inversion de la matrice (méthode itérative) --

    MPI_Barrier(MPI_COMM_WORLD);

    if (rank == 0)
    {
        std::cout << "\n=== Resolution systeme ===" << std::endl;
    }

    const double startResolutionCached = MPI_Wtime();

    Vecteur vect_p_cached = gradConjMatrixFree(maillage, vect_b, quadrature_maillage, tolGradConj, maxIterGradConj);

    MPI_Barrier(MPI_COMM_WORLD);
    const double endResolutionCached = MPI_Wtime();

    const double tempsResolutionCached = endResolutionCached - startResolutionCached;

    // -- Erreur avec le vecteur de p sur les milieux des bords --

    if (rank == 0)
    {
        Vecteur vect_p_ana(nbSegmentsMaillage);
        for (unsigned int i = 0; i < nbSegmentsMaillage; i++)
        {
            vect_p_ana[i] = p_analytique(maillage[i].milieu, idxTroncature);
        }

        double erreurPCached = (vect_p_cached - vect_p_ana).norm() / vect_p_ana.norm();

        std::cout << "\n=== Analyse de p ===" << std::endl;

        std::cout << "  Erreur p cached / p analytique L2 : "
                  << erreurPCached
                  << std::endl;
    }

    // -- Reconstruction de la solution --

    double tempsCached = 0.0;

    if (rank == 0)
    {

        std::cout << "\n=== Reconstruction de la solution ===" << std::endl;

        // -- Maillage pour la solution --

        vector<Point> pointsSolution;

        const double xmin = -L / 2.0;
        const double xmax = L / 2.0;
        const double ymin = -L / 2.0;
        const double ymax = L / 2.0;

        const double distanceMin = rayon + delta;

        for (double x = xmin; x <= xmax; x += pasSolution)
        {
            for (double y = ymin; y <= ymax; y += pasSolution)
            {
                const double distance =
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

        double erreurTotaleL2 = 0.0;
        double erreurTotaleMax = 0.0;

        for (const Point &Pj : pointsSolution)
        {
            complex<double> uCached = 0.0;

            // -- Solution analytique --

            const complex<double> uExact = u_N_plus_analytique(Pj, rayon, idxTroncature);

            // -- Interpolation cte --

            // Green cached
            const double startCached = MPI_Wtime();
            ;

            for (unsigned int i = 0; i < nbSegmentsMaillage; i++)
            {
                uCached += vect_p_cached[i] * integ_simple([&Pj](const Point &Q)
                                                           { return green_cached_vec(Pj, Q); },
                                                           quadrature_maillage[i]);
            }

            const double endCached = MPI_Wtime();

            tempsCached += endCached - startCached;

            // Erreurs locales
            const double erreurTotale = erreur_relative(uCached, uExact);

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

        erreurTotaleL2 = std::sqrt(erreurTotaleL2 / static_cast<double>(nbPointsSolution));

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
    const double end = MPI_Wtime();

    double seconds = end - start;

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

    MPI_Win_free(&green_cache_win);
    MPI_Finalize();

    return 0;
}