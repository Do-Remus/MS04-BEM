#include "../COCG.hpp"

Vecteur gradConjMatrixFree(const Maillage &maillage, const Vecteur &b, const std::vector<QuadratureSegment> &quadrature_maillage, Real tol, unsigned int maxIter)
{
    /* === Initialisation === */

    const std::size_t n = b.size();
    Vecteur x(n, 0.0);

    Vecteur r = b; // car x0 = 0
    Vecteur d = r;

    const Real bnorm = b.norm();

    if (bnorm == 0.0)
        return x;

    Real relativeResidual = 1.0;

    Complex rho = r.produitBilineaire(r);

    MPI_Barrier(MPI_COMM_WORLD);
    Real startBoucle = MPI_Wtime();

    /* === Boucle principale === */

    for (unsigned int iter = 0; iter < maxIter; ++iter)
    {
        Vecteur Ad = produit_A_cached_blocked(maillage, d, quadrature_maillage);
        const Complex denom = d.produitBilineaire(Ad);
        const Real scale = d.norm() * Ad.norm();

        if (scale == 0.0 || std::abs(denom) < 1e-20 * scale)
        {
            if (mpi_rank == 0)
            {
                std::cerr
                    << "  COCG : breakdown, "
                    << "denominateur nul "
                    << "a l'iteration "
                    << iter
                    << ", résidu relatif = "
                    << relativeResidual
                    << std::endl;
            }

            return x;
        }

        const Complex alpha = rho / denom;

        for (std::size_t i = 0; i < n; ++i)
        {
            x[i] += alpha * d[i];
            r[i] -= alpha * Ad[i];
        }

        relativeResidual = r.norm() / bnorm;

        if (relativeResidual < tol)
        {
            if (mpi_rank == 0)
            {
                std::cout
                    << "  COCG converge en "
                    << iter + 1
                    << " iterations, "
                    << "résidu relatif = "
                    << relativeResidual
                    << std::endl;
            }

            return x;
        }

        const Complex rhoNew = r.produitBilineaire(r);

        if (std::abs(rho) < 1e-30)
        {
            if (mpi_rank == 0)
            {
                std::cerr
                    << "  COCG : breakdown de rho "
                    << "a l'iteration "
                    << iter
                    << std::endl;
            }

            return x;
        }

        const Complex beta = rhoNew / rho;

        for (std::size_t i = 0; i < n; ++i)
        {
            d[i] = r[i] + beta * d[i];
        }

        rho = rhoNew;

        MPI_Barrier(MPI_COMM_WORLD);
        const Real endBoucle = MPI_Wtime();
        const Real tempsBoucle = endBoucle - startBoucle;
        startBoucle = endBoucle;

        if (mpi_rank == 0)
        {
            std::cout << "    Iter " << iter + 1 << " en " << tempsBoucle << " s, avec un résidu relatif de " << relativeResidual << endl;
        }
    }

    if (mpi_rank == 0)
    {
        std::cout
            << "  COCG : nombre maximal "
            << "d'iterations atteint."
            << std::endl;
    }

    return x;
}
