#include "../COCG.hpp"

Vecteur preconditionneur_diagonal(const Vecteur &diagonale)
{
    const std::size_t N = diagonale.size();
    Vecteur S(N, 0.0);

    for (std::size_t i = 0; i < N; ++i)
    {
        S[i] = 1.0 / std::sqrt(diagonale[i]);
    }

    return S;
}

Vecteur produit_A_preconditionne(const Vecteur &x, const Vecteur &S, const Vecteur &diagonale, const Maillage &maillage, const std::vector<QuadratureSegment> &quadrature_maillage)
{
    Vecteur Sx = x;

    for (std::size_t i = 0; i < x.size(); ++i)
    {
        Sx[i] *= S[i];
    }

    Vecteur ASx = produit_A_x(Sx, maillage, quadrature_maillage, diagonale);

    Vecteur y = ASx;

    for (std::size_t i = 0; i < x.size(); ++i)
    {
        y[i] *= S[i];
    }

    return y;
}

Vecteur COCG(const Maillage &maillage, const Vecteur &b, const std::vector<QuadratureSegment> &quadrature_maillage, Real tol, unsigned int maxIter)
{
    /* === Initialisation === */

    const std::size_t n = b.size();

    /*
     * Calcul et construction du preconditionneur.
     */
    Vecteur diagonale = diagonale_A(maillage, quadrature_maillage);
    // Vecteur S(n, 1.0); // cas sans preconditionnement
    Vecteur S = preconditionneur_diagonal(diagonale);

    /*
     * Second membre preconditionne :
     *
     *     b_preconditionne = S b
     */
    Vecteur b_preconditionne = b;

    for (std::size_t i = 0; i < n; ++i)
    {
        b_preconditionne[i] *= S[i];
    }

    /*
     * On resout :
     *
     *     S A S y = S b
     *
     * puis :
     *
     *     x = S y
     */
    Vecteur y(n, 0.0);

    Vecteur r = b_preconditionne;
    Vecteur d = r;

    const Real bnorm = b_preconditionne.norm();

    if (bnorm == 0.0)
        return y;

    Real relativeResidual = 1.0;

    Complex rho = r.produitBilineaire(r);

    MPI_Barrier(MPI_COMM_WORLD);
    Real startBoucle = MPI_Wtime();

    /* === Boucle principale === */

    for (unsigned int iter = 0; iter < maxIter; ++iter)
    {
        Vecteur Ad = produit_A_preconditionne(d, S, diagonale, maillage, quadrature_maillage);
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

            return y;
        }

        const Complex alpha = rho / denom;

        for (std::size_t i = 0; i < n; ++i)
        {
            y[i] += alpha * d[i];
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

            Vecteur x(n, 0.0);

            for (std::size_t i = 0; i < n; ++i)
            {
                x[i] = S[i] * y[i];
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

            return y;
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

    Vecteur x(n, 0.0);

    for (std::size_t i = 0; i < n; ++i)
    {
        x[i] = S[i] * y[i];
    }

    return x;
}