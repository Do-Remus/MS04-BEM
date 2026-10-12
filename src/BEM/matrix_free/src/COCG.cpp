#include "../COCG.hpp"

Vecteur preconditionneur_diagonal(const Vecteur &diagonale)
{
    const Index N = diagonale.size();
    Vecteur S = Vecteur::Zero(N);

    for (Index i = 0; i < N; ++i)
    {
        S[i] = 1.0 / std::sqrt(diagonale[i]);
    }

    return S;
}

Vecteur produit_A_preconditionne(const Vecteur &x, const Vecteur &S, const Vecteur &diagonale, const QuadratureSIMD &quadrature)
{
    Vecteur Sx = S.cwiseProduct(x);

    Vecteur ASx = produit_A_x(Sx, quadrature, diagonale);

    Vecteur y = S.cwiseProduct(ASx);

    return y;
}

Vecteur COCG(const Maillage &maillage, const Vecteur &b, const std::vector<QuadratureSegment> &quadrature_maillage, const QuadratureSIMD &quadrature_simd, Real tol, Index maxIter)
{
    /* === Initialisation === */

    const Index n = b.size();

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
    Vecteur b_preconditionne = S.cwiseProduct(b);
    /*
     * On resout :
     *
     *     S A S y = S b
     *
     * puis :
     *
     *     x = S y
     */
    Vecteur y = Vecteur::Zero(n);

    Vecteur r = b_preconditionne;
    Vecteur d = r;

    const Real bnorm = b_preconditionne.norm();

    if (bnorm == 0.0)
        return y;

    Real relativeResidual = 1.0;

    Complex rho = produitBilineaire(r, r);

    MPI_Barrier(MPI_COMM_WORLD);
    Real startBoucle = MPI_Wtime();

    /* === Boucle principale === */

    for (Index iter = 0; iter < maxIter; ++iter)
    {
        Vecteur Ad = produit_A_preconditionne(d, S, diagonale, quadrature_simd);
        const Complex denom = produitBilineaire(Ad, d);
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

        y += alpha * d;
        r -= alpha * Ad;

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

            Vecteur x = S.cwiseProduct(y);

            return x;
        }

        const Complex rhoNew = produitBilineaire(r, r);

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

            Vecteur x = S.cwiseProduct(y);

            return x;
        }

        const Complex beta = rhoNew / rho;

        d = r + beta * d;
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
            << "    COCG : nombre maximal "
            << "d'iterations atteint."
            << std::endl;
    }

    Vecteur x = S.cwiseProduct(y);

    return x;
}
