#include "../COCG.hpp"

Vecteur produit_H_preconditionne(const Vecteur &x, const Vecteur &S, const ArbreBlocs &arbre_blocs)
{
    const Vecteur Sx = S.cwiseProduct(x);

    const Vecteur HSx = appliquer_H(arbre_blocs, Sx);

    const Vecteur y = S.cwiseProduct(HSx);

    return y;
}

Vecteur COCG(const Maillage &maillage, const Vecteur &b, const std::vector<QuadratureSegment> &quadrature_maillage, const ArbreBlocs &arbre_blocs, Real tol, Index maxIter, bool utiliser_preconditionnement)
{
    /* === Initialisation === */

    const Index n = b.size();

    /*
     * Calcul et construction du preconditionneur.
     */
    Vecteur diagonale = diagonale_matrice_bem(maillage, quadrature_maillage);
    // Vecteur S(n, 1.0); // cas sans preconditionnement
    /* === Préconditionnement optionnel === */

    Vecteur S = Vecteur::Ones(n);

    if (utiliser_preconditionnement)
    {
        const Vecteur diagonale = diagonale_matrice_bem(maillage, quadrature_maillage);

        S = preconditionneur_diagonal(diagonale);
    }

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
    // Initialisation diagonale
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
        Vecteur Ad = produit_H_preconditionne(d, S, arbre_blocs);
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

            return S.cwiseProduct(y);
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

        if (mpi_rank == 0 && ((iter % 50) == 0))
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
