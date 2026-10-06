#include "../operator.hpp"

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

        Complex aii = integ_double_log(S);

        aii += integ_double([](const Point &Q, const Point &P)
                            { return green_reguliere_cached(Q, P); },
                            qS,
                            qS);

        y_local[i] += aii * x[i];

        for (unsigned int j = i + 1; j < N; ++j)
        {
            const Complex aij = integ_double([](const Point &Q, const Point &P)
                                             { return green_cached(Q, P); },
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

    MPI_Allreduce(y_local.data(), y.data(), static_cast<int>(N), mpi_complex_type(), MPI_SUM, MPI_COMM_WORLD);

    return y;
}

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
    std::vector<Complex> yi(BLOCK);
    std::vector<Complex> yj(BLOCK);

    for (std::size_t ib = rank * BLOCK; ib < N; ib += size * BLOCK)
    {
        const std::size_t i_end = std::min(ib + BLOCK, N);
        const std::size_t ni = i_end - ib;

        // Remise à zéro du bloc i
        std::fill(yi.begin(), yi.begin() + ni, Complex(0.0));

        // -- DIAGONALE --

        for (std::size_t i = ib; i < i_end; ++i)
        {
            const Segment &S = maillage[i];
            const QuadratureSegment &qS = quadrature_maillage[i];

            Complex aii = integ_double_log(S);

            aii += integ_double([](const Point &Q, const Point &P)
                                { return green_reguliere_cached(Q, P); },
                                qS,
                                qS);

            yi[i - ib] += aii * x[i];
        }

        // -- HORS-DIAGONALE : i < j --

        for (std::size_t jb = ib; jb < N; jb += BLOCK)
        {
            const std::size_t j_end = std::min(jb + BLOCK, N);
            const std::size_t nj = j_end - jb;

            // Remise à zéro du bloc j
            std::fill(yj.begin(), yj.begin() + nj, Complex(0.0));

            for (std::size_t i = ib; i < i_end; ++i)
            {
                const QuadratureSegment &qS = quadrature_maillage[i];

                const std::size_t j_start = std::max(jb, i + 1);

                for (std::size_t j = j_start; j < j_end; ++j)
                {
                    const Complex aij = integ_double([](const Point &Q, const Point &P)
                                                     { return green_cached(Q, P); },
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

    MPI_Allreduce(y_local.data(), y.data(), static_cast<int>(N), mpi_complex_type(), MPI_SUM, MPI_COMM_WORLD);

    return y;
}
