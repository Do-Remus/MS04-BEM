#include "../operator.hpp"

Vecteur produit_A_x(const Vecteur &x, const Maillage &maillage, const std::vector<QuadratureSegment> &quadrature_maillage, const Vecteur &diagonale, const std::size_t BLOCK)
{
    const std::size_t N = maillage.size();

    Vecteur y_local(N, 0.0);

    // Buffers réutilisés
    std::vector<Complex> yi(BLOCK);
    std::vector<Complex> yj(BLOCK);

    for (std::size_t ib = mpi_rank * BLOCK; ib < N; ib += mpi_size * BLOCK)
    {
        const std::size_t i_end = std::min(ib + BLOCK, N);
        const std::size_t ni = i_end - ib;

        // Remise à zéro du bloc i
        std::fill(yi.begin(), yi.begin() + ni, Complex(0.0));

        // -- DIAGONALE --

        for (std::size_t i = ib; i < i_end; ++i)
        {
            yi[i - ib] += diagonale[i] * x[i];
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
                    const Complex aij = integ_double(green_cached, qS, quadrature_maillage[j]);

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

Vecteur diagonale_A(const Maillage &maillage, const std::vector<QuadratureSegment> &quadrature_maillage)
{
    const std::size_t N = maillage.size();
    Vecteur diagonaleLocale(N, 0.0);

    /*
     * Chaque processus calcule une partie de la diagonale.
     */
    for (std::size_t i = mpi_rank; i < N; i += mpi_size)
    {
        const Segment &S = maillage[i];
        const QuadratureSegment &qS = quadrature_maillage[i];

        Complex aii = integ_double_log(S);
        aii += integ_double(green_reguliere_cached, qS, qS);

        diagonaleLocale[i] = aii;
    }

    /*
     * Chaque processus récupère la diagonale complète.
     */
    Vecteur diagonale(N, 0.0);

    MPI_Allreduce(diagonaleLocale.data(), diagonale.data(), static_cast<int>(N), mpi_complex_type(), MPI_SUM, MPI_COMM_WORLD);

    return diagonale;
}