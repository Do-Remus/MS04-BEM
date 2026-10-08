#include "../operator.hpp"

Vecteur produit_A_x(const Vecteur &x, const QuadratureSIMD &quadrature, const Vecteur &diagonale, const std::size_t NB_BLOCK)
{
    const std::size_t BLOCK = NB_BLOCK * SIMD_WIDTH;

    constexpr std::size_t IBLOCK = 4;

    const std::size_t N = quadrature.N;

    Vecteur y_local(N, 0.0);
    std::vector<Real> y_source_re(BLOCK, 0.0);
    std::vector<Real> y_source_im(BLOCK, 0.0);

    /*
     * Blocs cibles.
     */
    for (std::size_t debut_bloc_cible = mpi_rank * BLOCK; debut_bloc_cible < N; debut_bloc_cible += mpi_size * BLOCK)
    {
        const std::size_t fin_bloc_cible = std::min(debut_bloc_cible + BLOCK, N);
        const std::size_t nombre_lignes_cibles = fin_bloc_cible - debut_bloc_cible;

        std::vector<Real> y_cible_re(nombre_lignes_cibles, 0.0);
        std::vector<Real> y_cible_im(nombre_lignes_cibles, 0.0);

        /*
         * Diagonale.
         */
        for (std::size_t indice_cible = debut_bloc_cible; indice_cible < fin_bloc_cible; ++indice_cible)
        {
            const Complex contribution_diagonale = diagonale[indice_cible] * x[indice_cible];

            y_cible_re[indice_cible - debut_bloc_cible] = contribution_diagonale.real();
            y_cible_im[indice_cible - debut_bloc_cible] = contribution_diagonale.imag();
        }

        /*
         * Blocs sources.
         */
        for (std::size_t debut_bloc_source = debut_bloc_cible; debut_bloc_source < N; debut_bloc_source += BLOCK)
        {
            const std::size_t fin_bloc_source = std::min(debut_bloc_source + BLOCK, N);
            const std::size_t nombre_sources = fin_bloc_source - debut_bloc_source;

            std::fill(y_source_re.begin(), y_source_re.begin() + nombre_sources, 0.0);
            std::fill(y_source_im.begin(), y_source_im.begin() + nombre_sources, 0.0);

            /*
             * Micro-blocs cibles complets.
             */
            std::size_t debut_micro_bloc_cible = debut_bloc_cible;

            for (; debut_micro_bloc_cible + IBLOCK <= fin_bloc_cible; debut_micro_bloc_cible += IBLOCK)
            {
                calculer_micro_bloc_cible<IBLOCK>(debut_micro_bloc_cible, debut_bloc_cible, debut_bloc_source, fin_bloc_source, x, quadrature, y_cible_re, y_cible_im, y_source_re, y_source_im);
            }

            /*
             * Lignes cibles restantes.
             */
            for (; debut_micro_bloc_cible < fin_bloc_cible; ++debut_micro_bloc_cible)
            {
                calculer_ligne_cible(debut_micro_bloc_cible, debut_bloc_cible, debut_bloc_source, fin_bloc_source, x, quadrature, y_cible_re, y_cible_im, y_source_re, y_source_im);
            }

            /*
             * Écriture du bloc source.
             */
            for (std::size_t indice_source = debut_bloc_source; indice_source < fin_bloc_source; ++indice_source)
            {
                y_local[indice_source] += Complex(y_source_re[indice_source - debut_bloc_source], y_source_im[indice_source - debut_bloc_source]);
            }
        }

        /*
         * Écriture du bloc cible.
         */
        for (std::size_t indice_cible = debut_bloc_cible; indice_cible < fin_bloc_cible; ++indice_cible)
        {
            y_local[indice_cible] += Complex(y_cible_re[indice_cible - debut_bloc_cible], y_cible_im[indice_cible - debut_bloc_cible]);
        }
    }

    /*
     * Réduction MPI.
     */
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
