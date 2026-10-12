#include "../operator.hpp"

Vecteur produit_A_x(const Vecteur &x, const QuadratureSIMD &quadrature, const Vecteur &diagonale, const Index NB_BLOCK)
{
    const Index BLOCK = NB_BLOCK * SIMD_WIDTH;

    constexpr Index IBLOCK = 4;

    const Index N = quadrature.N;

    Vecteur y_local = Vecteur::Zero(N);
    std::vector<Real> y_source_re(BLOCK, 0.0);
    std::vector<Real> y_source_im(BLOCK, 0.0);

    /*
     * Blocs cibles.
     */
    for (Index debut_bloc_cible = mpi_rank * BLOCK; debut_bloc_cible < N; debut_bloc_cible += mpi_size * BLOCK)
    {
        const Index fin_bloc_cible = std::min(debut_bloc_cible + BLOCK, N);
        const Index nombre_lignes_cibles = fin_bloc_cible - debut_bloc_cible;

        std::vector<Real> y_cible_re(nombre_lignes_cibles, 0.0);
        std::vector<Real> y_cible_im(nombre_lignes_cibles, 0.0);

        /*
         * Diagonale.
         */
        for (Index indice_cible = debut_bloc_cible; indice_cible < fin_bloc_cible; ++indice_cible)
        {
            const Complex contribution_diagonale = diagonale[indice_cible] * x[indice_cible];

            y_cible_re[indice_cible - debut_bloc_cible] = contribution_diagonale.real();
            y_cible_im[indice_cible - debut_bloc_cible] = contribution_diagonale.imag();
        }

        /*
         * Blocs sources.
         */
        for (Index debut_bloc_source = debut_bloc_cible; debut_bloc_source < N; debut_bloc_source += BLOCK)
        {
            const Index fin_bloc_source = std::min(debut_bloc_source + BLOCK, N);
            const Index nombre_sources = fin_bloc_source - debut_bloc_source;

            std::fill(y_source_re.begin(), y_source_re.begin() + nombre_sources, 0.0);
            std::fill(y_source_im.begin(), y_source_im.begin() + nombre_sources, 0.0);

            /*
             * Micro-blocs cibles complets.
             */
            Index debut_micro_bloc_cible = debut_bloc_cible;

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
            for (Index indice_source = debut_bloc_source; indice_source < fin_bloc_source; ++indice_source)
            {
                y_local[indice_source] += Complex(y_source_re[indice_source - debut_bloc_source], y_source_im[indice_source - debut_bloc_source]);
            }
        }

        /*
         * Écriture du bloc cible.
         */
        for (Index indice_cible = debut_bloc_cible; indice_cible < fin_bloc_cible; ++indice_cible)
        {
            y_local[indice_cible] += Complex(y_cible_re[indice_cible - debut_bloc_cible], y_cible_im[indice_cible - debut_bloc_cible]);
        }
    }

    /*
     * Réduction MPI.
     */
    Vecteur y = Vecteur::Zero(N);

    MPI_Allreduce(y_local.data(), y.data(), static_cast<int>(N), mpi_complex_type(), MPI_SUM, MPI_COMM_WORLD);

    return y;
}

Vecteur diagonale_A(const Maillage &maillage, const std::vector<QuadratureSegment> &quadrature_maillage)
{
    const Index N = maillage.size();
    Vecteur diagonaleLocale = Vecteur::Zero(N);

    /*
     * Chaque processus calcule une partie de la diagonale.
     */
    for (Index i = mpi_rank; i < N; i += mpi_size)
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
    Vecteur diagonale = Vecteur::Zero(N);

    MPI_Allreduce(diagonaleLocale.data(), diagonale.data(), static_cast<int>(N), mpi_complex_type(), MPI_SUM, MPI_COMM_WORLD);

    return diagonale;
}
