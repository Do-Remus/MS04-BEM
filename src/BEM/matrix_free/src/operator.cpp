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

Vecteur produit_A_x_vec(
    const Vecteur &x,
    const Maillage &maillage,
    const std::vector<QuadratureSegment> &quadrature_maillage,
    const Vecteur &diagonale,
    const std::size_t BLOCK)
{
    const std::size_t N =
        maillage.size();

    const std::size_t nq =
        quadrature_maillage[0].size();

    Vecteur y_local(
        N,
        0.0);

    std::vector<Real> yj_re(
        BLOCK,
        0.0);

    std::vector<Real> yj_im(
        BLOCK,
        0.0);

    for (std::size_t ib = mpi_rank * BLOCK;
         ib < N;
         ib += mpi_size * BLOCK)
    {
        const std::size_t i_end =
            std::min(
                ib + BLOCK,
                N);

        const std::size_t ni =
            i_end - ib;

        std::vector<Real> yi_re(
            ni,
            0.0);

        std::vector<Real> yi_im(
            ni,
            0.0);

        /*
         * Diagonale.
         */
        for (std::size_t i = ib;
             i < i_end;
             ++i)
        {
            const Complex aii_x =
                diagonale[i] * x[i];

            yi_re[i - ib] =
                aii_x.real();

            yi_im[i - ib] =
                aii_x.imag();
        }

        /*
         * Blocs sources.
         */
        for (std::size_t jb = ib;
             jb < N;
             jb += BLOCK)
        {
            const std::size_t j_end =
                std::min(
                    jb + BLOCK,
                    N);

            const std::size_t nj =
                j_end - jb;

            std::fill(
                yj_re.begin(),
                yj_re.begin() + nj,
                0.0);

            std::fill(
                yj_im.begin(),
                yj_im.begin() + nj,
                0.0);

            /*
             * Bloc cible.
             */
            for (std::size_t i = ib;
                 i < i_end;
                 ++i)
            {
                const QuadratureSegment &qI =
                    quadrature_maillage[i];

                const std::size_t j_start =
                    std::max(
                        jb,
                        i + 1);

                std::size_t j =
                    j_start;

                /*
                 * ==================================================
                 * 4 segments sources simultanément.
                 * ==================================================
                 */
                for (; j + 3 < j_end; j += 4)
                {
                    Real aij_re0 = 0.0;
                    Real aij_im0 = 0.0;

                    Real aij_re1 = 0.0;
                    Real aij_im1 = 0.0;

                    Real aij_re2 = 0.0;
                    Real aij_im2 = 0.0;

                    Real aij_re3 = 0.0;
                    Real aij_im3 = 0.0;

                    const QuadratureSegment &qJ0 =
                        quadrature_maillage[j];

                    const QuadratureSegment &qJ1 =
                        quadrature_maillage[j + 1];

                    const QuadratureSegment &qJ2 =
                        quadrature_maillage[j + 2];

                    const QuadratureSegment &qJ3 =
                        quadrature_maillage[j + 3];

                    for (std::size_t a = 0;
                         a < nq;
                         ++a)
                    {
                        const Real Px =
                            qI[a].point.x;

                        const Real Py =
                            qI[a].point.y;

                        const Real wi =
                            qI[a].weight;

                        for (std::size_t b = 0;
                             b < nq;
                             ++b)
                        {
                            Real G_re0;
                            Real G_im0;

                            Real G_re1;
                            Real G_im1;

                            Real G_re2;
                            Real G_im2;

                            Real G_re3;
                            Real G_im3;

                            const __m256d Bx =
                                _mm256_set_pd(
                                    qJ3[b].point.x,
                                    qJ2[b].point.x,
                                    qJ1[b].point.x,
                                    qJ0[b].point.x);

                            const __m256d By =
                                _mm256_set_pd(
                                    qJ3[b].point.y,
                                    qJ2[b].point.y,
                                    qJ1[b].point.y,
                                    qJ0[b].point.y);

                            green_cached_4(
                                Px,
                                Py,
                                Bx,
                                By,

                                G_re0,
                                G_im0,

                                G_re1,
                                G_im1,

                                G_re2,
                                G_im2,

                                G_re3,
                                G_im3);

                            const Real facteur0 =
                                wi *
                                qJ0[b].weight;

                            const Real facteur1 =
                                wi *
                                qJ1[b].weight;

                            const Real facteur2 =
                                wi *
                                qJ2[b].weight;

                            const Real facteur3 =
                                wi *
                                qJ3[b].weight;

                            aij_re0 +=
                                facteur0 *
                                G_re0;

                            aij_im0 +=
                                facteur0 *
                                G_im0;

                            aij_re1 +=
                                facteur1 *
                                G_re1;

                            aij_im1 +=
                                facteur1 *
                                G_im1;

                            aij_re2 +=
                                facteur2 *
                                G_re2;

                            aij_im2 +=
                                facteur2 *
                                G_im2;

                            aij_re3 +=
                                facteur3 *
                                G_re3;

                            aij_im3 +=
                                facteur3 *
                                G_im3;
                        }
                    }

                    /*
                     * xi.
                     */
                    const Real xi_re =
                        x[i].real();

                    const Real xi_im =
                        x[i].imag();

                    /*
                     * Aij * xj
                     */
                    const Real x0_re =
                        x[j].real();

                    const Real x0_im =
                        x[j].imag();

                    yi_re[i - ib] +=
                        aij_re0 * x0_re -
                        aij_im0 * x0_im;

                    yi_im[i - ib] +=
                        aij_re0 * x0_im +
                        aij_im0 * x0_re;

                    const Real x1_re =
                        x[j + 1].real();

                    const Real x1_im =
                        x[j + 1].imag();

                    yi_re[i - ib] +=
                        aij_re1 * x1_re -
                        aij_im1 * x1_im;

                    yi_im[i - ib] +=
                        aij_re1 * x1_im +
                        aij_im1 * x1_re;

                    const Real x2_re =
                        x[j + 2].real();

                    const Real x2_im =
                        x[j + 2].imag();

                    yi_re[i - ib] +=
                        aij_re2 * x2_re -
                        aij_im2 * x2_im;

                    yi_im[i - ib] +=
                        aij_re2 * x2_im +
                        aij_im2 * x2_re;

                    const Real x3_re =
                        x[j + 3].real();

                    const Real x3_im =
                        x[j + 3].imag();

                    yi_re[i - ib] +=
                        aij_re3 * x3_re -
                        aij_im3 * x3_im;

                    yi_im[i - ib] +=
                        aij_re3 * x3_im +
                        aij_im3 * x3_re;

                    /*
                     * Aji * xi.
                     *
                     * Aji = Aij.
                     */
                    yj_re[j - jb] +=
                        aij_re0 * xi_re -
                        aij_im0 * xi_im;

                    yj_im[j - jb] +=
                        aij_re0 * xi_im +
                        aij_im0 * xi_re;

                    yj_re[j + 1 - jb] +=
                        aij_re1 * xi_re -
                        aij_im1 * xi_im;

                    yj_im[j + 1 - jb] +=
                        aij_re1 * xi_im +
                        aij_im1 * xi_re;

                    yj_re[j + 2 - jb] +=
                        aij_re2 * xi_re -
                        aij_im2 * xi_im;

                    yj_im[j + 2 - jb] +=
                        aij_re2 * xi_im +
                        aij_im2 * xi_re;

                    yj_re[j + 3 - jb] +=
                        aij_re3 * xi_re -
                        aij_im3 * xi_im;

                    yj_im[j + 3 - jb] +=
                        aij_re3 * xi_im +
                        aij_im3 * xi_re;
                }

                /*
                 * ==================================================
                 * Reste scalaire.
                 * ==================================================
                 */
                for (; j < j_end; ++j)
                {
                    const QuadratureSegment &qJ =
                        quadrature_maillage[j];

                    Real aij_re =
                        0.0;

                    Real aij_im =
                        0.0;

                    for (std::size_t a = 0;
                         a < nq;
                         ++a)
                    {
                        const Real Px =
                            qI[a].point.x;

                        const Real Py =
                            qI[a].point.y;

                        const Real wi =
                            qI[a].weight;

                        for (std::size_t b = 0;
                             b < nq;
                             ++b)
                        {
                            const Complex G =
                                green_cached(
                                    qI[a].point,
                                    qJ[b].point);

                            const Real facteur =
                                wi *
                                qJ[b].weight;

                            aij_re +=
                                facteur *
                                G.real();

                            aij_im +=
                                facteur *
                                G.imag();
                        }
                    }

                    const Real xj_re =
                        x[j].real();

                    const Real xj_im =
                        x[j].imag();

                    yi_re[i - ib] +=
                        aij_re * xj_re -
                        aij_im * xj_im;

                    yi_im[i - ib] +=
                        aij_re * xj_im +
                        aij_im * xj_re;

                    const Real xi_re =
                        x[i].real();

                    const Real xi_im =
                        x[i].imag();

                    yj_re[j - jb] +=
                        aij_re * xi_re -
                        aij_im * xi_im;

                    yj_im[j - jb] +=
                        aij_re * xi_im +
                        aij_im * xi_re;
                }
            }

            /*
             * Écriture du bloc source.
             */
            for (std::size_t j = jb;
                 j < j_end;
                 ++j)
            {
                y_local[j] +=
                    Complex(
                        yj_re[j - jb],
                        yj_im[j - jb]);
            }
        }

        /*
         * Écriture du bloc cible.
         */
        for (std::size_t i = ib;
             i < i_end;
             ++i)
        {
            y_local[i] +=
                Complex(
                    yi_re[i - ib],
                    yi_im[i - ib]);
        }
    }

    Vecteur y(
        N,
        0.0);

    MPI_Allreduce(
        y_local.data(),
        y.data(),
        static_cast<int>(N),
        mpi_complex_type(),
        MPI_SUM,
        MPI_COMM_WORLD);

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