#include "headers/utils.hpp"

void initialiser_green_cache()
{
    green_cache_step = 0.1 * min(pasSolution, pasMaillage) / k;
    green_cache_step_inv = 1.0 / green_cache_step;
    green_cache_step_inv2 = green_cache_step_inv * green_cache_step_inv;
    delta = green_cache_step;
    max_index = static_cast<unsigned int>(std::ceil(L * std::sqrt(2.0) / green_cache_step));

    const MPI_Aint n = static_cast<MPI_Aint>(max_index + 1);

    const MPI_Aint taille = n * static_cast<MPI_Aint>(sizeof(Complex));

    // Communicateur contenant les processus du même nœud.
    MPI_Comm shm_comm;
    MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &shm_comm);

    int shm_rank;
    MPI_Comm_rank(shm_comm, &shm_rank);

    // Un seul processus réserve physiquement la mémoire.
    // Les autres demandent 0 octet.
    MPI_Aint taille_locale = (shm_rank == 0) ? taille : 0;

    MPI_Win_allocate_shared(taille_locale, sizeof(Complex), MPI_INFO_NULL, shm_comm, &green_cache, &green_cache_win);

    // Récupération du pointeur vers la mémoire du processus 0.
    if (shm_rank != 0)
    {
        MPI_Aint taille_partagee;
        int disp_unit;
        void *ptr = nullptr;

        MPI_Win_shared_query(green_cache_win, 0, &taille_partagee, &disp_unit, &ptr);

        green_cache = static_cast<Complex *>(ptr);
    }

    // Le processus 0 initialise le cache.
    if (shm_rank == 0)
    {
        green_cache[0] = 0.0;

        for (std::size_t index = 1; index <= max_index; ++index)
        {
            const Real distance = index * green_cache_step;

            green_cache[index] = (I / Real(4.0)) * hankel_n(k * distance, 0);
        }
    }

    // Garantit que tous les processus voient
    // un cache complètement initialisé.
    MPI_Barrier(shm_comm);

    MPI_Comm_free(&shm_comm);
}

Complex u_N_plus_analytique(const Point &P1, Real radius, int N)
{
    Real theta = P1.theta();
    Real x = k * P1.norm();
    Real ka = k * radius;
    Complex iterative_i = 1;
    Complex partial_sum = -(boost::math::cyl_bessel_j(0, ka) / hankel_n(ka, 0)) * hankel_n(x, 0);

    for (int n = 1; n <= N; n++)
    {
        iterative_i *= -I;
        partial_sum -= Real(2.0) * iterative_i * (boost::math::cyl_bessel_j(n, ka) / hankel_n(ka, n)) * hankel_n(x, n) * cos(n * theta); // positive part of the sum
    }

    return partial_sum;
}

Complex u_inc(const Point &P1)
{
    Real theta = P1.theta();
    return exp(-I * k * P1.norm() * cos(theta));
}

Complex q_analytique(const Point &P1, int N)
{
    Real theta = P1.theta();
    Real ka = k * P1.norm();
    Complex iterative_i = 1;
    Complex partial_sum = k * boost::math::cyl_bessel_j(0, ka) * hankel_n(ka, 1) / hankel_n(ka, 0);

    for (int n = 1; n <= N; n++)
    {
        iterative_i *= -I;
        partial_sum -= k * iterative_i * (boost::math::cyl_bessel_j(n, ka) * ((hankel_n(ka, n - 1) - hankel_n(ka, n + 1))) / hankel_n(ka, n)) * cos(n * theta); // positive and negative part of the sum
    }

    return partial_sum;
}

Complex p_analytique(const Point &P1, int N)
{
    Real theta = P1.theta();
    Real ka = k * P1.norm();
    Complex iterative_i = 1;
    Complex partial_sum = k * boost::math::cyl_bessel_j(1, ka) - k * boost::math::cyl_bessel_j(0, ka) * hankel_n(ka, 1) / hankel_n(ka, 0);

    for (int n = 1; n <= N; n++)
    {
        iterative_i *= -I;
        partial_sum += k * iterative_i * cos(n * theta) * (boost::math::cyl_bessel_j(n, ka) * (hankel_n(ka, n - 1) - hankel_n(ka, n + 1)) / hankel_n(ka, n) - (boost::math::cyl_bessel_j(n - 1, ka) - boost::math::cyl_bessel_j(n + 1, ka)));
    }

    return partial_sum;
}

void export_obsctacles(const string &filename, vector<Cercle> cercles)
{
    ofstream f(filename);

    if (!f.is_open())
    {
        cout << "ERROR: Le fichier " << filename << " n'a pas pu être ouvert" << endl;
        exit(-1);
    }

    for (unsigned int i = 0; i < cercles.size(); i++)
    {
        const Cercle &s = cercles[i];
        f << s.centre.x << " " << s.centre.y << " ";
        f << s.rayon << endl;
    }
    f.close();

    return;
}

Real erreur_relative(const Complex &u, const Complex &reference, Real eps)
{
    const Real denom = std::max(std::abs(reference), eps);
    return std::abs(u - reference) / denom;
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

        Complex aii =
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
            const Complex aij = integ_double([](const Point &Q, const Point &P)
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

            Complex aii = integ_simple([&S](const Point &Q)
                                       { return integ_simple_log(Q, S); },
                                       qS);

            aii += integ_double([](const Point &Q, const Point &P)
                                { return green_reguliere_cached_vec(Q, P); },
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
                    // const Complex aij = integ_double([](const Point &Q,
                    //                                     const Point &P)
                    //                                  { return green_cached_vec(Q, P); },
                    //                                  qS,
                    //                                  quadrature_maillage[j]);

                    const Complex aij = integ_double_green_vec(qS, quadrature_maillage[j]);

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

Vecteur gradConjMatrixFree(const Maillage &maillage, const Vecteur &b, const std::vector<QuadratureSegment> &quadrature_maillage, Real tol, unsigned int maxIter)
{
    const std::size_t n = b.size();
    Vecteur x(n, 0.0);

    Vecteur r = b; // car x0 = 0
    Vecteur d = r;

    const Real bnorm = b.norm();

    if (bnorm == 0.0)
        return x;

    Real relativeResidual = 1.0;

    Complex rho = r.produitBilineaire(r);

    for (unsigned int iter = 0; iter < maxIter; ++iter)
    {
        Vecteur Ad = produit_A_cached_blocked(maillage, d, quadrature_maillage);
        const Complex denom = d.produitBilineaire(Ad);
        const Real scale = d.norm() * Ad.norm();

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

        const Complex alpha = rho / denom;

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

        const Complex rhoNew = r.produitBilineaire(r);

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

        const Complex beta = rhoNew / rho;

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
