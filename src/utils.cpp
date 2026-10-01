#include "headers/utils.hpp"

Complex *green_cache = nullptr;
MPI_Win green_cache_win = MPI_WIN_NULL;
Real green_cache_step = 0.0001;
Real green_cache_step_inv = 1 / green_cache_step;
Real green_cache_step_inv2 = green_cache_step * green_cache_step;
std::size_t max_index = 200000;
Real delta = green_cache_step;

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
