#include "../green.hpp"

// Cache de green (défaults)
double lambda = 2.0 * pi / k;
Complex *green_cache = nullptr;
MPI_Win green_cache_win = MPI_WIN_NULL;
Real green_cache_step_fine = 0.001;
Real green_cache_step_fine_inv = 1 / green_cache_step_fine;
Real green_cache_step_large = 0.05 * lambda;
Real green_cache_step_large_inv = 1 / green_cache_step_large;
Real green_cache_cutoff = 0.5 * lambda;
std::size_t green_cache_size = 0;
std::size_t green_cache_n_fine = 0;

void initialiser_green_cache()
{
    // -- Paramètres du cache --

    const double lambda = 2.0 * pi / k;
    const double distance_max = L * std::sqrt(2.0);
    double step_old = 0.1 * pasMaillage / k;
    std::size_t N_old = static_cast<std::size_t>(distance_max / step_old);

    double beta = 0.05; // multiplicateur old précision près de 0
    double alpha = 0.1; // multiplicateur de lambda (zone cache fin)

    green_cache_step_fine = beta * step_old;
    green_cache_step_fine_inv = 1 / green_cache_step_fine;

    green_cache_cutoff = alpha * lambda;

    green_cache_step_large = (distance_max - green_cache_cutoff) / (N_old - green_cache_cutoff / green_cache_step_fine);
    green_cache_step_large_inv = 1 / green_cache_step_large;

    delta = green_cache_step_fine;

    // -- Nombre de cases dans chaque zone --

    green_cache_n_fine = static_cast<std::size_t>(green_cache_cutoff / green_cache_step_fine + 0.5);
    std::size_t n_large = static_cast<std::size_t>((distance_max - green_cache_cutoff) / green_cache_step_large + 0.5);

    green_cache_size = green_cache_n_fine + n_large + 1;

    // -- Mémoire partagée MPI --

    const MPI_Aint taille = static_cast<MPI_Aint>(green_cache_size) * static_cast<MPI_Aint>(sizeof(Complex));

    // Communicateur contenant les processus du même nœud.
    MPI_Comm shm_comm;
    MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &shm_comm);

    int shm_rank;
    MPI_Comm_rank(shm_comm, &shm_rank);

    // Un seul processus réserve physiquement la mémoire.
    // Les autres demandent 0 octet.
    MPI_Aint taille_locale = (shm_rank == 0) ? taille : 0;

    MPI_Win_allocate_shared(taille_locale, sizeof(Complex), MPI_INFO_NULL, shm_comm, &green_cache, &green_cache_win);

    // -- Les autres processus récupèrent le pointeur partagé --

    if (shm_rank != 0)
    {
        MPI_Aint taille_partagee;
        int disp_unit;
        void *ptr = nullptr;

        MPI_Win_shared_query(green_cache_win, 0, &taille_partagee, &disp_unit, &ptr);

        green_cache = static_cast<Complex *>(ptr);
    }

    // -- Initialisation du cache pour le processus 0 --

    if (shm_rank == 0)
    {

        for (std::size_t index = 0; index < green_cache_size; ++index)
        {
            double distance;

            if (index <= green_cache_n_fine)
            {
                // Zone fine
                distance = index * green_cache_step_fine;
            }
            else
            {
                // Zone grossière
                distance = green_cache_cutoff + (index - green_cache_n_fine) * green_cache_step_large;
            }

            if (distance == 0.0)
            {
                green_cache[index] = 0.0;
            }
            else
            {
                green_cache[index] = (I / Real(4.0)) * hankel_n(k * distance, 0);
            }
        }
    }

    // -- Tous les processus attendent que le cache soit initialisé --

    MPI_Barrier(shm_comm);
    MPI_Comm_free(&shm_comm);

    green_cache_simd.cutoff = simd_set1(green_cache_cutoff);
    green_cache_simd.fine = simd_set1(green_cache_step_fine);
    green_cache_simd.fine_inv = simd_set1(green_cache_step_fine_inv);
    green_cache_simd.large = simd_set1(green_cache_step_large);
    green_cache_simd.large_inv = simd_set1(green_cache_step_large_inv);
    green_cache_simd.n_fine = simd_set1(static_cast<Real>(green_cache_n_fine));
}
