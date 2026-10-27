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
    const Real distance_max = L * std::sqrt(Real(2.0));

    // Paramètres initiaux de la grille
    const Real step_old = Real(0.1) * pasMaillage / k;
    const Index N_old = static_cast<Index>(distance_max / step_old);

    const Real beta = Real(0.05);
    const Real alpha = Real(0.1);

    // Pas fin
    green_cache_step_fine = beta * step_old;
    green_cache_step_fine_inv = Real(1.0) / green_cache_step_fine;

    // Nombre d'intervalles fins
    const Real cutoff_initial = alpha * lambda;

    green_cache_n_fine = static_cast<Index>(cutoff_initial / green_cache_step_fine + Real(0.5));

    // Cutoff effectif : exactement sur la grille fine
    green_cache_cutoff = static_cast<Real>(green_cache_n_fine) * green_cache_step_fine;

    // Pas large : les deux grilles partagent exactement le cutoff
    const Index n_large = (N_old > green_cache_n_fine) ? N_old - green_cache_n_fine : Index(1);

    green_cache_step_large = (distance_max - green_cache_cutoff) / static_cast<Real>(n_large);

    green_cache_step_large_inv = Real(1.0) / green_cache_step_large;

    delta = green_cache_step_fine;

    // Taille : indices 0 à n_fine+n_large inclus
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
            Real distance;

            if (index <= green_cache_n_fine)
            {
                // Zone fine
                distance = static_cast<Real>(index) * green_cache_step_fine;
            }
            else
            {
                // Zone grossière
                distance = green_cache_cutoff + static_cast<Real>(index - green_cache_n_fine) * green_cache_step_large;
            }

            if (distance == 0.0)
            {
                green_cache[index] = Complex(0.0, 0.0);
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
