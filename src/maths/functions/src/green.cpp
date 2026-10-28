#include "../green.hpp"

// Cache de green (défaults)
Complex *green_cache = nullptr;
MPI_Win green_cache_win = MPI_WIN_NULL;

Real green_cache_step = 0.0;
Real green_cache_step_inv = 0.0;

std::size_t green_cache_size = 0;

void initialiser_green_cache()
{
    // -- Paramètres du cache --

    const Real distance_max = L * std::sqrt(Real(2.0));

    if (!(facteur_step_green > Real(0.0)) || !std::isfinite(facteur_step_green))
    {
        throw std::invalid_argument("facteur_step_green doit être strictement positif et fini.");
    }

    if (!(pasMaillage > Real(0.0)) || !std::isfinite(pasMaillage))
    {
        throw std::invalid_argument("pasMaillage doit être strictement positif et fini.");
    }

    green_cache_step = facteur_step_green * pasMaillage;
    green_cache_step_inv = Real(1.0) / green_cache_step;
    delta = green_cache_step;

    // Deux points supplémentaires pour sécuriser l'interpolation.
    green_cache_size = static_cast<Index>(std::ceil(distance_max * green_cache_step_inv)) + 2;

    // -- Allocation de la mémoire partagée MPI --

    const MPI_Aint taille = static_cast<MPI_Aint>(green_cache_size) * static_cast<MPI_Aint>(sizeof(Complex));

    MPI_Comm shm_comm;
    MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &shm_comm);

    int shm_rank;
    MPI_Comm_rank(shm_comm, &shm_rank);

    const MPI_Aint taille_locale = (shm_rank == 0) ? taille : 0;

    MPI_Win_allocate_shared(taille_locale, sizeof(Complex), MPI_INFO_NULL, shm_comm, &green_cache, &green_cache_win);

    // -- Récupération du pointeur partagé par les autres processus --

    if (shm_rank != 0)
    {
        MPI_Aint taille_partagee;
        int disp_unit;
        void *ptr = nullptr;

        MPI_Win_shared_query(green_cache_win, 0, &taille_partagee, &disp_unit, &ptr);

        green_cache = static_cast<Complex *>(ptr);
    }

    // -- Initialisation par le processus 0 --

    if (shm_rank == 0)
    {
        for (Index index = 0; index < green_cache_size; ++index)
        {
            const Real distance = static_cast<Real>(index) * green_cache_step;

            if (distance == Real(0.0))
            {
                // Valeur artificielle : le noyau Green est singulier en r = 0.
                green_cache[index] = Complex(0.0, 0.0);
            }
            else
            {
                green_cache[index] = (I / Real(4.0)) * hankel_n(k * distance, 0);
            }
        }
    }

    // Synchronisation entre les processus partageant la mémoire.
    MPI_Barrier(shm_comm);
    MPI_Comm_free(&shm_comm);

    // -- Paramètres SIMD --

    green_cache_simd.step = simd_set1(green_cache_step);
    green_cache_simd.step_inv = simd_set1(green_cache_step_inv);
}