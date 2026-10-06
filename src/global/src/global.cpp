#include "../global.hpp"

// Cache de green
double lambda = 2.0 * pi / k;
Complex *green_cache = nullptr;
MPI_Win green_cache_win = MPI_WIN_NULL;
Real green_cache_step_fine = 0.001;
Real green_cache_step_fine_inv = 1 / green_cache_step_fine;
Real green_cache_step_large = 0.05 * lambda;
Real green_cache_step_large_inv = 1 / green_cache_step_large;
Real green_cache_cutoff = 0.5 * lambda;
Real delta = green_cache_step_fine;
std::size_t green_cache_size = 0;
