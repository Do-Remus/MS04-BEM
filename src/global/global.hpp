#ifndef VAR_GOLBAL_HPP_INCLUDED
#define VAR_GLOBAL_HPP_INCLUDED

#include "external.hpp"
#include "types.hpp"
#include "constantes.hpp"
#include "config.hpp"

// Cache pour green
extern Complex *green_cache;
extern MPI_Win green_cache_win;
extern Real green_cache_step_fine;
extern Real green_cache_step_fine_inv;
extern Real green_cache_step_large;
extern Real green_cache_step_large_inv;
extern Real green_cache_cutoff;
extern Real delta;
extern std::size_t green_cache_size;

#endif