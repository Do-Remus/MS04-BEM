#ifndef CONSTANTES_HPP_INCLUDED
#define CONSTANTES_HPP_INCLUDED

#include "external.hpp"

#define PRECISION_ZERO_DOUBLE 1e-10

extern const Real pi;
extern complex<Real> I;
extern const Real gamma_euler;

// Cache pour green
extern Complex *green_cache;
extern MPI_Win green_cache_win;
extern Real green_cache_step;
extern Real green_cache_step_inv;
extern Real green_cache_step_inv2;
extern std::size_t max_index;
extern Real delta;

#endif