#include "config/constantes.hpp"

/* Définition des constantes globales */

const Real pi = atan(1.) * 4;
complex<Real> I(0., 1.);
const Real gamma_euler = 0.57721566490153286060; // Approximated value

// Cache pour green
Complex *green_cache = nullptr;
MPI_Win green_cache_win = MPI_WIN_NULL;
Real green_cache_step = 0.0001;
Real green_cache_step_inv = 1 / green_cache_step;
Real green_cache_step_inv2 = green_cache_step * green_cache_step;
std::size_t max_index = 200000;
Real delta = green_cache_step;
