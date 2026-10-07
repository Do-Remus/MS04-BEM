#ifndef FONCTIONS_GREEN_HPP_INCLUDED
#define FONCTIONS_GREEN_HPP_INCLUDED

#include "../../global/commun.hpp"

#include "../geometrie/commun.hpp"

/* === Fonction de Hankel === */

// Hankel premiere espece ordre n
inline complex<Real> hankel_n(const Real x, const int n)
{
    const Real J = boost::math::cyl_bessel_j(n, x);
    const Real Y = boost::math::cyl_neumann(n, x);

    return {J, Y};
}

/* === Fonction de green === */

/* --- Version normale --- */

// Initialisateur du cache
void initialiser_green_cache();

// Green
inline complex<Real> green(const Point &p1, const Point &p2)
{
    Real r = (p1 - p2).norm();

    return (I / Real(4.0)) * hankel_n(k * r, 0);
}

// Partie reguliere de Green
inline complex<Real> green_reguliere(const Point &P1, const Point &P2)
{
    Real x = (P2 - P1).norm();

    if (x < PRECISION_ZERO_DOUBLE)
    {
        return (I / Real(4.0)) - (Real(1.0) / (Real(2.0) * pi)) * (gamma_euler + log(k / 2));
    }
    return green(P1, P2) + (Real(1.0) / (Real(2.0) * pi)) * log(x);
}

/* --- Version cached --- */

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

// Green cached
inline complex<Real> green_cached(const Point &p1, const Point &p2)
{
    const Real dx = p1.x - p2.x;
    const Real dy = p1.y - p2.y;

    const Real distance = std::sqrt(dx * dx + dy * dy);

    std::size_t index;

    const std::size_t n_fine = static_cast<std::size_t>(green_cache_cutoff / green_cache_step_fine + 0.5);

    if (distance < green_cache_cutoff)
    {
        index = static_cast<std::size_t>(distance / green_cache_step_fine + 0.5);
    }
    else
    {
        index = n_fine + static_cast<std::size_t>((distance - green_cache_cutoff) / green_cache_step_large + 0.5);
    }

    return green_cache[index];
}

// Partie reguliere de Green cached
inline complex<Real> green_reguliere_cached(const Point &P1, const Point &P2)
{
    Real x = (P2 - P1).norm();

    if (x < green_cache_step_fine)
    {
        return (I / Real(4.0)) - (Real(1.0) / (Real(2.0) * pi)) * (gamma_euler + log(k / 2));
    }
    return green_cached(P1, P2) + (Real(1.0) / (Real(2.0) * pi)) * log(x);
}

#endif