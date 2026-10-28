#ifndef FONCTIONS_GREEN_HPP_INCLUDED
#define FONCTIONS_GREEN_HPP_INCLUDED

#include "../../global/commun.hpp"
#include "../../SIMD/commun.hpp"

#include "../geometrie/commun.hpp"

/* === Fonction de Hankel === */

// Hankel premiere espece ordre n
inline Complex hankel_n(const Real x, const int n)
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
inline Complex green(const Point &p1, const Point &p2)
{
    Real r = (p1 - p2).norm();

    return (I / Real(4.0)) * hankel_n(k * r, 0);
}

// Partie reguliere de Green
inline Complex green_reguliere(const Point &P1, const Point &P2)
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

extern Real green_cache_step;
extern Real green_cache_step_inv;

extern Real delta;

extern Index green_cache_size;

// Avec interpolation linéaire (beaucoup mieux)
inline Complex green_cached(const Point &p1, const Point &p2)
{
    const Real dx = p1.x - p2.x;
    const Real dy = p1.y - p2.y;

    const Real distance = std::sqrt(dx * dx + dy * dy);

    const Real position = distance * green_cache_step_inv;

    Index index = static_cast<Index>(position);

    // Protection contre un accès à index + 1 hors du tableau.
    if (index >= green_cache_size - 1)
    {
        index = green_cache_size - 2;
    }

    const Real alpha = position - static_cast<Real>(index);

    return (Real(1.0) - alpha) * green_cache[index] + alpha * green_cache[index + 1];
}

// Version vectorisée
struct GreenCacheSIMD
{
    Simd step;
    Simd step_inv;
};

inline GreenCacheSIMD green_cache_simd;

inline void green_cached_SIMD(const Simd Ax, const Simd Ay, const Simd Bx, const Simd By, Simd &G_re, Simd &G_im)
{
    // Distances entre les points.
    const Simd dx = simd_sub(Ax, Bx);
    const Simd dy = simd_sub(Ay, By);

    const Simd distance = simd_sqrt(simd_add(simd_mul(dx, dx), simd_mul(dy, dy)));

    // Position réelle dans le cache.
    const Simd position = simd_mul(distance, green_cache_simd.step_inv);

    // Indice inférieur : troncature, comme dans la version scalaire.
    const SimdIndex index = simd_double_to_int(position);

    // Conversion de l'indice entier en SIMD réel.
    // Prévoir un helper de conversion si ton abstraction SIMD ne le fournit pas.
    const Simd index_real = simd_index_to_real(index);

    // Coefficient d'interpolation.
    const Simd alpha = simd_sub(position, index_real);
    const Simd one_minus_alpha = simd_sub(SIMD_UN, alpha);

    // Indices des deux valeurs à interpoler.
    const SimdIndex index_next = simd_index_add(index, SIMD_UN_ID);

    const SimdIndex index_re = simd_index_mul(index, SIMD_DEUX_ID);
    const SimdIndex index_im = simd_index_add(index_re, SIMD_UN_ID);

    const SimdIndex index_next_re = simd_index_mul(index_next, SIMD_DEUX_ID);
    const SimdIndex index_next_im = simd_index_add(index_next_re, SIMD_UN_ID);

    const Real *cache = reinterpret_cast<const Real *>(green_cache);

    // Lecture des deux valeurs réelles et imaginaires.
    const Simd re0 = simd_gather(cache, index_re);
    const Simd im0 = simd_gather(cache, index_im);
    const Simd re1 = simd_gather(cache, index_next_re);
    const Simd im1 = simd_gather(cache, index_next_im);

    // Interpolation linéaire.
    G_re = simd_add(simd_mul(one_minus_alpha, re0), simd_mul(alpha, re1));
    G_im = simd_add(simd_mul(one_minus_alpha, im0), simd_mul(alpha, im1));
}
// Partie reguliere de Green cached
inline Complex green_reguliere_cached(const Point &P1, const Point &P2)
{
    Real x = (P2 - P1).norm();

    if (x < green_cache_step)
    {
        return (I / Real(4.0)) - (Real(1.0) / (Real(2.0) * pi)) * (gamma_euler + log(k / 2));
    }
    return green_cached(P1, P2) + (Real(1.0) / (Real(2.0) * pi)) * log(x);
}

#endif