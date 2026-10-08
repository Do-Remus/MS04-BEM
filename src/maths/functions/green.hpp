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
extern Real green_cache_step_fine;
extern Real green_cache_step_fine_inv;
extern Real green_cache_step_large;
extern Real green_cache_step_large_inv;
extern Real green_cache_cutoff;
extern Real delta;
extern std::size_t green_cache_n_fine;
extern std::size_t green_cache_size;

// Green cached
inline Complex green_cached(const Point &p1, const Point &p2)
{
    const Real dx = p1.x - p2.x;
    const Real dy = p1.y - p2.y;

    const Real distance = std::sqrt(dx * dx + dy * dy);

    std::size_t index;

    if (distance < green_cache_cutoff)
    {
        index = static_cast<std::size_t>(distance / green_cache_step_fine + 0.5);
    }
    else
    {
        index = green_cache_n_fine + static_cast<std::size_t>((distance - green_cache_cutoff) / green_cache_step_large + 0.5);
    }

    return green_cache[index];
}

// Version vectorisée
struct GreenCacheSIMD
{
    Simd cutoff;
    Simd fine;
    Simd fine_inv;
    Simd large;
    Simd large_inv;
    Simd n_fine;
};

inline GreenCacheSIMD green_cache_simd;

inline void green_cached_SIMD(const Simd Ax, const Simd Ay, const Simd Bx, const Simd By, Simd &G_re, Simd &G_im)
{
    /*
     * Distances aux points sources.
     */
    const Simd dx = simd_sub(Ax, Bx);
    const Simd dy = simd_sub(Ay, By);
    const Simd distance = simd_sqrt(simd_add(simd_mul(dx, dx), simd_mul(dy, dy)));

    /*
     * Indices dans la partie fine.
     *
     * Même formule que green_cached().
     */
    const Simd index_fine = simd_add(simd_mul(distance, green_cache_simd.fine_inv), SIMD_HALF);

    /*
     * Indices dans la partie grossière.
     *
     * Même formule que green_cached().
     */
    const Simd distance_large = simd_sub(distance, green_cache_simd.cutoff);
    const Simd index_large = simd_add(green_cache_simd.n_fine, simd_add(simd_mul(distance_large, green_cache_simd.large_inv), SIMD_HALF));

    /*
     * Choix entre les deux parties du cache.
     */
    const SimdMask masque = simd_cmp_lt(distance, green_cache_simd.cutoff);
    const Simd index_selectionne = simd_blend(index_large, index_fine, masque);

    /*
     * Conversion des indices en entiers.
     *
     * Les indices sont positifs, donc la troncature
     * correspond au static_cast<std::size_t>() de
     * green_cached().
     */
    const SimdIndex index = simd_double_to_int(index_selectionne);

    /*
     * green_cache est un tableau de Complex :
     *
     *   green_cache[0] = {re0, im0}
     *   green_cache[1] = {re1, im1}
     *   ...
     *
     * Vu comme un tableau de double :
     *
     *   cache[2*i]     = real
     *   cache[2*i + 1] = imag
     */
    const SimdIndex index_re = simd_index_mul(index, SIMD_DEUX);
    const SimdIndex index_im = simd_index_add(index_re, SIMD_UN);

    const Real *cache = reinterpret_cast<const Real *>(green_cache);

    /*
     * Lecture SIMD des valeurs réelles et imaginaires.
     */
    G_re = simd_gather(cache, index_re);
    G_im = simd_gather(cache, index_im);
}

// Partie reguliere de Green cached
inline Complex green_reguliere_cached(const Point &P1, const Point &P2)
{
    Real x = (P2 - P1).norm();

    if (x < green_cache_step_fine)
    {
        return (I / Real(4.0)) - (Real(1.0) / (Real(2.0) * pi)) * (gamma_euler + log(k / 2));
    }
    return green_cached(P1, P2) + (Real(1.0) / (Real(2.0) * pi)) * log(x);
}

#endif