#ifndef FONCTIONS_GREEN_HPP_INCLUDED
#define FONCTIONS_GREEN_HPP_INCLUDED

#include "../../global/commun.hpp"

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
extern std::size_t green_cache_size;

// Green cached
inline void green_cached_4(
    const Real Ax,
    const Real Ay,
    const __m256d Bx,
    const __m256d By,
    Real &G0_re,
    Real &G0_im,
    Real &G1_re,
    Real &G1_im,
    Real &G2_re,
    Real &G2_im,
    Real &G3_re,
    Real &G3_im)
{
    const __m256d Ax_vec =
        _mm256_set1_pd(Ax);

    const __m256d Ay_vec =
        _mm256_set1_pd(Ay);

    /*
     * Distances aux quatre points sources.
     */
    const __m256d dx =
        _mm256_sub_pd(
            Ax_vec,
            Bx);

    const __m256d dy =
        _mm256_sub_pd(
            Ay_vec,
            By);

    const __m256d distance =
        _mm256_sqrt_pd(
            _mm256_add_pd(
                _mm256_mul_pd(dx, dx),
                _mm256_mul_pd(dy, dy)));

    /*
     * Nombre de cases dans la partie fine du cache.
     */
    const std::size_t n_fine =
        static_cast<std::size_t>(
            green_cache_cutoff /
                green_cache_step_fine +
            0.5);

    /*
     * Indices dans la partie fine.
     *
     * Même formule que green_cached().
     */
    const __m256d index_fine =
        _mm256_add_pd(
            _mm256_div_pd(
                distance,
                _mm256_set1_pd(
                    green_cache_step_fine)),
            _mm256_set1_pd(0.5));

    /*
     * Indices dans la partie grossière.
     *
     * Même formule que green_cached().
     */
    const __m256d index_large =
        _mm256_add_pd(
            _mm256_set1_pd(
                static_cast<Real>(n_fine)),
            _mm256_add_pd(
                _mm256_div_pd(
                    _mm256_sub_pd(
                        distance,
                        _mm256_set1_pd(
                            green_cache_cutoff)),
                    _mm256_set1_pd(
                        green_cache_step_large)),
                _mm256_set1_pd(0.5)));

    /*
     * Choix entre les deux parties du cache.
     */
    const __m256d masque =
        _mm256_cmp_pd(
            distance,
            _mm256_set1_pd(
                green_cache_cutoff),
            _CMP_LT_OQ);

    const __m256d index_selectionne =
        _mm256_blendv_pd(
            index_large,
            index_fine,
            masque);

    /*
     * Conversion des quatre indices en entiers.
     *
     * Les indices sont positifs, donc la troncature
     * correspond au static_cast<std::size_t>() de
     * green_cached().
     */
    const __m128i index =
        _mm256_cvttpd_epi32(
            index_selectionne);

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
    const __m128i deux =
        _mm_set1_epi32(2);

    const __m128i index_re =
        _mm_mullo_epi32(
            index,
            deux);

    const __m128i index_im =
        _mm_add_epi32(
            index_re,
            _mm_set1_epi32(1));

    const Real *cache =
        reinterpret_cast<const Real *>(
            green_cache);

    /*
     * Lecture SIMD des quatre valeurs réelles.
     */
    const __m256d G_re =
        _mm256_i32gather_pd(
            cache,
            index_re,
            8);

    /*
     * Lecture SIMD des quatre valeurs imaginaires.
     */
    const __m256d G_im =
        _mm256_i32gather_pd(
            cache,
            index_im,
            8);

    /*
     * Retour vers les quatre résultats scalaires.
     */
    Real re[4];
    Real im[4];

    _mm256_storeu_pd(
        re,
        G_re);

    _mm256_storeu_pd(
        im,
        G_im);

    G0_re = re[0];
    G0_im = im[0];

    G1_re = re[1];
    G1_im = im[1];

    G2_re = re[2];
    G2_im = im[2];

    G3_re = re[3];
    G3_im = im[3];
}

inline Complex green_cached(const Point &p1, const Point &p2)
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