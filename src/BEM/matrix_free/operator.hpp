#ifndef BEM_MF_OPERATOR_HPP_INCLUDED
#define BEM_MF_OPERATOR_HPP_INCLUDED

#include "../../global/commun.hpp"
#include "../../maths/commun.hpp"
#include "../../SIMD/commun.hpp"

inline Complex calculer_interaction_scalaire(const std::size_t indice_cible, const std::size_t indice_source, const QuadratureSIMD &quadrature)
{
    const std::size_t nombre_points = quadrature.N;
    const std::size_t nombre_points_quadrature = quadrature.nq;

    Complex coefficient = 0.0;

    for (std::size_t indice_quadrature_cible = 0; indice_quadrature_cible < nombre_points_quadrature; ++indice_quadrature_cible)
    {
        const std::size_t decalage_cible = indice_cible * nombre_points_quadrature + indice_quadrature_cible;

        const Point point_cible{quadrature.target_x[decalage_cible], quadrature.target_y[decalage_cible]};

        const Real poids_cible = quadrature.target_w[decalage_cible];

        for (std::size_t indice_quadrature_source = 0; indice_quadrature_source < nombre_points_quadrature; ++indice_quadrature_source)
        {
            const std::size_t decalage_source = indice_quadrature_source * nombre_points + indice_source;

            const Point point_source{quadrature.source_x[decalage_source], quadrature.source_y[decalage_source]};

            const Complex fonction_green = green_cached(point_cible, point_source);

            const Real facteur = poids_cible * quadrature.source_w[decalage_source];

            coefficient += facteur * fonction_green;
        }
    }

    return coefficient;
}

inline void accumuler_interaction_scalaire(const Complex coefficient, const std::size_t indice_cible, const std::size_t indice_source, const std::size_t debut_bloc_cible, const std::size_t debut_bloc_source, const Vecteur &x, std::vector<Real> &yi_re, std::vector<Real> &yi_im, std::vector<Real> &yj_re, std::vector<Real> &yj_im)
{
    const Real coefficient_re = coefficient.real();
    const Real coefficient_im = coefficient.imag();

    const Real x_source_re = x[indice_source].real();
    const Real x_source_im = x[indice_source].imag();

    yi_re[indice_cible - debut_bloc_cible] += coefficient_re * x_source_re - coefficient_im * x_source_im;
    yi_im[indice_cible - debut_bloc_cible] += coefficient_re * x_source_im + coefficient_im * x_source_re;

    const Real x_cible_re = x[indice_cible].real();
    const Real x_cible_im = x[indice_cible].imag();

    yj_re[indice_source - debut_bloc_source] += coefficient_re * x_cible_re - coefficient_im * x_cible_im;
    yj_im[indice_source - debut_bloc_source] += coefficient_re * x_cible_im + coefficient_im * x_cible_re;
}

template <std::size_t IBLOCK>
inline void calculer_micro_bloc_cible(const std::size_t debut_micro_bloc_cible, const std::size_t debut_bloc_cible, const std::size_t debut_bloc_source, const std::size_t fin_bloc_source, const Vecteur &x, const QuadratureSIMD &quadrature, std::vector<Real> &yi_re, std::vector<Real> &yi_im, std::vector<Real> &yj_re, std::vector<Real> &yj_im)
{
    const std::size_t nombre_points = quadrature.N;
    const std::size_t nombre_points_quadrature = quadrature.nq;

    /*
     * Xi des lignes du micro-bloc.
     */
    Real x_cible_re[IBLOCK];
    Real x_cible_im[IBLOCK];

    Simd x_cible_re_vec[IBLOCK];
    Simd x_cible_im_vec[IBLOCK];

    for (std::size_t indice_ligne = 0; indice_ligne < IBLOCK; ++indice_ligne)
    {
        x_cible_re[indice_ligne] = x[debut_micro_bloc_cible + indice_ligne].real();
        x_cible_im[indice_ligne] = x[debut_micro_bloc_cible + indice_ligne].imag();

        x_cible_re_vec[indice_ligne] = simd_set1(x_cible_re[indice_ligne]);
        x_cible_im_vec[indice_ligne] = simd_set1(x_cible_im[indice_ligne]);
    }

    /*
     * Partie proche de la diagonale.
     */
    const std::size_t debut_sources_simd = std::max(debut_bloc_source, debut_micro_bloc_cible + IBLOCK);

    for (std::size_t indice_cible = debut_micro_bloc_cible; indice_cible < debut_micro_bloc_cible + IBLOCK; ++indice_cible)
    {
        const std::size_t debut_source = std::max(debut_bloc_source, indice_cible + 1);
        const std::size_t fin_sources_scalaires = std::min(debut_sources_simd, fin_bloc_source);

        for (std::size_t indice_source = debut_source; indice_source < fin_sources_scalaires; ++indice_source)
        {
            const Complex coefficient = calculer_interaction_scalaire(indice_cible, indice_source, quadrature);

            accumuler_interaction_scalaire(coefficient, indice_cible, indice_source, debut_bloc_cible, debut_bloc_source, x, yi_re, yi_im, yj_re, yj_im);
        }
    }

    /*
     * Partie SIMD.
     */
    std::size_t indice_source = debut_sources_simd;

    for (; indice_source + SIMD_WIDTH - 1 < fin_bloc_source; indice_source += SIMD_WIDTH)
    {
        Simd coefficient_re_vec[IBLOCK];
        Simd coefficient_im_vec[IBLOCK];

        for (std::size_t indice_ligne = 0; indice_ligne < IBLOCK; ++indice_ligne)
        {
            coefficient_re_vec[indice_ligne] = simd_zero();
            coefficient_im_vec[indice_ligne] = simd_zero();
        }

        Simd y_source_re_acc = simd_load(&yj_re[indice_source - debut_bloc_source]);
        Simd y_source_im_acc = simd_load(&yj_im[indice_source - debut_bloc_source]);

        const Simd x_source_re_vec = simd_load_complex_re(x, indice_source);
        const Simd x_source_im_vec = simd_load_complex_im(x, indice_source);

        for (std::size_t indice_quadrature_cible = 0; indice_quadrature_cible < nombre_points_quadrature; ++indice_quadrature_cible)
        {
            Simd point_cible_x_vec[IBLOCK];
            Simd point_cible_y_vec[IBLOCK];
            Simd poids_cible_vec[IBLOCK];

            for (std::size_t indice_ligne = 0; indice_ligne < IBLOCK; ++indice_ligne)
            {
                const std::size_t decalage_cible = (debut_micro_bloc_cible + indice_ligne) * nombre_points_quadrature + indice_quadrature_cible;

                point_cible_x_vec[indice_ligne] = simd_set1(quadrature.target_x[decalage_cible]);
                point_cible_y_vec[indice_ligne] = simd_set1(quadrature.target_y[decalage_cible]);
                poids_cible_vec[indice_ligne] = simd_set1(quadrature.target_w[decalage_cible]);
            }

            for (std::size_t indice_quadrature_source = 0; indice_quadrature_source < nombre_points_quadrature; ++indice_quadrature_source)
            {
                const std::size_t decalage_source = indice_quadrature_source * nombre_points + indice_source;

                const Simd point_source_x_vec = simd_load(&quadrature.source_x[decalage_source]);
                const Simd point_source_y_vec = simd_load(&quadrature.source_y[decalage_source]);
                const Simd poids_source_vec = simd_load(&quadrature.source_w[decalage_source]);

                Simd fonction_green_re_vec[IBLOCK];
                Simd fonction_green_im_vec[IBLOCK];

                for (std::size_t indice_ligne = 0; indice_ligne < IBLOCK; ++indice_ligne)
                {
                    green_cached_SIMD(point_cible_x_vec[indice_ligne], point_cible_y_vec[indice_ligne], point_source_x_vec, point_source_y_vec, fonction_green_re_vec[indice_ligne], fonction_green_im_vec[indice_ligne]);
                }

                for (std::size_t indice_ligne = 0; indice_ligne < IBLOCK; ++indice_ligne)
                {
                    const Simd facteur_vec = simd_mul(poids_cible_vec[indice_ligne], poids_source_vec);

                    coefficient_re_vec[indice_ligne] = simd_add(coefficient_re_vec[indice_ligne], simd_mul(facteur_vec, fonction_green_re_vec[indice_ligne]));
                    coefficient_im_vec[indice_ligne] = simd_add(coefficient_im_vec[indice_ligne], simd_mul(facteur_vec, fonction_green_im_vec[indice_ligne]));
                }
            }
        }

        /*
         * Contributions Aij * xj et Aji * xi.
         */
        for (std::size_t indice_ligne = 0; indice_ligne < IBLOCK; ++indice_ligne)
        {
            const Simd y_cible_re_vec = simd_sub(simd_mul(coefficient_re_vec[indice_ligne], x_source_re_vec), simd_mul(coefficient_im_vec[indice_ligne], x_source_im_vec));
            const Simd y_cible_im_vec = simd_add(simd_mul(coefficient_re_vec[indice_ligne], x_source_im_vec), simd_mul(coefficient_im_vec[indice_ligne], x_source_re_vec));

            yi_re[debut_micro_bloc_cible + indice_ligne - debut_bloc_cible] += simd_horizontal_sum(y_cible_re_vec);
            yi_im[debut_micro_bloc_cible + indice_ligne - debut_bloc_cible] += simd_horizontal_sum(y_cible_im_vec);

            const Simd y_source_re_vec = simd_sub(simd_mul(coefficient_re_vec[indice_ligne], x_cible_re_vec[indice_ligne]), simd_mul(coefficient_im_vec[indice_ligne], x_cible_im_vec[indice_ligne]));
            const Simd y_source_im_vec = simd_add(simd_mul(coefficient_re_vec[indice_ligne], x_cible_im_vec[indice_ligne]), simd_mul(coefficient_im_vec[indice_ligne], x_cible_re_vec[indice_ligne]));

            y_source_re_acc = simd_add(y_source_re_acc, y_source_re_vec);
            y_source_im_acc = simd_add(y_source_im_acc, y_source_im_vec);
        }

        simd_store(&yj_re[indice_source - debut_bloc_source], y_source_re_acc);
        simd_store(&yj_im[indice_source - debut_bloc_source], y_source_im_acc);
    }

    /*
     * Reste scalaire du bloc SIMD.
     */
    for (; indice_source < fin_bloc_source; ++indice_source)
    {
        for (std::size_t indice_ligne = 0; indice_ligne < IBLOCK; ++indice_ligne)
        {
            const std::size_t indice_cible = debut_micro_bloc_cible + indice_ligne;

            if (indice_source <= indice_cible)
            {
                continue;
            }

            const Complex coefficient = calculer_interaction_scalaire(indice_cible, indice_source, quadrature);

            accumuler_interaction_scalaire(coefficient, indice_cible, indice_source, debut_bloc_cible, debut_bloc_source, x, yi_re, yi_im, yj_re, yj_im);
        }
    }
}

inline void calculer_ligne_cible(const std::size_t indice_cible, const std::size_t debut_bloc_cible, const std::size_t debut_bloc_source, const std::size_t fin_bloc_source, const Vecteur &x, const QuadratureSIMD &quadrature, std::vector<Real> &yi_re, std::vector<Real> &yi_im, std::vector<Real> &yj_re, std::vector<Real> &yj_im)
{
    const std::size_t nombre_points = quadrature.N;
    const std::size_t nombre_points_quadrature = quadrature.nq;

    const std::size_t decalage_cible = indice_cible * nombre_points_quadrature;

    const Real x_cible_re = x[indice_cible].real();
    const Real x_cible_im = x[indice_cible].imag();

    const Simd x_cible_re_vec = simd_set1(x_cible_re);
    const Simd x_cible_im_vec = simd_set1(x_cible_im);

    std::size_t indice_source = std::max(debut_bloc_source, indice_cible + 1);

    /*
     * SIMD sur les sources.
     */
    for (; indice_source + SIMD_WIDTH - 1 < fin_bloc_source; indice_source += SIMD_WIDTH)
    {
        Simd coefficient_re_vec = simd_zero();
        Simd coefficient_im_vec = simd_zero();

        for (std::size_t indice_quadrature_cible = 0; indice_quadrature_cible < nombre_points_quadrature; ++indice_quadrature_cible)
        {
            const Real point_cible_x = quadrature.target_x[decalage_cible + indice_quadrature_cible];
            const Real point_cible_y = quadrature.target_y[decalage_cible + indice_quadrature_cible];
            const Real poids_cible = quadrature.target_w[decalage_cible + indice_quadrature_cible];

            const Simd point_cible_x_vec = simd_set1(point_cible_x);
            const Simd point_cible_y_vec = simd_set1(point_cible_y);
            const Simd poids_cible_vec = simd_set1(poids_cible);

            for (std::size_t indice_quadrature_source = 0; indice_quadrature_source < nombre_points_quadrature; ++indice_quadrature_source)
            {
                const std::size_t decalage_source = indice_quadrature_source * nombre_points + indice_source;

                const Simd point_source_x_vec = simd_load(&quadrature.source_x[decalage_source]);
                const Simd point_source_y_vec = simd_load(&quadrature.source_y[decalage_source]);
                const Simd poids_source_vec = simd_load(&quadrature.source_w[decalage_source]);

                Simd fonction_green_re_vec;
                Simd fonction_green_im_vec;

                green_cached_SIMD(point_cible_x_vec, point_cible_y_vec, point_source_x_vec, point_source_y_vec, fonction_green_re_vec, fonction_green_im_vec);

                const Simd facteur_vec = simd_mul(poids_cible_vec, poids_source_vec);

                coefficient_re_vec = simd_add(coefficient_re_vec, simd_mul(facteur_vec, fonction_green_re_vec));
                coefficient_im_vec = simd_add(coefficient_im_vec, simd_mul(facteur_vec, fonction_green_im_vec));
            }
        }

        const Simd x_source_re_vec = simd_load_complex_re(x, indice_source);
        const Simd x_source_im_vec = simd_load_complex_im(x, indice_source);

        const Simd y_cible_re_vec = simd_sub(simd_mul(coefficient_re_vec, x_source_re_vec), simd_mul(coefficient_im_vec, x_source_im_vec));
        const Simd y_cible_im_vec = simd_add(simd_mul(coefficient_re_vec, x_source_im_vec), simd_mul(coefficient_im_vec, x_source_re_vec));

        yi_re[indice_cible - debut_bloc_cible] += simd_horizontal_sum(y_cible_re_vec);
        yi_im[indice_cible - debut_bloc_cible] += simd_horizontal_sum(y_cible_im_vec);

        const Simd y_source_re_vec = simd_sub(simd_mul(coefficient_re_vec, x_cible_re_vec), simd_mul(coefficient_im_vec, x_cible_im_vec));
        const Simd y_source_im_vec = simd_add(simd_mul(coefficient_re_vec, x_cible_im_vec), simd_mul(coefficient_im_vec, x_cible_re_vec));

        Simd y_source_re_acc = simd_load(&yj_re[indice_source - debut_bloc_source]);
        Simd y_source_im_acc = simd_load(&yj_im[indice_source - debut_bloc_source]);

        y_source_re_acc = simd_add(y_source_re_acc, y_source_re_vec);
        y_source_im_acc = simd_add(y_source_im_acc, y_source_im_vec);

        simd_store(&yj_re[indice_source - debut_bloc_source], y_source_re_acc);
        simd_store(&yj_im[indice_source - debut_bloc_source], y_source_im_acc);
    }

    /*
     * Reste scalaire.
     */
    for (; indice_source < fin_bloc_source; ++indice_source)
    {
        const Complex coefficient = calculer_interaction_scalaire(indice_cible, indice_source, quadrature);

        accumuler_interaction_scalaire(coefficient, indice_cible, indice_source, debut_bloc_cible, debut_bloc_source, x, yi_re, yi_im, yj_re, yj_im);
    }
}

// Calcul de Ax par blocks et parallelisé et vectorisé
Vecteur produit_A_x(const Vecteur &x, const QuadratureSIMD &quadrature, const Vecteur &diagonale, const std::size_t NB_BLOCK = 1);

// Calcul de la diagonal de l'opérateur A
Vecteur diagonale_A(const Maillage &maillage, const std::vector<QuadratureSegment> &quadrature_maillage);

#endif