#ifndef BEM_OPERATOR_HPP_INCLUDED
#define BEM_OPERATOR_HPP_INCLUDED

#include "../global/commun.hpp"
#include "../maths/commun.hpp"

Matrice get_matrice_bem(const Maillage &maillage, const std::vector<QuadratureSegment> &quadrature_maillage);

Vecteur ligne_bem(const Eigen::Index i, const Maillage &maillage, const std::vector<QuadratureSegment> &quadrature_maillage);

Vecteur colonne_bem(const Eigen::Index j, const Maillage &maillage, const std::vector<QuadratureSegment> &quadrature_maillage);

inline Complex matrice_bem(const Eigen::Index i, const Eigen::Index j, const Maillage &maillage, const std::vector<QuadratureSegment> &quadrature_maillage)
{
    const Segment &Si = maillage[i];
    const QuadratureSegment &qi = quadrature_maillage[i];

    /*
     * Diagonale régularisée.
     */
    if (i == j)
        return integ_double_log(Si) + integ_double(green_reguliere, qi, qi);

    /*
     * Hors-diagonale : calcul d'un seul triangle.
     */
    return integ_double(green, qi, quadrature_maillage[j]);
}

ACA_decomposition aca_bem(const Maillage &maillage, const std::vector<QuadratureSegment> &quadrature_maillage, const Eigen::Index rang_max, const double tolerance);

#endif