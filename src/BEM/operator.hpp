#ifndef BEM_OPERATOR_HPP_INCLUDED
#define BEM_OPERATOR_HPP_INCLUDED

#include "../global/commun.hpp"
#include "../maths/commun.hpp"

Matrice get_matrice_bem(const Maillage &maillage, const std::vector<QuadratureSegment> &quadrature_maillage);

Vecteur ligne_bem(const Index i, const Maillage &maillage, const std::vector<QuadratureSegment> &quadrature_maillage);

Vecteur colonne_bem(const Index j, const Maillage &maillage, const std::vector<QuadratureSegment> &quadrature_maillage);

Vecteur ligne_bem(const Index i, const std::vector<Index> &selection_j, const Maillage &maillage, const std::vector<QuadratureSegment> &quadrature_maillage);

Vecteur colonne_bem(const Index j, const std::vector<Index> &selection_i, const Maillage &maillage, const std::vector<QuadratureSegment> &quadrature_maillage);

inline Complex matrice_bem(const Index i, const Index j, const Maillage &maillage, const std::vector<QuadratureSegment> &quadrature_maillage)
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

#endif