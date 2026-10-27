#ifndef BEM_ACA_DECOMPOSITION_HPP_INCLUDED
#define BEM_ACA_DECOMPOSITION_HPP_INCLUDED

#include "../../global/commun.hpp"
#include "../../maths/commun.hpp"

#include "../operator.hpp"

ACA_decomposition aca_bem(const std::vector<Index> selection_i, const std::vector<Index> selection_j, const Maillage &maillage, const std::vector<QuadratureSegment> &quadrature_maillage, const Index rang_max, const double tolerance);

struct ComparaisonVecteurs
{
    double erreur_relative;
    double erreur_max;
};

inline ComparaisonVecteurs comparer_vecteurs(
    const Vecteur &direct,
    const Vecteur &cached)
{
    if (direct.size() != cached.size())
        throw std::invalid_argument("Tailles incompatibles");

    const double norme_direct = direct.norm();

    return {
        (cached - direct).norm() /
            std::max(norme_direct, std::numeric_limits<double>::min()),
        (cached - direct).cwiseAbs().maxCoeff()};
}

#endif