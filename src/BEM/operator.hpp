#ifndef BEM_OPERATOR_HPP_INCLUDED
#define BEM_OPERATOR_HPP_INCLUDED

#include "../global/commun.hpp"
#include "../maths/commun.hpp"

Matrice get_matrice_bem(const Maillage &maillage, const std::vector<QuadratureSegment> &quadrature_maillage);

#endif