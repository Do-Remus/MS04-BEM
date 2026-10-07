#ifndef BEM_MF_OPERATOR_HPP_INCLUDED
#define BEM_MF_OPERATOR_HPP_INCLUDED

#include "../../global/commun.hpp"
#include "../../maths/commun.hpp"

// Calcul de Ax par blocks et parallelisé
Vecteur produit_A_x(const Vecteur &x, const Maillage &maillage, const std::vector<QuadratureSegment> &quadrature_maillage, const Vecteur &diagonale, const std::size_t BLOCK = 8);

Vecteur diagonale_A(const Maillage &maillage, const std::vector<QuadratureSegment> &quadrature_maillage);

#endif