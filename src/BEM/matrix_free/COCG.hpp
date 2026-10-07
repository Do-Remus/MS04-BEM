#ifndef BEM_MF_COCG_HPP_INCLUDED
#define BEM_MF_COCG_HPP_INCLUDED

#include "../../global/commun.hpp"
#include "../../maths/commun.hpp"

#include "operator.hpp"

Vecteur preconditionneur_diagonal(const Vecteur &diagonale);

Vecteur produit_A_preconditionne(const Vecteur &x, const Vecteur &S, const Vecteur &diagonale, const Maillage &maillage, const std::vector<QuadratureSegment> &quadrature_maillage);

// Algorithme COCG pour resoudre un system pour A sym complexe de manière itérative sans caculer stocker A
Vecteur COCG(const Maillage &maillage, const Vecteur &b, const std::vector<QuadratureSegment> &quadrature_maillage, Real tol, unsigned int maxIter);

#endif