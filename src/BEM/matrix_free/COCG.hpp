#ifndef BEM_MF_COCG_HPP_INCLUDED
#define BEM_MF_COCG_HPP_INCLUDED

#include "../../global/commun.hpp"
#include "../../maths/commun.hpp"

#include "operator.hpp"

// Algorithme COCG pour resoudre un system pour A sym complexe de manière itérative sans caculer stocker A
Vecteur gradConjMatrixFree(const Maillage &maillage, const Vecteur &b, const std::vector<QuadratureSegment> &quadrature_maillage, Real tol, unsigned int maxIter);

#endif