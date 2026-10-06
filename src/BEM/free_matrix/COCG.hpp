#ifndef BEM_FM_COCG_HPP_INCLUDED
#define BEM_FM_COCG_HPP_INCLUDED

#include "../../global/global.hpp"
#include "../../geometrie/maillage.hpp"
#include "../../algebre/vecteur.hpp"
#include "operator.hpp"

// Algorithme COCG pour resoudre un system pour A sym complexe de manière itérative sans caculer stocker A
Vecteur gradConjMatrixFree(const Maillage &maillage, const Vecteur &b, const std::vector<QuadratureSegment> &quadrature_maillage, Real tol, unsigned int maxIter);

#endif