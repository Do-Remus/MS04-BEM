// Algorithme COCG pour resoudre un system pour A sym complexe de manière itérative sans caculer stocker A
#ifndef BEM_ACA_COCG_HPP_INCLUDED
#define BEM_ACA_COCG_HPP_INCLUDED

#include "../../global/commun.hpp"
#include "../../maths/commun.hpp"

#include "operator.hpp"
#include "../operator.hpp"

Vecteur COCG(const Maillage &maillage, const Vecteur &b, const std::vector<QuadratureSegment> &quadrature_maillage, const ArbreBlocs &arbre_blocs, Real tol, Index maxIter, bool utiliser_preconditionnement);

#endif