#ifndef BEM_FM_OPERATOR_HPP_INCLUDED
#define BEM_FM_OPERATOR_HPP_INCLUDED

#include "../../global/global.hpp"
#include "../../geometrie/maillage.hpp"
#include "../../algebre/vecteur.hpp"
#include "../../integration/integrale.hpp"
#include "../../functions/green.hpp"

// Calcul de Ax parallelisé
Vecteur produit_A_cached(const Maillage &maillage, const Vecteur &x, const std::vector<QuadratureSegment> &quadrature_maillage);

// Calcul de Ax par blocks et parallelisé
Vecteur produit_A_cached_blocked(const Maillage &maillage, const Vecteur &x, const std::vector<QuadratureSegment> &quadrature_maillage);

#endif