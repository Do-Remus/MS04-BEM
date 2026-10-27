#ifndef BEM_ACA_DECOMPOSITION_HPP_INCLUDED
#define BEM_ACA_DECOMPOSITION_HPP_INCLUDED

#include "../../global/commun.hpp"
#include "../../maths/commun.hpp"

#include "../operator.hpp"

ACA_decomposition aca_bem(const std::vector<Index> selection_i, const std::vector<Index> selection_j, const Maillage &maillage, const std::vector<QuadratureSegment> &quadrature_maillage, const Index rang_max, const double tolerance);

#endif