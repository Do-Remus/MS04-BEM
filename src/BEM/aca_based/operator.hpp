#ifndef BEM_ACA_OPERATOR_HPP_INCLUDED
#define BEM_ACA_OPERATOR_HPP_INCLUDED

#include "../../global/commun.hpp"
#include "../../maths/commun.hpp"

#include "../commun.hpp"

ACA_decomposition aca_bem(const std::vector<Eigen::Index> selection_i, const std::vector<Eigen::Index> selection_j, const Maillage &maillage, const std::vector<QuadratureSegment> &quadrature_maillage, const Eigen::Index rang_max, const double tolerance);

#endif