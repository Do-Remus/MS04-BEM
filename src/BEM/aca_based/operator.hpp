#ifndef BEM_ACA_OPERATOR_HPP_INCLUDED
#define BEM_ACA_OPERATOR_HPP_INCLUDED

#include "../../global/commun.hpp"
#include "../../maths/commun.hpp"

#include "../hmat/commun.hpp"
#include "../operator.hpp"
#include "../aca_based/decomposition.hpp"

void remplir_blocs_bem(Bloc &bloc, const Maillage &maillage, const std::vector<QuadratureSegment> &quadrature_maillage, const Index rang_max, const double tolerance);

#endif