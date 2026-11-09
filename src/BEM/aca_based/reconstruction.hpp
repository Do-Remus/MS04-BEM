#ifndef BEM_ACA_RECONSTRUCTION_HPP_INCLUDED
#define BEM_ACA_RECONSTRUCTION_HPP_INCLUDED

#include "../../global/commun.hpp"
#include "../../maths/commun.hpp"

#include "../hmat/commun.hpp"
#include "../reconstruction.hpp"
#include "../aca_based/decomposition.hpp"

void remplir_blocs_solution(Bloc &bloc, const std::vector<Point> &pointsSolution, const std::vector<QuadratureSegment> &quadrature_maillage, const Index rang_max, const double tolerance);


//uniquement pour etudier la SVD des blocs denses et admissibles
void construire_bloc_forcage_dense(Bloc &bloc,
                                   const std::vector<Point> &pointsSolution,
                                   const std::vector<QuadratureSegment> &quadrature_maillage,
                                   std::ostream &sortie,
                                   const double seuil_svd);
#endif