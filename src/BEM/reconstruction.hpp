#ifndef BEM_RECONSTRUCTION_HPP_INCLUDED
#define BEM_RECONSTRUCTION_HPP_INCLUDED

#include "../global/commun.hpp"
#include "../maths/commun.hpp"

Vecteur reconstruire_solution(const std::vector<Point> &pointsSolution, const Maillage &maillage, const std::vector<QuadratureSegment> &quadrature_maillage, const Vecteur &vect_p_cached);

Complex coefficient_reconstruction(const Index j, const Index i, const std::vector<Point> &pointsSolution, const std::vector<QuadratureSegment> &quadrature_maillage);

#endif