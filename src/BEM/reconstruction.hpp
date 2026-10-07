#ifndef BEM_RECONSTRUCTION_HPP_INCLUDED
#define BEM_RECONSTRUCTION_HPP_INCLUDED

#include "../global/commun.hpp"
#include "../maths/commun.hpp"

std::vector<Complex> reconstruire_solution(const std::vector<Point> &pointsSolution, const Maillage &maillage, const std::vector<QuadratureSegment> &quadrature_maillage, const Vecteur &vect_p_cached);

#endif