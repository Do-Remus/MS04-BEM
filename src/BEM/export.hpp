#ifndef BEM_EXPORT_HPP_INCLUDED
#define BEM_EXPORT_HPP_INCLUDED

#include "../global/commun.hpp"
#include "erreurs.hpp"

void ecrire_solution(const std::vector<Point> &pointsSolution, const std::vector<Complex> &u_cached, const ErreursSolution &erreurs);

#endif