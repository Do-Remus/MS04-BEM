#ifndef BEM_ERREURS_HPP_INCLUDED
#define BEM_ERREURS_HPP_INCLUDED

#include "../global/commun.hpp"
#include "../maths/commun.hpp"

struct ErreursSolution
{
    std::vector<Complex> u_analytique;
    std::vector<Real> erreurs;

    Real erreurL2;
    Real erreurMax;
};

ErreursSolution calculer_erreurs_solution(const std::vector<Point> &pointsSolution, const std::vector<Complex> &u_cached);

Real calculer_erreur_relative_p(const Maillage &maillage, const Vecteur &vect_p_cached);

#endif
