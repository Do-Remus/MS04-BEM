#ifndef BEM_ACA_OPERATOR_HPP_INCLUDED
#define BEM_ACA_OPERATOR_HPP_INCLUDED

#include "../../global/commun.hpp"
#include "../../maths/commun.hpp"

#include "../hmat/commun.hpp"

void appliquer_bloc(const Bloc &bloc, const Vecteur &x, Vecteur &y);

Vecteur appliquer_H(const ArbreBlocs &arbre_blocs, const Vecteur &x);

#endif