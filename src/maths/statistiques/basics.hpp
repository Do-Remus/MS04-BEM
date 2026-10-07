#ifndef MATHS_STATISTIQUES_BASICS_HPP_INCLUDED
#define MATHS_STATISTIQUES_BASICS_HPP_INCLUDED

#include "../../global/commun.hpp"

// Calcul d'une erreur relative en gérant les division par 0
Real erreur_relative(const Complex &u, const Complex &reference, Real eps = 1e-14);

#endif