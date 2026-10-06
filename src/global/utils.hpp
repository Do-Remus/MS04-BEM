#ifndef UTILS_HPP_INCLUDED
#define UTILS_HPP_INCLUDED

#include "types.hpp"

// Calcul d'une erreur relative en gérant les division par 0
Real erreur_relative(const Complex &u, const Complex &reference, Real eps = 1e-14);

// retire la fon des lignes de texte
string trim(const string &str);

#endif