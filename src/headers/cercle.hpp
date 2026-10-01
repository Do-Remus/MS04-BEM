#ifndef CERCLE_HPP_INCLUDED
#define CERCLE_HPP_INCLUDED
#include "../config/config.hpp"
#include "../config/external.hpp"
#include "../config/constantes.hpp"
#include "point.hpp"

/* Classe Cercle */

class Cercle
{
public:
    Point centre;
    Real rayon = 1.;
    Cercle() {};
    Cercle(const Real rayon, const Point &A);
};

/* Fonctions associées à la classe Cercle */

#endif