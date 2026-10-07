#ifndef SEGMENT_HPP_INCLUDED
#define SEGMENT_HPP_INCLUDED

#include "../../global/commun.hpp"

#include "point.hpp"

/* Classe Segment */

class Segment
{
public:
    Point P1;
    Point P2;
    Point milieu;
    Point vecteur_norm;
    Real norm;
    Segment() {};
    Segment(const Point &A, const Point &B);
    Point normale() const;
};

/* Fonctions associées à la classe Segment */

bool operator==(const Segment &AB, const Segment &CD);

#endif