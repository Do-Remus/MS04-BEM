#ifndef SEGMENT_HPP_INCLUDED
#define SEGMENT_HPP_INCLUDED
#include "../config/config.hpp"
#include "../config/external.hpp"
#include "../config/constantes.hpp"
#include "point.hpp"

/* Classe Segment */

class Segment
{
public:
    Point P1;
    Point P2;
    Segment() {};
    Segment(const Point &A, const Point &B);
    double norm() const;
    Point normale() const;
};

/* Fonctions associées à la classe Segment */

bool operator==(const Segment &AB, const Segment &CD);

#endif