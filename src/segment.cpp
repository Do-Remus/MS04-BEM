#include "headers/segment.hpp"

/* Fonctions de la classe Segment */

Segment::Segment(const Point &A, const Point &B)
{
    this->P1 = A;
    this->P2 = B;
    this->milieu = (A + B) * 0.5;
    this->vecteur_norm = B - A;
    this->norm = vecteur_norm.norm();
}

Point Segment::normale() const
{
    Point normale(this->P2.y - this->P1.y, this->P1.x - this->P2.x);
    return normale;
}

bool operator==(const Segment &AB, const Segment &CD)
{
    if (AB.P1 == CD.P1 && AB.P2 == CD.P2)
    {
        return true;
    }
    if (AB.P1 == CD.P2 && AB.P2 == CD.P1)
    {
        return true;
    }
    return false;
}
