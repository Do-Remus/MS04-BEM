#include "headers/segment.hpp"

/* Fonctions de la classe Segment */

Segment::Segment(const Point &A, const Point &B)
{
    this->P1 = A;
    this->P2 = B;
}

double Segment::norm() const
{
    Point C = this->P1 - this->P2;
    return C.norm();
}

Point Segment::normale() const
{
    Point normale(this->P2.y - this->P1.y, this->P1.x - this->P2.x);
    return normale;
}

Point Segment::milieu() const
{
    Point milieu((this->P2.x + this->P1.x) / 2, (this->P2.y + this->P1.y) / 2);
    return milieu;
}
/* Fonction associées à la classe Segment */

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
