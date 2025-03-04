#ifndef SEGMENT_HPP_INCLUDED
#define SEGMENT_HPP_INCLUDED
#include "point.hpp"

class Segment
{
public:
    Point P1;
    Point P2;
    Segment() {};
    Segment(const Point &A, const Point &B);
    double norm();
};

Segment::Segment(const Point &A, const Point &B)
{
    this->P1 = A;
    this->P2 = B;
}

double Segment::norm()
{   Point C=this->P1 - this->P2;
    return C.norm();
}

#endif