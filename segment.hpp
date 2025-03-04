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


bool operator==(const Segment AB, const Segment CD){
    if(AB.P1==CD.P1 && AB.P2==CD.P2){
        return true;
    }
    if(AB.P1==CD.P2 && AB.P2==CD.P1){
        return true;
    }
    return false;

}
#endif