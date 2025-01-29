#ifndef SEGMENT_HPP_INCLUDED
#define SEGMENT_HPP_INCLUDED
#include "point.hpp"

class Segment
{public:
    Point P1;
    Point P2;
    Segment(){};
    Segment(const Point& A, const Point& B);

};

Segment::Segment(const Point& A, const Point& B){
    this->P1=A;
    this->P2=B;
}



#endif