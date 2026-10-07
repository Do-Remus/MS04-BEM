#ifndef POINT_HPP_INCLUDED
#define POINT_HPP_INCLUDED

#include "../../global/commun.hpp"

/* Classe Point */

class Point
{
public:
    Real x;
    Real y;
    Point();
    Point(Real x, Real y);
    Point &operator+=(const Point &B);
    Point &operator-=(const Point &B);
    Point &operator*=(const Real a);
    Point &operator/=(const Real a);

    Real theta() const;
    Real norm() const;
};

/* Fontions associées */

ostream &operator<<(ostream &out, const Point &A);

Point operator+(const Point &A, const Point &B);

Point operator-(const Point &A, const Point &B);

Point operator*(const Point &A, const Real a);

Point operator*(const Real a, const Point &A);

Point operator/(const Point &A, const Real a);

Real operator*(const Point &A, const Point &B);

bool operator==(const Point &A, const Point &B);

Real operator|(const Point &A, const Point &B);

#endif
