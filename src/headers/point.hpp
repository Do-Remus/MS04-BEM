#ifndef POINT_HPP_INCLUDED
#define POINT_HPP_INCLUDED
#include "../config/config.hpp"
#include "../config/external.hpp"
#include "../config/constantes.hpp"

/* Classe Point */

class Point
{
public:
    double x;
    double y;
    Point();
    Point(double x, double y);
    Point &operator+=(const Point &B);
    Point &operator-=(const Point &B);
    Point &operator*=(const double a);
    Point &operator/=(const double a);
    double norm();
};

/* Fontions associées */

ostream &operator<<(ostream &out, const Point &A);

Point operator+(const Point &A, const Point &B);

Point operator-(const Point &A, const Point &B);

Point operator*(const Point &A, const double a);

Point operator*(const double a, const Point &A);

Point operator/(const Point &A, const double a);

double operator*(const Point &A, const Point &B);

bool operator==(const Point &A, const Point &B);

double operator|(const Point A, const Point B);

#endif
