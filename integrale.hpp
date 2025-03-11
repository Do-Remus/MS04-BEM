#ifndef INTEGRALE_HPP_INCLUDED
#define INTEGRALE_HPP_INCLUDED
#include <iostream>
#include <math.h>
#include "config.hpp"
#include "segment.hpp"

typedef complex<double> (*fun_d_P)(const Point &);
typedef complex<double> (*fun_d_P2)(const Point &, const Point &);

complex<double> integ_simple(Segment AB, fun_d_P f, int Nbpas1)
// permet de faire une intÈgrale simple sur le segment AB
// sur f:R^2->R, avec Nbpas1 intervalles sur AB
{
    complex<double> result = 0;
    double dAB = AB.norm() / Nbpas1;
    for (int i = 0; i < Nbpas1; i++)
    {
        result += dAB * f(AB.P1 + dAB * i * (AB.P2 - AB.P1));
    }
    return result;
}

complex<double> integ_double(Segment AB, Segment CD, fun_d_P2 f, int Nbpas1, int Nbpas2)
// permet de faire une intÈgrale double sur les segments AB et CD
// sur f:R^2xR^2->R, avec pas1 intervalles sur AB et pas2 intervalles sur CD
{
    complex<double> result = 0;
    double dAB = AB.norm() / Nbpas1;
    double dCD = CD.norm() / Nbpas2;
    double dS = abs((dAB * (AB.P2 - AB.P1)) * (dCD * (CD.P2 - CD.P1))); // produit vectoriel des deux petits vecteurs;
    for (int i = 0; i < Nbpas1; i++)
    {
        for (int j = 0; j < Nbpas2; j++)
        {
            result += dS * f(AB.P1 + dAB * (AB.P2 - AB.P1) * i, CD.P1 + dCD * (CD.P2 - CD.P1) * j);
        }
    }
    return result;
}

complex<double> log_norm(const Point &A, const Point &B)
{
    return log((A - B).norm());
}

complex<double> integ_log(Segment AB, Segment CD, int Nbpas1, int Nbpas2)
{
    if (AB == CD)
    {
        complex<double> a = AB.norm() * AB.norm() * log(AB.norm()) - 3 * AB.norm() * AB.norm() / 2;
        return a;
    }
    else if (AB.P1 == CD.P1 || AB.P1 == CD.P2 || AB.P2 == CD.P1 || AB.P2 == CD.P2)
    {
        complex<double> sp = AB.norm() / 2;
        complex<double> sm = -AB.norm() / 2;
        complex<double> tp = CD.norm() / 2;
        complex<double> tm = -CD.norm() / 2;
        complex<double> sigma = 1 / AB.norm() * (((CD.P1 + CD.P2) / 2 - (AB.P1 + AB.P2) / 2) | (AB.P1 - AB.P2));
        complex<double> sigmam = 1 / AB.norm() * ((CD.P2 - (AB.P1 + AB.P2) / 2) | (AB.P1 - AB.P2));
        complex<double> to = 1 / CD.norm() * (((AB.P1 + AB.P2) / 2 - (CD.P1 + CD.P2) / 2) | (CD.P1 - CD.P2));
        complex<double> top = 1 / CD.norm() * ((AB.P1 - (CD.P1 + CD.P2) / 2) | (CD.P1 - CD.P2));
        complex<double> cost = (AB.P1 - AB.P2) | (CD.P1 - CD.P2) / AB.norm() / CD.norm();
        complex<double> sint = (AB.P1 - AB.P2) * (CD.P1 - CD.P2) / AB.norm() / CD.norm();
        complex<double> dalphap = tp * sint - (to + sigma * cost) / sint;
        complex<double> dbetap = sp * sint - (to + sigma * cost) / sint;
        complex<double> a = AB.norm() / 2 * ((tp - top) * log((CD.P1 - AB.P1).norm()) - (tm - top) * log((CD.P2 - AB.P1).norm()));
        complex<double> b = CD.norm() / 2 * ((sp - sigmam) * log((AB.P1 - CD.P2).norm()) - (sm - sigmam) * log((AB.P2 - CD.P2).norm()));
        complex<double> c = dbetap * (atan((tp - top) / dbetap) - atan((tm - top) / dbetap));
        complex<double> d = dalphap * (atan((sp - sigmam) / dalphap) - atan((sm - sigmam) / dalphap));
        complex<double> e = -3 * AB.norm() * CD.norm() / 2;
        return a + b + c + d + e;
    }
    else
    {
        return (integ_double(AB, CD, log_norm, Nbpas1, Nbpas2));
    }
}

#endif