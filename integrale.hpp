#include <iostream>
#include <math.h>
#include "segment.hpp"

typedef double (*fun_d_P)(const Point &);
typedef double (*fun_d_P2)(const Point &, const Point &);

double integ_simple(Segment AB, fun_d_P f, int Nbpas1)
// permet de faire une intÈgrale simple sur le segment AB
// sur f:R^2->R, avec Nbpas1 intervalles sur AB
{
    double result = 0;
    double dAB = AB.norm() / Nbpas1;
    for (int i = 0; i < Nbpas1; i++)
    {
        result += dAB * f(AB.P1 + dAB * i * (AB.P2 - AB.P1));
    }
}

double integ_double(Segment AB, Segment CD, fun_d_P2 f, int Nbpas1, int Nbpas2)
// permet de faire une intÈgrale double sur les segments AB et CD
// sur f:R^2xR^2->R, avec pas1 intervalles sur AB et pas2 intervalles sur CD
{
    double result = 0;
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
}

double integ_log(Segment AB, Segment CD, int Nbpas)
{
}
