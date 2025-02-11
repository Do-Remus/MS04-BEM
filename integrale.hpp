#include <iostream>
#include "point.hpp"


double integ(Point A,Point B, Point C, Point D, std::function<double(Point,Point)> f,int pas1,int pas2)
//permet de faire une intÈgrale double sur les segments AB et CD
//sur f:R^2xR^2->R, avec pas1 intervalles sur AB et pas2 intervalles sur CD
{
    double result = 0;
    Point dAB=(B-A)/pas1;
    Point dCD=(D-C)/pas2;
    double dS=abs(dAB*dCD);//produit vect des deux petits vecteurs;
    for(int i=0; i<pas1;i++){
        for(int j=0;j<pas2;j++){
            result+=dS*f(A+dAB*i,C+dCD*j);
        }
    }
}