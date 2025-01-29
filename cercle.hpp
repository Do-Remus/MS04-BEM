#ifndef Cercle_HPP_INCLUDED
#define Cercle_HPP_INCLUDED

#include "point.hpp"

class Cercle  
{
public:
    Point centre;
    double rayon = 1.;
    Cercle(){};
    Cercle(const double rayon, const Point& A );
};



Cercle::Cercle(const double r, const Point& A){
    if(r<=0){
        cout<<"error : rayon cercle doit etre strictement positif, r= " <<r<<endl;
        exit(-1);
    }
    this->rayon=r;
    this->centre=A;
}


#endif 