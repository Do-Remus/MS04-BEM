#include "headers/cercle.hpp"

Cercle::Cercle(const double r, const Point &A)
{
    if (r <= 0)
    {
        cout << "error : rayon cercle doit etre strictement positif, r= " << r << endl;
        exit(-1);
    }
    this->rayon = r;
    this->centre = A;
}
