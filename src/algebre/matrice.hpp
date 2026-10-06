#ifndef MATRICE_HPP_INCLUDED
#define MATRICE_HPP_INCLUDED

#include "../global/global.hpp"
#include "vecteur.hpp"

class Matrice
{
protected:
    Vecteur coefs; // coeficients rangés par ligne
public:
    int n; // nombre de lignes
    int m; // nombre de colonnes

    Matrice(int n, int m, complex<Real> v = (complex<Real>)0);

    virtual complex<Real> &operator()(int i, int j);
    virtual complex<Real> operator()(int i, int j) const;
    virtual Vecteur operator*(const Vecteur &v) const;

    virtual ~Matrice() = default;
};

ostream &operator<<(ostream &out, const Matrice &A);

Matrice generateur_aleatoire_matrice(const int rang, const int taille);

#endif