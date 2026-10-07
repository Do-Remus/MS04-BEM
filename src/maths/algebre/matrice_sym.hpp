#ifndef MATRICE_SYM_HPP_INCLUDED
#define MATRICE_SYM_HPP_INCLUDED

#include "../../global/commun.hpp"

#include "vecteur.hpp"

class MatriceSym
{
protected:
    Vecteur coefs; // coefficients rangés par ligne
public:
    int n; // ordre de la MatriceSym
    MatriceSym(int m = 0, complex<Real> v = (complex<Real>)0);
    MatriceSym(const Vecteur &d);
    complex<Real> &operator()(int i, int j);
    Vecteur operator*(const Vecteur &v) const;
    complex<Real> operator()(int i, int j) const;
    void decomposition_LDL(MatriceSym &L, Vecteur &D) const;
};

/* Fonctions associées à la class Matrice Sym */

Vecteur resolution_systeme_lineaire(const MatriceSym &A, const Vecteur &P);

Vecteur COCG(const MatriceSym &A, const Vecteur &b, Real tol = 1e-10, unsigned int maxIter = 1000);

#endif
