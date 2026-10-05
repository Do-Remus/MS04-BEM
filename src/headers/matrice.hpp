#ifndef MATRICE_HPP_INCLUDED
#define MATRICE_HPP_INCLUDED
#include "../config/config.hpp"
#include "../config/external.hpp"
#include "../config/constantes.hpp"

/* Class Vecteur */

class Vecteur : public std::vector<Complex>
{
public:
    using std::vector<Complex>::vector; // récupère les constructeurs de vector

    Complex produitHermitien(const Vecteur &v) const;
    Complex produitBilineaire(const Vecteur &v) const;
    Vecteur operator+(const Vecteur &v) const;
    Vecteur operator-(const Vecteur &v) const;
    Vecteur operator*(const Complex &a) const;

    friend Vecteur operator*(const Complex &a, const Vecteur &v)
    {
        return v * a;
    }

    Real norm() const;
};

/* Fonctions associées au type Vecteur */

ostream &operator<<(ostream &out, const Vecteur &u);

/*Classe Matrice générale*/
class Matrice
{
protected:
    Vecteur coefs; // coeficients rangés par ligne
public:
    int n; // nombre de lignes
    int m; // nombre de colonnes
    Matrice(int n, int m, complex<Real> v = (complex<Real>)0);
    complex<Real> &operator()(int i, int j);
    complex<Real> operator()(int i, int j) const;
    Vecteur operator*(const Vecteur v) const;
};

ostream &operator<<(ostream &out, const Matrice &A);

/* Classe MatriceSym */

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

ostream &operator<<(ostream &out, const MatriceSym &A);

Vecteur resolution_systeme_lineaire(const MatriceSym &A, const Vecteur &P);

Vecteur gradConjMatSym(const MatriceSym &A, const Vecteur &b, Real tol = 1e-10, unsigned int maxIter = 1000);

Matrice generateur_matrice(const int rang, const int taille);


Vecteur generateur_aléatoire_vecteur(const int taille);


#endif
