#ifndef MATRICE_HPP_INCLUDED
#define MATRICE_HPP_INCLUDED
#include "../config/config.hpp"
#include "../config/external.hpp"
#include "../config/constantes.hpp"

/* Class Vecteur */

class Vecteur : public std::vector<std::complex<double>>
{
public:
    using Complex = std::complex<double>;
    using std::vector<Complex>::vector; // récupère les constructeurs de vector

    // Produit scalaire hermitien
    Complex operator*(const Vecteur &v) const
    {
        Complex resultat = 0.0;

        for (std::size_t i = 0; i < size(); ++i)
            resultat += std::conj((*this)[i]) * v[i];

        return resultat;
    }

    // Addition
    Vecteur operator+(const Vecteur &v) const
    {
        Vecteur resultat(size());

        for (std::size_t i = 0; i < size(); ++i)
            resultat[i] = (*this)[i] + v[i];

        return resultat;
    }

    // Soustraction
    Vecteur operator-(const Vecteur &v) const
    {
        Vecteur resultat(size());

        for (std::size_t i = 0; i < size(); ++i)
            resultat[i] = (*this)[i] - v[i];

        return resultat;
    }

    // Scalaire * vecteur
    Vecteur operator*(const Complex &a) const
    {
        Vecteur resultat(size());

        for (std::size_t i = 0; i < size(); ++i)
            resultat[i] = a * (*this)[i];

        return resultat;
    }

    // Vecteur * scalaire
    friend Vecteur operator*(const Complex &a, const Vecteur &v)
    {
        return v * a;
    }

    // Norme euclidienne
    double norm() const
    {
        return std::sqrt(std::real((*this) * (*this)));
    }
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
    Matrice(int n, int m, complex<double> v = (complex<double>)0);
    complex<double> &operator()(int i, int j);
    complex<double> operator()(int i, int j) const;
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
    MatriceSym(int m = 0, complex<double> v = (complex<double>)0);
    MatriceSym(const Vecteur &d);
    complex<double> &operator()(int i, int j);
    complex<double> operator()(int i, int j) const;
    void decomposition_LDL(MatriceSym &L, Vecteur &D) const;
};

ostream &operator<<(ostream &out, const MatriceSym &A);

Vecteur resolution_systeme_lineaire(const MatriceSym &A, const Vecteur &P);

Vecteur gradientConjugue(const Matrice &A, const Vecteur &b, double tol = 1e-10, int maxIter = 1000);

#endif
