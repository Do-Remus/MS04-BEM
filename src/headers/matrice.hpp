#ifndef MATRICE_HPP_INCLUDED
#define MATRICE_HPP_INCLUDED
#include "../config/config.hpp"
#include "../config/external.hpp"
#include "../config/constantes.hpp"

/* Typedef Vecteur */

typedef vector<complex<double>> Vecteur;

/* Fonctions associées au type Vecteur */

ostream &operator<<(ostream &out, const Vecteur &u);


/*Classe Matrice générale*/
class Matrice 
{
    protected:
        Vecteur coefs; //coeficients rangés par ligne
    public:
        int n; //nombre de lignes
        int m; //nombre de colonnes
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

#endif
