#ifndef MATRICE_HPP_INCLUDED
#define MATRICE_HPP_INCLUDED

#include "../../global/commun.hpp"

#include "vecteur.hpp"

using Matrice = Eigen::MatrixXcd;

// class Matrice
// {
// protected:
//     Vecteur coefs; // coeficients rangés par ligne
// public:
//     int n; // nombre de lignes
//     int m; // nombre de colonnes

//     Matrice(int n, int m, complex<Real> v = (complex<Real>)0);

//     virtual complex<Real> &operator()(int i, int j);
//     virtual complex<Real> operator()(int i, int j) const;
//     //lecture d'une ligne
//     virtual Vecteur operator()(int i) const;
//     virtual Vecteur operator*(const Vecteur &v) const;

//     virtual ~Matrice() = default;
// };

// ostream &operator<<(ostream &out, const Matrice &A);

Matrice matrice_aleatoire(const unsigned int rang, const Index taille);

vector<Vecteur> full_pivoted_cross_approximation(const Matrice& A, const Index rang, double tolerance);

vector<Vecteur> partial_pivoted_cross_approximation(const Matrice& A, const Index rang, double tolerance);

#endif