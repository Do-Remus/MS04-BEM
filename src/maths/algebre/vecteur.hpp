#ifndef VECTEUR_HPP_INCLUDED
#define VECTEUR_HPP_INCLUDED

#include "../../global/commun.hpp"

using Vecteur = Eigen::VectorXcd;

// class Vecteur : public std::vector<Complex>
// {
// public:
//     using std::vector<Complex>::vector; // récupère les constructeurs de vector

//     Vecteur operator+(const Vecteur &v) const;
//     Vecteur operator-(const Vecteur &v) const;
//     Vecteur operator*(const Complex &a) const;

//     friend Vecteur operator*(const Complex &a, const Vecteur &v)
//     {
//         return v * a;
//     }

//     Real norm() const;
// };

// // Fonctions associées au type Vecteur

// ostream &operator<<(ostream &out, const Vecteur &u);

// Vecteur generateur_aleatoire_vecteur(const int taille);

Complex produitBilineaire(const Vecteur &x, const Vecteur &y);

#endif