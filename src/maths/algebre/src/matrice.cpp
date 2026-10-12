#include "../matrice.hpp"

/* Fonctions associées à la classe Matrice */

Matrice matrice_aleatoire(const unsigned int rang, const Index taille)
{
    if (rang > taille)
        throw std::invalid_argument("Le rang ne peut pas dépasser la taille.");

    if (rang == 0 || taille == 0)
        return Matrice::Zero(taille, taille);

    Matrice U = Matrice::Random(taille, rang);
    Matrice V = Matrice::Random(rang, taille);

    return U * V;
}

// ostream &operator<<(ostream &out, const Matrice &A)
// {
//     for (int i = 0; i < A.n; i++)
//     {
//         for (int j = 0; j < A.m; j++)
//         {
//             out << A(i, j) << " ";
//         }
//         out << endl;
//     }
//     return out;
// }

// /* Méthodes de la classe Matrice */

// Matrice::Matrice(int N, int M, complex<Real> v)
// {
//     m = M;
//     n = N;
//     coefs.resize(n * m, v);
//     // constructeur matrice a valeurs constantes
// }

// complex<Real> &Matrice::operator()(int i, int j)
// {
//     if (i >= n || j >= m || i < 0 || j < 0)
//     {
//         cout << "overflow in matrix coordinates ; i= " << i << ", j= " << j << " matrix size (n,m) = (" << n << "," << m << ")" << endl;
//         std::exit(1);
//     }
//     return coefs[i * m + j];
// }

// complex<Real> Matrice::operator()(int i, int j) const
// {
//     if (i >= n || j >= m || i < 0 || j < 0)
//     {
//         cout << "overflow in matrix coordinates ; i= " << i << ", j= " << j << " matrix size (n,m) = (" << n << "," << m << ")" << endl;
//         std::exit(1);
//     }
//     return coefs[i * m + j];
// }

// Vecteur Matrice::operator()(int i) const
// {
//     if (i >= n || i < 0)
//     {
//         cout << "overflow in matrix coordinates ; i= " << i << " n lines from matrix size (n,m) = (" << n << "," << m << ")" << endl;
//         std::exit(1);
//     }
//     Vecteur ligne(m, 0);
//     for (int j = 0; j < m; j++)
//         ligne[j] = (*this)(i, j);
//     return ligne;
// }

// Vecteur Matrice::operator*(const Vecteur &v) const
// {
//     Vecteur b(n);
//     for (int i = 0; i < n; i++)
//     {
//         b[i] = 0;
//         for (int j = 0; j < m; j++)
//         {
//             b[i] += (*this)(i, j) * v[j];
//         }
//     }
//     return b;
// }
