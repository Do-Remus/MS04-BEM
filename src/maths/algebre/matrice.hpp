#ifndef MATRICE_HPP_INCLUDED
#define MATRICE_HPP_INCLUDED

#include "../../global/commun.hpp"

#include "vecteur.hpp"

constexpr double SEUIL_PIVOT = 1e-15;

using Matrice = Eigen::MatrixXcd;

struct ACA_decomposition
{
    Matrice U;
    Matrice V;
    Eigen::Index rang = 0;
};

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

ACA_decomposition full_pivoted_cross_approximation(const Matrice &A, const Eigen::Index rang_max, const double tolerance);

ACA_decomposition partial_pivoted_cross_approximation(const Matrice &A, const Eigen::Index rang_max, const double tolerance);

inline Matrice vecteurs_vers_matrice(const std::vector<Vecteur> &vecteurs)
{
    if (vecteurs.empty())
        return Matrice(0, 0);

    const Eigen::Index n = vecteurs[0].size();
    const Eigen::Index r = static_cast<Eigen::Index>(vecteurs.size());

    Matrice M(n, r);

    for (Eigen::Index k = 0; k < r; ++k)
    {
        if (vecteurs[k].size() != n)
            throw std::invalid_argument("Tailles de vecteurs incompatibles");

        M.col(k) = vecteurs[k];
    }

    return M;
}

inline Eigen::Index indice_max_absolu(const Vecteur &v)
{
    if (v.size() == 0)
        return -1;

    Eigen::Index indice_max = -1;
    double valeur_max = SEUIL_PIVOT;

    for (Eigen::Index i = 0; i < v.size(); ++i)
    {
        const double valeur = std::abs(v[i]);

        if (valeur > valeur_max)
        {
            valeur_max = valeur;
            indice_max = i;
        }
    }

    return indice_max;
}

inline Eigen::Index indice_max_absolu_non_utilise(const Vecteur &v, const std::vector<bool> &utilise)
{
    const Eigen::Index n = v.size();

    if (static_cast<Eigen::Index>(utilise.size()) != n)
        throw std::invalid_argument("Tailles incompatibles");

    Eigen::Index indice_max = -1;
    double valeur_max = SEUIL_PIVOT;

    for (Eigen::Index i = 0; i < n; ++i)
    {
        if (utilise[i])
            continue;

        const double valeur = std::abs(v[i]);

        if (valeur > valeur_max)
        {
            valeur_max = valeur;
            indice_max = i;
        }
    }

    return indice_max;
}

inline std::pair<Eigen::Index, Eigen::Index> indices_max_absolu(const Matrice &A)
{
    if (A.rows() == 0 || A.cols() == 0)
        return {-1, -1};

    Eigen::Index i_max = -1;
    Eigen::Index j_max = -1;
    double valeur_max = SEUIL_PIVOT;

    for (Eigen::Index i = 0; i < A.rows(); ++i)
    {
        for (Eigen::Index j = 0; j < A.cols(); ++j)
        {
            const double valeur = std::abs(A(i, j));

            if (valeur > valeur_max)
            {
                valeur_max = valeur;
                i_max = i;
                j_max = j;
            }
        }
    }

    return {i_max, j_max};
}

#endif