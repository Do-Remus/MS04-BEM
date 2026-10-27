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
    Index rang = 0;
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

ACA_decomposition full_pivoted_cross_approximation(const Matrice &A, const Index rang_max, const double tolerance);

ACA_decomposition partial_pivoted_cross_approximation(const Matrice &A, const Index rang_max, const double tolerance);

inline Matrice vecteurs_vers_matrice(const std::vector<Vecteur> &vecteurs)
{
    if (vecteurs.empty())
        return Matrice(0, 0);

    const Index n = vecteurs[0].size();
    const Index r = static_cast<Index>(vecteurs.size());

    Matrice M(n, r);

    for (Index k = 0; k < r; ++k)
    {
        if (vecteurs[k].size() != n)
            throw std::invalid_argument("Tailles de vecteurs incompatibles");

        M.col(k) = vecteurs[k];
    }

    return M;
}

inline std::optional<Index> indice_max_absolu(const Vecteur &v)
{
    if (v.size() == 0)
        return std::nullopt;

    std::optional<Index> indice_max = std::nullopt;
    double valeur_max = SEUIL_PIVOT;

    for (Index i = 0; i < v.size(); ++i)
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

inline std::optional<Index> indice_max_absolu_non_utilise(const Vecteur &v, const std::vector<bool> &utilise)
{
    const Index n = v.size();

    if (static_cast<Index>(utilise.size()) != n)
        throw std::invalid_argument("Tailles incompatibles");

    std::optional<Index> indice_max = std::nullopt;
    double valeur_max = SEUIL_PIVOT;

    for (Index i = 0; i < n; ++i)
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

inline std::pair<std::optional<Index>, std::optional<Index>> indices_max_absolu(const Matrice &A)
{
    if (A.rows() == 0 || A.cols() == 0)
        return {std::nullopt, std::nullopt};

    std::optional<Index> i_max = std::nullopt;
    std::optional<Index> j_max = std::nullopt;
    double valeur_max = SEUIL_PIVOT;

    for (Index i = 0; i < A.rows(); ++i)
    {
        for (Index j = 0; j < A.cols(); ++j)
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