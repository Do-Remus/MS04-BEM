#ifndef BEM_ACA_DECOMPOSITION_HPP_INCLUDED
#define BEM_ACA_DECOMPOSITION_HPP_INCLUDED

#include "../../global/commun.hpp"
#include "../../maths/commun.hpp"

#include "../operator.hpp"
#include "../reconstruction.hpp"

ACA_decomposition aca_bem(const std::vector<Index> &selection_i, const std::vector<Index> &selection_j, const Maillage &maillage, const std::vector<QuadratureSegment> &quadrature_maillage, const Index rang_max, const double tolerance);

ACA_decomposition aca_reconstruction(const std::vector<Index> &selection_i, const std::vector<Index> &selection_j, const std::vector<Point> &pointsSolution, const std::vector<QuadratureSegment> &quadrature_maillage, const Index rang_max, const double tolerance);

template <typename FonctionCoefficient>
ACA_decomposition aca_generique(const std::vector<Index> &selection_i, const std::vector<Index> &selection_j, FonctionCoefficient coefficient, const Index rang_max, const double tolerance)
{
    const Index m = static_cast<Index>(selection_i.size());
    const Index n = static_cast<Index>(selection_j.size());

    ACA_decomposition aca;

    const Index rang_max_reel = std::min(rang_max, std::min(m, n));

    if (m == 0 || n == 0 || rang_max_reel <= 0)
    {
        aca.U = Matrice(m, 0);
        aca.V = Matrice(n, 0);
        aca.rang = 0;
        return aca;
    }

    std::vector<Vecteur> U;
    std::vector<Vecteur> V;

    U.reserve(rang_max_reel);
    V.reserve(rang_max_reel);

    std::vector<bool> ligne_utilisee(m, false);

    Index k = 0;
    Index i_pivot = 0;

    double norm_Frobenius_A_k_squared = 0.0;
    const double tolerance_squared = tolerance * tolerance;

    while (k < rang_max_reel)
    {
        // Calcul de la ligne du résidu au pivot i_pivot.
        Vecteur ligne_residuelle(n);

        for (Index j = 0; j < n; ++j)
        {
            ligne_residuelle[j] = coefficient(selection_i[i_pivot], selection_j[j]);
        }

        for (Index j = 0; j < k; ++j)
        {
            ligne_residuelle -= U[j][i_pivot] * V[j];
        }

        const auto j_pivot_opt = indice_max_absolu(ligne_residuelle);

        if (!j_pivot_opt)
            break;

        const Index j_pivot = *j_pivot_opt;

        // Stockage de la ligne résiduelle.
        V.push_back(ligne_residuelle);

        // Calcul de la colonne du résidu au pivot j_pivot.
        Vecteur colonne_residuelle(m);

        for (Index i = 0; i < m; ++i)
        {
            colonne_residuelle[i] = coefficient(selection_i[i], selection_j[j_pivot]);
        }

        for (Index i = 0; i < k; ++i)
        {
            colonne_residuelle -= V[i][j_pivot] * U[i];
        }

        const Complex pivot = ligne_residuelle[j_pivot];

        colonne_residuelle /= pivot;

        U.push_back(colonne_residuelle);

        // Mise à jour de la norme de Frobenius de l'approximation.
        const double contribution_squared = colonne_residuelle.squaredNorm() * ligne_residuelle.squaredNorm();

        Complex cross_sum(0.0, 0.0);

        for (Index i = 0; i < k; ++i)
        {
            cross_sum +=
                U[i].dot(colonne_residuelle) * V[i].dot(ligne_residuelle);
        }

        norm_Frobenius_A_k_squared += 2.0 * cross_sum.real() + contribution_squared;
        norm_Frobenius_A_k_squared = std::max(0.0, norm_Frobenius_A_k_squared);

        ++k;

        // Critère d'arrêt par stagnation.
        if (contribution_squared <= tolerance_squared * norm_Frobenius_A_k_squared)
        {
            break;
        }

        ligne_utilisee[i_pivot] = true;

        const auto i_pivot_opt = indice_max_absolu_non_utilise(colonne_residuelle, ligne_utilisee);

        if (!i_pivot_opt)
            break;

        i_pivot = *i_pivot_opt;
    }

    aca.U = vecteurs_vers_matrice(U);
    aca.V = vecteurs_vers_matrice(V);
    aca.rang = k;

    return aca;
}

#endif