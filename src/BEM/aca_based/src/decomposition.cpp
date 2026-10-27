#include "../decomposition.hpp"

ACA_decomposition aca_bem(const std::vector<Index> selection_i, const std::vector<Index> selection_j, const Maillage &maillage, const std::vector<QuadratureSegment> &quadrature_maillage, const Index rang_max, const double tolerance)
{
    // Index
    const Index m = static_cast<Index>(selection_i.size());
    const Index n = static_cast<Index>(selection_j.size());

    // Structure de retour
    ACA_decomposition aca;

    // Le rang ne peut pas dépasser min(m,n)
    const Index rang_max_reel = std::min(rang_max, std::min(m, n));

    // Cas limites
    if (m == 0 || n == 0 || rang_max_reel <= 0)
    {
        aca.U = Matrice(m, 0);
        aca.V = Matrice(n, 0);
        aca.rang = 0;

        return aca;
    }

    // Initialisation des matrices
    vector<Vecteur> U;
    U.reserve(rang_max_reel);

    vector<Vecteur> V;
    V.reserve(rang_max_reel);

    std::vector<bool> ligne_utilisee(m, false);

    // Boucle principale
    Index k = 0;
    Index i_pivot = 0;

    // ||A_k||_F^2 et <A, A_k>_F
    double norm_Frobenius_A_k_squared = 0.0;
    double tolerance_squared = tolerance * tolerance;

    while (k < rang_max_reel)
    {
        // Calcul de la ligne i_pivot du résidu,
        // sans construire toute la matrice résiduelle.
        Vecteur ligne_residuelle = ligne_bem(selection_i[i_pivot], selection_j, maillage, quadrature_maillage);

        for (Index j = 0; j < k; ++j)
            ligne_residuelle -= U[j][i_pivot] * V[j];

        // Pivot maximal sur la ligne sélectionnée
        std::optional<Index> j_pivot_opt = indice_max_absolu(ligne_residuelle);

        // Pivot nul sur cette ligne, elle n'est pas utilisable
        if (!j_pivot_opt)
        {
            break;
        }

        Index j_pivot = *j_pivot_opt;

        V.push_back(ligne_residuelle);

        // Calcul de la colonne du résidu, normalisée
        Vecteur colonne_residuelle = colonne_bem(selection_j[j_pivot], selection_i, maillage, quadrature_maillage);

        for (Index i = 0; i < k; ++i)
            colonne_residuelle -= V[i][j_pivot] * U[i];

        colonne_residuelle /= ligne_residuelle[j_pivot];

        U.push_back(colonne_residuelle);

        // Contribution de rang 1 : ||u_k v_k^T||_F^2
        const double contribution_squared = colonne_residuelle.squaredNorm() * ligne_residuelle.squaredNorm();

        // Mise a jour de ||A_k||_F^2
        Complex cross_sum(0.0, 0.0);

        for (Index i = 0; i < k; ++i)
            cross_sum += U[i].dot(colonne_residuelle) * V[i].dot(ligne_residuelle);

        norm_Frobenius_A_k_squared += 2.0 * cross_sum.real() + contribution_squared;
        norm_Frobenius_A_k_squared = std::max(0.0, norm_Frobenius_A_k_squared);

        ++k;

        // Critere de stagnation:
        // ||u_k|| ||v_k|| <= tolerance * ||A_k||_F
        if (contribution_squared <= tolerance_squared * norm_Frobenius_A_k_squared)
        {
            break;
        }

        // Marquer la ligne de pivot comme utilisee
        ligne_utilisee[i_pivot] = true;

        // Choisir la ligne suivante avec max |u[i]|
        // parmi les lignes non encore utilisees.
        std::optional<Index> i_pivot_opt = indice_max_absolu_non_utilise(colonne_residuelle, ligne_utilisee);

        // Toutes les lignes ont été utilisées
        if (i_pivot < 0)
            break;

        i_pivot = *i_pivot_opt;
    }

    // Conversion des facteurs en matrices
    aca.U = vecteurs_vers_matrice(U);
    aca.V = vecteurs_vers_matrice(V);
    aca.rang = k;

    return aca;
}
