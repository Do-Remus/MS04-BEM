#include "../operator.hpp"

Matrice get_matrice_bem(const Maillage &maillage, const std::vector<QuadratureSegment> &quadrature_maillage)
{
    const Index N = maillage.size();
    Matrice A = Matrice::Zero(N, N);

    for (Index i = 0; i < N; ++i)
    {
        const Segment &Si = maillage[i];
        const QuadratureSegment &qi = quadrature_maillage[i];

        /*
         * Diagonale régularisée.
         */
        A(i, i) = integ_double_log(Si) + integ_double(green_reguliere, qi, qi);

        /*
         * Hors-diagonale : calcul d'un seul triangle.
         */
        for (Index j = i + 1; j < N; ++j)
        {
            const Complex aij = integ_double(green, qi, quadrature_maillage[j]);

            A(i, j) = aij;
            A(j, i) = aij;
        }
    }

    return A;
}

Vecteur ligne_bem(const Index i, const Maillage &maillage, const std::vector<QuadratureSegment> &quadrature_maillage)
{
    const Index N = static_cast<Index>(maillage.size());

    Vecteur ligne(N);

    for (Index j = 0; j < N; ++j)
    {
        ligne[j] = matrice_bem(i, j, maillage, quadrature_maillage);
    }

    return ligne;
}

Vecteur colonne_bem(const Index j, const Maillage &maillage, const std::vector<QuadratureSegment> &quadrature_maillage)
{
    const Index N = static_cast<Index>(maillage.size());

    Vecteur colonne(N);

    for (Index i = 0; i < N; ++i)
    {
        colonne[i] = matrice_bem(i, j, maillage, quadrature_maillage);
    }

    return colonne;
}

Vecteur ligne_bem(const Index i, const std::vector<Index> &selection_j, const Maillage &maillage, const std::vector<QuadratureSegment> &quadrature_maillage)
{
    const Index n = static_cast<Index>(selection_j.size());

    Vecteur ligne(n);

    for (Index j_local = 0; j_local < n; ++j_local)
    {
        const Index j_global = selection_j[j_local];

        ligne[j_local] = matrice_bem(i, j_global, maillage, quadrature_maillage);
    }

    return ligne;
}

Vecteur colonne_bem(const Index j, const std::vector<Index> &selection_i, const Maillage &maillage, const std::vector<QuadratureSegment> &quadrature_maillage)
{
    const Index m = static_cast<Index>(selection_i.size());

    Vecteur colonne(m);

    for (Index i_local = 0; i_local < m; ++i_local)
    {
        const Index i_global = selection_i[i_local];

        colonne[i_local] = matrice_bem(i_global, j, maillage, quadrature_maillage);
    }

    return colonne;
}

ACA_decomposition aca_bem(const Maillage &maillage, const std::vector<QuadratureSegment> &quadrature_maillage, const Index rang_max, const double tolerance)
{
    // Index
    const Index m = static_cast<Index>(maillage.size());
    const Index n = m;

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
        Vecteur ligne_residuelle = ligne_bem(i_pivot, maillage, quadrature_maillage);

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
        Vecteur colonne_residuelle = colonne_bem(j_pivot, maillage, quadrature_maillage);

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
        if (!i_pivot_opt)
            break;

        i_pivot = *i_pivot_opt;
    }

    // Conversion des facteurs en matrices
    aca.U = vecteurs_vers_matrice(U);
    aca.V = vecteurs_vers_matrice(V);
    aca.rang = k;

    return aca;
}
