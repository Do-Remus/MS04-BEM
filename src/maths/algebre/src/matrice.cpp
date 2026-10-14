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

ACA_decomposition full_pivoted_cross_approximation(const Matrice &A, const Eigen::Index rang_max, const double tolerance)
{
    // Index
    const Eigen::Index m = A.rows();
    const Eigen::Index n = A.cols();

    // Structure de retour
    ACA_decomposition aca;

    // Norme de A
    const double norm_A = A.norm();

    // Cas limites
    if (m == 0 || n == 0 || norm_A == 0.0 || rang_max <= 0 || tolerance >= 1.0)
    {
        aca.U = Matrice(m, 0);
        aca.V = Matrice(n, 0);
        aca.rang = 0;

        return aca;
    }

    // Le rang ne peut pas dépasser min(m,n)
    const Eigen::Index rang_max_reel = std::min(rang_max, std::min(m, n));

    // Initialisation des matrices
    Matrice R_k = A;

    vector<Vecteur> U;
    U.reserve(rang_max_reel);

    vector<Vecteur> V;
    V.reserve(rang_max_reel);

    // Boucle principale
    Eigen::Index k = 0;

    while (k < rang_max_reel && R_k.norm() > tolerance * norm_A)
    {
        // Trouver le pivot (valeur maximale en valeur absolue )
        const auto [i_pivot, j_pivot] = indices_max_absolu(R_k);

        // Aucun pivot non nul : le résidu est nul
        if (i_pivot < 0 || std::abs(R_k(i_pivot, j_pivot)) == 0.0)
            break;

        // Mise à jour des vecteurs
        V.push_back(R_k.row(i_pivot).transpose());
        U.push_back(R_k.col(j_pivot) / R_k(i_pivot, j_pivot));

        // Mise à jour du résidu
        R_k -= U.back() * V.back().transpose();

        // Mise à jour du rang
        ++k;
    }

    aca.U = vecteurs_vers_matrice(U);
    aca.V = vecteurs_vers_matrice(V);
    aca.rang = k;

    return aca;
}

ACA_decomposition partial_pivoted_cross_approximation(const Matrice &A, const Eigen::Index rang_max, const double tolerance)
{
    // Index
    const Eigen::Index m = A.rows();
    const Eigen::Index n = A.cols();

    // Structure de retour
    ACA_decomposition aca;

    // Le rang ne peut pas dépasser min(m,n)
    const Eigen::Index rang_max_reel = std::min(rang_max, std::min(m, n));

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
    Eigen::Index k = 0;
    Eigen::Index i_pivot = 0;

    // ||A_k||_F^2 et <A, A_k>_F
    double norm_Frobenius_A_k_squared = 0.0;
    double tolerance_squared = tolerance * tolerance;

    while (k < rang_max_reel)
    {
        // Calcul de la ligne i_pivot du résidu,
        // sans construire toute la matrice résiduelle.
        Vecteur ligne_residuelle = A.row(i_pivot).transpose();

        for (Eigen::Index j = 0; j < k; ++j)
            ligne_residuelle -= U[j][i_pivot] * V[j];

        // Pivot maximal sur la ligne sélectionnée
        Eigen::Index j_pivot = indice_max_absolu(ligne_residuelle);

        // Pivot nul sur cette ligne, elle n'est pas utilisable
        if (j_pivot < 0)
        {
            break;
        }

        V.push_back(ligne_residuelle);

        // Calcul de la colonne du résidu, normalisée
        Vecteur colonne_residuelle = A.col(j_pivot);

        for (Eigen::Index i = 0; i < k; ++i)
            colonne_residuelle -= V[i][j_pivot] * U[i];

        colonne_residuelle /= ligne_residuelle[j_pivot];

        U.push_back(colonne_residuelle);

        // Contribution de rang 1 : ||u_k v_k^T||_F^2
        const double contribution_squared = colonne_residuelle.squaredNorm() * ligne_residuelle.squaredNorm();

        // Mise a jour de ||A_k||_F^2
        Complex cross_sum(0.0, 0.0);

        for (Eigen::Index i = 0; i < k; ++i)
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
        i_pivot = indice_max_absolu_non_utilise(colonne_residuelle, ligne_utilisee);

        // Toutes les lignes ont été utilisées
        if (i_pivot < 0)
            break;
    }

    // Conversion des facteurs en matrices
    aca.U = vecteurs_vers_matrice(U);
    aca.V = vecteurs_vers_matrice(V);
    aca.rang = k;

    return aca;
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
