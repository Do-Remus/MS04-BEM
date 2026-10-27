#ifndef BEM_OPERATOR_HPP_INCLUDED
#define BEM_OPERATOR_HPP_INCLUDED

#include "../global/commun.hpp"
#include "../maths/commun.hpp"

Matrice get_matrice_bem(const Maillage &maillage, const std::vector<QuadratureSegment> &quadrature_maillage);

Vecteur ligne_bem(const Index i, const Maillage &maillage, const std::vector<QuadratureSegment> &quadrature_maillage);

Vecteur colonne_bem(const Index j, const Maillage &maillage, const std::vector<QuadratureSegment> &quadrature_maillage);

Vecteur ligne_bem(const Index i, const std::vector<Index> &selection_j, const Maillage &maillage, const std::vector<QuadratureSegment> &quadrature_maillage);

Vecteur colonne_bem(const Index j, const std::vector<Index> &selection_i, const Maillage &maillage, const std::vector<QuadratureSegment> &quadrature_maillage);

// CACHED
inline Complex matrice_bem(const Index i, const Index j, const Maillage &maillage, const std::vector<QuadratureSegment> &quadrature_maillage)
{
    const Segment &Si = maillage[i];
    const QuadratureSegment &qi = quadrature_maillage[i];

    /*
     * Diagonale régularisée.
     */
    if (i == j)
        return integ_double_log(Si) + integ_double(green_reguliere_cached, qi, qi);

    /*
     * Hors-diagonale : calcul d'un seul triangle.
     */
    return integ_double(green_cached, qi, quadrature_maillage[j]);
}

// Calcul de la diagonal de l'opérateur A
inline Vecteur diagonale_matrice_bem(const Maillage &maillage, const std::vector<QuadratureSegment> &quadrature_maillage)
{
    const Index N = maillage.size();
    Vecteur diagonaleLocale = Vecteur::Zero(N);

    /*
     * Chaque processus calcule une partie de la diagonale.
     */
    for (Index i = mpi_rank; i < N; i += mpi_size)
    {
        const Segment &S = maillage[i];
        const QuadratureSegment &qS = quadrature_maillage[i];

        Complex aii = integ_double_log(S);
        aii += integ_double(green_reguliere_cached, qS, qS);

        diagonaleLocale[i] = aii;
    }

    /*
     * Chaque processus récupère la diagonale complète.
     */
    Vecteur diagonale = Vecteur::Zero(N);

    MPI_Allreduce(diagonaleLocale.data(), diagonale.data(), static_cast<int>(N), mpi_complex_type(), MPI_SUM, MPI_COMM_WORLD);

    return diagonale;
}

inline Vecteur preconditionneur_diagonal(const Vecteur &diagonale)
{
    const Index N = diagonale.size();
    Vecteur S = Vecteur::Zero(N);

    for (Index i = 0; i < N; ++i)
    {
        S[i] = 1.0 / std::sqrt(diagonale[i]);
    }

    return S;
}

#endif