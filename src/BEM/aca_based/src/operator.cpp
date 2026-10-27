#include "../operator.hpp"

void appliquer_bloc(const Bloc &bloc, const Vecteur &x, Vecteur &y)
{
    // Aucun descendant appartenant à ce rang :
    // inutile de parcourir ce sous-arbre.
    if (bloc.nb_feuilles_rank_mpi == 0)
    {
        return;
    }

    // --------------------------------------------------
    // 1. Bloc interne : parcourir les enfants existants.
    // --------------------------------------------------
    if (bloc.type == TypeBloc::Interne)
    {
        for (Index i = 0; i < 2; ++i)
        {
            for (Index j = 0; j < 2; ++j)
            {
                if (bloc.fils[i][j])
                {
                    appliquer_bloc(*bloc.fils[i][j], x, y);
                }
            }
        }

        return;
    }

    // --------------------------------------------------
    // 2. Une feuille ne travaille que sur son rang MPI.
    // --------------------------------------------------
    if (bloc.proprietaire_mpi != mpi_rank)
        return;

    const auto &indices_i = bloc.cluster_lignes->indices;
    const auto &indices_j = bloc.cluster_colonnes->indices;

    const Index m = indices_i.size();
    const Index n = indices_j.size();

    // --------------------------------------------------
    // 3. Contribution directe : y_I += B * x_J.
    // --------------------------------------------------

    // Extraction du vecteur global.
    Vecteur x_j(n);

    for (Index j = 0; j < n; ++j)
        x_j[j] = x[indices_j[j]];

    Vecteur y_i(m);

    if (bloc.type == TypeBloc::Dense)
    {
        y_i.noalias() = bloc.dense * x_j;
    }
    else if (bloc.type == TypeBloc::Admissible)
    {
        const Vecteur temporaire = bloc.V.transpose() * x_j;
        y_i.noalias() = bloc.U * temporaire;
    }
    else
    {
        throw std::logic_error("Type de bloc inconnu.");
    }

    // Ajout de la contribution aux lignes globales.
    for (Index i = 0; i < m; ++i)
        y[indices_i[i]] += y_i[i];

    // --------------------------------------------------
    // 4. Contribution transposée :
    //    y_J += B^T * x_I.
    //
    //    Elle est nécessaire pour les blocs non diagonaux
    //    dont la branche symétrique n'a pas été construite.
    // --------------------------------------------------

    if (!bloc.est_diagonal)
    {
        // Extraction du vecteur global.
        Vecteur x_i(m);

        for (Index i = 0; i < m; ++i)
            x_i[i] = x[indices_i[i]];

        Vecteur y_j(n);

        if (bloc.type == TypeBloc::Dense)
        {
            y_j.noalias() = bloc.dense.transpose() * x_i;
        }
        else // TypeBloc::Admissible
        {
            const Vecteur temporaire =
                bloc.U.transpose() * x_i;

            y_j.noalias() = bloc.V * temporaire;
        }

        for (Index j = 0; j < n; ++j)
            y[indices_j[j]] += y_j[static_cast<Eigen::Index>(j)];
    }
}

Vecteur appliquer_H(const ArbreBlocs &arbre_blocs, const Vecteur &x)
{
    if (!arbre_blocs)
    {
        throw std::invalid_argument("Impossible d'appliquer une H-matrice vide.");
    }

    if (x.size() != static_cast<Eigen::Index>(arbre_blocs->cluster_colonnes->indices.size()))
    {
        throw std::invalid_argument("La dimension du vecteur x est incompatible.");
    }

    Vecteur y = Vecteur::Zero(x.size());

    appliquer_bloc(*arbre_blocs, x, y);

    // Combiner les contributions de tous les processus.
    if (mpi_size > 1)
    {
        if (y.size() > std::numeric_limits<int>::max())
        {
            throw std::overflow_error("Le vecteur est trop grand pour MPI_Allreduce.");
        }

        MPI_Allreduce(MPI_IN_PLACE, y.data(), static_cast<int>(y.size()), MPI_CXX_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
    }

    return y;
}
