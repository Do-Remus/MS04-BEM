#include "../operator.hpp"

void appliquer_bloc(const Bloc &bloc, const Vecteur &x, Vecteur &y)
{
    // Aucun descendant appartenant à ce rang :
    // inutile de parcourir ce sous-arbre.
    if (bloc.nb_feuilles_rank_mpi == 0)
    {
        return;
    }

    // Bloc interne : appliquer récursivement les quatre sous-blocs.
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

    // Une feuille ne contribue que sur son processus propriétaire.
    if (bloc.proprietaire_mpi != mpi_rank)
        return;

    const auto &indices_i = bloc.cluster_lignes->indices;
    const auto &indices_j = bloc.cluster_colonnes->indices;

    const Index m = indices_i.size();
    const Index n = indices_j.size();

    // Extraction du vecteur correspondant au cluster de colonnes.
    Vecteur x_local(n);

    for (Index j = 0; j < n; ++j)
    {
        x_local[j] = x[indices_j[j]];
    }

    Vecteur y_local(m);

    if (bloc.type == TypeBloc::Dense)
    {
        y_local.noalias() = bloc.dense * x_local;
    }
    else if (bloc.type == TypeBloc::Admissible)
    {
        const Vecteur temporaire = bloc.V.transpose() * x_local;
        y_local.noalias() = bloc.U * temporaire;
    }
    else
    {
        throw std::logic_error("Type de bloc inconnu.");
    }

    // Ajout de la contribution aux lignes globales.
    for (Index i = 0; i < m; ++i)
    {
        y[indices_i[i]] += y_local[i];
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
