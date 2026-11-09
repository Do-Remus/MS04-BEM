#include "../arbre_blocs.hpp"

namespace
{
    void ajouter_bloc_dense(const Bloc &bloc, Matrice &A, const bool est_symetrique)
    {
        // Sous-arbre sans feuille sur ce rang : rien n'a été calculé.
        if (bloc.nb_feuilles_rank_mpi == 0)
            return;

        if (bloc.type == TypeBloc::Interne)
        {
            for (const auto &ligne : bloc.fils)
                for (const auto &enfant : ligne)
                    if (enfant)
                        ajouter_bloc_dense(*enfant, A, est_symetrique);
            return;
        }

        const auto &I = bloc.cluster_lignes->indices;
        const auto &J = bloc.cluster_colonnes->indices;

        Matrice local;
        if (bloc.type == TypeBloc::Dense)
            local = bloc.dense;
        else // Admissible : A_loc ≈ U * V
            local = bloc.U * bloc.V.transpose();

        // Hors diagonale, le bloc n'est stocké que d'un côté si la matrice
        // est symétrique : on écrit aussi la transposée (pas le conjugué).
        const bool miroir = est_symetrique && (bloc.cluster_lignes != bloc.cluster_colonnes);

        for (std::size_t i = 0; i < I.size(); ++i)
        {
            for (std::size_t j = 0; j < J.size(); ++j)
            {
                A(I[i], J[j]) = local(i, j);
                if (miroir)
                    A(J[j], I[i]) = local(i, j);
            }
        }
    }

    // MPI_Allreduce sur un tampon complexe, découpé en morceaux
    // pour ne pas dépasser la limite int du paramètre count.
    void allreduce_complexe(Complex *data, const std::size_t taille, MPI_Comm comm)
    {
        const std::size_t morceau_max = static_cast<std::size_t>(1) << 28; // 2^28 complexes
        for (std::size_t debut = 0; debut < taille; debut += morceau_max)
        {
            const int count = static_cast<int>(std::min(morceau_max, taille - debut));
            MPI_Allreduce(MPI_IN_PLACE, data + debut, count, MPI_C_DOUBLE_COMPLEX, MPI_SUM, comm);
        }
    }
}

// Retourne la matrice complète, identique sur tous les rangs.
// Doit être appelée par TOUS les rangs (appel collectif).
Matrice convertion_aca_dense(const ArbreBlocs &arbre_blocs, const bool est_symetrique = true)
{
    const Index N = arbre_blocs->cluster_lignes->indices.size();
    const Index M = arbre_blocs->cluster_colonnes->indices.size();

    Matrice A = Matrice::Zero(N, M);
    ajouter_bloc_dense(*arbre_blocs, A, est_symetrique);

    allreduce_complexe(A.data(), static_cast<std::size_t>(A.size()), MPI_COMM_WORLD);

    return A;
}

void diagnostic_blocs(const Bloc &bloc, const Matrice &A_ref, const bool est_symetrique)
{
    if (bloc.nb_feuilles_rank_mpi == 0)
        return;

    if (bloc.type == TypeBloc::Interne)
    {
        for (const auto &ligne : bloc.fils)
            for (const auto &enfant : ligne)
                if (enfant)
                    diagnostic_blocs(*enfant, A_ref, est_symetrique);
        return;
    }

    const auto &I = bloc.cluster_lignes->indices;
    const auto &J = bloc.cluster_colonnes->indices;

    Matrice ref(I.size(), J.size());
    for (std::size_t i = 0; i < I.size(); ++i)
        for (std::size_t j = 0; j < J.size(); ++j)
            ref(i, j) = A_ref(I[i], J[j]);

    Matrice local;
    if (bloc.type == TypeBloc::Dense)
        local = bloc.dense;
    else
    {
        if (bloc.U.cols() != bloc.V.cols()){
            throw std::logic_error("U/V : rangs incompatibles");
        }
        local = bloc.U * bloc.V.transpose();
    }

    const double err = (local - ref).norm() / ref.norm();
    std::cout << (bloc.type == TypeBloc::Dense ? "DENSE " : "ADMIS ")
              << I.size() << "x" << J.size()
              << " diag=" << bloc.est_diagonal
              << " rang=" << (bloc.type == TypeBloc::Dense ? 0 : bloc.U.cols())
              << " err_rel=" << err << std::endl;
}

