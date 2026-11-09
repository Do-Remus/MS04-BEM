
#include "../reconstruction.hpp"


static const char *nom_type_bloc(TypeBloc t)
{
    switch (t)
    {
    case TypeBloc::Dense:      return "DENSE";
    case TypeBloc::Admissible: return "ADMISSIBLE";
    default:                   return "INTERNE";
    }
}

// Format d'une ligne du fichier (une ligne par feuille) :
// TYPE m n niveau_lignes niveau_colonnes est_diagonal nb_valeurs s_0 s_1 ... s_{k-1}
// Les valeurs singulières sont écrites jusqu'à ce que s_k / s_0 < seuil_svd
// (la première valeur sous le seuil est incluse).
void construire_bloc_forcage_dense(Bloc &bloc,
                                   const std::vector<Point> &pointsSolution,
                                   const std::vector<QuadratureSegment> &quadrature_maillage,
                                   std::ostream &sortie,
                                   const double seuil_svd)
{
    if (bloc.nb_feuilles_rank_mpi == 0)
        return;

    if (bloc.type == TypeBloc::Interne)
    {
        for (auto &ligne : bloc.fils)
            for (auto &enfant : ligne)
                if (enfant)
                    construire_bloc_forcage_dense(*enfant, pointsSolution, quadrature_maillage, sortie, seuil_svd);
        return;
    }

    const auto &indices_lignes = bloc.cluster_lignes->indices;
    const auto &indices_colonnes = bloc.cluster_colonnes->indices;

    const Index m = indices_lignes.size();
    const Index n = indices_colonnes.size();
    bloc.dense.resize(m, n);

    for (Index i = 0; i < m; ++i)
        for (Index j = 0; j < n; ++j)
            bloc.dense(i, j) = coefficient_reconstruction(indices_lignes[i], indices_colonnes[j],
                                                          pointsSolution, quadrature_maillage);

    // Valeurs singulières seulement (pas de U/V) : beaucoup plus rapide.
    Eigen::BDCSVD<Matrice> svd(bloc.dense);
    const auto &sigma = svd.singularValues(); // triées par ordre décroissant

    Index nb = 0;
    if (sigma.size() > 0)
    {
        const double s0 = sigma(0);
        for (Index k = 0; k < sigma.size(); ++k)
        {
            ++nb;
            if (s0 > 0.0 && sigma(k) / s0 < seuil_svd)
                break;
        }
    }

    sortie << nom_type_bloc(bloc.type) << ' '
           << m << ' ' << n << ' '
           << bloc.cluster_lignes->niveau << ' ' << bloc.cluster_colonnes->niveau << ' '
           << static_cast<int>(bloc.est_diagonal) << ' '
           << nb << std::scientific << std::setprecision(8);
    for (Index k = 0; k < nb; ++k)
        sortie << ' ' << sigma(k);
    sortie << '\n';

    // Aucun facteur ACA pour ce bloc.
    bloc.U.resize(0, 0);
    bloc.V.resize(0, 0);
}

void remplir_blocs_solution(Bloc &bloc, const std::vector<Point> &pointsSolution, const std::vector<QuadratureSegment> &quadrature_maillage, const Index rang_max, const double tolerance)
{
    // Aucun calcul si ce sous-arbre ne contient
    // aucune feuille appartenant au rang MPI courant.
    if (bloc.nb_feuilles_rank_mpi == 0)
        return;

    // 1. Bloc interne : parcourir les enfants existants.
    if (bloc.type == TypeBloc::Interne)
    {
        for (auto &ligne : bloc.fils)
        {
            for (auto &enfant : ligne)
            {
                if (enfant)
                {
                    remplir_blocs_solution(*enfant, pointsSolution, quadrature_maillage, rang_max, tolerance);
                }
            }
        }

        return;
    }

    // La feuille appartient nécessairement au rang MPI courant.
    const auto &indices_lignes = bloc.cluster_lignes->indices;
    const auto &indices_colonnes = bloc.cluster_colonnes->indices;

    const Index m = indices_lignes.size();
    const Index n = indices_colonnes.size();

    // 2. Bloc dense : calculer tous les coefficients.
    if (bloc.type == TypeBloc::Dense)
    {
        bloc.dense.resize(m, n);

        for (Index i = 0; i < m; ++i)
        {
            for (Index j = 0; j < n; ++j)
            {
                bloc.dense(i, j) = coefficient_reconstruction(indices_lignes[i], indices_colonnes[j], pointsSolution, quadrature_maillage);
            }
        }

        // Aucun facteur ACA pour ce bloc.
        bloc.U.resize(0, 0);
        bloc.V.resize(0, 0);

        return;
    }

    // 3. Bloc admissible : approximation ACA.
    if (bloc.type == TypeBloc::Admissible)
    {
        auto resultat = aca_reconstruction(indices_lignes, indices_colonnes, pointsSolution, quadrature_maillage, rang_max, tolerance);

        bloc.U = std::move(resultat.U);
        bloc.V = std::move(resultat.V);

        bloc.dense.resize(0, 0);

        return;
    }

    throw std::logic_error("Type de bloc inconnu dans remplir_blocs_solution.");
}
