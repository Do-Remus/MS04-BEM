#include "../operator.hpp"

void remplir_blocs_bem(Bloc &bloc, const Maillage &maillage, const std::vector<QuadratureSegment> &quadrature_maillage, const Index rang_max, const double tolerance)
{
    // Aucun calcul à effectuer si ce sous-arbre ne contient
    // aucune feuille appartenant au rang MPI courant.
    if (bloc.nb_feuilles_rank_mpi == 0)
    {
        return;
    }

    // --------------------------------------------------
    // 1. Bloc interne : remplir récursivement les fils.
    // --------------------------------------------------
    if (bloc.type == TypeBloc::Interne)
    {
        for (auto &ligne : bloc.fils)
        {
            for (auto &enfant : ligne)
            {
                if (enfant) // Pour un bloc diagonal fils[1][0] == nullptr et donc ne sera pas exploré
                    remplir_blocs_bem(*enfant, maillage, quadrature_maillage, rang_max, tolerance);
            }
        }

        return;
    }

    // À ce stade, le bloc est une feuille appartenant
    // nécessairement au rang MPI courant.
    // Récupérer les indices des clusters associés au bloc.
    const auto &indices_lignes = bloc.cluster_lignes->indices;
    const auto &indices_colonnes = bloc.cluster_colonnes->indices;

    // --------------------------------------------------
    // 2. Bloc dense : calculer et stocker les coefficients.
    // --------------------------------------------------
    if (bloc.type == TypeBloc::Dense)
    {
        const Index m = indices_lignes.size();
        const Index n = indices_colonnes.size();

        bloc.dense.resize(m, n);

        if (bloc.est_diagonal)
        {
            // Le bloc est carré : calcul de la moitié supérieure.
            for (Index i = 0; i < m; ++i)
            {
                for (Index j = i; j < n; ++j)
                {
                    const Complex aij = matrice_bem(indices_lignes[i], indices_colonnes[j], maillage, quadrature_maillage);

                    bloc.dense(i, j) = aij;

                    if (i != j)
                        bloc.dense(j, i) = aij;
                }
            }
        }
        else
        {
            // Bloc non diagonal : remplissage classique.
            for (Eigen::Index i = 0; i < m; ++i)
            {
                for (Eigen::Index j = 0; j < n; ++j)
                {
                    bloc.dense(i, j) = matrice_bem(indices_lignes[i], indices_colonnes[j], maillage, quadrature_maillage);
                }
            }
        }

        // Aucun facteur ACA n'est nécessaire pour ce bloc.
        bloc.U.resize(0, 0);
        bloc.V.resize(0, 0);

        return;
    }

    // --------------------------------------------------
    // 3. Bloc admissible : calculer l'approximation ACA.
    // --------------------------------------------------
    if (bloc.type == TypeBloc::Admissible)
    {
        auto resultat = aca_bem(indices_lignes, indices_colonnes, maillage, quadrature_maillage, rang_max, tolerance);

        // Transférer les facteurs ACA dans le bloc.
        bloc.U = std::move(resultat.U);
        bloc.V = std::move(resultat.V);

        // Ce bloc ne stocke pas sa matrice dense.
        bloc.dense.resize(0, 0);

        return;
    }

    std::cerr
        << "[Erreur] Type de bloc inconnu"
        << " | type = " << static_cast<int>(bloc.type)
        << " | rang MPI = " << mpi_rank
        << " | feuilles locales = " << bloc.nb_feuilles_rank_mpi
        << " | taille lignes = " << bloc.cluster_lignes->indices.size()
        << " | taille colonnes = " << bloc.cluster_colonnes->indices.size()
        << '\n';

    throw std::logic_error("Type de bloc inconnu.");
}
