#include "../arbre_blocs.hpp"

bool sont_admissibles(const BoiteEnglobante &boite_a, const BoiteEnglobante &boite_b, const double eta)
{
    return std::min(boite_a.diametre, boite_b.diametre) <= eta * distance(boite_a, boite_b);
}

ArbreBlocs construire_bloc_recursif(const Noeud &tau, const Noeud &sigma, const double eta, const double alpha, const bool est_symetrique)
{
    // Création du bloc associé aux deux clusters.
    auto bloc = std::make_unique<Bloc>();

    bloc->cluster_lignes = &tau;
    bloc->cluster_colonnes = &sigma;

    // Un bloc est diagonal si ses deux clusters sont le même nœud et la matrice est considéré comme symétrique.
    bloc->est_diagonal = (&tau == &sigma) && est_symetrique;

    // 0. Cluster vide trouvé
    if (tau.indices.empty() || sigma.indices.empty())
        throw std::logic_error("Un cluster vide a été rencontré dans l'arbre des blocs.");

    // 1. Le bloc est admissible : c'est une feuille. (ne peut fondamentalement pas être diagonal)
    if (sont_admissibles(tau.boite, sigma.boite, eta))
    {
        bloc->type = TypeBloc::Admissible;

        const double m = static_cast<double>(tau.indices.size());
        const double n = static_cast<double>(sigma.indices.size());

        bloc->poids = alpha * (m + n);

        return bloc;
    }
    // 2. Le bloc n'est pas admissible et ne peut plus
    // être subdivisé si l'un des clusters est une feuille.
    else if (tau.est_feuille() || sigma.est_feuille())
    {
        bloc->type = TypeBloc::Dense;

        const double m = static_cast<double>(tau.indices.size());
        const double n = static_cast<double>(sigma.indices.size());

        bloc->poids = m * n;

        return bloc;
    }
    // 3. Le bloc n'est pas admissible :
    // on le subdivise en quatre blocs fils.
    else
    {
        bloc->type = TypeBloc::Interne;

        bloc->fils[0][0] = construire_bloc_recursif(*tau.fils_gauche, *sigma.fils_gauche, eta, alpha, est_symetrique);
        bloc->fils[0][1] = construire_bloc_recursif(*tau.fils_gauche, *sigma.fils_droit, eta, alpha, est_symetrique);
        bloc->fils[1][1] = construire_bloc_recursif(*tau.fils_droit, *sigma.fils_droit, eta, alpha, est_symetrique);

        // En mode symétrique, on ne construit pas le bloc miroir
        // d'un bloc diagonal : il sera pris en compte par transposition
        if (!bloc->est_diagonal)
            bloc->fils[1][0] = construire_bloc_recursif(*tau.fils_droit, *sigma.fils_gauche, eta, alpha, est_symetrique);

        return bloc;
    }
}

void collecter_feuilles(Bloc &bloc, std::vector<Bloc *> &feuilles)
{
    if (bloc.type != TypeBloc::Interne)
    {
        feuilles.push_back(&bloc);
        return;
    }

    for (auto &ligne : bloc.fils)
    {
        for (auto &enfant : ligne)
        {
            if (enfant) // Pour un bloc diagonal fils[1][0] == nullptr et donc ne sera pas exploré
                collecter_feuilles(*enfant, feuilles);
        }
    }
}

void repartir_feuilles_ponderee(Bloc &racine, const int mpi_size)
{
    std::vector<Bloc *> feuilles;
    collecter_feuilles(racine, feuilles);

    // Les feuilles les plus coûteuses sont traitées en premier.
    std::stable_sort(feuilles.begin(), feuilles.end(), [](const Bloc *a, const Bloc *b)
                     { return a->poids > b->poids; });

    std::vector<double> charge(static_cast<std::size_t>(mpi_size), 0.0);

    for (Bloc *feuille : feuilles)
    {
        const auto it = std::min_element(charge.begin(), charge.end());
        const int rang = static_cast<int>(std::distance(charge.begin(), it));

        feuille->proprietaire_mpi = rang;
        charge[static_cast<std::size_t>(rang)] += feuille->poids;
    }
}

Index recalculer_feuilles_locales(Bloc &bloc, const int mpi_rank)
{
    if (bloc.type != TypeBloc::Interne)
    {
        bloc.nb_feuilles_rank_mpi = (bloc.proprietaire_mpi == mpi_rank) ? 1 : 0;

        return bloc.nb_feuilles_rank_mpi;
    }

    Index total = 0;

    for (auto &ligne : bloc.fils)
    {
        for (auto &enfant : ligne)
        {
            if (enfant)
                total += recalculer_feuilles_locales(*enfant, mpi_rank);
        }
    }

    bloc.nb_feuilles_rank_mpi = total;
    return total;
}

ArbreBlocs construire_arbre_blocs(const Arbre &arbre_clusters, const double eta, const double alpha, const bool est_symetrique)
{
    return construire_arbre_blocs(arbre_clusters, arbre_clusters, eta, alpha, est_symetrique);
}

ArbreBlocs construire_arbre_blocs(const Arbre &arbre_lignes, const Arbre &arbre_colonnes, const double eta, const double alpha, const bool est_symetrique)
{
    if (!arbre_lignes || !arbre_colonnes)
        return nullptr;

    // Étape 1 : construction de la structure de blocs.
    ArbreBlocs racine = construire_bloc_recursif(*arbre_lignes, *arbre_colonnes, eta, alpha, est_symetrique);

    // Étape 2 : attribution pondérée des propriétaires.
    repartir_feuilles_ponderee(*racine, mpi_size);

    // Étape 3 : compteurs correspondant au rang courant.
    recalculer_feuilles_locales(*racine, mpi_rank);

    return racine;
}

void ecrire_arbre_blocs(const ArbreBlocs &arbre, const std::string &fichier, const IntervallesClusters &intervalles_lignes, const IntervallesClusters &intervalles_colonnes, const bool est_symetrique)
{
    std::ofstream out(fichier);

    if (!out)
        throw std::runtime_error("Impossible d'ouvrir " + fichier);

    out << "id,parent,niveau,type,est_diagonal,est_symetrique,m,n,"
           "i_debut,i_fin,j_debut,j_fin,"
           "x_min_l,x_max_l,y_min_l,y_max_l,"
           "x_min_c,x_max_c,y_min_c,y_max_c\n";

    out << std::setprecision(17);

    if (!arbre)
        return;

    auto nom_type = [](TypeBloc type) -> const char *
    {
        switch (type)
        {
        case TypeBloc::Interne:
            return "Interne";
        case TypeBloc::Admissible:
            return "Admissible";
        case TypeBloc::Dense:
            return "Dense";
        }
        return "Inconnu";
    };

    std::size_t prochain_id = 0;

    auto visiter = [&](auto &&self, const Bloc &bloc, const long long parent_id, const std::size_t niveau) -> void
    {
        const std::size_t id = prochain_id++;

        const auto &l = bloc.cluster_lignes->boite;
        const auto &c = bloc.cluster_colonnes->boite;

        const auto [i_debut, i_fin] = intervalles_lignes.at(bloc.cluster_lignes);
        const auto [j_debut, j_fin] = intervalles_colonnes.at(bloc.cluster_colonnes);

        out << id << ','
            << parent_id << ','
            << niveau << ','
            << nom_type(bloc.type) << ','
            << (bloc.est_diagonal ? 1 : 0) << ','
            << (est_symetrique ? 1 : 0) << ','
            << bloc.cluster_lignes->indices.size() << ','
            << bloc.cluster_colonnes->indices.size() << ','
            << i_debut << ','
            << i_fin << ','
            << j_debut << ','
            << j_fin << ','
            << l.x_min << ','
            << l.x_max << ','
            << l.y_min << ','
            << l.y_max << ','
            << c.x_min << ','
            << c.x_max << ','
            << c.y_min << ','
            << c.y_max << '\n';

        for (const auto &ligne : bloc.fils)
        {
            for (const auto &enfant : ligne)
            {
                if (enfant)
                    self(self, *enfant, static_cast<long long>(id), niveau + 1);
            }
        }
    };

    visiter(visiter, *arbre, -1, 0);
}

void compter_blocs_locaux(const Bloc &bloc, const int rang, StatistiquesBlocs &stats)
{
    if (bloc.type == TypeBloc::Interne)
    {
        for (const auto &ligne : bloc.fils)
        {
            for (const auto &enfant : ligne)
            {
                if (enfant)
                {
                    compter_blocs_locaux(*enfant, rang, stats);
                }
            }
        }
        return;
    }

    if (bloc.proprietaire_mpi != rang)
        return;

    ++stats.feuilles;

    const std::uint64_t m = static_cast<std::uint64_t>(bloc.cluster_lignes->indices.size());
    const std::uint64_t n = static_cast<std::uint64_t>(bloc.cluster_colonnes->indices.size());

    if (bloc.type == TypeBloc::Dense)
    {
        ++stats.denses;
        stats.coefficients_denses += m * n;
    }
    else if (bloc.type == TypeBloc::Admissible)
    {
        ++stats.admissibles;

        const std::uint64_t r = static_cast<std::uint64_t>(bloc.U.cols());

        stats.coefficients_aca += r * (m + n);
    }
}

void appliquer_bloc(const Bloc &bloc, const Vecteur &x, Vecteur &y, const bool est_symetrique)
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
                    appliquer_bloc(*bloc.fils[i][j], x, y, est_symetrique);
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
    //    Si la matrice est symétrique.
    //    Elle est nécessaire pour les blocs non diagonaux
    //    dont la branche symétrique n'a pas été construite.
    // --------------------------------------------------

    if (!bloc.est_diagonal && est_symetrique)
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

Vecteur appliquer_H(const ArbreBlocs &arbre_blocs, const Vecteur &x, const bool est_symetrique)
{
    if (!arbre_blocs)
    {
        throw std::invalid_argument("Impossible d'appliquer une H-matrice vide.");
    }

    const Eigen::Index nb_lignes = static_cast<Eigen::Index>(arbre_blocs->cluster_lignes->indices.size());
    const Eigen::Index nb_colonnes = static_cast<Eigen::Index>(arbre_blocs->cluster_colonnes->indices.size());

    if (x.size() != nb_colonnes)
    {
        throw std::invalid_argument("La dimension du vecteur x est incompatible.");
    }

    Vecteur y = Vecteur::Zero(nb_lignes);

    appliquer_bloc(*arbre_blocs, x, y, est_symetrique);

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
