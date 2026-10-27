#include "../arbre_blocs.hpp"

bool sont_admissibles(const BoiteEnglobante &boite_a, const BoiteEnglobante &boite_b, const double eta)
{
    return std::min(boite_a.diametre, boite_b.diametre) <= eta * distance(boite_a, boite_b);
}

ArbreBlocs construire_bloc_recursif(const Noeud &tau, const Noeud &sigma, const double eta, const double alpha)
{
    // Création du bloc associé aux deux clusters.
    auto bloc = std::make_unique<Bloc>();

    bloc->cluster_lignes = &tau;
    bloc->cluster_colonnes = &sigma;

    // Un bloc est diagonal si ses deux clusters sont le même nœud.
    bloc->est_diagonal = (&tau == &sigma);

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

        bloc->fils[0][0] = construire_bloc_recursif(*tau.fils_gauche, *sigma.fils_gauche, eta, alpha);
        bloc->fils[0][1] = construire_bloc_recursif(*tau.fils_gauche, *sigma.fils_droit, eta, alpha);
        bloc->fils[1][1] = construire_bloc_recursif(*tau.fils_droit, *sigma.fils_droit, eta, alpha);
        if (!bloc->est_diagonal)
            bloc->fils[1][0] = construire_bloc_recursif(*tau.fils_droit, *sigma.fils_gauche, eta, alpha);

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

ArbreBlocs construire_arbre_blocs(const Arbre &arbre_clusters, const double eta, const double alpha)
{
    if (!arbre_clusters)
        return nullptr;

    // Étape 1 : construction de la structure de blocs.
    ArbreBlocs racine = construire_bloc_recursif(*arbre_clusters, *arbre_clusters, eta, alpha);

    // Étape 2 : attribution pondérée des propriétaires.
    repartir_feuilles_ponderee(*racine, mpi_size);

    // Étape 3 : compteurs correspondant au rang courant.
    recalculer_feuilles_locales(*racine, mpi_rank);

    return racine;
}

void remplir_blocs(Bloc &bloc, const Maillage &maillage, const std::vector<QuadratureSegment> &quadrature_maillage, const Index rang_max, const double tolerance)
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
                    remplir_blocs(*enfant, maillage, quadrature_maillage, rang_max, tolerance);
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

void ecrire_arbre_blocs(const ArbreBlocs &arbre, const std::string &fichier, const IntervallesClusters &intervalles)
{
    std::ofstream out(fichier);

    if (!out)
        throw std::runtime_error("Impossible d'ouvrir " + fichier);

    out << "id,parent,niveau,type,est_diagonal,m,n,"
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

    auto visiter = [&](auto &&self, const Bloc &bloc, long long parent_id, std::size_t niveau) -> void
    {
        const std::size_t id = prochain_id++;

        const auto &l = bloc.cluster_lignes->boite;
        const auto &c = bloc.cluster_colonnes->boite;

        const auto [i_debut, i_fin] = intervalles.at(bloc.cluster_lignes);
        const auto [j_debut, j_fin] = intervalles.at(bloc.cluster_colonnes);

        out << id << ','
            << parent_id << ','
            << niveau << ','
            << nom_type(bloc.type) << ','
            << (bloc.est_diagonal ? 1 : 0) << ','
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
