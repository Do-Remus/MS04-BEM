#include "../conditionners.hpp"

// LUT

std::pair<Permutation, Eigen::VectorXd> construire_ordre_et_echelles(const std::vector<Coordonnees> &centres)
{
    const Index N = centres.size();

    if (N == 0)
        throw std::invalid_argument("Le maillage est vide.");

    if (N > static_cast<Index>(std::numeric_limits<int>::max()))
        throw std::overflow_error("MPI_DOUBLE_INT ne peut pas indexer ce maillage.");

    Permutation permutation(N);
    Eigen::VectorXd echelles(N);

    if (N == 1)
    {
        permutation[0] = 0;
        echelles(0) = 0.0;
        return {std::move(permutation), std::move(echelles)};
    }

    // --------------------------------------------------
    // Répartition des indices entre les processus
    // --------------------------------------------------

    const Index debut = N * static_cast<Index>(mpi_rank) / static_cast<Index>(mpi_size);
    const Index fin = N * static_cast<Index>(mpi_rank + 1) / static_cast<Index>(mpi_size);

    const Index taille_locale = fin - debut;

    // --------------------------------------------------
    // Coordonnées
    // --------------------------------------------------

    std::vector<double> x(N);
    std::vector<double> y(N);

    double xmin = centres[0](0);
    double ymin = centres[0](1);

    for (Index i = 0; i < N; ++i)
    {
        x[i] = centres[i](0);
        y[i] = centres[i](1);

        xmin = std::min(xmin, x[i]);
        ymin = std::min(ymin, y[i]);
    }

    // --------------------------------------------------
    // Premier point
    // --------------------------------------------------

    Index premier = 0;
    double meilleure_distance = std::numeric_limits<double>::max();

    for (Index i = 0; i < N; ++i)
    {
        const double dx = x[i] - xmin;
        const double dy = y[i] - ymin;
        const double d2 = dx * dx + dy * dy;

        if (d2 < meilleure_distance)
        {
            meilleure_distance = d2;
            premier = i;
        }
    }

    permutation[N - 1] = premier;

    // --------------------------------------------------
    // Distances minimales locales
    // --------------------------------------------------

    std::vector<double> distances_min(taille_locale);

    double distance_voisin_min_locale = std::numeric_limits<double>::max();

    for (Index i = debut; i < fin; ++i)
    {
        const double dx = x[i] - x[premier];
        const double dy = y[i] - y[premier];
        const double d2 = dx * dx + dy * dy;

        distances_min[i - debut] = d2;

        if (i != premier)
            distance_voisin_min_locale = std::min(distance_voisin_min_locale, d2);
    }

    double distance_voisin_min = 0.0;

    MPI_Allreduce(&distance_voisin_min_locale, &distance_voisin_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

    echelles(N - 1) = std::sqrt(distance_voisin_min);

    // IMPORTANT : exclure le premier point des sélections suivantes.
    if (premier >= debut && premier < fin)
        distances_min[premier - debut] = -1.0;

    // --------------------------------------------------
    // Reverse max-min
    // --------------------------------------------------

    struct CoupleMaxLoc
    {
        double distance;
        int indice;
    };

    for (Index k = N - 1; k-- > 0;)
    {
        // Maximum local.
        CoupleMaxLoc local{-1.0, std::numeric_limits<int>::max()};

        for (Index i = 0; i < taille_locale; ++i)
        {
            if (distances_min[i] > local.distance)
            {
                local.distance = distances_min[i];
                local.indice = static_cast<int>(debut + i);
            }
        }

        // Maximum global, avec départage par indice minimal.
        CoupleMaxLoc global{};

        MPI_Allreduce(&local, &global, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);

        const Index meilleur = static_cast<Index>(global.indice);
        const double distance_max = global.distance;

        permutation[k] = meilleur;
        echelles(k) = std::sqrt(distance_max);

        const double xm = x[meilleur];
        const double ym = y[meilleur];

        // Mise à jour uniquement de la portion locale.
        for (Index i = 0; i < taille_locale; ++i)
        {
            const Index indice_global = debut + i;

            const double dx = x[indice_global] - xm;
            const double dy = y[indice_global] - ym;

            const double d2 = dx * dx + dy * dy;

            if (d2 < distances_min[i])
                distances_min[i] = d2;
        }

        // Le point sélectionné doit rester exclu des sélections suivantes.
        if (meilleur >= debut && meilleur < fin)
            distances_min[meilleur - debut] = -1.0;
    }

    std::vector<bool> indices_vus(N, false);

    return {std::move(permutation), std::move(echelles)};
}

SupportsDistribues construire_supports_sparse_KD(const std::vector<Coordonnees> &centres, const Permutation &permutation, const Eigen::VectorXd &echelles, const double rho)
{
    const Index N = centres.size();

    if (permutation.size() != N || echelles.size() != static_cast<Eigen::Index>(N))
        throw std::invalid_argument("Dimensions incompatibles.");

    if (rho <= 0.0)
        throw std::invalid_argument("rho doit être strictement positif.");

    const Index debut = N * static_cast<Index>(mpi_rank) / static_cast<Index>(mpi_size);
    const Index fin = N * static_cast<Index>(mpi_rank + 1) / static_cast<Index>(mpi_size);

    // KD-tree local à chaque processus.
    KDTree2D arbre(centres);

    // Inverse de la permutation.
    std::vector<Index> position(N);

    for (Index j = 0; j < N; ++j)
        position[permutation[j]] = j;

    SupportsDistribues resultat;
    resultat.debut = debut;
    resultat.supports.resize(fin - debut);

    std::vector<Index> candidats;
    candidats.reserve(32);

    for (Index j = debut; j < fin; ++j)
    {
        auto &support = resultat.supports[j - debut];

        const Index indice_j = permutation[j];
        const auto &yj = centres[indice_j];

        const double xj = yj(0);
        const double y_coord = yj(1);

        const double rayon = rho * echelles(j);

        support.push_back(j);

        if (rayon <= 0.0)
            continue;

        arbre.rechercher_dans_rayon(
            xj, y_coord, rayon, candidats);

        for (const Index indice_i : candidats)
        {
            const Index i = position[indice_i];

            if (i <= j)
                continue;

            const auto &yi = centres[indice_i];

            const double dx = yi(0) - xj;
            const double dy = yi(1) - y_coord;

            const double seuil = rho * std::min(echelles(i), echelles(j));

            if (dx * dx + dy * dy < seuil * seuil)
                support.push_back(i);
        }

        std::sort(support.begin() + 1, support.end());
    }

    return resultat;
}

FacteursSparseDistribues construire_facteurs_sparse(const AccesMatrice &acces_A, const Permutation &permutation, const SupportsDistribues &supports_distribues, const double tolerance_pivot)
{
    const Index N = permutation.size();

    const Index debut = supports_distribues.debut;
    const Index fin = debut + supports_distribues.supports.size();

    if (fin > N)
        throw std::runtime_error("Intervalle de supports invalide.");

    const Eigen::Index n = static_cast<Eigen::Index>(N);

    std::vector<TripletSparse> triplets_L;
    std::vector<TripletSparse> triplets_UT;

    // Estimation du nombre de coefficients locaux.
    Index nombre_coefficients = 0;

    for (const auto &support : supports_distribues.supports)
        nombre_coefficients += support.size();

    triplets_L.reserve(nombre_coefficients);
    triplets_UT.reserve(nombre_coefficients);

    for (Index j = debut; j < fin; ++j)
    {
        const auto &support = supports_distribues.supports[j - debut];

        const Index m = support.size();

        if (m == 0)
            continue;

        // Position de j dans son support.
        Index position_j = m;

        for (Index a = 0; a < m; ++a)
        {
            if (support[a] == j)
            {
                position_j = a;
                break;
            }
        }

        if (position_j == m)
            throw std::runtime_error("Le support ne contient pas son indice central.");

        // Extraction de A_{S_j,S_j}.
        Matrice A_local(m, m);

        for (Index a = 0; a < m; ++a)
        {
            const Eigen::Index ia = static_cast<Eigen::Index>(permutation[support[a]]);

            for (Index b = 0; b < m; ++b)
            {
                const Eigen::Index ib = static_cast<Eigen::Index>(permutation[support[b]]);

                A_local(a, b) = acces_A(ia, ib);
            }
        }

        // Résolution A_local q = e_j.

        // Vecteur e_j.
        Vecteur e = Vecteur::Zero(m);
        e(position_j) = Complex(1.0, 0.0);

        // Construction de L_S : A_local * qL = e_j
        Eigen::FullPivLU<Matrice> lu(A_local);
        lu.setThreshold(tolerance_pivot);

        if (lu.rank() != static_cast<Eigen::Index>(m))
        {
            throw std::runtime_error(
                "Sous-matrice singuliere dans la construction de L.");
        }

        const Vecteur qL = lu.solve(e);
        const Complex q0 = qL(position_j);

        if (std::abs(q0) <= tolerance_pivot * std::max(qL.norm(), 1e-30))
        {
            throw std::runtime_error(
                "Normalisation impossible dans la construction de L.");
        }

        // Construction de U_S^T : A_local^T * qUT = e_j
        Eigen::FullPivLU<Matrice> luT(A_local.transpose());
        luT.setThreshold(tolerance_pivot);

        if (luT.rank() != static_cast<Eigen::Index>(m))
        {
            throw std::runtime_error(
                "Sous-matrice singuliere dans la construction de UT.");
        }

        const Vecteur qUT = luT.solve(e);

        // Colonne j de L_S.
        // Ligne j de U_S^T.
        for (Index a = 0; a < m; ++a)
        {
            const Eigen::Index i = static_cast<Eigen::Index>(support[a]);
            const Eigen::Index jj = static_cast<Eigen::Index>(j);

            triplets_L.emplace_back(i, jj, qL(a) / q0);
            triplets_UT.emplace_back(jj, i, qUT(a));
        }
    }

    FacteursSparseDistribues facteurs;

    facteurs.debut = debut;
    facteurs.fin = fin;

    facteurs.L_locale.resize(n, n);
    facteurs.UT_locale.resize(n, n);

    facteurs.L_locale.setFromTriplets(triplets_L.begin(), triplets_L.end());
    facteurs.UT_locale.setFromTriplets(triplets_UT.begin(), triplets_UT.end());

    facteurs.L_locale.makeCompressed();
    facteurs.UT_locale.makeCompressed();

    return facteurs;
}

FacteursSparseDistribues construire_preconditionneur_sparse(const Maillage &maillage, const std::vector<QuadratureSegment> &quadrature_maillage, const double rho, const double tolerance_pivot)
{
    // Récupération des centres des éléments du maillage.
    std::vector<Coordonnees> centres;
    centres.reserve(maillage.size());

    for (const auto &element : maillage)
        centres.emplace_back(element.milieu.x, element.milieu.y);

    if (centres.empty())
        throw std::invalid_argument("Le maillage est vide.");

    // Construction de la permutation et des échelles.
    auto [permutation, echelles] = construire_ordre_et_echelles(centres);

    // Construction des supports locaux.
    SupportsDistribues supports = construire_supports_sparse_KD(centres, permutation, echelles, rho);

    // Accès aux coefficients de la matrice BEM.
    const AccesMatrice acces_A = [&maillage, &quadrature_maillage](Eigen::Index i, Eigen::Index j) -> Complex
    {
        return matrice_bem(i, j, maillage, quadrature_maillage);
    };

    return construire_facteurs_sparse(acces_A, permutation, supports, tolerance_pivot);
}

Vecteur appliquer_L(const FacteursSparseDistribues &facteurs, const Vecteur &x)
{
    const Eigen::Index N = x.size();

    Vecteur y_local = Vecteur::Zero(N);

    for (Eigen::Index j = facteurs.debut; j < facteurs.fin; ++j)
    {
        for (SparseMatrix::InnerIterator it(facteurs.L_locale, j); it; ++it)
        {
            y_local(it.row()) += it.value() * x(j);
        }
    }

    Vecteur y = Vecteur::Zero(N);

    MPI_Allreduce(y_local.data(), y.data(), static_cast<int>(N), MPI_CXX_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);

    return y;
}

Vecteur appliquer_UT(const FacteursSparseDistribues &facteurs, const Vecteur &x)
{
    const Eigen::Index N = x.size();

    Vecteur y_local = Vecteur::Zero(N);

    // Les lignes de UT sont distribuées entre les processus.
    for (Eigen::Index j = 0; j < facteurs.UT_locale.outerSize(); ++j)
    {
        for (SparseMatrix::InnerIterator it(facteurs.UT_locale, j); it; ++it)
        {
            const Eigen::Index i = it.row();

            if (i >= static_cast<Eigen::Index>(facteurs.debut) && i < static_cast<Eigen::Index>(facteurs.fin))
            {
                y_local(i) += it.value() * x(it.col());
            }
        }
    }

    Vecteur y = Vecteur::Zero(N);

    MPI_Allreduce(y_local.data(), y.data(), static_cast<int>(N), MPI_CXX_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);

    return y;
}

// BLOC

void extraire_blocs_diagonaux_locaux(const Bloc *bloc, std::vector<std::vector<Index>> &clusters)
{
    if (bloc == nullptr)
        return;

    // Bloc terminal.
    if (bloc->type != TypeBloc::Interne)
    {
        if (bloc->est_diagonal &&
            bloc->proprietaire_mpi == mpi_rank)
        {
            if (bloc->cluster_lignes == nullptr)
                throw std::runtime_error("Cluster lignes manquant.");

            clusters.push_back(bloc->cluster_lignes->indices);
        }

        return;
    }

    // Parcours des quatre enfants.
    for (int i = 0; i < 2; ++i)
    {
        for (int j = 0; j < 2; ++j)
        {
            extraire_blocs_diagonaux_locaux(
                bloc->fils[i][j].get(),
                clusters);
        }
    }
}

PreconditionneurBlocs construire_inverses_blocs_diagonaux(const ArbreBlocs &arbre_blocs, const Maillage &maillage, const std::vector<QuadratureSegment> &quadrature_maillage)
{
    if (!arbre_blocs)
        throw std::runtime_error("Arbre de blocs vide.");

    std::vector<std::vector<Index>> clusters;

    extraire_blocs_diagonaux_locaux(
        arbre_blocs.get(),
        clusters);

    PreconditionneurBlocs preconditionneur;
    preconditionneur.reserve(clusters.size());

    for (const auto &indices : clusters)
    {
        const Eigen::Index n =
            static_cast<Eigen::Index>(indices.size());

        if (n == 0)
            continue;

        Matrice D(n, n);

        for (Eigen::Index i = 0; i < n; ++i)
        {
            for (Eigen::Index j = 0; j < n; ++j)
            {
                D(i, j) = matrice_bem(
                    static_cast<Eigen::Index>(indices[i]),
                    static_cast<Eigen::Index>(indices[j]),
                    maillage,
                    quadrature_maillage);
            }
        }

        Eigen::FullPivLU<Matrice> lu(D);

        if (!lu.isInvertible())
        {
            throw std::runtime_error(
                "Un bloc diagonal est singulier.");
        }

        BlocDiagonalInverse bloc;
        bloc.indices = indices;
        bloc.inverse = lu.solve(Matrice::Identity(n, n));

        preconditionneur.push_back(std::move(bloc));
    }

    return preconditionneur;
}
