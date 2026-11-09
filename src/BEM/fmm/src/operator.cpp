#include "../operator.hpp"

Index copier_noeud_FMM(const Noeud *noeud, const Index parent, StructureFMM &fmm)
{
    if (noeud == nullptr)
        return INDICE_INVALIDE;

    const Index id = fmm.noeuds.size();

    NoeudFMM nouveau;

    nouveau.parent = parent;
    nouveau.fils_gauche = INDICE_INVALIDE;
    nouveau.fils_droit = INDICE_INVALIDE;

    nouveau.niveau = noeud->niveau;

    // Géométrie
    nouveau.boite = noeud->boite;
    nouveau.boite.calculer_diametre();

    nouveau.centre = Point(Real(0.5) * (nouveau.boite.x_min + nouveau.boite.x_max), Real(0.5) * (nouveau.boite.y_min + nouveau.boite.y_max));
    nouveau.rayon = Real(0.5) * nouveau.boite.diametre;

    // Indices uniquement dans les feuilles
    if (noeud->est_feuille())
        nouveau.indices = noeud->indices;

    // Initialisation des développements
    const Index taille = 2 * fmm.ordre + 1;

    nouveau.multipole = Vecteur::Zero(taille);
    nouveau.local = Vecteur::Zero(taille);

    fmm.noeuds.push_back(std::move(nouveau));
    fmm.nombre_niveaux = std::max(fmm.nombre_niveaux, noeud->niveau + 1);

    // Copie récursive
    const Index gauche = copier_noeud_FMM(noeud->fils_gauche.get(), id, fmm);
    const Index droit = copier_noeud_FMM(noeud->fils_droit.get(), id, fmm);

    fmm.noeuds[id].fils_gauche = gauche;
    fmm.noeuds[id].fils_droit = droit;

    return id;
}

void construire_interactions_FMM(const Index id_A, const Index id_B, StructureFMM &fmm)
{
    // Informations des cellules
    const Index niveau_A = fmm.noeuds[id_A].niveau;
    const Index niveau_B = fmm.noeuds[id_B].niveau;

    const bool feuille_A = fmm.noeuds[id_A].est_feuille();
    const bool feuille_B = fmm.noeuds[id_B].est_feuille();

    const Index gauche_A = fmm.noeuds[id_A].fils_gauche;
    const Index droit_A = fmm.noeuds[id_A].fils_droit;

    const Index gauche_B = fmm.noeuds[id_B].fils_gauche;
    const Index droit_B = fmm.noeuds[id_B].fils_droit;

    // 1. Interaction lointaine :
    // les deux cellules doivent être au même niveau.
    if (niveau_A == niveau_B && admissibles_M2L(fmm.noeuds[id_A], fmm.noeuds[id_B], fmm.eta))
    {
        fmm.interactions_M2L.push_back({id_A, id_B});
        return;
    }

    // 2. Interaction directe entre deux feuilles
    if (feuille_A && feuille_B)
    {
        fmm.interactions_proches.push_back({id_A, id_B});
        return;
    }

    // 3. Niveaux différents :
    // subdiviser prioritairement la cellule la moins profonde.
    if (niveau_A < niveau_B)
    {
        if (!feuille_A)
        {
            construire_interactions_FMM(gauche_A, id_B, fmm);
            construire_interactions_FMM(droit_A, id_B, fmm);
        }
        else
        {
            construire_interactions_FMM(id_A, gauche_B, fmm);
            construire_interactions_FMM(id_A, droit_B, fmm);
        }

        return;
    }

    if (niveau_B < niveau_A)
    {
        if (!feuille_B)
        {
            construire_interactions_FMM(id_A, gauche_B, fmm);
            construire_interactions_FMM(id_A, droit_B, fmm);
        }
        else
        {
            construire_interactions_FMM(gauche_A, id_B, fmm);
            construire_interactions_FMM(droit_A, id_B, fmm);
        }

        return;
    }

    // 4. Même niveau, deux cellules internes :
    // subdiviser les deux cellules.
    if (!feuille_A && !feuille_B)
    {
        construire_interactions_FMM(gauche_A, gauche_B, fmm);
        construire_interactions_FMM(gauche_A, droit_B, fmm);
        construire_interactions_FMM(droit_A, gauche_B, fmm);
        construire_interactions_FMM(droit_A, droit_B, fmm);

        return;
    }

    // 5. Même niveau, seule A est interne
    if (!feuille_A)
    {
        construire_interactions_FMM(gauche_A, id_B, fmm);
        construire_interactions_FMM(droit_A, id_B, fmm);

        return;
    }

    // 6. Même niveau, seule B est interne
    if (!feuille_B)
    {
        construire_interactions_FMM(id_A, gauche_B, fmm);
        construire_interactions_FMM(id_A, droit_B, fmm);
    }
}

void construire_blocs_denses_FMM(StructureFMM &fmm, const Maillage &maillage, const std::vector<QuadratureSegment> &quadratures)
{
    fmm.blocs_denses.clear();

    for (const InteractionFMM &interaction : fmm.interactions_proches)
    {
        const Index id_cible = interaction.cible;
        const Index id_source = interaction.source;

        // On ne stocke qu'une orientation de chaque interaction.
        if (id_cible > id_source)
            continue;

        const NoeudFMM &cible = fmm.noeuds[id_cible];
        const NoeudFMM &source = fmm.noeuds[id_source];

        const Index taille_cible = cible.indices.size();
        const Index taille_source = source.indices.size();

        BlocDenseFMM bloc;

        bloc.cible = id_cible;
        bloc.source = id_source;

        bloc.coefficients.resize(taille_cible, taille_source);

        // Calcul des coefficients du bloc.
        for (Index a = 0; a < taille_cible; ++a)
        {
            const Index i = cible.indices[a];

            for (Index b = 0; b < taille_source; ++b)
            {
                const Index j = source.indices[b];

                bloc.coefficients(a, b) = matrice_bem(i, j, maillage, quadratures);
            }
        }

        fmm.blocs_denses.push_back(std::move(bloc));
    }
}

void precalculer_bases_points_FMM(StructureFMM &fmm, const std::vector<QuadratureSegment> &quadratures)
{
    const int p = static_cast<int>(fmm.ordre);
    const Index N = quadratures.size();

    // Calcul du nombre total de points de quadrature.
    Index nombre_points = 0;

    fmm.offsets_quadratures.resize(N + 1);

    for (Index i = 0; i < N; ++i)
    {
        fmm.offsets_quadratures[i] = nombre_points;
        nombre_points += quadratures[i].size();
    }

    fmm.offsets_quadratures[N] = nombre_points;

    // Allocation des fonctions de base.
    fmm.bases_points.resize(nombre_points, 2 * p + 1);
    fmm.bases_points.setZero();

    // Calcul des fonctions de base pour les feuilles.
    for (const NoeudFMM &noeud : fmm.noeuds)
    {
        if (!noeud.est_feuille())
            continue;

        const Point &centre = noeud.centre;

        for (const Index i : noeud.indices)
        {
            const Index offset = fmm.offsets_quadratures[i];

            for (Index q = 0; q < quadratures[i].size(); ++q)
            {
                const Point &point = quadratures[i][q].point;

                const Real r = (point - centre).norm();
                const Real theta = angle_polaire(point, centre);

                const Real kr = fmm.k * r;

                for (int n = -p; n <= p; ++n)
                {
                    const Index colonne = indice_FMM(n, fmm.ordre);

                    const Real J = std::real(bessel_J(n, kr));

                    fmm.bases_points(offset + q, colonne) =
                        J * phase(n, theta);
                }
            }
        }
    }
}

void precalculer_M2M_FMM(StructureFMM &fmm)
{
    const int p = static_cast<int>(fmm.ordre);
    const Index taille = 4 * fmm.ordre + 1;
    const Real k = fmm.k;

    fmm.translations_M2M.clear();
    fmm.translations_M2M.reserve(2 * fmm.noeuds.size());

    for (Index id = 0; id < fmm.noeuds.size(); ++id)
    {
        const NoeudFMM &parent = fmm.noeuds[id];

        if (parent.est_feuille())
            continue;

        const Index fils[2] = {
            parent.fils_gauche,
            parent.fils_droit};

        for (const Index id_fils : fils)
        {
            if (id_fils == INDICE_INVALIDE)
                continue;

            const NoeudFMM &enfant = fmm.noeuds[id_fils];

            const Point d = enfant.centre - parent.centre;
            const Real r = d.norm();

            TranslationFMM translation;

            translation.cible = id;
            translation.source = id_fils;
            translation.coefficients.resize(taille);
            translation.coefficients.setZero();

            if (r < Real(1e-14))
            {
                // Translation identité.
                translation.coefficients[indice_translation(0, fmm.ordre)] = Complex(1.0, 0.0);
            }
            else
            {
                const Real theta = std::atan2(d.y, d.x);
                const Real kr = k * r;

                for (int q = -2 * p; q <= 2 * p; ++q)
                {
                    translation.coefficients[indice_translation(q, fmm.ordre)] = bessel_J(q, kr) * phase(-q, theta);
                }
            }

            fmm.translations_M2M.push_back(std::move(translation));
        }
    }
}

void precalculer_M2L_FMM(StructureFMM &fmm)
{
    const int p = static_cast<int>(fmm.ordre);
    const Index taille = 4 * fmm.ordre + 1;

    const Real k = fmm.k;
    const Complex facteur_green(0.0, 0.25);

    fmm.translations_M2L.clear();
    fmm.translations_M2L.reserve(fmm.interactions_M2L.size());

    for (const InteractionFMM &interaction : fmm.interactions_M2L)
    {
        const NoeudFMM &cible = fmm.noeuds[interaction.cible];
        const NoeudFMM &source = fmm.noeuds[interaction.source];

        const Point R = cible.centre - source.centre;
        const Real distance = R.norm();

        if (distance <= source.rayon + cible.rayon)
        {
            throw std::runtime_error(
                "Interaction M2L non convergente.");
        }

        const Real theta = std::atan2(R.y, R.x);
        const Real kr = k * distance;

        TranslationFMM translation;

        translation.cible = interaction.cible;
        translation.source = interaction.source;

        translation.coefficients.resize(taille);
        translation.coefficients.setZero();

        for (int q = -2 * p; q <= 2 * p; ++q)
        {
            translation.coefficients[indice_translation(q, fmm.ordre)] =
                facteur_green * hankel_H1(q, kr) * phase(q, theta);
        }

        fmm.translations_M2L.push_back(std::move(translation));
    }
}

void precalculer_L2L_FMM(StructureFMM &fmm)
{
    const int p = static_cast<int>(fmm.ordre);
    const Index taille = 4 * fmm.ordre + 1;

    const Real k = fmm.k;

    fmm.translations_L2L.clear();
    fmm.translations_L2L.reserve(2 * fmm.noeuds.size());

    for (Index id = 0; id < fmm.noeuds.size(); ++id)
    {
        const NoeudFMM &parent = fmm.noeuds[id];

        if (parent.est_feuille())
            continue;

        const Index fils[2] = {
            parent.fils_gauche,
            parent.fils_droit};

        for (const Index id_fils : fils)
        {
            if (id_fils == INDICE_INVALIDE)
                continue;

            const NoeudFMM &enfant = fmm.noeuds[id_fils];

            const Point d = enfant.centre - parent.centre;
            const Real r = d.norm();

            TranslationFMM translation;

            translation.cible = id_fils;
            translation.source = id;

            translation.coefficients.resize(taille);
            translation.coefficients.setZero();

            if (r < Real(1e-14))
            {
                // Translation identité.
                translation.coefficients[indice_translation(0, fmm.ordre)] = Complex(1.0, 0.0);
            }
            else
            {
                const Real theta = std::atan2(d.y, d.x);
                const Real kr = k * r;

                for (int q = -2 * p; q <= 2 * p; ++q)
                {
                    translation.coefficients[indice_translation(q, fmm.ordre)] =
                        bessel_J(q, kr) * phase(q, theta);
                }
            }

            fmm.translations_L2L.push_back(std::move(translation));
        }
    }
}

void precalculer_operateurs_FMM(StructureFMM &fmm, const std::vector<QuadratureSegment> &quadratures)
{
    precalculer_bases_points_FMM(fmm, quadratures);

    precalculer_M2M_FMM(fmm);

    precalculer_M2L_FMM(fmm);

    precalculer_L2L_FMM(fmm);
}

StructureFMM construire_structure_FMM(const Arbre &arbre_clusters, const Real k, const Index ordre, const Maillage &maillage, const std::vector<QuadratureSegment> &quadratures)
{
    StructureFMM fmm;

    fmm.k = k;
    fmm.ordre = ordre;
    fmm.eta = 1.5;

    fmm.nombre_niveaux = 0;
    fmm.racine = INDICE_INVALIDE;

    if (!arbre_clusters)
        return fmm;

    fmm.interactions_proches.clear();
    fmm.interactions_M2L.clear();

    fmm.racine = copier_noeud_FMM(arbre_clusters.get(), INDICE_INVALIDE, fmm);
    construire_interactions_FMM(fmm.racine, fmm.racine, fmm);
    construire_blocs_denses_FMM(fmm, maillage, quadratures);
    precalculer_operateurs_FMM(fmm, quadratures);

    return fmm;
}

void produit_blocs_denses_FMM(const StructureFMM &fmm, const Vecteur &x, Vecteur &y)
{
    for (const BlocDenseFMM &bloc : fmm.blocs_denses)
    {
        const NoeudFMM &cible = fmm.noeuds[bloc.cible];
        const NoeudFMM &source = fmm.noeuds[bloc.source];

        const Index taille_cible = cible.indices.size();
        const Index taille_source = source.indices.size();

        Vecteur x_source(taille_source);

        for (Index b = 0; b < taille_source; ++b)
            x_source(b) = x(source.indices[b]);

        // Interaction diagonale.
        if (bloc.cible == bloc.source)
        {
            Vecteur y_local(taille_cible);

            y_local.noalias() = bloc.coefficients * x_source;

            for (Index a = 0; a < taille_cible; ++a)
                y(cible.indices[a]) += y_local(a);
        }
        // Interaction hors diagonale.
        else
        {
            // Extraction des sous-vecteurs locaux.
            Vecteur x_cible(taille_cible);

            for (Index a = 0; a < taille_cible; ++a)
                x_cible(a) = x(cible.indices[a]);

            Vecteur y_cible(taille_cible);
            Vecteur y_source(taille_source);

            // Contribution A_XY * x_Y
            y_cible.noalias() = bloc.coefficients * x_source;

            // Contribution A_YX * x_X = A_XY^T * x_X
            y_source.noalias() = bloc.coefficients.transpose() * x_cible;

            for (Index a = 0; a < taille_cible; ++a)
                y(cible.indices[a]) += y_cible(a);

            for (Index b = 0; b < taille_source; ++b)
                y(source.indices[b]) += y_source(b);
        }
    }
}

void calculer_P2M(StructureFMM &fmm, const std::vector<QuadratureSegment> &quadratures, const Vecteur &x)
{
    const int p = static_cast<int>(fmm.ordre);

    for (NoeudFMM &noeud : fmm.noeuds)
    {
        if (!noeud.est_feuille())
            continue;

        for (const Index j : noeud.indices)
        {
            const Index offset = fmm.offsets_quadratures[j];

            for (Index q = 0; q < quadratures[j].size(); ++q)
            {
                const Real poids = quadratures[j][q].weight;

                const auto base = fmm.bases_points.row(offset + q);

                for (int n = -p; n <= p; ++n)
                {
                    const Index colonne = indice_FMM(n, fmm.ordre);

                    noeud.multipole(colonne) += x(j) * poids * std::conj(base(colonne));
                }
            }
        }
    }
}

void calculer_M2M(StructureFMM &fmm)
{
    const int p = static_cast<int>(fmm.ordre);

    for (Index id = fmm.translations_M2M.size(); id-- > 0;)
    {
        const TranslationFMM &translation = fmm.translations_M2M[id];

        NoeudFMM &parent = fmm.noeuds[translation.cible];
        const NoeudFMM &enfant = fmm.noeuds[translation.source];

        for (int n = -p; n <= p; ++n)
        {
            const Index indice_n = indice_FMM(n, fmm.ordre);
            Complex somme = Complex(0.0, 0.0);

            for (int m = -p; m <= p; ++m)
            {
                const int q = n - m;

                somme += translation.coefficients[indice_translation(q, fmm.ordre)] * enfant.multipole(indice_FMM(m, fmm.ordre));
            }

            parent.multipole(indice_n) += somme;
        }
    }
}

void calculer_M2L(StructureFMM &fmm)
{
    const int p = static_cast<int>(fmm.ordre);

    for (const TranslationFMM &translation : fmm.translations_M2L)
    {
        NoeudFMM &cible = fmm.noeuds[translation.cible];
        const NoeudFMM &source = fmm.noeuds[translation.source];

        for (int n = -p; n <= p; ++n)
        {
            Complex somme = Complex(0.0, 0.0);

            for (int m = -p; m <= p; ++m)
            {
                const int q = m - n;

                somme += translation.coefficients[indice_translation(q, fmm.ordre)] * source.multipole(indice_FMM(m, fmm.ordre));
            }

            cible.local(indice_FMM(n, fmm.ordre)) += somme;
        }
    }
}

void calculer_L2L(StructureFMM &fmm)
{
    const int p = static_cast<int>(fmm.ordre);

    for (const TranslationFMM &translation : fmm.translations_L2L)
    {
        NoeudFMM &enfant = fmm.noeuds[translation.cible];
        const NoeudFMM &parent = fmm.noeuds[translation.source];

        for (int n = -p; n <= p; ++n)
        {
            Complex somme = Complex(0.0, 0.0);

            for (int m = -p; m <= p; ++m)
            {
                const int q = m - n;

                somme += translation.coefficients[indice_translation(q, fmm.ordre)] * parent.local(indice_FMM(m, fmm.ordre));
            }

            enfant.local(indice_FMM(n, fmm.ordre)) += somme;
        }
    }
}

void calculer_L2P(const StructureFMM &fmm, const std::vector<QuadratureSegment> &quadratures, Vecteur &y)
{
    const int p = static_cast<int>(fmm.ordre);

    for (const NoeudFMM &noeud : fmm.noeuds)
    {
        if (!noeud.est_feuille())
            continue;

        for (const Index i : noeud.indices)
        {
            const Index offset = fmm.offsets_quadratures[i];

            for (Index q = 0; q < quadratures[i].size(); ++q)
            {
                Complex valeur_locale = Complex(0.0, 0.0);

                const auto base = fmm.bases_points.row(offset + q);

                for (int n = -p; n <= p; ++n)
                {
                    const Index colonne = indice_FMM(n, fmm.ordre);

                    valeur_locale += noeud.local(colonne) * base(colonne);
                }

                y(i) += quadratures[i][q].weight * valeur_locale;
            }
        }
    }
}

Vecteur produit_FMM_x(StructureFMM &fmm, const Maillage &maillage, const std::vector<QuadratureSegment> &quadratures, const Vecteur &x)
{
    const Index N = static_cast<Index>(x.size());

    if (quadratures.size() != N)
        throw std::runtime_error(
            "produit_FMM_x : nombre de quadratures incompatible avec x.");

    Vecteur y = Vecteur::Zero(N);

    // 0. Réinitialisation des développements
    reinitialiser_developpements_FMM(fmm);

    // 1. Interactions proches : blocs denses Eigen
    produit_blocs_denses_FMM(fmm, x, y);

    // 2. P2M : construction des multipôles des feuilles
    calculer_P2M(fmm, quadratures, x);

    // 3. M2M : remontée des multipôles
    calculer_M2M(fmm);

    // 4. M2L : interactions lointaines
    calculer_M2L(fmm);

    // 5. L2L : propagation des développements locaux
    calculer_L2L(fmm);

    // 6. L2P : évaluation des développements locaux
    calculer_L2P(fmm, quadratures, y);

    return y;
}
