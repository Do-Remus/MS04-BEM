#include "../arbre_clusters.hpp"

double distance(const BoiteEnglobante &boite_a, const BoiteEnglobante &boite_b)
{
    double distance_x = std::max({0.0, boite_a.x_min - boite_b.x_max, boite_b.x_min - boite_a.x_max});
    double distance_y = std::max({0.0, boite_a.y_min - boite_b.y_max, boite_b.y_min - boite_a.y_max});
    return std::sqrt(distance_x * distance_x + distance_y * distance_y);
}

void subdivision_arbre_recursive(Noeud &noeud, const Index N_leafs, const std::vector<Point> &points)
{

    if (noeud.indices.size() <= N_leafs)
        return;

    const double largeur = (noeud.boite.x_max - noeud.boite.x_min);
    const double hauteur = (noeud.boite.y_max - noeud.boite.y_min);

    // separation geometrique en 2 groupe divisée selon la plus longue direction
    // trouver quelle est la direction choisie (x ou y) puis definir nouvelles boites englobantes definies comme les moities de la boite englobante du parent
    const bool couper_x = (largeur >= hauteur);
    const double milieu = couper_x ? (noeud.boite.x_min + noeud.boite.x_max) * 0.5 : (noeud.boite.y_min + noeud.boite.y_max) * 0.5;

    // indices temporaires
    std::vector<Index> indices_droit;
    std::vector<Index> indices_gauche;

    indices_droit.reserve(noeud.indices.size());
    indices_gauche.reserve(noeud.indices.size());

    // iteration sur les points indexés dans le parent pour séparation sur les 2 fils
    for (const Index i : noeud.indices)
    {
        const double coord = couper_x ? points[i].x : points[i].y;

        if (coord <= milieu)
            indices_droit.push_back(i);
        else
            indices_gauche.push_back(i);
    }

    if (indices_droit.empty() || indices_gauche.empty())
    {

        if (indices_droit.empty())
        {
            if (couper_x)
                noeud.boite.x_min = milieu;
            else
                noeud.boite.y_min = milieu;
        }
        else
        {
            if (couper_x)
                noeud.boite.x_max = milieu;
            else
                noeud.boite.y_max = milieu;
        }

        noeud.boite.calculer_diametre();

        // Recommencer sur le même parent avec sa nouvelle boîte.
        subdivision_arbre_recursive(noeud, N_leafs, points);
        return;
    }

    // La coupure est valide : créer les deux fils.
    noeud.fils_droit = std::make_unique<Noeud>();
    noeud.fils_gauche = std::make_unique<Noeud>();

    Noeud &droit = *noeud.fils_droit;
    Noeud &gauche = *noeud.fils_gauche;

    droit.niveau = noeud.niveau + 1;
    gauche.niveau = noeud.niveau + 1;

    droit.indices = std::move(indices_droit);
    gauche.indices = std::move(indices_gauche);

    // Boîtes géométriques issues de la coupure.
    droit.boite.x_min = noeud.boite.x_min;
    droit.boite.y_min = noeud.boite.y_min;
    droit.boite.x_max = couper_x ? milieu : noeud.boite.x_max;
    droit.boite.y_max = couper_x ? noeud.boite.y_max : milieu;

    gauche.boite.x_min = couper_x ? milieu : noeud.boite.x_min;
    gauche.boite.y_min = couper_x ? noeud.boite.y_min : milieu;
    gauche.boite.x_max = noeud.boite.x_max;
    gauche.boite.y_max = noeud.boite.y_max;

    droit.boite.calculer_diametre();
    gauche.boite.calculer_diametre();

    subdivision_arbre_recursive(*noeud.fils_gauche, N_leafs, points);
    subdivision_arbre_recursive(*noeud.fils_droit, N_leafs, points);
}

Arbre construire_arbre_maillage(const Maillage &maillage, const Index N_leafs)
{
    Arbre arbre = std::make_unique<Noeud>();
    // -- Initialisation du 1er noeud --

    arbre->niveau = 0;

    const std::vector<Point> &points = maillage.PointsMilieuMaillage();
    arbre->indices.resize(points.size());

    arbre->boite.x_min = std::numeric_limits<double>::infinity();
    arbre->boite.x_max = -std::numeric_limits<double>::infinity();

    arbre->boite.y_min = std::numeric_limits<double>::infinity();
    arbre->boite.y_max = -std::numeric_limits<double>::infinity();

    for (Index i = 0; i < points.size(); ++i)
    {
        arbre->indices[i] = i;

        const Point &point = points[i];

        arbre->boite.x_min = std::min(arbre->boite.x_min, point.x);
        arbre->boite.x_max = std::max(arbre->boite.x_max, point.x);

        arbre->boite.y_min = std::min(arbre->boite.y_min, point.y);
        arbre->boite.y_max = std::max(arbre->boite.y_max, point.y);
    }

    arbre->boite.calculer_diametre();

    // -- Construction récursive de l'arbre --
    subdivision_arbre_recursive(*arbre, N_leafs, points);

    return arbre;
}

std::vector<Index> construire_permutation(const Arbre &arbre)
{
    std::vector<Index> permutation;

    if (!arbre)
        return permutation;

    // Parcours en profondeur : gauche puis droite.
    auto visiter = [&](auto &&self, const Noeud &noeud) -> void
    {
        if (noeud.est_feuille())
        {
            for (const Index i : noeud.indices)
                permutation.push_back(i);

            return;
        }

        if (noeud.fils_gauche)
            self(self, *noeud.fils_gauche);

        if (noeud.fils_droit)
            self(self, *noeud.fils_droit);
    };

    visiter(visiter, *arbre);

    return permutation;
}

std::vector<Index> construire_permutation_inverse(const std::vector<Index> &permutation)
{
    std::vector<Index> inverse(permutation.size());
    std::vector<bool> vu(permutation.size(), false);

    for (Index k = 0; k < permutation.size(); ++k)
    {
        const Index i = permutation[k];

        if (i >= permutation.size() || vu[i])
            throw std::runtime_error(
                "La permutation contient un indice invalide ou un doublon.");

        inverse[i] = k;
        vu[i] = true;
    }

    return inverse;
}

void construire_intervalles_recursifs(const Noeud &noeud, Index &position, IntervallesClusters &intervalles)
{
    const Index debut = position;

    if (noeud.est_feuille())
    {
        position += noeud.indices.size();
        intervalles[&noeud] = {debut, position};
        return;
    }

    if (noeud.fils_gauche)
    {
        construire_intervalles_recursifs(
            *noeud.fils_gauche, position, intervalles);
    }

    if (noeud.fils_droit)
    {
        construire_intervalles_recursifs(
            *noeud.fils_droit, position, intervalles);
    }

    intervalles[&noeud] = {debut, position};
}

IntervallesClusters construire_intervalles_clusters(const Arbre &arbre)
{
    IntervallesClusters intervalles;

    if (!arbre)
        return intervalles;

    Index position = 0;

    construire_intervalles_recursifs(*arbre, position, intervalles);

    return intervalles;
}

void ecrire_arbre_clusters(const Arbre &arbre, const std::string &fichier)
{
    std::ofstream out(fichier);

    if (!out)
        throw std::runtime_error("Impossible d'ouvrir " + fichier);

    out << "id,parent,niveau,taille,x_min,x_max,y_min,y_max,diametre\n";
    out << std::setprecision(17);

    if (!arbre)
        return;

    std::size_t prochain_id = 0;

    auto visiter = [&](auto &&self, const Noeud &noeud, long long parent_id) -> void
    {
        const std::size_t id = prochain_id++;

        out << id << ','
            << parent_id << ','
            << noeud.niveau << ','
            << noeud.indices.size() << ','
            << noeud.boite.x_min << ','
            << noeud.boite.x_max << ','
            << noeud.boite.y_min << ','
            << noeud.boite.y_max << ','
            << noeud.boite.diametre << '\n';

        if (noeud.fils_gauche)
            self(self, *noeud.fils_gauche, static_cast<long long>(id));

        if (noeud.fils_droit)
            self(self, *noeud.fils_droit, static_cast<long long>(id));
    };

    visiter(visiter, *arbre, -1);
}
