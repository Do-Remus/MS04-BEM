#ifndef BEM_ARBRE_CLUSTERS_HPP_INCLUDED
#define BEM_ARBRE_CLUSTERS_HPP_INCLUDED

#include "../../global/commun.hpp"
#include "../../maths/commun.hpp"

struct BoiteEnglobante
{
    double x_min;
    double x_max;
    double y_min;
    double y_max;
    double diametre;

    void calculer_diametre()
    {
        diametre = std::hypot(x_max - x_min, y_max - y_min);
    }
};

struct Noeud
{
    std::vector<std::size_t> indices;
    BoiteEnglobante boite;
    std::size_t niveau;

    std::unique_ptr<Noeud> fils_gauche;
    std::unique_ptr<Noeud> fils_droit;

    bool est_feuille() const
    {
        return !fils_gauche && !fils_droit;
    }
};

using Arbre = std::unique_ptr<Noeud>;

using Intervalle = std::pair<Index, Index>;

using IntervallesClusters = std::unordered_map<const Noeud *, Intervalle>;

double distance(const BoiteEnglobante &boite_a, const BoiteEnglobante &boite_b);

void subdivision_arbre_recursive(Noeud &neoud, const Index N_leafs, const std::vector<Point> &points);

Arbre construire_arbre_maillage(const Maillage &maillage, const Index N_leafs = 100);

std::vector<Index> construire_permutation(const Arbre &arbre);

std::vector<Index> construire_permutation_inverse(const std::vector<Index> &permutation);

void construire_intervalles_recursifs(const Noeud &noeud, Index &position, IntervallesClusters &intervalles);

IntervallesClusters construire_intervalles_clusters(const Arbre &arbre);

void ecrire_arbre_clusters(const Arbre &arbre, const std::string &fichier);

#endif
