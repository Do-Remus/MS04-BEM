#ifndef MAILLAGE_HPP_INCLUDED
#define MAILLAGE_HPP_INCLUDED

#include "../../global/commun.hpp"

#include "cercle.hpp"
#include "segment.hpp"

/* Structure pour la structuration du maillage */

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

enum class TypeBloc
{
    Interne,    // Possède quatre blocs fils
    Admissible, // Bloc admissible : compression ACA possible
    Dense,      // Bloc non admissible : stockage dense
    Vide        // Bloc avec aucun point
};

struct Bloc
{
    const Noeud *cluster_lignes;
    const Noeud *cluster_colonnes;

    TypeBloc type = TypeBloc::Interne;

    std::unique_ptr<Bloc> fils[2][2];
};

using ArbreBlocs = std::unique_ptr<Bloc>;

/* Classe Maillage */

class Maillage : public vector<Segment>
{
public:
    Maillage() {};
    Maillage(const vector<Cercle> &cercles, const Real pas_maillage);
    Maillage(const vector<Point> &points);

    void ajoute_cercle(const Real pas_maillage, const Cercle &Ob);

    void export_maillage(const string &filename) const;

    vector<Point> PointsMaillage() const;
    vector<Point> PointsMilieuMaillage() const;
};

/* Fonctions associées à la classe */

double distance(const BoiteEnglobante &boite_a, const BoiteEnglobante &boite_b);

bool sont_admissibles(const BoiteEnglobante &boite_a, const BoiteEnglobante &boite_b, const double eta);

std::vector<Point> construire_maillage_solution(const Maillage &maillage);

void subdivision_arbre_recursive(Noeud &neoud, const Index N_leafs, const std::vector<Point> &points);

Arbre construire_arbre_maillage(const Maillage &maillage, const Index N_leafs = 100);

ArbreBlocs construire_bloc_recursif(const Noeud &tau, const Noeud &sigma, const double eta);

std::vector<Index> construire_permutation(const Arbre &arbre);

std::vector<Index> construire_permutation_inverse(const std::vector<Index> &permutation);

void construire_intervalles_recursifs(const Noeud &noeud, Index &position, IntervallesClusters &intervalles);

IntervallesClusters construire_intervalles_clusters(const Arbre &arbre);

void ecrire_arbre_clusters(const Arbre &arbre, const std::string &fichier);

ArbreBlocs construire_arbre_blocs(const Arbre &arbre_clusters, const double eta);

void ecrire_arbre_blocs(const ArbreBlocs &arbre, const std::string &fichier, const IntervallesClusters &intervalles);

#endif