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
};

/* Fonctions associées à la classe */

std::vector<Point> construire_maillage_solution(const Maillage &maillage);

Arbre construire_arbre_maillage(const Maillage &maillage, const Index N_leafs = 100);

#endif