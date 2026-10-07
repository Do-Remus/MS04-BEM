#ifndef MAILLAGE_HPP_INCLUDED
#define MAILLAGE_HPP_INCLUDED

#include "../../global/commun.hpp"

#include "cercle.hpp"
#include "segment.hpp"

/* Classe Maillage */

class Maillage : public vector<Segment>
{
public:
    Maillage() {};
    Maillage(const vector<Cercle> &cercles, const Real pas_maillage);
    Maillage(const vector<Point> &points);

    void ajoute_cercle(const Real pas_maillage, const Cercle &Ob);

    void export_maillage(const string &filename) const;

    vector<Point> PointsMaillage();
};

/* Fonctions associées à la class */

std::vector<Point> construire_maillage_solution(const Maillage &maillage);

#endif