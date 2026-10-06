#ifndef STRUCT_MAILLAGE_HPP_INCLUDED
#define STRUCT_MAILLAGE_HPP_INCLUDED

#include "../global/global.hpp"
#include "../geometrie/cercle.hpp"
#include "../geometrie/segment.hpp"

/* Classe Maillage */

class Maillage : public vector<Segment>
{
public:
    Maillage() {};
    Maillage(const vector<Cercle> &cercles, const Real pas_maillage);
    Maillage(const vector<Point> &points);

    void ajoute_cercle(const Real pas_maillage, const Cercle &Ob);

    void export_maillage(const string &filename);

    vector<Point> PointsMaillage();
};

#endif