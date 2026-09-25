#ifndef MAILLAGE_HPP_INCLUDED
#define MAILLAGE_HPP_INCLUDED
#include "../config/config.hpp"
#include "../config/external.hpp"
#include "../config/constantes.hpp"
#include "cercle.hpp"
#include "segment.hpp"

/* Classe Maillage */

class Maillage : public vector<Segment>
{
public:
    Maillage() {};
    Maillage(const vector<Cercle> &cercles, const double pas_maillage);
    Maillage(const vector<Point> &points);
    

    void ajoute_cercle(const double pas_maillage, const Cercle &Ob);
    void export_maillage(const string &filename);
    vector<Point> PointsMaillage();
};

/* Fonctions associées à la classe Maillage */

#endif