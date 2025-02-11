#ifndef MAILLAGE_HPP_INCLUDED
#define MAILLAGE_HPP_INCLUDED
#include "cercle.hpp"
#include "segment.hpp"
#include <vector>
#include <iostream>
#include <fstream>
using namespace std;

const double pi = atan(1.) * 4;

class Maillage : public vector<Segment>
{
public:
    Maillage() {};
    Maillage(const vector<Cercle> cercles, const double pas_maillage);
    void ajoute_cercle(const double pas_mailage, const Cercle &Ob);
    void export_maillage(const string &filename);
};

void Maillage::ajoute_cercle(const double pas_maillage, const Cercle &Ob)
{
    // Verification pas_maillage > 0
    if (pas_maillage <= 0)
    {
        cout << "ERROR: Le pas du maillage doit-être strictement positif." << endl;
        exit(-1);
    }

    // Construction du cercle
    const int nb_segments = ceil(2 * pi * Ob.rayon / pas_maillage) + 1;
    const double pas_cercle = 2 * pi / nb_segments;
    for (int i = 0; i < nb_segments; i++)
    {
        Point A(Ob.centre.x + Ob.rayon * cos(pas_cercle * i), Ob.centre.y + Ob.rayon * sin(pas_cercle * i));
        Point B(Ob.centre.x + Ob.rayon * cos(pas_cercle * (i + 1)), Ob.centre.y + Ob.rayon * sin(pas_cercle * (i + 1)));
        this->push_back(Segment(A, B));
    }
    return;
}

void Maillage::export_maillage(const string &filename)
{
    ofstream f(filename);

    cout << filename << endl;
    cout << f.is_open() << endl;
    if (!f.is_open())
    {
        cout << "ERROR: Le fichier " << filename << " n'a pas pu être ouvert" << endl;
        exit(-1);
    }

    for (unsigned int i = 0; i < this->size(); i++)
    {
        const Segment &s = this->operator[](i);
        f << s.P1.x << " " << s.P1.y << " ";
        f << s.P2.x << " " << s.P2.y << endl;
    }
    f.close();

    return;
}

Maillage::Maillage(const vector<Cercle> cercles, const double pas_maillage)
{
    // Verification du strict positivité du pas du maillage fait dans ajoute cercle
    for (unsigned int i = 0; i < cercles.size(); i++)
    {
        ajoute_cercle(pas_maillage, cercles[i]);
    }
}

#endif