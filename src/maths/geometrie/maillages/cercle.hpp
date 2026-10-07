#ifndef MAILLAGE_CERCLE_HPP_INCLUDED
#define MAILLAGE_CERCLE_HPP_INCLUDED

#include "../../../global/commun.hpp"

#include "../maillage.hpp"

inline Maillage maillage_cercle()
{
    Point O(0, 0);
    Cercle cercle(rayon, O);

    Maillage maillage;
    maillage.ajoute_cercle(pasMaillage, cercle);

    return maillage;
}

#endif