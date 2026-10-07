#ifndef MUR_ACOUSTIQUE_HPP_INCLUDED
#define MUR_ACOUSTIQUE_HPP_INCLUDED

#include "../../../global/commun.hpp"

#include "../maillage.hpp"

/*
 * ============================================================
 * Frontière 2D d'un mur acoustique
 * ============================================================
 *
 * Le contour est un polygone exact dont le côté droit possède
 * des pointes triangulaires.
 *
 * Le maillage est construit directement à partir des sommets.
 *
 * ============================================================
 */
inline Maillage maillage_mur_acoustique(const double largeur = 0.5, const double hauteur = 4.0, const unsigned int nbPointes = 8, const double profondeur = 0.6, const double pas = 0.01, const bool centrer = true)
{
    /* === Vérification des paramètres === */
    if (largeur <= 0.0 || hauteur <= 0.0 || profondeur <= 0.0 || nbPointes == 0 || pas <= 0.0)
    {
        return Maillage();
    }

    /* === Construction du contour === */
    vector<Point> sommets;

    const double dy = hauteur / static_cast<double>(nbPointes);

    // Côté gauche
    sommets.emplace_back(0.0, 0.0);

    // Bas
    sommets.emplace_back(largeur, 0.0);

    /* === Pointes du mur === */
    for (unsigned int i = 0; i < nbPointes; ++i)
    {
        const double y0 = static_cast<double>(i) * dy;
        const double y1 = static_cast<double>(i + 1) * dy;
        const double ym = 0.5 * (y0 + y1);

        // Base de la pointe
        sommets.emplace_back(largeur, y0);

        // Extrémité de la pointe
        sommets.emplace_back(largeur + profondeur, ym);

        // Retour sur la base
        sommets.emplace_back(largeur, y1);
    }

    // Haut
    sommets.emplace_back(0.0, hauteur);

    /* === Centrage === */
    if (centrer)
    {
        Real xmin = sommets[0].x;
        Real xmax = sommets[0].x;

        for (const Point &P : sommets)
        {
            xmin = std::min(xmin, P.x);
            xmax = std::max(xmax, P.x);
        }

        const Real centreX = 0.5 * (xmin + xmax);

        for (Point &P : sommets)
        {
            P.x -= centreX;
        }
    }

    /* === Maillage initial === */
    Maillage maillageInitial(sommets);

    const unsigned int M = maillageInitial.size();

    if (M == 0)
    {
        return maillageInitial;
    }

    /* === Longueur totale === */
    Real longueur = 0.0;

    for (const Segment &S : maillageInitial)
    {
        longueur += S.norm;
    }

    if (longueur <= PRECISION_ZERO_DOUBLE)
    {
        return maillageInitial;
    }

    /* === Nombre de segments === */
    unsigned int N = static_cast<unsigned int>(std::round(longueur / pas));

    if (N < 3)
    {
        N = 3;
    }

    const Real h = longueur / static_cast<Real>(N);

    /* === Longueurs cumulées === */
    vector<Real> arc(M + 1, 0.0);

    for (unsigned int i = 0; i < M; ++i)
    {
        arc[i + 1] = arc[i] + maillageInitial[i].norm;
    }

    /* === Placement des points === */
    vector<Point> points;
    points.reserve(N);

    unsigned int j = 0;

    for (unsigned int i = 0; i < N; ++i)
    {
        const Real s = static_cast<Real>(i) * h;

        while (j + 1 < M && arc[j + 1] < s)
        {
            ++j;
        }

        const Segment &S = maillageInitial[j];

        const Real longueurSegment = S.norm;

        Real t = 0.0;

        if (longueurSegment > PRECISION_ZERO_DOUBLE)
        {
            t = (s - arc[j]) / longueurSegment;
        }

        points.emplace_back(
            S.P1.x + t * (S.P2.x - S.P1.x),
            S.P1.y + t * (S.P2.y - S.P1.y));
    }

    return Maillage(points);
}

#endif