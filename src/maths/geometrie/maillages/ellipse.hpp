#ifndef MAILLAGE_ELLIPSE_HPP_INCLUDED
#define MAILLAGE_ELLIPSE_HPP_INCLUDED

#include "../../../global/commun.hpp"
#include "../maillage.hpp"

/*
 * ============================================================
 * Frontière 2D d'une ellipse
 * ============================================================
 *
 * L'ellipse est centrée à l'origine et définie par :
 *
 *      x(theta) = a cos(theta)
 *      y(theta) = b sin(theta)
 *
 * où a et b sont les demi-axes horizontal et vertical.
 *
 * Le maillage final est construit directement sur la courbe
 * analytique. Les points sont rééchantillonnés pour obtenir
 * des segments de longueur quasi constante.
 *
 * ============================================================
 */

inline Maillage maillage_ellipse(
    const Real demi_largeur,
    const Real demi_hauteur,
    const Real pas = 0.02)
{
    // === Vérification des paramètres ===

    if (!(demi_largeur > Real(0.0)) ||
        !(demi_hauteur > Real(0.0)) ||
        !(pas > Real(0.0)) ||
        !std::isfinite(demi_largeur) ||
        !std::isfinite(demi_hauteur) ||
        !std::isfinite(pas))
    {
        return Maillage();
    }

    // === Approximation fine de l'ellipse ===

    const unsigned int nLongueur = 10000;

    vector<Point> pointsFins;
    pointsFins.reserve(nLongueur);

    for (unsigned int i = 0; i < nLongueur; ++i)
    {
        const Real theta =
            2.0 * pi * static_cast<Real>(i) /
            static_cast<Real>(nLongueur);

        pointsFins.emplace_back(
            demi_largeur * std::cos(theta),
            demi_hauteur * std::sin(theta));
    }

    // === Calcul de la longueur approchée ===

    Real longueur = 0.0;

    vector<Real> longueurs(nLongueur);

    for (unsigned int i = 0; i < nLongueur; ++i)
    {
        const unsigned int j = (i + 1) % nLongueur;

        longueurs[i] = (pointsFins[i] - pointsFins[j]).norm();
        longueur += longueurs[i];
    }

    if (longueur <= PRECISION_ZERO_DOUBLE)
    {
        return Maillage();
    }

    // === Nombre de segments ===

    unsigned int N =
        static_cast<unsigned int>(std::round(longueur / pas));

    if (N < 3)
    {
        N = 3;
    }

    const Real h = longueur / static_cast<Real>(N);

    // === Rééchantillonnage en longueur d'arc ===

    vector<Point> points;
    points.reserve(N);

    unsigned int j = 0;
    Real cumul = 0.0;

    for (unsigned int i = 0; i < N; ++i)
    {
        const Real cible = static_cast<Real>(i) * h;

        while (j < nLongueur - 1 &&
               cumul + longueurs[j] < cible)
        {
            cumul += longueurs[j];
            ++j;
        }

        const unsigned int j2 = (j + 1) % nLongueur;

        const Real longueurSegment = longueurs[j];

        Real t = 0.0;

        if (longueurSegment > PRECISION_ZERO_DOUBLE)
        {
            t = (cible - cumul) / longueurSegment;
        }

        points.emplace_back(
            pointsFins[j].x +
                t * (pointsFins[j2].x - pointsFins[j].x),

            pointsFins[j].y +
                t * (pointsFins[j2].y - pointsFins[j].y));
    }

    Maillage maillage(points);

    maillage.nom = "ellipse";

    return maillage;
}

#endif