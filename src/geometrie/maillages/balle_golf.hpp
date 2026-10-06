#ifndef BALLE_GOLF_HPP_INCLUDED
#define BALLE_GOLF_HPP_INCLUDED

#include "../../global/global.hpp"
#include "../maillage.hpp"

/*
 * ============================================================
 * Frontière 2D d'une balle de golf
 * ============================================================
 *
 * La frontière est définie en coordonnées polaires par
 *
 *      r(theta) = R - profondeur * creux(theta)
 *
 * avec des alvéoles régulièrement réparties.
 *
 * Le maillage final est construit directement sur la courbe
 * analytique. Les points sont choisis de manière à obtenir
 * des segments de longueur euclidienne quasi constante.
 *
 * ============================================================
 */
inline Maillage maillageBalleGolf(const double rayon, const unsigned int nbAlveoles = 12, const double profondeur = 0.12, const double remplissage = 0.8, const double pas = 0.02)
{
    /* === Vérification des paramètres === */
    if (rayon <= 0.0 || nbAlveoles == 0 || profondeur < 0.0 || pas <= 0.0)
    {
        return Maillage();
    }

    const double remplissageEffectif = std::max(0.0, std::min(1.0, remplissage));

    /* === Fonction de creux === */
    auto creux = [&](const Real theta) -> Real
    {
        const Real periode = 2.0 * pi / static_cast<Real>(nbAlveoles);

        Real thetaLocal = std::fmod(theta, periode);

        if (thetaLocal < 0.0)
        {
            thetaLocal += periode;
        }

        const Real centre = 0.5 * periode;
        const Real demiLargeur = 0.5 * remplissageEffectif * periode;

        const Real distanceCentre = std::abs(thetaLocal - centre);

        if (distanceCentre >= demiLargeur)
        {
            return 0.0;
        }

        const Real x = distanceCentre / demiLargeur;

        return 0.5 * (1.0 + std::cos(pi * x));
    };

    /* === Courbe analytique === */
    auto rayonTheta = [&](const Real theta) -> Real
    {
        return rayon - profondeur * creux(theta);
    };

    /* === Approximation fine de la courbe === */
    const unsigned int nLongueur = 10000;

    vector<Point> pointsFins;
    pointsFins.reserve(nLongueur);

    for (unsigned int i = 0; i < nLongueur; ++i)
    {
        const Real theta = 2.0 * pi * static_cast<Real>(i) / static_cast<Real>(nLongueur);

        const Real r = rayonTheta(theta);

        pointsFins.emplace_back(
            r * std::cos(theta),
            r * std::sin(theta));
    }

    /* === Longueur approchée === */
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

    /* === Nombre de segments === */
    unsigned int N = static_cast<unsigned int>(std::round(longueur / pas));

    if (N < 3)
    {
        N = 3;
    }

    const Real h = longueur / static_cast<Real>(N);

    /* === Rééchantillonnage en longueur d'arc === */
    vector<Point> points;
    points.reserve(N);

    unsigned int j = 0;
    Real cumul = 0.0;

    for (unsigned int i = 0; i < N; ++i)
    {
        const Real cible = static_cast<Real>(i) * h;

        while (j < nLongueur - 1 && cumul + longueurs[j] < cible)
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
            pointsFins[j].x + t * (pointsFins[j2].x - pointsFins[j].x),
            pointsFins[j].y + t * (pointsFins[j2].y - pointsFins[j].y));
    }

    return Maillage(points);
}

#endif