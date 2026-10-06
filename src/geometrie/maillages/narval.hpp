#ifndef NARVAL_HPP_INCLUDED
#define NARVAL_HPP_INCLUDED

#include "../../global/global.hpp"
#include "../maillage.hpp"

/*
 * ============================================================
 * Frontière 2D d'un narval
 * ============================================================
 *
 * La géométrie est définie par une spline de Catmull-Rom
 * centripète fermée.
 *
 * Le maillage final est construit directement sur la spline.
 *
 * ============================================================
 */
inline Maillage maillageNarval(const double echelle = 1.0, const double pas = 0.008, const bool centrer = true)
{
    /* === Points de contrôle === */
    vector<Point> controle = {
        Point(-1.00, 0.00),
        Point(-0.85, 0.45),
        Point(-0.45, 0.65),
        Point(0.00, 0.55),
        Point(0.45, 0.70),
        Point(0.85, 0.40),
        Point(1.00, 0.00),
        Point(0.85, -0.40),
        Point(0.45, -0.70),
        Point(0.00, -0.55),
        Point(-0.45, -0.65),
        Point(-0.85, -0.45)};

    /* === Vérification des paramètres === */
    if (echelle <= 0.0 || pas <= 0.0 || controle.size() < 4)
    {
        return Maillage();
    }

    /* === Spline Catmull-Rom === */
    auto catmullRom = [](const Point &P0, const Point &P1, const Point &P2, const Point &P3, const Real t) -> Point
    {
        const Real t2 = t * t;
        const Real t3 = t2 * t;

        const Real c0 = -0.5 * t3 + t2 - 0.5 * t;
        const Real c1 = 1.5 * t3 - 2.5 * t2 + 1.0;
        const Real c2 = -1.5 * t3 + 2.0 * t2 + 0.5 * t;
        const Real c3 = 0.5 * t3 - 0.5 * t2;

        return Point(
            c0 * P0.x + c1 * P1.x + c2 * P2.x + c3 * P3.x,
            c0 * P0.y + c1 * P1.y + c2 * P2.y + c3 * P3.y);
    };

    /* === Construction de la courbe fine === */
    const unsigned int nbControle = controle.size();

    const unsigned int nParSegment = 1000;

    vector<Point> pointsFins;

    pointsFins.reserve(nbControle * nParSegment);

    for (unsigned int i = 0; i < nbControle; ++i)
    {
        const unsigned int i0 = (i + nbControle - 1) % nbControle;
        const unsigned int i1 = i;
        const unsigned int i2 = (i + 1) % nbControle;
        const unsigned int i3 = (i + 2) % nbControle;

        for (unsigned int j = 0; j < nParSegment; ++j)
        {
            const Real t = static_cast<Real>(j) / static_cast<Real>(nParSegment);

            pointsFins.push_back(
                catmullRom(
                    controle[i0],
                    controle[i1],
                    controle[i2],
                    controle[i3],
                    t));
        }
    }

    /* === Mise à l'échelle === */
    for (Point &P : pointsFins)
    {
        P.x *= echelle;
        P.y *= echelle;
    }

    /* === Centrage horizontal === */
    if (centrer)
    {
        Real xmin = pointsFins[0].x;
        Real xmax = pointsFins[0].x;

        for (const Point &P : pointsFins)
        {
            xmin = std::min(xmin, P.x);
            xmax = std::max(xmax, P.x);
        }

        const Real centreX = 0.5 * (xmin + xmax);

        for (Point &P : pointsFins)
        {
            P.x -= centreX;
        }
    }

    /* === Longueur de la courbe === */
    const unsigned int M = pointsFins.size();

    Real longueur = 0.0;

    vector<Real> longueurs(M);

    for (unsigned int i = 0; i < M; ++i)
    {
        const unsigned int j = (i + 1) % M;

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

        while (j < M - 1 && cumul + longueurs[j] < cible)
        {
            cumul += longueurs[j];
            ++j;
        }

        const unsigned int j2 = (j + 1) % M;

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