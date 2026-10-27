#ifndef NARVAL_HPP_INCLUDED
#define NARVAL_HPP_INCLUDED

#include "../../../global/commun.hpp"

#include "../maillage.hpp"

/*
 * ============================================================
 * Frontière 2D d'un narval
 * ============================================================
 *
 * Géométrie construite avec des courbes de Bézier cubiques.
 * Narval orienté vers la gauche.
 *
 * - Corps arrondi et fuselé
 * - Défense longue et fine
 * - Queue bilobée
 * - Contour fermé et lisse
 *
 * Le contour est rééchantillonné uniformément en longueur d'arc.
 * ============================================================
 */
inline Maillage maillage_narval(
    const double echelle = 1.0,
    const double pas = 0.008,
    const bool centrer = true)
{
    if (!std::isfinite(echelle) ||
        !std::isfinite(pas) ||
        echelle <= 0.0 ||
        pas <= 0.0)
    {
        return Maillage();
    }

    // --------------------------------------------------------
    // 1. Outils géométriques
    // --------------------------------------------------------
    auto bezier = [](const Point &P0,
                     const Point &P1,
                     const Point &P2,
                     const Point &P3,
                     const Real t) -> Point
    {
        const Real s = 1.0 - t;
        const Real s2 = s * s;
        const Real t2 = t * t;

        const Real b0 = s2 * s;
        const Real b1 = 3.0 * s2 * t;
        const Real b2 = 3.0 * s * t2;
        const Real b3 = t2 * t;

        return Point(
            b0 * P0.x + b1 * P1.x + b2 * P2.x + b3 * P3.x,
            b0 * P0.y + b1 * P1.y + b2 * P2.y + b3 * P3.y);
    };

    // --------------------------------------------------------
    // 2. Construction du contour dans le sens antihoraire
    // --------------------------------------------------------
    //
    // Chaque courbe est définie par 4 points de Bézier.
    // La fin d'une courbe est le début de la suivante.
    //
    vector<Point> contour;
    const unsigned int echantillons = 300;

    auto ajouter_bezier = [&](const Point &P0,
                              const Point &P1,
                              const Point &P2,
                              const Point &P3)
    {
        for (unsigned int i = 0; i < echantillons; ++i)
        {
            const Real t =
                static_cast<Real>(i) /
                static_cast<Real>(echantillons);

            contour.push_back(bezier(P0, P1, P2, P3, t));
        }
    };

    /*
     * Le contour commence à la pointe de la défense.
     *
     * Coordonnées :
     * x négatif = avant du narval
     * x positif = arrière et queue
     */

    // A. Défense : bord supérieur, de la pointe vers la tête
    ajouter_bezier(
        Point(-2.20, 0.00),
        Point(-1.75, 0.005),
        Point(-1.15, 0.035),
        Point(-0.70, 0.075));

    // B. Front et dos
    ajouter_bezier(
        Point(-0.70, 0.11),
        Point(-0.60, 0.31),
        Point(-0.30, 0.43),
        Point(0.08, 0.43));

    // C. Dos vers l'arrière du corps
    ajouter_bezier(
        Point(0.08, 0.43),
        Point(0.55, 0.43),
        Point(0.95, 0.24),
        Point(1.23, 0.12));

    // D. Transition vers le lobe supérieur de la queue
    ajouter_bezier(
        Point(1.23, 0.12),
        Point(1.34, 0.10),
        Point(1.38, 0.24),
        Point(1.54, 0.34));

    // E. Extrémité arrondie du lobe supérieur
    ajouter_bezier(
        Point(1.54, 0.34),
        Point(1.58, 0.24),
        Point(1.51, 0.08),
        Point(1.41, 0.015));

    // F. Encoche centrale de la queue
    ajouter_bezier(
        Point(1.41, 0.015),
        Point(1.37, 0.00),
        Point(1.37, -0.015),
        Point(1.41, -0.035));

    // G. Bord extérieur du lobe inférieur
    ajouter_bezier(
        Point(1.41, -0.035),
        Point(1.52, -0.12),
        Point(1.59, -0.27),
        Point(1.54, -0.34));

    // H. Retour de la queue vers le ventre
    ajouter_bezier(
        Point(1.54, -0.34),
        Point(1.38, -0.25),
        Point(1.33, -0.11),
        Point(1.22, -0.12));

    // I. Ventre du corps
    ajouter_bezier(
        Point(1.22, -0.12),
        Point(0.92, -0.27),
        Point(0.50, -0.40),
        Point(0.08, -0.40));

    // J. Ventre et menton
    ajouter_bezier(
        Point(0.08, -0.40),
        Point(-0.28, -0.40),
        Point(-0.56, -0.29),
        Point(-0.70, -0.11));

    // K. Défense : bord inférieur vers la pointe
    ajouter_bezier(
        Point(-0.70, -0.075),
        Point(-1.15, -0.035),
        Point(-1.75, -0.005),
        Point(-2.20, 0.00));

    // --------------------------------------------------------
    // 3. Mise à l'échelle
    // --------------------------------------------------------
    for (Point &P : contour)
    {
        P.x *= echelle;
        P.y *= echelle;
    }

    // --------------------------------------------------------
    // 4. Centrage de la boîte englobante
    // --------------------------------------------------------
    if (centrer)
    {
        Real xmin = contour[0].x;
        Real xmax = contour[0].x;
        Real ymin = contour[0].y;
        Real ymax = contour[0].y;

        for (const Point &P : contour)
        {
            xmin = std::min(xmin, P.x);
            xmax = std::max(xmax, P.x);
            ymin = std::min(ymin, P.y);
            ymax = std::max(ymax, P.y);
        }

        const Real cx = 0.5 * (xmin + xmax);
        const Real cy = 0.5 * (ymin + ymax);

        for (Point &P : contour)
        {
            P.x -= cx;
            P.y -= cy;
        }
    }

    // --------------------------------------------------------
    // 5. Longueurs cumulées du contour fermé
    // --------------------------------------------------------
    const unsigned int M =
        static_cast<unsigned int>(contour.size());

    if (M < 3)
    {
        return Maillage();
    }

    vector<Real> longueurs(M);
    vector<Real> cumul(M + 1, 0.0);

    for (unsigned int i = 0; i < M; ++i)
    {
        const unsigned int j = (i + 1) % M;

        const Real dx = contour[j].x - contour[i].x;
        const Real dy = contour[j].y - contour[i].y;

        longueurs[i] = std::sqrt(dx * dx + dy * dy);
        cumul[i + 1] = cumul[i] + longueurs[i];
    }

    const Real longueurTotale = cumul[M];

    if (!std::isfinite(longueurTotale) ||
        longueurTotale <= PRECISION_ZERO_DOUBLE)
    {
        return Maillage();
    }

    // --------------------------------------------------------
    // 6. Nombre de segments du maillage final
    // --------------------------------------------------------
    const Real nombreReel = longueurTotale / pas;

    if (!std::isfinite(nombreReel) ||
        nombreReel > static_cast<Real>(
                         std::numeric_limits<unsigned int>::max()))
    {
        return Maillage();
    }

    unsigned int N = static_cast<unsigned int>(
        std::round(nombreReel));

    N = std::max(3u, N);

    const Real h = longueurTotale / static_cast<Real>(N);

    // --------------------------------------------------------
    // 7. Rééchantillonnage uniforme en longueur d'arc
    // --------------------------------------------------------
    vector<Point> points;
    points.reserve(N);

    unsigned int j = 0;

    for (unsigned int i = 0; i < N; ++i)
    {
        const Real cible = static_cast<Real>(i) * h;

        while (j + 1 < M && cumul[j + 1] < cible)
        {
            ++j;
        }

        const unsigned int j2 = (j + 1) % M;
        const Real L = longueurs[j];

        Real t = 0.0;

        if (L > PRECISION_ZERO_DOUBLE)
        {
            t = (cible - cumul[j]) / L;
        }

        t = std::max(0.0, std::min(1.0, t));

        points.emplace_back(
            contour[j].x + t * (contour[j2].x - contour[j].x),
            contour[j].y + t * (contour[j2].y - contour[j].y));
    }

    return Maillage(points);
}

#endif