#ifndef MUR_ACOUSTIQUE_HPP_INCLUDED
#define MUR_ACOUSTIQUE_HPP_INCLUDED

#include "../../../global/commun.hpp"
#include "../maillage.hpp"

inline Maillage maillage_mur_acoustique(const double largeur = 0.5, const double hauteur = 4.0, const unsigned int nbPointes = 8, const double profondeur = 0.6, const double pas = 0.01, const bool centrer = true)
{
    // Vérification des paramètres.
    if (largeur <= 0.0 || hauteur <= 0.0 ||
        profondeur <= 0.0 || nbPointes == 0 || pas <= 0.0)
    {
        return Maillage();
    }

    // =========================================================
    // 1. Construction du contour polygonal
    // =========================================================

    vector<Point> sommets;
    const double dy = hauteur / static_cast<double>(nbPointes);

    // Bas du mur, puis côté droit.
    sommets.emplace_back(0.0, 0.0);
    sommets.emplace_back(largeur, 0.0);
    sommets.emplace_back(largeur, hauteur);

    // Haut du mur.
    sommets.emplace_back(0.0, hauteur);

    // Pointes orientées vers la gauche :
    // l'onde incidente arrive de la gauche vers la droite.
    for (unsigned int i = 1; i <= nbPointes; ++i)
    {
        const double yHaut =
            hauteur - static_cast<double>(i - 1) * dy;

        const double yBas =
            hauteur - static_cast<double>(i) * dy;

        const double yMilieu = 0.5 * (yHaut + yBas);

        // Extrémité de la pointe.
        sommets.emplace_back(-profondeur, yMilieu);

        // Base de la pointe suivante.
        // On évite de répéter le sommet inférieur final :
        // le contour sera fermé par le constructeur Maillage.
        if (i < nbPointes)
        {
            sommets.emplace_back(0.0, yBas);
        }
    }

    // =========================================================
    // 2. Centrage horizontal et vertical
    // =========================================================

    if (centrer)
    {
        double xmin = sommets.front().x;
        double xmax = sommets.front().x;
        double ymin = sommets.front().y;
        double ymax = sommets.front().y;

        for (const Point &P : sommets)
        {
            xmin = std::min(xmin, P.x);
            xmax = std::max(xmax, P.x);
            ymin = std::min(ymin, P.y);
            ymax = std::max(ymax, P.y);
        }

        const double centreX = 0.5 * (xmin + xmax);
        const double centreY = 0.5 * (ymin + ymax);

        for (Point &P : sommets)
        {
            P.x -= centreX;
            P.y -= centreY;
        }
    }
    // =========================================================
    // 3. Construction du contour initial
    // =========================================================

    Maillage maillageInitial(sommets);
    const unsigned int M = maillageInitial.size();

    if (M == 0)
    {
        return maillageInitial;
    }

    // Longueur totale du contour.
    double longueur = 0.0;

    for (const Segment &S : maillageInitial)
    {
        longueur += S.norm;
    }

    if (longueur <= PRECISION_ZERO_DOUBLE)
    {
        return maillageInitial;
    }

    // =========================================================
    // 4. Détermination du nombre de segments
    // =========================================================

    unsigned int N = static_cast<unsigned int>(
        std::round(longueur / pas));

    N = std::max(3u, N);

    // Espacement uniforme le long du contour.
    const double h = longueur / static_cast<double>(N);

    // =========================================================
    // 5. Longueurs cumulées
    // =========================================================

    vector<double> arc(M + 1, 0.0);

    for (unsigned int i = 0; i < M; ++i)
    {
        arc[i + 1] = arc[i] + maillageInitial[i].norm;
    }

    // =========================================================
    // 6. Discrétisation uniforme du contour
    // =========================================================

    vector<Point> points;
    points.reserve(N);

    unsigned int j = 0;

    for (unsigned int i = 0; i < N; ++i)
    {
        const double s = static_cast<double>(i) * h;

        while (j + 1 < M && arc[j + 1] < s)
        {
            ++j;
        }

        const Segment &S = maillageInitial[j];
        const double longueurSegment = S.norm;

        double t = 0.0;

        if (longueurSegment > PRECISION_ZERO_DOUBLE)
        {
            t = (s - arc[j]) / longueurSegment;
        }

        points.emplace_back(
            S.P1.x + t * (S.P2.x - S.P1.x),
            S.P1.y + t * (S.P2.y - S.P1.y));
    }

    // =========================================================
    // 7. Maillage final
    // =========================================================

    Maillage maillage(points);
    maillage.nom = "mur_anechoique";

    return maillage;
}

#endif