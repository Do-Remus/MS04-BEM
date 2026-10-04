#ifndef NARVAL_HPP_INCLUDED
#define NARVAL_HPP_INCLUDED

#include <vector>
#include <array>
#include <cmath>
#include <algorithm>
#include "maillage.hpp" // Point, Maillage

/*
 * Frontière 2D d'un narval (vue de dessus, symétrique par rapport à y = 0) :
 * corps fuselé, tête arrondie, défense (tusk) pointant vers +x, nageoires
 * pectorales et lobes de la queue (flukes) vers -x.
 *
 * Construction :
 *   1) demi-contour supérieur défini par points de contrôle (unités : longueur
 *      du corps ~ 2.4, défense incluse ~ 4.1 avant mise à l'échelle) ;
 *   2) symétrie -> contour fermé, orienté anti-horaire ;
 *   3) spline de Catmull-Rom centripète fermée (pas de boucles ni d'overshoot
 *      gênant aux pointes) ;
 *   4) rééchantillonnage à pas constant en abscisse curviligne ;
 *   5) recentrage de la boîte englobante sur l'origine (optionnel).
 *
 * ATTENTION au pas : la défense fait ~0.05*echelle de large, les flukes et la
 * pointe des nageoires ~0.07*echelle. Prendre pas <= largeur/5 environ
 * (0.008 convient pour echelle = 1).
 *
 * Retourne les points SANS répéter le premier à la fin (comme pointsBalleGolf).
 */
inline std::vector<Point> pointsNarval(const double echelle = 1.0,
                                       const double pas = 0.008,
                                       const bool centrer = true)
{
    // Demi-contour supérieur, de la pointe de la défense vers l'encoche de la queue
    // (y > 0 sauf pointe et encoche, situées sur l'axe).
    const std::vector<std::array<double, 2>> haut = {
        {2.40, 0.000}, // pointe de la défense
        {2.30, 0.014},
        {2.00, 0.026},
        {1.40, 0.031},
        {1.14, 0.034}, // sortie de la défense
        {1.06, 0.100}, // museau
        {0.92, 0.190},
        {0.80, 0.250}, // base avant nageoire
        {0.74, 0.330},
        {0.66, 0.410}, // pointe nageoire pectorale
        {0.62, 0.350},
        {0.54, 0.305}, // base arrière nageoire
        {0.40, 0.310}, // largeur max du corps
        {0.00, 0.300},
        {-0.40, 0.250},
        {-0.80, 0.170},
        {-1.10, 0.100},
        {-1.30, 0.060}, // pédoncule
        {-1.42, 0.140},
        {-1.56, 0.300},
        {-1.68, 0.420}, // pointe du lobe de queue
        {-1.64, 0.280},
        {-1.57, 0.120},
        {-1.50, 0.000}  // encoche centrale
    };

    // Contour fermé anti-horaire : haut (droite -> gauche), puis bas (gauche -> droite)
    std::vector<std::array<double, 2>> ctrl = haut;
    for (int i = static_cast<int>(haut.size()) - 2; i >= 1; --i)
        ctrl.push_back({haut[i][0], -haut[i][1]});

    const int m = static_cast<int>(ctrl.size());
    auto C = [&](int i) -> const std::array<double, 2> &
    { return ctrl[((i % m) + m) % m]; };

    // Catmull-Rom centripète (Barry-Goldman) sur le segment P1->P2
    auto interp = [&](int i, double t) -> std::array<double, 2>
    {
        const auto &P0 = C(i - 1), &P1 = C(i), &P2 = C(i + 1), &P3 = C(i + 2);
        auto d = [](const std::array<double, 2> &a, const std::array<double, 2> &b)
        { return std::pow(std::hypot(a[0] - b[0], a[1] - b[1]), 0.5); };
        const double t0 = 0.0, t1 = t0 + d(P0, P1), t2 = t1 + d(P1, P2), t3 = t2 + d(P2, P3);
        const double u = t1 + t * (t2 - t1);
        std::array<double, 2> A1, A2, A3, B1, B2, R;
        for (int c = 0; c < 2; ++c)
        {
            A1[c] = (t1 - u) / (t1 - t0) * P0[c] + (u - t0) / (t1 - t0) * P1[c];
            A2[c] = (t2 - u) / (t2 - t1) * P1[c] + (u - t1) / (t2 - t1) * P2[c];
            A3[c] = (t3 - u) / (t3 - t2) * P2[c] + (u - t2) / (t3 - t2) * P3[c];
            B1[c] = (t2 - u) / (t2 - t0) * A1[c] + (u - t0) / (t2 - t0) * A2[c];
            B2[c] = (t3 - u) / (t3 - t1) * A2[c] + (u - t1) / (t3 - t1) * A3[c];
            R[c] = (t2 - u) / (t2 - t1) * B1[c] + (u - t1) / (t2 - t1) * B2[c];
        }
        return R;
    };

    // Échantillonnage fin de la spline + abscisse curviligne
    const int nSub = 400;
    std::vector<std::array<double, 2>> fin;
    fin.reserve(static_cast<size_t>(m) * nSub + 1);
    for (int i = 0; i < m; ++i)
        for (int j = 0; j < nSub; ++j)
            fin.push_back(interp(i, static_cast<double>(j) / nSub));
    fin.push_back(fin.front()); // fermeture

    std::vector<double> arc(fin.size(), 0.0);
    for (size_t i = 1; i < fin.size(); ++i)
        arc[i] = arc[i - 1] + std::hypot(fin[i][0] - fin[i - 1][0], fin[i][1] - fin[i - 1][1]);
    const double longueur = arc.back();

    // Rééchantillonnage à pas constant
    const unsigned int n = std::max(3u, static_cast<unsigned int>(std::ceil(longueur / (pas / echelle))));
    std::vector<std::array<double, 2>> res;
    res.reserve(n);
    size_t idx = 0;
    for (unsigned int i = 0; i < n; ++i)
    {
        const double s = longueur * i / n;
        while (idx + 2 < arc.size() && arc[idx + 1] < s)
            ++idx;
        const double f = (s - arc[idx]) / (arc[idx + 1] - arc[idx]);
        res.push_back({fin[idx][0] + f * (fin[idx + 1][0] - fin[idx][0]),
                       fin[idx][1] + f * (fin[idx + 1][1] - fin[idx][1])});
    }

    double cx = 0.0;
    if (centrer)
    {
        double xmin = 1e300, xmax = -1e300;
        for (const auto &p : res)
        {
            xmin = std::min(xmin, p[0]);
            xmax = std::max(xmax, p[0]);
        }
        cx = 0.5 * (xmin + xmax);
    }

    std::vector<Point> pts;
    pts.reserve(n);
    for (const auto &p : res)
        pts.emplace_back((p[0] - cx) * echelle, p[1] * echelle);
    return pts;
}

#endif