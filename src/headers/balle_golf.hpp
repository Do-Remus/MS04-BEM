#ifndef BALLE_GOLF_HPP_INCLUDED
#define BALLE_GOLF_HPP_INCLUDED

#include <vector>
#include <cmath>
#include "maillage.hpp" // Point, Maillage

/*
 * Frontière 2D d'une balle de golf (coupe / silhouette), alvéoles amplifiées.
 *
 * Courbe en coordonnées polaires, étoilée par rapport au centre (donc sans
 * auto-intersection) :
 *
 *      r(theta) = R - profondeur * sum_j bump( (theta - theta_j) / demiLargeur )
 *
 * avec theta_j = 2*pi*j/nbAlveoles + decalage et
 *      bump(u) = cos^2(pi*u/2) si |u| < 1, 0 sinon   (classe C^1, support compact).
 *
 * Le pas est ~uniforme en abscisse curviligne : dtheta = pas / sqrt(r^2 + r'^2),
 * puis le nombre de points est ajusté pour fermer exactement la courbe.
 *
 * Retourne des points ordonnés dans le sens trigonométrique, SANS répéter le
 * premier point à la fin : Maillage(points) doit refermer la courbe
 * (dernier point -> premier point), comme le fait déjà votre code P1 avec
 * (i+1) % N.
 *
 * Paramètres :
 *   rayon          R, rayon de la sphère lisse
 *   nbAlveoles     nombre d'alvéoles sur le pourtour
 *   profondeur     profondeur max d'une alvéole (en fraction de R : 0.12 = 12 %)
 *                  -> une vraie balle fait ~0.5 %, ici on amplifie volontairement
 *   remplissage    fraction de l'écart angulaire occupée par une alvéole, dans ]0,1[
 *                  (<1 : bords lisses entre alvéoles ; =1 : alvéoles jointives)
 *   pas            pas de maillage souhaité (longueur des segments)
 */
inline std::vector<Point> pointsBalleGolf(const double rayon,
                                          const unsigned int nbAlveoles = 12,
                                          const double profondeur = 0.12,
                                          const double remplissage = 0.8,
                                          const double pas = 0.02)
{
    const double pi_ = std::acos(-1.0);
    const double ecart = 2.0 * pi_ / nbAlveoles;
    const double demiLargeur = 0.5 * remplissage * ecart; // support de chaque alvéole
    const double decalage = 0.5 * ecart;                  // alvéole centrée entre 2 "crêtes" en theta=0
    const double creux = profondeur * rayon;

    // r(theta) et r'(theta)
    auto rEtDerivee = [&](double theta, double &r, double &dr)
    {
        // on se ramène à l'alvéole la plus proche
        double t = std::fmod(theta - decalage, ecart);
        if (t > 0.5 * ecart)
            t -= ecart;
        if (t < -0.5 * ecart)
            t += ecart;
        r = rayon;
        dr = 0.0;
        if (std::abs(t) < demiLargeur)
        {
            const double u = t / demiLargeur;
            const double c = std::cos(0.5 * pi_ * u);
            r -= creux * c * c;
            // d/dtheta[-creux*cos^2(pi u/2)] = creux*(pi/2)*sin(pi u)/demiLargeur
            dr = creux * (0.5 * pi_ / demiLargeur) * std::sin(pi_ * u);
        }
    };

    // 1) intégration de la longueur d'arc (pour fermer proprement)
    const unsigned int nFin = 200000;
    const double dth = 2.0 * pi_ / nFin;
    std::vector<double> arc(nFin + 1, 0.0);
    for (unsigned int i = 0; i < nFin; ++i)
    {
        double r, dr;
        rEtDerivee((i + 0.5) * dth, r, dr);
        arc[i + 1] = arc[i] + std::sqrt(r * r + dr * dr) * dth;
    }
    const double longueur = arc[nFin];
    const unsigned int n = std::max(3u, static_cast<unsigned int>(std::ceil(longueur / pas)));

    // 2) points équirépartis en abscisse curviligne
    std::vector<Point> pts;
    pts.reserve(n);
    unsigned int idx = 0;
    for (unsigned int i = 0; i < n; ++i)
    {
        const double s = longueur * i / n;
        while (idx + 1 < nFin && arc[idx + 1] < s)
            ++idx;
        const double f = (s - arc[idx]) / (arc[idx + 1] - arc[idx]);
        const double theta = (idx + f) * dth;
        double r, dr;
        rEtDerivee(theta, r, dr);
        pts.emplace_back(r * std::cos(theta), r * std::sin(theta));
    }
    return pts;
}

#endif