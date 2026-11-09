#ifndef MAILLAGE_RECTANGLE_HPP_INCLUDED
#define MAILLAGE_RECTANGLE_HPP_INCLUDED

#include "../../../global/commun.hpp"

#include "../maillage.hpp"

/*
 * ============================================================
 * Frontière 2D d'un rectangle
 * ============================================================
 *
 * Le rectangle est centré à l'origine et possède pour sommets :
 *
 *      (-largeur/2,  hauteur/2)    (largeur/2,  hauteur/2)
 *
 *      (-largeur/2, -hauteur/2)    (largeur/2, -hauteur/2)
 *
 * Les points sont répartis uniformément sur les quatre côtés.
 *
 * ============================================================
 */

inline Maillage maillage_rectangle(
    const Real largeur,
    const Real hauteur,
    const Real pas = 0.02)
{
    // === Vérification des paramètres ===

    if (!(largeur > Real(0.0)) ||
        !(hauteur > Real(0.0)) ||
        !(pas > Real(0.0)) ||
        !std::isfinite(largeur) ||
        !std::isfinite(hauteur) ||
        !std::isfinite(pas))
    {
        return Maillage();
    }

    // === Nombre de segments par dimension ===

    unsigned int N_largeur =
        static_cast<unsigned int>(std::round(largeur / pas));

    unsigned int N_hauteur =
        static_cast<unsigned int>(std::round(hauteur / pas));

    if (N_largeur < 1)
        N_largeur = 1;

    if (N_hauteur < 1)
        N_hauteur = 1;

    const Real h_largeur =
        largeur / static_cast<Real>(N_largeur);

    const Real h_hauteur =
        hauteur / static_cast<Real>(N_hauteur);

    // === Construction des points ===

    const unsigned int N = 2 * (N_largeur + N_hauteur);

    vector<Point> points;
    points.reserve(N);

    const Real demi_largeur = largeur / Real(2.0);
    const Real demi_hauteur = hauteur / Real(2.0);

    // Côté inférieur : gauche -> droite
    for (unsigned int i = 0; i < N_largeur; ++i)
    {
        const Real t = static_cast<Real>(i) * h_largeur;

        points.emplace_back(
            -demi_largeur + t,
            -demi_hauteur);
    }

    // Côté droit : bas -> haut
    for (unsigned int i = 0; i < N_hauteur; ++i)
    {
        const Real t = static_cast<Real>(i) * h_hauteur;

        points.emplace_back(
            demi_largeur,
            -demi_hauteur + t);
    }

    // Côté supérieur : droite -> gauche
    for (unsigned int i = 0; i < N_largeur; ++i)
    {
        const Real t = static_cast<Real>(i) * h_largeur;

        points.emplace_back(
            demi_largeur - t,
            demi_hauteur);
    }

    // Côté gauche : haut -> bas
    for (unsigned int i = 0; i < N_hauteur; ++i)
    {
        const Real t = static_cast<Real>(i) * h_hauteur;

        points.emplace_back(
            -demi_largeur,
            demi_hauteur - t);
    }

    Maillage maillage(points);

    maillage.nom = "rectangle";

    return maillage;
}

#endif