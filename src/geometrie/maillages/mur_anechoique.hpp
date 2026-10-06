#ifndef MUR_ACOUSTIQUE_HPP_INCLUDED
#define MUR_ACOUSTIQUE_HPP_INCLUDED

#include "../../global/global.hpp"
#include "../maillage.hpp"

/*
 * Frontière 2D d'un mur acoustique : un rectangle dont le côté DROIT est
 * hérissé de pointes triangulaires (type "wedges" d'une chambre anéchoïque),
 * pointant vers +x.
 *
 *            largeur
 *       <--------------> profondeur
 *       +--------------+\
 *       |              | >      <- pointe 3
 *       |              |/ \
 *       |              |   >    <- pointe 2
 *       |              |\ /
 *       |              | >      <- pointe 1
 *       +--------------+/
 *
 * Contour fermé, orienté anti-horaire (normales sortantes), sans répéter le
 * premier point. Les sommets (coins du rectangle, bases et pointes des
 * triangles) sont des points du maillage exacts : chaque arête est découpée en
 * ceil(longueur/pas) segments égaux, donc la longueur des segments est
 * <= pas partout (et proche de pas).
 *
 * Paramètres :
 *   largeur      épaisseur du rectangle (sans les pointes)
 *   hauteur      hauteur totale du mur (les pointes se partagent cette hauteur)
 *   nbPointes    nombre de triangles sur le côté droit
 *   profondeur   longueur d'une pointe (vers +x)
 *   pas          pas de maillage souhaité
 *   centrer      true : boîte englobante (pointes comprises) centrée sur l'origine
 */
inline std::vector<Point> pointsMurAcoustique(const double largeur = 0.5,
                                              const double hauteur = 4.0,
                                              const unsigned int nbPointes = 8,
                                              const double profondeur = 0.6,
                                              const double pas = 0.01,
                                              const bool centrer = true)
{
    const double x0 = 0.0, x1 = largeur; // face gauche, base des pointes
    const double y0 = -0.5 * hauteur, y1 = 0.5 * hauteur;
    const double hp = hauteur / std::max(1u, nbPointes);

    // Sommets du contour, anti-horaire, en partant du coin bas-gauche
    std::vector<Point> sommets;
    sommets.emplace_back(x0, y0); // bas-gauche
    sommets.emplace_back(x1, y0); // bas-droite (base de la 1re pointe)
    for (unsigned int j = 0; j < nbPointes; ++j)
    {
        const double yb = y0 + j * hp;
        sommets.emplace_back(x1 + profondeur, yb + 0.5 * hp); // pointe j
        sommets.emplace_back(x1, y0 + (j + 1) * hp);          // base suivante
    }
    // le dernier sommet ajouté est (x1, y1) = haut-droite
    sommets.emplace_back(x0, y1); // haut-gauche

    // Découpage de chaque arête (sommet i -> sommet i+1)
    std::vector<Point> pts;
    const size_t m = sommets.size();
    for (size_t i = 0; i < m; ++i)
    {
        const Point &A = sommets[i];
        const Point &B = sommets[(i + 1) % m];
        const double len = std::hypot(B.x - A.x, B.y - A.y);
        const unsigned int ns = std::max(1u, static_cast<unsigned int>(std::ceil(len / pas)));
        for (unsigned int s = 0; s < ns; ++s) // B exclu : sera le début de l'arête suivante
        {
            const double t = static_cast<double>(s) / ns;
            pts.emplace_back(A.x + t * (B.x - A.x), A.y + t * (B.y - A.y));
        }
    }

    if (centrer)
    {
        const double cx = 0.5 * (x0 + x1 + profondeur);
        for (Point &p : pts)
            p.x -= cx;
    }
    return pts;
}

#endif