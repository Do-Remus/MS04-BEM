#ifndef UTILS_HPP_INCLUDED
#define UTILS_HPP_INCLUDED
#include <time.h>
#include <iostream>
#include "cercle.hpp"
#include "maillage.hpp"

bool obstacle_valide(double rayonAleatoire, double xAleatoire, double yAleatoire, vector<Cercle> obstacles, double hauteur, double epaisseur)
{
    // Verification d'appartenance au rectangle
    if (xAleatoire - rayonAleatoire < 0 or xAleatoire + rayonAleatoire > epaisseur)
    {
        return false;
    }
    if (yAleatoire - rayonAleatoire < 0 or yAleatoire + rayonAleatoire > hauteur)
    {
        return false;
    }

    // Verification de non recouvrement
    for (unsigned int i = 0; i < obstacles.size(); i++)
    {
        double dx = obstacles[i].centre.x - xAleatoire;
        double dy = obstacles[i].centre.y - yAleatoire;
        double distance = sqrt(dx * dx + dy * dy);

        if (distance < obstacles[i].rayon + rayonAleatoire)
        {
            return false;
        }
    }

    return true;
}

Maillage genere_maillage_couche_diffusante(unsigned int nbObstacles, double hauteur, double epaisseur, double pas_maillage)
{
    double rayonLimitant = min(hauteur, epaisseur) / 2;
    double rayonAleatoire;
    double xAleatoire;
    double yAleatoire;
    vector<Cercle> obstacles;
    for (unsigned int i = 0; i < nbObstacles; i++)
    {
        do
        {
            rayonAleatoire = rayonLimitant * (double)rand() / RAND_MAX;
            xAleatoire = epaisseur * (double)rand() / RAND_MAX;
            yAleatoire = hauteur * (double)rand() / RAND_MAX;
        } while (!obstacle_valide(rayonAleatoire, xAleatoire, yAleatoire, obstacles, hauteur, epaisseur));

        Point nouvCentre(xAleatoire, yAleatoire);
        Cercle nouvObstacle(rayonAleatoire, nouvCentre);
        obstacles.push_back(nouvObstacle);
    }

    Maillage maillage(obstacles, pas_maillage);

    return maillage;
}

// int genere_coefficient_vecteur_P(Vecteur &P, const Maillage &maillage)
// {
//     for (int i = 0; i < maillage.size(); i++)
//     {
//         P[i] = integ_simple(maillage[i], pTheta);
//     }

//     return 0;
// }

// int genere_coefficient_matrice_A(Matrice &A, const Maillage &maillage)
// {
//     for (int i = 0; i < maillage.size(); i++)
//     {
//         for (int j = 0; j < maillage.size(); j++)
//         {
//             A(i, j) = integ_double(maillage[i], maillage[j], f) + integ_log(maillage[i], maillage[j]) / (2 * pi);
//         }
//     }

//     return 0;
// }

#endif