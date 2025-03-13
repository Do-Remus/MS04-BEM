#include "headers/utils.hpp"

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

complex<double> green_reguliere(const Point &P1, const Point &P2)
{
    return 0.0;
}

complex<double> p_theta(const Point &P)
{
    return p_theta_config(P.x, P.y);
}

int genere_coefficient_vecteur_P(Vecteur &P, const Maillage &maillage, double pas)
{
    for (unsigned int i = 0; i < maillage.size(); i++)
    {
        Segment seg = maillage[i];
        int nbPas = (int)(seg.norm() / pas);
        P[i] = integ_simple(maillage[i], p_theta, nbPas);
    }

    return 0;
}

int genere_coefficient_matrice_A(MatriceSym &A, const Maillage &maillage, double pas)
{
    for (unsigned int i = 0; i < maillage.size(); i++)
    {
        for (unsigned int j = 0; j < maillage.size(); j++)
        {
            Segment seg1 = maillage[i];
            Segment seg2 = maillage[j];
            int nbPas1 = (int)(seg1.norm() / pas);
            int nbPas2 = (int)(seg2.norm() / pas);
            A(i, j) = integ_double(seg1, seg2, green_reguliere, nbPas1, nbPas2) + integ_log(seg1, seg2, nbPas1, nbPas2) / (2 * pi);
        }
    }

    return 0;
}
