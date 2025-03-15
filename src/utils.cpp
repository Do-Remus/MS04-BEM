#include "headers/utils.hpp"

bool obstacle_valide(double rayonAleatoire, double xAleatoire, double yAleatoire, vector<Cercle> &obstacles, double hauteur, double epaisseur)
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

Maillage genere_maillage_couche_diffusante(unsigned int nbObstacles, double hauteur, double epaisseur, double pas_maillage, vector<Cercle> &obstacles)
{
    double rayonLimitant = min(hauteur, epaisseur) / 2;
    double rayonAleatoire;
    double xAleatoire;
    double yAleatoire;
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

complex<double> hankel(const Point &P1)
{
    Point P = P1;
    double x = k * P.norm();
    return j0(x) + I * y0(x);
}

complex<double> green_reguliere(const Point &P1, const Point &P2)
{
    double x = k * (P1 - P2).norm();
    if (x == 0)
    {
        return 1. + 2. * (I / pi) * (gamma_euler - log(2));
    }
    return j0(x) + I * y0(x) - (2. * I / pi) * log(x);
}

complex<double> hankel_derivate(const Point &P1)
{
    Point P = P1;
    double x = k * P.norm();
    return -k * (j1(x) + I * y1(x));
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
            A(i, j) = (1 / 4. * I) * integ_double(seg1, seg2, green_reguliere, nbPas1, nbPas2) + integ_log(seg1, seg2, nbPas1, nbPas2) / (2 * pi);
        }
    }

    return 0;
}

complex<double> p(const Maillage &maillage, double pas, const Vecteur &Q, const Vecteur P, const Point &x)
{
    complex<double> result = 0;
    for (unsigned int i = 0; i < maillage.size(); i++)
    {
        Segment seg = maillage[i];
        int nbPas = (int)(seg.norm() / pas);
        result = result - (1. / 4 * I) * (P[i] * integrale_pour_p(seg, hankel_derivate, nbPas, x) + Q[i] * integ_simple_pour_p(seg, hankel, nbPas, x));
    }
    return result;
}

bool position_dans_obstacles(Point x, const vector<Cercle> &obstacles)
{
    for (unsigned int i = 0; i < obstacles.size(); i++)
    {
        if ((x - obstacles[i].centre).norm() < obstacles[i].rayon)
        {
            return true;
        }
    }
    return false;
}

void exporte_solution(const string &filename, const Maillage &maillage, vector<Cercle> &obstacles, const Vecteur &solution, const Vecteur &P, const double e, const double h, const double pasMaillage, const unsigned int nbPasVisualisation)
{
    ofstream f(filename);

    if (!f.is_open())
    {
        cout << "ERROR: Le fichier " << filename << " n'a pas pu être ouvert" << endl;
        exit(-1);
    }

    for (unsigned int i = 0; i < nbPasVisualisation; i++)
    {
        for (unsigned int j = 0; j < nbPasVisualisation; j++)
        {
            Point IJ((i / (double)nbPasVisualisation) * 5 * e, (j / (double)nbPasVisualisation) * 5 * h);
            if (!position_dans_obstacles(IJ, obstacles))
            {
                complex<double> valeur = p(maillage, 0.01 * pasMaillage, solution, P, IJ);
                f << IJ.x << " " << IJ.y << " " << valeur.real() << " " << valeur.imag() << endl;
            }
        }
    }
    f.close();

    return;
}
