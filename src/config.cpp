#include "config/config.hpp"

/* Définition des varibles globales */

double theta = 0;
bool effectuerTests = true;
bool effectuerLaSimulation = true;

/* Déclaration des fonctions globales */

complex<double> p_theta_config(double x, double y)
{
    return (exp(I * (x * cos(theta) + y * sin(theta))));
}
