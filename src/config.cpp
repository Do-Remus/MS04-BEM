#include "config/config.hpp"

/* Définition des varibles globales */

double theta = 0;
bool effectuerTests = false;
bool effectuerLaSimulation = true;
double k = (4*pi*pi)+0.5; // frequence des ondes : il y a un nombre au plus denombrable de frequences pour lesquelles le pb n'est pas bien pose
//k doit etre strictement positif

/* Déclaration des fonctions globales */

complex<double> p_theta_config(double x, double y)
{
    return (exp(I * k*(x * cos(theta) + y * sin(theta))));
}
