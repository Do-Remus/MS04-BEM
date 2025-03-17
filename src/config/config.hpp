#ifndef CONFIG_HPP_INCLUDED
#define CONFIG_HPP_INCLUDED

#include "external.hpp"
#include "constantes.hpp"

extern double theta;
extern bool effectuerTests;
extern bool effectuerLaSimulation;
extern double k; // frequence des ondes : il y a un nombre au plus denombrable de frequences pour lesquelles le pb n'est pas bien pose
extern double h;
extern double e;
extern double pasMaillage;
extern unsigned int nbObstacles;
extern double pasIntegrale;
extern unsigned int nbPasExport;

complex<double> p_theta_config(double x, double y);
void get_config(const string &filename);

#endif