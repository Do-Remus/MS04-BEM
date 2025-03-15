#ifndef CONFIG_HPP_INCLUDED
#define CONFIG_HPP_INCLUDED

#include "external.hpp"
#include "constantes.hpp"

extern double theta;
extern bool effectuerTests;
extern bool effectuerLaSimulation;
extern double k; // frequence des ondes : il y a un nombre au plus denombrable de frequences pour lesquelles le pb n'est pas bien pose

complex<double> p_theta_config(double x, double y);

#endif