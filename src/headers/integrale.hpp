#ifndef INTEGRALE_HPP_INCLUDED
#define INTEGRALE_HPP_INCLUDED
#include "../config/config.hpp"
#include "../config/external.hpp"
#include "../config/constantes.hpp"
#include "segment.hpp"
#include "legendre.hpp"

/* Typedef : type pour passer des fonctions en argument */

typedef complex<double> (*fun_d_P)(const Point &);

typedef complex<double> (*fun_double)(const double &);

typedef complex<double> (*fun_d_P2)(const Point &, const Point &);

/* Fonctions d'integration */

// permet de faire une intègrale simple sur le segment AB
// sur f:R^2->R, avec Nbpas1 intervalles sur AB
complex<double> integ_simple(fun_double f, int n_ordre);
complex<double> integ_simple(fun_double f, double a, double b, int n_ordre);
complex<double> integ_simple(fun_d_P f, Segment &AB, int n_ordre);

#endif