#ifndef INTEGRALE_HPP_INCLUDED
#define INTEGRALE_HPP_INCLUDED
#include "../config/config.hpp"
#include "../config/external.hpp"
#include "../config/constantes.hpp"
#include "segment.hpp"
#include "legendre.hpp"

/* Typedef : type pour passer des fonctions en argument */

using fun_d_P = function<complex<double>(const Point &)>;
using fun_double = function<complex<double>(const double &)>;
using fun_d_P2 = function<complex<double>(const Point &, const Point &)>;

/* Fonctions d'integration */

// permet de faire une intègrale simple sur le segment AB
// sur f:R^2->R, avec Nbpas1 intervalles sur AB
complex<double> integ_simple(fun_double f, int n_ordre);
complex<double> integ_simple(fun_double f, double a, double b, int n_ordre);
complex<double> integ_simple(fun_d_P f, Segment &AB, int n_ordre);

#endif