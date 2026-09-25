#ifndef INTEGRALE_HPP_INCLUDED
#define INTEGRALE_HPP_INCLUDED
#include "../config/config.hpp"
#include "../config/external.hpp"
#include "../config/constantes.hpp"
#include "segment.hpp"

/* Typedef : type pour passer des fonctions en argument */

typedef complex<double> (*fun_d_P)(const Point &);

typedef complex<double> (*fun_double)(const double &);

typedef complex<double> (*fun_d_P2)(const Point &, const Point &);

/* Fonctions d'integration */

// permet de faire une intègrale simple sur le segment AB
// sur f:R^2->R, avec Nbpas1 intervalles sur AB
complex<double> integ_simple(Segment &AB, fun_d_P f, int Nbpas1);

complex<double> integ_simple_segsimple(fun_double f);

complex<double> integ_simple_segsimple_n(fun_double f, int n_ordre);
complex<double> integ_simple_segment_ab_n(fun_double f, int n_ordre, double a, double b);

// permet de faire une intègrale double sur les segments AB et CD
// sur f:R^2xR^2->R, avec pas1 intervalles sur AB et pas2 intervalles sur CD
complex<double> integ_double(Segment &AB, Segment &CD, fun_d_P2 f, int Nbpas1, int Nbpas2);

// fonction donnant le log de la norme de la difference entre 2 points
complex<double> log_norm(const Point &A, const Point &B);

// calcul de l'integral du log de la norme selon 2 segments
complex<double> integ_log(Segment &AB, Segment &CD, int Nbpas1, int Nbpas2);

// fonction integrale utilisé pour le calcul de p : incorporation normale au segment
complex<double> integrale_pour_p(Segment &AB, fun_d_P f, int Nbpas1, const Point &x);

// fonction integrale simple pour calcul de p
complex<double> integ_simple_pour_p(Segment &AB, fun_d_P f, int Nbpas1, const Point &x);

#endif