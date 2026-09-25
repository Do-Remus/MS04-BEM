#ifndef UTILS_HPP_INCLUDED
#define UTILS_HPP_INCLUDED
#include "../config/config.hpp"
#include "../config/external.hpp"
#include "../config/constantes.hpp"
#include "maillage.hpp"
#include "matrice.hpp"
#include "integrale.hpp"

// Hankel premiere espece ordre n
complex<double> hankel_n(const double x, const int n);

// solution approche (jusqu'au terme N de la somme) exterieure pour cas 1 disque de rayon radius evalué au point P = (r, theta)
complex<double> u_N_plus_analytique(const Point &P1, double radius, int N, const string &filename);

// dn(u+) sur le bord du disque Gamma (check P sur Gamma a faire en amont)
complex<double> q_analytique(const Point &P1, int N, const string &filename);

// - dn(u+) - dn(u_inc) sur le bord du disque Gamma (check P sur Gamma a faire en amont)
complex<double> p_analytique(const Point &P1, int N, const string &filename);

// solution u+ sur le domaine exterieur borné sur un carré de longueur totale = 2*longueur, stockage dans filename, approx à l'ordre N
void exporte_solution_analytique(const string &filename, const double radius_obstacle, const unsigned int nbPasVisualisation, const double longueur, const int N);

// Exporte les obstacles sous forme de cercles
void export_obsctacles(const string &filename, vector<Cercle> cercles);

#endif