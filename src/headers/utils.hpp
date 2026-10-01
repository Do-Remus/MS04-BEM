#ifndef UTILS_HPP_INCLUDED
#define UTILS_HPP_INCLUDED
#include "../config/config.hpp"
#include "../config/external.hpp"
#include "../config/constantes.hpp"
#include "maillage.hpp"
#include "matrice.hpp"
#include "integrale.hpp"

// Cache pour green
extern Complex *green_cache;
extern MPI_Win green_cache_win;
extern Real green_cache_step;
extern Real green_cache_step_inv;
extern std::size_t max_index;
extern Real delta;

// Initialisateur du cache
void initialiser_green_cache();

// Hankel premiere espece ordre n
complex<Real> hankel_n(const Real x, const int n);

// Green
complex<Real> green(const Point &p1, const Point &p2);

// Green cached with a vector (faster but fixed size at compile time)
complex<Real> green_cached_vec(const Point &p1, const Point &p2);

// Partie reguliere de Green
complex<Real> green_reguliere(const Point &P1, const Point &P2);

// Partie reguliere de Green cached
complex<Real> green_reguliere_cached_vec(const Point &P1, const Point &P2);

// solution approche (jusqu'au terme N de la somme) exterieure pour cas 1 disque de rayon radius evalué au point P = (r, theta)
complex<Real> u_N_plus_analytique(const Point &P1, Real radius, int N);

// u incident
complex<Real> u_inc(const Point &P1);

// dn(u+) sur le bord du disque Gamma (check P sur Gamma a faire en amont)
complex<Real> q_analytique(const Point &P1, int N);

// - dn(u+) - dn(u_inc) sur le bord du disque Gamma (check P sur Gamma a faire en amont)
complex<Real> p_analytique(const Point &P1, int N);

// Exporte les obstacles sous forme de cercles
void export_obsctacles(const string &filename, vector<Cercle> cercles);

// Calcul plus propre de l'erreur
Real erreur_relative(const Complex &u, const Complex &reference, Real eps = 1e-14);

#endif