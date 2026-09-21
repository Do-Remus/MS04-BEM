#ifndef UTILS_HPP_INCLUDED
#define UTILS_HPP_INCLUDED
#include "../config/config.hpp"
#include "../config/external.hpp"
#include "../config/constantes.hpp"
#include "maillage.hpp"
#include "matrice.hpp"
#include "integrale.hpp"

//Hankel premiere espece ordre n
complex<double> hankel_n(const double x, const int n);


//solution approche (jusqu'au terme N de la somme) exterieure pour cas 1 disque de rayon radius evalué au point P = (r, theta)
complex<double> u_N_plus(const Point &P1, double radius, int N,  const string& filename);

// dn(u+) sur le bord du disque Gamma (check P sur Gamma a faire en amont)
complex<double> q(const Point &P1, int N,  const string& filename);

// - dn(u+) - dn(u_inc) sur le bord du disque Gamma (check P sur Gamma a faire en amont)
complex<double> p(const Point &P1, int N,  const string& filename);


//solution u+ sur le domaine exterieur borné sur un carré de longueur totale = 2*longueur, stockage dans filename, approx à l'ordre N
void exporte_solution_analytique(const string &filename, const double radius_obstacle, const unsigned int nbPasVisualisation, const double longueur, const int N);




/*
// teste si un obstacle est dans une position valide (dans le carré et sans intersections avec les autres obstacles déjà générés)
bool obstacle_valide(double rayonAleatoire, double xAleatoire, double yAleatoire, vector<Cercle> &obstacles, double hauteur, double epaisseur);

// Genere le maillage associé à une couche diffusante de paramètres données (nombre d'obstacles, épaisseur, hauteur,...)
Maillage genere_maillage_couche_diffusante(unsigned int nbObstacles, double hauteur, double epaisseur, double pas_maillage, vector<Cercle> &obstacles);

// Fonction de green formatée pour prendre des points en entrée
complex<double> green_reguliere(const Point &P1, const Point &P2);

// Fonction p_theta de l'énoncé formatée pour prendre des points en entrée
complex<double> p_theta(const Point &P);

// Genere les coefficients du second membre P de l'équation
int genere_coefficient_vecteur_P(Vecteur &P, const Maillage &maillage, double pas);

// Genere les coefficients de la matrice A de l'équation
int genere_coefficient_matrice_A(MatriceSym &A, const Maillage &maillage, double pas);

// Hankel nu=0 1er ordre
complex<double> hankel(const Point &P1);


// derive de hankel nu=0, 1er ordre
complex<double> hankel_derivate(const Point &P1);

// calcule la valeur de p(x) pour x dans le rectangle h x e
complex<double> p(const Maillage &maillage, double pas, const Vecteur &Q, const Vecteur P, const Point &x);

// Test si le point est dans l'un des obstacles
bool position_dans_obstacles(Point x, const vector<Cercle> &obstacles);

// Reconstruit et exporte les champs p(x)
void exporte_solution(const string &filename, const Maillage &maillage, vector<Cercle> &obstacles, const Vecteur &solution, const Vecteur &P, const double e, const double h, const double pasMaillage, const unsigned int nbPasVisualisation);
*/

// Exporte les obstacles sous forme de cercles
void export_obsctacles(const string &filename, vector<Cercle> cercles);


#endif