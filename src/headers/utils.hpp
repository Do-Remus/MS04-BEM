#ifndef UTILS_HPP_INCLUDED
#define UTILS_HPP_INCLUDED
#include "maillage.hpp"
#include "matrice.hpp"
#include "integrale.hpp"

// teste si un obstacle est dans une position valide (dans le carré et sans intersections avec les autres obstacles déjà générés)
bool obstacle_valide(double rayonAleatoire, double xAleatoire, double yAleatoire, vector<Cercle> obstacles, double hauteur, double epaisseur);

// Genere le maillage associé à une couche diffusante de paramètres données (nombre d'obstacles, épaisseur, hauteur,...)
Maillage genere_maillage_couche_diffusante(unsigned int nbObstacles, double hauteur, double epaisseur, double pas_maillage);

// Fonction de green formatée pour prendre des points en entrée
complex<double> green_reguliere(const Point &P1, const Point &P2);

// Fonction p_theta de l'énoncé formatée pour prendre des points en entrée
complex<double> p_theta(const Point &P);

// Genere les coefficients du second membre P de l'équation
int genere_coefficient_vecteur_P(Vecteur &P, const Maillage &maillage, double pas);

// Genere les coefficients de la matrice A de l'équation
int genere_coefficient_matrice_A(MatriceSym &A, const Maillage &maillage, double pas);

//Hankel nu=0 1er ordre
complex<double> hankel(const Point& P1);

//derive de hankel nu=0, 1er ordre
complex<double> hankel_derivate(const Point& P1);

//calcule la valeur de p(x) pour x dans le rectangle h x e
complex<double> p(const Maillage &maillage, double pas, const Vecteur& Q,const Vecteur P,  const Point& x );

#endif