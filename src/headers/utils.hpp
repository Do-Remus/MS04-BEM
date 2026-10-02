#ifndef UTILS_HPP_INCLUDED
#define UTILS_HPP_INCLUDED
#include "../config/config.hpp"
#include "../config/external.hpp"
#include "../config/constantes.hpp"
#include "maillage.hpp"
#include "matrice.hpp"
#include "integrale.hpp"

// Initialisateur du cache
void initialiser_green_cache();

// Hankel premiere espece ordre n
inline complex<Real> hankel_n(const Real x, const int n)
{
    const Real J = boost::math::cyl_bessel_j(n, x);
    const Real Y = boost::math::cyl_neumann(n, x);

    return {J, Y};
}

// Green
inline complex<Real> green(const Point &p1, const Point &p2)
{
    Real r = (p1 - p2).norm();

    return (I / Real(4.0)) * hankel_n(k * r, 0);
}

// Green cached with a vector (faster but fixed size at compile time)
inline complex<Real> green_cached_vec(const Point &p1, const Point &p2)
{
    const Real dx = p1.x - p2.x;
    const Real dy = p1.y - p2.y;

    const Real d2 = dx * dx + dy * dy;

    const std::size_t index = static_cast<std::size_t>(std::sqrt(d2) * green_cache_step_inv + 0.5);

    return green_cache[index];
}

// Partie reguliere de Green
inline complex<Real> green_reguliere(const Point &P1, const Point &P2)
{
    Real x = (P2 - P1).norm();

    if (x < PRECISION_ZERO_DOUBLE)
    {
        return (I / Real(4.0)) - (Real(1.0) / (Real(2.0) * pi)) * (gamma_euler + log(k / 2));
    }
    return green(P1, P2) + (Real(1.0) / (Real(2.0) * pi)) * log(x);
}

// Partie reguliere de Green cached
inline complex<Real> green_reguliere_cached_vec(const Point &P1, const Point &P2)
{
    Real x = (P2 - P1).norm();

    if (x < green_cache_step)
    {
        return (I / Real(4.0)) - (Real(1.0) / (Real(2.0) * pi)) * (gamma_euler + log(k / 2));
    }
    return green_cached_vec(P1, P2) + (Real(1.0) / (Real(2.0) * pi)) * log(x);
}

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

// Calcul de Ax parallelisé
Vecteur produit_A_cached(const Maillage &maillage, const Vecteur &x, const std::vector<QuadratureSegment> &quadrature_maillage);

// Calcul de Ax par blocks et parallelisé
Vecteur produit_A_cached_blocked(const Maillage &maillage, const Vecteur &x, const std::vector<QuadratureSegment> &quadrature_maillage);

// Algorithme COCG pour resoudre un system pour A sym complexe de manière itérative sans caculer stocker A
Vecteur gradConjMatrixFree(const Maillage &maillage, const Vecteur &b, const std::vector<QuadratureSegment> &quadrature_maillage, Real tol, unsigned int maxIter);

#endif