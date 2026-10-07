#ifndef LEGENDRE_HPP_INCLUDED
#define LEGENDRE_HPP_INCLUDED

#include "../../global/commun.hpp"

#include "../geometrie/commun.hpp"

/* Datatype */
struct LegendreData
{
    std::vector<Real> roots;
    std::vector<Real> weights;
};

struct QuadraturePoint
{
    Point point;
    Real weight;
};

using QuadratureSegment = std::vector<QuadraturePoint>;

/* Fonctions sur les polynômes de Legendre avec mémoïsation */

// return the roots of Legendre polynomial of order n
const std::vector<Real> &get_legendre_roots(unsigned n_ordre);

// return roots and weights at root for Legendre polynomial of order n
const LegendreData &get_legendre_data(unsigned n_ordre);

// return the quadrature associated with the legendre polynomial
const std::vector<QuadratureSegment> get_quadrature_maillage(const Maillage &maillage, const LegendreData &data);

#endif