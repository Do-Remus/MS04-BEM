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

// Stockage optimisé pour SIMD
struct QuadratureSIMD
{
    std::size_t N;
    std::size_t nq;

    /*
     * Organisation cible :
     *
     * target_x[i * nq + a]
     * target_y[i * nq + a]
     * target_w[i * nq + a]
     */
    std::vector<Real> target_x;
    std::vector<Real> target_y;
    std::vector<Real> target_w;

    /*
     * Organisation source transposée :
     *
     * source_x[a * N + j]
     * source_y[a * N + j]
     * source_w[a * N + j]
     *
     * Ainsi j,j+1,j+2,j+3 sont contigus.
     */
    std::vector<Real> source_x;
    std::vector<Real> source_y;
    std::vector<Real> source_w;
};

/* Fonctions sur les polynômes de Legendre avec mémoïsation */

// return the roots of Legendre polynomial of order n
const std::vector<Real> &get_legendre_roots(unsigned n_ordre);

// return roots and weights at root for Legendre polynomial of order n
const LegendreData &get_legendre_data(unsigned n_ordre);

// return the quadrature associated with the legendre polynomial
const std::vector<QuadratureSegment> get_quadrature_maillage(const Maillage &maillage, const LegendreData &data);

// return the optimised for SIMD version of a the quadrature for all segments
QuadratureSIMD construire_quadrature_SIMD(const Maillage &maillage, const std::vector<QuadratureSegment> &quadrature_maillage);

#endif