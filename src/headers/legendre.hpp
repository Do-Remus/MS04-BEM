#ifndef LEGENDRE_HPP_INCLUDED
#define LEGENDRE_HPP_INCLUDED
#include "../config/config.hpp"
#include "../config/external.hpp"
#include "../config/constantes.hpp"

/* Datatype */
struct LegendreData
{
    std::vector<double> roots;
    std::vector<double> weights;
};

/* Fonctions sur les polynômes de Legendre avec mémoïsation */

// return the roots of Legendre polynomial of order n
const std::vector<double> &get_legendre_roots(unsigned n_ordre);

// return roots and weights at root for Legendre polynomial of order n
const LegendreData &get_legendre_data(unsigned n_ordre);

#endif