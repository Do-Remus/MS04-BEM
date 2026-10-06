#ifndef FONCTIONS_ANALYTIQUES_HPP_INCLUDED
#define FONCTIONS_ANALYTIQUES_HPP_INCLUDED

#include "../global/global.hpp"
#include "../geometrie/point.hpp"
#include "green.hpp"

Complex u_inc(const Point &P1);

Complex u_N_plus_analytique(const Point &P1, Real radius, int N);

Complex q_analytique(const Point &P1, int N);

Complex p_analytique(const Point &P1, int N);

#endif