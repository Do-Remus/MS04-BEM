#ifndef TEST_FUNCTIONS_HPP_INCLUDED
#define TEST_FUNCTIONS_HPP_INCLUDED
#include "../config/config.hpp"
#include "../config/external.hpp"
#include "../config/constantes.hpp"
#include "point.hpp"

void export_fd_q_p(const string &filename_q, const string &filename_p, double a, int N, unsigned int nbPoints, double h, bool centered = true);

#endif