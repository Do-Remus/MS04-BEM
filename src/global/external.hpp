#ifndef EXTERNAL_HPP_INCLUDED
#define EXTERNAL_HPP_INCLUDED

// General
#include <time.h>
#include <cstdlib>
#include <stdexcept>
#include <iomanip>
#include <unordered_map>
#include <optional>

// Maths
#include <cmath>
#include <complex>
#include <vector>
#include <functional>
#include <random>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/legendre.hpp>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <unsupported/Eigen/IterativeSolvers>

// Parallel work
#include <mpi.h>
extern int mpi_rank;
extern int mpi_size;

// Reading files
#include <algorithm>
#include <cctype>
#include <iostream>
#include <fstream>

using namespace std;

#endif