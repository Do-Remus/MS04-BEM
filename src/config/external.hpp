#ifndef EXTERNAL_HPP_INCLUDED
#define EXTERNAL_HPP_INCLUDED

// General
#include <time.h>
#include <cstdlib>
#include <iomanip>
#include <unordered_map>

// Maths
#include <cmath>
#include <complex>
#include <vector>
#include <functional>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/legendre.hpp>

// Parallel work
#include <mpi.h>

// Reading files
#include <algorithm>
#include <cctype>
#include <iostream>
#include <fstream>

using namespace std;

using Real = double;
using Complex = std::complex<Real>;

inline MPI_Datatype mpi_complex_type()
{
    if constexpr (std::is_same_v<Real, float>)
        return MPI_C_FLOAT_COMPLEX;
    else
        return MPI_C_DOUBLE_COMPLEX;
}

#endif