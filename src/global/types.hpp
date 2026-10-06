#ifndef TYPES_HPP_INCLUDED
#define TYPES_HPP_INCLUDED

#include "external.hpp"

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