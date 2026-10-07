#include "../second_membre.hpp"

Vecteur calculer_second_membre(const std::vector<QuadratureSegment> &quadrature_maillage)
{
    const unsigned int nbSegmentsMaillage = quadrature_maillage.size();

    Vecteur b_local(nbSegmentsMaillage, 0.0);

    for (unsigned int i = mpi_rank; i < nbSegmentsMaillage; i += mpi_size)
    {
        b_local[i] = integ_simple([](const Point &Q)
                                  { return -u_inc(Q); },
                                  quadrature_maillage[i]);
    }

    Vecteur vect_b(nbSegmentsMaillage, 0.0);

    MPI_Allreduce(b_local.data(), vect_b.data(), static_cast<int>(nbSegmentsMaillage), mpi_complex_type(), MPI_SUM, MPI_COMM_WORLD);

    return vect_b;
}
