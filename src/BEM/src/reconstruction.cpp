#include "../reconstruction.hpp"

std::vector<Complex> reconstruire_solution(const std::vector<Point> &pointsSolution, const Maillage &maillage, const std::vector<QuadratureSegment> &quadrature_maillage, const Vecteur &vect_p_cached)
{
    const unsigned int nbPoints = pointsSolution.size();

    // Chaque processus ne calcule qu'une partie des points.
    std::vector<Complex> u_local(nbPoints, Complex(0.0));

    for (unsigned int j = mpi_rank; j < nbPoints; j += mpi_size)
    {
        const Point &Pj = pointsSolution[j];

        for (unsigned int i = 0; i < maillage.size(); i++)
        {
            u_local[j] += vect_p_cached[i] * integ_simple([&Pj](const Point &Q)
                                                          { return green_cached(Pj, Q); },
                                                          quadrature_maillage[i]);
        }
    }

    // Tous les processus récupèrent le vecteur complet.
    std::vector<Complex> u_cached(nbPoints, Complex(0.0));

    MPI_Allreduce(u_local.data(), u_cached.data(), static_cast<int>(nbPoints), mpi_complex_type(), MPI_SUM, MPI_COMM_WORLD);

    return u_cached;
}
