#include "../reconstruction.hpp"

Vecteur reconstruire_solution(const std::vector<Point> &pointsSolution, const Maillage &maillage, const std::vector<QuadratureSegment> &quadrature_maillage, const Vecteur &vect_p_cached)
{
    const Index nbPoints = pointsSolution.size();

    // Chaque processus ne calcule qu'une partie des points.
    Vecteur u_local = Vecteur::Zero(nbPoints);

    for (Index j = mpi_rank; j < nbPoints; j += mpi_size)
    {
        const Point &Pj = pointsSolution[j];

        for (Index i = 0; i < maillage.size(); i++)
        {
            u_local[j] += vect_p_cached[i] * integ_simple([&Pj](const Point &Q)
                                                          { return green_cached(Pj, Q); },
                                                          quadrature_maillage[i]);
        }
    }

    // Tous les processus récupèrent le vecteur complet.
    Vecteur u_cached = Vecteur::Zero(nbPoints);

    MPI_Allreduce(u_local.data(), u_cached.data(), static_cast<int>(nbPoints), mpi_complex_type(), MPI_SUM, MPI_COMM_WORLD);

    return u_cached;
}
