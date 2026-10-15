#include "../operator.hpp"

Matrice get_matrice_bem(const Maillage &maillage, const std::vector<QuadratureSegment> &quadrature_maillage)
{
    const Index N = maillage.size();
    Matrice A = Matrice::Zero(N);

    for (Index i = 0; i < N; ++i)
    {
        const Segment &Si = maillage[i];
        const QuadratureSegment &qi = quadrature_maillage[i];

        /*
         * Diagonale régularisée.
         */
        A(i, i) = integ_double_log(Si) + integ_double(green_reguliere, qi, qi);

        /*
         * Hors-diagonale : calcul d'un seul triangle.
         */
        for (Index j = i + 1; j < N; ++j)
        {
            const Complex aij = integ_double(green, qi, quadrature_maillage[j]);

            A(i, j) = aij;
            A(j, i) = aij;
        }
    }

    return A;
}
