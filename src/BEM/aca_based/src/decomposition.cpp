#include "../decomposition.hpp"

ACA_decomposition aca_bem(const std::vector<Index> &selection_i, const std::vector<Index> &selection_j, const Maillage &maillage, const std::vector<QuadratureSegment> &quadrature_maillage, const Index rang_max, const double tolerance)
{
    auto coefficient = [&](const Index i, const Index j) -> Complex
    {
        return matrice_bem(i, j, maillage, quadrature_maillage);
    };

    return aca_generique(selection_i, selection_j, coefficient, rang_max, tolerance);
}

ACA_decomposition aca_reconstruction(const std::vector<Index> &selection_i, const std::vector<Index> &selection_j, const std::vector<Point> &pointsSolution, const std::vector<QuadratureSegment> &quadrature_maillage, const Index rang_max, const double tolerance)
{
    auto coefficient = [&](const Index i, const Index j) -> Complex
    {
        return coefficient_reconstruction(i, j, pointsSolution, quadrature_maillage);
    };

    return aca_generique(selection_i, selection_j, coefficient, rang_max, tolerance);
}
