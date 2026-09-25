#include "headers/integrale.hpp"

complex<double> integ_simple(fun_double f, int n_ordre)
{
    const LegendreData &data = get_legendre_data(n_ordre);

    complex<double> result = 0.0;

    for (std::size_t i = 0; i < data.roots.size(); ++i)
    {
        const double x = data.roots[i];
        const double w = data.weights[i];

        result += w * f(x);

        if (x != 0.0) // racine symétrique -x (la racine 0 existe une seule fois, n impair)
            result += w * f(-x);
    }

    return result;
}

complex<double> integ_simple(fun_double f, double a, double b, int n_ordre)
{
    const LegendreData &data = get_legendre_data(n_ordre);

    const double milieu = (a + b) / 2.0;
    const double demi_longueur = (b - a) / 2.0;

    complex<double> result = 0.0;

    for (std::size_t i = 0; i < data.roots.size(); ++i)
    {
        const double r = data.roots[i];
        const double w = data.weights[i];
        const double correction = r * demi_longueur;

        result += w * f(milieu + correction);

        if (r > 0.0) // racine symétrique -r
            result += w * f(milieu - correction);
    }

    return demi_longueur * result;
}

complex<double> integ_simple(fun_d_P f, Segment &AB, int n_ordre)
{
    const LegendreData &data = get_legendre_data(n_ordre);

    const Point milieu = (AB.P1 + AB.P2) * 0.5;
    const Point vecteur_demi_longueur = (AB.P1 - AB.P2) * 0.5;

    complex<double> result = 0.0;

    for (std::size_t i = 0; i < data.roots.size(); ++i)
    {
        const double r = data.roots[i];
        const double w = data.weights[i];
        const Point correction = r * vecteur_demi_longueur;

        result += w * f(milieu + correction);

        if (r > 0.0)
            result += w * f(milieu - correction);
    }

    return vecteur_demi_longueur.norm() * result;
}
