#ifndef INTEGRALE_HPP_INCLUDED
#define INTEGRALE_HPP_INCLUDED
#include "../config/config.hpp"
#include "../config/external.hpp"
#include "../config/constantes.hpp"
#include "segment.hpp"
#include "legendre.hpp"

/* Fonctions d'integration */

// permet de faire une intègrale simple sur le segment AB
// sur f:R^2->R, avec Nbpas1 intervalles sur AB
template <typename F>
complex<double> integ_simple(F &&f, const LegendreData &data)
{
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

template <typename F>
complex<double> integ_simple(F &&f, double a, double b, const LegendreData &data)
{
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

template <typename F>
std::complex<double> integ_simple(F &&f, const Segment &AB, const LegendreData &data)
{
    std::complex<double> result = 0.0;

    for (std::size_t i = 0; i < data.roots.size(); ++i)
    {
        const double r = data.roots[i];
        const double w = data.weights[i];
        const Point correction = (0.5 * r) * AB.vecteur_norm;

        result += w * f(AB.milieu + correction);

        if (r > 0.0)
            result += w * f(AB.milieu - correction);
    }

    return 0.5 * AB.norm * result;
}

#endif