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

        result += w * f(AB.milieu + correction, r);

        if (r > 0.0)
            result += w * f(AB.milieu - correction, -r);
    }

    return 0.5 * AB.norm * result;
}



template <typename F>
std::complex<double> integ_double(F &&f, const Segment &AB, const Segment &CD,
                                  const LegendreData &data)
{
    // Intégrale extérieure sur CD : y = point courant de CD (t = paramètre, ignoré ici)
    return integ_simple(
        [&f, &AB, &data](const Point &y, double /*t*/)
        {
            // Intégrale intérieure sur AB, y fixé : x = point courant de AB
            return integ_simple(
                [&f, &y](const Point &x, double /*s*/)
                {
                    return f(x, y);
                },
                AB, data);
        },
        CD, data);
}

std::complex<double> integ_simple_log(const Point &X, const Segment &AB){
    //a utiliser uniquement pour X sur le segment AB
    const Point DAe = (AB.P1 - X);
    const Point DBe = (AB.P2 - X);
    Point tau = (AB.P2 - AB.P1);
    //verifie X dans AB: DEBUGGAGE
    if(std::abs(tau*DAe)>1e-10 || DBe.norm() > tau.norm() || DAe.norm() > tau.norm()){
        cout<<"X ="<<X<<" n'est pas sur le segment AB, A="<<A<<", B ="<<B<<endl;
        exit(1);
    }
    tau = tau/tau.norm();
    return (1/2*M_PI)*((DBe|tau)*log(DBe.norm()) - (DAe|tau)*log(DAe.norm()) - tau.norm());

}


#endif