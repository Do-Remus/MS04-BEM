#ifndef INTEGRALE_HPP_INCLUDED
#define INTEGRALE_HPP_INCLUDED

#include "../global/global.hpp"
#include "../geometrie/segment.hpp"
#include "legendre.hpp"

/* Fonctions d'intégration */

template <typename F>
inline Complex integ_simple(F &&f, const QuadratureSegment &quadrature)
{
    Complex result = 0.0;

    for (const auto &q : quadrature)
    {
        result += q.weight * f(q.point);
    }

    return result;
}

template <typename F>
inline Complex integ_double(F &&f, const QuadratureSegment &AB, const QuadratureSegment &CD)
{
    Complex result = 0.0;

    for (const auto &qy : CD)
    {
        for (const auto &qx : AB)
        {
            result += qx.weight * qy.weight * f(qx.point, qy.point);
        }
    }

    return result;
}

/* Fonctions pré-intégrées */

inline Complex integ_double_log(const Segment &AB)
{
    const Real L = AB.norm;

    if (L <= PRECISION_ZERO_DOUBLE)
        return Complex(0.0, 0.0);

    return (L * L / (Real(2.0) * pi)) * (Real(1.5) - std::log(L));
}

#endif