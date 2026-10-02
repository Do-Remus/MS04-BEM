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

inline Complex integ_double_green_vec(const QuadratureSegment &AB, const QuadratureSegment &CD)
{
    Complex result = 0.0;

    for (const auto &qy : CD)
    {
        for (const auto &qx : AB)
        {

            const Real dx = qx.point.x - qy.point.x;
            const Real dy = qx.point.y - qy.point.y;

            const Real d2 = dx * dx + dy * dy;

            const std::size_t index = static_cast<std::size_t>(std::sqrt(d2) * green_cache_step_inv + 0.5);

            result += qy.weight * qx.weight * green_cache[index];
        }
    }

    return result;
}

Complex integ_simple_log(const Point &X, const Segment &AB);

#endif