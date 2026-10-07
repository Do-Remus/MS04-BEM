#include "../solutions_analytiques.hpp"

Complex u_inc(const Point &P1)
{
    Real theta = P1.theta();
    return exp(-I * k * P1.norm() * cos(theta));
}

Complex u_N_plus_analytique(const Point &P1, Real radius, int N)
{
    Real theta = P1.theta();
    Real x = k * P1.norm();
    Real ka = k * radius;
    Complex iterative_i = 1;
    Complex partial_sum = -(boost::math::cyl_bessel_j(0, ka) / hankel_n(ka, 0)) * hankel_n(x, 0);

    for (int n = 1; n <= N; n++)
    {
        iterative_i *= -I;
        partial_sum -= Real(2.0) * iterative_i * (boost::math::cyl_bessel_j(n, ka) / hankel_n(ka, n)) * hankel_n(x, n) * cos(n * theta); // positive part of the sum
    }

    return partial_sum;
}

Complex q_analytique(const Point &P1, int N)
{
    Real theta = P1.theta();
    Real ka = k * P1.norm();
    Complex iterative_i = 1;
    Complex partial_sum = k * boost::math::cyl_bessel_j(0, ka) * hankel_n(ka, 1) / hankel_n(ka, 0);

    for (int n = 1; n <= N; n++)
    {
        iterative_i *= -I;
        partial_sum -= k * iterative_i * (boost::math::cyl_bessel_j(n, ka) * ((hankel_n(ka, n - 1) - hankel_n(ka, n + 1))) / hankel_n(ka, n)) * cos(n * theta); // positive and negative part of the sum
    }

    return partial_sum;
}

Complex p_analytique(const Point &P1, int N)
{
    Real theta = P1.theta();
    Real ka = k * P1.norm();
    Complex iterative_i = 1;
    Complex partial_sum = k * boost::math::cyl_bessel_j(1, ka) - k * boost::math::cyl_bessel_j(0, ka) * hankel_n(ka, 1) / hankel_n(ka, 0);

    for (int n = 1; n <= N; n++)
    {
        iterative_i *= -I;
        partial_sum += k * iterative_i * cos(n * theta) * (boost::math::cyl_bessel_j(n, ka) * (hankel_n(ka, n - 1) - hankel_n(ka, n + 1)) / hankel_n(ka, n) - (boost::math::cyl_bessel_j(n - 1, ka) - boost::math::cyl_bessel_j(n + 1, ka)));
    }

    return partial_sum;
}