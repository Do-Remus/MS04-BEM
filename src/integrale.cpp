#include "headers/integrale.hpp"

std::complex<double> integ_simple_log(const Point &X, const Segment &AB)
{
    const Point DAe = AB.P1 - X;
    const Point DBe = AB.P2 - X;

    const double normDAe = DAe.norm();
    const double normDBe = DBe.norm();

    const Point tau = AB.vecteur_norm / AB.norm;

    if (std::abs(AB.vecteur_norm * DAe) > PRECISION_ZERO_DOUBLE || normDBe > AB.norm || normDAe > AB.norm)
    {
        cout << "X =" << X
             << " n'est pas sur le segment AB, A=" << AB.P1
             << ", B=" << AB.P2 << endl;
        exit(1);
    }

    auto x_log_x = [](double x)
    {
        if (x == 0.0)
            return 0.0;

        return x * std::log(x);
    };

    const double a = DAe | tau;
    const double b = DBe | tau;

    return -(1.0 / (2.0 * pi)) * (x_log_x(std::abs(b)) + x_log_x(std::abs(a)) - AB.norm);
}
