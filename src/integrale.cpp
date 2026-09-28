#include "headers/integrale.hpp"

std::complex<double> integ_simple_log(const Point &X, const Segment &AB)
{
    // a utiliser uniquement pour X sur le segment AB
    const Point DAe = (AB.P1 - X);
    const Point DBe = (AB.P2 - X);
    Point tau = (AB.P2 - AB.P1);
    // verifie X dans AB: DEBUGGAGE
    if (std::abs(tau * DAe) > 1e-10 || DBe.norm() > tau.norm() || DAe.norm() > tau.norm())
    {
        cout << "X =" << X << " n'est pas sur le segment AB, A=" << AB.P1 << ", B =" << AB.P2 << endl;
        exit(1);
    }
    tau = tau / tau.norm();
    return -(1.0 / (2.0 * pi)) * ((DBe | tau) * log(DBe.norm()) - (DAe | tau) * log(DAe.norm()) - tau.norm());
}
