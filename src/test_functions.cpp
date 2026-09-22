#include "headers/test_functions.hpp"
#include "headers/utils.hpp"

// ---------------------------------------------------------------------------
// Needs: #include <functional> #include <iomanip> #include <cmath>   (C++11 is enough)
// To paste in utils.cpp (and declare in utils.hpp):
//   void export_fd_q_p(const string &filename_q, const string &filename_p,
//                      double a, int N, unsigned int nbPoints, double h, bool centered = true);
// ---------------------------------------------------------------------------

// Finite-difference approximation of
//     q = d(u+)/dn            (n = e_r on the circle, so d/dn = d/dr)
//     p = -d(u+)/dn - d(uinc)/dn = -d(u+ + uinc)/dr
// on the boundary r = a, compared with the analytic q() and p().
// One line per angle:  theta  Re(fd) Im(fd)  Re(exact) Im(exact)  |fd - exact|
// (lines starting with '#' are a header: k, a, N, h, scheme).
//
//   centered = true : (f(a+h) - f(a-h)) / (2h)               O(h^2), uses r < a (the series
//                                                             is analytic there)
//   centered = false: (-3 f(a) + 4 f(a+h) - f(a+2h)) / (2h)  O(h^2), only uses r >= a
//
// Angles are the midpoints theta_j = pi (j + 1/2) / nbPoints in (0, pi): this avoids
// theta = 0 and pi, where acos((A|P)/r) can return NaN through rounding (ratio = 1 + eps).
void export_fd_q_p(const string &filename_q, const string &filename_p,
                   double a, int N, unsigned int nbPoints, double h, bool centered)
{
    ofstream fq(filename_q);
    ofstream fp(filename_p);
    if (!fq.is_open() || !fp.is_open())
    {
        cout << "ERROR: impossible d'ouvrir " << filename_q << " ou " << filename_p << endl;
        exit(-1);
    }
    const string bin = "outputs/bin.txt"; // scratch file required by u_N_plus, q and p

    for (ofstream *f : {&fq, &fp})
    {
        *f << setprecision(17);
        *f << "# k=" << k << " a=" << a << " N=" << N << " h=" << h
           << " scheme=" << (centered ? "centered" : "one-sided") << "\n";
        *f << "# theta Re(fd) Im(fd) Re(exact) Im(exact) abs_error\n";
    }

    double maxErrQ = 0, maxErrP = 0;
    for (unsigned int j = 0; j < nbPoints; j++)
    {
        const double th = 2 * pi * (j + 0.5) / nbPoints;
        const Point P(a * cos(th), a * sin(th));

        // u+ (scattered) and u+ + uinc (total) at radius r along the ray of angle theta
        auto u_scat = [&](double r) -> complex<double>
        { return u_N_plus(Point(r * cos(th), r * sin(th)), a, N, bin); };
        auto u_tot = [&](double r) -> complex<double>
        { return u_scat(r) + exp(-I * k * r * cos(th)); };
        auto deriv = [&](const function<complex<double>(double)> &f) -> complex<double>
        {
            if (centered)
                return (f(a + h) - f(a - h)) / (2. * h);
            return (f(a + h) - f(a)) / (h);
        };

        const complex<double> q_fd = deriv(u_scat);
        const complex<double> p_fd = -deriv(u_tot);
        const complex<double> q_ex = q(P, N, bin);
        const complex<double> p_ex = p(P, N, bin);

        fq << th << " " << q_fd.real() << " " << q_fd.imag() << " " << q_ex.real() << " "
           << q_ex.imag() << " " << abs(q_fd - q_ex) << "\n";
        fp << th << " " << p_fd.real() << " " << p_fd.imag() << " " << p_ex.real() << " "
           << p_ex.imag() << " " << abs(p_fd - p_ex) << "\n";
        maxErrQ = max(maxErrQ, abs(q_fd - q_ex));
        maxErrP = max(maxErrP, abs(p_fd - p_ex));
    }
    cout << "FD check (h=" << h << "): max |q_fd - q| = " << maxErrQ
         << ", max |p_fd - p| = " << maxErrP << endl;
}