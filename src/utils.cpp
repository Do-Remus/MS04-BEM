#include "headers/utils.hpp"

#include <stdexcept>
#include <limits>
/*
complex<double> hankel_n(const double x, const int n)
{
    const double J = boost::math::cyl_bessel_j(n, x);
    const double Y = boost::math::cyl_neumann(n, x);

    return {J, Y};
}
*/
complex<double> hankel_n(const double x, const int n)
{
    double Jn, Yn;
    try
    {
        Jn = jn(n, x);
    }
    catch (const std::exception &)
    {
        Jn = 0.0;   // J_n(x) underflows towards 0 for n large at fixed x; it never overflows
    }
    try
    {
        Yn = yn(n, x);
    }
    catch (const std::exception &)
    {
        // Y_n(x) -> -infinity as n grows past x, monotonically (no sign oscillation),
        // so -inf is always the right value to substitute here.
        Yn = -std::numeric_limits<double>::infinity();
    }
    return Jn + I * Yn;
}



complex<double> green(const Point &p1, const Point &p2)
{
    double r = (p1 - p2).norm();

    return (I / 4.0) * hankel_n(k * r, 0);
}

complex<double> green_cached_map(const Point &p1, const Point &p2)
{
    static std::unordered_map<std::size_t, complex<double>> green_cache;

    const double distance = (p1 - p2).norm();

    const std::size_t index = static_cast<std::size_t>(distance / green_cache_step + 0.5);

    auto it = green_cache.find(index);

    if (it != green_cache.end())
    {
        return it->second;
    }

    const double quantized_distance = index * green_cache_step;

    const complex<double> value =
        (I / 4.0) * hankel_n(k * quantized_distance, 0);

    green_cache[index] = value;

    return value;
}

complex<double> green_cached_vec(const Point &p1, const Point &p2)
{
    static std::vector<std::complex<double>> green_cache(max_index + 1);
    static std::vector<bool> computed(max_index + 1, false);

    const double distance = (p1 - p2).norm();

    const std::size_t index = static_cast<std::size_t>(distance / green_cache_step + 0.5);

    assert(index < max_index + 1);

    if (computed[index])
        return green_cache[index];

    const double quantized_distance = index * green_cache_step;

    const complex<double> value =
        (I / 4.0) * hankel_n(k * quantized_distance, 0);

    green_cache[index] = value;
    computed[index] = true;

    return value;
}

complex<double> green_reguliere(const Point &P1, const Point &P2)
{
    double x = (P2 - P1).norm();

    if (x < PRECISION_ZERO_DOUBLE)
    {
        return (I / 4.0) - (1.0 / (2.0 * pi)) * (gamma_euler + log(k / 2));
    }
    return green(P1, P2) + (1 / (2 * pi)) * log(x);
}

complex<double> green_reguliere_cached_vec(const Point &P1, const Point &P2)
{
    double x = (P2 - P1).norm();

    if (x < PRECISION_ZERO_DOUBLE)
    {
        return (I / 4.0) - (1.0 / (2.0 * pi)) * (gamma_euler + log(k / 2));
    }
    return green_cached_vec(P1, P2) + (1 / (2 * pi)) * log(x);
}

#ifdef TP0
complex<double> u_N_plus_analytique(const Point &P1, double radius, int N, const string &filename)
#else
complex<double> u_N_plus_analytique(const Point &P1, double radius, int N)
#endif
{
#ifdef TP0
    ofstream f(filename);

    if (!f.is_open())
    {
        cout << "ERROR: Le fichier " << filename << " n'a pas pu être ouvert" << endl;
        exit(-1);
    }
#endif

    double theta = P1.theta();
    double x = k * P1.norm();
    double ka = k * radius;
    complex<double> iterative_i = 1;
    complex<double> partial_sum = -(boost::math::cyl_bessel_j(0, ka) / hankel_n(ka, 0)) * hankel_n(x, 0);

#ifdef TP0
    f << 0 << " " << partial_sum.real() << " " << partial_sum.imag() << endl;
#endif

    for (int n = 1; n <= N; n++)
    {
        iterative_i *= -I;
        partial_sum -= 2. * iterative_i * (boost::math::cyl_bessel_j(n, ka) / hankel_n(ka, n)) * hankel_n(x, n) * cos(n * theta); // positive part of the sum
#ifdef TP0
        f << n << " " << partial_sum.real() << " " << partial_sum.imag() << endl;
#endif
    }

#ifdef TP0
    f.close();
#endif

    return partial_sum;
}

complex<double> u_inc(const Point &P1)
{
    double theta = P1.theta();
    return exp(-I * k * P1.norm() * cos(theta));
}

#ifdef TP0
complex<double> q_analytique(const Point &P1, int N, const string &filename)
#else
complex<double> q_analytique(const Point &P1, int N)
#endif
{
#ifdef TP0
    ofstream f(filename);

    if (!f.is_open())
    {
        cout << "ERROR: Le fichier " << filename << " n'a pas pu être ouvert" << endl;
        exit(-1);
    }
#endif

    double theta = P1.theta();
    double ka = k * P1.norm();
    complex<double> iterative_i = 1;
    complex<double> partial_sum = k * boost::math::cyl_bessel_j(0, ka) * hankel_n(ka, 1) / hankel_n(ka, 0);

#ifdef TP0
    f << 0 << " " << partial_sum.real() << " " << partial_sum.imag() << endl;
#endif

    for (int n = 1; n <= N; n++)
    {
        iterative_i *= -I;
        partial_sum -= k * iterative_i * (boost::math::cyl_bessel_j(n, ka) * ((hankel_n(ka, n - 1) - hankel_n(ka, n + 1))) / hankel_n(ka, n)) * cos(n * theta); // positive and negative part of the sum
#ifdef TP0
        f << n << " " << partial_sum.real() << " " << partial_sum.imag() << endl;
#endif
    }

#ifdef TP0
    f.close();
#endif

    return partial_sum;
}

#ifdef TP0
complex<double> p_analytique(const Point &P1, int N, const string &filename)
#else
complex<double> p_analytique(const Point &P1, int N)
#endif
{
#ifdef TP0
    ofstream f(filename);

    if (!f.is_open())
    {
        cout << "ERROR: Le fichier " << filename << " n'a pas pu être ouvert" << endl;
        exit(-1);
    }
#endif

    double theta = P1.theta();
    double ka = k * P1.norm();
    complex<double> iterative_i = 1;
    complex<double> partial_sum = k * boost::math::cyl_bessel_j(1, ka) - k * boost::math::cyl_bessel_j(0, ka) * hankel_n(ka, 1) / hankel_n(ka, 0);

#ifdef TP0
    f << 0 << " " << partial_sum.real() << " " << partial_sum.imag() << endl;
#endif

    for (int n = 1; n <= N; n++)
    {
        iterative_i *= -I;
        partial_sum += k * iterative_i * cos(n * theta) * (boost::math::cyl_bessel_j(n, ka) * (hankel_n(ka, n - 1) - hankel_n(ka, n + 1)) / hankel_n(ka, n) - (boost::math::cyl_bessel_j(n - 1, ka) - boost::math::cyl_bessel_j(n + 1, ka)));
#ifdef TP0
        f << n << " " << partial_sum.real() << " " << partial_sum.imag() << endl;
#endif
    }

#ifdef TP0
    f.close();
#endif

    return partial_sum;
}

#ifdef TP0
void exporte_solution_analytique(const string &filename, const double radius_obstacle, const unsigned int nbPasVisualisation, const double longueur, const int N)
{
    ofstream f(filename);

    if (!f.is_open())
    {
        cout << "ERROR: Le fichier " << filename << " n'a pas pu être ouvert" << endl;
        exit(-1);
    }
    string bin = "outputs/bin.txt";

    for (unsigned int i = 0; i < nbPasVisualisation; i++)
    {
        for (unsigned int j = 0; j < nbPasVisualisation; j++)
        {
            Point IJ(-1 * longueur + (i / (double)nbPasVisualisation) * 2 * longueur, -1 * longueur + (j / (double)nbPasVisualisation) * 2 * longueur);
            if (IJ.norm() >= radius_obstacle)
            {
                complex<double> valeur = u_N_plus_analytique(IJ, radius_obstacle, N, bin);
                f << IJ.x << " " << IJ.y << " " << valeur.real() << " " << valeur.imag() << endl;
                // cout<<IJ.x << " " << IJ.y << " " << valeur.real() << " " << valeur.imag() << endl;
            }
        }
    }
    f.close();

    return;
}
#endif

void export_obsctacles(const string &filename, vector<Cercle> cercles)
{
    ofstream f(filename);

    if (!f.is_open())
    {
        cout << "ERROR: Le fichier " << filename << " n'a pas pu être ouvert" << endl;
        exit(-1);
    }

    for (unsigned int i = 0; i < cercles.size(); i++)
    {
        const Cercle &s = cercles[i];
        f << s.centre.x << " " << s.centre.y << " ";
        f << s.rayon << endl;
    }
    f.close();

    return;
}

double erreur_relative(const std::complex<double> &u, const std::complex<double> &reference, double eps)
{
    const double denom = std::max(std::abs(reference), eps);
    return std::abs(u - reference) / denom;
}


#ifdef TP0
// ---------------------------------------------------------------------------
// Paste in utils.cpp (or test_functions.cpp) and declare in the matching header:
//
//   void etudie_convergence_airy(const vector<double> &kValues, double a, int Nmax,
//                                const vector<double> &tolValues, const string &outPrefix);
//
// Needs: #include <vector> #include <cmath> #include <iomanip> #include <algorithm>
//        #include <limits>  (hankel_n's overflow guard, if not already applied)
// ---------------------------------------------------------------------------

// Studies the truncation of the series for u+(a, theta=0) (the illuminated point, theta = 0,
// where convergence is slowest since all terms add constructively) as a function of k, for a
// fixed a, and for several target tolerances. For each k it:
//   1. builds the partial sums S_n, n = 0..Nmax, and the successive differences
//        diff_n = |S_n - S_{n-1}|
//      stopping as soon as a term becomes non-finite (the Jn/Yn overflow for n past ~ka),
//      which always happens well after the series has converged.
//   2. writes diff_n to one file per k (independent of tolerance):
//        outputs/airyResu_<outPrefix>_k=<k>_a=<a>.txt
//        header: "# k=.. a=.. ka=.. Neff=.."
//        rows  : n diff_n
//   3. for every tol in tolValues, finds n_star(tol) = smallest n such that diff_m < tol for
//      every m in [n, Neff] (so the series has truly settled below tol, not just dipped once),
//      and appends one row to a single summary file (long format: one row per (k, tol) pair):
//        outputs/airySummary_<outPrefix>.txt
//        header: "# k a ka tol n_star"
//        rows  : k a ka tol n_star        (n_star = -1 if never reached before Neff)
//
// The Airy-type estimate N_series ~= ka + C*(ka)^(1/3) and the best-fit C for each tolerance
// are NOT computed here: fit C in post-processing (plot_airy.py) from the (ka, n_star) pairs,
// since C should be fit from the data rather than assumed.
void etudie_convergence_airy(const vector<double> &kValues, double a, int Nmax,
                             const vector<double> &tolValues, const string &outPrefix)
{
    const string summaryName = "outputs/airySummary_" + outPrefix + ".txt";
    ofstream fsum(summaryName);
    if (!fsum.is_open())
    {
        cout << "ERROR: impossible d'ouvrir " << summaryName << endl;
        exit(-1);
    }
    fsum << setprecision(10);
    fsum << "# k a ka tol n_star\n";

    for (double k_local : kValues)
    {
        const double ka = k_local * a;
        const double theta = 0.0; // illuminated point: slowest convergence

        complex<double> iterative_i = 1;
        complex<double> S_prev = -(jn(0, ka) / hankel_n(ka, 0)) * hankel_n(ka, 0);
        complex<double> S = S_prev;

        const string fname = "outputs/airyResu_" + outPrefix + "_k=" + std::to_string(k_local)
                            + "_a=" + std::to_string(a) + ".txt";
        ofstream f(fname);
        if (!f.is_open())
        {
            cout << "ERROR: impossible d'ouvrir " << fname << endl;
            exit(-1);
        }
        f << setprecision(10);

        // Stop as soon as a term stops being finite: for n large enough (roughly n > ka),
        // Jn(ka) underflows towards 0 while Yn(ka) (inside hankel_n) overflows, so the term
        // becomes non-finite. This always happens well after the series has converged, but
        // well before Nmax if Nmax is set generously; once it happens every later diff would
        // spuriously fail the "stays converged" test below, so the search only looks at the
        // finite part, n = 1..Neff.
        vector<double> diffs(Nmax + 1, std::numeric_limits<double>::quiet_NaN()); // index 0 unused
        int Neff = 0;
        for (int n = 1; n <= Nmax; n++)
        {
            iterative_i *= -I;
            const complex<double> term = 2. * iterative_i * (jn(n, ka) / hankel_n(ka, n)) * hankel_n(ka, n) * cos(n * theta);
            if (!isfinite(term.real()) || !isfinite(term.imag()))
            {
                cout << "  (k=" << k_local << ": term becomes non-finite at n=" << n
                     << ", stopping the sum there)" << endl;
                break;
            }
            S -= term;
            diffs[n] = abs(S - S_prev);
            S_prev = S;
            Neff = n;
        }

        f << "# k=" << k_local << " a=" << a << " ka=" << ka << " Neff=" << Neff << "\n";
        for (int n = 1; n <= Neff; n++)
            f << n << " " << diffs[n] << "\n";
        f.close();

        for (double tol : tolValues)
        {
            int n_star = -1;
            for (int n = Neff; n >= 1; n--)
            {
                if (!(diffs[n] < tol))
                {
                    n_star = n + 1;
                    break;
                }
            }
            if (n_star == -1) n_star = 1;        // every diff was already below tol
            if (n_star > Neff) n_star = -1;        // never converged before the sum broke down

            cout << "k=" << k_local << " (ka=" << ka << "), tol=" << tol << ": n_star = " << n_star
                 << (n_star == -1 ? "  [not reached, increase Nmax]" : "") << endl;

            fsum << k_local << " " << a << " " << ka << " " << tol << " " << n_star << "\n";
        }
    }
    fsum.close();
}

#endif