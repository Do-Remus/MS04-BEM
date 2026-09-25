#include "headers/utils.hpp"

complex<double> hankel_n(const double x, const int n)
{
    return boost::math::cyl_bessel_j(n, x) + I * boost::math::cyl_neumann(n, x);
}

complex<double> green(const Point &p1, const Point &p2)
{
    return (I / 4.0) * hankel_n(k * (p1 - p2).norm(), 0);
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
