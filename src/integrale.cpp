#include "headers/integrale.hpp"

complex<double> integ_simple(Segment &AB, fun_d_P f, int Nbpas1)
{
    complex<double> result = 0;
    double dAB = AB.norm() / Nbpas1;
    for (int i = 0; i < Nbpas1; i++)
    {
        result += dAB * f(AB.P1 + dAB * i * (AB.P2 - AB.P1));
    }
    return result;
}

complex<double> integ_simple_legendre_n(fun_double f, int n_ordre)
{
    const LegendreData &data = get_legendre_data(n_ordre);

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

complex<double> integ_simple_legendre_n(fun_double f, double a, double b, int n_ordre)
{
    const LegendreData &data = get_legendre_data(n_ordre);

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

complex<double> integ_simple_legendre_n(fun_d_P f, Segment &AB, int n_ordre)
{
    const LegendreData &data = get_legendre_data(n_ordre);

    const Point milieu = (AB.P1 + AB.P2) * 0.5;
    const Point vecteur_demi_longueur = (AB.P1 - AB.P2) * 0.5;

    complex<double> result = 0.0;

    for (std::size_t i = 0; i < data.roots.size(); ++i)
    {
        const double r = data.roots[i];
        const double w = data.weights[i];
        const Point correction = r * vecteur_demi_longueur;

        result += w * f(milieu + correction);

        if (r > 0.0)
            result += w * f(milieu - correction);
    }

    return vecteur_demi_longueur.norm() * result;
}

complex<double> integ_double(Segment &AB, Segment &CD, fun_d_P2 f, int Nbpas1, int Nbpas2)
{
    complex<double> result = 0;
    double dAB = AB.norm() / Nbpas1;
    double dCD = CD.norm() / Nbpas2;
    double dS = dAB * dCD;
    // double dS = abs((dAB * (AB.P2 - AB.P1)) * (dCD * (CD.P2 - CD.P1))); // produit vectoriel des deux petits vecteurs;
    for (int i = 0; i < Nbpas1; i++)
    {
        for (int j = 0; j < Nbpas2; j++)
        {
            result += dS * f(AB.P1 + dAB * (AB.P2 - AB.P1) * i, CD.P1 + dCD * (CD.P2 - CD.P1) * j);
        }
    }
    return result;
}

complex<double> log_norm(const Point &A, const Point &B)
{
    return log(k * (A - B).norm());
}

complex<double> integ_log(Segment &AB, Segment &CD, int Nbpas1, int Nbpas2)
{
    if (AB == CD)
    {
        return (AB.norm() * AB.norm()) * log(k * AB.norm()) - 3 * AB.norm() * AB.norm() / 2;
    }
    else if (AB.P1 == CD.P1 || AB.P1 == CD.P2 || AB.P2 == CD.P1 || AB.P2 == CD.P2)
    // test si il y a un sommet commun entre les 2 segments
    //  disjonction de cas necessaire ensuite pour revenir aux notations du sujet
    {
        Point Amoins;
        Point Aplus;
        Point Bmoins;
        Point Bplus;

        if (AB.P1 == CD.P1)
        { // cas  A=C
            Amoins = AB.P1;
            Aplus = AB.P2;
            Bmoins = CD.P2;
            Bplus = CD.P1;
        }
        else if (AB.P2 == CD.P1)
        { // cas  B=C
            Amoins = AB.P2;
            Aplus = AB.P1;
            Bmoins = CD.P2;
            Bplus = CD.P1;
        }
        else if (AB.P1 == CD.P2)
        { // cas A=D
            Amoins = AB.P1;
            Aplus = AB.P2;
            Bmoins = CD.P1;
            Bplus = CD.P2;
        }
        else // cas B=D ie AB.P2 == CD.P2
        {
            Amoins = AB.P2;
            Aplus = AB.P1;
            Bmoins = CD.P1;
            Bplus = CD.P2;
        }

        complex<double> sp = AB.norm() / 2;
        complex<double> sm = -AB.norm() / 2;
        complex<double> tp = CD.norm() / 2;
        complex<double> tm = -CD.norm() / 2;
        complex<double> sigma = 1 / AB.norm() * (((CD.P1 + CD.P2) / 2 - (AB.P1 + AB.P2) / 2) | (AB.P1 - AB.P2));
        complex<double> sigmam = 1 / AB.norm() * ((CD.P2 - (AB.P1 + AB.P2) / 2) | (AB.P1 - AB.P2));
        complex<double> to = 1 / CD.norm() * (((AB.P1 + AB.P2) / 2 - (CD.P1 + CD.P2) / 2) | (CD.P1 - CD.P2));
        complex<double> top = 1 / CD.norm() * ((AB.P1 - (CD.P1 + CD.P2) / 2) | (CD.P1 - CD.P2));
        complex<double> cost = (AB.P1 - AB.P2) | (CD.P1 - CD.P2) / AB.norm() / CD.norm();
        complex<double> sint = (AB.P1 - AB.P2) * (CD.P1 - CD.P2) / AB.norm() / CD.norm();
        complex<double> dalphap = tp * sint - (to + sigma * cost) / sint;
        complex<double> dbetap = sp * sint - (to + sigma * cost) / sint;
        complex<double> terme1 = AB.norm() / 2 * ((tp - top) * log((Bplus - Aplus).norm()) - (tm - top) * log((Bmoins - Aplus).norm()));
        complex<double> terme2 = CD.norm() / 2 * ((sp - sigmam) * log((Aplus - Bmoins).norm()) - (sm - sigmam) * log((Amoins - Bmoins).norm()));
        complex<double> terme3 = dbetap * (atan((tp - top) / dbetap) - atan((tm - top) / dbetap));
        complex<double> terme4 = dalphap * (atan((sp - sigmam) / dalphap) - atan((sm - sigmam) / dalphap));
        complex<double> terme5 = (log(k) - 3. / 2) * AB.norm() * CD.norm(); // terme en k apparait ici
        return terme1 + terme2 + terme3 + terme4 + terme5;
    }
    else
    {
        return (integ_double(AB, CD, log_norm, Nbpas1, Nbpas2));
    }
}

complex<double> integrale_pour_p(Segment &AB, fun_d_P f, int Nbpas1, const Point &x)
{
    complex<double> result = 0;
    double dAB = AB.norm() / (double)Nbpas1;
    // cout<<"dAB"<<dAB<<endl;
    for (int i = 0; i < Nbpas1; i++)
    {
        Point y = AB.P1 + dAB * i * (AB.P2 - AB.P1);
        Point normale(AB.P2.y - AB.P1.y, AB.P1.x - AB.P2.x);
        double facteur_normale = (x - y) | normale;
        result += p_theta_config((x - y).x, (x - y).y) * f(x - y) * facteur_normale / (x - y).norm();
        // cout<<"resultat inter integ p"<<result<<endl;
    }
    return result;
}

complex<double> integ_simple_pour_p(Segment &AB, fun_d_P f, int Nbpas1, const Point &x)
{
    complex<double> result = 0;
    double dAB = AB.norm() / Nbpas1;
    for (int i = 0; i < Nbpas1; i++)
    {
        result += dAB * f(x - (AB.P1 + dAB * i * (AB.P2 - AB.P1)));
    }
    return result;
}
