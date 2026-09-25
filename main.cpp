#include <iostream>
#include <time.h>
#include "src/config/config.hpp"
#include "src/config/constantes.hpp"
#include "src/config/external.hpp"
#include "src/headers/utils.hpp"

int main()
{
    // initialisation du random
    srand(time(NULL));
    const string config = "config.txt";
    get_config(config);

    if (effectuerTests)
    {
        // tests Point
        Point M(1, 1);
        Point C(1, 2);
        std::cout << "A = (" << M.x << "," << M.y << ")" << endl;
        std::cout << "A =" << M << endl;
        std::cout << "M+C = " << M + C;
        std::cout << "M-C = " << M - C;
        std::cout << "3*C = " << 3 * C;
        std::cout << "C/5 = " << C / 5;
        std::cout << "norme C = " << C.norm() << endl;
        std::cout << "produit vectoriel composante z = " << M * C << endl;

        // tests Cercle
        Cercle B(10, M);
        std::cout << "obstacle point " << B.centre.x << "," << B.centre.y << endl;
        std::cout << " obstacle rayon = " << B.rayon << endl;

        // tests Segment
        Segment S(M, C);
        std::cout << "segment P1 " << S.P1.x << "," << S.P1.y << endl;
        std::cout << "segment P2 " << S.P2.x << "," << S.P2.y << endl;

        // tests Maillage
        vector<Cercle> obstaclesTest;
        Point O(0, 0);
        double a = 1.;
        Cercle Cercle0(a, O);
        obstaclesTest.push_back(Cercle0);
        Maillage monMaillage;
        std::cout << " k= " << k << ", a =" << a << " N =" << N << endl;
        monMaillage.ajoute_cercle(pasMaillage, Cercle0);
        monMaillage.export_maillage("outputs/test.txt");
        std::cout << "exported maillage" << endl;
    }

    if (effectuerLaSimulation)
    {
#ifdef TP0
        /* ----- TP0 ----- */

        // -- Creation Maillage --
        vector<Cercle> obstaclesTest;
        Point O(0, 0);
        double a = 1.;
        Cercle Cercle0(a, O);
        obstaclesTest.push_back(Cercle0);
        Maillage monMaillage;
        std::cout << "Info Générales: fréquence k = " << k << ", Taille du maillage a = " << a << ", Pas du maillage = " << pasMaillage << endl;
        monMaillage.ajoute_cercle(pasMaillage, Cercle0);
        std::cout << "Export du maillage..." << endl;
        monMaillage.export_maillage(cheminFichierMaillage);

        // -- Sommes partielles U_n^+ , q et p --
        std::cout << "Sommes partielles..." << endl;

        vector<Point> pointsCercle = monMaillage.PointsMaillage();
        std::cout << "Nombre de point du maillage = " << pointsCercle.size() << endl;

        std::cout << "Création des fichiers..." << endl;
        std::string filename_qConverged = "outputs/qResu_k=" + std::to_string(k) + "_a=" + std::to_string(a) + "_N=" + std::to_string(N) + ".txt";
        ofstream qC(filename_qConverged);
        std::string filename_pConverged = "outputs/pResu_k=" + std::to_string(k) + "_a=" + std::to_string(a) + "_N=" + std::to_string(N) + ".txt";
        ofstream pC(filename_pConverged);

        if (!qC.is_open())
        {
            std::cout << "ERROR: Le fichier " << filename_qConverged << " n'a pas pu être ouvert" << endl;
            exit(-1);
        }
        if (!pC.is_open())
        {
            std::cout << "ERROR: Le fichier " << filename_pConverged << " n'a pas pu être ouvert" << endl;
            exit(-1);
        }

        qC << "k =" << k << " a=" << a << endl;
        pC << "k =" << k << " a=" << a << endl;

        for (unsigned int j = 0; j < 10; j++)
        {
            int i = j * pointsCercle.size() / 10;
            Point P = pointsCercle[i];
            double theta = P.theta();
            std::string filename_q = "outputs/qResu_k=" + std::to_string(k) + "_a=" + std::to_string(a) + "_theta=" + std::to_string(theta) + "_N=" + std::to_string(N) + ".txt";
            std::string filename_p = "outputs/pResu_k=" + std::to_string(k) + "_a=" + std::to_string(a) + "_theta=" + std::to_string(theta) + "_N=" + std::to_string(N) + ".txt";
            std::string filename_solExt = "outputs/solExtResu_k=" + std::to_string(k) + "_a=" + std::to_string(a) + "_theta=" + std::to_string(theta) + "_N=" + std::to_string(N) + ".txt";
            complex<double> q_approche = q_analytique(P, N, filename_q);
            complex<double> p_approche = p_analytique(P, N, filename_p);

            u_N_plus_analytique(P, a, N, filename_solExt); // complex<double> sol_ext =
            qC << P << " " << theta << " " << q_approche.real() << " " << q_approche.imag() << endl;
            pC << P << " " << theta << " " << p_approche.real() << " " << p_approche.imag() << endl;
        }

        double err = 0;
        for (int t = 0; t < 10; t++)
        {
            double th = pi * t / 9.0;
            Point P(3. * a * cos(th), 3. * a * sin(th));
            std::string bin_theta = "outputs/solExtFarResu_k=" + std::to_string(k) + "_a=" + std::to_string(a) + "_theta=" + std::to_string(th) + "_N=" + std::to_string(N) + ".txt";
            complex<double> u = u_N_plus_analytique(P, a, N, bin_theta);
            err = max(err, abs(u + exp(-I * k * a * cos(th)))); // +: u = -uinc
        }
        std::cout << "N=" << N << " boundary error = " << err << endl;

        const string sol_externe = "outputs/u_N_plusResu_k=" + std::to_string(k) + "_a=" + std::to_string(a) + "_N=" + std::to_string(N) + ".txt";
        export_obsctacles(cheminFichierObstacles, obstaclesTest);
        exporte_solution_analytique(sol_externe, a, nbPasExport, 10, N);

        double Hdiff;
        const unsigned int nbPoints = 50;
        const bool centered = true;

        for (int i = 1; i <= 8; i++)
        {
            const string filename_q = "outputs/qFD_k=" + std::to_string(k) + "_a=" + std::to_string(a) + "_N=" + std::to_string(N) + "_h=" + std::to_string(Hdiff) + ".txt";
            const string filename_p = "outputs/pFD_k=" + std::to_string(k) + "_a=" + std::to_string(a) + "_N=" + std::to_string(N) + "_h=" + std::to_string(Hdiff) + ".txt";

            Hdiff = pow(10., -i);

            ofstream fq(filename_q);
            ofstream fp(filename_p);
            if (!fq.is_open() || !fp.is_open())
            {
                std::cout << "ERROR: impossible d'ouvrir " << filename_q << " ou " << filename_p << endl;
                exit(-1);
            }
            const string bin = "outputs/bin.txt"; // scratch file required by u_N_plus, q and p

            for (ofstream *f : {&fq, &fp})
            {
                *f << setprecision(17);
                *f << "# k=" << k << " a=" << a << " N=" << N << " h=" << Hdiff
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
                { return u_N_plus_analytique(Point(r * cos(th), r * sin(th)), a, N, bin); };
                auto u_tot = [&](double r) -> complex<double>
                { return u_scat(r) + exp(-I * k * r * cos(th)); };
                auto deriv = [&](const function<complex<double>(double)> &f) -> complex<double>
                {
                    if (centered)
                        return (f(a + Hdiff) - f(a - Hdiff)) / (2. * Hdiff);
                    return (f(a + Hdiff) - f(a)) / (Hdiff);
                };

                const complex<double> q_fd = deriv(u_scat);
                const complex<double> p_fd = -deriv(u_tot);
                const complex<double> q_ex = q_analytique(P, N, bin);
                const complex<double> p_ex = p_analytique(P, N, bin);

                fq << th << " " << q_fd.real() << " " << q_fd.imag() << " " << q_ex.real() << " "
                   << q_ex.imag() << " " << abs(q_fd - q_ex) << "\n";
                fp << th << " " << p_fd.real() << " " << p_fd.imag() << " " << p_ex.real() << " "
                   << p_ex.imag() << " " << abs(p_fd - p_ex) << "\n";
                maxErrQ = max(maxErrQ, abs(q_fd - q_ex));
                maxErrP = max(maxErrP, abs(p_fd - p_ex));
            }
            std::cout << "FD check (h=" << Hdiff << "): max |q_fd - q| = " << maxErrQ
                      << ", max |p_fd - p| = " << maxErrP << endl;
        }
#endif

#ifdef TP1
        /* ----- TP1 ----- */

        // -- Parametres --
        const double rayon = 1.;
        const double delta = 0.2;
        const double pas = 0.001;
        const unsigned int nbPointSolution = 10000;
        const unsigned int idxTroncature = 25;
        const unsigned int ordre = 4;

        // -- Création maillage --
        Point O(0, 0);
        std::cout << "Info Maillage: centre du cercle = " << O << ", rayon du cercle = " << rayon << ", pas du maillage = " << pas << endl;
        Cercle cercle(rayon, O);
        Maillage maillage;
        maillage.ajoute_cercle(pas, cercle);
        const unsigned int nbSegmentsMaillage = maillage.size();

        // -- Calcul du vecteur de p sur les milieux des bords --
        vector<complex<double>> vect_p;
        for (unsigned int i = 0; i < maillage.size(); i++)
        {
            vect_p.push_back(p_analytique(maillage[i].milieu(), idxTroncature));
        }

        // -- Créqtion du fichier de résultat --
        const string filename = "outputs/u_sol_k=" + std::to_string(k) + "_rayon=" + std::to_string(rayon) + +"_pas=" + std::to_string(pas) + "_nbPtSol=" + std::to_string(nbPointSolution) + "_idxTroncature=" + std::to_string(idxTroncature) + "_ordre=" + std::to_string(ordre) + "_Delta=" + std::to_string(delta) + ".txt";
        ofstream file(filename);

        if (!file.is_open())
        {
            std::cout << "ERROR: Le fichier " << filename << " n'a pas pu être ouvert" << endl;
            exit(-1);
        }

        // -- Construction solution approchée--
        double rayon_sol = rayon + delta;
        complex<double> uPj;
        double theta;

        for (unsigned int j = 0; j < nbPointSolution; j++)
        {
            uPj = 0.0;
            theta = 2.0 * pi * static_cast<double>(j) / static_cast<double>(nbPointSolution);

            const Point Pj{rayon_sol * std::cos(theta), rayon_sol * std::sin(theta)};

            fun_d_P f = [&Pj](const Point &Q)
            {
                return green(Pj, Q);
            };

            for (unsigned int i = 0; i < nbSegmentsMaillage; i++)
            {
                uPj += integ_simple(f, maillage[i], ordre) * vect_p[i];
            }

            file << Pj.x << " " << Pj.y << " " << uPj.real() << " " << uPj.imag() << endl;
        }

        file.close();

#endif
    }

    return 0;
}
