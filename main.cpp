#include <iostream>
#include <time.h>
#include "src/config/config.hpp"
#include "src/config/constantes.hpp"
#include "src/config/external.hpp"
#include "src/headers/utils.hpp"

Vecteur produit_A_cached(const Maillage &maillage, const Vecteur &x, const LegendreData &legendreData)
{
    const unsigned int N = maillage.size();
    Vecteur y(N, 0.0);

    for (unsigned int i = 0; i < N; ++i)
    {
        // Diagonale A_ii
        {
            const Segment S = maillage[i];

            std::complex<double> aii =
                integ_simple(
                    [&S](const Point &Q, double)
                    {
                        return integ_simple_log(Q, S);
                    },
                    S, legendreData);

            aii += integ_double(
                [](const Point &Q, const Point &P, double)
                {
                    return green_reguliere_cached_vec(Q, P);
                },
                S, S, legendreData);

            y[i] += aii * x[i];
        }

        // Triangle supérieur : A_ij = A_ji
        for (unsigned int j = i + 1; j < N; ++j)
        {
            const std::complex<double> aij =
                integ_double(
                    [](const Point &Q, const Point &P, double)
                    {
                        return green_cached_vec(Q, P);
                    },
                    maillage[i], maillage[j], legendreData);

            // Contribution de A_ij * x_j à y_i
            y[i] += aij * x[j];

            // Contribution de A_ji * x_i à y_j
            y[j] += aij * x[i];
        }
    }

    return y;
}

Vecteur gradConjMatrixFree(const Maillage &maillage, const Vecteur &b, const LegendreData &legendreData, double tol, unsigned int maxIter)
{
    const std::size_t n = b.size();

    Vecteur x(n, 0.0);

    Vecteur r = b;
    Vecteur d = r;

    const double bnorm = b.norm();
    double relativeResidual = 0.0;

    if (bnorm == 0.0)
        return x;

    std::complex<double> rho = r.produitBilineaire(r);

    for (unsigned int iter = 0; iter < maxIter; ++iter)
    {
        Vecteur Ad = produit_A_cached(maillage, d, legendreData);
        const std::complex<double> denom = d.produitBilineaire(Ad);
        const double scale = d.norm() * Ad.norm();

        if (scale == 0 || std::abs(denom) < 1e-20 * scale)
        {
            std::cerr << "GCMS : breakdown, denominateur nul"
                      << " a l'iteration " << iter
                      << ", avec un résidu relatif de : "
                      << relativeResidual << std::endl;
            return x;
        }

        const std::complex<double> alpha = rho / denom;
        x = x + alpha * d;
        r = r - alpha * Ad;
        relativeResidual = r.norm() / bnorm;

        if (relativeResidual < tol)
        {
            std::cout << "GCMS converge en "
                      << iter + 1
                      << " iterations avec un résidu relatif de : "
                      << relativeResidual << std::endl;

            return x;
        }

        const std::complex<double> rhoNew = r.produitBilineaire(r);

        if (std::abs(rho) < 1e-30)
        {

            std::cerr << "COCG : breakdown de rho, avec résidu relatif de : "
                      << relativeResidual
                      << std::endl;
            return x;
        }

        const std::complex<double> beta = rhoNew / rho;

        d = r + beta * d;
        rho = rhoNew;
    }

    std::cout << "GCMS : nombre maximal d'iterations atteint."
              << std::endl;

    return x;
}

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

        // Tests Matrice
        Matrice MMM(3, 3);
        cout << MMM << endl;
        MMM(1, 1) = 1;
        MMM(1, 2) = 2;
        MMM(3, 2) = 3;
        cout << MMM << endl;

        Matrice BBB(3, 2);
        cout << BBB << endl;
        BBB(3, 2) = 1;
        BBB(1, 2) = 2;
        BBB(2, 1) = 0;
        BBB(1, 1) = 3;
        cout << BBB << endl;

        Vecteur v(2);
        v[0] = 1;
        v[1] = 1;
        cout << BBB * v << endl;

        // ===================================================================
        // Test 1 : erreur de quadrature de Gauss-Legendre (ordre n=10) sur
        // les monomes de la base canonique x^p, p = 0..20, integres sur [-1,1]
        // (intervalle "naturel" des racines de Legendre)
        // ===================================================================
        {
            const string filename = string("outputs/convergence_quadrature_monome.txt");
            ofstream file(filename);

            if (!file.is_open())
            {
                std::cout << "ERROR: Le fichier " << filename << " n'a pas pu être ouvert" << endl;
                exit(-1);
            }

            const unsigned int ordreQuadrature = 10;
            const LegendreData &dataMonome = get_legendre_data(ordreQuadrature);

            file << "# Erreur de la quadrature de Gauss-Legendre (n=" << ordreQuadrature
                 << ") sur les monomes x^p integres sur [-1,1]" << endl;
            file << "# degre_p erreur_absolue" << endl;

            for (int p = 0; p <= 20; p++)
            {
                auto monome = [p](double x) -> complex<double>
                { return std::pow(x, p); };

                const complex<double> approx = integ_simple(monome, dataMonome);

                // integrale exacte de x^p sur [-1,1] : 0 si p impair, 2/(p+1) si p pair
                const double exact = (p % 2 == 0) ? (2.0 / (p + 1)) : 0.0;

                const double erreur = std::abs(approx - exact);

                file << p << " " << erreur << endl;
            }

            file.close();
            std::cout << "Fichier " << filename << " exporte (degre d'exactitude attendu = "
                      << 2 * ordreQuadrature - 1 << ")" << endl;
        }

        // ===================================================================
        // Test 2 : convergence de la quadrature de Gauss-Legendre vers
        // l'integrale de la fonction de Hankel H_0^(1)(x) = J_0(x) + i*Y_0(x)
        // sur trois segments a differentes distances de l'origine, afin
        // d'observer si la convergence est plus rapide en champ lointain.
        // ===================================================================
        {
            const string filename = string("outputs/convergence_quadrature_hankel.txt");
            ofstream file(filename);

            if (!file.is_open())
            {
                std::cout << "ERROR: Le fichier " << filename << " n'a pas pu être ouvert" << endl;
                exit(-1);
            }

            // fonction a integrer : H_0^(1)(x)
            auto H0 = [](double x) -> complex<double>
            { return hankel_n(x, 0); };

            struct SegmentTest
            {
                double a;
                double b;
            };

            const vector<SegmentTest> segments = {
                {0.01, 0.02},  // champ proche
                {1.01, 1.02},  // champ intermediaire
                {10.01, 10.02} // champ lointain
            };

            const unsigned int ordreMax = 20;       // ordres de quadrature testes : 1..ordreMax
            const unsigned int ordreReference = 40; // ordre pris comme reference "exacte"

            // valeurs de reference (quadrature d'ordre eleve)
            const LegendreData &dataRef = get_legendre_data(ordreReference);
            vector<complex<double>> references;
            for (const auto &s : segments)
                references.push_back(integ_simple(H0, s.a, s.b, dataRef));

            file << "# Convergence de la quadrature de Gauss-Legendre vers l'integrale de H_0^(1)"
                 << " sur differents segments (reference = ordre " << ordreReference << ")" << endl;
            file << "# ordre_n";
            for (const auto &s : segments)
                file << " erreur_[" << s.a << "," << s.b << "]";
            file << endl;

            for (unsigned int n = 1; n <= ordreMax; n++)
            {
                const LegendreData &data = get_legendre_data(n);

                file << n;

                for (std::size_t i = 0; i < segments.size(); i++)
                {
                    const complex<double> approx = integ_simple(H0, segments[i].a, segments[i].b, data);
                    const double erreur = std::abs(approx - references[i]);
                    file << " " << erreur;
                }
                file << endl;
            }

            file.close();
            std::cout << "Fichier " << filename << " exporte." << endl;
        }

        // ===================================================================
        // Test 3 : erreur de u+ (BEM, quadrature + Green non cachee) par
        // rapport a u_analytique, sur un cercle test de rayon superieur a
        // la frontiere, en fonction de l'angle theta, pour differents
        // ordres de quadrature n_q (pas de maillage FIXE).
        // ===================================================================
        {
            const string filename = string("outputs/erreur_u_theta_vs_nq.txt");
            ofstream file(filename);

            if (!file.is_open())
            {
                std::cout << "ERROR: Le fichier " << filename << " n'a pas pu être ouvert" << endl;
                exit(-1);
            }

            // -- Parametres fixes --
            const double kTest = 10.0;
            k = kTest; // k est une variable globale utilisee par green(), p_analytique(), u_N_plus_analytique()

            const double rayonTest = 1.0;       // rayon de la frontiere (cercle diffractant)
            const double pasMaillageTest = 0.1; // pas de maillage FIXE pour ce test
            const unsigned int idxTroncatureTest = 25;

            const double R_test = 2.0 * rayonTest; // rayon du cercle test > rayon de la frontiere
            const unsigned int nbThetaTest = 100;

            const vector<unsigned int> nq_values = {1, 2, 3, 4, 5, 6, 8, 10, 15, 20};

            // -- Construction du maillage frontiere (fixe pour ce test) --
            Point O(0, 0);
            Cercle cercleTest(rayonTest, O);
            Maillage maillageTest;
            maillageTest.ajoute_cercle(pasMaillageTest, cercleTest);
            const unsigned int nbSegments = maillageTest.size();

            // -- Donnee de Neumann p sur le maillage (independante de n_q) --
            vector<complex<double>> vect_p_test;
            for (unsigned int i = 0; i < nbSegments; i++)
                vect_p_test.push_back(p_analytique(maillageTest[i].milieu, idxTroncatureTest));

            // -- Points du cercle test et solution analytique associee --
            vector<Point> pointsCercleTest;
            vector<double> thetasTest;
            vector<complex<double>> uExactTest;
            for (unsigned int j = 0; j < nbThetaTest; j++)
            {
                const double theta = 2.0 * pi * j / nbThetaTest;
                Point Pj(R_test * cos(theta), R_test * sin(theta));
                pointsCercleTest.push_back(Pj);
                thetasTest.push_back(theta);
                uExactTest.push_back(u_N_plus_analytique(Pj, rayonTest, idxTroncatureTest));
            }

            // -- En-tete du fichier --
            file << "# Erreur de u+ (BEM) vs u_analytique sur cercle test de rayon R="
                 << R_test << ", en fonction de l'ordre de quadrature n_q" << endl;
            file << "# k=" << kTest << " rayon=" << rayonTest << " pasMaillage=" << pasMaillageTest
                 << " idxTroncature=" << idxTroncatureTest << " R_test=" << R_test
                 << " nbThetaTest=" << nbThetaTest << endl;
            file << "# n_q_values:";
            for (unsigned int nq : nq_values)
                file << " " << nq;
            file << endl;
            file << "# theta";
            for (unsigned int nq : nq_values)
                file << " erreur_nq" << nq;
            file << endl;

            // -- Calcul de l'erreur pour chaque n_q --
            for (unsigned int j = 0; j < nbThetaTest; j++)
            {
                const Point &Pj = pointsCercleTest[j];

                file << thetasTest[j];

                for (unsigned int nq : nq_values)
                {
                    const LegendreData &data = get_legendre_data(nq);

                    complex<double> uApprox = 0.0;
                    for (unsigned int i = 0; i < nbSegments; i++)
                    {
                        uApprox += vect_p_test[i] * integ_simple(
                                                        [&Pj](const Point &Q, double)
                                                        { return green(Pj, Q); },
                                                        maillageTest[i],
                                                        data);
                    }

                    const double erreur = std::abs(uApprox - uExactTest[j]);
                    file << " " << erreur;
                }

                file << endl;
            }

            file.close();
            std::cout << "Fichier " << filename << " exporte." << endl;
        }

        // ===================================================================
        // Test 4 : erreur de u+ (BEM) par rapport a u_analytique, sur le
        // meme type de cercle test, en fonction de l'angle theta, pour
        // differents pas de maillage de la frontiere (n_q FIXE).
        // ===================================================================
        {
            const string filename = string("outputs/erreur_u_theta_vs_pasMaillage.txt");
            ofstream file(filename);

            if (!file.is_open())
            {
                std::cout << "ERROR: Le fichier " << filename << " n'a pas pu être ouvert" << endl;
                exit(-1);
            }

            // -- Parametres fixes --
            const double kTest = 10.0;
            k = kTest;

            const double rayonTest = 1.0;
            const unsigned int idxTroncatureTest = 25;
            const unsigned int nqFixe = 8; // ordre de quadrature FIXE pour ce test

            const double R_test = 2.0 * rayonTest;
            const unsigned int nbThetaTest = 100;

            const vector<double> pas_values = {0.2, 0.1, 0.05, 0.025, 0.01, 0.001, 0.0001};

            const LegendreData &dataFixe = get_legendre_data(nqFixe);

            // -- Points du cercle test et solution analytique (independants du maillage) --
            vector<Point> pointsCercleTest;
            vector<double> thetasTest;
            vector<complex<double>> uExactTest;
            for (unsigned int j = 0; j < nbThetaTest; j++)
            {
                const double theta = 2.0 * pi * j / nbThetaTest;
                Point Pj(R_test * cos(theta), R_test * sin(theta));
                pointsCercleTest.push_back(Pj);
                thetasTest.push_back(theta);
                uExactTest.push_back(u_N_plus_analytique(Pj, rayonTest, idxTroncatureTest));
            }

            // -- En-tete du fichier --
            file << "# Erreur de u+ (BEM) vs u_analytique sur cercle test de rayon R="
                 << R_test << ", en fonction du pas de maillage de la frontiere" << endl;
            file << "# k=" << kTest << " rayon=" << rayonTest << " n_q=" << nqFixe
                 << " idxTroncature=" << idxTroncatureTest << " R_test=" << R_test
                 << " nbThetaTest=" << nbThetaTest << endl;
            file << "# pasMaillage_values:";
            for (double h : pas_values)
                file << " " << h;
            file << endl;
            file << "# theta";
            for (double h : pas_values)
                file << " erreur_h" << h;
            file << endl;

            // -- Pour chaque pas de maillage : construction du maillage + calcul de l'erreur --
            vector<vector<double>> erreurs(pas_values.size(), vector<double>(nbThetaTest, 0.0));

            for (std::size_t ih = 0; ih < pas_values.size(); ih++)
            {
                Point O(0, 0);
                Cercle cercleTest(rayonTest, O);
                Maillage maillageTest;
                maillageTest.ajoute_cercle(pas_values[ih], cercleTest);
                const unsigned int nbSegments = maillageTest.size();

                vector<complex<double>> vect_p_test;
                for (unsigned int i = 0; i < nbSegments; i++)
                    vect_p_test.push_back(p_analytique(maillageTest[i].milieu, idxTroncatureTest));

                for (unsigned int j = 0; j < nbThetaTest; j++)
                {
                    const Point &Pj = pointsCercleTest[j];

                    complex<double> uApprox = 0.0;
                    for (unsigned int i = 0; i < nbSegments; i++)
                    {
                        uApprox += vect_p_test[i] * integ_simple(
                                                        [&Pj](const Point &Q, double)
                                                        { return green(Pj, Q); },
                                                        maillageTest[i],
                                                        dataFixe);
                    }

                    erreurs[ih][j] = std::abs(uApprox - uExactTest[j]);
                }
            }

            // -- Ecriture --
            for (unsigned int j = 0; j < nbThetaTest; j++)
            {
                file << thetasTest[j];
                for (std::size_t ih = 0; ih < pas_values.size(); ih++)
                    file << " " << erreurs[ih][j];
                file << endl;
            }

            file.close();
            std::cout << "Fichier " << filename << " exporte." << endl;
        }

        // ===================================================================
        // Test 5 : erreur de quadrature PURE. Contrairement au Test 3, on ne
        // compare plus a u_analytique (qui melange erreur de maillage et
        // erreur de quadrature) mais a une reference calculee SUR LE MEME
        // MAILLAGE avec un ordre de quadrature tres eleve. Cela isole l'effet
        // de n_q seul, sans compensation fortuite avec l'erreur de maillage.
        // ===================================================================
        {
            const string filename = string("outputs/erreur_u_theta_vs_nq_quadrature_pure.txt");
            ofstream file(filename);

            if (!file.is_open())
            {
                std::cout << "ERROR: Le fichier " << filename << " n'a pas pu être ouvert" << endl;
                exit(-1);
            }

            // -- Parametres fixes (memes que le Test 3) --
            const double kTest = 10.0;
            k = kTest;

            const double rayonTest = 1.0;
            const double pasMaillageTest = 0.1; // pas de maillage FIXE pour ce test
            const unsigned int idxTroncatureTest = 25;

            const double R_test = 2.0 * rayonTest;
            const unsigned int nbThetaTest = 100;

            const vector<unsigned int> nq_values = {1, 2, 3, 4, 5, 6, 8, 10, 15, 20};
            const unsigned int nqReference = 30; // ordre "exact" pour ce meme maillage

            // -- Construction du maillage frontiere (fixe, identique pour tous les n_q) --
            Point O(0, 0);
            Cercle cercleTest(rayonTest, O);
            Maillage maillageTest;
            maillageTest.ajoute_cercle(pasMaillageTest, cercleTest);
            const unsigned int nbSegments = maillageTest.size();

            // -- Donnee de Neumann p sur le maillage (independante de n_q) --
            vector<complex<double>> vect_p_test;
            for (unsigned int i = 0; i < nbSegments; i++)
                vect_p_test.push_back(p_analytique(maillageTest[i].milieu, idxTroncatureTest));

            // -- Points du cercle test --
            vector<Point> pointsCercleTest;
            vector<double> thetasTest;
            for (unsigned int j = 0; j < nbThetaTest; j++)
            {
                const double theta = 2.0 * pi * j / nbThetaTest;
                pointsCercleTest.emplace_back(R_test * cos(theta), R_test * sin(theta));
                thetasTest.push_back(theta);
            }

            // -- Fonction utilitaire : calcule u+ (BEM) sur ce maillage pour un ordre n_q donne --
            auto calcule_u_BEM = [&](const LegendreData &data) -> vector<complex<double>>
            {
                vector<complex<double>> u(nbThetaTest, 0.0);
                for (unsigned int j = 0; j < nbThetaTest; j++)
                {
                    const Point &Pj = pointsCercleTest[j];
                    complex<double> uApprox = 0.0;
                    for (unsigned int i = 0; i < nbSegments; i++)
                    {
                        uApprox += vect_p_test[i] * integ_simple(
                                                        [&Pj](const Point &Q, double)
                                                        { return green(Pj, Q); },
                                                        maillageTest[i],
                                                        data);
                    }
                    u[j] = uApprox;
                }
                return u;
            };

            // -- Reference : meme maillage, ordre de quadrature tres eleve --
            const LegendreData &dataReference = get_legendre_data(nqReference);
            const vector<complex<double>> uReference = calcule_u_BEM(dataReference);

            // -- En-tete du fichier --
            file << "# Erreur de quadrature PURE : u+ (BEM, ordre n_q) vs u+ (BEM, ordre nqReference="
                 << nqReference << ") sur le MEME maillage (pasMaillage=" << pasMaillageTest << ")" << endl;
            file << "# k=" << kTest << " rayon=" << rayonTest << " pasMaillage=" << pasMaillageTest
                 << " idxTroncature=" << idxTroncatureTest << " R_test=" << R_test
                 << " nbThetaTest=" << nbThetaTest << " nqReference=" << nqReference << endl;
            file << "# n_q_values:";
            for (unsigned int nq : nq_values)
                file << " " << nq;
            file << endl;
            file << "# theta";
            for (unsigned int nq : nq_values)
                file << " erreur_nq" << nq;
            file << endl;

            // -- Calcul de l'erreur pure pour chaque n_q --
            vector<vector<double>> erreurs(nq_values.size(), vector<double>(nbThetaTest, 0.0));

            for (std::size_t iq = 0; iq < nq_values.size(); iq++)
            {
                const LegendreData &data = get_legendre_data(nq_values[iq]);
                const vector<complex<double>> uApprox = calcule_u_BEM(data);

                for (unsigned int j = 0; j < nbThetaTest; j++)
                    erreurs[iq][j] = std::abs(uApprox[j] - uReference[j]);
            }

            // -- Ecriture --
            for (unsigned int j = 0; j < nbThetaTest; j++)
            {
                file << thetasTest[j];
                for (std::size_t iq = 0; iq < nq_values.size(); iq++)
                    file << " " << erreurs[iq][j];
                file << endl;
            }

            file.close();
            std::cout << "Fichier " << filename << " exporte." << endl;
        }
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

        // -- Parametres du problème --
        k = 10;

        // -- Parametres Maillage --
        const double rayon = 1.; // rayon cercle du maillage
        pasMaillage = 0.005;     // pas du maillage

        // -- Parametres Solution --
        const double L = 4.;              // Domaine LxL pour le calcul de la solution
        const double pasSolution = 0.1;   // pas de la solution
        const double delta = pasSolution; // distance minimale entre les points de solution et le cercle

        // -- Parametre  d'approximation --
        const unsigned int idxTroncature = 25;
        const unsigned int ordre = 4;
        green_cache_step = 0.1 * min(pasMaillage, pasSolution) / k;
        max_index = static_cast<unsigned int>(std::ceil(L * std::sqrt(2.0) / green_cache_step));

        std::cout << "\n=== Parametres du probleme ===" << std::endl;
        std::cout << "  Nombre d'onde k          : " << k << std::endl;

        std::cout << "\n  --- Approximation ---" << std::endl;
        std::cout << "  Indice de troncature N   : " << idxTroncature << std::endl;
        std::cout << "  Ordre de quadrature      : " << ordre << std::endl;
        std::cout << "  Pas du cache de Green    : " << green_cache_step << std::endl;
        std::cout << "  Nombre max d'indices     : " << max_index << std::endl;

        // -- Start time --
        std::clock_t start = std::clock();

        // -- Création maillage --
        Point O(0, 0);
        Cercle cercle(rayon, O);
        Maillage maillage;
        maillage.ajoute_cercle(pasMaillage, cercle);
        const unsigned int nbSegmentsMaillage = maillage.size();

        std::cout << "\n=== Maillage du cercle ===" << std::endl;
        std::cout << "  Centre              : " << O << std::endl;
        std::cout << "  Rayon               : " << rayon << std::endl;
        std::cout << "  Pas                 : " << pasMaillage << std::endl;
        std::cout << "  Nombre de segments  : " << nbSegmentsMaillage << std::endl;

        // -- Maillage pour la solution --
        vector<Point> pointsSolution;

        const double xmin = -L / 2.0;
        const double xmax = L / 2.0;
        const double ymin = -L / 2.0;
        const double ymax = L / 2.0;

        const double distanceMin = rayon + delta;

        for (double x = xmin; x <= xmax; x += pasSolution)
        {
            for (double y = ymin; y <= ymax; y += pasSolution)
            {
                const double distance = std::sqrt(x * x + y * y);

                // On conserve uniquement les points
                // suffisamment éloignés du cercle.
                if (distance >= distanceMin)
                {
                    pointsSolution.emplace_back(x, y);
                }
            }
        }

        const unsigned int nbPointsSolution = pointsSolution.size();

        std::cout << "\n=== Maillage de solution ===" << std::endl;
        std::cout << "  Domaine             : [" << xmin << ", " << xmax << "] x ["
                  << ymin << ", " << ymax << "]" << std::endl;
        std::cout << "  Taille du domaine   : " << L << " x " << L << std::endl;
        std::cout << "  Pas                 : " << pasSolution << std::endl;
        std::cout << "  Distance à la frontière : " << delta << std::endl;
        std::cout << "  Rayon minimal autorisé   : " << distanceMin << std::endl;
        std::cout << "  Nombre de points    : " << nbPointsSolution << std::endl;

        // -- Calcul du vecteur de p sur les milieux des bords --
        Vecteur vect_p(nbSegmentsMaillage);
        for (unsigned int i = 0; i < nbSegmentsMaillage; i++)
        {
            vect_p[i] = p_analytique(maillage[i].milieu, idxTroncature);
        }

        // -- Calcul du vecteur de p aux noeuds
        Vecteur pNoeuds(nbSegmentsMaillage);
        for (unsigned int i = 0; i < nbSegmentsMaillage; ++i)
        {
            pNoeuds[i] = p_analytique(maillage[i].P1, idxTroncature);
        }

        // -- Création du fichier de résultat --
        const string filename = string("outputs/u") + "_k" + std::to_string(k) + "_R" + std::to_string(rayon) + "_L" + std::to_string(L) + "_hM" + std::to_string(pasMaillage) + "_hS" + std::to_string(pasSolution) + "_d" + std::to_string(delta) + "_N" + std::to_string(idxTroncature) + "_q" + std::to_string(ordre) + ".txt";
        ofstream file(filename);

        if (!file.is_open())
        {
            std::cout << "ERROR: Le fichier " << filename << " n'a pas pu être ouvert" << endl;
            exit(-1);
        }

        // -- Récupération des coefficients de Legendre --
        const LegendreData &legendreData = get_legendre_data(ordre);

        // -- Construction solution approchée --
        double erreurCacheL2 = 0.0;
        double erreurExacteL2 = 0.0;
        double erreurTotaleL2 = 0.0;

        double erreurCacheMax = 0.0;
        double erreurExacteMax = 0.0;
        double erreurTotaleMax = 0.0;

        double tempsGreen = 0.0;
        double tempsCached = 0.0;

        double erreurCacheL2Lin = 0.0;
        double erreurExacteL2Lin = 0.0;
        double erreurTotaleL2Lin = 0.0;

        double erreurCacheMaxLin = 0.0;
        double erreurExacteMaxLin = 0.0;
        double erreurTotaleMaxLin = 0.0;

        double tempsGreenLin = 0.0;
        double tempsCachedLin = 0.0;

        for (const Point &Pj : pointsSolution)
        {
            complex<double> uGreen = 0.0;
            complex<double> uCached = 0.0;

            complex<double> uGreenLin = 0.0;
            complex<double> uCachedLin = 0.0;

            // -- Solution analytique --

            const complex<double> uExact = u_N_plus_analytique(Pj, rayon, idxTroncature);

            // -- Interpolation cte --

            // Green exacte
            std::clock_t startGreen = std::clock();

            for (unsigned int i = 0; i < nbSegmentsMaillage; i++)
            {
                uGreen += vect_p[i] * integ_simple(
                                          [&Pj](const Point &Q, double)
                                          {
                                              return green(Pj, Q);
                                          },
                                          maillage[i],
                                          legendreData);
            }

            std::clock_t endGreen = std::clock();

            tempsGreen += static_cast<double>(endGreen - startGreen) / CLOCKS_PER_SEC;

            // Green cached
            std::clock_t startCached = std::clock();

            for (unsigned int i = 0; i < nbSegmentsMaillage; i++)
            {
                uCached += vect_p[i] * integ_simple(
                                           [&Pj](const Point &Q, double)
                                           {
                                               return green_cached_vec(Pj, Q);
                                           },
                                           maillage[i],
                                           legendreData);
            }

            std::clock_t endCached = std::clock();

            tempsCached += static_cast<double>(endCached - startCached) / CLOCKS_PER_SEC;

            // Erreurs locales
            const double erreurCache = erreur_relative(uCached, uGreen);
            const double erreurExacte = erreur_relative(uGreen, uExact);
            const double erreurTotale = erreur_relative(uCached, uExact);

            // Accumulation
            erreurCacheL2 += erreurCache * erreurCache;
            erreurExacteL2 += erreurExacte * erreurExacte;
            erreurTotaleL2 += erreurTotale * erreurTotale;

            erreurCacheMax = std::max(erreurCacheMax, erreurCache);
            erreurExacteMax = std::max(erreurExacteMax, erreurExacte);
            erreurTotaleMax = std::max(erreurTotaleMax, erreurTotale);

            // -- Interpolation linéaire --

            // Green exacte Lin
            std::clock_t startGreenLin = std::clock();

            for (unsigned int i = 0; i < nbSegmentsMaillage; i++)
            {
                const std::complex<double> pA = pNoeuds[i];
                const std::complex<double> pB = pNoeuds[(i + 1) % nbSegmentsMaillage];
                uGreenLin += integ_simple([&Pj, pA, pB](const Point &Q, double r)
                                          {
                                            const double N1 = 0.5 * (1.0 - r);
                                            const double N2 = 0.5 * (1.0 + r);

                                            const std::complex<double> p = N1 * pA + N2 * pB;

                                            return green(Pj, Q) * p; },
                                          maillage[i],
                                          legendreData);
            }

            std::clock_t endGreenLin = std::clock();

            tempsGreenLin += static_cast<double>(endGreenLin - startGreenLin) / CLOCKS_PER_SEC;

            // Green cached Lin
            std::clock_t startCachedLin = std::clock();

            for (unsigned int i = 0; i < nbSegmentsMaillage; i++)
            {
                const std::complex<double> pA = pNoeuds[i];
                const std::complex<double> pB = pNoeuds[(i + 1) % nbSegmentsMaillage];
                uCachedLin += integ_simple([&Pj, pA, pB](const Point &Q, double r)
                                           {
                                            const double N1 = 0.5 * (1.0 - r);
                                            const double N2 = 0.5 * (1.0 + r);

                                            const std::complex<double> p = N1 * pA + N2 * pB;

                                            return green_cached_vec(Pj, Q) * p; },
                                           maillage[i],
                                           legendreData);
            }

            std::clock_t endCachedLin = std::clock();

            tempsCachedLin += static_cast<double>(endCachedLin - startCachedLin) / CLOCKS_PER_SEC;

            // Erreurs locales
            const double erreurCacheLin = erreur_relative(uCachedLin, uGreenLin);
            const double erreurExacteLin = erreur_relative(uGreenLin, uExact);
            const double erreurTotaleLin = erreur_relative(uCachedLin, uExact);

            // Accumulation
            erreurCacheL2Lin += erreurCacheLin * erreurCacheLin;
            erreurExacteL2Lin += erreurExacteLin * erreurExacteLin;
            erreurTotaleL2Lin += erreurTotaleLin * erreurTotaleLin;

            erreurCacheMaxLin = std::max(erreurCacheMaxLin, erreurCacheLin);
            erreurExacteMaxLin = std::max(erreurExacteMaxLin, erreurExacteLin);
            erreurTotaleMaxLin = std::max(erreurTotaleMaxLin, erreurTotaleLin);

            // -- Ecriture Fichier --

            file
                << Pj.x << " "
                << Pj.y << " "

                << uExact.real() << " "
                << uExact.imag() << " "

                << uGreen.real() << " "
                << uGreen.imag() << " "

                << uCached.real() << " "
                << uCached.imag() << " "

                << uGreenLin.real() << " "
                << uGreenLin.imag() << " "

                << uCachedLin.real() << " "
                << uCachedLin.imag() << " "

                << erreurCacheLin << " "
                << erreurExacteLin << " "
                << erreurTotaleLin << " "

                << erreurCache << " "
                << erreurExacte << " "
                << erreurTotale

                << "\n";
        }

        file.close();

        // -- Erreurs L2 --

        erreurCacheL2 = std::sqrt(erreurCacheL2 / static_cast<double>(nbPointsSolution));
        erreurExacteL2 = std::sqrt(erreurExacteL2 / static_cast<double>(nbPointsSolution));
        erreurTotaleL2 = std::sqrt(erreurTotaleL2 / static_cast<double>(nbPointsSolution));

        erreurCacheL2Lin = std::sqrt(erreurCacheL2Lin / static_cast<double>(nbPointsSolution));
        erreurExacteL2Lin = std::sqrt(erreurExacteL2Lin / static_cast<double>(nbPointsSolution));
        erreurTotaleL2Lin = std::sqrt(erreurTotaleL2Lin / static_cast<double>(nbPointsSolution));

        // -- Prints --

        std::cout << "\n=== Erreurs ===" << std::endl;

        std::cout << "    P0 erreur L² : Green = "
                  << erreurExacteL2
                  << ", Cache = "
                  << erreurCacheL2
                  << ", Totale = "
                  << erreurTotaleL2
                  << endl;

        std::cout << "    P0 erreur L∞ : Green = "
                  << erreurExacteMax
                  << ", Cache = "
                  << erreurCacheMax
                  << ", Totale = "
                  << erreurTotaleMax
                  << endl;

        std::cout << "    P1 erreur L² : Green = "
                  << erreurExacteL2Lin
                  << ", Cache = "
                  << erreurCacheL2Lin
                  << ", Totale = "
                  << erreurTotaleL2Lin
                  << endl;

        std::cout << "    P1 erreur L∞ : Green = "
                  << erreurExacteMaxLin
                  << ", Cache = "
                  << erreurCacheMaxLin
                  << ", Totale = "
                  << erreurTotaleMaxLin
                  << endl;

        // -- End time --

        std::clock_t end = std::clock();
        double seconds = static_cast<double>(end - start) / CLOCKS_PER_SEC;

        std::cout << "\n=== Temps ===" << std::endl;
        std::cout << "    execution           : " << seconds << " seconds" << endl;
        std::cout << "    green               : " << tempsGreen << " seconds" << endl;
        std::cout << "    green cached        : " << tempsCached << " seconds" << endl;
        std::cout << "    speedup             : " << tempsGreen / tempsCached << " x" << std::endl;
        std::cout << "    gain de temps       : " << (1.0 - tempsCached / tempsGreen) * 100.0 << " %" << std::endl;
        std::cout << "    green lin           : " << tempsGreenLin << " seconds" << endl;
        std::cout << "    green cached lin    : " << tempsCachedLin << " seconds" << endl;
        std::cout << "    speedup lin         : " << tempsGreenLin / tempsCachedLin << " x" << std::endl;
        std::cout << "    gain de temps lin   : " << (1.0 - tempsCachedLin / tempsGreenLin) * 100.0 << " %" << std::endl;

#endif

#ifdef TP2
        /* ----- TP2 ----- */

        // -- Parametres du problème --

        k = k;

        // -- Parametres Maillage --

        const double rayon = 1.; // rayon cercle du maillage
        pasMaillage = 0.02;      // pas du maillage

        // -- Parametres Solution --

        const double L = 4.;              // Domaine LxL pour le calcul de la solution
        const double pasSolution = 0.1;   // pas de la solution
        const double delta = pasSolution; // distance minimale entre les points de solution et le cercle

        // -- Parametre  d'approximation --

        const unsigned int idxTroncature = 25;
        const unsigned int ordre = 4;
        green_cache_step = 0.1 * min(pasMaillage, pasSolution) / k;
        max_index = static_cast<unsigned int>(std::ceil(L * std::sqrt(2.0) / green_cache_step));
        double tolGradConj = 1e-8;
        unsigned int maxIterGradConj = 1000;

        // -- Affichage Parametres --

        std::cout << "\n=== Parametres du probleme ===" << std::endl;
        std::cout << "  Nombre d'onde k          : " << k << std::endl;

        std::cout << "\n  --- Approximation ---" << std::endl;
        std::cout << "  Indice de troncature N   : " << idxTroncature << std::endl;
        std::cout << "  Ordre de quadrature      : " << ordre << std::endl;
        std::cout << "  Pas du cache de Green    : " << green_cache_step << std::endl;
        std::cout << "  Nombre max d'indices     : " << max_index << std::endl;

        // -- Start time --

        std::clock_t start = std::clock();

        // -- Création maillage --

        Point O(0, 0);
        Cercle cercle(rayon, O);
        Maillage maillage;
        maillage.ajoute_cercle(pasMaillage, cercle);
        const unsigned int nbSegmentsMaillage = maillage.size();

        std::cout << "\n=== Maillage du cercle ===" << std::endl;
        std::cout << "  Centre              : " << O << std::endl;
        std::cout << "  Rayon               : " << rayon << std::endl;
        std::cout << "  Pas                 : " << pasMaillage << std::endl;
        std::cout << "  Nombre de segments  : " << nbSegmentsMaillage << std::endl;

        // -- Maillage pour la solution --

        vector<Point> pointsSolution;

        const double xmin = -L / 2.0;
        const double xmax = L / 2.0;
        const double ymin = -L / 2.0;
        const double ymax = L / 2.0;

        const double distanceMin = rayon + delta;

        for (double x = xmin; x <= xmax; x += pasSolution)
        {
            for (double y = ymin; y <= ymax; y += pasSolution)
            {
                const double distance =
                    std::sqrt(x * x + y * y);

                // On conserve uniquement les points
                // suffisamment éloignés du cercle.
                if (distance >= distanceMin)
                {
                    pointsSolution.emplace_back(x, y);
                }
            }
        }

        const unsigned int nbPointsSolution = pointsSolution.size();

        std::cout << "\n=== Maillage de solution ===" << std::endl;
        std::cout << "  Domaine             : [" << xmin << ", " << xmax << "] x ["
                  << ymin << ", " << ymax << "]" << std::endl;
        std::cout << "  Taille du domaine   : " << L << " x " << L << std::endl;
        std::cout << "  Pas                 : " << pasSolution << std::endl;
        std::cout << "  Distance à la frontière : " << delta << std::endl;
        std::cout << "  Rayon minimal autorisé   : " << distanceMin << std::endl;
        std::cout << "  Nombre de points    : " << nbPointsSolution << std::endl;

        // -- Récupération des coefficients de Legendre --

        const LegendreData &legendreData = get_legendre_data(ordre);

        // -- Calcul du vecteur de b de la FV de l'équation intégrale --

        Vecteur vect_b(nbSegmentsMaillage);
        for (unsigned int i = 0; i < nbSegmentsMaillage; i++)
        {
            vect_b[i] = integ_simple([](const Point &Q, double)
                                     { return -u_inc(Q); },
                                     maillage[i],
                                     legendreData);
        }

        // ================================================================
        //   Version exacte A
        // ================================================================

        // -- Calcul de la matrice A de la FV de l'équation intégrale --

        std::cout << "\n=== Construction matrice exacte ===" << std::endl;

        std::clock_t startConstructionGreen = std::clock();

        MatriceSym A(nbSegmentsMaillage);

        for (unsigned int i = 0; i < nbSegmentsMaillage; i++)
        {
            for (unsigned int j = 0; j <= i; j++)
            {
                if (i != j) // cas général
                {
                    A(i, j) = integ_double([](const Point &Q, const Point &P, double)
                                           { return green(Q, P); },
                                           maillage[i],
                                           maillage[j],
                                           legendreData);
                }
                else
                {
                    const Segment S = maillage[i];
                    A(i, j) = integ_simple([&S](const Point &Q, double)
                                           { return integ_simple_log(Q, S); },
                                           S,
                                           legendreData);

                    A(i, j) += integ_double([](const Point &Q, const Point &P, double)
                                            { return green_reguliere(Q, P); },
                                            maillage[i],
                                            maillage[j],
                                            legendreData);
                }
            }
        }

        std::clock_t endConstructionGreen = std::clock();

        const double tempsConstructionGreen = static_cast<double>(endConstructionGreen - startConstructionGreen) / CLOCKS_PER_SEC;

        // ================================================================
        //   Version cached A
        // ================================================================

        // -- Calcul de la matrice A de la FV de l'équation intégrale --

        std::cout << "\n=== Construction matrice cached ===" << std::endl;

        std::clock_t startConstructionCached = std::clock();

        MatriceSym A_cached(nbSegmentsMaillage);

        for (unsigned int i = 0; i < nbSegmentsMaillage; i++)
        {
            for (unsigned int j = 0; j <= i; j++)
            {
                if (i != j) // cas général
                {
                    A_cached(i, j) = integ_double([](const Point &Q, const Point &P, double)
                                                  { return green_cached_vec(Q, P); },
                                                  maillage[i],
                                                  maillage[j],
                                                  legendreData);
                }
                else
                {
                    const Segment S = maillage[i];
                    A_cached(i, j) = integ_simple([&S](const Point &Q, double)
                                                  { return integ_simple_log(Q, S); },
                                                  S,
                                                  legendreData);

                    A_cached(i, j) += integ_double([](const Point &Q, const Point &P, double)
                                                   { return green_reguliere_cached_vec(Q, P); },
                                                   maillage[i],
                                                   maillage[j],
                                                   legendreData);
                }
            }
        }

        std::clock_t endConstructionCached = std::clock();

        const double tempsConstructionCached = static_cast<double>(endConstructionCached - startConstructionCached) / CLOCKS_PER_SEC;

        // -- Inversion de la matrice (méthode itérative) --

        std::cout << "\n=== Resolution systeme exact ===" << std::endl;

        std::clock_t startResolutionGreen = std::clock();

        Vecteur vect_p = gradConjMatSym(A, vect_b, tolGradConj, maxIterGradConj);

        std::clock_t endResolutionGreen = std::clock();

        const double tempsResolutionGreen = static_cast<double>(endResolutionGreen - startResolutionGreen) / CLOCKS_PER_SEC;

        std::cout << "\n=== Resolution systeme cached ===" << std::endl;

        std::clock_t startResolutionCached = std::clock();

        Vecteur vect_p_cached = gradConjMatSym(A_cached, vect_b, tolGradConj, maxIterGradConj);

        std::clock_t endResolutionCached = std::clock();

        const double tempsResolutionCached = static_cast<double>(endResolutionCached - startResolutionCached) / CLOCKS_PER_SEC;

        // -- Calcul du vecteur de p sur les milieux des bords --
        Vecteur vect_p_ana(nbSegmentsMaillage);
        for (unsigned int i = 0; i < nbSegmentsMaillage; i++)
        {
            vect_p_ana[i] = p_analytique(maillage[i].milieu, idxTroncature);
        }

        // -- Erreur du cache --

        double erreurPCache = (vect_p_cached - vect_p).norm() / vect_p.norm();

        // -- Comparaison p et p_ana --

        double erreurPExacte = (vect_p - vect_p_ana).norm() / vect_p_ana.norm();
        double erreurPCached = (vect_p_cached - vect_p_ana).norm() / vect_p_ana.norm();

        // -- Affichage erreurs p --

        std::cout << "\n=== Analyse de p ===" << std::endl;

        std::cout << "  Erreur p exacte / p analytique L2 : "
                  << erreurPExacte
                  << std::endl;

        std::cout << "  Erreur p cached / p analytique L2 : "
                  << erreurPCached
                  << std::endl;

        std::cout << "  Difference p / p_cached L2 : "
                  << erreurPCache
                  << std::endl;

        // -- Affichage temps systeme Ap = b --

        std::cout << "\n=== Temps de construction ===" << std::endl;

        std::cout << "  Matrice Green exacte : "
                  << tempsConstructionGreen
                  << " s"
                  << std::endl;

        std::cout << "  Matrice Green cached : "
                  << tempsConstructionCached
                  << " s"
                  << std::endl;

        std::cout << "  Speedup construction : "
                  << tempsConstructionGreen / tempsConstructionCached
                  << " x"
                  << std::endl;

        std::cout << "  Gain construction : "
                  << (1.0 -
                      tempsConstructionCached /
                          tempsConstructionGreen) *
                         100.0
                  << " %"
                  << std::endl;

        std::cout << "\n=== Temps de resolution ===" << std::endl;

        std::cout << "  Resolution exacte : "
                  << tempsResolutionGreen
                  << " s"
                  << std::endl;

        std::cout << "  Resolution cached : "
                  << tempsResolutionCached
                  << " s"
                  << std::endl;

        // -- Création du fichier de résultat --
        const string filename = string("outputs/u") + "_k" + std::to_string(k) + "_R" + std::to_string(rayon) + "_L" + std::to_string(L) + "_hM" + std::to_string(pasMaillage) + "_hS" + std::to_string(pasSolution) + "_d" + std::to_string(delta) + "_N" + std::to_string(idxTroncature) + "_q" + std::to_string(ordre) + ".txt";
        ofstream file(filename);

        if (!file.is_open())
        {
            std::cout << "ERROR: Le fichier " << filename << " n'a pas pu être ouvert" << endl;
            exit(-1);
        }

        // -- Construction solution approchée --
        double erreurCacheL2 = 0.0;
        double erreurExacteL2 = 0.0;
        double erreurTotaleL2 = 0.0;

        double erreurCacheMax = 0.0;
        double erreurExacteMax = 0.0;
        double erreurTotaleMax = 0.0;

        double tempsGreen = 0.0;
        double tempsCached = 0.0;

        for (const Point &Pj : pointsSolution)
        {
            complex<double> uGreen = 0.0;
            complex<double> uCached = 0.0;

            // -- Solution analytique --

            const complex<double> uExact = u_N_plus_analytique(Pj, rayon, idxTroncature);

            // -- Interpolation cte --

            // Green exacte
            std::clock_t startGreen = std::clock();

            for (unsigned int i = 0; i < nbSegmentsMaillage; i++)
            {
                uGreen += vect_p[i] * integ_simple(
                                          [&Pj](const Point &Q, double)
                                          {
                                              return green(Pj, Q);
                                          },
                                          maillage[i],
                                          legendreData);
            }

            std::clock_t endGreen = std::clock();

            tempsGreen += static_cast<double>(endGreen - startGreen) / CLOCKS_PER_SEC;

            // Green cached
            std::clock_t startCached = std::clock();

            for (unsigned int i = 0; i < nbSegmentsMaillage; i++)
            {
                uCached += vect_p_cached[i] * integ_simple(
                                                  [&Pj](const Point &Q, double)
                                                  {
                                                      return green_cached_vec(Pj, Q);
                                                  },
                                                  maillage[i],
                                                  legendreData);
            }

            std::clock_t endCached = std::clock();

            tempsCached += static_cast<double>(endCached - startCached) / CLOCKS_PER_SEC;

            // Erreurs locales
            const double erreurCache = erreur_relative(uCached, uGreen);
            const double erreurExacte = erreur_relative(uGreen, uExact);
            const double erreurTotale = erreur_relative(uCached, uExact);

            // Accumulation
            erreurCacheL2 += erreurCache * erreurCache;
            erreurExacteL2 += erreurExacte * erreurExacte;
            erreurTotaleL2 += erreurTotale * erreurTotale;

            erreurCacheMax = std::max(erreurCacheMax, erreurCache);
            erreurExacteMax = std::max(erreurExacteMax, erreurExacte);
            erreurTotaleMax = std::max(erreurTotaleMax, erreurTotale);

            // -- Ecriture Fichier --

            file
                << Pj.x << " "
                << Pj.y << " "

                << uExact.real() << " "
                << uExact.imag() << " "

                << uGreen.real() << " "
                << uGreen.imag() << " "

                << uCached.real() << " "
                << uCached.imag() << " "

                << erreurCache << " "
                << erreurExacte << " "
                << erreurTotale

                << "\n";
        }

        file.close();

        // -- Erreurs L2 --

        erreurCacheL2 = std::sqrt(erreurCacheL2 / static_cast<double>(nbPointsSolution));
        erreurExacteL2 = std::sqrt(erreurExacteL2 / static_cast<double>(nbPointsSolution));
        erreurTotaleL2 = std::sqrt(erreurTotaleL2 / static_cast<double>(nbPointsSolution));

        // -- Prints --

        std::cout << "\n=== Erreurs ===" << std::endl;

        std::cout << "    P0 erreur L² : Green = "
                  << erreurExacteL2
                  << ", Cache = "
                  << erreurCacheL2
                  << ", Totale = "
                  << erreurTotaleL2
                  << endl;

        std::cout << "    P0 erreur L∞ : Green = "
                  << erreurExacteMax
                  << ", Cache = "
                  << erreurCacheMax
                  << ", Totale = "
                  << erreurTotaleMax
                  << endl;

        // -- End time --

        std::clock_t end = std::clock();
        double seconds = static_cast<double>(end - start) / CLOCKS_PER_SEC;

        std::cout << "\n=== Temps ===" << std::endl;
        std::cout << "    execution           : " << seconds << " seconds" << endl;
        std::cout << "    green               : " << tempsGreen << " seconds" << endl;
        std::cout << "    green cached        : " << tempsCached << " seconds" << endl;
        std::cout << "    speedup             : " << tempsGreen / tempsCached << " x" << std::endl;
        std::cout << "    gain de temps       : " << (1.0 - tempsCached / tempsGreen) * 100.0 << " %" << std::endl;

#endif

#ifdef TP2_etude_CG
        /* ----- TP2 (étude CG) : temps du gradient conjugué en fonction de N et nq -----
         *
         * Pour chaque couple (nombre de segments N, nombre de points de quadrature nq) :
         *   - construction de A (version "cached"), chronométrée pour information,
         *   - résolution A p = b par gradConjMatSym, chronométrée (minimum sur
         *     plusieurs répétitions pour limiter le bruit de mesure),
         *   - erreur relative sur p (pour vérifier que la résolution a convergé).
         *
         * Résultats dans outputs/etude_temps_CG.txt
         *
         * Remarque : gradConjMatSym ne renvoie pas son nombre d'itérations. Le temps
         * mesuré combine donc "coût d'une itération" (~N^2) et "nombre d'itérations"
         * (qui dépend du conditionnement de A, donc de N et nq).
         */
        {
            // -- Parametres du probleme (k vient de la config, comme dans TP2) --

            const double rayon = 1.;
            const double L = 4.; // sert uniquement a dimensionner le cache de Green
            const unsigned int idxTroncature = 25;
            const double tolGradConj = 1e-8;
            const unsigned int maxIterGradConj = 1000;

            // -- Parametres de l'etude --

            const vector<unsigned int> nbSegmentsValues = {10, 20, 50, 100, 200, 500,
                                                           1000, 2000, 5000, 10000};

            const unsigned int nqMin = 4;
            const unsigned int nqMax = 20;
            const unsigned int nqPas = 2; // 1 pour tous les entiers (2x plus long)

            // Le gradient conjugue est rejoue plusieurs fois et on garde le minimum.
            const unsigned int nbRepetitions = 3;

            // -- Cache de Green (regle une fois, sur le maillage le plus fin) --

            const unsigned int nbSegmentsMax =
                *std::max_element(nbSegmentsValues.begin(), nbSegmentsValues.end());
            const double pasMin = 2.0 * pi * rayon / nbSegmentsMax;
            green_cache_step = 0.1 * pasMin / k;
            max_index = static_cast<unsigned int>(std::ceil(L * std::sqrt(2.0) / green_cache_step));

            std::cout << "\n=== Etude temps du gradient conjugue ===" << std::endl;
            std::cout << "  Nombre d'onde k   : " << k << std::endl;
            std::cout << "  Segments          : " << nbSegmentsValues.front()
                      << " -> " << nbSegmentsMax << " (" << nbSegmentsValues.size() << " valeurs)" << std::endl;
            std::cout << "  Quadrature        : " << nqMin << " -> " << nqMax
                      << " (pas " << nqPas << ")" << std::endl;
            std::cout << "  Repetitions du CG : " << nbRepetitions << " (minimum retenu)" << std::endl;
            std::cout << "  Tolerance / iter. : " << tolGradConj << " / " << maxIterGradConj << std::endl;

            // -- Fichier de sortie --

            const string filename = "outputs/etude_temps_CG.txt";
            ofstream file(filename);

            if (!file.is_open())
            {
                std::cout << "ERROR: Le fichier " << filename << " n'a pas pu être ouvert" << endl;
                exit(-1);
            }

            file << "# Temps du gradient conjugue (gradConjMatSym) en fonction de N et nq" << endl;
            file << "# k=" << k << " rayon=" << rayon << " idxTroncature=" << idxTroncature
                 << " tolGradConj=" << tolGradConj << " maxIterGradConj=" << maxIterGradConj
                 << " nbRepetitions=" << nbRepetitions << endl;
            file << "# Temps en secondes (horloge murale). temps_CG = minimum sur les repetitions." << endl;
            file << "# nbSegments_demande nbSegments_reel pasMaillage nq "
                 << "temps_assemblage temps_CG erreur_p" << endl;

            const Point O(0, 0);
            const Cercle cercle(rayon, O);

            // -- Construction d'un maillage du cercle avec n segments (cf. TP2_etude) --

            auto construireMaillage = [&](unsigned int n, double &pasUtilise)
            {
                const double pasIdeal = 2.0 * pi * rayon / n;
                const double facteurs[] = {1.0, 1.0 - 1e-9, 1.0 + 1e-9, 1.0 - 1e-6, 1.0 + 1e-6};

                Maillage meilleur;
                int meilleurEcart = -1;

                for (double f : facteurs)
                {
                    Maillage m;
                    m.ajoute_cercle(pasIdeal * f, cercle);
                    const int ecart = std::abs(static_cast<int>(m.size()) - static_cast<int>(n));

                    if (meilleurEcart < 0 || ecart < meilleurEcart)
                    {
                        meilleur = m;
                        meilleurEcart = ecart;
                        pasUtilise = pasIdeal * f;
                    }
                    if (ecart == 0)
                        break;
                }
                return meilleur;
            };

            using Horloge = std::chrono::steady_clock;
            double puits = 0.0; // empeche le compilateur d'eliminer les resolutions repetees

            // -- Boucle principale --

            for (unsigned int nDemande : nbSegmentsValues)
            {
                double pas = 0.0;
                const Maillage maillage = construireMaillage(nDemande, pas);
                const unsigned int n = maillage.size();

                if (n != nDemande)
                    std::cout << "  [Attention] " << nDemande << " segments demandes, "
                              << n << " obtenus (pas = " << pas << ")" << std::endl;

                Vecteur vect_p_ana(n);
                for (unsigned int i = 0; i < n; i++)
                    vect_p_ana[i] = p_analytique(maillage[i].milieu, idxTroncature);

                for (unsigned int nq = nqMin; nq <= nqMax; nq += nqPas)
                {
                    const LegendreData &legendreData = get_legendre_data(nq);

                    // Second membre
                    Vecteur vect_b(n);
                    for (unsigned int i = 0; i < n; i++)
                    {
                        vect_b[i] = integ_simple([](const Point &Q, double)
                                                 { return -u_inc(Q); },
                                                 maillage[i],
                                                 legendreData);
                    }

                    // Assemblage de A (cached), chronometre pour information
                    const auto ta0 = Horloge::now();

                    MatriceSym A_cached(n);
                    for (unsigned int i = 0; i < n; i++)
                    {
                        for (unsigned int j = 0; j <= i; j++)
                        {
                            if (i != j)
                            {
                                A_cached(i, j) = integ_double([](const Point &Q, const Point &P, double)
                                                              { return green_cached_vec(Q, P); },
                                                              maillage[i],
                                                              maillage[j],
                                                              legendreData);
                            }
                            else
                            {
                                const Segment S = maillage[i];
                                A_cached(i, j) = integ_simple([&S](const Point &Q, double)
                                                              { return integ_simple_log(Q, S); },
                                                              S,
                                                              legendreData);

                                A_cached(i, j) += integ_double([](const Point &Q, const Point &P, double)
                                                               { return green_reguliere_cached_vec(Q, P); },
                                                               maillage[i],
                                                               maillage[j],
                                                               legendreData);
                            }
                        }
                    }

                    const auto ta1 = Horloge::now();
                    const double tempsAssemblage = std::chrono::duration<double>(ta1 - ta0).count();

                    // Gradient conjugue : 1ere resolution gardee pour l'erreur,
                    // les suivantes ne servent qu'a affiner la mesure du temps.
                    const auto tc0 = Horloge::now();
                    const Vecteur vect_p = gradConjMatSym(A_cached, vect_b, tolGradConj, maxIterGradConj);
                    const auto tc1 = Horloge::now();
                    double tempsCG = std::chrono::duration<double>(tc1 - tc0).count();

                    for (unsigned int r = 1; r < nbRepetitions; r++)
                    {
                        const auto t0 = Horloge::now();
                        const Vecteur vect_rep = gradConjMatSym(A_cached, vect_b, tolGradConj, maxIterGradConj);
                        const auto t1 = Horloge::now();

                        tempsCG = std::min(tempsCG, std::chrono::duration<double>(t1 - t0).count());
                        puits += vect_rep.norm();
                    }

                    const double erreur = (vect_p - vect_p_ana).norm() / vect_p_ana.norm();

                    // Ecriture (flush a chaque ligne : run long)
                    file << nDemande << " " << n << " " << pas << " " << nq << " "
                         << tempsAssemblage << " " << tempsCG << " " << erreur << endl;

                    std::cout << "  N=" << n << " nq=" << nq
                              << " | t_assemblage=" << tempsAssemblage << " s"
                              << " | t_CG=" << tempsCG << " s"
                              << " | err_p=" << erreur << std::endl;
                }
            }

            file.close();
            std::cout << "Fichier " << filename << " exporte. (controle: " << puits << ")" << std::endl;
        }
#endif

#ifdef TP2_bis
        /* ----- TP2 : matrice free ----- */

        // -- Parametres du problème --

        k = k;

        // -- Parametres Maillage --

        const double rayon = 1.; // rayon cercle du maillage
        pasMaillage = 0.001;     // pas du maillage

        // -- Parametres Solution --

        const double L = 4.;              // Domaine LxL pour le calcul de la solution
        const double pasSolution = 0.1;   // pas de la solution
        const double delta = pasSolution; // distance minimale entre les points de solution et le cercle

        // -- Parametre  d'approximation --

        const unsigned int idxTroncature = 25;
        const unsigned int ordre = 4;
        green_cache_step = 0.1 * min(pasMaillage, pasSolution) / k;
        max_index = static_cast<unsigned int>(std::ceil(L * std::sqrt(2.0) / green_cache_step));
        double tolGradConj = 1e-8;
        unsigned int maxIterGradConj = 1000;

        // -- Affichage Parametres --

        std::cout << "\n=== Parametres du probleme ===" << std::endl;
        std::cout << "  Nombre d'onde k          : " << k << std::endl;

        std::cout << "\n  --- Approximation ---" << std::endl;
        std::cout << "  Indice de troncature N   : " << idxTroncature << std::endl;
        std::cout << "  Ordre de quadrature      : " << ordre << std::endl;
        std::cout << "  Pas du cache de Green    : " << green_cache_step << std::endl;
        std::cout << "  Nombre max d'indices     : " << max_index << std::endl;

        // -- Start time --

        std::clock_t start = std::clock();

        // -- Création maillage --

        Point O(0, 0);
        Cercle cercle(rayon, O);
        Maillage maillage;
        maillage.ajoute_cercle(pasMaillage, cercle);
        const unsigned int nbSegmentsMaillage = maillage.size();

        std::cout << "\n=== Maillage du cercle ===" << std::endl;
        std::cout << "  Centre              : " << O << std::endl;
        std::cout << "  Rayon               : " << rayon << std::endl;
        std::cout << "  Pas                 : " << pasMaillage << std::endl;
        std::cout << "  Nombre de segments  : " << nbSegmentsMaillage << std::endl;

        // -- Maillage pour la solution --

        vector<Point> pointsSolution;

        const double xmin = -L / 2.0;
        const double xmax = L / 2.0;
        const double ymin = -L / 2.0;
        const double ymax = L / 2.0;

        const double distanceMin = rayon + delta;

        for (double x = xmin; x <= xmax; x += pasSolution)
        {
            for (double y = ymin; y <= ymax; y += pasSolution)
            {
                const double distance =
                    std::sqrt(x * x + y * y);

                // On conserve uniquement les points
                // suffisamment éloignés du cercle.
                if (distance >= distanceMin)
                {
                    pointsSolution.emplace_back(x, y);
                }
            }
        }

        const unsigned int nbPointsSolution = pointsSolution.size();

        std::cout << "\n=== Maillage de solution ===" << std::endl;
        std::cout << "  Domaine             : [" << xmin << ", " << xmax << "] x ["
                  << ymin << ", " << ymax << "]" << std::endl;
        std::cout << "  Taille du domaine   : " << L << " x " << L << std::endl;
        std::cout << "  Pas                 : " << pasSolution << std::endl;
        std::cout << "  Distance à la frontière : " << delta << std::endl;
        std::cout << "  Rayon minimal autorisé   : " << distanceMin << std::endl;
        std::cout << "  Nombre de points    : " << nbPointsSolution << std::endl;

        // -- Récupération des coefficients de Legendre --

        const LegendreData &legendreData = get_legendre_data(ordre);

        // -- Calcul du vecteur de b de la FV de l'équation intégrale --

        Vecteur vect_b(nbSegmentsMaillage);
        for (unsigned int i = 0; i < nbSegmentsMaillage; i++)
        {
            vect_b[i] = integ_simple([](const Point &Q, double)
                                     { return -u_inc(Q); },
                                     maillage[i],
                                     legendreData);
        }

        // -- Inversion de la matrice (méthode itérative) --

        std::cout << "\n=== Resolution systeme cached ===" << std::endl;

        std::clock_t startResolutionCached = std::clock();

        Vecteur vect_p_cached = gradConjMatrixFree(maillage, vect_b, legendreData, tolGradConj, maxIterGradConj);

        std::clock_t endResolutionCached = std::clock();

        const double tempsResolutionCached = static_cast<double>(endResolutionCached - startResolutionCached) / CLOCKS_PER_SEC;

        // -- Erreur avec le vecteur de p sur les milieux des bords --
        Vecteur vect_p_ana(nbSegmentsMaillage);
        for (unsigned int i = 0; i < nbSegmentsMaillage; i++)
        {
            vect_p_ana[i] = p_analytique(maillage[i].milieu, idxTroncature);
        }

        double erreurPCached = (vect_p_cached - vect_p_ana).norm() / vect_p_ana.norm();

        std::cout << "\n=== Analyse de p ===" << std::endl;

        std::cout << "  Erreur p cached / p analytique L2 : "
                  << erreurPCached
                  << std::endl;

        // -- Affichage temps systeme Ap = b --
        std::cout << "\n=== Temps de resolution ===" << std::endl;

        std::cout << "  Resolution cached : "
                  << tempsResolutionCached
                  << " s"
                  << std::endl;

        // -- Création du fichier de résultat --
        const string filename = string("outputs/u") + "_k" + std::to_string(k) + "_R" + std::to_string(rayon) + "_L" + std::to_string(L) + "_hM" + std::to_string(pasMaillage) + "_hS" + std::to_string(pasSolution) + "_d" + std::to_string(delta) + "_N" + std::to_string(idxTroncature) + "_q" + std::to_string(ordre) + ".txt";
        ofstream file(filename);

        if (!file.is_open())
        {
            std::cout << "ERROR: Le fichier " << filename << " n'a pas pu être ouvert" << endl;
            exit(-1);
        }

        // -- Construction solution approchée --

        double erreurTotaleL2 = 0.0;
        double erreurTotaleMax = 0.0;

        double tempsCached = 0.0;

        for (const Point &Pj : pointsSolution)
        {
            complex<double> uCached = 0.0;

            // -- Solution analytique --

            const complex<double> uExact = u_N_plus_analytique(Pj, rayon, idxTroncature);

            // -- Interpolation cte --

            // Green cached
            std::clock_t startCached = std::clock();

            for (unsigned int i = 0; i < nbSegmentsMaillage; i++)
            {
                uCached += vect_p_cached[i] * integ_simple(
                                                  [&Pj](const Point &Q, double)
                                                  {
                                                      return green_cached_vec(Pj, Q);
                                                  },
                                                  maillage[i],
                                                  legendreData);
            }

            std::clock_t endCached = std::clock();

            tempsCached += static_cast<double>(endCached - startCached) / CLOCKS_PER_SEC;

            // Erreurs locales
            const double erreurTotale = erreur_relative(uCached, uExact);

            // Accumulation
            erreurTotaleL2 += erreurTotale * erreurTotale;
            erreurTotaleMax = std::max(erreurTotaleMax, erreurTotale);

            // -- Ecriture Fichier --

            file
                << Pj.x << " "
                << Pj.y << " "

                << uExact.real() << " "
                << uExact.imag() << " "

                << uCached.real() << " "
                << uCached.imag() << " "

                << erreurTotale

                << "\n";
        }

        file.close();

        // -- Erreurs L2 --

        erreurTotaleL2 = std::sqrt(erreurTotaleL2 / static_cast<double>(nbPointsSolution));

        // -- Prints --

        std::cout << "\n=== Erreurs ===" << std::endl;

        std::cout << "    P0 erreur L² : Totale = "
                  << erreurTotaleL2
                  << endl;

        std::cout << "    P0 erreur L∞ : Totale = "
                  << erreurTotaleMax
                  << endl;

        // -- End time --

        std::clock_t end = std::clock();
        double seconds = static_cast<double>(end - start) / CLOCKS_PER_SEC;

        std::cout << "\n=== Temps ===" << std::endl;

        std::cout << "    execution        : "
                  << seconds << " s" << std::endl;

        std::cout << "    resolution GC    : "
                  << tempsResolutionCached << " s" << std::endl;

        std::cout << "    reconstruction u : "
                  << tempsCached << " s" << std::endl;

#endif
    }

    return 0;
}