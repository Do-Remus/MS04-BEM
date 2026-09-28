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
        pasMaillage = 0.01;      // pas du maillage

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
        std::cout << "  Distance minimale   : " << delta << std::endl;
        std::cout << "  Rayon exclu         : " << distanceMin << std::endl;
        std::cout << "  Nombre de points    : " << nbPointsSolution << std::endl;

        // -- Calcul du vecteur de p sur les milieux des bords --
        vector<complex<double>> vect_p;
        for (unsigned int i = 0; i < nbSegmentsMaillage; i++)
        {
            vect_p.push_back(p_analytique(maillage[i].milieu, idxTroncature));
        }

        // -- Calcul du vecteur de p aux noeuds
        std::vector<std::complex<double>> pNoeuds;
        for (unsigned int i = 0; i < nbSegmentsMaillage; ++i)
        {
            pNoeuds.push_back(p_analytique(maillage[i].P1, idxTroncature));
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
            const double erreurCache = std::abs(uGreen - uCached);
            const double erreurExacte = std::abs(uGreen - uExact);
            const double erreurTotale = std::abs(uExact - uCached);

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
            const double erreurCacheLin = std::abs(uGreenLin - uCachedLin);
            const double erreurExacteLin = std::abs(uGreenLin - uExact);
            const double erreurTotaleLin = std::abs(uExact - uCachedLin);

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
        /* ----- TP1 ----- */

        // -- Parametres du problème --
        k = 10;

        // -- Parametres Maillage --
        const double rayon = 1.; // rayon cercle du maillage
        pasMaillage = 0.01;      // pas du maillage

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
        std::cout << "  Distance minimale   : " << delta << std::endl;
        std::cout << "  Rayon exclu         : " << distanceMin << std::endl;
        std::cout << "  Nombre de points    : " << nbPointsSolution << std::endl;

#endif
    }

    return 0;
}
