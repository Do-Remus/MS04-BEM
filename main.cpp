#include <iostream>
#include <time.h>
#include "src/config/config.hpp"
#include "src/config/constantes.hpp"
#include "src/config/external.hpp"
#include "src/headers/utils.hpp"
#include "src/headers/test_functions.hpp"

complex<double> f(const Point &A)
{
    return A.x * A.x;
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
        cout << "A = (" << M.x << "," << M.y << ")" << endl;
        cout << "A =" << M << endl;
        cout << "M+C = " << M + C;
        cout << "M-C = " << M - C;
        cout << "3*C = " << 3 * C;
        cout << "C/5 = " << C / 5;
        cout << "norme C = " << C.norm() << endl;
        cout << "produit vectoriel composante z = " << M * C << endl;

        // tests Cercle
        Cercle B(10, M);
        cout << "obstacle point " << B.centre.x << "," << B.centre.y << endl;
        cout << " obstacle rayon = " << B.rayon << endl;

        // tests Segment
        Segment S(M, C);
        cout << "segment P1 " << S.P1.x << "," << S.P1.y << endl;
        cout << "segment P2 " << S.P2.x << "," << S.P2.y << endl;

        // tests Maillage
        vector<Cercle> obstaclesTest;
        Point O(0, 0);
        double a = 1.;
        Cercle Cercle0(a, O);
        obstaclesTest.push_back(Cercle0);
        Maillage monMaillage;
        cout << " k= " << k << ", a =" << a << " N =" << N << endl;
        monMaillage.ajoute_cercle(pasMaillage, Cercle0);
        monMaillage.export_maillage("outputs/test.txt");
        cout << "exported maillage" << endl;
        // for (unsigned int i = 0; i < obstaclesTest.size(); i++)
        //{
        //     cout << "Cercle " << i << ": r = " << obstaclesTest[i].rayon << " centre = " << obstaclesTest[i].centre;
        // }

        /*
        // tests sommes partielles U_n^+ , q et p
        cout<<"tests sommes partielles"<<endl;
        vector<Point> pointsCercle = monMaillage.PointsMaillage();
        cout<<"out of PointsMaillage method"<<endl;
        std::string filename_qConverged = "outputs/qResu_k=" + std::to_string(k) + "_a=" + std::to_string(a) +  "_N=" + std::to_string(N) + ".txt";
        ofstream qC(filename_qConverged);
        std::string filename_pConverged = "outputs/pResu_k=" + std::to_string(k)  + "_a=" + std::to_string(a) + "_N=" + std::to_string(N) +  ".txt";
        ofstream pC(filename_pConverged);

        if (!qC.is_open())
        {
        cout << "ERROR: Le fichier " << filename_qConverged << " n'a pas pu être ouvert" << endl;
        exit(-1);
        }
        if (!pC.is_open())
        {
        cout << "ERROR: Le fichier " << filename_pConverged << " n'a pas pu être ouvert" << endl;
        exit(-1);
        }
        qC << "k ="<<k<<" a="<<a<<endl;
        pC << "k ="<<k<<" a="<<a<<endl;
        cout<<"number of points in mesh ="<<pointsCercle.size()<<endl;
        for (unsigned int j =0; j<10;j++){
            int i = j*pointsCercle.size()/10;
            Point P = pointsCercle[i];
            double theta = P.theta();
            std::string filename_q = "outputs/qResu_k=" + std::to_string(k) + "_a=" + std::to_string(a) + "_theta=" + std::to_string(theta) +  "_N=" + std::to_string(N) + ".txt";
            std::string filename_p = "outputs/pResu_k=" + std::to_string(k) + "_a=" + std::to_string(a) + "_theta=" + std::to_string(theta) + "_N=" + std::to_string(N) +  ".txt";
            std::string filename_solExt = "outputs/solExtResu_k=" + std::to_string(k) + "_a=" + std::to_string(a) + "_theta=" + std::to_string(theta) + "_N=" + std::to_string(N) +  ".txt";
            complex<double> q_approche = q(P, N,filename_q);
            complex<double> p_approche = p(P, N, filename_p);

            complex<double> sol_ext= u_N_plus(P, a, N, filename_solExt);
            qC << P << " "<< theta<<" "<< q_approche.real() << " "<< q_approche.imag()<<endl;
            pC << P << " "<< theta<<" "<< p_approche.real() << " "<< p_approche.imag()<<endl;
        }

        double err = 0;
        for (int t = 0; t < 10; t++) {
        double th = M_PI * t / 9.0;
            Point P(3.*a*cos(th), 3.*a*sin(th));
            std::string bin_theta = "outputs/solExtFarResu_k=" + std::to_string(k) + "_a=" + std::to_string(a) + "_theta=" + std::to_string(th) +  "_N=" + std::to_string(N)+ ".txt";
            complex<double> u = u_N_plus(P, a, N, bin_theta);
            err = max(err, abs(u + exp(-I * k * a * cos(th))));   // +: u = -uinc
        }
        cout << "N=" << N << " boundary error = " << err << endl;



    const string sol_externe = "outputs/u_N_plusResu_k=" + std::to_string(k) + "_a=" + std::to_string(a) +  "_N=" + std::to_string(N) + ".txt";
    export_obsctacles(cheminFichierObstacles,obstaclesTest );
    exporte_solution_analytique(sol_externe, a, nbPasExport, 10, N);
    double Hdiff;
    for(int i=1; i<=8; i++){
        Hdiff = pow(10., -i);
        export_fd_q_p("outputs/qFD_k=" + std::to_string(k) + "_a=" + std::to_string(a) + "_N=" + std::to_string(N) + "_h=" +std::to_string(Hdiff) + ".txt",
              "outputs/pFD_k=" + std::to_string(k) + "_a=" + std::to_string(a) + "_N=" + std::to_string(N) + "_h=" +std::to_string(Hdiff)+ ".txt",
              a, N, 50, Hdiff);

 */
    }

    // test des quadratures (calculs d'intégrales)
    fun_double func = [](const double &x) -> complex<double>
    {
        return 3 * x + x * x + 2 - x * x * x * x;
    };

    fun_double x2 = [](const double &x) -> complex<double>
    {
        return x * x;
    };

    complex<double> result = integ_simple_segsimple(func);
    cout << "value =" << result << "error =" << result - 14. / 3. + 2. / 5. << endl;

    complex<double> result1 = integ_simple_segsimple_n(func, 2);
    complex<double> result2 = integ_simple_segsimple_n(func, 4);
    cout << "for integ_n, n=2 value =" << result1 << "error =" << result1 - 14. / 3. + 2. / 5. << endl;
    cout << "for integ_n, n=4 value =" << result2 << "error =" << result2 - 14. / 3. + 2. / 5. << endl;

    complex<double> result3 = integ_simple_segment_ab_n(x2, 2, 0, 4);
    cout << "x2 btw 0 and 4 =" << result3 << endl;

    return 0;
}