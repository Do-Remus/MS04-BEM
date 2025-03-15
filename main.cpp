#include <iostream>
#include <time.h>
#include "src/config/config.hpp"
#include "src/config/constantes.hpp"
#include "src/config/external.hpp"
#include "src/headers/utils.hpp"

complex<double> f(const Point &A)
{
    return A.x * A.x;
}

int main()
{
    // initialisation du random
    srand(time(NULL));

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
        Cercle B;
        cout << "obstacle point " << B.centre.x << "," << B.centre.y << endl;
        cout << " obstacle rayon = " << B.rayon << endl;

        // tests Segment
        Segment S(M, C);
        cout << "segment P1 " << S.P1.x << "," << S.P1.y << endl;
        cout << "segment P2 " << S.P2.x << "," << S.P2.y << endl;

        // tests Maillage
        vector<Cercle> obstaclesTest;
        Maillage monMaillage = genere_maillage_couche_diffusante(5, 10, 5, 0.01, obstaclesTest);
        monMaillage.export_maillage("outputs/test.txt");
        for (unsigned int i = 0; i < obstaclesTest.size(); i++)
        {
            cout << "Cercle " << i << ": r = " << obstaclesTest[i].rayon << " centre = " << obstaclesTest[i].centre;
        }

        // tests Vecteur
        Vecteur u = {1, 2, 3};
        cout << u;

        // tests Integrales
        Point P1(0, 0);
        Point P2(1, 0);
        Segment Seg(P1, P2);
        cout << "Integrale: " << integ_simple(Seg, f, 100) << endl;

        // test green reguliere
        for (int i = 10; i >= 0; i--)
        {
            Point R(i / 10., 0);
            cout << "G -log(" << i / 10. << ") =" << green_reguliere(R, P1);
        }
        cout << endl;

        // test integrale pour p
        cout << "integrale pour p" << integrale_pour_p(S, hankel_derivate, 10, P1) << endl;

        // tests MatriceSym
        MatriceSym L(5);
        Vecteur D(5);
        MatriceSym AAA(5);
        cout << AAA << endl;
        AAA(1, 1) = 4;
        AAA(2, 2) = 9;
        AAA(0, 0) = 3;
        AAA(3, 3) = 5;
        AAA(4, 4) = 1;
        AAA(1, 3) = -2;
        AAA(2, 4) = -0.5;
        cout << AAA << endl;
        // décomposition LDL
        AAA.decomposition_LDL(L, D);
        cout << L << endl;
        cout << D << endl;
        // résolution sytème linéaire
        Vecteur PY(5);
        PY[0] = 4;
        PY[1] = 1;
        PY[2] = -8;
        PY[3] = 6;
        PY[4] = 0;
        cout << PY << endl;
        Vecteur Y = resolution_systeme_lineaire(AAA, PY);
        cout << PY << endl;
        cout << Y << endl;
    }

    if (effectuerLaSimulation)
    {
        // programme
        cout << "Création du maillage..." << endl;
        double h = 10.;
        double e = 2.;
        double pas = 0.1;
        vector<Cercle> obstables;
        Maillage maillage = genere_maillage_couche_diffusante(1, h, e, pas, obstables);
        maillage.export_maillage("outputs/maillage.txt");

        cout << "Génération des matrices..." << endl;
        MatriceSym A(maillage.size(), 0);
        Vecteur P(maillage.size(), 0);
        genere_coefficient_matrice_A(A, maillage, pas * 0.1);
        genere_coefficient_vecteur_P(P, maillage, pas * 0.1);

        cout << "Résolution du système..." << endl;
        Vecteur solution = resolution_systeme_lineaire(A, P);

        cout << "La solution est:" << solution << endl;
        exporte_solution("outputs/output.txt", maillage, solution, P, e, h, pas, 1000);
    }

    return 0;
}
