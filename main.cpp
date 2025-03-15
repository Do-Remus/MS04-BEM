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
        Maillage monMaillage = genere_maillage_couche_diffusante(5, 10, 5, 0.01);
        monMaillage.export_maillage("outputs/test.txt");

        // tests Vecteur
        Vecteur u = {1, 2, 3};
        cout << u;

        // tests Integrales 
        Point P1(0, 0);
        Point P2(1, 0);
        Segment Seg(P1, P2);
        cout << "Integrale: " << integ_simple(Seg, f, 100) << endl;

        //test green reguliere
        for(int i=10; i>=0; i--){
            Point R(i/10.,0);
            cout<<"G -log("<<i/10.<<") ="<<green_reguliere(R,P1);
        }
        cout<<endl;

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
        Maillage maillageFinale = genere_maillage_couche_diffusante(1, h, e, pas);
        maillageFinale.export_maillage("outputs/maillage.txt");

        cout << "Génération des matrices..." << endl;
        MatriceSym A(maillageFinale.size(), 0);
        Vecteur P(maillageFinale.size(), 0);
        genere_coefficient_matrice_A(A, maillageFinale, pas*0.1);
        cout<<"A= "<<A<<endl;
        genere_coefficient_vecteur_P(P, maillageFinale, pas*0.1);
        cout<<"P = "<<P<<endl;

        cout << "Résolution du système..." << endl;
        Vecteur solution = resolution_systeme_lineaire(A, P);

        cout << "La solution est:" << solution << endl;
        for(int i=0; i<10;i++){
            for(int j=0; j<10; j++){
                Point IJ((i/10.)*e,(j/10.)*h);
                cout<< "p("<<(i/10.)*e<<","<<(j/10.)*h<<") ="<<p(maillageFinale,pas,solution, P,IJ)<<endl;
            }
        }

    }

    return 0;
}
