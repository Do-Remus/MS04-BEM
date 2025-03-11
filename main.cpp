#include <iostream>
#include <time.h>
#include "point.hpp"
#include "cercle.hpp"
#include "cercle.hpp"
#include "segment.hpp"
#include "utils.hpp"
#include "matrice.hpp"
#include "integrale.hpp"
using namespace std;

complex<double> f(const Point &A)
{
    return A.x * A.x;
}

int main()
{
    srand(time(NULL));
    Point M(1, 1);
    Point C(1, 2);
    cout << "A =(" << M.x << "," << M.y << ")" << endl;
    Cercle B;
    cout << "obstacle point " << B.centre.x << "," << B.centre.y << endl;
    cout << " obstacle rayon = " << B.rayon << endl;

    Segment S(M, C);
    cout << "segment P1" << S.P1.x << "," << S.P1.y << endl;
    cout << "segment P2" << S.P2.x << "," << S.P2.y << endl;
    // test ostream function
    cout << "A =" << M << endl;
    // test operations
    cout << "M+C" << M + C;
    cout << "M-C" << M - C;
    cout << "3*C" << 3 * C;
    cout << "C/5" << C / 5;
    cout << "norme C" << C.norm() << endl;
    cout << "produit vectoriel composante z =" << M * C;

    Maillage monMaillage = genere_maillage_couche_diffusante(5, 10, 5, 0.01);

    monMaillage.export_maillage("output.txt");

    // test Vecteur
    Vecteur u = {1, 2, 3};
    cout << u;

    // test Integrales

    Point P1(0, 0);
    Point P2(1, 0);
    Segment Seg(P1, P2);
    cout << "Integrale: " << integ_simple(Seg, f, 100) << endl;


    // test decomposition LDL
    MatriceSym L(5);
    Vecteur D(5);
    // test MatriceSym
    MatriceSym AAA(5);
    cout<<AAA<<endl;
    AAA(1,1)=4;
    AAA(2,2)=9;
    AAA(0,0)=3;
    AAA(3,3)=5;
    AAA(4,4)=1;
    AAA(1,3)=-2;
    AAA(2,4)=-0.5;

    cout<<AAA<<endl;

    AAA.decomposition_LDL(L,D);
    cout<<L<<endl;
    cout<<D<<endl;

    Vecteur PY(5);
    PY[0]=4;
    PY[1]=1;
    PY[2]=-8;
    PY[3]=6;
    PY[4]=0;
    cout<<PY<<endl;
    Vecteur Y=resolution_systeme_lineaire(AAA,PY);
    cout<<PY<<endl;
    cout<<Y<<endl;
    // programme
    // double h = 10.;
    // double e = 2.;
    // double pas = 0.001;
    // Maillage maillageFinale = genere_maillage_couche_diffusante(1, h, e, pas);
    // MatriceSym A;
    // Vecteur P;
    // genere_coefficient_matrice_A(A, maillageFinale, pas);
    // genere_coefficient_vecteur_P(P, maillageFinale, pas);
    // Vecteur solution = resolution_systeme_lineaire(A, P);
    return 0;
}
