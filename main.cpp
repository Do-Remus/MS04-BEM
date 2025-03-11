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

double f(const Point &A)
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

    // test Matrice

    // test MatriceSym

    // test decomposition LDL

    // programme
    double h = 10.;
    double e = 2.;
    double pas = 0.001;
    Maillage maillageFinale = genere_maillage_couche_diffusante(1, h, e, pas);
    MatriceSym A;
    Vecteur P;
    genere_coefficient_matrice_A(A, maillageFinale, pas);
    genere_coefficient_vecteur_P(P, maillageFinale, pas);
    Vecteur solution = resolution_systeme_lineaire(A, P);
    return 0;
}
