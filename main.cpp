#include <iostream>
#include <time.h>
#include "point.hpp"
#include "cercle.hpp"
#include "cercle.hpp"
#include "segment.hpp"
#include "utils.hpp"
using namespace std;

int main()
{
    srand(time(NULL));
    Point A(1, 1);
    Point C(1, 2);
    cout << "A =(" << A.x << "," << A.y << ")" << endl;
    Cercle B;
    cout << "obstacle point " << B.centre.x << "," << B.centre.y << endl;
    cout << " obstacle rayon = " << B.rayon << endl;

    Segment S(A, C);
    cout << "segment P1: " << S.P1.x << "," << S.P1.y << endl;
    cout << "segment P2: " << S.P2.x << "," << S.P2.y << endl;

    cout << "A =" << A << endl;

    Maillage monMaillage = genere_maillage_couche_diffusante(5, 10, 5, 0.5);

    monMaillage.export_maillage("output.txt");

    return 0;
}
