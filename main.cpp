#include <iostream>
#include <time.h>
#include "point.hpp"
#include "cercle.hpp"
#include "cercle.hpp"
#include "segment.hpp"
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
    cout << "segment P1" << S.P1.x << "," << S.P1.y << endl;
    cout << "segment P2" << S.P2.x << "," << S.P2.y << endl;
    //test ostream function
    cout<<"A ="<<A<<endl;
    //test operations
    cout<<"A+C"<<A+C;
    cout<<"A-C"<<A-C;
    cout<<"3*C"<<3*C;
    cout<<"C/5"<<C/5;
    cout<<"norme C"<<C.norm()<<endl;
    cout<<"produit vectoriel composante z ="<<A^C
    

    return 0;
}

