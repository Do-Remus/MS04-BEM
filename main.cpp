#include <iostream>
#include "point.hpp"
#include "obstacle.hpp"
#include "segment.hpp"
using namespace std;


int main(){
    Point A(1,1);
    Point C(1,2);
    cout<<"A =("<<A.x<<","<<A.y<<")"<<endl;
    Obstacle B;
    cout<<"obstacle point "<<B.centre.x<<","<<B.centre.y<<endl;
    cout<<" obstacle rayon = "<<B.rayon<<endl;

    Segment S(A,C);
    cout<<"segment P1"<<S.P1.x<<","<<S.P1.y<<endl;
    cout<<"segment P2"<<S.P2.x<<","<<S.P2.y<<endl;

    return 0;
}

