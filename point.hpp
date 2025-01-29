#ifndef POINT_HPP_INCLUDED
#define POINT_HPP_INCLUDED


class Point
{
public:
    double x;
    double y;
    Point();
    Point(double x, double y);

};

Point::Point(){
    //basic constuctor
    this->x=0;
    this->y=0;
}

Point::Point(double a, double b){
    this->x=a;
    this->y=b;
}


#endif 
