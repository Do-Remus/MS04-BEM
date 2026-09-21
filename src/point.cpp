#include "headers/point.hpp"

/* Fonction de la classe Point */

Point::Point()
{
    // basic constuctor
    this->x = 0;
    this->y = 0;
}

Point::Point(double a, double b)
{
    this->x = a;
    this->y = b;
}

Point &Point::operator+=(const Point &B)
{
    this->x += B.x;
    this->y += B.y;
    return *this;
}

Point &Point::operator-=(const Point &B)
{
    this->x -= B.x;
    this->y -= B.y;
    return *this;
}

Point &Point::operator*=(const double a)
{
    this->x *= a;
    this->y *= a;
    return *this;
}

Point &Point::operator/=(const double a)
{
    if (a == 0)
    {
        cout << "division par 0";
        exit(-1);
    }
    this->x /= a;
    this->y /= a;
    return *this;
}

double Point::norm() const
{
    return sqrt(this->x * this->x + this->y * this->y);
}


double Point::theta() const { return std::atan2(this->y, this->x); }

/* Fonctions associées à la classe Points */

ostream &operator<<(ostream &out, const Point &A)
{
    out << "(" << A.x << ";" << A.y << ")" << endl;
    return out;
}

Point operator+(const Point &A, const Point &B)
{
    Point R = A;
    R += B;
    return R;
}

Point operator-(const Point &A, const Point &B)
{
    Point R = A;
    R -= B;
    return R;
}

Point operator*(const Point &A, const double a)
{
    Point R = A;
    R *= a;
    return R;
}

Point operator*(const double a, const Point &A)
{
    Point R = A;
    R *= a;
    return R;
}

Point operator/(const Point &A, const double a)
{
    if(a==0){
        cout<<"division par 0 lors d'une division pour Point"<<endl;
        exit(-1);
    }
    Point R = A;
    R /= a;
    return R;
}

double operator*(const Point &A, const Point &B)
{
    return A.x * B.y - A.y * B.x;
}

bool operator==(const Point &A, const Point &B)
{
    if (A.x == B.x && A.y == B.y)
    {
        return true;
    }
    return false;
}

double operator|(const Point &A, const Point &B)
{
    return A.x * B.x + A.y * B.y;
}
