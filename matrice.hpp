#ifndef MATRICE_HPP_INCLUDED
#define MATRICE_HPP_INCLUDED
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
using namespace std;



typedef vector<double> Vecteur;

class MatriceSym
{private :
    Vecteur coefs;                       // coefficients rangés par ligne
public:
    int n;                               // ordre de la MatriceSym
    MatriceSym (int m=0, double v=0.);      // constructeur dimensions et coefs constants
    MatriceSym (const Vecteur& d);          // constructeur d'une MatriceSym diagonale
    MatriceSym (const initializer_list<Vecteur>& vs); // constructeur d'une MatriceSym à partir de vecteurs
    double& coef(int i,int j) ;          // acces au coef i,j (1->n)
    double  coef(int i, int j) const ;   // acces au coef i,j (1->n) const
    double& operator() (int i, int j);
    double operator() (int i,int j) const;
    Vecteur operator() (int i) const;


};

MatriceSym::MatriceSym (int m, double v)      // constructeur dimensions et coefs constants
{
    n=max(m,0);
    if(m==0) {return;}
    coefs.resize(n*n,v);
}
MatriceSym::MatriceSym (const Vecteur& d)     // constructeur d'une MatriceSym diagonale
{
    n=d.size();
    coefs.resize(n*n,0.);
    for(int i=0;i<n;i++)
        coefs[i*n+i]=d[i];
}
MatriceSym::MatriceSym(const initializer_list<Vecteur>& vs) // constructeur à partir de vecteurs
{
    n=vs.size();
    coefs.resize(n*n);
    int k=0;
    for(auto& li : vs)
    {
        if(li.size()!=size_t(n)) {cout<<"taille d'un vecteur incompatible"<<endl; exit(-1);}
        for(int i=0;i<n;i++,k++) coefs[k]=li[i];
    }
    for(int i=0; i<n;i++){
        for(int j=0;j<n;j++){
            if(this->coefs[(i-1)*this->n+j-1]!=this->coefs[(j-1)*this->n+i-1]){
                cout<<"matrice non symetrique"<<endl;
                exit(-1);
            }
        }
    }
    
}
double& MatriceSym::coef(int i,int j)       // acces au coef i,j (1->n)
{
     if(i<=0 || j<=0) {cout<<"coef(i,j) : i,j en dehors des bornes"<<endl; exit(-1);}
     return coefs[(i-1)*n+j-1];
}

double& MatriceSym::operator() (int i, int j){
      if(i<=0 || j<=0 || j>=n || i >=this->n) {cout<<"coef(i,j) : i,j en dehors des bornes"<<endl; exit(-1);}
     return this->coefs[(i-1)*this->n+j-1];
}




double  MatriceSym::coef(int i, int j) const // acces au coef i,j (1->n) const
{
     if(i<=0 || j<=0) {cout<<"coef(i,j) : i,j en dehors des bornes"<<endl; exit(-1);}
     return coefs[(i-1)*n+j-1];
}

double MatriceSym::operator() (int i, int j) const {
      if(i<=0 || j<=0 || j>=n || i >=this->n) {cout<<"coef(i,j) : i,j en dehors des bornes"<<endl; exit(-1);}
     return this->coefs[(i-1)*this->n+j-1];

}





#endif
