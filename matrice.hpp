#ifndef MATRICE_HPP_INCLUDED
#define MATRICE_HPP_INCLUDED
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <complex>
using namespace std;

typedef vector<complex<double>> Vecteur;
ostream &operator<<(ostream &out, const Vecteur &u)
{
    int n = u.size();
    out << '(' << u[0];
    for (int i = 1; i < n; i++)
    {
        out << "," << u[i];
    }
    out << ") \n";
    return out;
}

class MatriceSym
{protected :
    Vecteur coefs;                       // coefficients rangés par ligne
public:
    int n;                               // ordre de la MatriceSym
    MatriceSym (int m=0, complex<double> v= (complex<double>) 0);
    MatriceSym (const Vecteur& d);
    complex<double>& operator() (int i, int j);
    complex<double> operator() (int i,int j) const;
    void decomposition_LDL( MatriceSym& L, Vecteur& D) const; 
};

MatriceSym::MatriceSym (int m, complex<double> v)      // constructeur dimensions et coefs constants
{
    n=max(m,0);
    if(m==0) {return;}
    coefs.resize((n*(n+1))/2,v);
}

MatriceSym::MatriceSym (const Vecteur& d)     // constructeur d'une MatriceSym diagonale
{
    n=d.size();
    coefs.resize(n*(n+1)/2,0.);
    for(int i=0;i<n;i++){
        coefs[(n*(n+1))/2 - i*(i+1)/2]=d[i];
        }
}



complex<double>& MatriceSym::operator() (int i, int j){
    if(i<0 || j<0 || j>=n || i >=n) {cout<<"coef(i,j) :" <<i<<","<<j<<" en dehors des bornes"<<endl; exit(-1);}
    if(j<i){
        return coefs[i*(i+1)/2 + j];
     }
     return coefs[j*(j+1)/2 + i];
}




complex<double> MatriceSym::operator() (int i, int j) const {
    if(i<0 || j<0 || j>=n || i >=n) {cout<<"coef(i,j) :" <<i<<","<<j<<" en dehors des bornes"<<endl; exit(-1);}
    if(j<i){
        return coefs[i*(i+1)/2 + j];
     }
     return coefs[j*(j+1)/2 + i ];

}

ostream& operator <<(ostream& out, const MatriceSym& A){
    for(int i=0; i<A.n;i++){
        for(int j=0; j<A.n; j++){
            out<<A(i,j)<<" ";
        }
        out << endl;
    }
    return out;
}



//factorisation LDL
void MatriceSym::decomposition_LDL( MatriceSym& L, Vecteur& D) const {
    // init MatriceSym L et MatriceSym diagonale D
    if (L.n != n){
        cout<<"dimension problem with L.n ="<<L.n<<"and A.n ="<<n<<endl;
        exit(-1);
     } //bonne taille pour L (matrice triangulaire inferieure assimilée à matrice sym


    for (int j = 0; j <n; ++j) {
        // elements diagonaux de L sont 1
        L(j, j) =  1.0;
        // calcul D[i]
        D[j] = (*this)(j, j);
        for (int k = 0; k < j; ++k)
        {
            D[j] -= L(j, k) * L(j, k) * D[k];
        }
        // Verification caractère défini positif

        if (D[j] == (complex<double>) 0) {
            cout << "Element diagonal nul pour j="<<j<< endl;
            exit(-1);
        }

        // calcule L(i,j) pour i<j<n
        for (int i = j + 1; i < n; ++i)
        {
            L(i, j) = (*this)(i, j);
            for (int k = 0; k < j; ++k){
                L(i, j) -= L(j, k) * L(i, k) * D[k];
            }
            L(i, j) /= D[j];
        }
        


    }

    return;
}

Vecteur resolution_systeme_lineaire(const MatriceSym& A, const Vecteur& P){
    //va resoudre pour Q AQ=P en 3 parties 
    //1 -- A=LDL'
    //2 -- resolution LY=P pour Y
    //3 -- resolution DL'Q=Y

    int n=A.n;
    //decomposition de A=LDL'
    MatriceSym L(n, 0);
    Vecteur D(n); //represente la matrice diagonale D
    A.decomposition_LDL(L,D);

    // resolution LY=P pour Y par recurrence simple grace a la forme finale de la methode de Gauss du systeme lineaire
    //(substition)
    Vecteur Y(n);
    // initialisation
    Y[0] = P[0];
    // heredite
    for (int i = 1; i < n; i++)
    {
        Y[i] += P[i];
        for (int k = 0; k<i; k++)
        {
            Y[i] -= L(i, k) * Y[k];
        }
    }

    // resolution DL'Q=Y : recurrence simple decroissante sur l'indice
    Vecteur Q(n);
    // initialisation
    Q[0] = Y[0] / D[0];
    // heredite
    for (int i = n-1; i > 0; i--)
    {
        Q[i] += Y[i] / D[i];
        for (int k = i + 1; k < n; k++)
        {
            Q[i] -= L(k, i) * Q[k]; // attention L'(i,k)=L(k,i)
        }
    }
    // solution
    return Q;
}

#endif
