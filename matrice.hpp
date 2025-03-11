#ifndef MATRICE_HPP_INCLUDED
#define MATRICE_HPP_INCLUDED
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
using namespace std;

typedef vector<double> Vecteur;
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
    MatriceSym (int m=0, double v=0.);
    MatriceSym (const Vecteur& d);
    double& operator() (int i, int j);
    double operator() (int i,int j) const;
    void decomposition_LDL( MatriceSym& L, Vecteur& D); 
    
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
    coefs.resize(n*(n+1)/2,0.);
    for(int i=0;i<n;i++){
        coefs[n*(n+1)/2 - i*(i+1)/2]=d[i];
        }
}



double& MatriceSym::operator() (int i, int j){
    if(i<0 || j<0 || j>=n || i >=n) {cout<<"coef(i,j) : i,j en dehors des bornes"<<endl; exit(-1);}
    if(j<i){
        return coefs[n*(n+1)/2 -i*(i+1)/2 + j - i];
     }
     return coefs[n*(n+1)/2 -j*(j+1)/2 + i - j];
}




double MatriceSym::operator() (int i, int j) const {
    if(i<0 || j<0 || j>=n || i >=n) {cout<<"coef(i,j) : i,j en dehors des bornes"<<endl; exit(-1);}
    if(j<i){
        return coefs[n*(n+1)/2 -i*(i+1)/2 + j - i];
     }
     return coefs[n*(n+1)/2 -j*(j+1)/2 + i - j];

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
void MatriceSym::decomposition_LDL( MatriceSym& L, Vecteur& D) {
    // init MatriceSym L et MatriceSym diagonale D
    D.resize(n);

    for (int i = 0; i <n; ++i) {
        // calcul D[i]
        D[i] = A(i, i);
        for (int k = 0; k < i; ++k)
        {
            D[i] -= L(i, k) * L(i, k) * D[i];
        }

        // Verification caractère défini positif
        if (D[i] == 0) {
            cout << "Matrice non définie positive!" << endl;
            exit(-1);
        }

        // calcule L(i,j) pour j < i
        for (int j = i + 1; j < n; ++j)
        {
            L(j, i) = A(j, i);
            for (int k = 0; k < i; ++k){
                L(j, i) -= L(j, k) * L(i, k) * D[k];
            }
            L(j, i) /= D[i];
        }

        // elements diagonaux de L sont 1
        L(i, i) = 1.0;
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
    Matrice L(n, 0);
    Vecteur D(n); //represente la matrice diagonale D
    decomposition_LDL(A,L,D);

    // resolution LY=P pour Y par recurrence simple grace a la forme finale de la methode de Gauss du systeme lineaire
    //(substition)
    Vecteur Y(n);
    // initialisation
    Y[1] = P[1];
    // heredite
    for (int i = 2; i < n; i++)
    {
        Y[i] += P[i];
        for (int k = 1; i < k - 1; k++)
        {
            Y[i] -= L(i, k) * Y[k];
        }
    }

    // resolution DL'Q=Y : recurrence simple decroissante sur l'indice
    Vecteur Q(n);
    // initialisation
    Q[1] = Y[1] / D[1];
    // heredite
    for (int i = n; i > 0; i--)
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
