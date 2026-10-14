#include "../matrice.hpp"

/* Fonctions associées à la classe Matrice */

Matrice matrice_aleatoire(const unsigned int rang, const Index taille)
{
    if (rang > taille)
        throw std::invalid_argument("Le rang ne peut pas dépasser la taille.");

    if (rang == 0 || taille == 0)
        return Matrice::Zero(taille, taille);

    Matrice U = Matrice::Random(taille, rang);
    Matrice V = Matrice::Random(rang, taille);

    return U * V;
}



vector<Vecteur> full_pivoted_cross_approximation(const Matrice& A, const Index rang, double tolerance){
    //initialisation de la matrice R et A localement
    Matrice Aiter = Matrice::Zero(A.rows(),A.cols());
    Matrice Riter(A);

    // first [rang-1] vectors are the u_k, then v_k are from [rang] to [2*rang -1]
    vector<Vecteur> U_V(2*rang);

    //boucle principale
    int k=0;
    double Anorm = A.norm();
    while(k<rang && Riter.norm() >= tolerance*Anorm){
        //trouver le pivot (valeur maximale en valeur absolue )
        Index i_pivot =0;
        Index j_pivot=0;
        double abs_value_pivot=-1;
        for(Index i=0; i<A.rows(); i++){
            for(Index j=0; j<A.cols(); j++){
                double a = abs(Riter(i,j));
                if(a > abs_value_pivot){
                    abs_value_pivot = a;
                    i_pivot=i;
                    j_pivot =j;
                }
            }
        }

        //extraire u_k+1 et v_k+1, et iteration Aiter, Riter
       U_V[k+rang] = Riter.row(i_pivot);
       U_V[k] = Riter.col(j_pivot)/Riter(i_pivot, j_pivot);
       Aiter = Aiter + U_V[k]*U_V[k+rang].transpose();
       cout<<"||A -Aiter|| = "<<(A-Aiter).norm()<<endl;
       Riter = Riter - U_V[k]*U_V[k+rang].transpose();
       k++;
    }
    return U_V;

}


ACA_decomposition partial_pivoted_cross_approximation(const Matrice& A, const Index rang, double tolerance){
    // first [rang-1] vectors are the u_k, then v_k are from [rang] to [2*rang -1]
    ACA_decomposition sol;
    sol.U.reserve(rang);
    sol.V.reserve(rang);
    //vector<Vecteur> U_V(2*rang);
    //boucle principale
    Index k=0;
    Index next_ipivot_index =0;
    double Anorm = A.norm();
    double tol_sq = tolerance*tolerance;
    Index Arows = A.rows();
    Index Acols = A.cols();
    vector<int> row_used(Arows,0);
    row_used[0] =1; //first line will be used to as first iteration
    double Aiter_Frobenius_norm = 0.;
    double new_rank1_contribution = 1.;
    while(k<=rang && new_rank1_contribution >= tol_sq*Aiter_Frobenius_norm){
        //init pour les vecteurs u_k+1 et v_k+1
        Vecteur u = Vecteur::Zero(Arows);
        Vecteur v = Vecteur::Zero(Acols);
        //trouver le pivot  (valeur maximale en valeur absolue) !!!! sur la ligne k !!!!
        //i pivot is next_ipivot_index, j_pivot is searched
        Index j_pivot=0;
        double abs_value_pivot=-1;
        complex<double> value_pivot = 0.;
        for(Index j=0; j<Acols; j++){
            complex<double> Riter_ipivot_j = A(next_ipivot_index,j);
            for(Index l=0; l< k; l++){
                Riter_ipivot_j -= sol.U[l][next_ipivot_index] * sol.V[l][j];
            }
            //extraction du coeff Riter_kj = Riter_ipivot_j la k ieme ligne de Rk pour definir v_k+1
            v[j] = Riter_ipivot_j;
            
            //check maximalité de la norme du coeff
            double a = abs(Riter_ipivot_j);
            if(a > abs_value_pivot){
                abs_value_pivot = a;
                value_pivot = Riter_ipivot_j;
                j_pivot =j;
            }
        }
        if(abs_value_pivot < tol_sq){
            return sol;
        }
       //calcul u_k+1 sans construire Riter 
       for(Index i=0; i< Arows ; i++){
            complex<double> Riter_i_jpivot = A(i,j_pivot);
            for(Index l=0; l<k; l++){
                Riter_i_jpivot -=sol.U[l][i] * sol.V[l][j_pivot];
            }
            u[i] = Riter_i_jpivot / value_pivot;
       }

       //push de u et v dans le struct 
        sol.U.push_back(u);
        sol.V.push_back(v);
       //calcul iteratif de la norme de Frobenius pour Ak
       complex<double> sum=0.;
       for(Index l=0; l<k; l++){
            //complex<double> uTu = U_V[k].transpose()*U_V[l].conjugate();
            //complex<double> vTv = U_V[l+rang].transpose()*U_V[k+rang].conjugate();
            sum += sol.U[l].dot(sol.U[k]) * sol.V[l].dot(sol.V[k]);
       }
       new_rank1_contribution =  sol.U[k].squaredNorm()* sol.V[k].squaredNorm();
       cout<<"k ="<<k<<", new_rank1_contribution = "<<new_rank1_contribution<<endl;
       Aiter_Frobenius_norm += 2*sum.real() + new_rank1_contribution ;
       
       //cout<<"abs(||A|| -||Aiter||) = "<<abs(Anorm - Aiter_Frobenius_norm)<<endl;
       //Riter = Riter - U_V[k]*U_V[k+rang].transpose();
      
       //Find best index of line for next iteration
       double best=0.;
       for(Index i=0; i<Arows; i++){
        if(row_used[i]==0 && abs(sol.U[k][i]) > best){
            best = abs(sol.U[k][i]);
            next_ipivot_index = i;
        }
       }
       //add next_ipivot_index to used rows
       row_used[next_ipivot_index] =1;
       //k increment for while loop
        k++;
    }
    return sol;

}

// ostream &operator<<(ostream &out, const Matrice &A)
// {
//     for (int i = 0; i < A.n; i++)
//     {
//         for (int j = 0; j < A.m; j++)
//         {
//             out << A(i, j) << " ";
//         }
//         out << endl;
//     }
//     return out;
// }

// /* Méthodes de la classe Matrice */

// Matrice::Matrice(int N, int M, complex<Real> v)
// {
//     m = M;
//     n = N;
//     coefs.resize(n * m, v);
//     // constructeur matrice a valeurs constantes
// }

// complex<Real> &Matrice::operator()(int i, int j)
// {
//     if (i >= n || j >= m || i < 0 || j < 0)
//     {
//         cout << "overflow in matrix coordinates ; i= " << i << ", j= " << j << " matrix size (n,m) = (" << n << "," << m << ")" << endl;
//         std::exit(1);
//     }
//     return coefs[i * m + j];
// }

// complex<Real> Matrice::operator()(int i, int j) const
// {
//     if (i >= n || j >= m || i < 0 || j < 0)
//     {
//         cout << "overflow in matrix coordinates ; i= " << i << ", j= " << j << " matrix size (n,m) = (" << n << "," << m << ")" << endl;
//         std::exit(1);
//     }
//     return coefs[i * m + j];
// }

// Vecteur Matrice::operator()(int i) const
// {
//     if (i >= n || i < 0)
//     {
//         cout << "overflow in matrix coordinates ; i= " << i << " n lines from matrix size (n,m) = (" << n << "," << m << ")" << endl;
//         std::exit(1);
//     }
//     Vecteur ligne(m, 0);
//     for (int j = 0; j < m; j++)
//         ligne[j] = (*this)(i, j);
//     return ligne;
// }

// Vecteur Matrice::operator*(const Vecteur &v) const
// {
//     Vecteur b(n);
//     for (int i = 0; i < n; i++)
//     {
//         b[i] = 0;
//         for (int j = 0; j < m; j++)
//         {
//             b[i] += (*this)(i, j) * v[j];
//         }
//     }
//     return b;
// }
