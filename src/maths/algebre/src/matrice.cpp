#include "../matrice.hpp"

/* Fonctions associées à la classe Matrice */

ostream &operator<<(ostream &out, const Matrice &A)
{
    for (int i = 0; i < A.n; i++)
    {
        for (int j = 0; j < A.m; j++)
        {
            out << A(i, j) << " ";
        }
        out << endl;
    }
    return out;
}




Matrice generateur_aleatoire_matrice(const int rang, const int taille)
{
    if (rang > taille)
    {
        cout << "rang est plus grand que la taille" << endl;
        std::exit(1);
    }
    Matrice mat(taille, taille, 0);
    unsigned int i = 0;
    while (i < rang)
    {
        // generateur d'une ligne
        Vecteur v = generateur_aleatoire_vecteur(taille);

        // verification de liberté des vecteurs v1,..,v_r
        for(unsigned int j=0; j<i; j++){
            if( abs((v-mat(i)).produitBilineaire(mat(i))) < 1e-8){
                v[j] += 1;
            }
        }
        //ecriture de la ligne générée dans la matrice
        for(unsigned int j=0; j<taille; j++){
            mat(i,j) = v[j];
        }
        i++;
    }
    // generateur aléatoire de combinaisons linéaires pour completer la matrice
    for(unsigned int i=rang; i<taille; i++){
        Vecteur coefs_combinaison_lineaire = generateur_aleatoire_vecteur(rang);
        for(unsigned int j=0; j<rang; j++){
            for(unsigned int k =0; k<taille; k++){
                mat(i,k)+=  coefs_combinaison_lineaire[j]* mat(j,k);
            }
        }
    }
    // generateur aléatoire de permutation des colonnes (ou lignes)
    return mat;
}

/* Méthodes de la classe Matrice */

Matrice::Matrice(int N, int M, complex<Real> v)
{
    m = M;
    n = N;
    coefs.resize(n * m, v);
    // constructeur matrice a valeurs constantes
}

complex<Real> &Matrice::operator()(int i, int j)
{
    if (i >= n || j >= m || i < 0 || j < 0)
    {
        cout << "overflow in matrix coordinates ; i= " << i << ", j= " << j << " matrix size (n,m) = (" << n << "," << m << ")" << endl;
        std::exit(1);
    }
    return coefs[i * m + j];
}

complex<Real> Matrice::operator()(int i, int j) const
{
    if (i >= n || j >= m || i < 0 || j < 0)
    {
        cout << "overflow in matrix coordinates ; i= " << i << ", j= " << j << " matrix size (n,m) = (" << n << "," << m << ")" << endl;
        std::exit(1);
    }
    return coefs[i * m + j];
}

Vecteur Matrice::operator()(int i) const
{
    if (i >= n || i < 0)
    {
        cout << "overflow in matrix coordinates ; i= " << i << " n lines from matrix size (n,m) = (" << n << "," << m << ")" << endl;
        std::exit(1);
    }
    Vecteur ligne(m, 0);
    for (int j = 0; j < m; j++)
        ligne[j] = (*this)(i, j);
    return ligne;
}


Vecteur Matrice::operator*(const Vecteur &v) const
{
    Vecteur b(n);
    for (int i = 0; i < n; i++)
    {
        b[i] = 0;
        for (int j = 0; j < m; j++)
        {
            b[i] += (*this)(i, j) * v[j];
        }
    }
    return b;
}
