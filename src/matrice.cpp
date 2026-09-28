#include "headers/matrice.hpp"

/* Fonctions associées au type Vecteur */

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

Vecteur Vecteur::operator+(const Vecteur &v) const
{
    Vecteur resultat(size());

    for (std::size_t i = 0; i < size(); ++i)
        resultat[i] = (*this)[i] + v[i];

    return resultat;
}

Vecteur Vecteur::operator-(const Vecteur &v) const
{
    Vecteur resultat(size());

    for (std::size_t i = 0; i < size(); ++i)
        resultat[i] = (*this)[i] - v[i];

    return resultat;
}

Vecteur Vecteur::operator*(const std::complex<double> &a) const
{
    Vecteur resultat(size());

    for (std::size_t i = 0; i < size(); ++i)
        resultat[i] = a * (*this)[i];

    return resultat;
}

double Vecteur::norm() const
{
    return std::sqrt(
        std::real(produitHermitien(*this)));
}

std::complex<double> Vecteur::produitHermitien(const Vecteur &v) const
{
    if (size() != v.size())
        throw std::invalid_argument("Tailles incompatibles");

    std::complex<double> resultat = 0.0;

    for (std::size_t i = 0; i < size(); ++i)
        resultat += std::conj((*this)[i]) * v[i];

    return resultat;
}

std::complex<double> Vecteur::produitBilineaire(const Vecteur &v) const
{
    if (size() != v.size())
        throw std::invalid_argument("Tailles incompatibles");

    std::complex<double> resultat = 0.0;

    for (std::size_t i = 0; i < size(); ++i)
        resultat += (*this)[i] * v[i];

    return resultat;
}

/* Fonctions de la classe Matrice */

Matrice::Matrice(int N, int M, complex<double> v)
{
    m = M;
    n = N;
    coefs.resize(n * m, v);
    // constructeur matrice a valeurs constantes
}

complex<double> &Matrice::operator()(int i, int j)
{
    if (i > n || j > m || i < 1 || j < 1)
    {
        cout << "overflow in matrix coordinates ; i= " << i << ", j= " << j << " matrix size (n,m) = (" << n << "," << m << ")" << endl;
        exit(1);
    }
    return coefs[(i - 1) * m + j - 1];
}

complex<double> Matrice::operator()(int i, int j) const
{
    if (i > n || j > m || i < 1 || j < 1)
    {
        cout << "overflow in matrix coordinates ; i= " << i << ", j= " << j << " matrix size (n,m) = (" << n << "," << m << ")" << endl;
        exit(1);
    }
    return coefs[(i - 1) * m + j - 1];
}

Vecteur Matrice::operator*(const Vecteur v) const
{
    Vecteur b(n);
    for (int i = 0; i < n; i++)
    {
        b[i] = 0;
        for (int j = 0; j < m; j++)
        {
            b[i] += (*this)(i + 1, j + 1) * v[j];
            // cout<<" (*this)(i+1,j+1) * v[j] = "<<(*this)(i+1,j+1) <<"* "<<v[j]<<" = "<< b[i]<<endl;
        }
    }
    return b;
}

/*Fonctions associées à la classe Matrice*/

ostream &operator<<(ostream &out, const Matrice &A)
{
    for (int i = 1; i <= A.n; i++)
    {
        for (int j = 1; j <= A.m; j++)
        {
            out << A(i, j) << " ";
        }
        out << endl;
    }
    return out;
}

/* Fonctions de la classe MatriceSym */

MatriceSym::MatriceSym(int m, complex<double> v) // constructeur dimensions et coefs constants
{
    n = max(m, 0);
    if (m == 0)
    {
        return;
    }
    coefs.resize((n * (n + 1)) / 2, v);
}

MatriceSym::MatriceSym(const Vecteur &d)
{
    n = static_cast<int>(d.size());

    coefs.resize(n * (n + 1) / 2, 0.0);

    for (int i = 0; i < n; ++i)
    {
        (*this)(i, i) = d[i];
    }
}

complex<double> &MatriceSym::operator()(int i, int j)
{
    if (i < 0 || j < 0 || j >= n || i >= n)
    {
        cout << "coef(i,j) :" << i << "," << j << " en dehors des bornes" << endl;
        exit(-1);
    }
    if (j < i)
    {
        return coefs[i * (i + 1) / 2 + j];
    }
    return coefs[j * (j + 1) / 2 + i];
}

Vecteur MatriceSym::operator*(const Vecteur &v) const
{
    if (v.size() != static_cast<size_t>(n))
    {
        cout << "Erreur : dimensions incompatibles pour MatriceSym * Vecteur"
             << " : matrice = " << n << "x" << n
             << ", vecteur = " << v.size()
             << endl;
        exit(-1);
    }

    Vecteur b(n, 0.0);

    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            b[i] += (*this)(i, j) * v[j];
        }
    }

    return b;
}

complex<double> MatriceSym::operator()(int i, int j) const
{
    if (i < 0 || j < 0 || j >= n || i >= n)
    {
        cout << "coef(i,j) :" << i << "," << j << " en dehors des bornes" << endl;
        exit(-1);
    }
    if (j < i)
    {
        return coefs[i * (i + 1) / 2 + j];
    }
    return coefs[j * (j + 1) / 2 + i];
}

void MatriceSym::decomposition_LDL(MatriceSym &L, Vecteur &D) const
{
    // init MatriceSym L et MatriceSym diagonale D
    if (L.n != n)
    {
        cout << "dimension problem with L.n =" << L.n << "and A.n =" << n << endl;
        exit(-1);
    } // bonne taille pour L (matrice triangulaire inferieure assimilée à matrice sym

    for (int j = 0; j < n; ++j)
    {
        // elements diagonaux de L sont 1
        L(j, j) = 1.0;
        // calcul D[i]
        D[j] = (*this)(j, j);
        for (int k = 0; k < j; ++k)
        {
            D[j] -= L(j, k) * L(j, k) * D[k];
        }
        // Verification caractère défini positif

        if (D[j] == (complex<double>)0)
        {
            cout << "Element diagonal nul pour j=" << j << endl;
            exit(-1);
        }

        // calcule L(i,j) pour i<j<n
        for (int i = j + 1; i < n; ++i)
        {
            L(i, j) = (*this)(i, j);
            for (int k = 0; k < j; ++k)
            {
                L(i, j) -= L(j, k) * L(i, k) * D[k];
            }
            L(i, j) /= D[j];
        }
    }

    return;
}

/* Fonctions associées à la classe MatriceSym */

ostream &operator<<(ostream &out, const MatriceSym &A)
{
    for (int i = 0; i < A.n; i++)
    {
        for (int j = 0; j < A.n; j++)
        {
            out << A(i, j) << " ";
        }
        out << endl;
    }
    return out;
}

Vecteur resolution_systeme_lineaire(const MatriceSym &A, const Vecteur &P)
{
    // va resoudre pour Q AQ=P en 3 parties
    // 1 -- A=LDL'
    // 2 -- resolution LY=P pour Y
    // 3 -- resolution DL'Q=Y

    int n = A.n;
    // decomposition de A=LDL'
    MatriceSym L(n, 0);
    Vecteur D(n); // represente la matrice diagonale D
    A.decomposition_LDL(L, D);

    // resolution LY=P pour Y par recurrence simple grace a la forme finale de la methode de Gauss du systeme lineaire
    //(substition)
    Vecteur Y(n);
    // initialisation
    Y[0] = P[0];
    // heredite
    for (int i = 1; i < n; i++)
    {
        Y[i] += P[i];
        for (int k = 0; k < i; k++)
        {
            Y[i] -= L(i, k) * Y[k];
        }
    }

    // resolution DL'Q=Y : recurrence simple decroissante sur l'indice
    Vecteur Q(n);
    // initialisation
    Q[0] = Y[0] / D[0];
    // heredite
    for (int i = n - 1; i > 0; i--)
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

Vecteur gradConjMatSym(const MatriceSym &A, const Vecteur &b, double tol, unsigned int maxIter)
{
    const std::size_t n = b.size();

    Vecteur x(n, 0.0);

    Vecteur r = b - A * x;
    Vecteur d = r;

    const double bnorm = b.norm();
    double relativeResidual = 0.0;

    if (bnorm == 0.0)
        return x;

    std::complex<double> rho = r.produitBilineaire(r);

    for (unsigned int iter = 0; iter < maxIter; ++iter)
    {
        Vecteur Ad = A * d;
        const std::complex<double> denom = d.produitBilineaire(Ad);
        const double scale = d.norm() * Ad.norm();

        if (scale == 0 || std::abs(denom) < 1e-20 * scale)
        {
            std::cerr << "GCMS : breakdown, denominateur nul"
                      << " a l'iteration " << iter
                      << ", avec un résidu relatif de : "
                      << relativeResidual << std::endl;
            return x;
        }

        const std::complex<double> alpha = rho / denom;
        x = x + alpha * d;
        r = r - alpha * Ad;
        relativeResidual = r.norm() / bnorm;

        if (relativeResidual < tol)
        {
            std::cout << "GCMS converge en "
                      << iter + 1
                      << " iterations avec un résidu relatif de : "
                      << relativeResidual << std::endl;

            return x;
        }

        const std::complex<double> rhoNew = r.produitBilineaire(r);

        if (std::abs(rho) < 1e-30)
        {

            std::cerr << "COCG : breakdown de rho, avec résidu relatif de : "
                      << relativeResidual
                      << std::endl;
            return x;
        }

        const std::complex<double> beta = rhoNew / rho;

        d = r + beta * d;
        rho = rhoNew;
    }

    std::cout << "GCMS : nombre maximal d'iterations atteint."
              << std::endl;

    return x;
}
