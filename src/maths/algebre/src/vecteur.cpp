#include "../vecteur.hpp"

/* Fonctions de la classe Vecteur */

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

Vecteur generateur_aleatoire_vecteur(const int taille)
{
    Vecteur v(taille);
    for (unsigned int i = 0; i < taille; i++)
    {
        v[i] = 100 * rand() / RAND_MAX;
    }
    return v;
}

/* Méthodes de la classe Vecteur */

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

Vecteur Vecteur::operator*(const Complex &a) const
{
    Vecteur resultat(size());

    for (std::size_t i = 0; i < size(); ++i)
        resultat[i] = a * (*this)[i];

    return resultat;
}

Real Vecteur::norm() const
{
    return std::sqrt(
        std::real(produitHermitien(*this)));
}

Complex Vecteur::produitHermitien(const Vecteur &v) const
{
    if (size() != v.size())
        throw std::invalid_argument("Tailles incompatibles");

    Complex resultat = 0.0;

    for (std::size_t i = 0; i < size(); ++i)
        resultat += std::conj((*this)[i]) * v[i];

    return resultat;
}

Complex Vecteur::produitBilineaire(const Vecteur &v) const
{
    if (size() != v.size())
        throw std::invalid_argument("Tailles incompatibles");

    Complex resultat = 0.0;

    for (std::size_t i = 0; i < size(); ++i)
        resultat += (*this)[i] * v[i];

    return resultat;
}
