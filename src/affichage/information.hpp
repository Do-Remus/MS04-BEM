#ifndef AFFICHAGE_INFORMATION_HPP_INCLUDED
#define AFFICHAGE_INFORMATION_HPP_INCLUDED

#include "../global/commun.hpp"
#include "../maths/commun.hpp"

template <typename T>
void affichage_ligne(const std::string &nom, const T &valeur)
{
    std::cout << "    "
              << std::left << std::setw(LONGEUR_TEXT_AFFICHAGE)
              << nom
              << " : "
              << valeur
              << std::endl;
}

void affichage_parametres();

void affichage_approximations();

void affichage_maillage(const Point &centre, unsigned int nbSegments, const std::string &nomFichier);

void affichage_erreur_p(Real erreur);

void affichage_maillage_solution(unsigned int nbPointsGrille, unsigned int nbPointsSolution);

void affichage_erreurs_solution(Real erreurL2, Real erreurLinf);

#endif