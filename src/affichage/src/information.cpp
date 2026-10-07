#include "../information.hpp"

void affichage_parametres()
{
    affichage_ligne("Nombre d'onde k", k);
}

void affichage_approximations()
{
    affichage_ligne("Indice de troncature N", idxTroncature);
    affichage_ligne("Ordre de quadrature", ordre);
    affichage_ligne("Pas large du cache", green_cache_step_large);
    affichage_ligne("Pas fin du cache", green_cache_step_fine);
    affichage_ligne("Cutoff du cache", green_cache_cutoff);
    affichage_ligne("Taille du cache", green_cache_size);

    const Real memoire =
        green_cache_size * sizeof(Complex) / (1024.0 * 1024.0);

    affichage_ligne(
        "Memoire cache Green",
        std::to_string(memoire) + " MiB");
}

void affichage_maillage(const Point &centre, unsigned int nbSegments, const std::string &nomFichier)
{
    affichage_ligne("Centre", centre);
    affichage_ligne("Rayon", rayon);
    affichage_ligne("Pas", pasMaillage);
    affichage_ligne("Nombre de segments", nbSegments);
    affichage_ligne("Maillage exporte", nomFichier);
}

void affichage_erreur_p(Real erreur)
{
    affichage_ligne("Erreur L2 p", erreur);
}

void affichage_maillage_solution(unsigned int nbPointsGrille, unsigned int nbPointsSolution)
{
    affichage_ligne("Domaine", std::string("[") +
                                   std::to_string(-L / 2.0) + ", " +
                                   std::to_string(L / 2.0) + "] x [" +
                                   std::to_string(-L / 2.0) + ", " +
                                   std::to_string(L / 2.0) + "]");
    affichage_ligne("Taille du domaine", std::to_string(L) + " x " + std::to_string(L));
    affichage_ligne("Pas", pasSolution);
    affichage_ligne("Distance a la frontiere", delta);
    affichage_ligne("Rayon minimal autorise", rayon + delta);
    affichage_ligne("Points de grille", nbPointsGrille);
    affichage_ligne("Nombre de points", nbPointsSolution);
}

void affichage_erreurs_solution(Real erreurL2, Real erreurLinf)
{
    affichage_ligne("Erreur L2", erreurL2);
    affichage_ligne("Erreur Linf", erreurLinf);
}
