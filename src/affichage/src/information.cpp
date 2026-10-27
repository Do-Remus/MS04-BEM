#include "../information.hpp"

void affichage_parametres()
{
    affichage_ligne("Nombre d'onde k", k);
}

void affichage_approximations()
{
    affichage_ligne("Indice de troncature N", idxTroncature);
    affichage_ligne("Ordre de quadrature", ordre);
    affichage_ligne("Rang maximal ACA", rang_max);
    affichage_ligne("Tolerance ACA", toleranceAca);
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

void affichage_erreur_relative_matrice(Real erreur)
{
    affichage_ligne("Erreur relative matrice H", erreur);
}

void affichage_maillage(const Point &centre, unsigned int nbSegments, const std::string &nomFichier)
{
    affichage_ligne("Centre", centre);
    affichage_ligne("Rayon", rayon);
    affichage_ligne("Pas", pasMaillage);
    affichage_ligne("Nombre de segments", nbSegments);
    affichage_ligne("Maillage exporte", nomFichier);
}

void affichage_erreur_relative_p(Real erreur)
{
    affichage_ligne("Erreur relative L2 p", erreur);
}

void affichage_maillage_solution(unsigned int nbPointsSolution)
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
    affichage_ligne("Nombre de points", nbPointsSolution);
}

void affichage_erreurs_solution(Real erreurL2, Real erreurLinf)
{
    affichage_ligne("Erreur L2", erreurL2);
    affichage_ligne("Erreur Linf", erreurLinf);
}

void affichage_statistiques_matrice_H(const ArbreBlocs &arbre)
{
    if (!arbre)
    {
        affichage_ligne("Matrice H", "Arbre de blocs vide");
        return;
    }

    struct Statistiques
    {
        std::size_t nb_blocs_internes = 0;
        std::size_t nb_blocs_admissibles = 0;
        std::size_t nb_blocs_denses = 0;

        std::size_t coefficients_denses = 0;
        std::size_t coefficients_UV = 0;

        std::size_t rang_maximal = 0;
        std::size_t somme_rangs = 0;
    };

    Statistiques stats;

    std::function<void(const Bloc &)> parcourir = [&](const Bloc &bloc)
    {
        switch (bloc.type)
        {
        case TypeBloc::Interne:
            ++stats.nb_blocs_internes;

            for (const auto &ligne : bloc.fils)
            {
                for (const auto &fils : ligne)
                {
                    if (fils)
                        parcourir(*fils);
                }
            }
            break;

        case TypeBloc::Dense:
            ++stats.nb_blocs_denses;
            stats.coefficients_denses +=
                static_cast<std::size_t>(bloc.dense.size());
            break;

        case TypeBloc::Admissible:
        {
            ++stats.nb_blocs_admissibles;

            const std::size_t rang =
                static_cast<std::size_t>(bloc.U.cols());

            stats.rang_maximal =
                std::max(stats.rang_maximal, rang);

            stats.somme_rangs += rang;

            stats.coefficients_UV +=
                static_cast<std::size_t>(bloc.U.size()) +
                static_cast<std::size_t>(bloc.V.size());

            break;
        }
        }
    };

    parcourir(*arbre);

    const std::size_t N = arbre->cluster_lignes ? arbre->cluster_lignes->indices.size() : 0;

    const std::size_t coefficients_H = stats.coefficients_denses + stats.coefficients_UV;
    const long double coefficients_denses_complets = static_cast<long double>(N) * N;
    const long double octets_par_coefficient = sizeof(Matrice::Scalar);
    const long double memoire_H = coefficients_H * octets_par_coefficient;
    const long double memoire_dense = coefficients_denses_complets * octets_par_coefficient;
    const double memoire_H_MiB = static_cast<double>(memoire_H / (1024.0L * 1024.0L));
    const double memoire_dense_MiB = static_cast<double>(memoire_dense / (1024.0L * 1024.0L));
    const double pourcentage_stockage = coefficients_denses_complets > 0 ? 100.0 * static_cast<double>(coefficients_H) / static_cast<double>(coefficients_denses_complets) : 0.0;
    const double rang_moyen = stats.nb_blocs_admissibles > 0 ? static_cast<double>(stats.somme_rangs) / static_cast<double>(stats.nb_blocs_admissibles) : 0.0;

    affichage_ligne("Taille de la matrice", std::to_string(N) + " x " + std::to_string(N));

    affichage_ligne("Blocs internes", stats.nb_blocs_internes);
    affichage_ligne("Blocs admissibles", stats.nb_blocs_admissibles);
    affichage_ligne("Blocs denses", stats.nb_blocs_denses);

    affichage_ligne("Coefficients denses stockes", stats.coefficients_denses);
    affichage_ligne("Coefficients des facteurs U et V", stats.coefficients_UV);
    affichage_ligne("Nombre total de coefficients stockes", coefficients_H);

    affichage_ligne("Rang maximal ACA", stats.rang_maximal);
    affichage_ligne("Rang moyen ACA", std::to_string(rang_moyen));

    affichage_ligne("Taux de stockage", std::to_string(pourcentage_stockage) + " %");
    affichage_ligne("Memoire estimee matrice H", std::to_string(memoire_H_MiB) + " MiB");
    affichage_ligne("Memoire matrice dense", std::to_string(memoire_dense_MiB) + " MiB");
}
