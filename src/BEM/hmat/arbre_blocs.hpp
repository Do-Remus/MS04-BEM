#ifndef BEM_ARBRE_BLOCS_HPP_INCLUDED
#define BEM_ARBRE_BLOCS_HPP_INCLUDED

#include "../../global/commun.hpp"
#include "../../maths/commun.hpp"

#include "arbre_clusters.hpp"
#include "../operator.hpp"
#include "../aca_based/decomposition.hpp"

enum class TypeBloc
{
    Interne,    // Possède quatre blocs fils
    Admissible, // Bloc admissible : compression ACA possible
    Dense       // Bloc non admissible : stockage dense
};

struct Bloc
{
    const Noeud *cluster_lignes = nullptr;
    const Noeud *cluster_colonnes = nullptr;

    TypeBloc type = TypeBloc::Interne;

    int proprietaire_mpi = 0;
    Index nb_feuilles_rank_mpi = 0;

    double poids = 0.0;

    // Le bloc est-il sur la diagonal ?
    bool est_diagonal = false;

    // Pour le stockage de la matrice
    Matrice dense;
    Matrice U;
    Matrice V;

    std::unique_ptr<Bloc> fils[2][2];
};

using ArbreBlocs = std::unique_ptr<Bloc>;

struct StatistiquesBlocs
{
    Index feuilles = 0;
    Index denses = 0;
    Index admissibles = 0;

    std::uint64_t coefficients_denses = 0;
    std::uint64_t coefficients_aca = 0;
};

bool sont_admissibles(const BoiteEnglobante &boite_a, const BoiteEnglobante &boite_b, const double eta);

ArbreBlocs construire_bloc_recursif(const Noeud &tau, const Noeud &sigma, const double eta, const double alpha);

void collecter_feuilles(Bloc &bloc, std::vector<Bloc *> &feuilles);

void repartir_feuilles_ponderee(Bloc &racine, const int mpi_size);

Index recalculer_feuilles_locales(Bloc &bloc, const int mpi_rank);

ArbreBlocs construire_arbre_blocs(const Arbre &arbre_clusters, const double eta, const double alpha);

void remplir_blocs(Bloc &bloc, const Maillage &maillage, const std::vector<QuadratureSegment> &quadrature_maillage, const Index rang_max, const double tolerance);

void ecrire_arbre_blocs(const ArbreBlocs &arbre, const std::string &fichier, const IntervallesClusters &intervalles);

void compter_blocs_locaux(const Bloc &bloc, const int rang, StatistiquesBlocs &stats);

#endif