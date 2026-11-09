#ifndef BEM_FMM_STRUCTURE_HPP_INCLUDED
#define BEM_FMM_STRUCTURE_HPP_INCLUDED

#include "../../global/commun.hpp"
#include "../../maths/commun.hpp"

#include "../hmat/arbre_clusters.hpp"

constexpr Index INDICE_INVALIDE = std::numeric_limits<Index>::max();

struct InteractionFMM
{
    Index cible;
    Index source;
};

struct BlocDenseFMM
{
    Index cible;
    Index source;

    // Lignes : indices de la cible
    // Colonnes : indices de la source
    Matrice coefficients;
};

struct TranslationFMM
{
    Index cible;
    Index source;

    // Coefficients indexés par q dans [-2p, 2p]
    Vecteur coefficients;
};

struct NoeudFMM
{
    Index parent;
    Index fils_gauche;
    Index fils_droit;

    BoiteEnglobante boite;
    Point centre;
    Real rayon;
    Index niveau;

    std::vector<Index> indices;

    Vecteur multipole;
    Vecteur local;

    bool est_feuille() const
    {
        return fils_gauche == INDICE_INVALIDE && fils_droit == INDICE_INVALIDE;
    }
};

struct StructureFMM
{
    std::vector<NoeudFMM> noeuds;
    Index racine;

    std::vector<InteractionFMM> interactions_proches;
    std::vector<InteractionFMM> interactions_M2L;

    std::vector<BlocDenseFMM> blocs_denses;

    // Opérateurs de translation pré-calculés
    std::vector<TranslationFMM> translations_M2M;
    std::vector<TranslationFMM> translations_M2L;
    std::vector<TranslationFMM> translations_L2L;

    // Fonctions de base aux points de quadrature
    // Une ligne par point de quadrature
    // Une colonne par ordre n dans [-p,p]
    Matrice bases_points;

    // Décalage du premier point de quadrature de chaque segment
    std::vector<Index> offsets_quadratures;

    Real k;
    Index ordre;
    Real eta;

    Index nombre_niveaux;
};

#endif