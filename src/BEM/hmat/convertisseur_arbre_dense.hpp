#ifndef CONVERTISSEUR_ARBRE_DENSE_HPP_INCLUDED
#define CONVERTISSEUR_ARBRE_DENSE_HPP_INCLUDED
#include "../../global/commun.hpp"
#include "../../maths/commun.hpp"

#include "arbre_clusters.hpp"
#include "arbre_blocs.hpp"

void ajouter_bloc_dense(const Bloc &bloc, Matrice &A, const bool est_symetrique);

Matrice convertion_aca_dense(const ArbreBlocs &arbre_blocs, const bool est_symetrique = true);

void diagnostic_blocs(const Bloc &bloc, const Matrice &A_ref, const bool est_symetrique);
#endif