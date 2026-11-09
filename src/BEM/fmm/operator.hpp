#ifndef BEM_FMM_OPERATOR_HPP_INCLUDED
#define BEM_FMM_OPERATOR_HPP_INCLUDED

#include "../../global/commun.hpp"
#include "../../maths/commun.hpp"

#include "../operator.hpp"
#include "structure.hpp"

inline Index indice_FMM(const int n, const Index ordre)
{
    return static_cast<Index>(n + static_cast<int>(ordre));
}

inline Index indice_translation(const int q, const Index ordre)
{
    return static_cast<Index>(q + 2 * static_cast<int>(ordre));
}

inline void reinitialiser_developpements_FMM(StructureFMM &fmm)
{
    for (NoeudFMM &noeud : fmm.noeuds)
    {
        noeud.multipole.setZero();
        noeud.local.setZero();
    }
}

inline Complex coefficient_translation_M2M(const int n, const int m, const Real k, const Point &deplacement)
{
    const int q = n - m;

    const Real r = deplacement.norm();
    const Real theta = std::atan2(deplacement.y, deplacement.x);

    return bessel_J(q, k * r) * phase(-q, theta);
}

inline Complex coefficient_translation_M2L(const int n, const int m, const Real k, const Point &deplacement)
{
    const int q = m - n;

    const Real r = deplacement.norm();
    const Real theta = std::atan2(deplacement.y, deplacement.x);

    return Complex(0.0, 0.25) * hankel_H1(q, k * r) * phase(q, theta);
}

inline bool admissibles_M2L(const NoeudFMM &A, const NoeudFMM &B, const Real eta)
{
    const Real distance = (A.centre - B.centre).norm();

    return distance > eta * (A.rayon + B.rayon);
}

Index copier_noeud_FMM(const Noeud *noeud, const Index parent, StructureFMM &fmm);

void construire_interactions_FMM(const Index id_A, const Index id_B, StructureFMM &fmm);

void construire_blocs_denses_FMM(StructureFMM &fmm, const Maillage &maillage, const std::vector<QuadratureSegment> &quadratures);

StructureFMM construire_structure_FMM(const Arbre &arbre_clusters, const Real k, const Index ordre, const Maillage &maillage, const std::vector<QuadratureSegment> &quadratures);

void produit_blocs_denses_FMM(const StructureFMM &fmm, const Vecteur &x, Vecteur &y);

void calculer_P2M(StructureFMM &fmm, const std::vector<QuadratureSegment> &quadratures, const Vecteur &x);

void calculer_M2M(StructureFMM &fmm);

void calculer_M2L(StructureFMM &fmm);

void calculer_L2L(StructureFMM &fmm);

void calculer_L2P(const StructureFMM &fmm, const std::vector<QuadratureSegment> &quadratures, Vecteur &y);

Vecteur produit_FMM_x(StructureFMM &fmm, const Maillage &maillage, const std::vector<QuadratureSegment> &quadratures, const Vecteur &x);

#endif