
#ifndef BEM_CONDITIONNERS_HPP_INCLUDED
#define BEM_CONDITIONNERS_HPP_INCLUDED

#include "../global/commun.hpp"
#include "../maths/commun.hpp"

#include "operator.hpp"
#include "hmat/arbre_blocs.hpp"

// Condtionneur Diagaonal

// Calcul de la diagonal de l'opérateur A

inline Vecteur diagonale_matrice_bem(const Maillage &maillage, const std::vector<QuadratureSegment> &quadrature_maillage)
{
    const Index N = maillage.size();
    Vecteur diagonaleLocale = Vecteur::Zero(N);

    /*
     * Chaque processus calcule une partie de la diagonale.
     */
    for (Index i = mpi_rank; i < N; i += mpi_size)
    {
        const Segment &S = maillage[i];
        const QuadratureSegment &qS = quadrature_maillage[i];

        Complex aii = integ_double_log(S);
        aii += integ_double(green_reguliere_cached, qS, qS);

        diagonaleLocale[i] = aii;
    }

    /*
     * Chaque processus récupère la diagonale complète.
     */
    Vecteur diagonale = Vecteur::Zero(N);

    MPI_Allreduce(diagonaleLocale.data(), diagonale.data(), static_cast<int>(N), mpi_complex_type(), MPI_SUM, MPI_COMM_WORLD);

    return diagonale;
}

inline Vecteur preconditionneur_diagonal(const Vecteur &diagonale)
{
    const Index N = diagonale.size();
    Vecteur S = Vecteur::Zero(N);

    for (Index i = 0; i < N; ++i)
    {
        S[i] = 1.0 / std::sqrt(diagonale[i]);
    }

    return S;
}

// Conditionneur LUT

using TripletSparse = Eigen::Triplet<Complex, Eigen::Index>;

using Coordonnees = Eigen::Vector2d;

using AccesMatrice = std::function<Complex(Eigen::Index, Eigen::Index)>;

using Permutation = std::vector<Eigen::Index>;

using SupportsSparse = std::vector<std::vector<Eigen::Index>>;

struct SupportsDistribues
{
    Index debut = 0;
    SupportsSparse supports;
};

struct FacteursSparseDistribues
{
    Index debut = 0;
    Index fin = 0;

    SparseMatrix L_locale;
    SparseMatrix UT_locale;
};

class KDTree2D
{
private:
    struct NoeudKD
    {
        Index indice;
        Index gauche;
        Index droite;

        double xmin, xmax;
        double ymin, ymax;
    };

    const std::vector<Eigen::Vector2d> &centres;

    std::vector<NoeudKD> noeuds;

    Index racine = std::numeric_limits<Index>::max();

    // --------------------------------------------------
    // Construction récursive
    // --------------------------------------------------

    Index construire(std::vector<Index> &indices, const Index debut, const Index fin, const int axe)
    {
        if (debut >= fin)
            return std::numeric_limits<Index>::max();

        const Index milieu = debut + (fin - debut) / 2;

        std::nth_element(indices.begin() + debut, indices.begin() + milieu, indices.begin() + fin, [&](const Index i, const Index j)
                         { return centres[i](axe) < centres[j](axe); });

        const Index indice = indices[milieu];

        // Réserver la position du nœud avant les appels récursifs.
        const Index id = noeuds.size();

        noeuds.push_back({indice, std::numeric_limits<Index>::max(), std::numeric_limits<Index>::max(), centres[indice](0), centres[indice](0), centres[indice](1), centres[indice](1)});

        const int axe_suivant = 1 - axe;

        const Index gauche = construire(indices, debut, milieu, axe_suivant);
        const Index droite = construire(indices, milieu + 1, fin, axe_suivant);

        noeuds[id].gauche = gauche;
        noeuds[id].droite = droite;

        // Mise à jour de la boîte englobante.
        auto &n = noeuds[id];

        for (const Index fils : {gauche, droite})
        {
            if (fils == std::numeric_limits<Index>::max())
                continue;

            const auto &f = noeuds[fils];

            n.xmin = std::min(n.xmin, f.xmin);
            n.xmax = std::max(n.xmax, f.xmax);

            n.ymin = std::min(n.ymin, f.ymin);
            n.ymax = std::max(n.ymax, f.ymax);
        }

        return id;
    }

    // --------------------------------------------------
    // Distance minimale entre un point et une boîte
    // --------------------------------------------------

    double distance_boite2(const NoeudKD &n, const double x, const double y) const
    {
        const double dx = std::max({n.xmin - x, 0.0, x - n.xmax});

        const double dy = std::max({n.ymin - y, 0.0, y - n.ymax});

        return dx * dx + dy * dy;
    }

    // --------------------------------------------------
    // Recherche dans un rayon
    // --------------------------------------------------

    void rechercher(const Index id, const double x, const double y, const double rayon2, std::vector<Index> &resultats) const
    {
        if (id == std::numeric_limits<Index>::max())
            return;

        const auto &n = noeuds[id];

        // La boîte ne peut contenir aucun point recherché.
        if (distance_boite2(n, x, y) > rayon2)
            return;

        // Tester le point associé au nœud.
        const auto &p = centres[n.indice];

        const double dx = p(0) - x;
        const double dy = p(1) - y;

        if (dx * dx + dy * dy <= rayon2)
            resultats.push_back(n.indice);

        rechercher(n.gauche, x, y, rayon2, resultats);
        rechercher(n.droite, x, y, rayon2, resultats);
    }

public:
    // --------------------------------------------------
    // Constructeur
    // --------------------------------------------------

    explicit KDTree2D(const std::vector<Eigen::Vector2d> &centres_) : centres(centres_)
    {
        const Index N = centres.size();

        noeuds.reserve(N);

        std::vector<Index> indices(N);
        std::iota(indices.begin(), indices.end(), Index(0));

        racine = construire(indices, 0, N, 0);
    }

    // --------------------------------------------------
    // Recherche publique
    // --------------------------------------------------

    void rechercher_dans_rayon(const double x, const double y, const double rayon, std::vector<Index> &resultats) const
    {
        resultats.clear();

        if (rayon < 0.0)
            return;

        rechercher(racine, x, y, rayon * rayon, resultats);
    }
};

std::pair<Permutation, Eigen::VectorXd> construire_ordre_et_echelles(const std::vector<Coordonnees> &centres);

SupportsDistribues construire_supports_sparse_KD(const std::vector<Coordonnees> &centres, const Permutation &permutation, const Eigen::VectorXd &echelles, const double rho);

FacteursSparseDistribues construire_facteurs_sparse(const AccesMatrice &acces_A, const Permutation &permutation, const SupportsDistribues &supports_distribues, const double tolerance_pivot = 1e-12);

FacteursSparseDistribues construire_preconditionneur_sparse(const Maillage &maillage, const std::vector<QuadratureSegment> &quadrature_maillage, const double rho = 2.0, const double tolerance_pivot = 1e-12);

Vecteur appliquer_L(const FacteursSparseDistribues &facteurs, const Vecteur &x);

Vecteur appliquer_UT(const FacteursSparseDistribues &facteurs, const Vecteur &x);

// BLOC

struct BlocDiagonalInverse
{
    std::vector<Index> indices;
    Matrice inverse;
};

using PreconditionneurBlocs = std::vector<BlocDiagonalInverse>;

void extraire_blocs_diagonaux_locaux(const Bloc *bloc, std::vector<std::vector<Index>> &clusters);

PreconditionneurBlocs construire_inverses_blocs_diagonaux(const ArbreBlocs &arbre_blocs, const Maillage &maillage, const std::vector<QuadratureSegment> &quadrature_maillage);

#endif
