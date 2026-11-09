// Algorithme GMRES pour resoudre un system pour A sym complexe de manière itérative sans caculer stocker A
#ifndef BEM_ACA_GMRES_HPP_INCLUDED
#define BEM_ACA_GMRES_HPP_INCLUDED

#include "../../global/commun.hpp"
#include "../../maths/commun.hpp"

#include "operator.hpp"
#include "../conditionners.hpp"

struct ResultatGMRES
{
    Vecteur solution;
    Eigen::Index iterations;
    double erreur;
    Eigen::ComputationInfo info;
};

// Version vanilla

class OperateurH_Eigen;

namespace Eigen::internal
{
    template <>
    struct traits<OperateurH_Eigen>
        : public traits<Eigen::SparseMatrix<Complex>>
    {
    };
}

class OperateurH_Eigen : public Eigen::EigenBase<OperateurH_Eigen>
{
public:
    using Scalar = Complex;
    using RealScalar = double;
    using StorageIndex = int;
    using Index = Eigen::Index;

    enum
    {
        ColsAtCompileTime = Eigen::Dynamic,
        MaxColsAtCompileTime = Eigen::Dynamic,
        IsRowMajor = false
    };

    OperateurH_Eigen(
        const ArbreBlocs &arbre_blocs,
        const Index N,
        const bool est_symetrique = false)
        : arbre_blocs_(arbre_blocs),
          N_(static_cast<Eigen::Index>(N)),
          est_symetrique_(est_symetrique)
    {
    }

    Index rows() const { return N_; }
    Index cols() const { return N_; }

    template <typename Rhs>
    Eigen::Product<
        OperateurH_Eigen,
        Rhs,
        Eigen::AliasFreeProduct>
    operator*(const Eigen::MatrixBase<Rhs> &x) const
    {
        return Eigen::Product<
            OperateurH_Eigen,
            Rhs,
            Eigen::AliasFreeProduct>(*this, x.derived());
    }

    Vecteur appliquer(const Vecteur &x) const
    {
        return appliquer_H(
            arbre_blocs_,
            x,
            est_symetrique_);
    }

private:
    const ArbreBlocs &arbre_blocs_;
    Eigen::Index N_;
    bool est_symetrique_;
};

namespace Eigen::internal
{
    template <typename Rhs>
    struct generic_product_impl<
        OperateurH_Eigen,
        Rhs,
        SparseShape,
        DenseShape,
        GemvProduct>
        : generic_product_impl_base<
              OperateurH_Eigen,
              Rhs,
              generic_product_impl<
                  OperateurH_Eigen,
                  Rhs,
                  SparseShape,
                  DenseShape,
                  GemvProduct>>
    {
        using Scalar =
            typename Product<OperateurH_Eigen, Rhs>::Scalar;

        template <typename Dest>
        static void scaleAndAddTo(
            Dest &dst,
            const OperateurH_Eigen &lhs,
            const Rhs &rhs,
            const Scalar &alpha)
        {
            const Vecteur x = rhs;
            const Vecteur y = lhs.appliquer(x);

            dst += alpha * y;
        }
    };
}

inline ResultatGMRES resoudre_GMRES_Id(const ArbreBlocs &arbre_blocs, const Vecteur &b, const double tolerance, const Index max_iterations, const Index restart, const bool est_symetrique = false)
{
    OperateurH_Eigen A(arbre_blocs, b.size(), est_symetrique);

    { // TEST
        Vecteur x = Vecteur::Random(b.size());

        Vecteur Ax = appliquer_H(arbre_blocs, x, true);

        double erreur_preconditionneur =
            (Ax - x).norm() /
            std::max(x.norm(), 1e-30);

        if (mpi_rank == 0)
            cout << "Erreur approximation identité : " << erreur_preconditionneur << endl;
    }

    Eigen::GMRES<
        OperateurH_Eigen,
        Eigen::IdentityPreconditioner>
        gmres;

    gmres.set_restart(restart);
    gmres.setMaxIterations(max_iterations);
    gmres.setTolerance(tolerance);

    gmres.compute(A);

    ResultatGMRES resultat;

    resultat.solution = gmres.solve(b);
    resultat.iterations = gmres.iterations();
    resultat.erreur = gmres.error();
    resultat.info = gmres.info();

    return resultat;
}

// Version préconditionnée

class OperateurPreconditionneSparseLU;

namespace Eigen::internal
{
    template <>
    struct traits<OperateurPreconditionneSparseLU>
        : public traits<Eigen::SparseMatrix<Complex>>
    {
    };
}

class OperateurPreconditionneSparseLU : public Eigen::EigenBase<OperateurPreconditionneSparseLU>
{
public:
    using Scalar = Complex;
    using RealScalar = double;
    using StorageIndex = int;
    using Index = Eigen::Index;

    enum
    {
        ColsAtCompileTime = Eigen::Dynamic,
        MaxColsAtCompileTime = Eigen::Dynamic,
        IsRowMajor = false
    };

    OperateurPreconditionneSparseLU(const ArbreBlocs &arbre_blocs, const FacteursSparseDistribues &facteurs, const bool est_symetrique = false) : arbre_blocs_(arbre_blocs), facteurs_(facteurs), est_symetrique_(est_symetrique) {}

    Index rows() const
    {
        return static_cast<Index>(facteurs_.L_locale.rows());
    }

    Index cols() const
    {
        return rows();
    }

    template <typename Rhs>
    Eigen::Product<OperateurPreconditionneSparseLU, Rhs, Eigen::AliasFreeProduct> operator*(const Eigen::MatrixBase<Rhs> &x) const
    {
        return Eigen::Product<OperateurPreconditionneSparseLU, Rhs, Eigen::AliasFreeProduct>(*this, x.derived());
    }

    // Calcule UT * H * L * x.
    Vecteur appliquer(const Vecteur &x) const
    {
        const Vecteur Hx = appliquer_H(arbre_blocs_, x, est_symetrique_);

        const Vecteur UTHx = appliquer_UT(facteurs_, Hx);

        return appliquer_L(facteurs_, UTHx);
    }

private:
    const ArbreBlocs &arbre_blocs_;
    const FacteursSparseDistribues &facteurs_;
    bool est_symetrique_;
};

namespace Eigen::internal
{
    template <typename Rhs>
    struct generic_product_impl<OperateurPreconditionneSparseLU, Rhs, SparseShape, DenseShape, GemvProduct> : generic_product_impl_base<OperateurPreconditionneSparseLU, Rhs, generic_product_impl<OperateurPreconditionneSparseLU, Rhs>>
    {
        using Scalar = typename Product<OperateurPreconditionneSparseLU, Rhs>::Scalar;

        template <typename Dest>
        static void scaleAndAddTo(Dest &dst, const OperateurPreconditionneSparseLU &lhs, const Rhs &rhs, const Scalar &alpha)
        {
            const Vecteur x = rhs;
            const Vecteur y = lhs.appliquer(x);

            dst += alpha * y;
        }
    };
}

inline ResultatGMRES resoudre_GMRES_preconditionne(const ArbreBlocs &arbre_blocs, const FacteursSparseDistribues &facteurs, const Vecteur &b, const double tolerance, const Eigen::Index max_iterations, const bool est_symetrique = false)
{
    // Opérateur L * UT * H.
    OperateurPreconditionneSparseLU A(arbre_blocs, facteurs, est_symetrique);

    { // TEST : approximation de l'identité
        const Vecteur x = Vecteur::Random(b.size());
        const Vecteur Hx = appliquer_H(arbre_blocs, x, est_symetrique);
        const Vecteur UTHx = appliquer_UT(facteurs, Hx);
        const Vecteur LUTHx = appliquer_L(facteurs, UTHx);
        const double erreur = (LUTHx - x).norm() / std::max(x.norm(), 1e-30);

        if (mpi_rank == 0)
            cout << "Erreur approximation identite : "
                 << erreur << endl;
    }

    // Second membre préconditionné : L * UT * b.
    const Vecteur UTb = appliquer_UT(facteurs, b);
    const Vecteur b_preconditionne = appliquer_L(facteurs, UTb);

    // Construction du solveur GMRES.
    Eigen::GMRES<OperateurPreconditionneSparseLU, Eigen::IdentityPreconditioner> gmres;

    gmres.set_restart(40);
    gmres.setMaxIterations(max_iterations);
    gmres.setTolerance(tolerance);

    gmres.compute(A);

    // Résolution directe pour p :
    // (L * UT * H) p = L * UT * b.
    ResultatGMRES resultat;

    resultat.solution = gmres.solve(b_preconditionne);
    resultat.iterations = gmres.iterations();
    resultat.erreur = gmres.error();
    resultat.info = gmres.info();

    return resultat;
}

// preconditionnement par bloc

inline Vecteur appliquer_preconditionneur_blocs(const PreconditionneurBlocs &preconditionneur, const Vecteur &x)
{
    // Contribution calculée par le processus courant.
    Vecteur y_local = Vecteur::Zero(x.size());

    for (const auto &bloc : preconditionneur)
    {
        const Eigen::Index n = static_cast<Eigen::Index>(bloc.indices.size());

        Vecteur x_bloc(n);

        for (Eigen::Index i = 0; i < n; ++i)
        {
            x_bloc(i) = x(static_cast<Eigen::Index>(bloc.indices[i]));
        }

        const Vecteur y_bloc = bloc.inverse * x_bloc;

        for (Eigen::Index i = 0; i < n; ++i)
        {
            y_local(static_cast<Eigen::Index>(bloc.indices[i])) = y_bloc(i);
        }
    }

    // Assemblage des contributions de tous les processus.
    Vecteur y = Vecteur::Zero(x.size());

    MPI_Allreduce(y_local.data(), y.data(), static_cast<int>(y.size()), MPI_CXX_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);

    return y;
}

class OperateurPreconditionneBlocs;

namespace Eigen::internal
{
    template <>
    struct traits<OperateurPreconditionneBlocs> : public traits<Eigen::SparseMatrix<Complex>>
    {
    };
}

class OperateurPreconditionneBlocs : public Eigen::EigenBase<OperateurPreconditionneBlocs>
{
public:
    using Scalar = Complex;
    using RealScalar = double;
    using StorageIndex = int;
    using Index = Eigen::Index;

    enum
    {
        ColsAtCompileTime = Eigen::Dynamic,
        MaxColsAtCompileTime = Eigen::Dynamic,
        IsRowMajor = false
    };

    OperateurPreconditionneBlocs(const ArbreBlocs &arbre_blocs, const PreconditionneurBlocs &preconditionneur, const Index N, const bool est_symetrique = false) : arbre_blocs_(arbre_blocs), preconditionneur_(preconditionneur), N_(N), est_symetrique_(est_symetrique) {}

    Index rows() const { return N_; }
    Index cols() const { return N_; }

    template <typename Rhs>
    Eigen::Product<OperateurPreconditionneBlocs, Rhs, Eigen::AliasFreeProduct> operator*(const Eigen::MatrixBase<Rhs> &x) const
    {
        return Eigen::Product<OperateurPreconditionneBlocs, Rhs, Eigen::AliasFreeProduct>(*this, x.derived());
    }

    // Calcule M^-1 * H * x.
    Vecteur appliquer(const Vecteur &x) const
    {
        const Vecteur Ax = appliquer_H(arbre_blocs_, x, est_symetrique_);

        return appliquer_preconditionneur_blocs(preconditionneur_, Ax);
    }

private:
    const ArbreBlocs &arbre_blocs_;
    const PreconditionneurBlocs &preconditionneur_;
    Index N_;
    bool est_symetrique_;
};

namespace Eigen::internal
{
    template <typename Rhs>
    struct generic_product_impl<OperateurPreconditionneBlocs, Rhs, SparseShape, DenseShape, GemvProduct> : generic_product_impl_base<OperateurPreconditionneBlocs, Rhs, generic_product_impl<OperateurPreconditionneBlocs, Rhs, SparseShape, DenseShape, GemvProduct>>
    {
        using Scalar = typename Product<OperateurPreconditionneBlocs, Rhs>::Scalar;

        template <typename Dest>
        static void scaleAndAddTo(Dest &dst, const OperateurPreconditionneBlocs &lhs, const Rhs &rhs, const Scalar &alpha)
        {
            const Vecteur x = rhs;
            const Vecteur y = lhs.appliquer(x);

            dst += alpha * y;
        }
    };
}

inline ResultatGMRES resoudre_GMRES_BLOC(const ArbreBlocs &arbre_blocs, const PreconditionneurBlocs &preconditionneur, const Vecteur &b, const double tolerance, const Eigen::Index max_iterations, const Eigen::Index restart = 40, const bool est_symetrique = false)
{
    // Opérateur M^-1 * H.
    OperateurPreconditionneBlocs A(arbre_blocs, preconditionneur, b.size(), est_symetrique);

    // Test de l'approximation de l'identité.
    {
        Vecteur x = Vecteur::Zero(b.size());

        if (mpi_rank == 0)
        {
            x = Vecteur::Random(b.size());
        }

        MPI_Bcast(x.data(), static_cast<int>(x.size()), MPI_CXX_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);

        const Vecteur Hx = appliquer_H(arbre_blocs, x, est_symetrique);
        const Vecteur M_inv_Hx = appliquer_preconditionneur_blocs(preconditionneur, Hx);
        const double erreur = (M_inv_Hx - x).norm() / std::max(x.norm(), 1e-30);

        if (mpi_rank == 0)
        {
            cout << "Erreur approximation identite (bloc) : "
                 << erreur << endl;
        }
    }

    // Second membre préconditionné : M^-1 * b.
    const Vecteur b_preconditionne = appliquer_preconditionneur_blocs(preconditionneur, b);

    // Construction de GMRES.
    Eigen::GMRES<OperateurPreconditionneBlocs, Eigen::IdentityPreconditioner> gmres;

    gmres.set_restart(restart);
    gmres.setMaxIterations(max_iterations);
    gmres.setTolerance(tolerance);

    gmres.compute(A);

    // Résolution : M^-1 * H * p = M^-1 * b.
    ResultatGMRES resultat;

    resultat.solution = gmres.solve(b_preconditionne);
    resultat.iterations = gmres.iterations();
    resultat.erreur = gmres.error();
    resultat.info = gmres.info();

    return resultat;
}

#endif