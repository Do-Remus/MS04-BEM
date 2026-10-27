#ifndef BEM_ACA_OPERATOR_HPP_INCLUDED
#define BEM_ACA_OPERATOR_HPP_INCLUDED

#include "../../global/commun.hpp"
#include "../../maths/commun.hpp"

#include "../hmat/commun.hpp"

void appliquer_bloc(const Bloc &bloc, const Vecteur &x, Vecteur &y);

Vecteur appliquer_H(const ArbreBlocs &arbre_blocs, const Vecteur &x);

// Adaptateur matriciel : aucun stockage de la matrice H.
class OperateurH;

namespace Eigen
{
    namespace internal
    {

        template <>
        struct traits<OperateurH>
            : public Eigen::internal::traits<Eigen::SparseMatrix<::Complex>>
        {
        };

    } // namespace internal
} // namespace Eigen

class OperateurH : public Eigen::EigenBase<OperateurH>
{
public:
    using Scalar = ::Complex;
    using RealScalar = double;
    using StorageIndex = int;

    enum
    {
        ColsAtCompileTime = Eigen::Dynamic,
        MaxColsAtCompileTime = Eigen::Dynamic,
        IsRowMajor = false
    };

    OperateurH(const ArbreBlocs &arbre, Eigen::Index n)
        : arbre_(&arbre), n_(n) {}

    Eigen::Index rows() const { return n_; }
    Eigen::Index cols() const { return n_; }

    template <typename Rhs>
    Eigen::Product<OperateurH, Rhs, Eigen::AliasFreeProduct>
    operator*(const Eigen::MatrixBase<Rhs> &x) const
    {
        return Eigen::Product<OperateurH, Rhs, Eigen::AliasFreeProduct>(
            *this, x.derived());
    }

    // Produit H*x calculé directement par l'arbre de blocs.
    Vecteur appliquer(const Vecteur &x) const
    {
        return appliquer_H(*arbre_, x);
    }

private:
    const ArbreBlocs *arbre_;
    Eigen::Index n_;
};

namespace Eigen
{
    namespace internal
    {

        template <typename Rhs>
        struct generic_product_impl<
            OperateurH, Rhs, SparseShape, DenseShape, GemvProduct>
            : generic_product_impl_base<
                  OperateurH, Rhs,
                  generic_product_impl<OperateurH, Rhs>>
        {

            using Scalar = typename Product<OperateurH, Rhs>::Scalar;

            template <typename Dest>
            static void scaleAndAddTo(
                Dest &dst,
                const OperateurH &lhs,
                const Rhs &rhs,
                const Scalar &alpha)
            {

                Vecteur x = rhs;
                Vecteur y = lhs.appliquer(x);

                dst += alpha * y;
            }
        };

    } // namespace internal
} // namespace Eigen

#endif