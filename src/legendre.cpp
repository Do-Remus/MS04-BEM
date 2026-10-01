#include "headers/legendre.hpp"

const std::vector<Real> &get_legendre_roots(unsigned n_ordre)
{
    static std::unordered_map<unsigned, std::vector<Real>> cache;

    auto it = cache.find(n_ordre);

    if (it == cache.end())
    {
        auto roots = boost::math::legendre_p_zeros<Real>(n_ordre);
        it = cache.emplace(n_ordre, std::move(roots)).first;
    }

    return it->second;
}

const LegendreData &get_legendre_data(unsigned n_ordre)
{
    static std::unordered_map<unsigned, LegendreData> cache;

    auto it = cache.find(n_ordre);

    if (it == cache.end())
    {
        LegendreData data;

        data.roots = boost::math::legendre_p_zeros<Real>(n_ordre);

        data.weights.reserve(data.roots.size());

        for (Real x : data.roots)
        {
            Real dp = boost::math::legendre_p_prime(n_ordre, x);

            Real w = 2.0 / ((1.0 - x * x) * dp * dp);

            data.weights.push_back(w);
        }

        it = cache.emplace(n_ordre, std::move(data)).first;
    }

    return it->second;
}

const std::vector<QuadratureSegment> get_quadrature_maillage(const Maillage &maillage, const LegendreData &data)
{
    std::vector<QuadratureSegment> quadrature(maillage.size());

    for (std::size_t i = 0; i < maillage.size(); ++i)
    {
        const Segment &S = maillage[i];
        auto &q = quadrature[i];

        // Nombre de points de quadrature.
        // Pour Gauss-Legendre d'ordre pair : 2 * roots.size()
        // Pour un ordre impair contenant r = 0 : 2*n - 1
        q.reserve(2 * data.roots.size());

        for (std::size_t k = 0; k < data.roots.size(); ++k)
        {
            const Real r = data.roots[k];
            const Real w = data.weights[k];

            const Point correction = (0.5 * r) * S.vecteur_norm;
            const Real weight = 0.5 * S.norm * w;

            q.push_back({S.milieu + correction, weight});

            if (r > 0.0)
            {
                q.push_back({S.milieu - correction, weight});
            }
        }
    }

    return quadrature;
}