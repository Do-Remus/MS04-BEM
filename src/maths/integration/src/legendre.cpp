#include "../legendre.hpp"

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

QuadratureSIMD construire_quadrature_SIMD(const Maillage &maillage, const std::vector<QuadratureSegment> &quadrature_maillage)
{
    QuadratureSIMD q;

    q.N = maillage.size();
    q.nq = quadrature_maillage[0].size();

    const std::size_t N = q.N;
    const std::size_t nq = q.nq;

    /*
     * Allocation.
     */
    q.target_x.resize(N * nq);
    q.target_y.resize(N * nq);
    q.target_w.resize(N * nq);

    q.source_x.resize(N * nq);
    q.source_y.resize(N * nq);
    q.source_w.resize(N * nq);

    /*
     * Remplissage.
     */
    for (std::size_t i = 0; i < N; ++i)
    {
        for (std::size_t a = 0; a < nq; ++a)
        {
            const auto &qa = quadrature_maillage[i][a];

            /*
             * Layout cible :
             *
             * [i][a]
             */
            q.target_x[i * nq + a] = qa.point.x;
            q.target_y[i * nq + a] = qa.point.y;
            q.target_w[i * nq + a] = qa.weight;

            /*
             * Layout source transposé :
             *
             * [a][i]
             */
            q.source_x[a * N + i] = qa.point.x;
            q.source_y[a * N + i] = qa.point.y;
            q.source_w[a * N + i] = qa.weight;
        }
    }

    return q;
}
