#include "headers/legendre.hpp"

const std::vector<double> &get_legendre_roots(unsigned n_ordre)
{
    static std::unordered_map<unsigned, std::vector<double>> cache;

    auto it = cache.find(n_ordre);

    if (it == cache.end())
    {
        auto roots = boost::math::legendre_p_zeros<double>(n_ordre);
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

        data.roots = boost::math::legendre_p_zeros<double>(n_ordre);

        data.weights.reserve(data.roots.size());

        for (double x : data.roots)
        {
            double dp = boost::math::legendre_p_prime(n_ordre, x);

            double w = 2.0 / ((1.0 - x * x) * dp * dp);

            data.weights.push_back(w);
        }

        it = cache.emplace(n_ordre, std::move(data)).first;
    }

    return it->second;
}
