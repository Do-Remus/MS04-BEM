#include "../basics.hpp"

Real erreur_relative(const Complex &u, const Complex &reference, Real eps)
{
    const Real denom = std::max(std::abs(reference), eps);
    return std::abs(u - reference) / denom;
}
