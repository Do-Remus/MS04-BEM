#include "../utils.hpp"

Real erreur_relative(const Complex &u, const Complex &reference, Real eps)
{
    const Real denom = std::max(std::abs(reference), eps);
    return std::abs(u - reference) / denom;
}

string trim(const string &str)
{
    const size_t first = str.find_first_not_of(" \t\r\n");
    if (first == string::npos)
        return "";

    const size_t last = str.find_last_not_of(" \t\r\n");
    return str.substr(first, last - first + 1);
}
