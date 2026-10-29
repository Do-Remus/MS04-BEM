#include "../utils.hpp"

string trim(const string &str)
{
    const size_t first = str.find_first_not_of(" \t\r\n");
    if (first == string::npos)
        return "";

    const size_t last = str.find_last_not_of(" \t\r\n");
    return str.substr(first, last - first + 1);
}

std::string formater_nombre(double valeur)
{
    std::ostringstream oss;
    oss << std::defaultfloat << std::setprecision(6) << valeur;
    std::string resultat = oss.str();

    // Uniformiser la notation scientifique : 1e-07 -> 1e-7
    const auto pos = resultat.find('e');
    if (pos != std::string::npos && pos + 1 < resultat.size())
    {
        std::size_t debut = pos + 1;

        if (resultat[debut] == '+' || resultat[debut] == '-')
            ++debut;

        while (debut + 1 < resultat.size() && resultat[debut] == '0')
            resultat.erase(debut, 1);
    }

    return resultat;
}
