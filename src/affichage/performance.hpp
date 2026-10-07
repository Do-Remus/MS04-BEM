#ifndef AFFICHAGE_PERFORMANCE_HPP_INCLUDED
#define AFFICHAGE_PERFORMANCE_HPP_INCLUDED

#include "../global/commun.hpp"

class Timers
{
private:
    struct Timer
    {
        double start = 0.0;
        double elapsed = 0.0;
        bool running = false;
    };

    std::unordered_map<std::string, Timer> timers;

    // Ordre de création des timers
    std::vector<std::string> names;

public:
    // Démarre un timer.
    void start(const std::string &name);

    // Arrête un timer.
    void end(const std::string &name);

    // Retourne le temps mesuré.
    double get(const std::string &name) const;

    // Retourne les noms des timers dans leur ordre de création.
    const std::vector<std::string> &get_names() const;

    // Affiche le timer voulu
    void print(const std::string &name) const;

    // Affiche tous les timers.
    void print_all(const std::string &title = "=== Temps ===") const;
};

#endif
