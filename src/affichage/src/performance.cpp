#include "../performance.hpp"

void Timers::start(const std::string &name)
{
    MPI_Barrier(MPI_COMM_WORLD);

    auto it = timers.find(name);

    // Premier appel à start() pour ce nom
    if (it == timers.end())
    {
        timers[name] = Timer{};
        names.push_back(name);
        it = timers.find(name);
    }

    Timer &timer = it->second;

    if (timer.running)
    {
        throw std::runtime_error("Le timer '" + name + "' est deja en cours.");
    }

    timer.start = MPI_Wtime();
    timer.running = true;
}

void Timers::end(const std::string &name)
{
    MPI_Barrier(MPI_COMM_WORLD);

    auto it = timers.find(name);

    if (it == timers.end())
    {
        throw std::runtime_error(
            "Le timer '" + name + "' n'existe pas.");
    }

    Timer &timer = it->second;

    if (!timer.running)
    {
        throw std::runtime_error(
            "Le timer '" + name + "' n'est pas en cours.");
    }

    timer.elapsed = MPI_Wtime() - timer.start;
    timer.running = false;
}

double Timers::get(const std::string &name) const
{
    auto it = timers.find(name);

    if (it == timers.end())
    {
        throw std::runtime_error("Le timer '" + name + "' n'existe pas.");
    }

    return it->second.elapsed;
}

const std::vector<std::string> &Timers::get_names() const
{
    return names;
}

void Timers::print(const std::string &name) const
{
    std::cout << "    "
              << std::left
              << std::setw(LONGEUR_TEXT_AFFICHAGE)
              << ("Temps (" + name + ")")
              << " : "
              << get(name)
              << " s"
              << std::endl;
}

void Timers::print_all(const std::string &title) const
{
    std::cout << "\n"
              << title << std::endl;

    for (const std::string &name : names)
    {
        std::cout << "    "
                  << std::left
                  << std::setw(LONGEUR_TEXT_AFFICHAGE)
                  << name
                  << " : "
                  << get(name)
                  << " s"
                  << std::endl;
    }
}