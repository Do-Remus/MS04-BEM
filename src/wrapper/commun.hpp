#ifndef WRAPPER_COMMUN_HPP_INCLUDED
#define WRAPPER_COMMUN_HPP_INCLUDED

#include "../global/commun.hpp"
#include "../BEM/commun.hpp"
#include "../maths/commun.hpp"
#include "../affichage/commun.hpp"

/* Wrapper pour avoir un main plus propre avec une logique plus claire */

template <typename F>
decltype(auto) run(const std::string &titre, Timers &timers, const std::string &nomTimer, F &&fonction)
{
    if (mpi_rank == 0)
        std::cout << "\n=== " << titre << " ===" << std::endl;

    timers.start(nomTimer);

    if constexpr (std::is_void_v<std::invoke_result_t<F &>>)
    {
        std::forward<F>(fonction)();

        timers.end(nomTimer);

        if (mpi_rank == 0)
            timers.print(nomTimer);
    }
    else
    {
        decltype(auto) resultat =
            std::forward<F>(fonction)();

        timers.end(nomTimer);

        if (mpi_rank == 0)
            timers.print(nomTimer);

        return resultat;
    }
}

template <typename F, typename A>
decltype(auto) run(const std::string &titre, Timers &timers, const std::string &nomTimer, F &&fonction, A &&affichage)
{
    if (mpi_rank == 0)
        std::cout << "\n=== " << titre << " ===" << std::endl;

    timers.start(nomTimer);

    if constexpr (std::is_void_v<std::invoke_result_t<F &>>)
    {
        std::forward<F>(fonction)();

        timers.end(nomTimer);

        if (mpi_rank == 0)
        {
            timers.print(nomTimer);
            std::forward<A>(affichage)();
        }
    }
    else
    {
        decltype(auto) resultat =
            std::forward<F>(fonction)();

        timers.end(nomTimer);

        if (mpi_rank == 0)
        {
            timers.print(nomTimer);
            std::forward<A>(affichage)(resultat);
        }

        return resultat;
    }
}

#endif