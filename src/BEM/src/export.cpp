#include "../export.hpp"

void ecrire_solution(const std::vector<Point> &pointsSolution, const Vecteur &u_cached, const ErreursSolution &erreurs, const std::string &fichier)
{
    std::ofstream file(fichier);

    if (!file.is_open())
    {
        std::cout << "ERROR: Le fichier "
                  << fichier
                  << " n'a pas pu être ouvert"
                  << std::endl;

        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    for (Index j = 0; j < pointsSolution.size(); j++)
    {
        const Complex &uExact = erreurs.u_analytique[j];
        const Complex &uCached = u_cached[j];

        file << pointsSolution[j].x << " "
             << pointsSolution[j].y << " "
             << uExact.real() << " "
             << uExact.imag() << " "
             << uCached.real() << " "
             << uCached.imag() << " "
             << erreurs.erreurs[j]
             << "\n";
    }

    file.close();
}
