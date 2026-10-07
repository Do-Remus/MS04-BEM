#include "../export.hpp"

void ecrire_solution(const std::vector<Point> &pointsSolution, const std::vector<Complex> &u_cached, const ErreursSolution &erreurs)
{
    const string filename = string("outputs/cercle_u") +
                            "_k" + std::to_string(k) +
                            "_R" + std::to_string(rayon) +
                            "_L" + std::to_string(L) +
                            "_hM" + std::to_string(pasMaillage) +
                            "_hS" + std::to_string(pasSolution) +
                            "_d" + std::to_string(delta) +
                            "_N" + std::to_string(idxTroncature) +
                            "_q" + std::to_string(ordre) +
                            ".txt";

    std::ofstream file(filename);

    if (!file.is_open())
    {
        std::cout << "ERROR: Le fichier "
                  << filename
                  << " n'a pas pu être ouvert"
                  << std::endl;

        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    for (unsigned int j = 0; j < pointsSolution.size(); j++)
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