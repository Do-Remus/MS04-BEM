#include "config/config.hpp"

/* Définition des varibles globales par default */

Real k = 10; // frequence des ondes : il y a un nombre au plus denombrable de frequences pour lesquelles le pb n'est pas bien pose
Real rayon = 1.0;
Real L = 4.0;
Real pasMaillage = 0.001;
Real pasSolution = 0.1;
unsigned int idxTroncature = 25;
unsigned int ordre = 4;
unsigned int maxIterGradConj = 1000;
Real tolGradConj = 1e-8;

// -- Précalculé --

// Cache pour green
double lambda = 2.0 * pi / k;
Complex *green_cache = nullptr;
MPI_Win green_cache_win = MPI_WIN_NULL;
Real green_cache_step_fine = 0.001;
Real green_cache_step_fine_inv = 1 / green_cache_step_fine;
Real green_cache_step_large = 0.05 * lambda;
Real green_cache_step_large_inv = 1 / green_cache_step_large;
Real green_cache_cutoff = 0.5 * lambda;
Real delta = green_cache_step_fine;
std::size_t green_cache_size = 0;

/* Déclaration des fonctions globales */

void get_config(const string &filename)
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    ifstream f(filename);
    if (!f.is_open())
    {
        if (rank == 0)
        {
            cout << "ERROR: Le fichier " << filename << " n'a pas pu être ouvert" << endl;
        }
        exit(-1);
    }

    string ligne{""};
    while (getline(f, ligne))
    {
        stringstream streamLigne(ligne);
        string nom{""};
        string valeur{""};

        if (!getline(streamLigne, nom, '='))
        {
            continue;
        }

        if (!getline(streamLigne, valeur))
        {
            continue;
        }

        nom = trim(nom);
        valeur = trim(valeur);

        if (nom == "k")
        {
            k = atof(valeur.c_str());
            if (!(k > 0.000000001))
            {
                if (rank == 0)
                {
                    cout << "La valeur de k doit-être strictement supérieur à 10-9.";
                }
                exit(-1);
            }
        }

        if (nom == "rayon")
        {
            rayon = atof(valeur.c_str());
            if (!(rayon > 0))
            {
                if (rank == 0)
                {
                    cout << "La valeur de rayon doit-être strictement positive.";
                }
                exit(-1);
            }
        }

        if (nom == "L")
        {
            L = atof(valeur.c_str());
            if (!(L > 0))
            {
                if (rank == 0)
                {
                    cout << "La valeur de L doit-être strictement positive.";
                }
                exit(-1);
            }
        }

        if (nom == "pasMaillage")
        {
            pasMaillage = atof(valeur.c_str());
            if (!(pasMaillage > 0))
            {
                if (rank == 0)
                {
                    cout << "La valeur de pasMaillage doit-être strictement positive.";
                }
                exit(-1);
            }
        }

        if (nom == "pasSolution")
        {
            pasSolution = atof(valeur.c_str());
            if (!(pasSolution > 0))
            {
                if (rank == 0)
                {
                    cout << "La valeur de pasSolution doit-être strictement positive.";
                }
                exit(-1);
            }
        }

        if (nom == "idxTroncature")
        {
            idxTroncature = atoi(valeur.c_str());
            if (!(idxTroncature > 0))
            {
                if (rank == 0)
                {
                    cout << "La valeur de nbObstacles doit-être strictement positive.";
                }
                exit(-1);
            }
        }

        if (nom == "ordre")
        {
            ordre = atoi(valeur.c_str());
            if (!(ordre > 0))
            {
                if (rank == 0)
                {
                    cout << "La valeur de ordre doit-être strictement positive.";
                }
                exit(-1);
            }
        }

        if (nom == "maxIterGradConj")
        {
            maxIterGradConj = atoi(valeur.c_str());
            if (!(maxIterGradConj > 0))
            {
                if (rank == 0)
                {
                    cout << "La valeur de maxIterGradConj doit-être strictement positive.";
                }
                exit(-1);
            }
        }

        if (nom == "tolGradConj")
        {
            tolGradConj = atof(valeur.c_str());
            if (!(tolGradConj > 0))
            {
                if (rank == 0)
                {
                    cout << "La valeur de tolGradConj doit-être strictement positive.";
                }
                exit(-1);
            }
        }
    }

    return;
}

string trim(const string &str)
{
    const auto first = str.find_first_not_of(" \t\r\n");
    if (first == string::npos)
        return "";

    const auto last = str.find_last_not_of(" \t\r\n");
    return str.substr(first, last - first + 1);
}
