#include "config/config.hpp"

/* Définition des varibles globales par default */

double theta = 0;
bool effectuerTests = false;
bool effectuerLaSimulation = true;
double k = 2; // frequence des ondes : il y a un nombre au plus denombrable de frequences pour lesquelles le pb n'est pas bien pose
double h = 10;
double e = 2;
double pasMaillage = 0.1;
unsigned int nbObstacles = 1;
double pasIntegrale = 0.01;
unsigned int nbPasExport = 100;

/* Déclaration des fonctions globales */

void get_config(const string &filename)
{
    ifstream f(filename);
    if (!f.is_open())
    {
        cout << "ERROR: Le fichier " << filename << " n'a pas pu être ouvert" << endl;
        exit(-1);
    }

    string ligne{""};
    while (getline(f, ligne))
    {
        stringstream streamLigne(ligne);
        string nom{""};
        string valeur{""};

        if (!getline(streamLigne, nom, '=')) // Cas d'erreur de lecture (passage à la ligne suivante)
        {
            continue;
        }

        if (!getline(streamLigne, valeur)) // Cas d'erreur de lecture (passage à la ligne suivante)
        {
            continue;
        }

        if (nom == "theta")
        {
            cout << "Récupération de theta:" << endl;
            theta = atof(valeur.c_str());
            cout << "     theta = " << theta << endl;
        }

        if (nom == "k")
        {
            cout << "Récupération de k:" << endl;
            k = atof(valeur.c_str());
            if (!(k > 0.000000001))
            {
                cout << "La valeur de k doit-être strictement supérieur à 10-9.";
                exit(-1);
            }
            cout << "     k = " << k << endl;
        }

        if (nom == "h")
        {
            cout << "Récupération de h:" << endl;
            k = atof(valeur.c_str());
            if (!(k > 0))
            {
                cout << "La valeur de h doit-être strictement positive.";
                exit(-1);
            }
            cout << "     h = " << h << endl;
        }

        if (nom == "e")
        {
            cout << "Récupération de e:" << endl;
            k = atof(valeur.c_str());
            if (!(k > 0))
            {
                cout << "La valeur de e doit-être strictement positive.";
                exit(-1);
            }
            cout << "     e = " << e << endl;
        }

        if (nom == "pasMaillage")
        {
            cout << "Récupération de pasMaillage:" << endl;
            pasMaillage = atof(valeur.c_str());
            if (!(pasMaillage > 0))
            {
                cout << "La valeur de pasMaillage doit-être strictement positive.";
                exit(-1);
            }
            cout << "     pasMaillage = " << pasMaillage << endl;
        }

        if (nom == "pasIntegrale")
        {
            cout << "Récupération de pasIntegrale:" << endl;
            pasIntegrale = atof(valeur.c_str());
            if (!(pasIntegrale > 0))
            {
                cout << "La valeur de pasIntegrale doit-être strictement positive.";
                exit(-1);
            }
            cout << "     pasIntegrale = " << pasIntegrale << endl;
        }

        if (nom == "nbObstacles")
        {
            cout << "Récupération de nbObstacles:" << endl;
            pasIntegrale = atoi(valeur.c_str());
            if (!(nbObstacles > 0))
            {
                cout << "La valeur de nbObstacles doit-être strictement positive.";
                exit(-1);
            }
            cout << "     nbObstacles = " << nbObstacles << endl;
        }

        if (nom == "nbPasExport")
        {
            cout << "Récupération de nbPasExport:" << endl;
            nbPasExport = atoi(valeur.c_str());
            if (!(nbPasExport > 0))
            {
                cout << "La valeur de nbPasExport doit-être strictement positive.";
                exit(-1);
            }
            cout << "     nbPasExport = " << nbPasExport << endl;
        }

        if (nom == "effectuerTests")
        {
            cout << "Récupération de effectuerTests:" << endl;
            if (valeur == "true")
            {
                effectuerTests = true;
            }
            else if (valeur == "false")
            {
                effectuerTests = false;
            }
            else
            {
                cout << "La valeur de effectuerTests doit-être 'true' ou 'false'." << endl;
            }
            cout << "     effectuerTests = " << effectuerTests << endl;
        }

        if (nom == "effectuerLaSimulation")
        {
            cout << "Récupération de effectuerLaSimulation:" << endl;
            if (valeur == "true")
            {
                effectuerLaSimulation = true;
            }
            else if (valeur == "false")
            {
                effectuerLaSimulation = false;
            }
            else
            {
                cout << "La valeur de effectuerLaSimulation doit-être 'true' ou 'false'." << endl;
            }
            cout << "     effectuerLaSimulation = " << effectuerLaSimulation << endl;
        }
    }

    if (pasMaillage <= pasIntegrale)
    {
        cout << "Le pas du maillage doit-être inférieur au pas pour les intégrales." << endl;
        exit(-1);
    }

    return;
}

complex<double> p_theta_config(double x, double y)
{
    return (exp(I * (x * cos(theta) + y * sin(theta))));
}
