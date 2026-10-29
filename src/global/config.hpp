#ifndef CONFIG_HPP_INCLUDED
#define CONFIG_HPP_INCLUDED

#include "external.hpp"
#include "constantes.hpp"
#include "utils.hpp"

// -- Config --
extern Real k; // frequence des ondes : il y a un nombre au plus denombrable de frequences pour lesquelles le pb n'est pas bien pose
extern Real rayon;
extern Real L;
extern Real pasMaillage;
extern Real pasSolution;
extern Index N_leaves;
extern Real parametre_admissibilite;
extern Real facteur_admissibilite;
extern Index rang_max;
extern Real toleranceAca;
extern Real facteur_step_green;
extern unsigned int idxTroncature;
extern unsigned int ordre;
extern unsigned int maxIterGradConj;
extern Real tolGradConj;
extern std::string nom_maillage;
extern bool affichageErreur;

// -- Fonctions --

void get_config(const string &filename);

#endif