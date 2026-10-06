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
extern unsigned int idxTroncature;
extern unsigned int ordre;
extern unsigned int maxIterGradConj;
extern Real tolGradConj;

// -- Fonctions --

void get_config(const string &filename);

#endif