#ifndef CONFIG_HPP_INCLUDED
#define CONFIG_HPP_INCLUDED

#include "external.hpp"
#include "constantes.hpp"

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

// -- Precalculé --

// Cache pour green
extern Complex *green_cache;
extern MPI_Win green_cache_win;
extern Real green_cache_step_fine;
extern Real green_cache_step_fine_inv;
extern Real green_cache_step_large;
extern Real green_cache_step_large_inv;
extern Real green_cache_cutoff;
extern Real delta;
extern std::size_t green_cache_size;

// -- Fonctions --
string trim(const string &str);
void get_config(const string &filename);

#endif