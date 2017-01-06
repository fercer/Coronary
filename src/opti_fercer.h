/*
*******************************************************************************************************************
*    CENTRO DE INVESTIGACION EN MATEMATICAS
*    DOCTORADO EN CIENCIAS DE LA COMPUTACION
*
*    FERNANDO CERVANTES SANCHEZ
*    AGO - 2016
*******************************************************************************************************************
*/


#ifndef OPTI_FERCER_H_INCLUDED
#define OPTI_FERCER_H_INCLUDED

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include <math.h>


#include "rand_fercer.h"


/** STAUS:	Define una estructura que almacena las semillas requeridas por el generador de numeros pseudo-aleatorios Hybrid Taus    **/
typedef struct STAUS{
    unsigned int z1, z2, z3, lcg_seed;
} STAUS;





/** NODO_ESTRUCT:	Almacena la estructura de la cadena utilizada por el MIMIC y CHOWLIU y las probabilidades condicionales dado el valor del padre    **/
typedef struct NODO_ESTRUCT{
    int varID;
    double prob_cond_0, prob_cond_1, covar;
    // prob_cond_0: Probabilidad de que varID sea 0 dado que el nodo padre sea 0,
    // prob_cond_1: Probabilidad de que varID sea 0 dado que el nodo padre sea 1,
    // EN EL CASO CONTINUO:
    // prob_cond_0: Media de la variable varID,
    // prob_cond_1: Varianza de la variable varID,
    // covar: Covarianza entre la variable varID  y su padre.

    int n_hijos;
    struct NODO_ESTRUCT** nodos_hijos;

    double entropia; // La entropia marginal
    double info;    // La informacion mutua entre el nodo actual y su padre.
} NODO_ESTRUCT;



/** OPTIMET:             Puntero a funcion del metodo de optimizacion    **/
typedef INDIV* (*OPTIMET) (FUNPARS* mis_optipars, FUNPTR mi_fun, FUNPARS *mis_pars, struct STAUS *mi_semilla);




struct STAUS* ini_semilla(unsigned int semilla_i);

unsigned int lcg_r(unsigned int *mi_semilla);
unsigned int tausStep(unsigned int *z, const int S1, const int S2, const int S3, const unsigned int M);
double HybTaus(const double par1, const double par2, struct STAUS *mi_semilla);
double anorm_est(struct STAUS *mi_semilla);
double anorm(const double par1, const double par2, STAUS *mi_semilla);

void shuffle(int *indices, const int n0, const int n1, STAUS *mi_semilla);

void genIndiv(INDIV *indiv, struct NODO_ESTRUCT* nodo_padre, struct NODO_ESTRUCT* nodo_hijo, struct STAUS *mi_semilla);
void genPob(INDIV *pob, const int n_pob, struct NODO_ESTRUCT* estructura, FUNPTR mi_fun, FUNPARS *mis_pars, struct STAUS *mi_semilla);

void genIndivCont(INDIV *indiv, NODO_ESTRUCT* nodo_padre, NODO_ESTRUCT* nodo_hijo, FUNPARS *mis_pars, STAUS *mi_semilla);
void genPobCont(INDIV *pob, const int n_pob, NODO_ESTRUCT* estructura, FUNPTR mi_fun, FUNPARS *mis_pars, STAUS *mi_semilla);

void liberarNodos( struct NODO_ESTRUCT *estructura );

void calcFrecs( INDIV *pob, const int n_pob, const int n_vars, unsigned int **frecuencias );
void calcInfoMutua(double **info_mat, const int n_pob, const int n_vars, unsigned int **frecuencias );

void calcCovars(INDIV *pob, const int n_pob, const int n_vars, double **covars_mat, double *medias);
void calcInfoMutuaCont(double **info_mat, const int n_vars, double **covars_mat);

struct NODO_ESTRUCT *estructMIMIC(const int n_pob, const int n_vars, unsigned int **frecuencias, double **info_mat);
struct NODO_ESTRUCT *estructMIMICCont(const int n_vars, double **covars_mat, double *medias, double **info_mat);

struct NODO_ESTRUCT *estructCHOWLIU(const int n_pob, const int n_vars, unsigned int **frecuencias, double **info_mat);
struct NODO_ESTRUCT *estructCHOWLIUCont(const int n_vars, double **covars_mat, double *medias, double **info_mat);

INDIV *optiMIMIC(FUNPARS* mis_optipars, FUNPTR mi_fun, FUNPARS *mis_pars, struct STAUS *mi_semilla);
INDIV *optiMIMICCont(FUNPARS* mis_optipars, FUNPTR mi_fun, FUNPARS *mis_pars, struct STAUS *mi_semilla);

INDIV *optiCHOWLIU(FUNPARS* mis_optipars, FUNPTR mi_fun, FUNPARS *mis_pars, struct STAUS *mi_semilla);
INDIV *optiCHOWLIUCont( FUNPARS *mis_optipars, FUNPTR mi_fun, FUNPARS *mis_pars, STAUS *mi_semilla);

void selecPob(INDIV* sel, const INDIV* pob, const int n_pob, const int n_vars, struct STAUS *mi_semilla);
void crossOver(INDIV* pob, INDIV* sel, const int n_pob, struct STAUS *mi_semilla, int *indices, const double p_cruza, const double p_mutacion, FUNPTR mi_fun, FUNPARS *mis_pars);
INDIV *optiGA(FUNPARS* mis_optipars, FUNPTR mi_fun, FUNPARS *mis_pars, struct STAUS *mi_semilla);

double entropiaModelo(struct NODO_ESTRUCT *raiz);
double infomutuaModelo(struct NODO_ESTRUCT *raiz);

void mostrarModelo( struct NODO_ESTRUCT *nodo_padre, struct NODO_ESTRUCT *nodo_hijo );

#endif //OPTI_FERCER_H_INCLUDED
