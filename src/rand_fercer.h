/************************************************************************************************************
*                                                                                                           *
* CENTRO DE INVESTIGACION EN MATEMATICAS                                                                    *
* DOCTORADO EN CIENCIAS DE LA COMPUTACION                                                                   *
* FERNANDO CERVANTES SANCHEZ                                                                                *
*                                                                                                           *
* FILE NAME: rand_fercer.h                                                                                  *
*                                                                                                           *
* DEVELOPMENT HISTORY:                                                                                      *
* Date        Author        Change Id    Release    Description Of Change                                   *
* Nov/2016    Fernando C.   0            1.0        Creation                                                *
*                                                                                                           *
************************************************************************************************************/


#ifndef RAND_FERCER_H_INCLUDED
#define RAND_FERCER_H_INCLUDED

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>


/* STAUS is the structure that stores the information necesary for the random number generator Hybrid Taus */
typedef struct {
	unsigned int z1, z2, z3, lcg_seed; /* Seed required by the Hybrid Taustep method to generate random numbers */
} STAUS;


/***********************************************************************************************************/
STAUS* initSeed(unsigned int initial_seed);


/***********************************************************************************************************/
unsigned int lcgR(unsigned int *my_seed);
unsigned int tausStep(unsigned int *z, const int S1, const int S2, const int S3, const unsigned int M);
double HybTaus(const double par1, const double par2, STAUS *my_seed);

double anorm_est(STAUS *my_seed);
double anorm(const double par1, const double par2, STAUS *my_seed);



#endif //RAND_FERCER_H_INCLUDED
