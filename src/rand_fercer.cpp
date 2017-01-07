/************************************************************************************************************
*                                                                                                           *
* CENTRO DE INVESTIGACION EN MATEMATICAS                                                                    *
* DOCTORADO EN CIENCIAS DE LA COMPUTACION                                                                   *
* FERNANDO CERVANTES SANCHEZ                                                                                *
*                                                                                                           *
* FILE NAME: rand_fercer.cpp                                                                                *
*                                                                                                           *
* PURPOSE: Implementation of random number generators.                                                      *
*                                                                                                           *
* FILE REFERENCES:                                                                                          *
* Name        I/O        Description                                                                        *
* None        ----       ----------                                                                         *
*                                                                                                           *
* ABNORMAL TERMINATION CONDITIONS, ERROR AND WARNING MESSAGES:                                              *
* None                                                                                                      *
*                                                                                                           *
* NOTES:                                                                                                    *
*                                                                                                           *
* DEVELOPMENT HISTORY:                                                                                      *
* Date        Author        Change Id    Release    Description Of Change                                   *
* Nov/2016    Fernando C.   0            1.0        Creation                                                *
*                                                                                                           *
************************************************************************************************************/

#include "rand_fercer.h"




/************************************************************************************************************
*                                                                                                           *
* FUNCTION NAME: initSeed                                                                                   *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT        TYPE                I/O        DESCRIPTION                                                *
* initial_seed    unsigned int        input      The seed to initialize the random number generators        *
*                                                                                                           *
* RETURNS:                                                                                                  *
* An initialized seed for further use in random number generators                                           *
*                                                                                                           *
************************************************************************************************************/
STAUS* initSeed(unsigned int initial_seed) {
	if ( initial_seed == 0 ) {
#if defined(_WIN32) || defined(_WIN64)
		srand((unsigned int)time(NULL));
#else
		srand(clock());
#endif
		initial_seed = rand();
	}

	STAUS *my_seed = (STAUS*)malloc(sizeof(STAUS));
	my_seed->z1 = lcgR(&initial_seed);
	my_seed->z2 = lcgR(&initial_seed);
	my_seed->z3 = lcgR(&initial_seed);
	my_seed->lcg_seed = rand();

	return my_seed;
}






/***********************************************************************************************************/
/************************************************************************************************************
*                                                                                                           *
* FUNCTION NAME: lgcR                                                                                       *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT        TYPE                I/O        DESCRIPTION                                                *
* my_seed         unsigned int*       input      The seed to generate a random number                       *
*                                                                                                           *
* RETURNS:                                                                                                  *
* A random number between 0 and 2147483648                                                                  *
*                                                                                                           *
************************************************************************************************************/
unsigned int lcgR(unsigned int *my_seed) {
	return *my_seed = (1103515245u * *my_seed + 12345) % 2147483648;
}




/************************************************************************************************************
*                                                                                                           *
* FUNCTION NAME: tausStep                                                                                   *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT        TYPE                I/O        DESCRIPTION                                                *
* z               unsigned int*       input      Seed required to generate a random number                  *
* S1              unsigned int        input      Shift made to the seed 'z'                                 *
* S2              unsigned int        input      Shift made to the seed 'z'                                 *
* S3              unsigned int        input      Shift made to the seed 'z'                                 *
* M               unsigned int        input      Shift made to the seed 'z'                                 *
*                                                                                                           *
* RETURNS:                                                                                                  *
* A random number of type unsigned int                                                                      *
*                                                                                                           *
************************************************************************************************************/
unsigned int tausStep(unsigned int *z, const int S1, const int S2, const int S3, const unsigned int M) {
	const unsigned int b = (((*z << S1) ^ *z) >> S2);
	return *z = (((*z & M) << S3) ^ b);
}




/************************************************************************************************************
*                                                                                                           *
* FUNCTION NAME: HybTaus                                                                                    *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT        TYPE                I/O        DESCRIPTION                                                *
* par1            const double        input      Lower limit of the generated random number                 *
* par2            const double        input      Upper limit of the generated random number                 *
* my_seed         STAUS*              input      The seed used to generate the random number                *
*                                                                                                           *
* RETURNS:                                                                                                  *
* A random number between 'par1' and 'par2' of type double                                                  *
*                                                                                                           *
************************************************************************************************************/
double HybTaus(const double par1, const double par2, STAUS *my_seed) {
	double num = 2.3283064365387e-10 * (
		// Periods
		tausStep(&my_seed->z1, 13, 19, 12, 4294967294UL) ^ // p1=2^31-1
		tausStep(&my_seed->z2, 2, 25, 4, 4294967288UL) ^ // p2=2^30-1
		tausStep(&my_seed->z3, 3, 11, 17, 4294967280UL) ^ // p3=2^28-1
		lcgR(&my_seed->lcg_seed)	// p4=2^32
		);
	return (par2 - par1) * num + par1;
}






/************************************************************************************************************
*                                                                                                           *
* FUNCTION NAME: anorm_est                                                                                  *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT        TYPE                I/O        DESCRIPTION                                                *
* my_seed         STAUS*              input      The seed used to generate the random number                *
*                                                                                                           *
* RETURNS:                                                                                                  *
* A random from a normal standard distribution.                                                             *
*                                                                                                           *
************************************************************************************************************/
double anorm_est(STAUS *my_seed)
{
	double x1, x2, w, y1;//, y2;
	do {
		x1 = HybTaus(-1.0, 1.0, my_seed);
		x2 = HybTaus(-1.0, 1.0, my_seed);
		w = x1*x1 + x2*x2;
	} while (w >= 1.0);
	w = sqrt((-2.0 * log(w)) / w);
	y1 = x1 * w;
	//	y2 = x2 * w;
	return y1;
}






/************************************************************************************************************
*                                                                                                           *
* FUNCTION NAME: anorm_est                                                                                  *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT        TYPE                I/O        DESCRIPTION                                                *
* par1            const double         I         Parameter of localization of the normal distribution (mean)*
* par2            const double         I         Parameter of scale of the normal distribution (varianze)   *
* my_seed         STAUS*              input      The seed used to generate the random number                *
*                                                                                                           *
* RETURNS:                                                                                                  *
* A random from a normal distribution.                                                                      *
*                                                                                                           *
************************************************************************************************************/
double anorm(const double par1, const double par2, STAUS *my_seed)
{
	return sqrt(par2)*anorm_est(my_seed) + par1;
}