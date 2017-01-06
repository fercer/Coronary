/************************************************************************************************************
*                                                                                                           *
* CENTRO DE INVESTIGACION EN MATEMATICAS                                                                    *
* DOCTORADO EN CIENCIAS DE LA COMPUTACION                                                                   *
* FERNANDO CERVANTES SANCHEZ                                                                                *
*                                                                                                           *
* FILE NAME: performance_functions.h                                                                        *
*                                                                                                           *
* DEVELOPMENT HISTORY:                                                                                      *
* Date        Author        Change Id    Release    Description Of Change                                   *
* Nov/2016    Fernando C.   0            1.0        Creation                                                *
*                                                                                                           *
************************************************************************************************************/

#ifndef PERFORMANCE_FUNCTIONS_H_INCLUDED
#define PERFORMANCE_FUNCTIONS_H_INCLUDED

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <vector>

#if defined( _WIN32 ) || defined( _WIN64 )
	#include <Windows.h>
#else
	#include <sys/time.h>
#endif


#include "rand_fercer.h"
#include "filtros.h"
#include "IMGCONT.h"

#ifndef NDEBUG
#define MSGDEB( mensaje )   printf( "%s\n", mensaje)
#define MSGDEBPTR( mensaje, dato_ptr)   printf( "%s %p\n", mensaje, dato_ptr)
#define MSGDEBINT( mensaje, dato_int)   printf( "%s %i\n", mensaje, dato_int)
#define MSGDEBDBL( mensaje, dato_dbl)   printf( "%s %1.12f\n", mensaje, dato_dbl)
#define MSGDEBSTR( mensaje, dato_str)   printf( "%s %s", mensaje, dato_str);\
                                               printf("\n")
#define MSGDEBSTRCHR( mensaje, dato_str, n)   printf( "%s ", mensaje);\
                                               for( int i = 0; i < n; i++){\
                                                printf("%i ", (int)dato_str[i]);\
                                               }\
                                               printf("\n")
#define MSGDEBVEC( mensaje, dato_dbl, n)   printf( "%s {%f", mensaje, dato_dbl[0]);\
                                               for( int i = 1; i < n; i++){\
                                                   printf(", %f", dato_dbl[i]);\
                                               }\
                                               printf("}\n")
#else
#define MSGDEB( mensaje )
#define MSGDEBPTR( mensaje, dato_ptr)
#define MSGDEBINT( mensaje, dato_int)
#define MSGDEBDBL( mensaje, dato_dbl)
#define MSGDEBSTR( mensaje, dato_str)
#define MSGDEBSTRCHAR( mensaje, dato_str, n)
#define MSGDEBVEC( mensaje, dato_dbl, n)
#endif


#if defined( _WIN32 ) || defined( _WIN64 )
#define TIMERS LARGE_INTEGER t_ini, t_end, t_freq;\
				   QueryPerformanceFrequency(&t_freq)
#define GETTIME_INI QueryPerformanceCounter(&t_ini)
#define GETTIME_END QueryPerformanceCounter(&t_end)
#define DIFTIME (double)(t_end.QuadPart - t_ini.QuadPart)/(double)t_freq.QuadPart

#else
#define TIMERS struct timeval t_ini, t_end
#define GETTIME_INI gettimeofday( &t_ini, NULL)
#define GETTIME_END gettimeofday( &t_end, NULL)
#define DIFTIME ((t_end.tv_sec*1e6 + t_end.tv_usec) - (t_ini.tv_sec*1e6 + t_ini.tv_usec) )/ 1e6
#endif

#ifndef MY_INF
	#define MY_INF 179769313486231570814527423731704356798070567525844996598917476803157260780028538760589558632766878171540458953514382464234321326889464182768467546703537516986049910576551282076245490090389328944075868508455133942304583236903222948165808559332123348274797826204144723168738177180919299881250404026184124858368.0
#endif


class PERFORMANCE_FUNCTIONS {

private:
	std::vector< IMGCONT > * my_img_response;
	std::vector< IMGCONT > * my_img_response_mask;
	std::vector< IMGCONT > * my_img_thresholded_response;
	std::vector< IMGCONT > * my_img_groundtruth;

	unsigned int my_performance_imgs_count;

	double my_ROC_value;
	double my_ACC_value;

public:
	double calcROC();
	double calcCorCon();
	double calcAccuracy();

	void setInputResponse(std::vector<IMGCONT> * new_img_response );
	void setInputThreshold(std::vector<IMGCONT> * new_img_response_threshold);
	void setInputGroundtruth(std::vector<IMGCONT> * new_img_groundtruth);
	void setInputResponseMask(std::vector<IMGCONT>* new_img_response_mask);
};
#endif //PERFORMANCE_FUNCTIONS_H_INCLUDED
