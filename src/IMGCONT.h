/************************************************************************************************************
*                                                                                                           *
* CENTRO DE INVESTIGACION EN MATEMATICAS                                                                    *
* DOCTORADO EN CIENCIAS DE LA COMPUTACION                                                                   *
* FERNANDO CERVANTES SANCHEZ                                                                                *
*                                                                                                           *
* NAME:    IMGCONT.h                                                                                        *
*                                                                                                           *
* PURPOSE: Definition of the IMGCONT Class.                                                                 *
*                                                                                                           *
* GLOBAL VARIABLES:                                                                                         *
* Variable    Type    Description                                                                           *
* none        ---     ----------                                                                            *
*                                                                                                           *
* DEVELOPMENT HISTORY:                                                                                      *
* Date           Author        Change Id    Release        Description Of Change                            *
* 25/Dic/2016    Fernando C.   0            1.0            Creation                                         *
*                                                                                                           *
************************************************************************************************************/


#ifndef IMGCONT_H_INCLUDED
#define IMGCONT_H_INCLUDED

#ifdef BUILD_GDCM_VERSION
    #include <gdcmImageReader.h>
    #include <gdcmImage.h>
    #include <gdcmReader.h>
    #include <gdcmTag.h>
#endif

#if defined(_WIN32)
	#include MY_ZLIB
	#include MY_PNG
#else
	#include <png.h>
#endif

#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#include <iostream>

#include <omp.h>

#if defined(_WIN32) || defined(_WIN64)
	#include <time.h>
	#define COLOR_NORMAL
	#define COLOR_BOLD
	#define COLOR_UNDER
	#define COLOR_BLINK
	#define COLOR_INVERSE
	#define COLOR_BLACK
	#define COLOR_RED
	#define COLOR_GREEN
	#define COLOR_YELLOW
	#define COLOR_BLUE
	#define COLOR_MAGENTA
	#define COLOR_CYAN
	#define COLOR_WHITE
	#define COLOR_BACK_BLACK
	#define COLOR_BACK_RED
	#define COLOR_BACK_GREEN
	#define COLOR_BACK_YELLOW
	#define COLOR_BACK_BLUE
	#define COLOR_BACK_MAGENTA
	#define COLOR_BACK_CYAN
	#define COLOR_BACK_WHITE
	#define COLOR_RESET
#else
	#define COLOR_NORMAL            "\33[0m"
	#define COLOR_BOLD              "\33[1m"
	#define COLOR_UNDER             "\33[4m"
	#define COLOR_BLINK             "\33[5m"
	#define COLOR_INVERSE           "\33[7m"
	#define COLOR_BLACK             "\33[30m"
	#define COLOR_RED               "\33[31m"
	#define COLOR_GREEN             "\33[32m"
	#define COLOR_YELLOW            "\33[33m"
	#define COLOR_BLUE              "\33[34m"
	#define COLOR_MAGENTA           "\33[35m"
	#define COLOR_CYAN              "\33[36m"
	#define COLOR_WHITE             "\33[37m"
	#define COLOR_BACK_BLACK        "\33[40m"
	#define COLOR_BACK_RED          "\33[41m"
	#define COLOR_BACK_GREEN        "\33[42m"
	#define COLOR_BACK_YELLOW       "\33[43m"
	#define COLOR_BACK_BLUE         "\33[44m"
	#define COLOR_BACK_MAGENTA      "\33[45m"
	#define COLOR_BACK_CYAN         "\33[46m"
	#define COLOR_BACK_WHITE        "\33[47m"
	#define COLOR_RESET             "\x1b[0m"
#endif


#ifdef _OPENM
    #define TIMERS double t_ini, t_fin
    #define GETTIME_INI t_ini = omp_get_wtime()
    #define GETTIME_FIN t_fin = omp_get_wtime()
    #define DIFTIME (t_fin - t_ini)
    #define OMP_ENABLED true
#else
    #ifdef _WIN32
            #define TIMERS time_t t_ini, t_fin
            #define GETTIME_INI time (&t_ini)
            #define GETTIME_FIN time (&t_fin)
            #define DIFTIME difftime( t_fin, t_ini)
    #else
            #include <sys/time.h>
            #define TIMERS struct timeval t_ini, t_fin
            #define GETTIME_INI gettimeofday( &t_ini, NULL)
            #define GETTIME_FIN gettimeofday( &t_fin, NULL)
            #define DIFTIME ((t_fin.tv_sec*1e6 + t_fin.tv_usec) - (t_ini.tv_sec*1e6 + t_ini.tv_usec) )/ 1e6
    #endif
    #define omp_get_num_threads() 1
    #define omp_set_num_threads(cores)
    #define omp_get_thread_num() 0
    #define OMP_ENABLED false
#endif


#ifndef NDEBUG
	#define DEB_MSG(MENSAJE) std::cout << MENSAJE << std::endl;
#else
	#define DEB_MSG(MENSAJE)
#endif

#define MY_INF 179769313486231570814527423731704356798070567525844996598917476803157260780028538760589558632766878171540458953514382464234321326889464182768467546703537516986049910576551282076245490090389328944075868508455133942304583236903222948165808559332123348274797826204144723168738177180919299881250404026184124858368.0


class IMGCONT {
public:
	/*----------------------------------------------------------------------------- PUBLIC v ------------- */
	/** IMG_TYPE: **/
	typedef enum IMG_TYPE { IMGPNG, IMGPGM } IMG_TYPE;

	/** PIX_TYPE: **/
	typedef enum PIX_TYPE { PIX_SKL, PIX_END, PIX_BRANCH, PIX_CROSS } PIX_TYPE;

	/** THRESHOLD_TYPE: **/
	typedef enum THRESHOLD_TYPE { TRESH_LEVEL, TRESH_OTSU, TRESH_RIDLER_CALVARD } THRESHOLD_TYPE;

	/** CONNECTED_ALG: **/
	typedef enum CONNECTED_ALG { CONN_DYN, CONN_ITER } CONNECTED_ALG;

	IMGCONT();

	IMGCONT(const unsigned int new_height, const unsigned int new_width, const double init_val = 0.0);
	IMGCONT(const unsigned int new_height, const unsigned int new_width, const double * src_data);

	IMGCONT(const IMGCONT & img_src);

	~IMGCONT();

	IMGCONT & operator=(const IMGCONT & img_src);
	
	void Load(const char * src_path, const unsigned int level = 0);
	void Save(const char * out_path, const IMG_TYPE output_type = IMGPGM);

	double getPix(const unsigned int row_y, const unsigned int col_x);
	void setPix(const unsigned int row_y, const unsigned int col_x, const double new_val);
	void setDimensions(const unsigned int new_height, const unsigned int new_width, const double init_val = 0.0);


	int getHeight();
	int getWidth();

	void lengthFilter(const unsigned int threshold_length, CONNECTED_ALG  my_connected_algorithm = CONN_DYN);

	/*----------------------------------------------------------------------------- PUBLIC ^ ------------- */

protected:
	/*----------------------------------------------------------------------------- PROTECTED v ---------- */


	unsigned int my_height;       /* Height of the image */
	unsigned int my_width;        /* Width of the image */

	/************************************************************************/
	double *my_img_data;          /* The array where the image is contained */
	/************************************************************************/

	/* DICOM extracted information */
	double SID, SOD, DDP, DISO;
	double LAORAO, CRACAU;
	double WCenter, WWidth;
	double pixX, pixY, cenX, cenY;

	double * my_FOV_mask;

	void writeLog(const char *message);

	int LoadPNG(const char *src_path, const unsigned int level = 0);
	int LoadPGM(const char *src_path, const unsigned int level = 0);
	int LoadDICOM(const char *src_path, const unsigned int level = 0);

	void SavePGM(const char *out_path, const double my_min, const double my_max);
	void SavePNG(const char *out_path, const double my_min, const double my_max);

	void computeConnected(double * img_ptr, const int x, const int y, int *my_sets, unsigned int* number_of_labeled, bool* was_visited, const int number_of_labels);
	inline void increaseSetSize(int * my_labels, const int equiv_A, const int equiv_B, const int max_number_of_labels);
	unsigned int * connectedSets_Iterative(double * img_ptr, int * my_sets);
	unsigned int * connectedSets_Dynamic(double * img_ptr, int * my_sets);

	void lengthFilter(double * img_ptr, const unsigned int threshold_length, CONNECTED_ALG my_connected_algorithm = CONN_ITER);

	inline unsigned char erosionMask(double * ptr_tmp, const int x, const int y);
	void erode(double * img_src);

	void computeMaskFOV();
	void computeMask();
	/*----------------------------------------------------------------------------- PROTECTED ^ ---------- */
};

#endif // IMGCONT_H_INCLUDED
