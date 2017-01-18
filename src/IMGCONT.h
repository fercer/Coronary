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
#include <vector>

#include <omp.h>
#include <math.h>

#include "libpng_src/readpng.h"

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

	/** THRESHOLD_ALG: **/
	typedef enum THRESHOLD_ALG { THRESH_LEVEL, THRESH_OTSU, THRESH_RIDLER_CALVARD } THRESHOLD_ALG;

	/** CONNECTED_ALG: **/
	typedef enum CONNECTED_ALG { CONN_DYN, CONN_ITER } CONNECTED_ALG;
	
	/** PIX_PAIR:   **/
	typedef struct PIX_PAIR {
		double my_pos_x;
		double my_pos_y;
		double my_x_r;
		double my_y_r;
		double my_radious;
		double my_angle_alpha;
		int my_deep_level;
		int my_n_children;
		PIX_PAIR *my_children[3];
		PIX_TYPE my_pix_type;
	} PIX_PAIR;

	IMGCONT();

	IMGCONT(const unsigned int new_height, const unsigned int new_width, const double init_val = 0.0);
	IMGCONT(const unsigned int new_height, const unsigned int new_width, const double * src_data);

	IMGCONT(const IMGCONT & img_src);

	~IMGCONT();

	IMGCONT & operator=(const IMGCONT & img_src);

	void Load(const char * src_path, const unsigned int level = 0);
	void Save(const char * out_path, const IMG_TYPE output_type = IMGPGM);

	void setImageData(const unsigned int new_height, const unsigned int new_width, const double * src_data);
	double getPix(const unsigned int row_y, const unsigned int col_x);
	void setPix(const unsigned int row_y, const unsigned int col_x, const double new_val);
	void setDimensions(const unsigned int new_height, const unsigned int new_width, const double init_val = 0.0);

	unsigned int getHeight();
	unsigned int getWidth();

	void lengthFilter(const unsigned int threshold_length, CONNECTED_ALG  my_connected_algorithm = CONN_DYN);

	void regionFill();

	void threshold(const THRESHOLD_ALG my_threshold_alg = THRESH_LEVEL, const double threshold_value = 0.5);

	double getMaximum();
	double getMinimum();

	void normalize(const double fixed_min, const double fixed_max);
	void normalize();

	double * getMask();
	double * getDistancesMap();
	double * getBoundaries();
	double * getSkeleton();
	PIX_PAIR * getSkeletonFeatures();
	int getSkeletonFeaturesDeep();

	void Rotate(const double my_rotation_tetha);

private:

	unsigned int my_height;       /* Height of the image */
	unsigned int my_width;        /* Width of the image */

	/************************************************************************/
	double * my_img_data;          /* The array where the image is contained */
	double * my_FOV_mask;
	double * my_dist_map;
	double * my_boundaries;
	double * my_skeleton;
	/************************************************************************/

	/* DICOM extracted information */
	double SID, SOD, DDP, DISO;
	double LAORAO, CRACAU;
	double WCenter, WWidth;
	double pixX, pixY, cenX, cenY;
	
	int my_max_distance;
	int my_skeleton_graph_deep;
	PIX_PAIR * my_skeleton_features;

	char my_err_msg[512];
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

	inline unsigned char erosionMask(double * erode_ptr, const int pos_x, const int pos_y);
	void erode(double * img_src);
	
	inline unsigned char dilMask(double * mask_dil_ptr, const unsigned int pos_x, const unsigned int pos_y);
	void fillMask();

	void computeMaskFOV();
	void computeMask();


	void computeDistancesMap();
	void computeBoundaries();

	bool regionFilling9(const unsigned int pos_x, const unsigned int pos_y);
	bool regionFilling7(const unsigned int pos_x, const unsigned int pos_y);
	bool regionFilling5(const unsigned int pos_x, const unsigned int pos_y);
	bool regionFilling3(const unsigned int pos_x, const unsigned int pos_y);

	double threshold_by_Otsu(const double min, const double max);
	double threshold_by_Ridler_and_Calvard(const double min_intensity, const double max_intensity);

	PIX_PAIR * computeSkeletonGraph(double * skl_tmp, const unsigned int pos_x, const unsigned int pos_y, int *nivel, const unsigned char *lutabla, bool *was_visited);

	void deleteSkeletonGraph(PIX_PAIR *graph_root);
	void extractSkeletonFeatures();

	PIX_PAIR * copySkeletonGraph(PIX_PAIR * src_graph_root);
	void copyFrom(const IMGCONT & img_src);

	inline unsigned char sklMask(double * skl_temp, const unsigned int pos_x, const unsigned int pos_y);
	void computeSkeleton();

	inline double linearInterpolation(const int pos_i, const int pos_j, const double mapping_y, const double mapping_x);
};

#endif // IMGCONT_H_INCLUDED
