/******
    CENTRO DE INVESTIGACION EN MATEMATICAS
    MAESTRIA EN COMPUTACION Y MATEMATICAS INDUSTRIALES

    FERNANDO CERVANTES SANCHEZ
    ### - 20##
*****/

#ifndef IMGVTK_H_INCLUDED
#define IMGVTK_H_INCLUDED

// Librerias de uso comun con QT:
//#include <QPlainTextEdit>

#ifdef BUILD_VTK_VERSION
    // Librerias de uso comun:
    #include <vtkVersion.h>
    #include <vtkSmartPointer.h>

    // Librerias para trabajar con imagenes:
    #include <vtkImageData.h>
    #include <vtkImageReader2Factory.h>
    #include <vtkImageReader2.h>
    #include <vtkImageExtractComponents.h>
    #include <vtkPNGWriter.h>
#endif

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

#include <math.h>
#include <vector>
#include <algorithm>    // std::sort, std::fill_n
#include <iterator>

#include <iostream>

#include <omp.h>

#include "args_fercer.h"

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
	#define COLOR_INVERSE		"\33[7m"
	#define COLOR_BLACK             "\33[30m"
	#define COLOR_RED               "\33[31m"
	#define COLOR_GREEN             "\33[32m"
	#define COLOR_YELLOW		"\33[33m"
	#define COLOR_BLUE              "\33[34m"
	#define COLOR_MAGENTA		"\33[35m"
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
    #define DEB_MSG(MENSAJE) using namespace std;\
                             cout << MENSAJE << endl
#else
    #define DEB_MSG(MENSAJE)
#endif


#define Mi_PI 3.14159265
#define MY_INF 179769313486231570814527423731704356798070567525844996598917476803157260780028538760589558632766878171540458953514382464234321326889464182768467546703537516986049910576551282076245490090389328944075868508455133942304583236903222948165808559332123348274797826204144723168738177180919299881250404026184124858368.0


/* Source libpng writepng.h */
#define TEXT_TITLE    0x01
#define TEXT_AUTHOR   0x02
#define TEXT_DESC     0x04
#define TEXT_COPY     0x08
#define TEXT_EMAIL    0x10
#define TEXT_URL      0x20

#define TEXT_TITLE_OFFSET        0
#define TEXT_AUTHOR_OFFSET      72
#define TEXT_COPY_OFFSET     (2*72)
#define TEXT_EMAIL_OFFSET    (3*72)
#define TEXT_URL_OFFSET      (4*72)
#define TEXT_DESC_OFFSET (5*72)


// C L A S E: IMGVTK  ---------------------------------------------------------------------------------------------- v
class IMGVTK{
    public: //----------------------------------------------------------------------------- PUBLIC ------- v
        // T I P O S        D E     D A T O S       Y       E S T R U C T U R A S      P U B L I C A S        // IMG_IDX
        /** IMG_IDX:   **/
        typedef enum IMG_IDX { BASE, GROUNDTRUTH, MASK, SKELETON, SEGMENT, THRESHOLD, MAPDIST, BORDERS } IMG_IDX;

        /** TIPO_IMG:   **/
        typedef enum TIPO_IMG { PNG, PGM } TIPO_IMG;

        /** TIPO_CARACT:   **/
        typedef enum TIPO_CARACT { PIX_SKL, PIX_END, PIX_BRANCH, PIX_CROSS } TIPO_CARACT;

        /** TIPO_UMBRAL:    **/
        typedef enum TIPO_UMBRAL { NIVEL, OTSU, RIDLER_CALVARD} TIPO_UMBRAL;

        /** ALG_CONJUNTOS:    **/
        typedef enum ALG_CONJUNTOS { DINAMICO, ITERATIVO} ALG_CONJUNTOS;
		

		/* Source: libpng  writepng.c */
		typedef struct _mainprog_info {
			double gamma;
			long width;
			long height;
			time_t modtime;
			FILE *infile;
			FILE *outfile;
			void *png_ptr;
			void *info_ptr;
			unsigned char *image_data;
			unsigned char **row_pointers;
			char *title;
			char *author;
			char *desc;
			char *copyright;
			char *email;
			char *url;
			int filter;    /* command-line-filter flag, not PNG row filter! */
			int pnmtype;
			int sample_depth;
			int interlaced;
			int have_bg;
			int have_time;
			int have_text;
			jmp_buf jmpbuf;
			unsigned char bg_red;
			unsigned char bg_green;
			unsigned char bg_blue;
		} mainprog_info;



		/** IMGCONT: IMaGe CONTainer  **/
		typedef struct IMGCONT {
			unsigned int my_height;    /* Height of the image */
			unsigned int my_width;     /* Width of the image */
			unsigned int my_offset_x;  /* Offset in X axis */
			unsigned int my_offset_y;  /* Offset in Y axis */
			double *my_img_data;       /* The array where the image is contained */

			IMGCONT(const unsigned int new_height, const unsigned int new_width, const unsigned int new_off_x = 0,
				const unsigned int new_off_y = 0, const double init_val = 0.0) {
				my_height = new_height ;
				my_width = new_width;
				my_offset_x = new_off_x;
				my_offset_y = new_off_y;
				my_img_data = (double *)malloc((my_height + 2 * my_offset_y) * (my_width + 2 * my_offset_x)
					* sizeof(double));

				for (unsigned int xy = 0; xy < (my_height + 2 * my_offset_y)*(my_width + 2 * my_offset_x); xy++) {
					*(my_img_data + xy) = init_val;
				}
			}
			
			IMGCONT( const IMGCONT& img_src) {
				my_height = img_src.my_height;
				my_width = img_src.my_width;
				my_offset_x = img_src.my_offset_x;
				my_offset_y = img_src.my_offset_y;

				if (my_img_data) {
					free(my_img_data);
				}

				my_img_data = (double*)malloc(my_height * my_width * sizeof(double));

				memcpy(my_img_data, img_src.my_img_data, my_height * my_width * sizeof(double));
			}

			~IMGCONT() {
				free(my_img_data);
			}

			IMGCONT & operator= (const IMGCONT & img_src) {

				if ((my_height > img_src.my_height) || (my_width > img_src.my_width)) {

					my_offset_x = (my_width - img_src.my_width) / 2;
					my_offset_y = (my_height - img_src.my_height) / 2;

					for (unsigned int y = 0; y = img_src.my_height; y++) {
						for (unsigned int x = 0; x = img_src.my_width; x++) {
							*(my_img_data + (y + my_offset_y)*my_width + x + my_offset_x) = 
								*(img_src.my_img_data + y*img_src.my_width + x);
						}
					}

					my_offset_x += img_src.my_offset_x;
					my_offset_y += img_src.my_offset_y;

				}
				else if ((my_height < img_src.my_height) || (my_width < img_src.my_width)) {
					my_height = img_src.my_height;
					my_width = img_src.my_width;
					if (my_img_data) {
						free(my_img_data);
					}

					my_offset_x = img_src.my_offset_x;
					my_offset_y = img_src.my_offset_y;

					my_img_data = (double*)malloc(my_height * my_width * sizeof(double));
				}
				else {
					my_offset_x = img_src.my_offset_x;
					my_offset_y = img_src.my_offset_y;

					memcpy(my_img_data, img_src.my_img_data, my_height * my_width * sizeof(double));
				}

				return *this;
			}
		} IMGCONT;

        /** PIX_PAR:   **/
        typedef struct PIX_PAR {
            double x, y, x_r, y_r;
            double radio, alpha;
            int nivel;
            int n_hijos;
            PIX_PAR *hijos[3];
            TIPO_CARACT pix_tipo;
        } PIX_PAR;


        // M E T O D O S      P U B L I C O S
		IMGCONT * definirMask( IMGCONT *img_src);
        void skeletonization(IMG_IDX img_idx);
        void umbralizar(IMG_IDX img_idx, const TIPO_UMBRAL tipo_umb, const double nivel);

        void lengthFilter(IMG_IDX img_idx, const int min_length, ALG_CONJUNTOS mi_alg);
        void regionFill(IMG_IDX img_idx);
        void mapaDistancias(IMG_IDX img_idx);
        void detectarBorde(IMG_IDX img_idx);
        double medirExactitud();

        void Cargar(const IMG_IDX img_idx, const char *ruta_origen, const bool enmascarar, const int nivel);
		void Guardar(IMG_IDX img_idx, const char *ruta, const TIPO_IMG tipo_salida);

//        void setLog( QPlainTextEdit *log );

        IMGVTK();
		IMGVTK(const IMGVTK &origen );
		IMGVTK(const char *ruta_origen, const bool enmascarar, const int nivel);

        ~IMGVTK();

        // O P E R A D O R E S  S O B R E C A R G A D O S
        IMGVTK& operator= ( const IMGVTK &origen );

        // M I E M B R O S      P U B L I C O S
    private: //----------------------------------------------------------------------------- PRIVATE ----- v
        // T I P O S        D E     D A T O S      P R I V A D O S

        // M E T O D O S      P R I V A D O S
        void escribirLog( const char *mensaje );

		IMGCONT* CargarPNG(const char *ruta_origen);
		IMGCONT* CargarPGM(const char *ruta_origen);
		IMGCONT* CargarDICOM(const char *ruta_origen, const unsigned int nivel);

		void GuardarPGM(IMGCONT *img_src, const char *ruta, const double min, const double max);
		void GuardarPNG(IMGCONT *img_src, const char *ruta, const double min, const double max);

		IMGCONT* maskFOV(IMGCONT * img_src);
		void fillMask(IMGCONT *img_src, IMGCONT *mask_src);

        PIX_PAR *grafoSkeleton(double *skl_tmp, const int x, const int y, int *nivel, const unsigned char *lutabla, bool *visitados);
        void extraerCaract(IMG_IDX img_idx);
        void borrarSkeleton( PIX_PAR *raiz );

        // ---------------- Mascaras para region filling
		bool regionFilling9(IMGCONT *img_src, const int x, const int y);
        bool regionFilling7(IMGCONT *img_src, const int x, const int y);
        bool regionFilling5(IMGCONT *img_src, const int x, const int y);
        bool regionFilling3(IMGCONT *img_src, const int x, const int y);
		void regionFill(IMGCONT *img_src);

        double umbralizarOTSU(const double *img_ptr, const double min, const double max);
        double umbralizarRIDCAL( const double *img_ptr, const double min, const double max);

        inline unsigned char sklMask(IMGCONT *skl_ptr, const int x, const int y);

		inline unsigned char dilMask(IMGCONT *mask_dil, const int x, const int y);

        inline unsigned char erosionMask(IMGCONT *ptr_tmp, const int x, const int y);
        void erosionar(IMGCONT *img_src);

		void conexo(IMGCONT *img_src, const int x, const int y, int *conjuntos, unsigned int* n_etiquetados, bool* visitados, const int num_etiquetas);

		unsigned int *conjuntosConexosDinamico(IMGCONT *img_src, int *conjuntos);

        inline void ampliarConjunto(int *etiquetas, const int equiv_A, const int equiv_B, const int max_etiquetas);

        unsigned int *conjuntosConexos(IMGCONT *img_src, int *conjuntos);

		void lengthFilter(IMGCONT *img_src, const unsigned int min_length, ALG_CONJUNTOS mi_alg);

        char* setRuta( const char *ruta_input );
		
        // M I E M B R O S      P R I V A D O S
        //QPlainTextEdit *mi_log;
        int max_dist;
		char err_msg[512];

		int rows, cols, rows_cols;
		int n_niveles;
		IMGCONT *my_base;
		IMGCONT *my_groundtruth;
		IMGCONT *my_mask;
		IMGCONT *my_skeleton;
		IMGCONT *my_response;
		IMGCONT *my_segmented;
		IMGCONT *my_distmap;
		IMGCONT *my_boundaries;

		//// Datos extraidos del archivo DICOM:
		double SID, SOD, DDP, DISO;
		double LAORAO, CRACAU;
		double WCenter, WWidth;
		double pixX, pixY, cenX, cenY;

		PIX_PAR *pix_caract;
};
// C L A S E: IMGVTK  ---------------------------------------------------------------------------------------------- ^

#endif // IMG_PROC_H_INCLUDED
