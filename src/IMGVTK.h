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


#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#include <math.h>
#include <vector>
#include <algorithm>    // std::sort, std::fill_n
#include <iterator>

#include <iostream>

#include <omp.h>

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


#ifdef _WIN32
	#include <time.h>
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
#define INF 1e100;


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
        void definirMask(double *img_src, double *mask_src , const int mis_rens, const int mis_cols);
        void skeletonization(IMG_IDX img_idx);
        void umbralizar(IMG_IDX img_idx, const TIPO_UMBRAL tipo_umb, const double nivel);

        void lengthFilter(IMG_IDX img_idx, const int min_length, ALG_CONJUNTOS mi_alg);
        void regionFill(IMG_IDX img_idx);
        void mapaDistancias(IMG_IDX img_idx);
        void detectarBorde(IMG_IDX img_idx);
        double medirExactitud();

        void Cargar(const IMG_IDX img_idx, const char *ruta_origen, const bool enmascarar, const int nivel);
        void Cargar(const IMG_IDX img_idx, char **rutas , const int n_imgs, const bool enmascarar);

        void Guardar( IMG_IDX img_idx, const char *ruta, const TIPO_IMG tipo_salida );

//        void setLog( QPlainTextEdit *log );

        IMGVTK();
        IMGVTK(const IMGVTK &origen );
        IMGVTK(char **rutas_origen, const int n_imgs, const bool enmascarar);
        IMGVTK(const char *ruta_origen, const bool enmascarar, const int nivel);

        ~IMGVTK();

        // O P E R A D O R E S  S O B R E C A R G A D O S
        IMGVTK& operator= ( const IMGVTK &origen );

        // M I E M B R O S      P U B L I C O S
#ifdef BUILD_VTK_VERSION
        vtkSmartPointer<vtkImageData> base;
        vtkSmartPointer<vtkImageData> ground;
        vtkSmartPointer<vtkImageData> mask;
        vtkSmartPointer<vtkImageData> skeleton;
        vtkSmartPointer<vtkImageData> segment;
        vtkSmartPointer<vtkImageData> threshold;
        vtkSmartPointer<vtkImageData> mapa_dist;
        vtkSmartPointer<vtkImageData> borders;
#endif

        int rows, cols, rows_cols;
        int n_niveles;
        double *base_ptr;
        double *gt_ptr;
        double *mask_ptr;
        double *skl_ptr;
        double *segment_ptr;
        double *threshold_ptr;
        double *map_ptr;
        double *borders_ptr;

        //// Datos extraidos del archivo DICOM:
        double SID, SOD, DDP, DISO;
        double LAORAO, CRACAU;
        double WCenter, WWidth;
        double pixX, pixY, cenX, cenY;

        PIX_PAR *pix_caract;

    private: //----------------------------------------------------------------------------- PRIVATE ----- v
        // T I P O S        D E     D A T O S      P R I V A D O S
        /** ####:   **/
        // M E T O D O S      P R I V A D O S
        void escribirLog( const char *mensaje );
#ifdef BUILD_VTK_VERSION
        void Cargar(const char *ruta_origen, vtkSmartPointer<vtkImageData> img_src, vtkSmartPointer<vtkImageData> mask_src, const int nivel, const bool enmascarar);
        void Cargar(vtkSmartPointer<vtkImageData> img_src, vtkSmartPointer<vtkImageData> mask_src, char **rutas, const int n_imgs, const bool enmascarar);
#endif
        int *Cargar(const char *ruta_origen, double **img_src, double **mask_src, const int nivel, const bool enmascarar);
        int *Cargar(double **img_src, double **mask_src, char **rutas, const int n_imgs, const bool enmascarar);

        void maskFOV(double *img_tmp, double *mask_tmp, const int mis_cols, const int mis_rens);
        void fillMask(double *img_tmp, double *mask_tmp, const int mis_cols, const int mis_rens);

        PIX_PAR *grafoSkeleton(double *skl_tmp, const int x, const int y, int *nivel, const unsigned char *lutabla, bool *visitados);
        void extraerCaract(IMG_IDX img_idx);
        void borrarSkeleton( PIX_PAR *raiz );

        // ---------------- Mascaras para region filling
        bool regionFilling9( const double *ptr, const int x, const int y, const int mis_cols, const int mis_rens);
        bool regionFilling7( const double *ptr, const int x, const int y, const int mis_cols, const int mis_rens);
        bool regionFilling5( const double *ptr, const int x, const int y, const int mis_cols, const int mis_rens);
        bool regionFilling3( const double *ptr, const int x, const int y, const int mis_cols, const int mis_rens);
        void regionFill( double *ptr, const int mis_cols, const int mis_rens);

        double umbralizarOTSU(const double *img_ptr, const double min, const double max);
        double umbralizarRIDCAL( const double *img_ptr, const double min, const double max);

        inline unsigned char sklMask(const double *skl_ptr, const int x, const int y, const int mis_cols, const int mis_rens);

        inline unsigned char dilMask(const double *mask_dil, const int x, const int y , const int mis_cols, const int mis_rens);
        inline unsigned char erosionMask(const double *ptr_tmp, const int x, const int y, const int mis_cols, const int mis_rens);
        void erosionar(double *ptr, const int mis_cols, const int mis_rens);

        void conexo(const double *ptr, const int x, const int y, int *conjuntos, unsigned int* n_etiquetados, bool* visitados, const int num_etiquetas, const int mis_cols, const int mis_rens);
        unsigned int *conjuntosConexosDinamico(const double *ptr, int *conjuntos, const int mis_cols, const int mis_rens);

        inline void ampliarConjunto(int *etiquetas, const int equiv_A, const int equiv_B, const int max_etiquetas);
        unsigned int *conjuntosConexos(const double *ptr, int *conjuntos, const int mis_cols, const int mis_rens);
        void lengthFilter(double *ptr, const int min_length , const int mis_cols, const int mis_rens, ALG_CONJUNTOS mi_alg);

        char* setRuta( const char *ruta_input );

        // M I E M B R O S      P R I V A D O S
        //QPlainTextEdit *mi_log;
        int max_dist;
};
// C L A S E: IMGVTK  ---------------------------------------------------------------------------------------------- ^

#endif // IMG_PROC_H_INCLUDED
