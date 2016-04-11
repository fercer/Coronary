/******
    CENTRO DE INVESTIGACION EN MATEMATICAS
    MAESTRIA EN COMPUTACION Y MATEMATICAS INDUSTRIALES

    FERNANDO CERVANTES SANCHEZ
    ### - 20##
*****/

#ifndef IMGVTK_H_INCLUDED
#define IMGVTK_H_INCLUDED

// Librerias de uso comun:
#include <vtkVersion.h>
#include <vtkSmartPointer.h>

// Librerias para trabajar con imagenes:
#include <vtkImageData.h>
#include <vtkImageReader2Factory.h>
#include <vtkImageReader2.h>
#include <vtkImageExtractComponents.h>
#include <vtkPNGWriter.h>

#include <gdcmImageReader.h>
#include <gdcmImage.h>
#include <gdcmReader.h>
#include <gdcmTag.h>

#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#include <math.h>
#include <vector>
#include <algorithm>    // std::sort, std::fill_n
#include <iterator>

#include <iostream>

#include <omp.h>

#ifdef _OPENMP
    #define TIMERS double t_ini, t_fin
    #define GETTIME_INI t_ini = omp_get_wtime()
    #define GETTIME_FIN t_fin = omp_get_wtime()
    #define DIFTIME (t_fin - t_ini)
    #define OMP_ENABLED true
#else
    #include <sys/time.h>
    #define TIMERS struct timeval t_ini, t_fin
    #define GETTIME_INI gettimeofday( &t_ini, NULL)
    #define GETTIME_FIN gettimeofday( &t_fin, NULL)
    #define DIFTIME ((t_fin.tv_sec*1e6 + t_fin.tv_usec) - (t_ini.tv_sec*1e6 + t_ini.tv_usec) )/ 1e6
    #define omp_get_num_threads() 1
    #define omp_set_num_threads(cores)
    #define omp_get_thread_num() 0
    #define OMP_ENABLED false
#endif


#ifndef NDEBUG
    #define DEB_MSG(MENSAJE) using namespace std;\
                             cout << MENSAJE << endl;
#else
    #define DEB_MSG(MENSAJE)
#endif


// C L A S E: IMGVTK  ---------------------------------------------------------------------------------------------- v
class IMGVTK{
    public: //----------------------------------------------------------------------------- PUBLIC ------- v
        // T I P O S        D E     D A T O S       Y       E S T R U C T U R A S      P U B L I C A S        // IMG_IDX
        /** IMG_IDX:   **/
        typedef enum{ BASE, MASK, SKELETON, SEGMENT } IMG_IDX;


        /** TIPO_IMG:   **/
        typedef enum{ PNG, PGM } TIPO_IMG;

        /** TIPO_CARACT:   **/
        typedef enum{ PIX_END, PIX_BRANCH, PIX_CROSS } TIPO_CARACT;

        /** PIX_PAR:   **/
        typedef struct{
            unsigned int x, y;
            TIPO_CARACT pix_tipo;
        } PIX_PAR;


        // M E T O D O S      P U B L I C O S
        void definirMask( vtkSmartPointer<vtkImageData> img_src, vtkSmartPointer<vtkImageData> mask_src );
        void skeletonization(IMG_IDX img_idx);
        void umbralizar(IMG_IDX img_idx);
        void lengthFilter(const int min_length);
        void regionFill();

        void Cargar(const char *ruta_origen, const bool enmascarar, const int nivel);
        void Cargar(char **rutas , const int n_imgs, const bool enmascarar);

        void Guardar(const char *ruta, const TIPO_IMG tipo_salida );

        IMGVTK();
        IMGVTK(const IMGVTK &origen );
        IMGVTK(char **rutas_origen, const int n_imgs, const bool enmascarar);
        IMGVTK(const char *ruta_origen, const bool enmascarar, const int nivel);

        ~IMGVTK();

        // O P E R A D O R E S  S O B R E C A R G A D O S
        IMGVTK& operator= ( const IMGVTK &origen );

        // M I E M B R O S      P U B L I C O S
        vtkSmartPointer<vtkImageData> base;
        vtkSmartPointer<vtkImageData> mask;
        vtkSmartPointer<vtkImageData> skeleton;
        vtkSmartPointer<vtkImageData> segment;

        int rens, cols, rens_cols;
        int n_caracts;
        double *base_ptr;
        double *skl_ptr;
        double *mask_ptr;
        double *segment_ptr;

        //// Datos extraidos del archivo DICOM:
        double SID, SOD, DDP;
        double LAORAO, CRACAU;
        double WCenter, WWidth;
        double pixX, pixY, cenX, cenY;
        bool esDICOM;

        PIX_PAR *pix_caract;

    private: //----------------------------------------------------------------------------- PRIVATE ----- v
        // T I P O S        D E     D A T O S      P R I V A D O S
        /** ####:   **/
        // M E T O D O S      P R I V A D O S

        void Cargar(const char *ruta_origen, vtkSmartPointer<vtkImageData> img_src, vtkSmartPointer<vtkImageData> mask_src, const int nivel, const bool enmascarar);
        void Cargar(vtkSmartPointer<vtkImageData> img_src, vtkSmartPointer<vtkImageData> mask_src, char **rutas, const int n_imgs, const bool enmascarar);

        void mapaDistancias();

        void maskFOV(double *img_tmp, double *mask_tmp, const int mis_cols, const int mis_rens);
        void fillMask(double *img_tmp, double *mask_tmp, const int mis_cols, const int mis_rens);

        void extraerCaract();

        // ---------------- Mascaras para region filling
        bool regionFilling9( const double *ptr, const int x, const int y, const int mis_cols, const int mis_rens);
        bool regionFilling7( const double *ptr, const int x, const int y, const int mis_cols, const int mis_rens);
        bool regionFilling5( const double *ptr, const int x, const int y, const int mis_cols, const int mis_rens);
        bool regionFilling3( const double *ptr, const int x, const int y, const int mis_cols, const int mis_rens);
        void regionFill( double *ptr, const int mis_cols, const int mis_rens);

        inline unsigned char sklMask(const double *skl_ptr, const int x, const int y, const int mis_cols, const int mis_rens);

        inline unsigned char dilMask(const double *mask_dil, const int x, const int y , const int mis_cols, const int mis_rens);
        inline unsigned char erosionMask(const double *ptr_tmp, const int x, const int y, const int mis_cols, const int mis_rens);
        void erosionar(double *ptr, const int mis_cols, const int mis_rens);

        void conexo(const double *ptr, const int x, const int y, int *conjuntos, unsigned int* n_etiquetados, bool* visitados, const int num_etiquetas, const int mis_cols, const int mis_rens);
        unsigned int *conjuntosConexosDinamico(const double *ptr, int *conjuntos, const int mis_cols, const int mis_rens);

        unsigned int *conjuntosConexos(const double *ptr, int *conjuntos, const int mis_cols, const int mis_rens);
        void lengthFilter(double *ptr, const int min_length , const int mis_cols, const int mis_rens);

        char* setRuta( const char *ruta_input );

        // M I E M B R O S      P R I V A D O S
        vtkSmartPointer<vtkImageData> mapa_dist;
        double *map_ptr;


};
// C L A S E: IMGVTK  ---------------------------------------------------------------------------------------------- ^

#endif // IMG_PROC_H_INCLUDED
