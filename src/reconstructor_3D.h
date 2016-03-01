/******
    CENTRO DE INVESTIGACION EN MATEMATICAS
    MAESTRIA EN COMPUTACION Y MATEMATICAS INDUSTRIALES

    FERNANDO CERVANTES SANCHEZ
    SEP - 2015
*****/

#ifndef RECONSTRUCTOR_HPP_INCLUDED
#define RECONSTRUCTOR_HPP_INCLUDED


// Librerias de uso comun:
// Librerias de uso comun:
#include <vtkVersion.h>
#include <vtkSmartPointer.h>

// Librerias para trabajar con imagenes:
#include <vtkImageData.h>
#include <vtkImageReader2Factory.h>
#include <vtkImageReader2.h>
#include <vtkImageExtractComponents.h>
#include <vtkImageActor.h>
#include <vtkImageMapper3D.h>

// Librerias para generar mallas:
#include <vtkTriangle.h>
#include <vtkCellArray.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>

// Librerias para visualizacion
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkOutlineFilter.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkProperty.h>
#include <vtkSphereSource.h>



#include <gdcmImageReader.h>
#include <gdcmImage.h>
#include <gdcmReader.h>
#include <gdcmTag.h>

#include <assert.h>

#include <string.h>
#include <stdlib.h>

#include <iostream>
#include <fstream>

#include <omp.h>



#include "IMGVTK.h"
#include "filtros.h"


#ifdef _OPENMP
    #define TIMERS double t_ini, t_fin
    #define GETTIME_INI t_ini = omp_get_wtime()
    #define GETTIME_FIN t_fin = omp_get_wtime()
    #define DIFTIME (t_fin - t_ini)
#else
    #include <sys/time.h>
    #define TIMERS struct timeval t_ini, t_fin
    #define GETTIME_INI gettimeofday( &t_ini, NULL)
    #define GETTIME_FIN gettimeofday( &t_fin, NULL)
    #define DIFTIME ((t_fin.tv_sec*1e6 + t_fin.tv_usec) - (t_ini.tv_sec*1e6 + t_ini.tv_usec) )/ 1e6
    #define omp_get_num_threads() 1
    #define omp_set_num_threads(cores)
    #define omp_get_thread_num() 0
#endif


#define PI 3.14159265


#ifndef NDEBUG
    #define DEB_MSG(MENSAJE) using namespace std;\
                             cout << MENSAJE << endl;
#else
    #define DEB_MSG(MENSAJE)
#endif

// C L A S E: RECONS3D  ------------------------------------------------------------------------ v
class RECONS3D{

    private: //----------------------------------------------------------------------------- PRIVATE ----- v
    // M I E M B R O S      P R I V A D O S
        // Miembros para cargar las imagenes:
        IMGVTK img_base;
        IMGVTK img_delin;
        IMGVTK img_segment;

        bool esDICOM;

        // Miembros para visualizar la segmentacion 3D:
        vtkSmartPointer<vtkRenderer> mi_renderer;

        void renderizar();
        void mostrarImagen(vtkSmartPointer<vtkImageData> &imagen , const int nivel);
        void agregarEsfera( const double x, const double y, const double z, const double radio, double color[3] );


    // M E T O D O S       P R I V A D O S
    //-------------------------------------------------------------------------------------- PRIVATE ----- ^

    public: //------------------------------------------------------------------------------- PUBLIC----- v
    // M I E M B R O S      P U B L I C O S

    // C O N S T R U C T O R E S    /   D E S T R U C T O R E S
        RECONS3D(char **rutasbase_input , char **rutasground_input, const int n_imgs);
        RECONS3D(const char *rutabase_input, const char *rutadelin_input);

    // M E T O D O S        P U B L I C O S
        void segmentarImagenBase();
        void skeletonize();

    // O P E R A D O R E S  S O B R E C A R G A D O S
    //--------------------------------------------------------------------------------------- PUBLIC ----- ^
};
// C L A S E: RECONS3D  ------------------------------------------------------------------------- ^

#endif //RECONSTRUCTOR_HPP_INCLUDED
