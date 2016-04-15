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
#include <vtkTransform.h>
#include <vtkPointData.h>

// Arreglos de VTK:
#include <vtkUnsignedCharArray.h>
#include <vtkDoubleArray.h>

// Actores:
#include <vtkActor.h>
#include <vtkImageActor.h>
#include <vtkAxesActor.h>

// Mappers:
#include <vtkImageMapper3D.h>
#include <vtkDataSetMapper.h>
#include <vtkImageSliceMapper.h>
#include <vtkPolyDataMapper.h>

// Librerias para trabajar con imagenes:
#include <vtkImageData.h>
#include <vtkImageReader2Factory.h>
#include <vtkImageReader2.h>
#include <vtkImageExtractComponents.h>
#include <vtkImageSlice.h>
#include <vtkImageDataGeometryFilter.h>

// Librerias para generar mallas:
#include <vtkCellArray.h>
#include <vtkPolyData.h>
#include <vtkCellData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkStructuredGrid.h>
#include <vtkStructuredGridGeometryFilter.h>

// Librerias para visualizacion
#include <vtkCamera.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkOutlineFilter.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkProperty.h>

// Librerias para formas geometricas:
#include <vtkVertex.h>
#include <vtkLine.h>
#include <vtkTriangle.h>
#include <vtkSphereSource.h>
#include <vtkPyramid.h>
#include <vtkGlyph3D.h>


#include <gdcmImageReader.h>
#include <gdcmImage.h>
#include <gdcmReader.h>
#include <gdcmTag.h>
#include <gdcmPrivateTag.h>

#include <assert.h>

#include <string.h>
#include <stdlib.h>

#include <iostream>
#include <fstream>
#include <vector>

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
    // T I P O S        D E     D A T O S       Y       E S T R U C T U R A S       P R I V A D A S
        // POS
        typedef struct{
            double puntos[5][3];
        } POS;

    // M I E M B R O S      P R I V A D O S
        // Miembros para cargar las imagenes:
        std::vector< IMGVTK > imgs_base;
        std::vector< IMGVTK > imgs_delin;
        std::vector< bool > existe_ground;

        int n_angios;

        // Miembros para visualizar la segmentacion 3D:
        std::vector< vtkSmartPointer<vtkRenderer> > mis_renderers;
        vtkSmartPointer<vtkRenderer> renderer_global;
        std::vector< vtkSmartPointer<vtkPoints> > puntos;
        std::vector< vtkSmartPointer<vtkCellArray> > pixeles;

    // M E T O D O S       P R I V A D O S
        void renderizar(vtkSmartPointer<vtkRenderer> mi_renderer);
        void mostrarImagen(IMGVTK &imagen, IMGVTK::IMG_IDX img_idx, vtkSmartPointer<vtkRenderer> mi_renderer );
        void mostrarImagen(const int angio_ID, IMGVTK::IMG_IDX img_idx);
        void agregarEjes(vtkSmartPointer<vtkRenderer> mi_renderer);
        void agregarEsfera(const double x, const double y, const double z, const double radio, double color[3], vtkSmartPointer<vtkRenderer> mi_renderer );

        void mallarPuntos(const int angio_ID);
        POS posicionDefecto(const double ancho, const double alto, const double punta);
        void mostrarDetector(const int angio_ID);
    //-------------------------------------------------------------------------------------- PRIVATE ----- ^

    public: //------------------------------------------------------------------------------- PUBLIC----- v
    // M I E M B R O S      P U B L I C O S

    // C O N S T R U C T O R E S    /   D E S T R U C T O R E S
        RECONS3D();
        ~RECONS3D();

    // M E T O D O S        P U B L I C O S
        void agregarInput(char **rutasbase_input, char **rutasground_input, const int n_imgs);
        void agregarInput(const char *rutabase_input, const int nivel_l, const int nivel_u, const char *rutaground_input);

        void segmentarImagenBase();
        void segmentarImagenBase( const int angio_ID );
        void skeletonize();
        void skeletonize(const int angio_ID);

    // O P E R A D O R E S  S O B R E C A R G A D O S
    //--------------------------------------------------------------------------------------- PUBLIC ----- ^
};
// C L A S E: RECONS3D  ------------------------------------------------------------------------- ^

#endif //RECONSTRUCTOR_HPP_INCLUDED
