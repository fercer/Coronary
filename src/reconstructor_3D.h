/******
    CENTRO DE INVESTIGACION EN MATEMATICAS
    MAESTRIA EN COMPUTACION Y MATEMATICAS INDUSTRIALES

    FERNANDO CERVANTES SANCHEZ
    SEP - 2015
*****/

#ifndef RECONSTRUCTOR_HPP_INCLUDED
#define RECONSTRUCTOR_HPP_INCLUDED


// Librerias de uso comun:
#include <vtkVersion.h>
#include <vtkSmartPointer.h>

// Librerias para trabajar con imagenes:
#include <vtkImageData.h>
#include <vtkImageReader2Factory.h>
#include <vtkImageReader2.h>
#include <vtkImageActor.h>

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


#include <gdcmImageReader.h>
#include <gdcmImage.h>
#include <gdcmReader.h>
#include <gdcmTag.h>


#include "IMGVTK.h"
#include "filtros.h"

// C L A S E: RECONS3D  ------------------------------------------------------------------------ v
class RECONS3D{

    private: //----------------------------------------------------------------------------- PRIVATE ----- v
    // M I E M B R O S      P R I V A D O S
        // Miembros para cargar las imagenes:
        IMGVTK img_base;
        IMGVTK img_delin;
        IMGVTK img_segment;

        gdcm::DataSet &ds;

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
