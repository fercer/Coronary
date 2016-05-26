/******
    CENTRO DE INVESTIGACION EN MATEMATICAS
    MAESTRIA EN COMPUTACION Y MATEMATICAS INDUSTRIALES

    FERNANDO CERVANTES SANCHEZ
    SEP - 2015
*****/

#ifndef RECONSTRUCTOR_HPP_INCLUDED
#define RECONSTRUCTOR_HPP_INCLUDED


// Librerias de uso comun con QT:
#include <QPlainTextEdit>


// Librerias de uso comun:
#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkTransform.h>
#include <vtkPointData.h>

// Arreglos de VTK:
#include <vtkUnsignedCharArray.h>
#include <vtkDoubleArray.h>
#include <vtkMath.h>
#include <vtkTransform.h>
#include <vtkMatrix4x4.h>
#include <vtkTransformPolyDataFilter.h>

// Actores:
#include <vtkActor.h>
#include <vtkActor2D.h>
#include <vtkImageActor.h>
#include <vtkAxesActor.h>

// Mappers:
#include <vtkImageMapper.h>
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
#include <vtkQuad.h>
#include <vtkTriangle.h>
#include <vtkSphereSource.h>
#include <vtkCylinderSource.h>
#include <vtkPyramid.h>
#include <vtkGlyph3D.h>
#include <vtkArrowSource.h>
#include <vtkSurfaceReconstructionFilter.h>
#include <vtkContourFilter.h>
#include <vtkReverseSense.h>


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


#define PI 3.14159265


// C L A S E: RECONS3D  ------------------------------------------------------------------------ v
class RECONS3D{

    private: //----------------------------------------------------------------------------- PRIVATE ----- v
    // T I P O S        D E     D A T O S       Y       E S T R U C T U R A S       P R I V A D A S
        // POS define la estructura para almacenar la posicion del artefacto que emula el detector.
        typedef struct{
            double puntos[5][3];
        } POS;

        // NORISO define la estructura para almacenar la normal y el isocentro de una proyeccion.
        typedef struct{
            double direccion[3]; //Las componentes de la normal
            double origen[3]; // El cenntro de la imagen.
        } NORCEN;

    // M I E M B R O S      P R I V A D O S
        QPlainTextEdit *mi_log;

        // Miembros para cargar las imagenes:
        std::vector< IMGVTK > imgs_base;
        std::vector< IMGVTK > imgs_delin;
        std::vector< bool > existe_ground;

        int n_angios;
        int detalle;

        // Miembros para visualizar la segmentacion 3D:
        std::vector< vtkSmartPointer<vtkRenderer> > mis_renderers;
        vtkSmartPointer<vtkRenderer> renderer_global;
        std::vector< vtkSmartPointer<vtkPoints> > puntos;
        std::vector< vtkSmartPointer<vtkCellArray> > pixeles;
        std::vector< NORCEN > normal_centros;

        FILTROS filtro;

    // M E T O D O S       P R I V A D O S
        void escribirLog( const char *mensaje );

        void renderizar(vtkSmartPointer<vtkRenderer> mi_renderer);
        void mostrarImagen(IMGVTK &imagen, IMGVTK::IMG_IDX img_idx, vtkSmartPointer<vtkRenderer> mi_renderer );
        void mostrarImagen(const int angio_ID, IMGVTK::IMG_IDX img_idx);


        void agregarVector(NORCEN org_dir, const double t, double color[], vtkSmartPointer<vtkRenderer> &mi_renderer);
        void agregarEjes(vtkSmartPointer<vtkRenderer> &mi_renderer);
        void agregarEsfera(const double x, const double y, const double z, const double radio, double color[3], vtkSmartPointer<vtkRenderer> mi_renderer );

        void isoCentro( const int angio_ID );
        void mallarPuntos(const int angio_ID);
        POS posicionDefecto(const double ancho, const double alto, const double punta);
        void mostrarDetector(const int angio_ID);

        void mostrarRadios(vtkSmartPointer<vtkPoints> puntos, vtkSmartPointer<vtkCellArray> cilindros, int *n_pix, IMGVTK::PIX_PAR *grafo, const double DDP, const double crl, const double srl, const double ccc, const double scc, const int nivel_detalle, FILE *fp_cilindros);
        void mostrarRadios(vtkSmartPointer<vtkPoints> puntos, vtkSmartPointer<vtkCellArray> vert_skl, vtkSmartPointer<vtkUnsignedCharArray> grafo_nivel, int *n_pix, IMGVTK::PIX_PAR *grafo, const double DDP, const double crl, const double srl, const double ccc, const double scc, const int n_niveles);
    //-------------------------------------------------------------------------------------- PRIVATE ----- ^

    public: //------------------------------------------------------------------------------- PUBLIC----- v
    // M I E M B R O S      P U B L I C O S

    // C O N S T R U C T O R E S    /   D E S T R U C T O R E S
        RECONS3D();
        ~RECONS3D();

    // M E T O D O S        P U B L I C O S
        void agregarInput(char **rutasbase_input, char **rutasground_input, const int n_imgs);
        void agregarInput(const char *rutabase_input, const int nivel_l, const int nivel_u, const char *rutaground_input);
        void agregarInput(const char *rutabase_input);

        void agregarGroundtruth(const char *rutaground_input, const int angio_ID);

        void segmentarImagenBase();
        void segmentarImagenBase(const int angio_ID );
        void skeletonize();
        void skeletonize(const int angio_ID);

        void mostrarBase( const int angio_ID );
        void mostrarGroundtruth(  const int angio_ID  );

        void setFiltroLimites(const FILTROS::PARAMETRO par, const double inf, const double sup, const double delta);
        void setFiltroEntrenamiento( const FILTROS::EVO_MET evo_met, const int m_iters, const int pob);
        void setFiltroEval( const FILTROS::FITNESS fit_fun);
        void setFiltroMetodo( const FILTROS::SEG_FILTRO metodo_filtrado);
        void setFiltroParametros( const FILTROS::PARAMETRO par, const double val );
        void setFiltroParametros();

        vtkSmartPointer< vtkRenderer > getRenderer();
        vtkSmartPointer< vtkRenderer > getRenderer( const int angio_ID );


        void setLog( QPlainTextEdit *log );
    // O P E R A D O R E S  S O B R E C A R G A D O S
    //--------------------------------------------------------------------------------------- PUBLIC ----- ^
};
// C L A S E: RECONS3D  ------------------------------------------------------------------------- ^

#endif //RECONSTRUCTOR_HPP_INCLUDED
