/******
    CENTRO DE INVESTIGACION EN MATEMATICAS
    MAESTRIA EN COMPUTACION Y MATEMATICAS INDUSTRIALES

    FERNANDO CERVANTES SANCHEZ
    SEP - 2015
*****/

#ifndef RECONSTRUCTOR_HPP_INCLUDED
#define RECONSTRUCTOR_HPP_INCLUDED


// Librerias de uso comun:
#ifdef BUILD_VTK_VERSION
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
    #include <vtkChartXY.h>
    #include <vtkPlot.h>
    #include <vtkTable.h>
    #include <vtkIntArray.h>
    #include <vtkContextView.h>
    #include <vtkContextScene.h>
    #include <vtkAxis.h>

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
#endif

#ifdef BUILD_GCM_VERSION
    #include <gdcmImageReader.h>
    #include <gdcmImage.h>
    #include <gdcmReader.h>
    #include <gdcmTag.h>
    #include <gdcmPrivateTag.h>
#endif

#ifdef BUILD_GUI_VERSION
    #include "QVTKWidget.h"
    #include "QVTKInteractor.h"

    #include <QPlainTextEdit>
    #include <QTextEdit>
    #include <QProgressBar>
#endif

#include <assert.h>

#include <string.h>
#include <stdlib.h>

#include <iostream>
#include <fstream>
#include <vector>

#include <omp.h>

#include "args_fercer.h"
#include "IMGCONT.h"


#define MY_PI 3.14159265


// C L A S E: RECONS3D  ------------------------------------------------------------------------ v
class RECONS3D{
	public:
	typedef enum IMG_IDX { IMG_BASE, IMG_GROUNDTRUTH, IMG_MASK, IMG_BOUNDARIES, IMG_DIST_MAP, IMG_THRESHOLD, IMG_RESPONSE, IMG_SKELETON }IMG_IDX;

    private: //----------------------------------------------------------------------------- PRIVATE ----- v
    // T I P O S        D E     D A T O S       Y       E S T R U C T U R A S       P R I V A D A S

        // POS define la estructura para almacenar la posicion del artefacto que emula el detector.
        typedef struct POS {
            double puntos[5][3];
        } POS;

        // NORISO define la estructura para almacenar la normal y el isocentro de una proyeccion.
		typedef struct NORCEN {
			double direccion[3]; //Las componentes de la normal
			double origen[3]; // El cenntro de la imagen.
		} NORCEN;

		/* IMGCONT_NODE */
		typedef struct IMGCONT_NODE {
			IMGCONT *my_img_container;
			IMGCONT_NODE *my_next_node;
			IMGCONT_NODE *my_previous_node;
		} IMGCONT_NODE;

		/* IMGCONT_LIST */
		typedef struct IMGCONT_LIST {
			unsigned int my_nodes_count;
			IMGCONT_NODE *my_first_node;
			IMGCONT_NODE *my_last_node;
		} IMGCONT_LIST;

		/* IMGCONT_FEAT_SET */
		typedef struct IMGCONT_FEAT_SET {
			IMGCONT *my_base_ptr;
			IMGCONT *my_groundtruth_ptr;
			IMGCONT *my_mask_ptr;
			IMGCONT *my_response_ptr;
			IMGCONT *my_response_threshold_ptr;
			IMGCONT *my_thrs_response_map_dists_ptr;
			IMGCONT *my_thrs_response_skeleton_ptr;
			IMGCONT *my_thrs_response_boundaries_ptr;

			IMGCONT_FEAT_SET() {
				my_base_ptr = NULL;
				my_groundtruth_ptr = NULL;
				my_mask_ptr = NULL;
				my_response_ptr = NULL;
				my_response_threshold_ptr = NULL;
				my_thrs_response_map_dists_ptr = NULL;
				my_thrs_response_skeleton_ptr = NULL;
				my_thrs_response_boundaries_ptr = NULL;
			}

			IMGCONT_FEAT_SET(const IMGCONT_FEAT_SET & src_feat_set) {
				my_base_ptr = src_feat_set.my_base_ptr;
				my_groundtruth_ptr = src_feat_set.my_groundtruth_ptr;
				my_mask_ptr = src_feat_set.my_mask_ptr;
				my_response_ptr = src_feat_set.my_response_ptr;
				my_response_threshold_ptr = src_feat_set.my_response_threshold_ptr;
				my_thrs_response_map_dists_ptr = src_feat_set.my_thrs_response_map_dists_ptr;
				my_thrs_response_skeleton_ptr = src_feat_set.my_thrs_response_skeleton_ptr;
				my_thrs_response_boundaries_ptr = src_feat_set.my_thrs_response_boundaries_ptr;
			}

			IMGCONT_FEAT_SET & operator= (const IMGCONT_FEAT_SET & src_feat_set) {
				my_base_ptr = src_feat_set.my_base_ptr;
				my_groundtruth_ptr = src_feat_set.my_groundtruth_ptr;
				my_mask_ptr = src_feat_set.my_mask_ptr;
				my_response_ptr = src_feat_set.my_response_ptr;
				my_response_threshold_ptr = src_feat_set.my_response_threshold_ptr;
				my_thrs_response_map_dists_ptr = src_feat_set.my_thrs_response_map_dists_ptr;
				my_thrs_response_skeleton_ptr = src_feat_set.my_thrs_response_skeleton_ptr;
				my_thrs_response_boundaries_ptr = src_feat_set.my_thrs_response_boundaries_ptr;
			}


			~IMGCONT_FEAT_SET() {
				if (my_base_ptr) {
					delete my_base_ptr;
				}

				if (my_groundtruth_ptr) {
					delete my_groundtruth_ptr;
				}

				if (my_mask_ptr) {
					delete my_mask_ptr;
				}

				if (my_response_ptr) {
					delete my_response_ptr;
				}

				if (my_response_threshold_ptr) {
					delete my_response_threshold_ptr;
				}

				if (my_thrs_response_map_dists_ptr) {
					delete  my_thrs_response_map_dists_ptr;
				}

				if (my_thrs_response_skeleton_ptr) {
					delete my_thrs_response_skeleton_ptr;
				}

				if (my_thrs_response_boundaries_ptr) {
					delete my_thrs_response_boundaries_ptr;
				}
			}

		} IMGCONT_FEAT_SET;

    // M I E M B R O S      P R I V A D O S
        // Miembros para cargar las imagenes:
        IMGCONT_LIST *imgs_base;
		IMGCONT_LIST *imgs_groundtruth;
		IMGCONT_LIST *imgs_mask;
		IMGCONT_LIST *imgs_response;
		IMGCONT_LIST *imgs_response_threshold;
		IMGCONT_LIST *imgs_thrs_response_map_dists;
		IMGCONT_LIST *imgs_thrs_response_skeleton;
		IMGCONT_LIST *imgs_thrs_response_boundaries;

		std::vector< IMGCONT_FEAT_SET > my_img_features_ptr;

        std::vector< bool > existe_ground;

        int n_angios;
        int detalle;
        std::vector<int*> hist;
        std::vector<double> h_media, h_desvest, h_suma;

        // Miembros para visualizar la segmentacion 3D:
#ifdef BUILD_VTK_VERSION
        std::vector< vtkSmartPointer<vtkRenderer> > mis_renderers;
        vtkSmartPointer<vtkRenderer> renderer_global;
        std::vector< vtkSmartPointer<vtkPoints> > puntos;
        std::vector< vtkSmartPointer<vtkCellArray> > pixeles;
        std::vector< NORCEN > normal_centros;
        std::vector<  vtkSmartPointer<vtkContextView> > view;
#endif
		ARGUMENTS *my_args;
        FILE *fp_log;
        
#ifdef BUILD_GUI_VERSION
        QTextEdit *mi_txtLog;
        QProgressBar *mi_pBar;
#endif

    // M E T O D O S       P R I V A D O S
		IMGCONT * getIMGCONTPointer(const IMG_IDX my_img_idx, const unsigned int angios_ID);

        void escribirLog( const char *mensaje );
        void barraProgreso( const int avance, const int milestones );

#ifdef BUILD_VTK_VERSION
        void renderizar(vtkSmartPointer<vtkRenderer> mi_renderer);
        void mostrarImagen(IMGVTK::IMG_IDX img_idx, vtkSmartPointer<vtkRenderer> mi_renderer, const int angios_ID);
        void mostrarImagen(const int angios_ID, IMGVTK::IMG_IDX img_idx);

        void agregarVector(NORCEN org_dir, const double t, double color[], vtkSmartPointer<vtkRenderer> &mi_renderer);
        void agregarEjes(vtkSmartPointer<vtkRenderer> &mi_renderer);
        void agregarEsfera(const double x, const double y, const double z, const double radio, double color[3], vtkSmartPointer<vtkRenderer> mi_renderer );

        void isoCentro( const int angios_ID );
        void mallarPuntos(const int angios_ID);
        POS posicionDefecto(const double ancho, const double alto, const double punta);
        void mostrarDetector(const int angios_ID);

        void mostrarRadios(vtkSmartPointer<vtkPoints> puntos, vtkSmartPointer<vtkCellArray> cilindros, int *n_pix, IMGVTK::PIX_PAR *grafo, const double DDP, const double crl, const double srl, const double ccc, const double scc, const int nivel_detalle, FILE *fp_cilindros);
        void mostrarRadios(vtkSmartPointer<vtkPoints> puntos, vtkSmartPointer<vtkCellArray> vert_skl, vtkSmartPointer<vtkUnsignedCharArray> grafo_nivel, int *n_pix, IMGVTK::PIX_PAR *grafo, const double DDP, const double crl, const double srl, const double ccc, const double scc, const int n_niveles);
#endif
        void mostrarRadios(int *n_pix, IMGCONT::PIX_PAIR *grafo, const double DDP, const double crl, const double srl, const double ccc, const double scc, const int nivel_detalle, FILE *fp_cilindros);

    //-------------------------------------------------------------------------------------- PRIVATE ----- ^

    public: //------------------------------------------------------------------------------- PUBLIC----- v
    // M I E M B R O S      P U B L I C O S

    // C O N S T R U C T O R E S    /   D E S T R U C T O R E S
        RECONS3D();
		RECONS3D(ARGUMENTS *input_arguments);
        ~RECONS3D();

		void defineArguments();

    // M E T O D O S        P U B L I C O S
        void agregarInput(const char *rutabase_input, const int nivel_l, const int nivel_u, const char *rutaground_input, bool enmascarar);
        void agregarInput(const char *rutabase_input, bool enmascarar);

        void agregarGroundtruth(const char *rutaground_input, const int angios_ID);

        void leerConfiguracion( const char *ruta_conf);

        double medirExactitud(const int angios_ID);

		void clasAnchos(IMGCONT::PIX_PAIR *my_graph_root, const int angios_ID);

        void skeletonize(const int angios_ID);
        void skeletonize(const int angios_ID, const int nivel_detalle);

#ifdef BUILD_VTK_VERSION
        void mostrarBase( const int angios_ID );
        void mostrarGroundtruth(  const int angios_ID  );
#endif

        void setFiltroEntrenamiento(const FILTROS::EVO_MET evo_met);
        void setFiltroEntrenamientoPars(const FILTROS::EVO_MET_PAR evo_par, const double val);
        void setFiltroEval( const FILTROS::FITNESS fit_fun);
        void setFiltroMetodo( const FILTROS::SEG_FILTRO metodo_filtrado);

        void setFiltroParametros( const FILTROS::PARAMETRO par, const FILTROS::LIMITES lim, const double val);
        void setFiltroParametros( const FILTROS::PARAMETRO par, const double val );
        void setFiltroParametros();

        void Guardar(const char *ruta, IMG_IDX img_idx, IMGCONT::IMG_TYPE tipo_img, const int angios_ID);

        double getHist_desvest(const int angios_ID);
        double getHist_media(const int angios_ID);
        double getHist_suma(const int angios_ID);
        int getRows( const int angios_ID );
        int getCols( const int angios_ID );
        int getNangios();
        double *get_pixelData( const int angios_ID, IMG_IDX img_idx );

#ifdef BUILD_VTK_VERSION
        vtkRenderWindow *getHist(const int angios_ID);
        void setIteratorHist(vtkRenderWindowInteractor *interactor , const int angios_ID);
        vtkSmartPointer< vtkRenderer > getRenderer();
        vtkSmartPointer< vtkRenderer > getRenderer( const int angios_ID );
#endif

        void setLog( const char *ruta_log );
        void setFiltroLog( FILE *fplog );
        void setFiltroLog( const char* ruta_log );

#ifdef BUILD_GUI_VERSION
        void setLog( QTextEdit *txtLog );
        void setFiltroLog( QTextEdit *txtLog );
        void setProgressBar( QProgressBar *pBar );
        void setFiltroProgressBar( QProgressBar *pBar );
#endif

        void umbralizar(const IMGCONT::THRESHOLD_ALG tipo_umb, const double nivel, const int angios_ID);
        void lengthFilter(IMG_IDX my_img_idx, const int min_length, const int angios_ID);
	
    // O P E R A D O R E S  S O B R E C A R G A D O S
    //--------------------------------------------------------------------------------------- PUBLIC ----- ^
};
// C L A S E: RECONS3D  ------------------------------------------------------------------------- ^

#endif //RECONSTRUCTOR_HPP_INCLUDED
