/******
    CENTRO DE INVESTIGACION EN MATEMATICAS
    MAESTRIA EN COMPUTACION Y MATEMATICAS INDUSTRIALES

    FERNANDO CERVANTES SANCHEZ
    SEP - 2015
*****/


#include "reconstructor_3D.h"

// C L A S E: RECONS3D  ------------------------------------------------------------------------ v

//------------------------------------------------------------------------------------------------ PRIVATE ----- v
// M E T O D O S       P R I V A D O S



/*  Metodo: renderizar()

    Funcion: Renderiza los 'actores' contenidos en 'renderer' en una ventana de VTK.
*/
void RECONS3D::renderizar( vtkSmartPointer<vtkRenderer> mi_renderer ){

    mi_renderer->ResetCamera();

    // Crear una ventana temporal:
    vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow->AddRenderer(mi_renderer);

    // Crear interactuador temporal:
    vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();

    // Crear la ventana y mostrar la renderizacion:
    renderWindowInteractor->SetRenderWindow(renderWindow);
    //renderWindowInteractor->Initialize();

    renderWindowInteractor->Start();
}




/*  Metodo: mostrarImagen( vtkSmartPointer<vtkImageData> imagen )

    Funcion: Muestra la imagen en una ventana VTK.
*/
void RECONS3D::mostrarImagen( vtkSmartPointer<vtkImageData> &img_src, vtkSmartPointer<vtkRenderer> mi_renderer){
    // Crear un actor temporal
    vtkSmartPointer<vtkImageActor> actor = vtkSmartPointer<vtkImageActor>::New();

    #if VTK_MAJOR_VERSION <= 5
        actor->SetInput(img_src);
    #else
        actor->SetInputData(img_src);
    #endif

    // Agregar el actor de la imagen al renderizador.
    mi_renderer->AddActor(actor);
    mi_renderer->ResetCamera();
}




/*  Metodo: agregarEsfera

    Funcion: Agrega una esfera al renderizador.
*/
void RECONS3D::agregarEsfera(const double x, const double y, const double z, const double radio, double color[3], vtkSmartPointer<vtkRenderer> mi_renderer){
    vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();

    sphereSource->SetCenter(x, y, z);
    sphereSource->SetRadius(radio);

    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputConnection(sphereSource->GetOutputPort());

    // Crear un actor temporal
    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);

    actor->GetProperty()->SetColor(color);

    // Agregar el actor de la imagen al renderizador.
    mi_renderer->AddActor(actor);
    mi_renderer->ResetCamera();
}


//------------------------------------------------------------------------------------------------ PRIVATE ----- ^
//-------------------------------------------------------------------------------------------------- PUBLIC----- v

// M I E M B R O S      P U B L I C O S
/*  Metodo: agregarInput

    Funcion: Define las rutas de las imagenes que son usadas para reconstruir una arteria coronaria.
*/
void RECONS3D::agregarInput(char **rutasbase_input, char **rutasground_input, const int n_imgs){
    esDICOM.push_back(false);

    n_angios++;
    IMGVTK *imgs_temp = imgs_base;
    imgs_base = new IMGVTK [n_angios];
    memcpy(imgs_base, imgs_temp, (n_angios-1)*sizeof(IMGVTK));

    // Cargar la imagen base como la concatenacion de todas las imagenes base:
    imgs_base[n_angios-1].Cargar(rutasbase_input, n_imgs, true);

    imgs_temp = imgs_delin;
    imgs_delin = new IMGVTK [n_angios];
    memcpy(imgs_delin, imgs_temp, (n_angios-1)*sizeof(IMGVTK));
    // Cargar la imagen delineada como la concatenacion de todas las imagenes delineadas:
    imgs_delin[n_angios-1].Cargar(rutasground_input, n_imgs, false);

    imgs_temp = imgs_segment;
    imgs_segment = new IMGVTK [n_angios];
    memcpy(imgs_segment, imgs_temp, (n_angios-1)*sizeof(IMGVTK));

    // Mostrar la imagen en un renderizador
    mis_renderers.push_back(vtkSmartPointer<vtkRenderer>::New());
    mostrarImagen(imgs_base[n_angios-1].base, mis_renderers[n_angios-1]);
}



/*  Metodo: agregarInput

    Funcion: Define las rutas de las imagenes que son usadas para reconstruir una arteria coronaria.
*/
void RECONS3D::agregarInput(const char *rutabase_input, const char *rutaground_input, const int nivel){
    n_angios++;

    esDICOM.push_back(true);

    const int ruta_l = strlen(rutabase_input);
    DEB_MSG("Extension del archivo de entrada: " << (rutabase_input + ruta_l - 3));
    esDICOM[n_angios-1] = esDICOM[n_angios-1] * strcmp(rutabase_input + ruta_l - 3, "png");
    esDICOM[n_angios-1] = esDICOM[n_angios-1] * strcmp(rutabase_input + ruta_l - 3, "jpg");
    esDICOM[n_angios-1] = esDICOM[n_angios-1] * strcmp(rutabase_input + ruta_l - 4, "jpeg");
    esDICOM[n_angios-1] = esDICOM[n_angios-1] * strcmp(rutabase_input + ruta_l - 3, "bmp");

    IMGVTK *imgs_temp = imgs_base;
    imgs_base = new IMGVTK [n_angios];
    memcpy(imgs_base, imgs_temp, (n_angios-1)*sizeof(IMGVTK));

    imgs_temp = imgs_delin;
    imgs_delin = new IMGVTK [n_angios];
    memcpy(imgs_delin, imgs_temp, (n_angios-1)*sizeof(IMGVTK));

    imgs_temp = imgs_segment;
    imgs_segment = new IMGVTK [n_angios];
    memcpy(imgs_segment, imgs_temp, (n_angios-1)*sizeof(IMGVTK));

    if(esDICOM[n_angios-1]){

        // Abrir el archivo DICOM
        gdcm::ImageReader DICOMreader;
        DICOMreader.SetFileName( rutabase_input );
        DICOMreader.Read();

        gdcm::File &file = DICOMreader.GetFile();
        gdcm::DataSet &ds = file.GetDataSet();
        std::stringstream strm;
        strm.str("");
        DEB_MSG("TAG[" << (0x12) << "," << (0x456) << "]: ");
        if( ds.FindDataElement(gdcm::PrivateTag (0x12, 0x456)) ){
            ds.GetDataElement( gdcm::PrivateTag (0x12, 0x456) ).GetValue().Print(strm);
            DEB_MSG( strm.str() );
        }else{
            DEB_MSG("No funciona . . .");
        }
        const gdcm::Image &gimage = DICOMreader.GetImage();
        imgs_base[n_angios-1].Cargar(gimage, nivel);

    }else{
        imgs_base[n_angios-1].Cargar(rutabase_input, true);
    }

    imgs_delin[n_angios-1].Cargar(rutaground_input, false);

    // Mostrar la imagen en un renderizador
    mis_renderers.push_back(vtkSmartPointer<vtkRenderer>::New());
    mostrarImagen(imgs_base[n_angios-1].base, mis_renderers[n_angios-1]);
}





/*  Metodo: agregarPosicion

    Funcion: Define la posicion en que se encuentra alguna imagen base.
*/
void RECONS3D::agregarPosicion(const double RAO_LAO, const double CAU_CRA, const double Distance_source_to_patient, const double Distance_source_to_detector){

}




/*  Metodo: segmentarImagen()

    Funcion: Multiplica la imagen base por el filtro para extraer las intensidades de los pixeles segmentados
*/
void RECONS3D::segmentarImagenBase(){
    imgs_segment[0] = imgs_base[0];

    mostrarImagen(imgs_base[0].mask, mis_renderers[0]);

    renderizar( mis_renderers[0]);


    FILTROS filtro_gabor;
    filtro_gabor.setFiltro(FILTROS::SS_GABOR);
    filtro_gabor.setFitness(FILTROS::ROC);
    filtro_gabor.setEvoMet(FILTROS::EDA_BUMDA, 50, 30);
    filtro_gabor.setInput(imgs_base[0], imgs_delin[0]);
    filtro_gabor.setOutput(imgs_segment[0]);

//    filtro_gabor.setPar(FILTROS::PAR_K, 12);
//    filtro_gabor.setPar(FILTROS::PAR_DELTA, 1e-4);
//    filtro_gabor.setPar(FILTROS::PAR_SIGMA, 2.83);
//    filtro_gabor.setPar(FILTROS::PAR_L, 13);
//    filtro_gabor.setPar(FILTROS::PAR_T, 15);

//    filtro_gabor.filtrar();
    //using namespace std;
//    TIMERS;
//    GETTIME_INI;

    // Parametros a optimizar:
//    filtro_gabor.setLim(FILTROS::PAR_L, 8.0, 15.0, 0.0001);
//    filtro_gabor.setLim(FILTROS::PAR_T, 8.0, 15.0, 0.0001);
//    filtro_gabor.setLim(FILTROS::PAR_SIGMA, 1.0, 5.0, 0.0001);

    // Parametros fijos:
    filtro_gabor.setPar(FILTROS::PAR_T, 15);
    filtro_gabor.setPar(FILTROS::PAR_L, 2.65);
    filtro_gabor.setPar(FILTROS::PAR_K, 180);
    filtro_gabor.setPar(FILTROS::PAR_DELTA, 1e-4);

//    filtro_gabor.setPar();

//    GETTIME_FIN;
//    FILTROS::INDIV elite = filtro_gabor.getPars();
//    cout << "Tiempo: " << DIFTIME << " segundos" << endl << "Mejores parametros: T = " << elite.vars[0] << ", L = " << elite.vars[1] << ", K = " << elite.vars[2] << ", sigma = " << elite.vars[3] << ", delta = " << elite.vars[4] << endl;

    filtro_gabor.filtrar();
    imgs_segment[0].umbralizar();
    mostrarImagen(imgs_segment[0].base, mis_renderers[0]);
    renderizar(mis_renderers[0]);

}



/*  Metodo: skeletonize()

    Funcion: Obtiene el esqueleto de la imagen y muestra los puntos de interes.
*/
void RECONS3D::skeletonize(){
    imgs_delin[0].skeletonization();
    int n_caracts = imgs_delin[0].n_caracts;

    for( int c = 0; c < n_caracts; c++ ){
        double color[3];
        switch( imgs_delin[0].pix_caract[c].pix_tipo){
            case IMGVTK::PIX_END:
                color[0] = 1.0;
                color[1] = 0.0;
                color[2] = 0.0;
                break;
            case IMGVTK::PIX_BRANCH:
                color[0] = 0.0;
                color[1] = 1.0;
                color[2] = 0.0;
                break;
            case IMGVTK::PIX_CROSS:
                color[0] = 0.0;
                color[1] = 0.0;
                color[2] = 1.0;
                break;
        }

        //agregarEsfera( imgs_delin[0].pix_caract[c].x, imgs_delin[0].pix_caract[c].y, 0.0, 1.5, color, 0 );
    }

//    const int rens = img_delin.rens;
//    const int cols = img_delin.cols;
//    double color_skl[] = {1.0, 1.0, 1.0};

//    unsigned char *skl_ptr = img_delin.skl_ptr;

//    for( int y = 0; y < rens; y++){
//        for( int x = 0; x < cols; x++){
//            if( skl_ptr[(x+1)+(y+1)*(cols+1)] > 0 ){
//                agregarEsfera( x, y, 0.0, 0.5, color_skl );
//            }
//        }
//    }

    mostrarImagen( imgs_delin[0].skeleton, mis_renderers[0] );
    renderizar(mis_renderers[0]);
}



// C O N S T R U C T O R E S    /   D E S T R U C T O R E S
/*  Constructor ()
    Funcion: Constructor por default.
*/
RECONS3D::RECONS3D(){
    esDICOM.push_back(false);
    mis_renderers.push_back( vtkSmartPointer<vtkRenderer>::New() );
    double color[] = {1.0, 1.0, 1.0};
    agregarEsfera(0.0, 0.0, 0.0, 5.0, color, renderer_global);
    renderizar(renderer_global);

    n_angios = 0;
}


/*  Constructor (char **rutas_input, char **rutasground_input , char **rutasmask_input, const int n_imgs)
    Funcion: Recibe las rutas de las imagenes usadas como conjunto de entrenamiento.
*/
RECONS3D::RECONS3D(char **rutasbase_input, char **rutasground_input, const int n_imgs){
    n_angios = 0;
    agregarInput(rutasbase_input, rutasground_input, n_imgs);

    // Preparar el renderer Global:
    renderer_global = vtkSmartPointer<vtkRenderer>::New();
    double color[] = {1.0, 1.0, 1.0};
    agregarEsfera(0.0, 0.0, 0.0, 5.0, color, renderer_global);
    renderizar(renderer_global);
}


/*  Constructor ( const char *rutas_input )
    Funcion: Recibe las rutas de las imagenes usadas como conjunto de entrenamiento.
*/
RECONS3D::RECONS3D(const char *rutabase_input, const char *rutaground_input, const int nivel){
    n_angios = 0;
    agregarInput(rutabase_input, rutaground_input, nivel);

    // Preparar el renderer Global:
    renderer_global = vtkSmartPointer<vtkRenderer>::New();
    double color[] = {1.0, 1.0, 1.0};
    agregarEsfera(0.0, 0.0, 0.0, 5.0, color, renderer_global);
    renderizar(renderer_global);
}



/*  Destructor
    Funcion: Libera la memoria utilizada para almacenar las imagenes.
*/
RECONS3D::~RECONS3D(){
    if( n_angios ){
        delete [] imgs_base;
        delete [] imgs_delin;
        delete [] imgs_segment;
    }
}

//-------------------------------------------------------------------------------------------------- PUBLIC----- ^

// C L A S E: RECONS3D  ------------------------------------------------------------------------ ^
