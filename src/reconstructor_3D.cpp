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
void RECONS3D::renderizar(){

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
void RECONS3D::mostrarImagen( vtkSmartPointer<vtkImageData> &img_src, const int nivel ){

    int dimensions[3];
    img_src->GetDimensions(dimensions);
    const int mis_cols = dimensions[0];
    const int mis_rens = dimensions[1];
    const int mis_rens_cols = mis_cols*mis_rens;


    // Elegir el nivel que se quiere mostar:
    vtkSmartPointer<vtkImageData> img_temp = vtkSmartPointer<vtkImageData>::New();

    // Alojar memoria para la imagen:
    img_temp->SetExtent(0, mis_cols-1, 0, mis_rens-1, nivel, nivel);
    img_temp->AllocateScalars( VTK_UNSIGNED_CHAR, 1);
    img_temp->SetOrigin(0.0, 0.0, (double)nivel);
    img_temp->SetSpacing(1.0, 1.0, 1.0);

    unsigned char *img_ptr = static_cast<unsigned char*>(img_src->GetScalarPointer(0,0,nivel));
    unsigned char *temp_ptr = static_cast<unsigned char*>(img_temp->GetScalarPointer(0,0,nivel));
    memcpy(temp_ptr, img_ptr, mis_rens_cols*sizeof(unsigned char));

    // Crear un actor temporal
    vtkSmartPointer<vtkImageActor> actor = vtkSmartPointer<vtkImageActor>::New();

    DEB_MSG("Actor input antes: " << actor->GetInput() << ", " << img_temp);
    #if VTK_MAJOR_VERSION <= 5
        actor->SetInput(img_temp);
    #else
        actor->SetInputData(img_temp);
    #endif
    DEB_MSG("Actor input despues: " << actor->GetInput() << ", " << img_temp);
    actor->GetInput()->GetDimensions(dimensions);

#ifndef NDEBUG
    double bounds[6];
    actor->GetDisplayBounds(bounds);
    DEB_MSG("Bounds: {x = [" << bounds[0] << "," << bounds[1] << "]; y = [" << bounds[2] << "," << bounds[3] << "]; z = [" << bounds[4] << "," << bounds[5] << "]}" );
    DEB_MSG("Slice number: [" << actor->GetWholeZMin() << "," << actor->GetWholeZMax() << "]");
#endif

    // Agregar el actor de la imagen al renderizador.
    mi_renderer->AddActor(actor);
    mi_renderer->ResetCamera();
}




/*  Metodo: agregarEsfera

    Funcion: Agrega una esfera al renderizador.
*/
void RECONS3D::agregarEsfera( const double x, const double y, const double z, const double radio, double color[3] ){
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
/*  Metodo: segmentarImagen()

    Funcion: Multiplica la imagen base por el filtro para extraer las intensidades de los pixeles segmentados
*/
void RECONS3D::segmentarImagenBase(){
    img_segment = img_base;

    for( int i = 0; i < 57; i++){
        mostrarImagen(img_base.mask, i);
    }

    renderizar();

    /*
    FILTROS filtro_gabor;
    filtro_gabor.setFiltro(FILTROS::GMF);
    filtro_gabor.setFitness(FILTROS::ROC);
    filtro_gabor.setEvoMet(FILTROS::EDA_BUMDA, 50, 30);
    filtro_gabor.setInput(img_base, img_delin);
    filtro_gabor.setOutput(img_segment);
    */
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
    /*
    // Parametros fijos:
    filtro_gabor.setPar(FILTROS::PAR_T, 15);
    filtro_gabor.setPar(FILTROS::PAR_L, 11);
    filtro_gabor.setPar(FILTROS::PAR_SIGMA, 2.82);
    filtro_gabor.setPar(FILTROS::PAR_K, 12);
    filtro_gabor.setPar(FILTROS::PAR_DELTA, 1e-4);
    */
//    filtro_gabor.setPar();

//    GETTIME_FIN;
//    FILTROS::INDIV elite = filtro_gabor.getPars();
//    cout << "Tiempo: " << DIFTIME << " segundos" << endl << "Mejores parametros: T = " << elite.vars[0] << ", L = " << elite.vars[1] << ", K = " << elite.vars[2] << ", sigma = " << elite.vars[3] << ", delta = " << elite.vars[4] << endl;
    /*
    filtro_gabor.filtrar();
    img_segment.umbralizar();
    mostrarImagen(img_segment.base);
    renderizar();
    */
}



/*  Metodo: skeletonize()

    Funcion: Obtiene el esqueleto de la imagen y muestra los puntos de interes.
*/
void RECONS3D::skeletonize(){
    img_delin.skeletonization();
    int n_caracts = img_delin.n_caracts;

    for( int c = 0; c < n_caracts; c++ ){
        double color[3];
        switch( img_delin.pix_caract[c].pix_tipo){
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

        agregarEsfera( img_delin.pix_caract[c].x, img_delin.pix_caract[c].y, 0.0, 1.5, color );
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

    mostrarImagen( img_delin.skeleton, 0 );
    renderizar();
}



// C O N S T R U C T O R E S    /   D E S T R U C T O R E S
/*  Constructor (char **rutas_input, char **rutasground_input , char **rutasmask_input, const int n_imgs)
    Funcion: Recibe las rutas de las imagenes usadas como conjunto de entrenamiento.
*/
RECONS3D::RECONS3D(char **rutasbase_input, char **rutasground_input, const int n_imgs){
    // Cargar la imagen base como la concatenacion de todas las imagenes base:
    img_base.Cargar(rutasbase_input, n_imgs, true);

    // Cargar la imagen delineada como la concatenacion de todas las imagenes delineadas:
    img_delin.Cargar(rutasground_input, n_imgs, false);

    mi_renderer = vtkSmartPointer<vtkRenderer>::New();
}


/*  Constructor ( const char *rutas_input )
    Funcion: Recibe las rutas de las imagenes usadas como conjunto de entrenamiento.
*/
RECONS3D::RECONS3D(const char *rutabase_input, const char *rutaground_input){
    img_base.Cargar( rutabase_input, true);
    img_delin.Cargar( rutaground_input, false);
    mi_renderer = vtkSmartPointer<vtkRenderer>::New();

    if( img_base.esDICOM ){
        gdcm::ImageReader DICOMreader;
        DEB_MSG("ruta base: " << rutabase_input);
        DICOMreader.SetFileName( rutabase_input );
        DICOMreader.Read();
        gdcm::File &file = DICOMreader.GetFile();
        gdcm::FileMetaInformation &fmi = file.GetHeader();
        std::stringstream strm;
        strm.str("");
        const gdcm::Image &gimage = DICOMreader.GetImage();
        DEB_MSG(" Buffer length: " << gimage.GetBufferLength());
        if( fmi.FindDataElement(gdcm::Tag (0x0002, 0x0002)) ){
            fmi.GetDataElement( gdcm::Tag (0x0002, 0x0002) ).GetValue().Print(strm);
            DEB_MSG("TAG[" << (0x0002) << "," << (0x0002) << "]: " << strm.str() );
        }else{
            DEB_MSG("No funciona . . .");
        }
    }
}

//-------------------------------------------------------------------------------------------------- PUBLIC----- ^

// C L A S E: RECONS3D  ------------------------------------------------------------------------ ^
