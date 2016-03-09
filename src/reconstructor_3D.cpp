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



void RECONS3D::agregarEjes(vtkSmartPointer<vtkRenderer> mi_renderer){

    // iamgen temporal para colocar los vectores:
    vtkSmartPointer<vtkImageData> imageX = vtkSmartPointer<vtkImageData>::New();
    imageX->SetDimensions(1,1,1);

    vtkSmartPointer<vtkImageData> imageY= vtkSmartPointer<vtkImageData>::New();
    imageY->SetDimensions(1,1,1);

    vtkSmartPointer<vtkImageData> imageZ = vtkSmartPointer<vtkImageData>::New();
    imageZ->SetDimensions(1,1,1);

  #if VTK_MAJOR_VERSION <= 5
    imageX->SetNumberOfScalarComponents(3);
    imageX->SetScalarTypeToFloat();
    imageX->AllocateScalars();

    imageY->SetNumberOfScalarComponents(3);
    imageY->SetScalarTypeToFloat();
    imageY->AllocateScalars();

    imageZ->SetNumberOfScalarComponents(3);
    imageZ->SetScalarTypeToFloat();
    imageZ->AllocateScalars();
  #else
    imageX->AllocateScalars(VTK_FLOAT,3);
    imageY->AllocateScalars(VTK_FLOAT,3);
    imageZ->AllocateScalars(VTK_FLOAT,3);
  #endif

    {
    float* pixelX = static_cast<float*>(imageX->GetScalarPointer(0,0,0));
    pixelX[0] = 700.0;
    pixelX[1] = 0.0;
    pixelX[2] = 0.0;
    }

    {
    float* pixelY = static_cast<float*>(imageY->GetScalarPointer(0,0,0));
    pixelY[0] = 0.0;
    pixelY[1] = 700.0;
    pixelY[2] = 0.0;
    }

    {
    float* pixelZ = static_cast<float*>(imageZ->GetScalarPointer(0,0,0));
    pixelZ[0] = 0.0;
    pixelZ[1] = 0.0;
    pixelZ[2] = 700.0;
    }

    imageX->GetPointData()->SetActiveVectors("ImageScalars");
    imageY->GetPointData()->SetActiveVectors("ImageScalars");
    imageZ->GetPointData()->SetActiveVectors("ImageScalars");

    vtkSmartPointer<vtkArrowSource> arrowX = vtkSmartPointer<vtkArrowSource>::New();
    arrowX->Update();

    vtkSmartPointer<vtkArrowSource> arrowY = vtkSmartPointer<vtkArrowSource>::New();
    arrowY->Update();

    vtkSmartPointer<vtkArrowSource> arrowZ = vtkSmartPointer<vtkArrowSource>::New();
    arrowZ->Update();



    vtkSmartPointer<vtkGlyph3D> glyphX = vtkSmartPointer<vtkGlyph3D>::New();
    vtkSmartPointer<vtkGlyph3D> glyphY = vtkSmartPointer<vtkGlyph3D>::New();
    vtkSmartPointer<vtkGlyph3D> glyphZ = vtkSmartPointer<vtkGlyph3D>::New();



    glyphX->SetSourceConnection(arrowX->GetOutputPort());
    glyphX->OrientOn();
    glyphX->SetVectorModeToUseVector();

    glyphY->SetSourceConnection(arrowY->GetOutputPort());
    glyphY->OrientOn();
    glyphY->SetVectorModeToUseVector();

    glyphZ->SetSourceConnection(arrowZ->GetOutputPort());
    glyphZ->OrientOn();
    glyphZ->SetVectorModeToUseVector();

#if VTK_MAJOR_VERSION <= 5
    glyphX->SetInputConnection(imageX->GetProducerPort());
    glyphY->SetInputConnection(imageY->GetProducerPort());
    glyphZ->SetInputConnection(imageZ->GetProducerPort());
#else
    glyphX->SetInputData(imageX);
    glyphY->SetInputData(imageY);
    glyphZ->SetInputData(imageZ);
#endif

    glyphX->Update();
    glyphY->Update();
    glyphZ->Update();


    // Create actors
    vtkSmartPointer<vtkImageSliceMapper> mapperX = vtkSmartPointer<vtkImageSliceMapper>::New();
    vtkSmartPointer<vtkImageSliceMapper> mapperY = vtkSmartPointer<vtkImageSliceMapper>::New();
    vtkSmartPointer<vtkImageSliceMapper> mapperZ = vtkSmartPointer<vtkImageSliceMapper>::New();

#if VTK_MAJOR_VERSION <= 5
    mapperX->SetInputConnection(imageX->GetProducerPort());
    mapperY->SetInputConnection(imageY->GetProducerPort());
    mapperZ->SetInputConnection(imageZ->GetProducerPort());
#else
    mapperX->SetInputData(imageX);
    mapperY->SetInputData(imageY);
    mapperZ->SetInputData(imageZ);
#endif

    vtkSmartPointer<vtkImageSlice> sliceX = vtkSmartPointer<vtkImageSlice>::New();
    sliceX->SetMapper(mapperX);

    vtkSmartPointer<vtkImageSlice> sliceY = vtkSmartPointer<vtkImageSlice>::New();
    sliceY->SetMapper(mapperY);

    vtkSmartPointer<vtkImageSlice> sliceZ = vtkSmartPointer<vtkImageSlice>::New();
    sliceZ->SetMapper(mapperZ);



    vtkSmartPointer<vtkPolyDataMapper> ejeXmapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    ejeXmapper->SetInputConnection(glyphX->GetOutputPort());

    vtkSmartPointer<vtkPolyDataMapper> ejeYmapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    ejeYmapper->SetInputConnection(glyphY->GetOutputPort());

    vtkSmartPointer<vtkPolyDataMapper> ejeZmapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    ejeZmapper->SetInputConnection(glyphZ->GetOutputPort());

    vtkSmartPointer<vtkActor> actorX = vtkSmartPointer<vtkActor>::New();
    vtkSmartPointer<vtkActor> actorY = vtkSmartPointer<vtkActor>::New();
    vtkSmartPointer<vtkActor> actorZ = vtkSmartPointer<vtkActor>::New();

    actorX->SetMapper(ejeXmapper);
    actorY->SetMapper(ejeYmapper);
    actorZ->SetMapper(ejeZmapper);

//    mi_renderer->AddActor(actorX);
    mi_renderer->AddActor(actorY);
//    mi_renderer->AddActor(actorZ);
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

    // Agregar el detector y la fuente en posociones por defecto:
    detector.push_back( posicionDefecto( imgs_base[n_angios-1].cols, imgs_base[n_angios-1].rens, imgs_base[n_angios-1].rens/2 ) );
    fuente.push_back( posicionDefecto( imgs_base[n_angios-1].cols, imgs_base[n_angios-1].rens, -imgs_base[n_angios-1].rens/2 ) );
}



/*  Metodo: posicionDefecto

    Funcion: Retorna la posicion por defecto del detector y fuente:
*/
RECONS3D::POS RECONS3D::posicionDefecto( const double ancho, const double alto, const double punta){
    POS defecto;
    /* Base 1
     *      X-------0
     *      |   0   |
     *      0-------0
    */
    defecto.puntos[0][0] = -ancho/2; defecto.puntos[0][1] =  alto/2; defecto.puntos[0][2] = 0.0;

    /* Base 2
     *      0-------X
     *      |   0   |
     *      0-------0
    */
    defecto.puntos[1][0] =  ancho/2; defecto.puntos[1][1] =  alto/2; defecto.puntos[1][2] = 0.0;

    /* Base 3
     *      0-------0
     *      |   0   |
     *      0-------X
    */
    defecto.puntos[2][0] =  ancho/2; defecto.puntos[2][1] = -alto/2; defecto.puntos[2][2] = 0.0;

    /* Base 4
     *      0-------0
     *      |   0   |
     *      X-------0
    */
    defecto.puntos[3][0] = -ancho/2; defecto.puntos[3][1] = -alto/2; defecto.puntos[3][2] = 0.0;

    /* Altura
     *      0-------0
     *      |   X   |
     *      0-------0
    */
    defecto.puntos[4][0] =  0.0; defecto.puntos[4][1] =  0.0; defecto.puntos[4][2] =  punta;

    return defecto;
}




/*  Metodo: agregarPosicion

    Funcion: Define la posicion en que se encuentra alguna imagen base.
*/
void RECONS3D::moverPosicion(const int angio_ID, const double RAO_LAO, const double CAU_CRA, const double Distance_source_to_patient, const double Distance_source_to_detector){
    /// Mover el detector y fuente a las posiciones definidas:
    // Mover el detector a la posicion indicada como SID - SOD:
    const double Distance_patient_to_detector = Distance_source_to_detector - Distance_source_to_patient;

    DEB_MSG("SID: " << Distance_source_to_patient << ", SOD: " << Distance_source_to_detector << ", DDP:" << Distance_patient_to_detector);

    POS det_pos = detector[angio_ID];
    POS fnt_pos = fuente[angio_ID];

    for( int i = 0; i < 5; i++){
        det_pos.puntos[i][2] += Distance_patient_to_detector;
        fnt_pos.puntos[i][2] -= Distance_source_to_patient;
    }

    // Realizar los cambios en los puntos del detector y fuente:
    //// Generar una malla no estructurada:
    vtkSmartPointer<vtkUnstructuredGrid> det_ugrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtkSmartPointer<vtkUnstructuredGrid> fnt_ugrid = vtkSmartPointer<vtkUnstructuredGrid>::New();

    // Definir los puntos del detector y puente
    vtkSmartPointer<vtkPoints> det_points = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkPoints> fnt_points = vtkSmartPointer<vtkPoints>::New();

    for(int i = 0; i < 5; i++){
        det_points->InsertNextPoint( det_pos.puntos[i] );
        fnt_points->InsertNextPoint( fnt_pos.puntos[i] );
    }

    //// Definir la forma del detector y fuente:
    vtkSmartPointer<vtkPyramid> det_pyr = vtkSmartPointer<vtkPyramid>::New();
    vtkSmartPointer<vtkPyramid> fnt_pyr = vtkSmartPointer<vtkPyramid>::New();
    for(int i = 0; i < 5; ++i){
        det_pyr->GetPointIds()->SetId(i, i);
        fnt_pyr->GetPointIds()->SetId(i, i);
    }

    det_ugrid->SetPoints(det_points);
    det_ugrid->InsertNextCell(det_pyr->GetCellType(),det_pyr->GetPointIds());

    fnt_ugrid->SetPoints(fnt_points);
    fnt_ugrid->InsertNextCell(fnt_pyr->GetCellType(),fnt_pyr->GetPointIds());

    vtkSmartPointer<vtkDataSetMapper> det_mapper = vtkSmartPointer<vtkDataSetMapper>::New();
    vtkSmartPointer<vtkActor> det_actor = vtkSmartPointer<vtkActor>::New();

    vtkSmartPointer<vtkDataSetMapper> fnt_mapper = vtkSmartPointer<vtkDataSetMapper>::New();
    vtkSmartPointer<vtkActor> fnt_actor = vtkSmartPointer<vtkActor>::New();

#if VTK_MAJOR_VERSION <= 5
    det_mapper->SetInputConnection(det_ugrid->GetProducerPort());
    fnt_mapper->SetInputConnection(fnt_ugrid->GetProducerPort());
#else
    det_mapper->SetInputData(det_ugrid);
    fnt_mapper->SetInputData(fnt_ugrid);
#endif

    det_actor->SetMapper(det_mapper);
    fnt_actor->SetMapper(fnt_mapper);

    double det_col[] = {0.5, 0.0, 0.5};// Morado
    double fnt_col[] = {1.0, 0.546875, 0.0};// Naranja

    det_actor->GetProperty()->SetColor(det_col);
    fnt_actor->GetProperty()->SetColor(fnt_col);

    renderer_global->AddActor(det_actor);
    renderer_global->AddActor(fnt_actor);

    renderizar(renderer_global);
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
    //agregarEsfera(0.0, 0.0, 0.0, 256.0, color, renderer_global);
    agregarEjes(renderer_global);
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
    //agregarEsfera(0.0, 0.0, 0.0, 256.0, color, renderer_global);
    agregarEjes(renderer_global);
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
    //agregarEsfera(0.0, 0.0, 0.0, 256.0, color, renderer_global);
    agregarEjes(renderer_global);
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
