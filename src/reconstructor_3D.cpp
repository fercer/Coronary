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

    vtkSmartPointer<vtkPolyData> ejesPolyData = vtkSmartPointer<vtkPolyData>::New();

    // Posicion de los puntos en el espacio:
    double origen[3] = { 0.0, 0.0, 0.0 };
    double X[3] = { 2096.0, 0.0, 0.0 };
    double Y[3] = { 0.0, 2096.0, 0.0 };
    double Z[3] = { 0.0, 0.0, 2096.0 };

    // Create a vtkPoints container and store the points in it
    vtkSmartPointer<vtkPoints> ejesPts = vtkSmartPointer<vtkPoints>::New();
    ejesPts->InsertNextPoint(origen);
    ejesPts->InsertNextPoint(X);
    ejesPts->InsertNextPoint(Y);
    ejesPts->InsertNextPoint(Z);

    ejesPolyData->SetPoints(ejesPts);

    // Crear las lineas de los ejes:
    vtkSmartPointer<vtkLine> ejeX = vtkSmartPointer<vtkLine>::New();
    ejeX->GetPointIds()->SetId(0, 0);
    ejeX->GetPointIds()->SetId(1, 1);

    vtkSmartPointer<vtkLine> ejeY = vtkSmartPointer<vtkLine>::New();
    ejeY->GetPointIds()->SetId(0, 0);
    ejeY->GetPointIds()->SetId(1, 2);

    vtkSmartPointer<vtkLine> ejeZ = vtkSmartPointer<vtkLine>::New();
    ejeZ->GetPointIds()->SetId(0, 0);
    ejeZ->GetPointIds()->SetId(1, 3);

    vtkSmartPointer<vtkCellArray> ejes = vtkSmartPointer<vtkCellArray>::New();
    ejes->InsertNextCell(ejeX);
    ejes->InsertNextCell(ejeY);
    ejes->InsertNextCell(ejeZ);

    ejesPolyData->SetLines(ejes);

    unsigned char colX[3] = { 255, 0, 0 };
    unsigned char colY[3] = { 0, 255, 0 };
    unsigned char colZ[3] = { 0, 0, 255 };

    vtkSmartPointer<vtkUnsignedCharArray> colores = vtkSmartPointer<vtkUnsignedCharArray>::New();
    colores->SetNumberOfComponents(3);
    colores->InsertNextTupleValue(colX);
    colores->InsertNextTupleValue(colY);
    colores->InsertNextTupleValue(colZ);
    ejesPolyData->GetCellData()->SetScalars(colores);

    // Setup the visualization pipeline
    vtkSmartPointer<vtkPolyDataMapper> ejesMapper = vtkSmartPointer<vtkPolyDataMapper>::New();

#if VTK_MAJOR_VERSION <= 5
    ejesMapper->SetInput(ejesPolyData);
#else
    ejesMapper->SetInputData(ejesPolyData);
#endif

    vtkSmartPointer<vtkActor> ejesActor = vtkSmartPointer<vtkActor>::New();
    ejesActor->SetMapper(ejesMapper);

    mi_renderer->AddActor( ejesActor );
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

    imgs_base[n_angios-1].Cargar(rutabase_input, true, nivel);
    imgs_delin[n_angios-1].Cargar(rutaground_input, false, 0);

    // Mostrar la imagen en un renderizador
    mis_renderers.push_back(vtkSmartPointer<vtkRenderer>::New());
    mostrarImagen(imgs_base[n_angios-1].base, mis_renderers[n_angios-1]);

    // Agregar el detector y la fuente en posociones por defecto:
    detector.push_back( posicionDefecto( imgs_base[n_angios-1].cols, imgs_base[n_angios-1].rens, imgs_base[n_angios-1].rens/2 ) );

    // Mover el detector a su posicion definida por el archivo DICOM:
    moverPosicion(n_angios-1);
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
void RECONS3D::moverPosicion(const int angio_ID){
    /// Mover el detector y fuente a las posiciones definidas:
    POS det_pos = detector[angio_ID];

    // Mover los puntos segun SID y SOD:
    for( int i = 0; i < 5; i++){
        det_pos.puntos[i][2] += imgs_base[angio_ID].DDP;
    }

    // Rotar los puntos segun LAO/RAO y CAU/CRA:
    const double crl = cos(imgs_base[angio_ID].LAORAO/180.0 * PI);
    const double srl = sin(imgs_base[angio_ID].LAORAO/180.0 * PI);
    const double ccc = cos(imgs_base[angio_ID].CRACAU/180.0 * PI);
    const double scc = sin(imgs_base[angio_ID].CRACAU/180.0 * PI);

    //// Rotacion usando el eje x como base:
    for( int i = 0; i < 5; i++){
        const double old_y = det_pos.puntos[i][1];
        det_pos.puntos[i][1] = crl*old_y - srl*det_pos.puntos[i][2];
        det_pos.puntos[i][2] = srl*old_y + crl*det_pos.puntos[i][2];
    }

    //// Rotacion usando el eje y como base:
    for( int i = 0; i < 5; i++){
        const double old_x = det_pos.puntos[i][0];
        det_pos.puntos[i][0] = ccc*old_x - scc*det_pos.puntos[i][2];
        det_pos.puntos[i][2] = scc*old_x + ccc*det_pos.puntos[i][2];
    }

    // Realizar los cambios en los puntos del detector y fuente:
    //// Generar una malla no estructurada:
    vtkSmartPointer<vtkUnstructuredGrid> det_ugrid = vtkSmartPointer<vtkUnstructuredGrid>::New();

    // Definir los puntos del detector y puente
    vtkSmartPointer<vtkPoints> det_points = vtkSmartPointer<vtkPoints>::New();

    for(int i = 0; i < 5; i++){
        det_points->InsertNextPoint( det_pos.puntos[i] );
    }

    //// Definir la forma del detector y fuente:
    vtkSmartPointer<vtkPyramid> det_pyr = vtkSmartPointer<vtkPyramid>::New();
    for(int i = 0; i < 5; ++i){
        det_pyr->GetPointIds()->SetId(i, i);
    }

    det_ugrid->SetPoints(det_points);
    det_ugrid->InsertNextCell(det_pyr->GetCellType(),det_pyr->GetPointIds());

    vtkSmartPointer<vtkDataSetMapper> det_mapper = vtkSmartPointer<vtkDataSetMapper>::New();
    vtkSmartPointer<vtkActor> det_actor = vtkSmartPointer<vtkActor>::New();

#if VTK_MAJOR_VERSION <= 5
    det_mapper->SetInputConnection(det_ugrid->GetProducerPort());
#else
    det_mapper->SetInputData(det_ugrid);
#endif

    det_actor->SetMapper(det_mapper);

    double det_col[] = {0.5, 0.0, 0.5};// Morado

    det_actor->GetProperty()->SetColor(det_col);

    renderer_global->AddActor(det_actor);

    renderizar(renderer_global);
}




/*  Metodo: segmentarImagen()

    Funcion: Multiplica la imagen base por el filtro para extraer las intensidades de los pixeles segmentados
*/
void RECONS3D::segmentarImagenBase(){
    imgs_segment[0] = imgs_base[0];

    mostrarImagen(imgs_base[0].base, mis_renderers[0]);

    renderizar( mis_renderers[0]);

/*
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
*/
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
    mis_renderers.push_back( vtkSmartPointer<vtkRenderer>::New() );
    double color[] = {1.0, 1.0, 1.0};
    agregarEsfera(0.0, 0.0, 0.0, 100.0, color, renderer_global);
    agregarEjes(renderer_global);

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
    agregarEsfera(0.0, 0.0, 0.0, 100.0, color, renderer_global);
    agregarEjes(renderer_global);
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
    agregarEsfera(0.0, 0.0, 0.0, 100.0, color, renderer_global);
    agregarEjes(renderer_global);
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
