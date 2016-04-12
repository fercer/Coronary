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




/*  Metodo: mallaPuntos

    Funcion: Genera una malla de puntos para representar la imagen en un renderizador global.
*/
void RECONS3D::mallarPuntos( const int angio_ID ){


    // Rotar los puntos segun LAO/RAO y CAU/CRA:
    const double crl = cos(imgs_base[angio_ID].LAORAO/180.0 * PI);
    const double srl = sin(imgs_base[angio_ID].LAORAO/180.0 * PI);
    const double ccc = cos(imgs_base[angio_ID].CRACAU/180.0 * PI);
    const double scc = sin(imgs_base[angio_ID].CRACAU/180.0 * PI);

    const int mis_cols = imgs_base[angio_ID].cols+1;
    const int mis_rens = imgs_base[angio_ID].rens+1;

   // puntos[angio_ID]->SetNumberOfPoints((mis_cols+1)*(mis_rens+1));

    const double mi_pixX = imgs_base[angio_ID].pixX;
    const double mi_pixY = imgs_base[angio_ID].pixY;

    const double mi_DDP = imgs_base[angio_ID].DDP;

    double yy = - mis_rens/2 * mi_pixY;

    for( int y = 0; y < mis_rens; y++){
        double xx = - mis_cols/2 * mi_pixX;
        for(int x = 0; x < mis_cols; x++){

            double xx_3D, yy_3D, zz_3D;

            // Mover los puntos segun indica el SID y SOD:

            //// Rotacion usando el eje 'x' como base:
            yy_3D = crl*yy - srl*mi_DDP;
            zz_3D = srl*yy + crl*mi_DDP;


            //// Rotacion usando el eje 'y' como base:
            xx_3D = ccc*xx - scc*zz_3D;
            zz_3D = scc*xx + ccc*zz_3D;

            puntos[angio_ID]->InsertNextPoint(xx_3D, yy_3D, zz_3D);
            xx += mi_pixX;
        }
        yy += mi_pixY;
    }

}



/*  Metodo: mostrarImagen

    Funcion: Muestra la imagen en una ventana VTK.
*/
void RECONS3D::mostrarImagen( IMGVTK &img_src, IMGVTK::IMG_IDX img_idx, vtkSmartPointer<vtkRenderer> mi_renderer){
    // Crear un actor temporal
    vtkSmartPointer<vtkImageActor> actor = vtkSmartPointer<vtkImageActor>::New();

    int mis_cols = img_src.cols;
    int mis_rens = img_src.rens;

    vtkSmartPointer< vtkImageData> img_ptr = NULL;

    switch(img_idx){
        case IMGVTK::BASE:
            img_ptr = img_src.base;
            break;
        case IMGVTK::MASK:
            img_ptr = img_src.mask;
            break;
        case IMGVTK::SKELETON:
            img_ptr = img_src.skeleton;
            mis_rens+=2;
            mis_cols+=2;
            break;
        case IMGVTK::SEGMENT:
            img_ptr = img_src.segment;
            break;
    }

    const int mis_rens_cols = mis_rens*mis_cols;

    /// Crear una imagen temporal:
    vtkSmartPointer<vtkImageData> img_tmp = vtkSmartPointer<vtkImageData>::New();
    img_tmp->SetExtent(0, mis_cols-1, 0, mis_rens-1, 0, 0);
    img_tmp->AllocateScalars( VTK_UNSIGNED_CHAR, 1);
    img_tmp->SetOrigin(0.0, 0.0, 0.0);
    img_tmp->SetSpacing(1.0, 1.0, 1.0);

    unsigned char *img_tmp_ptr = static_cast<unsigned char*>(img_tmp->GetScalarPointer(0,0,0));
    double *img_src_ptr = static_cast<double*>(img_ptr->GetScalarPointer(0,0,0));

    double max =-1e100;
    double min = 1e100;

    for( int xy = 0; xy < mis_rens_cols; xy++){
        if(*(img_src_ptr + xy) < min){
            min = *(img_src_ptr + xy);
        }
        if(*(img_src_ptr + xy) > max){
            max = *(img_src_ptr + xy);
        }
    }

    const double rango = max - min;
    for( int xy = 0; xy < mis_rens_cols; xy++){
        *(img_tmp_ptr + xy) = (unsigned char)(255.0 * (*(img_src_ptr + xy) - min) / rango);
    }

    actor->SetInputData(img_tmp);

    // Agregar el actor de la imagen al renderizador.
    mi_renderer->AddActor(actor);
    mi_renderer->ResetCamera();
}





/*  Metodo: mostrarImagen

    Funcion: Muestra la imagen en una ventana VTK usando el renderizador global.
*/
void RECONS3D::mostrarImagen( const int angio_ID, IMGVTK::IMG_IDX img_idx){

    int mis_cols = imgs_base[angio_ID].cols;
    int mis_rens = imgs_base[angio_ID].rens;

    double *img_ptr = NULL;

    switch(img_idx){
        case IMGVTK::BASE:
            img_ptr = imgs_base[angio_ID].base_ptr;
            break;
        case IMGVTK::MASK:
            img_ptr = imgs_base[angio_ID].mask_ptr;
            break;
        case IMGVTK::SEGMENT:
            img_ptr = imgs_base[angio_ID].segment_ptr;
            break;
    }

    vtkSmartPointer<vtkCellArray> quads = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkDoubleArray> intensidades = vtkSmartPointer<vtkDoubleArray>::New();

    intensidades->SetNumberOfComponents(1);
    //intensidades->SetNumberOfTuples(mis_cols*mis_rens);
    intensidades->SetName("Intensidad");

    for( int i = 0; i < mis_rens; i++){
        for(int j = 0; j < mis_cols; j++){
            vtkSmartPointer<vtkQuad> quad = vtkSmartPointer<vtkQuad>::New();
            quad->GetPointIds()->SetId(0,( j )+( i )*(mis_cols+1));
            quad->GetPointIds()->SetId(1,(j+1)+( i )*(mis_cols+1));
            quad->GetPointIds()->SetId(2,(j+1)+(i+1)*(mis_cols+1));
            quad->GetPointIds()->SetId(3,( j )+(i+1)*(mis_cols+1));

            quads->InsertNextCell(quad);
            //const double intensidad = img_ptr[i*mis_cols + j];

            //intensidades->SetTuple3(i*mis_cols+j, intensidad, intensidad, intensidad);
            intensidades->InsertNextValue(img_ptr[i*mis_cols + j]);
        }
    }

    vtkSmartPointer<vtkPolyData> rejilla = vtkSmartPointer<vtkPolyData>::New();
    rejilla->SetPoints( puntos[angio_ID] );
    rejilla->SetPolys(quads);

    rejilla->GetCellData()->SetScalars(intensidades);


    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputData(rejilla);

    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);

    // Setup render window, renderer, and interactor
    renderer_global->AddActor(actor);
}




/*  Metodo: agregarEjes

    Funcion: Muestra los ejes en que se encuentra la escena global respecto al paciente.
*/
void RECONS3D::agregarEjes(vtkSmartPointer<vtkRenderer> mi_renderer){

    vtkSmartPointer<vtkPolyData> ejesPolyData = vtkSmartPointer<vtkPolyData>::New();

    // Posicion de los puntos en el espacio:
    double origen[3] = { 0.0, 0.0, 0.0 };
    double X[3] = { 512.0, 0.0, 0.0 };
    double Y[3] = { 0.0, 512.0, 0.0 };
    double Z[3] = { 0.0, 0.0, 512.0 };

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
    ejesMapper->SetInputData(ejesPolyData);

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



/*  Metodo: posicionDefecto

    Funcion: Retorna la posicion por defecto del detector y fuente:
*/
inline RECONS3D::POS RECONS3D::posicionDefecto( const double ancho, const double alto, const double punta){
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




/*  Metodo: moverPosicion

    Funcion: Mueve los puntos de la imagen base hasta la orientacion definida por los angulos (en el archivo DICOM)
*/
void RECONS3D::mostrarDetector(const int angio_ID){

    // Rotar los puntos segun LAO/RAO y CAU/CRA:
    const double crl = cos(imgs_base[angio_ID].LAORAO/180.0 * PI);
    const double srl = sin(imgs_base[angio_ID].LAORAO/180.0 * PI);
    const double ccc = cos(imgs_base[angio_ID].CRACAU/180.0 * PI);
    const double scc = sin(imgs_base[angio_ID].CRACAU/180.0 * PI);

    const int mis_rens = imgs_base[n_angios].rens;
    const int mis_cols = imgs_base[n_angios].cols;

    const double mi_pixX = imgs_base[angio_ID].pixX;
    const double mi_pixY = imgs_base[angio_ID].pixY;

    //------------------------------------------------------------------------------------------------------------------ FIGURA DEL DETECTOR
    /// Mover el detector y fuente a las posiciones definidas:
    POS det_pos = posicionDefecto( mis_cols*mi_pixX, mis_rens*mi_pixY, mis_rens*mi_pixX/2 );

    // Mover los puntos segun SID y SOD:
    for( int i = 0; i < 5; i++){
        det_pos.puntos[i][2] += imgs_base[angio_ID].DDP;
    }

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
    det_mapper->SetInputData(det_ugrid);

    vtkSmartPointer<vtkActor> det_actor = vtkSmartPointer<vtkActor>::New();
    det_actor->SetMapper(det_mapper);

    double det_col[] = {0.5, 0.0, 0.5};// Morado

    det_actor->GetProperty()->SetColor(det_col);

    renderer_global->AddActor(det_actor);
}



//------------------------------------------------------------------------------------------------ PRIVATE ----- ^
//-------------------------------------------------------------------------------------------------- PUBLIC----- v

// M I E M B R O S      P U B L I C O S
/*  Metodo: agregarInput

    Funcion: Define las rutas de las imagenes que son usadas para reconstruir una arteria coronaria.
*/
void RECONS3D::agregarInput(char **rutasbase_input, char **rutasground_input, const int n_imgs){
    imgs_base.push_back(IMGVTK(rutasbase_input, n_imgs, true));
    imgs_delin.push_back(IMGVTK(rutasground_input, n_imgs, false));
    imgs_base.push_back(IMGVTK());

    // Mostrar la imagen en un renderizador
    mis_renderers.push_back(vtkSmartPointer<vtkRenderer>::New());
    mostrarImagen(imgs_base[n_angios], IMGVTK::BASE, mis_renderers[n_angios]);

    n_angios++;
}



/*  Metodo: agregarInput

    Funcion: Define las rutas de las imagenes que son usadas para reconstruir una arteria coronaria.
*/
void RECONS3D::agregarInput(const char *rutabase_input, const char *rutaground_input, const int nivel){

    imgs_base.push_back(IMGVTK(rutabase_input, true, nivel));
    imgs_delin.push_back(IMGVTK(rutaground_input, false, 0));
    imgs_base.push_back(IMGVTK());

    // Mostrar la imagen en un renderizador
    mis_renderers.push_back(vtkSmartPointer<vtkRenderer>::New());
    puntos.push_back(vtkSmartPointer<vtkPoints>::New());


    // Mover el detector a su posicion definida por el archivo DICOM:
    mallarPuntos(n_angios);
    //mostrarDetector(n_angios);

    mostrarImagen(n_angios, IMGVTK::BASE);

    renderizar( renderer_global );


    n_angios++;
}




/*  Metodo: segmentarImagen()

    Funcion: Multiplica la imagen base por el filtro para extraer las intensidades de los pixeles segmentados
*/
void RECONS3D::segmentarImagenBase( const int angio_ID ){
    FILTROS filtro;
    filtro.setFiltro(FILTROS::SS_GABOR);
    filtro.setFitness(FILTROS::ROC);
    filtro.setEvoMet(FILTROS::EDA_BUMDA, 50, 30);

    filtro.setInput(imgs_base[angio_ID]);
    filtro.setInputGround(imgs_delin[angio_ID]);

    // Parametros fijos:
    filtro.setPar(FILTROS::PAR_L, 2.65);
    filtro.setPar(FILTROS::PAR_T, 15);
    filtro.setPar(FILTROS::PAR_K, 180);
    filtro.setPar(FILTROS::PAR_DELTA, 1e-4);
    filtro.filtrar();

#ifndef NDEBUG
    mostrarImagen(imgs_base[angio_ID], IMGVTK::SEGMENT, mis_renderers[angio_ID]);
    renderizar(mis_renderers[angio_ID]);
#endif

    imgs_base[angio_ID].umbralizar(IMGVTK::SEGMENT);

    /// Falta el linking broken vessels ..................

#ifndef NDEBUG
    mostrarImagen(imgs_base[angio_ID], IMGVTK::SEGMENT, mis_renderers[angio_ID]);
    renderizar(mis_renderers[angio_ID]);
#endif
    imgs_base[angio_ID].lengthFilter(imgs_base[angio_ID].cols * 5);

#ifndef NDEBUG
    mostrarImagen(imgs_base[angio_ID], IMGVTK::SEGMENT, mis_renderers[angio_ID]);
    renderizar(mis_renderers[angio_ID]);
#endif
    imgs_base[angio_ID].regionFill();

#ifndef NDEBUG
    mostrarImagen(imgs_base[angio_ID], IMGVTK::SEGMENT, mis_renderers[angio_ID]);
    renderizar(mis_renderers[angio_ID]);
#endif

}



/*  Metodo: skeletonize()

    Funcion: Obtiene el esqueleto de la imagen y muestra los puntos de interes.
*/
void RECONS3D::skeletonize(const int angio_ID){
    imgs_base[angio_ID].skeletonization(IMGVTK::SEGMENT);

    int n_caracts = imgs_base[angio_ID].n_caracts;

    for( int c = 0; c < n_caracts; c++ ){
        double color[3];
        switch( imgs_base[angio_ID].pix_caract[c].pix_tipo){
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

        agregarEsfera( imgs_base[angio_ID].pix_caract[c].x, imgs_base[angio_ID].pix_caract[c].y, 0.0, 1.5, color,  mis_renderers[angio_ID] );
    }


    /// Obtener el esqueleto de la imagen segmentada y graficar esferas en cada punto caracteristico
    //mostrarImagen( imgs_segment[angio_ID].skeleton, mis_renderers[angio_ID] );
    //renderizar( mis_renderers[angio_ID] );
}



// C O N S T R U C T O R E S    /   D E S T R U C T O R E S
/*  Constructor ()
    Funcion: Constructor por default.
*/
RECONS3D::RECONS3D(){
    renderer_global = vtkSmartPointer<vtkRenderer>::New();

    double color[] = {1.0, 1.0, 1.0};
    agregarEsfera(0.0, 0.0, 0.0, 50.0, color, renderer_global);
    agregarEjes(renderer_global);

    n_angios = 0;
}



/*  Destructor
    Funcion: Libera la memoria utilizada para almacenar las imagenes.
*/
RECONS3D::~RECONS3D(){

}

//-------------------------------------------------------------------------------------------------- PUBLIC----- ^

// C L A S E: RECONS3D  ------------------------------------------------------------------------ ^
