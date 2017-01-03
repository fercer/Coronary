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


/************************************************************************************************************
* RECONS3D::PRIVATE                                                                                         *
*                                                                                                           *
* FUNCTION NAME: getIMGCONTPointer                                                                          *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT                  TYPE                      I/O  DESCRIPTION                                      *
* my_img_id                 const IMG_IDX              I   The image identificator to return its pointer.   *
*                                                                                                           *
* RETURNS:                                                                                                  *
* The image container pointer                                                                               *
*                                                                                                           *
************************************************************************************************************/
IMGCONT * RECONS3D::getIMGCONTPointer(const IMG_IDX my_img_idx, const unsigned int angios_ID)
{
	switch (my_img_idx) {
	case IMG_BASE:
		return my_img_features_ptr[angios_ID].my_base_ptr;
	case IMG_BOUNDARIES:
		return my_img_features_ptr[angios_ID].my_thrs_response_boundaries_ptr;
	case IMG_DIST_MAP:
		return my_img_features_ptr[angios_ID].my_thrs_response_map_dists_ptr;
	case IMG_GROUNDTRUTH:
		return my_img_features_ptr[angios_ID].my_groundtruth_ptr;
	case IMG_MASK:
		return my_img_features_ptr[angios_ID].my_mask_ptr;
	case IMG_RESPONSE:
		return my_img_features_ptr[angios_ID].my_response_ptr;
	case IMG_SKELETON:
		return my_img_features_ptr[angios_ID].my_thrs_response_skeleton_ptr;
	case IMG_THRESHOLD:
		return my_img_features_ptr[angios_ID].my_response_threshold_ptr;
	}
}



/*  Metodo: escribirLog

    Funcion: Escribe un mensaje en el log.
*/
void RECONS3D::escribirLog( const char *mensaje ){
    if(fp_log){
        fprintf(fp_log, "%s", mensaje);

#ifdef BUILD_GUI_VERSION
    }else if(mi_txtLog){
        mi_txtLog->append( mensaje );
#endif

    }else{
        std::cout << mensaje;
    }
}



/*  Metodo: barraProgreso

    Funcion: Actualiza la barra de progreso.
*/
void RECONS3D::barraProgreso( const int progress, const int max_progress ){

#ifdef BUILD_GUI_VERSION
    if( mi_pBar ){
        mi_pBar->setMaximum( max_progress );
        mi_pBar->setValue( progress );
        return;
    }
#endif

    /// Limpiar el resto de la linea:
    int max_ancho = 100;
    for( int i = 0; i < max_ancho; i++){
        printf("\r");
    }
    printf(COLOR_BACK_RED "[");
    int progress_completed = (int)((double)max_ancho * (double)progress / (double)max_progress);
    for( int i = 0; i < progress_completed; i++){
        printf(COLOR_BACK_GREEN " ");
    }
    for( int i = progress_completed; i < max_ancho; i++){
        printf(COLOR_BACK_CYAN " ");
    }
    printf(COLOR_BACK_RED "]" COLOR_NORMAL);
    fflush(stdout);
}



/*  Metodo: renderizar()

    Funcion: Renderiza los 'actores' contenidos en 'renderer' en una ventana de VTK.
*/
#ifdef BUILD_VTK_VERSION
void RECONS3D::renderizar( vtkSmartPointer<vtkRenderer> mi_renderer ){
    mi_renderer->ResetCamera();
DEB_MSG("Camara reset");
    // Crear una ventana temporal:
    vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow->AddRenderer(mi_renderer);
DEB_MSG("Ventana generada");
    // Crear interactuador temporal:
    vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();

    // Crear la ventana y mostrar la renderizacion:
    renderWindowInteractor->SetRenderWindow(renderWindow);
    renderWindowInteractor->Initialize();
DEB_MSG("Iterador para la ventana lista");

    renderWindowInteractor->Start();
DEB_MSG("Retornando...");
}
#endif


/*  Metodo: isoCentro

    Funcion: Calcula el isocentro de la imagen 'angios_ID' como la linea que va desde el centro de la imagen al centro de la fuente.
*/
#ifdef BUILD_VTK_VERSION
void RECONS3D::isoCentro( const int angios_ID ){
    // Calcular la ecuacion del plano:
    const double crl = cos((*(imgs_base + angios_ID))->LAORAO/180.0 * Mi_PI);
    const double srl = sin((*(imgs_base + angios_ID))->LAORAO/180.0 * Mi_PI);
    const double ccc = cos((*(imgs_base + angios_ID))->CRACAU/180.0 * Mi_PI);
    const double scc = sin((*(imgs_base + angios_ID))->CRACAU/180.0 * Mi_PI);

    const int mis_cols = (*(imgs_base + angios_ID))->cols;
    const int mis_rens = (*(imgs_base + angios_ID))->rows;

    const double mi_DDP = (*(imgs_base + angios_ID))->DDP;

    const double mi_pixX = (*(imgs_base + angios_ID))->pixX;
    const double mi_pixY = (*(imgs_base + angios_ID))->pixY;

    // Primer punto (0,0, DDP) ------------------------------------------------------------------------------
    double xx = 0.0, yy = 0.0;
    double Px, Py, Pz;

    //// Rotacion usando el eje 'x' como base:
    Py = crl*yy - srl*mi_DDP;
    Pz = srl*yy + crl*mi_DDP;
    //// Rotacion usando el eje 'y' como base:
    Px = ccc*xx - scc*Pz;
    Pz = scc*xx + ccc*Pz;


    // Segundo punto (mis_cols,0, DDP) ------------------------------------------------------------------------------
    xx = mis_cols*mi_pixX/2;
    double Qx, Qy, Qz;

    //// Rotacion usando el eje 'x' como base:
    Qy = crl*yy - srl*mi_DDP;
    Qz = srl*yy + crl*mi_DDP;
    //// Rotacion usando el eje 'y' como base:
    Qx = ccc*xx - scc*Qz;
    Qz = scc*xx + ccc*Qz;



    // Tercer punto (0, mis_rens, DDP) ------------------------------------------------------------------------------
    xx = 0.0;
    yy = mis_rens*mi_pixY/2;
    double Rx, Ry, Rz;

    //// Rotacion usando el eje 'x' como base:
    Ry = crl*yy - srl*mi_DDP;
    Rz = srl*yy + crl*mi_DDP;
    //// Rotacion usando el eje 'y' como base:
    Rx = ccc*xx - scc*Rz;
    Rz = scc*xx + ccc*Rz;


    //// Calcular la normal al plano:
    std::vector< NORCEN >::iterator mi_normal = normal_centros.begin() + angios_ID;
    mi_normal->direccion[0] = (Qy - Py)*(Rz - Pz) - (Qz - Pz)*(Ry - Py);
    mi_normal->direccion[1] = (Qz - Pz)*(Rx - Px) - (Qx - Px)*(Rz - Pz);
    mi_normal->direccion[2] = (Qx - Px)*(Ry - Py) - (Qy - Py)*(Rx - Px);

    /// Definir el centro de la imagen:
    mi_normal->origen[0] = Px;
    mi_normal->origen[1] = Py;
    mi_normal->origen[2] = Pz;
}
#endif


/*  Metodo: mallaPuntos

    Funcion: Genera una malla de puntos para representar la imagen en un renderizador global.
*/

#ifdef BUILD_VTK_VERSION
void RECONS3D::mallarPuntos( const int angios_ID ){
    // Rotar los puntos segun LAO/RAO y CAU/CRA:
    const double crl = cos((*(imgs_base + angios_ID))->LAORAO/180.0 * Mi_PI);
    const double srl = sin((*(imgs_base + angios_ID))->LAORAO/180.0 * Mi_PI);
    const double ccc = cos((*(imgs_base + angios_ID))->CRACAU/180.0 * Mi_PI);
    const double scc = sin((*(imgs_base + angios_ID))->CRACAU/180.0 * Mi_PI);

    const int mis_cols = (*(imgs_base + angios_ID))->cols;
    const int mis_rens = (*(imgs_base + angios_ID))->rows;
    const double mi_pixX = (*(imgs_base + angios_ID))->pixX;
    const double mi_pixY = (*(imgs_base + angios_ID))->pixY;

    const double mi_DDP = (*(imgs_base + angios_ID))->DDP;

    double yy = (1 - mis_rens)*mi_pixY/2.0;

    for( int y = 0; y < mis_rens; y++){
        double xx = (1 - mis_cols)*mi_pixX/2.0;
        for(int x = 0; x < mis_cols; x++){
            double xx_3D, yy_3D, zz_3D;

            // Mover los puntos segun indica el SID y SOD:

            //// Rotacion usando el eje 'x' como base:
            yy_3D = crl*yy - srl*mi_DDP;
            zz_3D = srl*yy + crl*mi_DDP;


            //// Rotacion usando el eje 'y' como base:
            xx_3D = ccc*xx - scc*zz_3D;
            zz_3D = scc*xx + ccc*zz_3D;

            puntos[angios_ID]->InsertNextPoint(xx_3D, yy_3D, zz_3D);

            vtkSmartPointer< vtkVertex > pix = vtkSmartPointer< vtkVertex >::New();
            pix->GetPointIds()->SetId(0, x + y*mis_cols);

            pixeles[angios_ID]->InsertNextCell(pix);
            xx += mi_pixX;
        }
        yy += mi_pixY;
    }
}
#endif



/*  Metodo: mostrarImagen

    Funcion: Muestra la imagen en una ventana VTK.
*/

#ifdef BUILD_VTK_VERSION
void RECONS3D::mostrarImagen(IMGVTK::IMG_IDX img_idx, vtkSmartPointer<vtkRenderer> mi_renderer, const int angios_ID){

    int mis_cols = (*(imgs_base + angios_ID))->cols;
    int mis_rens = (*(imgs_base + angios_ID))->rows;

    vtkSmartPointer< vtkImageData> img_ptr = NULL;

    switch(img_idx){
    case IMGVTK::BASE:
        img_ptr = (*(imgs_base + angios_ID))->base;
        break;
    case IMGVTK::GROUNDTRUTH:
        img_ptr = (*(imgs_base + angios_ID))->ground;
        break;
    case IMGVTK::MASK:
        img_ptr = (*(imgs_base + angios_ID))->mask;
        break;
    case IMGVTK::SKELETON:
        img_ptr = (*(imgs_base + angios_ID))->skeleton;
        mis_rens+=2;
        mis_cols+=2;
        break;
    case IMGVTK::SEGMENT:
        img_ptr = (*(imgs_base + angios_ID))->segment;
        break;
    case IMGVTK::THRESHOLD:
        img_ptr = (*(imgs_base + angios_ID))->threshold;
        break;
    case IMGVTK::MAPDIST:
        img_ptr = (*(imgs_base + angios_ID))->mapa_dist;
        break;
    case IMGVTK::BORDERS:
        img_ptr = (*(imgs_base + angios_ID))->borders;
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

    double max =-MY_INF;
    double min = MY_INF;

    for( int xy = 0; xy < mis_rens_cols; xy++){
        if(*(img_src_ptr + xy) < min){
            min = *(img_src_ptr + xy);
        }
        if(*(img_src_ptr + xy) > max){
            max = *(img_src_ptr + xy);
        }
    }

DEB_MSG("min: " << min << ", max: " << max);
    const double rango = max - min;
    for( int xy = 0; xy < mis_rens_cols; xy++){
        *(img_tmp_ptr + xy) = (unsigned char)(255.0 * (*(img_src_ptr + xy) - min) / rango);
    }

    // Crear un mapper y actor temporales
    vtkSmartPointer<vtkImageMapper> mapper = vtkSmartPointer<vtkImageMapper>::New();
    mapper->SetInputData(img_tmp);
    mapper->SetColorWindow(255);
    mapper->SetColorLevel(127.5);

    vtkSmartPointer<vtkActor2D> actor = vtkSmartPointer<vtkActor2D>::New();
    actor->SetMapper(mapper);

    // Agregar el actor de la imagen al renderizador.
    mi_renderer->Clear();
    mi_renderer->AddActor2D(actor);
}
#endif



/*  Metodo: mostrarImagen

    Funcion: Muestra la imagen en una ventana VTK usando el renderizador global.
*/

#ifdef BUILD_VTK_VERSION
void RECONS3D::mostrarImagen( const int angios_ID, IMGVTK::IMG_IDX img_idx){
    const int mis_cols = (*(imgs_base + angios_ID))->cols;
    const int mis_rens = (*(imgs_base + angios_ID))->rows;

    int offset_y = 0, offset_x = 0;
    double *img_ptr = NULL;

    switch(img_idx){
        case IMGVTK::BASE:
            img_ptr = (*(imgs_base + angios_ID))->base_ptr;
            break;
        case IMGVTK::MASK:
            img_ptr = (*(imgs_base + angios_ID))->mask_ptr;
            break;
        case IMGVTK::SKELETON:
            img_ptr = (*(imgs_base + angios_ID))->skl_ptr;
            offset_x=1;
            offset_y=1;
            break;
        case IMGVTK::SEGMENT:
            img_ptr = (*(imgs_base + angios_ID))->segment_ptr;
            break;
        case IMGVTK::THRESHOLD:
            img_ptr = (*(imgs_base + angios_ID))->threshold_ptr;
            break;
        case IMGVTK::MAPDIST:
            img_ptr = (*(imgs_base + angios_ID))->map_ptr;
            break;
        case IMGVTK::BORDERS:
                img_ptr = (*(imgs_base + angios_ID))->borders_ptr;
                break;
    }

    vtkSmartPointer<vtkUnsignedCharArray> intensidades = vtkSmartPointer<vtkUnsignedCharArray>::New();

    intensidades->SetNumberOfComponents(3);
    intensidades->SetName("Intensidadsdasdsads");
    unsigned char intensidad_xy[3];

    for( int i = 0; i < mis_rens; i++){
        for(int j = 0; j < mis_cols; j++){
            intensidad_xy[0] = (unsigned char)(255.0 * img_ptr[(i+offset_y)*(mis_cols+offset_x*2) + (j+offset_x)]);
            intensidad_xy[1] = (unsigned char)(255.0 * img_ptr[(i+offset_y)*(mis_cols+offset_x*2) + (j+offset_x)]);
            intensidad_xy[2] = (unsigned char)(255.0 * img_ptr[(i+offset_y)*(mis_cols+offset_x*2) + (j+offset_x)]);
            intensidades->InsertNextTupleValue( intensidad_xy );
        }
    }

    vtkSmartPointer<vtkPolyData> rejilla = vtkSmartPointer<vtkPolyData>::New();
    rejilla->SetPoints( puntos[angios_ID] );
    rejilla->SetVerts( pixeles[angios_ID] );
    rejilla->GetCellData()->SetScalars(intensidades);

    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputData(rejilla);

    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    renderer_global->AddActor(actor);

    // Mostrar el isocentro:

    /// Posicionar una esfera en el lugar del isocentro:
    double norma = normal_centros[angios_ID].direccion[0]*normal_centros[angios_ID].direccion[0] +
                   normal_centros[angios_ID].direccion[1]*normal_centros[angios_ID].direccion[1] +
                   normal_centros[angios_ID].direccion[2]*normal_centros[angios_ID].direccion[2];

    norma = sqrt(norma);

    const double t_iso = (*(imgs_base + angios_ID))->DISO / norma;

    NORCEN mi_isocen;
    mi_isocen.origen[0] = normal_centros[angios_ID].origen[0];
    mi_isocen.origen[1] = normal_centros[angios_ID].origen[1];
    mi_isocen.origen[2] = normal_centros[angios_ID].origen[2];
    mi_isocen.direccion[0] = - t_iso * normal_centros[angios_ID].direccion[0];
    mi_isocen.direccion[1] = - t_iso * normal_centros[angios_ID].direccion[1];
    mi_isocen.direccion[2] = - t_iso * normal_centros[angios_ID].direccion[2];

    double color[] = {0.0, 0.0, 1.0};
    agregarEsfera(normal_centros[angios_ID].origen[0], normal_centros[angios_ID].origen[1], normal_centros[angios_ID].origen[2], 17.0, color, renderer_global);

    color[0] = 0.0; color[1] = 1.0; color[2] = 0.0;
    agregarVector(mi_isocen, 1.0, color, renderer_global);
    agregarEsfera(mi_isocen.origen[0] + mi_isocen.direccion[0], mi_isocen.origen[1] + mi_isocen.direccion[1], mi_isocen.origen[2] + mi_isocen.direccion[2], 20.0, color, renderer_global);
}
#endif




/*  Metodo: agregarVector

    Funcion: Agrega un vector al renderizador definido.
*/
#ifdef BUILD_VTK_VERSION
void RECONS3D::agregarVector(NORCEN org_dir, const double t, double color[3], vtkSmartPointer<vtkRenderer> &mi_renderer){

    double fin[3];
    fin[0] = org_dir.origen[0] + t*org_dir.direccion[0];
    fin[1] = org_dir.origen[1] + t*org_dir.direccion[1];
    fin[2] = org_dir.origen[2] + t*org_dir.direccion[2];


    double normalizedX[3];
    double normalizedY[3];
    double normalizedZ[3];
    vtkMath::Subtract( fin, org_dir.origen, normalizedX);
    double length = vtkMath::Norm(normalizedX);
DEB_MSG("length: " << length);
    vtkMath::Normalize(normalizedX);

    // Z es el producto cruz entre un vector cualquiera y el vector X:
    double arbitrary[3] = {1.0, 1.0, 1.0};
    vtkMath::Cross(normalizedX, arbitrary, normalizedZ);

    // Y es el producto cruz entre X y Z
    vtkMath::Cross(normalizedZ, normalizedX, normalizedY);

    // Generar la matriz base donde se genera el vector:
    vtkSmartPointer<vtkMatrix4x4> matrix_base = vtkSmartPointer<vtkMatrix4x4>::New();
    matrix_base->Identity();
    matrix_base->SetElement(0, 0, normalizedX[0]); matrix_base->SetElement(0, 1, normalizedY[0]); matrix_base->SetElement(0, 2, normalizedZ[0]);
    matrix_base->SetElement(1, 0, normalizedX[1]); matrix_base->SetElement(1, 1, normalizedY[1]); matrix_base->SetElement(1, 2, normalizedZ[1]);
    matrix_base->SetElement(2, 0, normalizedX[2]); matrix_base->SetElement(2, 1, normalizedY[2]); matrix_base->SetElement(2, 2, normalizedZ[2]);

    // Definir la transformacion:
    vtkSmartPointer<vtkTransform> transformacion = vtkSmartPointer<vtkTransform>::New();
    transformacion->Translate( org_dir.origen ); // Mover la flecha al origen
    transformacion->Concatenate( matrix_base ); // rotar la flecha en la direccion dada
    transformacion->Scale(length, length, length); // dimensionar la flecha segun 't'

    // Definir la transformacion del polydata de la flecha:
    vtkSmartPointer<vtkArrowSource> vec_flecha = vtkSmartPointer<vtkArrowSource>::New();
    vec_flecha->SetShaftRadius( 0.005 );
    vec_flecha->SetShaftResolution( 3 );

    vec_flecha->SetTipLength( 0.025 );
    vec_flecha->SetTipRadius( 0.010 );
    vec_flecha->SetTipResolution( 3 );

    vtkSmartPointer<vtkTransformPolyDataFilter> transformar_flecha = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
    transformar_flecha->SetTransform(transformacion);
    transformar_flecha->SetInputConnection(vec_flecha->GetOutputPort());

    // Generar el mapper y el actor:
    vtkSmartPointer<vtkPolyDataMapper> vec_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    vec_mapper->SetInputConnection(transformar_flecha->GetOutputPort());

    vtkSmartPointer<vtkActor> vec_actor = vtkSmartPointer<vtkActor>::New();
    vec_actor->SetMapper(vec_mapper);
    vec_actor->GetProperty()->SetColor( color );

    // Agregar el actor al renderer:
    mi_renderer->AddActor( vec_actor );
}
#endif



/*  Metodo: agregarEjes

    Funcion: Muestra los ejes en que se encuentra la escena global respecto al paciente.
*/
#ifdef BUILD_VTK_VERSION
void RECONS3D::agregarEjes(vtkSmartPointer<vtkRenderer> &mi_renderer){
    double color[3];
    NORCEN eje_dir;
    eje_dir.origen[0] = 0.0;
    eje_dir.origen[1] = 0.0;
    eje_dir.origen[2] = 0.0;

    // Eje X:
    color[0] = 1.0; color[1] = 0.0; color[2] = 0.0;
    eje_dir.direccion[0] = 1.0; eje_dir.direccion[1] = 0.0; eje_dir.direccion[2] = 0.0;
    agregarVector( eje_dir, 500.0, color, mi_renderer);

    // Eje Y:
    color[0] = 0.0; color[1] = 1.0; color[2] = 0.0;
    eje_dir.direccion[0] = 0.0; eje_dir.direccion[1] = 1.0; eje_dir.direccion[2] = 0.0;
    agregarVector( eje_dir, 500.0, color, mi_renderer);

    // Eje Z:
    color[0] = 0.0; color[1] = 0.0; color[2] = 1.0;
    eje_dir.direccion[0] = 0.0; eje_dir.direccion[1] = 0.0; eje_dir.direccion[2] = 1.0;
    agregarVector( eje_dir, 500.0, color, mi_renderer);

}
#endif





/*  Metodo: agregarEsfera

    Funcion: Agrega una esfera al renderizador.
*/
#ifdef BUILD_VTK_VERSION
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
#endif



/*  Metodo: posicionDefecto

    Funcion: Retorna la posicion por defecto del detector y fuente:
*/
#ifdef BUILD_VTK_VERSION
inline RECONS3D::POS RECONS3D::posicionDefecto( const double ancho, const double alto, const double punta){
    POS defecto;
    defecto.puntos[0][0] = -ancho/2; defecto.puntos[0][1] =  alto/2; defecto.puntos[0][2] = 0.0;
    defecto.puntos[1][0] =  ancho/2; defecto.puntos[1][1] =  alto/2; defecto.puntos[1][2] = 0.0;
    defecto.puntos[2][0] =  ancho/2; defecto.puntos[2][1] = -alto/2; defecto.puntos[2][2] = 0.0;
    defecto.puntos[3][0] = -ancho/2; defecto.puntos[3][1] = -alto/2; defecto.puntos[3][2] = 0.0;
    defecto.puntos[4][0] =  0.0; defecto.puntos[4][1] =  0.0; defecto.puntos[4][2] =  punta;
    return defecto;
}
#endif




/*  Metodo: moverPosicion

    Funcion: Mueve los puntos de la imagen base hasta la orientacion definida por los angulos (en el archivo DICOM)
*/
#ifdef BUILD_VTK_VERSION
void RECONS3D::mostrarDetector(const int angios_ID){
    // Rotar los puntos segun LAO/RAO y CAU/CRA:
    const double crl = cos((*(imgs_base + angios_ID))->LAORAO/180.0 * Mi_PI);
    const double srl = sin((*(imgs_base + angios_ID))->LAORAO/180.0 * Mi_PI);
    const double ccc = cos((*(imgs_base + angios_ID))->CRACAU/180.0 * Mi_PI);
    const double scc = sin((*(imgs_base + angios_ID))->CRACAU/180.0 * Mi_PI);

DEB_MSG("Mostrando detector para: " << angios_ID << ", LAORAO: " << (*(imgs_base + angios_ID))->LAORAO << ", CRACAU: " << (*(imgs_base + angios_ID))->CRACAU << ", SID: " << (*(imgs_base + angios_ID))->SID << ", SOD: " << (*(imgs_base + angios_ID))->SOD );

const int mis_rens = (*(imgs_base + angios_ID))->rows;
const int mis_cols = (*(imgs_base + angios_ID))->cols;

    const double mi_pixX = (*(imgs_base + angios_ID))->pixX;
    const double mi_pixY = (*(imgs_base + angios_ID))->pixY;

    //------------------------------------------------------------------------------------------------------------------ FIGURA DEL DETECTOR
    /// Mover el detector y fuente a las posiciones definidas:
    POS det_pos = posicionDefecto( mis_cols*mi_pixX, mis_rens*mi_pixY, mis_rens*mi_pixX/2 );

    // Mover los puntos segun SID y SOD:
    for( int i = 0; i < 5; i++){
        det_pos.puntos[i][2] += (*(imgs_base + angios_ID))->DDP;
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
    det_points->InsertNextPoint( det_pos.puntos[0] );
    det_points->InsertNextPoint( det_pos.puntos[1] );
    det_points->InsertNextPoint( det_pos.puntos[2] );
    det_points->InsertNextPoint( det_pos.puntos[3] );
    det_points->InsertNextPoint( det_pos.puntos[4] );

    //// Definir la forma del detector y fuente:
    vtkSmartPointer<vtkPyramid> det_pyr = vtkSmartPointer<vtkPyramid>::New();
    det_pyr->GetPointIds()->SetId(0, 0);
    det_pyr->GetPointIds()->SetId(1, 1);
    det_pyr->GetPointIds()->SetId(2, 2);
    det_pyr->GetPointIds()->SetId(3, 3);
    det_pyr->GetPointIds()->SetId(4, 4);

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
#endif


//------------------------------------------------------------------------------------------------ PRIVATE ----- ^
//-------------------------------------------------------------------------------------------------- PUBLIC----- v

// M I E M B R O S      P U B L I C O S
void RECONS3D::lengthFilter(IMG_IDX my_img_idx, const int min_length, const int angios_ID)
{
	this->getIMGCONTPointer(my_img_idx, angios_ID)->lengthFilter(min_length);
}





/*  Metodo: medirExactitud

    Funcion: Mide la exactitud entre la imagen segmentada y el ground-truth
*/
double RECONS3D::medirExactitud(const int angios_ID)
{
	/*
    DEB_MSG("midiendo exactitud ...");

    double acc = (*(imgs_base + angios_ID))->medirExactitud();

    char mensaje[512] = "X.XXXXXXXXXXXXXXXX \n";
#if defined(_WIN32) || defined(_WIN64)
	sprintf_s(mensaje, 512 * sizeof(char), "%1.16f\n", acc);
#else
    sprintf(mensaje, "%1.16f\n", acc);
#endif
    escribirLog(mensaje);

    return acc;
	*/
	return 0.0;
}





/*  Metodo: getNangios

    Funcion: Retorna el numero de angiogramas cargados dentro del reconstructor.
*/
int RECONS3D::getNangios(){
    return (n_angios+1);
}



/*  Metodo: agregarInput

    Funcion: Define las rutas de las imagenes que son usadas para reconstruir una arteria coronaria.
*/
void RECONS3D::agregarInput(const char *rutabase_input, const int nivel_l, const int nivel_u, const char *rutaground_input, bool enmascarar) {
	for (int curr_level = nivel_l; curr_level <= nivel_u; curr_level++) {
		n_angios++;

		IMGCONT_FEAT_SET void_img_features_ptr;
		my_img_features_ptr.push_back(void_img_features_ptr);

		my_img_features_ptr[n_angios].my_base_ptr = new IMGCONT();
		my_img_features_ptr[n_angios].my_groundtruth_ptr = new IMGCONT();
		my_img_features_ptr[n_angios].my_response_ptr = new IMGCONT();
		my_img_features_ptr[n_angios].my_response_threshold_ptr = new IMGCONT();
		my_img_features_ptr[n_angios].my_thrs_response_boundaries_ptr = new IMGCONT();
		my_img_features_ptr[n_angios].my_thrs_response_map_dists_ptr = new IMGCONT();
		my_img_features_ptr[n_angios].my_thrs_response_skeleton_ptr = new IMGCONT();


		/* Load the images from their corresponding paths */
		my_img_features_ptr[n_angios].my_base_ptr->Load(rutabase_input, curr_level);
		if (rutaground_input) {
			my_img_features_ptr[n_angios].my_groundtruth_ptr->Load(rutaground_input);
			existe_ground.push_back(true);
		}
		else {
			existe_ground.push_back(false);
		}

		if (enmascarar) {
			my_img_features_ptr[n_angios].my_mask_ptr = new IMGCONT(my_img_features_ptr[n_angios].my_base_ptr->getHeight(),
				my_img_features_ptr[n_angios].my_base_ptr->getWidth(),
				my_img_features_ptr[n_angios].my_base_ptr->getMask());
		}
		else {
			my_img_features_ptr[n_angios].my_mask_ptr = new IMGCONT();
		}
		

		hist.push_back(NULL);
		h_suma.push_back(0.0);
		h_media.push_back(0.0);
		h_desvest.push_back(0.0);
	}

}



/*  Metodo: agregarInput

    Funcion: Define las rutas de las imagenes que son usadas para reconstruir una arteria coronaria.
*/
void RECONS3D::agregarInput( const char *rutabase_input, bool enmascarar ){

    n_angios++;


	my_img_features_ptr[n_angios].my_base_ptr = new IMGCONT();
	my_img_features_ptr[n_angios].my_groundtruth_ptr = new IMGCONT();
	my_img_features_ptr[n_angios].my_response_ptr = new IMGCONT();
	my_img_features_ptr[n_angios].my_response_threshold_ptr = new IMGCONT();
	my_img_features_ptr[n_angios].my_thrs_response_boundaries_ptr = new IMGCONT();
	my_img_features_ptr[n_angios].my_thrs_response_map_dists_ptr = new IMGCONT();
	my_img_features_ptr[n_angios].my_thrs_response_skeleton_ptr = new IMGCONT();


	/* Load the images from their corresponding paths */
	my_img_features_ptr[n_angios].my_base_ptr->Load(rutabase_input);

	existe_ground.push_back(false);

	if (enmascarar) {
		my_img_features_ptr[n_angios].my_mask_ptr = new IMGCONT(my_img_features_ptr[n_angios].my_base_ptr->getHeight(),
			my_img_features_ptr[n_angios].my_base_ptr->getWidth(),
			my_img_features_ptr[n_angios].my_base_ptr->getMask());
	}
	else {
		my_img_features_ptr[n_angios].my_mask_ptr = new IMGCONT(my_img_features_ptr[n_angios].my_base_ptr->getHeight(),
			my_img_features_ptr[n_angios].my_base_ptr->getWidth(), 1.0);
	}

	
    hist.push_back( NULL );
    h_suma.push_back(0.0);
    h_media.push_back(0.0);
    h_desvest.push_back(0.0);

	//mostrarImagen(IMGVTK::BASE, mis_renderers.at(n_angios), n_angios);

}





/*  Metodo: agregarGroundtruth

    Funcion: Define las rutas de las imagenes que son usadas para reconstruir una arteria coronaria.
*/
void RECONS3D::agregarGroundtruth(const char *rutaground_input, const int angios_ID ){
    if( angios_ID > n_angios ){
        char mensaje[512] = "\n<<Error: El angiograma XXX no ha sido agregado al reconstructor>>\n\n";
#if defined(_WIN32) || defined(_WIN64)
		sprintf_s(mensaje, 512 * sizeof(char), "\n<<Error: El angiograma %i no ha sido agregado al reconstructor>>\n\n", angios_ID);
#else
        sprintf(mensaje, "\n<<Error: El angiograma %i no ha sido agregado al reconstructor>>\n\n", angios_ID);
#endif
        escribirLog( mensaje );
    }else{
		my_img_features_ptr[angios_ID].my_groundtruth_ptr->Load(rutaground_input);

        existe_ground[ angios_ID ] = true;
    }
}




/*  Metodo: leerConfiguracion
    Funcion: Carga una configuracion predefinida para el filtro de deteccion
*/
void RECONS3D::leerConfiguracion(const char *ruta_conf){
#if defined(_WIN32) || defined(_WIN64)
	FILE *fp_config;
	fopen_s(&fp_config, ruta_conf, "r");
#else
	FILE *fp_config = fopen(ruta_conf, "r");
#endif

    if( !fp_config ){
        char mensaje_err[512] = "\n<< Error: Cannot open file:\'%s\'\n\n";
#if defined(_WIN32) || defined(_WIN64)
		sprintf_s(mensaje_err, 544 *sizeof(char), "\n<< Error: Cannot open file:\'%s\'\n\n", ruta_conf);
#else
        sprintf(mensaje_err, "\n<< Error: Cannot open file:\'%s\'\n\n", ruta_conf);
#endif
        escribirLog( mensaje_err );
        return;
    }

    char tmp = fgetc( fp_config );
    char tmp_par[128] = "", tmp_val[128] = "";
    char *tmp_str_ptr = tmp_par;

    bool fixed[4] = {false, false, false, false};

    do{
        tmp = fgetc( fp_config );
        if( tmp == '\n' || tmp == EOF ){
            *(tmp_str_ptr) = '\0';

            DEB_MSG( tmp_par << ":" << tmp_val);

            if( strcmp( tmp_par, "GMF" ) == 0 ){
                if( atoi(tmp_val) ){
                    filtro.setFiltro( FILTROS::GMF );
                }
            }else if( strcmp( tmp_par, "SSG") == 0 ){
                if( atoi(tmp_val) ){
                    filtro.setFiltro( FILTROS::SS_GABOR );
                }
            }else if( strcmp( tmp_par, "T") == 0 ){
                if( tmp_val[0] != 'o' ){
                    const double par_t = atof(tmp_val);
                    filtro.setPar( FILTROS::PAR_T, par_t );
                    fixed[0] = true;
                }
            }else if( strcmp( tmp_par, "L") == 0 ){
                if( tmp_val[0] != 'o' ){
                    const double par_l = atof(tmp_val);
                    filtro.setPar( FILTROS::PAR_L, par_l );
                    fixed[1] = true;
                }
            }else if( strcmp( tmp_par, "K") == 0 ){
                if( tmp_val[0] != 'o' ){
                    const double par_k = atof(tmp_val);
                    filtro.setPar( FILTROS::PAR_K, par_k);
                    fixed[2] = true;
                }
            }else if( strcmp( tmp_par, "sig") == 0 ){
                if( tmp_val[0] != 'o' ){
                    const double par_s = atof(tmp_val);
                    filtro.setPar( FILTROS::PAR_SIGMA, par_s);
                    fixed[3] = true;
                }
            }else if( strcmp( tmp_par, "Tinf") == 0 && !fixed[0] ){
                const double inf_t = atof(tmp_val);
                filtro.setLim( FILTROS::PAR_T, FILTROS::INFERIOR, inf_t );
            }else if( strcmp( tmp_par, "Tsup") == 0 && !fixed[0] ){
                const double sup_t = atof(tmp_val);
                filtro.setLim( FILTROS::PAR_T, FILTROS::SUPERIOR, sup_t );
            }else if( strcmp( tmp_par, "Tdel") == 0 && !fixed[0] ){
                const double del_t = atof(tmp_val);
                filtro.setLim( FILTROS::PAR_T, FILTROS::DELTA, del_t );
            }else if( strcmp( tmp_par, "Linf") == 0 && !fixed[1] ){
                const double inf_l = atof(tmp_val);
                filtro.setLim( FILTROS::PAR_L, FILTROS::INFERIOR, inf_l );
            }else if( strcmp( tmp_par, "Lsup") == 0 && !fixed[1] ){
                const double sup_l = atof(tmp_val);
                filtro.setLim( FILTROS::PAR_L, FILTROS::SUPERIOR, sup_l );
            }else if( strcmp( tmp_par, "Ldel") == 0 && !fixed[1] ){
                const double del_l = atof(tmp_val);
                filtro.setLim( FILTROS::PAR_L, FILTROS::DELTA, del_l );
            }else if( strcmp( tmp_par, "Kinf") == 0 && !fixed[2] ){
                const double inf_k = atof(tmp_val);
                filtro.setLim( FILTROS::PAR_K, FILTROS::INFERIOR, inf_k );
            }else if( strcmp( tmp_par, "Ksup") == 0 && !fixed[2] ){
                const double sup_k = atof(tmp_val);
                filtro.setLim( FILTROS::PAR_K, FILTROS::SUPERIOR, sup_k );
            }else if( strcmp( tmp_par, "Kdel") == 0 && !fixed[2] ){
                const double del_k = atof(tmp_val);
                filtro.setLim( FILTROS::PAR_K, FILTROS::DELTA, del_k );
            }else if( strcmp( tmp_par, "Sinf") == 0 && !fixed[3] ){
                const double inf_s = atof(tmp_val);
                filtro.setLim( FILTROS::PAR_SIGMA, FILTROS::INFERIOR, inf_s );
            }else if( strcmp( tmp_par, "Ssup") == 0 && !fixed[3] ){
                const double sup_s = atof(tmp_val);
                filtro.setLim( FILTROS::PAR_SIGMA, FILTROS::SUPERIOR, sup_s );
            }else if( strcmp( tmp_par, "Sdel") == 0 && !fixed[3] ){
                const double del_s = atof(tmp_val);
                filtro.setLim( FILTROS::PAR_SIGMA, FILTROS::DELTA, del_s );
            }else if( strcmp( tmp_par, "ROC") == 0 ){
                if( atoi(tmp_val) ){
                    filtro.setFitness( FILTROS::ROC );
                }
            }else if( strcmp( tmp_par, "CC") == 0 ){
                if( atoi(tmp_val) ){
                    filtro.setFitness( FILTROS::CORCON );
                }
            }else if( strcmp( tmp_par, "Unset") == 0 ){
                if( atoi(tmp_val) ){
                    filtro.setEvoMet( FILTROS::EVO_UNSET );
                }
            }else if( strcmp( tmp_par, "Exh") == 0 ){
                if( atoi(tmp_val) ){
                    filtro.setEvoMet( FILTROS::EXHAUSTIVA );
                }
            }else if( strcmp( tmp_par, "GA") == 0 ){
                if( atoi(tmp_val) ){
                    filtro.setEvoMet( FILTROS::EA_GA );
                }
            }else if( strcmp( tmp_par, "UMDA") == 0 ){
                if( atoi(tmp_val) ){
                    filtro.setEvoMet( FILTROS::EDA_UMDA );
                }
            }else if( strcmp( tmp_par, "BUMDA") == 0 ){
                if( atoi(tmp_val) ){
                    filtro.setEvoMet( FILTROS::EDA_BUMDA );
                }
            }else if( strcmp( tmp_par, "DE") == 0 ){
                if( atoi(tmp_val) ){
                    filtro.setEvoMet( FILTROS::EA_DE );
                }
            }else if( strcmp( tmp_par, "pop_size") == 0 ){
                const int pop_size = atoi(tmp_val);
                filtro.setEvoMetPar( FILTROS::POPSIZE, (double)pop_size);
            }else if( strcmp( tmp_par, "max_gen") == 0 ){
                const int max_gen = atoi(tmp_val);
                filtro.setEvoMetPar( FILTROS::MAXGEN, (double)max_gen);
            }else if( strcmp( tmp_par, "CR") == 0 ){
                const double cr = atof(tmp_val);
                DEB_MSG("CR: " << cr << " / " << tmp_val);
                filtro.setEvoMetPar( FILTROS::CR, cr);
            }else if( strcmp( tmp_par, "MR") == 0 ){
                const double mr = atof(tmp_val);
                filtro.setEvoMetPar( FILTROS::MR, mr);
            }

            tmp_str_ptr = tmp_par;

        }else if( tmp == '\t'){
            *(tmp_str_ptr) = '\0';
            tmp_str_ptr = tmp_val;
        }else{
            *(tmp_str_ptr) = tmp;
            tmp_str_ptr++;
        }
    }while ( tmp != EOF);

    fclose( fp_config );
}




/*  Metodo: setFiltroEntrenamiento

    Funcion:
*/
void RECONS3D::setFiltroEntrenamiento(const FILTROS::EVO_MET evo_met){
    filtro.setEvoMet(evo_met);
}



/*  Metodo: setFiltroEntrenamientoPars

    Funcion:
*/
void RECONS3D::setFiltroEntrenamientoPars(const FILTROS::EVO_MET_PAR evo_par, const double val){
    filtro.setEvoMetPar(evo_par, val);
}


/*  Metodo: setFiltroEval

    Funcion:
*/
void RECONS3D::setFiltroEval(const FILTROS::FITNESS fit_fun){
    filtro.setFitness( fit_fun );
}



/*  Metodo: setFiltroMetodo

    Funcion:
*/
void RECONS3D::setFiltroMetodo(const FILTROS::SEG_FILTRO metodo_filtrado){
    filtro.setFiltro( metodo_filtrado );
}




/*  Metodo: setFiltroParametros

    Funcion:
*/
void RECONS3D::setFiltroParametros(const FILTROS::PARAMETRO par, const FILTROS::LIMITES lim, const double val){
    filtro.setLim( par, lim, val );
}



/*  Metodo: setFiltroParametros

    Funcion:
*/
void RECONS3D::setFiltroParametros(const FILTROS::PARAMETRO par, const double val){
    filtro.setPar(par, val);
}



/*  Metodo: setFiltroParametros

    Funcion:
*/
void RECONS3D::setFiltroParametros(){
    filtro.setPar();
}




/*  Metodo: Guardar

    Funcion: Almacena una imagen del reconstructor a la ruta especificada, en formato PNG o PGM.
*/
void RECONS3D::Guardar(const char *ruta, IMG_IDX img_idx, IMGCONT::IMG_TYPE tipo_img, const int angios_ID){
	this->getIMGCONTPointer(img_idx, angios_ID)->Save(ruta, tipo_img);
}


/*  Metodo: getRows

    Funcion:
*/
int RECONS3D::getRows( const int angios_ID ){
	return this->getIMGCONTPointer(IMG_BASE, angios_ID)->getHeight();
}



/*  Metodo: getCols

    Funcion:
*/
int RECONS3D::getCols( const int angios_ID ){
	return this->getIMGCONTPointer(IMG_BASE, angios_ID)->getWidth();
}




/*  Metodo: get_pixelData

    Funcion: Retorna el apuntador a la informacion de la imagen seleccionada con img_idx
*/
double* RECONS3D::get_pixelData(const int angios_ID, IMG_IDX img_idx) {
	/*
	double *img_ptr = NULL;

	switch(img_idx){
	case IMGVTK::BASE:
		img_ptr = (*(imgs_base + angios_ID))->base_ptr;
		break;
	case IMGVTK::GROUNDTRUTH:
		img_ptr = (*(imgs_base + angios_ID))->gt_ptr;
		break;
	case IMGVTK::MASK:
		img_ptr = (*(imgs_base + angios_ID))->mask_ptr;
		break;
	case IMGVTK::SKELETON:
		img_ptr = (*(imgs_base + angios_ID))->skl_ptr;
		break;
	case IMGVTK::SEGMENT:
		img_ptr = (*(imgs_base + angios_ID))->segment_ptr;
		break;
	case IMGVTK::THRESHOLD:
		img_ptr = (*(imgs_base + angios_ID))->threshold_ptr;
		break;
	case IMGVTK::MAPDIST:
		img_ptr = (*(imgs_base + angios_ID))->map_ptr;
		break;
	case IMGVTK::BORDERS:
		img_ptr = (*(imgs_base + angios_ID))->borders_ptr;
		break;
	}

	return img_ptr;
	*/
	return NULL;
}



/*  Metodo: setLog

    Funcion: Define la ruta donde se guarda el log del proceso de reconstruccion.
*/
void RECONS3D::setLog(const char *ruta_log){
    if(fp_log){
        fclose(fp_log);
        fp_log = NULL;
    }
#if defined(_WIN32) || defined(_WIN64)
	fopen_s(&fp_log, ruta_log, "w");
#else
    fp_log = fopen(ruta_log, "w");
#endif
}


/*  Metodo: setLog
 * 
 *    Funcion: Define el objeto tipo QPlainTextEdit donde se guarda el log del proceso de reconstruccion.
 */
#ifdef BUILD_GUI_VERSION
void RECONS3D::setLog(QTextEdit *txtLog)
{
    mi_txtLog =  txtLog;
}
#endif



/*  Metodo: setFiltroLog

    Funcion: Define el objeto tipo QPlainTextEdit donde se guarda el log del proceso de reconstruccion.
*/
#ifdef BUILD_GUI_VERSION
void RECONS3D::setFiltroLog(QTextEdit *txtLog)
{
    filtro.setLog( txtLog );
}
#endif



/*  Metodo: segmentarImagenBase

    Funcion: Aplica el filtro con los parametros definidos
*/
void RECONS3D::segmentarImagenBase( const int angios_ID ){

    filtro.setInput(*(imgs_base + angios_ID), angios_ID, angios_ID);

    if( filtro.getParametrosOptimizar() > 0 ){

        DEB_MSG("Comenzando entrenamiento de parametros ...");

        filtro.setPar();
    }

    TIMERS;
    GETTIME_INI;
    filtro.filtrar();
    GETTIME_FIN;
    char mensaje[512] = "XXX.XXXXXXXX\n";
#if defined(_WIN32) || defined(_WIN64)
	sprintf_s(mensaje, 512 *sizeof(char), "%3.8f\n", DIFTIME);
#else
    sprintf(mensaje, "%3.8f\n", DIFTIME);
#endif
    escribirLog( mensaje );
}



/*  Metodo: mostrarRadios

    Funcion: Muestra el radio de cada pixel del esqueleto.
*/
#ifdef BUILD_VTK_VERSION
void RECONS3D::mostrarRadios(vtkSmartPointer<vtkPoints> puntos, vtkSmartPointer<vtkCellArray> cilindros, int *n_pix, IMGVTK::PIX_PAR *grafo, const double DDP, const double crl, const double srl, const double ccc, const double scc, const int nivel_detalle, FILE *fp_cilindros){

    const double theta_inc = 2 * Mi_PI / (double)detalle;

    const double xx = grafo->x;
    const double yy = grafo->y;

    for( int h = 0; h < grafo->n_hijos; h++){

        double radio = grafo->radio;

        IMGVTK::PIX_PAR *sig_hijo = grafo->hijos[h];
        int n_prof = 0;
        while( 1 ){
            radio += sig_hijo->radio;
            n_prof++;

            if( sig_hijo->pix_tipo != IMGVTK::PIX_SKL || n_prof >= nivel_detalle){
                break;
            }

            sig_hijo = sig_hijo->hijos[0];
        }

        /// Promedio de los radios de esta seccion:
        radio /= (double)n_prof;

        double x_fin = sig_hijo->x;
        double y_fin = sig_hijo->y;

        const double alpha = atan2(y_fin - yy, x_fin - xx) + Mi_PI / 2.0;

        const double cal = cos(alpha);
        const double sal = sin(alpha);

        double r_temp;

        double theta = 0.0;


        /// Primera iteracion:
        {
            const double cth = cos(theta);
            const double sth = sin(theta);

            /// Puntos inicio:
            double xx_3D_ini = xx + cal*cth*radio;
            double yy_3D_ini = yy + sal*cth*radio;
            double zz_3D_ini = DDP - sth*radio;

            // Mover los puntos segun indica el LAO/RAO y CRA/CAU:
            //// Rotacion usando el eje 'x' como base:
            r_temp = crl*yy_3D_ini - srl*zz_3D_ini;
            zz_3D_ini = srl*yy_3D_ini + crl*zz_3D_ini;
            yy_3D_ini = r_temp;


            //// Rotacion usando el eje 'y' como base:
            r_temp = ccc*xx_3D_ini - scc*zz_3D_ini;
            zz_3D_ini = scc*xx_3D_ini + ccc*zz_3D_ini;
            xx_3D_ini = r_temp;

            puntos->InsertNextPoint(xx_3D_ini, yy_3D_ini, zz_3D_ini);

            /// Puntos final:
            double xx_3D_fin = x_fin + cal*cth*radio;
            double yy_3D_fin = y_fin + sal*cth*radio;
            double zz_3D_fin = DDP - sth*radio;

            // rotar el punto para que se dirija hacia el punto con el que se midio el radio.

            // Mover los puntos segun indica el SID y SOD:
            //// Rotacion usando el eje 'x' como base:
            r_temp = crl*yy_3D_fin - srl*zz_3D_fin;
            zz_3D_fin = srl*yy_3D_fin + crl*zz_3D_fin;
            yy_3D_fin = r_temp;


            //// Rotacion usando el eje 'y' como base:
            r_temp = ccc*xx_3D_fin - scc*zz_3D_fin;
            zz_3D_fin = scc*xx_3D_fin + ccc*zz_3D_fin;
            xx_3D_fin = r_temp;

            puntos->InsertNextPoint(xx_3D_fin, yy_3D_fin, zz_3D_fin);

            theta += theta_inc;
        }

        for( int i = 1; i <= detalle; i++){
            const double cth = cos(theta);
            const double sth = sin(theta);

            /// Puntos inicio:
            double xx_3D_ini = xx + cal*cth*radio;
            double yy_3D_ini = yy + sal*cth*radio;
            double zz_3D_ini = DDP - sth*radio;

            // Mover los puntos segun indica el LAO/RAO y CRA/CAU:
            //// Rotacion usando el eje 'x' como base:
            r_temp = crl*yy_3D_ini - srl*zz_3D_ini;
            zz_3D_ini = srl*yy_3D_ini + crl*zz_3D_ini;
            yy_3D_ini = r_temp;


            //// Rotacion usando el eje 'y' como base:
            r_temp = ccc*xx_3D_ini - scc*zz_3D_ini;
            zz_3D_ini = scc*xx_3D_ini + ccc*zz_3D_ini;
            xx_3D_ini = r_temp;

            puntos->InsertNextPoint(xx_3D_ini, yy_3D_ini, zz_3D_ini);

            /// Puntos final:
            double xx_3D_fin = x_fin + cal*cth*radio;
            double yy_3D_fin = y_fin + sal*cth*radio;
            double zz_3D_fin = DDP - sth*radio;

            // rotar el punto para que se dirija hacia el punto con el que se midio el radio.

            // Mover los puntos segun indica el SID y SOD:
            //// Rotacion usando el eje 'x' como base:
            r_temp = crl*yy_3D_fin - srl*zz_3D_fin;
            zz_3D_fin = srl*yy_3D_fin + crl*zz_3D_fin;
            yy_3D_fin = r_temp;


            //// Rotacion usando el eje 'y' como base:
            r_temp = ccc*xx_3D_fin - scc*zz_3D_fin;
            zz_3D_fin = scc*xx_3D_fin + ccc*zz_3D_fin;
            xx_3D_fin = r_temp;

            puntos->InsertNextPoint(xx_3D_fin, yy_3D_fin, zz_3D_fin);

            /// Generar un nuevoQuad con estos dos puntos y los anteriores:
            vtkSmartPointer<vtkQuad> nuevo_quad = vtkSmartPointer<vtkQuad>::New();

            nuevo_quad->GetPointIds()->SetId(0, *(n_pix) + 2*(i-1));
            nuevo_quad->GetPointIds()->SetId(1, *(n_pix) + 2*(i-1)+1);
            nuevo_quad->GetPointIds()->SetId(2, *(n_pix) + 2*i+1);
            nuevo_quad->GetPointIds()->SetId(3, *(n_pix) + 2*i);

            cilindros->InsertNextCell(nuevo_quad);

            theta += theta_inc;
        }

        if( fp_cilindros ){
            double xx_3D_ini = xx;
            double yy_3D_ini = yy;
            double zz_3D_ini = DDP;

            //// Rotacion usando el eje 'x' como base:
            r_temp = crl*yy_3D_ini - srl*zz_3D_ini;
            zz_3D_ini = srl*yy_3D_ini + crl*zz_3D_ini;
            yy_3D_ini = r_temp;


            //// Rotacion usando el eje 'y' como base:
            r_temp = ccc*xx_3D_ini - scc*zz_3D_ini;
            zz_3D_ini = scc*xx_3D_ini + ccc*zz_3D_ini;
            xx_3D_ini = r_temp;

            double xx_3D_fin = x_fin;
            double yy_3D_fin = y_fin;
            double zz_3D_fin = DDP;

            //// Rotacion usando el eje 'x' como base:
            r_temp = crl*yy_3D_fin - srl*zz_3D_fin;
            zz_3D_fin = srl*yy_3D_fin + crl*zz_3D_fin;
            yy_3D_fin = r_temp;


            //// Rotacion usando el eje 'y' como base:
            r_temp = ccc*xx_3D_fin - scc*zz_3D_fin;
            zz_3D_fin = scc*xx_3D_fin + ccc*zz_3D_fin;
            xx_3D_fin = r_temp;

            fprintf(fp_cilindros, "%f %f %f %f %f %f %f\n",xx_3D_ini, yy_3D_ini, zz_3D_ini, radio, xx_3D_fin, yy_3D_fin, zz_3D_fin);
        }

        *(n_pix) = *(n_pix) + (detalle+1)*2;
        mostrarRadios(puntos, cilindros, n_pix, sig_hijo, DDP, crl, srl, ccc, scc, nivel_detalle, fp_cilindros);
    }
}
#endif







/*  Metodo: mostrarRadios

    Funcion: Muestra el radio de cada pixel del esqueleto.
*/
void RECONS3D::mostrarRadios(int *n_pix, IMGCONT::PIX_PAIR *grafo, const double DDP, const double crl, const double srl, const double ccc, const double scc, const int nivel_detalle, FILE *fp_cilindros){

    const double xx = grafo->my_pos_x;
    const double yy = grafo->my_pos_y;

    for( int h = 0; h < grafo->my_n_children; h++){

        double radio = grafo->my_radious;

        IMGCONT::PIX_PAIR *sig_hijo = grafo->my_children[h];
        int n_prof = 0;
        while( 1 ){
            radio += sig_hijo->my_radious;
            n_prof++;

            if( sig_hijo->my_pix_type != IMGCONT::PIX_SKL || n_prof >= nivel_detalle){
                break;
            }

            sig_hijo = sig_hijo->my_children[0];
        }

        /// Promedio de los radios de esta seccion:
        radio /= (double)n_prof;

        double x_fin = sig_hijo->my_pos_x;
        double y_fin = sig_hijo->my_pos_y;
        double r_temp;

        if( fp_cilindros ){
            double xx_3D_ini = xx;
            double yy_3D_ini = yy;
            double zz_3D_ini = DDP;

            //// Rotacion usando el eje 'x' como base:
            r_temp = crl*yy_3D_ini - srl*zz_3D_ini;
            zz_3D_ini = srl*yy_3D_ini + crl*zz_3D_ini;
            yy_3D_ini = r_temp;


            //// Rotacion usando el eje 'y' como base:
            r_temp = ccc*xx_3D_ini - scc*zz_3D_ini;
            zz_3D_ini = scc*xx_3D_ini + ccc*zz_3D_ini;
            xx_3D_ini = r_temp;

            double xx_3D_fin = x_fin;
            double yy_3D_fin = y_fin;
            double zz_3D_fin = DDP;

            //// Rotacion usando el eje 'x' como base:
            r_temp = crl*yy_3D_fin - srl*zz_3D_fin;
            zz_3D_fin = srl*yy_3D_fin + crl*zz_3D_fin;
            yy_3D_fin = r_temp;


            //// Rotacion usando el eje 'y' como base:
            r_temp = ccc*xx_3D_fin - scc*zz_3D_fin;
            zz_3D_fin = scc*xx_3D_fin + ccc*zz_3D_fin;
            xx_3D_fin = r_temp;

            fprintf(fp_cilindros, "%f %f %f %f %f %f %f\n",xx_3D_ini, yy_3D_ini, zz_3D_ini, radio, xx_3D_fin, yy_3D_fin, zz_3D_fin);
        }

        *(n_pix) = *(n_pix) + (detalle+1)*2;
        mostrarRadios(n_pix, sig_hijo, DDP, crl, srl, ccc, scc, nivel_detalle, fp_cilindros);
    }
}





/*  Metodo: mostrarRadios

    Funcion: Muestra el radio de cada pixel del esqueleto.
*/
#ifdef BUILD_VTK_VERSION
void RECONS3D::mostrarRadios(vtkSmartPointer<vtkPoints> puntos, vtkSmartPointer<vtkCellArray> vert_skl, vtkSmartPointer<vtkUnsignedCharArray> grafo_nivel, int *n_pix, IMGVTK::PIX_PAR *grafo, const double DDP, const double crl, const double srl, const double ccc, const double scc, const int n_niveles){

    unsigned char nivel[3];

    nivel[0] = (unsigned char)(255.0 * (double)(grafo->nivel + 1) / (double)n_niveles);
    nivel[1] = (unsigned char)(255.0 * (double)(grafo->nivel + 1) / (double)n_niveles);
    nivel[2] = (unsigned char)(255.0 * (double)(grafo->nivel + 1) / (double)n_niveles);

    double theta = 0.0;
    const double theta_inc = 2 * Mi_PI / (double)detalle;
    const double xx = grafo->x;
    const double yy = grafo->y;
    const double radio = grafo->radio;
    const double alpha = grafo->alpha;

    const double cal = cos(alpha);
    const double sal = sin(alpha);
    double r_temp;

    for( int i = 0; i < detalle; i++){
        const double cth = cos(theta);
        const double sth = sin(theta);
        double xx_3D = xx + cal*cth*radio;
        double yy_3D = yy + sal*cth*radio;
        double zz_3D = DDP - sth*radio;

        // rotar el punto para que se dirija hacia el punto con el que se midio el radio.

        // Mover los puntos segun indica el SID y SOD:
        //// Rotacion usando el eje 'x' como base:
        r_temp = crl*yy_3D - srl*zz_3D;
        zz_3D = srl*yy_3D + crl*zz_3D;
        yy_3D = r_temp;


        //// Rotacion usando el eje 'y' como base:
        r_temp = ccc*xx_3D - scc*zz_3D;
        zz_3D = scc*xx_3D + ccc*zz_3D;
        xx_3D = r_temp;

        puntos->InsertNextPoint(xx_3D, yy_3D, zz_3D);
        vtkSmartPointer< vtkVertex > pix = vtkSmartPointer< vtkVertex >::New();
        pix->GetPointIds()->SetId(0, *n_pix);

        *n_pix = *n_pix + 1;
        vert_skl->InsertNextCell(pix);

        grafo_nivel->InsertNextTupleValue( nivel );
        theta += theta_inc;
    }

    double xx_3D = xx;
    double yy_3D = yy;
    double zz_3D = DDP;

    // Mover los puntos segun indica el SID y SOD:
    //// Rotacion usando el eje 'x' como base:
    r_temp = crl*yy_3D - srl*zz_3D;
    zz_3D = srl*yy_3D + crl*zz_3D;
    yy_3D = r_temp;

    //// Rotacion usando el eje 'y' como base:
    r_temp = ccc*xx_3D - scc*zz_3D;
    zz_3D = scc*xx_3D + ccc*zz_3D;
    xx_3D = r_temp;

    puntos->InsertNextPoint(xx_3D, yy_3D, zz_3D);
    vtkSmartPointer< vtkVertex > pix = vtkSmartPointer< vtkVertex >::New();
    pix->GetPointIds()->SetId(0, *n_pix);

    *n_pix = *n_pix + 1;
    vert_skl->InsertNextCell(pix);

    for( int i = 0; i < grafo->n_hijos; i++){
        mostrarRadios(puntos, vert_skl, grafo_nivel, n_pix, grafo->hijos[i], DDP, crl, srl, ccc, scc, n_niveles);
    }
}
#endif




/*  Metodo: maxminAnchos
    Funcion: Obtiene el maximo y minimo de los anchos estimados.
*/
void maxminAnchos( IMGCONT::PIX_PAIR *grafo, double *max, double *min){
    if( (2 * grafo->my_radious) > *max ){
        *max = 2 * grafo->my_radious;
    }


    if( (2 * grafo->my_radious) < *min ){
        *min = 2 * grafo->my_radious;
    }


    for( int i = 0; i < grafo->my_n_children; i++){
        maxminAnchos( grafo->my_children[i], max, min );
    }
}




/*  Metodo: clasAnchos
    Funcion: Clasifica los anchos en las clases del histograma ( las clases tienen un rango interclase de 0.5's)
*/
void RECONS3D::clasAnchos(IMGCONT::PIX_PAIR *my_graph_root, const int angios_ID){
    const double ancho = my_graph_root->my_radious * 2.0;

    const int clase = (int)floor(ancho);
    const double residuo = ancho - floor(ancho);
    const int clase_med = (residuo >= 0.75 ) ? 2 : ((residuo >= 0.25) ? 1 : 0);

    hist[angios_ID][2*clase + clase_med]++;

    for( int i = 0; i < my_graph_root->my_n_children; i++){
        clasAnchos(my_graph_root->my_children[i], angios_ID );
    }

}


/*  Metodo: skeletonize

    Funcion: Obtiene el esqueleto de la imagen y genera el histograma de los anchos
*/
void RECONS3D::skeletonize(const int angios_ID){
	/*
	DEB_MSG("Extrayendo esquelto a " << angios_ID);
    if( !(*(imgs_base + angios_ID))->skl_ptr ){

		(*(imgs_base + angios_ID))->skeletonization(IMGVTK::THRESHOLD);

        DEB_MSG("Caracteristicas puntero: " << (*(imgs_base + angios_ID))->pix_caract);

        // Generar el histograma:
        //// Recorrer los anchos y detectar el minimo y el maximo
        if ( !(*(imgs_base + angios_ID))->pix_caract ){
            DEB_MSG("No existe grafo alguno...");
            return;
        }else{
            DEB_MSG("Generando Histograma...");
        }

        double max = -MY_INF;
        double min =  MY_INF;

        maxminAnchos((*(imgs_base + angios_ID))->pix_caract, &max, &min );

        DEB_MSG("Max ancho: " << max << ", min ancho: " << min);

        // Armar el histograma:

        const double residuo = max - floor(max);
        double inc = 0.0;
        if( residuo > 0.5 ){
            inc = 1.0;
        }else if( residuo > 0.0 ){
            inc = 0.5;
        }

        max = floor(max) + inc;
        const int n_clases = (int)(max * 2.0) + 1;

        DEB_MSG("n clases: " << n_clases);
        hist.at( angios_ID ) = new int [n_clases];
        memset(hist[ angios_ID ], 0, n_clases*sizeof(int) );
        clasAnchos((*(imgs_base + angios_ID))->pix_caract, angios_ID);

        for(int i = 0; i < n_clases; i++){
            h_suma[ angios_ID ] = h_suma[ angios_ID ] + (double)i * 0.5 * hist[angios_ID][ i ];
        }
        h_media[angios_ID] = h_suma[angios_ID] / (double)n_clases;

        for(int i = 0; i < n_clases; i++){
            double temp = (double)i * 0.5 * hist[angios_ID][ i ] - h_media[angios_ID];
            h_desvest[angios_ID] = h_desvest[angios_ID] + temp * temp;
        }

        h_desvest[angios_ID] = sqrt( h_desvest[angios_ID] / (double)n_clases );

#ifdef BUILD_VTK_VERSION
        // Inicializar los objetos para la visualizacion con VTK
        view[angios_ID]->GetRenderer()->SetBackground(1.0, 1.0, 1.0);
        view[angios_ID]->GetRenderWindow()->SetSize(400, 300);

        vtkSmartPointer< vtkChartXY > hist_chart = vtkSmartPointer< vtkChartXY >::New();
        view[angios_ID]->GetScene()->AddItem(hist_chart);

        vtkSmartPointer< vtkTable > table =  vtkSmartPointer< vtkTable >::New();

        vtkSmartPointer<vtkDoubleArray> arrAnchos = vtkSmartPointer<vtkDoubleArray>::New();
        arrAnchos->SetName("Vessel width estimation");
        table->AddColumn(arrAnchos);

        vtkSmartPointer<vtkIntArray> arrFrecs = vtkSmartPointer<vtkIntArray>::New();
        arrFrecs->SetName("Frecuency");
        table->AddColumn(arrFrecs);
        table->SetNumberOfRows(n_clases);

        for(int i = 0; i < n_clases; i++){
            table->SetValue(i,0,i*0.5);
            table->SetValue(i,1,hist[angios_ID][i]);
        }

        vtkPlot *line = NULL;
        line = hist_chart->AddPlot(vtkChart::BAR);
        line->SetColor(0.0, 0.0, 1.0);

        hist_chart->GetAxis(vtkAxis::LEFT)->SetTitle("Frecuency");
        hist_chart->GetAxis(vtkAxis::BOTTOM)->SetTitle("Vessel width estimation");

        line->SetInputData(table, 0, 1);
#endif
    }
	*/
}





/*  Metodo: skeletonize

    Funcion: Obtiene el esqueleto de la imagen y muestra los puntos de interes.
*/
void RECONS3D::skeletonize(const int angios_ID, const int nivel_detalle){
	/*
    skeletonize(angios_ID);

    if( !(*(imgs_base + angios_ID))->pix_caract ){
        DEB_MSG("No existe grafo alguno...");
        return;
    }

    // Rotar los puntos segun LAO/RAO y CAU/CRA:
    const double crl = cos((*(imgs_base + angios_ID))->LAORAO/180.0 * Mi_PI);
    const double srl = sin((*(imgs_base + angios_ID))->LAORAO/180.0 * Mi_PI);
    const double ccc = cos((*(imgs_base + angios_ID))->CRACAU/180.0 * Mi_PI);
    const double scc = sin((*(imgs_base + angios_ID))->CRACAU/180.0 * Mi_PI);

    const double DDP = (*(imgs_base + angios_ID))->DDP;
    const int n_niveles = (*(imgs_base + angios_ID))->n_niveles;

DEB_MSG("numero de niveles: " << n_niveles);

    // Recorrer el grafo y generar en 3D una burda reconstruccion.
DEB_MSG("Mostrando los radios de cada " << nivel_detalle << " pixeles del esqueleto...");

#ifdef BUILD_VTK_VERSION
    vtkSmartPointer< vtkPoints > puntos = vtkSmartPointer< vtkPoints >::New();
    vtkSmartPointer< vtkCellArray > cilindros = vtkSmartPointer< vtkCellArray >::New();
#endif

    int n_pix = 0;

    char nom_cilindros[] = "cilindros_XXX_YYY.dat";
#if defined(_WIN32) || defined(_WIN64)
	sprintf_s(nom_cilindros, "cilindros_%i_%i.dat", angios_ID, nivel_detalle);
	FILE *fp_cilindros;
	fopen_s(&fp_cilindros, nom_cilindros, "w");
#else
    sprintf( nom_cilindros, "cilindros_%i_%i.dat", angios_ID, nivel_detalle);
	FILE *fp_cilindros = fopen(nom_cilindros, "w");
#endif

    fprintf( fp_cilindros, "X1 Y1 Z1 RADIO X2 Y2 Z2\n");

#ifdef BUILD_VTK_VERSION
    mostrarRadios(puntos, cilindros, &n_pix, (*(imgs_base + angios_ID))->pix_caract, DDP, crl, srl, ccc, scc, nivel_detalle, fp_cilindros);
#else
    mostrarRadios(&n_pix, (*(imgs_base + angios_ID))->pix_caract, DDP, crl, srl, ccc, scc, nivel_detalle, fp_cilindros);
#endif

    fclose( fp_cilindros );

#ifdef BUILD_VTK_VERSION
    /// Generar los cilindros:
    vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
    polydata->SetPoints(puntos);
    polydata->SetPolys(cilindros);

    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputData(polydata);
    DEB_MSG("Generando Mapper");

    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    DEB_MSG("Generando Actor");

    double color[] = {0.0, 1.0, 0.0};
    actor->GetProperty()->SetColor(color);

    DEB_MSG("Agregando al renderizador");
    renderer_global->AddActor( actor );

    renderizar( renderer_global );
#endif
*/
}



/*  Metodo: getHist

    Funcion: Retorna la ventana del visualizador donde se muestra el histograma.
*/
#ifdef BUILD_VTK_VERSION
vtkRenderWindow* RECONS3D::getHist( const int angios_ID ){
    return view[angios_ID]->GetRenderWindow();
}
#endif



double RECONS3D::getHist_desvest( const int angios_ID )
{
    return h_desvest[angios_ID];
}

double RECONS3D::getHist_media( const int angios_ID )
{
    return h_media[angios_ID];
}

double RECONS3D::getHist_suma( const int angios_ID )
{
    return h_suma[angios_ID];
}



/*  Metodo: setIteratorHist

    Funcion: Define el interactor para visualizar el histograma.
*/
#ifdef BUILD_VTK_VERSION
void RECONS3D::setIteratorHist( vtkRenderWindowInteractor *interactor, const int angios_ID ){
    view[angios_ID]->SetInteractor( interactor );
}
#endif



/*  Metodo: mostrarBase

    Funcion: Muestra la imagen base en el renderizador del input 'angios_ID'
*/
#ifdef BUILD_VTK_VERSION
void RECONS3D::mostrarBase( const int angios_ID ){
    if( angios_ID > n_angios ){
        char mensaje_err[] = "\n<<Error: El angiograma XXX no ha sido agregado al reconstructor>>\n\n";
        sprintf(mensaje_err, "\n<<Error: El angiograma %i no ha sido agregado al reconstructor>>\n\n", angios_ID);
        DEB_MSG(mensaje_err);
        escribirLog( mensaje_err );
    }else{
        mostrarImagen( IMGVTK::BASE, mis_renderers[ angios_ID ], angios_ID);
    }
}
#endif



/*  Metodo: mostrarGroundtruth

    Funcion: Muestra la imagen ground-truth en el renderizador del input 'angios_ID'
*/
#ifdef BUILD_VTK_VERSION
void RECONS3D::mostrarGroundtruth( const int angios_ID ){
    if( (angios_ID > n_angios) || (!existe_ground.at(angios_ID)) ){
        char mensaje[512] = "\n<<Error: El ground-truth para el angiograma XXX no ha sido agregado al reconstructor>>\n\n";
#if defined(_WIN32) || defined(_WIN64)
		sprintf_s(mensaje, 512 * sizeof(char), "\n<<Error: El ground-truth para el angiograma %i no ha sido agregado al reconstructor>>\n\n", angios_ID);
#else
        sprintf(mensaje, "\n<<Error: El ground-truth para el angiograma %i no ha sido agregado al reconstructor>>\n\n", angios_ID);
#endif
        DEB_MSG(mensaje);
        escribirLog( mensaje );
    }else{
        mostrarImagen( IMGVTK::GROUNDTRUTH, mis_renderers[ angios_ID ], angios_ID);
    }
}
#endif


/*  Metodo: getRenderer

    Funcion: Retorna el renderizador global.
*/
#ifdef BUILD_VTK_VERSION
vtkSmartPointer< vtkRenderer > RECONS3D::getRenderer(){
    return renderer_global;
}
#endif


/*  Metodo: getRenderer

    Funcion: Retorna el renderizador correspondiente al input 'angios_ID'.
*/
#ifdef BUILD_VTK_VERSION
vtkSmartPointer< vtkRenderer > RECONS3D::getRenderer( const int angios_ID ){
    if( angios_ID > n_angios ){
        char mensaje[512] = "\n<<Error: El angiograma XXX no ha sido agregado al reconstructor>>\n\n";
#if defined(_WIN32) || defined(_WIN64)
		sprintf_s(mensaje, 512 * sizeof(char), "\n<<Error: El angiograma %i no ha sido agregado al reconstructor>>\n\n", angios_ID);
#else
        sprintf(mensaje, "\n<<Error: El angiograma %i no ha sido agregado al reconstructor>>\n\n", angios_ID);
#endif
        DEB_MSG(mensaje);
        escribirLog( mensaje );
        return NULL;
    }else{
        return mis_renderers[ angios_ID ];
    }
}
#endif





/*  Metodo: setFiltroLog

    Funcion: Define el editor donde se escribiran todos los logs del filtrado
*/
void RECONS3D::setFiltroLog( FILE *fplog ){
    filtro.setLog( fplog );
}



/*  Metodo: setFiltroLog

    Funcion: Define el editor donde se escribiran todos los logs del filtrado
*/
void RECONS3D::setFiltroLog( const char *ruta_log ){
    filtro.setLog( ruta_log );
}





/*  Metodo: setProgressLog

    Funcion: Define el objeto donde se visualiza el progreso del proceso de filtrado.
*/
#ifdef BUILD_GUI_VERSION
void RECONS3D::setProgressBar(QProgressBar *pBar){
    mi_pBar = pBar;
}
#endif


/*  Metodo: setFiltroLog

    Funcion: Define el objeto donde se visualiza el progreso del proceso de filtrado.
*/
#ifdef BUILD_GUI_VERSION
void RECONS3D::setFiltroProgressBar(QProgressBar *pBar){
    filtro.setProgressBar( pBar );
}
#endif




/*  Metodo: umbralizar

    Funcion: Umbraliza la imagen definida con 'img_idx'
*/
void RECONS3D::umbralizar(const IMGCONT::THRESHOLD_ALG tipo_umb, const double nivel, const int angios_ID)
{
	IMGCONT *responses_threshold_temp = this->getIMGCONTPointer(IMG_THRESHOLD, angios_ID);
	*responses_threshold_temp = *this->getIMGCONTPointer(IMG_RESPONSE, angios_ID);

	responses_threshold_temp->threshold(tipo_umb, nivel);
}




// C O N S T R U C T O R E S    /   D E S T R U C T O R E S
/*  Constructor ()
    Funcion: Constructor por default.
*/
RECONS3D::RECONS3D(){
#ifdef BUILD_VTK_VERSION
    renderer_global = vtkSmartPointer<vtkRenderer>::New();
    double color[] = {1.0, 1.0, 1.0};
//    agregarEsfera(0.0, 0.0, 0.0, 10.0, color, renderer_global);
//    agregarEjes(renderer_global);
#endif

    fp_log = NULL;
    detalle = 180;

#ifdef BUILD_GUI_VERSION
    mi_txtLog = NULL;
    mi_pBar = NULL;
#endif

	imgs_base = new IMGCONT_LIST;
	imgs_base->my_first_node = NULL;
	imgs_base->my_last_node = NULL;
	imgs_base->my_nodes_count = 0;

	imgs_groundtruth = new IMGCONT_LIST;
	imgs_groundtruth->my_first_node = NULL;
	imgs_groundtruth->my_last_node = NULL;
	imgs_groundtruth->my_nodes_count = 0;

	imgs_mask = new IMGCONT_LIST;
	imgs_mask->my_first_node = NULL;
	imgs_mask->my_last_node = NULL;
	imgs_mask->my_nodes_count = 0;

	imgs_response = new IMGCONT_LIST;
	imgs_response->my_first_node = NULL;
	imgs_response->my_last_node = NULL;
	imgs_response->my_nodes_count = 0;

	imgs_response_threshold = new IMGCONT_LIST;
	imgs_response_threshold->my_first_node = NULL;
	imgs_response_threshold->my_last_node = NULL;
	imgs_response_threshold->my_nodes_count = 0;

	imgs_thrs_response_map_dists = new IMGCONT_LIST;
	imgs_thrs_response_map_dists->my_first_node = NULL;
	imgs_thrs_response_map_dists->my_last_node = NULL;
	imgs_thrs_response_map_dists->my_nodes_count = 0;

	imgs_thrs_response_skeleton = new IMGCONT_LIST;
	imgs_thrs_response_skeleton->my_first_node = NULL;
	imgs_thrs_response_skeleton->my_last_node = NULL;
	imgs_thrs_response_skeleton->my_nodes_count = 0;

	imgs_thrs_response_boundaries = new IMGCONT_LIST;
	imgs_thrs_response_boundaries->my_first_node = NULL;
	imgs_thrs_response_boundaries->my_last_node = NULL;
	imgs_thrs_response_boundaries->my_nodes_count = 0;
	
    n_angios = -1;
}







/************************************************************************************************************
*                                                                                                           *
* FUNCTION NAME: RECONS3D                                                                                   *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT        TYPE                I/O        DESCRIPTION                                                *
* input_arguments ARGUMENTS*          input      The arguments object                                       *
*                                                                                                           *
* RETURNS:                                                                                                  *
* Constructor for RECONS3D that defines the reconstruction arguments                                        *
*                                                                                                           *
************************************************************************************************************/
RECONS3D::RECONS3D(ARGUMENTS *input_arguments)
{
	my_args = input_arguments;
	defineArguments();
	
#ifdef BUILD_VTK_VERSION
	renderer_global = vtkSmartPointer<vtkRenderer>::New();
	double color[] = { 1.0, 1.0, 1.0 };
	//    agregarEsfera(0.0, 0.0, 0.0, 10.0, color, renderer_global);
	//    agregarEjes(renderer_global);
#endif

	fp_log = NULL;
	detalle = 180;

#ifdef BUILD_GUI_VERSION
	mi_txtLog = NULL;
	mi_pBar = NULL;
#endif


	imgs_base = new IMGCONT_LIST;
	imgs_base->my_first_node = NULL;
	imgs_base->my_last_node = NULL;
	imgs_base->my_nodes_count = 0;
	
	imgs_groundtruth = new IMGCONT_LIST;
	imgs_groundtruth->my_first_node = NULL;
	imgs_groundtruth->my_last_node = NULL;
	imgs_groundtruth->my_nodes_count = 0;
	
	imgs_mask = new IMGCONT_LIST;
	imgs_mask->my_first_node = NULL;
	imgs_mask->my_last_node = NULL;
	imgs_mask->my_nodes_count = 0;


	imgs_response = new IMGCONT_LIST;
	imgs_response->my_first_node = NULL;
	imgs_response->my_last_node = NULL;
	imgs_response->my_nodes_count = 0;

	imgs_response_threshold = new IMGCONT_LIST;
	imgs_response_threshold->my_first_node = NULL;
	imgs_response_threshold->my_last_node = NULL;
	imgs_response_threshold->my_nodes_count = 0;

	imgs_thrs_response_map_dists = new IMGCONT_LIST;
	imgs_thrs_response_map_dists->my_first_node = NULL;
	imgs_thrs_response_map_dists->my_last_node = NULL;
	imgs_thrs_response_map_dists->my_nodes_count = 0;

	imgs_thrs_response_skeleton = new IMGCONT_LIST;
	imgs_thrs_response_skeleton->my_first_node = NULL;
	imgs_thrs_response_skeleton->my_last_node = NULL;
	imgs_thrs_response_skeleton->my_nodes_count = 0;

	imgs_thrs_response_boundaries = new IMGCONT_LIST;
	imgs_thrs_response_boundaries->my_first_node = NULL;
	imgs_thrs_response_boundaries->my_last_node = NULL;
	imgs_thrs_response_boundaries->my_nodes_count = 0;


	n_angios = -1;
}



/*  Destructor
    Funcion: Libera la memoria utilizada para almacenar las imagenes.
*/
RECONS3D::~RECONS3D() {
	if (fp_log) {
		fclose(fp_log);
		fp_log = NULL;
	}
	
    for( std::vector<int*>::iterator it = hist.begin(); it != hist.end(); it++){
        if( *it ){
            delete *it;
        }
    }
}





void RECONS3D::defineArguments()
{
	my_args->newArgument(
		"Input angiography to be used as base image (.PNG, .BMP, .JPEG/.JPG, o DICOM file)",
		"b", "base", "NULL", true);

	my_args->newArgument(
		"Lower level to extract the images (from a DICOM file only)",
		"bl", "lower", 0, true);

	my_args->newArgument(
		"Upper level to extract the images (from a DICOM file only)",
		"bu", "upper", 0, true);

	my_args->newArgument(
		"Input ground-truth image for the base image (.PNG, .BMP, .JPEG/.JPG)",
		"g", "ground-truth", "NULL", true);
	
	my_args->newArgument(
		"Paths file to dataset base images used for training", "dsb", "dataset-base", "NULL", true);

	my_args->newArgument(
		"Paths file to dataset ground-truth images used for training", "dsg", "dataset-ground", "NULL", true);

	my_args->newArgument(
		"Path to save segmentation output image", "oseg", "output-segment", "segment.pgm", true);

	my_args->newArgument(
		"Path to save thresholding output image", "othr", "output-threshold", "threshold.pgm", true);

	my_args->newArgument(
		"Reconstruction process' log", "lr", "log-reconstruction", "NULL", true);

	my_args->newArgument(
		"Configuration file with the segmentation instructions", "c", "config", "NULL", true);

	my_args->newArgument(
		"Filtering process' log", "lf", "log-filter", "NULL", true);

	my_args->newArgument(
		"Treshold segmented image?", "thr", "threshold", "no", true);

}

//---- ----------------------------------------------------------------------------- PUBLIC----- ^

// C L A S E: RECONS3D  ------------------------------------------------------------------------ ^
