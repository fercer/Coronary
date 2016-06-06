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

/*  Metodo: escribirLog

    Funcion: Escribe un mensaje en el log.
*/
void RECONS3D::escribirLog( const char *mensaje ){
    if( mi_log ){
        mi_log->appendPlainText( mensaje );
    }else{
        std:cout << mensaje;
    }
}



/*  Metodo: renderizar()

    Funcion: Renderiza los 'actores' contenidos en 'renderer' en una ventana de VTK.
*/
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



/*  Metodo: isoCentro

    Funcion: Calcula el isocentro de la imagen 'angio_ID' como la linea que va desde el centro de la imagen al centro de la fuente.
*/
void RECONS3D::isoCentro( const int angio_ID ){
    // Calcular la ecuacion del plano:
    const double crl = cos(imgs_base[angio_ID].LAORAO/180.0 * PI);
    const double srl = sin(imgs_base[angio_ID].LAORAO/180.0 * PI);
    const double ccc = cos(imgs_base[angio_ID].CRACAU/180.0 * PI);
    const double scc = sin(imgs_base[angio_ID].CRACAU/180.0 * PI);

    const int mis_cols = imgs_base[angio_ID].cols;
    const int mis_rens = imgs_base[angio_ID].rows;

    const double mi_DDP = imgs_base[angio_ID].DDP;

    const double mi_pixX = imgs_base[angio_ID].pixX;
    const double mi_pixY = imgs_base[angio_ID].pixY;

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
    std::vector< NORCEN >::iterator mi_normal = normal_centros.begin() + angio_ID;
    mi_normal->direccion[0] = (Qy - Py)*(Rz - Pz) - (Qz - Pz)*(Ry - Py);
    mi_normal->direccion[1] = (Qz - Pz)*(Rx - Px) - (Qx - Px)*(Rz - Pz);
    mi_normal->direccion[2] = (Qx - Px)*(Ry - Py) - (Qy - Py)*(Rx - Px);

    /// Definir el centro de la imagen:
    mi_normal->origen[0] = Px;
    mi_normal->origen[1] = Py;
    mi_normal->origen[2] = Pz;
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

    const int mis_cols = imgs_base[angio_ID].cols;
    const int mis_rens = imgs_base[angio_ID].rows;
    const double mi_pixX = imgs_base[angio_ID].pixX;
    const double mi_pixY = imgs_base[angio_ID].pixY;

    const double mi_DDP = imgs_base[angio_ID].DDP;

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

            puntos[angio_ID]->InsertNextPoint(xx_3D, yy_3D, zz_3D);

            vtkSmartPointer< vtkVertex > pix = vtkSmartPointer< vtkVertex >::New();
            pix->GetPointIds()->SetId(0, x + y*mis_cols);

            pixeles[angio_ID]->InsertNextCell(pix);
            xx += mi_pixX;
        }
        yy += mi_pixY;
    }

}



/*  Metodo: mostrarImagen

    Funcion: Muestra la imagen en una ventana VTK.
*/
void RECONS3D::mostrarImagen(IMGVTK::IMG_IDX img_idx, vtkSmartPointer<vtkRenderer> mi_renderer, const int angio_ID){

    int mis_cols = imgs_base[angio_ID].cols;
    int mis_rens = imgs_base[angio_ID].rows;

    vtkSmartPointer< vtkImageData> img_ptr = NULL;

    switch(img_idx){
    case IMGVTK::BASE:
        img_ptr = imgs_base[angio_ID].base;
        break;
    case IMGVTK::GROUNDTRUTH:
        img_ptr = imgs_base[angio_ID].ground;
        break;
    case IMGVTK::MASK:
        img_ptr = imgs_base[angio_ID].mask;
        break;
    case IMGVTK::SKELETON:
        img_ptr = imgs_base[angio_ID].skeleton;
        mis_rens+=2;
        mis_cols+=2;
        break;
    case IMGVTK::SEGMENT:
        img_ptr = imgs_base[angio_ID].segment;
        break;
    case IMGVTK::THRESHOLD:
        img_ptr = imgs_base[angio_ID].threshold;
        break;
    case IMGVTK::MAPDIST:
        img_ptr = imgs_base[angio_ID].mapa_dist;
        break;
    case IMGVTK::BORDERS:
        img_ptr = imgs_base[angio_ID].borders;
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

    double max =-INF;
    double min = INF;

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
    actor->SetPosition(25,0);

    // Agregar el actor de la imagen al renderizador.
    mi_renderer->Clear();
    mi_renderer->AddActor2D(actor);
}





/*  Metodo: mostrarImagen

    Funcion: Muestra la imagen en una ventana VTK usando el renderizador global.
*/
void RECONS3D::mostrarImagen( const int angio_ID, IMGVTK::IMG_IDX img_idx){
    const int mis_cols = imgs_base[angio_ID].cols;
    const int mis_rens = imgs_base[angio_ID].rows;

    int offset_y = 0, offset_x = 0;
    double *img_ptr = NULL;

    switch(img_idx){
        case IMGVTK::BASE:
            img_ptr = imgs_base[angio_ID].base_ptr;
            break;
        case IMGVTK::MASK:
            img_ptr = imgs_base[angio_ID].mask_ptr;
            break;
        case IMGVTK::SKELETON:
            img_ptr = imgs_base[angio_ID].skl_ptr;
            offset_x=1;
            offset_y=1;
            break;
        case IMGVTK::SEGMENT:
            img_ptr = imgs_base[angio_ID].segment_ptr;
            break;
        case IMGVTK::THRESHOLD:
            img_ptr = imgs_base[angio_ID].threshold_ptr;
            break;
        case IMGVTK::MAPDIST:
            img_ptr = imgs_base[angio_ID].map_ptr;
            break;
        case IMGVTK::BORDERS:
                img_ptr = imgs_base[angio_ID].borders_ptr;
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
    rejilla->SetPoints( puntos[angio_ID] );
    rejilla->SetVerts( pixeles[angio_ID] );
    rejilla->GetCellData()->SetScalars(intensidades);

    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputData(rejilla);

    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    renderer_global->AddActor(actor);

    // Mostrar el isocentro:

    /// Posicionar una esfera en el lugar del isocentro:
    double norma = normal_centros[angio_ID].direccion[0]*normal_centros[angio_ID].direccion[0] +
                   normal_centros[angio_ID].direccion[1]*normal_centros[angio_ID].direccion[1] +
                   normal_centros[angio_ID].direccion[2]*normal_centros[angio_ID].direccion[2];

    norma = sqrt(norma);

    const double t_iso = imgs_base[angio_ID].DISO / norma;

    NORCEN mi_isocen;
    mi_isocen.origen[0] = normal_centros[angio_ID].origen[0];
    mi_isocen.origen[1] = normal_centros[angio_ID].origen[1];
    mi_isocen.origen[2] = normal_centros[angio_ID].origen[2];
    mi_isocen.direccion[0] = - t_iso * normal_centros[angio_ID].direccion[0];
    mi_isocen.direccion[1] = - t_iso * normal_centros[angio_ID].direccion[1];
    mi_isocen.direccion[2] = - t_iso * normal_centros[angio_ID].direccion[2];

    double color[] = {0.0, 0.0, 1.0};
    agregarEsfera(normal_centros[angio_ID].origen[0], normal_centros[angio_ID].origen[1], normal_centros[angio_ID].origen[2], 17.0, color, renderer_global);

    color[0] = 0.0; color[1] = 1.0; color[2] = 0.0;
    agregarVector(mi_isocen, 1.0, color, renderer_global);
    agregarEsfera(mi_isocen.origen[0] + mi_isocen.direccion[0], mi_isocen.origen[1] + mi_isocen.direccion[1], mi_isocen.origen[2] + mi_isocen.direccion[2], 20.0, color, renderer_global);
}



/*  Metodo: agregarVector

    Funcion: Agrega un vector al renderizador definido.
*/
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





/*  Metodo: agregarEjes

    Funcion: Muestra los ejes en que se encuentra la escena global respecto al paciente.
*/
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

DEB_MSG("Mostrando detector para: " << angio_ID << ", LAORAO: " << imgs_base[angio_ID].LAORAO << ", CRACAU: " << imgs_base[angio_ID].CRACAU << ", SID: " << imgs_base[angio_ID].SID << ", SOD: " << imgs_base[angio_ID].SOD );

    const int mis_rens = imgs_base[n_angios].rows;
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



//------------------------------------------------------------------------------------------------ PRIVATE ----- ^
//-------------------------------------------------------------------------------------------------- PUBLIC----- v

// M I E M B R O S      P U B L I C O S


/*  Metodo: umbralizar

    Funcion:
*/
void RECONS3D::umbralizar( const int angio_ID, const IMGVTK::TIPO_UMBRAL mi_umbral, const double umbral ){
    imgs_base[ angio_ID ].umbralizar( IMGVTK::SEGMENT, mi_umbral, umbral );
}


/*  Metodo: agregarInput

    Funcion: Define las rutas de las imagenes que son usadas para reconstruir una arteria coronaria.
*/
void RECONS3D::agregarInput(const char *rutabase_input, const int nivel_l, const int nivel_u, const char *rutaground_input){

DEB_MSG("Ruta ground: " << rutaground_input);
    NORCEN norcen_temp;

    if( nivel_u > nivel_l ){ // Si hay varios niveles:

        for( int i = nivel_l; i <= nivel_u; i++){
            n_angios++;

            imgs_base.push_back(IMGVTK(rutabase_input, true, i));
            mis_renderers.push_back(vtkSmartPointer<vtkRenderer>::New());
            puntos.push_back(vtkSmartPointer<vtkPoints>::New());
            pixeles.push_back(vtkSmartPointer<vtkCellArray>::New());
            normal_centros.push_back( norcen_temp );

            /// Mover el detector a su posicion definida por el archivo DICOM:
            mallarPuntos(n_angios);
            isoCentro(n_angios);
            existe_ground.push_back( false );
        }

    }else{
        n_angios++;

        imgs_base.push_back(IMGVTK(rutabase_input, true, nivel_l));
        imgs_base[n_angios].setLog( mi_log );
        existe_ground.push_back( false );

        mis_renderers.push_back(vtkSmartPointer<vtkRenderer>::New());
        puntos.push_back(vtkSmartPointer<vtkPoints>::New());
        pixeles.push_back(vtkSmartPointer<vtkCellArray>::New());
        normal_centros.push_back( norcen_temp );


        /// Mover el detector a su posicion definida por el archivo DICOM:
        mallarPuntos(n_angios);
        isoCentro(n_angios);

        if( strcmp(rutaground_input, "NULL") ){
            DEB_MSG("Abriendo " << rutaground_input << " como ground truth");
            imgs_base[ n_angios ].Cargar(IMGVTK::GROUNDTRUTH, rutaground_input, false, 0);
            existe_ground[ n_angios ] = true;
        }

    }
}



/*  Metodo: agregarInput

    Funcion: Define las rutas de las imagenes que son usadas para reconstruir una arteria coronaria.
*/
void RECONS3D::agregarInput( const char *rutabase_input ){

    n_angios++;

    NORCEN norcen_temp;

    imgs_base.push_back(IMGVTK(rutabase_input, true, 0));
    imgs_base[n_angios].setLog( mi_log );
    existe_ground.push_back( false );

    mis_renderers.push_back(vtkSmartPointer<vtkRenderer>::New());
    puntos.push_back(vtkSmartPointer<vtkPoints>::New());
    pixeles.push_back(vtkSmartPointer<vtkCellArray>::New());
    normal_centros.push_back( norcen_temp );

    /// Mover el detector a su posicion definida por el archivo DICOM:
    mallarPuntos(n_angios);
    isoCentro(n_angios);

}





/*  Metodo: agregarInput

    Funcion: Define las rutas de las imagenes que son usadas para reconstruir una arteria coronaria.
*/
void RECONS3D::agregarInput( char **rutasbase_input, const int n_imgs){


    n_angios++;

    NORCEN norcen_temp;

    imgs_base.push_back(IMGVTK(rutasbase_input, n_imgs, true));
    imgs_base[n_angios].setLog( mi_log );
    existe_ground.push_back( false );

    // Mostrar la imagen en un renderizador
    mis_renderers.push_back(vtkSmartPointer<vtkRenderer>::New());

    puntos.push_back(vtkSmartPointer<vtkPoints>::New());
    pixeles.push_back(vtkSmartPointer<vtkCellArray>::New());
    normal_centros.push_back( norcen_temp );

    /// Mover el detector a su posicion definida por el archivo DICOM:
    mallarPuntos(n_angios);
    isoCentro(n_angios);

}





/*  Metodo: agregarGroundtruth

    Funcion: Define las rutas de las imagenes que son usadas para reconstruir una arteria coronaria.
*/
void RECONS3D::agregarGroundtruth(const char *rutaground_input, const int angio_ID ){
    if( angio_ID > n_angios ){
        char mensaje[] = "\n<<Error: El angiograma XXX no ha sido agregado al reconstructor>>\n\n";
        sprintf(mensaje, "\n<<Error: El angiograma %i no ha sido agregado al reconstructor>>\n\n", angio_ID);
        escribirLog( mensaje );
    }else{
        imgs_base[ angio_ID ].Cargar(IMGVTK::GROUNDTRUTH, rutaground_input, false, 0);
        existe_ground[ angio_ID ] = true;
    }
}





/*  Metodo: agregarGroundtruth

    Funcion: Define las rutas de las imagenes que son usadas para reconstruir una arteria coronaria.
*/
void RECONS3D::agregarGroundtruth(char **rutasground_input, const int n_imgs, const int angio_ID ){
    if( angio_ID > n_angios ){
        char mensaje[] = "\n<<Error: El angiograma XXX no ha sido agregado al reconstructor>>\n\n";
        sprintf(mensaje, "\n<<Error: El angiograma %i no ha sido agregado al reconstructor>>\n\n", angio_ID);
        escribirLog( mensaje );
    }else{
        imgs_base[ angio_ID ].Cargar(IMGVTK::GROUNDTRUTH, rutasground_input, n_imgs, false);
        existe_ground[ angio_ID ] = true;
    }
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




/*  Metodo: getRows

    Funcion:
*/
int RECONS3D::getRows( const int angio_ID ){
    return imgs_base[angio_ID].rows;
}



/*  Metodo: getCols

    Funcion:
*/
int RECONS3D::getCols( const int angio_ID ){
    return imgs_base[angio_ID].cols;
}




/*  Metodo: get_pixelData

    Funcion: Retorna el apuntador a la informacion de la imagen seleccionada con img_idx
*/
double* RECONS3D::get_pixelData(const int angio_ID, IMGVTK::IMG_IDX img_idx){

    double *img_ptr = NULL;

    switch(img_idx){
    case IMGVTK::BASE:
        img_ptr = imgs_base[angio_ID].base_ptr;
        break;
    case IMGVTK::GROUNDTRUTH:
        img_ptr = imgs_base[angio_ID].gt_ptr;
        break;
    case IMGVTK::MASK:
        img_ptr = imgs_base[angio_ID].mask_ptr;
        break;
    case IMGVTK::SKELETON:
        img_ptr = imgs_base[angio_ID].skl_ptr;
        break;
    case IMGVTK::SEGMENT:
        img_ptr = imgs_base[angio_ID].segment_ptr;
        break;
    case IMGVTK::THRESHOLD:
        img_ptr = imgs_base[angio_ID].threshold_ptr;
        break;
    case IMGVTK::MAPDIST:
        img_ptr = imgs_base[angio_ID].map_ptr;
        break;
    case IMGVTK::BORDERS:
        img_ptr = imgs_base[angio_ID].borders_ptr;
        break;
    }

    return img_ptr;
}


/*  Metodo: segmentarImagen

    Funcion: Aplica el filtro con los parametros definidos
*/
void RECONS3D::segmentarImagenBase( const int angio_ID ){

    filtro.setInput(imgs_base[angio_ID]);

    if( filtro.getParametrosOptimizar() >= 1 ){
        filtro.setPar();
    }

    filtro.filtrar();
    escribirLog( "\nEl filtrado de la imagen base termino con exito\n" );
}



/*  Metodo: mostrarRadios

    Funcion: Muestra el radio de cada pixel del esqueleto.
*/
void RECONS3D::mostrarRadios(vtkSmartPointer<vtkPoints> puntos, vtkSmartPointer<vtkCellArray> cilindros, int *n_pix, IMGVTK::PIX_PAR *grafo, const double DDP, const double crl, const double srl, const double ccc, const double scc, const int nivel_detalle, FILE *fp_cilindros){

    const double theta_inc = 2 * PI / (double)detalle;

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

        const double alpha = atan2(y_fin - yy, x_fin - xx) + PI / 2.0;

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








/*  Metodo: mostrarRadios

    Funcion: Muestra el radio de cada pixel del esqueleto.
*/
void RECONS3D::mostrarRadios(vtkSmartPointer<vtkPoints> puntos, vtkSmartPointer<vtkCellArray> vert_skl, vtkSmartPointer<vtkUnsignedCharArray> grafo_nivel, int *n_pix, IMGVTK::PIX_PAR *grafo, const double DDP, const double crl, const double srl, const double ccc, const double scc, const int n_niveles){

    unsigned char nivel[3];

    nivel[0] = (unsigned char)(255.0 * (double)(grafo->nivel + 1) / (double)n_niveles);
    nivel[1] = (unsigned char)(255.0 * (double)(grafo->nivel + 1) / (double)n_niveles);
    nivel[2] = (unsigned char)(255.0 * (double)(grafo->nivel + 1) / (double)n_niveles);

    double theta = 0.0;
    const double theta_inc = 2 * PI / (double)detalle;
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




/*  Metodo: skeletonize

    Funcion: Obtiene el esqueleto de la imagen y muestra los puntos de interes.
*/
void RECONS3D::skeletonize(const int angio_ID){

    DEB_MSG("Extrayendo esquelto a " << angio_ID);
    imgs_base[angio_ID].skeletonization(IMGVTK::GROUNDTRUTH);


    //if( !imgs_base[angio_ID].pix_caract ){
    if( !imgs_base[angio_ID].pix_caract ){
        DEB_MSG("No existe grafo alguno...");
        return;
    }

    // Rotar los puntos segun LAO/RAO y CAU/CRA:
    const double crl = cos(imgs_base[angio_ID].LAORAO/180.0 * PI);
    const double srl = sin(imgs_base[angio_ID].LAORAO/180.0 * PI);
    const double ccc = cos(imgs_base[angio_ID].CRACAU/180.0 * PI);
    const double scc = sin(imgs_base[angio_ID].CRACAU/180.0 * PI);

    const double DDP = 0.0;//imgs_base[angio_ID].DDP;
    const int n_niveles = imgs_base[angio_ID].n_niveles;

    DEB_MSG("numero de niveles: " << n_niveles);

    {
        // Recorrer el grafo y generar en 3D una burda reconstruccion.

        DEB_MSG("Mostrando los radios de cada 100 pixeles del esqueleto...");
        vtkSmartPointer< vtkPoints > puntos = vtkSmartPointer< vtkPoints >::New();
        vtkSmartPointer< vtkCellArray > cilindros = vtkSmartPointer< vtkCellArray >::New();
        int n_pix = 0;

        char nom_cilindros[] = "cilindros_XXX.dat";
        sprintf( nom_cilindros, "cilindros_%i.dat", angio_ID);
        FILE *fp_cilindros = fopen(nom_cilindros, "w");
        fprintf( fp_cilindros, "X1 Y1 Z1 RADIO X2 Y2 Z2\n");
        mostrarRadios(puntos, cilindros, &n_pix, imgs_base[angio_ID].pix_caract, DDP, crl, srl, ccc, scc, 100, fp_cilindros);

        fclose( fp_cilindros );

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
    }
}





/*  Metodo: mostrarBase

    Funcion: Muestra la imagen base en el renderizador del input 'angio_ID'
*/
void RECONS3D::mostrarBase( const int angio_ID ){
    if( angio_ID > n_angios ){
        char mensaje_err[] = "\n<<Error: El angiograma XXX no ha sido agregado al reconstructor>>\n\n";
        sprintf(mensaje_err, "\n<<Error: El angiograma %i no ha sido agregado al reconstructor>>\n\n", angio_ID);
        DEB_MSG(mensaje_err);
        escribirLog( mensaje_err );
    }else{
        mostrarImagen( IMGVTK::BASE, mis_renderers[ angio_ID ], angio_ID);
        renderizar( mis_renderers[angio_ID] );
    }
}





/*  Metodo: mostrarGroundtruth

    Funcion: Muestra la imagen ground-truth en el renderizador del input 'angio_ID'
*/
void RECONS3D::mostrarGroundtruth( const int angio_ID ){
    if( (angio_ID > n_angios) || (!existe_ground.at(angio_ID)) ){
        char mensaje[] = "\n<<Error: El ground-truth para el angiograma XXX no ha sido agregado al reconstructor>>\n\n";
        sprintf(mensaje, "\n<<Error: El ground-truth para el angiograma %i no ha sido agregado al reconstructor>>\n\n", angio_ID);
        DEB_MSG(mensaje);
        escribirLog( mensaje );
    }else{
        mostrarImagen( IMGVTK::GROUNDTRUTH, mis_renderers[ angio_ID ], angio_ID);
    }
}








/*  Metodo: getRenderer

    Funcion: Retorna el renderizador global.
*/
vtkSmartPointer< vtkRenderer > RECONS3D::getRenderer(){
    return renderer_global;
}




/*  Metodo: getRenderer

    Funcion: Retorna el renderizador correspondiente al input 'angio_ID'.
*/
vtkSmartPointer< vtkRenderer > RECONS3D::getRenderer( const int angio_ID ){
    if( angio_ID > n_angios ){
        char mensaje[] = "\n<<Error: El angiograma XXX no ha sido agregado al reconstructor>>\n\n";
        sprintf(mensaje, "\n<<Error: El angiograma %i no ha sido agregado al reconstructor>>\n\n", angio_ID);
        DEB_MSG(mensaje);
        escribirLog( mensaje );
        return NULL;
    }else{
        return mis_renderers[ angio_ID ];
    }
}



/*  Metodo: setLog

    Funcion: Define el editor donde se escribiran todos los logs del sistema
*/
void RECONS3D::setLog(QPlainTextEdit *log ){
    mi_log = log;

    for( std::vector< IMGVTK >::iterator it = imgs_base.begin(); it != imgs_base.end(); it++){
        it->setLog( mi_log );
    }

    escribirLog( "Reconstructor 3D" );
    escribirLog( "by: Fernando Cervantes-Sanchez\n" );
}




// C O N S T R U C T O R E S    /   D E S T R U C T O R E S
/*  Constructor ()
    Funcion: Constructor por default.
*/
RECONS3D::RECONS3D(){
    mi_log = NULL;
    renderer_global = vtkSmartPointer<vtkRenderer>::New();

    detalle = 180;

    double color[] = {1.0, 1.0, 1.0};
//    agregarEsfera(0.0, 0.0, 0.0, 10.0, color, renderer_global);
//    agregarEjes(renderer_global);

    n_angios = -1;
}



/*  Destructor
    Funcion: Libera la memoria utilizada para almacenar las imagenes.
*/
RECONS3D::~RECONS3D(){

}

//-------------------------------------------------------------------------------------------------- PUBLIC----- ^

// C L A S E: RECONS3D  ------------------------------------------------------------------------ ^
