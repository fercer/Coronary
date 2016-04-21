﻿/******
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




/*  Metodo: isoCentro

    Funcion: Calcula un isocentro como la media de los isocentros.
*/
void RECONS3D::isoCentro( int *angios_ID ){

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
    const int mis_rens = imgs_base[angio_ID].rens;

    const double mi_DDP = imgs_base[angio_ID].DDP;

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
    xx = mis_cols;
    double Qx, Qy, Qz;

    //// Rotacion usando el eje 'x' como base:
    Qy = crl*yy - srl*mi_DDP;
    Qz = srl*yy + crl*mi_DDP;
    //// Rotacion usando el eje 'y' como base:
    Qx = ccc*xx - scc*Qz;
    Qz = scc*xx + ccc*Qz;



    // Tercer punto (0, mis_rens, DDP) ------------------------------------------------------------------------------
    xx = 0.0;
    yy = mis_rens;
    double Rx, Ry, Rz;

    //// Rotacion usando el eje 'x' como base:
    Ry = crl*yy - srl*mi_DDP;
    Rz = srl*yy + crl*mi_DDP;
    //// Rotacion usando el eje 'y' como base:
    Rx = ccc*xx - scc*Rz;
    Rz = scc*xx + ccc*Rz;


    //// Calcular la normal al plano:
    std::vector< NORCEN >::iterator mi_normal = normal_centros.begin() + angio_ID;
    mi_normal->Nx = (Qy - Py)*(Rz - Pz) - (Qz - Pz)*(Ry - Py);
    mi_normal->Ny = (Qx - Px)*(Rz - Pz) - (Qz - Pz)*(Rx - Px);
    mi_normal->Nz = (Qx - Px)*(Ry - Py) - (Qy - Py)*(Rx - Px);

    /// Definir el centro de la imagen:
    const double mi_pixX = imgs_base[angio_ID].pixX;
    const double mi_pixY = imgs_base[angio_ID].pixY;
    xx = (mis_cols+1)*mi_pixX / 2;
    yy = (mis_rens+1)*mi_pixY / 2;

    double Cx, Cy, Cz;

    //// Rotacion usando el eje 'x' como base:
    Cy = crl*yy - srl*mi_DDP;
    Cz = srl*yy + crl*mi_DDP;
    //// Rotacion usando el eje 'y' como base:
    Cx = ccc*xx - scc*Cz;
    Cz = scc*xx + ccc*Cz;

    mi_normal->Cx = Cx;
    mi_normal->Cy = Cy;
    mi_normal->Cz = Cz;
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
    const int mis_rens = imgs_base[angio_ID].rens;
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
            pix->GetPointIds()->SetId(0, xx + yy*mis_cols);

            pixeles[angio_ID]->InsertNextCell(pix);
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
        case IMGVTK::THRESHOLD:
            img_ptr = img_src.threshold;
            break;
    }

DEB_MSG("("<< img_idx <<"/" << IMGVTK::BASE << ") img_ptr: " << img_ptr);
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
    TIMERS;
    GETTIME_INI;
    const int mis_cols = imgs_base[angio_ID].cols;
    const int mis_rens = imgs_base[angio_ID].rens;
    const int mi_pixX = imgs_base[angio_ID].pixX;

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
    }

    vtkSmartPointer<vtkDoubleArray> intensidades = vtkSmartPointer<vtkDoubleArray>::New();

    intensidades->SetNumberOfComponents(1);
    intensidades->SetName("Intensidadsdasdsads");

    for( int i = 0; i < mis_rens; i++){
        for(int j = 0; j < mis_cols; j++){
            intensidades->InsertNextValue(img_ptr[(i+offset_y)*(mis_cols+offset_x*2) + (j+offset_x)]);
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
    const double ICx = normal_centros[angio_ID].Cx + imgs_base[angio_ID].DISO * normal_centros[angio_ID].Nx;
    const double ICy = normal_centros[angio_ID].Cy + imgs_base[angio_ID].DISO * normal_centros[angio_ID].Ny;
    const double ICz = normal_centros[angio_ID].Cz + imgs_base[angio_ID].DISO * normal_centros[angio_ID].Nz;
    GETTIME_FIN;

DEB_MSG("Tiempo en segundos para extraer un dato desde el vector de IMGVTK's: " << DIFTIME);
    // Setup render window, renderer, and interactor
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

DEB_MSG("Mostrando detector para: " << angio_ID << ", LAORAO: " << imgs_base[angio_ID].LAORAO << ", CRACAU: " << imgs_base[angio_ID].CRACAU << ", SID: " << imgs_base[angio_ID].SID << ", SOD: " << imgs_base[angio_ID].SOD );

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
void RECONS3D::agregarInput(const char *rutabase_input, const int nivel_l, const int nivel_u, const char *rutaground_input){

DEB_MSG("Ruta ground: " << rutaground_input);
    NORCEN norcen_temp;

    if( nivel_u > nivel_l ){ // Si hay varios niveles:

        for( int i = nivel_l; i <= nivel_u; i++){
            imgs_base.push_back(IMGVTK(rutabase_input, true, i));
            mis_renderers.push_back(vtkSmartPointer<vtkRenderer>::New());
            puntos.push_back(vtkSmartPointer<vtkPoints>::New());
            pixeles.push_back(vtkSmartPointer<vtkCellArray>::New());
            normal_centros.push_back( norcen_temp );

            // Mover el detector a su posicion definida por el archivo DICOM:
            mallarPuntos(n_angios);
            isoCentro(n_angios);

            mostrarImagen(imgs_base[n_angios], IMGVTK::BASE, mis_renderers[n_angios]);
            renderizar(mis_renderers[n_angios]);

            existe_ground.push_back( false );
            n_angios++;
        }

    }else{

        imgs_base.push_back(IMGVTK(rutabase_input, true, nivel_l));
        mis_renderers.push_back(vtkSmartPointer<vtkRenderer>::New());
        puntos.push_back(vtkSmartPointer<vtkPoints>::New());
        pixeles.push_back(vtkSmartPointer<vtkCellArray>::New());
        normal_centros.push_back( norcen_temp );


        // Mover el detector a su posicion definida por el archivo DICOM:
        mallarPuntos(n_angios);
        isoCentro(n_angios);

        mostrarImagen( n_angios, IMGVTK::BASE);
        renderizar(renderer_global);
//        mostrarImagen(imgs_base[n_angios], IMGVTK::BASE, mis_renderers[n_angios]);
//        renderizar(mis_renderers[n_angios]);

        if( strcmp(rutaground_input, "NULL") ){
            imgs_delin.push_back(IMGVTK(rutaground_input, false, 0));
            existe_ground.push_back( true );
        }else{
            existe_ground.push_back( false );
        }
        n_angios++;
    }
}




/*  Metodo: segmentarImagen()

    Funcion: Aplica el filtro a todas las imagenes.
*/
void RECONS3D::segmentarImagenBase(){
    for( int i = 0; i < n_angios; i++){
        segmentarImagenBase( i );
//        mostrarImagen(i, IMGVTK::SEGMENT);
    }
//    renderizar(renderer_global);
}




/*  Metodo: segmentarImagen()

    Funcion: Aplica el filtro de Gabor con los parametros definidos como mejors( tras la optimizacion con un algoritmo de evolucion computacional)
*/
void RECONS3D::segmentarImagenBase( const int angio_ID ){
    FILTROS filtro;
    filtro.setFiltro(FILTROS::SS_GABOR);
//    filtro.setFiltro(FILTROS::GMF);
    filtro.setFitness(FILTROS::ROC);
    filtro.setEvoMet(FILTROS::EDA_BUMDA, 50, 30);

    filtro.setInput(imgs_base[angio_ID]);
    if( existe_ground[angio_ID] ){
        filtro.setInputGround(imgs_delin[angio_ID]);
    }

    // Parametros fijos (SS Gabor):
    filtro.setPar(FILTROS::PAR_L, 2.65);
    filtro.setPar(FILTROS::PAR_T, 15);
    filtro.setPar(FILTROS::PAR_K, 180);
    filtro.setPar(FILTROS::PAR_DELTA, 1e-4);

    // Parametros fijos (GMF):
//    filtro.setPar(FILTROS::PAR_L, 13);
//    filtro.setPar(FILTROS::PAR_T, 15);
//    filtro.setPar(FILTROS::PAR_K, 12);
//    filtro.setPar(FILTROS::PAR_SIGMA, 2.82);
//    filtro.setPar(FILTROS::PAR_DELTA, 1e-4);
    filtro.filtrar();

    mostrarImagen(imgs_base[angio_ID], IMGVTK::SEGMENT, mis_renderers[angio_ID]);
    renderizar(mis_renderers[angio_ID]);


    imgs_base[angio_ID].umbralizar(IMGVTK::SEGMENT);

    mostrarImagen(imgs_base[angio_ID], IMGVTK::THRESHOLD, mis_renderers[angio_ID]);
    renderizar(mis_renderers[angio_ID]);

    /// Falta el linking broken vessels ..................
    ///
    imgs_base[angio_ID].lengthFilter(IMGVTK::THRESHOLD, 1000);
    //imgs_base[angio_ID].lengthFilter(IMGVTK::THRESHOLD, imgs_base[angio_ID].cols * 5);
    imgs_base[angio_ID].regionFill(IMGVTK::THRESHOLD);
    mostrarImagen(imgs_base[angio_ID], IMGVTK::THRESHOLD, mis_renderers[angio_ID]);
    renderizar(mis_renderers[angio_ID]);
}




/*  Metodo: skeletonize

    Funcion: Obtiene el esqueleto de todas las imagenes:
*/
void RECONS3D::skeletonize(){
    for( int i = 0; i < n_angios; i++){
        skeletonize( i );
    }
    renderizar(renderer_global);
}



/*  Metodo: skeletonize

    Funcion: Obtiene el esqueleto de la imagen y muestra los puntos de interes.
*/
void RECONS3D::skeletonize(const int angio_ID){
    imgs_base[angio_ID].skeletonization(IMGVTK::THRESHOLD);

    int n_caracts = imgs_base[angio_ID].n_caracts;

    if( !n_caracts ){
        return;
    }

    vtkSmartPointer<vtkCellArray> verticesSkeleton = vtkSmartPointer<vtkCellArray>::New();

    const double mis_cols = imgs_base[angio_ID].cols;
    const double mis_rens = imgs_base[angio_ID].rens;

    for( int c = 0; c < n_caracts; c++ ){
        const int xx = imgs_base[angio_ID].pix_caract[c].x - 1;
        const int yy = imgs_base[angio_ID].pix_caract[c].y - 1;

        vtkSmartPointer<vtkVertex> pix = vtkSmartPointer<vtkVertex>::New();
        pix->GetPointIds()->SetId(0,  (xx + yy*mis_cols) );

        verticesSkeleton->InsertNextCell( pix );
    }

    vtkSmartPointer<vtkPolyData> polyDataSkeleton = vtkSmartPointer<vtkPolyData>::New();
    polyDataSkeleton->SetPoints( puntos[angio_ID] );
    polyDataSkeleton->SetVerts( verticesSkeleton );

    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputData(polyDataSkeleton);

    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);

    double color[] = {1.0, 0.0, 0.0};
    actor->GetProperty()->SetColor(color);

    renderer_global->AddActor( actor );
}



// C O N S T R U C T O R E S    /   D E S T R U C T O R E S
/*  Constructor ()
    Funcion: Constructor por default.
*/
RECONS3D::RECONS3D(){
    renderer_global = vtkSmartPointer<vtkRenderer>::New();

    double color[] = {1.0, 1.0, 1.0};
    agregarEsfera(0.0, 0.0, 0.0, 5.0, color, renderer_global);
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
