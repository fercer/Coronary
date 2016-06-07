/******
    CENTRO DE INVESTIGACION EN MATEMATICAS
    MAESTRIA EN COMPUTACION Y MATEMATICAS INDUSTRIALES

    FERNANDO CERVANTES SANCHEZ
    ### - 20##
*****/


#include "IMGVTK.h"

// C L A S E: IMGVTK  ---------------------------------------------------------------------------------------------- v
//----------------------------------------------------------------------------- PRIVATE ------- v
    // M E T O D O S      P R I V A D O S


/*  Metodo: escribirLog

    Funcion: Escribe un mensaje en el log.
*/
void IMGVTK::escribirLog( const char *mensaje ){
    DEB_MSG( mi_log );
    if( mi_log ){
        mi_log->appendPlainText( mensaje );
    }else{
        std::cout << mensaje;
    }
}



//****************************************************************************************************************************************
//                              H E R R A M I E N T A S     D E     P R O C E S A M I E N T O       D E     I M A G E N E S
//****************************************************************************************************************************************
/*  Metodo: mapaDistancias
    Funcion: Obtiene el mapa de distancias de los pixeles a los bordes.
*/
void IMGVTK::mapaDistancias( IMG_IDX img_idx ){
    double *img_ptr = NULL;

    switch( img_idx ){
    case BASE:
        img_ptr = base_ptr;
        break;
    case GROUNDTRUTH:
        img_ptr = gt_ptr;
        break;
    case MASK:
        img_ptr = mask_ptr;
        break;
    case SKELETON:
        img_ptr = skl_ptr;
        break;
    case SEGMENT:
        img_ptr = segment_ptr;
        break;
    case THRESHOLD:
        img_ptr = threshold_ptr;
        break;
    }


    if( !map_ptr ){
        mapa_dist = vtkSmartPointer<vtkImageData>::New();
        mapa_dist->SetExtent(0, cols-1, 0, rows-1, 0, 0);
        mapa_dist->AllocateScalars( VTK_DOUBLE, 1);
        mapa_dist->SetOrigin(0.0, 0.0, 0.0);
        mapa_dist->SetSpacing(1.0, 1.0, 1.0);

        map_ptr = static_cast<double*>(mapa_dist->GetScalarPointer(0,0,0));
    }

    for( int xy = 0 ; xy < rows_cols; xy++ ){
        *(map_ptr + xy) = (*(img_ptr + xy) < 1.0) ? 0.0 : INF;
    }

    double *f = new double[ rows > cols ? rows : cols  ];
    double *dh = new double[ rows ];
    int *vh= new int[ rows ];
    double *zh = new double[ rows+1 ];

    // transform along columns
    for (int x = 0; x < cols; x++){
        for (int y = 0; y < rows; y++){
            f[y] = *(map_ptr + y*cols + x);
        }

        int k = 0;
        vh[0] = 0;
        zh[0] = -INF;
        zh[1] = +INF;

        for (int q = 1; q < rows; q++){
            double s  = ((f[q]+(q*q))-(f[vh[k]]+(vh[k]*vh[k])))/(2*q-2*vh[k]);

            while(s <= zh[k]){
                k--;
                s  = ((f[q]+(q*q))-(f[vh[k]]+(vh[k]*vh[k])))/(2*q-2*vh[k]);
            }

            k++;
            vh[k] = q;
            zh[k] = s;
            zh[k+1] = +INF;
        }

        k = 0;
        for (int y = 0; y < rows; y++){
            while (zh[k+1] < y){
                k++;
            }
            *(map_ptr + y*cols + x) = ((y-vh[k])*(y-vh[k])) + f[vh[k]];
        }

    }
    delete [] dh;
    delete [] vh;
    delete [] zh;

    double *dw = new double[cols];
    int *vw = new int[cols];
    double *zw = new double[cols+1];

    // transform along rows
    for (int y = 0; y < rows; y++){
        for (int x = 0; x < cols; x++){
            f[x] = *(map_ptr + y*cols + x);
        }
        int k = 0;
        vw[0] = 0;
        zw[0] = -INF;
        zw[1] = +INF;

        for (int q = 1; q < cols; q++){
            double s  = ((f[q]+(q*q))-(f[vw[k]]+(vw[k]*vw[k])))/(2*q-2*vw[k]);
            while (s <= zw[k]){
                k--;
                s  = ((f[q]+(q*q))-(f[vw[k]]+(vw[k]*vw[k])))/(2*q-2*vw[k]);
            }
            k++;
            vw[k] = q;
            zw[k] = s;
            zw[k+1] = +INF;
        }

        k = 0;
        for( int x = 0; x < cols; x++){
            while (zw[k+1] < x){
                k++;
            }
            map_ptr[y*cols + x] = sqrt( ((x-vw[k])*(x-vw[k])) + f[vw[k]] );

            if( (int)(*(map_ptr + y*cols + x) + 1) > max_dist){
                max_dist = (int)(*(map_ptr + y*cols + x) + 1);
            }
        }
    }

    DEB_MSG("Max dist: " << max_dist);


    delete [] dw;
    delete [] vw;
    delete [] zw;
    delete [] f;
}



/*  Metodo: detectarBorde
	Funcion: Detecta los bordes de la imagen en base a la transformada de la distancia.
*/
void IMGVTK::detectarBorde( IMG_IDX img_idx ){
    if (!map_ptr) {
        mapaDistancias( img_idx );
	}

	if (!borders) {
		borders = vtkSmartPointer<vtkImageData>::New();

        borders->SetExtent(0, cols - 1, 0, rows - 1, 0, 0);
        borders->AllocateScalars(VTK_DOUBLE, 1);
        borders->SetOrigin(0.0, 0.0, 0.0);
        borders->SetSpacing(1.0, 1.0, 1.0);

        borders_ptr = static_cast<double*>(borders->GetScalarPointer(0, 0, 0));
    }

    memset(borders_ptr, 0, rows_cols*sizeof(double));

	for (int xy = 0; xy < rows_cols; xy++) {
		if ( (*(base_ptr + xy) > 0.5) && (*(map_ptr + xy) < 2.0) ) {
			*(borders_ptr + xy) = 1.0;
		}
	}
}



/*  Metodo: regionFilling9
    Funcion: Rellena espacios vacios dentro del conjunto de pixeles resaltado.
*/
bool IMGVTK::regionFilling9( const double *ptr, const int x, const int y, const int mis_cols, const int mis_rens){
    int n_hits = 0;

    for( int m = 0; m < 9; m++){
        // Arriba:
        if( (y-4) > 0 ){
            if( ((x-4+m) > 0) && ((x-4+m) < mis_cols) ){
                n_hits += (ptr[(x-4+m) + (y-4)*mis_cols] > 0);
            }else{
                n_hits++;
            }
        }else{
            n_hits++;
        }

        // Abajo:
        if( (y+4) < mis_rens ){
            if( ((x-4+m) > 0) && ((x-4+m) < mis_cols) ){
                n_hits += (ptr[(x-4+m) + (y+4)*mis_cols] > 0);
            }else{
                n_hits++;
            }
        }else{
            n_hits++;
        }

        // Izquierda:
        if( (x-4) > 0 ){
            if( ((y-4+m) > 0) && ((y-4+m) < mis_rens) ){
                n_hits += (ptr[(x-4) + (y-4+m)*mis_cols] > 0);
            }else{
                n_hits++;
            }
        }else{
            n_hits++;
        }

        // Derecha:
        if( (x+4) < mis_cols ){
            if( ((y-4+m) > 0) && ((y-4+m) < mis_rens) ){
                n_hits += (ptr[(x+4) + (y-4+m)*mis_cols] > 0);
            }else{
                n_hits++;
            }
        }else{
            n_hits++;
        }
    }
    return (n_hits == 36);
}



/*  Metodo: regionFilling7
    Funcion: Rellena espacios vacios dentro del conjunto de pixeles resaltado.
*/
bool IMGVTK::regionFilling7( const double*ptr, const int x, const int y, const int mis_cols, const int mis_rens){

    int n_hits = 0;

    for( int m = 0; m < 7; m++){
        // Arriba:
        if( (y-3) > 0 ){
            if( ((x-3+m) > 0) && ((x-3+m) < mis_cols) ){
                n_hits += (ptr[(x-3+m) + (y-3)*mis_cols] > 0);
            }else{
                n_hits++;
            }
        }else{
            n_hits++;
        }

        // Abajo:
        if( (y+3) < mis_rens ){
            if( ((x-3+m) > 0) && ((x-3+m) < mis_cols) ){
                n_hits += (ptr[(x-3+m) + (y+3)*mis_cols] > 0);
            }else{
                n_hits++;
            }
        }else{
            n_hits++;
        }

        // Izquierda:
        if( (x-3) > 0 ){
            if( ((y-3+m) > 0) && ((y-3+m) < mis_rens) ){
                n_hits += (ptr[(x-3) + (y-3+m)*mis_cols] > 0);
            }else{
                n_hits++;
            }
        }else{
            n_hits++;
        }

        // Derecha:
        if( (x+3) < mis_cols ){
            if( ((y-3+m) > 0) && ((y-3+m) < mis_rens) ){
                n_hits += (ptr[(x+3) + (y-3+m)*mis_cols] > 0);
            }else{
                n_hits++;
            }
        }else{
            n_hits++;
        }
    }

    return (n_hits == 28);
}



/*  Metodo: regionFilling5
    Funcion: Rellena espacios vacios dentro del conjunto de pixeles resaltado.
*/
bool IMGVTK::regionFilling5( const double *ptr, const int x, const int y, const int mis_cols, const int mis_rens){
    int n_hits = 0;

    for( int m = 0; m < 5; m++){
        // Arriba:
        if( (y-2) > 0 ){
            if( ((x-2+m) > 0) && ((x-2+m) < mis_cols) ){
                n_hits += (ptr[(x-2+m) + (y-2)*mis_cols] > 0);
            }else{
                n_hits++;
            }
        }else{
            n_hits++;
        }

        // Abajo:
        if( (y+2) < mis_rens ){
            if( ((x-2+m) > 0) && ((x-2+m) < mis_cols) ){
                n_hits += (ptr[(x-2+m) + (y+2)*mis_cols] > 0);
            }else{
                n_hits++;
            }
        }else{
            n_hits++;
        }

        // Izquierda:
        if( (x-2) > 0 ){
            if( ((y-2+m) > 0) && ((y-2+m) < mis_rens) ){
                n_hits += (ptr[(x-2) + (y-2+m)*mis_cols] > 0);
            }else{
                n_hits++;
            }
        }else{
            n_hits++;
        }

        // Derecha:
        if( (x+2) < mis_cols ){
            if( ((y-2+m) > 0) && ((y-2+m) < mis_rens) ){
                n_hits += (ptr[(x+2) + (y-2+m)*mis_cols] > 0);
            }else{
                n_hits++;
            }
        }else{
            n_hits++;
        }
    }

    return (n_hits == 20);
}



/*  Metodo: regionFilling3
    Funcion: Rellena espacios vacios dentro del conjunto de pixeles resaltado.
*/
bool IMGVTK::regionFilling3( const double *ptr, const int x, const int y, const int mis_cols, const int mis_rens ){
    int n_hits = 0;

    for( int m = 0; m < 3; m++){
        // Arriba:
        if( (y-1) > 0 ){
            if( ((x-1+m) > 0) && ((x-1+m) < mis_cols) ){
                n_hits += (ptr[(x-1+m) + (y-1)*mis_cols] > 0);
            }else{
                n_hits++;
            }
        }else{
            n_hits++;
        }

        // Abajo:
        if( (y+1) < mis_rens ){
            if( ((x-1+m) > 0) && ((x-1+m) < mis_cols) ){
                n_hits += (ptr[(x-1+m) + (y+1)*mis_cols] > 0);
            }else{
                n_hits++;
            }
        }else{
            n_hits++;
        }

        // Izquierda:
        if( (x-1) > 0 ){
            if( ((y-1+m) > 0) && ((y-1+m) < mis_rens) ){
                n_hits += (ptr[(x-1) + (y-1+m)*mis_cols] > 0);
            }else{
                n_hits++;
            }
        }else{
            n_hits++;
        }

        // Derecha:
        if( (x+1) < mis_cols ){
            if( ((y-1+m) > 0) && ((y-1+m) < mis_rens) ){
                n_hits += (ptr[(x+1) + (y-1+m)*mis_cols] > 0);
            }else{
                n_hits++;
            }
        }else{
            n_hits++;
        }
    }

    return (n_hits == 12);
}



/*  Metodo: conexo
    Funcion: Metodo recursivo (dinamico) para encontrar los conjuntos conexos utilizando la conectividad 8.
*/
void IMGVTK::conexo(const double *ptr, const int x, const int y, int *conjuntos, unsigned int* n_etiquetados, bool* visitados, const int num_etiquetas, const int mis_cols, const int mis_rens){
    visitados[x + y*mis_cols] = true;

    // Si el pixel es parte del fondo, no se hace nada:
    if( ptr[x + y*mis_cols] < 1.0){
        return;
    }else{
        conjuntos[x + y*mis_cols] = num_etiquetas;
        n_etiquetados[ num_etiquetas ] ++;
    }

    // De lo contrario, se revisan los pixeles vecinos:
    /// Superiores:
    if( y > 0 ){
        if( x > 0 ){
            /// Ariba - Izquierda
            if( ! visitados[x-1 + (y-1)*mis_cols] ){ // Si el pixel arriba-izquierda no ha sido visitado:
                conexo( ptr, x-1, y-1, conjuntos, n_etiquetados, visitados, num_etiquetas, mis_cols, mis_rens );
            }
        }

        /// Arriba
        if( ! visitados[x + (y-1)*mis_cols] ){ // Si el pixel arriba no ha sido visitado:
            conexo( ptr, x, y-1, conjuntos, n_etiquetados, visitados, num_etiquetas, mis_cols, mis_rens );
        }

        /// Ariba - Derecha
        if( x < (mis_cols-1) ){
            if( ! visitados[x+1 + (y-1)*mis_cols] ){ // Si el pixel arriba-derecha no ha sido visitado:
                conexo( ptr, x+1, y-1, conjuntos, n_etiquetados, visitados, num_etiquetas, mis_cols, mis_rens );
            }
        }
    }

    /// Centrales:
    if( x > 0 ){
        /// Centro - Izquierda
        if( ! visitados[x-1 + y*mis_cols] ){ // Si el pixel centro-izquierda no ha sido visitado:
            conexo( ptr, x-1, y, conjuntos, n_etiquetados, visitados, num_etiquetas, mis_cols, mis_rens );
        }
    }

    /// Centro - Derecha
    if( x < (mis_cols-1) ){
        if( ! visitados[x+1 + y*mis_cols] ){ // Si el pixel arriba-izquierda no ha sido visitado:
            conexo( ptr, x+1, y, conjuntos, n_etiquetados, visitados, num_etiquetas, mis_cols, mis_rens );
        }
    }

    /// Inferiores:
    if( y < (mis_rens-1) ){
        if( x > 0 ){
            /// Abajo - Izquierda
            if( ! visitados[x-1 + (y+1)*mis_cols] ){ // Si el pixel arriba-izquierda no ha sido visitado:
                conexo( ptr, x-1, y+1, conjuntos, n_etiquetados, visitados, num_etiquetas, mis_cols, mis_rens );
            }
        }

        /// Abajo
        if( ! visitados[x + (y+1)*mis_cols] ){ // Si el pixel arriba no ha sido visitado:
            conexo( ptr, x, y+1, conjuntos, n_etiquetados, visitados, num_etiquetas, mis_cols, mis_rens );
        }

        /// Abajo - Derecha
        if( (x < (mis_cols-1)) ){
            if( ! visitados[x+1 + (y+1)*mis_cols] ){ // Si el pixel arriba-derecha no ha sido visitado:
                conexo( ptr, x+1, y+1, conjuntos, n_etiquetados, visitados, num_etiquetas, mis_cols, mis_rens );
            }
        }
    }
}



/*  Metodo: conjuntosConexos
    Funcion: Obtiene los conjuntos conexos de la imagen usando programacion dinamica (Llega a generar stackoverflow ... ).
*/
unsigned int* IMGVTK::conjuntosConexosDinamico(const double *ptr, int *conjuntos, const int mis_cols, const int mis_rens){
    int num_etiquetas = 0;

    const int mis_rens_cols = mis_cols*mis_rens;

    bool *visitados = new bool [mis_rens_cols];
    memset(visitados, 0, mis_rens_cols*sizeof(bool));

    unsigned int *tmp = new unsigned int [mis_rens_cols];
    memset(tmp, 0, sizeof(unsigned int)*mis_rens_cols);


    memset(conjuntos, -1, sizeof(int) * mis_rens_cols);

    for( int y = 0; y< mis_rens; y++){
        for( int x = 0; x < mis_cols; x++){
            /// Se mueve la etiqueta si el pixel es de un nuevo conjunto:
            if( !visitados[x + y*mis_cols]) {
                conexo( ptr, x, y, conjuntos, tmp, visitados, num_etiquetas, mis_cols, mis_rens );
                if( ptr[x + y*mis_cols] > 0 ){
                    num_etiquetas++;
                }
            }
        }
    }

    unsigned int *n_etiquetados = new unsigned int [num_etiquetas+1];
    memcpy(n_etiquetados + 1, tmp, num_etiquetas*sizeof(unsigned int) );
    n_etiquetados[0] = num_etiquetas;

    delete [] tmp;
    delete [] visitados;

    DEB_MSG( "Numero de etiquetas puestas: " << num_etiquetas);

    return n_etiquetados;
}




/*  Metodo: ampliarConjunto
    Funcion: Amplia la memoria requerida para las etiqeutas de conjuntos equivalentes.
*/
inline void IMGVTK::ampliarConjunto( int **etiquetas, const int equiv_A, const int equiv_B){
    int base, nueva_eq;

    if( equiv_A < equiv_B ){
        base = equiv_A;
        nueva_eq = equiv_B;
    }else{
        base = equiv_B;
        nueva_eq = equiv_A;
    }


    int n_elementos = *(*(etiquetas + base));
    int n_espacios = *(*(etiquetas + base) + 1);

    if( (n_elementos + 1) > n_espacios ){

        int *swap = *(etiquetas + base);
        *(etiquetas + base) = new int [ 2 + n_espacios + 50 ];
        memcpy( *(etiquetas + base) + 2, swap + 2, n_espacios*sizeof(int));
        delete [] swap;

        *(*(etiquetas + base)) = n_elementos;
        *(*(etiquetas + base) + 1) = n_espacios + 50;
    }

    bool encontrado = false;

    for( int i = 1; i < n_elementos; i++){
        encontrado = ( *(*(etiquetas + base) + 2 + i) == nueva_eq );
        if( encontrado ){
            break;
        }
    }

    if( !encontrado ){
        *(*(etiquetas + base)) = n_elementos + 1;
        *(*(etiquetas + base) + 2 + n_elementos) = nueva_eq;
    }
}




/*  Metodo: conjuntosConexos
    Funcion: Obtiene los conjuntos conexos de la imagen.
*/
unsigned int* IMGVTK::conjuntosConexos(const double *ptr, int *conjuntos, const int mis_cols, const int mis_rens){

    const int mis_rens_cols = mis_cols*mis_rens;

    int max_etiquetas = mis_rens;
    int **val_etiquetas = new int* [max_etiquetas];

    for( int i = 0; i < max_etiquetas; i++){
        *(val_etiquetas + i) = new int [2 + mis_cols];
        *(*(val_etiquetas + i)) = 1;
        *(*(val_etiquetas + i) + 1) = mis_cols;
        *(*(val_etiquetas + i) + 2) = -1;
    }

    int *pix_etq = new int [(mis_rens+2)*(mis_cols+2)];
    memset( pix_etq, -1, (mis_rens+2)*(mis_cols+2) * sizeof(int));

    int n_conjunto = -1;

    for( int y = 1; y <= mis_rens; y++){
        for( int x = 1; x <= mis_cols; x++){
            if( *(ptr + (y-1)*mis_cols + (x-1)) > 0.0 ){

                const int NO = *(pix_etq + (y-1)*(mis_cols+2) + (x-1));
                const int N  = *(pix_etq + (y-1)*(mis_cols+2) + x);
                const int NE = *(pix_etq + (y-1)*(mis_cols+2) + (x+1));
                const int E  = *(pix_etq + y*(mis_cols+2) + (x+1));
                const int SE = *(pix_etq + (y+1)*(mis_cols+2) + (x+1));
                const int S  = *(pix_etq + (y+1)*(mis_cols+2) + x);
                const int SO = *(pix_etq + (y+1)*(mis_cols+2) + (x-1));
                const int O  = *(pix_etq + y*(mis_cols+2) + (x-1));


                int base = -1;

                /// NO
                if( NO >= 0 ){
                    base = NO;
                }

                /// N
                if( N >= 0 ){
                    if( (base >= 0) && (base != N) ){
                        ampliarConjunto( val_etiquetas, base, N);
                    }else{
                        base = N;
                    }
                }

                /// NE
                if( NE >= 0 ){
                    if( (base >= 0) && (base != NE) ){
                        ampliarConjunto( val_etiquetas, base, NE);
                    }else{
                        base = NE;
                    }
                }

                /// E
                if( E >= 0){
                    if( (base >= 0) && (base != E) ){
                        ampliarConjunto( val_etiquetas, base, E);
                    }else{
                        base = E;
                    }
                }

                /// SE
                if( SE >= 0 ){
                    if( (base >= 0) && (base != SE) ){
                        ampliarConjunto( val_etiquetas, base, SE);
                    }else{
                        base = SE;
                    }
                }

                /// S
                if( S >= 0 ){
                    if( (base >= 0) && (base != S) ){
                        ampliarConjunto( val_etiquetas, base, S);
                    }else{
                        base = S;
                    }
                }

                /// SO
                if( SO >= 0 ){
                    if( (base >= 0) && (base != SO) ){
                        ampliarConjunto( val_etiquetas, base, SO);
                    }else{
                        base = SO;
                    }
                }

                /// O
                if( O >= 0 ){
                    if( (base >= 0) && (base != O) ){
                        ampliarConjunto( val_etiquetas, base, O);
                    }else{
                        base = O;
                    }
                }
                if( base >= 0 ){
                    *(pix_etq + y*(mis_cols+2) + x) = base;
                }else{
                    if( (++n_conjunto) > max_etiquetas ){
                        int **swap = val_etiquetas;
                        val_etiquetas = new int* [max_etiquetas + mis_cols];
                        memcpy( val_etiquetas, swap, max_etiquetas*sizeof(int*) );
                        delete [] swap;

                        for( int i = max_etiquetas; i < (max_etiquetas + mis_cols); i++){
                            *(val_etiquetas + i) = new int [mis_cols + 2];
                            *(*(val_etiquetas + i)) = 1;
                            *(*(val_etiquetas + i) + 1) = mis_cols;
                            *(*(val_etiquetas + i) + 2) = -1;
                        }

                        max_etiquetas += mis_cols;
                    }
                    *(*(val_etiquetas + n_conjunto) + 2) = n_conjunto;
                    *(pix_etq + y*(mis_cols+2) + x) = n_conjunto;
                }
            }
        }
    }

    DEB_MSG("Quedaron " << n_conjunto << " activados");

    /// Determinar cuantos conjuntos diferentes hay y encontrar las etiquetas equivalentes:
    int etq_activas = n_conjunto+1;
    n_conjunto = 0;

    int *diferentes = new int [etq_activas];
    memset(diferentes, -1, etq_activas*sizeof(int) );

    for( int i = 0; i < etq_activas; i++){
        bool encontrado = false;
        const int base = *(*(val_etiquetas + i) + 2);

        for( int j = 0; j < n_conjunto; j++){
            if( base == *(diferentes + j) ){
                encontrado = true;
                break;
            }
        }

        if( !encontrado ){
            *(diferentes + n_conjunto) = base;
            n_conjunto++;
        }

        for( int j = 0; j < *(*(val_etiquetas + i)); j++){
            const int equivalente = *(*(val_etiquetas + i) + 2 + j);
            *(*(val_etiquetas + equivalente) + 2) = base;
        }
    }

    for( int i = 0; i < etq_activas; i++){
        const int base = *(*(val_etiquetas + i) + 2);

        int nueva_base;
        for( int j = 0; j < n_conjunto; j++){
            if( base == *(diferentes + j) ){
                nueva_base = j;
                break;
            }
        }
        *(*(val_etiquetas + i) + 2) = nueva_base;
    }

    memset( conjuntos, -1, mis_rens_cols*sizeof(int) );

    // Indicar a que conjunto pertenece cada pixel:
    for( int y = 0; y < mis_rens; y++){
        for( int x = 0; x < mis_cols; x++){
            const int etiqueta = *(pix_etq + (y+1)*(mis_cols+2) + (x+1));
            if( etiqueta >= 0 ){
                *(conjuntos + x + y*mis_cols) = *(*(val_etiquetas + etiqueta) + 2);
            }
        }
    }

    /// Contar cuantos elementos hay en cada grupo:
    unsigned int *n_etiquetados = new unsigned int [n_conjunto+1];
    memset(n_etiquetados + 1, 0, n_conjunto*sizeof(unsigned int));
    n_etiquetados[0] = n_conjunto;


    for( int xy = 0; xy < mis_rens_cols; xy++){
        if( *(conjuntos + xy) >= 0 ){
            *(n_etiquetados + 1 + *(conjuntos + xy)) = *(n_etiquetados + 1 + *(conjuntos + xy)) + 1;
        }
    }

    for( int i = 0; i < max_etiquetas; i++){
        delete [] *(val_etiquetas + i);
    }

    delete [] diferentes;
    delete [] val_etiquetas;
    delete [] pix_etq;

    return n_etiquetados;

}




/*  Metodo:  sklMask
    Funcion: Mascara usada para la extraccion del esqueleto.
*/
inline unsigned char IMGVTK::sklMask( const double *skl_ptr, const int x, const int y, const int mis_cols, const int mis_rens ){
    return   1*(skl_ptr[(x-1) + (y-1)*mis_cols] > 0) + /* P2 */
             2*(skl_ptr[  x   + (y-1)*mis_cols] > 0) + /* P3 */
             4*(skl_ptr[(x+1) + (y-1)*mis_cols] > 0) + /* P4 */
             8*(skl_ptr[(x+1) +   y  *mis_cols] > 0) + /* P5 */
            16*(skl_ptr[(x+1) + (y+1)*mis_cols] > 0) + /* P6 */
            32*(skl_ptr[  x   + (y+1)*mis_cols] > 0) + /* P7*/
            64*(skl_ptr[(x-1) + (y+1)*mis_cols] > 0) + /* P8 */
           128*(skl_ptr[(x-1) +   y  *mis_cols] > 0);  /* P9 */
}




/*  Metodo: lengthFiltering
    Funcion: Filtra los conjuntos con numero de pixeles menor a 'min_length'
*/
void IMGVTK::lengthFilter(double *ptr, const int min_length, const int mis_cols, const int mis_rens ){

    const int mis_rens_cols = mis_cols*mis_rens;

    int *mis_conjuntos = new int [mis_rens_cols];

    unsigned int *mis_n_etiquetados = conjuntosConexos(ptr, mis_conjuntos, mis_cols, mis_rens);
DEB_MSG("n etiquetas: " << mis_n_etiquetados[0] );

    for( int xy = 0; xy < mis_rens_cols; xy++){
        if( (mis_conjuntos[ xy ] >= 0) && (mis_n_etiquetados[ mis_conjuntos[ xy ]+1 ] < min_length)){
            ptr[ xy ] = 0.0;
        }
    }

    delete [] mis_conjuntos;
    delete [] mis_n_etiquetados;
}



/*  Metodo: lengthFiltering (Publica)
    Funcion: Filtra los conjuntos con numero de pixeles menor a 'min_length'
*/
void IMGVTK::lengthFilter(IMG_IDX img_idx, const int min_length){

    double *img_ptr = NULL;
    switch( img_idx ){
        case BASE:
            img_ptr = base_ptr;
            break;
        case MASK:
            img_ptr = mask_ptr;
            break;
        case SEGMENT:
            img_ptr = segment_ptr;
            break;
        case THRESHOLD:
            img_ptr = threshold_ptr;
            break;
    }

    lengthFilter(img_ptr, min_length, cols, rows);
}


/*  Metodo: regionFill
    Funcion: Rellena vacios dentro del cuerpo de la arteria segmentada.
*/
void IMGVTK::regionFill( double *ptr, const int mis_cols, const int mis_rens  ){
    // -------- Mascara 9x9
    for( int y = 0; y < mis_rens; y++){
        for( int x = 0; x < mis_cols; x++){
            if( ptr[x+y*mis_cols] < 1.0 ){
                ptr[x+y*mis_cols] = (regionFilling9(ptr, x, y, mis_cols, mis_rens) ? 1.0 : 0.0);
            }
        }
    }
    // -------- Mascara 7x7
    for( int y = 0; y < mis_rens; y++){
        for( int x = 0; x < mis_cols; x++){
            if( ptr[x+y*mis_cols] < 1.0 ){
                ptr[x+y*mis_cols] = (regionFilling7(ptr, x, y, mis_cols, mis_rens) ? 1.0 : 0.0);
            }
        }
    }
    // -------- Mascara 5x5
    for( int y = 0; y < mis_rens; y++){
        for( int x = 0; x < mis_cols; x++){
            if( ptr[x+y*mis_cols] < 1.0 ){
                ptr[x+y*mis_cols] = (regionFilling5(ptr, x, y, mis_cols, mis_rens) ? 1.0 : 0.0);
            }
        }
    }
    // -------- Mascara 3x3
    for( int y = 0; y < mis_rens; y++){
        for( int x = 0; x < mis_cols; x++){
            if( ptr[x+y*mis_cols] < 1.0 ){
                ptr[x+y*mis_cols] = (regionFilling3(ptr, x, y, mis_cols, mis_rens) ? 1.0 : 0.0);
            }
        }
    }
}




/*  Metodo: regionFill (Publica)
    Funcion: Rellena vacios dentro del cuerpo de la arteria segmentada.
*/
void IMGVTK::regionFill(IMG_IDX img_idx){

    double *img_ptr = NULL;
    switch( img_idx ){
        case BASE:
            img_ptr = base_ptr;
            break;
        case MASK:
            img_ptr = mask_ptr;
            break;
        case SEGMENT:
            img_ptr = segment_ptr;
            break;
        case THRESHOLD:
            img_ptr = threshold_ptr;
            break;
    }

    regionFill(img_ptr, cols, rows);
}



/*  Metodo:  dilMask
    Funcion: Mascara usada para dilatar unaimagen.
*/
inline unsigned char IMGVTK::dilMask( const double *mask_dil, const int x, const int y, const int mis_cols, const int mis_rens){
    return (mask_dil[( x ) + (y-1)*mis_cols] > 0) +
           (mask_dil[(x+1) + ( y )*mis_cols] > 0) +
           (mask_dil[( x ) + (y+1)*mis_cols] > 0) +
           (mask_dil[(x-1) + ( y )*mis_cols] > 0);
}



/*  Metodo:  erosionMask
    Funcion: Mascara para erosion usando un disco de radio 5
*/
inline unsigned char IMGVTK::erosionMask( const double *ptr_tmp, const int x, const int y, const int mis_cols, const int mis_rens ){
    return                                                                                   (ptr_tmp[(x-2) + (y-4)*mis_cols] > 0) + (ptr_tmp[(x-1) + (y-4)*mis_cols] > 0) + (ptr_tmp[( x ) + (y-4)*mis_cols] > 0) + (ptr_tmp[(x+1) + (y-4)*mis_cols] > 0) + (ptr_tmp[(x+2) + (y-1)*mis_cols] > 0) +
                                                     (ptr_tmp[(x-3) + (y-3)*mis_cols] > 0) + (ptr_tmp[(x-2) + (y-3)*mis_cols] > 0) + (ptr_tmp[(x-1) + (y-3)*mis_cols] > 0) + (ptr_tmp[( x ) + (y-3)*mis_cols] > 0) + (ptr_tmp[(x+1) + (y-3)*mis_cols] > 0) + (ptr_tmp[(x+2) + (y-2)*mis_cols] > 0) + (ptr_tmp[(x+3) + (y-3)*mis_cols] > 0) +
             (ptr_tmp[(x-4) + (y-2)*mis_cols] > 0) + (ptr_tmp[(x-3) + (y-2)*mis_cols] > 0) + (ptr_tmp[(x-2) + (y-2)*mis_cols] > 0) + (ptr_tmp[(x-1) + (y-2)*mis_cols] > 0) + (ptr_tmp[( x ) + (y-2)*mis_cols] > 0) + (ptr_tmp[(x+1) + (y-2)*mis_cols] > 0) + (ptr_tmp[(x+2) + (y-3)*mis_cols] > 0) + (ptr_tmp[(x+3) + (y-2)*mis_cols] > 0) + (ptr_tmp[(x+4) + (y-2)*mis_cols] > 0) +
             (ptr_tmp[(x-4) + (y-1)*mis_cols] > 0) + (ptr_tmp[(x-3) + (y-1)*mis_cols] > 0) + (ptr_tmp[(x-2) + (y-1)*mis_cols] > 0) + (ptr_tmp[(x-1) + (y-1)*mis_cols] > 0) + (ptr_tmp[( x ) + (y-1)*mis_cols] > 0) + (ptr_tmp[(x+1) + (y-1)*mis_cols] > 0) + (ptr_tmp[(x+2) + (y-4)*mis_cols] > 0) + (ptr_tmp[(x+3) + (y-1)*mis_cols] > 0) + (ptr_tmp[(x+4) + (y-1)*mis_cols] > 0) +
             (ptr_tmp[(x-4) + ( y )*mis_cols] > 0) + (ptr_tmp[(x-3) + ( y )*mis_cols] > 0) + (ptr_tmp[(x-2) + ( y )*mis_cols] > 0) + (ptr_tmp[(x-1) + ( y )*mis_cols] > 0)                  +                        (ptr_tmp[(x+1) + ( y )*mis_cols] > 0) + (ptr_tmp[(x+2) + ( y )*mis_cols] > 0) + (ptr_tmp[(x+3) + ( y )*mis_cols] > 0) + (ptr_tmp[(x+4) + ( y )*mis_cols] > 0) +
             (ptr_tmp[(x-4) + (y+1)*mis_cols] > 0) + (ptr_tmp[(x-3) + (y+1)*mis_cols] > 0) + (ptr_tmp[(x-2) + (y+1)*mis_cols] > 0) + (ptr_tmp[(x-1) + (y+1)*mis_cols] > 0) + (ptr_tmp[( x ) + (y+1)*mis_cols] > 0) + (ptr_tmp[(x+1) + (y+1)*mis_cols] > 0) + (ptr_tmp[(x+2) + (y+1)*mis_cols] > 0) + (ptr_tmp[(x+3) + (y+1)*mis_cols] > 0) + (ptr_tmp[(x+4) + (y+1)*mis_cols] > 0) +
             (ptr_tmp[(x-4) + (y+2)*mis_cols] > 0) + (ptr_tmp[(x-3) + (y+2)*mis_cols] > 0) + (ptr_tmp[(x-2) + (y+2)*mis_cols] > 0) + (ptr_tmp[(x-1) + (y+2)*mis_cols] > 0) + (ptr_tmp[( x ) + (y+2)*mis_cols] > 0) + (ptr_tmp[(x+1) + (y+2)*mis_cols] > 0) + (ptr_tmp[(x+2) + (y+2)*mis_cols] > 0) + (ptr_tmp[(x+3) + (y+2)*mis_cols] > 0) + (ptr_tmp[(x+4) + (y+2)*mis_cols] > 0) +
                                                     (ptr_tmp[(x-3) + (y+3)*mis_cols] > 0) + (ptr_tmp[(x-2) + (y+3)*mis_cols] > 0) + (ptr_tmp[(x-1) + (y+3)*mis_cols] > 0) + (ptr_tmp[( x ) + (y+3)*mis_cols] > 0) + (ptr_tmp[(x+1) + (y+3)*mis_cols] > 0) + (ptr_tmp[(x+2) + (y+3)*mis_cols] > 0) + (ptr_tmp[(x+3) + (y+3)*mis_cols] > 0) +
                                                                                             (ptr_tmp[(x-2) + (y+4)*mis_cols] > 0) + (ptr_tmp[(x-1) + (y+4)*mis_cols] > 0) + (ptr_tmp[( x ) + (y+4)*mis_cols] > 0) + (ptr_tmp[(x+1) + (y+4)*mis_cols] > 0) + (ptr_tmp[(x+2) + (y+4)*mis_cols] > 0);
}



/*  Metodo: erosionar
    Funcion: Rellena vacios dentro del cuerpo de la arteria segmentada.
*/
void IMGVTK::erosionar( double *ptr , const int mis_cols, const int mis_rens){
    double *ptr_tmp = new double [(mis_rens+8)*(mis_cols+8)];

    for(int xy = 0; xy < (mis_rens+8)*(mis_cols+8); xy++){
        ptr_tmp[xy] = 1.0;
    }

    for(int y = 0; y < mis_rens; y++){
        memcpy( ptr_tmp + (y+4)*(mis_cols+8) + 4, ptr + y*mis_cols, mis_cols*sizeof(double));
    }

    for( int y = 0; y < mis_rens; y++){
        for( int x = 0; x < mis_cols; x++){
            if(ptr_tmp[(x+4)+(y+4)*(mis_cols+8)] > 0){
                const unsigned char resp = erosionMask(ptr_tmp, x+4, y+4, mis_cols+8, mis_rens);
                if( (0 < resp) && (resp < 68) ){
                    ptr[x+y*mis_cols] = 0.0;
                }
            }
        }
    }

    delete [] ptr_tmp;
}



/*  Metodo: maskFOV
    Funcion: Obtiene la mascara del contorno
*/
void IMGVTK::maskFOV( double * img_tmp, double *mask_tmp, const int mis_cols, const int mis_rens){

    const int mis_rens_cols = mis_cols * mis_rens;

    // Se umbraliza al 0.1 la imagen original
    for( int xy = 0; xy < mis_rens_cols; xy++){
        mask_tmp[xy] = (img_tmp[xy] < 0.1) ? 0.0 : 1.0;
    }

    // Se eliminan los conjuntos pequeños
    lengthFilter(mask_tmp, 1000, mis_cols, mis_rens);

    // Se erosiona la mascara:
    erosionar(mask_tmp, mis_cols, mis_rens);

    // Se eliminan los conjuntos grandes que no esten en las esquinas:
    // Se extraen las etiquetas de los conjuntos que se conectan a las esquinas:
    for( int xy = 0; xy < mis_rens_cols; xy++){
        *(mask_tmp + xy) = 1.0 - *(mask_tmp + xy);
    }

    int *mis_conjuntos = new int [mis_rens_cols];

    // Se buscan los conjuntos que no esten en las esquinas para eliminarlos
    unsigned int *mis_n_etiquetados = conjuntosConexos(mask_tmp, mis_conjuntos, mis_cols, mis_rens);
    delete [] mis_n_etiquetados;


    //// etiqueta de los conjuntos donde existe una esquina
    const int NO = mis_conjuntos[0], NE = mis_conjuntos[mis_cols-1], SO = mis_conjuntos[(mis_rens-1)*mis_cols], SE = mis_conjuntos[mis_rens_cols-1];
    DEB_MSG("NO: " << NO);
    DEB_MSG("NE: " << NE);
    DEB_MSG("SO: " << SO);
    DEB_MSG("SE: " << SE);

    for( int xy = 0; xy < mis_rens_cols; xy++){
        if( (mis_conjuntos[xy] >= 0) && ((mis_conjuntos[xy] == NO) || (mis_conjuntos[xy] == NE) || (mis_conjuntos[xy] == SO) || (mis_conjuntos[xy] == SE)) ){
            mask_tmp[xy] = 0.0;
        }else{
            mask_tmp[xy] = 1.0;
        }
    }
    delete [] mis_conjuntos;
}



/*  Metodo: fillMask
    Funcion: Se llenan los espacios de la mascara con la media de la imagen circundante.
*/
void IMGVTK::fillMask( double *img_tmp, double *mask_tmp, const int mis_cols, const int mis_rens ){

    PIX_PAR par_tmp;
    par_tmp.pix_tipo = PIX_CROSS;

    double *mask_dil = new double [(mis_rens+2)*(mis_cols+2)];
    for( int xy = 0; xy < (mis_rens+2)*(mis_cols+2); xy++){
        mask_dil[xy] = 0.0  ;
    }

    for( int y = 0; y < mis_rens; y++){
        memcpy( mask_dil + (y+1)*(mis_cols+2) + 1, mask_tmp + y*mis_cols, mis_cols*sizeof(double));
    }

    int iter = 0;
    while( iter < mis_cols ){
        // Obtener el borde de la mascara realizando una dilatacion:
        std::vector< PIX_PAR > borde;

        // Si es la primer pasada, se verifican todos los pixeles para encontrar el borde
        for( int y = 0; y < mis_rens; y++){
            for( int x = 0; x < mis_cols; x++){
                if( mask_dil[(x+1) + (y+1)*(mis_cols+2)] < 1.0 ){
                    // Si es HIT se considera como borde, si es FIT no:
                    const unsigned char resp = dilMask(mask_dil, x+1, y+1, mis_cols+2, mis_rens);
                    if( resp > 0){
                        // Se guardan las coordenadas de los pixeles que forman el borde:
                        par_tmp.x = x;
                        par_tmp.y = y;
                        borde.push_back( par_tmp );
                    }
                }
            }
        }

        const int n_borde = borde.size();
        DEB_MSG("Iter: " << iter << ", pixeles en el borde: " <<  n_borde);

        //// Si ya no existen pixeles en el borde, se termina el ciclo:
        if( n_borde == 0){
            break;
        }

        // Para cada pixel en el borde, se calcula la media usando una ventana de 21 x 21 pixeles.
        for( int b = 0; b < n_borde; b++ ){
            const int x_act = borde[b].x;
            const int y_act = borde[b].y;
            const int offset_x_izq = (x_act < 10) ? 0 : (x_act - 10);
            const int offset_x_der = (x_act >= (mis_cols - 10)) ? (mis_cols-1) : (x_act + 10);
            const int offset_y_sup = (y_act < 10) ? 0 : (y_act - 10);
            const int offset_y_inf = (y_act >= (mis_rens - 10)) ? (mis_rens-1) : (y_act + 10);

            double suma = 0.0;
            int n_vecinos = 0;

            for( int y = offset_y_sup; y <= offset_y_inf; y++){
                for( int x = offset_x_izq; x <= offset_x_der; x++){
                    if( mask_dil[(x+1) + (y+1)*(mis_cols+2)] > 0 ){
                        suma += (double)img_tmp[ x + y*mis_cols ];
                        n_vecinos++;
                    }
                }
            }

            img_tmp[x_act + y_act*mis_cols] = suma / n_vecinos;
            mask_dil[(x_act+1) + (y_act+1)*(mis_cols+2)] = 1.0;
        }
        iter++;
    }
    delete [] mask_dil;
}





/*  Metodo: grafoSkeleton
    Funcion: Genera un grafo a partir del esqueleto.
*/
IMGVTK::PIX_PAR* IMGVTK::grafoSkeleton(double *skl_tmp, const int x, const int y, int *nivel, const unsigned char *lutabla, bool *visitados){
    if( *(visitados + x + y*cols) ){
        return NULL;
    }

    PIX_PAR *temp = new PIX_PAR;

    const unsigned char resp = sklMask( skl_tmp, x, y, cols+2, rows);

    temp->x = (x-1 - (double)cols/2)*pixX;
    temp->y = (y-1 - (double)rows/2)*pixY;

    temp->n_hijos = 0;

    temp->hijos[0] = NULL;
    temp->hijos[1] = NULL;
    temp->hijos[2] = NULL;

    temp->nivel = *nivel;

    /// Calcular el radio de la arteria en el pixel actual:
    const int min_x = ((x - max_dist - 1) < 0) ? 0 : (x - max_dist - 1);
    const int max_x = ((x + max_dist - 1) >= cols) ? cols : (x + max_dist - 1);
    const int min_y = ((y - max_dist - 1) < 0) ? 0 : (y - max_dist - 1);
    const int max_y = ((y + max_dist - 1) >= rows) ? rows : (y + max_dist - 1);

    double dist, radio = (max_dist + 0.5) * (max_dist + 0.5);
    double x_r, y_r;
    for( int yy = min_y; yy < max_y; yy++){
        for( int xx = min_x; xx < max_x; xx++){
            dist = (yy - y + 1.5)*(yy - y + 1.5) + (xx - x + 1.5)*(xx - x + 1.5);
            if( (borders_ptr[xx + yy*cols] > 0.0) && (dist < radio) ){
                radio = dist;
                x_r = (double)xx;
                y_r = (double)yy;
            }
        }
    }

    temp->radio = sqrt(radio) * pixX;
    temp->y_r = (y_r - (double)rows/2)*pixY;
    temp->x_r = (x_r - (double)cols/2)*pixX;
    temp->alpha = atan2(temp->y_r - temp->y, temp->x_r - temp->x);// + PI / 2.0;

    switch( lutabla[ resp ] ){
        case (unsigned char)1:{ /* END point*/
            temp->pix_tipo = PIX_END;
            break;
        }
        case (unsigned char)2:{ /* BRANCH point */
            temp->pix_tipo = PIX_BRANCH;
            break;
        }
        case (unsigned char)3:{ /* CROSS point */
            temp->pix_tipo = PIX_CROSS;
            break;
        }
        default:{
            temp->pix_tipo = PIX_SKL;
            break;
        }
    }

    *(visitados + x + y*cols) = true;

    /// NO
    if( (resp & (unsigned char)1) && (*(skl_tmp + (x-1) + (y-1)*(cols+2)) < 2.0)){
        if( temp->pix_tipo == PIX_CROSS || temp->pix_tipo == PIX_BRANCH ){
            *nivel = *nivel + 1;
        }
        temp->hijos[temp->n_hijos] = grafoSkeleton(skl_tmp, x-1, y-1, nivel, lutabla, visitados);

        if( temp->hijos[temp->n_hijos] ){
            temp->n_hijos++;
        }
    }

    /// N
    if( (resp & (unsigned char)2) && (*(skl_tmp + x + (y-1)*(cols+2)) < 2.0) ){
        if( temp->pix_tipo == PIX_CROSS || temp->pix_tipo == PIX_BRANCH ){
            *nivel = *nivel + 1;
        }
        temp->hijos[temp->n_hijos] = grafoSkeleton(skl_tmp, x, y-1, nivel, lutabla, visitados);
        if( temp->hijos[temp->n_hijos] ){
            temp->n_hijos++;
        }
    }

    /// NE
    if( (resp & (unsigned char)4) && (*(skl_tmp + (x+1) + (y-1)*(cols+2)) < 2.0) ){
        if( temp->pix_tipo == PIX_CROSS || temp->pix_tipo == PIX_BRANCH ){
            *nivel = *nivel + 1;
        }
        temp->hijos[temp->n_hijos] = grafoSkeleton(skl_tmp, x+1, y-1, nivel, lutabla, visitados);

        if( temp->hijos[temp->n_hijos] ){
            temp->n_hijos++;
        }
    }

    /// E
    if( (resp & (unsigned char)8) && (*(skl_tmp + (x+1) + y*(cols+2)) < 2.0) ){
        if( temp->pix_tipo == PIX_CROSS || temp->pix_tipo == PIX_BRANCH ){
            *nivel = *nivel + 1;
        }
        temp->hijos[temp->n_hijos] = grafoSkeleton(skl_tmp, x+1, y, nivel, lutabla, visitados);

        if( temp->hijos[temp->n_hijos] ){
            temp->n_hijos++;
        }
    }

    /// SE
    if( (resp & (unsigned char)16) && (*(skl_tmp + (x+1) + (y+1)*(cols+2)) < 2.0) ){
        if( temp->pix_tipo == PIX_CROSS || temp->pix_tipo == PIX_BRANCH ){
            *nivel = *nivel + 1;
        }
        temp->hijos[temp->n_hijos] = grafoSkeleton(skl_tmp, x+1, y+1, nivel, lutabla, visitados);

        if( temp->hijos[temp->n_hijos] ){
            temp->n_hijos++;
        }
    }

    /// S
    if( (resp & (unsigned char)32) && (*(skl_tmp + x + (y+1)*(cols+2)) < 2.0) ){
        if( temp->pix_tipo == PIX_CROSS || temp->pix_tipo == PIX_BRANCH ){
            *nivel = *nivel + 1;
        }
        temp->hijos[temp->n_hijos] = grafoSkeleton(skl_tmp, x, y+1, nivel, lutabla, visitados);

        if( temp->hijos[temp->n_hijos] ){
            temp->n_hijos++;
        }
    }

    /// SO
    if( (resp & (unsigned char)64) && (*(skl_tmp + (x-1) + (y+1)*(cols+2)) < 2.0) ){
        if( temp->pix_tipo == PIX_CROSS || temp->pix_tipo == PIX_BRANCH ){
            *nivel = *nivel + 1;
        }
        temp->hijos[temp->n_hijos] = grafoSkeleton(skl_tmp, x-1, y+1, nivel, lutabla, visitados);

        if( temp->hijos[temp->n_hijos] ){
            temp->n_hijos++;
        }
    }

    /// O
    if( (resp & (unsigned char)128) && (*(skl_tmp + (x-1) + y*(cols+2)) < 2.0)){
        if( temp->pix_tipo == PIX_CROSS || temp->pix_tipo == PIX_BRANCH ){
            *nivel = *nivel + 1;
        }
        temp->hijos[temp->n_hijos] = grafoSkeleton(skl_tmp, x-1, y, nivel, lutabla, visitados);

        if( temp->hijos[temp->n_hijos] ){
            temp->n_hijos++;
        }
    }

    /// Revisar si en verdad sigue siendo del tipo que creia ser:
    switch( temp->n_hijos ){
        case 0: temp->pix_tipo = PIX_END;
                break;
        case 1: temp->pix_tipo = PIX_SKL;
                break;
        case 2: temp->pix_tipo = PIX_BRANCH;
                break;
        case 3: temp->pix_tipo = PIX_CROSS;
                break;
    }

    return temp;
}




/*  Metodo: extraerCaract
    Funcion: Extrae los pixeles caracteristicos (end y branch points) a partir del esqueleot de la imagen.
*/
void IMGVTK::extraerCaract( IMG_IDX img_idx ){
    if( !borders_ptr ){
        detectarBorde( img_idx );
    }

    const unsigned char tabla[] = {
        0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        1, 0, 0, 0, 0, 2, 0, 0, 0, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        1, 0, 0, 0, 0, 2, 0, 0, 0, 2, 2, 0, 0, 0, 0, 0, 0, 0, 2, 0, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 2, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 2, 0, 2, 0, 0, 0, 2, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    double *skl_tmp = new double [(rows+2)*(cols+2)];
    memcpy( skl_tmp, skl_ptr, (rows+2)*(cols+2)*sizeof(double));


    /// Buscar un punto 'end' del esqueleto y empezar a generar el grafo a aprtir de ahi.
    bool *visitados = new bool [rows_cols];
    memset(visitados, 0, rows_cols*sizeof(bool));

    int x_ini, y_ini, xy = cols+2;
    unsigned char resp;
    do{
        xy++;
        x_ini = xy % (cols+2);
        y_ini = (int) (double)xy / ((double)cols+2);
        resp = sklMask( skl_ptr, x_ini, y_ini, cols+2, rows) * (*(skl_ptr + xy) > 0.0);
    }while( tabla[resp] != (unsigned char)1 );

    *(skl_tmp + xy ) = 2.0;
    int nivel = 0;

    pix_caract = grafoSkeleton(skl_tmp, x_ini, y_ini, &nivel, tabla, visitados);

    DEB_MSG("Encontrados " << nivel << " niveles");

    n_niveles = nivel;

    delete [] visitados;
    delete [] skl_tmp;
}



/*  Metodo: borrarSkeleton
    Funcion: Borra todos los nodos hijo de este nodo.
*/
void IMGVTK::borrarSkeleton( PIX_PAR *raiz ){

    if( raiz->n_hijos ){
        for( int i = 0; i < raiz->n_hijos; i++ ){
            borrarSkeleton( raiz->hijos[i] );
            delete raiz->hijos[i];
        }
    }
}




// C L A S E: IMGVTK  ---------------------------------------------------------------------------------------------- v
//----------------------------------------------------------------------------- PUBLIC ------- v
    // M E T O D O S      P U B L I C O S
/*
 *  Fuente: A Fast Parallel Algorithm for thinning digital patterns
 *  De:    T. Y. ZHANG and C. Y. SUEN
 *
*/
/*  Metodo: skeletonization
    Funcion: Obtiene el esqueleto de la imagen segmentada.
*/
void IMGVTK::skeletonization(IMG_IDX img_idx){

    double *img_ptr = NULL;

    switch( img_idx ){
        case BASE:
            img_ptr = base_ptr;
            break;
        case MASK:
            img_ptr = mask_ptr;
            break;
        case SKELETON:
            img_ptr = skl_ptr;
            break;
        case SEGMENT:
            img_ptr = segment_ptr;
            break;
        case THRESHOLD:
            img_ptr = threshold_ptr;
            break;
    }

    if( !skl_ptr ){
        skeleton = vtkSmartPointer<vtkImageData>::New();

        skeleton->SetExtent(0, cols+1, 0, rows+1, 0, 0);
        skeleton->AllocateScalars( VTK_DOUBLE, 1);
        skeleton->SetOrigin(0.0, 0.0, 0.0);
        skeleton->SetSpacing(1.0, 1.0, 1.0);

        skl_ptr = static_cast<double*>(skeleton->GetScalarPointer(0,0,0));
    }

    memset(skl_ptr, 0, (rows+2)*(cols+2)*sizeof(double));

    for( int y = 0; y < rows; y++ ){
        memcpy( skl_ptr + (y+1)*(cols+2)+1, img_ptr + y*cols, cols*sizeof(double) );
    }

    double *skl_mark = new double [(rows+2)*(cols+2)];
    memcpy(skl_mark, skl_ptr, (rows+2)*(cols+2)*sizeof(double));

    const unsigned char tabla[] = {
        0, 0, 0, 1, 0, 0, 1, 3, 0, 0, 3, 1, 1, 0, 1, 3, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 2, 0, 3, 0, 3, 3,
        0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 3, 0, 2, 2,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        2, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 2, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 3, 0, 2, 0,
        0, 0, 3, 1, 0, 0, 1, 3, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
        3, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        2, 3, 1, 3, 0, 0, 1, 3, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        2, 3, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 3, 3, 0, 1, 0, 0, 0, 0, 2, 2, 0, 0, 2, 0, 0, 0};
    int n_borrado;

    do{
        n_borrado = 0;
        // Primer paso:
        for(int y = 1; y <= rows; y++){
            for(int x = 1; x <= cols; x++){
                if(skl_ptr[x+y*(cols+2)] > 0){
                    const unsigned char resp = tabla[ sklMask(skl_ptr, x, y, cols+2, rows) ];
                    if(resp == 1 || resp == 3 ){
                        skl_mark[x+y*(cols+2)] = 0.0;
                        n_borrado++;
                    }
                }
            }
        }

        memcpy(skl_ptr, skl_mark, (rows+2)*(cols+2)*sizeof(double));

        // Segundo paso:
        for(int y = 1; y <= rows; y++){
            for(int x = 1; x <= cols; x++){
                if(skl_ptr[x+y*(cols+2)] > 0){
                    unsigned char resp = tabla[ sklMask(skl_ptr, x, y, cols+2, rows) ];
                    if(resp == 2 || resp == 3 ){
                        skl_mark[x+y*(cols+2)] = 0.0;
                        n_borrado++;
                    }
                }
            }
        }

        memcpy(skl_ptr, skl_mark, (rows+1)*(cols+1)*sizeof(double));

        DEB_MSG("Borrados: " << n_borrado);

    }while( n_borrado > 0 );

    delete [] skl_mark;



    extraerCaract( img_idx );
}



/*  Metodo: definirMask
    Funcion: Define una mascara para normalizar los pixeles de la imagen.
*/
void IMGVTK::definirMask( vtkSmartPointer<vtkImageData> img_src, vtkSmartPointer<vtkImageData> mask_src ){

    double *img_tmp = static_cast<double*>(img_src->GetScalarPointer(0,0,0));

    int dims[3];

    img_src->GetDimensions( dims );

    const int mis_cols = dims[0];
    const int mis_rens = dims[1];
    mask_src->SetExtent(0, mis_cols-1, 0, mis_rens-1, 0, 0);
    mask_src->AllocateScalars( VTK_DOUBLE, 1);
    mask_src->SetOrigin(0.0, 0.0, 0.0);
    mask_src->SetSpacing(1.0, 1.0, 1.0);

    double *mask_tmp = static_cast<double*>(mask_src->GetScalarPointer(0,0,0));
    maskFOV( img_tmp, mask_tmp, mis_cols, mis_rens );
    fillMask( img_tmp, mask_tmp, mis_cols, mis_rens );

    char mensaje[] = "Mascara generada exitosamente.\n";
    escribirLog( mensaje );
}





/*  Metodo: umbralizarOTSU
    Funcion: Utiliza el metodo de Otsu para umbralizar la imagen y separar el fondo y el primer plano de la imagen.
*/
double IMGVTK::umbralizarOTSU( const double *img_ptr, const double min, const double max){
    const int clases = 256;//(int)(1.0 + 3.322*log10(rens_cols));
    double *histograma_frecuencias = (double*) calloc(clases, sizeof(double));

    // Obtener el histograma de frecuencias a cada nivel de grises de la imagen y las frecuencias acumuladas:
    double suma = 0.0;
    const double fraccion = 1.0 / (double)rows_cols;

    for( int xy = 0; xy < rows_cols; xy ++){
        const int clase_i = (int)((clases-1) * (*(img_ptr + xy) - min)/(max - min));
        histograma_frecuencias[ clase_i ] += fraccion;
        suma += (double)(clase_i+1);
    }

    suma *= fraccion;

    double suma_back = 0;
    double peso_back = 0.0;
    double varianza_entre, max_varianza_entre = -1.0;
    double umbral;

    for( int k = 0; k < clases; k++){
        // Calcular el peso del back y fore-ground:
        peso_back += (double)histograma_frecuencias[k];

        // Calcular la media del back y fore-ground:
        suma_back += (double)(histograma_frecuencias[k]*(k+1));

        // Calcular la varianza entre el fore y back ground:
        varianza_entre = (suma*peso_back - suma_back);
        varianza_entre *= varianza_entre;
        varianza_entre /= (peso_back*(1.0 - peso_back));

        if( varianza_entre > max_varianza_entre ){
            max_varianza_entre = varianza_entre;
            umbral = ((double)k / (double)(clases-1));
        }
    }

    free( histograma_frecuencias );
    DEB_MSG("Umbral: " << umbral);

    return umbral;
}





/*  Metodo: umbralizarRIDCAL
    Funcion: Utiliza el metodo de Otsu para umbralizar la imagen y separar el fondo y el primer plano de la imagen.
*/
double IMGVTK::umbralizarRIDCAL( const double *img_ptr, const double min, const double max){
    double umbral_nuevo = 0.0;
    // Obtener la media de la intensidad de la imagen como primera aproximacion del umbral
    const double rango = 1.0 / (max - min);
    const double fraccion = 1.0 / (double)rows_cols;
    for( int xy = 0; xy < rows_cols; xy ++){
        umbral_nuevo += *(img_ptr + xy) * rango;
    }
    umbral_nuevo *= fraccion;

    double umbral;
    while( 1 ){
        const double umbral_previo = umbral_nuevo;
        umbral = umbral_previo;

        double m_abajo = 0.0;
        double m_arriba = 0.0;

        for( int xy = 0; xy < rows_cols; xy++){
            if( *(img_ptr + xy) <= umbral ){
                m_abajo += *(img_ptr + xy);
            }else{
                m_arriba += *(img_ptr + xy);
            }
        }

        m_abajo *= fraccion;
        m_arriba *= fraccion;

        umbral_nuevo = (m_arriba + m_abajo) / 2.0;

        if( fabs(umbral_nuevo - umbral_previo) < 1e-5 ){
            break;
        }
    }

    DEB_MSG("Umbral: " << umbral);

    return umbral;
}




/*  Metodo: umbralizar
    Funcion: Utiliza el metodo de Otsu o Ridler & Calvard para umbralizar la imagen y separar el fondo y el primer plano de la imagen.
*/
void IMGVTK::umbralizar(IMG_IDX img_idx, const TIPO_UMBRAL tipo_umb, const double nivel){

    if( !threshold_ptr ){
        threshold = vtkSmartPointer<vtkImageData>::New();

        threshold->SetExtent(0, cols-1, 0, rows-1, 0, 0);
        threshold->AllocateScalars( VTK_DOUBLE, 1);
        threshold->SetOrigin(0.0, 0.0, 0.0);
        threshold->SetSpacing(1.0, 1.0, 1.0);

        threshold_ptr = static_cast<double*>(threshold->GetScalarPointer(0,0,0));
    }

    double *img_ptr = NULL;
    switch( img_idx ){
        case BASE:
            img_ptr = base_ptr;
            break;
        case MASK:
            img_ptr = mask_ptr;
            break;
        case SEGMENT:
            img_ptr = segment_ptr;
            break;
    }

    double max =-1e100;
    double min = 1e100;

    for( int xy = 0; xy < rows_cols; xy++){
        if( *(img_ptr + xy) < min){
            min = *(img_ptr + xy);
        }
        if( *(img_ptr + xy) > max){
            max = *(img_ptr + xy);
        }
    }

    double umbral = 0.0;

    switch( tipo_umb ){
        case NIVEL:
            umbral = nivel;
            break;
        case OTSU:
            umbral = umbralizarOTSU( img_ptr, min, max);
            break;
        case RIDLER_CALVARD:
            umbral = umbralizarRIDCAL( img_ptr, min, max);
            break;
    }

    // Se umbraliza la imagen con el valor optimo encontrado:
    for( int xy = 0; xy < rows_cols; xy++){
        *(threshold_ptr + xy) = ( ((*(img_ptr + xy) - min) / (max - min)) >= umbral) ? 1.0 : 0.0;
    }
}



/*  Metodo: Cargar
    Funcion: Cargar desde un archivo DICOM/JPEG/PNG/BMP la imagen a formato VTK.
*/
void IMGVTK::Cargar(const char *ruta_origen, vtkSmartPointer<vtkImageData> img_src, vtkSmartPointer<vtkImageData> mask_src, const int nivel, const bool enmascarar){

    const int ruta_l = strlen(ruta_origen);

    if( ruta_l <= 0 ){
        escribirLog( "\n<<Error: Hay un problema con la ruta de entrada, tiene longitud 0 >>\n\n" );
        return;
    }

DEB_MSG("Extension del archivo de entrada: " << (ruta_origen + ruta_l - 3));
    esDICOM = true;
    esDICOM *= strcmp(ruta_origen + ruta_l - 3, "png");
    esDICOM *= strcmp(ruta_origen + ruta_l - 3, "jpg");
    esDICOM *= strcmp(ruta_origen + ruta_l - 4, "jpeg");
    esDICOM *= strcmp(ruta_origen + ruta_l - 3, "bmp");

    if( esDICOM ){
        gdcm::ImageReader DICOMreader;
        DICOMreader.SetFileName( ruta_origen );

        DICOMreader.Read();

        gdcm::File &file = DICOMreader.GetFile();
        gdcm::DataSet &ds = file.GetDataSet();

////---------- Extraer SID (Distancia de la fuente al detector): -----------------------------------------
        {
            const gdcm::DataElement &de = ds.GetDataElement( gdcm::Tag(0x18, 0x1110) );
            const gdcm::ByteValue *bv = de.GetByteValue();
            if( bv ){
                std::string strm(bv->GetPointer(), bv->GetLength());
                SID = atof( strm.c_str() );
            }
        }
////---------- Extraer SOD (Distancia de la fuente al paciente): -----------------------------------------
        {
            const gdcm::DataElement &de = ds.GetDataElement( gdcm::Tag(0x18, 0x1111) );
            const gdcm::ByteValue *bv = de.GetByteValue();
            if( bv ){
                std::string strm(bv->GetPointer(), bv->GetLength());
                SOD = atof( strm.c_str() );
            }
        }
////---------- Calcular el DDP (Distancia del detector al paciente): ---------------------------------------------------------------------
        DDP = SID - SOD;
////---------- Calcular el Magnification Factor: -----------------------------------------
        const double Magnification = SID / SOD;

////---------- Extraer pixY\pixX (Distancia entre el centro de los pixeles en el eje Y\X): -----------------------------------------
        {
            const gdcm::DataElement &de = ds.GetDataElement( gdcm::Tag(0x18, 0x1164) );
            const gdcm::ByteValue *bv = de.GetByteValue();
            if( bv ){
                std::string strm(bv->GetPointer(), bv->GetLength());
DEB_MSG("pixXY: " << strm);
                char *pixXYstr = new char [bv->GetLength()];
                memcpy(pixXYstr, strm.c_str(), bv->GetLength() * sizeof(char ));
                char *tmp = strchr(pixXYstr,'\\');
                pixXYstr[ tmp - pixXYstr ] = '\0';
                pixY = atof(pixXYstr);// / Magnification;
                pixX = atof(tmp+1);// / Magnification;
DEB_MSG("pixY: " << pixY << ", pixX: " << pixX);
                delete [] pixXYstr;
            }
        }
////---------- Extraer LAO/RAO (Angulo del detector en direccion izquierda(-) a derecha(+)): -----------------------------------------
        {
            const gdcm::DataElement &de = ds.GetDataElement( gdcm::Tag(0x18, 0x1510) );
            const gdcm::ByteValue *bv = de.GetByteValue();
            if( bv ){
                std::string strm(bv->GetPointer(), bv->GetLength());
                LAORAO = atof( strm.c_str() );
DEB_MSG("LAO: " << LAORAO);
            }
        }
////---------- Extraer CRA/CAU (Angulo del detector en direccion cranial(-) a caudal(+)): -----------------------------------------
        {
            const gdcm::DataElement &de = ds.GetDataElement( gdcm::Tag(0x18, 0x1511) );
            const gdcm::ByteValue *bv = de.GetByteValue();
            if( bv ){
                std::string strm(bv->GetPointer(), bv->GetLength());
                CRACAU = atof( strm.c_str() );
DEB_MSG("CRA: " << CRACAU);
            }
        }

////---------- Extraer Window Width (Ancho de la ventana donde se encuentran los valores de interes): -----------------------------------------
        {
            const gdcm::DataElement &de = ds.GetDataElement( gdcm::Tag(0x18, 0x7030) );
            const gdcm::ByteValue *bv = de.GetByteValue();
            if( bv ){
                std::string strm(bv->GetPointer(), bv->GetLength());
DEB_MSG("Field of view origin: " << strm);
            }
        }

////---------- Extraer Source to Isocenter (Distancia de la fuente al isocentro): -----------------------------------------
        {
            const gdcm::DataElement &de = ds.GetDataElement( gdcm::Tag(0x21, 0x1017) );
            const gdcm::ByteValue *bv = de.GetByteValue();
            if( bv ){
                gdcm::Element<gdcm::VR::SL, gdcm::VM::VM1_n> el;
                el.Set( de.GetValue() );
                const double SISO = el.GetValue();
                // Restar de la distancia de la fuente al detector (SID) la distancia del detector al isocentro para obtener la distancia del detector al isocentro (DISO).
                DISO = SID - SISO;
            }
        }
////---------- Extraer Window Center (Centro de los vlaores de interes): -----------------------------------------
        {
            const gdcm::DataElement &de = ds.GetDataElement( gdcm::Tag(0x28, 0x1050) );
            const gdcm::ByteValue *bv = de.GetByteValue();
            if( bv ){
                std::string strm(bv->GetPointer(), bv->GetLength());
                WCenter = atof( strm.c_str() );
DEB_MSG("Window Center: " << WCenter);
            }
        }
////---------- Extraer Window Width (Ancho de la ventana donde se encuentran los valores de interes): -----------------------------------------
        {
            const gdcm::DataElement &de = ds.GetDataElement( gdcm::Tag(0x28, 0x1051) );
            const gdcm::ByteValue *bv = de.GetByteValue();
            if( bv ){
                std::string strm(bv->GetPointer(), bv->GetLength());
                WWidth = atof( strm.c_str() );
DEB_MSG("Window Width: " << WWidth);
            }
        }
////---------- Extraer el ECG: -----------------------------------------
        int ecg_dim = 0, ecg_np = 0;
        {
            const gdcm::DataElement &de = ds.GetDataElement( gdcm::Tag(0x5000, 0x0005) );
            const gdcm::ByteValue *bv = de.GetByteValue();
            if( bv ){
                gdcm::Element<gdcm::VR::US, gdcm::VM::VM1_n> el;
                el.Set( de.GetValue() );
                ecg_dim = el.GetValue();
DEB_MSG("[ECG] Dimensiones: " << ecg_dim );
            }
        }
        {
            const gdcm::DataElement &de = ds.GetDataElement( gdcm::Tag(0x5000, 0x0010) );
            const gdcm::ByteValue *bv = de.GetByteValue();
            if( bv ){
                gdcm::Element<gdcm::VR::US, gdcm::VM::VM1_n> el;
                el.Set( de.GetValue() );
                ecg_np = el.GetValue();
DEB_MSG("[ECG] Numero de puntos: " << ecg_np );
            }
        }
        int n_niveles = 1;
        {
            const gdcm::DataElement &de = ds.GetDataElement( gdcm::Tag(0x8, 0x2143) );
            const gdcm::ByteValue *bv = de.GetByteValue();
            if( bv ){
                std::string strm(bv->GetPointer(), bv->GetLength());
                n_niveles = atof( strm.c_str() );
DEB_MSG("Numero de fotogramas: " << n_niveles);
            }
        }
        {
            const gdcm::DataElement &de = ds.GetDataElement( gdcm::Tag(0x5000, 0x3000) );
            const gdcm::ByteValue *bv = de.GetByteValue();
            if( bv ){
                gdcm::Element<gdcm::VR::OW, gdcm::VM::VM1_n> el;
                el.Set( de.GetValue() );
                char *nombre_ecg = new char [18];
                strcpy(nombre_ecg, ruta_origen + ruta_l - 7);
                sprintf( nombre_ecg, "%s_%i.dat", nombre_ecg, n_niveles);
                FILE *fp = fopen(nombre_ecg, "w");

                for( int i = 0; i < ecg_np; i++){
                    fprintf(fp, "%i %i\n", i, el.GetValue(i));
                }
                fclose( fp );
                delete [] nombre_ecg;
            }
        }
/*
            std::string num_points_str = file.GetEntryString(0x5000,0x0010);
            unsigned short num_points;
            convert.clear();
            convert.str(num_points_str);
            convert >> num_points;
DEB_MSG("Number of Points: " << num_points);

            std::string data_type = file.GetEntryString(0x5000,0x0020);
DEB_MSG("Type of Data: " << data_type);

            std::string curve_desc = file.GetEntryString(0x5000,0x0022);
DEB_MSG("Curve Description: " << curve_desc);

            std::string data_rep_str = file.GetEntryString(0x5000,0x0103);
            unsigned short data_rep;
            convert.clear();
            convert.str(data_rep_str);
            convert >> data_rep;

            gdcm::DocEntry *pCurveDataDoc = file.GetDocEntry(0x5000, 0x3000);
            gdcm::DataEntry *pCurveData = dynamic_cast<gdcm::DataEntry *>(pCurveDataDoc);
            uint8_t *curve_data = pCurveData->GetBinArea();

            for(int i = 0; i < num_points; i++){
DEB_MSG("Pt(" << i <<  ") = " << ((unsigned short*)curve_data)[i]);
            }
*/

        ///----------------------------------------------------- Leer imagenes
        DICOMreader.Read();
        const gdcm::Image &gimage = DICOMreader.GetImage();
DEB_MSG("Buffer length: " << gimage.GetBufferLength());
        char *buffer = new char[gimage.GetBufferLength()];
        gimage.GetBuffer(buffer);

        const unsigned int* dimension = gimage.GetDimensions();
        const int mis_cols = dimension[0];
        const int mis_rens = dimension[1];
        const int mis_niveles = dimension[2];
        const int mis_rens_cols = mis_rens*mis_cols;

DEB_MSG("Cols: " << mis_cols << ", rens: " << mis_rens  << ", niveles: " << mis_niveles);

        // Alojar memoria para la imagen:
        img_src->SetExtent(0, mis_cols-1, 0, mis_rens-1, 0, 0);
        img_src->AllocateScalars( VTK_DOUBLE, 1);
        img_src->SetOrigin(0.0, 0.0, 0.0);
        img_src->SetSpacing(1.0, 1.0, 1.0);

        double *img_tmp = static_cast<double*>(img_src->GetScalarPointer(0,0,0));

        mask_src->SetExtent(0, mis_cols-1, 0, mis_rens-1, 0,0);
        mask_src->AllocateScalars( VTK_DOUBLE, 1);
        mask_src->SetOrigin(0.0, 0.0, 0.0);
        mask_src->SetSpacing(1.0, 1.0, 1.0);

        gdcm::PhotometricInterpretation scl_comps = gimage.GetPhotometricInterpretation();
        gdcm::PixelFormat pix_format = gimage.GetPixelFormat();

        switch(scl_comps){
            case gdcm::PhotometricInterpretation::RGB:{
DEB_MSG("Imagen DICOM en RGB...");
                    if( pix_format == gdcm::PixelFormat::UINT8 ){

                        img_tmp = static_cast<double*>(img_src->GetScalarPointer(0,0,0));
                        for( int y = 0; y < mis_rens; y++){
                            for( int x = 0; x < mis_cols; x++){
                                const double pixR = (double)(unsigned char)*(buffer + 3*x   + y*mis_cols*3 + nivel*mis_rens_cols*3) - WCenter + 0.5;
                                const double pixG = (double)(unsigned char)*(buffer + 3*x+1 + y*mis_cols*3 + nivel*mis_rens_cols*3) - WCenter + 0.5;
                                const double pixB = (double)(unsigned char)*(buffer + 3*x+2 + y*mis_cols*3 + nivel*mis_rens_cols*3) - WCenter + 0.5;
                                double pix = (0.297)*pixR + (0.589)*pixG + (0.114)*pixB;
                                if( pix <= -((WWidth-1) / 2)){
                                    pix = 0.0;
                                }else if(pix > ((WWidth-1) / 2)){
                                    pix = 0.0;
                                }else{
                                    pix = pix / (WWidth -1) + 0.5;
                                }
                                *(img_tmp + (mis_rens-y-1)*mis_rens + x ) = pix; // 255.0;
                            }
                        }
                    }else{
                        escribirLog( COLOR_BACK_BLACK COLOR_RED "< ERROR AL LEER ARCHIVO DICOM: Formato de imagen RGB no soportado >\n" COLOR_NORMAL);
                    }
                    break;
                }
            case gdcm::PhotometricInterpretation::MONOCHROME1:
            case gdcm::PhotometricInterpretation::MONOCHROME2:{
DEB_MSG("Imagen DICOM en escala de grises...");
                    if( pix_format == gdcm::PixelFormat::UINT8 ){
DEB_MSG("Tipo UINT8");

                        img_tmp = static_cast<double*>(img_src->GetScalarPointer(0,0,0));
                        for( int y = 0; y < mis_rens; y++){
                            for( int x = 0; x < mis_cols; x++){
                                double pix = (double)(unsigned char)*(buffer + nivel*mis_rens_cols + x + y*mis_cols) - WCenter + 0.5;

                                if( pix <= -((WWidth-1) / 2)){
                                    pix = 0.0;
                                }else if(pix > ((WWidth-1) / 2)){
                                    pix = 1.0;
                                }else{
                                    pix = pix / (WWidth -1) + 0.5;
                                }
                                *(img_tmp + (mis_rens-y-1)*mis_rens + x) = pix; // 255.0;
                            }
                        }

                    }else if( pix_format == gdcm::PixelFormat::UINT16 ){
DEB_MSG("Tipo UINT16");
                        unsigned short *buffer16 = (unsigned short*)buffer;

                        img_tmp = static_cast<double*>(img_src->GetScalarPointer(0,0,0));
                        for( int y = 0; y < mis_rens; y++){
                            for( int x = 0; x < mis_cols; x++){
                                const double pixR = (double)((unsigned char)*(buffer16 + 3*x   + y*mis_cols*3 + nivel*mis_rens_cols*3) / 16)  - WCenter + 0.5;
                                const double pixG = (double)((unsigned char)*(buffer16 + 3*x+1 + y*mis_cols*3 + nivel*mis_rens_cols*3) / 16)  - WCenter + 0.5;
                                const double pixB = (double)((unsigned char)*(buffer16 + 3*x+2 + y*mis_cols*3 + nivel*mis_rens_cols*3) / 16)  - WCenter + 0.5;
                                double pix = (0.297)*pixR + (0.589)*pixG + (0.114)*pixB;
                                if( pix <= -((WWidth-1) / 2)){
                                    pix = 0.0;
                                }else if(pix > ((WWidth-1) / 2)){
                                    pix = 0.0;
                                }else{
                                    pix = pix / (WWidth -1) + 0.5;
                                }
                                *(img_tmp + (mis_rens-y-1)*mis_rens + x) =  pix; // 255.0;
                            }
                        }

                    }else{
                        escribirLog( COLOR_BACK_BLACK COLOR_RED "< ERROR AL LEER ARCHIVO DICOM: Formato de imagen RGB no soportado >\n" COLOR_NORMAL);
                    }
                    break;
                }
        }
        delete [] buffer;


    }else{ ///------------------------------------------------------------------------------------------------



        // Leer la imagen en RBG o Escala de Grises:
        vtkSmartPointer<vtkImageReader2Factory> readerFactory = vtkSmartPointer<vtkImageReader2Factory>::New();
        vtkSmartPointer<vtkImageReader2> imgReader = readerFactory->CreateImageReader2( ruta_origen );

        imgReader->SetFileName( ruta_origen );
        imgReader->Update();

        // Pasar temporalmente la imagen de RGB a escala de grises en la imagen de entrada.
        int dims[3];
        imgReader->GetOutput()->GetDimensions( dims );

        int scl_comps = imgReader->GetOutput()->GetNumberOfScalarComponents();

        const int mis_cols = dims[0];
        const int mis_rens = dims[1];
        const int mis_rens_cols = mis_rens*mis_cols;

        // Alojar memoria para la imagen:
        img_src->SetExtent(0, mis_cols-1, 0, mis_rens-1, 0, 0);
        img_src->AllocateScalars( VTK_DOUBLE, 1);
        img_src->SetOrigin(0.0, 0.0, 0.0);
        img_src->SetSpacing(1.0, 1.0, 1.0);

        double *img_tmp = static_cast<double*>(img_src->GetScalarPointer(0,0,0));


        switch(scl_comps){
            // La imagen esta en escala de grises:
            case 1:
            case 2:{
				DEB_MSG("La imagen esta en escala de grises");
                vtkSmartPointer<vtkImageExtractComponents> extractGreyFilter = vtkSmartPointer<vtkImageExtractComponents>::New();
                extractGreyFilter->SetInputConnection(imgReader->GetOutputPort());
                extractGreyFilter->SetComponents(0);
                extractGreyFilter->Update();

                unsigned char *gris = static_cast<unsigned char*>(extractGreyFilter->GetOutput()->GetScalarPointer(0,0,0));

                for( int y = 0; y < mis_rens; y++){
                    for( int x = 0; x < mis_cols; x++){
                        *(img_tmp + (mis_rens - 1 - y)*mis_cols + x ) = (double)*(gris + y*mis_cols + x) / 255.0;
                    }
                }
                break;
            }
            // La imagen esta en RGB:
            case 3:{
				DEB_MSG("La imagen esta en RGB");
                vtkSmartPointer<vtkImageExtractComponents> extractRedFilter = vtkSmartPointer<vtkImageExtractComponents>::New();
                extractRedFilter->SetInputConnection(imgReader->GetOutputPort());
                extractRedFilter->SetComponents(0);
                extractRedFilter->Update();

                vtkSmartPointer<vtkImageExtractComponents> extractGreenFilter = vtkSmartPointer<vtkImageExtractComponents>::New();
                extractGreenFilter->SetInputConnection(imgReader->GetOutputPort());
                extractGreenFilter->SetComponents(1);
                extractGreenFilter->Update();

                vtkSmartPointer<vtkImageExtractComponents> extractBlueFilter = vtkSmartPointer<vtkImageExtractComponents>::New();
                extractBlueFilter->SetInputConnection(imgReader->GetOutputPort());
                extractBlueFilter->SetComponents(2);
                extractBlueFilter->Update();

                vtkImageData *ptR = extractRedFilter->GetOutput();
                vtkImageData *ptG = extractGreenFilter->GetOutput();
                vtkImageData *ptB = extractBlueFilter->GetOutput();

                double color;

                for( int y = 0 ; y < mis_rens; y++){
                    for( int x = 0; x < mis_cols; x++ ){
                        color = ((0.297)*ptR->GetScalarComponentAsDouble(x,y,0,0) +
                                 (0.589)*ptG->GetScalarComponentAsDouble(x,y,0,0) +
                                 (0.114)*ptB->GetScalarComponentAsDouble(x,y,0,0));

                        *(img_tmp + (mis_rens - 1 - y)*mis_cols + x) = color / 255.0;
                    }
                }
                break;
            }
            // La imagen esta en RGBA:
            case 4:{
				DEB_MSG("La Imagen esta en RGB{A}");
                vtkSmartPointer<vtkImageExtractComponents> extractRedFilter = vtkSmartPointer<vtkImageExtractComponents>::New();
                extractRedFilter->SetInputConnection(imgReader->GetOutputPort());
                extractRedFilter->SetComponents(0);
                extractRedFilter->Update();

                vtkSmartPointer<vtkImageExtractComponents> extractGreenFilter = vtkSmartPointer<vtkImageExtractComponents>::New();
                extractGreenFilter->SetInputConnection(imgReader->GetOutputPort());
                extractGreenFilter->SetComponents(1);
                extractGreenFilter->Update();

                vtkSmartPointer<vtkImageExtractComponents> extractBlueFilter = vtkSmartPointer<vtkImageExtractComponents>::New();
                extractBlueFilter->SetInputConnection(imgReader->GetOutputPort());
                extractBlueFilter->SetComponents(2);
                extractBlueFilter->Update();

                vtkSmartPointer<vtkImageExtractComponents> extractAlphaFilter = vtkSmartPointer<vtkImageExtractComponents>::New();
                extractAlphaFilter->SetInputConnection(imgReader->GetOutputPort());
                extractAlphaFilter->SetComponents(3);
                extractAlphaFilter->Update();

                vtkImageData *ptR = extractRedFilter->GetOutput();
                vtkImageData *ptG = extractGreenFilter->GetOutput();
                vtkImageData *ptB = extractBlueFilter->GetOutput();
                //vtkImageData *ptA = extractAlphaFilter->GetOutput();
                double color;

                for( int y = 0 ; y < mis_rens; y++){
                    for( int x = 0; x < mis_cols; x++ ){
                        color = ((0.297)*ptR->GetScalarComponentAsDouble(x,y,0,0) +
                                 (0.589)*ptG->GetScalarComponentAsDouble(x,y,0,0) +
                                 (0.114)*ptB->GetScalarComponentAsDouble(x,y,0,0));

                        *(img_tmp + (mis_rens - 1 - y)*mis_cols + x) = color / 255.0;
                    }
                }
                break;
            }
        }

        imgReader->Delete();
    }

    char mensaje[512];
    sprintf( mensaje, "Imagen " COLOR_GREEN "'%s'" COLOR_NORMAL " cargada exitosamente.\n", ruta_origen );
    escribirLog( mensaje );

    if(enmascarar){
        definirMask(img_src, mask_src);
    }
}



/*  Metodo: Cargar
    Funcion: Cargar la imagen a formato VTK desde varios archivos (Solo para imagenes JPEG, BMP o PNG).
*/
void IMGVTK::Cargar(vtkSmartPointer<vtkImageData> img_src, vtkSmartPointer<vtkImageData> mask_src, char **rutas, const int n_imgs, const bool enmascarar){
    vtkSmartPointer<vtkImageData> auxiliar = vtkSmartPointer<vtkImageData>::New();
    vtkSmartPointer<vtkImageData> mask_auxiliar = vtkSmartPointer<vtkImageData>::New();

    Cargar(rutas[0], auxiliar, mask_auxiliar, 0, enmascarar);

    double *aux_tmp = static_cast<double*>(auxiliar->GetScalarPointer(0, 0, 0));
    double *mskaux_tmp = static_cast<double*>(mask_auxiliar->GetScalarPointer(0, 0, 0));

    int dims[3];
    auxiliar->GetDimensions(dims);

    const int mis_cols = dims[0];
    const int mis_rens = dims[1];

    const int n_cols = mis_cols * n_imgs;

    // Extract data:
    img_src->SetExtent(0, n_cols-1, 0, mis_rens-1, 0, 0);
    img_src->AllocateScalars( VTK_DOUBLE, 1);
    img_src->SetOrigin(0.0, 0.0, 0.0);
    img_src->SetSpacing(1.0, 1.0, 1.0);

    mask_src->SetExtent(0, n_cols-1, 0, mis_rens-1, 0, 0);
    mask_src->AllocateScalars( VTK_DOUBLE, 1);
    mask_src->SetOrigin(0.0, 0.0, 0.0);
    mask_src->SetSpacing(1.0, 1.0, 1.0);

    double *base_ptr_tmp = static_cast<double*>(img_src->GetScalarPointer(0, 0, 0));
    double *mask_ptr_tmp = static_cast<double*>(mask_src->GetScalarPointer(0, 0, 0));
    for( int y = 0; y < mis_rens; y++){
        memcpy(base_ptr_tmp + y*n_cols, aux_tmp + y*mis_cols, mis_cols*sizeof(double));
        if(enmascarar){
            memcpy(mask_ptr_tmp + y*n_cols, mskaux_tmp + y*mis_cols, mis_cols*sizeof(double));
        }
    }

    for( int img_i = 1; img_i < n_imgs; img_i++){
        Cargar(rutas[img_i], auxiliar, mask_auxiliar, 0, enmascarar);
        aux_tmp = static_cast<double*>(auxiliar->GetScalarPointer(0, 0, 0));
        mskaux_tmp = static_cast<double*>(mask_auxiliar->GetScalarPointer(0, 0, 0));

        for( int y = 0; y < mis_rens; y++){
            memcpy(base_ptr_tmp + y*n_cols + img_i*mis_cols, aux_tmp+ y*mis_cols, mis_cols*sizeof(double));
            if(enmascarar){
                memcpy(mask_ptr_tmp + y*n_cols + img_i*mis_cols, mskaux_tmp + y*mis_cols, mis_cols*sizeof(double));
            }
        }
    }
}



/*  Metodo: Cargar
    Funcion: Cargar la imagen a formato VTK desde un archivo.
*/
void IMGVTK::Cargar(const IMG_IDX img_idx, const char *ruta_origen, const bool enmascarar, const int nivel){

    DEB_MSG("Cargando imagen desde: " COLOR_GREEN << ruta_origen << COLOR_NORMAL);

    vtkSmartPointer<vtkImageData> img_temp = NULL;

    switch( img_idx){
    case BASE:
        if( !base ){
            base = vtkSmartPointer<vtkImageData>::New();
        }
        img_temp = base;

        if( !mask ){
            mask = vtkSmartPointer<vtkImageData>::New();
        }
        break;

    case GROUNDTRUTH:
        if( !ground ){
            ground = vtkSmartPointer<vtkImageData>::New();
        }
        img_temp = ground;
        break;
    }


    Cargar(ruta_origen, img_temp, mask, nivel, enmascarar);

    switch( img_idx ){
    case BASE:
        base_ptr = static_cast<double*>(base->GetScalarPointer(0, 0, 0));
        mask_ptr = static_cast<double*>(mask->GetScalarPointer(0, 0, 0));
        break;

    case GROUNDTRUTH:
        gt_ptr = static_cast<double*>(ground->GetScalarPointer(0, 0, 0));
        break;
    }

    if( !cols ){
        int dims[3];
        img_temp->GetDimensions( dims );

        cols = dims[0];
        rows = dims[1];
        rows_cols = rows*cols;
    }

    if( !segment ){
        segment = vtkSmartPointer<vtkImageData>::New();
        segment->SetExtent(0, cols-1, 0, rows-1, 0, 0);
        segment->AllocateScalars( VTK_DOUBLE, 1);
        segment->SetOrigin(0.0, 0.0, 0.0);
        segment->SetSpacing(1.0, 1.0, 1.0);

        segment_ptr = static_cast<double*>(segment->GetScalarPointer(0, 0, 0));
    }
}



/*  Metodo: Cargar
    Funcion: Cargar la imagen a formato VTK desde varios archivos.
*/
void IMGVTK::Cargar(const IMG_IDX img_idx, char **rutas_origen, const int n_imgs, const bool enmascarar){

    DEB_MSG("Cargando imagen desde: " COLOR_GREEN << rutas_origen << COLOR_NORMAL);

    vtkSmartPointer<vtkImageData> img_temp = NULL;

    switch( img_idx){
    case BASE:
        if( !base ){
            base = vtkSmartPointer<vtkImageData>::New();
        }
        img_temp = base;

        if( !mask ){
            mask = vtkSmartPointer<vtkImageData>::New();
        }
        break;

    case GROUNDTRUTH:
        if( !ground ){
            ground = vtkSmartPointer<vtkImageData>::New();
        }
        img_temp = ground;
        break;
    }

    Cargar( img_temp, mask, rutas_origen, n_imgs, enmascarar );

    switch( img_idx ){
    case BASE:
        base_ptr = static_cast<double*>(base->GetScalarPointer(0, 0, 0));
        mask_ptr = static_cast<double*>(mask->GetScalarPointer(0, 0, 0));
        break;

    case GROUNDTRUTH:
        gt_ptr = static_cast<double*>(ground->GetScalarPointer(0, 0, 0));
    }


    if( !cols ){
        int dims[3];
        img_temp->GetDimensions( dims );

        cols = dims[0];
        rows = dims[1];
        rows_cols = rows*cols;
    }

    if( !segment ){
        segment = vtkSmartPointer<vtkImageData>::New();
        segment->SetExtent(0, cols-1, 0, rows-1, 0, 0);
        segment->AllocateScalars( VTK_DOUBLE, 1);
        segment->SetOrigin(0.0, 0.0, 0.0);
        segment->SetSpacing(1.0, 1.0, 1.0);

        segment_ptr = static_cast<double*>(segment->GetScalarPointer(0, 0, 0));
    }
}




/*  Metodo: Guardar
    Funcion: Guarda la imagen en la ruta especificada con la extension especificada.
*/
void IMGVTK::Guardar(IMG_IDX img_idx, const char *ruta, const TIPO_IMG tipo_salida ){
    double *img_ptr = NULL;
    int offset_x = 0, offset_y = 0;

    switch( img_idx ){
    case BASE:
        img_ptr = base_ptr;
        break;
    case GROUNDTRUTH:
        img_ptr = gt_ptr;
        break;
    case MASK:
        img_ptr = mask_ptr;
        break;
    case SKELETON:
        img_ptr = skl_ptr;
        offset_x = 1;
        offset_y = 1;
        break;
    case SEGMENT:
        img_ptr = segment_ptr;
        break;
    case MAPDIST:
        img_ptr = map_ptr;
        break;
    case THRESHOLD:
        img_ptr = threshold_ptr;
        break;
    case BORDERS:
        img_ptr = borders_ptr;
        break;
    }

    // Alojar memoria para la imagen:
    vtkSmartPointer<vtkImageData> img_out = vtkSmartPointer<vtkImageData>::New();
    img_out->SetExtent(0, cols-1, 0, rows-1, 0, 0);
    img_out->AllocateScalars( VTK_UNSIGNED_CHAR, 1);
    img_out->SetOrigin(0.0, 0.0, 0.0);
    img_out->SetSpacing(1.0, 1.0, 1.0);
    unsigned char *img_out_ptr = static_cast<unsigned char*>(img_out->GetScalarPointer(0,0,0));

    double min = INF;
    double max =-INF;

    for( int y = 0; y < rows; y++ ){
        for( int x = 0; x < cols; x++ ){
            if( min > *(img_ptr + (x + offset_x) + (y + offset_y)*(cols + offset_x*2)) ){
                min = *(img_ptr + (x + offset_x) + (y + offset_y)*(cols + offset_x*2));
            }
            if( max < *(img_ptr + (x + offset_x) + (y + offset_y)*(cols + offset_x*2)) ){
                max = *(img_ptr + (x + offset_x) + (y + offset_y)*(cols + offset_x*2));
            }
        }
    }

    if( fabs(max - min) < 1e-12 ){
        max = 1.0;
        min = 0.0;
    }

    DEB_MSG(COLOR_BACK_GREEN COLOR_BACK_BLACK "GUARDAR" COLOR_NORMAL " min: " << min << ", " << max);

    for( int y = 0; y < rows; y++ ){
        for( int x = 0; x < cols; x++ ){
            *(img_out_ptr + x + y*cols) = (unsigned char)(255.0 * (*(img_ptr + (x + offset_x) + (y + offset_y)*(cols + offset_x*2)) - min) / (max - min));
        }
    }

    switch(tipo_salida){
        case PGM:{
            break;
        }

        case PNG:{
            vtkSmartPointer<vtkPNGWriter> png_output = vtkSmartPointer<vtkPNGWriter>::New();
            png_output->SetFileName( ruta );
            png_output->SetInputData( img_out );
            png_output->Write();
            break;
        }
    }
}



/*  Metodo: setLog

    Funcion: Define el editor donde se escribiran todos los logs del sistema
*/
void IMGVTK::setLog( QPlainTextEdit *log ){
    DEB_MSG( "Agregando log: " << log );
    mi_log = log;
}



/* CONSTRUCTORES */
IMGVTK::IMGVTK(){
    mi_log = NULL;

    cols = 0;
    rows = 0;
    max_dist = 0;

    base_ptr = NULL;
    gt_ptr = NULL;
    skl_ptr = NULL;
    mask_ptr = NULL;
    map_ptr = NULL;
    borders_ptr = NULL;
    pix_caract = NULL;
    segment_ptr = NULL;
    threshold_ptr = NULL;

    // Defaults:
    SID = 1100.0;
    SOD = 400.0;
    DDP = SID - SOD;
    DISO = SID / 2;
    LAORAO = 20.0;
    CRACAU = 20.0;
    pixX = 0.308;
    pixY = 0.308;
    WCenter = 127.5;
    WWidth = 255.0;
}



IMGVTK::IMGVTK( const IMGVTK &origen ){
    mi_log = origen.mi_log;
    max_dist = origen.max_dist;

    base_ptr = NULL;
    gt_ptr = NULL;
    skl_ptr = NULL;
    mask_ptr = NULL;
    map_ptr = NULL;
    borders_ptr = NULL;
    pix_caract = NULL;
    segment_ptr = NULL;
    threshold_ptr = NULL;

    // Defaults:
    SID = origen.SID;
    SOD = origen.SOD;
    DDP = SID - SOD;
    DISO = origen.DISO;
    LAORAO = origen.LAORAO;
    CRACAU = origen.CRACAU;
    pixX = origen.pixX;
    pixY = origen.pixY;
    WCenter = origen.WCenter;
    WWidth = origen.WWidth;


    cols = origen.cols;
    rows = origen.rows;
    rows_cols = rows * cols;

    if(origen.base_ptr){

        base = vtkSmartPointer<vtkImageData>::New();
        base->SetExtent(0, cols-1, 0, rows-1, 0, 0);
        base->AllocateScalars( VTK_DOUBLE, 1);
        base->SetOrigin(0.0, 0.0, 0.0);
        base->SetSpacing(1.0, 1.0, 1.0);

        base_ptr = static_cast<double*>(base->GetScalarPointer(0, 0, 0));
        memcpy(base_ptr, origen.base_ptr, rows_cols*sizeof(double));

        if(origen.gt_ptr){
            ground = vtkSmartPointer<vtkImageData>::New();
            ground->SetExtent(0, cols-1, 0, rows-1, 0, 0);
            ground->AllocateScalars( VTK_DOUBLE, 1);
            ground->SetOrigin(0.0, 0.0, 0.0);
            ground->SetSpacing(1.0, 1.0, 1.0);
            gt_ptr = static_cast<double*>(ground->GetScalarPointer(0, 0, 0));

            memcpy(gt_ptr, origen.gt_ptr, rows_cols*sizeof(double));
        }

        if(origen.mask_ptr){
            mask = vtkSmartPointer<vtkImageData>::New();
            mask->SetExtent(0, cols-1, 0, rows-1, 0, 0);
            mask->AllocateScalars( VTK_DOUBLE, 1);
            mask->SetOrigin(0.0, 0.0, 0.0);
            mask->SetSpacing(1.0, 1.0, 1.0);
            mask_ptr = static_cast<double*>(mask->GetScalarPointer(0, 0, 0));

            memcpy(mask_ptr, origen.mask_ptr, rows_cols*sizeof(double));
        }

        if(origen.segment_ptr){
            segment = vtkSmartPointer<vtkImageData>::New();
            segment->SetExtent(0, cols-1, 0, rows-1, 0, 0);
            segment->AllocateScalars( VTK_DOUBLE, 1);
            segment->SetOrigin(0.0, 0.0, 0.0);
            segment->SetSpacing(1.0, 1.0, 1.0);
            segment_ptr = static_cast<double*>(segment->GetScalarPointer(0, 0, 0));

            memcpy(segment_ptr, origen.segment_ptr, rows_cols*sizeof(double));
        }

        if(origen.skl_ptr){
            skeleton = vtkSmartPointer<vtkImageData>::New();
            skeleton->SetExtent(0, cols+1, 0, rows+1, 0, 0);
            skeleton->AllocateScalars( VTK_DOUBLE, 1);
            skeleton->SetOrigin(0.0, 0.0, 0.0);
            skeleton->SetSpacing(1.0, 1.0, 1.0);
            skl_ptr = static_cast<double*>(skeleton->GetScalarPointer(0, 0, 0));

            memcpy(skl_ptr, origen.skl_ptr, rows_cols*sizeof(double));
        }

        if(origen.map_ptr){
            mapa_dist = vtkSmartPointer<vtkImageData>::New();
            mapa_dist->SetExtent(0, cols-1, 0, rows-1, 0, 0);
            mapa_dist->AllocateScalars( VTK_DOUBLE, 1);
            mapa_dist->SetOrigin(0.0, 0.0, 0.0);
            mapa_dist->SetSpacing(1.0, 1.0, 1.0);
            map_ptr = static_cast<double*>(mapa_dist->GetScalarPointer(0, 0, 0));

            memcpy(map_ptr, origen.map_ptr, rows_cols*sizeof(double));
        }

        if(origen.borders_ptr){
            borders = vtkSmartPointer<vtkImageData>::New();
            borders->SetExtent(0, cols-1, 0, rows-1, 0, 0);
            borders->AllocateScalars( VTK_DOUBLE, 1);
            borders->SetOrigin(0.0, 0.0, 0.0);
            borders->SetSpacing(1.0, 1.0, 1.0);
            borders_ptr = static_cast<double*>(borders->GetScalarPointer(0, 0, 0));

            memcpy(borders_ptr, origen.borders_ptr, rows_cols*sizeof(double));
        }
    }
}



IMGVTK::IMGVTK( char **rutas_origen, const int n_imgs, const bool enmascarar){
    mi_log = NULL;
    max_dist = 0;

    base_ptr = NULL;
    gt_ptr = NULL;
    skl_ptr = NULL;
    mask_ptr = NULL;
    map_ptr = NULL;
    borders_ptr = NULL;
    pix_caract = NULL;
    segment_ptr = NULL;
    threshold_ptr = NULL;

    // Defaults:
    SID = 1100.0;
    SOD = 400.0;
    DDP = SID - SOD;
    DISO = SID / 2;
    LAORAO = 20.0;
    CRACAU = 20.0;
    pixX = 0.308;
    pixY = 0.308;
    WCenter = 127.5;
    WWidth = 255.0;


    cols = 0;
    rows = 0;
    max_dist = 0;

    Cargar(BASE, rutas_origen, n_imgs, enmascarar );
}


IMGVTK::IMGVTK( const char *ruta_origen, const bool enmascarar, const int nivel){
    mi_log = NULL;
    max_dist = 0;

    base_ptr = NULL;
    gt_ptr = NULL;
    skl_ptr = NULL;
    mask_ptr = NULL;
    map_ptr = NULL;
    borders_ptr = NULL;
    pix_caract = NULL;
    segment_ptr = NULL;
    threshold_ptr = NULL;

    // Defaults:
    SID = 1100.0;
    SOD = 400.0;
    DDP = SID - SOD;
    DISO = SID / 2;
    LAORAO = 20.0;
    CRACAU = 20.0;
    pixX = 0.308;
    pixY = 0.308;
    WCenter = 127.5;
    WWidth = 255.0;

    cols = 0;
    rows = 0;
    max_dist = 0;

    DEB_MSG("Cargando imagen desde: " COLOR_GREEN << ruta_origen << COLOR_NORMAL);

    Cargar(BASE, ruta_origen, enmascarar, nivel);
}




/* DESTRUCTOR */
IMGVTK::~IMGVTK(){
    if(pix_caract){
        /// Liberar memoria recursivamente:
        borrarSkeleton( pix_caract );
        delete pix_caract;
    }
}




// O P E R A D O R E S  S O B R E C A R G A D O S
// El operador de copia extrae unicamente el contenido de la imagen original
IMGVTK& IMGVTK::operator= ( const IMGVTK &origen ){
    mi_log = origen.mi_log;
    max_dist = origen.max_dist;

    base_ptr = NULL;
    gt_ptr = NULL;
    skl_ptr = NULL;
    mask_ptr = NULL;
    map_ptr = NULL;
    borders_ptr = NULL;
    pix_caract = NULL;
    segment_ptr = NULL;
    threshold_ptr = NULL;


    // Defaults:
    SID = origen.SID;
    SOD = origen.SOD;
    DDP = SID - SOD;
    DISO = origen.DISO;
    LAORAO = origen.LAORAO;
    CRACAU = origen.CRACAU;
    pixX = origen.pixX;
    pixY = origen.pixY;
    WCenter = origen.WCenter;
    WWidth = origen.WWidth;


    cols = origen.cols;
    rows = origen.rows;
    rows_cols = rows * cols;

    if(origen.base_ptr){

        base = vtkSmartPointer<vtkImageData>::New();
        base->SetExtent(0, cols-1, 0, rows-1, 0, 0);
        base->AllocateScalars( VTK_DOUBLE, 1);
        base->SetOrigin(0.0, 0.0, 0.0);
        base->SetSpacing(1.0, 1.0, 1.0);

        base_ptr = static_cast<double*>(base->GetScalarPointer(0, 0, 0));
        memcpy(base_ptr, origen.base_ptr, rows_cols*sizeof(double));

        if(origen.gt_ptr){
            mask = vtkSmartPointer<vtkImageData>::New();
            ground->SetExtent(0, cols-1, 0, rows-1, 0, 0);
            ground->AllocateScalars( VTK_DOUBLE, 1);
            ground->SetOrigin(0.0, 0.0, 0.0);
            ground->SetSpacing(1.0, 1.0, 1.0);
            gt_ptr = static_cast<double*>(ground->GetScalarPointer(0, 0, 0));

            memcpy(gt_ptr, origen.gt_ptr, rows_cols*sizeof(double));
        }


        if(origen.mask_ptr){
            mask = vtkSmartPointer<vtkImageData>::New();
            mask->SetExtent(0, cols-1, 0, rows-1, 0, 0);
            mask->AllocateScalars( VTK_DOUBLE, 1);
            mask->SetOrigin(0.0, 0.0, 0.0);
            mask->SetSpacing(1.0, 1.0, 1.0);
            mask_ptr = static_cast<double*>(mask->GetScalarPointer(0, 0, 0));

            memcpy(mask_ptr, origen.mask_ptr, rows_cols*sizeof(double));
        }

        if(origen.segment_ptr){
            segment = vtkSmartPointer<vtkImageData>::New();
            segment->SetExtent(0, cols-1, 0, rows-1, 0, 0);
            segment->AllocateScalars( VTK_DOUBLE, 1);
            segment->SetOrigin(0.0, 0.0, 0.0);
            segment->SetSpacing(1.0, 1.0, 1.0);
            segment_ptr = static_cast<double*>(segment->GetScalarPointer(0, 0, 0));

            memcpy(segment_ptr, origen.segment_ptr, rows_cols*sizeof(double));
        }

        if(origen.skl_ptr){
            skeleton = vtkSmartPointer<vtkImageData>::New();
            skeleton->SetExtent(0, cols+1, 0, rows+1, 0, 0);
            skeleton->AllocateScalars( VTK_DOUBLE, 1);
            skeleton->SetOrigin(0.0, 0.0, 0.0);
            skeleton->SetSpacing(1.0, 1.0, 1.0);
            skl_ptr = static_cast<double*>(skeleton->GetScalarPointer(0, 0, 0));

            memcpy(skl_ptr, origen.skl_ptr, rows_cols*sizeof(double));
        }

        if(origen.map_ptr){
            mapa_dist = vtkSmartPointer<vtkImageData>::New();
            mapa_dist->SetExtent(0, cols-1, 0, rows-1, 0, 0);
            mapa_dist->AllocateScalars( VTK_DOUBLE, 1);
            mapa_dist->SetOrigin(0.0, 0.0, 0.0);
            mapa_dist->SetSpacing(1.0, 1.0, 1.0);
            map_ptr = static_cast<double*>(mapa_dist->GetScalarPointer(0, 0, 0));

            memcpy(map_ptr, origen.map_ptr, rows_cols*sizeof(double));
        }

        if(origen.borders_ptr){
            borders = vtkSmartPointer<vtkImageData>::New();
            borders->SetExtent(0, cols-1, 0, rows-1, 0, 0);
            borders->AllocateScalars( VTK_DOUBLE, 1);
            borders->SetOrigin(0.0, 0.0, 0.0);
            borders->SetSpacing(1.0, 1.0, 1.0);
            borders_ptr = static_cast<double*>(borders->GetScalarPointer(0, 0, 0));

            memcpy(borders_ptr, origen.borders_ptr, rows_cols*sizeof(double));
        }
    }

	return *this;
}





// M E T O D O S      P R I V A D O S
/*  Metodo: setRuta
    Funcion: Copia al argumento a una variable local.
*/
char* IMGVTK::setRuta( const char *ruta_input ){
    const int l_src = strlen( ruta_input ) + 1;
    char *mi_ruta = new char [l_src];
    strcpy( mi_ruta, ruta_input );
    return mi_ruta;
}


// C L A S E: IMGVTK  ------------------------------------------------------------------------------------------------- ^
