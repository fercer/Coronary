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




//****************************************************************************************************************************************
//                              H E R R A M I E N T A S     D E     P R O C E S A M I E N T O       D E     I M A G E N E S
//****************************************************************************************************************************************
/*  Metodo: mapaDistancias
    Funcion: Obtiene el mapa de distancias de los pixeles a los bordes.
*/
void IMGVTK::mapaDistancias(){
    float INF = 1e9;

    float *map_tmp = new float [rens_cols];

    for( int xy = 0 ; xy < rens_cols; xy++ ){
        map_tmp[xy] = (base_ptr[xy] == 0) ? 0.0 : INF;
    }

    float *f = new float[ rens > cols ? rens : cols  ];

    // transform along columns
    for (int x = 0; x < cols; x++){
        for (int y = 0; y < rens; y++){
            f[y] = *(map_tmp + y*cols + x);
        }

        float *dh = new float[ rens ];
        int *vh= new int[ rens ];
        float *zh = new float[ rens+1 ];
        int k = 0;
        vh[0] = 0;
        zh[0] = -INF;
        zh[1] = +INF;

        for (int q = 1; q < rens; q++){
            float s  = ((f[q]+(q*q))-(f[vh[k]]+(vh[k]*vh[k])))/(2*q-2*vh[k]);

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
        for (int y = 0; y < rens; y++){
            while (zh[k+1] < y){
                k++;
            }
            *(map_tmp + y*cols + x) = ((y-vh[k])*(y-vh[k])) + f[vh[k]];
        }

        delete [] dh;delete [] vh;delete [] zh;
    }

    //normalizarlas a [0,255] y mostrar base
    float mini = INF, maxi = -INF;

    // transform along rows
    for (int y = 0; y < rens; y++){
        for (int x = 0; x < cols; x++){
            f[x] = map_tmp[y*cols+ x];
        }

        float *dw = new float[cols];
        int *vw = new int[cols];
        float *zw = new float[cols+1];
        int k = 0;
        vw[0] = 0;
        zw[0] = -INF;
        zw[1] = +INF;

        for (int q = 1; q < cols; q++){
            float s  = ((f[q]+(q*q))-(f[vw[k]]+(vw[k]*vw[k])))/(2*q-2*vw[k]);
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
            map_tmp[y*cols + x] = sqrt( ((x-vw[k])*(x-vw[k])) + f[vw[k]] );

            if( map_tmp[y*cols + x] > maxi){
                maxi = map_tmp[y*cols + x];
            }
            if( map_tmp[y*cols + x] < mini ){
                mini = map_tmp[y*cols + x];
            }
        }

        delete [] dw;
        delete [] vw;
        delete [] zw;
    }

    delete [] f;

    if( !map_ptr ){
        mapa_dist = vtkSmartPointer<vtkImageData>::New();
    }

    mapa_dist->SetExtent(0, cols-1, 0, rens-1, 0, 0);
    mapa_dist->AllocateScalars( VTK_UNSIGNED_CHAR, 1);
    mapa_dist->SetOrigin(0.0, 0.0, 0.0);
    mapa_dist->SetSpacing(1.0, 1.0, 1.0);

    map_ptr = static_cast<unsigned char *>(mapa_dist->GetScalarPointer(0,0,0));

    for(int xy=0; xy < rens_cols; xy++){
        map_ptr[xy] = (unsigned char)((255.0/maxi-mini)*( map_tmp[xy] - mini));
    }

    delete [] map_tmp;
}



/*  Metodo: regionFilling9
    Funcion: Rellena espacios vacios dentro del conjunto de pixeles resaltado.
*/
bool IMGVTK::regionFilling9( const unsigned char *ptr, const int x, const int y, const int mis_cols, const int mis_rens){
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

    DEB_MSG("n_hits 9 [" << x << "," << y << "]:" << n_hits);

    return (n_hits == 36);
}



/*  Metodo: regionFilling7
    Funcion: Rellena espacios vacios dentro del conjunto de pixeles resaltado.
*/
bool IMGVTK::regionFilling7( const unsigned char *ptr, const int x, const int y, const int mis_cols, const int mis_rens){

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

    DEB_MSG("n_hits 7 [" << x << "," << y << "]:" << n_hits);

    return (n_hits == 28);
}



/*  Metodo: regionFilling5
    Funcion: Rellena espacios vacios dentro del conjunto de pixeles resaltado.
*/
bool IMGVTK::regionFilling5( const unsigned char *ptr, const int x, const int y, const int mis_cols, const int mis_rens){
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

    DEB_MSG("n_hits 5 [" << x << "," << y << "]:" << n_hits);

    return (n_hits == 20);
}



/*  Metodo: regionFilling3
    Funcion: Rellena espacios vacios dentro del conjunto de pixeles resaltado.
*/
bool IMGVTK::regionFilling3( const unsigned char *ptr, const int x, const int y, const int mis_cols, const int mis_rens ){
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

    DEB_MSG("n_hits 3 [" << x << "," << y << "]:" << n_hits);

    return (n_hits == 12);
}



/*  Metodo: conexo
    Funcion: Metodo recursivo (dinamico) para encontrar los conjuntos conexos utilizando la conectividad 8.
*/
void IMGVTK::conexo(const unsigned char *ptr, const int x, const int y, int *conjuntos, unsigned int* n_etiquetados, bool* visitados, const int num_etiquetas, const int mis_cols, const int mis_rens){
    visitados[x + y*mis_cols] = true;

    // Si el pixel es parte del fondo, no se hace nada:
    if( ptr[x + y*mis_cols] < 255){
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
unsigned int* IMGVTK::conjuntosConexosDinamico(const unsigned char *ptr, int *conjuntos, const int mis_cols, const int mis_rens){
    int num_etiquetas = 0;

    const int mis_rens_cols = mis_cols*mis_rens;

    bool *visitados = new bool [mis_rens_cols];
    memset(visitados, 0, mis_rens_cols*sizeof(unsigned char));
DEB_MSG("visitados listo !");
    unsigned int *tmp = new unsigned int [mis_rens_cols];
    memset(tmp, 0, sizeof(unsigned int)*mis_rens_cols);
DEB_MSG("tmp listo !");

    memset(conjuntos, -1, sizeof(int) * mis_rens_cols);
DEB_MSG("conjuntos listo !");


DEB_MSG("pixeles: " << mis_rens_cols);
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




/*  Metodo: conjuntosConexos
    Funcion: Obtiene los conjuntos conexos de la imagen.
*/
unsigned int* IMGVTK::conjuntosConexos(unsigned char *ptr, int *conjuntos, const int mis_cols, const int mis_rens){


    const int mis_rens_cols = mis_cols*mis_rens;

    memset(conjuntos, -1, sizeof(int) * mis_rens_cols);

    //// Primer escaneo para determinar las etiquetas:
    // Verificar el primer pixel:
    if( ptr[0] == 255 ){
        conjuntos[0] = 0;
    }

    int eq, ep, ex;
    std::vector<int> etiquetas;
    etiquetas.push_back( 0 );
    int num_etiquetas = 0;

    // Verificar primer fila:
    for( int x = 1; x < mis_cols; x++){
        if( ptr[x] == 255 ){
            eq = conjuntos[x-1];
            ex = eq;
            if( eq == -1 ){
                num_etiquetas++;
                etiquetas.push_back(num_etiquetas);
                ex = num_etiquetas;
            }
            conjuntos[x] = ex;
        }
    }

    // Verifiar primer columna:
    for( int y = 1; y < mis_rens; y++){
        if( ptr[y*mis_cols] == 255 ){
            ep = conjuntos[(y-1)*mis_cols];
            ex = ep;
            if( ep == -1){
                num_etiquetas++;
                etiquetas.push_back( num_etiquetas );
                ex = num_etiquetas;
            }
            conjuntos[y*mis_cols] = ex;
        }
    }

    // Verificar el resto de la imagen:
    for( int y = 1; y < mis_rens; y++){
        for( int x = 1; x < mis_cols; x++){
            if( ptr[x + y*mis_cols] == 255 ){
                ep = conjuntos[x + (y-1)*mis_cols];
                eq = conjuntos[(x-1) + y*mis_cols];
                ex = ep;
                if( (ep == -1) && (eq == -1) ){
                    num_etiquetas++;
                    etiquetas.push_back( num_etiquetas );
                    ex = num_etiquetas;
                }else if( (ep != -1) && (eq != -1) && (etiquetas[ep] != etiquetas[eq]) ){
                    // Las etiquetas ep y eq son equivalentes, y por lo tanto se igualan:
                    const int k_eq = etiquetas[eq];
                    const int k_ep = etiquetas[ep];
                    for( std::vector<int>::iterator k = etiquetas.begin(); k < etiquetas.end(); k++){
                        if( *k == k_eq ){
                            *k = k_ep;
                        }
                    }

                    ex = ep;
                }else if( eq != -1){
                    ex = eq;
                }
                conjuntos[x+y*mis_cols] = ex;
            }
        }
    }
    DEB_MSG( "Numero de etiquetas puestas: " << etiquetas.size());


    // Acomodar las etiquetas:
    std::sort(etiquetas.begin(), etiquetas.end());
    ex = -1;
    num_etiquetas = 0;

    std::vector<int> diferentes;
    for( std::vector<int>::iterator k = etiquetas.begin(); k < etiquetas.end(); k++){
        if( *k != ex ){
            num_etiquetas++;
            ex = *k;
            diferentes.push_back(ex);
        }
    }

    DEB_MSG("Numero de etiquetas diferentes: " << (diferentes.size()));
    unsigned int *n_etiquetados = new unsigned int [num_etiquetas+1];
    memset( n_etiquetados, 0, (num_etiquetas+1)*sizeof(unsigned int) );
    n_etiquetados[0] = num_etiquetas;
    num_etiquetas = -1;
    ex = -1;
    for( std::vector<int>::iterator k = diferentes.begin(); k < diferentes.end(); k++){
        if( *k != ex ){
            ex = *k;
            num_etiquetas++;
        }
        for( int xy = 0; xy < mis_rens_cols; xy++){
            if( (conjuntos[xy]>=0) && ( etiquetas[conjuntos[xy]] == *k ) ){
                conjuntos[xy] = num_etiquetas;
                n_etiquetados[num_etiquetas+1]++;
            }
        }
    }

    return n_etiquetados;
}




/*  Metodo:  sklMask
    Funcion: Mascara usada para la extraccion del esqueleto.
*/
inline unsigned char IMGVTK::sklMask( const unsigned char *skl_ptr, const int x, const int y, const int mis_cols, const int mis_rens ){
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
void IMGVTK::lengthFilter(unsigned char *ptr, const int min_length, const int mis_cols, const int mis_rens ){

    const int mis_rens_cols = mis_cols*mis_rens;

    int *mis_conjuntos = new int [mis_rens_cols];

    DEB_MSG("cols: " << mis_cols << ", " << mis_rens);
    unsigned int *mis_n_etiquetados = conjuntosConexos(ptr, mis_conjuntos, mis_cols, mis_rens);
    DEB_MSG("n etiquetas: " << mis_n_etiquetados[0] );

    for( int xy = 0; xy < mis_rens_cols; xy++){
        if( (mis_conjuntos[ xy ] >= 0) && (mis_n_etiquetados[ mis_conjuntos[ xy ]+1 ] < min_length) ){
            ptr[ xy ] = 0;
        }
    }

    delete [] mis_conjuntos;
    delete [] mis_n_etiquetados;
}



/*  Metodo: regionFill
    Funcion: Rellena vacios dentro del cuerpo de la arteria segmentada.
*/
void IMGVTK::regionFill( unsigned char *ptr, const int mis_cols, const int mis_rens  ){
    // -------- Mascara 9x9
    for( int y = 0; y < mis_rens; y++){
        for( int x = 0; x < mis_cols; x++){
            if( ptr[x+y*mis_cols] == 0 ){
                ptr[x+y*mis_cols] = (regionFilling9(ptr, x, y, mis_cols, mis_rens) ? 255 : 0);
            }
        }
    }
    // -------- Mascara 7x7
    for( int y = 0; y < mis_rens; y++){
        for( int x = 0; x < mis_cols; x++){
            if( ptr[x+y*mis_cols] == 0 ){
                ptr[x+y*mis_cols] = (regionFilling7(ptr, x, y, mis_cols, mis_rens) ? 255 : 0);
            }
        }
    }
    // -------- Mascara 5x5
    for( int y = 0; y < mis_rens; y++){
        for( int x = 0; x < mis_cols; x++){
            if( ptr[x+y*mis_cols] == 0 ){
                ptr[x+y*mis_cols] = (regionFilling5(ptr, x, y, mis_cols, mis_rens) ? 255 : 0);
            }
        }
    }
    // -------- Mascara 3x3
    for( int y = 0; y < mis_rens; y++){
        for( int x = 0; x < mis_cols; x++){
            if( ptr[x+y*mis_cols] == 0 ){
                ptr[x+y*mis_cols] = (regionFilling3(ptr, x, y, mis_cols, mis_rens) ? 255 : 0);
            }
        }
    }
}





/*  Metodo:  dilMask
    Funcion: Mascara usada para dilatar unaimagen.
*/
inline unsigned char IMGVTK::dilMask( const unsigned char *mask_dil, const int x, const int y, const int mis_cols, const int mis_rens){
    return (mask_dil[( x ) + (y-1)*mis_cols] > 0) +
           (mask_dil[(x+1) + ( y )*mis_cols] > 0) +
           (mask_dil[( x ) + (y+1)*mis_cols] > 0) +
           (mask_dil[(x-1) + ( y )*mis_cols] > 0);
}



/*  Metodo:  erosionMask
    Funcion: Mascara para erosion usando un disco de radio 5
*/
inline unsigned char IMGVTK::erosionMask( const unsigned char *ptr_tmp, const int x, const int y, const int mis_cols, const int mis_rens ){
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
void IMGVTK::erosionar( unsigned char *ptr , const int mis_cols, const int mis_rens){
    unsigned char *ptr_tmp = new unsigned char [(mis_rens+8)*(mis_cols+8)];
    memset(ptr_tmp, 255, (mis_rens+8)*(mis_cols+8)*sizeof(unsigned char));

    for(int y = 0; y < mis_rens; y++){
        memcpy( ptr_tmp + (y+4)*(mis_cols+8) + 4, ptr + y*mis_cols, mis_cols*sizeof(unsigned char));
    }

    for( int y = 0; y < mis_rens; y++){
        for( int x = 0; x < mis_cols; x++){
            if(ptr_tmp[(x+4)+(y+4)*(mis_cols+8)] > 0){
                const unsigned char resp = erosionMask(ptr_tmp, x+4, y+4, mis_cols+8, mis_rens);
                if( (0 < resp) && (resp < 68) ){
                    ptr[x+y*mis_cols] = 0;
                }
            }
        }
    }

    delete [] ptr_tmp;
}



/*  Metodo: maskFOV
    Funcion: Obtiene la mascara del contorno
*/
void IMGVTK::maskFOV( unsigned char * img_tmp, unsigned char *mask_tmp, const int mis_cols, const int mis_rens){

    const int mis_rens_cols = mis_cols * mis_rens;

    // Se umbraliza al 0.1 la imagen original
    for( int xy = 0; xy < mis_rens_cols; xy++){
        mask_tmp[xy] = ((double)img_tmp[xy] / 255.0 < 0.1) ? 0 : 255;
    }

    // Se eliminan los conjuntos pequeños
    //lengthFilter(mask_tmp, 1000, mis_cols, mis_rens);

    // Se erosiona la mascara:
    erosionar(mask_tmp, mis_cols, mis_rens);

    // Se eliminan los conjuntos grandes que no esten en las esquinas:
    // Se extraen las etiquetas de los conjuntos que se conectan a las esquinas:
    for( int xy = 0; xy < mis_rens_cols; xy++){
        mask_tmp[xy] = 255*!mask_tmp[xy];
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
            mask_tmp[xy] = 0;
        }else{
            mask_tmp[xy] = 255;
        }
    }

    delete [] mis_conjuntos;
}



/*  Metodo: fillMask
    Funcion: Se llenan los espacios de la mascara con la media de la imagen circundante.
*/
void IMGVTK::fillMask( unsigned char *img_tmp, unsigned char *mask_tmp, const int mis_cols, const int mis_rens ){

    PIX_PAR par_tmp;
    par_tmp.pix_tipo = PIX_CROSS;

    unsigned char *mask_dil = new unsigned char [(mis_rens+2)*(mis_cols+2)];
    memset(mask_dil, 255, (mis_rens+2)*(mis_cols+2)*sizeof( unsigned char ));
    for( int y = 0; y < mis_rens; y++){
        memcpy( mask_dil + (y+1)*(mis_cols+2) + 1, mask_tmp + y*mis_cols, mis_cols*sizeof(unsigned char));
    }

    int iter = 0;
    while( iter < 50 ){
        // Obtener el borde de la mascara realizando una dilatacion:
        std::vector< PIX_PAR > borde;

        // Si es la primer pasada, se verifican todos los pixeles para encontrar el borde
        for( int y = 0; y < mis_rens; y++){
            for( int x = 0; x < mis_cols; x++){
                if( mask_dil[(x+1) + (y+1)*(mis_cols+2)] < 255 ){
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

            img_tmp[x_act + y_act*mis_cols] = (unsigned char)( suma / n_vecinos );
            mask_dil[(x_act+1) + (y_act+1)*(mis_cols+2)] = 255;
        }
        iter++;
    }
    delete [] mask_dil;
}



/*  Metodo: extraerCaract
    Funcion: Extrae los pixeles caracteristicos (end y branch points) a partir del esqueleot de la imagen.
*/
void IMGVTK:: extraerCaract(){
    std::vector<PIX_PAR> pares_caract;
    PIX_PAR temp;

    const unsigned char tabla[] = {
        0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        1, 0, 0, 0, 0, 2, 0, 0, 0, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        1, 0, 0, 0, 0, 2, 0, 0, 0, 2, 2, 0, 0, 0, 0, 0, 0, 0, 2, 0, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 2, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 2, 0, 2, 0, 0, 0, 2, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    // Buscar "End points", "Branch points" y "Cross points"
    for( int y = 1; y <= rens; y++){
        for( int x = 1; x <= cols; x++){
            if( skl_ptr[x+y*(cols+1)] > 0 ){
                const unsigned char resp = tabla[ sklMask( skl_ptr, x, y, cols+1, rens) ];
                switch( resp ){
                    case 1:{ /* END point*/
                        temp.x = x;
                        temp.y = y;
                        temp.pix_tipo = PIX_END;
                        pares_caract.push_back( temp );
                        break;
                    }
                    case 2:{ /* BRANCH point */
                        temp.x = x;
                        temp.y = y;
                        temp.pix_tipo = PIX_BRANCH;
                        pares_caract.push_back( temp );
                        break;
                    }
                    case 3:{ /* CROSS point */
                        temp.x = x;
                        temp.y = y;
                        temp.pix_tipo = PIX_CROSS;
                        pares_caract.push_back( temp );
                        break;
                    }
                }
            }
        }
    }

    n_caracts = pares_caract.size();
    if( n_caracts > 0 ){
        pix_caract = new PIX_PAR [ n_caracts ];
        for( int c = 0; c < n_caracts; c++){
            pix_caract[c].x = pares_caract[c].x;
            pix_caract[c].y = pares_caract[c].y;
            pix_caract[c].pix_tipo = pares_caract[c].pix_tipo;
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
void IMGVTK::skeletonization(){

    if( !skl_ptr ){
        skeleton = vtkSmartPointer<vtkImageData>::New();
    }

    skeleton->SetExtent(0, cols, 0, rens, 0, 0);
    skeleton->AllocateScalars( VTK_UNSIGNED_CHAR, 1);
    skeleton->SetOrigin(0.0, 0.0, 0.0);
    skeleton->SetSpacing(1.0, 1.0, 1.0);

    skl_ptr = static_cast<unsigned char *>(skeleton->GetScalarPointer(0,0,0));
    memset(skl_ptr, 0, (rens+1)*(cols+1)*sizeof(unsigned char));

    for( int y = 0; y < rens; y++ ){
        memcpy( skl_ptr + (y+1)*(cols+1)+1, base_ptr + y*cols, cols*sizeof(unsigned char) );
    }

    unsigned char *skl_mark = new unsigned char [(rens+1)*(cols+1)];
    memcpy(skl_mark, skl_ptr, (rens+1)*(cols+1)*sizeof(unsigned char));

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
        for(int y = 1; y <= rens; y++){
            for(int x = 1; x <= cols; x++){
                if(skl_ptr[x+y*(cols+1)] > 0){
                    const unsigned char resp = tabla[ sklMask(skl_ptr, x, y, cols+1, rens) ];
                    if(resp == 1 || resp == 3 ){
                        skl_mark[x+y*(cols+1)] = 0;
                        n_borrado++;
                    }
                }
            }
        }

        memcpy(skl_ptr, skl_mark, (rens+1)*(cols+1)*sizeof(unsigned char));

        // Segundo paso:
        for(int y = 1; y <= rens; y++){
            for(int x = 1; x <= cols; x++){
                if(skl_ptr[x+y*(cols+1)] > 0){
                    unsigned char resp = tabla[ sklMask(skl_ptr, x, y, cols, rens) ];
                    if(resp == 2 || resp == 3 ){
                        skl_mark[x+y*(cols+1)] = 0;
                        n_borrado++;
                    }
                }
            }
        }

        memcpy(skl_ptr, skl_mark, (rens+1)*(cols+1)*sizeof(unsigned char));

        DEB_MSG("Borrados: " << n_borrado);

    }while( n_borrado > 0 );

//    for( int y = 0; y < rens; y++ ){
//        memcpy( base_ptr + y*cols, skl_ptr + (y+1)*(cols+1)+1, cols*sizeof(unsigned char) );
//    }

    delete [] skl_mark;

    extraerCaract();
}




/*  Metodo: definirMask
    Funcion: Define una mascara para normalizar los pixeles de la imagen.
*/
void IMGVTK::definirMask( vtkSmartPointer<vtkImageData> img_src, vtkSmartPointer<vtkImageData> mask_src ){

    unsigned char *img_tmp = static_cast<unsigned char *>(img_src->GetScalarPointer(0,0,0));

    int dims[3];

    img_src->GetDimensions( dims );

    const int mis_cols = dims[0];
    const int mis_rens = dims[1];

    DEB_MSG("img cols: " << mis_cols << ", img rens" << mis_rens);

    mask_src->SetExtent(0, mis_cols-1, 0, mis_rens-1, 0, 0);
    mask_src->AllocateScalars( VTK_UNSIGNED_CHAR, 1);
    mask_src->SetOrigin(0.0, 0.0, 0.0);
    mask_src->SetSpacing(1.0, 1.0, 1.0);

    unsigned char *mask_tmp = static_cast<unsigned char *>(mask_src->GetScalarPointer(0,0,0));
    maskFOV( img_tmp, mask_tmp, mis_cols, mis_rens );
    fillMask( img_tmp, mask_tmp, mis_cols, mis_rens );
}








/*  Metodo: umbralizar
    Funcion: Utiliza el metodo de Otsu para umbralizar la imagen y separar el fondo y el primer plano de la imagen.
*/
void IMGVTK::umbralizar(){

    const int clases = 256;

    unsigned int *histograma = new unsigned int [clases];
    memset(histograma, 0, clases*sizeof(unsigned int));

    double max =-1e100;
    double min = 1e100;

    for( int xy = 0; xy < rens_cols; xy++){
        if(base_ptr[xy] < min){
            min = base_ptr[xy];
        }
        if(base_ptr[xy] > max){
            max = base_ptr[xy];
        }
    }

    // Obtener el histograma de frecuencias a cada nivel de grises de la imagen:
    for( int xy = 0; xy < rens_cols; xy ++){
        histograma[ (int)(clases*(base_ptr[xy]-min) / max) ]++;
    }


    // Obtener las frecuencias acumuladas:
    double suma = 0.0;
    for( int k = 1; k < clases; k++){
        suma += (double)(histograma[k] * k);
    }

    double suma_back = 0;
    double media_fore, media_back;
    double peso_fore, peso_back = 0.0;
    double varianza_entre, max_varianza_entre = -1.0;
    unsigned char umbral;

    for( int k = 0; k < clases; k++){
        // Calcular el peso del back y fore-ground:
        peso_back += (double)histograma[k];
        peso_fore = (double)rens_cols - peso_back;

        // Calcular la media del back y fore-ground:
        suma_back += (double)(histograma[k]*k);
        media_back = suma_back / peso_back;
        media_fore = (suma - suma_back) / peso_fore;

        // Calcular la varianza entre el fore y back ground:
        varianza_entre = peso_back*peso_fore*(media_back - media_fore)*(media_back - media_fore);

        if( varianza_entre > max_varianza_entre ){
            max_varianza_entre = varianza_entre;
            umbral = (unsigned char)k;
        }
    }

    // Se umbraliza la imagen:
    for( int xy = 0; xy < rens_cols; xy++){
        base_ptr[xy] = (base_ptr[xy] >= umbral) ? 255 : 0;
    }

    delete [] histograma;
}



/*  Metodo: Cargar
    Funcion: Cargar desde un archivo DICOM/JPEG/PNG/BMP la imagen a formato VTK.
*/
void IMGVTK::Cargar(const char *ruta_origen, vtkSmartPointer<vtkImageData> img_src, vtkSmartPointer<vtkImageData> mask_src, const int nivel, const bool enmascarar){

    const int ruta_l = strlen(ruta_origen);
DEB_MSG("Extension del archivo de entrada: " << (ruta_origen + ruta_l - 3));
    esDICOM = true;
    esDICOM *= strcmp(ruta_origen + ruta_l - 3, "png");
    esDICOM *= strcmp(ruta_origen + ruta_l - 3, "jpg");
    esDICOM *= strcmp(ruta_origen + ruta_l - 4, "jpeg");
    esDICOM *= strcmp(ruta_origen + ruta_l - 3, "bmp");

    if( esDICOM ){
        gdcm::ImageReader DICOMreader;
        DICOMreader.SetFileName( ruta_origen );

/*
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
////---------- Extraer DDP (Distancia del detector al paciente): ---------------------------------------------------------------------
        DDP = SID - SOD;
////---------- Extraer pixY\pixX (Distancia entre el centro de los pixeles en el eje Y\X): -----------------------------------------
        {
            const gdcm::DataElement &de = ds.GetDataElement( gdcm::Tag(0x18, 0x1164) );
            const gdcm::ByteValue *bv = de.GetByteValue();
            if( bv ){
                std::string strm(bv->GetPointer(), bv->GetLength());
                char *pixXYstr = new char [bv->GetLength()];
                memcpy(pixXYstr, strm.c_str(), bv->GetLength() * sizeof(char ));
                char *tmp = strchr(pixXYstr,'\\');
DEB_MSG("original: " << pixXYstr);
                pixXYstr[ tmp - pixXYstr ] = '\0';
DEB_MSG("trimmed: Y = " << pixXYstr << ", X = " << (tmp+1));

                pixY = atof(pixXYstr);
                pixX = atof(tmp+1);
            }
        }
////---------- Extraer LAO/RAO (Angulo del detector en direccion izquierda(-) a derecha(+)): -----------------------------------------
        {
            const gdcm::DataElement &de = ds.GetDataElement( gdcm::Tag(0x18, 0x1510) );
            const gdcm::ByteValue *bv = de.GetByteValue();
            if( bv ){
                std::string strm(bv->GetPointer(), bv->GetLength());
                LAORAO = atof( strm.c_str() );
            }
        }
////---------- Extraer CRA/CAU (Angulo del detector en direccion cranial(-) a caudal(+)): -----------------------------------------
        {
            const gdcm::DataElement &de = ds.GetDataElement( gdcm::Tag(0x18, 0x1511) );
            const gdcm::ByteValue *bv = de.GetByteValue();
            if( bv ){
                std::string strm(bv->GetPointer(), bv->GetLength());
                CRACAU = atof( strm.c_str() );
            }
        }
*/

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
        img_src->AllocateScalars( VTK_UNSIGNED_CHAR, 1);
        img_src->SetOrigin(0.0, 0.0, 0.0);
        img_src->SetSpacing(1.0, 1.0, 1.0);

        unsigned char *img_tmp = static_cast<unsigned char*>(img_src->GetScalarPointer(0,0,0));

        mask_src->SetExtent(0, mis_cols-1, 0, mis_rens-1, 0,0);
        mask_src->AllocateScalars( VTK_UNSIGNED_CHAR, 1);
        mask_src->SetOrigin(0.0, 0.0, 0.0);
        mask_src->SetSpacing(1.0, 1.0, 1.0);

        gdcm::PhotometricInterpretation scl_comps = gimage.GetPhotometricInterpretation();
        gdcm::PixelFormat pix_format = gimage.GetPixelFormat();

        switch(scl_comps){
            case gdcm::PhotometricInterpretation::RGB:{
DEB_MSG("Imagen DICOM en RGB...");
                    if( pix_format == gdcm::PixelFormat::UINT8 ){

                        img_tmp = static_cast<unsigned char*>(img_src->GetScalarPointer(0,0,0));
                        for( int xy = 0; xy < mis_rens_cols*3; xy+=3){
                            img_tmp[xy/3] = (unsigned char)buffer[xy + nivel*mis_rens_cols];
                        }
                    }else{
                        using namespace std;
                        cout << "============ ERROR AL CARGAR ARCHIVO DICOM ===========\nFormato de imagen RGB no soportado." << endl;
                    }
                    break;
                }
            case gdcm::PhotometricInterpretation::MONOCHROME2:{
DEB_MSG("Imagen DICOM en escala de grises...");
                    if( pix_format == gdcm::PixelFormat::UINT8 ){
DEB_MSG("Tipo UINT8");

                        img_tmp = static_cast<unsigned char*>(img_src->GetScalarPointer(0,0,0));
                        memcpy( img_tmp, buffer + nivel*mis_rens_cols, mis_rens_cols*sizeof(unsigned char));

                    }else if( pix_format == gdcm::PixelFormat::UINT16 ){
DEB_MSG("Tipo UINT16");
                        unsigned short *buffer16 = (unsigned short*)buffer;

                        img_tmp = static_cast<unsigned char*>(img_src->GetScalarPointer(0,0,0));
                        for( int xy = 0; xy < mis_rens_cols; xy+=3){
                            img_tmp[xy] = (unsigned char)buffer16[xy + nivel*mis_rens_cols*3] / 16;
                        }

                    }else{
                        using namespace std;
                        cout << "============ ERROR AL CARGAR ARCHIVO DICOM ===========\nFormato de imagen RGB no soportado." << endl;
                    }
                    break;
                }
        }
        delete [] buffer;


    }else{ ///------------------------------------------------------------------------------------------------



        // Leer la imagen en RBG o Escala de Grises:
        vtkSmartPointer<vtkImageReader2Factory> readerFactory = vtkSmartPointer<vtkImageReader2Factory>::New();
        vtkSmartPointer<vtkImageReader2> imgReader = readerFactory->CreateImageReader2( ruta_origen );

DEB_MSG("cargando: " << ruta_origen << " (PNG, BMP, JPEG)")

#ifndef NDEBUG
    FILE *fp = fopen(ruta_origen,"r");
    DEB_MSG("Existe el archivo?: " << (fp?"si":"no") );
    if(fp) fclose( fp );
#endif

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
        img_src->AllocateScalars( VTK_UNSIGNED_CHAR, 1);
        img_src->SetOrigin(0.0, 0.0, 0.0);
        img_src->SetSpacing(1.0, 1.0, 1.0);


        switch(scl_comps){
            // La imagen esta en escala de grises:
            case 1:{
DEB_MSG("La imagen esta en escala de grises")
                vtkSmartPointer<vtkImageExtractComponents> extractGreyFilter = vtkSmartPointer<vtkImageExtractComponents>::New();
                extractGreyFilter->SetInputConnection(imgReader->GetOutputPort());
                extractGreyFilter->SetComponents(0);
                extractGreyFilter->Update();

                unsigned char *gris = static_cast<unsigned char*>(extractGreyFilter->GetOutput()->GetScalarPointer(0,0,0));
                unsigned char *img_tmp = static_cast<unsigned char*>(img_src->GetScalarPointer(0,0,0));

                memcpy(img_tmp, gris, mis_rens_cols*sizeof(unsigned char));
                break;
            }
            // La imagen esta en RGB:
            case 3:{
DEB_MSG("La imagen esta en RGB")
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

                for( int i = 0 ; i < mis_rens; i++){
                    for( int j = 0; j < mis_cols; j++ ){
                        color = ((0.297)*ptR->GetScalarComponentAsDouble(j,i,0,0) +
                                 (0.589)*ptG->GetScalarComponentAsDouble(j,i,0,0) +
                                 (0.114)*ptB->GetScalarComponentAsDouble(j,i,0,0));

                        img_src->SetScalarComponentFromDouble(j,i,0,0,color);
                    }
                }
                break;
            }
            // La imagen esta en RGBA:
            case 4:{
DEB_MSG("La Imagen esta en RGB{A}")
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

                for( int i = 0 ; i < mis_rens; i++){
                    for( int j = 0; j < mis_cols; j++ ){
                        color = ((0.297)*ptR->GetScalarComponentAsDouble(j,i,0,0) +
                                 (0.589)*ptG->GetScalarComponentAsDouble(j,i,0,0) +
                                 (0.114)*ptB->GetScalarComponentAsDouble(j,i,0,0));

                        img_src->SetScalarComponentFromDouble(j,i,0,0,color);// * ptA->GetScalarComponentAsDouble(j,i,0,0));
                    }
                }
                break;
            }
        }

        imgReader->Delete();
    }

DEB_MSG("Imagen cargada");

    if(enmascarar){
        definirMask(img_src, mask_src);
DEB_MSG("Mascara obtenida");
    }
}



/*  Metodo: Cargar
    Funcion: Cargar la imagen a formato VTK desde varios archivos (Solo para imagenes JPEG, BMP o PNG).
*/
void IMGVTK::Cargar(vtkSmartPointer<vtkImageData> img_src, vtkSmartPointer<vtkImageData> mask_src, char **rutas, const int n_imgs, const bool enmascarar){
    vtkSmartPointer<vtkImageData> auxiliar = vtkSmartPointer<vtkImageData>::New();
    vtkSmartPointer<vtkImageData> mask_auxiliar = vtkSmartPointer<vtkImageData>::New();

    Cargar(rutas[0], auxiliar, mask_auxiliar, 0, enmascarar);

    unsigned char *aux_tmp = static_cast<unsigned char*>(auxiliar->GetScalarPointer(0, 0, 0));
    unsigned char *mskaux_tmp = static_cast<unsigned char*>(mask_auxiliar->GetScalarPointer(0, 0, 0));

    int dims[3];
    auxiliar->GetDimensions(dims);

    const int mis_cols = dims[0];
    const int mis_rens = dims[1];

    const int n_cols = mis_cols * n_imgs;

    // Extract data:
    img_src->SetExtent(0, n_cols-1, 0, mis_rens-1, 0, 0);
    img_src->AllocateScalars( VTK_UNSIGNED_CHAR, 1);
    img_src->SetOrigin(0.0, 0.0, 0.0);
    img_src->SetSpacing(1.0, 1.0, 1.0);

    mask_src->SetExtent(0, n_cols-1, 0, mis_rens-1, 0, 0);
    mask_src->AllocateScalars( VTK_UNSIGNED_CHAR, 1);
    mask_src->SetOrigin(0.0, 0.0, 0.0);
    mask_src->SetSpacing(1.0, 1.0, 1.0);

    unsigned char *base_ptr_tmp = static_cast<unsigned char*>(img_src->GetScalarPointer(0, 0, 0));
    unsigned char *mask_ptr_tmp = static_cast<unsigned char*>(mask_src->GetScalarPointer(0, 0, 0));
    for( int y = 0; y < mis_rens; y++){
        memcpy(base_ptr_tmp + y*n_cols, aux_tmp + y*mis_cols, mis_cols*sizeof(unsigned char));
        if(enmascarar){
            memcpy(mask_ptr_tmp + y*n_cols, mskaux_tmp + y*mis_cols, mis_cols*sizeof(unsigned char));
        }
    }

    for( int img_i = 1; img_i < n_imgs; img_i++){
        Cargar(rutas[img_i], auxiliar, mask_auxiliar, 0, enmascarar);
        aux_tmp = static_cast<unsigned char*>(auxiliar->GetScalarPointer(0, 0, 0));
        mskaux_tmp = static_cast<unsigned char*>(mask_auxiliar->GetScalarPointer(0, 0, 0));

        for( int y = 0; y < mis_rens; y++){
            memcpy(base_ptr_tmp + y*n_cols + img_i*mis_cols, aux_tmp+ y*mis_cols, mis_cols*sizeof(unsigned char));
            if(enmascarar){
                memcpy(mask_ptr_tmp + y*n_cols + img_i*mis_cols, mskaux_tmp + y*mis_cols, mis_cols*sizeof(unsigned char));
            }
        }
    }
}



/*  Metodo: Cargar
    Funcion: Cargar la imagen a formato VTK desde un archivo.
*/
void IMGVTK::Cargar(const char *ruta_origen, const bool enmascarar, const int nivel){
    base = vtkSmartPointer<vtkImageData>::New();
    mask = vtkSmartPointer<vtkImageData>::New();

    Cargar(ruta_origen, base, mask, nivel, enmascarar);

    base_ptr = static_cast<unsigned char*>(base->GetScalarPointer(0, 0, 0));
    mask_ptr = static_cast<unsigned char*>(mask->GetScalarPointer(0, 0, 0));

    int dims[3];
    base->GetDimensions( dims );

    cols = dims[0];
    rens = dims[1];
    rens_cols = rens*cols;
}



/*  Metodo: Cargar
    Funcion: Cargar la imagen a formato VTK desde varios archivos.
*/
void IMGVTK::Cargar(char **rutas , const int n_imgs, const bool enmascarar){
    base = vtkSmartPointer<vtkImageData>::New();
    mask = vtkSmartPointer<vtkImageData>::New();

    Cargar( base, mask, rutas, n_imgs, enmascarar );

    base_ptr = static_cast<unsigned char*>(base->GetScalarPointer(0, 0, 0));
    mask_ptr = static_cast<unsigned char*>(mask->GetScalarPointer(0, 0, 0));

    int dims[3];
    base->GetDimensions( dims );

    cols = dims[0];
    rens = dims[1];
    rens_cols = rens*cols;
}



/*  Metodo: Guardar
    Funcion: Guarda la imagen.
*/
void IMGVTK::Guardar(){
    switch(tipo_salida){
        case PGM:{
            break;
        }

        case PNG:{
            vtkSmartPointer<vtkPNGWriter> png_output = vtkSmartPointer<vtkPNGWriter>::New();
            png_output->SetFileName( ruta_salida );          
            png_output->SetInputData( base );
            png_output->Write();
            break;
        }
    }
}



/*  Metodo: Guardar
    Funcion: Guarda la imagen en la ruta especificada con la extension especificada.
*/
void IMGVTK::Guardar( const char *ruta, const TIPO_IMG tipo_salida_src ){
    if(ruta_salida){
        delete [] ruta_salida;
    }
    ruta_salida = setRuta(ruta);
    tipo_salida = tipo_salida_src;
    Guardar();
}



/* CONSTRUCTORES */
IMGVTK::IMGVTK(){
    ruta_salida = NULL;
    tipo_salida = PNG;

    cols = 0;
    rens = 0;

    map_ptr = NULL;
    base_ptr = NULL;
    skl_ptr = NULL;
    pix_caract = NULL;
    mask_ptr = NULL;

    // Defaults:
    SID = 1100.0;
    SOD = 400.0;
    DDP = SID - SOD;
    LAORAO = 20.0;
    CRACAU = 20.0;
    pixX = 0.308;
    pixY = 0.308;
}



IMGVTK::IMGVTK( const IMGVTK &original ){
    ruta_salida = NULL;
    tipo_salida = PNG;

    map_ptr = NULL;
    base_ptr = NULL;
    skl_ptr = NULL;
    pix_caract = NULL;
    mask_ptr = NULL;

    // Defaults:
    SID = original.SID;
    SOD = original.SOD;
    DDP = SID - SOD;
    LAORAO = original.LAORAO;
    CRACAU = original.CRACAU;
    pixX = original.pixX;
    pixY = original.pixY;

    if(original.base_ptr){

        if(original.ruta_salida){
            if(ruta_salida){
                delete [] ruta_salida;
            }
            ruta_salida = setRuta(original.ruta_salida);
        }

        tipo_salida = original.tipo_salida;

        cols = original.cols;
        rens = original.rens;

DEB_MSG("cols = " << cols);
DEB_MSG("rens = " << rens);

        base = vtkSmartPointer<vtkImageData>::New();
        base->SetExtent(0, cols-1, 0, rens-1, 0, 0);
        base->AllocateScalars( VTK_UNSIGNED_CHAR, 1);
        base->SetOrigin(0.0, 0.0, 0.0);
        base->SetSpacing(1.0, 1.0, 1.0);

        base_ptr = static_cast<unsigned char*>(base->GetScalarPointer(0, 0, 0));

        mask = vtkSmartPointer<vtkImageData>::New();
        mask->SetExtent(0, cols-1, 0, rens-1, 0, 0);
        mask->AllocateScalars( VTK_UNSIGNED_CHAR, 1);
        mask->SetOrigin(0.0, 0.0, 0.0);
        mask->SetSpacing(1.0, 1.0, 1.0);
        mask_ptr = static_cast<unsigned char*>(mask->GetScalarPointer(0, 0, 0));
        rens_cols = rens * cols;

        memcpy(base_ptr, original.base_ptr, rens_cols*sizeof(unsigned char));

        if(original.mask_ptr){
            memcpy(mask_ptr, original.mask_ptr, rens_cols*sizeof(unsigned char));
        }else{
            memset(mask_ptr, 0, rens_cols*sizeof(unsigned char));
        }
    }
}



IMGVTK::IMGVTK( char **rutas_origen, const int n_imgs, const bool enmascarar){
    ruta_salida = NULL;
    tipo_salida = PNG;

    map_ptr = NULL;
    base_ptr = NULL;
    skl_ptr = NULL;
    pix_caract = NULL;
    mask_ptr = NULL;

    // Defaults:
    SID = 1100.0;
    SOD = 400.0;
    DDP = SID - SOD;
    LAORAO = 20.0;
    CRACAU = 20.0;
    pixX = 0.308;
    pixY = 0.308;

    Cargar(rutas_origen, n_imgs, enmascarar );
}


IMGVTK::IMGVTK( const char *ruta_origen, const bool enmascarar, const int nivel){
    ruta_salida = NULL;
    tipo_salida = PNG;

    map_ptr = NULL;
    base_ptr = NULL;
    skl_ptr = NULL;
    pix_caract = NULL;
    mask_ptr = NULL;

    // Defaults:
    SID = 1100.0;
    SOD = 400.0;
    DDP = SID - SOD;
    LAORAO = 20.0;
    CRACAU = 20.0;
    pixX = 0.308;
    pixY = 0.308;

    Cargar(ruta_origen, enmascarar, nivel);
}


/* DESTRUCTOR */
IMGVTK::~IMGVTK(){
    if(ruta_salida){
        delete [] ruta_salida;
    }
    if(pix_caract){
        delete [] pix_caract;
    }
}




// O P E R A D O R E S  S O B R E C A R G A D O S
// El operador de copia extrae unicamente el contenido de la imagen original
IMGVTK& IMGVTK::operator= ( const IMGVTK &origen ){

    // Defaults:
    SID = origen.SID;
    SOD = origen.SOD;
    DDP = SID - SOD;
    LAORAO = origen.LAORAO;
    CRACAU = origen.CRACAU;
    pixX = origen.pixX;
    pixY = origen.pixY;

    if(origen.base_ptr){

        if(origen.ruta_salida){
            if(ruta_salida){
                delete [] ruta_salida;
            }
            ruta_salida = setRuta(origen.ruta_salida);
        }

        tipo_salida = origen.tipo_salida;

        cols = origen.cols;
        rens = origen.rens;

        base = vtkSmartPointer<vtkImageData>::New();
        base->SetExtent(0, cols-1, 0, rens-1, 0, 0);
        base->AllocateScalars( VTK_UNSIGNED_CHAR, 1);
        base->SetOrigin(0.0, 0.0, 0.0);
        base->SetSpacing(1.0, 1.0, 1.0);

        base_ptr = static_cast<unsigned char*>(base->GetScalarPointer(0, 0, 0));


        mask = vtkSmartPointer<vtkImageData>::New();
        mask->SetExtent(0, cols-1, 0, rens-1, 0, 0);
        mask->AllocateScalars( VTK_UNSIGNED_CHAR, 1);
        mask->SetOrigin(0.0, 0.0, 0.0);
        mask->SetSpacing(1.0, 1.0, 1.0);
        mask_ptr = static_cast<unsigned char*>(mask->GetScalarPointer(0, 0, 0));
        rens_cols = rens * cols;

        if(origen.base_ptr){
            memcpy(base_ptr, origen.base_ptr, rens_cols*sizeof(unsigned char));
        }else{
            memset(mask_ptr, 0, rens_cols*sizeof(unsigned char));
        }

        if(origen.mask_ptr){
            memcpy(mask_ptr, origen.mask_ptr, rens_cols*sizeof(unsigned char));
        }else{
            memset(mask_ptr, 0, rens_cols*sizeof(unsigned char));
        }
    }
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
