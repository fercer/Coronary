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
    std::cout << mensaje;
}



//****************************************************************************************************************************************
//                              H E R R A M I E N T A S     D E     P R O C E S A M I E N T O       D E     I M A G E N E S
//****************************************************************************************************************************************
/*  Metodo: mapaDistancias
    Funcion: Obtiene el mapa de distancias de los pixeles a los bordes.
*/
void IMGVTK::mapaDistancias(IMG_IDX img_idx) {
	IMGCONT *img_tmp = NULL;

	switch (img_idx) {
	case BASE:
		img_tmp = my_base;
		break;
	case GROUNDTRUTH:
		img_tmp = my_groundtruth;
		break;
	case MASK:
		img_tmp = my_mask;
		break;
	case SKELETON:
		img_tmp = my_skeleton;
		break;
	case SEGMENT:
		img_tmp = my_response;
		break;
	case THRESHOLD:
		img_tmp = my_segmented;
		break;
	}


	if (!my_distmap) {
		my_distmap = new IMGCONT(img_tmp->my_height, img_tmp->my_width);
	}

	for (int xy = 0; xy < (my_distmap->my_height * my_distmap->my_width); xy++) {
		*(my_distmap->my_img_data + xy) = (*(img_tmp->my_img_data + xy) < 1.0) ? 0.0 : MY_INF;
	}

	double *f = new double[my_distmap->my_height > my_distmap->my_width ?
		my_distmap->my_height :
		my_distmap->my_width];

	double *dh = new double[my_distmap->my_height];
	int *vh = new int[my_distmap->my_height];
	double *zh = new double[my_distmap->my_height + 1];

	// transform along columns
	for (int x = 0; x < cols; x++) {
		for (int y = 0; y < my_distmap->my_height; y++) {
			f[y] = *(my_distmap->my_img_data + y*my_distmap->my_width + x);
		}

		int k = 0;
		vh[0] = 0;
		zh[0] = -MY_INF;
		zh[1] = MY_INF;

		for (int q = 1; q < my_distmap->my_height; q++) {
			double s = ((f[q] + (q*q)) - (f[vh[k]] + (vh[k] * vh[k]))) / (2 * q - 2 * vh[k]);

			while (s <= zh[k]) {
				k--;
				s = ((f[q] + (q*q)) - (f[vh[k]] + (vh[k] * vh[k]))) / (2 * q - 2 * vh[k]);
			}

			k++;
			vh[k] = q;
			zh[k] = s;
			zh[k + 1] = MY_INF;
		}

		k = 0;
		for (int y = 0; y < my_distmap->my_height; y++) {
			while (zh[k + 1] < y) {
				k++;
			}
			*(my_distmap->my_img_data + y*my_distmap->my_width + x) = ((y - vh[k])*(y - vh[k])) + f[vh[k]];
		}

	}
	delete[] dh;
	delete[] vh;
	delete[] zh;

	double *dw = new double[my_distmap->my_width];
	int *vw = new int[my_distmap->my_width];
	double *zw = new double[my_distmap->my_width + 1];

	// transform along rows
	for (int y = 0; y < rows; y++) {
		for (int x = 0; x < my_distmap->my_width; x++) {
			f[x] = *(my_distmap->my_img_data + y*my_distmap->my_width + x);
		}
		int k = 0;
		vw[0] = 0;
		zw[0] = -MY_INF;
		zw[1] = +MY_INF;

		for (int q = 1; q < my_distmap->my_width; q++) {
			double s = ((f[q] + (q*q)) - (f[vw[k]] + (vw[k] * vw[k]))) / (2 * q - 2 * vw[k]);
			while (s <= zw[k]) {
				k--;
				s = ((f[q] + (q*q)) - (f[vw[k]] + (vw[k] * vw[k]))) / (2 * q - 2 * vw[k]);
			}
			k++;
			vw[k] = q;
			zw[k] = s;
			zw[k + 1] = +MY_INF;
		}

		k = 0;
		for (int x = 0; x < my_distmap->my_width; x++) {
			while (zw[k + 1] < x) {
				k++;
			}
			*(my_distmap->my_img_data + y*my_distmap->my_width + x) = sqrt(((x - vw[k])*(x - vw[k])) + f[vw[k]]);

			if (((int)*(my_distmap->my_img_data + y*my_distmap->my_width + x) + 1) > max_dist) {
				max_dist = (int)*(my_distmap->my_img_data + y*my_distmap->my_width + x) + 1;
			}
		}
	}

	delete[] dw;
	delete[] vw;
	delete[] zw;
	delete[] f;
}



/*  Metodo: detectarBorde
	Funcion: Detecta los bordes de la imagen en base a la transformada de la distancia.
*/
void IMGVTK::detectarBorde(IMG_IDX img_idx)
{

	if (!my_distmap) {
		mapaDistancias(img_idx);
	}

	IMGCONT *img_tmp = NULL;

	switch (img_idx) {
	case BASE:
		img_tmp = my_base;
		break;
	case GROUNDTRUTH:
		img_tmp = my_groundtruth;
		break;
	case MASK:
		img_tmp = my_mask;
		break;
	case SKELETON:
		img_tmp = my_skeleton;
		break;
	case SEGMENT:
		img_tmp = my_response;
		break;
	case THRESHOLD:
		img_tmp = my_segmented;
		break;
	}

	if (!my_boundaries) {
		my_boundaries = new IMGCONT(img_tmp->my_height, img_tmp->my_width);
	}

	for (int xy = 0; xy < (img_tmp->my_height * img_tmp->my_width); xy++) {
		if ((*(img_tmp->my_img_data + xy) > 0.5) && (*(my_distmap->my_img_data + xy) < 2.0)) {
			*(my_boundaries->my_img_data + xy) = 1.0;
		}
	}
}



/*  Metodo: regionFilling9
    Funcion: Rellena espacios vacios dentro del conjunto de pixeles resaltado.
*/
bool IMGVTK::regionFilling9(IMGCONT *img_src, const int x, const int y)
{
	int n_hits = 0;

	for (int m = 0; m < 9; m++) {
		// Arriba:
		if ((y - 4) > 0) {
			if (((x - 4 + m) > 0) && ((x - 4 + m) < img_src->my_width)) {
				n_hits += (*(img_src->my_img_data + (x - 4 + m) + (y - 4)*img_src->my_width) > 0);
			}
			else {
				n_hits++;
			}
		}
		else {
			n_hits++;
		}

		// Abajo:
		if ((y + 4) < img_src->my_height) {
			if (((x - 4 + m) > 0) && ((x - 4 + m) < img_src->my_width)) {
				n_hits += (*(img_src->my_img_data + (x - 4 + m) + (y + 4)*img_src->my_width) > 0);
			}
			else {
				n_hits++;
			}
		}
		else {
			n_hits++;
		}

		// Izquierda:
		if ((x - 4) > 0) {
			if (((y - 4 + m) > 0) && ((y - 4 + m) < img_src->my_height)) {
				n_hits += (*(img_src->my_img_data + (x - 4) + (y - 4 + m)*img_src->my_width) > 0);
			}
			else {
				n_hits++;
			}
		}
		else {
			n_hits++;
		}

		// Derecha:
		if ((x + 4) < img_src->my_width) {
			if (((y - 4 + m) > 0) && ((y - 4 + m) < img_src->my_height)) {
				n_hits += (*(img_src->my_img_data + (x + 4) + (y - 4 + m)*img_src->my_width) > 0);
			}
			else {
				n_hits++;
			}
		}
		else {
			n_hits++;
		}
	}
	return (n_hits == 36);
}



/*  Metodo: regionFilling7
    Funcion: Rellena espacios vacios dentro del conjunto de pixeles resaltado.
*/
bool IMGVTK::regionFilling7(IMGCONT *img_src, const int x, const int y) 
{

	int n_hits = 0;

	for (int m = 0; m < 7; m++) {
		// Arriba:
		if ((y - 3) > 0) {
			if (((x - 3 + m) > 0) && ((x - 3 + m) < img_src->my_width)) {
				n_hits += (*(img_src->my_img_data + (x - 3 + m) + (y - 3)*img_src->my_width) > 0);
			}
			else {
				n_hits++;
			}
		}
		else {
			n_hits++;
		}

		// Abajo:
		if ((y + 3) < img_src->my_height) {
			if (((x - 3 + m) > 0) && ((x - 3 + m) < img_src->my_width)) {
				n_hits += (*(img_src->my_img_data + (x - 3 + m) + (y + 3)*img_src->my_width) > 0);
			}
			else {
				n_hits++;
			}
		}
		else {
			n_hits++;
		}

		// Izquierda:
		if ((x - 3) > 0) {
			if (((y - 3 + m) > 0) && ((y - 3 + m) < img_src->my_height)) {
				n_hits += (*(img_src->my_img_data + (x - 3) + (y - 3 + m)*img_src->my_width) > 0);
			}
			else {
				n_hits++;
			}
		}
		else {
			n_hits++;
		}

		// Derecha:
		if ((x + 3) < img_src->my_width) {
			if (((y - 3 + m) > 0) && ((y - 3 + m) < img_src->my_height)) {
				n_hits += (*(img_src->my_img_data + (x + 3) + (y - 3 + m)*img_src->my_width) > 0);
			}
			else {
				n_hits++;
			}
		}
		else {
			n_hits++;
		}
	}

	return (n_hits == 28);
}



/*  Metodo: regionFilling5
    Funcion: Rellena espacios vacios dentro del conjunto de pixeles resaltado.
*/
bool IMGVTK::regionFilling5(IMGCONT *img_src, const int x, const int y)
{
	int n_hits = 0;

	for (int m = 0; m < 5; m++) {
		// Arriba:
		if ((y - 2) > 0) {
			if (((x - 2 + m) > 0) && ((x - 2 + m) < img_src->my_width)) {
				n_hits += (*(img_src->my_img_data + (x - 2 + m) + (y - 2)*img_src->my_width) > 0);
			}
			else {
				n_hits++;
			}
		}
		else {
			n_hits++;
		}

		// Abajo:
		if ((y + 2) < img_src->my_height) {
			if (((x - 2 + m) > 0) && ((x - 2 + m) < img_src->my_width)) {
				n_hits += (*(img_src->my_img_data + (x - 2 + m) + (y + 2)*img_src->my_width) > 0);
			}
			else {
				n_hits++;
			}
		}
		else {
			n_hits++;
		}

		// Izquierda:
		if ((x - 2) > 0) {
			if (((y - 2 + m) > 0) && ((y - 2 + m) < img_src->my_height)) {
				n_hits += (*(img_src->my_img_data + (x - 2) + (y - 2 + m)*img_src->my_width) > 0);
			}
			else {
				n_hits++;
			}
		}
		else {
			n_hits++;
		}

		// Derecha:
		if ((x + 2) < img_src->my_width) {
			if (((y - 2 + m) > 0) && ((y - 2 + m) < img_src->my_height)) {
				n_hits += (*(img_src->my_img_data + (x + 2) + (y - 2 + m)*img_src->my_width) > 0);
			}
			else {
				n_hits++;
			}
		}
		else {
			n_hits++;
		}
	}

	return (n_hits == 20);
}



/*  Metodo: regionFilling3
    Funcion: Rellena espacios vacios dentro del conjunto de pixeles resaltado.
*/
bool IMGVTK::regionFilling3(IMGCONT *img_src, const int x, const int y)
{
	int n_hits = 0;

	for (int m = 0; m < 3; m++) {
		// Arriba:
		if ((y - 1) > 0) {
			if (((x - 1 + m) > 0) && ((x - 1 + m) < img_src->my_width)) {
				n_hits += (*(img_src->my_img_data + (x - 1 + m) + (y - 1)*img_src->my_width) > 0);
			}
			else {
				n_hits++;
			}
		}
		else {
			n_hits++;
		}

		// Abajo:
		if ((y + 1) < img_src->my_height) {
			if (((x - 1 + m) > 0) && ((x - 1 + m) < img_src->my_width)) {
				n_hits += (*(img_src->my_img_data + (x - 1 + m) + (y + 1)*img_src->my_width) > 0);
			}
			else {
				n_hits++;
			}
		}
		else {
			n_hits++;
		}

		// Izquierda:
		if ((x - 1) > 0) {
			if (((y - 1 + m) > 0) && ((y - 1 + m) < img_src->my_height)) {
				n_hits += (*(img_src->my_img_data + (x - 1) + (y - 1 + m)*img_src->my_width) > 0);
			}
			else {
				n_hits++;
			}
		}
		else {
			n_hits++;
		}

		// Derecha:
		if ((x + 1) < img_src->my_width) {
			if (((y - 1 + m) > 0) && ((y - 1 + m) < img_src->my_height)) {
				n_hits += (*(img_src->my_img_data + (x + 1) + (y - 1 + m)*img_src->my_width) > 0);
			}
			else {
				n_hits++;
			}
		}
		else {
			n_hits++;
		}
	}

	return (n_hits == 12);
}



/*  Metodo: conexo
    Funcion: Metodo recursivo (dinamico) para encontrar los conjuntos conexos utilizando la conectividad 8.
*/
void IMGVTK::conexo(IMGCONT *img_src, const int x, const int y, int *conjuntos, unsigned int* n_etiquetados, bool* visitados, const int num_etiquetas)
{
	*(visitados + x + y*img_src->my_width) = true;

	// Si el pixel es parte del fondo, no se hace nada:
	if (*(img_src->my_img_data + x + y*img_src->my_width) < 1.0) {
		return;
	}
	else {
		*(conjuntos + x + y*img_src->my_width) = num_etiquetas;
		*(n_etiquetados + num_etiquetas) = *(n_etiquetados + num_etiquetas) + 1;
	}

	// De lo contrario, se revisan los pixeles vecinos:
	/// Superiores:
	if (y > 0) {
		if (x > 0) {
			/// Ariba - Izquierda
			if (!*(visitados + x - 1 + (y - 1)*img_src->my_width)) { // Si el pixel arriba-izquierda no ha sido visitado:
				conexo(img_src, x - 1, y - 1, conjuntos, n_etiquetados, visitados, num_etiquetas);
			}
		}

		/// Arriba
		if (!*(visitados + x + (y - 1)*img_src->my_width)) { // Si el pixel arriba no ha sido visitado:
			conexo(img_src, x, y - 1, conjuntos, n_etiquetados, visitados, num_etiquetas);
		}

		/// Ariba - Derecha
		if (x < (img_src->my_width - 1)) {
			if (!*(visitados + x + 1 + (y - 1)*img_src->my_width)) { // Si el pixel arriba-derecha no ha sido visitado:
				conexo(img_src, x + 1, y - 1, conjuntos, n_etiquetados, visitados, num_etiquetas);
			}
		}
	}

	/// Centrales:
	if (x > 0) {
		/// Centro - Izquierda
		if (!*(visitados + x - 1 + y*img_src->my_width)) { // Si el pixel centro-izquierda no ha sido visitado:
			conexo(img_src, x - 1, y, conjuntos, n_etiquetados, visitados, num_etiquetas);
		}
	}

	/// Centro - Derecha
	if (x < (img_src->my_width - 1)) {
		if (!*(visitados + x + 1 + y*img_src->my_width)) { // Si el pixel arriba-izquierda no ha sido visitado:
			conexo(img_src, x + 1, y, conjuntos, n_etiquetados, visitados, num_etiquetas);
		}
	}

	/// Inferiores:
	if (y < (img_src->my_width - 1)) {
		if (x > 0) {
			/// Abajo - Izquierda
			if (!*(visitados + x - 1 + (y + 1)*img_src->my_width)) { // Si el pixel arriba-izquierda no ha sido visitado:
				conexo(img_src, x - 1, y + 1, conjuntos, n_etiquetados, visitados, num_etiquetas);
			}
		}

		/// Abajo
		if (!*(visitados + x + (y + 1)*img_src->my_width)) { // Si el pixel arriba no ha sido visitado:
			conexo(img_src, x, y + 1, conjuntos, n_etiquetados, visitados, num_etiquetas);
		}

		/// Abajo - Derecha
		if ((x < (img_src->my_width - 1))) {
			if (!*(visitados + x + 1 + (y + 1)*img_src->my_width)) { // Si el pixel arriba-derecha no ha sido visitado:
				conexo(img_src, x + 1, y + 1, conjuntos, n_etiquetados, visitados, num_etiquetas);
			}
		}
	}
}



/*  Metodo: conjuntosConexos
    Funcion: Obtiene los conjuntos conexos de la imagen usando programacion dinamica (Llega a generar stackoverflow ... ).
*/
unsigned int* IMGVTK::conjuntosConexosDinamico(IMGCONT *img_src, int *conjuntos){
    int num_etiquetas = 0;

	const int mis_rows_cols = img_src->my_height * img_src->my_width;

    bool *visitados = new bool [mis_rows_cols];
    memset(visitados, 0, mis_rows_cols*sizeof(bool));

    unsigned int *tmp = new unsigned int [mis_rows_cols];
    memset(tmp, 0, sizeof(unsigned int)*mis_rows_cols);

    memset(conjuntos, -1, sizeof(int) * mis_rows_cols);

    for( int y = 0; y< img_src->my_height; y++){
        for( int x = 0; x < img_src->my_width; x++){
            /// Se mueve la etiqueta si el pixel es de un nuevo conjunto:
            if( !visitados[x + y*img_src->my_width]) {
                conexo( img_src, x, y, conjuntos, tmp, visitados, num_etiquetas);
                if( *(img_src->my_img_data + x + y*img_src->my_width) > 0 ){
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
inline void IMGVTK::ampliarConjunto( int *etiquetas, const int equiv_A, const int equiv_B, const int max_etiquetas){
    int base, base_ant;

    if( equiv_A > equiv_B ){
        base = equiv_A;
        base_ant = equiv_B;
    }else{
        base = equiv_B;
        base_ant = equiv_A;
    }

    for(int i = 0; i < max_etiquetas; i++){
        if(*(etiquetas + i) == -1){
            break;
        }else if( *(etiquetas + i) == base_ant ){
            *(etiquetas + i) = base;
        }
    }
}




/*  Metodo: conjuntosConexos
    Funcion: Obtiene los conjuntos conexos de la imagen.
*/
unsigned int* IMGVTK::conjuntosConexos(IMGCONT *img_src, int *conjuntos){

    static int veces = 0;

    const int mis_rows_cols = img_src->my_height*img_src->my_width;

    int max_etiquetas = mis_rows_cols;
    int *val_etiquetas = new int [max_etiquetas];
    memset( val_etiquetas, -1, max_etiquetas * sizeof(int));

    int *pix_etq = new int [(img_src->my_height +2)*(img_src->my_width +2)];
    memset( pix_etq, -1, (img_src->my_height +2)*(img_src->my_width +2) * sizeof(int));

    int n_conjunto = -1;

    for( int y = 1; y <= img_src->my_height; y++){
        for( int x = 1; x <= img_src->my_width; x++){
            if( *(img_src->my_img_data + (y-1)*img_src->my_width + (x-1)) > 0.0 ){

                const int NO = *(pix_etq + (y-1)*(img_src->my_width +2) + (x-1));
                const int N  = *(pix_etq + (y-1)*(img_src->my_width +2) +   x  );
                const int NE = *(pix_etq + (y-1)*(img_src->my_width +2) + (x+1));
                const int E  = *(pix_etq +   y  *(img_src->my_width +2) + (x+1));
                const int SE = *(pix_etq + (y+1)*(img_src->my_width +2) + (x+1));
                const int S  = *(pix_etq + (y+1)*(img_src->my_width +2) +   x  );
                const int SO = *(pix_etq + (y+1)*(img_src->my_width +2) + (x-1));
                const int O  = *(pix_etq +   y  *(img_src->my_width +2) + (x-1));

                int base = -1;

                /// NO
                if( NO >= 0 ){
                    base = NO;
                }

                /// N
                if( N >= 0 ){
                    if( (base >= 0) && (base != N) ){
                        ampliarConjunto( val_etiquetas, base, N, max_etiquetas);
                    }else{
                        base = N;
                    }
                }

                /// NE
                if( NE >= 0 ){
                    if( (base >= 0) && (base != NE) ){
                        ampliarConjunto( val_etiquetas, base, NE, max_etiquetas);
                    }else{
                        base = NE;
                    }
                }

                /// E
                if( E >= 0){
                    if( (base >= 0) && (base != E) ){
                        ampliarConjunto( val_etiquetas, base, E, max_etiquetas);
                    }else{
                        base = E;
                    }
                }

                /// SE
                if( SE >= 0 ){
                    if( (base >= 0) && (base != SE) ){
                        ampliarConjunto( val_etiquetas, base, SE, max_etiquetas);
                    }else{
                        base = SE;
                    }
                }

                /// S
                if( S >= 0 ){
                    if( (base >= 0) && (base != S) ){
                        ampliarConjunto( val_etiquetas, base, S, max_etiquetas);
                    }else{
                        base = S;
                    }
                }

                /// SO
                if( SO >= 0 ){
                    if( (base >= 0) && (base != SO) ){
                        ampliarConjunto( val_etiquetas, base, SO, max_etiquetas);
                    }else{
                        base = SO;
                    }
                }

                /// O
                if( O >= 0 ){
                    if( (base >= 0) && (base != O) ){
                        ampliarConjunto( val_etiquetas, base, O, max_etiquetas);
                    }else{
                        base = O;
                    }
                }

                if( base >= 0 ){
                    *(pix_etq + y*(img_src->my_width +2) + x) = base;
                }else{
                    n_conjunto++;
                    *(val_etiquetas + n_conjunto) = n_conjunto;
                    *(pix_etq + y*(img_src->my_width +2) + x) = n_conjunto;
                }
            }
        }
    }

    DEB_MSG("Quedaron " << n_conjunto << " activados");
    n_conjunto++;

    // Indicar a que conjunto pertenece cada pixel:
    memset( conjuntos, -1, mis_rows_cols*sizeof(int));
    for( int y = 0; y < img_src->my_height; y++){
        for( int x = 0; x < img_src->my_width; x++){
            const int etiqueta = *(pix_etq + (y+1)*(img_src->my_width +2) + (x+1));
            if( etiqueta >= 0 ){
                *(conjuntos + x + y*img_src->my_width) = *(val_etiquetas + etiqueta);
            }
        }
    }
    delete [] val_etiquetas;
    delete [] pix_etq;


    /// Contar cuantos elementos hay en cada grupo:
    unsigned int *n_etiquetados = new unsigned int [n_conjunto+1];
    memset(n_etiquetados + 1, 0, n_conjunto*sizeof(unsigned int));
    n_etiquetados[0] = n_conjunto;

    for( int xy = 0; xy < mis_rows_cols; xy++){
        if( *(conjuntos + xy) >= 0 ){
            *(n_etiquetados + 1 + *(conjuntos + xy)) = *(n_etiquetados + 1 + *(conjuntos + xy)) + 1;
        }
    }

    return n_etiquetados;
}




/*  Metodo:  sklMask
    Funcion: Mascara usada para la extraccion del esqueleto.
*/
inline unsigned char IMGVTK::sklMask(IMGCONT *img_src, const int x, const int y) {
	return   1 * (*(img_src->my_img_data + (x - 1) + (y - 1)*img_src->my_width) > 0.0) + /* P2 */
		2 * (*(img_src->my_img_data + x + (y - 1)*img_src->my_width) > 0.0) + /* P3 */
		4 * (*(img_src->my_img_data + (x + 1) + (y - 1)*img_src->my_width) > 0.0) + /* P4 */
		8 * (*(img_src->my_img_data + (x + 1) + y  *img_src->my_width) > 0.0) + /* P5 */
		16 * (*(img_src->my_img_data + (x + 1) + (y + 1)*img_src->my_width) > 0.0) + /* P6 */
		32 * (*(img_src->my_img_data + x + (y + 1)*img_src->my_width) > 0.0) + /* P7*/
		64 * (*(img_src->my_img_data + (x - 1) + (y + 1)*img_src->my_width) > 0.0) + /* P8 */
		128 * (*(img_src->my_img_data + (x - 1) + y  *img_src->my_width) > 0.0);  /* P9 */
}




/*  Metodo: lengthFiltering
    Funcion: Filtra los conjuntos con numero de pixeles menor a 'min_length'
*/
void IMGVTK::lengthFilter(IMGCONT *img_src, const unsigned int min_length, ALG_CONJUNTOS mi_alg)
{
    int *mis_conjuntos = new int [img_src->my_height * img_src->my_width];

    unsigned int *mis_n_etiquetados = NULL;

    switch (mi_alg) {
    case DINAMICO:
        mis_n_etiquetados = conjuntosConexosDinamico(img_src, mis_conjuntos);
        break;
    case ITERATIVO:
        mis_n_etiquetados = conjuntosConexos(img_src, mis_conjuntos);
        break;
    }

    for( int xy = 0; xy < (img_src->my_height * img_src->my_width); xy++){
        if( ( *(mis_conjuntos + xy) >= 0) && ( *(mis_n_etiquetados + *(mis_conjuntos + xy ) + 1 ) < min_length)){
            *(img_src->my_img_data + xy) = 0.0;
        }
    }

    delete [] mis_conjuntos;
    delete [] mis_n_etiquetados;
}



/*  Metodo: lengthFiltering (Publica)
    Funcion: Filtra los conjuntos con numero de pixeles menor a 'min_length'
*/
void IMGVTK::lengthFilter(IMG_IDX img_idx, const int min_length, ALG_CONJUNTOS mi_alg){

    IMGCONT *img_tmp = NULL;

	switch (img_idx) {
	case BASE:
		img_tmp = my_base;
		break;
	case GROUNDTRUTH:
		img_tmp = my_groundtruth;
		break;
	case MASK:
		img_tmp = my_mask;
		break;
	case SKELETON:
		img_tmp = my_skeleton;
		break;
	case SEGMENT:
		img_tmp = my_response;
		break;
	case THRESHOLD:
		img_tmp = my_segmented;
		break;
	case BORDERS:
		img_tmp = my_boundaries;
		break;
	case MAPDIST:
		img_tmp = my_distmap;
		break;
	}

    lengthFilter(img_tmp, min_length, mi_alg);
}


/*  Metodo: regionFill
    Funcion: Rellena vacios dentro del cuerpo de la arteria segmentada.
*/
void IMGVTK::regionFill( IMGCONT *img_src )
{
    // -------- Mascara 9x9
    for( int y = 0; y < img_src->my_height; y++){
        for( int x = 0; x < img_src->my_width; x++){
            if( *(img_src->my_img_data + x + y*img_src->my_width) < 1.0 ){
				*(img_src->my_img_data + x + y*img_src->my_width) = (regionFilling9(img_src, x, y) ? 1.0 : 0.0);
            }
        }
    }
    // -------- Mascara 7x7
    for( int y = 0; y < img_src->my_height; y++){
        for( int x = 0; x < img_src->my_width; x++){
			if (*(img_src->my_img_data + x + y*img_src->my_width) < 1.0) {
				*(img_src->my_img_data + x + y*img_src->my_width) = (regionFilling7(img_src, x, y) ? 1.0 : 0.0);
			}
        }
    }
    // -------- Mascara 5x5
    for( int y = 0; y < img_src->my_height; y++){
        for( int x = 0; x < img_src->my_width; x++){
			if (*(img_src->my_img_data + x + y*img_src->my_width) < 1.0) {
				*(img_src->my_img_data + x + y*img_src->my_width) = (regionFilling5(img_src, x, y) ? 1.0 : 0.0);
			}
        }
    }
    // -------- Mascara 3x3
    for( int y = 0; y < img_src->my_height; y++){
        for( int x = 0; x < img_src->my_width; x++){
			if (*(img_src->my_img_data + x + y*img_src->my_width) < 1.0) {
				*(img_src->my_img_data + x + y*img_src->my_width) = (regionFilling3(img_src, x, y) ? 1.0 : 0.0);
			}
        }
    }
}




/*  Metodo: regionFill (Publica)
    Funcion: Rellena vacios dentro del cuerpo de la arteria segmentada.
*/
void IMGVTK::regionFill(IMG_IDX img_idx) {

	IMGCONT *img_tmp = NULL;

	switch (img_idx) {
	case BASE:
		img_tmp = my_base;
		break;
	case GROUNDTRUTH:
		img_tmp = my_groundtruth;
		break;
	case MASK:
		img_tmp = my_mask;
		break;
	case SKELETON:
		img_tmp = my_skeleton;
		break;
	case SEGMENT:
		img_tmp = my_response;
		break;
	case THRESHOLD:
		img_tmp = my_segmented;
		break;
	case BORDERS:
		img_tmp = my_boundaries;
		break;
	case MAPDIST:
		img_tmp = my_distmap;
		break;
	}

	regionFill(img_tmp);
}



/*  Metodo:  dilMask
    Funcion: Mascara usada para dilatar unaimagen.
*/
inline unsigned char IMGVTK::dilMask( IMGCONT *mask_dil, const int x, const int y){
    return (*(mask_dil->my_img_data + ( x ) + (y-1)*mask_dil->my_width) > 0.0) +
           (*(mask_dil->my_img_data + (x+1) + ( y )*mask_dil->my_width) > 0.0) +
           (*(mask_dil->my_img_data + ( x ) + (y+1)*mask_dil->my_width) > 0.0) +
           (*(mask_dil->my_img_data + (x-1) + ( y )*mask_dil->my_width) > 0.0);
}



/*  Metodo:  erosionMask
    Funcion: Mascara para erosion usando un disco de radio 5
*/
inline unsigned char IMGVTK::erosionMask(IMGCONT *ptr_tmp, const int x, const int y)
{
	return (*(ptr_tmp->my_img_data + (x - 2) + (y - 4)*ptr_tmp->my_width) > 0.0) +
		(*(ptr_tmp->my_img_data + (x - 1) + (y - 4)*ptr_tmp->my_width) > 0.0) +
		(*(ptr_tmp->my_img_data + (x)+(y - 4)*ptr_tmp->my_width) > 0.0) +
		(*(ptr_tmp->my_img_data + (x + 1) + (y - 4)*ptr_tmp->my_width) > 0.0) +
		(*(ptr_tmp->my_img_data + (x + 2) + (y - 1)*ptr_tmp->my_width) > 0.0) +
		(*(ptr_tmp->my_img_data + (x - 3) + (y - 3)*ptr_tmp->my_width) > 0.0) +
		(*(ptr_tmp->my_img_data + (x - 2) + (y - 3)*ptr_tmp->my_width) > 0.0) +
		(*(ptr_tmp->my_img_data + (x - 1) + (y - 3)*ptr_tmp->my_width) > 0.0) +
		(*(ptr_tmp->my_img_data + (x)+(y - 3)*ptr_tmp->my_width) > 0.0) +
		(*(ptr_tmp->my_img_data + (x + 1) + (y - 3)*ptr_tmp->my_width) > 0.0) +
		(*(ptr_tmp->my_img_data + (x + 2) + (y - 2)*ptr_tmp->my_width) > 0.0) +
		(*(ptr_tmp->my_img_data + (x + 3) + (y - 3)*ptr_tmp->my_width) > 0.0) +
		(*(ptr_tmp->my_img_data + (x - 4) + (y - 2)*ptr_tmp->my_width) > 0.0) +
		(*(ptr_tmp->my_img_data + (x - 3) + (y - 2)*ptr_tmp->my_width) > 0.0) +
		(*(ptr_tmp->my_img_data + (x - 2) + (y - 2)*ptr_tmp->my_width) > 0.0) +
		(*(ptr_tmp->my_img_data + (x - 1) + (y - 2)*ptr_tmp->my_width) > 0.0) +
		(*(ptr_tmp->my_img_data + (x)+(y - 2)*ptr_tmp->my_width) > 0.0) +
		(*(ptr_tmp->my_img_data + (x + 1) + (y - 2)*ptr_tmp->my_width) > 0.0) +
		(*(ptr_tmp->my_img_data + (x + 2) + (y - 3)*ptr_tmp->my_width) > 0.0) +
		(*(ptr_tmp->my_img_data + (x + 3) + (y - 2)*ptr_tmp->my_width) > 0.0) +
		(*(ptr_tmp->my_img_data + (x + 4) + (y - 2)*ptr_tmp->my_width) > 0.0) +
		(*(ptr_tmp->my_img_data + (x - 4) + (y - 1)*ptr_tmp->my_width) > 0.0) +
		(*(ptr_tmp->my_img_data + (x - 3) + (y - 1)*ptr_tmp->my_width) > 0.0) +
		(*(ptr_tmp->my_img_data + (x - 2) + (y - 1)*ptr_tmp->my_width) > 0.0) +
		(*(ptr_tmp->my_img_data + (x - 1) + (y - 1)*ptr_tmp->my_width) > 0.0) +
		(*(ptr_tmp->my_img_data + (x)+(y - 1)*ptr_tmp->my_width) > 0.0) +
		(*(ptr_tmp->my_img_data + (x + 1) + (y - 1)*ptr_tmp->my_width) > 0.0) +
		(*(ptr_tmp->my_img_data + (x + 2) + (y - 4)*ptr_tmp->my_width) > 0.0) +
		(*(ptr_tmp->my_img_data + (x + 3) + (y - 1)*ptr_tmp->my_width) > 0.0) +
		(*(ptr_tmp->my_img_data + (x + 4) + (y - 1)*ptr_tmp->my_width) > 0.0) +
		(*(ptr_tmp->my_img_data + (x - 4) + (y)*ptr_tmp->my_width) > 0.0) +
		(*(ptr_tmp->my_img_data + (x - 3) + (y)*ptr_tmp->my_width) > 0.0) +
		(*(ptr_tmp->my_img_data + (x - 2) + (y)*ptr_tmp->my_width) > 0.0) +
		(*(ptr_tmp->my_img_data + (x - 1) + (y)*ptr_tmp->my_width) > 0.0) +
		(*(ptr_tmp->my_img_data + (x + 1) + (y)*ptr_tmp->my_width) > 0.0) +
		(*(ptr_tmp->my_img_data + (x + 2) + (y)*ptr_tmp->my_width) > 0.0) +
		(*(ptr_tmp->my_img_data + (x + 3) + (y)*ptr_tmp->my_width) > 0.0) +
		(*(ptr_tmp->my_img_data + (x + 4) + (y)*ptr_tmp->my_width) > 0.0) +
		(*(ptr_tmp->my_img_data + (x - 4) + (y + 1)*ptr_tmp->my_width) > 0.0) +
		(*(ptr_tmp->my_img_data + (x - 3) + (y + 1)*ptr_tmp->my_width) > 0.0) +
		(*(ptr_tmp->my_img_data + (x - 2) + (y + 1)*ptr_tmp->my_width) > 0.0) +
		(*(ptr_tmp->my_img_data + (x - 1) + (y + 1)*ptr_tmp->my_width) > 0.0) +
		(*(ptr_tmp->my_img_data + (x)+(y + 1)*ptr_tmp->my_width) > 0.0) +
		(*(ptr_tmp->my_img_data + (x + 1) + (y + 1)*ptr_tmp->my_width) > 0.0) +
		(*(ptr_tmp->my_img_data + (x + 2) + (y + 1)*ptr_tmp->my_width) > 0.0) +
		(*(ptr_tmp->my_img_data + (x + 3) + (y + 1)*ptr_tmp->my_width) > 0.0) +
		(*(ptr_tmp->my_img_data + (x + 4) + (y + 1)*ptr_tmp->my_width) > 0.0) +
		(*(ptr_tmp->my_img_data + (x - 4) + (y + 2)*ptr_tmp->my_width) > 0.0) +
		(*(ptr_tmp->my_img_data + (x - 3) + (y + 2)*ptr_tmp->my_width) > 0.0) +
		(*(ptr_tmp->my_img_data + (x - 2) + (y + 2)*ptr_tmp->my_width) > 0.0) +
		(*(ptr_tmp->my_img_data + (x - 1) + (y + 2)*ptr_tmp->my_width) > 0.0) +
		(*(ptr_tmp->my_img_data + (x)+(y + 2)*ptr_tmp->my_width) > 0.0) +
		(*(ptr_tmp->my_img_data + (x + 1) + (y + 2)*ptr_tmp->my_width) > 0.0) +
		(*(ptr_tmp->my_img_data + (x + 2) + (y + 2)*ptr_tmp->my_width) > 0.0) +
		(*(ptr_tmp->my_img_data + (x + 3) + (y + 2)*ptr_tmp->my_width) > 0.0) +
		(*(ptr_tmp->my_img_data + (x + 4) + (y + 2)*ptr_tmp->my_width) > 0.0) +
		(*(ptr_tmp->my_img_data + (x - 3) + (y + 3)*ptr_tmp->my_width) > 0.0) +
		(*(ptr_tmp->my_img_data + (x - 2) + (y + 3)*ptr_tmp->my_width) > 0.0) +
		(*(ptr_tmp->my_img_data + (x - 1) + (y + 3)*ptr_tmp->my_width) > 0.0) +
		(*(ptr_tmp->my_img_data + (x)+(y + 3)*ptr_tmp->my_width) > 0.0) +
		(*(ptr_tmp->my_img_data + (x + 1) + (y + 3)*ptr_tmp->my_width) > 0.0) +
		(*(ptr_tmp->my_img_data + (x + 2) + (y + 3)*ptr_tmp->my_width) > 0.0) +
		(*(ptr_tmp->my_img_data + (x + 3) + (y + 3)*ptr_tmp->my_width) > 0.0) +
		(*(ptr_tmp->my_img_data + (x - 2) + (y + 4)*ptr_tmp->my_width) > 0.0) +
		(*(ptr_tmp->my_img_data + (x - 1) + (y + 4)*ptr_tmp->my_width) > 0.0) +
		(*(ptr_tmp->my_img_data + (x)+(y + 4)*ptr_tmp->my_width) > 0.0) +
		(*(ptr_tmp->my_img_data + (x + 1) + (y + 4)*ptr_tmp->my_width) > 0.0) +
		(*(ptr_tmp->my_img_data + (x + 2) + (y + 4)*ptr_tmp->my_width) > 0.0);
}



/*  Metodo: erosionar
    Funcion: Rellena vacios dentro del cuerpo de la arteria segmentada.
*/
void IMGVTK::erosionar( IMGCONT *img_src )
{
	IMGCONT *ptr_tmp = new IMGCONT(img_src->my_height, img_src->my_width, 4, 4, 1.0);
	*ptr_tmp = *img_src;

    for( int y = 0; y < img_src->my_height; y++){
        for( int x = 0; x < img_src->my_width; x++){
            if (*(ptr_tmp->my_img_data + (x+4) + (y+4)*(ptr_tmp->my_width)) > 0){
                const unsigned char resp = erosionMask(ptr_tmp, x+4, y+4);
                if( (0 < resp) && (resp < 68) ){
                    *(img_src->my_img_data + x + y*img_src->my_width) = 0.0;
                }
            }
        }
    }

    delete ptr_tmp;
}



/*  Metodo: maskFOV
    Funcion: Obtiene la mascara del contorno
*/
IMGVTK::IMGCONT * IMGVTK::maskFOV( IMGCONT *img_src )
{

	IMGCONT *new_mask = new IMGCONT(*img_src);

    // Se umbraliza al 0.1 la imagen original
    for( int xy = 0; xy < img_src->my_height * img_src->my_width; xy++){
		*(new_mask->my_img_data + xy) = (*(new_mask->my_img_data + xy) < 0.1) ? 0.0 : 1.0;
    }

    // Se eliminan los conjuntos pequeños
    lengthFilter(new_mask, 1000, ITERATIVO);

    // Se erosiona la mascara:
    erosionar(new_mask);

    // Se eliminan los conjuntos grandes que no esten en las esquinas:
    // Se extraen las etiquetas de los conjuntos que se conectan a las esquinas:
    for( int xy = 0; xy < new_mask->my_width * new_mask->my_height; xy++){
        *(new_mask->my_img_data + xy) = 1.0 - *(new_mask->my_img_data + xy);
    }

    int *mis_conjuntos = new int [new_mask->my_width * new_mask->my_height];

    // Se buscan los conjuntos que no esten en las esquinas para eliminarlos
	unsigned int *mis_n_etiquetados = conjuntosConexos(new_mask, mis_conjuntos);
    delete [] mis_n_etiquetados;


    //// etiqueta de los conjuntos donde existe una esquina
	const int NO = *(mis_conjuntos);
	const int NE = *(mis_conjuntos + new_mask->my_width - 1);
	const int SO = *(mis_conjuntos + (new_mask->my_height - 1)*new_mask->my_width);
	const int SE = *(mis_conjuntos + new_mask->my_height * new_mask->my_width - 1);

    for( int xy = 0; xy < (new_mask->my_height * new_mask->my_width); xy++){
        if( (*(mis_conjuntos + xy) >= 0) && ((mis_conjuntos[xy] == NO) || (mis_conjuntos[xy] == NE) || 
			(mis_conjuntos[xy] == SO) || (mis_conjuntos[xy] == SE)) ){
            *(new_mask->my_img_data + xy) = 0.0;
        }else{
			*(new_mask->my_img_data + xy) = 1.0;
        }
    }

    delete [] mis_conjuntos;

	return new_mask;
}



/*  Metodo: fillMask
    Funcion: Se llenan los espacios de la mascara con la media de la imagen circundante.
*/
void IMGVTK::fillMask( IMGCONT *img_src, IMGCONT *mask_src ){

    PIX_PAR par_tmp;
    par_tmp.pix_tipo = PIX_CROSS;

    IMGCONT *mask_dil = new IMGCONT(img_src->my_height+2, img_src->my_width+2, 1, 1);
	*mask_dil = *mask_src;

    int iter = 0;
    while( iter < mask_src->my_width ){
        // Obtener el borde de la mascara realizando una dilatacion:
        std::vector< PIX_PAR > borde;

        // Si es la primer pasada, se verifican todos los pixeles para encontrar el borde
        for( int y = 0; y < mask_src->my_height; y++){
            for( int x = 0; x < mask_src->my_width; x++){
                if( *(mask_dil->my_img_data + (x+1) + (y+1)*(mask_dil->my_width)) < 1.0 ){
                    // Si es HIT se considera como borde, si es FIT no:
					const unsigned char resp = dilMask(mask_dil, x + 1, y + 1);
                    if( resp > 0){
                        // Se guardan las coordenadas de los pixeles que forman el borde:
                        par_tmp.x = x;
                        par_tmp.y = y;
                        borde.push_back( par_tmp );
                    }
                }
            }
        }

        const int n_borde = (int)borde.size();
        DEB_MSG("Iter: " << iter << ", pixeles en el borde: " <<  n_borde);

        //// Si ya no existen pixeles en el borde, se termina el ciclo:
        if( n_borde == 0){
            break;
        }

        // Para cada pixel en el borde, se calcula la media usando una ventana de 21 x 21 pixeles.
        for( int b = 0; b < n_borde; b++ ){
            const int x_act = (int)borde[b].x;
            const int y_act = (int)borde[b].y;

            const int offset_x_izq = (x_act < 10) ?
				0 :
				(x_act - 10);

            const int offset_x_der = (x_act >= (mask_src->my_width - 10)) ? 
				(mask_src->my_width-1) : 
				(x_act + 10);

            const int offset_y_sup = (y_act < 10) ?
				0 : 
				(y_act - 10);

            const int offset_y_inf = (y_act >= (mask_src->my_height - 10)) ?
				(mask_src->my_height-1) :
				(y_act + 10);

            double suma = 0.0;
            int n_vecinos = 0;

            for( int y = offset_y_sup; y <= offset_y_inf; y++){
                for( int x = offset_x_izq; x <= offset_x_der; x++){
                    if( *(mask_dil->my_img_data + (x+1) + (y+1)*(mask_dil->my_width)) > 0.0 ){
                        suma += *(img_src->my_img_data + x + y*img_src->my_width);
                        n_vecinos++;
                    }
                }
            }

            *(img_src->my_img_data + x_act + y_act*img_src->my_width) = suma / n_vecinos;
            *(mask_dil->my_img_data + (x_act+1) + (y_act+1)*(mask_dil->my_width)) = 1.0;
        }
        iter++;
    }

    delete mask_dil;
}





/*  Metodo: grafoSkeleton
    Funcion: Genera un grafo a partir del esqueleto.
*/
IMGVTK::PIX_PAR* IMGVTK::grafoSkeleton(double *skl_tmp, const int x, const int y, int *nivel, const unsigned char *lutabla, bool *visitados) {
	/*
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
	const int max_x = ((x + max_dist - 1) > cols) ? cols : (x + max_dist - 1);
	const int min_y = ((y - max_dist - 1) < 0) ? 0 : (y - max_dist - 1);
	const int max_y = ((y + max_dist - 1) > rows) ? rows : (y + max_dist - 1);

	double dist, radio = MY_INF;
	bool visitado = false;
	double x_r, y_r;
	for( int yy = min_y; yy < max_y; yy++){
		for( int xx = min_x; xx < max_x; xx++){
			dist = (double)(yy - y + 0.5)*(double)(yy - y + 0.5) + (double)(xx - x + 0.5)*(double)(xx - x + 0.5);
			if( (borders_ptr[xx + yy*cols] > 0.0) && (dist < radio) ){
				radio = dist;
				x_r = (double)xx;
				y_r = (double)yy;
				visitado = true;
			}
		}
	}

	if( !visitado ){
		char mensaje[512] = "\nEl pixel [XXX, YYY] no tiene vecinos en la imagen de bordes a: DDD pixeles a la redonda\n";
#if defined(_WIN32) || defined(_WIN64)
		sprintf_s(mensaje, 512 * sizeof(char), "\nEl pixel [%i, %i] no tiene vecinos en la imagen de bordes a: %i pixeles a la redonda\n", x, y, max_dist);
#else
		sprintf( mensaje, "\nEl pixel [%i, %i] no tiene vecinos en la imagen de bordes a: %i pixeles a la redonda\n", x, y, max_dist);
#endif
		escribirLog(mensaje);
	}

	temp->radio = sqrt(radio) * pixX;
	temp->y_r = (y_r - (double)rows/2)*pixY;
	temp->x_r = (x_r - (double)cols/2)*pixX;
	temp->alpha = atan2(temp->y_r - temp->y, temp->x_r - temp->x);// + MY_PI / 2.0;

	switch( lutabla[ resp ] ){
		case (unsigned char)1:{ /* END point*
			temp->pix_tipo = PIX_END;
			break;
		}
		case (unsigned char)2:{ /* BRANCH point *
			temp->pix_tipo = PIX_BRANCH;
			break;
		}
		case (unsigned char)3:{ /* CROSS point *
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
	*/
	return NULL;
}




/*  Metodo: extraerCaract
    Funcion: Extrae los pixeles caracteristicos (end y branch points) a partir del esqueleot de la imagen.
*/
void IMGVTK::extraerCaract( IMG_IDX img_idx ){
	/*
    if( !my_boundaries ){
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
    memcpy( skl_tmp, img_src, (rows+2)*(cols+2)*sizeof(double));


    /// Buscar un punto 'end' del esqueleto y empezar a generar el grafo a aprtir de ahi.
    bool *visitados = new bool [rows_cols];
    memset(visitados, 0, rows_cols*sizeof(bool));

    int x_ini, y_ini, xy = cols+2;
    unsigned char resp;
    do{
        xy++;
        x_ini = xy % (cols+2);
        y_ini = (int)((double)xy / ((double)cols+2));
        resp = sklMask( img_src, x_ini, y_ini, cols+2, rows) * (*(img_src + xy) > 0.0);
    }while( tabla[resp] != (unsigned char)1 );

    *(skl_tmp + xy ) = 2.0;
    int nivel = 0;

    pix_caract = grafoSkeleton(skl_tmp, x_ini, y_ini, &nivel, tabla, visitados);

    n_niveles = nivel;

    DEB_MSG("Encontrados " << n_niveles << " niveles");

    delete [] visitados;
    delete [] skl_tmp;
	*/
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
void IMGVTK::skeletonization(IMG_IDX img_idx) {

	IMGCONT *img_tmp = NULL;

	switch (img_idx) {
	case BASE:
		img_tmp = my_base;
		break;
	case GROUNDTRUTH:
		img_tmp = my_groundtruth;
		break;
	case MASK:
		img_tmp = my_mask;
		break;
	case SEGMENT:
		img_tmp = my_response;
		break;
	case MAPDIST:
		img_tmp = my_distmap;
		break;
	case THRESHOLD:
		img_tmp = my_segmented;
		break;
	case BORDERS:
		img_tmp = my_boundaries;
		break;
	}

	if (!my_skeleton) {
		my_skeleton = new IMGCONT(img_tmp->my_height, img_tmp->my_width, 1, 1);
	}
	*my_skeleton = *img_tmp;

	IMGCONT *skl_mark = new IMGCONT(*my_skeleton);

	const unsigned char tabla[] = {
		0, 0, 0, 1, 0, 0, 1, 3, 0, 0, 3, 1, 1, 0, 1, 3, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 2, 0, 3, 0, 3, 3,
		0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 3, 0, 2, 2,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		2, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 2, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 3, 0, 2, 0,
		0, 0, 3, 1, 0, 0, 1, 3, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
		3, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		2, 3, 1, 3, 0, 0, 1, 3, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		2, 3, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 3, 3, 0, 1, 0, 0, 0, 0, 2, 2, 0, 0, 2, 0, 0, 0 };
	int n_borrado;

	do {
		n_borrado = 0;
		// Primer paso:
		for (int y = 1; y <= my_skeleton->my_height; y++) {
			for (int x = 1; x <= my_skeleton->my_width; x++) {
				if (*(my_skeleton->my_img_data + x + y*my_skeleton->my_width) > 0.0) {
					const unsigned char resp = *(tabla + sklMask(my_skeleton, x, y));
					if ((resp == 1) || (resp == 3)) {
						*(skl_mark->my_img_data + x + y * skl_mark->my_width) = 0.0;
						n_borrado++;
					}
				}
			}
		}

		*my_skeleton = *skl_mark;

		// Segundo paso:
		for (int y = 1; y <= rows; y++) {
			for (int x = 1; x <= cols; x++) {
				if (*(my_skeleton->my_img_data + x + y*my_skeleton->my_width) > 0.0) {
					unsigned char resp = *(tabla + sklMask(my_skeleton, x, y));
					if (resp == 2 || resp == 3) {
						*(skl_mark->my_img_data + x + y*skl_mark->my_width) = 0.0;
						n_borrado++;
					}
				}
			}
		}


		*my_skeleton = *skl_mark;

	} while (n_borrado > 0);

	delete skl_mark;
	extraerCaract(img_idx);
}




/*  Metodo: definirMask
    Funcion: Define una mascara para normalizar los pixeles de la imagen.
*/
IMGVTK::IMGCONT * IMGVTK::definirMask( IMGCONT * img_src )
{

    IMGCONT *new_mask = maskFOV( img_src );
    fillMask( img_src, new_mask );

    escribirLog( "Mascara generada exitosamente.\n" );

	return new_mask;
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
        const int clase_i = (int)((clases-1) * (*(img_ptr + xy) - min)/(max - min + 1e-12));
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

    return umbral;
}





/*  Metodo: umbralizarRIDCAL
    Funcion: Utiliza el metodo de Otsu para umbralizar la imagen y separar el fondo y el primer plano de la imagen.
*/
double IMGVTK::umbralizarRIDCAL( const double *img_ptr, const double min, const double max){
    double umbral_nuevo = 0.0;
    // Obtener la media de la read_intensity de la imagen como primera aproximacion del umbral
    const double rango = 1.0 / (max - min);
    const double fraccion = 1.0 / (double)rows_cols;

    for( int xy = 0; xy < rows_cols; xy ++){
        umbral_nuevo += (*(img_ptr + xy) - min) * rango;
    }
    umbral_nuevo *= fraccion;

    double umbral;
    while( 1 ){
        const double umbral_previo = umbral_nuevo;
        umbral = umbral_previo;

        double m_abajo = 0.0;
        double m_arriba = 0.0;
        int n_abajo = 0, n_arriba = 0;

        for( int xy = 0; xy < rows_cols; xy++){

            if( (*(img_ptr + xy) - min) * rango <= umbral ){
                m_abajo += (*(img_ptr + xy) - min) * rango;
                n_abajo++;
            }else{
                m_arriba += (*(img_ptr + xy) - min) * rango;
                n_arriba++;
            }
        }

        m_abajo /= (double)n_abajo;
        m_arriba /= (double)n_arriba;

        umbral_nuevo = (m_arriba + m_abajo) / 2.0;

        if( fabs(umbral_nuevo - umbral_previo) < 1e-5 ){
            break;
        }
    }

    return umbral;
}




/*  Metodo: umbralizar
    Funcion: Utiliza el metodo de Otsu o Ridler & Calvard para umbralizar la imagen y separar el fondo y el primer plano de la imagen.
*/
void IMGVTK::umbralizar(IMG_IDX img_idx, const TIPO_UMBRAL tipo_umb, const double nivel){
	/*
    if( !threshold_ptr ){
#ifdef BUILD_VTK_VERSION
        threshold = vtkSmartPointer< vtkImageData >::New();
        threshold->SetExtent(0, cols - 1, 0, rows - 1, 0, 0);
        threshold->AllocateScalars(VTK_DOUBLE, 1);
        threshold->SetOrigin(0.0, 0.0, 0.0);
        threshold->SetSpacing(1.0, 1.0, 1.0);

        threshold_ptr = static_cast< double* >(threshold->GetScalarPointer(0, 0, 0));
#else
        threshold_ptr = new double [rows_cols];
#endif
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

    DEB_MSG("Umbral: " << umbral);

    // Se umbraliza la imagen con el valor optimo encontrado:
    for( int xy = 0; xy < rows_cols; xy++){
        *(threshold_ptr + xy) = ( ((*(img_ptr + xy) - min) / (max - min)) >= umbral) ? 1.0 : 0.0;
    }
	*/
}





/*  Metodo: Cargar
    Funcion: Retorna la exactitud del clasificador resultante contra el ground-truth
*/
double IMGVTK::medirExactitud(){
	/*
    if( !gt_ptr ){
        char mensaje_error[] = "<<ERROR: No se cargo la imagen ground-truth\n";
        escribirLog( mensaje_error);
        return 0.0;
    }

    int FP = 0, FN = 0, TP = 0, TN = 0;

    for( int xy = 0; xy < rows_cols; xy++){
        if( *(mask_ptr + xy) > 0.5 ){
            if( *(gt_ptr + xy) > 0.5 ){
                if( *(threshold_ptr + xy) > 0.5 ){
                    TP++;
                }else{
                    FN++;
                }
            }else{
                if( *(threshold_ptr + xy) < 0.5 ){
                    TN++;
                }else{
                    FP++;
                }
            }
        }
    }

    return (double)(TP + TN) / (double)(TP + TN + FP + FN);
	*/
	return 0.0;
}


/*
*  Source: readpng.c form the libpng documentation
*  
*/
IMGVTK::IMGCONT* IMGVTK::CargarPNG(const char *ruta_origen)
{
	return NULL;
}




IMGVTK::IMGCONT * IMGVTK::CargarPGM(const char *ruta_origen)
{
#if defined(_WIN32) || defined(_WIN64)
	FILE *img_file;
	fopen_s(&img_file, ruta_origen, "r");
#else
	FILE *img_file = fopen(ruta_origen, "r");
#endif

	char temp_str[512] = "";
	fgets(temp_str, 512, img_file);  // Leer 'Magic number'
	fgets(temp_str, 512, img_file);  // Leer Comentario

	double max_intensity;
	unsigned int height;
	unsigned int width;

#if defined(_WIN32) || defined(_WIN64)
	fscanf_s(img_file, "%i", &width);    // Leer ancho
	fscanf_s(img_file, "%i", &height);    // Leer alto
	fscanf_s(img_file, "%lf", &max_intensity);  // Leer la escala de read_intensity maxima
#else
	fscanf(img_file, "%i", &width);    // Leer ancho
	fscanf(img_file, "%i", &height);    // Leer alto
	fscanf(img_file, "%lf", &max_intensity);  // Leer la escala de read_intensity maxima
#endif
	
	IMGCONT *new_img = new IMGCONT(height, width);

	int read_intensity;
	for (int xy = 0; xy < (height* width); xy++) {
#if defined(_WIN32) || defined(_WIN64)
		fscanf_s(img_file, "%i", &read_intensity);
#else
		fscanf(img_file, "%i", &read_intensity);
#endif
		*(new_img->my_img_data + xy) = (double)read_intensity / max_intensity;
	}

	fclose(img_file);

	return new_img;
}



/*  Metodo: CargarDICOM   */
IMGVTK::IMGCONT * IMGVTK::CargarDICOM(const char *ruta_origen, const unsigned int nivel)
{
#ifdef BUILD_GDCM_VERSION
        gdcm::ImageReader DICOMreader;
        DICOMreader.SetFileName( ruta_origen );

        DICOMreader.Read();

        gdcm::File &file = DICOMreader.GetFile();
        gdcm::DataSet &ds = file.GetDataSet();

////---------- Extraer SID (Distancia de la fuente al detector):
        {
            const gdcm::DataElement &de = ds.GetDataElement( gdcm::Tag(0x18, 0x1110) );
            const gdcm::ByteValue *bv = de.GetByteValue();
            if( bv ){
                std::string strm(bv->GetPointer(), bv->GetLength());
                SID = atof( strm.c_str() );
            }
        }
////---------- Extraer SOD (Distancia de la fuente al paciente):
        {
            const gdcm::DataElement &de = ds.GetDataElement( gdcm::Tag(0x18, 0x1111) );
            const gdcm::ByteValue *bv = de.GetByteValue();
            if( bv ){
                std::string strm(bv->GetPointer(), bv->GetLength());
                SOD = atof( strm.c_str() );
            }
        }
////---------- Calcular el DDP (Distancia del detector al paciente):
        DDP = SID - SOD;
////---------- Calcular el Magnification Factor:
        const double Magnification = SID / SOD;

////---------- Extraer pixY\pixX (Distancia entre el centro de los pixeles en el eje Y\X): 
        {
            const gdcm::DataElement &de = ds.GetDataElement( gdcm::Tag(0x18, 0x1164) );
            const gdcm::ByteValue *bv = de.GetByteValue();
            if( bv ){
                std::string strm(bv->GetPointer(), bv->GetLength());
                char *pixXYstr = new char [bv->GetLength()];
                memcpy(pixXYstr, strm.c_str(), bv->GetLength() * sizeof(char ));
                char *tmp = strchr(pixXYstr,'\\');
                pixXYstr[ tmp - pixXYstr ] = '\0';
                pixY = atof(pixXYstr);// / Magnification;
                pixX = atof(tmp+1);// / Magnification;
                delete [] pixXYstr;
            }
        }
////---------- Extraer LAO/RAO (Angulo del detector en direccion izquierda(-) a derecha(+)):
        {
            const gdcm::DataElement &de = ds.GetDataElement( gdcm::Tag(0x18, 0x1510) );
            const gdcm::ByteValue *bv = de.GetByteValue();
            if( bv ){
                std::string strm(bv->GetPointer(), bv->GetLength());
                LAORAO = atof( strm.c_str() );
            }
        }
////---------- Extraer CRA/CAU (Angulo del detector en direccion cranial(-) a caudal(+)):
        {
            const gdcm::DataElement &de = ds.GetDataElement( gdcm::Tag(0x18, 0x1511) );
            const gdcm::ByteValue *bv = de.GetByteValue();
            if( bv ){
                std::string strm(bv->GetPointer(), bv->GetLength());
                CRACAU = atof( strm.c_str() );
            }
        }

////---------- Extraer Window Width (Ancho de la ventana donde se encuentran los valores de interes):
        {
            const gdcm::DataElement &de = ds.GetDataElement( gdcm::Tag(0x18, 0x7030) );
            const gdcm::ByteValue *bv = de.GetByteValue();
            if( bv ){
                std::string strm(bv->GetPointer(), bv->GetLength());
            }
        }

////---------- Extraer Source to Isocenter (Distancia de la fuente al isocentro):
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
            }
        }
////---------- Extraer Window Width (Ancho de la ventana donde se encuentran los valores de interes):
        {
            const gdcm::DataElement &de = ds.GetDataElement( gdcm::Tag(0x28, 0x1051) );
            const gdcm::ByteValue *bv = de.GetByteValue();
            if( bv ){
                std::string strm(bv->GetPointer(), bv->GetLength());
                WWidth = atof( strm.c_str() );
            }
        }
////---------- Extraer el ECG:
        int ecg_dim = 0, ecg_np = 0;
        {
            const gdcm::DataElement &de = ds.GetDataElement( gdcm::Tag(0x5000, 0x0005) );
            const gdcm::ByteValue *bv = de.GetByteValue();
            if( bv ){
                gdcm::Element<gdcm::VR::US, gdcm::VM::VM1_n> el;
                el.Set( de.GetValue() );
                ecg_dim = el.GetValue();
            }
        }
        {
            const gdcm::DataElement &de = ds.GetDataElement( gdcm::Tag(0x5000, 0x0010) );
            const gdcm::ByteValue *bv = de.GetByteValue();
            if( bv ){
                gdcm::Element<gdcm::VR::US, gdcm::VM::VM1_n> el;
                el.Set( de.GetValue() );
                ecg_np = el.GetValue();
            }
        }
        int n_niveles = 1;
        {
            const gdcm::DataElement &de = ds.GetDataElement( gdcm::Tag(0x8, 0x2143) );
            const gdcm::ByteValue *bv = de.GetByteValue();
            if( bv ){
                std::string strm(bv->GetPointer(), bv->GetLength());
                n_niveles = atof( strm.c_str() );
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
        char *buffer = new char[gimage.GetBufferLength()];
        gimage.GetBuffer(buffer);

        const unsigned int* dimension = gimage.GetDimensions();
        unsigned int width = dimension[0];
        unsigned int height = dimension[1];
        const int mis_niveles = dimension[2];

        // Alojar memoria para la imagen:
		IMGCONT *new_img = new IMGCONT(height, width);

        gdcm::PhotometricInterpretation scl_comps = gimage.GetPhotometricInterpretation();
        gdcm::PixelFormat pix_format = gimage.GetPixelFormat();

        switch(scl_comps){
            case gdcm::PhotometricInterpretation::RGB:{
                    if( pix_format == gdcm::PixelFormat::UINT8 ){
                        for( int y = 0; y < my_height; y++){
                            for( int x = 0; x < my_width; x++){
                                const double pixR = (double)(unsigned char)*(buffer + 3*x   + y*my_width*3 + nivel*mis_rows_cols*3) - WCenter + 0.5;
                                const double pixG = (double)(unsigned char)*(buffer + 3*x+1 + y*my_width*3 + nivel*mis_rows_cols*3) - WCenter + 0.5;
                                const double pixB = (double)(unsigned char)*(buffer + 3*x+2 + y*my_width*3 + nivel*mis_rows_cols*3) - WCenter + 0.5;
                                double pix = (0.297)*pixR + (0.589)*pixG + (0.114)*pixB;
                                if( pix <= -((WWidth-1) / 2)){
                                    pix = 0.0;
                                }else if(pix > ((WWidth-1) / 2)){
                                    pix = 0.0;
                                }else{
                                    pix = pix / (WWidth -1) + 0.5;
                                }
                                *(new_img->my_img_data + (my_height-y-1)*my_height + x ) = pix; // 255.0;
                            }
                        }
                    }else{
                        escribirLog(COLOR_RED"<<ERROR AL LEER ARCHIVO DICOM:" COLOR_BACK_YELLOW " Formato de imagen RGB no soportado>>" COLOR_RESET "\n");
						return NULL;
                    }
                    break;
                }
            case gdcm::PhotometricInterpretation::MONOCHROME1:
            case gdcm::PhotometricInterpretation::MONOCHROME2:{
                    if( pix_format == gdcm::PixelFormat::UINT8 ){
                        for( int y = 0; y < my_height; y++){
                            for( int x = 0; x < my_width; x++){
                                double pix = (double)(unsigned char)*(buffer + nivel*mis_rows_cols + x + y*my_width) - WCenter + 0.5;

                                if( pix <= -((WWidth-1) / 2)){
                                    pix = 0.0;
                                }else if(pix > ((WWidth-1) / 2)){
                                    pix = 1.0;
                                }else{
                                    pix = pix / (WWidth -1) + 0.5;
                                }
                                *(new_img->my_img_data + (my_height-y-1)*my_height + x) = pix; // 255.0;
                            }
                        }

                    }else if( pix_format == gdcm::PixelFormat::UINT16 ){
                        unsigned short *buffer16 = (unsigned short*)buffer;
                        for( int y = 0; y < my_height; y++){
                            for( int x = 0; x < my_width; x++){
                                const double pixR = (double)((unsigned char)*(buffer16 + 3*x   + y*my_width*3 + nivel*mis_rows_cols*3) / 16)  - WCenter + 0.5;
                                const double pixG = (double)((unsigned char)*(buffer16 + 3*x+1 + y*my_width*3 + nivel*mis_rows_cols*3) / 16)  - WCenter + 0.5;
                                const double pixB = (double)((unsigned char)*(buffer16 + 3*x+2 + y*my_width*3 + nivel*mis_rows_cols*3) / 16)  - WCenter + 0.5;
                                double pix = (0.297)*pixR + (0.589)*pixG + (0.114)*pixB;
                                if( pix <= -((WWidth-1) / 2)){
                                    pix = 0.0;
                                }else if(pix > ((WWidth-1) / 2)){
                                    pix = 0.0;
                                }else{
                                    pix = pix / (WWidth -1) + 0.5;
                                }
                                *(new_img->my_img_data + (my_height-y-1)*my_height + x) =  pix; // 255.0;
                            }
                        }

                    }else{
                        escribirLog(COLOR_RED"<<ERROR AL LEER ARCHIVO DICOM:" COLOR_BACK_YELLOW " Formato de imagen RGB no soportado>>" COLOR_RESET "\n");
                    }
                    break;
                }
        }
        delete [] buffer;
#else
		return NULL;
#endif
}




/*  Metodo: Cargar
    Funcion: Cargar la imagen a formato VTK desde un archivo.
*/
void IMGVTK::Cargar(const IMG_IDX img_idx, const char *ruta_origen, const bool enmascarar, const int nivel) {

	IMGCONT *img_tmp = NULL;

	/* Verify the type of image to be read */
	const unsigned int len_path = (unsigned int)strlen(ruta_origen);
	if (strcmp(ruta_origen + len_path - 4, ".png")) {
		img_tmp = CargarPGM(ruta_origen);
	}
	else if (strcmp(ruta_origen + len_path - 4, ".pgm")) {
		img_tmp = CargarPGM(ruta_origen);
	}
	else if (strcmp(ruta_origen + len_path - 4, ".dcm")) {
		img_tmp = CargarDICOM(ruta_origen, nivel);
	}
	else {
#if defined(_WIN32) || defined(_WIN64)
		sprintf_s(err_msg, 512, "<<Error: Image format not supported in: \'%s\'>>\n", ruta_origen);
#else
		sprintf(err_msg,  COLOR_BACK_RED "<<" COLOR_WHITE "Error: Image format not supported in: \'%s\'>>" COLOR_RESET "\n", ruta_origen);
#endif
		escribirLog(err_msg);
	}

	switch (img_idx) {
	case BASE:
		my_base = new IMGCONT(img_tmp->my_height, img_tmp->my_width);
		*my_base = *img_tmp;

		if (enmascarar) {
			my_mask = definirMask(my_base);
		}

		break;
	case GROUNDTRUTH:
		my_groundtruth = new IMGCONT(img_tmp->my_height, img_tmp->my_width);
		*my_groundtruth = *img_tmp;
		break;
	case MASK:
		my_mask = new IMGCONT(img_tmp->my_height, img_tmp->my_width);
		*my_mask = *img_tmp;
		break;
	case SKELETON:
		my_skeleton = new IMGCONT(img_tmp->my_height, img_tmp->my_width, 1, 1);
		*my_skeleton = *img_tmp;
		break;
	case SEGMENT:
		my_response = new IMGCONT(img_tmp->my_height, img_tmp->my_width);
		*my_response = *img_tmp;
		break;
	case MAPDIST:
		my_distmap = new IMGCONT(img_tmp->my_height, img_tmp->my_width);
		*my_distmap = *img_tmp;
		break;
	case THRESHOLD:
		my_segmented = new IMGCONT(img_tmp->my_height, img_tmp->my_width);
		*my_segmented = *img_tmp;
		break;
	case BORDERS:
		my_boundaries = new IMGCONT(img_tmp->my_height, img_tmp->my_width);
		*my_boundaries = *img_tmp;
		break;
	}

	delete img_tmp;

}




/*  Metodo: GuardarPGM */
void IMGVTK::GuardarPGM(IMGCONT *img_src, const char *ruta, const double min, const double max)
{
#if defined(_WIN32) || defined(_WIN64)
	FILE *img_file;
	fopen_s(&img_file, ruta, "w");
#else
	FILE *img_file = fopen(ruta, "w");
#endif

	fprintf(img_file, "P2\n");
	fprintf(img_file, "# by FerCer\n");

	fprintf(img_file, "%i %i\n", img_src->my_width, img_src->my_height);
	fprintf(img_file, "255\n");

	int intensity;
	for (unsigned int y = 0; y < img_src->my_height; y++) {
		for (unsigned int x = 0; x < img_src->my_width; x++) {
			intensity = (unsigned char)(255.0 * (*(img_src->my_img_data + (x + img_src->my_offset_x) + (y + img_src->my_offset_y)*img_src->my_width) - min) / (max - min));
			fprintf(img_file, "%i\n", intensity);
		}
	}

	fclose(img_file);
}




/*  Metodo: GuardarPNG */
void IMGVTK::GuardarPNG(IMGCONT *img_src, const char *ruta, const double min, const double max)
{	
}




/*  Metodo: Guardar
    Funcion: Guarda la imagen en la ruta especificada con la extension especificada.
*/
void IMGVTK::Guardar(IMG_IDX img_idx, const char *ruta, const TIPO_IMG tipo_salida )
{

    IMGCONT *img_tmp = NULL;
	
    switch( img_idx ){
    case BASE:
        img_tmp = my_base;
        break;
    case GROUNDTRUTH:
        img_tmp = my_groundtruth;
        break;
    case MASK:
        img_tmp = my_mask;
        break;
    case SKELETON:
        img_tmp = my_skeleton;
        break;
    case SEGMENT:
        img_tmp = my_response;
        break;
    case MAPDIST:
        img_tmp = my_distmap;
        break;
    case THRESHOLD:
        img_tmp = my_segmented;
        break;
    case BORDERS:
        img_tmp = my_boundaries;
        break;
    }

    double min = MY_INF;
    double max =-MY_INF;

	for (int y = 0; y < img_tmp->my_height; y++) {
		for (int x = 0; x < img_tmp->my_width; x++) {
			if (min > *(img_tmp->my_img_data + (x + img_tmp->my_offset_x) + (y + img_tmp->my_offset_y)*(img_tmp->my_width))) {
				min = *(img_tmp->my_img_data + (x + img_tmp->my_offset_x) + (y + img_tmp->my_offset_y)*(img_tmp->my_width));
			}
			if (max < *(img_tmp->my_img_data + (x + img_tmp->my_offset_x) + (y + img_tmp->my_offset_y)*(img_tmp->my_width))) {
				max = *(img_tmp->my_img_data + (x + img_tmp->my_offset_x) + (y + img_tmp->my_offset_y)*(img_tmp->my_width));
			}
		}
	}

    if( fabs(max - min) < 1e-12 ){
        max = 1.0;
        min = 0.0;
    }

    switch(tipo_salida){
	case PGM: {
		GuardarPGM(img_tmp, ruta, min, max);
		break;
	}
        case PNG:{

			break;
        }
    }
}






/* CONSTRUCTORES */
IMGVTK::IMGVTK(){
    max_dist = 0;

    my_base = NULL;
    my_groundtruth = NULL;
    my_skeleton = NULL;
    my_mask = NULL;
    my_distmap = NULL;
    my_boundaries = NULL;
    my_response = NULL;
    my_segmented = NULL;
	
    pix_caract = NULL;

    // Defaults:
    SID = 0.0;
    SOD = 0.0;
    DDP = SID - SOD;
    DISO = SID / 2;
    LAORAO = 0.0;
    CRACAU = 0.0;
    pixX = 1.0;//0.308;
    pixY = 1.0;//0.308;
    WCenter = 127.5;
    WWidth = 255.0;
}



IMGVTK::IMGVTK(const IMGVTK &origen) {
	max_dist = origen.max_dist;

	my_base = NULL;
	my_groundtruth = NULL;
	my_skeleton = NULL;
	my_mask = NULL;
	my_distmap = NULL;
	my_boundaries = NULL;
	my_response = NULL;
	my_segmented = NULL;

	pix_caract = NULL;

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

	if (origen.my_base) {
		*my_base = *origen.my_base;
	}

	if (origen.my_groundtruth) {
		*my_groundtruth = *origen.my_groundtruth;
	}

	if (origen.my_mask) {
		*my_mask = *origen.my_mask;
	}

	if (origen.my_response) {
		*my_response = *origen.my_response;
	}

	if (origen.my_skeleton) {
		*my_skeleton = *origen.my_skeleton;
	}

	if (origen.my_distmap) {
		*my_distmap = *origen.my_distmap;
	}

	if (origen.my_boundaries) {
		*my_boundaries = *origen.my_boundaries;
	}

	if (origen.my_segmented) {
		*my_segmented = *origen.my_segmented;
	}
}




IMGVTK::IMGVTK( const char *ruta_origen, const bool enmascarar, const int nivel){
    max_dist = 0;

	my_base = NULL;
	my_groundtruth = NULL;
	my_skeleton = NULL;
	my_mask = NULL;
	my_distmap = NULL;
	my_boundaries = NULL;
	my_response = NULL;
	my_segmented = NULL;

	pix_caract = NULL;

    // Defaults:
    SID = 0.0;
    SOD = 0.0;
    DDP = SID - SOD;
    DISO = SID / 2;
    LAORAO = 0.0;
    CRACAU = 0.0;
    pixX = 1.0;//0.308;
    pixY = 1.0;//0.308;
    WCenter = 127.5;
    WWidth = 255.0;

    cols = 0;
    rows = 0;
    max_dist = 0;

    DEB_MSG("Cargando imagen desde: " COLOR_GREEN << ruta_origen COLOR_NORMAL);

    Cargar(BASE, ruta_origen, enmascarar, nivel);
}




/* DESTRUCTOR */
IMGVTK::~IMGVTK() {
	if (pix_caract) {
		/// Liberar memoria recursivamente:
		borrarSkeleton(pix_caract);
		delete pix_caract;
	}


	if (my_base) {
		delete my_base;
	}

	if (my_groundtruth) {
		delete my_groundtruth;
	}

	if (my_mask) {
		delete my_mask;
	}

	if (my_response) {
		delete my_response;
	}

	if (my_skeleton) {
		delete my_skeleton;
	}

	if (my_distmap) {
		delete my_distmap;
	}

	if (my_boundaries) {
		delete my_boundaries;
	}

	if (my_segmented) {
		delete my_segmented;
	}
}




// O P E R A D O R E S  S O B R E C A R G A D O S
// El operador de copia extrae unicamente el contenido de la imagen original
IMGVTK& IMGVTK::operator= ( const IMGVTK &origen ){
    max_dist = origen.max_dist;

	my_base = NULL;
	my_groundtruth = NULL;
	my_skeleton = NULL;
	my_mask = NULL;
	my_distmap = NULL;
	my_boundaries = NULL;
	my_response = NULL;
	my_segmented = NULL;

	pix_caract = NULL;


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
	
	if (origen.my_base) {
		*my_base = *origen.my_base;
	}

	if (origen.my_groundtruth) {
		*my_groundtruth = *origen.my_groundtruth;
	}

	if (origen.my_mask) {
		*my_mask = *origen.my_mask;
	}

	if (origen.my_response) {
		*my_response = *origen.my_response;
	}

	if (origen.my_skeleton) {
		*my_skeleton = *origen.my_skeleton;
	}

	if (origen.my_distmap) {
		*my_distmap = *origen.my_distmap;
	}

	if (origen.my_boundaries) {
		*my_boundaries = *origen.my_boundaries;
	}

	if (origen.my_segmented) {
		*my_segmented = *origen.my_segmented;
	}

	return *this;
}





// M E T O D O S      P R I V A D O S
/*  Metodo: setRuta
    Funcion: Copia al argumento a una variable local.
*/
char* IMGVTK::setRuta( const char *ruta_input ){
    const int l_src = (int)strlen( ruta_input ) + 1;
    char *mi_ruta = new char [512];
#if defined(_WIN32) || defined(_WIN64)
	sprintf_s(mi_ruta, 512 * sizeof(char), "%s" , ruta_input);
#else
	strcpy(mi_ruta, ruta_input);
#endif
    return mi_ruta;
}

// C L A S E: IMGVTK  ------------------------------------------------------------------------------------------------- ^
