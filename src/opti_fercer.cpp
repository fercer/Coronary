/*******************************************************************************************************************
*    CENTRO DE INVESTIGACION EN MATEMATICAS
*    DOCTORADO EN CIENCIAS DE LA COMPUTACION
*
*    FERNANDO CERVANTES SANCHEZ
*    AGO - 2016
*******************************************************************************************************************/


#include "opti_fercer.h"





/*
    FUNCION:        shuffle
    DESCRIPCION:    Realiza la permutacion aleatoria de los indices desde 0 hasta n-1, por medio del reordenamiento (barajamiento).

    FUENTE: http://stackoverflow.com/questions/15961119/creating-random-permutation-without-repeating
*/
void shuffle(int *indices, const int n0, const int n1, STAUS *mi_semilla){
        // Se asignen inicialmente los indices.
        for(int i = n0; i < n1; i++){
                indices[i - n0] = i;
        }

        // Se reordenan aleatoriamente de atras para adelante.
        for(int i = n1-1; i >= n0; i--){
                const int pos = (int)(HybTaus(0.0, 1.0, mi_semilla) * (double)(i + 1 - n0));
                const int tmp = indices[i - n0];
                indices[i - n0] = indices[pos];
                indices[pos] = tmp;
        }
}






/*  Metodo:  genIndiv
 *   Funcion: Genera un individuo en base a la estructura y parametros para su generacion.
 */
void genIndiv( INDIV *indiv, NODO_ESTRUCT* nodo_padre, NODO_ESTRUCT* nodo_hijo, STAUS *mi_semilla ){
    if( !nodo_hijo ){
        return;
    }

    if( *((char*)indiv->my_features + *(indiv->ptr_starts + nodo_padre->varID))){
        *((char*)indiv->my_features + *(indiv->ptr_starts + nodo_hijo->varID)) = (HybTaus(0.0, 1.0, mi_semilla) <= nodo_hijo->prob_cond_1) ? 0 : 1;
    }else{
        *((char*)indiv->my_features + *(indiv->ptr_starts + nodo_hijo->varID)) = (HybTaus(0.0, 1.0, mi_semilla) <= nodo_hijo->prob_cond_0) ? 0 : 1;
    }

    for( int k = 0; k < nodo_hijo->n_hijos; k++){
        genIndiv( indiv, nodo_hijo, nodo_hijo->nodos_hijos[k], mi_semilla );
    }
}







/*	Metodo:  genPob
 *   Funcion: Genera una poblacion dada una estructura y parametros para su generacion.
*/
void genPob(INDIV *pob, const int n_pob, NODO_ESTRUCT* estructura, FUNPTR mi_fun, FUNPARS *mis_pars, STAUS *mi_semilla){
    for( int i = 0; i < n_pob; i++ ){
        // Calcular el valor de la variable raiz definida por la esturctura:
        *((char*)(pob + i)->my_features + *( (pob + i)->ptr_starts + estructura->varID) ) = (HybTaus(0.0, 1.0, mi_semilla) <= estructura->prob_cond_0) ? 0 : 1;

        for(int k = 0; k < estructura->n_hijos; k++){
            genIndiv( pob + i, estructura, estructura->nodos_hijos[k], mi_semilla);
        }

        (pob + i)->fitness = mi_fun( mis_pars, pob + i );
    }
}





/*  Metodo:  genIndivCont
 *   Funcion: Genera un individuo en base a la estructura y parametros para su generacion.
 */
void genIndivCont( INDIV *indiv, NODO_ESTRUCT* nodo_padre, NODO_ESTRUCT* nodo_hijo, FUNPARS *mis_pars, STAUS *mi_semilla ){
    if( !nodo_hijo ){
        return;
    }

    // Generar la varianza condicional:
    const double var_cond = nodo_hijo->prob_cond_1 - (nodo_hijo->covar*nodo_hijo->covar) / nodo_padre->prob_cond_1;
    const double med_cond = nodo_hijo->prob_cond_0 + nodo_hijo->covar / nodo_padre->prob_cond_1 * ( *(double*)((char*)indiv->my_features + *(indiv->ptr_starts + nodo_padre->varID)) - nodo_padre->prob_cond_0 );

    double min = *(double*)((char*)mis_pars->parameters + *(mis_pars->ptr_starts + 2*nodo_hijo->varID  ));
    double max = *(double*)((char*)mis_pars->parameters + *(mis_pars->ptr_starts + 2*nodo_hijo->varID+1));

    // Generar un valor dentro de la region debusqueda:
    double valor;

    if( var_cond > 1e-12 ){
        valor = anorm( med_cond, var_cond, mi_semilla );
        valor = (min > valor) ? min : (valor > max) ? max : valor;
    }else{
        valor = med_cond;
    }

    *(double*)((char*)indiv->my_features + *(indiv->ptr_starts + nodo_hijo->varID)) = valor;
    for( int k = 0; k < nodo_hijo->n_hijos; k++){
        genIndivCont( indiv, nodo_hijo, nodo_hijo->nodos_hijos[k], mis_pars, mi_semilla );
    }
}





/*	Metodo:  genPobCont
 *   Funcion: Genera una poblacion dada una estructura y parametros para su generacion con variables continuas.
*/
void genPobCont(INDIV *pob, const int n_pob, NODO_ESTRUCT* estructura, FUNPTR mi_fun, FUNPARS *mis_pars, STAUS *mi_semilla){
    for( int i = 0; i < n_pob; i++ ){

        double min = *(double*)((char*)mis_pars->parameters + *(mis_pars->ptr_starts  ));
        double max = *(double*)((char*)mis_pars->parameters + *(mis_pars->ptr_starts+1));

        // Calcular el valor de la variable raiz definida por la esturctura:
        double valor;
        if( estructura->prob_cond_1 > 1e-12 ){
            valor = anorm( estructura->prob_cond_0, estructura->prob_cond_1, mi_semilla );
            valor = (min > valor) ? min : (valor > max) ? max : valor;
        }else{
            valor = estructura->prob_cond_0;
        }
        *((char*)(pob + i)->my_features + *( (pob + i)->ptr_starts + estructura->varID) ) = valor;

        for(int k = 0; k < estructura->n_hijos; k++){
            genIndivCont( pob + i, estructura, estructura->nodos_hijos[k], mis_pars, mi_semilla);
        }

        (pob + i)->fitness = mi_fun( mis_pars, pob + i );
    }
}



/*	Metodo:  liberarNodos
    Funcion: Libera la memoria solicitada para la estructura de datos empleada por el MIMIC y CHOWLIU.
*/
void liberarNodos( NODO_ESTRUCT *estructura ){
    if( estructura->n_hijos > 0 ){
        for( int i = 0; i < estructura->n_hijos; i++){
            liberarNodos( *(estructura->nodos_hijos + i) );
        }
        free( estructura->nodos_hijos );
    }
    free( estructura );
}





/*	Metodo:  calcFrecs
    Funcion: Cuenta el numero de 1's para cada variable y sus condicionales.
*/
void calcFrecs( INDIV *pob, const int n_pob, const int n_vars, unsigned int **frecuencias ){
    for(int i = 0; i < n_vars; i++){
        memset(*(frecuencias + i), 0, (n_vars-i)*sizeof(unsigned int));
    }
    //Se realiza un conteo de los valores de la poblacion
    for( int i = 0; i < n_pob; i++){
        for( int j = 0; j < n_vars; j++){
            if( *((char*)(pob + i)->my_features + *((pob + i)->ptr_starts + j)) ){
                for( int k = 0; k < (n_vars-j); k++){
                    frecuencias[j][k] += (unsigned int) *((char*)(pob + i)->my_features + *((pob + i)->ptr_starts+k+j));
                }
            }
        }
    }
}




/*	Metodo:  calcInfoMutua
    Funcion: Calcula las entropias marginales y condicionales en base a la matriz de frecuencias.
*/
void calcInfoMutua(double **info_mat, const int n_pob, const int n_vars, unsigned int **frecuencias ){
    // Calcular las entropias marginales a partir de los datos.
    for(int i = 0; i < n_vars; i++){
        const double P1 = (double)frecuencias[i][0] / (double)n_pob;
        const double P0 = 1.0 - P1;
        // Revisar P(x) = 0, para evitar indeterminaciones con 0.log(0), utilizando lim x->0   x * log(x) = 0:
        double logP0 = 0.0;
        if( P0 > 1e-12 ){
            logP0 = log( P0 );
        }
        double logP1 = 0.0;
        if( P1 > 1e-12 ){
            logP1 = log( P1 );
        }
        info_mat[i][0] = - (P0 * logP0 + P1 * logP1);
    }

    for( int i = 0; i < (n_vars-1); i++){
        // Frecuencia marginal de la variable A = 1:
        const double A1_M = (double)frecuencias[i][0] / (double)n_pob;
        
        for( int j = 1; j < (n_vars-i); j++){
            // Frecuencia de A = 1 y B = 1
            const double A1_B1 = (double)frecuencias[i][j] / (double)n_pob;
            
            // Frecuencia de A = 1 y B = 0:
            const double A1_B0 = A1_M - A1_B1;
            
            // Frecuencia marginal de la variable B = 1:
            const double B1_M = (double)frecuencias[i+j][0] / (double)n_pob;
            
            // Frecuencia marginal de la variable B = 0:
            const double B0_M = 1.0 - B1_M;
            
            // Frecuencia de A = 0 y B = 1:
            const double A0_B1 = B1_M - A1_B1;
            
            // Frecuencia de A = 0 y B = 0:
            const double A0_B0 = B0_M - A1_B0;
            

            // Revisar P(x) = 0, para evitar indeterminaciones con 0.log(0), utilizando lim x->0   x * log(x) = 0:
            double logA0B0 = 0.0;
            if( A0_B0 > 1e-12 ){
                logA0B0 = log( A0_B0 );
            }
            double logA0B1 = 0.0;
            if( A0_B1 > 1e-12 ){
                log( A0_B1 );
            }
            double logA1B0 = 0.0;
            if( A1_B0 > 1e-12 ){
                logA1B0 = log( A1_B0 );
            }
            double logA1B1 = 0.0;
            if( A1_B1 > 1e-12 ){
                logA1B1 = log( A1_B1 );
            }

            // calcular la informacion muta entre A y B:
            info_mat[i][j] = info_mat[i][0] + info_mat[i+j][0] + A0_B0*logA0B0 + A0_B1*logA0B1 + A1_B0*logA1B0 + A1_B1*logA1B1;
        }
    }
}










/*	Metodo:  calcCovars
    Funcion: Estima los parametros de media, varianza y covarianza de las distribuciones de las variables.
*/
void calcCovars( INDIV *pob, const int n_pob, const int n_vars, double **covars_mat, double *medias){

    for(int i = 0; i < n_vars; i++){
        memset(*(covars_mat + i), 0, (n_vars-i)*sizeof(double));
        *(medias + i) = 0.0;
    }

    // Se estiman las medias a partir de la poblacion:
    for( int i = 0; i < n_pob; i++){
        for( int j = 0; j < n_vars; j++){
            const double x = *(double*)((char*)(pob + i)->my_features + *(pob->ptr_starts + j));
            medias[j] += x / (double)n_pob;
        }
    }

    // Estimar las covarianzas:
    for( int i = 0; i < n_pob; i++){
        for( int j  = 0; j < n_vars; j++){
            const double x = *(double*)((char*)(pob + i)->my_features + *(pob->ptr_starts + j));
            for( int k = j; k < n_vars; k++){
                const double y = *(double*)((char*)(pob + i)->my_features + *(pob->ptr_starts + k));
                covars_mat[j][k-j] += (x - *(medias+j)) * (y - *(medias+k)) / (double)n_pob;
            }
        }
    }

}




/*	Metodo:  calcInfoMutuaCont
    Funcion: Calcula las entropias marginales y condicionales en base a la informacion de las distribuciones continuas de las variables.
*/
void calcInfoMutuaCont(double **info_mat, const int n_vars, double **covars_mat){

    // Calcular las entropias marginales a partir de los datos.
    for(int i = 0; i < n_vars; i++){
        info_mat[i][0] = log( 2.0 * Mi_PI * Mi_EXP * sqrt( covars_mat[i][0] ) );
    }

    for( int i = 0; i < (n_vars-1); i++){
        const double desv_est_1 = sqrt( covars_mat[i][0] );

        for( int j = 1; j < (n_vars-i); j++){
            const double desv_est_2 = sqrt( covars_mat[j][0] );

            const double corr_1_2 = covars_mat[i][j] / (desv_est_1 * desv_est_2);

            info_mat[i][j] = -log( 1.0 - corr_1_2*corr_1_2 )/ 2.0 ;
        }
    }
}





/*	Metodo:  estructMIMIC
    Funcion: Genera una estructura de cadena tipo MIMIC a partir de las entropias condicionales obtenidas de los datos.
*/
NODO_ESTRUCT *estructMIMIC(const int n_pob, const int n_vars, unsigned int **frecuencias, double **info_mat){

    NODO_ESTRUCT *estructura = (NODO_ESTRUCT*) malloc(sizeof(NODO_ESTRUCT));
    
    unsigned char *vars_asignadas = (unsigned char*) calloc(n_vars, sizeof(unsigned char));
    int varID = 0;
    double min_entro_marg = info_mat[0][0];
    for( int i = 1; i < n_vars; i++){
        if( min_entro_marg > info_mat[i][0] ){
            varID = i;
            min_entro_marg = info_mat[i][0];
        }
    }
    estructura->varID = varID;
    estructura->prob_cond_0 = 1.0 - (double)frecuencias[varID][0] / (double)n_pob;
    vars_asignadas[ varID ] = 1;

    estructura->entropia = min_entro_marg;
    estructura->info = 0;

    // Asignar los siguientes nodos dependiendo de la informacion mutua:
    NODO_ESTRUCT *nodo_padre = estructura;
    int n_asignados = 1;

    while( n_asignados < n_vars ){
        nodo_padre->n_hijos = 1;
        nodo_padre->nodos_hijos = (NODO_ESTRUCT**) malloc(sizeof(NODO_ESTRUCT*));
        nodo_padre->nodos_hijos[0] = (NODO_ESTRUCT*) malloc(sizeof(NODO_ESTRUCT));

        double max_info_mutua =-1e100;
        for( int i = 0; i < n_vars; i++){
            if( !vars_asignadas[i] ){
                int A =  (nodo_padre->varID < i) ? nodo_padre->varID : i;
                int B = ((nodo_padre->varID < i) ? i : nodo_padre->varID ) - A;
                
                const double info_mutua = info_mat[A][B];

                if( max_info_mutua < info_mutua ){
                    max_info_mutua = info_mutua;
                    varID = i;
                }
            }
        }

        nodo_padre->nodos_hijos[0]->varID = varID;
        nodo_padre->nodos_hijos[0]->info = max_info_mutua;
        nodo_padre->nodos_hijos[0]->entropia = info_mat[varID][0] - max_info_mutua; // Entropia condicional
        vars_asignadas[varID] = 1;
        n_asignados++;

        // Calcular las probabilidades condicionales a partir de los datos:
        int A =  (nodo_padre->varID < varID) ? nodo_padre->varID : varID;
        int B = ((nodo_padre->varID < varID) ? varID : nodo_padre->varID) - A;
        
        // Frecuencia marginal de la variable A = 1:
        const double A1_M = (double)frecuencias[A][0];
        
        // Frecuencia marginal de la variable A = 0:
        const double A0_M = (double)n_pob - A1_M;
        
        // Frecuencia de A = 1 y B = 1
        const double A1_B1 = (double)frecuencias[A][B];
        
        // Frecuencia de A = 1 y B = 0:
        const double A1_B0 = A1_M - A1_B1;
        
        // Frecuencia marginal de la variable B = 1:
        const double B1_M = (double)frecuencias[B+A][0];
        
        // Frecuencia marginal de la variable B = 0:
        const double B0_M = (double)n_pob - B1_M;

        // Frecuencia de A = 0 y B = 0:
        const double A0_B0 = B0_M - A1_B0;
        
        nodo_padre->nodos_hijos[0]->prob_cond_0 = A0_B0 / A0_M;
        nodo_padre->nodos_hijos[0]->prob_cond_1 = A1_B0 / A1_M;

        nodo_padre = nodo_padre->nodos_hijos[0];
    }
    free(vars_asignadas);

    nodo_padre->nodos_hijos = NULL;
    nodo_padre->n_hijos = 0;
    
    return estructura;
}





/*	Metodo:  estructMIMICCont
    Funcion: Genera una estructura de cadena tipo MIMIC a partir de las entropias condicionales obtenidas de los datos para el caso continuo.
*/
NODO_ESTRUCT *estructMIMICCont(const int n_vars, double **covars_mat, double *medias, double **info_mat){

    NODO_ESTRUCT *estructura = (NODO_ESTRUCT*) malloc(sizeof(NODO_ESTRUCT));

    unsigned char *vars_asignadas = (unsigned char*) calloc(n_vars, sizeof(unsigned char));
    int varID = 0;
    double min_entro_marg = info_mat[0][0];
    for( int i = 1; i < n_vars; i++){
        if( min_entro_marg > info_mat[i][0] ){
            varID = i;
            min_entro_marg = info_mat[i][0];
        }
    }
    estructura->varID = varID;
    // Media de la variable raiz:
    estructura->prob_cond_0 = medias[ varID ];
    // Varianza de la variable raiz:
    estructura->prob_cond_1 = covars_mat[ varID ][ 0 ];
    vars_asignadas[ varID ] = 1;

    estructura->entropia = min_entro_marg;
    estructura->info = 0;

    // Asignar los siguientes nodos dependiendo de la informacion mutua:
    NODO_ESTRUCT *nodo_padre = estructura;
    int n_asignados = 1;

    while( n_asignados < n_vars ){
        nodo_padre->n_hijos = 1;
        nodo_padre->nodos_hijos = (NODO_ESTRUCT**) malloc(sizeof(NODO_ESTRUCT*));
        nodo_padre->nodos_hijos[0] = (NODO_ESTRUCT*) malloc(sizeof(NODO_ESTRUCT));

        double max_info_mutua =-1e100;
        for( int i = 0; i < n_vars; i++){
            if( !vars_asignadas[i] ){
                int A =  (nodo_padre->varID < i) ? nodo_padre->varID : i;
                int B = ((nodo_padre->varID < i) ? i : nodo_padre->varID ) - A;

                const double info_mutua = info_mat[A][B];

                if( max_info_mutua < info_mutua ){
                    max_info_mutua = info_mutua;
                    varID = i;
                }
            }
        }

        nodo_padre->nodos_hijos[0]->varID = varID;
        nodo_padre->nodos_hijos[0]->info = max_info_mutua;
        nodo_padre->nodos_hijos[0]->entropia = info_mat[varID][0] - max_info_mutua; // Entropia condicional
        vars_asignadas[varID] = 1;
        n_asignados++;


        const int A =  (nodo_padre->varID < varID) ? nodo_padre->varID : varID;
        const int B = ((nodo_padre->varID < varID) ? varID : nodo_padre->varID ) - A;

        //Asignar media, varianza y covarianza:
        nodo_padre->nodos_hijos[0]->prob_cond_0 = medias[ varID ];
        nodo_padre->nodos_hijos[0]->prob_cond_1 = covars_mat[ varID ][ 0 ];
        nodo_padre->nodos_hijos[0]->covar = covars_mat[ A ][ B ];

        nodo_padre = nodo_padre->nodos_hijos[0];
    }
    free(vars_asignadas);

    nodo_padre->nodos_hijos = NULL;
    nodo_padre->n_hijos = 0;

    return estructura;
}





/*  Metodo:  estructCHOWLIU
*   Funcion: Genera una estructura de arbol de dependencia tipo CHOWLIU a partir de las entropias condicionales obtenidas de los datos.
*/
NODO_ESTRUCT *estructCHOWLIU(const int n_pob, const int n_vars, unsigned int **frecuencias, double **info_mat){
    // Definir la raiz como el nodo on menor entropia:
    int varID = 0;
    double min_entro_marg = info_mat[0][0];
    for( int i = 1; i < n_vars; i++){
        if( min_entro_marg > info_mat[i][0] ){
            varID = i;
            min_entro_marg = info_mat[i][0];
        }
    }

    // Asignar los siguientes nodos dependiendo de la informacion mutua:
    //// Algoritmo Prim:
    double *informaciones = (double*) malloc( n_vars * sizeof(double) );

    // El arreglo padres indica cual es el padre de cada variable (nodo):
    int *padres = (int*) malloc(n_vars * sizeof(int));
    unsigned char *vars_asignadas = (unsigned char*) calloc(n_vars, sizeof(unsigned char));
    for(int k = 0; k < n_vars; k++){
        informaciones[k] = -2;
        padres[k] = -1;
    }

    informaciones[ varID ] = 0.0;
    int n_asignados = 0;
    int var_actual = 0;
    while( n_asignados < n_vars ){
        double max_info_mutua = -1.0;
        
        for(int k = 0; k < n_vars; k++){
            if( vars_asignadas[k] ){
                continue;
            }
            if( max_info_mutua < informaciones[k] ){
                var_actual = k;
                max_info_mutua = informaciones[k];
            }
        }

        // Verificar las variables que dada la variable elejida, su informacion mutua se incrementa:
        for( int k = 0; k < var_actual; k++){
            if( vars_asignadas[ k ] ){
                continue;
            }
            if( informaciones[k] < info_mat[k][var_actual-k] ){
                informaciones[k] = info_mat[k][var_actual-k];
                padres[k] = var_actual;
            }
        }
        for( int k = (var_actual+1); k < n_vars; k++){
            if( vars_asignadas[ k ] ){
                continue;
            }
            if( informaciones[k] < info_mat[var_actual][k-var_actual] ){
                informaciones[k] = info_mat[var_actual][k-var_actual];
                padres[k] = var_actual;
            }
        }

        vars_asignadas[ var_actual ] = 1;
        n_asignados++;
    }
        
    free( vars_asignadas );
    free( informaciones );

    // Generar todos los nodos dentro de un arreglo temporal:
    NODO_ESTRUCT **nodos_temp = (NODO_ESTRUCT**) malloc( n_vars * sizeof(NODO_ESTRUCT*) );
    for( int i = 0; i < n_vars; i++){
        nodos_temp[i] = (NODO_ESTRUCT*) malloc( sizeof(NODO_ESTRUCT) );
        nodos_temp[i]->n_hijos = 0;
    }
    
    // Crear las dependencias segun el vector que indica cuantos hijos tiene cada uno de los nodos:
    for( int i = 0; i < n_vars; i++){
        if( padres[i] >= 0 ){ // Si padres[i] == -1 es el nodo raiz, de otra forma el nodo es hijo de alguno otro:
            nodos_temp[ padres[i] ]->n_hijos++;
        }
    }
    
    // Se requiere memoria para cada nodo para el numero de hijos que tiene:
    for( int i = 0; i < n_vars; i++){
        const int n_hijos = nodos_temp[ i ]->n_hijos;
        if( n_hijos > 0 ){
            nodos_temp[ i ]->nodos_hijos = (NODO_ESTRUCT**) calloc( n_hijos, sizeof(NODO_ESTRUCT*) );
            nodos_temp[ i ]->n_hijos = 0;
        }else{
            nodos_temp[ i ]->nodos_hijos = NULL;
        }
    }
    NODO_ESTRUCT *raiz = NULL;
    for( int i = 0; i < n_vars; i++){
        if( padres[i] >= 0 ){
            nodos_temp[ padres[i] ]->nodos_hijos[ nodos_temp[ padres[i] ]->n_hijos ] = nodos_temp[ i ];
            nodos_temp[ padres[i] ]->n_hijos++;
            nodos_temp[ i ]->varID = i;
            
            // Calcular las probabilidades condicionales a partir de los datos:
            int A =  (padres[i] < i) ? padres[i] : i;
            int B = ((padres[i] < i) ? i : padres[i]) - A;
            
            // Frecuencia marginal de la variable A = 1:
            const double A1_M = (double)frecuencias[A][0];
            
            // Frecuencia marginal de la variable A = 0:
            const double A0_M = (double)n_pob - A1_M;
            
            // Frecuencia de A = 1 y B = 1
            const double A1_B1 = (double)frecuencias[A][B];
            
            // Frecuencia de A = 1 y B = 0:
            const double A1_B0 = A1_M - A1_B1;
            
            // Frecuencia marginal de la variable B = 1:
            const double B1_M = (double)frecuencias[B+A][0];
            
            // Frecuencia marginal de la variable B = 0:
            const double B0_M = (double)n_pob - B1_M;
            
            // Frecuencia de A = 0 y B = 1:
            const double A0_B1 = B1_M - A1_B1;
            
            // Frecuencia de A = 0 y B = 0:
            const double A0_B0 = B0_M - A1_B0;
            
            nodos_temp[ i ]->prob_cond_0 = A0_B0 / A0_M;
            nodos_temp[ i ]->prob_cond_1 = A1_B0 / A1_M;
            
            nodos_temp[ i ]->info = info_mat[A][B];
            nodos_temp[ i ]->entropia = info_mat[i][0] - info_mat[A][B]; // Entropia condicional
            
        }else{
            // Calcular la probabilidad marginal para el nodo raiz:
            nodos_temp[ i ]->varID = i;
            nodos_temp[ i ]->prob_cond_0 = 1.0 - (double)frecuencias[i][0] / (double)n_pob;
            
            nodos_temp[ i ]->entropia = info_mat[i][0];
            nodos_temp[ i ]->info = 0;
            
            raiz = nodos_temp[ i ];
        }
    }

    free( padres );
    free( nodos_temp );
    return raiz;
}




/*  Metodo:  estructCHOWLIUCont
*   Funcion: Genera una estructura de arbol de dependencia tipo CHOWLIU a partir de las entropias condicionales obtenidas de los datos para el caso continuo.
*/
NODO_ESTRUCT *estructCHOWLIUCont(const int n_vars, double **covars_mat, double *medias, double **info_mat){
    // Definir la raiz como el nodo on menor entropia:
    int varID = 0;
    double min_entro_marg = info_mat[0][0];
    for( int i = 1; i < n_vars; i++){
        if( min_entro_marg > info_mat[i][0] ){
            varID = i;
            min_entro_marg = info_mat[i][0];
        }
    }

    // Asignar los siguientes nodos dependiendo de la informacion mutua:
    //// Algoritmo Prim:
    double *informaciones = (double*) malloc( n_vars * sizeof(double) );

    // El arreglo padres indica cual es el padre de cada variable (nodo):
    int *padres = (int*) malloc(n_vars * sizeof(int));
    unsigned char *vars_asignadas = (unsigned char*) calloc(n_vars, sizeof(unsigned char));
    for(int k = 0; k < n_vars; k++){
        informaciones[k] = -2;
        padres[k] = -1;
    }

    informaciones[ varID ] = 0.0;
    int n_asignados = 0;
    int var_actual = 0;
    while( n_asignados < n_vars ){
        double max_info_mutua = -1.0;

        for(int k = 0; k < n_vars; k++){
            if( vars_asignadas[k] ){
                continue;
            }
            if( max_info_mutua < informaciones[k] ){
                var_actual = k;
                max_info_mutua = informaciones[k];
            }
        }

        // Verificar las variables que dada la variable elejida, su informacion mutua se incrementa:
        for( int k = 0; k < var_actual; k++){
            if( vars_asignadas[ k ] ){
                continue;
            }
            if( informaciones[k] < info_mat[k][var_actual-k] ){
                informaciones[k] = info_mat[k][var_actual-k];
                padres[k] = var_actual;
            }
        }
        for( int k = (var_actual+1); k < n_vars; k++){
            if( vars_asignadas[ k ] ){
                continue;
            }
            if( informaciones[k] < info_mat[var_actual][k-var_actual] ){
                informaciones[k] = info_mat[var_actual][k-var_actual];
                padres[k] = var_actual;
            }
        }

        vars_asignadas[ var_actual ] = 1;
        n_asignados++;
    }

    free( vars_asignadas );
    free( informaciones );

    // Generar todos los nodos dentro de un arreglo temporal:
    NODO_ESTRUCT **nodos_temp = (NODO_ESTRUCT**) malloc( n_vars * sizeof(NODO_ESTRUCT*) );
    for( int i = 0; i < n_vars; i++){
        nodos_temp[i] = (NODO_ESTRUCT*) malloc( sizeof(NODO_ESTRUCT) );
        nodos_temp[i]->n_hijos = 0;
    }

    // Crear las dependencias segun el vector que indica cuantos hijos tiene cada uno de los nodos:
    for( int i = 0; i < n_vars; i++){
        if( padres[i] >= 0 ){ // Si padres[i] == -1 es el nodo raiz, de otra forma el nodo es hijo de alguno otro:
            nodos_temp[ padres[i] ]->n_hijos++;
        }
    }

    // Se requiere memoria para cada nodo para el numero de hijos que tiene:
    for( int i = 0; i < n_vars; i++){
        const int n_hijos = nodos_temp[ i ]->n_hijos;
        if( n_hijos > 0 ){
            nodos_temp[ i ]->nodos_hijos = (NODO_ESTRUCT**) calloc( n_hijos, sizeof(NODO_ESTRUCT*) );
            nodos_temp[ i ]->n_hijos = 0;
        }else{
            nodos_temp[ i ]->nodos_hijos = NULL;
        }
    }
    NODO_ESTRUCT *raiz = NULL;
    for( int i = 0; i < n_vars; i++){
        if( padres[i] >= 0 ){
            nodos_temp[ padres[i] ]->nodos_hijos[ nodos_temp[ padres[i] ]->n_hijos ] = nodos_temp[ i ];
            nodos_temp[ padres[i] ]->n_hijos++;
            nodos_temp[ i ]->varID = i;

            // Asignar las medias y varianzas correspondientes a cada nodo:
            int A =  (padres[i] < i) ? padres[i] : i;
            int B = ((padres[i] < i) ? i : padres[i]) - A;

            nodos_temp[ i ]->prob_cond_0 = medias[ i ];
            nodos_temp[ i ]->prob_cond_1 = covars_mat[ i ][ 0 ];
            nodos_temp[ i ]->covar = covars_mat[ A ][ B ];

            nodos_temp[ i ]->info = info_mat[A][B];
            nodos_temp[ i ]->entropia = info_mat[i][0] - info_mat[A][B]; // Entropia condicional

        }else{
            // Calcular la probabilidad marginal para el nodo raiz:
            nodos_temp[ i ]->varID = i;
            nodos_temp[ i ]->prob_cond_0 = medias[ i ];
            nodos_temp[ i ]->prob_cond_1 = covars_mat[ i ][ 0 ];

            nodos_temp[ i ]->entropia = info_mat[i][0];
            nodos_temp[ i ]->info = 0;

            raiz = nodos_temp[ i ];
        }
    }

    free( padres );
    free( nodos_temp );
    return raiz;
}





/*  Metodo:  entropiaModelo
 *   Funcion: Calcula la entropia condicional del modelo.
 */
double entropiaModelo(NODO_ESTRUCT *raiz){
    if( !raiz ){
        return 0.0;
    }
    
    double entropia = 0.0;
    for( int i = 0; i < raiz->n_hijos; i++){
        entropia += entropiaModelo( raiz->nodos_hijos[i] );
    }
    
    return entropia + raiz->entropia;
}


/*  Metodo:  infomutuaModelo
 *   Funcion: Calcula la informacion mutua del modelo.
 */
double infomutuaModelo(NODO_ESTRUCT *raiz){
    if( !raiz ){
        return 0.0;
    }
    
    double info_mutua = 0.0;
    for( int i = 0; i < raiz->n_hijos; i++){
        info_mutua += infomutuaModelo( raiz->nodos_hijos[i] );
    }
    
    return info_mutua + raiz->info;
}



/*  Metodo:  mostrarModelo
 *   Funcion: Muestra las dependencias del modelo (las aristas del arbol/cadena).
 */
void mostrarModelo( NODO_ESTRUCT *nodo_padre, NODO_ESTRUCT *nodo_hijo ){
    if( !nodo_hijo ){
        return;
    }

    if( !nodo_padre ){
        printf("Raiz[%i]\n", nodo_hijo->varID);
    }else{
        printf("Padre[%i]->Hijo[%i]\n", nodo_padre->varID, nodo_hijo->varID);
    }
    
    for( int i = 0; i < nodo_hijo->n_hijos; i++){
        mostrarModelo(nodo_hijo, nodo_hijo->nodos_hijos[i]);
    }
    
}



/*  Metodo:  optiMIMIC
 *   Funcion: Encuentra el optimo global de la funcion 'mi_fun' mediante el algoritmo MIMIC de optimizacion.
 */
INDIV *optiMIMIC(FUNPARS *mis_optipars, FUNPTR mi_fun, FUNPARS *mis_pars, STAUS *mi_semilla)
{

    const unsigned int n_iters = *(unsigned int*)((char*)mis_optipars->parameters + *(mis_optipars->ptr_starts));
    const unsigned int n_pob =   *(unsigned int*)((char*)mis_optipars->parameters + *(mis_optipars->ptr_starts + 1));
    const unsigned int n_seleccion = *(unsigned int*)((char*)mis_optipars->parameters + *(mis_optipars->ptr_starts + 2));
    const unsigned int n_vars = *(unsigned int*)((char*)mis_optipars->parameters + *(mis_optipars->ptr_starts + 3));


    //// Requerir memoria para las tablas de frecuencias e informacion mutua:
    unsigned int **frecuencias = (unsigned int**) malloc( n_vars * sizeof(unsigned int*) );
    double **info_mat = (double**) malloc( n_vars * sizeof(double*) );

    for( unsigned int i = 0; i < n_vars; i++){
        frecuencias[i] = (unsigned int*) malloc( (n_vars-i)* sizeof(unsigned int) );
        info_mat[i] =  (double*) malloc( (n_vars-i)* sizeof(double) );
    }

    //// Requerir memoria e inicializar la poblacion 0:
    INDIV *poblacion = (INDIV *) malloc( n_pob * sizeof( INDIV ) );
    for( unsigned int i = 0; i < n_pob; i++){
        (poblacion + i)->dim = n_vars;
        (poblacion + i)->ptr_starts = (size_t*) malloc((n_vars+1) * sizeof(size_t));
        *( (poblacion + i)->ptr_starts ) = 0;
        for(unsigned int k = 1; k <= n_vars; k++){
            *( (poblacion + i)->ptr_starts + k ) = *( (poblacion + i)->ptr_starts + k - 1) + sizeof(char);
        }
        (poblacion + i)->my_features = malloc( *((poblacion + i)->ptr_starts + n_vars) );

        // Generar la cadena de cada individuo de manera aleatoria.
        for(unsigned int k = 0; k < n_vars; k++){
            *((char*)(poblacion + i)->my_features + *((poblacion + i)->ptr_starts + k) ) = HybTaus( 0.0, 1.0, mi_semilla) <= 0.5 ? 0 : 1;
        }

        // Evaluar individuo en la funcion objetivo:
        (poblacion + i)->fitness = mi_fun( mis_pars, poblacion + i );
    }

    // Ordernar la poblacion segun sus fitness
    qsort((void*)poblacion, n_pob, sizeof(INDIV), compIndiv); // El mejor fitness queda en la posicion 0

    // Proponer la solucion al problema de optimizacion como el individuo con mejor fitness
    INDIV *mi_elite = (INDIV*) malloc( sizeof(INDIV) );
    mi_elite->dim = n_vars;
    mi_elite->ptr_starts = (size_t*) malloc( (n_vars+1)*sizeof(size_t) );
    mi_elite->my_features = malloc( *(poblacion->ptr_starts + n_vars) );

    mi_elite->fitness = poblacion->fitness;
    memcpy(mi_elite->ptr_starts, poblacion->ptr_starts, (n_vars+1)*sizeof(size_t));
    memcpy(mi_elite->my_features, poblacion->my_features, *(poblacion->ptr_starts + n_vars));

    /***************************************************************************************
     * Estimar las estructuras y parametros a partir de la poblacion '0':
     ***************************************************************************************/
    // Obtener las frecuencias de aparicion de 1's para cada variable y su condicional:
    calcFrecs( poblacion, n_seleccion, n_vars, frecuencias );

    // Calcular las Entropias marginales para todas las varibales y determinar el minimo como el nodo raiz de la cadena.
    calcInfoMutua(info_mat, n_seleccion, n_vars, frecuencias );

    // Generar el modelo a partir de la informacion obtenida dela poblacion en el momento '0':
    NODO_ESTRUCT *estructura_MIMIC = estructMIMIC(n_seleccion, n_vars, frecuencias, info_mat);

    printf("generacion\tevaluaciones\tevaluacion_elite\n");
    printf("%i\t%i\t%f\n", 0, n_pob, mi_elite->fitness);
    for(unsigned int t = 1; t < n_iters; t++){

        /***************************************************************************************
         * Generar una poblacion 't':
         ***************************************************************************************/
        genPob( poblacion, n_pob, estructura_MIMIC, mi_fun, mis_pars, mi_semilla );

        // Ordernar la poblacion segun sus fitness
        qsort((void*)poblacion, n_pob, sizeof(INDIV), compIndiv); // El mejor fitness queda en la posicion 0

        // Si un nuevo candidato mejora la solucion, se actualiza la solucion propuesta.
        if( mi_elite->fitness < poblacion->fitness ){
            mi_elite->fitness = poblacion->fitness;
            memcpy(mi_elite->my_features, poblacion->my_features, *(poblacion->ptr_starts + n_vars));
        }
        printf("%i\t%i\t%f\n", t, n_pob * (t+1), mi_elite->fitness);

        /***************************************************************************************
         * Estimar las estructuras y parametros a partir de la poblacion 't':
         ***************************************************************************************/
        // Obtener las frecuencias de aparicion de 1's para cada variable y su condicional:
        calcFrecs( poblacion, n_seleccion, n_vars, frecuencias );

        // Calcular las Entropias marginales para todas las varibales y determinar el minimo como el nodo raiz de la cadena.
        calcInfoMutua(info_mat, n_seleccion, n_vars, frecuencias );

        // Generar el modelo a partir de la informacion obtenida dela poblacion en el momento 't':
        liberarNodos( estructura_MIMIC );
        estructura_MIMIC = estructMIMIC(n_seleccion, n_vars, frecuencias, info_mat);
    }

    printf("%i\t%i\t%f\n", n_iters, n_pob * (n_iters+1), mi_elite->fitness);
    /***************************************************************************************
     *                                          Liberar memoria
     ***************************************************************************************/
    liberarNodos( estructura_MIMIC );

    for(unsigned int i = 0; i < n_vars; i++){
        free( frecuencias[i] );
        free( info_mat[i] );
    }
    free(frecuencias);
    free( info_mat );

    freeIndivs( poblacion, n_pob );

    return mi_elite;
}







/*  Metodo:  optiMIMICCont
 *   Funcion: Encuentra el optimo global de la funcion 'mi_fun' mediante el algoritmo MIMIC de optimizacion para variables continuas.
 */
INDIV *optiMIMICCont(FUNPARS *mis_optipars, FUNPTR mi_fun, FUNPARS *mis_pars, STAUS *mi_semilla)
{

    const unsigned int n_iters = *(unsigned int*)((char*)mis_optipars->parameters + *(mis_optipars->ptr_starts));
    const unsigned int n_pob =   *(unsigned int*)((char*)mis_optipars->parameters + *(mis_optipars->ptr_starts + 1));
    const unsigned int n_seleccion = *(unsigned int*)((char*)mis_optipars->parameters + *(mis_optipars->ptr_starts + 2));
    const unsigned int n_vars = *(unsigned int*)((char*)mis_optipars->parameters + *(mis_optipars->ptr_starts + 3));

    //// Requerir memoria para las tablas de frecuencias e informacion mutua:
    double **covs_mat = (double**) malloc( n_vars * sizeof(double*) );
    double *medias = (double*) malloc( n_vars*sizeof(double) );

    double **info_mat = (double**) malloc( n_vars * sizeof(double*) );

    for( unsigned int i = 0; i < n_vars; i++){
        *(covs_mat + i) = (double*) malloc( (n_vars-i)* sizeof(double) );
        *(info_mat + i) = (double*) malloc( (n_vars-i)* sizeof(double) );
    }

    //// Requerir memoria e inicializar la poblacion 0:
    INDIV *poblacion = (INDIV *) malloc( n_pob * sizeof( INDIV ) );

    for( unsigned int i = 0; i < n_pob; i++){

        (poblacion + i)->dim = n_vars;
        (poblacion + i)->ptr_starts = (size_t*) malloc((n_vars+1) * sizeof(size_t));
        *((poblacion + i)->ptr_starts) = 0;
        for(unsigned int k = 1; k <= n_vars; k++){
            *( (poblacion + i)->ptr_starts + k ) = *( (poblacion + i)->ptr_starts + k - 1) + sizeof(double);
        }
        (poblacion + i)->my_features = malloc( *((poblacion + i)->ptr_starts + n_vars) );

        // Generar la cadena de cada individuo de manera aleatoria.
        for(unsigned int k = 0; k < n_vars; k++){
            const double min = *(double*)((char*)mis_pars->parameters + *(mis_pars->ptr_starts + 2*k  ));
            const double max = *(double*)((char*)mis_pars->parameters + *(mis_pars->ptr_starts + 2*k+1));

            *(double*)((char*)(poblacion + i)->my_features + *((poblacion + i)->ptr_starts + k) ) = HybTaus( min, max, mi_semilla);
        }

        // Evaluar individuo en la funcion objetivo:
        (poblacion + i)->fitness = mi_fun( mis_pars, poblacion + i );
    }

    // Ordernar la poblacion segun sus fitness
    qsort((void*)poblacion, n_pob, sizeof(INDIV), compIndiv); // El mejor fitness queda en la posicion 0

    // Proponer la solucion al problema de optimizacion como el individuo con mejor fitness
    INDIV *mi_elite = (INDIV*) malloc( sizeof(INDIV) );
    mi_elite->dim = n_vars;
    mi_elite->ptr_starts = (size_t*) malloc( (n_vars+1)*sizeof(size_t) );
    mi_elite->my_features = malloc( *(poblacion->ptr_starts + n_vars) );

    mi_elite->fitness = poblacion->fitness;
    memcpy(mi_elite->ptr_starts, poblacion->ptr_starts, (n_vars+1)*sizeof(size_t));
    memcpy(mi_elite->my_features, poblacion->my_features, *(poblacion->ptr_starts + n_vars));

    /***************************************************************************************
     * Estimar las estructuras y parametros a partir de la poblacion '0':
     ***************************************************************************************/
    // Obtener las frecuencias de aparicion de 1's para cada variable y su condicional:
    calcCovars(poblacion, n_seleccion, n_vars, covs_mat, medias);

    // Calcular las Entropias marginales para todas las varibales y determinar el minimo como el nodo raiz de la cadena.
    calcInfoMutuaCont(info_mat, n_vars, covs_mat);

    // Generar el modelo a partir de la informacion obtenida dela poblacion en el momento '0':
    NODO_ESTRUCT *estructura_MIMIC = estructMIMICCont(n_vars, covs_mat, medias, info_mat);

    printf("generacion\tevaluaciones\tevaluacion_elite\n");
    printf("%i\t%i\t%f\n", 0, n_pob, mi_elite->fitness);
    for(unsigned int t = 1; t < n_iters; t++){

        /***************************************************************************************
         * Generar una poblacion 't':
         ***************************************************************************************/
        genPobCont( poblacion, n_pob, estructura_MIMIC, mi_fun, mis_pars, mi_semilla );

        // Ordernar la poblacion segun sus fitness
        qsort((void*)poblacion, n_pob, sizeof(INDIV), compIndiv); // El mejor fitness queda en la posicion 0

        // Si un nuevo candidato mejora la solucion, se actualiza la solucion propuesta.
        if( mi_elite->fitness < poblacion->fitness ){
            mi_elite->fitness = poblacion->fitness;
            memcpy(mi_elite->my_features, poblacion->my_features, *(poblacion->ptr_starts + n_vars));
        }
        printf("%i\t%i\t%f\n", t, n_pob * (t+1), mi_elite->fitness);

        /***************************************************************************************
         * Estimar las estructuras y parametros a partir de la poblacion 't':
         ***************************************************************************************/
        // Obtener las frecuencias de aparicion de 1's para cada variable y su condicional:
        calcCovars(poblacion, n_seleccion, n_vars, covs_mat, medias);

        // Calcular las Entropias marginales para todas las varibales y determinar el minimo como el nodo raiz de la cadena.
        calcInfoMutuaCont(info_mat, n_vars, covs_mat);

        // Generar el modelo a partir de la informacion obtenida dela poblacion en el momento 't':
        liberarNodos( estructura_MIMIC );
        estructura_MIMIC = estructMIMICCont(n_vars, covs_mat, medias, info_mat);
    }

    printf("%i\t%i\t%f\n", n_iters, n_pob * (n_iters+1), mi_elite->fitness);
    /***************************************************************************************
     *                                          Liberar memoria
     ***************************************************************************************/
    liberarNodos( estructura_MIMIC );

    for(unsigned int i = 0; i < n_vars; i++){
        free( covs_mat[i] );
        free( info_mat[i] );
    }
    free( covs_mat );
    free( info_mat );
    free( medias );

    freeIndivs( poblacion, n_pob );

    return mi_elite;
}




/*  Metodo:  optiCHOWLIU
 *   Funcion: Encuentra el optimo global de la funcion 'mi_fun' mediante el algoritmo CHOWLIU de optimizacion.
 */
INDIV *optiCHOWLIU( FUNPARS *mis_optipars, FUNPTR mi_fun, FUNPARS *mis_pars, STAUS *mi_semilla)
{

    const unsigned int n_iters = *(unsigned int*)((char*)mis_optipars->parameters + *(mis_optipars->ptr_starts));
    const unsigned int n_pob = *(unsigned int*)((char*)mis_optipars->parameters + *(mis_optipars->ptr_starts + 1));
    const unsigned int n_seleccion = *(unsigned int*)((char*)mis_optipars->parameters + *(mis_optipars->ptr_starts + 2));
    const unsigned int n_vars = *(unsigned int*)((char*)mis_optipars->parameters + *(mis_optipars->ptr_starts + 3));


    //// Requerir memoria para las tablas de frecuencias e informacion mutua:
    unsigned int **frecuencias = (unsigned int**) malloc( n_vars * sizeof(unsigned int*) );
    double **info_mat = (double**) malloc( n_vars * sizeof(double*) );
    for(unsigned int i = 0; i < n_vars; i++){
        frecuencias[i] = (unsigned int*) malloc( (n_vars-i)* sizeof(unsigned int) );
        info_mat[i] =  (double*) malloc( (n_vars-i)* sizeof(double) );
    }

    //// Requerir memoria e inicializar la poblacion 0:
    INDIV *poblacion = (INDIV *) malloc( n_pob * sizeof( INDIV ) );
    for(unsigned int i = 0; i < n_pob; i++){
        (poblacion + i)->dim = n_vars;
        (poblacion + i)->ptr_starts = (size_t*) malloc((n_vars+1) * sizeof(size_t));
        *( (poblacion + i)->ptr_starts ) = 0;
        for(unsigned int k = 1; k <= n_vars; k++){
            *( (poblacion + i)->ptr_starts + k ) = *( (poblacion + i)->ptr_starts + k - 1) + sizeof(char);
        }
        (poblacion + i)->my_features = malloc( *((poblacion + i)->ptr_starts + n_vars) );

        // Generar la cadena de cada individuo de manera aleatoria.
        for(unsigned int k = 0; k < n_vars; k++){
            *((char*)(poblacion + i)->my_features + *((poblacion + i)->ptr_starts + k) ) = HybTaus( 0.0, 1.0, mi_semilla) <= 0.5 ? 0 : 1;
        }

        // Evaluar individuo en la funcion objetivo:
        (poblacion + i)->fitness = mi_fun( mis_pars, poblacion + i );
    }

    // Ordernar la poblacion segun sus fitness
    qsort((void*)poblacion, n_pob, sizeof(INDIV), compIndiv); // El mejor fitness queda en la posicion 0

    // Proponer la solucion al problema de optimizacion como el individuo con mejor fitness
    INDIV *mi_elite = (INDIV*) malloc( sizeof(INDIV) );
    mi_elite->dim = n_vars;
    mi_elite->ptr_starts = (size_t*) malloc( (n_vars+1)*sizeof(size_t) );
    mi_elite->my_features = malloc( *(poblacion->ptr_starts + n_vars) );

    mi_elite->fitness = poblacion->fitness;
    memcpy(mi_elite->ptr_starts, poblacion->ptr_starts, (n_vars+1)*sizeof(size_t));
    memcpy(mi_elite->my_features, poblacion->my_features, *(poblacion->ptr_starts + n_vars));

    /***************************************************************************************
     * Estimar las estructuras y parametros a partir de la poblacion '0':
     ***************************************************************************************/
    // Obtener las frecuencias de aparicion de 1's para cada variable y su condicional:
    calcFrecs( poblacion, n_seleccion, n_vars, frecuencias );

    // Calcular las Entropias marginales para todas las varibales y determinar el minimo como el nodo raiz de la cadena.
    calcInfoMutua(info_mat, n_seleccion, n_vars, frecuencias );

    // Generar el modelo a partir de la informacion obtenida dela poblacion en el momento '0':
    NODO_ESTRUCT *estructura_CHOWLIU = estructCHOWLIU(n_seleccion, n_vars, frecuencias, info_mat);

    printf("generacion\tevaluaciones\tevaluacion_elite\n");
    printf("%i\t%i\t%f\n", 0, n_pob, mi_elite->fitness);
    for(unsigned int t = 1; t < n_iters; t++){

        /***************************************************************************************
         * Generar una poblacion 't':
         ***************************************************************************************/
        genPob( poblacion, n_pob, estructura_CHOWLIU, mi_fun, mis_pars, mi_semilla );

        // Ordernar la poblacion segun sus fitness
        qsort((void*)poblacion, n_pob, sizeof(INDIV), compIndiv); // El mejor fitness queda en la posicion 0

        // Si un nuevo candidato mejora la solucion, se actualiza la solucion propuesta.
        if( mi_elite->fitness < poblacion->fitness ){
            mi_elite->fitness = poblacion->fitness;
            memcpy(mi_elite->my_features, poblacion->my_features, *(poblacion->ptr_starts + n_vars));
        }
        printf("%i\t%i\t%f\n", t, n_pob * (t+1), mi_elite->fitness);

        /***************************************************************************************
         * Estimar las estructuras y parametros a partir de la poblacion 't':
         ***************************************************************************************/
        // Obtener las frecuencias de aparicion de 1's para cada variable y su condicional:
        calcFrecs( poblacion, n_seleccion, n_vars, frecuencias );

        // Calcular las Entropias marginales para todas las varibales y determinar el minimo como el nodo raiz de la cadena.
        calcInfoMutua(info_mat, n_seleccion, n_vars, frecuencias );

        // Generar el modelo a partir de la informacion obtenida dela poblacion en el momento 't':
        liberarNodos( estructura_CHOWLIU );
        estructura_CHOWLIU = estructCHOWLIU(n_seleccion, n_vars, frecuencias, info_mat);
    }
    printf("%i\t%i\t%f\n", n_iters, n_pob * (n_iters+1), mi_elite->fitness);

    /***************************************************************************************
     *                                          Liberar memoria
     ***************************************************************************************/
    liberarNodos( estructura_CHOWLIU );

    for(unsigned int i = 0; i < n_vars; i++){
        free( frecuencias[i] );
        free( info_mat[i] );
    }
    free(frecuencias);
    free( info_mat );

    freeIndivs( poblacion, n_pob );

    return mi_elite;
}







/*  Metodo:  optiCHOWLIUCont
 *   Funcion: Encuentra el optimo global de la funcion 'mi_fun' mediante el algoritmo CHOWLIU de optimizacion para variables continuas.
 */
INDIV *optiCHOWLIUCont( FUNPARS *mis_optipars, FUNPTR mi_fun, FUNPARS *mis_pars, STAUS *mi_semilla)
{

    const unsigned int n_iters = *(unsigned int*)((char*)mis_optipars->parameters + *(mis_optipars->ptr_starts));
    const unsigned int n_pob = *(unsigned int*)((char*)mis_optipars->parameters + *(mis_optipars->ptr_starts + 1));
    const unsigned int n_seleccion = *(unsigned int*)((char*)mis_optipars->parameters + *(mis_optipars->ptr_starts + 2));
    const unsigned int n_vars = *(unsigned int*)((char*)mis_optipars->parameters + *(mis_optipars->ptr_starts + 3));


    //// Requerir memoria para las tablas de frecuencias e informacion mutua:
    double **covs_mat = (double**) malloc( n_vars * sizeof(double*) );
    double *medias= (double*) malloc( n_vars * sizeof(double) );
    double **info_mat = (double**) malloc( n_vars * sizeof(double*) );
    for(unsigned int i = 0; i < n_vars; i++){
        covs_mat[i] = (double*) malloc( (n_vars-i)* sizeof(double) );
        info_mat[i] = (double*) malloc( (n_vars-i)* sizeof(double) );
    }

    //// Requerir memoria e inicializar la poblacion 0:
    INDIV *poblacion = (INDIV *) malloc( n_pob * sizeof( INDIV ) );
    const double min = *(double*)((char*)mis_pars->parameters + *(mis_pars->ptr_starts + 1));
    const double max = *(double*)((char*)mis_pars->parameters + *(mis_pars->ptr_starts + 2));

    for(unsigned int i = 0; i < n_pob; i++){
        (poblacion + i)->dim = n_vars;
        (poblacion + i)->ptr_starts = (size_t*) malloc((n_vars+1) * sizeof(size_t));
        *( (poblacion + i)->ptr_starts ) = 0;
        for(unsigned int k = 1; k <= n_vars; k++){
            *( (poblacion + i)->ptr_starts + k ) = *( (poblacion + i)->ptr_starts + k - 1) + sizeof(double);
        }
        (poblacion + i)->my_features = malloc( *((poblacion + i)->ptr_starts + n_vars) );

        // Generar la cadena de cada individuo de manera aleatoria.
        for(unsigned int k = 0; k < n_vars; k++){
            *(double*)((char*)(poblacion + i)->my_features + *((poblacion + i)->ptr_starts + k) ) = HybTaus( min, max, mi_semilla);
        }

        // Evaluar individuo en la funcion objetivo:
        (poblacion + i)->fitness = mi_fun( mis_pars, poblacion + i );
    }

    // Ordernar la poblacion segun sus fitness
    qsort((void*)poblacion, n_pob, sizeof(INDIV), compIndiv); // El mejor fitness queda en la posicion 0

    // Proponer la solucion al problema de optimizacion como el individuo con mejor fitness
    INDIV *mi_elite = (INDIV*) malloc( sizeof(INDIV) );
    mi_elite->dim = n_vars;
    mi_elite->ptr_starts = (size_t*) malloc( (n_vars+1)*sizeof(size_t) );
    mi_elite->my_features = malloc( *(poblacion->ptr_starts + n_vars) );

    mi_elite->fitness = poblacion->fitness;
    memcpy(mi_elite->ptr_starts, poblacion->ptr_starts, (n_vars+1)*sizeof(size_t));
    memcpy(mi_elite->my_features, poblacion->my_features, *(poblacion->ptr_starts + n_vars));

    /***************************************************************************************
     * Estimar las estructuras y parametros a partir de la poblacion '0':
     ***************************************************************************************/
    // Obtener las frecuencias de aparicion de 1's para cada variable y su condicional:
    calcCovars( poblacion, n_seleccion, n_vars, covs_mat, medias);

    // Calcular las Entropias marginales para todas las varibales y determinar el minimo como el nodo raiz de la cadena.
    calcInfoMutuaCont(info_mat, n_vars, covs_mat );

    // Generar el modelo a partir de la informacion obtenida dela poblacion en el momento '0':
    NODO_ESTRUCT *estructura_CHOWLIU = estructCHOWLIUCont(n_vars, covs_mat, medias, info_mat);

    printf("generacion\tevaluaciones\tevaluacion_elite\n");
    printf("%i\t%i\t%f\n", 0, n_pob, mi_elite->fitness);
    for(unsigned int t = 1; t < n_iters; t++){

        /***************************************************************************************
         * Generar una poblacion 't':
         ***************************************************************************************/
        genPobCont( poblacion, n_pob, estructura_CHOWLIU, mi_fun, mis_pars, mi_semilla );

        // Ordernar la poblacion segun sus fitness
        qsort((void*)poblacion, n_pob, sizeof(INDIV), compIndiv); // El mejor fitness queda en la posicion 0

        // Si un nuevo candidato mejora la solucion, se actualiza la solucion propuesta.
        if( mi_elite->fitness < poblacion->fitness ){
            mi_elite->fitness = poblacion->fitness;
            memcpy(mi_elite->my_features, poblacion->my_features, *(poblacion->ptr_starts + n_vars));
        }
        printf("%i\t%i\t%f\n", t, n_pob * (t+1), mi_elite->fitness);

        /***************************************************************************************
         * Estimar las estructuras y parametros a partir de la poblacion 't':
         ***************************************************************************************/
        // Obtener las frecuencias de aparicion de 1's para cada variable y su condicional:
        calcCovars( poblacion, n_seleccion, n_vars, covs_mat, medias);

        // Calcular las Entropias marginales para todas las varibales y determinar el minimo como el nodo raiz de la cadena.
        calcInfoMutuaCont(info_mat, n_vars, covs_mat);

        // Generar el modelo a partir de la informacion obtenida dela poblacion en el momento 't':
        liberarNodos( estructura_CHOWLIU );
        estructura_CHOWLIU = estructCHOWLIUCont(n_vars, covs_mat, medias, info_mat);
    }
    printf("%i\t%i\t%f\n", n_iters, n_pob * (n_iters+1), mi_elite->fitness);

    /***************************************************************************************
     *                                          Liberar memoria
     ***************************************************************************************/
    liberarNodos( estructura_CHOWLIU );

    for(unsigned int i = 0; i < n_vars; i++){
        free( covs_mat[i] );
        free( info_mat[i] );
    }
    free( covs_mat );
    free( info_mat );
    free( medias );

    freeIndivs( poblacion, n_pob );

    return mi_elite;
}







/*
    Metodo:     selecPob
    Funcion:    Genera el grupo de individuos seleccionados.
*/
void selecPob(INDIV* sel, const INDIV* pob, const int n_pob, const int n_vars, STAUS *mi_semilla){

    //// Suma de los fintess solo de los mejores 'n_seleccion' individuos de la poblacion:
    double suma_fitness = 0.0;
    for( int i = 0; i < n_pob; i++){
        suma_fitness += (pob + i)->fitness;
    }

    // Seleccionar el numero de individuos correspondiente phacia la nueva poblacion:
    int sel_pos = 0;
    for( int i = 0; i < n_pob; i++){
        // Se copia a la seleccion el individuo i el proporcionalmente a lo que aporta al fitness global:
        unsigned int n_sel = (unsigned int) round( (pob+i)->fitness / suma_fitness * n_pob );
		unsigned int j = 1;
        while( ( sel_pos < n_pob ) && (j <= n_sel) ){
            memcpy( (sel + sel_pos)->my_features, (pob + i)->my_features, *(pob->ptr_starts + n_vars) );
            (sel + sel_pos)->fitness = (pob + sel_pos)->fitness;
            sel_pos++;
            j++;
        }

        //si se completo la seleccion, se termina de tomar individuos:
        if( sel_pos == n_pob ){
            break;
        }
    }
}







/*  Metodo:  crossOver
 *   Funcion: Genera un conjunto de nuevos indiviudos a partir de la cruza de los individuos en la seleccion.
 */
void crossOver(INDIV* pob, INDIV* sel, const int n_pob, STAUS *mi_semilla, int *indices, const double p_cruza, const double p_mutacion, FUNPTR mi_fun, FUNPARS *mis_pars){
    const int n_bits = pob->dim;

    //// Crear una permutacion aleatoria de los indices para realizar la cruza:
    shuffle(indices, 0, n_pob, mi_semilla);

    // Generar dos nuevos individuos a partir de dos padres, el proceso se repite hasta completar la fraccion 'sel'.
    for(int i = 0; i < n_pob; i+=2){
        // Si el individuo no se cruza:
        if( HybTaus(0.0, 1.0, mi_semilla) > p_cruza ){
            memcpy( (pob + i  )->my_features, (sel + indices[i  ])->my_features, *(pob->ptr_starts + n_bits) );
            memcpy( (pob + i+1)->my_features, (sel + indices[i+1])->my_features, *(pob->ptr_starts + n_bits) );

            (pob + i  )->fitness = (sel + indices[i  ])->fitness;
            (pob + i+1)->fitness = (sel + indices[i+1])->fitness;

        }else{ // Si la cruza se lleva a cabo:
            int padre_1 = indices[i  ];
            int padre_2 = indices[i+1];
            int bits_ini = 0;

            //// Cruza uniforme
            // Realizar cortes hasta terminar con la secuencia de los individuos:
            while(bits_ini < n_bits){
                const unsigned int n_bits_cpy = (unsigned int) round( HybTaus(0.5, (double)(n_bits - bits_ini) + 0.49, mi_semilla) );
                const size_t v_bits = *(pob->ptr_starts + n_bits_cpy + bits_ini) - *(pob->ptr_starts + bits_ini);
                memcpy((char*)(pob + i  )->my_features + *(pob->ptr_starts + bits_ini), (char*)(sel + padre_1)->my_features + *(pob->ptr_starts + bits_ini), v_bits);
                memcpy((char*)(pob + i+1)->my_features + *(pob->ptr_starts + bits_ini), (char*)(sel + padre_2)->my_features + *(pob->ptr_starts + bits_ini), v_bits);

                int padre_swap = padre_1;
                padre_1 = padre_2;
                padre_2 = padre_swap;

                bits_ini += n_bits_cpy;
            }
        }

        // Realizar mutacion sobre los hijos generados con cierta probabilidad: 'prob_mutacion'
        // Mutar el hijo 1:
        if( HybTaus(0.0, 1.0, mi_semilla) <= p_mutacion ){
            // Seleccionar un gen a mutar:
            const unsigned int gen = (unsigned int) round( HybTaus(-0.5, (double)n_bits - 0.51, mi_semilla) );
            *((char*)(pob + i  )->my_features + *(pob->ptr_starts + gen)) = 1 - *((char*)(pob + i  )->my_features + *(pob->ptr_starts + gen));
        }
        //// Evaluar los individuos generados en la funcion objetivo:
        (pob + i  )->fitness = mi_fun( mis_pars, pob + i  );


        // Mutar el hijo 2:
        if( HybTaus(0.0, 1.0, mi_semilla) <= p_mutacion ){
            // Seleccionar un gen a mutar:
            const unsigned int gen = (unsigned int) round( HybTaus(-0.5, (double)n_bits - 0.51, mi_semilla) );
            *((char*)(pob + i+1)->my_features + *(pob->ptr_starts + gen)) = 1 - *((char*)(pob + i+1)->my_features + *(pob->ptr_starts + gen));
        }
        //// Evaluar los individuos generados en la funcion objetivo:
        (pob + i+1)->fitness = mi_fun( mis_pars, pob + i+1);
    }
}






/*  Metodo:  optiGA
 *   Funcion: Encuentra el optimo global de la funcion 'mi_fun' mediante el algoritmo genetico de optimizacion.
 */
INDIV *optiGA(FUNPARS *mis_optipars, FUNPTR mi_fun, FUNPARS *mis_pars, STAUS *mi_semilla)
{
    const unsigned int n_iters = *(unsigned int *)((char*)mis_optipars->parameters + *(mis_optipars->ptr_starts));
    const unsigned int n_pob = *(unsigned int *)((char*)mis_optipars->parameters + *(mis_optipars->ptr_starts + 1));
    const unsigned int n_vars = *(unsigned int *)((char*)mis_optipars->parameters + *(mis_optipars->ptr_starts + 3));
    const double p_cruza = *(double*)((char*)mis_optipars->parameters + *(mis_optipars->ptr_starts + 4));
    const double p_mutacion = *(double*)((char*)mis_optipars->parameters + *(mis_optipars->ptr_starts + 5));

    //// Requerir memoria e inicializar la poblacion 0:
    INDIV *poblacion = (INDIV *) malloc( n_pob * sizeof( INDIV ) );
    INDIV *seleccion = (INDIV *) malloc( n_pob * sizeof( INDIV ) );

    for(unsigned int i = 0; i < n_pob; i++){
        (poblacion + i)->dim = n_vars;
        (seleccion + i)->dim = n_vars;

        (poblacion + i)->ptr_starts = (size_t*) malloc((n_vars+1) * sizeof(size_t));
        (seleccion + i)->ptr_starts = (size_t*) malloc((n_vars+1) * sizeof(size_t));

        *( (poblacion + i)->ptr_starts ) = 0;
        *( (seleccion + i)->ptr_starts ) = 0;
        for(unsigned int k = 1; k <= n_vars; k++){
            *( (poblacion + i)->ptr_starts + k ) = *( (poblacion + i)->ptr_starts + k - 1) + sizeof(char);
            *( (seleccion + i)->ptr_starts + k ) = *( (seleccion + i)->ptr_starts + k - 1) + sizeof(char);
        }
        (poblacion + i)->my_features = malloc( *((poblacion + i)->ptr_starts + n_vars) );
        (seleccion + i)->my_features = malloc( *((seleccion + i)->ptr_starts + n_vars) );

        // Generar la cadena de cada individuo de manera aleatoria.
        for(unsigned int k = 0; k < n_vars; k++){
            *((char*)(poblacion + i)->my_features + *((poblacion + i)->ptr_starts + k) ) = HybTaus( 0.0, 1.0, mi_semilla) <= 0.5 ? 0 : 1;
        }

        // Evaluar individuo en la funcion objetivo:
        (poblacion + i)->fitness = mi_fun( mis_pars, poblacion + i );
    }

    // Ordernar la poblacion segun sus fitness
    qsort((void*)poblacion, n_pob, sizeof(INDIV), compIndiv); // El mejor fitness queda en la posicion 0

    // Proponer la solucion al problema de optimizacion como el individuo con mejor fitness
    INDIV *mi_elite = (INDIV*) malloc( sizeof(INDIV) );
    mi_elite->dim = n_vars;
    mi_elite->ptr_starts = (size_t*) malloc( (n_vars+1)*sizeof(size_t) );
    mi_elite->my_features = malloc( *(poblacion->ptr_starts + n_vars) );

    mi_elite->fitness = poblacion->fitness;
    memcpy(mi_elite->ptr_starts, poblacion->ptr_starts, (n_vars+1)*sizeof(size_t));
    memcpy(mi_elite->my_features, poblacion->my_features, *(poblacion->ptr_starts + n_vars));

    int* indices = (int*) malloc(n_pob*sizeof(int));

    printf("generacion\tevaluaciones\tevaluacion_elite\n");
    printf("%i\t%i\t%f\n", 0, n_pob, mi_elite->fitness);
    for(unsigned int t = 1; t < n_iters; t++){

        /***************************************************************************************
         Elegir el conjunto de individuos seleccionados a partir de la poblacion 't':
         ***************************************************************************************/
        selecPob( seleccion, poblacion, n_pob, n_vars, mi_semilla );

        /***************************************************************************************
         * Realizar la cruza de los individuos a partir de la seleccion 't':
         ***************************************************************************************/
        crossOver( poblacion, seleccion, n_pob, mi_semilla, indices, p_cruza, p_mutacion, mi_fun, mis_pars );

        // Ordernar la poblacion segun sus fitness
        qsort((void*)poblacion, n_pob, sizeof(INDIV), compIndiv); // El mejor fitness queda en la posicion 0

        // Si un nuevo candidato mejora la solucion, se actualiza la solucion propuesta.
        if( mi_elite->fitness < poblacion->fitness ){
            mi_elite->fitness = poblacion->fitness;
            memcpy(mi_elite->my_features, poblacion->my_features, *(poblacion->ptr_starts + n_vars));
        }
        printf("%i\t%i\t%f\n", t, n_pob * (t+1), mi_elite->fitness);

    }
    printf("%i\t%i\t%f\n", n_iters, n_pob * (n_iters+1), mi_elite->fitness);

    /***************************************************************************************
     *                                          Liberar memoria
     ***************************************************************************************/
    free(indices);
    freeIndivs( poblacion, n_pob );
    freeIndivs( seleccion, n_pob );

    return mi_elite;
}
