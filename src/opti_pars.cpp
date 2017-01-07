/************************************************************************************************************
*                                                                                                           *
* CENTRO DE INVESTIGACION EN MATEMATICAS                                                                    *
* DOCTORADO EN CIENCIAS DE LA COMPUTACION                                                                   *
* FERNANDO CERVANTES SANCHEZ                                                                                *
*                                                                                                           *
* FILE NAME: IMGCONT.cpp                                                                                    *
*                                                                                                           *
* PURPOSE: Implementation of the OPTI_PARS class for optimization of the parameters of the GMF and SSG.     *
*                                                                                                           *
* DEVELOPMENT HISTORY:                                                                                      *
* Date           Author        Change Id    Release    Description Of Change                                *
* 06/Jan/2017    Fernando C.   0            1.0        Creation                                             *
*                                                                                                           *
************************************************************************************************************/



#include "opti_pars.h"


/*  Metodo: setEvoMet
    Funcion: Define el metodo evolutivo que se usara para establecer los mejores aprametros del filtro.
*/
void OPTI_PARS::setEvoMet(const EVO_MET evo_met){
    metodo_elegido = evo_met;
}



/*  Metodo: setEvoMetPar
    Funcion: Define elparametro para el metodo de optimizacion.
*/
void OPTI_PARS::setEvoMetPar(const EVO_MET_PAR evo_par, const double val){
    switch( evo_par ){
    case POPSIZE:
        n_pob = (int)val;
        break;
    case MAXGEN:
        max_iters = (int)val;
        break;
    case CR:
        prob_cruza = val;
        seleccion = (int)(val * (double)n_pob);
        seleccion += (seleccion%2);
        DEB_MSG("Selection size: " << seleccion);
        break;
    case MR:
        prob_mutacion = val;
        break;
    }
}


/*  Metodo: setFitness
    Funcion: Define que funcion de evaluacion se usara en el metodo de optimizacion como funcion objetivo.
*/
void OPTI_PARS::setFitness( const FITNESS fit_fun){
    fitness_elegido = fit_fun;
}



/*  Metodo: OPTI_PARS

    Funcion: Inicializa los parametros en valores por defecto.
*/
OPTI_PARS::OPTI_PARS(){
#ifdef BUILD_GUI_VERSION
    mi_pBar = NULL;
    mi_txtLog = NULL;
#endif

    n_pars = 4;
    for( int p = 0; p < 4; p++){
        // Por defecto, se optimizan todos los parametros.
        pars_optim[p] = true;
    }

    mi_elite = new INDIV;
    mi_elite->eval = 0.0;
    memset( mi_elite->vars, 0, 4*sizeof(double) );
	memset( mi_elite->cadena, 0, 128 * sizeof(unsigned char));

    semilla = NULL;

    metodo_elegido = EVO_UNSET;
    fitness_elegido = FIT_UNSET;
	
    n_pob = 0;
    prob_mutacion = 0.0;
    seleccion = 0;
    max_iters = 0;
}




/*  Metodo: ~OPTI_PARS

    Funcion: Libera la memoria requerida por el objeto.
*/
OPTI_PARS::~OPTI_PARS(){
    delete mi_elite;
    if ( semilla ){
        delete semilla;
    }
}



/*  Metodo: setPar
    Funcion: Ejecuta el metodo de optimizacion seleccionado para encontrar los mejores parametros para el filtro sobre la imagen de origen.
*/
void OPTI_PARS::setPar(){
    switch( metodo_elegido ){
        case EA_DE:
            DE();
            break;

        case EA_GA:
            GA();
            break;

        case EDA_BUMDA:
            BUMDA();
            break;

        case EDA_UMDA:
            UMDA();
            break;

        case EXHAUSTIVA:
            busquedaExhaustiva();
            break;
    }
}


/*  Metodo: setPar
    Funcion: Establece los parametros del filtro.
*/
void OPTI_PARS::setPar( const PARAMETRO par, const double val ){
    mi_elite->vars[ par ] = val;
    // Si se fija el parametro, se omite este en la busqueda exhaustiva/metaheuristica
    pars_optim[ par ] = false;
    min_vars[ par ] = 0.0;
    lim_inf[ par ] = val;
}





/*  Metodo: getPars
    Funcion: Retorna los parametros utilzados para el filtro.
*/
OPTI_PARS::INDIV OPTI_PARS::getPars(){
    INDIV out_pars;
    memcpy(out_pars.vars, mi_elite->vars, 5*sizeof(double));
    out_pars.eval = mi_elite->eval;
    return out_pars;
}





/*  Metodo: getPars
    Funcion: Retorna los parametros utilzados para el filtro.
*/
int OPTI_PARS::getParametrosOptimizar(){
    int n_pars = 0;
    for( int i = 0; i < 4; i++){
        DEB_MSG("Par[" << i << "]: " << pars_optim[i] );
        if( pars_optim[i] ){
            n_pars++;
        }
    }

    DEB_MSG("Optimizar " << n_pars << " parametros.");

    return n_pars;
}




/*  Metodo: setLim
    Funcion: Establece los limites de busqueda de los parametros, si se utiliza para metodos de codificacion binaria, var_delta indica cuantos bits se utilizaran.
*/
void OPTI_PARS::setLim( const PARAMETRO par, const LIMITES lim, const double val){
    switch( lim ){
    case INFERIOR:
        lim_inf[ par ] = val;
        break;
    case SUPERIOR:
        lim_sup[ par ] = val;
        break;
    case DELTA:
        min_vars[ par ] = val;
        break;
    }
    pars_optim[ par ] = true;
}



/*  Metodo: fitnessROC
    Funcion: Evalua la eficacia de un detector utilizando el area bajo la curva de ROC:
*/
double OPTI_PARS::fitnessROC( INDIV *test )
{
	this->filter();
	return this->calcROC();
}




/*  Metodo: fitnessCorCon
    Funcion: Evalua los parametros 'L', 'T' y //'K'// para el filtro de Gabor y el la correlacion y contraste:
*/
double OPTI_PARS::fitnessCorCon( INDIV *test )
{
	this->filter();
    return this->calcCorCon();
}




/*
    Metodo:        compIndiv
    Funcion:    Compara el indidividuo A contra el B.
        Si la evaluacion del individuo A es menor a la del B, se retorna +,
                             si es igual, se retorna 0,
                             si es mayor, se retorna -.
*/
int OPTI_PARS::compIndiv(const void* A, const void* B){
    INDIV *indA = (INDIV*)A;
    INDIV *indB = (INDIV*)B;
    if(indA->eval < indB->eval) return  1;
    if(indA->eval == indB->eval) return  0;
    return -1;
}




// BUMDA ----- BUMDA ----- BUMDA ----- BUMDA ----- BUMDA ----- BUMDA ----- BUMDA ----- BUMDA ----- BUMDA ----- BUMDA ----- BUMDA ----- BUMDA ----- BUMDA ----- BUMDA ----- BUMDA ----- BUMDA ----- BUMDA ----- BUMDA ----- BUMDA ----- BUMDA ----- BUMDA ----- BUMDA ----- BUMDA ----- BUMDA ----- BUMDA ----- BUMDA ----- BUMDA ----- BUMDA ----- BUMDA ----- BUMDA ----- BUMDA
/*
    Metodo:        generarPobInicial
    Funcion:    Genera una poblacion inicial dentro de los limites establecidos.

*/
void OPTI_PARS::generarPobInicial(INDIV *poblacion){
    // VERIFICAR CUALES PARAMETROS SE BUSCAN:
    using namespace std;
    for(int i = 0; i < n_pob; i++){
        memcpy( poblacion[i].vars, mi_elite->vars, 5*sizeof(double) );

        for( unsigned int j = 0; j < n_pars; j++){
            const unsigned int k = idx_pars[j];
            (poblacion + i)->vars[k] = HybTaus(lim_inf[k], lim_sup[k], semilla);
        }

        switch( fitness_elegido ){
            case ROC:
                (poblacion + i)->eval = fitnessROC( poblacion + i);
                break;
            case CORCON:
                (poblacion + i)->eval = fitnessCorCon( poblacion + i);
                break;
        }
    }
}





/*
    Metodo:        generarPob
    Funcion:    Genera una poblacion de tamaño 'n_pob' individuos con una distribucion normal con medias y varianzas apra cada parametro los indicados.

*/
void OPTI_PARS::generarPob(double * medias, double * varianzas, INDIV *poblacion){
    double val_gen;

    //Se generan n_pob - 1 individuos, y permanece el elite que se tiene.
    // VERIFICAR CUALES PARAMETROS SE BUSCAN:
    using namespace std;
    for(int i = 0; i < n_pob; i++){
        for(unsigned int j = 0; j < n_pars; j++){
            const unsigned int k = idx_pars[j];
            val_gen =  anorm( medias[k], varianzas[k], semilla);
            val_gen = (val_gen <= lim_sup[ k ]) ? ((val_gen >= lim_inf[ k ]) ? val_gen : lim_inf[ k ]) : lim_sup[ k ];
            (poblacion + i)->vars[ k ] = val_gen;
        }
        switch( fitness_elegido ){
            case ROC:
                (poblacion + i)->eval = fitnessROC( poblacion + i);
                break;

            case CORCON:
                (poblacion + i)->eval = fitnessCorCon( poblacion + i);
                break;
        }
    }
}



/*
    Metodo:        calcularPars
    Funcion:    Calula la media y la varianza para cada variable de la poblacion, utilizando solo a los individuos de la seleccion.
*/
void OPTI_PARS::calcularPars(const INDIV *poblacion, const int truncamiento, double *medias, double *varianzas){

    double sum_evals = 0.0;
    const double gx_sel = (poblacion + truncamiento)->eval - 1e-6;
        //Se calculan las medias primero:
        memset(medias, 0, 4*sizeof(double));

        for( int i = 0; i < truncamiento; i++){//Este es el ciclo para toda la seleccion.
            const double g_bar = (poblacion + i)->eval - gx_sel;
            sum_evals += g_bar;
            for(unsigned int j = 0; j < n_pars; j++){
                const unsigned int k = idx_pars[j];
                *(medias + k) = *(medias + k) + (poblacion + i)->vars[ k ] * g_bar;
            }
        }

        for(unsigned int j = 0; j < n_pars; j++){
            *( medias + idx_pars[j] ) = *( medias + idx_pars[j] ) / sum_evals;
        }
        //Se calculan las varianzas:
        double tmp;
        memset(varianzas, 0, 4*sizeof(double));
        for( int i = 0; i < truncamiento; i++){
            const double g_bar = (poblacion + i)->eval - gx_sel;
            for(unsigned int j = 0; j < n_pars; j++){

                const unsigned int k = idx_pars[j];
                tmp = (poblacion+i)->vars[k] - *(medias + k);
                *(varianzas + k) = *(varianzas + k) + tmp*tmp * g_bar;
            }
        }

        for(unsigned int j = 0; j < n_pars; j++){
            *(varianzas + idx_pars[j] ) = *(varianzas + idx_pars[j] ) / (sum_evals);
        }
}



/*
    Metodo:     seleccionarPob
    Funcion:    Realiza la seleccion por truncamiento, redefiniendo el tetha para esta iteracion dentro del mismo apuntador de tetha.
*/
int OPTI_PARS::seleccionarPob(double *theta_t, const INDIV *poblacion){
    int truncamiento;

    for( truncamiento = (n_pob-1); truncamiento >= 0; truncamiento--){
        // Se trunca hasta el individuo con evaluacion en la funcion objetivo encima del theta anterior.
        if((poblacion + truncamiento)->eval > (*theta_t)){
            break;
        }
    }

    if( truncamiento == (n_pob-1) ){
        truncamiento = n_pob / 2;
    }
    //Se asigna el nuevo theta como la evaluacion del individuo en la posicion 'truncamiento'
    *theta_t = (poblacion + truncamiento)->eval;

    return truncamiento;
}



/*
    Metodo:        BUMDA
    Funcion:    Utiliza el algoritmo Boltzmann Univariate Marginal Distribution para encotnrar los parametros automaticamente.
*/
void OPTI_PARS::BUMDA(){
    INDIV *poblacion = new INDIV [n_pob + 1];
    // Poner los valores del elite por defecto:
    for( int i = 0; i < n_pob; i++){
        memcpy( poblacion + i, mi_elite, sizeof(INDIV));
    }

    double medias[4] = {0.0, 0.0, 0.0, 0.0};
    double varianzas[4] = {0.0, 0.0, 0.0, 0.0};

    // Se cuantan los parametros activos:
    n_pars = 0;
    for( int p = 0; p < 4; p++){
        if( pars_optim[p] ){
            idx_pars[n_pars] = p;
            n_pars ++;
        }
    }

    if(semilla){
        delete semilla;
    }
    semilla = initSeed(0);

    // Inicia el algoritmo BUMDA
    int k = 1, truncamiento = n_pob-1;
    bool procesar = true;

    //// Generar la primer poblacion:
    generarPobInicial(poblacion);
    qsort((void*)poblacion, n_pob, sizeof(INDIV), compIndiv); // El mejor fitness queda en la posicion 0

    //Se define el theta en el tiempo 0 como el minimo de la poblacion inicial.
    double tetha_t = (poblacion + n_pob - 1)->eval;
    char mensaje_iter[512] = "Best fit: X.XXXXXXXX, L: XX.XXX , T: XX.XXX, Sigma: XX.XXX, K: XXX.X :: Laste eval: X.XXXX, L: XX.XXX , T: XX.XXX, Sigma: XX.XXX, K: XXX.X :: en XXXXX.XXXX s.\n";
    escribirLog( "iteration\tbest_fit\tT\tL\tSigma\tK\telapsed_time\n" );

    do{
        // Guardar el elite como el mejor de los individuos en la posicion 'n_pob' del arreglo:
        memcpy(poblacion + n_pob, poblacion, sizeof(INDIV));

#if defined(_WIN32) || defined(_WIN64)
		sprintf_s(mensaje_iter, 512 * sizeof(char), "%i\t%1.8f\t%1.3f\t%1.3f\t%2.3f\t%3.0f\n", k, poblacion->eval, poblacion->vars[PAR_T], poblacion->vars[PAR_L], poblacion->vars[PAR_SIGMA], poblacion->vars[PAR_K]);
#else
        sprintf( mensaje_iter, "%i\t%1.8f\t%1.3f\t%1.3f\t%2.3f\t%3.0f\n", k, poblacion->eval, poblacion->vars[PAR_T], poblacion->vars[PAR_L], poblacion->vars[PAR_SIGMA], poblacion->vars[PAR_K]);
#endif
        escribirLog( mensaje_iter );
        barraProgreso( k, max_iters);

        calcularPars(poblacion, truncamiento, medias, varianzas);
        generarPob(medias, varianzas, poblacion);
        qsort((void*)poblacion, n_pob+1, sizeof(INDIV), compIndiv); // El mejor fitness queda en la posicion 0 ya comparando el elite

        //En truncamiento se guarda el indice del individuo a partir del cual se trunca la poblacion.
        truncamiento = seleccionarPob(&tetha_t, poblacion);
        k++;

        //// Verificar condiciones de paro:
		unsigned int condicion = 0;
        for(unsigned int j = 0; j < n_pars; j++){
            const unsigned int v = idx_pars[j];
            DEB_MSG("var[" << v << "] = " << *(varianzas + v) << "/" << *(min_vars + v));
            condicion += ( *(varianzas + v) <= *(min_vars + v));
        }
        procesar = (k < max_iters) && (condicion < n_pars);

    }while(procesar);

    memcpy(mi_elite, poblacion, sizeof(INDIV));

#if defined(_WIN32) || defined(_WIN64)
	sprintf_s(mensaje_iter, 512 * sizeof(char), "%i\t%1.8f\t%1.3f\t%1.3f\t%2.3f\t%3.0f\n", k, mi_elite->eval, mi_elite->vars[PAR_T], mi_elite->vars[PAR_L], mi_elite->vars[PAR_SIGMA], mi_elite->vars[PAR_K]);
#else
    sprintf( mensaje_iter, "%i\t%1.8f\t%1.3f\t%1.3f\t%2.3f\t%3.0f\n", k, mi_elite->eval, mi_elite->vars[PAR_T], mi_elite->vars[PAR_L], mi_elite->vars[PAR_SIGMA], mi_elite->vars[PAR_K]);
#endif
    escribirLog( mensaje_iter );
    barraProgreso( max_iters, max_iters);

    delete [] poblacion;
}




// UMDA ----- UMDA ----- UMDA ----- UMDA ----- UMDA ----- UMDA ----- UMDA ----- UMDA ----- UMDA ----- UMDA ----- UMDA ----- UMDA ----- UMDA ----- UMDA ----- UMDA ----- UMDA ----- UMDA ----- UMDA ----- UMDA ----- UMDA ----- UMDA ----- UMDA ----- UMDA ----- UMDA ----- UMDA ----- UMDA ----- UMDA ----- UMDA ----- UMDA ----- UMDA ----- UMDA ----- UMDA ----- UMDA ----- UMDA

/*
    Metodo:     generar_pob
    Funcion:   Genera una nueva poblacion en bsae al muestreo de cada atributo.
*/
void OPTI_PARS::generarPob(INDIV *poblacion, const double *probs, const double *deltas_var){
    for(int i = 0; i < n_pob; i++){
        // Se muestrean los genes para cada parametro:
        unsigned int bits_recorridos = 0;
        for(unsigned int j = 0; j < n_pars; j++){
            const unsigned int k = idx_pars[j];
            double cadena_val = 0;
            double pow_2 = 1.0;
            for( int b = 0; b < (int)min_vars[k]; b++, bits_recorridos++, pow_2*=2.0){
                (poblacion + i)->cadena[ bits_recorridos ] = (HybTaus(0.0, 1.0, semilla) <= *(probs + bits_recorridos)) ? 1 : 0;
                cadena_val += ((poblacion + i)->cadena[ bits_recorridos ] ? 1.0 : 0.0) * (double)pow_2;
            }
            // Asignar el valor segun la cadena formada para el atributo 'k':
            (poblacion + i)->vars[k] = *(deltas_var + k) * cadena_val + lim_inf[k];
        }

        switch( fitness_elegido ){
            case ROC:
                (poblacion + i)->eval = fitnessROC( poblacion + i);
                break;
            case CORCON:
                (poblacion + i)->eval = fitnessCorCon( poblacion + i);
                break;
        }
    }
}




/*
    Metodo:     calcular_pars
    Funcion:    Calcula las probabilidades de cada distribucion marginal.
*/
void OPTI_PARS::calcularPars(const INDIV *poblacion, const int n_bits, const int truncamiento, double *probs){
    // Se estiman las probabilidades para cada bit de todos los parametro:
    memset(probs, 0, n_bits*sizeof(double));
    const double fraccion = 1.0 / (double)truncamiento;

    for(int i = 0; i < truncamiento; i++){
        for( int b = 0; b < n_bits; b++){
            *(probs + b) = *(probs + b) + ((poblacion + i)->cadena[b] ? fraccion : 0.0);
        }
    }
}





/*
    Metodo:     UMDA
    Funcion:    Utiliza el algoritmo Univariate Marginal Distribution para encotnrar los parametros automaticamente.
*/
void OPTI_PARS::UMDA(){

    INDIV *poblacion = new INDIV [n_pob];
    // Poner los valores del elite por defecto:
    for( int i = 0; i < n_pob; i++){
        memcpy( poblacion + i, mi_elite, sizeof(INDIV));
    }

    // Se cuentan los parametros activos:
    n_pars = 0;
    double deltas_vars[4] = {0.0, 0.0, 0.0, 0.0};

    unsigned int n_bits = 0;
    for( int p = 0; p < 4; p++){
        if( pars_optim[p] ){
            idx_pars[ n_pars ] = p;

            // Calcular el tamano de paso
            const double max_bit_val = pow(2, min_vars[p]) - 1;
            deltas_vars[ p ] = (lim_sup[ p ] - lim_inf[ p ]) / max_bit_val;

			DEB_MSG("[" << p << "] delta: " << deltas_vars[p]);


            n_bits += (unsigned int)min_vars[ p ];

            n_pars ++;
        }
    }

    double *probs = new double [n_bits];
    // Inicializar las probabilidades en 0.5:
    for(unsigned int b = 0; b < n_bits; b++){
        probs[b] = 0.5;
    }

    if(semilla){
        delete semilla;
    }
    semilla = initSeed(0);

    // Inicia el algoritmo UMDA
    int k = 0, truncamiento = (int)( (double)n_pob * 0.6);
    bool procesar = true;


    char mensaje_iter[512] = "Best fit: X.XXXXXXXX, L: XX.XXX , T: XX.XXX, Sigma: XX.XXX, K: XXX.X :: Laste eval: X.XXXX, L: XX.XXX , T: XX.XXX, Sigma: XX.XXX, K: XXX.X :: en XXXXX.XXXX s.\n";
    escribirLog( "iteration\tbest_fit\tT\tL\tSigma\tK\telapsed_time\n" );
	
    mi_elite->eval = 0.0;

    do{
        generarPob(poblacion, probs, deltas_vars);
        qsort((void*)poblacion, n_pob, sizeof(INDIV), compIndiv);
        calcularPars(poblacion, n_bits, truncamiento, probs);

        // Guardar el elite como el mejor de los individuos en la posicion 'n_pob' del arreglo:
        if( mi_elite->eval < poblacion->eval ){
            memcpy(mi_elite, poblacion, sizeof(INDIV));
        }
        k++;

        //// Verificar condiciones de paro:
        procesar = (k < max_iters);

#if defined(_WIN32) || defined(_WIN64)
		sprintf_s(mensaje_iter, 512 * sizeof(char), "%i\t%1.8f\t%1.3f\t%1.3f\t%2.3f\t%3.0f\n", k, mi_elite->eval, mi_elite->vars[PAR_T], mi_elite->vars[PAR_L], mi_elite->vars[PAR_SIGMA], mi_elite->vars[PAR_K]);
#else
        sprintf( mensaje_iter, "%i\t%1.8f\t%1.3f\t%1.3f\t%2.3f\t%3.0f\n", k, mi_elite->eval, mi_elite->vars[PAR_T], mi_elite->vars[PAR_L], mi_elite->vars[PAR_SIGMA], mi_elite->vars[PAR_K]);
#endif
        escribirLog( mensaje_iter );
        barraProgreso( k, max_iters);
    }while(procesar);

#if defined(_WIN32) || defined(_WIN64)
	sprintf_s(mensaje_iter, 512 * sizeof(char), "%i\t%1.8f\t%1.3f\t%1.3f\t%2.3f\t%3.0f\n", k, mi_elite->eval, mi_elite->vars[PAR_T], mi_elite->vars[PAR_L], mi_elite->vars[PAR_SIGMA], mi_elite->vars[PAR_K]);
#else
    sprintf( mensaje_iter, "%i\t%1.8f\t%1.3f\t%1.3f\t%2.3f\t%3.0f\n", k, mi_elite->eval, mi_elite->vars[PAR_T], mi_elite->vars[PAR_L], mi_elite->vars[PAR_SIGMA], mi_elite->vars[PAR_K]);
#endif
    escribirLog( mensaje_iter );
    barraProgreso( k, max_iters);

    delete [] probs;
    delete [] poblacion;
}




// GA ----- GA ----- GA ----- GA ----- GA ----- GA ----- GA ----- GA ----- GA ----- GA ----- GA ----- GA ----- GA ----- GA ----- GA ----- GA ----- GA ----- GA ----- GA ----- GA ----- GA ----- GA ----- GA ----- GA ----- GA ----- GA ----- GA ----- GA ----- GA ----- GA ----- GA
/*
    Metodo:     generar_pob
    Funcion:   Genera una nueva poblacion en bsae al muestreo de cada atributo.
*/
double OPTI_PARS::generarPobInicial(INDIV *poblacion, const double *deltas_var){
    double sum_fitness = 0.0;

    for(int i = 0; i < n_pob; i++){
        // Se genera cada bit con la misma probabilidad de ser 0 o 1:

        unsigned int bits_recorridos = 0;
		for (unsigned int j = 0; j < n_pars; j++) {
			const unsigned int k = idx_pars[j];
			double pow_2 = 1.0;
			double cadena_val = 0.0;
			for (int b = 0; b < (int)min_vars[k]; b++, bits_recorridos++, pow_2*=2.0) {
				(poblacion + i)->cadena[bits_recorridos] = (HybTaus(0.0, 1.0, semilla) <= 0.5) ? 0 : 1;
				cadena_val += ((poblacion + i)->cadena[bits_recorridos] ? 1.0 : 0.0) * pow_2;
			}

			// Asignar el valor segun la cadena formada para el atributo 'k':
			(poblacion + i)->vars[k] = *(deltas_var + k) * cadena_val + lim_inf[k];
        }

        switch( fitness_elegido ){
            case ROC:
                (poblacion + i)->eval = fitnessROC( poblacion + i);
                break;
        }
        sum_fitness += (poblacion + i)->eval;
    }

    return sum_fitness;
}



/*
    Metodo:     selecPob
    Funcion:    Genera el grupo de individuos seleccionados.
*/
void OPTI_PARS::selecPob(INDIV* sel_grp, const INDIV* poblacion, double *fitness_acum, const double suma_fitness){
    // Calcular la funcion de distribucion de probabilidad F(x):
    *(fitness_acum + n_pob - 1) = (poblacion + n_pob - 1)->eval / suma_fitness;
    for( int i = n_pob-2; i >= 0; i--){
        *(fitness_acum + i) = *(fitness_acum + i + 1) + (poblacion + i)->eval / suma_fitness;
    }

    // Armar el grupo de seleccion:
    for( int i = 0; i < seleccion; i++){
        const double ruleta = HybTaus(0.0, 1.0, semilla);
        int sel_pos = n_pob;
        while( sel_pos > 0 ){
            sel_pos--;
            if( *(fitness_acum + sel_pos) >= ruleta ){
                break;
            }
        }
		memcpy( sel_grp + i, poblacion + sel_pos, sizeof(INDIV) );
    }
}



/*
    Metodo:     cruzaPob
    Funcion:    Calcula la siguiente generacion a partir de la cruza multipunto de los padres seleccionados, ademas se realiza una mutacion sobre un porcentage de los individuos.
*/
void OPTI_PARS::cruzaPob(INDIV* cruza, const INDIV* sel_grp, const unsigned int n_bits){

    // Generar dos nuevos individuos a partir de dos padres, el proceso se repite hasta completar la fraccion 'seleccion'.
    for(int i = 0; i < seleccion; i+=2){

        int padre_1 = i;
        int padre_2 = i + 1;
		unsigned int bits_ini = 0;
		
        // Realizar cortes hasta terminar con la secuencia de los individuos:
        while(bits_ini < n_bits){
            const unsigned int n_bits_cpy = (unsigned int)HybTaus(1.0, (double)(n_bits - bits_ini), semilla);
            memcpy( (cruza + i  )->cadena + bits_ini, (sel_grp + padre_1)->cadena + bits_ini, n_bits_cpy*sizeof(unsigned char));
            memcpy( (cruza + i+1)->cadena + bits_ini, (sel_grp + padre_2)->cadena + bits_ini, n_bits_cpy*sizeof(unsigned char));

            int padre_swap = padre_1;
            padre_1 = padre_2;
            padre_2 = padre_swap;

            bits_ini += n_bits_cpy;
        }

        // Realizar mutacion sobre los hijos generados con cierta probabilidad: 'prob_mutacion'
        // Mutar el hijo 1:
        if( HybTaus(0.0, 1.0, semilla) <= prob_mutacion ){
            // Seleccionar un gen a mutar:
            const unsigned int gen = (unsigned int)(HybTaus(0.0, (double)n_bits - 1e-8, semilla));
            (cruza + i)->cadena[gen] = 1 - (cruza + i)->cadena[gen];
        }

        // Mutar el hijo 2:
        if( HybTaus(0.0, 1.0, semilla) <= prob_mutacion ){
            // Seleccionar un gen a mutar:
            const unsigned int gen = (unsigned int)(HybTaus(0.0, (double)n_bits - 1e-8, semilla));
            (cruza + i + 1)->cadena[gen] = 1 - (cruza + i + 1)->cadena[gen];
        }
    }
}




/*
    Metodo:     generar_pob
    Funcion:   Genera una nueva poblacion en base al muestreo de cada atributo.
*/
double OPTI_PARS::generarPob(INDIV *poblacion, const INDIV *cruza, const INDIV *sel_grp, const double *deltas_var){
    // Construir la nueva poblacion a partir de los hijos generados y la seleccion previa:
    double suma_fitness = 0.0;

    for( int i = 0; i < seleccion; i++){
        memcpy( poblacion + i, cruza + i, sizeof(INDIV) );

        unsigned int bits_recorridos = 0;
        for(unsigned int j = 0; j < n_pars; j++){
            const unsigned int k = idx_pars[j];
            double pow_2 = 1.0;
            double cadena_val = 0.0;
            for( int b = 0; b < (int)min_vars[k]; b++, bits_recorridos++, pow_2*=2.0){
                (poblacion + i)->cadena[bits_recorridos] = (HybTaus(0.0, 1.0, semilla) <= 0.5) ? 0 : 1;
                cadena_val += ((poblacion + i)->cadena[bits_recorridos] ? 1.0 : 0.0) * (double)pow_2;
            }
            // Asignar el valor segun la cadena formada para el atributo 'k':
            (poblacion + i)->vars[ k ] = *(deltas_var + k) * cadena_val + lim_inf[k];
        }
        switch( fitness_elegido ){
        case ROC:
            (poblacion + i)->eval = fitnessROC( poblacion + i);
            break;
        case CORCON:
            (poblacion + i)->eval = fitnessCorCon( poblacion + i);
            break;
        }
        suma_fitness += (poblacion + i)->eval;
    }

    for( int i = seleccion; i < n_pob; i++){
        const unsigned int sel_pos = (unsigned int)(HybTaus(0.0, (double)seleccion - 1e-8, semilla));
        memcpy( poblacion + i, sel_grp + sel_pos, sizeof(INDIV) );
        suma_fitness += (poblacion + i)->eval;
    }

    return suma_fitness;
}




/*
    Metodo:     GA
    Funcion:    Utiliza el algoritmo genetico para encotnrar los parametros automaticamente.
*/
void OPTI_PARS::GA(){
    INDIV *poblacion = new INDIV [n_pob];
    INDIV *sel_grp = new INDIV [seleccion];
    INDIV *cruza = new INDIV [seleccion];

    // Poner los valores del elite por defecto:
    for( int i = 0; i < seleccion; i++){
        memcpy(poblacion + i, mi_elite, sizeof(INDIV));
        memcpy(cruza + i, mi_elite, sizeof(INDIV));
    }
    for( int i = seleccion; i < n_pob; i++){
        memcpy(poblacion + i, mi_elite, sizeof(INDIV));
    }

    // Se cuentan los parametros activos:
    n_pars = 0;
    double deltas_var[4] = {0.0, 0.0, 0.0, 0.0};

    unsigned int n_bits = 0;
    for( int p = 0; p < 4; p++){
        if( pars_optim[p] ){
            idx_pars[ n_pars ] = p;

            // Calcular el tamano de paso
            const double max_bit_val = pow(2, min_vars[p]) - 1;
            deltas_var[ p ] = (lim_sup[ p ] - lim_inf[ p ]) / max_bit_val;
            n_bits += (unsigned int)min_vars[ p ];
            n_pars ++;
        }
    }

    if(semilla){
        delete semilla;
    }
    semilla = initSeed(0);

    // Inicia el algoritmo GA
    int k = 1;
    bool procesar = true;

    //// Generar la primer poblacion:
    double suma_fitness = generarPobInicial(poblacion, deltas_var);
    qsort((void*)poblacion, n_pob, sizeof(INDIV), compIndiv);

    double *fitness_acum = new double [n_pob];

    char mensaje_iter[512] = "Best fit: X.XXXXXXXX, L: XX.XXX , T: XX.XXX, Sigma: XX.XXX, K: XXX.X :: Laste eval: X.XXXX, L: XX.XXX , T: XX.XXX, Sigma: XX.XXX, K: XXX.X :: en XXXXX.XXXX s.\n";
    escribirLog( "iteration\tbest_fit\tT\tL\tSigma\tK\telapsed_time\n" );
    mi_elite->eval = 0.0;

    do{
        selecPob(sel_grp, poblacion, fitness_acum, suma_fitness);
        cruzaPob(cruza, sel_grp, n_bits);
        suma_fitness = generarPob(poblacion, cruza, sel_grp, deltas_var);

		qsort((void*)poblacion, n_pob, sizeof(INDIV), compIndiv);
        // Guardar el elite como el mejor de los individuos en la posicion 'n_pob' del arreglo:
        if( mi_elite->eval < poblacion->eval ){
            memcpy(mi_elite, poblacion, sizeof(INDIV));
        }

#if defined(_WIN32) || defined(_WIN64)
		sprintf_s(mensaje_iter, 512 * sizeof(char), "%i\t%1.8f\t%1.3f\t%1.3f\t%2.3f\t%3.0f\n", k, mi_elite->eval, mi_elite->vars[PAR_T], mi_elite->vars[PAR_L], mi_elite->vars[PAR_SIGMA], mi_elite->vars[PAR_K]);
#else
        sprintf( mensaje_iter, "%i\t%1.8f\t%1.3f\t%1.3f\t%2.3f\t%3.0f\n", k, mi_elite->eval, mi_elite->vars[PAR_T], mi_elite->vars[PAR_L], mi_elite->vars[PAR_SIGMA], mi_elite->vars[PAR_K]);
#endif
        escribirLog( mensaje_iter );
        barraProgreso( k, max_iters);
        k++;

        //// Verificar condiciones de paro:
        procesar = (k < max_iters);


    }while(procesar);

#if defined(_WIN32) || defined(_WIN64)
	sprintf_s(mensaje_iter, 512 * sizeof(char), "%i\t%1.8f\t%1.3f\t%1.3f\t%2.3f\t%3.0f\n", k, mi_elite->eval, mi_elite->vars[PAR_T], mi_elite->vars[PAR_L], mi_elite->vars[PAR_SIGMA], mi_elite->vars[PAR_K]);
#else
    sprintf( mensaje_iter, "%i\t%1.8f\t%1.3f\t%1.3f\t%2.3f\t%3.0f\n", k, mi_elite->eval, mi_elite->vars[PAR_T], mi_elite->vars[PAR_L], mi_elite->vars[PAR_SIGMA], mi_elite->vars[PAR_K]);
#endif
    escribirLog( mensaje_iter );
    barraProgreso( k, max_iters);

    delete [] cruza;
    delete [] sel_grp;
    delete [] fitness_acum;
    delete [] poblacion;
}



// DE ----- DE ----- DE ----- DE ----- DE ----- DE ----- DE ----- DE ----- DE ----- DE ----- DE ----- DE ----- DE ----- DE ----- DE ----- DE ----- DE ----- DE ----- DE ----- DE ----- DE ----- DE ----- DE ----- DE ----- DE ----- DE ----- DE ----- DE ----- DE ----- DE ----- DE
/*
    Metodo:     diferenciarPoblacion
    Funcion:    Generar una nueva poblacion a partir de la poblacion anterior.

*/
void OPTI_PARS::diferenciarPoblacion( INDIV* poblacion, INDIV* pob_base ){

    // Crear una copia de la poblacion anterior a partir de la cual se genera la nueva poblacion:
    memcpy( pob_base, poblacion, n_pob * sizeof(INDIV) );

    for(int i = 0; i < n_pob; i++){
        // Diferenciar el individuo unicamente si se cumple la probabilidad de cruza:
        if( HybTaus(0.0, 1.0, semilla) <= prob_cruza ){

            // Generar tres indices aleatorios diferentes (excluyendo el indice actual):
            int idx_1 = (int)round(HybTaus(-0.5, (double)(n_pob-2) + 0.49, semilla));
            idx_1 += (idx_1 >= i) ? 1 : 0;

            int idx_2 = idx_1;
            do{
                idx_2 = (int)round(HybTaus(-0.5, (double)(n_pob-2) + 0.49, semilla));
                idx_2 += (idx_2 >= i) ? 1 : 0;
            }while( idx_2 == idx_1 );

            int idx_3 = idx_1;
            do{
                idx_3 = (int)round(HybTaus(-0.5, (double)(n_pob-2) + 0.49, semilla));
                idx_3 += (idx_3 >= i) ? 1 : 0;
            }while( (idx_3 == idx_1) || (idx_3 == idx_2) );

            DEB_MSG("[" << i << "]: {" << idx_1 << ", " << idx_2 << ", " << idx_3 << "}");

            // Diferenciar todas las variables del individuo:
            for(unsigned int j = 0; j < n_pars; j++){
                const unsigned int k = idx_pars[j];
                double val_gen = (pob_base + idx_1)->vars[k] + prob_mutacion*((pob_base + idx_2)->vars[k] - (pob_base + idx_3)->vars[k]);

                val_gen = (val_gen <= lim_sup[ k ]) ? ((val_gen >= lim_inf[ k ]) ? val_gen : lim_inf[ k ]) : lim_sup[ k ];
                (poblacion + i)->vars[ k ] = val_gen;
            }

            switch( fitness_elegido ){
                case ROC:
                    (poblacion + i)->eval = fitnessROC( poblacion + i);
                    break;

                case CORCON:
                    (poblacion + i)->eval = fitnessCorCon( poblacion + i);
                    break;
            }

            // Verificar si el individuo nuevo mejora con respecto del anterior:
            if( (poblacion + i)->eval <= (pob_base + i)->eval ){ // El individuo nuevo NO mejora la solucion del anterior:
                memcpy( poblacion + i, pob_base + i, sizeof(INDIV) );
            }
        }
    }

}


/*
    Metodo:     DE
    Funcion:    Differential Evolution algorithm for filter parameters optimization.

*/
void OPTI_PARS::DE(){
    INDIV *poblacion = new INDIV [n_pob];
    INDIV *pob_anterior = new INDIV [n_pob];

    // Poner los valores del elite por defecto:
    for( int i = 0; i < n_pob; i++){
        memcpy( poblacion + i, mi_elite, sizeof(INDIV));
    }

    // Se cuentan los parametros activos:
    n_pars = 0;
    for( int p = 0; p < 4; p++){
        if( pars_optim[p] ){
            idx_pars[n_pars] = p;
            n_pars ++;
        }
    }

    if(semilla){
        delete semilla;
    }
    semilla = initSeed(0);

    // Inicia el algoritmo DE
    int k = 1;
    bool procesar = true;

    //// Generar la primer poblacion:
    generarPobInicial(poblacion);
    qsort((void*)poblacion, n_pob, sizeof(INDIV), compIndiv); // El mejor fitness queda en la posicion 0

    //Se define el theta en el tiempo 0 como el minimo de la poblacion inicial.
    char mensaje_iter[512] = "Best fit: X.XXXXXXXX, L: XX.XXX , T: XX.XXX, Sigma: XX.XXX, K: XXX.X :: Laste eval: X.XXXX, L: XX.XXX , T: XX.XXX, Sigma: XX.XXX, K: XXX.X :: en XXXXX.XXXX s.\n";
    escribirLog( "iteration\tbest_fit\tT\tL\tSigma\tK\telapsed_time\n" );

    do{
#if defined(_WIN32) || defined(_WIN64)
		sprintf_s(mensaje_iter, 512 * sizeof(char), "%i\t%1.8f\t%1.3f\t%1.3f\t%2.3f\t%3.0f\n", k, poblacion->eval, poblacion->vars[PAR_T], poblacion->vars[PAR_L], poblacion->vars[PAR_SIGMA], poblacion->vars[PAR_K]);
#else
        sprintf( mensaje_iter, "%i\t%1.8f\t%1.3f\t%1.3f\t%2.3f\t%3.0f\n", k, poblacion->eval, poblacion->vars[PAR_T], poblacion->vars[PAR_L], poblacion->vars[PAR_SIGMA], poblacion->vars[PAR_K]);
#endif
        escribirLog( mensaje_iter );
        barraProgreso( k, max_iters);

        diferenciarPoblacion( poblacion, pob_anterior );
        qsort((void*)poblacion, n_pob, sizeof(INDIV), compIndiv); // El mejor fitness queda en la posicion 0 ya comparando el elite

        k++;

        //// Verificar condiciones de paro:
        procesar = (k < max_iters);

    }while(procesar);

    memcpy(mi_elite, poblacion, sizeof(INDIV));

#if defined(_WIN32) || defined(_WIN64)
	sprintf_s(mensaje_iter, 512 * sizeof(char), "%i\t%1.8f\t%1.3f\t%1.3f\t%2.3f\t%3.0f\n", k, mi_elite->eval, mi_elite->vars[PAR_T], mi_elite->vars[PAR_L], mi_elite->vars[PAR_SIGMA], mi_elite->vars[PAR_K]);
#else
    sprintf( mensaje_iter, "%i\t%1.8f\t%1.3f\t%1.3f\t%2.3f\t%3.0f\n", k, mi_elite->eval, mi_elite->vars[PAR_T], mi_elite->vars[PAR_L], mi_elite->vars[PAR_SIGMA], mi_elite->vars[PAR_K]);
#endif
    escribirLog( mensaje_iter );
    barraProgreso( max_iters, max_iters);

    delete [] poblacion;
    delete [] pob_anterior;
}





/*
    Metodo:     busquedaExhaustiva
    Funcion:    Realiza una busqueda exhaustiva en los limites de los parametros
*/
void OPTI_PARS::busquedaExhaustiva(){
    // Determinar cuales parametros se van a optimizar:
    n_pars = 0;

    INDIV *test = new INDIV;
    memcpy( test, mi_elite, sizeof(INDIV));

    unsigned int max_bit_val[4] = {0, 0, 0, 0};

    unsigned int n_bits = 1;
    for( int p = 0; p < 4; p++){
        if( pars_optim[p] ){
            idx_pars[ n_pars ] = p;

            // Determinar cuantas particiones tiene cada variable:
            max_bit_val[ p ] = (unsigned int) ((lim_sup[ p ] - lim_inf[ p ] + 1e-12) / min_vars[p]);
            n_bits *= max_bit_val[ p ] + 1;
            n_pars ++;
        }
    }


    int idx = 0;


    char mensaje_iter[512] = "Best fit: X.XXXXXXXX, L: XX.XXX , T: XX.XXX, Sigma: XX.XXX, K: XXX.X :: Laste eval: X.XXXX, L: XX.XXX , T: XX.XXX, Sigma: XX.XXX, K: XXX.X :: en XXXXX.XXXX s.\n";
    escribirLog( "iteration\tbest_fit\tT\tL\tSigma\tK\telapsed_time\n" );

    for(unsigned int i_L = 0; i_L <= max_bit_val[PAR_L]; i_L++){
        test->vars[PAR_L] = (double)i_L * min_vars[PAR_L] + lim_inf[PAR_L];

        for(unsigned int i_T = 0; i_T <= max_bit_val[PAR_T]; i_T++){
            test->vars[PAR_T] = (double)i_T * min_vars[PAR_T] + lim_inf[PAR_T];

            for(unsigned int i_K = 0; i_K <= max_bit_val[PAR_K]; i_K++){
                test->vars[PAR_K] = (double)i_K * min_vars[PAR_K] + lim_inf[PAR_K];

                for(unsigned int i_Sigma = 0; i_Sigma <= max_bit_val[PAR_SIGMA]; i_Sigma++){
                    test->vars[PAR_SIGMA] = (double)i_Sigma * min_vars[PAR_SIGMA] + lim_inf[PAR_SIGMA];

                    switch( fitness_elegido ){
                        case ROC:
                            test->eval = fitnessROC( test);
                            break;
                        case CORCON:
                            test->eval = fitnessCorCon( test);
                            break;
                    }

                    if( (test->eval) > (mi_elite->eval) ){
                        memcpy( mi_elite, test, sizeof(INDIV));
                    }
                    idx++;


#if defined(_WIN32) || defined(_WIN64)
					sprintf_s(mensaje_iter, 512 * sizeof(char), "%i\t%1.8f\t%1.3f\t%1.3f\t%2.3f\t%3.0f\n", idx, mi_elite->eval, mi_elite->vars[PAR_T], mi_elite->vars[PAR_L], mi_elite->vars[PAR_SIGMA], mi_elite->vars[PAR_K]);
#else
                    sprintf( mensaje_iter, "%i\t%1.8f\t%1.3f\t%1.3f\t%2.3f\t%3.0f\t%5.4f\n", idx, mi_elite->eval, mi_elite->vars[PAR_T], mi_elite->vars[PAR_L], mi_elite->vars[PAR_SIGMA], mi_elite->vars[PAR_K], DIFTIME);
#endif
                    escribirLog( mensaje_iter );
                    barraProgreso( idx, n_bits);
                }
            }
        }
    }

#if defined(_WIN32) || defined(_WIN64)
	sprintf_s(mensaje_iter, 512 * sizeof(char), "%i\t%1.8f\t%1.3f\t%1.3f\t%2.3f\t%3.0f\n", idx, mi_elite->eval, mi_elite->vars[PAR_T], mi_elite->vars[PAR_L], mi_elite->vars[PAR_SIGMA], mi_elite->vars[PAR_K]);
#else
    sprintf( mensaje_iter, "%i\t%1.8f\t%1.3f\t%1.3f\t%2.3f\t%3.0f\n", idx, mi_elite->eval, mi_elite->vars[PAR_T], mi_elite->vars[PAR_L], mi_elite->vars[PAR_SIGMA], mi_elite->vars[PAR_K]);
#endif
    escribirLog( mensaje_iter );
    barraProgreso( idx, n_bits);

    delete test;
}

















/************************************************************************************************************
* OPTI_PARS::PUBLIC                                                                                         *
*                                                                                                           *
* FUNCTION NAME: setInputBase                                                                               *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT                  TYPE                      I/O  DESCRIPTION                                      *
* new_img_base              std::vector<IMGCONT>*      I   The pointer to a vector of IMGCONT objects.      *
*                                                                                                           *
* RETURNS:                                                                                                  *
* Nothing.                                                                                                  *
*                                                                                                           *
************************************************************************************************************/
void OPTI_PARS::setInputBase(std::vector<IMGCONT>* new_img_base)
{
	this->setInputFilterBase(new_img_base);
}













/************************************************************************************************************
* OPTI_PARS::PUBLIC                                                                                         *
*                                                                                                           *
* FUNCTION NAME: setInputMask                                                                               *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT                  TYPE                      I/O  DESCRIPTION                                      *
* new_img_base_mask         std::vector<IMGCONT>*      I   The pointer to a vector of IMGCONT objects.      *
*                                                                                                           *
* RETURNS:                                                                                                  *
* Nothing.                                                                                                  *
*                                                                                                           *
************************************************************************************************************/
void OPTI_PARS::setInputMask(std::vector<IMGCONT>* new_img_mask)
{
	this->setInputFilterMask(new_img_mask);
	this->setInputPerformanceMask(new_img_mask);
}













/************************************************************************************************************
* OPTI_PARS::PUBLIC                                                                                         *
*                                                                                                           *
* FUNCTION NAME: setInputGroundtruth                                                                        *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT                  TYPE                      I/O  DESCRIPTION                                      *
* new_img_groundtruth       std::vector<IMGCONT>*      I   The pointer to a vector of IMGCONT objects.      *
*                                                                                                           *
* RETURNS:                                                                                                  *
* Nothing.                                                                                                  *
*                                                                                                           *
************************************************************************************************************/
void OPTI_PARS::setInputGroundtruth(std::vector<IMGCONT>* new_img_groundtruth)
{
	this->setInputFilterGroundtruth(new_img_groundtruth);
	this->setInputPerformanceGroundtruth(new_img_groundtruth);
}













/************************************************************************************************************
* OPTI_PARS::PUBLIC                                                                                         *
*                                                                                                           *
* FUNCTION NAME: setInputResponse                                                                           *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT                  TYPE                      I/O  DESCRIPTION                                      *
* new_img_response          std::vector<IMGCONT>*      I   The pointer to a vector of IMGCONT objects.      *
*                                                                                                           *
* RETURNS:                                                                                                  *
* Nothing.                                                                                                  *
*                                                                                                           *
************************************************************************************************************/
void OPTI_PARS::setInputResponse(std::vector<IMGCONT>* new_img_response)
{
	this->setInputFilterResponse(new_img_response);
	this->setInputPerformanceResponse(new_img_response);
}












/************************************************************************************************************
* OPTI_PARS::PUBLIC                                                                                         *
*                                                                                                           *
* FUNCTION NAME: setInputThreshold                                                                          *
*                                                                                                           *
* ARGUMENTS:                                                                                                *
* ARGUMENT                   TYPE                      I/O  DESCRIPTION                                     *
* new_img_response_threshold std::vector<IMGCONT>*      I   The pointer to a vector of IMGCONT objects.     *
*                                                                                                           *
* RETURNS:                                                                                                  *
* Nothing.                                                                                                  *
*                                                                                                           *
************************************************************************************************************/
void OPTI_PARS::setInputThreshold(std::vector<IMGCONT>* new_img_response_threshold)
{
	this->setInputPerformanceThreshold(new_img_response_threshold);
}