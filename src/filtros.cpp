/******
    CENTRO DE INVESTIGACION EN MATEMATICAS
    MAESTRIA EN COMPUTACION Y MATEMATICAS INDUSTRIALES

    FERNANDO CERVANTES SANCHEZ
   NOV - 2015
*****/


#include "filtros.h"

/*  Metodo: escribirLog

    Funcion: Escribe un mensaje en el log.
*/
void FILTROS::escribirLog( const char *mensaje ){
    if(mi_fplog){
        fprintf(mi_fplog, "%s", mensaje);
    }else{
        std::cout << mensaje << std::endl;
        fflush(stdout);
    }
}


/*  Metodo: barraProgreso

    Funcion: Actualiza la barra de progreso.
*/
void FILTROS::barraProgreso( const int avance, const int milestones ){
    /// Limpiar el resto de la linea:
    int max_ancho = 100;
    for( int i = 0; i <= max_ancho; i++){
        printf("\r");
    }
    printf(COLOR_BACK_RED "[");
    int avance_milestones = (int)((double)max_ancho * (double)avance / (double)milestones);
    for( int i = 0; i < avance_milestones; i++){
        printf(COLOR_BACK_GREEN " ");
    }
    for( int i = avance_milestones; i < max_ancho; i++){
        printf(COLOR_BACK_CYAN " ");
    }
    printf(COLOR_BACK_RED "]" COLOR_NORMAL);
    fflush(stdout);
}




/*  Metodo: interpolacion (Lineal)
    Funcion: Interpola el valor del pixel destino en base a los 4 pixeles origen.
*/
inline double FILTROS::interpolacion(const double *pix, const int j, const int i, const double x, const double y, const int mis_rens, const int mis_cols){
    double intensidad = 0.0;

    if( j >= 0 && j < mis_cols ){
        if( i >= 0 && i < mis_rens ){
            intensidad += *(pix +   i  *mis_cols +   j  ) * (1.0 - x) * (1.0 - y);
        }
        if( i >= -1  && i < (mis_rens-1) ){
            intensidad += *(pix + (i+1)*mis_cols +   j  ) * (1.0 - x) * (   y   );
        }
    }
    if( j >= -1 && j < (mis_cols-1)){
        if( i >=  0 && i < mis_rens ){
            intensidad += *(pix +   i  *mis_cols + (j+1)) * (   x   ) * (1.0 - y);
        }
        if( i >= -1 && i < (mis_rens-1) ){
            intensidad += *(pix + (i+1)*mis_cols + (j+1)) * (   x   ) * (   y   );
        }
    }

    return intensidad;
}




/*  Metodo: rotarImg
    Funcion: Rota una imagen almacenada como intensidad de 0 a 1 en un angulo theta.
*/
void FILTROS::rotarImg( const double *org, double *rot, const double ctheta, const double stheta, const int mis_rens, const int mis_cols, const int org_rens, const int org_cols){

    const double mitad_x = (double)(mis_cols - 1) / 2.0;
    const double mitad_y = (double)(mis_rens - 1) / 2.0;
    const double mitad_orgx = (double)(org_cols - 1) / 2.0;
    const double mitad_orgy = (double)(org_rens - 1) / 2.0;


    double x, y;
    for(int i = 0; i < (mis_rens-1); i++){
        for(int j = 0; j < (mis_cols-1); j++){
            x = ((double)j-mitad_x)*ctheta + ((double)i-mitad_y)*stheta + mitad_orgx;
            y =-((double)j-mitad_x)*stheta + ((double)i-mitad_y)*ctheta + mitad_orgy;
            const int flx = (int)x;
            const int fly = (int)y;
            const double delta_x = x - (double)flx;
            const double delta_y = y - (double)fly;

            *(rot + i*mis_cols + j) = interpolacion(org, flx, fly, delta_x, delta_y, org_rens, org_cols);
        }
    }
}




/*//        //      //      //      //      //      //        //      //      //      //      //      //        //      //      //      //      //      //        //      //      //      //      //      //
 *
 *              M E T O D O S       P U B L I C O S     C O N S T R U C T O R E S / D E S T R U C T O R E S / S E T / G E T
 *
 * //        //      //      //      //      //      //        //      //      //      //      //      //        //      //      //      //      //      //        //      //      //      //      //      //
*/
/*  Metodo: setEvoMet
    Funcion: Define el metodo evolutivo que se usara para establecer los mejores aprametros del filtro.
*/
void FILTROS::setEvoMet(const EVO_MET evo_met){
    metodo_elegido = evo_met;
}



/*  Metodo: setEvoMetPar
    Funcion: Define elparametro para el metodo de optimizacion.
*/
void FILTROS::setEvoMetPar(const EVO_MET_PAR evo_par, const double val){
    switch( evo_par ){
    case POPSIZE:
        n_pob = (int)val;
        break;
    case MAXGEN:
        max_iters = (int)val;
        break;
    case CR:
        seleccion = (int)(val * (double)n_pob);
        seleccion += (seleccion%2);
        DEB_MSG("Selection size: " << seleccion);
        break;
    case MR:
        prob_mutacion = val;
        break;
    }
}



/*  Metodo: setFiltro
    Funcion: Define el metodo de fitrado que se usara sobre la imagen de origen.
*/
void FILTROS::setFiltro( const SEG_FILTRO seg_fil){
    filtro_elegido = seg_fil;
    switch( filtro_elegido ){
        case SS_GABOR:
            pars_optim[3] = false;
            mi_elite->vars[3] = 0.0;
            lim_inf[3] = 0.0;
            min_vars[3] = 0.0;
            fftImgOrigen();
            break;
    }
}


/*  Metodo: setFitness
    Funcion: Define que funcion de evaluacion se usara en el metodo de optimizacion como funcion objetivo.
*/
void FILTROS::setFitness( const FITNESS fit_fun){
    fitness_elegido = fit_fun;
}



/*  Metodo: FILTROS

    Funcion: Inicializa los parametros en valores por defecto.
*/
FILTROS::FILTROS(){
    mi_ruta_log = NULL;
    mi_fplog = NULL;
    resp = NULL;

    rows = 0;
    cols = 0;
    rows_cols = 0;

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

    filtro_elegido = SEG_UNSET;
    metodo_elegido = EVO_UNSET;
    fitness_elegido = FIT_UNSET;

    transformada = false;

    n_pob = 0;
    prob_mutacion = 0.0;
    seleccion = 0;
    max_iters = 0;
}


/*  Metodo: ~FILTROS

    Funcion: Libera la memoria requerida por el objeto.
*/
FILTROS::~FILTROS(){
    delete mi_elite;
    if( semilla ){
        delete semilla;
    }

    if( resp ){
        delete [] resp;
    }

    if( mi_ruta_log ){
        fclose( mi_fplog );
        delete [] mi_ruta_log;
        mi_ruta_log = NULL;
    }

    TIMERS;
    GETTIME_INI;
    if(transformada){
        fftw_free(Img_fft);
        fftw_free(Img_fft_HPF);
    }
    GETTIME_FIN;
    std::cout << COLOR_BACK_GREEN COLOR_BLACK "Tiempo para liberar memoria del FFT: " << DIFTIME << " s." COLOR_NORMAL << std::endl;
}



/*  Metodo: setInputOriginal
    Funcion: Establece las imagenes origen yu ground truth.
*/
void FILTROS::setInput(IMGVTK &img_org){
    org = img_org.base_ptr;
    DEB_MSG("ORG ptr: " << img_org.base_ptr);
    mask = img_org.mask_ptr;
    dest = img_org.segment_ptr;
    ground_truth = img_org.gt_ptr;

    // Obtener las dimensiones de la imagen:
    cols = img_org.cols;
    rows = img_org.rows;
    rows_cols = img_org.rows_cols;

    DEB_MSG("Dimensiones para el filtro: " << cols << "x" << rows);

    resp = new double [rows_cols];

    switch( filtro_elegido ){
        case SS_GABOR:
            pars_optim[3] = false;
            mi_elite->vars[3] = 0.0;
            fftImgOrigen();
            break;
    }
}



/*  Metodo: setPar
    Funcion: Ejecuta el metodo de optimizacion seleccionado para encontrar los mejores parametros para el filtro sobre la imagen de origen.
*/
void FILTROS::setPar(){
    switch( metodo_elegido ){
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
        default: /* EVO_UNSET */
            char mensaje_error[] = COLOR_BACK_BLACK COLOR_RED "<<ERROR: " COLOR_YELLOW "No se ha definido el metodo de optimizacion" COLOR_NORMAL "\n";
            escribirLog( mensaje_error);
    }
}


/*  Metodo: setPar
    Funcion: Establece los parametros del filtro.
*/
void FILTROS::setPar( const PARAMETRO par, const double val ){
    mi_elite->vars[ par ] = val;
    // Si se fija el parametro, se omite este en la busqueda exhaustiva/metaheuristica
    pars_optim[ par ] = false;
    min_vars[ par ] = 0.0;
    lim_inf[ par ] = val;
}


/*  Metodo: getPars
    Funcion: Retorna los parametros utilzados para el filtro.
*/
FILTROS::INDIV FILTROS::getPars(){
    INDIV out_pars;
    memcpy(out_pars.vars, mi_elite->vars, 5*sizeof(double));
    out_pars.eval = mi_elite->eval;
    return out_pars;
}


/*  Metodo: getPars
    Funcion: Retorna los parametros utilzados para el filtro.
*/
int FILTROS::getParametrosOptimizar(){
    int n_pars = 0;
    for( int i = 0; i < 4; i++){
        DEB_MSG("Par[" << i << "]: " << pars_optim[i] );
        if( pars_optim[i] ){
            n_pars++;
        }
    }
    return n_pars;
}


/*  Metodo: setLim
    Funcion: Establece los limites de busqueda de los parametros, si se utiliza para metodos de codificacion binaria, var_delta indica cuantos bits se utilizaran.
*/
void FILTROS::setLim( const PARAMETRO par, const LIMITES lim, const double val){
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




/*  Metodo: filtrar
    Funcion: Aplica el filtro sobre la imagen de origen y la almacena sobre la imagen destino.
*/
void FILTROS::filtrar(){

    memset( resp, 0, rows_cols * sizeof(double));

    double fitness = 0.0;

    switch( fitness_elegido ){
        case ROC:
            fitness = fitnessROC(mi_elite, resp);
            DEB_MSG("Fitnes ROC: " << fitness << " :: L: " << mi_elite->vars[PAR_L] <<  ", T: " << mi_elite->vars[PAR_T] << ", K: " << mi_elite->vars[PAR_K] << ", sigma: " << mi_elite->vars[PAR_SIGMA]);
            break;
        case CORCON:
            fitness = fitnessCorCon(mi_elite, resp);
            DEB_MSG("Fitnes Cor Con: " << fitness << " :: L: " << mi_elite->vars[PAR_L] <<  ", T: " << mi_elite->vars[PAR_T] << ", K: " << mi_elite->vars[PAR_K] << ", sigma: " << mi_elite->vars[PAR_SIGMA]);
            break;
        default: /* FIT_UNSET */
            char mensaje_error[] = COLOR_BACK_BLACK COLOR_RED "<<ERROR: " COLOR_YELLOW "No se ha definido la funcion de evaluacion (fitness)" COLOR_NORMAL "\n";
            escribirLog( mensaje_error);
    }

    memcpy(dest, resp, rows_cols*sizeof(double));
}


/*  Metodo: setLog
    Funcion: Define el stream donde se imprimen los logs del sistema.
*/
void FILTROS::setLog(FILE *fplog){
    mi_fplog = fplog;
}




/*  Metodo: setLog
    Funcion: Define el stream donde se imprimen los logs del sistema.
*/
void FILTROS::setLog( const char *ruta_log){
    mi_ruta_log = new char [(int)strlen(ruta_log) + 1];
    sprintf(mi_ruta_log, "%s", ruta_log);
    if(mi_fplog){
        fclose( mi_fplog );
    }
    mi_fplog = fopen( mi_ruta_log, "w" );
}



/*//        //      //      //      //      //      //        //      //      //      //      //      //        //      //      //      //      //      //        //      //      //      //      //      //
 *
 *              M E T O D O S       D E     F I L T R A D O     D E     I M A G E N E S
 *
 * //        //      //      //      //      //      //        //      //      //      //      //      //        //      //      //      //      //      //        //      //      //      //      //      //
*/
/*  Metodo: respGMF

    Funcion: Obtiene la respuesta del filtro Gaussiano Multiescala, con desviacion estandar: sigma, largo del template: L, ancho del template: T, y numero de rotaciones que se hacen al filtro entre 0 y 180°: K.
*/
void FILTROS::respGMF(INDIV *test, double *resp){
	TIMERS;

	GETTIME_INI;

    const int L = round(test->vars[0]);
    const int T = round(test->vars[1]);
    const int K = round(test->vars[2]);
    const double sigma = test->vars[3];
    const int temp_dims = 1.5 * ((T > L) ? T : L);

    double **templates = new double*[K];

    //// Se calcula el template Gaussiano en rotacion 0°.
    double *gauss_0 = new double[T];
	double *gauss_ptr = gauss_0;

    ////// Se calcula una linea de la gaussiana para el template, luego se copia hacia abajo:
    const double sig_2 = 2.0 * sigma * sigma;
    double sum = 0.0;
    for( double x = -floor((double)T/2.0) + (double)(1 - T%2)/2.0; x <= floor((double)T/2.0) - (double)(1 - T%2)/2.0; x+=1.0){
        *(gauss_ptr) = 1.0 - exp( -(x*x / sig_2) );
        sum += *(gauss_ptr);
        gauss_ptr++;
    }
    sum *= L;

    //// Se resta la media a todo el template, y se divide entre la suma:
    gauss_ptr = gauss_0;
    const double media = sum / ((double)T * (double)L);

    for( double x = -floor((double)T/2.0) + (double)(1 - T%2)/2.0; x <= floor((double)T/2.0) - (double)(1 - T%2)/2.0; x+=1.0){
        *(gauss_ptr) = (*(gauss_ptr) - media) / sum;
         gauss_ptr++;
    }


    //// Se termina de construir el template a 0°:
    templates[0] = new double [temp_dims * temp_dims];
    memset(templates[0], 0, temp_dims*temp_dims*sizeof(double));
    for( int y = 0; y < L; y++ ){
        memcpy( templates[0] + (temp_dims/2 - L/2 + y)*temp_dims + (temp_dims/2 - T/2), gauss_0, T*sizeof(double) );
//        for( int x = 0; x < T; x++){
//             templates[0][(temp_dims/2 - L/2 + y)*temp_dims + (temp_dims/2 - T/2+x)] = 1.0;
//        }
    }


    // Rotar el template segun el numero de rotaciones 'K':
    const double theta_inc = 180.0 / (double)K;
    double theta = 0.0;

    for( int k = 1; k < K; k++){
        theta += theta_inc;
        const double ctheta = cos( -theta * PI/180.0 );
        const double stheta = sin( -theta * PI/180.0 );

        templates[k] = new double [temp_dims * temp_dims];
        memset(templates[k], 0, temp_dims*temp_dims*sizeof(double));
        rotarImg( templates[0], templates[k], ctheta, stheta, temp_dims, temp_dims, temp_dims, temp_dims);
    }

	delete[] gauss_0;
    ////--------------------------------------------------------- Aplicacion del filtro:

    double *resp_tmp = new double [ K * (cols + temp_dims - 1) * (rows + temp_dims - 1) ];
    for( int xy = 0; xy < K * (cols + temp_dims - 1) * (rows + temp_dims - 1); xy++ ){
        *(resp_tmp + xy) = -INF;
    }

	const int offset = (int)(temp_dims/2);
    const int mis_cols = cols;
    const int mis_rens = rows;
    const double *mi_org = org;

    #pragma omp parallel for shared(resp_tmp, mi_org, templates) firstprivate(mis_rens, mis_cols, temp_dims, K)
    for( int k = 0; k < K; k++ ){
        for( int yR = 0; yR < (mis_rens + temp_dims - 1); yR++){
			// Definir los limites en el eje y que pueden recorrerse del template:
			const int min_y = (yR > (temp_dims - 1)) ? (yR - temp_dims + 1) : 0;
            const int max_y = (yR < mis_rens) ? (yR + 1) : mis_rens;
            for( int xR = 0; xR < (mis_cols + temp_dims - 1); xR++){

				// Definir los limites en el eje x que pueden recorrerse del template:
                const int min_x = (xR > (temp_dims - 1)) ? (xR - temp_dims + 1) : 0;
                const int max_x = (xR < mis_cols) ? (xR + 1) : mis_cols;

				double resp_k = 0.0;
                // Convolucionar el template con la vecindad de pixeles de la imagen:
                for( int y = min_y; y < max_y; y++){
                    for( int x = min_x; x < max_x; x++){
                        resp_k += *(templates[k] + (max_y - y - 1)*temp_dims + (max_x - x - 1)) * *(mi_org + y*mis_cols + x);
                    }
                }

                if (resp_k > *(resp_tmp + yR*(mis_cols + temp_dims - 1) + xR)) {
                    *(resp_tmp + yR*(mis_cols + temp_dims - 1) + xR + k*(mis_cols + temp_dims - 1) * (mis_rens + temp_dims - 1)) = resp_k;
                }
            }
        }
    }

    // Liberar la matriz de templates:
    for( int xy = 0; xy < rows_cols; xy++){
        *(resp +xy) =-INF;
    }

    for (int y = 0; y < rows; y++) {
        for (int x = 0; x < cols; x++) {
            for( int k = 0; k < K; k++){
                if( *(resp + x + y*cols) < *(resp_tmp + (x + offset) + (y + offset)*(cols + temp_dims - 1) + k*(mis_cols + temp_dims - 1) * (mis_rens + temp_dims - 1)) ){
                    *(resp + x + y*cols) = *(resp_tmp + (x + offset) + (y + offset)*(cols + temp_dims - 1) + k*(mis_cols + temp_dims - 1) * (mis_rens + temp_dims - 1));
                }
            }
        }
    }

    for( int k = 0; k < K; k++){
        delete [] templates[k];
    }
    delete [] templates;

    delete [] resp_tmp;

	GETTIME_FIN;

	DEB_MSG("Filtrado en " << DIFTIME << " segundos.");
}



/*  Metodo: fftImgOrigen

    Funcion: Obtiene la transformada de Fourier de la imagen original.
*/
void FILTROS::fftImgOrigen(){
    if( rows_cols ){

		TIMERS;

		GETTIME_INI;

        if(!transformada){
            Img_fft = (fftw_complex*) fftw_malloc(rows*(cols/2+1)*sizeof(fftw_complex));
            Img_fft_HPF = (fftw_complex*) fftw_malloc(rows*(cols/2+1)*sizeof(fftw_complex));
            transformada = true;
        }

        double *Img_org = (double*) malloc(rows_cols * sizeof(double));
        for(int xy = 0; xy < rows_cols; xy++){
            *(Img_org + xy) = 1.0 - *(org + xy);
        }

        fftw_plan p_r2c = fftw_plan_dft_r2c_2d(rows, cols, Img_org, Img_fft, FFTW_ESTIMATE);
        fftw_execute(p_r2c);
        fftw_destroy_plan(p_r2c);
        free(Img_org);

		GETTIME_FIN;

		DEB_MSG("Tiempo obtener DFT: " << DIFTIME << " s.");

    }
}



/*  Metodo: respGabor

    Funcion: Obtiene la respuesta del filtro de escala simple de Gabor, largo del template: L, ancho del template: T, y numero de rotaciones que se hacen al filtro entre 0 y 180°: K.
*/
void FILTROS::respGabor(INDIV *test, double *resp){
    const double L = test->vars[0];
    const int T = round(test->vars[1]);
    const int K = round(test->vars[2]);

    // Calculate sx y sy:
    double sx2 = (double)T / (2.0*sqrt(2.0 * log(2.0)));
    sx2 *= sx2;
    const double sy2 = L*L * sx2;

    DEB_MSG("L: " << L << ", T: " << T << ", K: " << K << " :: sx2: " << sx2 << ", sy2: " << sy2);

    const double fx = 1.0 / (double)T;
    const double fy = 0.0;

    double *u = (double*) malloc( cols * sizeof(double));
    double *v = (double*) malloc( rows * sizeof(double));

    for( int x = 0; x < cols; x++){
        *(u + x) = (x - (double)cols/2.0) * (2.0 * PI / (double)cols);
    }
    for( int y = 0; y < rows; y++){
        *(v + y) = (y - (double)rows/2.0) * (2.0 * PI / (double)rows);
    }

    //// Generate the high-pass template HPF
    double *HPF = (double*) malloc( rows*(cols/2 + 1) * sizeof(double));

    DEB_MSG("rows: " << rows << ", cols: " << cols);

    for( int y = 0; y < rows/2; y++ ){
        const double v_y2 = *(v + y) * *(v + y);
        {
            const double u_x2 = *(u) * *(u);
            *(HPF + (y+rows/2)*(cols/2 + 1) + cols/2) = 1.0 - exp( -( sy2*(u_x2 + v_y2) )/2.0 );
        }
        for( int x = cols/2; x < cols; x++){
            const double u_x2 = *(u + x) * *(u + x);
            *(HPF + (y+rows/2)*(cols/2 + 1) + x-cols/2) = 1.0 - exp( -( sy2*(u_x2 + v_y2) )/2.0 );
        }
    }
    for( int y = rows/2; y < rows; y++ ){
        const double v_y2 = *(v + y) * *(v + y);
        {
            const double u_x2 = *(u) * *(u);
            *(HPF + (y-rows/2)*(cols/2 + 1) + cols/2) = 1.0 - exp( -( sy2*(u_x2 + v_y2) )/2.0 );
        }
        for( int x = cols/2; x < cols; x++){
            const double u_x2 = *(u + x) * *(u + x);
            *(HPF + (y-rows/2)*(cols/2 + 1) + x-cols/2) = 1.0 - exp( -( sy2*(u_x2 + v_y2) )/2.0 );
        }
    }

    //// Apply the high-pass filter to the image in the frequencies domain:
    for(int y = 0; y < rows; y++){
		for( int x = 0; x <= cols/2; x++){
            *(*(Img_fft_HPF + x+y*(cols/2+1))  ) = *(*(Img_fft + x+y*(cols/2+1))  ) * *(HPF + x+y*(cols/2 + 1));
            *(*(Img_fft_HPF + x+y*(cols/2+1))+1) = *(*(Img_fft + x+y*(cols/2+1))+1) * *(HPF + x+y*(cols/2 + 1));
        }
    }
    free( HPF );

    // 'Img_filter' is the temporal filtered image in the frequencies domain:
    fftw_complex *Img_filter = (fftw_complex*) fftw_malloc(rows*(cols/2+1)*sizeof(fftw_complex));

    // 'Img_resp' contains the response of the gabor filter at certain orientation 'theta':
    double *Img_resp = (double*) calloc(rows_cols, sizeof(double));

    // 'max_resp' saves the highest response for every pixel of the input:
    double *max_resp = (double*) malloc(rows_cols * sizeof(double));
    for( int xy = 0; xy < rows_cols; xy++){
        *(max_resp + xy) = -INF;
    }

    fftw_plan p_c2r;
    const double theta_increment = 180.0 / (double)K;

    double Gabor_xy, Vr, Ur;

    double *max_resp_angles = (double*) malloc( rows_cols * sizeof(double) );

    for( double theta = 0.0; theta < 180.0; theta+=theta_increment){
        const double stheta = sin(theta*PI/180.0);
        const double ctheta = cos(theta*PI/180.0);

        //// The Gabor filter is calculated for the rotated base at 'theta' degrees:
        for( int y = 0; y < rows/2; y++){
            const double v_y = *(v+y+rows/2);
            for( int x = 0; x < cols/2; x++){
                Ur = *(u + x+cols/2)*ctheta + v_y*stheta;
                Vr =-*(u + x+cols/2)*stheta + v_y*ctheta;
                Gabor_xy = exp(-(0.5)*(sx2*(Ur*Ur + 4.0*PI*fx*PI*fx) + sy2*(Vr*Vr + 4.0*PI*fy*PI*fy))) * cosh(2.0*PI*(sx2*fx*Ur + sy2*fy*Vr));
                *(*(Img_filter + y*(cols/2+1) + x)  ) = *(*(Img_fft_HPF + y*(cols/2+1) + x)  ) * Gabor_xy;
                *(*(Img_filter + y*(cols/2+1) + x)+1) = *(*(Img_fft_HPF + y*(cols/2+1) + x)+1) * Gabor_xy;
            }
            Ur = *(u)*ctheta + v_y*stheta;
            Vr =-*(u)*stheta + v_y*ctheta;
            Gabor_xy = exp(-(0.5)*(sx2*(Ur*Ur + 4.0*PI*fx*PI*fx) + sy2*(Vr*Vr + 4.0*PI*fy*PI*fy))) * cosh(2.0*PI*(sx2*fx*Ur + sy2*fy*Vr));
            *(*(Img_filter + y*(cols/2+1) + cols/2)  ) = *(*(Img_fft_HPF + y*(cols/2+1) + cols/2)  ) * Gabor_xy;
            *(*(Img_filter + y*(cols/2+1) + cols/2)+1) = *(*(Img_fft_HPF + y*(cols/2+1) + cols/2)+1) * Gabor_xy;
        }

        for( int y = rows/2; y < rows; y++){
            const double v_y = *(v + y-rows/2);
            for( int x = 0; x < cols/2; x++){
                Ur = *(u + x+cols/2)*ctheta + v_y*stheta;
                Vr =-*(u + x+cols/2)*stheta + v_y*ctheta;
                Gabor_xy = exp(-(0.5)*(sx2*(Ur*Ur + 4.0*PI*fx*PI*fx) + sy2*(Vr*Vr + 4.0*PI*fy*PI*fy))) * cosh(2.0*PI*(sx2*fx*Ur + sy2*fy*Vr));
                *(*(Img_filter + y*(cols/2+1) + x)  ) = *(*(Img_fft_HPF + y*(cols/2+1) + x)  ) * Gabor_xy;
                *(*(Img_filter + y*(cols/2+1) + x)+1) = *(*(Img_fft_HPF + y*(cols/2+1) + x)+1) * Gabor_xy;

            }
            Ur = *(u)*ctheta + v_y*stheta;
            Vr =-*(u)*stheta + v_y*ctheta;
            Gabor_xy = exp(-(0.5)*(sx2*(Ur*Ur + 4.0*PI*fx*PI*fx) + sy2*(Vr*Vr + 4.0*PI*fy*PI*fy))) * cosh(2.0*PI*(sx2*fx*Ur + sy2*fy*Vr));
            *(*(Img_filter + y*(cols/2+1) + cols/2)  ) = *(*(Img_fft_HPF + y*(cols/2+1) + cols/2)  ) * Gabor_xy;
            *(*(Img_filter + y*(cols/2+1) + cols/2)+1) = *(*(Img_fft_HPF + y*(cols/2+1) + cols/2)+1) * Gabor_xy;
        }

        //// Transform the response to the original domain:
        p_c2r = fftw_plan_dft_c2r_2d(rows, cols, Img_filter, Img_resp, FFTW_ESTIMATE);
        fftw_execute(p_c2r);
        fftw_destroy_plan(p_c2r);

        //// Update the highest response for every pixel:
        for( int xy = 0; xy < rows_cols; xy++){
            if( *(max_resp + xy) < *(Img_resp + xy)){
                *(max_resp + xy) = *(Img_resp + xy);
                *(max_resp_angles + xy) = -theta * PI / 180.0 + PI/2.0;
            }
        }
    }

    fftw_free(Img_filter);
    free(Img_resp);
    free(u);
    free(v);

    for( int xy = 0; xy < rows_cols; xy++){
        *(resp + xy) = (*(mask + xy) > 0.5) ? (*(max_resp + xy) / rows_cols) : 0.0;
    }

    free(max_resp);
    free(max_resp_angles);

}





/*  Metodo: comp_resp

    Funcion: Compara la respuesta de dos pixeles.
*/
int comp_resp( const void *a, const void *b){
    double resp_a = *(double*)a;
    double resp_b = *(double*)b;
    if( resp_a > resp_b ){
        return 1;
    }else if(resp_a < resp_b){
        return -1;
    }
    return 0;
}





/*  Metodo: calcROC

    Funcion: Calculates the area under the Receiver Operating Characteristic curve.
*/
double FILTROS::calcROC( double *resp ){

    if( !ground_truth ){
        char mensaje_error[] = COLOR_BACK_BLACK COLOR_RED "<<ERROR: " COLOR_YELLOW "No se cargo la imagen ground-truth" COLOR_NORMAL "\n";
        escribirLog( mensaje_error);
        return 0.0;
    }

    TIMERS;

    GETTIME_INI;

    const double delta = 1e-4;

    //// Count the number of positive and negative individuals on the ground-truth that belong to the masked area:
    int n_positive = 0, n_negative = 0;
    for( int xy = 0; xy < rows_cols; xy++){
        if( *(ground_truth + xy) > 0.5 && *(mask + xy) > 0.5){
            n_positive++;
        }else{
            n_negative++;
        }
    }

    double *positive_array = (double*) malloc( n_positive * sizeof(double));
    double *negative_array = (double*) malloc( n_negative * sizeof(double));
    n_positive = 0;
    n_negative = 0;

    //// Store the values of the response in the corresponding array according to the ground-truth:
    double max_resp =-INF;
    double min_resp = INF;

    for( int xy = 0; xy < rows_cols; xy++){
        //// Store the value only if belongs to the masked area:
        if( *(mask + xy) > 0.5 ){
            if( *(ground_truth + xy) >= 0.5 ){
                positive_array[n_positive] = resp[xy];
                n_positive++;
            }else{
                negative_array[n_negative] = resp[xy];
                n_negative++;
            }

            //// Determine minima and maxima in the response:
            if(max_resp < resp[xy]){
                max_resp = resp[xy];
            }
            if(min_resp > resp[xy]){
                min_resp = resp[xy];
            }
        }
    }

    DEB_MSG("min resp:" << COLOR_BLINK << min_resp << COLOR_NORMAL << ", max resp: " << COLOR_BLINK << max_resp << COLOR_NORMAL);

    //// Sort both arrays:
    qsort(positive_array, n_positive, sizeof(double), comp_resp);
    qsort(negative_array, n_negative, sizeof(double), comp_resp);

    //// Start with the most discriminant threshold
    double threshold = max_resp;
    int i_positive = (n_positive-1), i_negative = (n_negative-1), i_threshold = 0;
    double TP = 0.0, TN = (double)n_negative;

    double TPF_old = 0.0;
    double FPF_old = 0.0;

    double TPF_new, FPF_new;
    double Az = 0.0;

#ifndef NDEBUG
    FILE *fp_ROC = fopen("ROC_curve.fcs", "w");
#endif

    min_resp -= delta;
    while( threshold > min_resp ){
        //// Count the positive and negative individuals at this threshold:

        while( i_positive > 0 ){
            if(positive_array[i_positive] >= threshold){
                TP += 1.0;
                i_positive--;
            }else{
                break;
            }
        }

        while( i_negative > 0){
            if( negative_array[i_negative] > threshold){
                TN -= 1.0;
                i_negative--;
            }else{
                break;
            }
        }

        //// Compute the True Positive Fraction and the False Positive Fraction:
        TPF_new = TP/(double)n_positive;
        FPF_new = 1.0 - TN/(double)n_negative;

        //// Approximate the area under the ROC curve with the Trapezoid Rule:
        Az += (FPF_new - FPF_old)*(TPF_new + TPF_old)/2.0;

        //// Update the threshold:
        threshold -= delta;
        i_threshold++;

        // Mover los valores de la curva:
        TPF_old = TPF_new;
        FPF_old = FPF_new;

#ifndef NDEBUG
        fprintf( fp_ROC , "%i %f %f %f %f\n", i_threshold, threshold, TPF_old, FPF_old, Az);
#endif
    }

#ifndef NDEBUG
    fclose(fp_ROC);
#endif
    free(positive_array);
    free(negative_array);

    GETTIME_FIN;

    DEB_MSG("Az: " << Az << " en " << DIFTIME << " s.");

    return Az;
}



/*  Metodo: fitnessROC
    Funcion: Evalua la eficacia de un detector utilizando el area bajo la curva de ROC:
*/
double FILTROS::fitnessROC( INDIV *test, double *mi_resp ){
    //// Generar la respuesta del filtro de establecido para los parametros dados:
    switch(filtro_elegido){
        case SS_GABOR:
            respGabor(test, mi_resp);
            break;
        case GMF:
            respGMF(test, mi_resp);
            break;
    }

    return calcROC(mi_resp);
}








/*  Metodo: calcCorCon

    Funcion: Calculates the correlation and contrast of the filtered response.
*/
double FILTROS::calcCorCon(double *resp){

    //// Scale the response to 8 levels:
    const int levels = 8;

    double max_resp =-INF;
    double min_resp = INF;

    for( int xy = 0; xy< rows_cols; xy++){
        if( max_resp < *(resp + xy)){
            max_resp = *(resp + xy);
        }
        if( min_resp > *(resp + xy)){
            min_resp = *(resp + xy);
        }
    }

    const double range = max_resp - min_resp;
    const double scalar = (double)(levels - 1) / range;

    int *scaled_resp = (int*) malloc(rows_cols * sizeof(int));

    for( int xy = 0; xy < rows_cols; xy++){
        *(scaled_resp + xy) = (int)round( scalar * (*(resp + xy) - min_resp) );
    }

    //// Compute the GLCM matrix using an intensity at 8 levels, for 4 orientations:
    double *GLCM_mat = (double*) calloc((levels * levels * 4), sizeof(double));

    const double fraction_1 = 1.0 / (double)((cols-1) *   rows  );
    const double fraction_2 = 1.0 / (double)((cols-1) * (rows-1));
    const double fraction_3 = 1.0 / (double)(  cols   * (rows-1));
    const double fraction_4 = 1.0 / (double)((cols-1) * (rows-1));

    double contrast_1 = 0.0, contrast_2 = 0.0, contrast_3 = 0.0, contrast_4 = 0.0;

    double m_i[4] = {0.0, 0.0, 0.0, 0.0};
    double m_j[4] = {0.0, 0.0, 0.0, 0.0};

    for( int x = 0; x < cols-1; x++){
        // SI(x + 1, y    )
        const int i = *(scaled_resp +   x  );
        const int j = *(scaled_resp + (x+1));

        *(GLCM_mat + (levels*4)*i + 4*j  ) = *(GLCM_mat + (levels*4)*i + 4*j  ) + fraction_1;
        contrast_1 += (double)((i-j)*(i-j)) * fraction_1;

        *(m_i  ) += (double)(i+1) * fraction_1;
        *(m_j  ) += (double)(j+1) * fraction_1;
    }

    for( int y = 1; y < rows; y++){
        // SI(x   , y - 1)
        const int i = *(scaled_resp + (y+1)*cols - 1);
        const int j = *(scaled_resp +   y  *cols - 1);
        *(GLCM_mat + (levels*4)*i + 4*j+2) = *(GLCM_mat + (levels*4)*i + 4*j+2) + fraction_3;
        contrast_3 += (double)((i-j)*(i-j)) * fraction_3;

        *(m_i+2) += (double)(i+1) * fraction_3;
        *(m_j+2) += (double)(j+1) * fraction_3;
    }


    for( int y = 1; y < rows; y++){
        for( int x = 0; x < cols-1; x++){
            {// SI(x + 1, y    )
                const int i = *(scaled_resp + ( x ) + ( y )*cols);
                const int j = *(scaled_resp + (x+1) + ( y )*cols);
                *(GLCM_mat + (levels*4)*i + 4*j  ) = *(GLCM_mat + (levels*4)*i + 4*j  ) + fraction_1;
                contrast_1 += (double)((i-j)*(i-j)) * fraction_1;

                *(m_i  ) += (double)(i+1) * fraction_1;
                *(m_j  ) += (double)(j+1) * fraction_1;

            }
            {// SI(x + 1, y - 1)
                const int i = *(scaled_resp + ( x ) + ( y )*cols);
                const int j = *(scaled_resp + (x+1) + (y-1)*cols);
                *(GLCM_mat + (levels*4)*i + 4*j+1) = *(GLCM_mat + (levels*4)*i + 4*j+1) + fraction_2;
                contrast_2 += (double)((i-j)*(i-j)) * fraction_2;

                *(m_i+1) += (double)(i+1) * fraction_2;
                *(m_j+1) += (double)(j+1) * fraction_2;
            }
            {// SI(x    , y - 1)
                const int i = *(scaled_resp + ( x ) + ( y )*cols);
                const int j = *(scaled_resp + ( x ) + (y-1)*cols);
                *(GLCM_mat + (levels*4)*i + 4*j+2) = *(GLCM_mat + (levels*4)*i + 4*j+2) + fraction_3;
                contrast_3 += (double)((i-j)*(i-j)) * fraction_3;

                *(m_i+2) += (double)(i+1) * fraction_3;
                *(m_j+2) += (double)(j+1) * fraction_3;
            }
            {// SI(x - 1, y - 1)
                const int i = *(scaled_resp + (x+1) + ( y )*cols);
                const int j = *(scaled_resp + ( x ) + (y-1)*cols);
                *(GLCM_mat + (levels*4)*i + 4*j+3) = *(GLCM_mat + (levels*4)*i + 4*j+3) + fraction_4;
                contrast_4 += (double)((i-j)*(i-j)) * fraction_4;

                *(m_i+3) += (double)(i+1) * fraction_4;
                *(m_j+3) += (double)(j+1) * fraction_4;
            }
        }
    }

    free(scaled_resp);
    DEB_MSG( "Contrast: " << contrast_1 << ", " << contrast_2 << ", " << contrast_3 << ", " << contrast_4 << ", mean: " << (contrast_1 + contrast_2 + contrast_3 + contrast_4)/ 4.0);

    double s_i[4] = {0.0, 0.0, 0.0, 0.0};
    double s_j[4] = {0.0, 0.0, 0.0, 0.0};

    //// Calculate variance for i and j form GLCM:
    for( int i = 0; i < levels; i++){
        const double i_m_1 = (double)(i+1) - *(m_i  );
        const double i_m_2 = (double)(i+1) - *(m_i+1);
        const double i_m_3 = (double)(i+1) - *(m_i+2);
        const double i_m_4 = (double)(i+1) - *(m_i+3);

        for( int j = 0; j < levels; j++){
            const double j_m_1 = (double)(j+1) - *(m_j  );
            const double j_m_2 = (double)(j+1) - *(m_j+1);
            const double j_m_3 = (double)(j+1) - *(m_j+2);
            const double j_m_4 = (double)(j+1) - *(m_j+3);

            *(s_i  ) += i_m_1 * i_m_1 * *(GLCM_mat + (levels*4)*i + 4*j  );
            *(s_j  ) += j_m_1 * j_m_1 * *(GLCM_mat + (levels*4)*i + 4*j  );

            *(s_i+1) += i_m_2 * i_m_2 * *(GLCM_mat + (levels*4)*i + 4*j+1);
            *(s_j+1) += j_m_2 * j_m_2 * *(GLCM_mat + (levels*4)*i + 4*j+1);

            *(s_i+2) += i_m_3 * i_m_3 * *(GLCM_mat + (levels*4)*i + 4*j+2);
            *(s_j+2) += j_m_3 * j_m_3 * *(GLCM_mat + (levels*4)*i + 4*j+2);

            *(s_i+3) += i_m_4 * i_m_4 * *(GLCM_mat + (levels*4)*i + 4*j+3);
            *(s_j+3) += j_m_4 * j_m_4 * *(GLCM_mat + (levels*4)*i + 4*j+3);
        }
    }

    //// Calculate standard deviation for i and j form GLCM:
    *(s_i  ) = sqrt( *(s_i  ) );
    *(s_j  ) = sqrt( *(s_j  ) );

    *(s_i+1) = sqrt( *(s_i+1) );
    *(s_j+1) = sqrt( *(s_j+1) );

    *(s_i+2) = sqrt( *(s_i+2) );
    *(s_j+2) = sqrt( *(s_j+2) );

    *(s_i+3) = sqrt( *(s_i+3) );
    *(s_j+3) = sqrt( *(s_j+3) );

    //// Calculate correlation form GLCM:
    double correlation_1 = 0.0, correlation_2 = 0.0, correlation_3 = 0.0, correlation_4 = 0.0;
    for( int i = 0; i < levels; i++){
        const double i_m_1 = (double)(i+1) - *(m_i  );
        const double i_m_2 = (double)(i+1) - *(m_i+1);
        const double i_m_3 = (double)(i+1) - *(m_i+2);
        const double i_m_4 = (double)(i+1) - *(m_i+3);

        const double i_s_1 = *(s_i  );
        const double i_s_2 = *(s_i+1);
        const double i_s_3 = *(s_i+2);
        const double i_s_4 = *(s_i+3);

        for( int j = 0; j < levels; j++){
            const double j_m_1 = (double)(j+1) - *(m_j  );
            const double j_m_2 = (double)(j+1) - *(m_j+1);
            const double j_m_3 = (double)(j+1) - *(m_j+2);
            const double j_m_4 = (double)(j+1) - *(m_j+3);

            correlation_1 += *(GLCM_mat + (levels*4)*i + 4*j  ) * i_m_1*j_m_1 / (i_s_1 * *(s_j  ));
            correlation_2 += *(GLCM_mat + (levels*4)*i + 4*j+1) * i_m_2*j_m_2 / (i_s_2 * *(s_j+1));
            correlation_3 += *(GLCM_mat + (levels*4)*i + 4*j+2) * i_m_3*j_m_3 / (i_s_3 * *(s_j+2));
            correlation_4 += *(GLCM_mat + (levels*4)*i + 4*j+3) * i_m_4*j_m_4 / (i_s_4 * *(s_j+3));
        }
    }
    free(GLCM_mat);

    DEB_MSG( "Correlation: " << correlation_1 << ", " << correlation_2 << ", " << correlation_3 << ", " << correlation_4 << ", mean: " << (correlation_1 + correlation_2 + correlation_3 + correlation_4)/ 4.0);

    return (correlation_1 + correlation_2 + correlation_3 + correlation_4) / 4.0 + (double)((levels-1)*(levels-1)) - (contrast_1 + contrast_2 + contrast_3 + contrast_4)/4.0;
}



/*  Metodo: fitnessCorCon
    Funcion: Evalua los parametros 'L', 'T' y //'K'// para el filtro de Gabor y el la correlacion y contraste:
*/
double FILTROS::fitnessCorCon( INDIV *test, double *resp){
    //// Generar la respuesta del filtro de establecido para los parametros dados:
    switch(filtro_elegido){
        case SS_GABOR:
            respGabor(test, resp);
            break;
        case GMF:
            respGMF(test, resp);
            break;
    }

    return calcCorCon(resp);
}





/*//        //      //      //      //      //      //        //      //      //      //      //      //        //      //      //      //      //      //        //      //      //      //      //      //
 *
 *              G E N E R A D O R E S       D E     N U M E R O S       A L E A T O R I O S
 *
 * //        //      //      //      //      //      //        //      //      //      //      //      //        //      //      //      //      //      //        //      //      //      //      //      //
*/




/*	Metodo:        ini_semilla
    Funcion:: Inicializa una semilla para el generador HybTaus. Si la semilla dada por el usuario es 0, se genera una semilla aleatoriamente y se utliza, de otro modo se usa la semilla dada por el usuario.
*/
FILTROS::STAUS* FILTROS::ini_semilla(unsigned int semilla_i){
    if(!semilla_i){
        srand(clock());
        semilla_i = rand();
    }

    STAUS *mi_semilla = new STAUS;
    mi_semilla->z1 = lcg_r(&semilla_i);
    mi_semilla->z2 = lcg_r(&semilla_i);
    mi_semilla->z3 = lcg_r(&semilla_i);

    return mi_semilla;
}




/*	Metodo:         lcg_s
    Funcion: Genera un numero pseudo-aleatorio por medio del metodo Linear Congruential Generator. El generador utiliza una semilla static y los parametros que utiliza son los mismos que utiliza gcc.
*/
unsigned int FILTROS::lcg_s(){
    static unsigned int st_semilla = 0;
    if(!st_semilla){
        st_semilla = rand();
    }
    const unsigned int m = pow(2,31);
    return (1103515245u*st_semilla + 12345) % m;
}




/* Metodo:        lcg_r
    Funcion: Genera un numero pseudo-aleatorio por medio del metodo Linear Congruential Generator. El generador utiliza una semilla enviada por referencia, la cual actualiza y los parametros que utiliza son los mismos que utiliza gcc.
*/
unsigned int FILTROS::lcg_r(unsigned int *mi_semilla){
    const unsigned int m = pow(2,31);
    return *mi_semilla= (1103515245u* *mi_semilla + 12345) % m;
}




/*	Metodo:        tausStep
    Funcion: Genera un numero pseudo-aleatorio por medio del metodo Taus Step.
*/
unsigned int FILTROS::tausStep(unsigned int *z, const int S1, const int S2, const int S3, const unsigned int M){
    const unsigned int b = (((*z << S1) ^ *z) >> S2);
    return *z = (((*z & M) << S3) ^ b);
}




/*	Metodo:        HybTaus
    Funcion:: Genera un numero pseudo-aleatorio por medio del metodo Hybrid Taus Step entre par1 y par2.
*/
double FILTROS::HybTaus(const double par1, const double par2){
// Combined period is lcm(p1,p2,p3,p4)~ 2^121
    double num = 2.3283064365387e-10 * (
        // Periods
        tausStep(&semilla->z1, 13, 19, 12, 4294967294UL) ^ // p1=2^31-1
        tausStep(&semilla->z2, 2, 25, 4, 4294967288UL) ^ // p2=2^30-1
        tausStep(&semilla->z3, 3, 11, 17, 4294967280UL) ^ // p3=2^28-1
        lcg_s()	// p4=2^32
    );

    return (par2 - par1) * num + par1;
}




/*	Metodo         anorm_est
    Funcion: Genera un número aleatorio de la distribución normal estándar, por el método de transformación Box-Muller.
*/
double FILTROS::anorm_est(){
    double x1, x2, w, y1;//, y2;
    do{
        x1 = HybTaus(-1.0, 1.0);
        x2 = HybTaus(-1.0, 1.0);
        w = x1*x1 + x2*x2;
    }while( w >= 1.0);
    w = sqrt( (-2.0 * log(w)) / w);
    y1 = x1 * w;
//	y2 = x2 * w;
    return y1;
}




/*//        //      //      //      //      //      //        //      //      //      //      //      //        //      //      //      //      //      //        //      //      //      //      //      //
 *
 *              A L G O R I T M O S     E V O L U T I V O S
 *
 * //        //      //      //      //      //      //        //      //      //      //      //      //        //      //      //      //      //      //        //      //      //      //      //      //
*/

/*
    Metodo:        compIndiv
    Funcion:    Compara el indidividuo A contra el B.
        Si la evaluacion del individuo A es menor a la del B, se retorna +,
                             si es igual, se retorna 0,
                             si es mayor, se retorna -.
*/
int FILTROS::compIndiv(const void* A, const void* B){
    INDIV *indA = (INDIV*)A;
    INDIV *indB = (INDIV*)B;
    if(indA->eval < indB->eval) return  1;
    if(indA->eval == indB->eval) return  0;
    if(indA->eval > indB->eval) return -1;
}




// BUMDA ----- BUMDA ----- BUMDA ----- BUMDA ----- BUMDA ----- BUMDA ----- BUMDA ----- BUMDA ----- BUMDA ----- BUMDA ----- BUMDA ----- BUMDA ----- BUMDA ----- BUMDA ----- BUMDA ----- BUMDA ----- BUMDA ----- BUMDA ----- BUMDA ----- BUMDA ----- BUMDA ----- BUMDA ----- BUMDA ----- BUMDA ----- BUMDA ----- BUMDA ----- BUMDA ----- BUMDA ----- BUMDA ----- BUMDA ----- BUMDA
/*
    Metodo:        generarPobInicial
    Funcion:    Genera una poblacion inicial dentro de los limites establecidos.

*/
void FILTROS::generarPobInicial(INDIV *poblacion){
    // VERIFICAR CUALES PARAMETROS SE BUSCAN:
    using namespace std;
    for(int i = 0; i < n_pob; i++){
        memcpy( poblacion[i].vars, mi_elite->vars, 5*sizeof(double) );

        for( int j = 0; j < n_pars; j++){
            const unsigned int k = idx_pars[j];
            (poblacion + i)->vars[k] = HybTaus(lim_inf[k], lim_sup[k]);
        }

        switch( fitness_elegido ){
            case ROC:
                (poblacion + i)->eval = fitnessROC( poblacion + i, resp );
                break;
            case CORCON:
                (poblacion + i)->eval = fitnessCorCon( poblacion + i, resp );
                break;
        }
    }
}\





/*
    Metodo:        generarPob
    Funcion:    Genera una poblacion de tamaño 'n_pob' individuos con una distribucion normal con medias y varianzas apra cada parametro los indicados.

*/
void FILTROS::generarPob(double medias[4], double varianzas[4], INDIV *poblacion){
    double val_gen;

    //Se generan n_pob - 1 individuos, y permanece el elite que se tiene.
    // VERIFICAR CUALES PARAMETROS SE BUSCAN:
    using namespace std;
    for(register int i = 0; i < n_pob; i++){
        for( int j = 0; j < n_pars; j++){
            const unsigned int k = idx_pars[j];
            val_gen =  anorm( medias[k], varianzas[k]);
            val_gen = (val_gen <= lim_sup[ k ]) ? ((val_gen >= lim_inf[ k ]) ? val_gen : lim_inf[ k ]) : lim_sup[ k ];
            (poblacion + i)->vars[ k ] = val_gen;
        }
        switch( fitness_elegido ){
            case ROC:
                (poblacion + i)->eval = fitnessROC( poblacion + i, resp );
                break;

            case CORCON:
                (poblacion + i)->eval = fitnessCorCon( poblacion + i, resp );
                break;
        }
    }
}



/*
    Metodo:        calcularPars
    Funcion:    Calula la media y la varianza para cada variable de la poblacion, utilizando solo a los individuos de la seleccion.
*/
void FILTROS::calcularPars(const INDIV *poblacion, const int truncamiento, double *medias, double *varianzas){

    double sum_evals = 0.0;
    const double gx_sel = (poblacion + truncamiento)->eval - 1e-6;
        //Se calculan las medias primero:
        memset(medias, 0, 4*sizeof(double));

        for( int i = 0; i < truncamiento; i++){//Este es el ciclo para toda la seleccion.
            const double g_bar = (poblacion + i)->eval - gx_sel;
            sum_evals += g_bar;
            for( int j = 0; j < n_pars; j++){
                const unsigned int k = idx_pars[j];
                *(medias + k) = *(medias + k) + (poblacion + i)->vars[ k ] * g_bar;
            }
        }

        for( int j = 0; j < n_pars; j++){
            *( medias + idx_pars[j] ) = *( medias + idx_pars[j] ) / sum_evals;
        }
        //Se calculan las varianzas:
        double tmp;
        memset(varianzas, 0, 4*sizeof(double));
        for( int i = 0; i < truncamiento; i++){
            const double g_bar = (poblacion + i)->eval - gx_sel;
            for( int j = 0; j < n_pars; j++){

                const unsigned int k = idx_pars[j];
                tmp = (poblacion+i)->vars[k] - *(medias + k);
                *(varianzas + k) = *(varianzas + k) + tmp*tmp * g_bar;
            }
        }

        for( int j = 0; j < n_pars; j++){
            *(varianzas + idx_pars[j] ) = *(varianzas + idx_pars[j] ) / (sum_evals);
        }
}



/*
    Metodo:     seleccionarPob
    Funcion:    Realiza la seleccion por truncamiento, redefiniendo el tetha para esta iteracion dentro del mismo apuntador de tetha.
*/
int FILTROS::seleccionarPob(double *theta_t, const INDIV *poblacion){
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
void FILTROS::BUMDA(){
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
    semilla = ini_semilla(0);

    // Inicia el algoritmo BUMDA
    int k = 1, truncamiento = n_pob-1;
    bool procesar = true;

    //// Generar la primer poblacion:
    TIMERS;
    GETTIME_INI;
    generarPobInicial(poblacion);
    qsort((void*)poblacion, n_pob, sizeof(INDIV), compIndiv); // El mejor fitness queda en la posicion 0

    //Se define el theta en el tiempo 0 como el minimo de la poblacion inicial.
    double tetha_t = (poblacion + n_pob - 1)->eval;
    char mensaje_iter[] = "Best fit: X.XXXXXXXX, L: XX.XXX , T: XX.XXX, Sigma: XX.XXX, K: XXX.X :: Laste eval: X.XXXX, L: XX.XXX , T: XX.XXX, Sigma: XX.XXX, K: XXX.X :: en XXXXX.XXXX s.\n";
    escribirLog( "iteration\tbest_fit\tT\tL\tSigma\tK\telapsed_time\n" );

    do{
        // Guardar el elite como el mejor de los individuos en la posicion 'n_pob' del arreglo:
        memcpy(poblacion + n_pob, poblacion, sizeof(INDIV));

        GETTIME_FIN;
        sprintf( mensaje_iter, "%i\t%1.8f\t%1.3f\t%1.3f\t%2.3f\t%3.0f\t%5.4f\n", k, poblacion->eval, poblacion->vars[PAR_T], poblacion->vars[PAR_L], poblacion->vars[PAR_SIGMA], poblacion->vars[PAR_K], DIFTIME);
        escribirLog( mensaje_iter );
        barraProgreso( k, max_iters);

        calcularPars(poblacion, truncamiento, medias, varianzas);
        generarPob(medias, varianzas, poblacion);
        qsort((void*)poblacion, n_pob+1, sizeof(INDIV), compIndiv); // El mejor fitness queda en la posicion 0 ya comparando el elite

        //En truncamiento se guarda el indice del individuo a partir del cual se trunca la poblacion.
        truncamiento = seleccionarPob(&tetha_t, poblacion);
        k++;

        //// Verificar condiciones de paro:
        int condicion = 0;
        for(int j = 0; j < n_pars; j++){
            const unsigned int v = idx_pars[j];
            DEB_MSG("var[" << v << "] = " << *(varianzas + v) << "/" << *(min_vars + v));
            condicion += ( *(varianzas + v) <= *(min_vars + v));
        }
        procesar = (k < max_iters) && (condicion < n_pars);

    }while(procesar);

    memcpy(mi_elite, poblacion, sizeof(INDIV));

    GETTIME_FIN;
    sprintf( mensaje_iter, "%i\t%1.8f\t%1.3f\t%1.3f\t%2.3f\t%3.0f\t%5.4f\n", k, mi_elite->eval, mi_elite->vars[PAR_T], mi_elite->vars[PAR_L], mi_elite->vars[PAR_SIGMA], mi_elite->vars[PAR_K], DIFTIME);
    escribirLog( mensaje_iter );
    barraProgreso( max_iters, max_iters);

    delete [] poblacion;
}




// UMDA ----- UMDA ----- UMDA ----- UMDA ----- UMDA ----- UMDA ----- UMDA ----- UMDA ----- UMDA ----- UMDA ----- UMDA ----- UMDA ----- UMDA ----- UMDA ----- UMDA ----- UMDA ----- UMDA ----- UMDA ----- UMDA ----- UMDA ----- UMDA ----- UMDA ----- UMDA ----- UMDA ----- UMDA ----- UMDA ----- UMDA ----- UMDA ----- UMDA ----- UMDA ----- UMDA ----- UMDA ----- UMDA ----- UMDA

/*
    Metodo:     generar_pob
    Funcion:   Genera una nueva poblacion en bsae al muestreo de cada atributo.
*/
void FILTROS::generarPob(INDIV *poblacion, const double *probs, const double *deltas_var){
    for(int i = 0; i < n_pob; i++){
        // Se muestrean los genes para cada parametro:
        unsigned int bits_recorridos = 0;
        for( int j = 0; j < n_pars; j++){
            const unsigned int k = idx_pars[j];
            double cadena_val = 0;
            double pow_2 = 1.0;
            for( int b = 0; b < (int)min_vars[k]; b++, bits_recorridos++, pow_2*=2.0){
                (poblacion + i)->cadena[ bits_recorridos ] = (HybTaus(0.0, 1.0) <= *(probs + bits_recorridos)) ? 1 : 0;
                cadena_val += ((poblacion + i)->cadena[ bits_recorridos ] ? 1.0 : 0.0) * (double)pow_2;
            }
            // Asignar el valor segun la cadena formada para el atributo 'k':
            (poblacion + i)->vars[k] = *(deltas_var + k) * cadena_val + lim_inf[k];
			DEB_MSG(COLOR_GREEN "[" COLOR_BLUE << k << COLOR_GREEN "] cadena: " COLOR_BACK_WHITE COLOR_BLACK << cadena_val << "/" << (poblacion + i)->vars[k] << " :: " COLOR_BACK_BLACK COLOR_BLUE << "delta: " << *(deltas_var + k) << " :: lim_inf: " << lim_inf[k] << " :: lim_sup: " << lim_sup[k] << COLOR_NORMAL);
        }

        switch( fitness_elegido ){
            case ROC:
                (poblacion + i)->eval = fitnessROC( poblacion + i, resp );
                break;
            case CORCON:
                (poblacion + i)->eval = fitnessCorCon( poblacion + i, resp );
                break;
        }
    }
}




/*
    Metodo:     calcular_pars
    Funcion:    Calcula las probabilidades de cada distribucion marginal.
*/
void FILTROS::calcularPars(const INDIV *poblacion, const int n_bits, const int truncamiento, double *probs){
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
void FILTROS::UMDA(){

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
    for( int b = 0; b < n_bits; b++){
        probs[b] = 0.5;
    }

    if(semilla){
        delete semilla;
    }
    semilla = ini_semilla(0);

    // Inicia el algoritmo UMDA
    int k = 0, truncamiento = (int)( (double)n_pob * 0.6);
    bool procesar = true;


    char mensaje_iter[] = "Best fit: X.XXXXXXXX, L: XX.XXX , T: XX.XXX, Sigma: XX.XXX, K: XXX.X :: Laste eval: X.XXXX, L: XX.XXX , T: XX.XXX, Sigma: XX.XXX, K: XXX.X :: en XXXXX.XXXX s.\n";
    escribirLog( "iteration\tbest_fit\tT\tL\tSigma\tK\telapsed_time\n" );

    TIMERS;
    GETTIME_INI;

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

        GETTIME_FIN;
        sprintf( mensaje_iter, "%i\t%1.8f\t%1.3f\t%1.3f\t%2.3f\t%3.0f\t%5.4f\n", k, mi_elite->eval, mi_elite->vars[PAR_T], mi_elite->vars[PAR_L], mi_elite->vars[PAR_SIGMA], mi_elite->vars[PAR_K], DIFTIME);
        escribirLog( mensaje_iter );
        barraProgreso( k, max_iters);
    }while(procesar);

    GETTIME_FIN;
    sprintf( mensaje_iter, "%i\t%1.8f\t%1.3f\t%1.3f\t%2.3f\t%3.0f\t%5.4f\n", k, mi_elite->eval, mi_elite->vars[PAR_T], mi_elite->vars[PAR_L], mi_elite->vars[PAR_SIGMA], mi_elite->vars[PAR_K], DIFTIME);
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
double FILTROS::generarPobInicial(INDIV *poblacion, const double *deltas_var){
    double sum_fitness = 0.0;

    for(int i = 0; i < n_pob; i++){
        // Se genera cada bit con la misma probabilidad de ser 0 o 1:

        unsigned int bits_recorridos = 0;
		for (int j = 0; j < n_pars; j++) {
			const unsigned int k = idx_pars[j];
			double pow_2 = 1.0;
			double cadena_val = 0.0;
			for (int b = 0; b < (int)min_vars[k]; b++, bits_recorridos++, pow_2*=2.0) {
				(poblacion + i)->cadena[bits_recorridos] = (HybTaus(0.0, 1.0) <= 0.5) ? 0 : 1;
				cadena_val += ((poblacion + i)->cadena[bits_recorridos] ? 1.0 : 0.0) * pow_2;
			}

			// Asignar el valor segun la cadena formada para el atributo 'k':
			(poblacion + i)->vars[k] = *(deltas_var + k) * cadena_val + lim_inf[k];
        }

        switch( fitness_elegido ){
            case ROC:
                (poblacion + i)->eval = fitnessROC( poblacion + i, resp );
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
void FILTROS::selecPob(INDIV* sel_grp, const INDIV* poblacion, double *fitness_acum, const double suma_fitness){
    // Calcular la funcion de distribucion de probabilidad F(x):
    *(fitness_acum + n_pob - 1) = (poblacion + n_pob - 1)->eval / suma_fitness;
    for( int i = n_pob-2; i >= 0; i--){
        *(fitness_acum + i) = *(fitness_acum + i + 1) + (poblacion + i)->eval / suma_fitness;
    }

    // Armar el grupo de seleccion:
    for( int i = 0; i < seleccion; i++){
        const double ruleta = HybTaus(0.0, 1.0);
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
void FILTROS::cruzaPob(INDIV* cruza, const INDIV* sel_grp, const unsigned int n_bits){

    // Generar dos nuevos individuos a partir de dos padres, el proceso se repite hasta completar la fraccion 'seleccion'.
    for(int i = 0; i < seleccion; i+=2){

        int padre_1 = i;
        int padre_2 = i + 1;
        int bits_ini = 0;
		
        // Realizar cortes hasta terminar con la secuencia de los individuos:
        while(bits_ini < n_bits){
            const unsigned int n_bits_cpy = (unsigned int)HybTaus(1.0, (double)(n_bits - bits_ini));
            memcpy( (cruza + i  )->cadena + bits_ini, (sel_grp + padre_1)->cadena + bits_ini, n_bits_cpy*sizeof(unsigned char));
            memcpy( (cruza + i+1)->cadena + bits_ini, (sel_grp + padre_2)->cadena + bits_ini, n_bits_cpy*sizeof(unsigned char));

            int padre_swap = padre_1;
            padre_1 = padre_2;
            padre_2 = padre_swap;

            bits_ini += n_bits_cpy;
        }

        // Realizar mutacion sobre los hijos generados con cierta probabilidad: 'prob_mutacion'
        // Mutar el hijo 1:
        if( HybTaus(0.0, 1.0) <= prob_mutacion ){
            // Seleccionar un gen a mutar:
            const unsigned int gen = (unsigned int)(HybTaus(0.0, (double)n_bits - 1e-8));
            (cruza + i)->cadena[gen] = 1 - (cruza + i)->cadena[gen];
        }

        // Mutar el hijo 2:
        if( HybTaus(0.0, 1.0) <= prob_mutacion ){
            // Seleccionar un gen a mutar:
            const unsigned int gen = (unsigned int)(HybTaus(0.0, (double)n_bits - 1e-8));
            (cruza + i + 1)->cadena[gen] = 1 - (cruza + i + 1)->cadena[gen];
        }
    }
}




/*
    Metodo:     generar_pob
    Funcion:   Genera una nueva poblacion en base al muestreo de cada atributo.
*/
double FILTROS::generarPob(INDIV *poblacion, const INDIV *cruza, const INDIV *sel_grp, const double *deltas_var){
    // Construir la nueva poblacion a partir de los hijos generados y la seleccion previa:
    double suma_fitness = 0.0;

    for( int i = 0; i < seleccion; i++){
        memcpy( poblacion + i, cruza + i, sizeof(INDIV) );

        unsigned int bits_recorridos = 0;
        for( int j = 0; j < n_pars; j++){
            const unsigned int k = idx_pars[j];
            double pow_2 = 1.0;
            double cadena_val = 0.0;
            for( int b = 0; b < (int)min_vars[k]; b++, bits_recorridos++, pow_2*=2.0){
                (poblacion + i)->cadena[bits_recorridos] = (HybTaus(0.0, 1.0) <= 0.5) ? 0 : 1;
                cadena_val += ((poblacion + i)->cadena[bits_recorridos] ? 1.0 : 0.0) * (double)pow_2;
            }
            // Asignar el valor segun la cadena formada para el atributo 'k':
            (poblacion + i)->vars[ k ] = *(deltas_var + k) * cadena_val + lim_inf[k];
			DEB_MSG(COLOR_GREEN "[" COLOR_BLUE << k << COLOR_GREEN "] cadena: " COLOR_BACK_WHITE COLOR_BLACK << cadena_val << "/" << (poblacion + i)->vars[k] << " :: " COLOR_BACK_BLACK COLOR_BLUE << "delta: " << *(deltas_var + k) << " :: lim_inf: " << lim_inf[k] << " :: lim_sup: " << lim_sup[k] << COLOR_NORMAL);
        }
        switch( fitness_elegido ){
        case ROC:
            (poblacion + i)->eval = fitnessROC( poblacion + i, resp );
            break;
        case CORCON:
            (poblacion + i)->eval = fitnessCorCon( poblacion + i, resp );
            break;
        }
        suma_fitness += (poblacion + i)->eval;
    }

    for( int i = seleccion; i < n_pob; i++){
        const unsigned int sel_pos = (unsigned int)(HybTaus(0.0, (double)seleccion - 1e-8));
        memcpy( poblacion + i, sel_grp + sel_pos, sizeof(INDIV) );
        suma_fitness += (poblacion + i)->eval;
    }

    return suma_fitness;
}




/*
    Metodo:     GA
    Funcion:    Utiliza el algoritmo genetico para encotnrar los parametros automaticamente.
*/
void FILTROS::GA(){
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
    semilla = ini_semilla(0);

    // Inicia el algoritmo GA
    int k = 1;
    bool procesar = true;

    //// Generar la primer poblacion:
    TIMERS;
    GETTIME_INI;
    double suma_fitness = generarPobInicial(poblacion, deltas_var);
    qsort((void*)poblacion, n_pob, sizeof(INDIV), compIndiv);

    double *fitness_acum = new double [n_pob];

    char mensaje_iter[] = "Best fit: X.XXXXXXXX, L: XX.XXX , T: XX.XXX, Sigma: XX.XXX, K: XXX.X :: Laste eval: X.XXXX, L: XX.XXX , T: XX.XXX, Sigma: XX.XXX, K: XXX.X :: en XXXXX.XXXX s.\n";
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

        GETTIME_FIN;
        sprintf( mensaje_iter, "%i\t%1.8f\t%1.3f\t%1.3f\t%2.3f\t%3.0f\t%5.4f\n", k, mi_elite->eval, mi_elite->vars[PAR_T], mi_elite->vars[PAR_L], mi_elite->vars[PAR_SIGMA], mi_elite->vars[PAR_K], DIFTIME);
        escribirLog( mensaje_iter );
        barraProgreso( k, max_iters);
        k++;

        //// Verificar condiciones de paro:
        procesar = (k < max_iters);


    }while(procesar);

    GETTIME_FIN;
    sprintf( mensaje_iter, "%i\t%1.8f\t%1.3f\t%1.3f\t%2.3f\t%3.0f\t%5.4f\n", k, mi_elite->eval, mi_elite->vars[PAR_T], mi_elite->vars[PAR_L], mi_elite->vars[PAR_SIGMA], mi_elite->vars[PAR_K], DIFTIME);
    escribirLog( mensaje_iter );
    barraProgreso( k, max_iters);

    delete [] cruza;
    delete [] sel_grp;
    delete [] fitness_acum;
    delete [] poblacion;
}




/*
    Metodo:     busquedaExhaustiva
    Funcion:    Realiza una busqueda exhaustiva en los limites de los parametros
*/
void FILTROS::busquedaExhaustiva(){
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
            max_bit_val[ p ] = (lim_sup[ p ] - lim_inf[ p ] + 1e-12) / min_vars[p];
            n_bits *= max_bit_val[ p ] + 1;
            n_pars ++;
        }
    }


    int idx = 0;


    char mensaje_iter[] = "Best fit: X.XXXXXXXX, L: XX.XXX , T: XX.XXX, Sigma: XX.XXX, K: XXX.X :: Laste eval: X.XXXX, L: XX.XXX , T: XX.XXX, Sigma: XX.XXX, K: XXX.X :: en XXXXX.XXXX s.\n";
    escribirLog( "iteration\tbest_fit\tT\tL\tSigma\tK\telapsed_time\n" );


    TIMERS;

    GETTIME_INI;
    for( int i_L = 0; i_L <= max_bit_val[PAR_L]; i_L++){
        test->vars[PAR_L] = (double)i_L * min_vars[PAR_L] + lim_inf[PAR_L];

        for( int i_T = 0; i_T <= max_bit_val[PAR_T]; i_T++){
            test->vars[PAR_T] = (double)i_T * min_vars[PAR_T] + lim_inf[PAR_T];

            for( int i_K = 0; i_K <= max_bit_val[PAR_K]; i_K++){
                test->vars[PAR_K] = (double)i_K * min_vars[PAR_K] + lim_inf[PAR_K];

                for( int i_Sigma = 0; i_Sigma <= max_bit_val[PAR_SIGMA]; i_Sigma++){
                    test->vars[PAR_SIGMA] = (double)i_Sigma * min_vars[PAR_SIGMA] + lim_inf[PAR_SIGMA];

                    switch( fitness_elegido ){
                        case ROC:
                            test->eval = fitnessROC( test, resp );
                            break;
                        case CORCON:
                            test->eval = fitnessCorCon( test, resp );
                            break;
                    }

                    if( (test->eval) > (mi_elite->eval) ){
                        memcpy( mi_elite, test, sizeof(INDIV));
                    }
                    idx++;

                    GETTIME_FIN;
                    sprintf( mensaje_iter, "%i\t%1.8f\t%1.3f\t%1.3f\t%2.3f\t%3.0f\t%5.4f\n", idx, mi_elite->eval, mi_elite->vars[PAR_T], mi_elite->vars[PAR_L], mi_elite->vars[PAR_SIGMA], mi_elite->vars[PAR_K], DIFTIME);
                    escribirLog( mensaje_iter );
                    barraProgreso( idx, n_bits);
                }
            }
        }
    }

    GETTIME_FIN;
    sprintf( mensaje_iter, "%i\t%1.8f\t%1.3f\t%1.3f\t%2.3f\t%3.0f\t%5.4f\n", idx, mi_elite->eval, mi_elite->vars[PAR_T], mi_elite->vars[PAR_L], mi_elite->vars[PAR_SIGMA], mi_elite->vars[PAR_K], DIFTIME);
    escribirLog( mensaje_iter );
    barraProgreso( idx, n_bits);

    delete test;
}
